function out = swe1d_run(cfg)
% Solve SWE1D problem with given config cfg
% Time integration, initialization, recording at times tsave

% Check required fields
required = {'N','xl','xr','bed','init','tend'};
for f=required
  if ~isfield(cfg,f{1}), error('swe1d_run: missing field %s', f{1}); end
end

% Setup problem P
P = swe1d_setup(cfg);

% Spatial grid and bed
N = P.N; xc = P.xc; Bc = P.Bc;
x = xc; % physical coords centers
B = Bc(3:end-2); % physical cell centers

% Initial condition given as init handle: [h0,u0] = init(x,Bc)
[h0, u0] = cfg.init(x,B);

% Initialize conserved state: U = [h; hu]
U = zeros(2,N);
U(1,:) = max(0,h0); % h>=0
U(2,:) = U(1,:) .* u0;

% Time stepping parameters
tend = cfg.tend;
if isfield(cfg,'tsave')
  tsave_user = cfg.tsave;
else
  tsave_user = tend;
end
if ~isvector(tsave_user)
  error('swe1d_run: tsave must be vector of output times');
end
tsave_user = tsave_user(:).'; % row vector sorted
if any(tsave_user < 0) || any(diff(tsave_user)<0)
  error('swe1d_run: tsave must be nonnegative ascending');
end

% Time integrator
rk = P.rk;
if ~any(strcmp(rk,{'euler','ssprk2','ssprk3','rk4'}))
  error('swe1d_run: unknown rk integrator');
end

cfl = P.cfl;
dx = P.dx; g = P.g;
scheme = P.scheme;
order = P.order;
bc = P.bc;

% Prepare outputs
nt = length(tsave_user);
Hout = cell(1,nt);
HUout = cell(1,nt);
Tout = tsave_user;

% Store initial state separately (as out.H0 / out.HU0 internally)
H0 = U(1,:);
HU0 = U(2,:);

% Monitor outputs after every step if requested
monflag = isfield(cfg,'monitor') && cfg.monitor;

if monflag
  % initialize monitor arrays with doubling strategy
  mcap = 1000; % initial capacity
  mon.t = zeros(1,mcap);
  mon.umax = zeros(1,mcap);
  mon.mass = zeros(1,mcap);
  moncount = 1;
  mon.t(1) = 0;
  wet = U(1,:) > 1e-10;
  uvel = zeros(size(U(1,:)));
  uvel(wet) = abs(U(2,wet)./U(1,wet));
  mon.umax(1) = max(uvel);
  mon.mass(1) = sum(U(1,:));
end

% Set maxsteps param with override
maxsteps = 2e6;
if isfield(cfg,'maxsteps') && ~isempty(cfg.maxsteps)
    maxsteps = cfg.maxsteps;
end

tcur = 0;
nsteps = 0;
nupdates = 0;
eps = 1e-12;

nout = 1; % Index into tsave_user, next output time to hit

% --- Immediately record outputs for all tsave <= 0 (or tcur), no step ---
while nout <= nt && tsave_user(nout) <= tcur + eps
    Hout{nout} = U(1,:);
    HUout{nout} = U(2,:);
    nout = nout + 1;
end

while tcur < tend - eps

  % Find next output time to hit exactly
  if nout <= nt
    tnxt = tsave_user(nout);
  else
    tnxt = tend; % no more output times
  end

  % Compute max wave speed
  hmax = max(U(1,:));
  uvel = zeros(size(U(1,:)));
  wet = U(1,:) > 1e-10;
  uvel(wet) = abs(U(2,wet)./U(1,wet));
  smax = max(uvel + sqrt(g*U(1,:)));

  if smax == 0
    dt_cfl = cfl * dx; % no flow, arbitrary dt
  else
    dt_cfl = cfl * dx / smax;
  end

  % Check blow-up guard on dt_cfl BEFORE clipping
  if dt_cfl < 1e-14
      error('swe1d_run:blowup','timestep too small at t=%g (step %d)', tcur, nsteps);
  end

  dt = dt_cfl;

  % Adjust dt to land on output time exactly, only if tnxt >= tcur
  if nout <= nt && (tcur + dt > tnxt) && (tnxt >= tcur)
    dt = tnxt - tcur;
  end

  % Never step with dt <= 0 (can happen if tnxt == tcur)
  if dt <= 0
    % Instead record output times equal or less than current time immediately
    while nout <= nt && tsave_user(nout) <= tcur + eps
      Hout{nout} = U(1,:);
      HUout{nout} = U(2,:);
      nout = nout + 1;
    end
    % If after recording outputs we are at or beyond tend, break loop
    if tcur >= tend - eps
      break
    end
    % Now continue to next loop iteration without stepping
    continue
  end

  % RK step
  Unew = rk_step(U, P, dt, rk);

  % Diagnostics check: NaN in h indicates failure
  if any(isnan(Unew(1,:)))
    error('swe1d_run: negative or NaN depth encountered, abort');
  end

  U = Unew;
  tcur = tcur + dt;
  nsteps = nsteps + 1;

  % Blow-up guard checks
  if ~all(isfinite(U(:))) || max(U(1,:)) > 1e4 || min(U(1,:)) < -1e-8
      error('swe1d_run:blowup','solution blew up at t=%g (step %d)', tcur, nsteps);
  end
  % Removed dt < 1e-14 check here (moved above on dt_cfl before clipping)

  if nsteps > maxsteps
      error('swe1d_run:maxsteps','exceeded maximum steps %d at t=%g', maxsteps, tcur);
  end

  % Count updates: nupdates = nsteps * stages * N
  stages = 1;
  switch rk
    case 'euler', stages=1;
    case 'ssprk2', stages=2;
    case 'ssprk3', stages=3;
    case 'rk4', stages=4;
  end
  nupdates = nsteps * stages * N;

  % Record monitor outputs after EVERY step, if requested
  if monflag
    moncount = moncount + 1;
    if moncount > length(mon.t)
      % expand monitor arrays if needed
      mon.t = [mon.t zeros(1,length(mon.t))];
      mon.umax = [mon.umax zeros(1,length(mon.umax))];
      mon.mass = [mon.mass zeros(1,length(mon.mass))];
    end
    mon.t(moncount) = tcur;
    wet = U(1,:) > 1e-10;
    uvel = zeros(size(U(1,:)));
    uvel(wet) = abs(U(2,wet)./U(1,wet));
    mon.umax(moncount) = max(uvel);
    mon.mass(moncount) = sum(U(1,:));
  end

  % Record outputs only if current time exactly matches output time (within eps)
  if nout <= nt && abs(tcur - tnxt) < eps
    Hout{nout} = U(1,:);
    HUout{nout} = U(2,:);
    nout = nout + 1;

    % Additionally, if multiple tsave times equal to tcur, record all of them now
    % This covers repeated save times equal to tcur
    while nout <= nt && tsave_user(nout) <= tcur + eps
      Hout{nout} = U(1,:);
      HUout{nout} = U(2,:);
      nout = nout + 1;
    end
  end

  if tcur >= tend - eps
    break
  end

end

% Trim monitor fields
if monflag
  mon.t = mon.t(1:moncount);
  mon.umax = mon.umax(1:moncount);
  mon.mass = mon.mass(1:moncount);
else
  mon = [];
end

% If any output slots remain unset because tsave times not hit exactly, fill with last U
for k = nout:nt
  Hout{k} = U(1,:);
  HUout{k} = U(2,:);
end

out.x = x;
out.B = B;
out.t = tsave_user; % exactly what user requested (no prepended zero)
out.H = Hout;
out.HU = HUout;
out.nsteps = nsteps;
out.nupdates = nupdates;
out.mon = mon;
out.P = P;
out.H0 = H0; % initial state stored separately
out.HU0 = HU0;

end

function Unew = rk_step(U, P, dt, rk)
% single Runge-Kutta time step using scheme in rk

switch rk
  case 'euler'
    dU = swe1d_rhs(U, P);
    Unew = U + dt*dU;

  case 'ssprk2'
    % Heun (2-stage SSP RK)
    dU1 = swe1d_rhs(U, P);
    U1 = U + dt*dU1;
    dU2 = swe1d_rhs(U1, P);
    Unew = 0.5*U + 0.5*(U1 + dt*dU2);

  case 'ssprk3'
    % Shu-Osher 3-stage
    dU1 = swe1d_rhs(U, P);
    U1 = U + dt*dU1;
    dU2 = swe1d_rhs(U1, P);
    U2 = 0.75*U + 0.25*(U1 + dt*dU2);
    dU3 = swe1d_rhs(U2, P);
    Unew = (1/3)*U + (2/3)*(U2 + dt*dU3);

  case 'rk4'
    k1 = swe1d_rhs(U, P);
    k2 = swe1d_rhs(U + dt/2*k1, P);
    k3 = swe1d_rhs(U + dt/2*k2, P);
    k4 = swe1d_rhs(U + dt*k3, P);
    Unew = U + dt/6*(k1 + 2*k2 + 2*k3 + k4);

  otherwise
    error('rk_step: unknown integrator');
end
end
