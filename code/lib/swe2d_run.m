function out = swe2d_run(cfg)
% Set defaults
if ~isfield(cfg,'g'), cfg.g=9.81; end
if ~isfield(cfg,'cfl')
  if strcmp(cfg.method,'plain'), cfg.cfl=0.45; else cfg.cfl=0.9; end
end

Nx = cfg.Nx; Ny = cfg.Ny;
xl = cfg.xl; xr = cfg.xr; yl = cfg.yl; yr = cfg.yr;
dx = (xr - xl)/Nx; dy = (yr - yl)/Ny;

x = xl + ((1:Nx)'-0.5)*dx;
y = yl + ((1:Ny)'-0.5)*dy;

[X,Y] = ndgrid(x,y);
B = cfg.bed2(X,Y);

init_data = cfg.init(X,Y,B);
h = max(0,init_data(:,:,1));
u = init_data(:,:,2);
v = init_data(:,:,3);
hu = h.*u;
hv = h.*v;

U = zeros(3,Nx,Ny);
U(1,:,:) = h;
U(2,:,:) = hu;
U(3,:,:) = hv;

tend = cfg.tend;
if ~isfield(cfg,'tsave') || isempty(cfg.tsave)
  tsave = tend;
else
  tsave = cfg.tsave;
end
tsave = tsave(:).';
tsave = tsave(tsave <= tend & tsave >= 0);
if isempty(tsave), tsave = 0; end

g = cfg.g;
bc = cfg.bc;
method = cfg.method;
cfl = cfg.cfl;

t = 0;
nsteps = 0;

out.x = x;
out.y = y;
out.B = B;
out.method = method;

out.H = {};
out.HU = {};
out.HV = {};
out.t = [];
out.umax = [];

next_save_idx = 1;
nextsave = tsave(next_save_idx);

% Record initial state if first tsave is 0
if abs(t - nextsave) < 1e-14
  out.H{end+1} = reshape(U(1,:,:),Nx,Ny);
  out.HU{end+1} = reshape(U(2,:,:),Nx,Ny);
  out.HV{end+1} = reshape(U(3,:,:),Nx,Ny);
  out.t(end+1) = t;
  % Compute umax at initial state
  hh = reshape(U(1,:,:),Nx,Ny);
  uu = reshape(U(2,:,:),Nx,Ny);
  vv = reshape(U(3,:,:),Nx,Ny);
  wet = hh>1e-10;
  if any(wet(:))
    umax0 = max([0; abs(uu(wet))./hh(wet); abs(vv(wet))./hh(wet)]);
  else
    umax0 = 0;
  end
  out.umax(end+1) = umax0;
  next_save_idx = next_save_idx + 1;
  if next_save_idx <= length(tsave)
    nextsave = tsave(next_save_idx);
  else
    nextsave = tend+1; % past tend
  end
end

max_steps = 2e5;

while t < tend && next_save_idx <= length(tsave)
  h = reshape(U(1,:,:),Nx,Ny);
  hu = reshape(U(2,:,:),Nx,Ny);
  hv = reshape(U(3,:,:),Nx,Ny);
  hpos = h > 1e-10;
  c = sqrt(g*h);
  u = zeros(Nx,Ny);
  v = zeros(Nx,Ny);
  u(hpos) = hu(hpos)./h(hpos);
  v(hpos) = hv(hpos)./h(hpos);
  umax_local = max([0; abs(u(hpos)); abs(v(hpos))]);
  dt1 = dx./(umax_local + max(c(:)) + 1e-16);
  dt2 = dy./(umax_local + max(c(:)) + 1e-16);
  dt = cfl*min(dt1,dt2);

  % Clip dt to hit nextsave exactly, if possible
  if t + dt > nextsave
    dt = nextsave - t;
  end
  if ~(isfinite(dt) && dt > 0)
    error('Non-finite dt encountered');
  end
  
  L = rhs_pair(U,B,dx,dy,g,bc,method,dt);

  Unew = U + dt*L;
  if strcmp(method,'ctu')
    hnew = reshape(Unew(1,:,:),Nx,Ny);
    Unew(1,:,:) = max(0,hnew);
  end

  U = Unew;
  t = t + dt;
  nsteps = nsteps + 1;
  if nsteps > max_steps
    error('Maximum number of steps (%d) exceeded', max_steps);
  end

  % Check for finite values
  if any(~isfinite(U(:)))
    error('Non-finite state values encountered at t = %g', t);
  end

  % After stepping, record states at all tsave that have been reached
  while next_save_idx <= length(tsave) && (abs(t - tsave(next_save_idx)) < 1e-12 || t > tsave(next_save_idx))
    out.H{end+1} = reshape(U(1,:,:),Nx,Ny);
    out.HU{end+1} = reshape(U(2,:,:),Nx,Ny);
    out.HV{end+1} = reshape(U(3,:,:),Nx,Ny);
    out.t(end+1) = tsave(next_save_idx);

    % Compute umax for current state using wet mask
    hh = reshape(U(1,:,:),Nx,Ny);
    uu = reshape(U(2,:,:),Nx,Ny);
    vv = reshape(U(3,:,:),Nx,Ny);
    wet = hh > 1e-10;
    if any(wet(:))
      umax = max([0; abs(uu(wet))./hh(wet); abs(vv(wet))./hh(wet)]);
    else
      umax = 0;
    end
    out.umax(end+1) = umax;

    next_save_idx = next_save_idx + 1;
    if next_save_idx <= length(tsave)
      nextsave = tsave(next_save_idx);
    else
      nextsave = tend + 1; % so loop ends
    end
  end
end

out.nsteps = nsteps;
end

function Lx = xsweep(U,B,dx,g,bc)
% Inputs: U 3 x Nx x Ny, B Nx x Ny, returns Lx 3 x Nx x Ny

[~,Nx,Ny] = size(U);
Up = zeros(3,Nx+2,Ny);
Bp = zeros(Nx+2,Ny);
Up(:,2:Nx+1,:) = U; 
Bp(2:Nx+1,:) = B;

switch bc
  case 'wall'
    Up(:,1,:) = Up(:,2,:);
    Up(2,1,:) = -Up(2,1,:);
    Up(:,Nx+2,:) = Up(:,Nx+1,:);
    Up(2,Nx+2,:) = -Up(2,Nx+2,:);
    Bp(1,:) = Bp(2,:);
    Bp(Nx+2,:) = Bp(Nx+1,:);
  case 'open'
    Up(:,1,:) = Up(:,2,:);
    Up(:,Nx+2,:) = Up(:,Nx+1,:);
    Bp(1,:) = Bp(2,:);
    Bp(Nx+2,:) = Bp(Nx+1,:);
  case 'periodic'
    Up(:,1,:) = Up(:,Nx+1,:);
    Up(:,Nx+2,:) = Up(:,2,:);
    Bp(1,:) = Bp(Nx+1,:);
    Bp(Nx+2,:) = Bp(2,:);
  otherwise
    error('Unknown bc');
end

h = reshape(Up(1,:,:), Nx+2, Ny);
hu = reshape(Up(2,:,:), Nx+2, Ny);
hv = reshape(Up(3,:,:), Nx+2, Ny);
h(h<=1e-10) = 0;
u = zeros(size(h));
v = zeros(size(h));
wet = h > 1e-10;
u(wet) = hu(wet)./h(wet);
v(wet) = hv(wet)./h(wet);

eta = h + Bp;

etaL = eta(1:end-1,:); etaR = eta(2:end,:);
uL = u(1:end-1,:); uR = u(2:end,:);
vL = v(1:end-1,:); vR = v(2:end,:);
BL = Bp(1:end-1,:);
BR = Bp(2:end,:);

Bstar = max(BL,BR);
hL0 = max(0, etaL - BL);
hR0 = max(0, etaR - BR);
hLs = max(0, etaL - Bstar);
hRs = max(0, etaR - Bstar);

F = hll2d(hLs(:).', uL(:).', vL(:).', hRs(:).', uR(:).', vR(:).', g);
F = reshape(F,3,Nx+1,Ny);
FL = F;
FR = F;
FL(2,:,:) = FL(2,:,:) + reshape(g/2*(hL0.^2 - hLs.^2),1,Nx+1,Ny);
FR(2,:,:) = FR(2,:,:) + reshape(g/2*(hR0.^2 - hRs.^2),1,Nx+1,Ny);

% Flux difference
Lx = -(FL(:,2:end,:) - FR(:,1:end-1,:))/dx;
end

function Ly = ysweep(U,B,dy,g,bc)
% Rotate and permute data to reuse xsweep for y-direction
U2 = permute(U,[1 3 2]);
U2 = U2([1 3 2],:,:);
B2 = B.';
Lt = xsweep(U2,B2,dy,g,bc);
Lt = Lt([1 3 2],:,:);
Ly = permute(Lt,[1 3 2]);
end

function F = hll2d(hL,uL,vL,hR,uR,vR,g)
tol = 1e-10;
n = numel(hL);

UL = [hL; hL.*uL; hL.*vL];
UR = [hR; hR.*uR; hR.*vR];

FL = [hL.*uL; hL.*uL.^2 + g*hL.^2/2; hL.*uL.*vL];
FR = [hR.*uR; hR.*uR.^2 + g*hR.^2/2; hR.*uR.*vR];

cL = sqrt(g.*hL);
cR = sqrt(g.*hR);

sqrt_hL = sqrt(max(hL,0));
sqrt_hR = sqrt(max(hR,0));
denom = sqrt_hL + sqrt_hR;
u_hat = zeros(1,n);
mask_dryL = hL <= tol;
mask_dryR = hR <= tol;
both_wet = ~mask_dryL & ~mask_dryR & (denom > 0);
u_hat(both_wet) = (sqrt_hL(both_wet).*uL(both_wet) + sqrt_hR(both_wet).*uR(both_wet)) ./ denom(both_wet);

c_hat = sqrt(g.*(hL+hR)/2);

sL = min(uL - cL, u_hat - c_hat);
sR = max(uR + cR, u_hat + c_hat);

% Adjust speed for dry states
sL(mask_dryL & ~mask_dryR) = uR(mask_dryL & ~mask_dryR) - 2*cR(mask_dryL & ~mask_dryR);
sR(mask_dryL & ~mask_dryR) = uR(mask_dryL & ~mask_dryR) + cR(mask_dryL & ~mask_dryR);
sL(mask_dryR & ~mask_dryL) = uL(mask_dryR & ~mask_dryL) - cL(mask_dryR & ~mask_dryL);
sR(mask_dryR & ~mask_dryL) = uL(mask_dryR & ~mask_dryL) + 2*cL(mask_dryR & ~mask_dryL);

sLm = min(sL,0);
sRp = max(sR,0);
den = sRp - sLm;
den(den == 0) = 1;

F = (FL .* sRp - FR .* sLm + (UR - UL) .* (sLm .* sRp)) ./ den;

both_dry = (hL <= tol) & (hR <= tol);
F(:,both_dry) = 0;
end

function L = rhs_pair(U,B,dx,dy,g,bc,method,dt)
if strcmp(method,'plain')
  L = xsweep(U,B,dx,g,bc) + ysweep(U,B,dy,g,bc);
elseif strcmp(method,'ctu')
  LyU = ysweep(U,B,dy,g,bc);
  Uy = U + dt/2*LyU;
  hUy = reshape(Uy(1,:,:),size(U,2),size(U,3));
  Uy(1,:,:) = max(0,hUy);

  LxUy = xsweep(Uy,B,dx,g,bc);

  LxU = xsweep(U,B,dx,g,bc);
  Ux = U + dt/2*LxU;
  hUx = reshape(Ux(1,:,:),size(U,2),size(U,3));
  Ux(1,:,:) = max(0,hUx);

  LyUx = ysweep(Ux,B,dy,g,bc);

  L = LxUy + LyUx;
else
  error('Unknown method');
end
end
