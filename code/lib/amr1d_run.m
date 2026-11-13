function out = amr1d_run(cfg)
% 1D adaptive mesh refinement for well-balanced WENO-3 scheme swe1d hr order 3
% Levels l=0..Lmax, refinement ratio=2, nested adaptive hierarchy.
%
% Main function + local functions below

% Defaults and config checks
if ~isfield(cfg,'N0'), error('amr1d_run: N0 required'); end
N0 = cfg.N0;
if ~isfield(cfg,'xl') || ~isfield(cfg,'xr'), error('amr1d_run: xl,xr required'); end
xl = cfg.xl; xr = cfg.xr;
if ~isfield(cfg,'bed'), error('amr1d_run: bed required'); end
bed = cfg.bed;
if ~isfield(cfg,'init'), error('amr1d_run: init required'); end
initf = cfg.init;
if ~isfield(cfg,'tend'), error('amr1d_run: tend required'); end
tend = cfg.tend;
if isfield(cfg,'tsave'), tsave = cfg.tsave; else tsave = tend; end
if ~isvector(tsave), error('amr1d_run: tsave must be vector'); end
tsave_user = tsave(:).';
drop_initial = tsave_user(1) > 0;
if drop_initial
  tsave = [0 tsave_user];
else
  tsave = tsave_user;
end
if any(tsave<0) || any(diff(tsave)<0), error('amr1d_run: tsave must be ascending nonneg'); end

Lmax = 2;
if isfield(cfg,'Lmax'), Lmax = cfg.Lmax; end
nregrid = 5;
if isfield(cfg,'nregrid'), nregrid = cfg.nregrid; end
thr_rel = 2e-3;
if isfield(cfg,'thr'), thr_rel = cfg.thr; end
cfl = 0.4;
if isfield(cfg,'cfl'), cfl = cfg.cfl; end
g = 9.81;
if isfield(cfg,'g'), g = cfg.g; end
bc = 'wall';
if isfield(cfg,'bc'), bc = cfg.bc; end

mode = 'adaptive';
if isfield(cfg,'mode'), mode = cfg.mode; end

fixedflagfunc = [];
if isfield(cfg,'fixedflag'), fixedflagfunc = cfg.fixedflag; end

% No refined cells allowed within 4 level-0 cells of domain boundaries
no_refine_zone_l = 4;
no_refine_zone_r = N0 - 4 + 1; % index from 1

% Setup grids and geometry per level
U = cell(1,Lmax+1);
exists = cell(1,Lmax+1);
covered = cell(1,Lmax+1);
active = cell(1,Lmax+1);
dx = zeros(1,Lmax+1);
Nlev = zeros(1,Lmax+1);
xc = cell(1,Lmax+1);
Bc = cell(1,Lmax+1);
Bm = cell(1,Lmax+1);
Bp = cell(1,Lmax+1);

Nlev(1) = N0;
dx(1) = (xr - xl) / N0;
xc{1} = xl + ((1:N0) - 0.5)*dx(1);

% Compute grid and Bc for all levels by averaging max level cells
Nlev(Lmax+1) = N0 * 2^Lmax;
dx(Lmax+1) = dx(1)/2^Lmax;
xc{Lmax+1} = xl + ((1:Nlev(Lmax+1)) - 0.5)*dx(Lmax+1);

% Bed cell centers at finest level: average of bed.f at fine cell centers
Bc{Lmax+1} = average_bed(bed, xc{Lmax+1});

% For coarser levels average from children
for l = Lmax:-1:1
  Nlev(l) = N0*2^(l-1);
  dx(l) = dx(1)/2^(l-1);
  xc{l} = xl + ((1:Nlev(l)) - 0.5)*dx(l);
  % Average Bc from children
  Bc{l} = 0.5*(Bc{l+1}(1:2:end) + Bc{l+1}(2:2:end));
end

% Compute face beds for all levels (using analytic function at faces)
for l = 1:Lmax+1
  Nl = Nlev(l);
  dxl = dx(l);
  xf = xl + (0:Nl)*dxl;
  % Bm(k) = bed.f(xf(k),+1), Bp(k) = bed.f(xf(k),-1), length Nl+1
  Bm{l} = bed.f(xf, +1);
  Bp{l} = bed.f(xf, -1);
  % Bm,L and Bp,R aligned with faces
end

% --- build padded bed arrays Bcp, Mpad, Ppad for all levels ---
% Bcp: cell centered padded beds length Nl+4 (2 ghost cells each side)
% Bmp: left face bed padded length Nl+4 at faces k=-1..Nl+2 (for internal level_rhs)
% Bpp: right face bed padded length Nl+4 at faces k=-1..Nl+2

Bcp = cell(1,Lmax+1);
Bmp = cell(1,Lmax+1);
Bpp = cell(1,Lmax+1);

for l = 1:Lmax+1
  Nl = Nlev(l);
  dxl = dx(l);
  xl_l = xl;
  
  % Padded cell indices k = -1..Nl+2 (total Nl+4)
  k_padded = -1:(Nl+2);
  
  % Left and right faces of these padded cells (k=-1..Nl+2)
  xfL = xl_l + (k_padded - 1)*dxl;
  xfR = xl_l + k_padded*dxl;
  
  % Evaluate face bed arrays padded length Nl+4
  Bmp{l} = bed.f(xfL, +1);
  Bpp{l} = bed.f(xfR, -1);
  
  % For padded cell centers: average of finest bed at children for cell centers
  Bc_pad = zeros(1, Nl + 4);
  
  % Number of finest cells per padded cell at level l
  refinement_rate = 2^(Lmax - (l - 1));
  
  % Finest cell centers and bed (finest grid level Lmax+1)
  xcfine = xc{Lmax+1};
  Bfine = Bc{Lmax+1};
  
  % For each padded cell k, find finest cells inside and average
  for idx = 1:length(k_padded)
    kcell = k_padded(idx);
    % Cell center of level l padded cell k
    xc_cell = xl_l + (kcell - 0.5)*dxl;
    
    % Finest cell centers inside this coarse cell:
    % Finest cell indices covering cell k: from
    % leftmost finest cell with center > (cell left face)
    % to rightmost finest cell with center <= (cell right face)
    leftf = xl_l + (kcell - 1)*dxl;
    rightf = xl_l + kcell*dxl;
    
    % Find finest cells inside [leftf, rightf)
    % Indices of fine cells inside this cell:
    fine_inds = find(xcfine >= leftf & xcfine < rightf);
    if isempty(fine_inds)
      % Possible at domain boundaries: pick closest finest cell index (nearest)
      [~, minind] = min(abs(xcfine - xc_cell));
      fine_inds = minind;
    end
    
    % Average the finest bed values over fine_inds
    Bc_pad(idx) = mean(Bfine(fine_inds));
  end
  
  Bcp{l} = Bc_pad;
  
  % Interior cells slice indices relative to padding: 3 : Nl+2
  % Check that padded interior cell-centered bed equals unpadded Bc to 1e-15
  diff_check = max(abs(Bcp{l}(3:Nl+2) - Bc{l}));
  assert(diff_check < 1e-15, ...
    sprintf('Padded Bc interior mismatch at level %d max diff=%g', l, diff_check));
end

% Initialize existence mask, level 0 all true, others empty (no refined cells)
exists{1} = true(1,Nlev(1));
for l = 2:Lmax+1
  exists{l} = false(1,Nlev(l));
end

% Level 0 initial state on finest grid first
[hi,ui] = initf(xc{Lmax+1},Bc{Lmax+1});
h_fine = max(0, hi);
hu_fine = h_fine .* ui;
U{Lmax+1} = zeros(2,Nlev(Lmax+1));
U{Lmax+1}(1,:) = h_fine;
U{Lmax+1}(2,:) = hu_fine;

% Average down to coarser levels
for l = Lmax:-1:1
  U{l} = zeros(2,Nlev(l));
  U{l}(1,:) = 0.5*(U{l+1}(1,1:2:end) + U{l+1}(1,2:2:end));
  U{l}(2,:) = 0.5*(U{l+1}(2,1:2:end) + U{l+1}(2,2:2:end));
  exists{l}(:) = true; % Initially coarse all exists true
end

% Composite marking: initially no refined active cells except level0 fully active
for l = 1:Lmax
  covered{l} = false(1,Nlev(l));
end
covered{Lmax+1} = false;
for l = 1:Lmax
  covered{l} = exists{l+1}(1:2:end) & exists{l+1}(2:2:end);
end

for l = 1:Lmax+1
  active{l} = exists{l};
  if l <= Lmax
    active{l}(covered{l}) = false;
  end
end

% Determine initial refinement flags and regrid at t=0 (except if fixed 'fixedflag')
if strcmp(mode,'adaptive') || ~isempty(fixedflagfunc)
  % Compute threshold from LEVEL-1 (coarsest) data:
  eta0 = U{1}(1,:) + Bc{1};
  d_eta = max(eta0) - min(eta0);
  if d_eta < 1e-14
    thr = max(thr_rel, 1e-3);
  else
    thr = thr_rel * d_eta;
  end

  % Use fixedflagfunc on level0 xc if provided (logical mask)
  fixed_flag0 = false(1,Nlev(1));
  if ~isempty(fixedflagfunc)
    fixed_flag0 = fixedflagfunc(xc{1});
  end

  % Keep a copy of initial U before first regrid for sharp initial data overwrite
  U_init = U;

  % Call new regrid function (no refine_flag call)
  [U, exists] = regrid_hierarchy(U, exists, Bc, thr, fixed_flag0, no_refine_zone_l, no_refine_zone_r);

  % Overwrite cells that exist at levels l>=2 with exact finest-grid initial data averaged down
  for l = 2:Lmax+1
    exists_l = exists{l};
    if any(exists_l)
      U{l}(:,exists_l) = U_init{l}(:,exists_l);
    end
  end

  % Update covered and active masks after regrid
  for l = 1:Lmax
    covered{l} = exists{l+1}(1:2:end) & exists{l+1}(2:2:end);
  end
  for l = 1:Lmax+1
    active{l} = exists{l};
    if l <= Lmax, active{l}(covered{l}) = false; end
  end

  % Re-average down fine->coarse after updating covered/active masks
  U = average_down(U, covered, Lmax);
end

% Time stepping parameters
tcur = 0;
nsteps = 0;
nupdates = 0;
tout = 1; % next tsave index

ntsave = length(tsave);
Hout = cell(1,ntsave);
HUout = cell(1,ntsave);
levelmapout = cell(1,ntsave);

% Save initial output (prolong composed solution to finest)
[Hf, HUf, levmap] = prolong_composite(U, exists);
Hout{1} = Hf;
HUout{1} = HUf;
levelmapout{1} = levmap;
tout = 2;

% SSP-RK3 constants removed (not used anymore)

max_step = 1e5; % safety limit

ncells_hist = zeros(1,0); t_hist = zeros(1,0);

while tcur < tend && nsteps < max_step
  % Compute common dt from active cells on all levels
  dx_finest = dx(Lmax+1);
  [dt, maxwave] = compute_dt(U, exists, active, dx, g, cfl, Lmax);

  % Adjust dt to hit next output time exactly
  if tout <= ntsave && (tcur + dt > tsave(tout))
    dt = tsave(tout) - tcur;
  end

  % SSP-RK3 step on the hierarchy
  [U, nact] = ssprk3_hier(U, dt, exists, active, covered, Bcp, Bmp, Bpp, dx, g, bc, Lmax);
  nupdates = nupdates + 3 * nact;

  % After step update time and steps
  tcur = tcur + dt;
  nsteps = nsteps + 1;
  ncells_hist(end+1) = nact; t_hist(end+1) = tcur;

  if tout <= ntsave && abs(tcur - tsave(tout)) < 1e-12
    [Hf, HUf, levmap] = prolong_composite(U, exists);
    Hout{tout} = Hf;
    HUout{tout} = HUf;
    levelmapout{tout} = levmap;
    tout = tout + 1;
  end

  % Regrid every nregrid steps and if adaptive
  if strcmp(mode,'adaptive') && mod(nsteps,nregrid) == 0
    % Compute threshold from LEVEL-1 (coarsest) data:
    eta0 = U{1}(1,:) + Bc{1};
    d_eta = max(eta0) - min(eta0);
    if d_eta < 1e-14
      thr = max(thr_rel, 1e-3);
    else
      thr = thr_rel * d_eta;
    end

    % Override flags with fixedflag if provided
    fixed_flag0 = false(1,Nlev(1));
    if ~isempty(fixedflagfunc)
      fixed_flag0 = fixedflagfunc(xc{1});
    end

    % Call new regrid function (no refine_flag call)
    [U, exists] = regrid_hierarchy(U, exists, Bc, thr, fixed_flag0, no_refine_zone_l, no_refine_zone_r);

    % Update masks covered, active
    for l = 1:Lmax
      covered{l} = exists{l+1}(1:2:end) & exists{l+1}(2:2:end);
    end
    for l = 1:Lmax+1
      active{l} = exists{l};
      if l <= Lmax, active{l}(covered{l}) = false; end
    end
  end

  if tcur >= tend
    break
  end

end

% Final outputs: prolong solution and level map for all tsave filled
for k = 1:ntsave
  if isempty(Hout{k})
    [Hf, HUf, levmap] = prolong_composite(U, exists);
    Hout{k} = Hf;
    HUout{k} = HUf;
    levelmapout{k} = levmap;
  end
end

% Remove initial output if drop_initial true
if drop_initial
  Hout(1) = [];
  HUout(1) = [];
  levelmapout(1) = [];
  out.t = tsave_user;
else
  out.t = tsave_user;
end

% Output struct
out.x = xc{Lmax+1};
out.H = Hout;
out.HU = HUout;
out.level = levelmapout;
out.ncells_hist = ncells_hist; out.t_hist = t_hist;
out.nupdates = nupdates;
out.nsteps = nsteps;

end

function Bc_l = average_bed(bed, xcl)
% cell centers Bc at level points = average of finest bed at children (if finest)
% For finest level, evaluate directly
Bc_l = bed.f(xcl,0);
end

function flags = refine_flag(Ul, Bc_l, thr, fixflag)
% Flag cells on level l to refine if
% Either fixed flag or gradient of eta or hu above threshold thr
% Ul 2xN: [h; hu], Bc_l cell bed centers, thr scalar
% fixflag logical 1xN forcing refinement

h = Ul(1,:); hu = Ul(2,:);
eta = h + Bc_l;
N = length(h);
flags = false(1,N);
% Compute internal diffs for eta: centred 1-cell
if N>2
  eta_grad = zeros(1,N);
  hu_grad = zeros(1,N);
  eta_grad(2:end-1) = abs((eta(3:end)-eta(1:end-2))/2);
  hu_grad(2:end-1) = abs((hu(3:end)-hu(1:end-2))/2);
else
  eta_grad = zeros(size(eta));
  hu_grad = zeros(size(hu));
end

flags(2:end-1) = (eta_grad(2:end-1) > thr) | (hu_grad(2:end-1) > thr);
% Set fixed flag forced if provided
if ~isempty(fixflag)
  flags = flags | fixflag;
end
% Dilate flags by 3 cells (logical dilation)
for k=1:3
  flags = flags | [false flags(1:end-1)];
  flags = flags | [flags(2:end) false];
end
% Also disallow refinement in first 4 and last 4 on level 0 by caller code

end

function [U_new, exists_new] = regrid_hierarchy(U_old, exists_old, Bc, thr, fixed_flag0, no_refine_zone_l, no_refine_zone_r)
% Regrid hierarchy with correct flags, prolongation from current coarse solution,
% and threshold refinement properly applied.
%
% Inputs:
%   U_old        - cell array {1..Lmax+1} of 2xNl [h;hu], old solution data
%   exists_old   - cell array {1..Lmax+1} logical mask of existing cells
%   Bc           - cell array {1..Lmax+1} cell bed arrays 1xNl (unpadded)
%   thr          - scalar refinement threshold on one-cell difference at level 1
%   fixed_flag0  - logical 1xN0 refinement mask forced on level 1 (may be all false)
%   no_refine_zone_l - integer: left no-refine zone limit (flags(1:no_refine_zone_l)=false)
%   no_refine_zone_r - integer: right no-refine zone start index (flags(no_refine_zone_r:end)=false)
%
% Outputs:
%   U_new        - cell array {1..Lmax+1} 2xNl new solution data after refinement
%   exists_new   - cell array {1..Lmax+1} logical mask of existing cells after refinement

Lmax = length(U_old) - 1;
N0 = length(exists_old{1});

U_new = cell(1,Lmax+1);
exists_new = cell(1,Lmax+1);

% Initialize level 1 (coarsest) solution copy and exists (all true)
U_new{1} = U_old{1};
exists_new{1} = true(1,N0);

for l = 1:Lmax
    Nl = length(exists_new{l});
    Nlf = 2 * Nl;
    eta = zeros(1,Nl);
    hu = zeros(1,Nl);
    eta(:) = U_new{l}(1,:) + Bc{l}; % eta = h+Bc at level l current solution
    hu(:) = U_new{l}(2,:);
    
    % Compute centered one-cell differences only where both neighbours exist
    de = zeros(1,Nl);
    dq = zeros(1,Nl);
    if Nl>2
        de(2:end-1) = abs(eta(3:end) - eta(1:end-2)) / 2;
        dq(2:end-1) = abs(hu(3:end) - hu(1:end-2)) / 2;
    end
    
    % Validity mask for differences (centered diff): neighbors must exist
    ok = exists_new{l};
    if Nl>2
        ok(2:end-1) = ok(2:end-1) & exists_new{l}(1:end-2) & exists_new{l}(3:end);
    end
    ok([1 end]) = false;
    
    thr_l = thr / 2^(l-1);
    
    % Flags for refinement on level l from threshold check and existence
    flags = ok & ((de > thr_l) | (dq > thr_l));
    
    % For level 1 apply fixed_flag0 and no-refinement zones
    if l == 1
        flags = flags | fixed_flag0;
        flags(1:no_refine_zone_l) = false;
        flags(no_refine_zone_r:end) = false;
    end
    
    % Dilate flags by 3 cells, repeat 3 times
    for i=1:3
        flags = flags | [false flags(1:end-1)] | [flags(2:end) false];
    end
    
    % Keep only existing cells
    flags = flags & exists_new{l};
    
    % Reapply no-refine zones for level 1 after dilation
    if l == 1
        flags(1:no_refine_zone_l) = false;
        flags(no_refine_zone_r:end) = false;
    end
    
    % Proper nested refinement: parents and neighbors exist
    allowed = flags;
    if Nl>2
        allowed(2:end-1) = allowed(2:end-1) & exists_new{l}(1:end-2) & exists_new{l}(2:end-1) & exists_new{l}(3:end);
    end
    allowed([1 end]) = false;
    
    % Set existence of children at fine level l+1 from parents in allowed
    exists_new{l+1} = false(1, Nlf);
    idxp = find(allowed);
    exists_new{l+1}(2*idxp - 1) = true;
    exists_new{l+1}(2*idxp) = true;
    
    % Initialize U_new at level l+1
    U_new{l+1} = zeros(2, Nlf);
    
    % Copy existing fine cells from old data if both old and new exist at fine level
    fine_exists_old = exists_old{l+1};
    fine_exists_new = exists_new{l+1};
    common_exists = fine_exists_old & fine_exists_new;
    U_new{l+1}(:,common_exists) = U_old{l+1}(:,common_exists);
    
    % Find new fine cells that need prolongation from coarse level l
    new_fine = fine_exists_new & ~fine_exists_old;
    new_indices = find(new_fine);
    
    % Prepare for prolongation
    % For local minmod function
    mm = @(a,b) ( (a .* b <= 0) .* 0 ) + ( (a .* b > 0) .* ( abs(a) < abs(b)) .* a ) + ...
                 ( (a .* b > 0) .* ( abs(b) <= abs(a)) .* b );
             
    for i = new_indices
        p = ceil(i/2); % parent coarse cell index
        sgn = -1 + 2*((mod(i,2)==0)); % -1 if odd i(left child), +1 if even i(right child)
        
        % Compute limited slope for eta
        if p > 1 && p < Nl && exists_new{l}(p-1) && exists_new{l}(p+1)
            etaL = U_new{l}(1,p-1) + Bc{l}(p-1);
            etaC = U_new{l}(1,p)   + Bc{l}(p);
            etaR = U_new{l}(1,p+1) + Bc{l}(p+1);
            slope_eta = mm(etaR - etaC, etaC - etaL);
            
            huL = U_new{l}(2,p-1);
            huC = U_new{l}(2,p);
            huR = U_new{l}(2,p+1);
            slope_hu = mm(huR - huC, huC - huL);
        else
            % At edges or missing neighbors slope = 0
            etaC = U_new{l}(1,p) + Bc{l}(p);
            huC = U_new{l}(2,p);
            slope_eta = 0;
            slope_hu = 0;
        end
        
        eta_p = U_new{l}(1,p) + Bc{l}(p);
        hu_p = U_new{l}(2,p);
        
        % Fine eta and hu: linear limited reconstruction at fine cell center
        eta_f = eta_p + sgn * slope_eta * 0.25;
        hu_f  = hu_p  + sgn * slope_hu  * 0.25;
        
        % Use fine bed
        h_f = max(0, eta_f - Bc{l+1}(i));
        
        % Set fine cell data (h, hu)
        U_new{l+1}(:,i) = [h_f; hu_f];
    end
end

% Average down covered coarse cells to ensure consistency
for l = Lmax:-1:1
    children_exist = exists_new{l+1};
    covered = children_exist(1:2:end) & children_exist(2:2:end);
    cidx = find(covered);
    if ~isempty(cidx)
        U_new{l}(:,cidx) = 0.5*(U_new{l+1}(:,2*cidx-1) + U_new{l+1}(:,2*cidx));
    end
end

% Zero cells that do not exist, all levels
for l = 1:Lmax+1
    U_new{l}(:, ~exists_new{l}) = 0;
end

end

function [dt, smax] = compute_dt(U, exists, active, dx, g, cfl, Lmax)
% Compute stable time step dt with cfl on active cells of all levels on finest dx

% Active cells are only progressed, but dt based on finest dx active cells
smax = 0;
for l=1:Lmax+1
  active_l = active{l};
  if ~any(active_l), continue; end
  U_l = U{l}(:,active_l);
  if isempty(U_l), continue; end
  h = U_l(1,:);
  hu = U_l(2,:);
  tol=1e-10;
  u = zeros(size(h));
  wet = (h > tol);
  u(wet) = hu(wet)./h(wet);
  smax_l = max(abs(u) + sqrt(g*h));
  if smax_l > smax, smax = smax_l; end
end
if smax==0, smax=1; end
dt = cfl*dx(Lmax+1)/smax;
end

function dU_hier = amr_rhs(U, exists, active, covered, Bc, Bm, Bp, dx, g, bc, Lmax)
% Compute RHS dU on all levels padded with 2 ghosts, hydrostatic reconstruction hr order 3,
% ghost prolongation from coarse level for ghosts inside domain, lower level bc for coarsest,
% reflux corrections. Recomputes source only once in level_rhs and uses full arrays.
%
% NOTE: Bc, Bm, Bp are PADDED (2 ghosts per side): padded index = physical index + 2.

dU_hier = cell(1,Lmax+1);
FL_hier = cell(1,Lmax+1);
FR_hier = cell(1,Lmax+1);
S_hier = cell(1,Lmax+1);

% Pad arrays, fill ghosts with prolongation on refined boundaries or physical bc at level0
Upad = cell(1,Lmax+1);

for l=1:Lmax+1
  Nl = length(exists{l});
  Ucurr = U{l};
  existmask = exists{l};
  activemask = active{l};
  % Construct padded arrays with 2 ghosts each side for existing cells only
  if l==1
    % Use physical BC for level0 ghosts
    Uvalid = Ucurr(:,existmask);
    % Pad with bc
    Upad{l} = swe1d_pad(Uvalid, bc);
  else
    % Internal BC: prolong from coarse level below at boundaries of refined patch
    Uvalid = Ucurr(:,existmask);
    % Prepare padded array with 2 ghosts per side on the refined patch domain only
    % ghosts at boundary of exists true patches to be set by prolongation

    % Reconstruct full array with all exists false cells set zero
    Ufull = zeros(2,Nl);
    Ufull(:,existmask) = Ucurr(:,existmask);

    % Pad with zeros first: 2 ghost cells outside internal domain ok
    Upzero = zeros(2,Nl+4);
    Upzero(:,3:end-2) = Ufull;

    % Find boundaries of refined region (exists true adjacent to false)
    existpad = [false false existmask false false];
    i1 = find(diff(existpad)==1) - 1;    % first physical cell of each patch
    i2 = find(diff(existpad)==-1) - 2;   % last physical cell of each patch

    % For each refined patch in level l, prolong from coarse level l-1
    % Prolong level l-1 states at fine boundary ghosts
    % Note: assume contiguous single patch; multiple patches rare
    if isempty(i1) || isempty(i2)
      Upad{l} = Upzero; % no refined cells, ghost padding irrelevant
    else
      % Prolong ghosts on left boundary if near left patch boundary
      for p=1:length(i1)
        leftcell = i1(p);
        rightcell = i2(p);
        % Prolong left two ghosts if exist
        for ig=1:2
          ighost = leftcell + 2 - ig;
          if ighost>=1 && ighost<=Nl+4
            % fine ghost at cell ighost padded array index (relative to level l fine cells)
            % corresponding coarse cell at l-1:
            % fine level l cell index = ighost-2 corresponds to coarse cell ceil((ighost-2)/2)
            i_coarse = ceil((ighost-2)/2);
            if i_coarse >= 1 && i_coarse <= length(U{l-1})
              % Prolong coarse u, eta to fine ghost by linear slope limiting minmod
              % Coarse states at (eta, hu): eta = h + Bc
              hc = U{l-1}(1,i_coarse);
              huc = U{l-1}(2,i_coarse);
              etac = hc + Bc{l-1}(i_coarse+2);
              % Compute coarse slopes
              if i_coarse==1
                etaL = Bc{l-1}(i_coarse+2);
                etaR = Bc{l-1}(i_coarse+3);
                slope_eta = 0;
                slope_hu = 0;
              elseif i_coarse==length(U{l-1}(1,:))
                slope_eta=0;
                slope_hu=0;
              else
                etaL = U{l-1}(1,i_coarse-1) + Bc{l-1}(i_coarse-1+2);
                etaC = etac;
                etaR = U{l-1}(1,i_coarse+1) + Bc{l-1}(i_coarse+1+2);
                slope_eta = minmod((etaR - etaC),(etaC - etaL));
                huL = U{l-1}(2,i_coarse-1);
                huC = huc;
                huR = U{l-1}(2,i_coarse+1);
                slope_hu = minmod((huR - huC),(huC - huL));
              end

              % Prolong fine child's Bc, Bm, Bp
              % Estimate fine cell center and bed at ghost
              if ighost >= 3 && ighost <= Nl+2
                Bc_fine_ghost = Bc{l}(ighost);
                % eta at fine ghost children = coarse etac +- slope_eta/4 (corrected)
                eta_f = etac + (2*((ighost-2)-2*i_coarse)+1)*slope_eta/4;
                hu_f = huc + (2*((ighost-2)-2*i_coarse)+1)*slope_hu/4;
                h_f = max(0, eta_f - Bc_fine_ghost);
                Upzero(:,ighost) = [h_f; hu_f];
              end
            end
          end
        end
      end
      % Repeat same for right boundary patches for 2 right ghosts
      for p=1:length(i1)
        leftcell = i1(p);
        rightcell = i2(p);
        for ig=1:2
          ighost = rightcell + 2 + ig;
          if ighost >= 1 && ighost <= Nl+4
            i_coarse = ceil((ighost-2)/2);
            if i_coarse >= 1 && i_coarse <= length(U{l-1})
              hc = U{l-1}(1,i_coarse);
              huc = U{l-1}(2,i_coarse);
              etac = hc + Bc{l-1}(i_coarse+2);
              if i_coarse==1
                slope_eta = 0;
                slope_hu = 0;
              elseif i_coarse==length(U{l-1}(1,:))
                slope_eta=0; slope_hu=0;
              else
                etaL = U{l-1}(1,i_coarse-1) + Bc{l-1}(i_coarse-1+2);
                etaC = etac;
                etaR = U{l-1}(1,i_coarse+1) + Bc{l-1}(i_coarse+1+2);
                slope_eta = minmod((etaR - etaC),(etaC - etaL));
                huL = U{l-1}(2,i_coarse-1);
                huC = huc;
                huR = U{l-1}(2,i_coarse+1);
                slope_hu = minmod((huR - huC),(huC - huL));
              end
              if ighost >= 3 && ighost <= Nl+2
                Bc_fine_ghost = Bc{l}(ighost);
                eta_f = etac + (2*((ighost-2)-2*i_coarse)+1)*slope_eta/4;
                hu_f = huc + (2*((ighost-2)-2*i_coarse)+1)*slope_hu/4;
                h_f = max(0, eta_f - Bc_fine_ghost);
                Upzero(:,ighost) = [h_f; hu_f];
              end
            end
          end
        end
      end
      % For remaining ghosts not filled, copy nearest level l cells directly
      Upzero(:,1) = Upzero(:,3); 
      Upzero(:,2) = Upzero(:,3);
      Upzero(:,end) = Upzero(:,end-2); 
      Upzero(:,end-1) = Upzero(:,end-2);

      Upad{l} = Upzero;
    end
  end
end

% Compute fluxes and RHS for all levels (order 3, scheme 'hr')
for l=1:Lmax+1
  Up = Upad{l};
  Ufull = U{l};
  ex = exists{l};
  % Extract beds and dx
  Bc_l = Bc{l};
  Bm_l = Bm{l};
  Bp_l = Bp{l};
  dxl = dx(l);

  % Compute rhs using local level_rhs (copies swe1d_rhs hr 3 logic)
  % Modified level_rhs to return also source S
  [dUl, FLl, FRl, Sl] = level_rhs(Up, Bc_l, Bm_l, Bp_l, dxl, g);

  S_hier{l} = Sl;
  FL_hier{l} = FLl;
  FR_hier{l} = FRl;
end

% Reflux corrections: for l=0..Lmax-1, replace coarse face flux at active cell interfaces
for l=1:Lmax
  Nl = length(exists{l});
  % Coarse level fluxes FL,FR are size 2 x (Nl+1)
  FLl = FL_hier{l};
  FRl = FR_hier{l};
  FLf = FL_hier{l+1};
  FRf = FR_hier{l+1};

  ex_coarse = exists{l};
  active_coarse = active{l};
  ex_fine = exists{l+1};

  ia = find(active_coarse);
  % Vectorized reflux corrections
  % Right face: replace FLl(:,i+1) if i<Nl and covered{l}(i+1)
  ir = ia(ia < Nl & covered{l}(min(ia+1,Nl)));
  if ~isempty(ir)
    FLl(:,ir+1) = FLf(:,2*ir+1);
  end
  % Left face: replace FRl(:,i) if i>1 and covered{l}(i-1)
  il = ia(ia > 1 & covered{l}(max(ia-1,1)));
  if ~isempty(il)
    FRl(:,il) = FRf(:,2*il-1);
  end

  % Store corrected
  FL_hier{l} = FLl;
  FR_hier{l} = FRl;
end

% Compute final RHS with corrected fluxes and stored source terms S_hier
for l=1:Lmax+1
  Nl = length(exists{l});
  ex = exists{l};

  FLl = FL_hier{l};
  FRl = FR_hier{l};
  
  dU = zeros(2,Nl);
  % Use full length arrays for every level; then zero out non-active cells
  dU = -(FLl(:,2:Nl+1) - FRl(:,1:Nl))/dx(l) + S_hier{l};

  % Zero out cells not active
  dU(:,~active{l}) = 0;

  dU_hier{l} = dU;
end

% Note: remove average down step from here (caller must do average_down after RK stage)

end

function Uout = combine_levels(Uold, Ustage, dU, dt, a, b)
% RK convex combination of U with update: Unew = a*Uold + b*(Ustage + dt*dU)
% Only active cells updated, other cells copied unchanged
Lmaxp1 = length(Uold);
Uout = cell(1,Lmaxp1);
for l=1:Lmaxp1
  Uo = Uold{l};
  Us = Ustage{l};
  dUl = dU{l};
  Uout{l} = a*Uo + b*(Us + dt*dUl);
end
end

function mm = minmod(a,b)
% minmod limiter for vectors
mm = zeros(size(a));
idx = (a.*b)>0;
mm(idx) = sign(a(idx)) .* min(abs(a(idx)), abs(b(idx)));
end

function [dU, FL, FR, S] = level_rhs(Up, Bc, Bm, Bp, dx, g)
% Local level RHS function, copy of swe1d_rhs hr order 3
% Inputs:
%  Up: 2 x (N+4) array padded with 2 ghosts each side
%  Bc: cell centered padded bed 1 x (N+4)
%  Bm,Bp: face bed arrays padded length N+4 at faces k=-1..N+2
% Outputs:
%  dU: 2xN RHS array
%  FL, FR: fluxes 2 x (N+1) arrays
%  S: source term 2xN array

N = length(Bc) - 4; % padded length, interior length N=length(Bc)-4
h = Up(1,:);
hu = Up(2,:);
tol=1e-10;
u = zeros(size(h));
wet = (h > tol);
u(wet) = hu(wet)./h(wet);

% Hydrostatic reconstruction (order 3 with WENO3 faces)
eta = h + Bc;

% Left and right face states eta,u (padded 2 ghosts)
[etaL, etaR] = weno3_faces(eta);
[uL, uR] = weno3_faces(u);

% Beds at faces as seen from left and right cells
BL = Bp(2:N+2);
BR = Bm(3:N+3);

hL0 = max(0, etaL - BL);
hR0 = max(0, etaR - BR);
Bstar = max(BL, BR);
hL = max(0, etaL - Bstar);
hR = max(0, etaR - Bstar);

F = hll_flux(hL, uL, hR, uR, g);

corrL = g/2*(hL0.^2 - hL.^2);
corrR = g/2*(hR0.^2 - hR.^2);
FL = F;
FR = F;
FL(2,:) = FL(2,:) + corrL;
FR(2,:) = FR(2,:) + corrR;

h_m = eta(3:N+2) - Bm(3:N+2);
h_p = eta(3:N+2) - Bp(3:N+2);

S = zeros(2,N);
S(2,:) = -g * 0.5 * (h_m + h_p) .* (Bp(3:N+2) - Bm(3:N+2))/dx;

dU = -(FL(:,2:N+1) - FR(:,1:N))/dx + S;

end

function [Hf, HUf, levmap] = prolong_composite(U, exists)
% Prolong composite solution from hierarchy active arrays to finest grid by 
% piecewise constant copying active cells data and level map.
Lmaxp1 = length(U);
Nfine = length(U{Lmaxp1});
Hf = zeros(1,Nfine);
HUf = zeros(1,Nfine);
levmap = zeros(1,Nfine);

covered = cell(1,Lmaxp1);
covered{Lmaxp1} = false;
for l = 1:Lmaxp1-1
  covered{l} = exists{l+1}(1:2:end) & exists{l+1}(2:2:end);
end

covered{Lmaxp1} = false(1,Nfine);

% From coarsest to finest level, for each level find indices of cells that correspond
% to fine cells, assign data where exists and not covered by finer level cells
for l = 1:Lmaxp1
  idx = ceil((1:Nfine)/2^(Lmaxp1 - l));
  actcells = exists{l}(idx) & ~covered{l}(idx);
  % Assign only where actcells true
  Hf(actcells) = U{l}(1,idx(actcells));
  HUf(actcells) = U{l}(2,idx(actcells));
  levmap(actcells) = l - 1;
end

end

function [Unew, nact] = ssprk3_hier(U, dt, exists, active, covered, Bcp, Bmp, Bpp, dx, g, bc, Lmax)
  nact = sum(cellfun(@nnz, active));
  L0 = amr_rhs(U , exists, active, covered, Bcp, Bmp, Bpp, dx, g, bc, Lmax);
  U1 = comb(1,   U, 0, U , 1, dt, L0);                 % U1 = U + dt*L(U)
  U1 = average_down(U1, covered, Lmax);
  L1 = amr_rhs(U1, exists, active, covered, Bcp, Bmp, Bpp, dx, g, bc, Lmax);
  U2 = comb(3/4, U, 1/4, U1, 1/4, dt, L1);             % U2 = 3/4 U + 1/4 (U1 + dt L(U1))
  U2 = average_down(U2, covered, Lmax);
  L2 = amr_rhs(U2, exists, active, covered, Bcp, Bmp, Bpp, dx, g, bc, Lmax);
  Unew = comb(1/3, U, 2/3, U2, 2/3, dt, L2);           % Unew = 1/3 U + 2/3 (U2 + dt L(U2))
  Unew = average_down(Unew, covered, Lmax);
end

function W = comb(a, A, b, B, c, dt, dU)
  Lmaxp1 = length(A);
  W = cell(1,Lmaxp1);
  for l=1:Lmaxp1
    W{l} = a*A{l} + b*B{l} + c*dt*dU{l};
  end
end

function U = average_down(U, covered, Lmax)
% Average from fine to coarse for covered cells only (vectorized)
% U{l}(:,covered{l}) = 0.5*(U{l+1}(:,2*idx-1)+U{l+1}(:,2*idx)) with idx = find(covered{l})
  for l = Lmax:-1:1
    idx = find(covered{l});
    if isempty(idx)
      continue;
    end
    Ufine = U{l+1};
    Ucoarse = U{l};
    Ucoarse(:, idx) = 0.5*(Ufine(:, 2*idx-1) + Ufine(:, 2*idx));
    U{l} = Ucoarse;
  end
  % Clamp negative depths at all levels after avg down
  for l = 1:(Lmax+1)
    U{l}(1,:) = max(U{l}(1,:),0);
  end
end
