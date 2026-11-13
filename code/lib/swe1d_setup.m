function P = swe1d_setup(cfg)
% Setup grid and padded bed arrays from the configuration struct

% Defaults
if ~isfield(cfg,'g'), cfg.g=9.81; end
if ~isfield(cfg,'scheme'), cfg.scheme='hr'; end
if ~isfield(cfg,'order'), cfg.order=3; end
if ~isfield(cfg,'rk'), cfg.rk='ssprk3'; end
if ~isfield(cfg,'bc'), cfg.bc='wall'; end
if ~isfield(cfg,'N'), error('swe1d_setup: must specify N'); end
if ~isfield(cfg,'xl') || ~isfield(cfg,'xr'), error('swe1d_setup: must specify xl,xr'); end
N=cfg.N; xl=cfg.xl; xr=cfg.xr; order=cfg.order; scheme=cfg.scheme; bc=cfg.bc;

dx=(xr-xl)/N;
xc = xl + ((1:N) - 0.5)*dx; % cell centers 1xN
xf = xl + (0:N)*dx;         % faces 1x(N+1)

% Bed function handle f from bed struct
bed = cfg.bed;

% Allocate padded beds: size 1x(N+4)
P.Bc=zeros(1,N+4); % cell center values
P.Bm=zeros(1,N+4); % bed at left face in cell (hydrostatic left)
P.Bp=zeros(1,N+4); % bed at right face in cell (hydrostatic right)

% Compute interior physical cell beds first (indices 3:(N+2))
P.Bc(3:N+2) = bed.f(xc,0);

if order==3
  % faces for cells: left face = xf(i), right face=xf(i+1), i=1..N
  % for i=1:N: Bm(i) = B(xf(i), +1), Bp(i) = B(xf(i+1), -1)
  % But stored padded, so P.Bm(3:N+2) = B(xf(1:N), +1)
  %            P.Bp(3:N+2) = B(xf(2:N+1), -1)
  P.Bm(3:N+2) = bed.f(xf(1:N), +1);
  P.Bp(3:N+2) = bed.f(xf(2:N+1), -1);
elseif order==1 && strcmp(scheme,'hr')
  % order 1 hydrostatic reconstruction uses Bm=Bp=Bc (cell center) to get zero source term
  P.Bm(3:N+2) = P.Bc(3:N+2);
  P.Bp(3:N+2) = P.Bc(3:N+2);
else
  % For other cases (including 'central' scheme order 1)
  % no spec given; safest to default to order 3 formulation:
  P.Bm(3:N+2) = bed.f(xf(1:N), +1);
  P.Bp(3:N+2) = bed.f(xf(2:N+1), -1);
end

% Pad ghost cells: 2 cells each side, indices 1,2 and N+3,N+4
% Apply the same boundary conditions to all arrays and pads
switch bc
  case 'wall'
    % mirror cells; for bed fields:
    % Bc: mirror P.Bc(3,:), P.Bc(4,:) ... Bc(3) <-> Bc(2) mirrored over boundary, but swapping Bm<->Bp for ghost cells
    % We pad Bc by reflection:
    % Left ghost cells (1,2):
    P.Bc(2) = P.Bc(3);
    P.Bc(1) = P.Bc(4);
    % Right ghost cells (N+3, N+4):
    P.Bc(N+3) = P.Bc(N+2);
    P.Bc(N+4) = P.Bc(N+1);

    % For Bm and Bp in ghost cells swap Bm<->Bp
    % Left ghost cells: 
    P.Bm(2) = P.Bp(3); % swap inside ghost cell
    P.Bp(2) = P.Bm(3);
    P.Bm(1) = P.Bp(4);
    P.Bp(1) = P.Bm(4);

    % Right ghost cells:
    P.Bm(N+3)= P.Bp(N+2);
    P.Bp(N+3)= P.Bm(N+2);
    P.Bm(N+4)= P.Bp(N+1);
    P.Bp(N+4)= P.Bm(N+1);

  case 'periodic'
    % cyclic for all arrays
    P.Bc(1:2) = P.Bc(N+1:N+2);
    P.Bc(N+3:N+4) = P.Bc(3:4);
    P.Bm(1:2) = P.Bm(N+1:N+2);
    P.Bm(N+3:N+4) = P.Bm(3:4);
    P.Bp(1:2) = P.Bp(N+1:N+2);
    P.Bp(N+3:N+4) = P.Bp(3:4);

  case 'open'
    % copy edge values:
    P.Bc(1:2) = P.Bc(3);
    P.Bc(N+3:N+4) = P.Bc(N+2);
    P.Bm(1:2) = P.Bm(3);
    P.Bm(N+3:N+4) = P.Bm(N+2);
    P.Bp(1:2) = P.Bp(3);
    P.Bp(N+3:N+4) = P.Bp(N+2);

  otherwise
    error('swe1d_setup: unknown bc');
end

% Output, include inputs unchanged plus computed arrays
P.N = N;
P.xl = xl;
P.xr = xr;
P.dx = dx;
P.xc = xc;
P.xf = xf;
P.bed = bed;
P.scheme = scheme;
P.order = order;
P.bc = bc;
P.g = cfg.g;
P.rk = cfg.rk;
if isfield(cfg,'cfl')
  P.cfl = cfg.cfl;
else
  % default CFL: 0.4 for order 3, 0.9 for order 1
  if order==3, P.cfl=0.4; else P.cfl=0.9; end
end

end
