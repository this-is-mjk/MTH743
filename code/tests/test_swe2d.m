% test_swe2d: tests of the 2D solver (lake at rest, mass conservation, symmetry)
g = 9.81;

% Common parameters
tol_vel = 1e-12;
tol_eta = 1e-12;
tol_mass = 1e-12;
Nx1 = 40; Ny1 = 40;
Nx2 = 60; Ny2 = 60;
xl = 0; xr = 1; yl = 0; yr = 1;

% (a) Lake at rest, wall BC, both methods
bed = @(X,Y) 0.4*exp(-((X-0.5).^2 + (Y-0.5).^2)/0.02);
eta0 = 1;
init_lake = @(X,Y,B) cat(3, eta0 - B, zeros(size(X)), zeros(size(X)));

cfg1 = struct;
cfg1.Nx = Nx1; cfg1.Ny = Ny1; cfg1.xl = xl; cfg1.xr = xr; cfg1.yl = yl; cfg1.yr = yr;
cfg1.bed2 = bed;
cfg1.init = init_lake;
cfg1.tend = 1;
cfg1.tsave = 1;
cfg1.bc = 'wall';
cfg1.g = g;

methods = {'plain','ctu'};
for m=methods
  cfg1.method = m{1};
  out = swe2d_run(cfg1);
  H = out.H{end};
  HU = out.HU{end};
  HV = out.HV{end};
  eta = H + bed(out.x', out.y);
  max_u = max(abs(HU(:)./max(H(:),1e-16)));
  max_v = max(abs(HV(:)./max(H(:),1e-16)));
  max_eta_err = max(abs(eta(:)-eta0));
  assert(max_u<tol_vel,'Lake at rest u violation');
  assert(max_v<tol_vel,'Lake at rest v violation');
  assert(max_eta_err<tol_eta,'Lake at rest eta violation');
end

% (b) Mass conservation wall BC circular dam break flat bed
bed2 = @(X,Y) zeros(size(X));
r0 = 0.2;
xc = 0.5; yc = 0.5;
init_dam = @(X,Y,B) cat(3, double(((X-xc).^2+(Y-yc).^2) < r0^2)*2 + 1*double(((X-xc).^2+(Y-yc).^2) >= r0^2), ...
  zeros(size(X)), zeros(size(X)));
cfg2 = struct;
cfg2.Nx = Nx2; cfg2.Ny = Ny2; cfg2.xl = xl; cfg2.xr = xr; cfg2.yl = yl; cfg2.yr = yr;
cfg2.bed2 = bed2;
cfg2.init = init_dam;
cfg2.tend = 0.1;
cfg2.tsave = 0.1;
cfg2.bc = 'wall';
cfg2.g = g;

dx = (cfg2.xr - cfg2.xl) / cfg2.Nx;
dy = (cfg2.yr - cfg2.yl) / cfg2.Ny;

for m=methods
  cfg2.method = m{1};
  out = swe2d_run(cfg2);

  % Build grid for init_dam evaluation
  [X,Y] = ndgrid(out.x, out.y);
  B = zeros(size(X)); % flat bed
  d0 = init_dam(X, Y, B);
  h0 = d0(:,:,1);
  mass_init = sum(h0(:)) * dx * dy;

  Hend = out.H{end};
  mass_final = sum(Hend(:)) * dx * dy;

  rel_mass_err = abs(mass_final - mass_init)/mass_init;
  assert(rel_mass_err < tol_mass, 'Mass conservation violated');
end

% (c) Symmetry test for dam break
for m=methods
  cfg2.method = m{1};
  out = swe2d_run(cfg2);
  H = out.H{end};
  % Check symmetry: H(i,j) ≈ H(j,i)
  err_sym = max(abs(H - H.'), [], 'all');
  assert(err_sym < tol_eta, 'Symmetry violated');
end

fprintf('swe2d OK\n');
