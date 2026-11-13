% test_lib_b: tests of the 1D solver (lake at rest for every scheme and order, mass conservation)

tol_balanced = 1e-12;
tol_unbalanced = 1e-4;
tol_mass = 1e-12;
g = 9.81;

% (a) Lake at rest tests: gaussian bed [0.4 5 0.5], N=100, [0,10], 'wall', tend=2
N = 100;
xl = 0; xr = 10;
bed = make_bed('gaussian', [0.4 5 0.5]);
init_lake = @(x,Bc) deal(max(0,1-Bc), zeros(size(x)));

tend = 2;
% orders and schemes to test:
schemes = {'hr','central'};
orders = [1 3];
% Storage of results for printing
fprintf('Lake at rest test results:\n');
for sch = schemes
  for ord = orders
    if strcmp(sch{1},'central') && ord == 3
      continue % skip central order3 (not specified)
    end
    cfg.N = N;
    cfg.xl = xl; cfg.xr = xr;
    cfg.bed = bed;
    cfg.order = ord;
    cfg.scheme = sch{1};
    cfg.rk = 'ssprk3';
    cfg.bc = 'wall';
    cfg.cfl = 0.4; % safe
    cfg.init = init_lake;
    cfg.tend = tend;
    cfg.tsave = tend;

    % Run simulation
    out = swe1d_run(cfg);

    % Check max |hu| and max |h+B - eta0|
    H = out.H{end};
    HU = out.HU{end};
    eta = H + out.B;
    hu_max = max(abs(HU));
    eta0 = 1;
    err_eta = max(abs(eta - eta0));

    fprintf('Scheme %-8s order %d: max|hu|=%.2e max|h+B-1|=%.2e\n', sch{1}, ord, hu_max, err_eta);

    if strcmp(sch{1}, 'hr')
      assert(hu_max < tol_balanced, 'hr scheme failed velocity balance');
      assert(err_eta < tol_balanced, 'hr scheme failed free surface balance');
    elseif strcmp(sch{1}, 'central') && ord == 1
      fprintf('central order 1 max|hu|=%.2e (threshold %.2e)\n', hu_max, tol_unbalanced);
      assert(hu_max > tol_unbalanced, 'central scheme unexpectedly well balanced');
    end
  end
end

% (b) Mass conservation in dam break on flat bed (B=0), 'wall', order 3
bed_flat = make_bed('flat', []);
N = 100;
xl = 0; xr = 10;
init_dambreak = @(x,Bc) deal( (x < 5)*2 + (x >= 5)*1 , zeros(size(x)) );

cfg.N = N;
cfg.xl = xl; cfg.xr = xr;
cfg.bed = bed_flat;
cfg.order = 3;
cfg.scheme = 'hr';
cfg.rk = 'ssprk3';
cfg.bc = 'wall';
cfg.cfl = 0.4;
cfg.init = init_dambreak;
cfg.tend = 2;
cfg.tsave = 2;

out = swe1d_run(cfg);

Hfinal = out.H{end};
[h0chk,~] = init_dambreak(out.x, out.B);
mass_initial = sum(h0chk);
mass_final = sum(Hfinal);

mass_err = abs(mass_final - mass_initial);
rel_err = mass_err / mass_initial;
fprintf('Dam break mass conservation absolute error = %.2e, relative error = %.2e\n', mass_err, rel_err);
assert(rel_err < tol_mass, 'Dam break mass not conserved');

fprintf('lib_b OK\n');
