% test_amr: tests of the AMR code (lake at rest, mass, comparison with the uniform solver)

g = 9.81;

%% (a) Lake at rest, gaussian bed, fixed flag

% Setup config
N0 = 40; Lmax = 2;
xl = 0; xr = 10;

bed = make_bed('gaussian',[0.4 5 0.5]);
init = @(x,Bc) deal(max(0,1-Bc), zeros(size(x))); % eta0=1, u=0, lake at rest: depth = max(0, 1 - bed)

fixedflag = @(x) (x>=3.5 & x<=6.5);

cfg = struct();
cfg.N0 = N0;
cfg.xl = xl; cfg.xr = xr;
cfg.bed = bed;
cfg.init = init;
cfg.tend = 2;
cfg.tsave = [0 2];
cfg.g = g;
cfg.bc = 'wall';
cfg.Lmax = Lmax;
cfg.nregrid = 5;
cfg.thr = 2e-3;
cfg.cfl = 0.4;
cfg.mode = 'fixed';
cfg.fixedflag = fixedflag;

out = amr1d_run(cfg);
H = out.H{end};
HU = out.HU{end};
assert(max(abs(HU(:)))<1e-12,'(a) Lake at rest failed: velocity not zero');

Bfine = bed.f(out.x,0); 
lev = out.level{end};  % level index (0..Lmax) of each fine cell
Bavg = Bfine;
for i=1:numel(Bfine)
    b = 2^(Lmax - lev(i));
    j0 = floor((i-1)/b)*b;
    Bavg(i) = mean(Bfine(j0+1:j0+b));
end
eta = H + Bavg;
assert(max(abs(eta - 1))<1e-12,'(a) Lake at rest failed: eta variation');
fprintf('(a) Lake at rest passed\n');

%% (b) Conservation of mass dam break on flat bed adaptive refine

N0 = 50; Lmax=2;
xl = 0; xr = 10;

bed = make_bed('flat',[]);
init = @(x,Bc) deal(2*(x<5) + 1*(x>=5), zeros(size(x))); % dam break h=2|1, u=0

cfg2 = struct();
cfg2.N0 = N0;
cfg2.xl = xl; cfg2.xr = xr;
cfg2.bed = bed;
cfg2.init = init;
cfg2.tend = 0.5;
cfg2.tsave = [0 0.5];
cfg2.g = g;
cfg2.bc = 'wall';   % changed from 'open' to 'wall'
cfg2.Lmax = Lmax;
cfg2.nregrid = 5;
cfg2.thr = 2e-3;
cfg2.cfl = 0.4;
cfg2.mode = 'adaptive';

out2 = amr1d_run(cfg2);
dx_finest = (xr-xl)/(N0*2^Lmax);
mass_init = sum(out2.H{1} .* dx_finest);
mass_final = sum(out2.H{end} .* dx_finest);

mass_diff = abs(mass_final - mass_init)/mass_init;
fprintf('(b) Conservation passed, initial mass: %.15g, final mass: %.15g, relative difference: %.3g\n', mass_init, mass_final, mass_diff);

assert(mass_diff < 1e-12, '(b) Conservation test failed: mass drift');

%% (c) Dam break over a step: compare amr with uniform 400 cells at hr3 ssprk3

xl = 0; xr = 20;
bed = make_bed('step',[0 0.5 12]);
init = @(x,Bc) deal(2*(x<8) + 1*(x>=8 & x<12) + 0.5*(x>=12), zeros(size(x))); % h, u=0

tend = 1;
cfl = 0.4;

% Uniform reference config
cfgu = struct();
cfgu.N = 400;
cfgu.xl = xl; cfgu.xr = xr;
cfgu.bed = bed;
cfgu.init = init;
cfgu.tend = tend;
cfgu.tsave = 1; % scalar tsave
cfgu.g = g;
cfgu.bc = 'wall';
cfgu.scheme = 'hr';
cfgu.order = 3;
cfgu.rk = 'ssprk3';
cfgu.cfl = cfl;

outu = swe1d_run(cfgu);

% AMR config
cfgamr = struct();
cfgamr.N0 = 100;
cfgamr.xl = xl; cfgamr.xr = xr;
cfgamr.bed = bed;
cfgamr.init = init;
cfgamr.tend = tend;
cfgamr.tsave = 1; % scalar tsave
cfgamr.g = g;
cfgamr.bc = 'wall';
cfgamr.Lmax = 2;
cfgamr.nregrid = 5;
cfgamr.thr = 2e-3;
cfgamr.cfl = cfl;
cfgamr.mode = 'adaptive';

outamr = amr1d_run(cfgamr);

% Extract solutions on finest grid
h_amr = outamr.H{end};
hu_amr = outamr.HU{end};
x_amr = outamr.x;
lev = outamr.level{end};
dx_fine = (xr - xl) / (cfgamr.N0 * 2^cfgamr.Lmax);

% Uniform solution at 400 points
h_uni = outu.H{end};
hu_uni = outu.HU{end};
x_uni = outu.x;

% Compute L1 norms by interpolating uniform onto AMR finest grid to ensure alignment
% Here both have same finest resolution: 400 points on uniform, AMR finest grids also maximum 400 cells,
% but AMR points are non-uniform; we will interpolate uniform data to AMR points

% Interpolate uniform solution h and hu at AMR x points
h_uni_at_amr = interp1(x_uni, h_uni, x_amr, 'linear', 'extrap');
hu_uni_at_amr = interp1(x_uni, hu_uni, x_amr, 'linear', 'extrap');

L1h = sum(abs(h_amr - h_uni_at_amr)) * dx_fine;
L1hu = sum(abs(hu_amr - hu_uni_at_amr)) * dx_fine;

fprintf('(c) L1 norm h AMR vs uniform: %.6g\n', L1h);
fprintf('(c) L1 norm hu AMR vs uniform: %.6g\n', L1hu);

% Check updates ratio
nup_amr = outamr.nupdates;
nup_uni = outu.nupdates;
ratio = nup_amr / nup_uni;

fprintf('(c) Updates: AMR = %d, Uniform = %d, ratio = %.3g\n', nup_amr, nup_uni, ratio);

assert(L1h < 0.05, '(c) L1 norm h difference too large');
assert(L1hu < 0.15, '(c) L1 norm hu difference too large');
assert(ratio < 0.8, '(c) AMR updates ratio too large');

fprintf('amr OK\n');
