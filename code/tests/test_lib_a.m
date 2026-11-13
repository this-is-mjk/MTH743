% test_lib_a: unit tests of the basic library (bed profiles, HLL flux, WENO-3 faces)

% Test hll_flux with h,u same left and right: flux equals physical flux
g = 9.81;
h = [0, 0.1, 1, 2]; % row vector
u = [0, 1, -0.5, 2];
F = hll_flux(h,u,h,u,g);
F_exact = [h.*u; h.*u.^2 + g*h.^2/2];
assert(all(abs(F(:)-F_exact(:))<1e-14));

% Test weno3_faces on linear ramp and constant
N=10;
vfull = (0:N+3)/(N+3); % padded ramp v(1..N+4)
% Remove tiny floating errors in input (ensure exactly increasing):
[vL,vR] = weno3_faces(vfull);
% Exact face values for linear function: face k lies between padded cells k+1 and k+2,
% so exact face value is (k+0.5)/(N+3), k=1..N+1
xface = ((1:N+1) + 0.5)/(N+3);
assert(all(abs(vL - xface) < 1e-13));
assert(all(abs(vR - xface) < 1e-13));

% Const test: constant input stays constant to within 1 ulp
vconst = ones(1,N+4)*7.13;
[vL,vR] = weno3_faces(vconst);
assert(all(abs(vL - 7.13) < 1e-14));
assert(all(abs(vR - 7.13) < 1e-14));

% make_bed step test: at x=5 with side=-1 returns B0=0, side=+1 returns B1=1
bed = make_bed('step',[0 1 5]);
assert(bed.f(5,-1) == 0);
assert(bed.f(5,1) == 1);

fprintf('lib_a OK\n')
