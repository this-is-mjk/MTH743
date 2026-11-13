# code/

```
code/
├── setup_paths.m            addpath helper         (setup_paths())
├── run_tests.m              runs all self-tests    (ok = run_tests())
├── run_all.sh               runs every experiment in its own MATLAB process
├── lib/                     numerical library (see table below)
├── tests/                   self-tests (lake at rest, conservation, AMR vs uniform, 2D)
└── experiments/             one script per study; writes figures to outputs/figures,
                             tables to outputs/results (logs of run_all.sh: outputs/logs),
                             all created next to code/ on the first run
```

## Library (`lib/`)

| file | role |
|---|---|
| `make_bed.m` | bed profiles (`flat`, `gaussian`, `cosbump`, `step`, `sin2`, `slope`) with one-sided limits for jumps |
| `hll_flux.m` | HLL flux, Einfeldt speeds, dry-safe, branch-free |
| `weno3_faces.m` | WENO-3 face values (left- and right-biased) |
| `swe1d_setup.m`, `swe1d_pad.m` | grid, padded bed arrays, ghost cells (`wall`, `periodic`, `open`) |
| `swe1d_rhs.m` | spatial operator: schemes `central` (naive) / `hr` (hydrostatic reconstruction), orders 1 and 3 |
| `swe1d_run.m` | driver: CFL time loop, `euler`/`ssprk2`/`ssprk3`/`rk4`, monitors, blow-up guards |
| `amr1d_run.m` | 1D block AMR (ratio 2, any number of levels), flux replacement, nested regridding |
| `swe2d_run.m` | 2D unsplit solver, `plain` and `ctu` (transverse-corrected) |

## Quick start

```matlab
cd code
setup_paths; ok = run_tests();          % ~3 s, must print PASS for four tests
cfg = struct('N',200,'xl',0,'xr',10,'bed',make_bed('gaussian',[0.4 5 0.5]), ...
             'init',@(x,Bc) deal(max(0,1-Bc),0*x),'tend',10,'tsave',10, ...
             'scheme','hr','order',3,'rk','ssprk3','bc','wall','monitor',true);
out = swe1d_run(cfg);  max(abs(out.HU{end}))      % ~1e-15
```

```sh
./run_all.sh            # all experiments (figures + tables)  |  ./run_all.sh p1 amr   for a subset
```