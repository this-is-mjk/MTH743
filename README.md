# Well-balanced high-resolution finite volume schemes for the shallow water equations

MTH743 research project - Manas Jain Kuniya (230626), Department of Mathematics, IIT Kanpur - Instructor: Abhijit Biswas.

Reports and presentation (PDF, in [`docs/report/`](docs/report/)):
[proposal](docs/report/MTH743-proposal.pdf), [progress report](docs/report/MTH743-progress-report.pdf) (warm-ups and Phases 1-3),
[final report](docs/report/MTH743-final-report.pdf) (all four phases and conclusion), [presentation](docs/report/MTH743-presentation.pdf).

## The problem
The 1D shallow water equations over a bed $B(x)$,
$h_t + (hu)_x = 0$, $(hu)_t + (hu^2 + \tfrac12 g h^2)_x = -g h B_x$,
have a stationary *lake at rest* ($u=0$, $h+B=\text{const}$) in which the flux gradient and the bed source cancel exactly. Standard finite volume schemes discretise the two terms separately and do not preserve this
equilibrium: they create spurious waves of the size of the truncation error - larger than a small tsunami.

In this project I build, test and compare schemes that preserve the equilibrium to machine precision
(**hydrostatic reconstruction**), and add high order (**WENO-3 + SSP Runge-Kutta**), local mesh refinement (**AMR**) and a 2D extension with **transverse (CTU) corrections**, all in MATLAB.

## Results (details in [`docs/results.md`](docs/results.md))

| Measurement | baseline | final scheme |
|---|---|---|
| Lake at rest, Gaussian bump, N = 200, t = 10: max\|u\| [m/s] | 6.3e-2 | 1.8e-14 |
| Smooth test, observed order at N = 400 | 0.88 (order 1 scheme) | 2.99 |
| Perturbation 1e-3 on a lake over a steep bump: L1 error | 0.628 | 0.0182 |
| Dam break over a step, N = 400: L1 error of h | 0.145 (hr-1) | 0.0415 |
| Tsunami over a shelf: RMS deviation from Green's law | - | 0.84 % |
| Dam break, uniform N = 800 vs AMR (3 levels): L1(h), cell updates | 0.0186, 2.44e6 | 0.0199, 7.89e5 |
| 2D circular dam break: largest stable CFL, unsplit vs CTU | 0.6 | 1.2 |

## Repository map

```
README.md                  this file
docs/
  theory.md                the mathematics, with proofs (start with sections 4-6)
  results.md               all results and inference (numbers come from the tables the experiments write)
  report/                  proposal, progress report, final report, presentation (PDF)
code/
  lib/ tests/ experiments/ the solver library, self-tests and one script per study
  run_tests.m, run_all.sh  run the tests / all experiments
```
Running the experiments creates `outputs/` (figures, tables, logs); it is not part of the submission, everything in it can be regenerated.
See [`code/README.md`](code/README.md) for the file-level description of the code.

## Run it

Requires MATLAB (developed on R2026b).

```matlab
cd code
setup_paths; ok = run_tests();      % four self-tests, ~3 s, must print ALL OK = 1
```
```sh
cd code && ./run_all.sh      # every experiment: figures -> outputs/figures, tables -> outputs/results (created on first run)
./run_all.sh p1 amr                 # a subset (names: beds b1_advection b2_burgers b3_dam_flat b4_balance_law p1 p2 p3 pert dam rk tsunami amr amr_lake 2d 2d_stability)
```
## References Used

LeVeque, *Finite Volume Methods for Hyperbolic Problems* (2002); Mishra, *Numerical methods for hyperbolic PDEs*; Audusse-Bouchut-Bristeau-Klein-Perthame, SIAM J. Sci. Comput. 25 (2004);
Jiang-Shu, J. Comput. Phys. 126 (1996); Gottlieb-Shu-Tadmor, SIAM Rev. 43 (2001); Berger-Colella, J. Comput. Phys. 82 (1989).

## Use of AI assistance

I used Claude (Anthropic) as a helper for writing the MATLAB code and for writing the reports (LaTeX reports, presentation and Markdown notes). The topic, scope and phases follow my own proposal (`docs/report/MTH743-proposal.pdf`).
The code is checked by the self-tests in `code/tests`, against exact or fine-grid reference solutions, and by inspecting the figures; the defects found on the way are listed in the progress and final reports.
