# Results and inference

All numbers below are copied from the tables (`outputs/results/*.txt`) that the experiment scripts write, or from their logs (`outputs/logs/`); figures are written to `outputs/figures/`. These folders are created when the experiments are run and are not part of the submission. Mathematical background: [`theory.md`](theory.md). Reproduce everything with
`cd code && ./run_all.sh` (a few minutes per experiment on an Apple-silicon laptop; AMR and 2D are the longest).

Setup common to all runs: $g=9.81$, hydrostatic reconstruction ("hr") unless stated; orders: 1 = piecewise constant + forward Euler, 3 = WENO-3 + SSP-RK3.

| Question | Answer |
|---|---|
| Does the baseline scheme fail? | Yes: on N = 200 cells spurious velocity 3.8 cm/s (peak 6.3 cm/s) after t = 10; converges only like $\Delta x$ (velocity) / $\Delta x^{0.7}$ (surface) |
| Does hydrostatic reconstruction fix it? | Yes: $\max\lvert u\rvert = 10^{-15}\ldots10^{-13}$ (orders 1 and 3, three bathymetries, N = 50-400, t = 10); with a refined patch (AMR) and in 2D too |
| Is high order achieved? | WENO-3 + SSP-RK3: observed order **2.99** (N = 400) on a smooth test, 48x-206x smaller error than first order |
| Does it matter for small waves? | Yes: for a perturbation of 1e-3 of the depth the naive scheme's spurious signal is larger than the wave itself (scaled max 3.0 against ~0.5 expected), the well-balanced third-order scheme resolves it (34x smaller L1 error) |
| Tsunami shoaling | crest amplitude follows Green's law within **0.84%** RMS (x in 35-80 m) with equilibrium preserved to 3e-14 |
| Dam break over a step, AMR | AMR with 3 levels is within **7%** of the uniform-fine-grid error with **3.1x fewer cell updates** (but 4.5x *slower* in MATLAB wall time - see caveats) |
| Time integrators | max CFL 0.80 / 1.88 / 1.98 / 2.08 for Euler / SSP-RK2 / SSP-RK3 / RK4; SSP-RK3 kept (SSP, 2 registers) |
| 2D transverse corrections | CTU reduces the L1 error of a circular dam break by 27-37%, is 1.5-2x more isotropic, and doubles the stable time step (largest stable CFL 0.6 -> 1.2; same cost at the limit) |

---

## R0. Warm-up experiments (basic ideas before the main problem)

Four small experiments, each isolating one idea (scripts `exp_b1_advection`, `exp_b2_burgers`, `exp_b3_dam_flat`, `exp_b4_balance_law`; tables in `results/b1.txt` ... `b4.txt`).

**b1: linear advection** ($u_t + u_x = 0$, periodic, one period, exact cell averages as initial data). $L^1$ error of a smooth pulse and observed order from $N=320$ to $640$:

| scheme | $L^1$ error, $N=640$ | order |
|---|---|---|
| upwind | $5.20\times10^{-3}$ | 0.96 |
| Lax-Wendroff | $1.11\times10^{-4}$ | 2.00 |
| MUSCL-minmod + SSP-RK2 | $7.71\times10^{-4}$ | 1.82 |
| WENO-3 + SSP-RK3 | $1.69\times10^{-4}$ | 2.74 |
| WENO-3 + forward Euler | $1.03\times10^{-2}$ | 0.93 |

For a Gaussian plus a square wave ($N=200$, one period) the extreme values are: upwind $[0.000, 0.999]$ (smears), Lax-Wendroff $[-0.195, 1.195]$ (oscillates), MUSCL $[0.000, 0.998]$, WENO-3 + SSP-RK3 $[-0.000, 1.000]$, WENO-3 + Euler $[-0.007, 1.007]$ (total variation grows from 3.997 to 4.733 after five periods). **Inference.** Smooth accuracy and monotonicity are different goals; WENO-3 gives both, but only with an SSP integrator (forward Euler makes it first order and mildly unstable). Figures: `b1_advection_profiles.png`, `b1_advection_convergence.png`.

**b2: Burgers equation** ($u_L=1$, $u_R=0$, exact shock at $x=0.5$ at $t=0.4$). Godunov and Rusanov (conservative) give $0.5000$ and $0.5001$ at $N=800$ (errors shrinking with $N$) and conserve mass ($0.5000$); the non-conservative upwind scheme leaves the shock at $x=0.3000$ for every $N$ and the mass at $0.3000$. **Inference.** Only conservative schemes give the right shock speed (Lax-Wendroff theorem), hence all schemes here are in conservation form. Figure: `b2_burgers.png`.

**b3: dam break on a flat bed** (exact Riemann solution; wet $h_L=2$, $h_R=1$, $t=0.8$; dry $h_L=1$, $h_R=0$, $t=0.6$). Star state $h^*=1.4538$, $u^*=1.3058$, shock speed $4.183$. $L^1$ error of $h$ at $N=800$ (order between $N=400$ and $800$):

| scheme | wet | dry |
|---|---|---|
| Rusanov | $5.18\times10^{-2}$ (0.80) | $4.25\times10^{-2}$ (0.73) |
| HLL | $4.04\times10^{-2}$ (0.84) | $2.54\times10^{-2}$ (0.72) |
| WENO-3 + SSP-RK3 | $8.92\times10^{-3}$ (1.16) | blows up at every $N$ |

HLL is 22% (wet) and 40% (dry) more accurate than Rusanov. First-order HLL is positive (smallest depth 0) and stable up to CFL 1.1 (blow-up at 1.5). **Inference.** Errors with shocks and rarefactions converge slower than order 1 for all schemes; **WENO-3 blows up on a dry bed** (reconstruction of $u=q/h$ next to $h\approx0$): wet/dry fronts need a positivity-preserving treatment (open item). Figures: `b3_dam_wet.png`, `b3_dam_dry.png`, `b3_dam_convergence.png`.

**b4: scalar balance law** ($u_t + u_x = -u$, steady state $e^{-x}$; naive source versus a source integrated along the transport step). Starting from the exact steady state and running to $t=5$, the naive scheme drifts with error $7.4\times10^{-3}$ ($N=25$) to $2.3\times10^{-4}$ ($N=800$), order 1.00, the balanced scheme keeps it to $3\times10^{-15}\ldots1.6\times10^{-13}$. A pulse of $10^{-3}$ on the steady state ($N=100$, $t=0.4$, amplitude $6.7\times10^{-4}$) is reproduced by the balanced scheme with error 0.21 of its size, and drowned by the naive scheme (error $1.31\times10^{-3}$, twice the signal). **Inference.** The lake-at-rest problem in miniature: the cure is a scheme with the right steady state, not a smaller $\Delta x$. Figure: `b4_balance_law.png`.

---

## R1. Phase 1 - the baseline scheme destroys the lake at rest (`exp_p1`)

Lake at rest $\eta_0=1$ over $B = 0.4\,e^{-(x-5)^2/0.5}$, walls, first-order HLL + central-difference source, $t=10$.

| N | $\max\lvert u\rvert$ (end) | $\max\lvert u\rvert$ over time | $L^1$ error of $\eta$ (end) | order ($\eta$) |
|---|---|---|---|---|
| 50 | 1.60e-1 | 2.15e-1 | 2.70e-2 | - |
| 100 | 8.08e-2 | 1.20e-1 | 1.35e-2 | 1.00 |
| 200 | 3.77e-2 | 6.33e-2 | 9.41e-3 | 0.52 |
| 400 | 1.67e-2 | 3.25e-2 | 6.08e-3 | 0.63 |

**Inference.** The error is *consistent* but converges slowly: the final-time velocity drops by orders 0.99 / 1.10 / 1.17, the L1 error of the surface with mean order 0.72.
On a resolution that is already generous for a 10 m domain (N = 200, dx = 5 cm), spurious velocities of 3.8-6.3 cm/s remain, i.e. the same order as the particle velocity of a real tsunami in the deep ocean (about 5 cm/s for a 1 m wave on 4000 m, $u = \eta\sqrt{g/h}$). The failure is exactly the one predicted in `theory.md` section 4.2.
Figures: `p1_maxu_vs_t.png`, `p1_u_times.png`, `p1_eta_times.png`, `p1_animation.mp4`.

## R2. Phase 2 - hydrostatic reconstruction is well balanced (`exp_p2`)

Lake at rest, $t=10$, three beds (Gaussian, cosine bump, narrow tall Gaussian), N = 50...400.

| scheme | worst $\max\lvert u\rvert$ over all beds and $N$ | worst $\max\lvert\eta-\eta_0\rvert$ |
|---|---|---|
| central-1 (naive) | 6.3e-1 | 4.6e-2 |
| hr-1 | 1.1e-14 | 8e-16 |
| hr-3 (WENO-3 + SSP-RK3) | 1.3e-13 | 1.3e-13 |

**Inference.** The well-balanced schemes are stationary to round-off (the proposal's expected result 1, "machine-precision equilibrium, about 1e-15"): hr-1 stays within 1e-14, hr-3 within 1.3e-13. The slow growth of the hr-3 residual with N (about 1e-14 -> 1e-13 for N = 50 -> 400) is the accumulation of round-off over $O(N)$ more time steps and more cells, not a balance error: the structure of the proof in `theory.md` 6.2 leaves nothing that could scale with a discretisation parameter.
The naive scheme's error is 13-14 orders of magnitude larger and decays only slowly with N (`p2_maxu_vs_N.png`).

## R3. Phase 3 - WENO-3 + SSP-RK3 (`exp_p3`)

Smooth test (periodic, $B=\sin^2(\pi x)$, $h = 5 + e^{\cos 2\pi x}$, $hu = \sin(\cos 2\pi x)$, $t=0.1$), L1 error against the same scheme on N = 3200:

| N | order 1: L1(h) | obs. order | order 3: L1(h) | obs. order |
|---|---|---|---|---|
| 25 | 1.19e-1 | - | 5.13e-2 | - |
| 50 | 8.09e-2 | 0.55 | 1.77e-2 | 1.54 |
| 100 | 5.04e-2 | 0.68 | 4.02e-3 | 2.14 |
| 200 | 2.88e-2 | 0.81 | 6.05e-4 | 2.73 |
| 400 | 1.57e-2 | 0.88 | 7.60e-5 | **2.99** |

($hu$: 3.04 at N = 400.) **Inference.** The high-order scheme reaches its design order 3 asymptotically; at the same N the third-order scheme is 48x more accurate at N = 200 (6.05e-4 vs 2.88e-2) and 206x at N = 400, for 7.6x the cost per run (nupdates 1.58e6 vs 2.09e5 at N = 200).
The well-balanced property is not lost by the high-order machinery: lake at rest $\max\lvert u\rvert = 7.4\times10^{-15}$, while a *naive WENO-3* (reconstruct $h$, $u$, central source) gives $\max\lvert u\rvert = 0.118$. So WENO accuracy alone does not give well-balancing - the reconstruction variable ($\eta$) and the face-bed treatment do.
Numerical diffusion (`p3_pulse_diffusion.png`, narrow pulse, N = 200, t = 0.3): peak of $\eta-1$ is 0.113 (order 1), 0.147 (order 3) and 0.139 (reference): first order loses 19% of the peak, WENO-3 overshoots it by 5.6%.
*Note:* we have not investigated the 5.6% overshoot (probably the nonlinear steepening of a pulse of relative height 0.2 plus the lack of a limiter/positivity fix); it is smaller than the first-order loss and is listed as an open item.

## R4. Small perturbations of a lake at rest (`exp_pert`, LeVeque's quasi-steady test)

$x\in[0,1]$, $B = 0.25\bigl(\cos(10\pi(x-0.5)) + 1\bigr)$ on $\lvert x-0.5\rvert<0.1$, $\eta = 1+\varepsilon$ on $0.1<x<0.2$, open boundaries, N = 100, $t=0.2$; error of $(\eta-1)/\varepsilon$ against a fine hr-3 reference.

| eps | scheme | L1 error | $\max\lvert(\eta-1)/\varepsilon\rvert$ (numerical solution) |
|---|---|---|---|
| 0.2 | central-1 | 4.47e-2 | 0.29 |
| 0.2 | hr-1 | 4.37e-2 | 0.29 |
| 0.2 | hr-3 | 2.71e-2 | 0.40 |
| 0.001 | central-1 | **6.28e-1** | **3.04** |
| 0.001 | hr-1 | 4.47e-2 | 0.30 |
| 0.001 | hr-3 | **1.82e-2** | 0.54 |

The reference is the hr-3 solution on N = 3200, restricted to the 100 coarse cells by the block mean of the **surface** $h+B$ (see the correction note below). For a small-amplitude pulse linear theory predicts two half-pulses of height 0.5; the reference indeed shows a flat-topped pulse of height 0.5 (`pert_eps0010.png`), so the last column can be read against 0.5.

**Inference.** For a large perturbation (eps = 0.2) the three schemes are comparable (the physical signal swamps the lake-at-rest error; the naive and the first-order well-balanced scheme even have the same error). For eps = 0.001 - the realistic tsunami regime - the naive scheme's maximum scaled deviation is 3.04, six times the 0.5 expected for a half-pulse, i.e. its spurious signal is *larger than the wave itself*, and its L1 error is **34x** that of hr-3 (6.28e-1 vs 1.82e-2) and 14x that of hr-1. The well-balanced schemes resolve the wave; the third-order one is 2.5x more accurate than the first-order one. This is the central practical argument for the project's method.

hr-3 convergence for eps = 0.001 (N = 50 -> 400): errors 2.66e-2, 1.82e-2, 1.10e-2, 6.92e-3, observed orders 0.55, 0.72, 0.67. These are low because the solution contains discontinuities (square pulse): for a captured discontinuity a scheme of order p converges in L1 only like $\Delta x^{p/(p+1)}$, i.e. about 0.75 for p = 3, which is what is observed. The accuracy is limited by the data, not by the scheme.

**Correction note (2026-10-14).** An earlier version of this experiment formed the reference surface as the block mean of the *depth* plus the bed at the coarse cell centre. On the sharply curved bump $\overline{B} - B(x_c) \sim B''\,\Delta x^2/24 \sim 5\times10^{-4}$, which at eps = 0.001 is half of the signal: it appeared as a spurious half-height pulse at x = 0.5 in the reference and inflated the errors of the well-balanced schemes (the hr-3 error was 8.3e-2 instead of 1.8e-2, the ratio to the naive scheme 6.8 instead of 34, and the orders 1.7/1.6/1.3 were an artefact). It was found while inspecting `pert_eps0010.png` and fixed in `exp_pert.m`. The eps = 0.2 numbers are unchanged to 3 digits because the artefact is 0.2% of that signal.
Figures: `pert_eps0010.png`, `pert_eps2000.png`, `pert_eps0001.mp4`.

## R5. Dam break over a step (`exp_dam`, proposal expected result 2)

$[0,20]$, depths 2 | 1 | 0.5 above a step of height 0.5 at $x=12$ (surface 2 | 1 | 1), $t=2$; reference N = 6400.

| method | N | L1(h) | L1(hu) | nupdates |
|---|---|---|---|---|
| hr-1 | 400 | 1.45e-1 | 5.35e-1 | 1.64e5 |
| hr-3 | 400 | 4.15e-2 | 1.47e-1 | 6.13e5 |
| hr-3 | 100 | 1.76e-1 | 6.49e-1 | 3.84e4 |

Mass drift: 0 (hr-1), 2e-14 (hr-3, N = 400), 3e-13 (reference): conservation to round-off. **Inference.** WENO-3 gives 3.5x smaller error than first order on the same grid and sharpens the bore front and the rarefaction head (`dam_eta_t2.png`);
hr-3 with N = 100 has an error comparable to hr-1 with N = 400 (1.76e-1 vs 1.45e-1, 21% larger) at 23% of the cost. Figures: `dam_eta_t2.png`, `dam_u_t2.png`, `dam_eta_xt_stacked.png`, `dam_animation.mp4`.

## R6. Time integrators (`exp_rk`)

Moving pulse over a cosine bump, WENO-3 + hr, N = 100, $t=0.5$, error = *temporal* error against the same scheme with SSP-RK3 at CFL 0.05.

| integrator | stages | order | SSP | registers | largest stable CFL | error at CFL 0.4 | error at 0.8 x limit |
|---|---|---|---|---|---|---|---|
| Euler | 1 | 1 | yes | 1 | 0.80 | 1.36e-2 | 5.4e-2 |
| SSP-RK2 | 2 | 2 | yes | 2 | 1.88 | 2.9e-4 | 1.5e-2 |
| SSP-RK3 | 3 | 3 | yes | 2 | 1.98 | 8.5e-6 | 1.2e-2 |
| RK4 | 4 | 4 | no | 3 | 2.08 | 1.2e-6 | 1.3e-2 |

(CFL is relative to the *global* maximum wave speed, so values above 1 are possible.) Error versus cost (`rk_efficiency.png`): at roughly equal cost (87k vs 93k updates) SSP-RK3 has error 2.9e-5 and RK4 9.2e-6; at ~116k updates SSP-RK3 1.2e-5 vs RK4 3.9e-6.
**Inference.** Order matters most: going from Euler to SSP-RK2 gains 47x in accuracy at the same CFL; SSP-RK3 another 34x. In the smooth regime RK4 is about 3x more accurate than SSP-RK3 *per unit cost*, but it is not SSP (no positivity/TVD guarantee for the nonlinear shallow-water problem with dry fronts or shocks) and needs 3 registers; the project therefore keeps SSP-RK3 as the production integrator
(2 registers, SSP coefficient 1, third order), as the proposal intended. Near the stability limit all schemes lose accuracy (errors of 1e-2 at 0.8 x limit): the useful range is CFL 0.4-1.
Forward Euler cannot be recommended with WENO-3 (limit 0.8 and 1.4e-2 error at CFL 0.4).

## R7. Tsunami over a shelf (`exp_tsunami`)

Depth 10 m for $x\le30$, linear shoaling to 1 m at $x=90$, $\eta_0 = 0.05\,e^{-((x-15)/3)^2}$, hr-3, N = 1000. Lake-at-rest sanity run on the same grid: $\max\lvert u\rvert = 3.4\times10^{-14}\ \mathrm{m/s}$.
The right-going crest (initial amplitude $A_0 = 0.02498$) is tracked and compared with Green's law $A(x) = A_0\,\bigl(h_0/h(x)\bigr)^{1/4}$:

| x [m] | measured [m] | Green [m] | difference |
|---|---|---|---|
| 40 | 0.02611 | 0.02601 | +0.40% |
| 50 | 0.02738 | 0.02731 | +0.25% |
| 60 | 0.02899 | 0.02897 | +0.05% |
| 70 | 0.03121 | 0.03136 | -0.48% |
| 80 | 0.03470 | 0.03535 | -1.83% |
| 90 | 0.04102 | 0.04433 | -7.47% |

RMS relative deviation for $35<x<80$: **0.84%**. **Inference.** The scheme reproduces the linear shoaling law with 0.84% RMS deviation (worst 1.8% at x = 80) while the crest grows by a factor 1.6 between x = 35 and x = 90, while preserving the equilibrium ($3.4\times10^{-14}$) - the "small-amplitude signal over complex bathymetry is not destroyed by numerical noise" claim of the proposal. Note that this is a model problem (the pulse is narrower than the depth, so real water would be dispersive; the shallow water *equations* have no dispersion and are what is being tested). The departure close to the shore (-7.5% at 1 m depth, where $\eta/h = 4\%$) is the expected effect of nonlinearity and steepening, which Green's law (linear theory) neglects.
Figures: `tsunami_amplitude_comparison.png`, `tsunami_eta_times.png`, `tsunami_xt_contour_eta.png`, `tsunami_animation.mp4`.

## R8. Adaptive mesh refinement (`exp_amr`, `test_amr`; proposal expected result 2)

Library tests (`code/tests/test_amr.m`): lake at rest with a fixed refined patch preserved to **1e-15**; with a $10^{-3}$ pulse and 1, 2, 3 levels all populated the far field (not yet reached by the wave) stays within 8e-11, 1e-12, 7e-14 of $\eta=1$ (it falls as levels are added; `outputs/results/amr_lake.txt`), mass drift **5e-15** with walls, and for the dam break over a step AMR (N0 = 100, 2 levels) differs from the uniform 400-cell grid by L1(h) = **3.1e-5** (L1(hu) = 1.2e-4) using **40.5%** of the cell updates. The threshold controls refinement as it should: `thr` = 2e-2 / 2e-3 / 2e-4 give L1(h) = 8.0e-5 / 3.1e-5 / 6.1e-6 with 39 / 41 / 44% of the updates.

Experiment (dam break over the step, $t=2$, errors against the N = 6400 reference, compared on the finest 800-cell grid):

| method | L1(h) | L1(hu) | cell updates | error / uniform-800 | saving vs uniform-800 |
|---|---|---|---|---|---|
| uniform N = 100 | 1.76e-1 | 6.49e-1 | 3.81e4 | 9.48 | 64x |
| uniform N = 800 | 1.86e-2 | 6.59e-2 | 2.44e6 | 1.00 | 1x |
| **AMR N0 = 100, 3 levels** | **1.99e-2** | 1.12e-1 | 7.89e5 | **1.07** | **3.1x** |
| AMR N0 = 200, 2 levels | 1.88e-2 | 7.67e-2 | 9.93e5 | 1.01 | 2.5x |

**Inference.** AMR reproduces the accuracy of the uniform fine grid: its $L^1$ error in $h$ is 1.01-1.07 times that of the fine grid (the proposal's "WENO + AMR accurately resolves hydraulic jumps while saving computational overhead") with 3.1x fewer cell updates; the momentum error is somewhat larger (1.7x for the 3-level run); we did not diagnose why (candidates: the 3-cell flag dilation lets a bore leave the finest patch between regrids, and the coarse-fine flux replacement acts on the momentum flux of one side only).
**Caveats, stated plainly:** (i) the saving is in *cell updates*, not wall time - in this MATLAB implementation (full-length arrays per level, ghost filling in loops) the AMR run is about 4.4x *slower* in seconds (2.24 s vs 0.51 s in the last run; timings vary between runs); wall-clock gains would need a compact data structure or compiled code. (ii) There is no sub-cycling: the coarse levels use the fine time step, which limits the achievable saving to the ratio of active cells. (iii) The gradient indicator refines the whole smooth rarefaction fan to the finest level (`amr_level_map.png`), which is wasteful; a curvature or error-estimator (e.g. second differences of $\eta$) would keep it coarse.
*Threshold sweep* (`amr_sensitivity_thr.png`, 3 levels, `thr` from 1e-1 to 1e-4 of the surface range): with `thr = 1e-1` the error is 3.9e-2 (refinement too late); for every `thr <= 3e-2` the error is the same 1.99e-2, equal to the fine-grid limit, while the work grows from about 7.4e5 to 9.2e5 updates. So here the accuracy is saturated by the finest level and the cheapest threshold that reaches it (about 3e-2) is the right choice; tightening it only costs updates.
Figures: `amr_eta_compare.png`, `amr_eta_zoom.png`, `amr_level_map.png`, `amr_active_cells_time.png`, `amr_sensitivity_thr.png`, `amr_animation.mp4`.

## R9. 2D and transverse corrections (`exp_2d`, `exp_2d_stability`; proposal expected result 3, "if time permits")

First-order hydrostatic reconstruction on $[0,1]^2$, walls, two variants: `plain` (unsplit, $U + \Delta t\,(L_x + L_y)$) and `ctu` (transverse predictor, theory section 8.3).

*Lake at rest*, 2D Gaussian bump, 60x60, $t=2$: plain $\max\lvert hu\rvert=\max\lvert hv\rvert = 5.6\times10^{-16}$, $\max\lvert\eta-1\rvert = 4.4\times10^{-16}$; CTU $9.8\times10^{-16}$, $3.3\times10^{-16}$. Both are balanced to round-off.

*Circular dam break* ($h=2$ for $r<0.2$, else 1, $t=0.1$), L1 error of $h$ against the 400x400 CTU solution restricted by block averaging:

| method | 50x50 | 100x100 |
|---|---|---|
| plain | 4.04e-2 | 2.49e-2 |
| CTU | 2.93e-2 | 1.57e-2 |

CTU removes 27% (50x50) and 37% (100x100) of the error, and the error drops by a factor 1.6 (plain) and 1.9 (CTU) when the grid is refined, i.e. observed order 0.7-0.9, as expected for a shock with first-order schemes.

*Isotropy* (the proposal's "diagonal wave propagation"): a circular dam break must stay radially symmetric. At $t=0.05$, when the shock (radius about 0.4) has not yet reached any wall (distance 0.5), we compare the profile along an axis with the profile along a diagonal at equal radii (`2d_circ_dam_break_radprof.png`) and report $\dfrac{\overline{\lvert h_{\mathrm{diag}}(r) - h_{\mathrm{axis}}(r)\rvert}}{\max h - \min h}$ for $r\le0.45$:

| method | 50x50 | 100x100 |
|---|---|---|
| plain | 2.66e-2 | 1.97e-2 |
| CTU | 1.32e-2 | 1.33e-2 |
| reference 400x400 CTU | 8.3e-3 | |

CTU is 2.0x (50x50) and 1.5x (100x100) more isotropic than the unsplit scheme. The 400x400 reference is itself not perfectly isotropic (0.8%), so the metric has a floor of that size. *Why $t=0.05$:* at $t=0.1$ (used for the L1 table) the shock has already been reflected on the axes (walls at $r=0.5$) but not on the diagonals (corner at $r=0.71$), so an axis-versus-diagonal comparison would measure the wall, not the scheme; an earlier version of this experiment did exactly that (and, in addition, indexed the diagonal wrongly), which is why its anisotropy numbers and radial-profile figure were meaningless and were redone.

*Plane wave at 45 degrees* ($\eta = 1 + 0.05\,e^{-((x+y-0.5)/0.08)^2}$, i.e. centred on the diagonal coordinate $s=(x+y)/\sqrt2=0.5$, open boundaries, 100x100, $t=0.12$): the pulse splits into two travelling in opposite directions (peaks near $s=0.12$ and $s=0.88$, `2d_plane_wave_diag.png`). Peak of $\eta-1$ (the larger, left pulse, which is close to the open boundary) and its full width at half maximum, against the 400x400 CTU reference:

| | peak | FWHM |
|---|---|---|
| plain | 0.0217 | 0.947 |
| CTU | 0.0232 | 0.929 |
| reference | 0.0249 | 0.914 |

CTU loses 6.7% of the peak (plain 12.8%) and broadens the wave by 1.6% (plain 3.6%). The right pulse shows the same ordering (heights about 1.0193, 1.0210 and 1.0237 for plain, CTU and the reference): transverse terms reduce the numerical diffusion for waves not aligned with the grid.

*Time-step limit* (`exp_2d_stability`, 50x50 circular dam break to $t=0.3$, CFL based on the global maximum speed; a run is called *unstable* when the deviation from the method's own CFL = 0.3 solution grows by more than a factor 3 between consecutive CFL values, or $\max h$, $\min h$ leave a band around the baseline):

| method | largest stable CFL | first unstable | first blow-up | sweeps at the largest stable CFL |
|---|---|---|---|---|
| plain | 0.6 | 0.7 | 1.2 | 234 |
| CTU | 1.2 | 1.4 | 1.8 | 244 |

**Inference.** The transverse correction doubles the usable time step (ratio 2.0 on the scanned CFL grid with spacing 0.1, so the true ratio lies between 1.7 and 2.3), as theory predicts ($\nu\le1$ per direction for CTU, $\nu_x+\nu_y\le1$ for the unsplit scheme). Because a CTU step costs twice as many sweeps (4 vs 2), the *cost at each method's stability limit is the same* (234 vs 244 sweeps), so the gain of CTU is accuracy at equal cost and the
absence of grid imprinting for oblique waves, not speed. Only the first-order scheme is available in 2D.
Figures: `2d_circ_dam_break_surfs.png`, `2d_circ_dam_break_radprof.png`, `2d_plane_wave_diag.png`, `2d_stability.png`, `2d_dam_break.mp4`.

---

## What the project verified (against the proposal's expected results)

1. **Machine-precision equilibrium** - confirmed: $\max\lvert u\rvert$ = 1e-15...1e-13 for orders 1 and 3, on three bathymetries, with AMR levels and in 2D (5e-16...1e-15); the naive scheme gives 0.03-0.6 m/s (R1, R2).
2. **Shock-capturing on adaptive grids** - confirmed for the dam break over a step: AMR has 1.01-1.07 times the $L^1$ error (h) of the uniform fine solution with up to 3.1x fewer updates (R5, R8).
3. **Multidimensional validation** - implemented and verified on lake at rest (5e-16), mass conservation and symmetry; CTU reduces the error by 27-37%, broadens a diagonal wave less, and doubles the stable CFL (R9). The proposal's "correctly model complex diagonal wave propagation" is shown for a plane wave and a circular dam break on flat/smooth beds only.

## Limitations and open items

* WENO-3 pulse overshoot of 5.6% at N = 200 (R3) is unexplained.
* Wet cells only: dry-front handling (positivity-preserving reconstruction, "wet/dry interface") is *not* implemented although the proposal mentions dry/wet interfaces; the HLL flux is dry-safe but no runup test is included.
* The AMR work counter ignores the cost of regridding/ghost filling; wall time is worse than uniform in MATLAB (see R8).
* The AMR time step has no sub-cycling; the refinement indicator is gradient-based.
* The 2D scheme is first order in space (no WENO in 2D) and only tested on flat/smooth beds with walls (isotropy and plane-wave tests use flat beds).
* The code is verified by the four self-tests, independent measurements (exact Riemann solution, Green's law, fine-grid references) and inspection of every figure; the defects found and fixed on the way are listed in the progress and final reports. The code has not been audited line by line.
