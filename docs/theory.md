# Mathematics Used

Well-balanced high-resolution finite volume schemes for the shallow water equations (SWE).
Here I list the formulas that are implemented in `code/lib/`, with proofs wherever possible, to understand why each scheme does (or does not) preserve the lake at rest.

References:

1. **[L]** = LeVeque, *Finite Volume Methods for Hyperbolic Problems* (2002)
2. **[M]** = Mishra, *Numerical methods for hyperbolic PDEs*, lecture notes
3. **[A]** = Audusse, Bouchut, Bristeau, Klein, Perthame (2004), *A fast and stable well-balanced scheme with hydrostatic reconstruction for shallow water flows*, SIAM J. Sci. Comput. 25.

**Notation**

| symbol | meaning |
|---|---|
| $h$ | water depth |
| $u$ | velocity |
| $q = hu$ | discharge |
| $B(x)$ | bed elevation |
| $\eta = h + B$ | free-surface elevation |
| $g = 9.81\ \mathrm{m/s^2}$ | gravity |
| $c = \sqrt{g h}$ | gravity-wave speed |
| $C_i = [x_{i-1/2},\,x_{i+1/2}]$ | cell $i$ |
| $\Delta x$ | cell width |
| $U_i^n \approx \dfrac{1}{\Delta x}\displaystyle\int_{C_i} U(x,t^n)\,dx$ | cell average of $U = (h,q)^T$ |

---

## 1. The model

### 1.1 Balance law (mass and momentum)

$$
\begin{aligned}
h_t + q_x &= 0,\\
q_t + \Big(\frac{q^2}{h} + \frac{g h^2}{2}\Big)_x &= -\,g\,h\,B_x .
\end{aligned}
\tag{1}
$$

In compact form $U_t + F(U)_x = S(U,x)$ with

$$
U = \begin{pmatrix} h \\ q \end{pmatrix},\qquad
F(U) = \begin{pmatrix} q \\ q^2/h + g h^2/2 \end{pmatrix},\qquad
S(U,x) = \begin{pmatrix} 0 \\ -g h B_x \end{pmatrix}.
$$

See [L, 13.1, 18.7].

### 1.2 Hyperbolicity

The Jacobian is

$$
F'(U) = \begin{pmatrix} 0 & 1 \\ g h - u^2 & 2u \end{pmatrix},
\qquad
\lambda_{1,2} = u \mp c,\quad c = \sqrt{g h},
$$

real and distinct for $h > 0$, so the system is strictly hyperbolic. Both fields are genuinely nonlinear, so solutions contain shocks (hydraulic jumps, bores) and rarefactions; this is why conservative schemes and Riemann solvers are used ([L, 13.3-13.8]).

The CFL restriction is

$$
\Delta t \le \mathrm{CFL}\;\frac{\Delta x}{\max_i\,(|u_i| + c_i)} .
$$

### 1.3 The steady state that matters: lake at rest

For $u = 0$ the momentum equation reduces to $\big(\tfrac12 g h^2\big)_x = -g h B_x$. Since $\big(\tfrac12 g h^2\big)_x = g h\,h_x$,

$$
g\,h\,(h + B)_x = 0
\quad\Longrightarrow\quad
\eta = h + B = \text{const on wet cells.}
\tag{2}
$$

So the stationary solution is $u = 0$, $\eta = \eta_0$. The *physical* pressure term $g h h_x$ and the *source* $-g h B_x$ cancel exactly, pointwise.

A scheme is **well balanced** if its discrete flux difference and discrete source cancel exactly for discrete data satisfying (2), so that the discrete lake at rest is an exact fixed point of the scheme (up to round-off).

Why this matters: a tsunami in the open ocean has $\eta - \eta_0 \sim 1\ \mathrm{m}$ on a depth of $4000\ \mathrm{m}$, a relative amplitude of about $2.5\times10^{-4}$. A scheme whose lake-at-rest error is $O(\Delta x)$ or $O(\Delta x^2)$ times the bed variation produces false waves larger than the signal.

---

## 2. Finite volume framework

Integrating (1) over $C_i \times [t^n, t^{n+1}]$ gives the exact identity

$$
U_i^{n+1} = U_i^n - \frac{\Delta t}{\Delta x}\big(F_{i+1/2} - F_{i-1/2}\big) + \Delta t\, S_i ,
\tag{3}
$$

where $F_{i+1/2}$ is the time-averaged flux through the interface $x_{i+1/2}$. Summing over $i$ the flux terms telescope, so $\sum_i U_i\,\Delta x$ changes only through the boundary fluxes: **mass is conserved to round-off** with reflective walls (relative drift about $10^{-14}$ in all tests).

Method of lines:

$$
\frac{dU_i}{dt} = L_i(U) = -\frac{F_{i+1/2} - F_{i-1/2}}{\Delta x} + S_i ,
$$

advanced with an SSP Runge-Kutta scheme (section 6.3).

---

## 3. Numerical flux: HLL with Einfeldt wave speeds ([L, 15.3.5], [M, Ch. 9])

For left and right states $U_L, U_R$ at an interface choose speeds $s_L \le s_R$ that bound the Riemann fan:

$$
\hat u = \frac{\sqrt{h_L}\,u_L + \sqrt{h_R}\,u_R}{\sqrt{h_L} + \sqrt{h_R}},\qquad
\hat c = \sqrt{\tfrac12 g\,(h_L + h_R)},
$$

$$
s_L = \min\big(u_L - c_L,\ \hat u - \hat c\big),\qquad
s_R = \max\big(u_R + c_R,\ \hat u + \hat c\big).
\tag{4}
$$

(Einfeldt / Roe-averaged speeds. For a dry left state use $s_L = u_R - 2c_R$, $s_R = u_R + c_R$, and the mirror image for a dry right state.) The HLL flux is

$$
F_{\mathrm{HLL}} =
\begin{cases}
F(U_L), & s_L \ge 0,\\[4pt]
\dfrac{s_R F(U_L) - s_L F(U_R) + s_L s_R\,(U_R - U_L)}{s_R - s_L}, & s_L < 0 < s_R,\\[10pt]
F(U_R), & s_R \le 0.
\end{cases}
\tag{5}
$$

The implementation uses the branch-free form with $s_L^- = \min(s_L, 0)$ and $s_R^+ = \max(s_R, 0)$:

$$
F = \frac{s_R^+ F_L - s_L^- F_R + s_L^- s_R^+\,(U_R - U_L)}{s_R^+ - s_L^-}.
$$

Properties used below: **consistency** $F(U,U) = F(U)$, **conservation**, and for $U_L = U_R$ the numerical-diffusion term vanishes.

The sign matters: the term is $+\,s_L s_R\,(U_R - U_L)$. With the opposite sign the scheme is *anti*-diffusive and every perturbation grows exponentially (this exact bug was found and fixed in the 2D code, see the Phase 4 problems table in the final report).

---

## 4. Why the baseline scheme fails (Phase 1)

### 4.1 The naive scheme

Cell-centred FVM with the HLL flux on $(h,q)$, ignoring $B$ in the flux, and a central-difference source:

$$
S_i = \begin{pmatrix} 0 \\[2pt] -\,g\,h_i\,\dfrac{B_{i+1} - B_{i-1}}{2\,\Delta x} \end{pmatrix}.
\tag{6}
$$

### 4.2 Lake-at-rest residual (exact calculation)

Take $u = 0$, $\eta = \eta_0$ and $h_i = \eta_0 - B_i$, so $B_{i+1} - B_{i-1} = -(h_{i+1} - h_{i-1})$. With $q_L = q_R = 0$ the speeds are $s_L = -\hat c$, $s_R = +\hat c$ and (5) gives

$$
F^{h}_{i+1/2} = -\frac{\hat c}{2}\,(h_{i+1} - h_i),
\qquad
F^{q}_{i+1/2} = \frac{g}{4}\,\big(h_i^2 + h_{i+1}^2\big).
$$

The mass flux is pure numerical diffusion, non-zero because $h$ varies with the bed. Hence

$$
\frac{dh_i}{dt} = \frac{\hat c}{2\Delta x}\,\big(h_{i+1} - 2h_i + h_{i-1}\big) = \frac{\hat c\,\Delta x}{2}\,h'' + O(\Delta x^3) = O(\Delta x),
$$

and the depth is not stationary. For the momentum, using (6),

$$
\begin{aligned}
\frac{dq_i}{dt}
&= -\frac{g}{4\Delta x}\big(h_{i+1}^2 - h_{i-1}^2\big) + \frac{g\,h_i}{2\Delta x}\big(h_{i+1} - h_{i-1}\big)\\
&= -\frac{g}{2\Delta x}\,\big(h_{i+1} - h_{i-1}\big)\Big(\frac{h_{i+1} + h_{i-1}}{2} - h_i\Big)\\
&= -\frac{g\,\Delta x^2}{2}\,h'(x_i)\,h''(x_i) + O(\Delta x^4).
\end{aligned}
$$

So the *momentum* residual is $O(\Delta x^2)$ but the *mass* leak is $O(\Delta x)$. The leak changes $h$, destroys $\eta = \text{const}$, and the pressure imbalance then accelerates the water, giving spurious velocities of size $O(\Delta x)$ (plus waves bouncing between the walls). Both residuals are proportional to $h''$, the bed curvature: they live exactly where the bed bends.

Measured (Gaussian bump, $t = 10$, `outputs/results/p1.txt`): $\max|u|$ at the end of the run decreases with observed orders $0.99,\ 1.10,\ 1.17$ ($N = 50 \to 400$), only about first order, and the $L^1$ error of $\eta$ with a mean order of $0.72$. On $N = 200$ cells the spurious speed is $3.8\ \mathrm{cm/s}$ at the end of the run ($6.3\ \mathrm{cm/s}$ at its peak), the same order as the particle velocity of a real tsunami in the open ocean ($u = \eta\sqrt{g/h}$; for a $1\ \mathrm{m}$ wave on $4000\ \mathrm{m}$ depth about $5\ \mathrm{cm/s}$).

---

## 5. Hydrostatic reconstruction (Phase 2) [A]

### 5.1 Definition

At interface $i+\tfrac12$ with left and right data $(\eta_L, u_L, B_L)$ and $(\eta_R, u_R, B_R)$ define

$$
\begin{aligned}
B^* &= \max(B_L, B_R),\\
h_L &= \max(0,\ \eta_L - B_L), & h_R &= \max(0,\ \eta_R - B_R) &&\text{(unmodified depths)},\\
h_L^* &= \max(0,\ \eta_L - B^*), & h_R^* &= \max(0,\ \eta_R - B^*) &&\text{(reconstructed depths)},
\end{aligned}
\tag{7}
$$

$$
F = F_{\mathrm{HLL}}\Big(\big(h_L^*,\ h_L^* u_L\big),\ \big(h_R^*,\ h_R^* u_R\big)\Big),
$$

$$
F_L = F + \begin{pmatrix} 0 \\ \tfrac{g}{2}\big(h_L^2 - h_L^{*2}\big) \end{pmatrix}
\ \ \text{(flux seen by the cell on the left)},
\qquad
F_R = F + \begin{pmatrix} 0 \\ \tfrac{g}{2}\big(h_R^2 - h_R^{*2}\big) \end{pmatrix}
\ \ \text{(flux seen by the cell on the right)}.
$$

The update is

$$
U_i^{n+1} = U_i - \frac{\Delta t}{\Delta x}\Big(F_L\big(x_{i+1/2}\big) - F_R\big(x_{i-1/2}\big)\Big) + \Delta t\, S_i .
$$

For the first-order scheme the bed is piecewise constant, $B_L = B_i$, $B_R = B_{i+1}$, and $S_i = 0$.

### 5.2 Theorem (well-balancing)

*If $u_i = 0$ and $\eta_i = \eta_0$ for all $i$, the first-order scheme leaves the data unchanged exactly.*

*Proof.* Since $\eta_L = \eta_R = \eta_0$, we get $h_L^* = h_R^* = \max(0, \eta_0 - B^*) =: h^*$, and $u_L = u_R = 0$ gives equal reconstructed states $U_L^* = U_R^*$. By consistency of the HLL flux, $F = F(U^*) = \big(0,\ \tfrac12 g h^{*2}\big)^T$. Therefore

$$
F_L\big(x_{i+1/2}\big) = \begin{pmatrix} 0 \\ \tfrac g2 h^{*2} + \tfrac g2\big(h_i^2 - h^{*2}\big) \end{pmatrix} = \begin{pmatrix} 0 \\ \tfrac g2 h_i^2 \end{pmatrix},
\qquad
F_R\big(x_{i-1/2}\big) = \begin{pmatrix} 0 \\ \tfrac g2 h_i^2 \end{pmatrix}.
$$

(Cell $i$ is the *left* cell of face $i+\tfrac12$ and the *right* cell of face $i-\tfrac12$, so both corrections involve only $h_i$.) The flux difference is $\tfrac g2 h_i^2 - \tfrac g2 h_i^2 = 0$ and $S_i = 0$, so $U_i^{n+1} = U_i^n$. $\blacksquare$

The proof uses only: (a) $\eta$ is the reconstructed variable; (b) $F(U,U) = F(U)$; (c) the left and right corrections cancel inside the cell. It also holds with *dry* cells ($h^* = 0$), and the scheme preserves positivity of $h$ under the CFL condition because $h^*$ is built from physical (non-negative) depths.

Measured: $\max|u| \approx 2\times10^{-15}\ldots10^{-14}$ for three different bumps and $N = 50\ldots400$ (`p2.txt`, figure `p2_maxu_vs_N.png`).

Moving flows are *not* altered at leading order: when the bed is smooth, $B_R - B_L = O(\Delta x)$, so $h_L - h_L^* = O(\Delta x)$ and the correction is a consistent $O(\Delta x)$ discretisation of the source term (first-order accurate overall), as the convergence tables in `p3_smooth_accuracy.txt` show.

---

## 6. High-order reconstruction and time stepping (Phase 3)

### 6.1 WENO-3 ([M, 6.3-6.4], Jiang-Shu)

The face value from the left, $v^-_{i+1/2}$, uses cells $i-1, i, i+1$:

$$
\begin{aligned}
p_0 &= -\tfrac12 v_{i-1} + \tfrac32 v_i, & p_1 &= \tfrac12 v_i + \tfrac12 v_{i+1}, & d_0 &= \tfrac13,\ d_1 = \tfrac23,\\
\beta_0 &= (v_i - v_{i-1})^2, & \beta_1 &= (v_{i+1} - v_i)^2,\\
\alpha_k &= \frac{d_k}{(\varepsilon + \beta_k)^2}, & \omega_k &= \frac{\alpha_k}{\alpha_0 + \alpha_1}, & v^-_{i+1/2} &= \omega_0 p_0 + \omega_1 p_1 .
\end{aligned}
\tag{8}
$$

The ideal weights give the third-order central stencil $-\tfrac16 v_{i-1} + \tfrac56 v_i + \tfrac13 v_{i+1}$. Near a discontinuity the stencil that crosses it gets a weight of order $\beta^{-2}\to 0$, so the reconstruction is essentially non-oscillatory (second order at a smooth extremum, third order in smooth monotone regions). The right-biased value at the same face is the mirror image (cells $i+2, i+1, i$).

Property used below: **a constant is reconstructed exactly** ($\beta_k = 0$, $\omega_0 + \omega_1 = 1$, $p_k = \text{const}$).

### 6.2 Well-balanced high-order scheme

Reconstruct the *surface* $\eta$ and the *velocity* $u$ (not $h$ and $q$), obtaining $\eta^{L,R}_{i+1/2}$ and $u^{L,R}_{i+1/2}$. Take the bed at the faces *analytically*:

$$
B_i^+ = B\big(x_{i+1/2}^-\big)\ \ (\text{left limit}),\qquad
B_{i+1}^- = B\big(x_{i+1/2}^+\big)\ \ (\text{right limit}),
$$

and apply (7) with $B_L = B_i^+$, $B_R = B_{i+1}^-$. Add the cell source

$$
S_i = \begin{pmatrix} 0 \\[2pt] -\,g\,\dfrac{h_i^- + h_i^+}{2}\;\dfrac{B_i^+ - B_i^-}{\Delta x} \end{pmatrix},
\qquad
h_i^{\mp} = \eta_i^{\mp} - B_i^{\mp}.
\tag{9}
$$

**Theorem.** *If $\eta_i = \eta_0$ and $u_i = 0$, the scheme is stationary.*

*Proof.* WENO reproduces constants, so $\eta^{L,R} = \eta_0$ and $u^{L,R} = 0$ at every face. As in 5.2, $F_L(x_{i+1/2}) = \big(0,\ \tfrac g2 (h_i^+)^2\big)^T$ with $h_i^+ = \eta_0 - B_i^+$, and $F_R(x_{i-1/2}) = \big(0,\ \tfrac g2 (h_i^-)^2\big)^T$. The momentum flux difference divided by $\Delta x$ is

$$
-\frac{g}{2\Delta x}\Big((h_i^+)^2 - (h_i^-)^2\Big) = -\frac{g}{2\Delta x}\,\big(h_i^+ + h_i^-\big)\big(h_i^+ - h_i^-\big).
$$

Since $h_i^+ - h_i^- = -(B_i^+ - B_i^-)$, this equals $+\,g\,\dfrac{h_i^+ + h_i^-}{2}\,\dfrac{B_i^+ - B_i^-}{\Delta x}$, which cancels the source (9) exactly. $\blacksquare$

The identity is algebraic: it needs no smoothness of $B$, so it also works across a bed step placed at a face.

Measured: $\max|u| \approx 3\times10^{-15}\ldots10^{-13}$ for orders 1 and 3 over $t = 10$.

*Remark (why not reconstruct $h$).* Reconstructing $h$ and $B$ separately gives $h + B \ne \text{const}$ at the faces and destroys the balance at $O(\Delta x^3)$ times bed derivatives: tiny but not zero. Reconstructing $\eta$ makes the balance exact. A *naive WENO scheme* (reconstruct $h$, $u$, central-difference source) is not balanced: measured $\max|u| = 0.12$ in `p3`.

### 6.3 SSP Runge-Kutta ([M, Ch. 5-6], Gottlieb-Shu-Tadmor, [L, 10.4])

If forward Euler is stable (TVD or positivity-preserving) for $\Delta t \le \Delta t_{\mathrm{FE}}$, every *strong stability preserving* scheme that is a convex combination of Euler steps inherits that stability for $\Delta t \le c\,\Delta t_{\mathrm{FE}}$, with $c$ the SSP coefficient.

SSP-RK2 (Heun), $c = 1$, order 2:

$$
U^{(1)} = U + \Delta t\,L(U),\qquad
U^{+} = \tfrac12 U + \tfrac12\big(U^{(1)} + \Delta t\,L(U^{(1)})\big).
$$

SSP-RK3 (Shu-Osher), $c = 1$, order 3, two registers:

$$
\begin{aligned}
U^{(1)} &= U + \Delta t\,L(U),\\
U^{(2)} &= \tfrac34 U + \tfrac14\big(U^{(1)} + \Delta t\,L(U^{(1)})\big),\\
U^{+} &= \tfrac13 U + \tfrac23\big(U^{(2)} + \Delta t\,L(U^{(2)})\big).
\end{aligned}
$$

Well-balancing is preserved because every stage is a convex combination of Euler steps and $L(U_{\text{rest}}) = 0$.

Measured (`exp_rk.txt`): largest stable CFL $0.80$ (Euler), $1.88$ (SSP-RK2), $1.98$ (SSP-RK3), $2.08$ (classical RK4, not SSP) for WENO-3 with hydrostatic reconstruction on a moving pulse. The CFL here is relative to the *global* maximum speed, so the observed limits exceed 1.

---

## 7. Adaptive mesh refinement (Phase 4a)

### 7.1 Hierarchy

Levels $\ell = 0,\dots,L$ with ratio 2, $\Delta x_\ell = \Delta x_0 / 2^\ell$. A fine cell *exists* inside refined regions; a coarse cell is *covered* if both its children exist. Only *active* cells (existing and not covered) are advanced, and the composite solution is the set of active cells. One global $\Delta t$ (set by the finest level, no sub-cycling) is used, so Runge-Kutta stages need no time interpolation.

### 7.2 Ingredients

1. **Ghost cells** of a fine patch: prolongation of $(\eta, q)$ from the parent level with minmod-limited linear slopes, $\eta_{\text{child}} = \eta_p \pm s/4$ and $h = \max(0,\ \eta - B_{\text{fine}})$. A linear reconstruction with symmetric weights is conservative: the mean of the children equals the parent.
2. **Flux replacement (reflux)**: at a coarse-fine face the coarse cell uses the *fine* flux (same $\Delta t$, same face), so what leaves the coarse cell enters the fine cells exactly, which gives **conservation** (measured drift $5\times10^{-15}$; $10^{-9}$ only when waves leak through *open* boundaries, equal to the uniform-grid value).
3. **Average-down** after every RK stage: the covered coarse value is the mean of its children (conservative restriction).
4. **Regridding** every $n$ steps: flag cells with $\tfrac12|\eta_{i+1} - \eta_{i-1}| > \mathrm{thr}/2^\ell$ (or the same for $q$), dilate by 3 cells, enforce *proper nesting* (parent and both neighbours exist). Existing fine cells keep their data; new cells are prolonged from the *current* coarse solution (never from stale data).

### 7.3 Well-balancing across levels

Choose the cell bed on level $\ell$ as the average of the finest-level cell-centre bed values inside the cell, so that exactly

$$
B_c^{(\ell)}(i) = \tfrac12\Big(B_c^{(\ell+1)}(2i-1) + B_c^{(\ell+1)}(2i)\Big),
$$

and use the **analytic** face beds on every level. Then:

- (i) the $\eta$-prolongation of a constant is exact;
- (ii) the average-down of $h = \eta_0 - B_c^{(\ell+1)}$ gives $\eta_0 - B_c^{(\ell)}$;
- (iii) at a coarse-fine face both levels use the same $h^{\pm} = \eta_0 - B(x_{\text{face}})$, so the replaced flux equals what the coarse cell would compute itself.

By the proof in 6.2 the balance then holds cell by cell: the lake at rest is preserved to about $10^{-15}$ with 1, 2 or 3 levels (test (a) in `test_amr.m`). A naive index error that shifted the bed by two cells produced errors of $10^{-5}$ and was caught only by this test.

---

## 8. Two dimensions and transverse corrections (Phase 4b)

### 8.1 Equations

$$
\begin{aligned}
h_t + (hu)_x + (hv)_y &= 0,\\
(hu)_t + \big(hu^2 + \tfrac12 g h^2\big)_x + (huv)_y &= -\,g h B_x,\\
(hv)_t + (huv)_x + \big(hv^2 + \tfrac12 g h^2\big)_y &= -\,g h B_y .
\end{aligned}
$$

### 8.2 Unsplit scheme

$$
U^{n+1} = U^n + \Delta t\,\big(L_x(U) + L_y(U)\big),
$$

with $L_x$, $L_y$ the 1D hydrostatic-reconstruction operators applied in each direction (the normal momentum plays the role of "the momentum", the tangential momentum $h\,u_n v_t$ is advected). This "donor-cell" combination is exact for waves aligned with the grid but treats a wave travelling at angle $\theta$ by splitting its *flux* into the two axis directions. The cross-derivative term is missing, which gives a direction-dependent numerical diffusion (**grid imprinting**) and a restricted time step $\Delta t\,(|u| + c)\big(\tfrac{1}{\Delta x} + \tfrac{1}{\Delta y}\big) \le 1$.

### 8.3 Corner Transport Upwind ([L, 20.2-20.3, 21.6-21.7])

The CTU method adds the transverse contribution by updating the states used for the normal fluxes with a half-step transverse operator:

$$
U^{y} = U + \tfrac{\Delta t}{2}\,L_y(U),\qquad
U^{x} = U + \tfrac{\Delta t}{2}\,L_x(U),\qquad
U^{n+1} = U + \Delta t\,\Big(L_x\big(U^{y}\big) + L_y\big(U^{x}\big)\Big).
\tag{10}
$$

In LeVeque's wave-propagation language this is the transverse Riemann solve that splits the fluctuation $\mathcal{A}^+\Delta Q$ into up and down contributions. A Taylor expansion shows that (10) contains cross-term corrections of the form $\tfrac12\Delta t^2\,(AB + BA)\,U_{xy}$, so the leading anisotropic error of the donor-cell scheme cancels, and the stability limit rises to $\max(\nu_x, \nu_y) \le 1$.

Here $L_x$ and $L_y$ **include** the hydrostatic corrections and the source, so for the lake at rest $L_x(U) = L_y(U) = 0$ exactly (section 5.2 applied direction by direction), hence $U^x = U^y = U$ and the corrected scheme is well balanced too. The $y$-operator is implemented by permuting the data so that *the same verified $x$-sweep* is reused (no duplicated direction-dependent code).

---

## 9. Verification metrics used in the experiments

- $\max|u|$ and $\max|\eta - \eta_0|$ after a long time for the lake at rest (target about $10^{-15}$).
- $L^1$ error $\sum_i |U_i - U_i^{\mathrm{ref}}|\,\Delta x$ against the same scheme on a much finer grid, restricted to the coarse grid by **exact cell averaging of the compared quantity** (for a surface, the block mean of $h + B$). Observed order $p = \log_2(e_N / e_{2N})$.
- Conservation: relative drift of $\sum_i h_i \Delta x$ (walls).
- Efficiency: $n_{\text{updates}} = n_{\text{steps}} \cdot n_{\text{stages}} \cdot n_{\text{active cells}}$; AMR saving $= n_{\text{updates}}^{\text{uniform fine}} / n_{\text{updates}}^{\text{AMR}}$.
- Stability: largest CFL for which a run finishes with finite, bounded data, judged by a jump in the deviation from the low-CFL solution.
- 2D isotropy: for a radially symmetric initial condition, at a time before the shock reaches a wall, the mean absolute difference between the diagonal profile $h_{\mathrm{diag}}(r)$ and the axis profile $h_{\mathrm{axis}}(r)$ interpolated to the same radii, divided by $\max h - \min h$.
