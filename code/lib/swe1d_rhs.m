function [dU, FL, FR] = swe1d_rhs(U, P)
% Compute dU/dt for 1D SWE with given state U and problem struct P
% Supports schemes: 'central', 'hr'
% Supports order 1,3 + WENO3 faces
% Also returns face fluxes FL,FR for possible later use

g = P.g; N = P.N; dx = P.dx; bc = P.bc; scheme = P.scheme; order = P.order;
Bc = P.Bc; Bm = P.Bm; Bp = P.Bp;
tol = 1e-10;

% Pad state U: 2xN -> 2x(N+4)
Up = swe1d_pad(U, bc);

% Primitive variables:
% For hydrostatic reconstruction, reconstruct eta = h + Bc, u = hu ./ h with u=0 if h<=tol
h = Up(1,:);
hu = Up(2,:);
u = zeros(size(h)); 
wetmask = (h > tol);
u(wetmask) = hu(wetmask) ./ h(wetmask);

switch scheme
  case 'hr'
    % Hydrostatic reconstruction (use eta = h + Bc):
    eta = h + Bc;

    if order == 1
      % piecewise constant states: primitive states constant per cell, so face states constant from cell (left or right)
      % Face states: reconstruct (eta,u) either as constant or WENO-3
      % order 1 special: face arrays uniform, left is cell center of left cell
      % face k between padded cells k+1 and k+2 (k=1..N+1)
      % Left face state from padded cell k+1: etaL = eta(k+1), uL = u(k+1)
      % Right face state from padded cell k+2: etaR= eta(k+2), uR = u(k+2)
      etaL = eta(2:N+2); % faces 1..N+1 left state
      uL = u(2:N+2);
      etaR = eta(3:N+3);
      uR = u(3:N+3);
    elseif order == 3
      % WENO3 face values of eta and u, input vectors padded 1x(N+4)
      [etaL, etaR] = weno3_faces(eta);
      [uL, uR] = weno3_faces(u);
      % Length of etaL,R is (N+1)
      % etaL(k) from cell k+1, etaR(k) from cell k+2
    else
      error('swe1d_rhs: unsupported order');
    end

    % Bed at faces as seen from left (left) and right (right) cell:
    % For face k=1..N+1, from padded arrays Bp and Bm (which have length N+4,
    % face k lies between cell k+1 and k+2, ie beds at columns k+1 (left cell) and k+2 (right cell))
    BL = Bp(2:N+2); % left cell right face bed, Bp aligned with cell centers plus padding
    BR = Bm(3:N+3); % right cell left face bed

    % Hydrostatic Reconstruction:
    % hL0 = max(0, etaL - BL), hR0=max(0, etaR - BR)
    hL0 = max(0, etaL - BL);
    hR0 = max(0, etaR - BR);
    Bstar = max(BL, BR);
    hL = max(0, etaL - Bstar);
    hR = max(0, etaR - Bstar);

    % Velocities unchanged uL, uR
    % Compute HLL fluxes for hL,uL,hR,uR
    F = hll_flux(hL, uL, hR, uR, g); % 2x(N+1)

    % Compute hydrostatic flux corrections (per face)
    % FL = F + [0; g/2*(hL0.^2 - hL.^2)]
    % FR = F + [0; g/2*(hR0.^2 - hR.^2)]
    corrL = g/2 * (hL0.^2 - hL.^2);
    corrR = g/2 * (hR0.^2 - hR.^2);
    FL = F;
    FR = F;
    FL(2,:) = FL(2,:) + corrL;
    FR(2,:) = FR(2,:) + corrR;

    % Compute cell source term S
    % S = [0; -g*(h_m + h_p)/2*(Bp - Bm)/dx] for each cell i=1..N, where
    % h_m = eta(i) - Bm(i+2), h_p = eta(i) - Bp(i+2), Bm,Bp padded aligned with cells
    h_m = eta(3:N+2) - Bm(3:N+2);
    h_p = eta(3:N+2) - Bp(3:N+2);
    S2 = -g * (h_m + h_p)/2 .* (Bp(3:N+2) - Bm(3:N+2))/dx;
    S = zeros(2,N);
    S(2,:) = S2;

  case 'central'
    % Scheme with h,u states reconstructed

    if order == 1
      % Piecewise constant states (h,u) at cells:
      hL = Up(1,2:N+2);
      uL = zeros(1,N+1);
      % face number: k=1..N+1, left state from cell k and right state from cell k+1:
      % For order 1 central, reconstruct h,u constant per cell, so faces left/right states are:
      % left face k left state: U(:,k)
      % right face k right state: U(:,k+1)
      hL = Up(1,2:N+2); % cell centers for left cells
      uL = zeros(1,N+1); % init face values
      hR = Up(1,3:N+3);
      uR = zeros(1,N+1);

      uL = Up(2,2:N+2) ./ (hL + (hL<=tol));
      uR = Up(2,3:N+3) ./ (hR + (hR<=tol));
      uL(hL<=tol) = 0;
      uR(hR<=tol) = 0;
    elseif order == 3
      % order 3, reconstruct h,u with WENO3
      htot = Up(1,:);
      hL = zeros(1,N+1);
      hR = zeros(1,N+1);
      uvec = zeros(1,N+4);
      % compute velocity u = hu/h with u=0 if h<=tol
      uvec(:) = 0;
      wet = (htot > tol);
      uvec(wet) = Up(2,wet) ./ htot(wet);

      [hL, hR] = weno3_faces(htot);
      [uL, uR] = weno3_faces(uvec);
      % u zero where h small (face hL,R)
      uL(hL<=tol) = 0;
      uR(hR<=tol) = 0;
    else
      error('swe1d_rhs: unknown order for central scheme');
    end

    % HLL flux (no hydrostatic reconstruction)
    F = hll_flux(hL,uL,hR,uR,g);
    FL=F; FR=F;

    if order == 1
      % Central scheme uses source: S=[0; -g*h_i*(Bc(i+1)-Bc(i-1))/(2 dx)], central difference of bed at cell centers
      % Bed cell centers (N): P.Bc(3:N+2)
      Bc_interior = Bc(3:N+2);
      dBdx = (Bc(4:N+3) - Bc(2:N+1))/(2*dx); % 1xN central difference
      mass_h = U(1,:);
      S = zeros(2,N);
      S(2,:) = -g * mass_h .* dBdx;
    else
      % central scheme of order 3 is not used: no source term
      S = zeros(2,N);
    end

  otherwise
    error('swe1d_rhs: unknown scheme');
end

% Flux difference and source for interior cells i=1..N:
% dU/dt = -(FL(:,i+1)-FR(:,i))/dx + S(:,i)
dU = zeros(2,N);
dU = -(FL(:,2:N+1) - FR(:,1:N))/dx + S;

end
