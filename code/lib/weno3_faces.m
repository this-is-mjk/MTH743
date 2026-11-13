function [vL,vR] = weno3_faces(v)
  % v: 1 x (N+4) padded vector (2 ghosts each side), output vL,vR: 1 x (N+1)
  % Left-biased from cell k+1: faces 1 to N+1 use cells 0-based indexing shifted:
  % vL(k) = left-biased value at face k from padded cell k+1 = using v(k) to v(k+2)
  % Right-biased from cell k+2: vR(k) from padded cell k+2 = uses v(k+3), v(k+2), v(k+1)
  % (mirror the stencil as specified)

  eps=1e-6;
  N = length(v) - 4;
  vL = zeros(1,N+1);
  vR = zeros(1,N+1);

  % Left biased WENO3 at i+1/2 from cells i-1,i,i+1, i=1..N+1:
  % p0 = -v(i-1)/2 + 3/2 v(i)
  % p1 = 1/2 v(i) + 1/2 v(i+1)
  % b0 = (v(i)-v(i-1))^2
  % b1 = (v(i+1)-v(i))^2
  i = 1:(N+1);
  v_im1 = v(i);     % v(i-1) in formula with offset since indexing shifted: here v(i)=v indexed i, i=1..N+1 corresponds to v index from 1..N+1
  v_i   = v(i+1);
  v_ip1 = v(i+2);

  p0 = (-0.5)*v_im1 + 1.5*v_i;
  p1 = 0.5*v_i + 0.5*v_ip1;
  b0 = (v_i - v_im1).^2;
  b1 = (v_ip1 - v_i).^2;
  a0 = (1/3)./(eps + b0).^2;
  a1 = (2/3)./(eps + b1).^2;
  wsum = a0 + a1;
  vL(:) = (a0 .* p0 + a1 .* p1) ./ wsum;

  % Right biased: mirror image using cells b+1,b,a for face k=1..N+1
  % Indices: a = 2:N+2, b = 3:N+3, b+1=4:N+4, a-1=1:N+1
  % For vR(k), use stencil cells: v(b+1), v(b), v(a)
  % Smoothness indicators and nonlinear weights
  k = 1:(N+1);
  a = k + 1; % 2:N+2
  b = k + 2; % 3:N+3

  v_a   = v(a);
  v_b   = v(b);
  v_bp1 = v(b+1);

  q0 = (-0.5)*v_bp1 + 1.5*v_b;
  q1 = 0.5*v_b + 0.5*v_a;
  c0 = (v_bp1 - v_b).^2;
  c1 = (v_b - v_a).^2;

  alpha0 = (1/3) ./ (eps + c0).^2;
  alpha1 = (2/3) ./ (eps + c1).^2;
  wsum = alpha0 + alpha1;

  vR(:) = (alpha0 .* q0 + alpha1 .* q1) ./ wsum;
end
