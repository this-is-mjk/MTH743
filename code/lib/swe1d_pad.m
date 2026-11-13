function Up = swe1d_pad(U, bc)
% Pads 2xN state U with 2 ghost cells each side by boundary condition bc
% U: 2 x N, output Up: 2 x (N+4)

[~,N] = size(U);
Up = zeros(2,N+4);
Up(:,3:N+2) = U;

switch bc
  case 'wall'
    % mirror with hu inverted sign
    % Left ghost cells: indexes 2 (first ghost), 1 (second ghost)
    Up(1,2) = U(1,1);
    Up(2,2) = - U(2,1);
    Up(1,1) = U(1,2);
    Up(2,1) = - U(2,2);

    % Right ghost cells: indexes N+3, N+4
    Up(1,N+3) = U(1,N);
    Up(2,N+3) = - U(2,N);
    Up(1,N+4) = U(1,N-1);
    Up(2,N+4) = - U(2,N-1);

  case 'periodic'
    % left ghost copies from right end
    Up(:,1:2) = U(:,N-1:N);
    % right ghost copies from left end
    Up(:,N+3:N+4) = U(:,1:2);

  case 'open'
    % ghost cells copy the adjacent interior cell (edge cell)
    Up(:,1:2) = repmat(U(:,1),1,2);
    Up(:,N+3:N+4) = repmat(U(:,N),1,2);

  otherwise
    error('swe1d_pad: unknown bc');
end

end
