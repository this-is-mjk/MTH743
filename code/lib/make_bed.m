function bed=make_bed(type,params)
  % Return struct with fields f, type, params
  bed.type=type;
  bed.params=params;
  switch type
    case 'flat'
      bed.f=@(x,~) zeros(size(x));
    case 'gaussian'
      a=params(1); x0=params(2); s=params(3);
      bed.f=@(x,~) a*exp(-(x - x0).^2/(2*s^2));
    case 'cosbump'
      a=params(1); x0=params(2); w=params(3);
      bed.f=@(x,~) (abs(x - x0)<w).*(a*(cos(pi*(x - x0)/w)+1)/2);
    case 'step'
      B0=params(1); B1=params(2); xs=params(3);
      bed.f=@(x,side) step_bed(x,side,B0,B1,xs);
    case 'sin2'
      bed.f=@(x,~) sin(pi*x).^2;
    case 'slope'
      x1=params(1); x2=params(2); b1=params(3); b2=params(4);
      bed.f=@(x,~) slope_bed(x,x1,x2,b1,b2);
    otherwise
      error('make_bed: unknown bed type');
  end
end

function B=step_bed(x,side,B0,B1,xs)
  % side: -1 left limit, +1 right limit; only matters at x==xs
  B=zeros(size(x));
  leftmask = (x < xs) | ((x == xs) & (side==-1));
  rightmask = (x > xs) | ((x == xs) & (side==+1));
  B(leftmask)=B0;
  B(rightmask)=B1;
end

function B=slope_bed(x,x1,x2,b1,b2)
  B = zeros(size(x));
  idx1 = (x <= x1);
  idx2 = (x >= x2);
  idxlin = (~idx1) & (~idx2);
  B(idx1) = b1;
  B(idx2) = b2;
  % linear interpolation for x in (x1,x2)
  if any(idxlin)
    B(idxlin) = b1 + (b2 - b1).*(x(idxlin)-x1)./(x2 - x1);
  end
end
