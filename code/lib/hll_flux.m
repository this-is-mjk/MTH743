function F=hll_flux(hL,uL,hR,uR,g)
  tol=1e-10;
  % Same size inputs, output 2 x numel array
  sz = size(hL);
  N=numel(hL);
  hL = reshape(hL,1,N);
  uL = reshape(uL,1,N);
  hR = reshape(hR,1,N);
  uR = reshape(uR,1,N);

  cL = sqrt(g .* hL);
  cR = sqrt(g .* hR);
  sqrt_hL = sqrt(max(hL,0));
  sqrt_hR = sqrt(max(hR,0));
  denom = sqrt_hL + sqrt_hR;
  u_hat = zeros(1,N);
  mask_dryL = (hL <= tol);
  mask_dryR = (hR <= tol);
  % Compute weighted average velocity if both wet:
  both_wet = ~mask_dryL & ~mask_dryR & (denom > 0);
  u_hat(both_wet) = (sqrt_hL(both_wet).*uL(both_wet) + sqrt_hR(both_wet).*uR(both_wet)) ./ denom(both_wet);

  c_hat = sqrt(g .* (hL + hR)/2);

  % Compute wavespeeds sL,sR:
  sL = min(uL - cL, u_hat - c_hat);
  sR = max(uR + cR, u_hat + c_hat);

  % Adjust wavespeeds for dry states:
  % dry left only
  sL(mask_dryL & ~mask_dryR) = uR(mask_dryL & ~mask_dryR) - 2*cR(mask_dryL & ~mask_dryR);
  sR(mask_dryL & ~mask_dryR) = uR(mask_dryL & ~mask_dryR) + cR(mask_dryL & ~mask_dryR);
  % dry right only
  sL(mask_dryR & ~mask_dryL) = uL(mask_dryR & ~mask_dryL) - cL(mask_dryR & ~mask_dryL);
  sR(mask_dryR & ~mask_dryL) = uL(mask_dryR & ~mask_dryL) + 2*cL(mask_dryR & ~mask_dryL);
  % both dry: sL,sR not used, flux=0

  % Physical fluxes on left and right states:
  FL1 = hL .* uL;
  FL2 = hL .* uL.^2 + g .* hL.^2/2;
  FR1 = hR .* uR;
  FR2 = hR .* uR.^2 + g .* hR.^2/2;

  % Reshape to row vectors to vectorize easily
  sLm = min(sL,0);
  sRp = max(sR,0);
  denom_s = sRp - sLm;
  denom_s(denom_s==0) = NaN; % prevent division by zero but handle below

  F1 = (sRp .* FL1 - sLm .* FR1 + sLm .* sRp .* (hR - hL)) ./ denom_s;
  F2 = (sRp .* FL2 - sLm .* FR2 + sLm .* sRp .* (hR .* uR - hL .* uL)) ./ denom_s;

  % Handle dry/dry or denom=0 case: flux=0
  both_dry = (hL <= tol) & (hR <= tol);
  zero_flux_mask = both_dry | isnan(F1);
  F1(zero_flux_mask) = 0;
  F2(zero_flux_mask) = 0;

  F = [F1; F2];
end
