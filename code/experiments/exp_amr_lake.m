% exp_amr_lake: AMR with 1, 2 and 3 levels actually present: a pulse of size 1e-3 on a lake at rest over a Gaussian bed.
%   Far field (x < 2 or x > 6.5, not yet reached by the wave) must stay at eta = 1, hu = 0 across all coarse-fine interfaces.
root   = fullfile(fileparts(mfilename('fullpath')),'..','..');
resdir = fullfile(root,'outputs','results');
if ~isfolder(resdir), mkdir(resdir); end

bed = make_bed('gaussian',[0.4 5 0.5]); eps0 = 1e-3; tend = 0.3;
init = @(x,Bc) deal(max(0,1+eps0*exp(-((x-4.2)/0.1).^2)-Bc),0*x);
fid = fopen(fullfile(resdir,'amr_lake.txt'),'w');
fprintf(fid,'AMR, N0 = 100, pulse %.0e at x = 4.2 on a lake at rest (Gaussian bed), t = %.1f, thr = 1e-2\n',eps0,tend);
fprintf(fid,'far field: x < 2 or x > 6.5. L1 and updates ratio are against the uniform hr-3 solver on the finest grid.\n');
fprintf(fid,'%5s %-12s %16s %16s %14s %14s\n','Lmax','levels','max|eta-1| far','max|hu| far','L1(h) vs unif','updates ratio');
for L = 1:3
  c = struct('N0',100,'xl',0,'xr',10,'bed',bed,'init',init,'tend',tend,'tsave',tend,'bc','wall', ...
             'Lmax',L,'cfl',0.4,'nregrid',5,'thr',1e-2,'mode','adaptive');
  o = amr1d_run(c); lv = o.level{end}; x = o.x; Nf = numel(x);
  Bf = bed.f(x,0); Bavg = zeros(1,Nf);                      % cell bed of each (nested) cell = mean of the finest-level bed values
  for i = 1:Nf, b = 2^(L-lv(i)); j0 = floor((i-1)/b)*b; Bavg(i) = mean(Bf(j0+1:j0+b)); end
  eta = o.H{end} + Bavg; far = (x < 2.0) | (x > 6.5);
  u = struct('N',Nf,'xl',0,'xr',10,'bed',bed,'init',init,'tend',tend,'tsave',tend,'bc','wall', ...
             'scheme','hr','order',3,'rk','ssprk3','cfl',0.4);
  ru = swe1d_run(u);
  fprintf(fid,'%5d %-12s %16.2e %16.2e %14.2e %14.2f\n',L,mat2str(unique(lv)),max(abs(eta(far)-1)),max(abs(o.HU{end}(far))), ...
          sum(abs(o.H{end}-ru.H{end}))*10/Nf,o.nupdates/ru.nupdates);
end
fclose(fid); fprintf('EXP amr_lake DONE\n');
