% exp_2d: Phase 4, 2D unsplit (plain) versus CTU: lake at rest, circular dam break, isotropy, plane wave
set(0,'DefaultFigureVisible','off'); % headless for scripts and tests
s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

g = 9.81;

%% 1. Lake at rest in 2D
Nx = 60; Ny = 60;
xl = 0; xr = 1; yl = 0; yr = 1;

bedfn = @(X,Y) 0.4*exp(-((X-0.5).^2+(Y-0.5).^2)/0.02);

eta0 = 1;

initlake = @(X,Y,B) cat(3,max(0,eta0-B), zeros(size(X)), zeros(size(X)));

% timesave only at tend, to save memory
tsave1 = 2;

cfg1plain = struct('Nx',Nx,'Ny',Ny,'xl',xl,'xr',xr,'yl',yl,'yr',yr, ...
  'bed2',bedfn,'init',initlake, 'tend',tsave1, 'tsave',tsave1, ...
  'method','plain','bc','wall','cfl',0.45,'g',g);
cfg1ctu = cfg1plain; cfg1ctu.method = 'ctu'; cfg1ctu.cfl = 0.9;

out1plain = swe2d_run(cfg1plain);
out1ctu = swe2d_run(cfg1ctu);

hplain = out1plain.H{end};
huplain = out1plain.HU{end};
hvplain = out1plain.HV{end};
etaplain = hplain + out1plain.B;
maxhuplain = max(abs(huplain(:)));
maxhvplain = max(abs(hvplain(:)));
maxetaplain = max(abs(etaplain(:)-eta0));

hctu = out1ctu.H{end};
huctu = out1ctu.HU{end};
hvctu = out1ctu.HV{end};
etactu = hctu + out1ctu.B;
maxhuctu = max(abs(huctu(:)));
maxhvctu = max(abs(hvctu(:)));
maxetactu = max(abs(etactu(:) - eta0));

% Table 1 output
t1_names = {'max|hu|','max|hv|','max|eta - 1|'};
vals_plain = [maxhuplain maxhvplain maxetaplain];
vals_ctu = [maxhuctu maxhvctu maxetactu];

fid = fopen(fullfile(resdir,'2d.txt'),'w');
fprintf(fid,'EXP 2D: Lake at rest max|hu|, max|hv|, max|eta-1|\n');
fprintf(fid,'Method   | max|hu|      max|hv|      max|eta-1|\n');
fprintf(fid,'---------|-------------------------------\n');
fprintf(fid,'plain    | %12.5e %12.5e %12.5e\n',vals_plain);
fprintf(fid,'ctu      | %12.5e %12.5e %12.5e\n',vals_ctu);


%% 2. Circular dam break on flat bed

bedflat = @(X,Y) zeros(size(X));
center = 0.5;
r = @(X,Y) sqrt((X-center).^2+(Y-center).^2);

h0 = @(X,Y) double(r(X,Y)<0.2)*2 + double(r(X,Y)>=0.2)*1;
% Initial velocity zero everywhere
initdam = @(X,Y,B) cat(3, max(0,h0(X,Y)-B), zeros(size(X)), zeros(size(X)));

tend2 = 0.1;

% Configurations to run:
reslist = [50 100];
methlist = {'plain','ctu'};

% High-res reference 400x400 ctu
cfg_ref = struct('Nx',400,'Ny',400,'xl',0,'xr',1,'yl',0,'yr',1, ...
  'bed2',bedflat,'init',initdam,'tend',tend2,'tsave',tend2, ...
  'method','ctu','bc','wall','cfl',0.9,'g',g);
outref = swe2d_run(cfg_ref);

% Coarse ref by exact block averaging (no interp) for 100x100
Nrf = 400;
Nrc = 100;
ratio = Nrf/Nrc;
assert(ratio==4);
href = outref.H{end};
% Block averaging over 4x4 blocks
hx = reshape(href, ratio, Nrc, ratio, Nrc);
h_ref100 = squeeze(mean(mean(hx,1),3));
x_ref = outref.x(1:ratio:end);
y_ref = outref.y(1:ratio:end);
[X_ref,Y_ref] = ndgrid(x_ref,y_ref);

L1_err = zeros(length(methlist),length(reslist));

for im = 1:length(methlist)
  method_i = methlist{im};
  for ir = 1:length(reslist)
    N = reslist(ir);
    dx = 1/N;
    cfg_i = struct('Nx',N,'Ny',N,'xl',0,'xr',1,'yl',0,'yr',1, ...
      'bed2',bedflat,'init',initdam,'tend',tend2,'tsave',tend2, ...
      'method',method_i,'bc','wall','cfl', 0.45*(method_i=="plain")+0.9*(method_i=="ctu"),'g',g);
    outi = swe2d_run(cfg_i);
    h_i = outi.H{end};
    % Project reference to coarse grid shape by exact block averaging:
    % use exact block average from href 400x400 to N x N coarse
    blk = Nrf/N;
    hx = reshape(href, blk, N, blk, N);
    hR = squeeze(mean(mean(hx,1),3));
    % Compute L1 error by summation (no interp)
    L1_err(im,ir) = sum(abs(h_i - hR), 'all')*dx^2;

    % Save solutions for plotting for res=100 only for fig 1, 2
    if N==100 && ismember(method_i,{'plain','ctu'})
      if method_i=="plain"
        hplain100 = h_i;
      else
        hctu100 = h_i;
      end
    end
  end
end


% Write table for 2
fprintf(fid,'\nCircular dam break L1 error (vs 400x400 ctu ref)\n');
fprintf(fid,'Method/Res | L1 error   50x50      100x100\n');
fprintf(fid,'-----------|-----------------------------\n');
for im=1:length(methlist)
  fprintf(fid,'%9s | %10.3e %10.3e\n',methlist{im},...
    L1_err(im,1),L1_err(im,2));
end

% ---- Isotropy test at t = 0.05, BEFORE the shock reaches any wall ------------------------------------------------
% The shock front is at r ~ 0.2 + 4*t ~ 0.4 < 0.5 (distance to the nearest wall). At t = 0.1 the front has already
% been reflected on the axes (walls at r = 0.5) but not on the diagonals (corner at r = 0.71), which would show up as a
% fake anisotropy. The metric is evaluated only for r <= 0.45.
tiso = 0.05; riso = 0.45;
cfg_ri = cfg_ref; cfg_ri.tend = tiso; cfg_ri.tsave = tiso;
href_iso = swe2d_run(cfg_ri).H{end};
aniso = zeros(length(methlist),length(reslist));
for im = 1:length(methlist)
  for ir = 1:length(reslist)
    N = reslist(ir);
    cfg_i = struct('Nx',N,'Ny',N,'xl',0,'xr',1,'yl',0,'yr',1, ...
      'bed2',bedflat,'init',initdam,'tend',tiso,'tsave',tiso, ...
      'method',methlist{im},'bc','wall','cfl', 0.45*strcmp(methlist{im},'plain')+0.9*strcmp(methlist{im},'ctu'),'g',g);
    h_iso = swe2d_run(cfg_i).H{end};
    aniso(im,ir) = anisotropy_metric(h_iso, N, riso);
    if N==100, hiso100.(methlist{im}) = h_iso; end
  end
end
aniso_ref = anisotropy_metric(href_iso, 400, riso);
fprintf(fid,'\nIsotropy of the circular dam break at t = %.2f (shock still inside the domain), r <= %.2f:\n',tiso,riso);
fprintf(fid,'mean |h_diag(r) - h_axis(r)| / (max h - min h); axis profile interpolated to the diagonal radii; 0 = isotropic\n');
fprintf(fid,'Method/Res |   50x50     100x100\n');
for im=1:length(methlist)
  fprintf(fid,'%9s | %9.3e %9.3e\n',methlist{im},aniso(im,1),aniso(im,2));
end
fprintf(fid,'Reference 400x400 CTU: %9.3e\n',aniso_ref);
fprintf(fid,'Ratio aniso(plain)/aniso(ctu): 50x50 %.2f, 100x100 %.2f\n',aniso(1,1)/aniso(2,1),aniso(1,2)/aniso(2,2));

% Create coordinates for 100x100
xc100 = ((1:100)-0.5)/100;

f1 = figure('Color','w','Position',[100 100 900 600]);
subplot(1,3,1);
surf(xc100, xc100, hplain100'); shading interp; view(35,40); colorbar;
title('plain 100x100 h at t=0.1');
zlim([0.9 2.1]);

subplot(1,3,2);
surf(xc100, xc100, hctu100'); shading interp; view(35,40); colorbar;
title('ctu 100x100 h at t=0.1');
zlim([0.9 2.1]);

subplot(1,3,3);
surf(xc100, xc100, h_ref100'); shading interp; view(35,40); colorbar;
title('reference 400x400 ctu avg 100x100');
zlim([0.9 2.1]);
exportgraphics(f1, fullfile(imgdir,'2d_circ_dam_break_surfs.png'), 'Resolution',150);
close(f1);

% Radial profiles figure - axis and diagonal profiles of plain, ctu, reference vs r
f2 = figure('Color','w','Position',[100 100 900 600]);
hold on;

% Compute axis and diagonal profiles for reference 400x400 
[h_axis_ref, r_axis_ref, h_diag_ref, r_diag_ref] = extract_profiles(href_iso, 400);

plot(r_axis_ref, h_axis_ref, 'Color',[0 0 0], 'LineStyle','-', 'LineWidth', 2, 'DisplayName', 'Reference axis');
plot(r_diag_ref, h_diag_ref, 'Color',[0 0 0], 'LineStyle','--', 'LineWidth', 2, 'DisplayName', 'Reference diag');

% Compute axis and diagonal profiles for plain 100x100
[h_axis_plain, r_axis_plain, h_diag_plain, r_diag_plain] = extract_profiles(hiso100.plain, 100);
plot(r_axis_plain, h_axis_plain, 'Color',[0.85 0.325 0.098], 'LineStyle','-', 'LineWidth', 1.8, 'DisplayName', 'Plain axis');
plot(r_diag_plain, h_diag_plain, 'Color',[0.85 0.325 0.098], 'LineStyle','--', 'LineWidth', 1.8, 'DisplayName', 'Plain diag');

% Compute axis and diagonal profiles for ctu 100x100
[h_axis_ctu, r_axis_ctu, h_diag_ctu, r_diag_ctu] = extract_profiles(hiso100.ctu, 100);
plot(r_axis_ctu, h_axis_ctu, 'Color',[0 0.447 0.741], 'LineStyle','-', 'LineWidth', 1.8, 'DisplayName', 'CTU axis');
plot(r_diag_ctu, h_diag_ctu, 'Color',[0 0.447 0.741], 'LineStyle','--', 'LineWidth', 1.8, 'DisplayName', 'CTU diag');

xlabel('Radius r');
ylabel('Mean h(r)');
legend('Location','NorthEastOutside');
xlim([0 0.5]);
title('h(r) at t = 0.05 along an axis (solid) and a diagonal (dashed): isotropy test');

exportgraphics(f2, fullfile(imgdir,'2d_circ_dam_break_radprof.png'), 'Resolution',150);
close(f2);

%% 3. Plane wave at 45deg

Nxw = 100; Nyw = 100;
xlw = 0; xrw = 1; ylw = 0; yrw = 1;

bedflat0 = @(X,Y) zeros(size(X));
% centered along diagonal s = (X+Y)/sqrt(2)
s = @(X,Y) (X+Y)/sqrt(2);
w = 0.08;
eta_init = @(X,Y,B) 1 + 0.05*exp(-((s(X,Y) - 0.5)/w).^2);

% velocity zero:
initwave = @(X,Y,B) cat(3,max(0,eta_init(X,Y,B)-B), zeros(size(X)), zeros(size(X)));

tendw = 0.12;

cfgplaneplain = struct('Nx',Nxw,'Ny',Nyw,'xl',xlw,'xr',xrw,'yl',ylw,'yr',yrw, ...
  'bed2',bedflat0, 'init',initwave,'tend',tendw, 'tsave', tendw, ...
  'method','plain', 'bc', 'open', 'cfl',0.45 ,'g',g);

cfgplanectu = cfgplaneplain; cfgplanectu.method = 'ctu'; cfgplanectu.cfl = 0.9;

outwref = swe2d_run(struct('Nx',400,'Ny',400,'xl',xlw,'xr',xrw,'yl',ylw,'yr',yrw,...
  'bed2',bedflat0, 'init',initwave,'tend',tendw,'tsave', tendw, ...
  'method','ctu','bc','open','cfl',0.9,'g',g));

outwplain = swe2d_run(cfgplaneplain);
outwctu = swe2d_run(cfgplanectu);

Hplain = outwplain.H{end} + outwplain.B;
Hctu = outwctu.H{end} + outwctu.B;
Href = outwref.H{end} + outwref.B;
xw = outwplain.x; yw = outwplain.y;
[Xw,Yw] = ndgrid(xw,yw);

% Sample along diagonal s = (X+Y)/sqrt(2)
sd = (Xw + Yw)/sqrt(2);
sd_lin = linspace(0, sqrt(2), Nxw);

% Instead of interpolation from 400x400, use exact block averaging from Href onto 100x100 first:
% Compute block averaging ratio 4 just like before
blk = 400/100; % 4
hxw = reshape(Href, blk, 100, blk, 100);
H_ref100 = squeeze(mean(mean(hxw,1),3));

% Sample diagonal indices i=j in coarse 100x100 grid
diag_inds = 1:100;

eta_diag_plain = diag(Hplain(diag_inds,diag_inds));
eta_diag_ctu = diag(Hctu(diag_inds,diag_inds));
eta_diag_ref = diag(H_ref100(diag_inds,diag_inds));

% Peak amp and FWHM (full width at half max)
ampplain = max(eta_diag_plain) - 1;
ampctu = max(eta_diag_ctu) - 1;
ampref = max(eta_diag_ref) - 1;

% Background 1, so FWHM at value 1 + amp/2
halfval_plain = 1 + ampplain/2;
halfval_ctu = 1 + ampctu/2;
halfval_ref = 1 + ampref/2;

% Using the local function fwhm_width for better estimate:

fwhm_plain = fwhm_width(eta_diag_plain, sd_lin);
fwhm_ctu = fwhm_width(eta_diag_ctu, sd_lin);
fwhm_ref = fwhm_width(eta_diag_ref, sd_lin);

fprintf(fid,'\nPlane wave at 45deg: max amplitude and FWHM (along diagonal s)\n');
fprintf(fid,'Method | max eta-1      | FWHM (diagonal direction)\n');
fprintf(fid,'-------|----------------|-------------------------\n');
fprintf(fid,'plain  | %12.5e | %12.5e\n', ampplain, fwhm_plain);
fprintf(fid,'ctu    | %12.5e | %12.5e\n', ampctu, fwhm_ctu);
fprintf(fid,'ref400 | %12.5e | %12.5e\n', ampref, fwhm_ref);

% Plot diagonal profiles
fig3 = figure('Color','w','Position',[100 100 900 600]);
plot(sd_lin, eta_diag_plain, 'Color',[0.85 0.325 0.098], 'LineStyle','-', 'LineWidth',1.8, 'DisplayName','Plain 100x100');
hold on
plot(sd_lin, eta_diag_ctu, 'Color',[0 0.447 0.741], 'LineStyle','-', 'LineWidth',1.8, 'DisplayName','CTU 100x100');
plot(sd_lin, eta_diag_ref, 'Color',[0 0 0], 'LineStyle','-', 'LineWidth',2, 'DisplayName','Reference 400x400 CTU');
xlabel('Diagonal coordinate s (m)');
ylabel('\eta');
title('Plane wave at 45°: Free surface \eta along diagonal');
legend('Location','Best');
exportgraphics(fig3, fullfile(imgdir,'2d_plane_wave_diag.png'), 'Resolution',150);
close(fig3);


%% 5. Make MP4 video 2d_dam_break.mp4: ctu 100x100, tsave 0:0.005:0.3, wall BC

tsavev = linspace(0,0.3,61);
Nxv = 100; Nyv = 100;

cfgv = struct('Nx',Nxv,'Ny',Nyv,'xl',0,'xr',1,'yl',0,'yr',1, 'bed2',bedflat,...
  'init',initdam,'tend',0.3,'tsave',tsavev,'method','ctu','bc','wall','cfl',0.9,'g',g);

outv = swe2d_run(cfgv);

figv = figure('Color','w','Position',[100 100 900 600]);
vw = VideoWriter(fullfile(imgdir,'2d_dam_break.mp4'),'MPEG-4');
open(vw);
axis tight manual;
ax = gca;
ax.ZLim = [0.9 1.9];
view(35,40);
colorbar;

% Coordinates for 100x100
xv = ((1:100)-0.5)/100;
yv = xv;
[Xv,Yv] = ndgrid(xv,yv);

for k = 1:length(outv.t)
  clf(figv);
  surf(Xv,Yv,outv.H{k}'); 
  shading interp; 
  zlim([0.9 1.9]);
  title(sprintf('CTU dam break 100x100, t = %.3f s', outv.t(k)));
  colorbar;
  view(35,40);
  drawnow;
  frame = getframe(figv);
  writeVideo(vw,frame);
end
close(vw);
close(figv);

fclose(fid);

fprintf('EXP 2D DONE\n');

%% Local functions defined after script code

function a = anisotropy_metric(h, N, rmax)
% Isotropy of a radially symmetric field: mean |h_diag(r) - h_axis(r)| over r, relative to the range of h.
  [h_axis, r_axis, h_diag, r_diag] = extract_profiles(h, N);
  ha = interp1(r_axis, h_axis, r_diag, 'linear');
  ok = ~isnan(ha) & (r_diag <= rmax);
  d = abs(h_diag(ok) - ha(ok));
  a = trapz(r_diag(ok), d) / ((r_diag(find(ok,1,'last')) - r_diag(find(ok,1,'first'))) * (max(h(:)) - min(h(:))));
end

function [h_axis, r_axis, h_diag, r_diag] = extract_profiles(h, N)
% Axis and diagonal profiles of h(r) for plotting (column vectors)
  % Cell center coords:
  xc = ((1:N)-0.5)/N;
  [X,Y] = ndgrid(xc,xc);
  ic = N/2; % center corner between cells ic, ic+1 in each direction

  kvals = (1:(N/2-1))';
  % Axis ray right direction
  h_axis = 0.5*(h(ic+kvals, ic) + h(ic+kvals, ic+1));
  r_axis = 0.5*( hypot(X(ic+kvals, ic)-0.5, Y(ic+kvals, ic)-0.5) + ...
                 hypot(X(ic+kvals, ic+1)-0.5, Y(ic+kvals, ic+1)-0.5) );

  % Diagonal ray up-right direction:
  idx = sub2ind(size(h), ic+kvals, ic+kvals);     % true diagonal cells (not the whole sub-matrix)
  h_diag = h(idx);
  r_diag = hypot(X(idx)-0.5, Y(idx)-0.5);
end

function w = fwhm_width(v, s)
% Compute Full Width at Half Maximum of profile v sampled at positions s (monotonically increasing)
% v: vector values of profile
% s: vector positions corresponding to v
% Returns: width at half maximum using linear interpolation.

minv = min(v);
maxv = max(v);
half = (maxv + minv)/2;

% Find indices where v crosses half level
ind1 = find(v >= half, 1, 'first');
ind2 = find(v >= half, 1, 'last');

if isempty(ind1) || isempty(ind2)
    % No half crossings found
    w = NaN;
    return;
end

% Handle edge cases for ind1 and ind2 interpolation
if ind1 == 1
    x_lower = s(1);
else
    x1 = s(ind1-1); x2 = s(ind1);
    y1 = v(ind1-1); y2 = v(ind1);
    % Linear interpolation to find crossing
    x_lower = x1 + (half - y1) * (x2 - x1) / (y2 - y1);
end

if ind2 == length(v)
    x_upper = s(end);
else
    x3 = s(ind2); x4 = s(ind2+1);
    y3 = v(ind2); y4 = v(ind2+1);
    % Linear interpolation to find crossing
    x_upper = x3 + (half - y3) * (x4 - x3) / (y4 - y3);
end

w = x_upper - x_lower;
end
