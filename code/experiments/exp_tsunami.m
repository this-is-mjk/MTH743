% exp_tsunami: small tsunami over a shelf, crest amplitude against Green's law
set(0,'DefaultFigureVisible','off'); 
s=settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

% PARAMETERS
name = 'tsunami';
% Domain
xl=0; xr=100; N=1000;
g=9.81;
cfl=0.4;
scheme='hr'; order=3; rk='ssprk3';
bc='open';

% Bed setup: depth 10 for x<=30; linear rise to depth 1 at x=90
% bed B(x) negative of depth, bed elevation B= -depth
% make_bed('slope',[30 90 -10 -1]) the b_i params are bed elevation, so negative depths
bed = make_bed('slope',[30 90 -10 -1]);

% Grid setup
dx = (xr-xl)/N;
xc = xl + ((1:N)-0.5)*dx;

% Initial conditions eta, u=0
% eta = 0.05*exp(-((x-15)/3)^2), still water level 0 => h = eta - B(xc)
% Bc = bed at cell centers
Bc = bed.f(xc,0);
eta0 = 0.05*exp(-((xc-15)/3).^2);
init = @(x,Bc_) deal(max(0,eta0 - Bc_), zeros(size(x))); % [h,u]

% Save times
tend=20;
tsave = linspace(0,20,401);

% CONFIG structure
cfg.N=N; cfg.xl=xl; cfg.xr=xr; cfg.bed=bed;
cfg.init=init;
cfg.tend=tend;
cfg.tsave=tsave;
cfg.scheme=scheme; cfg.order=order; cfg.rk=rk; cfg.bc=bc; cfg.cfl=cfl;
cfg.monitor=true;

% RUN simulation
out = swe1d_run(cfg);

x = out.x; B = out.B; % bed at cell centers
H = out.H;    % cell depths h
HU = out.HU;  % hu momenta
nt = numel(H);

%% 1. Figure A: eta(x) at 5 times with filled bed in two stacked panels (full and zoom eta)
figA = figure('Color','w','Position',[100 100 900 600]);

times_plot = [5 10 15 18 20]; % chosen snapshot times (seconds)
[~,ind_plot] = arrayfun(@(t) min(abs(tsave - t)), times_plot);

% Compute eta = h + B at selected times
eta_plot = cellfun(@(h) h + B, H(ind_plot),'UniformOutput',false);

% Plot
c = lines(numel(ind_plot));
lw = 1.8; fs=12; ms=7;

tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact')

% Panel 1: full domain free surface and bed
ax1 = nexttile;
hold on
for i=1:numel(ind_plot)
    plot(x, eta_plot{i}, 'Color', c(i,:), 'LineWidth', lw, 'DisplayName',sprintf('t = %.1f s', tsave(ind_plot(i))));
end
plot(x, B, 'k', 'LineWidth', 1.5, 'DisplayName','Bed elevation');
hold off
xlabel('x [m]','FontSize',fs)
ylabel('\eta [m]','FontSize',fs)
legend('Location','southeast','FontSize',fs)
title('Free surface (amplitude < 5 cm, see zoom below) and bed','FontSize',fs)
grid on

% Panel 2: zoom on free surface eta (y-lims chosen)
ax2 = nexttile;
hold on
for i=1:numel(ind_plot)
    plot(x, eta_plot{i}, 'Color', c(i,:), 'LineWidth', lw);
end
plot(x, B, 'k', 'LineWidth', 1.5);
hold off
xlabel('x [m]','FontSize',fs)
ylabel('\eta [m] (zoom)','FontSize',fs)
ylim([-0.02 0.1])
title('Zoom on \eta(x)','FontSize',fs)
grid on

exportgraphics(figA, fullfile(imgdir,[name '_eta_times.png']), 'Resolution',150)

%% 2. Right-going pulse amplitude vs x and Green's law comparison

% Green's law amplitude scale:
% A(x) = A0 * (h0/h(x))^(1/4)
% h0 = depth in constant region ~ 10 ; at x in [20,30] constant
% find first snapshot where pulse center reaches x=28 in [20,30], find A0 numerically
% Then for snapshots find max eta over x>35 (exclude left pulse) and x location.

% Depth h(x) = eta(x) - B(x), but for Green's law depth is -B(x) (still water)
h0_region = [20 30];
h0_mask = (x>=h0_region(1)) & (x<=h0_region(2));
h0_true = -B(h0_mask);
h0 = mean(h0_true); % mean depth in constant region

% Find time when pulse center is near 28 in [20,30]. Use max eta in [20,30] per time
pulse_center_x = 28;
tol_find_center = 0.2; % m tolerance to consider center reached

% Find max eta and location in [20,30] for each time
max_eta_20_30 = zeros(1,nt);
idx_eta_20_30 = zeros(1,nt);
mask_20_30 = (x>=20) & (x<=30);
for k=1:nt
    [max_eta_20_30(k), idx] = max(H{k}(mask_20_30) + B(mask_20_30));
    idx_eta_20_30(k) = find(mask_20_30,1,'first') - 1 + idx;
end

% Find first time pulse center is near 28
xc_at_max_eta = x(idx_eta_20_30);
dist_center = abs(xc_at_max_eta - pulse_center_x);
time_idx_center = find(dist_center<=tol_find_center,1,'first');
if isempty(time_idx_center)
    error('No snapshot found with pulse center near x=28');
end

A0 = max_eta_20_30(time_idx_center) - 0; % initial amplitude, subtract still water level (0)
% Amplitude baseline is 0, no background eta

% Now for each snapshot find right-going pulse max eta over x>35, location, and amplitude above 0
mask_right = x>35;
measured_amp_all = zeros(1,nt);
measured_x_all = zeros(1,nt);

for k=1:nt
    eta_k = H{k}+B;
    [amp, idxmax] = max(eta_k(mask_right));
    measured_amp_all(k) = amp; % total eta
    xmask = x(mask_right);
    measured_x_all(k) = xmask(idxmax);
end

% Select the crest points of the right-going pulse:
%   start = first k with amp(k) >= 0.9*A0;
%   then keep k = start, start+1, ... while ALL hold: 
%       pos(k) >= pos(k-1) - 1e-9, pos(k) <= 95, amp(k) >= 0.4*A0  
%   (do NOT apply the position-increasing test to k=start).
amp_09A0_mask = measured_amp_all >= 0.9*A0;
start_idx = find(amp_09A0_mask,1,'first');

if isempty(start_idx)
    error('No snapshot found with amplitude >= 0.9*A0 to start tracking.');
end

filtered_indices = start_idx; % initialize with start_idx

for k = start_idx+1:nt
    ampk = measured_amp_all(k);
    posk = measured_x_all(k);
    posk_prev = measured_x_all(k-1);
    if (posk >= posk_prev - 1e-9) && (posk <= 95) && (ampk >= 0.4*A0)
        filtered_indices(end+1) = k; %#ok<AGROW>
    else
        break; % stop at first failure
    end
end

measured_x = measured_x_all(filtered_indices);
measured_amp = measured_amp_all(filtered_indices);

% Check minimum points count
keep_mask_35_80 = (measured_x > 35) & (measured_x < 80);
num_points_for_rms = sum(keep_mask_35_80);

if num_points_for_rms < 10
    error('Fewer than 10 points in filtered data for RMS calculation.');
end

% Ensure monotonicity and sort by position (should already be monotone)
[measured_x, sort_idx] = sort(measured_x);
measured_amp = measured_amp(sort_idx);

% Compute Green's law amplitudes at these positions
hbed = -interp1(x, B, measured_x);
green = A0 * (h0 ./ hbed) .^ (1/4);

% Compute RMS relative deviation for 35 < x < 80 using only filtered points in range
rms_mask = (measured_x > 35) & (measured_x < 80);
rel_diff = (measured_amp(rms_mask) - green(rms_mask)) ./ green(rms_mask);
rms_rel_dev = sqrt(mean(rel_diff.^2));

fprintf('Number of points used for RMS deviation and analysis: %d\n', num_points_for_rms);

% Plot measured and Green's law amplitudes vs position
figB = figure('Color','w','Position',[100 100 900 600]);
plot(measured_x, measured_amp, 'b.-', 'LineWidth', 1.8, 'MarkerSize', 10, 'DisplayName', 'Measured amplitude');
hold on
plot(measured_x, green, 'k-', 'LineWidth', 2, 'DisplayName', 'Green''s law A(x)');
xlabel('x [m]', 'FontSize', fs)
ylabel('right-going crest amplitude [m]', 'FontSize', fs)
xlim([30 95])
title(sprintf('Right-going pulse amplitude vs position and Green''s law (RMS rel dev = %.2f%%)', rms_rel_dev*100), 'FontSize', fs)
legend('Location', 'northwest', 'FontSize', fs)
grid on
hold off

exportgraphics(figB, fullfile(imgdir,[name '_amplitude_comparison.png']), 'Resolution', 150);

% Output table at x = 40, 50, 60, 70, 80, 90 (nearest tracked point)
output_x_targets = [40 50 60 70 80 90];
output_n = numel(output_x_targets);

fid = fopen(fullfile(resdir,[name '_amplitude.txt']),'w');
fprintf(fid,'%% x[m]   Measured[m]   Green[m]   RelDiff[%%]\n');
fprintf(fid,'%% RMS relative deviation (35<x<80): %.3f%%\n', rms_rel_dev*100);

for i = 1:output_n
    xt = output_x_targets(i);
    % Find index of nearest measured_x
    [~, idx_near] = min(abs(measured_x - xt));
    xm = measured_x(idx_near);
    Mm = measured_amp(idx_near);
    Gm = green(idx_near);
    rel_perc = 100*(Mm - Gm)/Gm;
    fprintf(fid, '%6.1f %12.6f %12.6f %12.3f\n', xm, Mm, Gm, rel_perc);
end

fclose(fid);

fprintf('Amplitude data saved with filtered pulse tracking and Green''s law comparison:\n');
fprintf('x[m]   Measured[m]   Green[m]   RelDiff[%%]\n');
for i = 1:output_n
    xt = output_x_targets(i);
    [~, idx_near] = min(abs(measured_x - xt));
    xm = measured_x(idx_near);
    Mm = measured_amp(idx_near);
    Gm = green(idx_near);
    rel_perc = 100*(Mm - Gm)/Gm;
    fprintf('%6.1f %12.6f %12.6f %12.3f\n', xm, Mm, Gm, rel_perc);
end
fprintf('RMS relative deviation (35 < x < 80) = %.3f%%\n', rms_rel_dev*100);

%% 3. x-t contourf of eta with bed profile below

% Prepare eta matrix for all times (nt x N)
eta_mat = zeros(nt,N);
for k=1:nt
    eta_mat(k,:) = H{k} + B;
end

% Figure
figC = figure('Color','w','Position',[100 100 900 600]);

% Use diverging colormap with white at zero
% Center at 0 by shifting color limits:
maxabs_eta = max(abs(eta_mat(:)));
clim = maxabs_eta * [-1 1];

% contourf x time 
contourf(x, tsave, eta_mat, 40, 'LineColor','none')
colormap(gca, redbluecmap);
caxis(clim)
colorbar('Ticks', linspace(clim(1), clim(2),5),...
         'TickLabels',num2cell(linspace(clim(1), clim(2),5),1));

xlabel('x [m]','FontSize',fs)
ylabel('Time [s]','FontSize',fs)
title('x-t contour of free surface elevation \eta with time','FontSize',fs)
hold on

% Bed profile below domain: plot at bottom with scaled height
yy = ylim;
bedscale = 0.3 * maxabs_eta; % bed scaled down
plot(x, B*bedscale + yy(1), 'k-', 'LineWidth', 2); 
hold off
grid on

exportgraphics(figC, fullfile(imgdir,[name '_xt_contour_eta.png']), 'Resolution',150)

%% 4. Lake-at-rest sanity run max |u| reported

% lake at rest: eta=0 everywhere, u=0, same grid
eta0_lake = zeros(size(xc));
init_lake = @(x,Bc_) deal(max(0,eta0_lake-Bc_), zeros(size(x)));
cfg_lake = cfg;
cfg_lake.init = init_lake;
cfg_lake.monitor = true;
cfg_lake.tend = 1; % one second run for sanity
cfg_lake.tsave = 0:0.1:1;

out_lake = swe1d_run(cfg_lake);
% max |u| during simulation from monitor output
max_u_lake = max(out_lake.mon.umax);

% Save and print sanity check max_u
fid = fopen(fullfile(resdir,[name '_lakeatrest_maxu.txt']),'w');
fprintf(fid, 'max|u| during lake-at-rest run: %.2e m/s\n', max_u_lake);
fclose(fid);
fprintf('Lake-at-rest max|u| = %.2e m/s\n', max_u_lake);

%% 5. Create MP4 animation of eta with bed (eta and bed plotted)
figM = figure('Color','w','Position',[100 100 900 600]);
ax = axes(figM);

% Prepare colormap and clim for animation
clim_anim = clim; % from above

% Create video object
mp4name = fullfile(imgdir, 'tsunami_animation.mp4'); % save MP4 in imgdir, NOT resdir
v = VideoWriter(mp4name,'MPEG-4');
v.FrameRate = 20;
open(v);

for k = 1:nt
    cla(ax)
    % plot eta and bed
    plot(ax,x, H{k} + B, 'b-', 'LineWidth',1.8); hold(ax,'on')
    plot(ax,x, B, 'k-', 'LineWidth',1.5)
    ylim([-0.02 0.12])
    xlim([xl xr])
    xlabel('x [m]','FontSize',fs)
    ylabel('\eta [m]','FontSize',fs)
    title(sprintf('Tsunami propagation t = %.2f s', tsave(k)),'FontSize',fs)
    grid(ax,'on')
    legend(ax,{'Free surface \eta','Bed'},'Location','northeast','FontSize',fs)
    drawnow
    frame = getframe(figM);
    writeVideo(v,frame);
end

close(v);
close(figM);

fprintf('Animation saved as %s\n', mp4name);

%% Save all figs + report ready prints
% Text output for main experiment info
txtname = fullfile(resdir,[name '.txt']);
fid = fopen(txtname,'w');
fprintf(fid,'EXP %s RESULTS\n', name);
fprintf(fid,'N=%d, xl=%.1f, xr=%.1f, scheme=%s, order=%d, rk=%s, cfl=%.2f\n', ...
    N, xl, xr, scheme, order, rk, cfl);
fprintf(fid,'Time snapshots saved: %d, tend=%.1f s\n', nt, tend);
fprintf(fid,'Amplitude start A0=%.6f at t=%.2f s (center x=%.2f)\n', A0, tsave(time_idx_center), pulse_center_x);
fprintf(fid,'Lake-at-rest max|u|=%.2e\n', max_u_lake);
fprintf(fid,'RMS relative deviation (35 < x < 80) = %.3f%%\n', rms_rel_dev*100);
fclose(fid);

fprintf('EXP %s DONE\n', name);

%% Local helper: redblue colormap, 64 points (white=zero)
function map = redbluecmap(n)
if nargin<1, n=64; end
half = ceil(n/2);
r = [linspace(0,1,half) ones(1,n-half)];
g = [linspace(0,1,half) linspace(1,0,n-half)];
b = [ones(1,half) linspace(1,0,n-half)];

map = [r(:) g(:) b(:)];
end
