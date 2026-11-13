% exp_p1: Phase 1, the baseline scheme (HLL + central source + Euler) on a lake at rest
function exp_p1
set(0,'DefaultFigureVisible','off');
s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';

% Setup directories
root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

% Common config
cfg.xl = 0; cfg.xr = 10;
cfg.bc = 'wall';
cfg.bed = make_bed('gaussian',[0.4 5 0.5]);
eta0 = 1;
cfg.init = @(x,Bc) deal(max(0,eta0-Bc), zeros(size(x))); % lake at rest
cfg.scheme = 'central';
cfg.order = 1;
cfg.rk = 'euler';
cfg.cfl = 0.5;
cfg.tend = 10;
cfg.monitor = true;

Ns = [50 100 200 400];

maxu_end = zeros(size(Ns));
maxu_time = zeros(size(Ns));
L1_eta_end = zeros(size(Ns));

for iN = 1:numel(Ns)
    cfg.N = Ns(iN);
    out = swe1d_run(cfg);
    % max|u| at end
    u_end = abs(out.HU{end} ./ max(out.H{end},1e-14));
    maxu_end(iN) = max(u_end);
    % max over time of max|u|
    maxu_time(iN) = max(out.mon.umax);
    % L1 error eta vs 1 at end
    dx = (cfg.xr - cfg.xl)/cfg.N;
    eta_end = out.H{end} + out.B;
    L1_eta_end(iN) = sum(abs(eta_end - 1))*dx;
end

% observed order: log2(eN/e2N)
ord_eta = nan(1,numel(Ns)-1);
for i=1:numel(Ns)-1
    if L1_eta_end(i+1)>0 && L1_eta_end(i)>0
        ord_eta(i) = log2(L1_eta_end(i)/L1_eta_end(i+1));
    end
end

% new colors for line plots (replace lines(7) and avoid yellows)
c = [0 0.447 0.741;
     0.85 0.325 0.098;
     0.466 0.674 0.188;
     0.494 0.184 0.556;
     0.301 0.745 0.933;
     0.2 0.2 0.2];

% Plot figure 1: max|u|(t) semilogy for all Ns
fig1 = figure('Color','w','Position',[100 100 900 600]);
hold on
for iN=1:numel(Ns)
    tmon = out_mon_get_tmon(Ns(iN), cfg, imgdir, cfg.bed, eta0); % helper function defined below for reuse
    plot_color = c(mod(iN-1,size(c,1))+1,:);
    semilogy(tmon{1}, tmon{2}, 'Color', plot_color, 'LineWidth',1.5, 'DisplayName',sprintf('N=%d',Ns(iN)));
end
% Add epsilon or reference line in black, LineWidth >=1.5 and in legend 
% (No such line originally, so no additions here)
xlabel('Time t');
ylabel('max|u|(t) (m/s)');
title('Max velocity vs time for scheme ''central'', order 1, Euler RK');
legend('Location','best');
grid on
set(findall(fig1,'-property','FontSize'),'FontSize',11);
exportgraphics(fig1, fullfile(imgdir,'p1_maxu_vs_t.png'),'Resolution',150);

% Figure 2 (N=200): eta(x,t=0,2,10), vs exact = 1 flat, subplot 3 panels stacked
Nref = 200;
cfg.N = Nref;
cfg.monitor = false;
cfg.tsave = [0 2 10];
out = swe1d_run(cfg);
x = out.x; B = out.B; dx = (cfg.xr-cfg.xl)/cfg.N;

fig2 = figure('Color','w','Position',[100 100 900 900]);
tplot = cfg.tsave;
for k=1:3
    eta = out.H{k}+B;
    subplot(3,1,k);
    plot(x, eta, 'Color', c(1,:), 'LineWidth',1.5);
    hold on
    plot(x, ones(size(x)), 'k--', 'LineWidth',1.5);
    title(sprintf('Free surface \\eta(x) at t=%.1f s (N=%d)', tplot(k), Nref));
    xlabel('x');
    ylabel('\eta (m)');
    grid on;
    axis([cfg.xl cfg.xr 0.6 1.6])
    set(gca,'FontSize',11)
end
exportgraphics(fig2, fullfile(imgdir,'p1_eta_times.png'),'Resolution',150);

% Figure 3: same t, u(x)
fig3 = figure('Color','w','Position',[100 100 900 900]);
for k=1:3
    h = out.H{k};
    hu = out.HU{k};
    u = zeros(size(h));
    wet = h>1e-12;
    u(wet) = hu(wet)./h(wet);
    subplot(3,1,k);
    plot(x, u, 'Color', c(2,:), 'LineWidth',1.5);
    title(sprintf('Velocity u(x) at t=%.1f s (N=%d)', tplot(k), Nref));
    xlabel('x');
    ylabel('u (m/s)');
    grid on;
    axis([cfg.xl cfg.xr -0.15 0.15])
    set(gca,'FontSize',11)
end
exportgraphics(fig3, fullfile(imgdir,'p1_u_times.png'),'Resolution',150);

% Figure 4: bed + water at t = 10 for N=200 (reuse figure 2 bottom)? new figure stacked 1 panel:
fig4 = figure('Color','w','Position',[100 100 900 300]);
eta = out.H{end} + out.B;
fill_x = [x(1)-dx/2 x+dx/2 x(end)+dx/2];
bed_fill = [out.B(1) out.B out.B(end)];
water_fill = [eta(1) eta eta(end)];
fill(fill_x, bed_fill, [0.8 0.5 0.2], 'EdgeColor','none'); hold on
fill(fill_x, water_fill, [0 0.3 0.7], 'EdgeColor','none', 'FaceAlpha',0.7);
plot(x, out.B, 'Color',[0.5 0.3 0.1],'LineWidth',2)
plot(x, eta, 'Color', c(1,:), 'LineWidth',1.5);
xlabel('x'); ylabel('Elevation (m)');
title('Bed B(x) and free surface \eta(x) at t=10s');
ylim([min(out.B)-0.1, max(eta)+0.1])
grid on
set(gca,'FontSize',11)
exportgraphics(fig4, fullfile(imgdir,'p1_bed_water.png'),'Resolution',150);

% Animation N=200 over t=[0,10] with tsave=201 frames
cfg.N = 200;
cfg.tsave = linspace(0,10,201);
cfg.monitor = false;
out = swe1d_run(cfg);

figA = figure('Color','w','Position',[100 100 900 600],'Visible','off');
ax1 = subplot(2,1,1);
hplt = plot(x, out.B, 'Color',[0.5 0.3 0.1],'LineWidth',2);
hold on
etaplt = plot(x, out.H{1}+out.B, 'Color', c(1,:), 'LineWidth',1.5);
ylabel('\eta, B (m)')
ylim([min(out.B)-0.2, eta0+0.1])
grid on
title('Free surface \eta and Bed B at time t=0.00s');
set(ax1,'FontSize',11)

ax2 = subplot(2,1,2);
h_u = plot(x, zeros(size(x)), 'Color', c(2,:), 'LineWidth',1.5);
ylabel('u (m/s)')
xlabel('x')
ylim([-0.15 0.15])
grid on
title('Velocity u at time t=0.00s');
set(ax2,'FontSize',11)

% Prepare video writer
vfile = fullfile(imgdir,'p1_animation.mp4');
vw = VideoWriter(vfile,'MPEG-4');
vw.FrameRate = 20;
open(vw);

for k=1:numel(cfg.tsave)
    tcur = cfg.tsave(k);
    eta_k = out.H{k} + out.B;
    u_k = zeros(size(out.H{k}));
    wet = out.H{k} > 1e-12;
    u_k(wet) = out.HU{k}(wet)./out.H{k}(wet);

    % update plots
    set(etaplt,'YData', eta_k)
    set(hplt,'YData', out.B)
    title(ax1, sprintf('Free surface \\eta and Bed B at time t=%.2f s', tcur))

    set(h_u,'YData', u_k)
    title(ax2, sprintf('Velocity u at time t=%.2f s', tcur))

    drawnow
    fr = getframe(figA);
    writeVideo(vw, fr);
end
close(vw);
close(figA);

% Write results table and conclusion
fname = fullfile(resdir,'p1.txt');
fid = fopen(fname,'wt');
fprintf(fid,'N\tmax|u|(end)\tmax|u|(all)\tL1 error eta (end)\torder L1 eta\n');
fprintf('N\tmax|u|(end)\tmax|u|(all)\tL1 error eta (end)\torder L1 eta\n');
for i=1:numel(Ns)
    ord_str = '----';
    if i>1 && ~isnan(ord_eta(i-1))
        ord_str = sprintf('%.2f', ord_eta(i-1));
    end
    fprintf(fid,'%d\t%e\t%e\t%e\t%s\n', Ns(i), maxu_end(i), maxu_time(i), L1_eta_end(i), ord_str);
    fprintf('%d\t%e\t%e\t%e\t%s\n', Ns(i), maxu_end(i), maxu_time(i), L1_eta_end(i), ord_str);
end

% Conclusion from orders
valid_ord = ord_eta(~isnan(ord_eta));
if isempty(valid_ord)
    concl = 'No convergence order computed.';
else
    avg_ord = mean(valid_ord);
    concl = sprintf('Error decreases with N at an observed rate of %.2f (log2 eN/e2N) in L1 eta error.', avg_ord);
end
fprintf(fid,'\nConclusion:\n%s\n', concl);
fprintf('\nConclusion:\n%s\n', concl);
fclose(fid);

fprintf('EXP p1 DONE\n');

end

function mon_out = out_mon_get_tmon(N, cfg_base, imgdir, bed, eta0)
% Helper to get t and max|u| monitor data for figure 1, run with monitor true (cache if needed)
cfg = cfg_base;
cfg.N = N; cfg.monitor = true;
out = swe1d_run(cfg);
mon_out = {out.mon.t, out.mon.umax};
end
