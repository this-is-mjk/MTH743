% exp_2d_stability: largest stable CFL of the plain and CTU 2D schemes
set(0,'DefaultFigureVisible','off');
s = settings; s.matlab.appearance.figure.GraphicsTheme.TemporaryValue='light';

root = fullfile(fileparts(mfilename('fullpath')),'..','..');
imgdir = fullfile(root,'outputs','figures');
resdir = fullfile(root,'outputs','results');
if ~isfolder(imgdir), mkdir(imgdir); end
if ~isfolder(resdir), mkdir(resdir); end

% Parameters
Nx = 50; Ny = 50;
xl = 0; xr = 1; yl = 0; yr = 1;
tend = 0.3; 
g = 9.81;
bc = 'wall';

% Bed and init function
bed2 = @(X,Y) zeros(size(X));
init = @(X,Y,B) cat(3, 2*(hypot(X-0.5,Y-0.5)<0.2) + 1.*(hypot(X-0.5,Y-0.5)>=0.2), ...
                 zeros(size(X)), zeros(size(X)));

methods = {'plain','ctu'};
cfls = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.4 1.6 1.8 2.0 2.4];
sweeps_per_step = struct('plain',2,'ctu',4);

% Run baseline reference solutions at CFL=0.3 for each method
h_ref = struct;
for m = 1:numel(methods)
  method = methods{m};
  cfg = struct;
  cfg.Nx = Nx; cfg.Ny = Ny; cfg.xl = xl; cfg.xr = xr; cfg.yl = yl; cfg.yr = yr;
  cfg.bed2 = bed2; cfg.init = init; cfg.tend = tend;
  cfg.bc = bc; 
  cfg.method = char(method); 
  cfg.cfl = 0.3; 
  cfg.g = g;
  try
    out_ref = swe2d_run(cfg);
    h_ref.(method) = out_ref.H{end}; % Nx x Ny at tend
  catch
    error('Reference run failed for method %s at CFL=0.3', method);
  end
end

% Storage for results:
% columns: method (char), cfl, nsteps, cost, max_h, min_h, dev, class (string)
data = {};

for m = 1:numel(methods)
  method = methods{m};
  for ci = 1:numel(cfls)
    cfl = cfls(ci);

    cfg = struct;
    cfg.Nx = Nx; cfg.Ny = Ny; cfg.xl = xl; cfg.xr = xr; cfg.yl = yl; cfg.yr = yr;
    cfg.bed2 = bed2; cfg.init = init; cfg.tend = tend;
    cfg.bc = bc; 
    cfg.method = char(method); 
    cfg.cfl = cfl; 
    cfg.g = g;

    finished = false; nsteps = NaN; maxh = NaN; minh = NaN; dev = NaN; class = "blow-up";
    try
      out = swe2d_run(cfg);
      h = out.H{end};
      if any(~isfinite(h(:)))
        class = "blow-up";
      else
        finished = true;
        nsteps = out.nsteps;
        maxh = max(h(:));
        minh = min(h(:));
        d = max(abs(h(:)-h_ref.(method)(:)));
        dev = d;
        % No classification here now, done later
        class = "unknown";
      end
    catch
      class = "blow-up";
    end
    cost = nsteps * sweeps_per_step.(method);

    data(end+1,:) = {char(method), cfl, nsteps, cost, maxh, minh, dev, class};
    fprintf('Method=%s CFL=%.2f Nsteps=%d cost=%d maxh=%.4f minh=%.4f dev=%.4e class=%s\n',...
      method,cfl,nsteps,cost,maxh,minh,dev,class);
  end
end

% Fix classifications according to new criterion:
% For each method, order by increasing CFL, first run is always stable by design (CFL=0.3)
% For k=2, use only max_h/min_h criteria; for k>=3 use growth and max_h/min_h criteria.
% class = 'blow-up' if failed or non-finite data
% class = 'unstable' if finished and (growth>3 or max_h > 1.4*max_h_0.3 or min_h<0.9*min_h_0.3), and also all later CFL runs finished become unstable
% class = 'stable' otherwise

methods_char = unique(char(data(:,1)),'rows'); % just the unique method names
% We can loop by the cell methods list again for convenience:
for m = 1:numel(methods)
  method = methods{m};
  rows = strcmp(string(data(:,1)), string(method));
  cfl_method = cell2mat(data(rows,2));
  [cfl_sorted, idx_sort] = sort(cfl_method);
  keys = find(rows);
  
  % Reorder indices by sorted CFL for method
  sorted_indices = keys(idx_sort);
  
  % Extract dev, max_h, min_h, class for these runs
  devs = cell2mat(data(sorted_indices,7));
  maxhs = cell2mat(data(sorted_indices,5));
  minhs = cell2mat(data(sorted_indices,6));
  classes = string(data(sorted_indices,8));
  nstepss = cell2mat(data(sorted_indices,3));
  
  % Find base reference values at CFL=0.3 (assumed first run)
  base_idx = find(abs(cfl_sorted - 0.3) < 1e-8,1);
  if isempty(base_idx)
    error('Reference CFL=0.3 run missing for method %s', method);
  end
  maxh_ref = maxhs(base_idx);
  minh_ref = minhs(base_idx);
  
  unstable_found = false;
  blowup_found = false;
  
  for k = 1:numel(cfl_sorted)
    % If blow-up or failed:
    if classes(k) == "blow-up"
      blowup_found = true;
      unstable_found = true;
      data{sorted_indices(k),8} = "blow-up";
    elseif unstable_found && nstepss(k)>0 && isfinite(maxhs(k)) && isfinite(minhs(k)) && isfinite(devs(k))
      % After instability found, all finished runs marked unstable
      data{sorted_indices(k),8} = "unstable";
    elseif classes(k) == "unknown"
      if k == 1
        % First run at CFL=0.3 always stable
        data{sorted_indices(k),8} = "stable";
      elseif k == 2
        % For second run, check only max_h/min_h criteria for instability
        condition_unstable = (maxhs(k) > 1.4*maxh_ref) || ...
                             (minhs(k) < 0.9*minh_ref);
        if condition_unstable
          data{sorted_indices(k),8} = "unstable";
          unstable_found = true;
        else
          data{sorted_indices(k),8} = "stable";
        end
      else
        % For k >= 3, check growth and max_h/min_h criteria
        growth = devs(k)/devs(k-1);
        condition_unstable = (growth > 3) || ...
                             (maxhs(k) > 1.4*maxh_ref) || ...
                             (minhs(k) < 0.9*minh_ref);
        if condition_unstable
          data{sorted_indices(k),8} = "unstable";
          unstable_found = true;
        else
          data{sorted_indices(k),8} = "stable";
        end
      end
    end
  end
end

% Write table file with new classifications
fname = fullfile(resdir,'2d_stability.txt');
fid = fopen(fname,'w');
fprintf(fid,'method cfl nsteps cost max_h min_h dev class\n');
for i=1:size(data,1)
  fprintf(fid,'%s %.2f %d %d %.6f %.6f %.6e %s\n', ...
    data{i,1}, data{i,2}, data{i,3}, data{i,4}, data{i,5}, data{i,6}, data{i,7}, data{i,8});
end
fclose(fid);

% Summaries for each method
summary_lines = strings(numel(methods)*4+1,1);
idx = 1;
largest_stable_CFL = NaN(numel(methods),1);
first_unstable_CFL = NaN(numel(methods),1);
first_blowup_CFL = NaN(numel(methods),1);
cost_at_largest_stable_CFL = NaN(numel(methods),1);

for m = 1:numel(methods)
  method = methods{m};
  rows = strcmp(string(data(:,1)), string(method));
  cfl_method = cell2mat(data(rows,2));
  class_method = string(data(rows,8));

  stable_cfls = cfl_method(class_method == "stable");
  if isempty(stable_cfls)
    largest_stable_cfl = NaN;
  else
    largest_stable_cfl = max(stable_cfls);
  end

  unstable_idx = find(class_method == "unstable",1,'first');
  if isempty(unstable_idx)
    first_unstable_cfl = NaN;
  else
    first_unstable_cfl = cfl_method(unstable_idx);
  end

  blowup_idx = find(class_method == "blow-up",1,'first');
  if isempty(blowup_idx)
    first_blowup_cfl = NaN;
  else
    first_blowup_cfl = cfl_method(blowup_idx);
  end

  % Cost at largest stable CFL
  cost_stable = NaN;
  if ~isnan(largest_stable_cfl)
    indc = find(cfl_method==largest_stable_cfl,1);
    cost_stable = data(rows,:);
    cost_stable = cell2mat(cost_stable(indc,4));
  end

  largest_stable_CFL(m) = largest_stable_cfl;
  first_unstable_CFL(m) = first_unstable_cfl;
  first_blowup_CFL(m) = first_blowup_cfl;
  cost_at_largest_stable_CFL(m) = cost_stable;

  summary_lines(idx) = sprintf('%s largest stable CFL = %.2f', method, largest_stable_cfl);
  idx = idx + 1;
  summary_lines(idx) = sprintf('%s first unstable CFL = %.2f', method, first_unstable_cfl);
  idx = idx + 1;
  summary_lines(idx) = sprintf('%s first blow-up CFL = %.2f', method, first_blowup_cfl);
  idx = idx + 1;
  summary_lines(idx) = sprintf('%s cost at largest stable CFL = %.0f sweeps', method, cost_stable);
  idx = idx + 1;
end

% Ratio of largest stable CFL (ctu/plain)
ctu_idx = find(strcmp(methods,'ctu'),1);
plain_idx = find(strcmp(methods,'plain'),1);
ctu_clean = largest_stable_CFL(ctu_idx);
plain_clean = largest_stable_CFL(plain_idx);
if ~isnan(ctu_clean) && ~isnan(plain_clean)
  ratio = ctu_clean/plain_clean;
  summary_lines(idx) = sprintf('Ratio of largest stable CFL (ctu/plain) = %.3f', ratio);
else
  summary_lines(idx) = "Ratio of largest stable CFL (ctu/plain) = NaN";
end

% Print all summary lines
fid = fopen(fullfile(resdir,'2d_stability_summary.txt'),'w');
for i=1:numel(summary_lines)
  fprintf('%s\n', summary_lines(i));
  fprintf(fid,'%s\n', summary_lines(i));
end
fclose(fid);

% === FIGURE ===
fig=figure('Color','w','Position',[100 100 900 600]);
hold on;
clrs = struct('plain',[0.85 0.325 0.098], 'ctu',[0 0.447 0.741]);

for m = 1:numel(methods)
  method = methods{m};
  rows = strcmp(string(data(:,1)), string(method));
  cfls_m_all = cell2mat(data(rows,2));
  devs_m_all = cell2mat(data(rows,7));
  class_m_all = string(data(rows,8));
  
  % Filter out blow-up points for plotting
  mask = (class_m_all ~= "blow-up") & (devs_m_all > 0);   % the baseline run (dev = 0) cannot be shown on a log axis
  cfls_m = cfls_m_all(mask);
  devs_m = devs_m_all(mask);

  % Fix small dev values for log scale
  devs_m(~isfinite(devs_m)) = 1e-16;
  devs_m(devs_m<1e-16) = 1e-16;

  % semilogy plot
  semilogy(cfls_m, devs_m, '-o', 'Color', clrs.(method), ...
    'DisplayName', sprintf('%s, largest stable CFL=%.2f', method, largest_stable_CFL(m)), ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', clrs.(method));

  % Vertical dashed line at largest stable CFL
  if ~isnan(largest_stable_CFL(m))
    xline(largest_stable_CFL(m), '--', 'Color', clrs.(method), 'LineWidth', 1.5, 'HandleVisibility','off');
  end
end

grid on;
set(gca,'YScale','log');          % hold on from a linear axes keeps the axes linear, so set it explicitly
xlabel('CFL number');
ylabel('max |h - h_{baseline, CFL 0.3}|');
title('Deviation vs CFL number for plain and CTU 2D schemes');
legend('Location','southeast');
set(gca,'FontSize',12);
ylim([5e-3 1]);          % data span 1e-2 ... 7e-1
hold off;

exportgraphics(fig,fullfile(imgdir,'2d_stability.png'),'Resolution',150);

fprintf('EXP 2d_stability DONE\n');
