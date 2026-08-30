% SIMULATE_LORENZ_CHAOS
% Runs the lorenz_sim Simulink model twice: nominal and a tiny perturbed
% initial condition, then animates the deterministic chaos (butterfly effect):
%   subplot 1 -> 3D attractor (nominal vs perturbed)
%   subplot 2 -> phase portrait (x-z projection) to see the trajectories split
%   subplot 3 -> ||Delta||(t) on log scale (exponential divergence)
% Requires lorenz_sim.slx in the same directory

% resolve project paths (handles src/ and src/build) and ensure all is on path
thisFile = mfilename('fullpath');
thisDir = fileparts(thisFile);
projRoot = thisDir;
for iter = 1:3
    if exist(fullfile(projRoot, 'models'), 'dir'), break; end
    projRoot = fileparts(projRoot);
end
if ~exist(fullfile(projRoot, 'models'), 'dir'), projRoot = pwd; end
addpath(fullfile(projRoot, 'src'));
addpath(fullfile(projRoot, 'src', 'build'));
addpath(fullfile(projRoot, 'src', 'analysis'));
addpath(fullfile(projRoot, 'src', 'visualization'));
addpath(fullfile(projRoot, 'models'));
mdlFile = fullfile(projRoot, 'models', 'lorenz_sim.slx');
if ~exist(mdlFile, 'file')
    build_lorenz_sim();
end
mdl = 'lorenz_sim';

if ~bdIsLoaded(mdl)
    load_system(mdl);
end

delta = 1e-6;   % tiny perturbation applied to x0
Tend  = 40;

% Ask whether to record the animation as an MP4 (e.g. for an Instagram reel).
% Default is NO, keeping the on-screen animation unchanged. The answer can be
% pre-answered headlessly through the base-workspace variable lorenz_rec_ask
% ('y'/'n'), which avoids blocking in automated/MCP runs.

if evalin('base', 'exist(''lorenz_rec_ask'', ''var'') && ~isempty(lorenz_rec_ask)')
    recChoice = lower(evalin('base', 'lorenz_rec_ask'));
else
    recChoice = lower(strtrim(input('Record animation as Video (MP4 format)? [y/N]: ', 's')));
end
doRecord = ~isempty(recChoice) && recChoice(1) == 'y';
if doRecord
    recFile = fullfile(projRoot, 'media', 'lorenz_chaos_reel.mp4');
    if ~exist(fullfile(projRoot, 'media'), 'dir'), mkdir(fullfile(projRoot, 'media')); end
    fprintf('Recording animation to %s (1080x1920, 30 fps, ~12 s) ...\n', recFile);
end

fprintf('Simulating: nominal (x0=1) and perturbed (x0=1+%g) ...\n', delta);
% NOTE: setVariable/setModelParameter on Simulink.SimulationInput silently fail
% in this setup, so we set the base-workspace variables the model reads and the
% model StopTime directly via set_param.
prevStopFcn = get_param(mdl, 'StopFcn');
if isempty(prevStopFcn), prevStopFcn = 'plot_lorenz_attractor'; end   % never leave it blank
set_param(mdl, 'StopFcn', '');
set_param(mdl, 'StopTime', num2str(Tend));   % run only Tend seconds

assignin('base', 'x0', 1.0); assignin('base', 'y0', 1.0); assignin('base', 'z0', 1.0);
sim(mdl);   % nominal run   -> xsim/ysim/zsim in base
t1 = xsim.Time; x1 = xsim.Data; y1 = ysim.Data; z1 = zsim.Data;

assignin('base', 'x0', 1.0 + delta);   % tiny perturbation on x0
sim(mdl);   % perturbed run -> overwrites xsim/ysim/zsim in base
t2 = xsim.Time; x2 = xsim.Data; y2 = ysim.Data; z2 = zsim.Data;

% restore base workspace + model settings used by standalone runs
assignin('base', 'x0', 1.0);
set_param(mdl, 'StopTime', '100');          % restore full stop time
set_param(mdl, 'StopFcn', prevStopFcn);   % restore for standalone model runs

% both runs produce different solver time grids -> resample onto a common grid
tq = linspace(0, Tend, 5000)';
x1 = interp1(t1, x1, tq); y1 = interp1(t1, y1, tq); z1 = interp1(t1, z1, tq);
x2 = interp1(t2, x2, tq); y2 = interp1(t2, y2, tq); z2 = interp1(t2, z2, tq);
t1 = tq;

d = sqrt((x1-x2).^2 + (y1-y2).^2 + (z1-z2).^2);

% Full Lyapunov spectrum via the Benettin method with Gram-Schmidt
% reorthonormalization of the tangent basis: [lambda_1; lambda_2; lambda_3].
% The method integrates the variational equations along a single orbit for
% 250 s (with tight tolerances), so it is robust against attractor saturation
% and does not depend on a single pair of trajectories.
lambda = lorenz_lyapunov_spectrum();
bar = '+---------------------------------------------+';
fprintf('%s\n', bar);
fprintf('| %-26s | %14s |\n', 'Quantity', 'Value');
fprintf('%s\n', bar);
fprintf('| %-26s | %14.3e |\n', 'Initial distance', d(1));
fprintf('| %-26s | %14.3f |\n', 'Final distance', d(end));
fprintf('| %-26s | %14.3f |\n', 'lambda_1 (Benettin)', lambda(1));
fprintf('| %-26s | %14.3f |\n', 'lambda_2 (Benettin)', lambda(2));
fprintf('| %-26s | %14.3f |\n', 'lambda_3 (Benettin)', lambda(3));
fprintf('| %-26s | %14.3f |\n', 'sum lambda (trace check)', sum(lambda));
fprintf('%s\n', bar);

% --- clean exponential-phase window (after incubation, before saturation)
iA = find(d > 1e-5, 1, 'first');
iB = find(d > 0.8, 1, 'first');
if isempty(iA), iA = 1; end
if isempty(iB), iB = numel(d); end
tA = t1(iA); tB = t1(iB);

X = [x1; x2]; Y = [y1; y2]; Z = [z1; z2];

n = numel(t1);
if doRecord
    nFrames = 360;   % 12 s at 30 fps
else
    nFrames = 180;   % as before
end
idx = round(linspace(1, n, nFrames));

close all;   % release figures left open by previous runs
if doRecord
    % 9:16 portrait that fits on the 1536x864 screen; exportgraphics renders
    % it at 1080x1920 regardless of the on-screen size.
    fig = figure('Name', 'Lorenz deterministic chaos', 'Color', [0.04 0.04 0.07], ...
        'Position', [60 40 450 800]);
else
    fig = figure('Name', 'Lorenz deterministic chaos', 'Color', [0.04 0.04 0.07], ...
        'Position', [200 50 820 920]);
end

cNom = [0.10 0.95 1.00];   % neon cyan    (nominal)
cPer = [1.00 0.22 0.85];   % neon magenta  (perturbed)
cSep = [0.45 1.00 0.55];   % neon green    (separation)
cRef = [1.00 0.92 0.25];   % neon yellow   (reference)
cAx  = [0.66 0.70 0.85];   % axis / label ink
cGrid= [0.18 0.42 0.50];   % faint grid
cBg  = [0.04 0.04 0.07];   % panel background
cInk = [0.90 0.95 1.00];   % title / text ink

% --- subplot 1: 3D attractor (cyberpunk dark theme) ---
ax1 = subplot(3,1,1); hold on; grid on; box on; az0 = -27; el0 = 22; view(ax1, az0, el0);
h1g = animatedline('Color', cNom*0.4, 'LineWidth', 2.5);   % subtle neon halo
h1  = animatedline('Color', cNom, 'LineWidth', 1.1);          % bright core
h2g = animatedline('Color', cPer*0.4, 'LineWidth', 2.5);
h2  = animatedline('Color', cPer, 'LineWidth', 1.1);
xlim([min(X) max(X)]); ylim([min(Y) max(Y)]); zlim([min(Z) max(Z)]);
xlabel('x', 'Color', cAx); ylabel('y', 'Color', cAx); zlabel('z', 'Color', cAx);
title('Lorenz attractor: nominal vs. perturbed trajectories', 'Color', cInk);
legend([h1 h2], {'nominal', 'perturbed'}, 'Location', 'northwest', 'Box', 'off', ...
       'TextColor', cInk, 'Color', [0.06 0.06 0.10]);
set(ax1, 'Color', cBg, 'XColor', cAx, 'YColor', cAx, 'ZColor', cAx, ...
         'GridColor', cGrid, 'GridAlpha', 0.45, 'FontSize', 11, ...
         'TitleFontSizeMultiplier', 1.1);

% --- subplot 2: phase portrait (x-z projection) ---
ax2 = subplot(3,1,2); hold on; grid on; box on;
h3g = animatedline('Color', cNom*0.4, 'LineWidth', 2.5);
h3  = animatedline('Color', cNom, 'LineWidth', 1.0);
h4g = animatedline('Color', cPer*0.4, 'LineWidth', 2.5);
h4  = animatedline('Color', cPer, 'LineWidth', 1.0);
xlim([min([x1;x2]) max([x1;x2])]); ylim([min([z1;z2]) max([z1;z2])]);
xlabel('x', 'Color', cAx); ylabel('z', 'Color', cAx);
title('Phase portrait (x-z): trajectories overlap, then split', 'Color', cInk);
legend([h3 h4], {'nominal', 'perturbed'}, 'Location', 'northwest', 'Box', 'off', ...
       'TextColor', cInk, 'Color', [0.06 0.06 0.10]);
set(ax2, 'Color', cBg, 'XColor', cAx, 'YColor', cAx, ...
         'GridColor', cGrid, 'GridAlpha', 0.45, 'FontSize', 11, ...
         'TitleFontSizeMultiplier', 1.1);

% --- subplot 3: separation (log) - dark theme; keeps log scale + reference + regions ---
ax3 = subplot(3,1,3);
hold on;
hdg = animatedline('Color', cSep*0.4, 'LineWidth', 2.5);
hd  = animatedline('Color', cSep, 'LineWidth', 1.1);
set(ax3, 'YScale', 'log'); grid on;
xlabel('t (s)', 'Color', cAx); ylabel('||\Delta|| (log)', 'Color', cAx);
title(sprintf('Separation between trajectories (Benettin \\lambda_1 \\approx %.2f, %.0f s)', lambda(1), 500), ...
      'Color', cInk);
xlim([0 40]); ylim([1e-6 30]);
% reference line: ideal exponential growth at the Lyapunov rate, anchored at the
% end of the exponential phase (d~1) and projected backwards in time
tref = linspace(tA, tB, 50);
dref = exp(-lambda(1)*(tB - tref));   % = 1 at tB, decays backwards
href = plot(tref, dref, 'Color', cRef, 'LineStyle', '--', 'LineWidth', 1.2);
xline(tA, 'Color', [0.40 0.60 0.70], 'LineStyle', ':', 'LineWidth', 0.8);
xline(tB, 'Color', [0.40 0.60 0.70], 'LineStyle', ':', 'LineWidth', 0.8);
text(mean([0 tA]), 15, 'Incubation', 'HorizontalAlignment', 'center', 'Color', cInk, 'BackgroundColor', [0.08 0.08 0.13]);
text(mean([tA tB]), 15, 'Exponential growth', 'HorizontalAlignment', 'center', 'Color', cInk, 'BackgroundColor', [0.08 0.08 0.13]);
text(mean([tB 40]), 15, 'Saturation', 'HorizontalAlignment', 'center', 'Color', cInk, 'BackgroundColor', [0.08 0.08 0.13]);
legend([hd href], {'||\Delta|| sim', sprintf('rate \\lambda\\approx%.2f', lambda(1))}, ...
       'Location', 'northwest', 'Box', 'off', 'TextColor', cInk, 'Color', [0.06 0.06 0.10]);
set(ax3, 'Color', cBg, 'XColor', cAx, 'YColor', cAx, ...
         'GridColor', cGrid, 'GridAlpha', 0.45, 'FontSize', 11, ...
         'TitleFontSizeMultiplier', 1.1);
hold off;

% --- optional MP4 recording setup ---
if doRecord
    vid = VideoWriter(recFile, 'MPEG-4');
    vid.FrameRate = 30;
    vid.Quality = 75;   % NOTE: ignored by the MPEG-4 encoder (fixed ~10 Mbps -> ~15 MB for 12 s)
    open(vid);
    % guarantee the writer is closed even if the script exits early
    cleanup = onCleanup(@() cleanupRecording(vid));
    % hide the animating figure during capture (faster, avoids the slow
    % on-screen redraw) and show a progress bar instead
    set(fig, 'Visible', 'off');
    wb = waitbar(0, sprintf('Recording %s ... 0/%d', recFile, nFrames), ...
        'Name', 'Lorenz reel', 'CreateCancelBtn', 'setappdata(gcbf,''canceling'',true)');
    setappdata(wb, 'canceling', false);
    cleanupWb = onCleanup(@() cleanupWaitbar(wb));
end

prev = 0;
for k = 1:nFrames
    i = idx(k);
    addpoints(h1g, x1(prev+1:i), y1(prev+1:i), z1(prev+1:i));
    addpoints(h1, x1(prev+1:i), y1(prev+1:i), z1(prev+1:i));
    addpoints(h2g, x2(prev+1:i), y2(prev+1:i), z2(prev+1:i));
    addpoints(h2, x2(prev+1:i), y2(prev+1:i), z2(prev+1:i));
    addpoints(h3g, x1(prev+1:i), z1(prev+1:i));
    addpoints(h3, x1(prev+1:i), z1(prev+1:i));
    addpoints(h4g, x2(prev+1:i), z2(prev+1:i));
    addpoints(h4, x2(prev+1:i), z2(prev+1:i));
    addpoints(hdg, t1(prev+1:i), d(prev+1:i));
    addpoints(hd, t1(prev+1:i), d(prev+1:i));
    view(ax1, az0 + k*0.9, el0);   % slow rotation about z while animating
    prev = i;
    if doRecord
        % capture offscreen and upscale to 1080x1920 (fast: the on-screen
        % figure is small, but display scaling keeps the capture ~1.9x bigger)
        % waitbar gives feedback instead of the slow on-screen animation
        if isgraphics(wb) && getappdata(wb, 'canceling')
            fprintf('Recording canceled at frame %d/%d.\n', k, nFrames);
            break;
        end
        writeVideo(vid, imresize(getframe(fig).cdata, [1920 1080]));
        if isgraphics(wb)
            waitbar(k/nFrames, wb, sprintf('Recording %s ... %d/%d', recFile, k, nFrames));
        end
    else
        drawnow limitrate
        pause(0.01)
    end
end

if doRecord
    if exist('wb','var'), try delete(wb); catch, end; end
    try delete(cleanupWb); catch, end
    set(fig, 'Visible', 'on'); drawnow
    close(vid);
    recInfo = dir(recFile);
    fprintf('Saved %s (%.1f x %.1f px, ~%.1f s, %.2f MB)\n', ...
        recFile, 1080, 1920, nFrames/vid.FrameRate, recInfo.bytes/1e6);
end
% plot_lorenz_attractor(x1, y1, z1, t1);   % 3D Lorenz attractor kept out of this script (StopFcn-only on standalone model runs)
disp('Animation complete.');

% ---------------------------------------------------------------------
function cleanupRecording(vid)
% Close the video writer if the script exits early (error, Ctrl-C).
% No-op if already closed. VideoWriter has no Status property, so just try.
if ~isempty(vid) && isvalid(vid)
    try close(vid); catch, end
end
end
function cleanupWaitbar(wb)
% Close the waitbar if the script exits early (keeps the desktop tidy).
if ~isempty(wb) && isgraphics(wb)
    try delete(wb); catch, end
end
end
