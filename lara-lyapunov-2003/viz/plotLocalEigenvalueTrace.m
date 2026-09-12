function fig = plotLocalEigenvalueTrace(results, rhoChaotic, rhoRegular, savePath)
%PLOTLOCALEIGENVALUETRACE  Replicate paper Figs. 3a / 3b: local eigenvalue traces.
%
%   fig = plotLocalEigenvalueTrace(results, rhoChaotic, rhoRegular, savePath)
%
%   WHY these figures (paper Sec. 5, Figs. 3a/3b and Sec. 6 didactic observation):
%   Near the chaotic regime, at least one eigenvalue of J_i (the frozen Jacobian
%   per window, Sec. 3 Eq. 5) shows strongly oscillatory, aperiodic behavior
%   with variable sign — the time trace alternates irregularly between local
%   stability and instability. The paper illustrates this with the real part
%   of a local eigenvalue at r=28 (chaotic, irregular) vs r=23 (regular,
%   periodic). Sec. 6 elevates this to a didactic signature of chaos:
%   "a lo largo de la trayectoria, la misma debe alternar irregularmente su
%   comportamiento de estabilidad e inestabilidad local."
%
%   This function plots, for two rho values from the sweep, the sampled
%   ReHistory (real parts of eig(J_i) sorted descending) versus tSampled.
%   When a sweep is not available for a rho (e.g. smoke tests), it simulates
%   a short trajectory on demand for that rho.
%
%   Inputs:
%     results     - struct array from runLorenzSweep (may be missing one of the two rhos)
%     rhoChaotic  - rho for chaotic trace (default 28)
%     rhoRegular  - rho for regular trace (default 23)
%     savePath    - optional export path (single figure with two subplots)
%
%   Output:
%     fig - figure handle (contains 2 subplots)

arguments
    results (:,1) struct
    rhoChaotic (1,1) double = 28
    rhoRegular (1,1) double = 23
    savePath string = ""
end

% Helper to fetch or synthesize ReHistory for a given rho
    function rec = fetchOrSimulate(rhoVal)
        idx = find([results.rho] == rhoVal, 1);
        if ~isempty(idx)
            rec = results(idx);
            return;
        end
        % On-demand short simulation (so Figs. 3a/3b still work in smoke mode)
        fprintf('[plotLocalEigenvalueTrace] rho=%g not in sweep — simulating short orbit (horizon 2000)...\n', rhoVal);
        addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'lara-lyapunov-2003','core'));
        addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'lara-lyapunov-2003','sim'));
        % Use paper IC and sigma/beta
        [t, X] = runSimulinkTrajectory(rhoVal, 2000, 10, 8/3, 0, 1, 0);
        eigRes = eigenvalueSpectrum(t, X, 10, rhoVal, 8/3, 50);
        rec = struct('rho',rhoVal,'tSampled',eigRes.tSampled,'ReHistory',eigRes.ReHistory,'Lambda',eigRes.Lambda);
    end

recC = fetchOrSimulate(rhoChaotic);
recR = fetchOrSimulate(rhoRegular);

fig = figure('Name', sprintf('Local Re(eig) traces: r=%g vs r=%g (Figs. 3a/3b)', rhoChaotic, rhoRegular), ...
    'Color','w','Position',[100 80 960 420]);

% Panel 1: chaotic (Fig. 3a)
ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
% Show all three real parts (paper shows one representative; showing three makes
% the sign-variable oscillation unmistakable). Sorted descending -> first trace
% is the most unstable direction.
for k = 1:size(recC.ReHistory,1)
    plot(ax1, recC.tSampled, recC.ReHistory(k,:), 'LineWidth', 1.1);
end
yline(ax1, 0, '--','Color',[0.5 0.5 0.5]);
xlabel(ax1,'t (s)'); ylabel(ax1,'Re(\lambda_{local}) (s^{-1})');
title(ax1, sprintf('Fig. 3a  r=%g (chaotic) — aperiodic, sign-variable', recC.rho));

% Panel 2: regular (Fig. 3b)
ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
for k = 1:size(recR.ReHistory,1)
    plot(ax2, recR.tSampled, recR.ReHistory(k,:), 'LineWidth', 1.1);
end
yline(ax2, 0, '--','Color',[0.5 0.5 0.5]);
xlabel(ax2,'t (s)'); ylabel(ax2,'Re(\lambda_{local}) (s^{-1})');
title(ax2, sprintf('Fig. 3b  r=%g (regular) — periodic', recR.rho));

sgtitle('Local Jacobian eigenvalue traces — didactic chaotic signature (paper Sec. 6)');

if strlength(savePath) > 0
    outDir = fileparts(savePath);
    if strlength(outDir) > 0 && ~exist(outDir,'dir'), mkdir(outDir); end
    exportgraphics(fig, savePath, 'Resolution', 300);
    fprintf('Saved local-eigenvalue traces to %s\n', savePath);
end
end
