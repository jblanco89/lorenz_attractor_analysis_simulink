function fig = plotMaxExponentConvergence(results, rhoTarget, savePath)
%PLOTMAXEXPONENTCONVERGENCE  Replicate paper Fig. 2: Lambda1 convergence in time.
%
%   fig = plotMaxExponentConvergence(results, rhoTarget, savePath)
%
%   WHY this figure (paper Sec. 5, Fig. 2):
%   As integration proceeds at fixed rho, the largest coefficient converges to
%   a quasi-stationary value (paper Sec. 3, Eq. 7b/8b: limit as n->inf / t->inf).
%   Fig. 2 overlays the Gram-Schmidt max exponent and the eigenvalue-average
%   max at rho=28, showing that both stabilize but at slightly different
%   values — sign/description equivalent, not digit-identical (Sec. 5, Sec. 6).
%   The time axis runs to tf=10000 in the paper; here it runs to the sweep
%   horizon (paper fidelity 10000, smoke horizons shorter still show the trend).
%
%   This function plots, for the rho closest to rhoTarget, the cumulative
%   Lambda_1(t) from eigenvalueSpectrum (LambdaHistory(1,:)) vs tSampled.
%
%   Inputs:
%     results   - struct array from runLorenzSweep
%     rhoTarget - scalar rho to plot (default 28)
%     savePath  - optional export path
%
%   Output:
%     fig - figure handle

arguments
    results (:,1) struct
    rhoTarget (1,1) double = 28
    savePath string = ""
end

% Find closest rho entry
[~, idx] = min(abs([results.rho] - rhoTarget));
rec = results(idx);

tSamp = rec.tSampled;
hist = rec.LambdaHistory;  % 3 x nSamples
if isempty(tSamp) || isempty(hist)
    error('plotMaxExponentConvergence:missingHistory','No LambdaHistory at rho=%g (was stride too large or sweep missing?)', rec.rho);
end

fig = figure('Name', sprintf('Lambda1 convergence at r=%g (Fig. 2)', rec.rho), 'Color','w','Position',[120 120 760 420]);
hold on; grid on; box on;

plot(tSamp, hist(1,:), 'LineWidth',1.8, 'Color',[0 0.45 0.74], 'DisplayName','Lambda_1 (eig avg)');
if size(hist,1) >= 2
    plot(tSamp, hist(2,:), 'LineWidth',1.2, 'Color',[0.85 0.33 0.10], 'LineStyle','--', 'DisplayName','Lambda_2');
    plot(tSamp, hist(3,:), 'LineWidth',1.2, 'Color',[0 0.62 0.45], 'LineStyle','--', 'DisplayName','Lambda_3');
end

% Horizontal asymptote at final value (quasi-stationary limit, Sec. 3)
yline(hist(1,end),'--','Color',[0.4 0.4 0.4],'DisplayName',sprintf('Lambda_1(\\infty)\\approx%.3f',hist(1,end)));

xlabel('t (s)'); ylabel('Lambda (s^{-1})');
title(sprintf('Convergence of spectrum coefficients at r=%.0f (paper Fig. 2, horizon %.0f s)', rec.rho, tSamp(end)));
legend('Location','best');
hold off;

if strlength(savePath) > 0
    outDir = fileparts(savePath);
    if strlength(outDir) > 0 && ~exist(outDir,'dir'), mkdir(outDir); end
    exportgraphics(fig, savePath, 'Resolution', 300);
    fprintf('Saved convergence figure to %s\n', savePath);
end
end
