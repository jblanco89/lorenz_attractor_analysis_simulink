function fig = plotSpectraVsR(results, savePath)
%PLOTSPECTRAVSR  Replicate paper Figs. 1a / 1b: spectra vs rho (=r).
%
%   fig = plotSpectraVsR(results, savePath)
%
%   WHY these figures (paper Sec. 5, Table 2 + Figs. 1a/1b):
%   As rho=r grows from 20 to 30, the largest exponent crosses zero between
%   24 and 25, marking the transition from fixed point (all Lambda<0) to
%   strange attractor (Lambda1>0). Figs. 1a (Gram-Schmidt) and 1b (eigenvalue
%   method) show the same qualitative jump; Table 2 quantifies it. The two
%   methods are not digit-identical, only sign/description-equivalent.
%
%   This plot shows the eigenvalue-method spectra from the sweep; if subset
%   Gram-Schmidt values are present (results.gsLambda), they are overlaid
%   with a second style for comparison. Divergence -13.666 is annotated as
%   the control (Sec. 5).
%
%   Inputs:
%     results  - struct array from runLorenzSweep
%     savePath - optional char/string to export figure (e.g. 'results/fig1.png')
%
%   Output:
%     fig - figure handle

arguments
    results (:,1) struct
    savePath string = ""
end

rhos = [results.rho];
L = cat(2, results.Lambda);  % 3 x nR (row 1=largest, due to sorted eig in eigenvalueSpectrum)
hasGS = isfield(results,'gsLambda') && any(arrayfun(@(s) all(isfinite(s.gsLambda)), results));
if hasGS
    Lgs = cat(2, results.gsLambda);
end

fig = figure('Name','Spectra vs rho (paper Figs. 1a/1b)','Color','w','Position',[100 100 840 520]);
hold on; grid on; box on;

cols = [0 0.45 0.74; 0.85 0.33 0.10; 0 0.62 0.45]; % blue, red, green
markers = {'o','s','d'};

% Eigenvalue-method lines (solid)
for k = 1:3
    plot(rhos, L(k,:), '-','Color',cols(k,:),'LineWidth',1.8,'Marker',markers{k},'MarkerFaceColor',cols(k,:),'DisplayName',sprintf('Lambda_%d (eig)',k));
end

% Gram-Schmidt subset overlaid (dashed, if available)
if hasGS
    for k = 1:3
        % Only plot finite entries
        ygs = Lgs(k,:);
        ygs(~isfinite(ygs)) = NaN;
        plot(rhos, ygs, '--','Color',cols(k,:)*0.7 + 0.3,'LineWidth',1.4,'Marker',markers{k},'MarkerEdgeColor','k','DisplayName',sprintf('lambda_%d (GS subset)',k));
    end
end

% Transition line between r=24 and r=25 (paper Sec. 5)
xline(24.5,'--','Color',[0.4 0.4 0.4],'LineWidth',1.1,'DisplayName','transition 24->25');
yline(0,'-','Color',[0.6 0.6 0.6],'LineWidth',0.8,'DisplayName','zero');

xlabel('r (= rho)'); ylabel('Lambda (s^{-1})');
title('Lyapunov-like spectra vs r — eigenvalue method (solid) and GS subset (dashed)');
legend('Location','best');
xlim([min(rhos)-0.5, max(rhos)+0.5]);

% Divergence annotation
if all(isfield(results,{'expectedDiv'}))
    expDiv = results(1).expectedDiv;
    text(mean(rhos), min(ylim)+0.05*range(ylim), sprintf('Divergence control: sum Lambda \\approx %.3f (paper Sec. 5)', expDiv), ...
        'FontSize',8,'Color',[0.3 0.3 0.3],'BackgroundColor',[1 1 0.9]);
end

hold off;

if strlength(savePath) > 0
    outDir = fileparts(savePath);
    if strlength(outDir) > 0 && ~exist(outDir,'dir'), mkdir(outDir); end
    exportgraphics(fig, savePath, 'Resolution', 300);
    fprintf('Saved spectra figure to %s\n', savePath);
end
end
