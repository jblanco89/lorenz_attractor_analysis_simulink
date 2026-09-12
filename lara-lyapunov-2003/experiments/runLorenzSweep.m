function results = runLorenzSweep(varargin)
%RUNLORENZSWEEP  Replicate paper Sec. 5 (Table 2, Figs. 1a/1b, 2, 3a/3b).
%
%   results = runLorenzSweep('Name',Value,...)
%
%   WHY this sweep (paper Sec. 5, Eq. 5.1 and Table 2):
%   The Lorenz equation possesses a rich phase space; for sigma=10, beta=8/3
%   and rho=r swept from 20 to 30 (step 1) with IC (0,1,0) and tf=10000, the
%   system transitions from a stable fixed point (three negative exponents,
%   r<=24) to a strange attractor (one positive exponent, rest negative,
%   negative divergence, r>=25). The paper reports Table 2 with two rows per r:
%   Gram-Schmidt (standard/Benettin) and the proposed eigenvalue-average method
%   (Sec. 3). Both give the same sign/description and the same divergence
%   -13.666, confirming the alternative spectrum.
%
%   This function replicates the eigenvalue half of Table 2 via the
%   simulink_lorenz model (hybrid Simulink + MATLAB), validates via the
%   divergence control (§3 property, §5 "variable de control"), classifies
%   each r, and optionally cross-checks a subset against the existing
%   Gram-Schmidt implementation src/analysis/lorenz_lyapunov_spectrum.m
%   (called as-is on that subset, with sigma/rho/beta set via assignin, per
%   user decision 3).
%
%   Name-Value options:
%     'RhoRange'  - vector of rho (=r) values (default 20:30)
%     'Sigma'     - scalar sigma (default 10)
%     'Beta'      - scalar beta  (default 8/3)
%     'Horizon'   - StopTime tf in s (default 10000, paper fidelity)
%     'Stride'    - eigenvalue subsampling factor (default 100 => 0.1 s)
%     'X0'        - 1x3 IC row [x0 y0 z0] (default [0 1 0] per paper §5)
%     'GsSubset'  - subset of RhoRange on which to run Gram-Schmidt
%                  cross-check (default [23 24 25 28]). Use [] to skip.
%     'Verbose'   - logical, print per-r progress (default true)
%     'SavePath'  - char/string path to save results .mat (default ''=no save)
%
%   Output:
%     results - struct array with fields per rho:
%       rho, Lambda (3x1), divergence, expectedDiv, attractorLabel,
%       tSampled, LambdaHistory, ReHistory, gsLambda (3x1 or NaN if not in subset)
%
%   Example:
%     % Quick smoke: two r values bracketing the transition (seconds, not hours)
%     res = runLorenzSweep('RhoRange',[24 25],'Horizon',500,'Stride',200);
%     % Full paper: 11 runs at tf=10000 (slow: ~10M steps each)
%     res = runLorenzSweep('SavePath','lara-lyapunov-2003/results/lorenz_sweep.mat');

% Parse options
p = inputParser;
addParameter(p, 'RhoRange', 20:30);
addParameter(p, 'Sigma', 10);
addParameter(p, 'Beta', 8/3);
addParameter(p, 'Horizon', 10000);
addParameter(p, 'Stride', 100);
addParameter(p, 'X0', [0 1 0]);
addParameter(p, 'GsSubset', [23 24 25 28]);
addParameter(p, 'Verbose', true);
addParameter(p, 'SavePath', "");
parse(p, varargin{:});
opts = p.Results;

rhoVals   = opts.RhoRange;
sigmaVal  = opts.Sigma;
betaVal   = opts.Beta;
horizon   = opts.Horizon;
stride    = opts.Stride;
x0Vec     = opts.X0;
gsSubset  = opts.GsSubset;
verbose   = opts.Verbose;
savePath  = string(opts.SavePath);

% Ensure core/sim on path (when called as script)
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
projRoot = fileparts(thisDir); projRoot = fileparts(projRoot);
addpath(fullfile(projRoot, 'lara-lyapunov-2003','core'));
addpath(fullfile(projRoot, 'lara-lyapunov-2003','sim'));
addpath(fullfile(projRoot, 'src','analysis'));

expectedDiv = lorenzDivergence(sigmaVal, betaVal);
if verbose
    fprintf('=== Lorenz sweep (paper Sec. 5) === rho=%s horizon=%.0f stride=%d ===\n', ...
        mat2str(rhoVals), horizon, stride);
    fprintf('Expected divergence sum Lambda = %.6f\n', expectedDiv);
end

% Preallocate results
nR = numel(rhoVals);
results = struct('rho',{},'Lambda',{},'divergence',{},'expectedDiv',{}, ...
    'attractorLabel',{},'tSampled',{},'LambdaHistory',{},'ReHistory',{}, ...
    'gsLambda',{},'nSamples',{});

for i = 1:nR
    rhoVal = rhoVals(i);
    if verbose
        fprintf('[%d/%d] rho=%g -> simulating (tf=%.0f)...\n', i, nR, rhoVal, horizon);
    end

    % Hybrid: Simulink provides the orbit; MATLAB averages eigenvalues (§3)
    [t, X] = runSimulinkTrajectory(rhoVal, horizon, sigmaVal, betaVal, x0Vec(1), x0Vec(2), x0Vec(3));

    eigRes = eigenvalueSpectrum(t, X, sigmaVal, rhoVal, betaVal, stride);

    Lambda = eigRes.Lambda;
    divVal = eigRes.divergence;

    % Paper §3 + §5 classification
    label = classifyAttractor(Lambda, divVal);

    % Divergence control (Sec. 5): warn if drift beyond tolerance
    divErr = abs(divVal - expectedDiv);
    if divErr > 0.05
        warning('runLorenzSweep:divergenceDrift', ...
            'rho=%g divergence %.4f deviates %.4f from expected %.4f (check solver/stride)', ...
            rhoVal, divVal, divErr, expectedDiv);
    end

    % Optional Gram-Schmidt cross-check on subset (user decision 3: use existing
    % lorenz_lyapunov_spectrum as-is, Tben=500, ode45, dT=0.5, 25% transient).
    gsLambda = nan(3,1);
    if ismember(rhoVal, gsSubset)
        if verbose, fprintf('  -> Gram-Schmidt cross-check (subset, Tben=500)...\n'); end
        % The spectrum function reads sigma/rho/beta from base workspace
        assignin('base','sigma', sigmaVal);
        assignin('base','rho',   rhoVal);
        assignin('base','beta',  betaVal);
        try
            gsLambda = lorenz_lyapunov_spectrum();
        catch ME
            warning('runLorenzSweep:gsFailed','GS failed at rho=%g: %s', rhoVal, ME.message);
        end
        % Paper only claims sign/class equivalence, not digit equality
        if verbose && all(isfinite(gsLambda))
            gsLabel = classifyAttractor(gsLambda, sum(gsLambda));
            match = (label == gsLabel);
            fprintf('  Lambda_eig=[%.4f %.4f %.4f] (%s)  GS=[%.4f %.4f %.4f] (%s)  classMatch=%d\n', ...
                Lambda(1),Lambda(2),Lambda(3), label, gsLambda(1),gsLambda(2),gsLambda(3), gsLabel, match);
            if ~match
                warning('runLorenzSweep:classMismatch','Class mismatch at rho=%g: eig=%s GS=%s', ...
                    rhoVal, label, gsLabel);
            end
        end
    end

    results(i).rho = rhoVal;
    results(i).Lambda = Lambda;
    results(i).divergence = divVal;
    results(i).expectedDiv = expectedDiv;
    results(i).attractorLabel = label;
    results(i).tSampled = eigRes.tSampled;
    results(i).LambdaHistory = eigRes.LambdaHistory;
    results(i).ReHistory = eigRes.ReHistory;
    results(i).gsLambda = gsLambda;
    results(i).nSamples = eigRes.nSamples;

    if verbose
        fprintf('  Lambda=[%.5f %.5f %.5f] sum=%.4f divErr=%.2e class=%s n=%d\n', ...
            Lambda(1),Lambda(2),Lambda(3), divVal, divErr, label, eigRes.nSamples);
    end

    % Free large orbit arrays before next rho iteration
    clear t X eigRes
end

% Print Table 2-like summary (Sec. 5)
if verbose
    printTable2(results);
end

% Optionally save
if strlength(savePath) > 0
    outDir = fileparts(savePath);
    if strlength(outDir) > 0 && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    save(savePath, 'results', 'opts');
    if verbose, fprintf('Saved sweep results to %s\n', savePath); end
end
end

function printTable2(results)
% Pretty-print a Table 2 style summary (eigenvalue method only;
% Gram-Schmidt column shown where subset was evaluated).
fprintf('\n--- Replicated Table 2 (eigenvalue method, Sec. 5) ---\n');
fprintf('%4s | %10s %10s %10s | %10s | %-12s | %s\n', ...
    'r','Lambda1','Lambda2','Lambda3','sumLambda','Class','GS (subset)');
fprintf('%s\n', repmat('-',1,82));
for k = 1:numel(results)
    r = results(k);
    if all(isfinite(r.gsLambda))
        gsStr = sprintf('[%.3f %.3f %.3f]', r.gsLambda(1),r.gsLambda(2),r.gsLambda(3));
    else
        gsStr = '-';
    end
    fprintf('%4g | %10.5f %10.5f %10.5f | %10.4f | %-12s | %s\n', ...
        r.rho, r.Lambda(1), r.Lambda(2), r.Lambda(3), r.divergence, r.attractorLabel, gsStr);
end
fprintf('%s\n\n', repmat('-',1,82));
end
