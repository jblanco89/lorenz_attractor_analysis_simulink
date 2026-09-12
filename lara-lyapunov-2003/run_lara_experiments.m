%% RUN_LARA_EXPERIMENTS  Entry point replicating Lara et al. (2003) ENIEF.
%
%   This script reproduces the paper "Estimacion de los exponentes de Lyapunov"
%   (Lara, Stoico, Machado, Castagnino, Mecanica Computacional Vol. XXII,
%   ENIEF 2003). It uses a hybrid design: the Simulink model
%   models/lorenz_sim.slx (fixed-step ode4, 0.001 s) supplies the orbit; pure
%   MATLAB averaging of Jacobian eigenvalues (Sec. 3, Eq. 7/8) supplies the
%   alternative spectrum. The script is modular and incremental: runLorenzSweep,
%   runCpuBenchmark and the viz layer are independent flags.
%
%   Paper -> code mapping (section refs in comments of each called function):
%   - Sec. 3 (Metodo de los autovalores): core/eigenvalueSpectrum.m
%       freezes J(t) -> J_i (Eq. 5), diagonalizes conceptually via eig,
%       averages Re(eig) (Eq. 7/7b, 8a/8b). SIM path: sim/runSimulinkTrajectory.m
%   - Sec. 4 (Resultados numericos, Table 1): experiments/runCpuBenchmark.m
%       pure-MATLAB timing of GS vs eigenvalue method on linear oscillators
%   - Sec. 5 (Ecuacion de Lorenz, Table 2, Figs. 1-3): experiments/runLorenzSweep.m
%       sweeps r=20..30 (here rho), tf=10000, IC (0,1,0), sigma=10, beta=8/3,
%       checks divergence -13.666, classifies attractors (Sec. 3 properties),
%       and viz/plot* functions reproduce Figs. 1a/1b, 2, 3a/3b.
%
%   User decisions applied:
%     1. Table 1 -> pure-MATLAB timing (no Simulink oscillator models)
%     2. Full paper fidelity tf=10000 by default (10M steps at 0.001; slow)
%     3. Gram-Schmidt cross-check on a subset with the existing
%        src/analysis/lorenz_lyapunov_spectrum.m as-is (Tben=500)
%
%   SMOKE vs FULL:
%     This file defaults to paper-fidelity. For quick iteration comment the
%     FULL block and uncomment the SMOKE block below. Smoke uses a short
%     horizon and only the two bracketing rhos (24->25) to verify the
%     fixed-point / strange transition in seconds instead of many minutes.
%
%   Restore invariant: sim/ helpers restore StopTime=100, StopFcn and
%   base-workspace x0=1,y0=1,z0=1, sigma=10,rho=28,beta=8/3 after each batch.

clear; clc; close all;

% --- Resolve paths (so the script works from any cwd) ---
thisFile = mfilename('fullpath');
if isempty(thisFile), thisFile = which('run_lara_experiments'); end
thisDir  = fileparts(thisFile);
projRoot = fileparts(thisDir);
addpath(fullfile(projRoot, 'src'));
addpath(fullfile(projRoot, 'src','build'));
addpath(fullfile(projRoot, 'src','analysis'));
addpath(fullfile(projRoot, 'models'));
addpath(fullfile(thisDir, 'core'));
addpath(fullfile(thisDir, 'sim'));
addpath(fullfile(thisDir, 'experiments'));
addpath(fullfile(thisDir, 'viz'));

% Ensure results output tree exists
resultsDir = fullfile(thisDir, 'results');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

% Verify Simulink model invariants before the long sweep (paper + AGENTS.md)
try
    ensureModelReady();
catch ME
    warning('run_lara_experiments:modelCheck','Model invariant check failed: %s', ME.message);
end

% ============================================================
% CONFIG — toggle SMOKE vs FULL here
% ============================================================

% -- FULL PAPER FIDELITY (tf=10000, 11 rhos, stride 100 -> ~100k eig samples/run) --
DO_SMOKE        = false;  % set true for fast smoke (overrides below if true)
% Common defaults
CFG.rhoRange    = 20:30;
CFG.horizon     = 10000;        % paper Sec. 5 tf
CFG.stride      = 100;          % effective epsilon = 0.1 s
CFG.gsSubset    = [23 24 25 28]; % cross-check subset (user decision 3)
CFG.runSweep    = true;
CFG.runBenchmark= true;
CFG.plotFigs    = true;
CFG.verbose     = true;

% -- SMOKE (uncomment to run fast, e.g. during development/CI) --
% CFG.rhoRange    = [24 25];
% CFG.horizon     = 500;
% CFG.stride      = 200;
% CFG.gsSubset    = [];          % skip GS to stay fast
% CFG.runSweep    = true;
% CFG.runBenchmark= false;
% CFG.plotFigs    = true;
% CFG.verbose     = true;

% If DO_SMOKE flag is set, force smoke values
if DO_SMOKE %#ok<UNRCH>
    fprintf('*** SMOKE MODE (quick verification) ***\n');
    CFG.rhoRange = [24 25];
    CFG.horizon  = 500;
    CFG.stride   = 200;
    CFG.gsSubset = [];
end

% ============================================================
% 1) LORENZ SWEEP (Sec. 5: Table 2 + attractor classes)
% ============================================================
lorenzResults = [];
if CFG.runSweep
    fprintf('\n========== 1) Lorenz sweep (Sec. 5, Table 2) ==========\n');
    lorenzResults = runLorenzSweep( ...
        'RhoRange', CFG.rhoRange, ...
        'Sigma', 10, 'Beta', 8/3, ...
        'Horizon', CFG.horizon, ...
        'Stride', CFG.stride, ...
        'X0', [0 1 0], ...
        'GsSubset', CFG.gsSubset, ...
        'Verbose', CFG.verbose, ...
        'SavePath', fullfile(resultsDir, 'lorenz_sweep.mat'));
end

% ============================================================
% 2) CPU BENCHMARK (Sec. 4: Table 1 — pure MATLAB timing)
% ============================================================
benchResults = [];
if CFG.runBenchmark
    fprintf('\n========== 2) CPU benchmark (Sec. 4, Table 1) ==========\n');
    benchResults = runCpuBenchmark( ...
        'Dims', 2:2:20, ...
        'NIter', 2000, ...
        'Dt', 0.5, ...
        'NRepeats', 3, ...
        'Verbose', CFG.verbose, ...
        'SavePath', fullfile(resultsDir, 'cpu_benchmark.mat'));
end

% ============================================================
% 3) FIGURES (Sec. 5: Figs. 1a/1b, 2, 3a/3b)
% ============================================================
if CFG.plotFigs && ~isempty(lorenzResults)
    fprintf('\n========== 3) Figures (Sec. 5) ==========\n');
    try
        fig1 = plotSpectraVsR(lorenzResults, fullfile(resultsDir, 'fig1_spectra_vs_r.png'));
    catch ME
        warning('run_lara_experiments:fig1','Fig. 1 failed: %s', ME.message);
    end
    try
        % Fig. 2 always available if sweep included r=28 or near
        fig2 = plotMaxExponentConvergence(lorenzResults, 28, fullfile(resultsDir, 'fig2_convergence.png'));
    catch ME
        warning('run_lara_experiments:fig2','Fig. 2 failed: %s', ME.message);
    end
    try
        fig3 = plotLocalEigenvalueTrace(lorenzResults, 28, 23, fullfile(resultsDir, 'fig3_local_eig_traces.png'));
    catch ME
        warning('run_lara_experiments:fig3','Fig. 3 failed: %s', ME.message);
    end
    if verboseExists()
        fprintf('Figures exported to %s\n', resultsDir);
    end
end

fprintf('\nDone. Results in %s\n', resultsDir);

% --- helper ---
function ensureModelReady()
    mdl = 'lorenz_sim';
    % Resolve repo root from this script location (works for scripts and functions)
    scriptPath = which('run_lara_experiments');
    if isempty(scriptPath), scriptPath = mfilename('fullpath'); end
    mdlFile = fullfile(fileparts(fileparts(scriptPath)), 'models', 'lorenz_sim.slx');
    % If file missing, build_lorenz_sim restores all invariants
    if ~exist(mdlFile,'file')
        build_lorenz_sim();
    end
    if ~bdIsLoaded(mdl), load_system(mdl); end
    % Check invariants (AGENTS.md § Model invariants)
    assert(strcmp(get_param(mdl,'Solver'),'ode4'), 'Solver must be ode4');
    assert(strcmp(get_param(mdl,'SolverType'),'Fixed-step'), 'SolverType must be Fixed-step');
    assert(strcmp(get_param(mdl,'FixedStep'),'0.001'), 'FixedStep must be 0.001');
    assert(strcmp(get_param(mdl,'ReturnWorkspaceOutputs'),'off'), 'ReturnWorkspaceOutputs must be off');
    % StopTime and StopFcn are intentionally mutated during sweeps but must be
    % 100 / plot_lorenz_attractor at rest. Restore if drifted.
    if ~strcmp(get_param(mdl,'StopTime'),'100')
        warning('ensureModelReady:stopTime','StopTime is %s (expected 100) — restoring.', get_param(mdl,'StopTime'));
        set_param(mdl,'StopTime','100'); save_system(mdl);
    end
    scf = get_param(mdl,'StopFcn');
    if isempty(scf) || ~strcmp(scf,'plot_lorenz_attractor')
        warning('ensureModelReady:stopFcn','StopFcn is ''%s'' (expected plot_lorenz_attractor) — restoring.', scf);
        set_param(mdl,'StopFcn','plot_lorenz_attractor'); save_system(mdl);
    end
    fprintf('[ensureModelReady] Model invariants OK: ode4 0.001 ReturnWorkspaceOutputs off StopFcn plot_lorenz_attractor\n');
end

function tf = verboseExists()
    tf = true;
end
