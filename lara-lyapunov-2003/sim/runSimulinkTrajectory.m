function [t, X] = runSimulinkTrajectory(rhoVal, horizon, sigmaVal, betaVal, x0Val, y0Val, z0Val)
%RUNSIMULINKTRAJECTORY  Simulate models/lorenz_sim.slx with controlled parameters.
%
%   [t, X] = runSimulinkTrajectory(rhoVal, horizon, sigmaVal, betaVal, x0Val, y0Val, z0Val)
%
%   Wraps the Simulink model that is the Lorenz system (paper Sec. 5, Eq. 5.1):
%       dx/dt = sigma*(y - x)  (see models/lorenz_sim.slx / src/build/build_lorenz_sim.m)
%       dy/dt = x*(rho - z) - y   (paper notation: rho == r)
%       dz/dt = x*y - beta*z
%   Paper Sec. 5 parameters: sigma=10, beta=8/3, rho=r in 20..30, IC (0,1,0),
%   tf=10000, divergence -(1+sigma+beta) as control.
%
%   WHY this wrapper (paper + repo gotchas, see AGENTS.md):
%   - Paper Sec. 5 studies the Lorenz orbit at long horizons (tf=10000). The repo
%     default is StopTime=100 for the chaos demo. This wrapper temporarily sets
%     StopTime to the requested horizon, as the paper did.
%   - Repo invariant: Integrator ICs are base-workspace variables x0/y0/z0.
%     Gotcha: SimulationInput.setVariable silently fails in this MCP setup, so
%     we use assignin('base', ...) + set_param(StopTime) before sim(model),
%     exactly like src/simulate_lorenz_chaos.m does.
%   - The model StopFcn='plot_lorenz_attractor' auto-plots on every sim and would
%     spawn 11 figures during a sweep. This wrapper blanks StopFcn during the
%     batch (and restores it after) — same pattern as simulate_lorenz_chaos.
%   - Model is script-generated (build_lorenz_sim.m). If missing, it is rebuilt.
%     Solver stays Fixed-step ode4, FixedStep 0.001, ReturnWorkspaceOutputs off
%     (timeseries xsim/ysim/zsim in base workspace) — never changed here.
%
%   Inputs (all optional, defaults = paper Sec. 5 sweep or repo defaults):
%     rhoVal  - rho (= r) value to assign to base variable 'rho' (default 28)
%     horizon - StopTime in seconds (default 100, paper sweep uses 10000)
%     sigmaVal,betaVal - sigma, beta (defaults 10, 8/3)
%     x0Val,y0Val,z0Val - initial conditions (paper sweep: 0,1,0; repo demo: 1,1,1)
%
%   Outputs:
%     t - N x 1 time vector (from xsim.Time)
%     X - N x 3 matrix [x y z] (from xsim.Data, ysim.Data, zsim.Data)

arguments
    rhoVal  (1,1) double = 28
    horizon (1,1) double = 100
    sigmaVal (1,1) double = 10
    betaVal  (1,1) double = 8/3
    x0Val (1,1) double = 0
    y0Val (1,1) double = 1
    z0Val (1,1) double = 0
end

% Resolve project root (same traversal as src/simulate_lorenz_chaos.m)
thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);            % .../lara-lyapunov-2003/sim
projRoot = fileparts(thisDir);             % .../lara-lyapunov-2003
projRoot = fileparts(projRoot);            % repo root
if ~exist(fullfile(projRoot, 'models'), 'dir')
    p = fileparts(mfilename('fullpath'));
    for iter = 1:4
        if exist(fullfile(p, 'models'), 'dir')
            projRoot = p; break;
        end
        p = fileparts(p);
    end
end
if ~exist(fullfile(projRoot, 'models'), 'dir')
    projRoot = pwd;
end

addpath(fullfile(projRoot, 'src'));
addpath(fullfile(projRoot, 'src', 'build'));
addpath(fullfile(projRoot, 'src', 'analysis'));
addpath(fullfile(projRoot, 'models'));

mdl = 'lorenz_sim';
mdlFile = fullfile(projRoot, 'models', [mdl '.slx']);

if ~exist(mdlFile, 'file')
    fprintf('[runSimulinkTrajectory] Model file missing — rebuilding via build_lorenz_sim()...\n');
    build_lorenz_sim();
end

if ~bdIsLoaded(mdl)
    load_system(mdl);
end

% Snapshot model state (must restore even on Ctrl-C / error)
prevStopTime = get_param(mdl, 'StopTime');
prevStopFcn  = get_param(mdl, 'StopFcn');
prevSaveTime = get_param(mdl, 'SaveTime');
prevSaveOutput = get_param(mdl, 'SaveOutput');
if isempty(prevStopFcn)
    prevStopFcn = 'plot_lorenz_attractor';  % AGENTS.md recovery
end

% Snapshot base workspace scalars
prevSigma = 10; prevRho = 28; prevBeta = 8/3;
prevX0 = 1; prevY0 = 1; prevZ0 = 1;
try prevSigma = evalin('base','sigma'); catch, end
try prevRho   = evalin('base','rho');   catch, end
try prevBeta  = evalin('base','beta');  catch, end
try prevX0    = evalin('base','x0');    catch, end
try prevY0    = evalin('base','y0');    catch, end
try prevZ0    = evalin('base','z0');    catch, end

% Guarantee restoration via onCleanup (captures snapshot variables by closure)
cleanupObj = onCleanup(@() restoreAll());

    function restoreAll()
        try set_param(mdl,'StopTime', prevStopTime); catch, end
        try set_param(mdl,'StopFcn',  prevStopFcn);  catch, end
        try set_param(mdl,'SaveTime', prevSaveTime); catch, end
        try set_param(mdl,'SaveOutput', prevSaveOutput); catch, end
        try assignin('base','sigma', prevSigma); catch, end
        try assignin('base','rho',   prevRho);   catch, end
        try assignin('base','beta',  prevBeta);  catch, end
        try assignin('base','x0',    prevX0);    catch, end
        try assignin('base','y0',    prevY0);    catch, end
        try assignin('base','z0',    prevZ0);    catch, end
    end

% Assign requested parameters to base workspace
assignin('base', 'sigma', sigmaVal);
assignin('base', 'rho',   rhoVal);
assignin('base', 'beta',  betaVal);
assignin('base', 'x0',    x0Val);
assignin('base', 'y0',    y0Val);
assignin('base', 'z0',    z0Val);

% Override horizon and silence StopFcn for batch
% Also disable tout/yout logging to avoid workspace conflict when called from a
% function (AGENTS.md notes StopFcn handling; this avoids "Cannot create tout")
set_param(mdl, 'StopFcn', '');
set_param(mdl, 'StopTime', num2str(horizon));
set_param(mdl, 'SaveTime', 'off');
set_param(mdl, 'SaveOutput', 'off');

% Run Simulink: To Workspace blocks log timeseries xsim/ysim/zsim to base.
% Use evalin('base') so the logging lands in the base workspace even when
% this helper is called from a function (repo gotcha: To Workspace targets base).
evalin('base', sprintf('sim(''%s'');', mdl));

% Retrieve logs
try
    xsim = evalin('base', 'xsim');
    ysim = evalin('base', 'ysim');
    zsim = evalin('base', 'zsim');
catch ME
    error('runSimulinkTrajectory:missingLog', ...
        'Expected timeseries xsim/ysim/zsim not found after sim: %s', ME.message);
end

t = xsim.Time;
X = [xsim.Data, ysim.Data, zsim.Data];

% Explicit restore (onCleanup also guards interruption)
restoreAll();

% Dismiss the onCleanup so it does not fire again at function exit
% (we already restored; clearing it avoids double-restore warning)
try delete(cleanupObj); catch, end
end
