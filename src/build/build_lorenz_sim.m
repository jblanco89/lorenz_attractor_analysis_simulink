function build_lorenz_sim()
%BUILD_LORENZ_SIM  Build a Simulink model of the Lorenz attractor
%using simple blocks (integrators, sums, gains, products).
%   dx/dt = sigma*(y - x)
%   dy/dt = x*(rho - z) - y
%   dz/dt = x*y - beta*z
%with sigma=10, rho=28, beta=8/3 (chaotic regime).
%
%Reproduces the model exactly as manually edited:
%  - block positions/sizes snap to the same coordinates
%  - orthogonal (90-degree) line routing identical to the edited model
%    (all 32 line segments replayed as point polylines)
%  - fixed-step RK4 (ode4) solver with 0.001 s step, StopTime 100
%  - color palette + TeX annotations (via decorate_lorenz_sim)
%
%Integrator initial conditions are taken from workspace variables
%x0, y0, z0 so they can be varied per-simulation via
%Simulink.SimulationInput.setVariable to study the butterfly effect.

% ensure project folders are on path (handles src/build location)
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
mdl = 'lorenz_sim';
mdlFile = fullfile(projRoot, 'models', [mdl '.slx']);

% --- reset state: close if loaded, delete stale file ---
if bdIsLoaded(mdl)
    close_system(mdl, 0);
end
if exist(mdlFile, 'file')
    delete(mdlFile);
end

% --- parameters and initial conditions in base workspace ---
sigma = 10; rho = 28; beta = 8/3;
x0 = 1; y0 = 1; z0 = 1;
assignin('base', 'sigma', sigma);
assignin('base', 'rho',   rho);
assignin('base', 'beta',  beta);
assignin('base', 'x0',    x0);
assignin('base', 'y0',    y0);
assignin('base', 'z0',    z0);

% --- create and open the model ---
new_system(mdl);
open_system(mdl);

% --- integrators (state = x, y, z) ---
add_block('simulink/Continuous/Integrator', [mdl '/Integrator x'], ...
    'InitialCondition', 'x0', 'Position', [405 185 445 225]);
add_block('simulink/Continuous/Integrator', [mdl '/Integrator y'], ...
    'InitialCondition', 'y0', 'Position', [405  75 445 115]);
add_block('simulink/Continuous/Integrator', [mdl '/Integrator z'], ...
    'InitialCondition', 'z0', 'Position', [405 305 445 345]);

% --- dx/dt = sigma*(y - x) ---
add_block('simulink/Math Operations/Sum', [mdl '/Sum dx'], ...
    'Inputs', '+-', 'Position', [170 195 190 215]);
add_block('simulink/Math Operations/Gain', [mdl '/Gain sigma'], ...
    'Gain', 'sigma', 'Position', [240 190 280 220]);

% --- dy/dt = x*(rho - z) - y ---
add_block('simulink/Sources/Constant', [mdl '/Const rho'], ...
    'Value', 'rho', 'Position', [85 30 135 60]);
add_block('simulink/Math Operations/Sum', [mdl '/Sum rhoz'], ...
    'Inputs', '+-', 'Position', [160 70 180 90]);
add_block('simulink/Math Operations/Product', [mdl '/Product x(rho-z)'], ...
    'Inputs', '**', 'Position', [260 57 300 88]);
add_block('simulink/Math Operations/Sum', [mdl '/Sum dy'], ...
    'Inputs', '+-', 'Position', [330 85 350 105]);

% --- dz/dt = x*y - beta*z ---
add_block('simulink/Math Operations/Product', [mdl '/Product xy'], ...
    'Inputs', '**', 'Position', [130 277 170 308]);
add_block('simulink/Math Operations/Gain', [mdl '/Gain beta'], ...
    'Gain', 'beta', 'Position', [260 355 300 385]);
add_block('simulink/Math Operations/Sum', [mdl '/Sum dz'], ...
    'Inputs', '+-', 'Position', [330 315 350 335]);

% --- sinks: log x,y,z as timeseries + scope for live viewing ---
add_block('simulink/Sinks/To Workspace', [mdl '/To Workspace x'], ...
    'VariableName', 'xsim', 'SaveFormat', 'Timeseries', 'Position', [540 185 600 215]);
add_block('simulink/Sinks/To Workspace', [mdl '/To Workspace y'], ...
    'VariableName', 'ysim', 'SaveFormat', 'Timeseries', 'Position', [540 80 600 110]);
add_block('simulink/Sinks/To Workspace', [mdl '/To Workspace z'], ...
    'VariableName', 'zsim', 'SaveFormat', 'Timeseries', 'Position', [540 285 600 315]);
add_block('simulink/Sinks/Scope', [mdl '/Scope'], ...
    'NumInputPorts', '3', 'Position', [535 385 605 435]);

% --- wiring: replay the exact orthogonal (90-degree) line routing ---
drawExactLines(mdl);

% --- solver: fixed-step RK4 (ode4) with 0.001 s step ---
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode4', ...
    'FixedStep', '0.001', 'StopTime', '100', 'StartTime', '0.0');

% Write trajectories to the base workspace (xsim/ysim/zsim) instead of a single
% SimulationOutput object, so the model StopFcn can read and plot them.
set_param(mdl, 'ReturnWorkspaceOutputs', 'off');

% --- StopFcn: plot the Lorenz attractor automatically when a simulation ends ---
set_param(mdl, 'StopFcn', 'plot_lorenz_attractor');

save_system(mdl, mdlFile);

decorate_lorenz_sim();   % colors + title/equations TeX annotations

fprintf('Model %s built: %d blocks, %d lines (RK4, step 0.001, t=0..100)\n', ...
    mdl, numel(find_system(mdl, 'LookUnderMasks', 'all', 'Type', 'block')), ...
    numel(get_param(mdl, 'Lines')));
fprintf('Saved to %s\n', mdlFile);
end

% ---------------------------------------------------------------------
function drawExactLines(mdl)
% Replays every line segment captured from the edited model as a point
% polyline via add_line(sys, points). Shared endpoints auto-merge into
% continuous polylines; segments starting on an interior point become
% branches, reproducing the manual 90-degree routing exactly.
%
% Draw order matters: trunks/singles first (source port start), then
% continuation segments, then branches (start on a previously drawn line).
TRUNKS = {
    [355 325;400 325], ...            % Sum dz   -> Integrator z
    [305 370;340 370;340 340], ...    % Gain beta-> Sum dz
    [175 295;340 295;340 310], ...    % Product xy -> Sum dz
    [355 95;400 95], ...              % Sum dy   -> Integrator y
    [305 75;340 75;340 80], ...       % Product x(rho-z) -> Sum dy
    [185 80;255 80], ...              % Sum rhoz -> Product x(rho-z)
    [450 325;475 325], ...            % Integrator z  trunk
    [140 45;170 45;170 65], ...       % Const rho -> Sum rhoz
    [285 205;400 205], ...            % Gain sigma -> Integrator x
    [195 205;235 205], ...            % Sum dx   -> Gain sigma
    [450 205;475 205;475 200], ...    % Integrator x  trunk
    [450 95;475 95], ...              % Integrator y  trunk
    };
CONT = {
    [475 325;475 330], ...            % Integrator z  cont
    [475 325;475 300], ...            % Integrator z  cont
    [475 200;487 200], ...            % Integrator x  cont
    [475 200;475 158;205 158], ...    % Integrator x  cont
    [475 95;475 135;340 135], ...     % Integrator y  cont
    [475 95;475 80], ...              % Integrator y  cont
    };
BRANCHES = {
    [475 330;475 425;530 425], ...               % -> Scope 3
    [475 300;535 300], ...                       % -> To Workspace z
    [475 300;236 300;236 370;255 370], ...       % -> Gain beta
    [475 330;499 330;499 142;170 142;170 95], ...% -> Sum rhoz 2
    [487 200;487 395;530 395], ...               % -> Scope 1
    [487 200;535 200], ...                       % -> To Workspace x
    [205 158;106 158;106 285;125 285], ...       % -> Product xy 1
    [205 158;205 65;255 65], ...                 % -> Product x(rho-z) 1
    [475 200;473 200;473 249;180 249;180 220], ...% -> Sum dx 2
    [340 135;340 150;82 150;82 300;125 300], ... % -> Product xy 2
    [340 135;340 110], ...                       % -> Sum dy 2
    [475 80;517 80;517 410;530 410], ...         % -> Scope 2
    [475 95;535 95], ...                         % -> To Workspace y
    [475 80;461 80;461 180;180 180;180 190], ... % -> Sum dx 1
    };
for k = 1:numel(TRUNKS),    add_line(mdl, TRUNKS{k});    end
for k = 1:numel(CONT),      add_line(mdl, CONT{k});      end
for k = 1:numel(BRANCHES),  add_line(mdl, BRANCHES{k});  end
end
