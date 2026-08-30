function decorate_lorenz_sim()
%DECORATE_LORENZ_SIM  Add LaTeX-rendered annotations (title, equations)
%and a professional color palette to the lorenz_sim model.
%Does NOT change block positions, sizes, or signal lines.
%
% Palette (block BackgroundColor + TeX RGB):
%   lightBlue -> dx/dt part    (Sum dx, Gain sigma, Integrator x)
%   orange [255 128 21] -> dy/dt part (Const rho, Sum rhoz, Product x(rho-z), Sum dy, Integrator y)
%   purple [204 140 242] -> dz/dt part (Product xy, Gain beta, Sum dz, Integrator z)
%   gray [224 224 224] -> sinks    (To Workspace x/y/z, Scope)
% Requires R2023b+ for the annotation 'tex' Interpreter.

mdl = 'lorenz_sim';
if ~bdIsLoaded(mdl)
    open_system(mdl);
end

% --- block colors grouped by equation part ---
dxPart = {'Sum dx', 'Gain sigma', 'Integrator x'};
dyPart = {'Const rho', 'Sum rhoz', 'Product x(rho-z)', 'Sum dy', 'Integrator y'};
dzPart = {'Product xy', 'Gain beta', 'Sum dz', 'Integrator z'};
sinks  = {'To Workspace x', 'To Workspace y', 'To Workspace z', 'Scope'};

for k = 1:numel(dxPart), set_param([mdl '/' dxPart{k}], 'BackgroundColor', 'lightBlue'); end
for k = 1:numel(dyPart), set_param([mdl '/' dyPart{k}], 'BackgroundColor', '[1 0.502 0.082]'); end
for k = 1:numel(dzPart), set_param([mdl '/' dzPart{k}], 'BackgroundColor', '[0.80 0.55 0.95]'); end
for k = 1:numel(sinks),  set_param([mdl '/' sinks{k}],  'BackgroundColor', '[0.88 0.88 0.88]'); end

% --- annotations rendered with the TeX (LaTeX) interpreter ---
% Simulink annotations support only the MATLAB "TeX" subset: \bf, \rm,
% \fontsize, \color[rgb]{...}, Greek letters, sub/superscripts, \newline.
% NOT supported (LaTeX-only): \dot, \frac, \sqrt, \,. dx/dt notation is used
% for the time derivatives instead.
titleTxt = ['\fontsize{16}\bf Lorenz Attractor\newline\fontsize{11}\rm\color[rgb]{0.25 0.25 0.25}' ...
    '\sigma=10, \rho=28, \beta=8/3   |   nominal (x_0=1) vs perturbed (x_0=1+10^{-6})'];
eqTxt = ['\fontsize{13}\bf Equations\newline\fontsize{14}' ...
    '\color[rgb]{0.55 0.74 0.96}dx/dt = \sigma(y-x)\newline' ...
    '\color[rgb]{1.0 0.502 0.082}dy/dt = x(\rho - z) - y\newline' ...
    '\color[rgb]{0.80 0.55 0.95}dz/dt = xy - \beta z'];

% remove any previously created annotations
anns = find_system(mdl, 'FindAll', 'on', 'Type', 'annotation');
for i = 1:numel(anns)
    a = get_param(anns(i), 'Object');
    if contains(char(a.Text), {'Atractor de Lorenz', 'Lorenz Attractor', 'Ecuaciones', 'Equations', 'Leyenda', 'Legend', 'Simulation parameters'})
        delete(a);
    end
end

a1 = Simulink.Annotation([mdl '/Lorenz Attractor']);
a1.Text = titleTxt; a1.Interpreter = 'tex'; a1.TeXMode = 'on';
a1.Position = [-265 28 48 78];

a2 = Simulink.Annotation([mdl '/Equations']);
a2.Text = eqTxt; a2.Interpreter = 'tex'; a2.TeXMode = 'on';
a2.Position = [-265 95 -148 168];

% simulation-parameter box (below the equations)
stp = get_param(mdl, 'StopTime');
fst = get_param(mdl, 'FixedStep');
paramTxt = sprintf('\\fontsize{13}\\bf Simulation parameters\\newline\\fontsize{10}\\rmSolver: ODE4 (Runge-Kutta, fixed-step)\\newlineStep size: %s s\\newlineStop time: %s s', fst, stp);
a3 = Simulink.Annotation([mdl '/Simulation parameters']);
a3.Text = paramTxt; a3.Interpreter = 'tex'; a3.TeXMode = 'on';
a3.Position = [-265 185 -88 248];

save_system(mdl);
fprintf('Decorated: %d blocks colored and 3 LaTeX annotations in %s\n', ...
    numel(dxPart)+numel(dyPart)+numel(dzPart)+numel(sinks), mdl);
end
