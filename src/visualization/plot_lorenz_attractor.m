function plot_lorenz_attractor(varargin)
%PLOT_LORENZ_ATTRACTOR  Animated, rotating 3D Lorenz attractor.
%   Called automatically as the lorenz_sim model StopFcn (reads xsim/ysim/zsim
%   from the base workspace), with a Simulink.SimulationOutput:
%       plot_lorenz_attractor(out)        % out = sim(...)  (or out(1) for run 1)
%   or with raw trajectory arrays (e.g. from simulate_lorenz_chaos):
%       plot_lorenz_attractor(x, y, z, t)
%
%The trajectory is revealed progressively while the view rotates around the
%vertical (z) axis. Single-hue palette (light -> deep blue) -- not rainbow.

x = []; y = []; z = []; t = [];

if nargin >= 4
    x = varargin{1}; y = varargin{2}; z = varargin{3}; t = varargin{4};
elseif nargin >= 1 && isa(varargin{1}, 'Simulink.SimulationOutput')
    so = varargin{1};
    try
        x = so.get('xsim').Data;
        y = so.get('ysim').Data;
        z = so.get('zsim').Data;
        t = so.get('xsim').Time;
    catch
        return
    end
elseif evalin('base', 'exist(''xsim'', ''var'')')
    try
        x = evalin('base', 'xsim.Data');
        y = evalin('base', 'ysim.Data');
        z = evalin('base', 'zsim.Data');
        t = evalin('base', 'xsim.Time');
    catch
        return
    end
end

if isempty(x) || ~isnumeric(x), return; end

% Smooth, dense resampling (pchip avoids overshoot) so the curve is not
% faceted/straight -- ~8000 points regardless of the input density.
n  = numel(t);
if n >= 10
    tq = linspace(t(1), t(end), 8000)';
    xi = pchip(t, x, tq); yi = pchip(t, y, tq); zi = pchip(t, z, tq);
else
    xi = x; yi = y; zi = z;
end
N  = numel(xi);
if N < 2, return; end

% single-hue gradient (light -> deep blue), split into a few hundred colored
% animated polylines so the rotation stays smooth (few graphics objects).
K  = 2000;
ck = [linspace(0.62, 0.05, K)' , linspace(0.83, 0.34, K)' , linspace(1.00, 0.72, K)'];
chunkOf = min(K, max(1, ceil((1:N)'/N*K)));

fig = figure('Name', 'Lorenz Attractor', 'Color', 'w', 'Position', [200 120 760 620]);
ax  = axes('Parent', fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
axis(ax, 'vis3d');
view(ax, -27, 22);

xlabel(ax, 'x'); ylabel(ax, 'y'); zlabel(ax, 'z');
title(ax, 'Lorenz attractor  (\sigma=10, \rho=28, \beta=8/3)', 'FontSize', 12);
set(ax, 'FontSize', 11, 'GridColor', [0.85 0.85 0.85], 'LineWidth', 0.8, ...
        'XColor', [0.25 0.25 0.25], 'YColor', [0.25 0.25 0.25], 'ZColor', [0.25 0.25 0.25]);

% --- animated reveal + slow, continuous rotation about z ---
% One animated line per gradient chunk -> only ~K graphics objects, so the
% rotation is smooth instead of repainting thousands of segments each frame.
nF        = 1500;     % total frames (higher = the animation takes longer)
revealEnd = 0.55;    % fraction of frames spent drawing the trajectory
az0 = -27; el0 = 22;

hs = gobjects(K, 1);
for k = 1:K
    hs(k) = animatedline(ax, 'Color', ck(k, :), 'LineWidth', 1.2);
end

shown = 0; prevChunk = 0;
for fr = 1:nF
    target = min(N, ceil((fr / (nF*revealEnd)) * N));
    while shown < target
        shown = shown + 1;
        c = chunkOf(shown);
        addpoints(hs(c), xi(shown), yi(shown), zi(shown));
        if shown > 1 && c ~= prevChunk
            % overlap the boundary point into the previous chunk to avoid gaps
            addpoints(hs(prevChunk), xi(shown), yi(shown), zi(shown));
        end
        prevChunk = c;
    end
    view(ax, az0 + fr*0.1, el0);   % slow rotation, continues after the reveal
    drawnow limitrate
end
drawnow
end
