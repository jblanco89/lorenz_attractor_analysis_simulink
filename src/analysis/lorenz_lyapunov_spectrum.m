function lambda = lorenz_lyapunov_spectrum()
%LORENZ_LYAPUNOV_SPECTRUM  Full Lyapunov spectrum of the Lorenz (1963) system.
%   lambda = lorenz_lyapunov_spectrum()
%
%   Returns the three Lyapunov exponents of the Lorenz system
%       dx/dt = sigma*(y - x)
%       dy/dt = x*(rho - z) - y
%       dz/dt = x*y - beta*z
%   as a 3x1 vector [lambda_1; lambda_2; lambda_3] with lambda_1 >= lambda_2 >=
%   lambda_3. The parameters sigma, rho, beta are read from the base workspace
%   and fall back to the chaotic-regime defaults sigma=10, rho=28, beta=8/3.
%
%   Method (Benettin et al., Gram-Schmidt variant): the 3x3 variational system
%   dQ/dt = J(x)*Q is integrated alongside the orbit, i.e. a 12-state coupled
%   ODE. After every window of dT seconds the tangent basis Q is reorthonormalized
%   by Gram-Schmidt; the log of each column norm, divided by dT, is one sample of
%   the corresponding exponent. The first 25% of the samples (transient needed for
%   the basis to align with the expanding directions) is discarded before averaging.
%
%   Accuracy check: sum(lambda) must equal trace(J) = -(sigma + 1 + beta)
%   (the contraction rate of phase-space volume).

% --- parameters from the base workspace (chaotic regime defaults) ---
try
    sigma = evalin('base', 'sigma'); rho = evalin('base', 'rho'); beta = evalin('base', 'beta');
catch
    sigma = 10; rho = 28; beta = 8/3;
end

% --- Lorenz vector field (states 1:3) ---
fLorenz = @(x) [sigma*(x(2)-x(1)); ...
                x(1)*(rho-x(3))-x(2); ...
                x(1)*x(2)-beta*x(3)];

% --- Jacobian of the vector field, evaluated along the orbit ---
J = @(x) [-sigma, sigma, 0; ...
          rho-x(3), -1, -x(1); ...
          x(2), x(1), -beta];

% --- augmented system: orbit + tangent basis (Q stored as 9 states) ---
fBen = @(t, y) [fLorenz(y(1:3)); ...
                reshape(J(y(1:3))*reshape(y(4:12), 3, 3), 9, 1)];

% tight tolerances: the variational equations are linear and their growth must
% be resolved accurately over a long (500 s) horizon
opts = odeset('AbsTol', 1e-10, 'RelTol', 1e-11);

% --- integration parameters ---
x0  = 1.0; y0 = 1.0; z0 = 1.0;   % initial condition (same as the nominal run)
Q0  = eye(3);                    % orthonormal tangent basis (unit columns)
dT  = 0.5;                       % reorthonormalization window (s)
Tben = 500;                      % total integration horizon (s); long enough for
                                 % lambda_1 to converge (see AGENTS.md check)

% --- Benettin loop: integrate window-by-window, reorthonormalize each dT ---
nWin   = ceil(Tben/dT);           % total number of windows
lamB   = zeros(3, nWin);          % per-window log-growth samples (3 x nWin)
state  = [x0; y0; z0];
Q      = Q0;
tB     = 0;
iw     = 0;
while tB < Tben
    Y = ode45(fBen, [tB tB+dT], [state; Q(:)], opts);
    s = Y.y(:, end);
    state = s(1:3);              % continue the orbit
    V = reshape(s(4:12), 3, 3);  % evolved tangent basis (columns)
    iw = iw + 1;

    % Gram-Schmidt orthonormalization, recording each column's growth
    for k = 1:3
        r = norm(V(:, k));
        if ~(r > 0) || ~isfinite(r)
            error('lorenz_lyapunov_spectrum:nonfinite', ...
                'Tangent norm is not positive and finite at t = %.1f s.', tB);
        end
        lamB(k, iw) = log(r)/dT;
        V(:, k) = V(:, k)/r;                     % normalize q_k
        for j = (k+1):3                          % remove q_k from the rest
            V(:, j) = V(:, j) - (V(:, j).'*V(:, k))*V(:, k);
        end
    end
    Q = V;                       % reorthonormalized basis for the next window
    tB = tB + dT;
end

% --- discard the initial transient, then average the remaining samples ---
k = round(0.25*nWin);
lambda = mean(lamB(:, (k+1):nWin), 2);

end