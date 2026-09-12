function divVal = lorenzDivergence(sigma, beta)
%LORENZDIVERGENCE  Phase-space divergence (trace of the Jacobian) for Lorenz.
%
%   divVal = lorenzDivergence(sigma, beta)
%
%   WHY this check (paper Sec. 3 properties of Lambda_k and Sec. 5 control):
%   The sum of the coefficient spectrum equals the flow divergence:
%       sum_k Lambda_k = trace(J) = div(f)
%   For Lorenz, div = -(sigma + 1 + beta) = -(1 + sigma + b) in paper notation.
%   Paper Sec. 5 uses this as the numerical control variable: every spectrum
%   estimate, by either method, must recover -13.666... for the canonical
%   parameters sigma=10, beta=8/3 (rho/r cancels). Any deviation signals an
%   integration or averaging bug before you claim new physics.
%
%   This is the same trace check used in src/analysis/lorenz_lyapunov_spectrum.m
%   and in src/simulate_lorenz_chaos.m table ("sum lambda (trace check)").
%
%   Inputs:
%     sigma, beta - Lorenz parameters
%   Output:
%     divVal - scalar divergence (negative for the dissipative Lorenz flow)

arguments
    sigma (1,1) double
    beta  (1,1) double
end

% Paper Sec. 5: "Como variable de control numerico, usamos la divergencia
% del sistema que es -(1+sigma+b) = -13.666" (with b == beta).
divVal = -(sigma + 1 + beta);
end
