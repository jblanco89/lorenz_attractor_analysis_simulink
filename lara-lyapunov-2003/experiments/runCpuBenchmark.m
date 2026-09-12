function results = runCpuBenchmark(varargin)
%RUNCPUBENCHMARK  Pure-MATLAB timing harness replicating paper Sec. 4 (Table 1).
%
%   results = runCpuBenchmark('Name',Value,...)
%
%   WHY this benchmark (paper Sec. 4, Table 1):
%   Sec. 4 compares CPU time of the standard Gram-Schmidt method (Benettin,
%   ref. 6, 8) vs the proposed eigenvalue-average method (Sec. 3) on the same
%   equations, parameters and stopping conditions. Linear oscillators of
%   dimension m=2..20 (even, step 2) are used. Finding: methods tie at m=2
%   and the proposed method wins increasingly with m (up to ~7.5x at m=18),
%   with identical convergence in number of iterations — the gain is cost per
%   iteration (Sec. 4 last paragraph).
%
%   This harness is pure MATLAB (no Simulink models), as required by user
%   decision 1, reusing src/analysis/gram_schmidt.m for the GS leg.
%   It is intentionally solver-agnostic and focuses on the per-window cost
%   that distinguishes the two estimators, in the spirit of the paper's
%   "estandarizar los codigos numericos" (Sec. 4).
%
%   Method for each m (all via numeric eig/QR, per Sec. 4 "restringir todos
%   las operaciones a calculos numericos"):
%     Linear oscillator matrix A_m (2x2 blocks [0 1; -k -c]) built by
%     buildOscillatorMatrix(m). For the benchmark it is the Jacobian itself
%     (linear dynamics => J is the system matrix, Sec. 3).
%     Eigenvalue leg: nIter windows, each computes eig(A_m) and accumulates
%       Re(eig) (the analogue of Lambda_k^n, Eq. 7). No QR, no tangent basis.
%     Gram-Schmidt leg: nIter windows, each evolves the m x m tangent basis
%       as V = expm(A_m * dT) * Q and reorthonormalizes with [Q,R]=gram_schmidt(V),
%       accumulating lamGS = log(diag(R))/dT (same formula used in
%       src/analysis/lorenz_lyapunov_spectrum.m). This mirrors the Benettin
%       window without per-window ODE integration, exposing the QR cost.
%
%   Options:
%     'Dims'      - vector of dimensions (default 2:2:20)
%     'NIter'     - windows per timing (default 2000; increase for stabler clock)
%     'Dt'        - effective window length in s for GS leg (default 0.5, matches Tben loop)
%     'NRepeats'  - repeats per method for median timing (default 3)
%     'Verbose'   - print per-m progress (default true)
%     'SavePath'  - path to save results .mat (default '' = no save)
%
%   Output:
%     results - struct array with fields: m, timeEig, timeGS, speedup, nIter, dt
%
%   Note: Absolute seconds reflect the modern machine (not the 2003 Pentium
%   MMX/32MB from Table 1); the ratio/speedup trend is the replicated signal.

p = inputParser;
addParameter(p, 'Dims', 2:2:20);
addParameter(p, 'NIter', 2000);
addParameter(p, 'Dt', 0.5);
addParameter(p, 'NRepeats', 3);
addParameter(p, 'Verbose', true);
addParameter(p, 'SavePath', "");
parse(p, varargin{:});
opts = p.Results;

dims     = opts.Dims;
nIter    = opts.NIter;
dt       = opts.Dt;
nRepeats = opts.NRepeats;
verbose  = opts.Verbose;
savePath = string(opts.SavePath);

thisFile = mfilename('fullpath');
thisDir  = fileparts(thisFile);
projRoot = fileparts(thisDir); projRoot = fileparts(projRoot);
addpath(fullfile(projRoot, 'src','analysis'));

if verbose
    fprintf('=== CPU benchmark (paper Sec. 4, Table 1) nIter=%d dt=%.2f ===\n', nIter, dt);
    fprintf('%4s | %10s %10s | %7s\n', 'm', 't_eig(s)', 't_GS(s)', 'speedup');
    fprintf('%s\n', repmat('-',1,40));
end

results = struct('m',{},'timeEig',{},'timeGS',{},'speedup',{},'nIter',{},'dt',{});

for idx = 1:numel(dims)
    m = dims(idx);

    % Paper Sec. 4: "La dinamica elegida corresponde a osciladores lineales."
    % Each 2D block is a linear oscillator m*x_ddot + c*x_dot + k*x = 0
    % with c=0.3, k=1+0.2*i (deterministic, stable).
    A = buildOscillatorMatrix(m);

    % Precompute the per-window propagator for the GS leg (same matrix every window)
    Phi = expm(A * dt);

    % --- time the eigenvalue leg (median of nRepeats) ---
    tEigReps = zeros(1, nRepeats);
    for rep = 1:nRepeats
        tic;
        % Paper Sec. 4: "restringir todos las operaciones a calculos numericos"
        % so eig is called numerically every window even though A is constant and
        % analytic formulas would exist for m<=4.
        s = zeros(m,1);
        for k = 1:nIter
            lam = eig(A);
            s = s + real(lam); %#ok<NASGU> prevent elimination
        end
        % Touch s to prevent dead-code elimination (defensive; not strictly needed)
        if ~isfinite(s(1)), error('nonfinite'); end
        tEigReps(rep) = toc;
    end
    timeEig = median(tEigReps);

    % --- time the Gram-Schmidt leg (median) ---
    tGsReps = zeros(1, nRepeats);
    for rep = 1:nRepeats
        tic;
        Q = eye(m);
        acc = zeros(m,1);
        for k = 1:nIter
            V = Phi * Q;               % exact window propagation
            [Q, R] = gram_schmidt(V);  % MGS QR, diag(R)>0
            lamWin = log(diag(R))/dt;  % Benettin per-window sample (m x 1)
            acc = acc + lamWin; %#ok<NASGU>
        end
        if ~isfinite(acc(1)), error('nonfinite'); end
        tGsReps(rep) = toc;
    end
    timeGS = median(tGsReps);

    speedup = timeGS / max(timeEig, eps);

    results(idx).m = m;
    results(idx).timeEig = timeEig;
    results(idx).timeGS  = timeGS;
    results(idx).speedup = speedup;
    results(idx).nIter   = nIter;
    results(idx).dt      = dt;

    if verbose
        fprintf('%4d | %10.4f %10.4f | %7.2f\n', m, timeEig, timeGS, speedup);
    end
end

if verbose
    fprintf('%s\n', repmat('-',1,40));
    fprintf('Interpretation (paper Sec. 4): tie at m=2, growing advantage of eigenvalue method with m.\n');
    fprintf('Absolute seconds are machine-dependent; the speedup trend is the replicated finding.\n\n');
end

if strlength(savePath) > 0
    outDir = fileparts(savePath);
    if strlength(outDir) > 0 && ~exist(outDir,'dir')
        mkdir(outDir);
    end
    save(savePath, 'results', 'opts');
    if verbose, fprintf('Saved benchmark results to %s\n', savePath); end
end
end

function A = buildOscillatorMatrix(m)
%BUILDOSCILLATORMATRIX  Block-diagonal linear oscillator system of dimension m.
%
%   Each oscillator i contributes a 2x2 block:
%       [ 0     1
%        -k_i  -c ]
%   with c=0.30, k_i = 1 + 0.2*i  (stable, underdamped, deterministic).
%   If m is odd, the last row/col is a single damped mode [-c].
%
%   WHY this shape (paper Sec. 4): Table 1 studies "osciladores lineales"
%   from 2 to 20 equations. A natural canonical form is block-diagonal
%   oscillators; trace and stability are controlled and reproducible.
%   A dense random A would hide the oscillator nature, while this keeps
%   the dimension sweep interpretable.
%
%   Equipped to handle any m, but the paper only needs even m.

A = zeros(m);
nOsc = floor(m/2);
c = 0.30;
for i = 1:nOsc
    k = 1 + 0.2*i;
    r = 2*i-1;
    A(r, r)   = 0;
    A(r, r+1) = 1;
    A(r+1, r) = -k;
    A(r+1, r+1) = -c;
end
if mod(m,2)==1
    A(m,m) = -c;
end
end
