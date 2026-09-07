function [Q, R] = gram_schmidt(V)
%GRAM_SCHMIDT  Modified Gram-Schmidt QR factorisation with positive diagonal.
%   [Q,R] = gram_schmidt(V) computes the thin QR factorisation of the matrix
%   V (M-by-N, M >= N) using the Modified Gram-Schmidt (MGS) process.
%
%   The factorisation satisfies:
%       V = Q*R,  Q'*Q = I_N,  R upper-triangular with diag(R) > 0.
%
%   This is the thin/economy variant (Q is M-by-N, R is N-by-N), mathematically
%   equivalent to MATLAB's [Q,R] = qr(V,0) up to the sign convention on the
%   diagonal. Unlike qr, the diagonal is forced positive so that log(diag(R))
%   is well-defined without abs(), as required by the Benettin method for
%   Lyapunov exponents (see lorenz_lyapunov_spectrum).
%
%   Method
%   ------
%   Right-looking MGS (Golub & Van Loan, Matrix Computations, 4th ed.,
%   Alg. 5.2.8; Bjorck, Numerics of Gram-Schmidt orthogonalization):
%
%       for k = 1:N
%           R(k,k) = norm(V(:,k))          % norm of the k-th column after
%                                          % removal of span{q1..q_{k-1}}
%           Q(:,k) = V(:,k) / R(k,k)       % normalise
%           for j = k+1:N
%               R(k,j) = Q(:,k)' * V(:,j)  % projection (V(:,j) is already
%                                          % cleaned from q1..q_{k-1})
%               V(:,j) = V(:,j) - R(k,j)*Q(:,k)
%           end
%       end
%
%   The inner vector V(:,j) is the progressively updated column, not the
%   original A(:,j). This is the defining difference between MGS and Classical
%   Gram-Schmidt (CGS): CGS projects from the original column
%   (R(i,j)=q_i'*A(:,j)), while MGS projects from the already-orthogonalised
%   intermediate vector. In exact arithmetic both coincide; in floating point
%   MGS satisfies ||Q'*Q - I|| = O(eps*kappa(V)) versus O(eps*kappa(V)^2) for
%   CGS (Bjorck 1967, 1994). For Benettin windows (Lorenz, dT=0.5) the tangent
%   basis is exponentially ill-conditioned (kappa ~ exp((lambda1-lambda3)*dT)
%   ~ 1.8e3), so CGS loses orthogonality within a few hundred windows and
%   corrupts log(R(k,k))/dT, while MGS remains accurate over the full
%   Tben=500 s horizon (1000 reorthonormalisations).
%
%   Numerical details
%   -----------------
%   - Diagonal forced positive: R(k,k) = norm(...) > 0, so no sign ambiguity.
%     MATLAB's qr may return R(k,k) < 0; taking log(abs(diag(R))) would mask
%     the convention difference. This function guarantees log(diag(R)) valid.
%   - Real and complex matrices are supported (uses Q(:,k)' for Hermitian
%     transpose; for real data this is identical to dot product).
%   - Input is not modified in the caller (copy-on-write in MATLAB); the
%     internal working copy V is overwritten column-by-column.
%   - Complexity O(M*N^2) flops, O(M*N) memory.
%
%   Input
%   -----
%   V : (M x N) double, M >= N, finite (no NaN/Inf), non-empty.
%
%   Outputs
%   -------
%   Q : (M x N) double, columns orthonormal to machine precision.
%   R : (N x N) double, upper-triangular, diag(R) > 0.
%
%   Errors
%   ------
%   - gram_schmidt:invalidInput       - V is not 2-D, empty, or M < N.
%   - gram_schmidt:nonfinite          - column norm is non-positive or non-finite
%                                       (rank deficiency or overflow).
%
%   Example
%   -------
%   V0 = [1 1 0; 1 0 1; 0 1 1];
%   [Q,R] = gram_schmidt(V0);
%   norm(V0 - Q*R)          % < 1e-15
%   norm(Q'*Q - eye(3))     % < 1e-15
%   all(diag(R) > 0)        % true
%
%   % Benettin usage (cf. lorenz_lyapunov_spectrum):
%   %   [Q,R] = gram_schmidt(V);  lamB(:,iw) = log(diag(R))/dT;
%
%   See also QR, LORENZ_LYAPUNOV_SPECTRUM.
%
%   References
%   ----------
%   [1] G. Benettin et al., Meccanica 15, 9-20, 1980.
%   [2] G. H. Golub & C. F. Van Loan, Matrix Computations, 4th ed., 2013.
%   [3] A. Bjorck, Numer. Math. 9, 1967; BIT 34, 1994.
%   [4] C. Skokos, Lect. Notes Phys. 790, 63-135, 2010.

arguments
    V (:,:) double
end

% --- validation: 2-D, non-empty, M >= N, finite ---
if isempty(V)
    error('gram_schmidt:invalidInput', 'Input V must be non-empty.');
end
[m, n] = size(V);
if m < n
    error('gram_schmidt:invalidInput', ...
        'Input V must satisfy M >= N (got %d-by-%d). Thin QR requires tall matrix.', m, n);
end
if ~all(isfinite(V(:)))
    error('gram_schmidt:invalidInput', ...
        'Input V must be finite (no NaN/Inf).');
end

% --- allocate outputs ---
Q = zeros(m, n);
R = zeros(n, n);

% Work on a copy to avoid aliasing the caller's variable
W = V;

for k = 1:n
    % Norm of the k-th column after removal of span{q1..q_{k-1}}.
    % Using 2-norm (Frobenius for vector); MATLAB's norm is numerically
    % stable (scaling via LAPACK's DNRM2, avoids overflow for large V).
    rk = norm(W(:, k));

    % Guard: rank deficiency, overflow, or NaN propagation.
    % Mirrors the check previously in lorenz_lyapunov_spectrum:73 but with
    % column index for diagnostics; identifier changed from
    % lorenz_lyapunov_spectrum:nonfinite to gram_schmidt:nonfinite.
    if ~(rk > 0) || ~isfinite(rk)
        error('gram_schmidt:nonfinite', ...
            'Column %d norm is not positive and finite (rk=%g). Input may be rank-deficient or contain overflow.', k, rk);
    end

    R(k, k) = rk;
    Q(:, k) = W(:, k) / rk;

    % Orthogonalise the remaining columns against q_k.
    % W(:,j) already lacks components in span{q1..q_{k-1}} (MGS property),
    % so R(k,j) = q_k' * W(:,j) equals q_k' * A(:,j) in exact arithmetic
    % but is more accurate in floating point.
    for j = (k+1):n
        rkj = Q(:, k)' * W(:, j);
        R(k, j) = rkj;
        W(:, j) = W(:, j) - rkj * Q(:, k);
    end
end

end
