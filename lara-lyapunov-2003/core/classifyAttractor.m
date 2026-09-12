function label = classifyAttractor(Lambda, divergence)
%CLASSIFYATTRACTOR  Classify asymptotic regime from the coefficient spectrum.
%
%   label = classifyAttractor(Lambda, divergence)
%
%   WHY these rules (paper Sec. 3, properties of Lambda_k, and Sec. 5, Table 2
%   discussion):
%   Lambda_k are time-averages of Re(eig(J)) along the orbit (§3, Eq. 7/8).
%   Their signs carry the stability verdict (same contract as the Lyapunov
%   spectrum: the paper only claims sign/coarse-description equivalence, not
%   digit equality — Sec. 3 last paragraphs and Sec. 6). Divergence = sum Lambda
%   gives the volume-contraction check (§3: "la suma es la divergencia del flujo
%   vectorial de la dinamica").
%
%   Decision table (paper §3 + §5 "Clase de Atractor"):
%     All Lambda < 0                        -> asymptotically stable orbit,
%                                            converges to a fixed point.
%                                            In Lorenz this is "Pto. Fijo" for
%                                            r = 20..24.
%     Exactly one Lambda > 0, rest < 0,
%       AND divergence < 0                  -> unstable but volume-contracting.
%                                            Paper conjectures chaotic dynamics
%                                            ("creemos que la dinamica es caotica
%                                            ya que al menos existe una direccion
%                                            donde hay inestabilidad y cualquier
%                                            elemento de volumen ... tiende a cero").
%                                            In Lorenz this is "Extraño" (strange
%                                            attractor) for r = 25..30.
%     divergence > 0                        -> strongly unstable, no attractor
%                                            ("fuerte inestabilidad", §3).
%
%   Inputs:
%     Lambda     - m x 1 vector of spectrum coefficients (Lambda_k)
%     divergence - scalar sum(Lambda) or theoretical div(f); if omitted,
%                  sum(Lambda) is used.
%
%   Output:
%     label - string: 'FixedPoint', 'Strange', or 'StrongInstability'
%             (Spanish equivalents in comments: Pto. Fijo / Extraño)

arguments
    Lambda (:,1) double
    divergence (1,1) double = sum(Lambda)
end

% Tolerance for sign: numerical Lambda ~ 0 should count as non-positive.
% Paper Table 2: Gram-Schmidt lambda2 ~ 1e-4 ("nulo") vs eigenvalue method
% lambda2 ~ -2.07 — the verdict still matches, so a small epsilon is appropriate.
tol = 1e-9;

nPos = sum(Lambda > tol);
nNeg = sum(Lambda < -tol);

if all(Lambda < -tol)
    % Paper §3: "Si todos los coeficientes Lambda_k son negativos, entonces la
    % orbita es asintoticamente estable, y por tanto la misma converge a un
    % conjunto estable." (§5: Pto. Fijo for r<=24).
    label = "FixedPoint";  % Pto. Fijo
elseif divergence > tol
    % Paper §3: "Cuando la divergencia es positiva, estamos en un caso de
    % fuerte inestabilidad."
    label = "StrongInstability";
elseif nPos == 1 && nNeg == numel(Lambda) - 1 && divergence < -tol
    % Paper §3: one positive, rest negative, divergence negative -> candidate
    % chaotic ("un caso de inestabilidad ... creemos que la dinamica es caotica").
    % Paper §5 validates this as "Extraño" for Lorenz r>=25.
    label = "Strange";  % Extraño / strange attractor
elseif nPos >= 1 && divergence < -tol
    % Generalization: any positive direction with overall contraction still
    % unstable / strange-like in the paper's conjecture.
    label = "Strange";
else
    % Fallback: unstable but not matching the clean cases above
    label = "StrongInstability";
end
end
