function J = lorenzJacobian(state, sigma, rho, beta)
%LORENZJACOBIAN  Jacobian of the Lorenz (1963) vector field at a state point.
%
%   J = lorenzJacobian(state, sigma, rho, beta)
%
%   Lorenz system (paper Eq. 5.1, Sec. 5):
%       dx/dt = sigma*(y - x)
%       dy/dt = x*(rho - z) - y
%       dz/dt = x*y - beta*z
%
%   WHY this Jacobian (paper Sec. 3, Eq. 4 and Sec. 5):
%   The variational equation is deltaDot = J(t) * delta (Eq. 4). J(t) is the
%   Jacobian evaluated along the orbit. The eigenvalue method freezes J at the
%   left endpoint of each integration window (J_i = J(t_i), Eq. 5) and
%   diagonalizes it. That frozen J_i is exactly what this function returns
%   for a single state X(t_i). In linear dynamics J is the system matrix
%   itself, so the computation simplifies (§3 paragraph after Eq. 6).
%
%   Inputs:
%     state  - [x y z] row or column vector (1x3 or 3x1)
%     sigma, rho, beta - Lorenz parameters (paper uses sigma, r, b;
%                        repo uses sigma, rho, beta; rho == r)
%
%   Output:
%     J - 3x3 Jacobian matrix

arguments
    state (1,3) double
    sigma (1,1) double
    rho   (1,1) double
    beta  (1,1) double
end

x = state(1); y = state(2); z = state(3);

% Paper §5, Eq. 5.1 Jacobian derived analytically:
%   df1/dx = -sigma, df1/dy = sigma, df1/dz = 0
%   df2/dx = rho - z, df2/dy = -1, df2/dz = -x
%   df3/dx = y, df3/dy = x, df3/dz = -beta
J = [-sigma, sigma, 0; ...
     rho - z, -1,   -x; ...
     y,        x,  -beta];
end
