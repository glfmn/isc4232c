function [U, xs] = fd_periodic(f, a, b, n)
%FD_PERIODIC Solve the preiodic BVP with dirichelet boundary conditions.
%   Solves using finite difference. Modified for periodic functions with no
%   explicit boundary condition but where u(a) = u(b).
%
%   [U, xs] = FD_PERIODIC(f, a, b, n) will return the approximate solution
%   U, and the x inputs where it was evaluated xs where f is the function
%   which defines the PDE, the boundary is at a,b, and the domain is (a,b),
%   and n is the number of discretization points

% Create mesh of x values
xs = linspace(a,b,n+2)';
xs = xs(2:end-1);
dx = (a-b)/(n+1); % mesh spacing

% Calculate the right-hand side of the matrix solve
rhs = -dx^2*f(xs);

% From the diagonals construct the matrix which represents the finite
% difference
e = ones(n,1);
L = spdiags([e -(2+dx^2)*e e], -1:1, n,n);
L(1,end) = 1;
L(end,1) = 1;

% Find the solution by solving the system using \
U = L\rhs;

end

