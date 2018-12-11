function [U, xs] = fd_sparse(f, a, b, ua, ub, n)
% FD_SPARSE Solve the BVP with dirichelet conditions with finite difference

% Create mesh of x values
xs = linspace(a,b,n+2)';
xs = xs(2:end-1);
dx = (a-b)/(n+1); % mesh spacing

% Calculate the right-hand side of the matrix solve
rhs = -dx^2*f(xs);
rhs(1) = rhs(1) - ua;
rhs(end) = rhs(end) - ub;

% From the diagonals construct the matrix which represents the finite
% difference
e = ones(n,1);
L = spdiags([e -2*e e], -1:1, n,n);

% Find the solution by solving the system using \
U = L\rhs;

end

