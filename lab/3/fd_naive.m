function [U, xs] = fd_naive(f, a, b, ua, ub, n)
% FD_NAIVE Solve the BVP with dirichelet conditions using finite difference

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
da = ones(n-1,1);
db = -2*ones(n,1);
dc = ones(n-1,1);

L = diag(da,-1) + diag(db,0) + diag(dc,1);

% Find the solution by solving the system using \
U = L\rhs;

end

