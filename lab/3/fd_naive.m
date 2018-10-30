function [U, xs] = fd_naive(f, a, b, A, B, N)
%FD_NAIVE Solve the BVP with dirichelet conditions using finite difference

% Create mesh of x values
xs = linspace(a,b,N+2)';
xs = xs(2:end-1);
dx = 1/(N+1); % mesh spacing

% Calculate the right-hand side of the matrix solve
rhs = -dx^2*f(xs);
rhs(1) = rhs(1) - A;
rhs(end) = rhs(end) - B;

% From the diagonals construct the matrix which represents the finite
% difference
da = ones(N-1,1);
db = -2*ones(N,1);
dc = ones(N-1,1);

L = diag(da,1) + diag(db,0) + diag(dc,-1);

% Find the solution by solving the system using \
U = L\rhs;

end

