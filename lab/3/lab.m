%%
% Solve the 1D equation -u'' = f in [0,2pi] with Dirichlet boundary
% conditions u(0) = A, u(2pi) - B

f = @(x) -6*x - sin(x);
A = 5;
B = 8*pi^3 + 5;

x0 = 0;
xn = 2*pi;

%%
% The exact solution of the BVP
uexact = @(x) x.^3 - sin(x) + 5;

%%
% Set up the grid of the numerical solution

N = 100;

xs = linspace(x0,xn,N+2)';
xs = xs(2:end-1);
U = fd_naive(f, x0, xn, A, B, N);

figure(1); clf; hold on
plot([x0;xs;xn], [A;U;B]);
plot([x0;xs;xn], [A;uexact(x);B], 'r--');

%%
%

% xs = linspace(a, b,N+2)';
% xs = xs(2:end-1);
% dx = 1/(N+1); % mesh spacing
% 
% rhs = -dx^2*f(xs);
% rhs(1) = rhs(1) - A;
% rhs(end) = rhs(end) - B;
% 
% e = ones(N,1);
% Lsparse = spdiags([e -2*e e], -1:1, N,N);
% U = Lsparse\rhs;
% 
% figure(2); clf; hold on
% plot([x0;x;xn], [A;U;B]);
% plot([x0;x;xn], [A;uexact(x);B], 'r--');
