%%
% Solve the 1D Poisson equation -u'' = f in [0,1] with Dirichlet boundary
% conditions u(0) = A, u(1) - B

f = @(x) pi^2*sin(pi*x);
A = 0;
B = 0;

%%
% The exact solution of the BVP
uexact = @(x) sin(pi*x);

%%
% Set up the grid of the numerical solution

N = 10000;

% Include the start and end points, and then remove them
x = linspace(0,1,N+2)';
x = x(2:end-1);
dx = 1/(N+1); % mesh spacing

rhs = -dx^2*f(x);
rhs(1) = rhs(1) - A;
rhs(end) = rhs(end) - B;

% Set up the diagonals of the matrix which correspond to the coefficients
% of the second centered difference
a = ones(N,1);
b = -2*ones(N,1);
c = ones(N,1);

L = diag(a(1:end-1),1) + diag(b,0) + diag(c(1:end-1),-1);

e = ones(N,1);
Lsparse = spdiags([e -2*e e], -1:1, N,N);

tic
U = L\rhs;
fprintf('Speed with backslash %.5e\n', toc);

tic
U = Lsparse\rhs;
fprintf('Speed with sparse backslash %.5e\n', toc);


figure(1); clf; hold on
plot([0;x;1], [A;U;B]);
plot([0;x;1], [A;uexact(x);B], 'r--');

%%
%

tic
U = thomas(a,b,c,rhs);
fprintf('Speed with thomas algorithm %.5e\n', toc);

figure(2); clf; hold on
plot([0;x;1], [A;U;B]);
plot([0;x;1], [A;uexact(x);B], 'r--');
clear(U,A,a,b,c);
