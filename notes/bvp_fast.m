%%
% Solve the 1D equation -u'' = f in [0,2pi] with Dirichlet boundary
% conditions u(0) = A, u(2pi) - B

f = @(x) -6*x + sin(x);
A = 5;
B = 8*pi^3 + 5;

x0 = 0;
xn = 2*pi;

%%
% The exact solution of the BVP
uexact = @(x) x.^3 - sin(x) + 5;

%%
% Set up the grid of the numerical solution

N = 64;

% Include the start and end points, and then remove them
x = linspace(x0,xn,N+2)';
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

%e = ones(N,1);
%Lsparse = spdiags([e -2*e e], -1:1, N,N);
%U = Lsparse\rhs;

U = L\rhs;


figure(1); clf; hold on
plot([x0;x;xn], [A;U;B]);
plot([x0;x;xn], [A;uexact(x);B], 'r--');

%%
%
U = thomas(a,b,c,rhs);

figure(2); clf; hold on
plot([x0;x;xn], [A;U;B]);
plot([x0;x;xn], [A;uexact(x);B], 'r--');
