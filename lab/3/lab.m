%%
% Gwen Lofman, Lab 3, ISC4232C
%%
%%
% Given a the function $f$ where $-u'' = f$ in $[0,2pi]$, with Dirichlet
% boundary conditions: $u(0) = A$, $u(2pi) = B$ solve the BVP where $f$ is
% defined as:
%
% $$
% f = -6x - \sin(x)
% $$

f = @(x) -6*x - sin(x);
ua = 5;
ub = 8*(pi^3) + 5;

a = 0;
b = 2*pi;

%%
% The exact solution of the BVP for the given boundary conditions is:
%
% $$
% u(x) = x^3 - \sin(x) + 5
% $$

uexact = @(x) x.^3 - sin(x) + 5;

% Set up the grid of the numerical solution
N = 64;

%%
% Using the naiive matrix solve to approximate the solution to the BVP
% yeilds the following plot:

[U, xs] = fd_naive(f, a, b, ua, ub, N);

figure(1); hold on
plot([a;xs;b], [ua;U;ub]);
plot([a;xs;b], [ua;uexact(xs);ub], 'r--');
title("Full matrix solution");
legend("Approximation", "Exact Solution");
hold off

%%
% Using the sparse matrix solve to approximate the solution to the BVP
% yeilds the following plot:


[U, xs] = fd_sparse(f, a, b, ua, ub, N);

figure(2); hold on
plot([a;xs;b], [ua;U;ub]);
plot([a;xs;b], [ua;uexact(xs);ub], 'r--');
title("Sparse matrix solution");
legend("Approximation", "Exact Solution");
hold off

%%
% Using the $L_{\infty}$ norm and the $L_2$ norm to calculate the error, we
% can perform a convergence study.  For both, since we used a second order
% method, the ratio of the errors should be 4, and the slope of the error
% should be -2 in loglog scale.

% Determine the powers of two to use for the error
exps = 3:12;

% Expression for the exact solution
dt = @(a,b,N) (b-a)/(N+1);
exact = @(N) uexact(linspace(a+dt(a,b,N), b-dt(a,b,N),N))';

% Expression for the sparse solution
sparse = @(N) fd_sparse(f,a,b,ua,ub,N);

% Infinity Norm
l_inf = @(N) max(exact(N) - sparse(N));

% L2 Norm
l_2 = @(N) 1/sqrt(N) * sqrt(sum((exact(N)-sparse(N)).^2));

% Calculate the errors
N = repmat(2,[1,numel(exps)]).^exps;
E_inf = arrayfun(l_inf, N);
E_2 = arrayfun(l_2, N);

% Plot the error in loglog scale
figure(3);
loglog(N, E_inf);
hold on;
loglog(N, E_2);
title("Convergence Study");
legend("L_{\infty} norm", "L_2 norm");
xlabel("N_{elements}");
ylabel("Error");

fprintf("\\begin{array}{r|rr|rr}\n");
fprintf("N & L_{\\infty} & ratio & L_2 & ratio \\\\\\hline\n");
for e = 2:numel(N);
    einf = E_inf(e-1)/E_inf(e);
    e2 = E_2(e-1)/E_2(e);
    fprintf("%4i & %.3e & %.3f & %.3e & %.3f \\\\\n", N(e), E_inf(e), einf, E_2(e), e2);
end
fprintf("\\end{array}\n");