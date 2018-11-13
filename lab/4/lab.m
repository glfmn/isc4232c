%% Lab 4: 2D Finite Difference
%%%
% *Gwen Zapata*
%%%
% *ISC4220C*
%%
%%
clear

%% Poisson Equation
% As a prototypical equation for studying the properties of 2D finite
% difference equations, we can use the Poisson Equation:
% $-\delta u(x,y) = f(x,y)$ where $\delta u$ is shorthand for
% $ u_{xx} + u_{yy} $.  Putting all of this together, we have:
%
% $$
% -u_{xx}-u_{yy} = f
% $$

%%
% Using this to build an example problem, we establish that, where 
% $x \in (0,1), y \in (0,1)$:
%
% $$
% -\delta u = -2\pi^2 \sin(\pi x) \cos(\pi y)
% $$

f = @(x, y) 2*pi^2 .* sin(pi*x) .* cos(pi*y);

%%
% And we have the boundary conditions:
%
%
% $$
% u(x,0) = sin(\pi x), u(x, 1) = -sin(\pi x)
% $$
%
% $$
% u(0, y) = 0, u(1, y) = 0,
% $$

% Conditions where y = 0 and 1 respectively
br = @(x) 0.*x;
bl = @(x) 0.*x;

% Conditions where x = 0 and 1 respectively
bb = @(y) sin(pi*y);
bt = @(y) -sin(pi*y);

%%
% We can use the sparse matrix solver to apply a finite difference stencil
% to the 2D PDE.
N = 2.^(2:7);

%%
% Performing a convergence study on this approach with
% $ N {4,8,16,32,64,128} $ and comparing the error with the $l_2$ norm
% against the exact solution $ u(x,y) = \sin(\pi x)\cos(\pi y) $.

u = @(x,y) sin(pi.*x) .* cos(pi.*y);

% Expression for the sparse solution
approximate = @(N) finitediff2d(f, N, br, bl, bb, bt);

% L2 Norm
l_2 = @(E, A) 1/sqrt(numel(A)) * sqrt(sum((E-A).^2));

E_2 = 1:numel(N);
for e = E_2
    [approx, x, y] = approximate(N(e));
    exact = u(x,y);
    E_2(e) = l_2(exact, approx);
end

% Plot the error in loglog scale
figure(1);
loglog(N, E_2);
title("Convergence Study");
legend("L_2 norm");
xlabel("N_{elements}");
ylabel("Error");

fprintf("\\begin{array}{r|rr}\n");
fprintf("   N & L_2   & ratio \\\\\\hline\n");
for e = 2:numel(N)
    e2 = E_2(e-1)/E_2(e);
    fprintf("%4i & %.3f & %.3f \\\\\n", N(e), E_2(e), e2);
end
fprintf("\\end{array}\n");