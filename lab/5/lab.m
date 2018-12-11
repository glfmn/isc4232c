%%
% Gwen Lofman, Lab 5, ISC4232C
%%
%% Mean Weighted Residuals: Exercise 1
%
% First, we define the constants of our problem.

alpha = 0;
N = 10;

%%
% Then we define the function f and the function u_h which approximates our
% exact solution u.

w = @(N, x) cos(x(:) * N);

u_h = @(c, x) (alpha * (0:numel(c)-1).^2 .* w((0:numel(c)-1), x)...
    + w((0:numel(c)-1), x)) * c(:);

f = @(x) exp(cos(x));

%%
% Set up the right-hand side of the problem and the A matrix to find the
% coefficients.

I = @(n) besseli(n, 1);
rhs = 2*pi*I(0:N-1)';
d1 = @(n) pi*ones(size(n));
edge = @(n) sin(2*pi*n)./n;

A = spdiags(d1(1:N)', 0, N, N);
A(1,1) = 2*pi;
A(2:N,1) = edge(1:N-1);
A(1,2:N) = edge(1:N-1);
A(1:N,2:N) = (1 + alpha*(1:N-1).^2) .* A(1:N,2:N);

%%
% Calculate coefficients of the residual and plot the residual.

c = A\rhs;

x = linspace(0,2*pi,100);

figure;
plot(x, u_h(c,x), x, f(x), 'r--');
title(['Yukawa equation with \alpha = ' num2str(alpha,'%d')]);
legend("u_h", "f");

figure;
plot(x, u_h(c,x) - f(x)');
title("residual");
legend("u_h - f");

%% Exercise 2
% Solve the heat equation which is a special case of the yukawa equation
% once Rothe's method is applied.

%%
% Redefine constants according to the problem:
% $$
% \Delta t = 0.1
% $$
%
% $$
% \alpha = \Delta t
% $$

dt = 0.1;
alpha = dt;
N = 10;
nx = 1000;

%%
% We update our $u_h$ to fit the problem:
%
% $$
% u_h(x) = \sum_{i=1}^N c_iw_i(x)
% $$

u_h = @(c, x) w((0:numel(c)-1), x) * c(:);

%%
% Create quadrature points for updating the rhs

x = linspace(0,2*pi, nx);

%%
% Create A with $\alpha = dt$

A = spdiags(d1(1:N)', 0, N, N);
A(1,1) = 2*pi;
A(2:N,1) = edge(1:N-1);
A(1,2:N) = edge(1:N-1);
A(1:N,2:N) = (1 + alpha*(1:N-1).^2) .* A(1:N,2:N);

rhs = 2*pi*I(0:N-1)';

%%
% Complete the exercise by plotting the change in u_h for each integral t.

figure;
hold on;
legend('Location', 'NorthEastOutside');
for t = 0:dt:5
    c = A\rhs;

    if abs(floor(t) - t) < 0.01
        plot(x/pi, u_h(c, x), 'DisplayName', ['u_h at t=' num2str(t, '%i')]);
    end

    rhs = sum(u_h(c, x) .* w((0:N-1), x))' * 2*pi/nx;
end
title("Heat equation IBVP");
xlabel("\pi*x");
ylabel("temperature");