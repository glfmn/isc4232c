% Forward in space, centered in time

% Solve heat equation in 1D with D = 0.4, domain is (0,1), time horizon is
% T = 1, M = 100 time steps, N = 9 interior spacial points, and the initial
% condition is sin(pi*x)

a = 0;
b = 1;
D = 0.4;
N = 20;
M = 400;
T = 1;
dx = (b-a)/(N+1);
dt = T/M;

lambda = D*dt/dx^2

x = (a:dx:b); % Spacial Domain
U = @(x) sin(pi*x); % Initial Condition

plot(x, U(x));
xlabel('x');
ylabel('u');
ylim([a b]);
title("Initial condition of u(x,t)")

U = U(x);
for k = 1:M
    U(1) = 0;
    U(end) = 0;
    U(2:end-1) ...
        = lambda * U(1:end-2) ...
        + (1-2*lambda) * U(2:end-1) ...
        + lambda * U(3:end);

    plot(x, U);
    xlabel('x');
    ylabel('u');
    ylim([a b]);
    title(['Time ' num2str(k*dt, '%4.2e' )]);
    pause(dt);
end