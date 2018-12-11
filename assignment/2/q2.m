%% Define the BVP as a system of IVPs

% Define the system of first-order IVPs
 
% w' = u''
w = @(y, t) y;
% v' = u'
v = @(y, t) 20*y*(y^2-1);

A = [0 1;
     1 0];
f = @(y, t) A*[w(y(1), t); v(y(2), t)];

% Boundary condition of the original problem
A = -1; % u(0)
B = 1; % u(1)

% Initial condition as a function of our guess for $w_0$ where $u(0)
Y0 = @(s) [s; A];

dt = 0.01;
s = 1;
[y, t] = forward_euler(f, Y0(s), dt, 0, 1);

%% Test plot

figure(1);
plot(t, y(1,:), t, y(2,:));
title(['Initial condition of [' num2str(s, '%d') '; -1]']);
legend('u_x', 'u');

%% Perform the shooting method

% Get the boundary condition
boundary = @(A, r) A(r, end);

% Get guessed boundary from the shooting method
shoot = @(s) boundary(forward_euler(f, Y0(s), dt, 0, 1), 2) - B;

s = bisection(shoot, 0.4, 1, @(a,b) abs(shoot(b)) < 1e-8);

%% Plot with final guess

[y, t] = forward_euler(f, Y0(s), dt, 0, 1); 
[y1, t1] = forward_euler(f, Y0(1), dt, 0, 1);
[y2, t2] = forward_euler(f, Y0(0.4), dt, 0, 1);

figure(2);
plot(t, y(2,:), t1, y1(2,:), '--', t2, y2(2,:), '--', 0, A, 'ro', 1, B, 'bo');
title(['Converged initial condition of [' num2str(s, '%d') '; -1]']);
legend('u_{converged}', 'u_{s = 1}', 'u_{s = 0.4}', 'u(0)', 'u(1)',...
       'Location', 'northwest');
xlabel('t');
ylabel('u(t)');
