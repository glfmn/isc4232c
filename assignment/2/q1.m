%% Visual verification

% Define the problem for the exact solution $u(x) = e^{cos(x)}$
u = @(x) exp(cos(x));
u_xx = @(x) exp(cos(x)).*(sin(x).^2 - cos(x));
f = @(x) u(x) - u_xx(x);
a = 0;
b = 2*pi;
n = 100;

[U, xs] = fd_periodic(f, a, b, n);

figure(1);
plot(xs, U, xs, u(xs), '--');
title(['Solution using finite difference for N = ' num2str(n, '%d')]);
legend("approximate solution", "exact solution");

%% Convergence Study

% Exponents of 2 for the number of discretization points
exps = 3:12;

% Expression for the approximate solution
approx = @(N) fd_periodic(f,a,b,N);

% Expression for the exact solution
dt = @(a,b,N) (b-a)/(N+1);
exact = @(N) u(linspace(a+dt(a,b,N), b-dt(a,b,N),N))';

% L2 Norm
l_2 = @(N) 1/sqrt(N) * sqrt(sum((exact(N)-approx(N)).^2));

% Calculate the errors
N = repmat(2,[1,numel(exps)]).^exps;
E_2 = arrayfun(l_2, N);

% Plot the error in loglog scale
figure(2);
loglog(N, E_2);
title("Convergence Study");
legend("L_2 norm");
xlabel("N_{elements}");
ylabel("Error");