x = linspace(0, pi, 1000);

phi1 = @(x) sin(x);
phi2 = @(x) sin(2*x);
phi3 = @(x) sin(3*x);
phi = [phi1(x) ; phi2(x) ; phi3(x)];

exact = @(x) x.*(pi-x);
approx = @(x, c1, c2, c3) phi1(x)*c1 + phi2(x)*c2 + phi3(x)*c3;

%% Plot the basis functions

c1 = 8/pi;
c2 = 0;
c3 = 8/(27*pi);

figure;
plot(x/pi, phi);
title("Basis functions of u^h");
legend("\phi_1", "\phi_2", "\phi_3", "Location", "southwest");
ylabel("y");
xlabel("x\pi");

figure;
plot(x/pi, [c1; c2; c3].*phi);
title("c_i\phi_i");
legend("c_1\phi_1", "c_2\phi_2", "c_3\phi_3");
ylabel("y");
xlabel("x\pi");


%% Compare the approximate $u^h$ to $u$

figure;
plot(x/pi, [exact(x); approx(x, c1, c2, c3)]);
title("u vs u^h");
legend("u", "u^h");
ylabel("y");
xlabel("x\pi");