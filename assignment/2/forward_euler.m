function [y, t] = forward_euler(f, y, dt, t, T)
%FORWARD_EULER Solve an IVP or system of IVPs
%   y must be a scalar or column vector; if y is a column vector, then
%   f(y, t) must return a vector of the same size.
%
%   [y, t] = FORWARD_EULER(...) returns a row vector of t values used, and
%   y is a matrix with the same number of columns as the original y0; thus,
%   if y0 is scalar then y is a row vector.

n = 1;
while t < T-dt
    t(n+1) = t(n) + dt;
    y(:,n+1) = y(:,n) + dt * f(y(:,n), t(n));
    n = n + 1;
end

end