function [x] = thomas(a,b,c,d)
%THOMAS Solve a tri-diagonal system in O(N) time

N = numel(b);
x = zeros(N,1);

c(1) = c(1)/b(1);
c(2:N-1) = arrayfun(@(c2, c1, b, a) c2/(b - a*c1), c(2:N-1), c(1:N-2), b(2:N-1), a(2:N-1));

d(1) = d(1)/b(1);
for k = 2:N
    d(k) = (d(k) - a(k) * d(k-1))/(b(k) - a(k)*c(k-1));
end

x(N) = d(N);
for k = N-1:-1:1
    x(k) = d(k) - c(k)*x(k+1);
end

end

