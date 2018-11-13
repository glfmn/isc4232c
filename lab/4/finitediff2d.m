function [U,x,y] = finitediff2d(f, N, br, bl, bt, bb)
%FINITEDIFF2D Summary of this function goes here
%   Detailed explanation goes here

% Create our set of discretization points
h = 1/(N+1);
x = linspace(h,1-h,N);
y = x;
[X, Y] = meshgrid(x,y);

rhs = -h^2*f(X,Y);
rhs(1:N,1) = rhs(1:N,1) + bl(y');
rhs(1:N,N) = rhs(1:N,N) + br(y');
rhs(1,1:N) = rhs(1,1:N) + bt(x);
rhs(N,1:N) = rhs(N,1:N) + bb(x);

D = [-ones(N^2,1) -d(N,N) 4*ones(N^2,1) -d(N,1) -ones(N^2,1)];
A = spdiags(D, [-N -1 0 1 N], N^2, N^2);

rhs = reshape(rhs, [N^2 1]);
x = reshape(X, [N^2 1]);
y = reshape(Y, [N^2 1]);

U = A\rhs;

end

function [d] = d(N,n)
% Create diagonals with a pattern where the nth element is a 0 in a
% repeating sequence of N vectors of N elements.  The resulting number of
% elements will be N^2

d = ones(N,1);
d(n) = 0;
d = repmat(d,N,1);

end
