function [ x, err, its ] = bisection( f, a, b, converge )
%BISECTION finds roots using the bisection method.
%   Evaluates function f over closed interval [a,b] to find the root,
%   requires that [a,b] change signs, otherwise a root may or may not
%   exist.
%
%   x = BISECTION(f,a,b) finds a root x of f in the interval [a,b]
%   BISECTION(f,a,b,converge) finishes its iterations when converge(a,b)
%   evaluates to true.
%   [ x, err, its ] = BISECTION(...) returns the most recent absolute error
%   term err, along with the number of iterations its.
%
%   Convergence criteria default to:
%
%       @(a,b) abs(a-b) < 10^-8
%
%   The function will also terminate if a-b == 0 evaluates to true, meaning
%   that the answer is known beyond machine epsilon (assuming a valid
%   interval [a,b]).

if nargin < 4, converge = @(a,b) abs(a-b) < 10^-8; end

if f(a)*f(b) > 0, error('f(a), f(b) must change signs.'); end

next = @(a,b) (a+b)/2;
its  = 0;

while ~converge(a,b) % Iterate until interval is sufficiently small
    % determine approximate error to help check convergence
    if its > 0, err = abs( (next(a,b)-x) / x ); end

    its = its + 1;

    x = next(a,b); % Claculate midpoint of the interval
    choice = f(a)*f(x); % Determine sign difference at a and xst

    if (choice <  0), b = x; end % Root is in [a, xst]
    if (choice == 0), break; end % Exit, root was found in macheps
    if (choice >  0), a = x; end % Root is in [xst, b]
end

end