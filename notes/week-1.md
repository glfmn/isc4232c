# Floating Point Arithmetic

\$$
-1^{sign}\(\sum{52}{i=1}b_{52 - i}^{-i}\)2^{exponent - 1023}
\$$

Given a register with 1 sign bit, 4 exponent bits, and 5 fraction bits, can we exactly represent 20?

Sign is 0 because 20 is positive

How do we determine fraction?
fraction can only represent  $[0,1)$ so we want to find out how to divide $20$ by a power of $2$ such that it fits that interval.

$$
\frac{20}{2^5} = \frac{20}{32} = 0.625 = \frac{1}{2} + \frac{1}{2^3}
$$

The corresponding binary is $10100$

## Error

### Absolute error

If $p^*$ approximates $p$:

$$
\|p - p^*\|
$$

### Relative error

If $p^*$ approximates $p$:

Tells you if you are doing a better job approximating the sun than the door, but be careful not to divide by zero!

$$
\frac{\|p - p^*\|}{\|p\|}
$$

### Roundoff error

Relative error between the number we could represent by the number we wanted to represent.

#### Subtraction

Suppose a machine with 16 bits for the fraction:

$$x = 1011011011011011$$
$$y = 1011011011011001$$

The roundoff error for each is:

$$x: 1.15x10^{-6}$$
$$y: 3.38x10^{-6}$$

But the relative error of subtracting these two values is $0.0523$, which means the subtraction has only two digits accuracy.  If we repeatedly take an action and loose a digit of accuracy each time, we very quickly loose all meaning in the numbers.

# Taylor Series

Suppose $\f : \[a,b\] \arrow \BB{R}$ is $n$-times differentiable at
$x_0 \in \[a,b\]$.  The Taylor series, or $n^{th}$-order polynomial, is:

$$
P(x) = T_n(x) = f(x_0) + f^'(x_0)(x - x_0)
              + \frac{f^{''}(x_0)}{\factorial{2}}(x - x_0)^2
              + ...
              + \frac{f^{\[n\]}{\factorial{n}}(x - x_0)^n
$$

# Newton's Method

Derived from the taylor expansion:

```
$$
g(x) = g(x_n) + g^'(x_n)(x-x_n) + g^{''}(x_n)\frac{(x-x_n)^2}{2!}
$$
```

We can create a formula for x* such that $g(x^*) = 0$ through substitution and a
little bit of algebra.

```
$$
x^* = x_n - \frac{g(x_n)}{'(x_n)}
$$
```

Note: we can iteratively approximate $x^*$ using the formula:

```
$$
x_{n+1} = x_n - \frac{g(x_n)}{g^'(x_n)}
$$
```

Graphically, you could interpret this as following the tangent line to the
x-axis, and using that x-coordinate as ${x_n+1}$

- Newton's method is sensitive to the initial guess
- Newton's method can oscillate between a few points which are not the roots.
- Newton's method does extend to higher dimensions
- Computing derivatives (or Jacobians in higher dimensions) is not always
  practical.

# Secant Method

Useful when the derivatives or Jacobian are not practical to calculate.


```
$$
x_{n+1} = x_n - \frac{x_n - x_{n-1}}{g(x_n)-g(x_{n-1})}
$$
```

# Lagrange Interpolation

The Taylor expansion is only accurate around the center of expansion; this means
we can potentially introduce error as we approximate roots away from the center
of interpolation.

If we use a set of points $x_1,...,x_n\in\[a,b\]$

```
$$
L_k(x) := \product^n_{i=1 where L \ne k}\frac{x-x_i}{x_k-x_i} for k=1,...,n
$$

$$
L_k(x_i) = 0
L_k(x_j) = 1
$$
```

Our set of points do not need to be spaced equally, but they must be unique.

The location of these nodes plays a central role in minimizing the error between
the lagrange approximation and the function.

Runge phenomenon: large oscillation in error as we approach the boundaries of
the function.

We can match not just the location of the points, but we can also match the

# Numerical Differentiation

Using Lagrange interpolation:

```
$$
p(x) = f(x_1)L_1(x)+f(x_2)L_2(x)
     = -\frac{f(x_0)}{h}(x-x_h-h)+\frac{f(x_0+h)}{h}(x-x_0)
$$
```
The error is the same as the previous Taylor interpolation so the error is
$O(h)$, or the smaller h is, the smaller the error is.

```
$$
f=1=>\frac{f(x_0+h)-f(x_0)}{h} = \frac{1-1}{h} = 0
f=x=>\frac{f(x_0+h)-f(x_0)}{h} = \frac{x_0+h-x_0}{h} = \frac{h}{h} = 0
f=x^2=>{f(x_0+h)-f(x_0)}{h} = \frac{(x_0+h)^2-x_0^2}{h} = 2x_0+h
$$
```

The approximation is always exact for linear polynomials but may include error
for nonlinear polynomials.
