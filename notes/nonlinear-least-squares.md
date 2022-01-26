@def reeval=true
# Nonlinear Least-Squares

```julia:setup
# hideall
using Random; Random.seed!(12)
```

\blurb{The nonlinear least-squares problem generalizes the linear least-squares problem to include nonlinear residual functions. The Gauss-Newton solves this problem as a sequence of least-squares subproblems.}

### Required Reading

- Sections [3.5](https://doi.org/10.1137/1.9781611973655.ch3) and [4.5](https://epubs.siam.org/doi/abs/10.1137/1.9781611973655.ch4) of [Beck](https://epubs.siam.org/doi/book/10.1137/1.9781611973655)


## Nonlinear residuals

The nonlinear least squares problem seeks to minimize the sum of squares

\begin{equation}\label{eq:nonlinear-ls}
  \min_{x\in\R^n}\ \half\|r(x)\|^2,
\end{equation}

where the $m$-vector of residuals
\begin{equation*}
  r(x)=\bmat{r_1(x)\\\vdots\\r_m(x)}
\end{equation*}
is composed of differentiable nonlinear functions $r_{i}:\Rn→\R$. The non-linear least squares problem reduces to linear least-squares when $r$ is affine, i.e., $r(x) = Ax-b$ for some $m$-by-$n$ matrix $A$ and $m$-vector of observations $b$.

## Linearizing the residual

We can solve non-linear least squares problem \eqref{eq:nonlinear-ls} by solving a sequence of linear least-squares problem, which result from linearization of $r$ at the current estimate of the ground truth $x$. 

The linear approximation of $r$ at a point $x^k \in \Rn$ is defined by

\begin{align*}
r^k(d) 
    &:= \bmat{
    r_1(x^k)+\nabla r_1(x^k)\T d \\ \vdots \\ r_m(x^k)+\nabla r_m(x^k)\T d}
\\  &= r(x^k) + J(x^k)d 
\end{align*}
where
$$
J(x)
= \bmat{
    \nabla r_1(x)^T\\ \vdots \\ \nabla r_m(x)^T}
$$
is the Jacobian of $r$ evaluated at $x$

\algorithm{Gauss-Newton method}{
\\ Choose a starting $x^{(0)}$\\
for $k=0,1,2,\ldots$
1. linearize $r$ at $x^{k}$ to define $r^k(d):=r(x^k)+J(x^k)d$
2. solve linear least-squares problem
   \begin{equation*}
     d^k = \argmin_d\ \half\norm{r^k(d)}^2
   \end{equation*} 
3. update iterate: $x^{k+1} = x^k + αd^k$ for some $α∈(0,1]$
}
If the linesearch parameter $\alpha$ is held fixed at 1, then this is a "pure" Gauss-Newton iteration. We'll later discuss options for when it's advisable to use smaller steps.

Here's a basic version of the method that takes as arguments the residual and Jacobian functions `r` and `J`, and a starting vector `x`.
```!
function gauss_newton(r, J, x; ε=1e-4, maxits=30)
	err = []
	for i = 1:maxits
        rk, Jk = r(x), J(x)
		push!(err, norm(Jk'rk))
		err[end] ≤ ε && break
        x = x - Jk\rk 
    end
    return x, err
end
; # hide
```
The condition `err[end] ≤ ε` causes the iterations to terminate when the latest residual `rk` is nearly orthogonal to all the gradients of the residual, i.e., the residual is orthogonal to the tangent space of the level-set for the nonlinear least-squares objective $\half\|r(x)\|^2$. This is exactly analogous to the optimality condition for linear least-squares.

## Example: Position estimation from ranges

\fig{./position-distances}

Let $x \in \R^2$ represent the unknown position of an object. Our aim is to find the position of this object relative to a set of $m$ beacons placed in fixed known positions $b_{i} \in \R^2$, $i = 1,\dots,m$. The only data available are range measurements $δ_i$ that give an estimate of the distance between each beacon $b_{i}$ and the object at $x$, i.e.,
$$
δ_i := \norm{x-b_i} + ν_i
$$
where each scalar $ν_i$ represents the error between the true (and unkown) distance $\norm{x-b_i}$ and the measurement $δ_i$.

We can obtain an estimate of the object's position $x$ by solving the nonlinear least-squares problem \eqref{eq:nonlinear-ls} where we define the $i$th residual between the reported distance $δ_i$ and the distance between the position $b_i$ of the $i$th beacon and $x$:

\begin{equation*}
  r_i(x) := δ_{i} - \|x-b_i\|.
\end{equation*}

Here's the residual function, which takes a vector of ranges $δ$ and a vector of positions `b`:
```!
r(x, δ, b) = [ δ[i] - norm(x - b[:,i]) for i in 1:length(δ) ]
; #hide
```

The following Julia packages for this example.
```!
using Random
using Statistics
using LinearAlgebra
using ForwardDiff
using Plots
```

We simulate data by placing `m` beacons in random locations:
```!
m = 13                 # number of beacons
b = 2rand(2, m) .- 1   # place beacons uniformly in the unit box 
x = zeros(2)           # true position 
x0 = 2rand(2) .- 1     # initial guess of the unknown position
ν = .5*rand(m)
δ = [ norm(x - b[:,i]) + ν[i] for i in 1:m]
; # hide
```

Place these items on a map (we'll wrap this in a function so that we can reuse it below):

```!
function plotmap(b, x, x0)
scatter(xlim=(-1,+1), ylim=(-1,+1), leg=:outertopright,frame=:box, aspect_ratio=:equal)
scatter!(b[1,:], b[2,:], label="beacons", shape=:square, ms=7)
scatter!(x[1,:], x[2,:], label="true position", shape=:xcross, ms=7)
scatter!(x0[1,:], x0[2,:], label="initial guess", c="yellow", ms=7)
end
plotmap(b, x, x0)
savefig(joinpath(@OUTPUT,"map")) # hide
```
\fig{map}

```!
J(x) = ForwardDiff.jacobian(x->r(x, δ, b), x)
xs, err = gauss_newton(x->r(x, δ, b), J, x0)
plot(err, yaxis=:log)
savefig(joinpath(@OUTPUT,"convergence")) # hide
```
\fig{convergence}

Plot the original map and overlay the obtained solution:
```!
plotmap(b, x, x0)
scatter!(xs[1,:], xs[2,:], label="solution", shape=:star7, c="green", ms=7)
savefig(joinpath(@OUTPUT,"map-soln")) # hide
```
\fig{map-soln}
