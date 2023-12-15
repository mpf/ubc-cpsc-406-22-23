@def title = "Descent Methods"

```julia:setup
# hideall
using Random; Random.seed!(12)
```

# Descent Methods

\blurb{Descent methods for unconstrained optimization follow at each iteration a search direction related to the current gradient. There are a variety of methods for generating valid descent directions, including steepest descent, scaled descent, coordinate descent, and variations of Newton descent.}

The methods we consider here minimize a smooth function $f:\R^n\to\R$ by generating a sequence of iterates
$$
x_{k+1} = x_k + α_k d_k
$$
that monotonically reduces the function value, i.e.,
$$
f(x_{k+1}) < f(x_k) \quad\text{for all}\quad k\in\N.
$$
The positive scalars $\alpha_k$ are **steplengths** and the vectors $d_k\in\R^n$ are [descent directions](/notes/gradients/#descent_directions). Our aim is to ensure that the generated sequence converges to a stationary point:
$$
\lim_{k\in\N}\|\nabla f(x_k)\| = 0.
$$

Below we describe methods for computing the descent directions and the steplengths.

## Running example

In our discussion we'll often use the [Rosenbrock test function](https://en.wikipedia.org/wiki/Rosenbrock_function) as a convenient example. Here's the function definition and its gradient and Hessian: 

```!
using LinearAlgebra
import ForwardDiff:gradient,hessian

f(x) = (1 - x[1])^2 + 100.0(x[2] - x[1]^2)^2
∇f(x) = gradient(f, x)
∇²f(x) = hessian(f, x)
;#hide
```

The minimizer is $x^*=(1,1)$, which we verify by checking that the gradient is zero and the Hessian is positive definite at that point:

```!
x̄ = [1., 1.]
@show norm(∇f(x̄))
@show isposdef(∇²f(x̄))
```

Here's a visualization of the function and the minimizer marked with a star: 
```julia:./code/descent
using Plots; gr()

x = y = -0.25:0.001:1.5
z = Surface((x,y)->f([x,y]), x, y)
contourf(x,y,z, levels=collect([10^k for k ∈ -4.0:3.0]), colorbar=:none, leg=false)
scatter!([1.],[1.], c="gold", ms=10, shape=:star5)
savefig(joinpath(@OUTPUT,"surface")) # hide
```
\fig{surface}

## Gradient-related directions

Observe that for any sequence of positive-definite matrices $\{B_k\}_{k\in\N}$, the vectors
$$
d_k = -B_k\nabla f(x_k)
$$
are descent directions for $f$ at points $x_k$ because the directional derivative
$$
f'(x_k;d_k) = ∇f(x_k)^T d_k = -∇f(x_k)\T B_k ∇f(x_k)<0
$$
whenever $f$ is not stationary at $x_k$, i.e., $\nabla f(x_k)\ne0$. This fact means that we can generate a whole family of descent directions through different choices the scaling matrices $B_k$. 

Stronger conditions are needed, however, to ensure that the directional derivative doesn't get arbitrarily close to zero. In particular, we need additional conditions to ensure that
$$
\limsup_{k\in\N}f'(x_k;d_k)<0.
$$
One approach is to require that the eigenvalues of the matrices $B_k$ are uniformly positive and bounded. Thus, we assume that there exist positive constants $\mu_1$ and $\mu_2$ such that
$$
0 < \mu_1 < \lambda_{\text{min}}(B_k) \le \lambda_{\text{max}}(B_k) < \mu_2 \quad \forall k\in\N,
$$
where $\lambda_{\text{min}}(B_k)$ and $\lambda_{\text{max}}(B_k)$, respectively, are the minimum and maximum eigenvalues of $B_k$. This condition implies that there are bounds on the directional gradient and the norm of the direction directly related to the gradient:
\begin{align}
  f'(x_k;d_k) &= -∇f(x_k)^T B_k∇f(x_k) < -μ_1∥∇f(x_k)∥^2 \\
  ∥d_k∥ &= ∥{-}B_k∇f(x_k)∥≤∥B_k∥⋅∥∇f(x_k)∥<μ_2∥∇f(x_k)∥.
\end{align}
Thus, the directional derivative of a gradient-related direction vanishes only when at stationary point, and the norm of a gradient-related direction is bounded by a multiple of the norm of the gradient.

Below we discuss the most-used search directions.


### Steepest descent direction

The negative gradient $d:=-∇f(x)$ is evidently a descent direction, because
$$
f'(x_k;-∇f(x_k)) = -∇f(x_k)^T ∇f(x_k) = -\|∇f(x_k)\|^2 < 0
$$
whenever $∇f(x_k)\ne0$. The negative gradient is often called the **steepest descent** direction because it's aligned with the unit-norm direction that makes the directional derivative most negative. In particular, use the 
[Cauchy-Scharz inequality](https://en.wikipedia.org/wiki/Cauchy%E2%80%93Schwarz_inequality) to deduce the bound
$$
  f'(x_k;v) = ∇f(x_k)^T v \ge -∥∇f(x_k)∥_2⋅∥v∥_2.
$$
for any direction $v$. Then observe that the direction $v=-∇f(x_k)$ achieves the lower bound.

The package [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl/blob/master/docs/src/index.md) implements the `gradient` function, which we negate and assign to `steepestDescent`:

```!
steepestDescent(f, x) = -gradient(f, x)
;#hide
```

### Newton direction

The Newton direction