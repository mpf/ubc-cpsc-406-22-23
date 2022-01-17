# Regularized Least Squares

```julia:setup
# hideall
using Random; Random.seed!(1234)
```

\blurb{Regularization introduces to an optimization problem prior information or assumptions about properties of an optimal solution. A regularization parameter governs the tradeoff between the fidelity to the unmodified problem and the desired solution properties.}

## The Pareto frontier

The standard [least-squares problem](least-squares) can be interpreted as a method for recovering some underlying object (say, a signal or image) encoded in the vector $x_0$ using noisy measurements of the form

\begin{equation} \label{eq:standard-ls}
   b := Ax_0 + w_1,
\end{equation}

where $w$ is random noise and the matrix $A$ describes the measurement process. The least-squares problem

\begin{equation*}
  \min_{x\in\R^n}\ \|Ax-b\|^2.
\end{equation*}

seeks a solution $x$ under the implicit assumption that the noise term $w_1$ is small.


But what if $x$ also needs to satisfy some other competing objective?  Suppose, for example, that we have also available other set of measurements of $x_0$ of the form $d = Cx_0+w_2$, where $w_2$ is also random noise. The "best" choice for $x$ is not necessarily \eqref{eq:standard-ls}, and instead we wish to choose $x$ to balance the two objective values

\begin{equation*}
f_1(x) = \|Ax - b\|^2 \quad\text{and}\quad f_2(x) = \|Cx - d\|^2.
\end{equation*}

Generally, we can make $f_1(x)$ or $f_2(x)$ small, but not both. The figure below sketches the relationship between the pair $\{f_1(x), f_2(x)\}$ for all values of $x$.

\fig{./pareto-curve}

The objective pairs on the boundary of two regions is the [Pareto frontier](https://en.wikipedia.org/wiki/Pareto_efficiency). We can compute these Pareto optimal solutions $x$ by minimizing the weighted sum objective

\begin{equation} \label{Regularized_LS_weights}
f_1(x)+\gamma f_2(x) = \|Ax-g\|^2+ λ\|Cx-d\|^2,
\end{equation}

where the positive parameter $\lambda$ provides the relative weight between the objectives. 

For fixed scalars $λ$ and $α$, the set 

$$
\ell(\lambda,\alpha) = \{\ (f_1(x),f_2(x)) \mid f_1(x) + \lambda f_2(x) = \alpha,\ x \in \R^n\ \}
$$ 

forms the graph of a line with slope of $-\lambda$. We may visualize the optimization problem \eqref{Regularized_LS_weights} as the problem of finding the smallest value $\alpha$ such that the line is tangent to the Pareto frontier.

\fig{./pareto-curve-levels}

## Tikhonov regularization

A particularly common example of regularized least-squares is [Tikhonov](https://en.wikipedia.org/wiki/Tikhonov_regularization), which has the form

\begin{equation*}
\min_{x} \half\|Ax-b\|^2 + \lambda\half\|Dx\|^2,
\end{equation*}

for which the objective can be expressed as

$$
\|Ax-b\|^2+\lambda\|Dx\|^2
 = \bigg\|\bmat{A\\\sqrt{\lambda}D} x - \bmat{b\\ 0}\bigg\|^2.
$$

If $D$ has full column rank, then the stacked matrix
$$
\bmat{A\\\sqrt{\lambda}D}
$$
necessarily also has full rank for any positive $\lambda$, which implies that the regularized problem always has a well-defined unique solution.

## Example: Signal denoising

Consider a noisy measurement

$$
b := x^\natural + w
$$

of a signal $x^\natural$, where the vector $w$ represents unknown noise. Here's a simple 1-dimensional noisy signal:

```julia:ls-reg-noisy1
using LinearAlgebra, Plots
n = 300
t = LinRange(0, 4, n)
x = @. sin(t) + t*cos(t)^2
w = 0.1*randn(n)
b = x + w
plot(b, leg=:topleft, label="noisy signal")
savefig(joinpath(@OUTPUT,"ls-reg-noisy1")) # hide
```
\fig{ls-reg-noisy1}

 The obvious least-squares problem

$$
\min_{x}\ \half\|x-b\|^2
$$

isn't useful because the optimal solution is simply the noisy measurements $b$, and so it doesn't yield any new information. But suppose we believe that the signal is "smooth" in the sense that the consecutive elements of $$ x=(x_1,\ldots,x_i,x_{i+1},\ldots,x_n)$$ change relatively little, i.e., the difference $|x_i-x_{i+1}|$ is small relative to $x$. In this case, we might balance the least-squares fit against the smoothness of the solution by instead solving the regularized least-squares problem 

\begin{equation}\label{Regularized_LS_identity}
  \min_{x}\ \half\underbrace{\vphantom{\sum_{i=1}}\|x-b\|^2}_{f_1(x)} + λ\half \underbrace{\sum_{i=1}^{n-1}(x_i - x_{i+1})^2}_{f_2(x)}.
\end{equation}

The role of the regularizer $f_2(x)$ is to promote smooth changes in the elements of $x$.

### Matrix notation

We can alternatively write the above minimization program in matrix notation. Define the $(n-1)$-by-$n$ finite difference matrix

$$
D = \bmat{ 1 & -1 & \phantom+0 & \cdots & \phantom+0 & \phantom+0\\
           0 & \phantom+1 & -1 & 0 & \cdots & \phantom+0\\
             &\ddots  & \ddots & \ddots &  &  \\
           \vdots &  & \phantom+0 & 1 & -1 & \phantom+0\\
            0 & \cdots &  & 0 & \phantom+1 & -1},
$$
which when applied to a vector $x$, yields a vector of the differences:

$$
Dx = \bmat{ x_1 - x_2 \\ x_2 - x_3 \\ \vdots \\ x_{n-1} - x_n }.
$$

Then we can rephrase the regularization objective as $$ f_2(x) = \sum_{i=1}^{n-1}(x_i - x_{i+1})^2 = \|Dx\|^2. $$ This allows for a reformulation of the weighted leas squares objective into a familiar least squares objective:

$$
\|x-b\|^2+\lambda\|Dx\|^2 = \bigg\|\underbrace{\bmat{I\\\sqrt{\lambda}D}}_{\hat{A}}x - 
\underbrace{\bmat{b\\ 0}}_{\hat{b}}\bigg\|^2.
$$

So the solution to the weighted least squares minimization program \eqref{Regularized_LS_identity} satisfies the normal equation $\hat{A}\T\hat{A}x = \hat{A}\T\hat{b}$, which simplifies to 

$$
(I + \lambda D\T D)x = b.
$$

```julia:ls-reg-noisy2
finiteDiff(n) = diagm(ones(n)) - diagm(+1 => ones(n-1))
finiteDiff(4)
```
\show{ls-reg-noisy2}

```julia:ls-reg-makeD
D = finiteDiff(n)
```

```julia:ls-reg-noisy3
λ = 100
plot(t, b, leg =:topleft, label="noisy data")
b̂ = [ b; zeros(n)]
for λ in LinRange(100,5000,5) 
    Â = [ I; √λ*D ]
    xLS = Â \ b̂
    plot!(t, xLS, label="λ = $(λ)")
end
savefig(joinpath(@OUTPUT,"ls-reg-noisy3")) # hide
```
\fig{ls-reg-noisy3}
