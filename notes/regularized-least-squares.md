# Regularized Least Squares

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

tries to remove the noise term $w$.


But what if $x$ also needs to satisfy some other competing objective?  Suppose, for example, that we have also available other set of measurements of $x_0$ of the form $d = Cx_0+w_2$, where $w_2$ is also random noise. The "best" choice for $x$ is not necessarily \eqref{eq:standard-ls}, and instead we wish to choose $x$ to balance the two objective values

\begin{equation*}
f_1(x) = \|Ax - b\|^2 \quad\text{and}\quad f_2(x) = \|Cx - d\|^2.
\end{equation*}

Generally, we can make $f_1(x)$ or $f_2(x)$ small, but not both. The figure below sketches the relationship between the pair $\{f_1(x), f_2(x)\}$ for all values of $x$.

\fig{./pareto-curve}

The objective pairs on the boundary of two regions is the [Pareto frontier](https://en.wikipedia.org/wiki/Pareto_efficiency). We can compute these Pareto optimal solutions $x$ by minimizing the weighted sum objective

\begin{equation} \label{Regularized_LS_weights}
f_1(x)+\gamma f_2(x) = \|Ax-g\|_2^2+ λ\|Fx-g\|_2^2,
\end{equation}

where the parameter positive parameter $\gamma$ provides the relative weight between the objectives. 

For fixed parameters $λ$ and $α$, the set 

$$
\ell(\gamma,\alpha) = \{(f_1(x),f_2(x)):f_1(x) +\gamma f_2(x) = \alpha, x \in \R^n\}
$$ 

forms the graph of a line with slope of $-\gamma$. One way to visualize the optimization problem \eqref{Regularized_LS_weights} is that finds the smallest value $\alpha$ such that the line is tangent to the Pareto frontier.

\fig{./pareto-curve-levels}

## Example: Signal denoising

Consider a noisy measurement

$$
b = \hat{x} + w,
$$

of a signal $\hat x$.

```julia:noisy
using LinearAlgebra, Plots
n = 300
t = LinRange(0, 4, n)
x̂ = @. sin(t) + t*cos(t)^2
w = 0.1*randn(n)
b = x̂ + w
plot(b, label="b")
savefig(joinpath(@OUTPUT,"noisy")) # hide
```
\fig{noisy}

 The vector $w$ represents unknown noise. The obvious least-squares problem

$$
\min_{x\in\R^n} \frac{1}{2}\|x-b\|_2^2.
$$

isn't useful because the optimal solution just gives us back the measurements $b$. But suppose we believe that the signal is "smooth" in the sense that the elements of $$ x=(x_1,\ldots,x_i,x_{i+1},\ldots,x_n)$$ change relatively little, i.e., the difference $|x_i-x_{i+1}|$ is small relative to $x$. In this case, we might balance the least-squares fit against the smoothness of the solution by instead solving the regularized least-squares problem 

\begin{equation}\label{Regularized_LS_identity}
  \min_{x}\ \half\underbrace{\vphantom{\sum_{i=1}}\|x-b\|^2}_{f_1(x)} + λ\half \underbrace{\sum_{i=1}^{n-1}(x_i - x_{i+1})^2}_{f_2(x)}.
\end{equation}

The role of the regularizer $f_2(x)$ is to promote smooth changes in the elements of $x$. We can alternatively write the above minimization program in matrix notation. Define the $(n-1)$-by-$n$ finite difference matrix

### Matrix notation

If we let

$$
D = \bmat{ 1 & -1 & \phantom+0 & \cdots & \phantom+0 & \phantom+0\\
           0 & \phantom+1 & -1 & 0 & \cdots & \phantom+0\\
             &\ddots  & \ddots & \ddots &  &  \\
           \vdots &  & \phantom+0 & 1 & -1 & \phantom+0\\
            0 & \cdots &  & 0 & \phantom+1 & -1},
$$

we can then rephrase the regularization objective as $$ f_2(x) = \sum_{i=1}^{n-1}(x_i - x_{i+1})^2 = \|Dx\|_2^2. $$ This allows for a reformulation of the weighted leas squares objective into a familiar least squares objective:

$$
\|x-b\|_2^2+\gamma\|Dx\|_2^2 = \bigg\|\underbrace{\bmat{I\\\sqrt{\gamma}D}}_{\hat{A}}x - 
\underbrace{\bmat{b\\ 0}}_{\hat{b}}\bigg\|^2.
$$

So the solution to the weighted least squares minimization program \eqref{Regularized_LS_identity} satisfies the normal equation $\hat{A}\T\hat{A}x = \hat{A}\T\hat{b}$, which simplifies to 

$$
(I + \gamma D\T D)x = b.
$$

## **Regularized least squares (aka Tikhonov)**

We now generalize the result to noisy linear observations of a signal. In this case, the model is

$$
b = Ax + w,
$$

where we added the measurement matrix $A \in \R^{m\times n}$. The  corresponding weighted-sum least squares
program is 

\begin{equation}
\min_{x\in\R^n} \frac{1}{2}\|Ax-b\|_2^2 + \frac{\lambda}{2}\|Dx\|_2^2,
\end{equation}

<!-- where $\|Dx\|_2^2$ is called the regularization penalty and $ \gamma $ is called the regularization 
parameter. The objective function can be reformulated as an least squares objective -->

$$
\|Ax-b\|_2^2+\gamma\|Dx\|_2^2
 = \bigg\|\underbrace{\bmat{A\\\sqrt{\gamma}D}}_{\hat{A}}x - 
\underbrace{\bmat{b\\ 0} }_{\hat{b}}\bigg\|_2^2.
$$

and the corresponding normal equations is

$$
(A\T A + \lambda D\T D)x = A\T b.
$$
