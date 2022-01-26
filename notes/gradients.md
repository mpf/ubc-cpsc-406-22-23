<!-- @def reeval=true -->
# Gradients

```julia:setup
# hideall
using Random; Random.seed!(1234)
```

\blurb{Gradients provide information on a function's sensitivity to perturbations in the input.}


### Required Reading

- [Section 1.5.1](https://doi.org/10.1137/1.9781611973655.ch3) of Beck
- [Chapter 2](https://algorithmsbook.com/optimization/files/optimization.pdf) of Kochenderfer and Wheeler

## Directional derivatives

The behavior of a function $f:\R^n\to\R$ along the ray $\{x+αd\mid α\in\R_+\}$, where $x$ and $d$ are $n$-vectors, is given by the univariate function
\begin{equation*}
 \phi(\alpha) = f(x+αd).
\end{equation*}
From standard calculus, the derivative of $\phi$ at the origin, when it exists, is the limit
\begin{equation*}
  \phi'(0) = \lim_{α\to0^+}\frac{\phi(α)-\phi(0)}{α}.
\end{equation*}
We thus arrive at the following definition.
\definition{directional derivative}{
  The directional derivative of a function $f:\R^n\to\R$ at a point $x\in\R^n$, along a direction $d\in\R^n$, is the limit
  \begin{equation*}
    f'(x;d) = \lim_{α\to0^+}\frac{f(x+αd)-f(x)}{α}.
  \end{equation*}
}

It follows immediately from this definition that the partial derivatives of $f$ are simply the directional derivatives of $f$ along the each of the canonical unit direction $e_1,\ldots,e_n$, i.e.,
\begin{equation} \label{eq:partial-derivative}
  \frac{\partial f}{\partial x_i}(x) \equiv f'(x;e_i). 
\end{equation}

### Descent directions

A nonzero vector $d$ is a **descent direction** of $f$ at $x$ if the directional derivative is negative:
\begin{equation}\label{eq:descent-dir}
  f'(x;d) < 0.
\end{equation}
It follows directly from the definition of the directional derivative that $f(x+αd) < f(x)$ for all positive $\alpha$ small enough.

## Gradient vector

The gradient of the function $f$ is the collection of all the partial derivatives:
\begin{equation*}
  \nabla f(x) = \bmat{
	  \frac{\partial f}{\partial x_1}(x)
	\\\vdots
	\\\frac{\partial f}{\partial x_n}(x)}.
\end{equation*}
The gradient and directional derivative are related via the formula
\begin{equation*}
  f'(x;d) = \nabla f(x)\T d.
\end{equation*}
If, for example, the direction $d$ to be the canonical unit direction $e_i$, then this formula reduces to
\begin{equation*}
  f'(x;e_i) = \nabla f(x)\T e_i = [\nabla f(x))]_i = \frac{\partial f}{\partial x_i}(x),
\end{equation*}
which confirms the identity \eqref{eq:partial-derivative}.

### Linear approximation

The gradient of a continuously differentiable function $f$ (i.e., $f$ is differentiable at all $x$ and $\nabla f$ is continuous) provides a local linear approximation of $f$ in the following sense:
\begin{equation*}
   f(x+d) = f(x) + \nabla f(x)\T d + \omicron(\norm{d}),
\end{equation*}
The residual $\omicron:\R_+\to\R$ of the approximation decays faster than $\norm{d}$, i.e., \begin{equation*}\lim_{α\to0+}\omicron(α)/α=0.\end{equation*} 

### Example

Fortunately, there are [good computational](https://juliadiff.org/) tools that automatically produce reliable gradients. Consider the 2-dimensional [Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function) and its gradient:
\begin{align*}
     f(x)        &= (a - x_1)^2 + b(x_2 - x_1^2)^2
  \\ \nabla f(x) &= \bmat{-2(a-x_1)-4b(x_2-x_1^2)x_1 \\ 2b(x_2-x_1^2)}
\end{align*}

Here is the code for $f$ and its gradient:
```!
a, b = 1, 100
f(x) = (a - x[1])^2 + b*(x[2] - x[1]^2)^2
∇f(x) = [-2(a - x[1]) - 4b*(x[2] - x[1]^2)*x[1] , 2b*(x[2] - x[1]^2) ]
; #hide
```

Instead of computing gradients by hand, as we did above, we can use [automatic differentiation](https://duckduckgo.com/?t=ffab&q=automatic+differentiation&atb=v274-1&ia=web), such as implemented in the package [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl), to compute these.
```!
using ForwardDiff
∇fad(x) = ForwardDiff.gradient(f, x)
; #hide
```
The function `∇fad` returns the value of the gradient at `x`. Let's compare the hand-computed gradient `∇f` against that automatically-computed gradient `∇fad` at a random point:
```!
x = rand(2) 
∇f(x) == ∇fad(x)
```


<!-- \fig{rosenbrock} -->
<!-- x = y = -5:01:5
contour(x,y,(x,y)->f([x,y]))
savefig(joinpath(@OUTPUT,"rosenbrock")) # hide -->

## Calculus rules

We derive calculus rules for linear and quadratic functions, which appear often in optimization.

### Linear functions

 Let $a\in\R^n$. The linear function
\begin{equation*}
  f(x) = a\T x = \sum_i^n a_i x_i
\end{equation*}
has the gradient $\nabla f(x) = a$, and so the gradient is constant. Here's a small example:
```!
a = collect(1:5)
ForwardDiff.gradient(x->a'x, rand(5))
; #hide
```

### Quadratic functions

 Let $A\in\R^{n\times n}$ be a square matrix. Consider the quadratic function
\begin{equation}\label{eq:quadratic-fn}
  f(x) = \half x\T A x.
\end{equation}
One way to derive the gradient of this function is to write out the quadratic function making explicit all of the coefficients in $A$. Here's another approach that uses the product rule:
\begin{align*}
    ∇f(x) = \half ∇(x\T a) + \half ∇(b\T x),
\end{align*}
where $a = Ax$ and $b:=A\T x$ are held fixed when applying the gradient. Because each of the functions in the right-hand side of this sum is a linear function, we can apply the calculus rule for linear functions to deduce that
\begin{equation*}
  ∇f(x) = \half Ax + \half A\T x = \half (A+A\T)x.
\end{equation*}
(Recall that $A$ is square.) The matrix $\half(A+A\T)$ is the _symmetric part_ of $A$. If $A$ is symmetric, i.e., $A = A\T$, then the gradient reduces to 
\begin{equation*}
  ∇f(x) = Ax.
\end{equation*}

But in optimization, we almost always assume that the matrix that defines the quadratic in \eqref{eq:quadratic-fn} is symmetric because always $x\T Ax = \half x\T(A+A\T)x$, and therefore we can instead work with the symmetric part of $A$.

**Example** (2-norm). Consider the two functions
\begin{equation*}
f_1(x) = \norm{x}_2 \quad\text{and}\quad f_2(x) = \half\norm{x}_2^2.
\end{equation*}
The function $f_2$ is of the form \eqref{eq:quadratic-fn} with $A=I$, and so $\nabla f_2(x) = x$. Use the chain rule to obtain the gradient of the $f_1$:
\begin{equation*}
  \nabla f_1(x) = \nabla (x\T x)^\half = \half (x\T x)^{-\half}\nabla (x\T x) = \frac{x}{\norm{x}_2},
\end{equation*}
which isn't differentiable at the origin.

## Visualizing gradients

Gradients can be understood geometrically in relation to the level-set of the function. The $α$-level set of a function $f:\R^n\to\R$ is the set of points that have an equal or lower value at $x$:
\begin{equation*}
  [f≤α] = \{x\in\Rn\mid f(x)≤α\}.
\end{equation*}
Fix any $x$ and consider the level set $[f≤f(x)]$. For any direction $d$ that's either a descent direction for $f$ or a tangent direction for $[f≤f(x)]$,
\begin{equation*}
   f'(x;d) = ∇f(x)\T d ≤ 0,
\end{equation*}
which implies that the gradient $\nabla f(x)$ is the outward normal to the level set.

\fig{./gradient-normal}