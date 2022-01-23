# QR Factorization

\blurb{The QR factorization of a matrix constructs an orthgonal basis for its columnspace. It's also one of the best computational tools for solving least-squares problems.}

```!
# hideall
import Random; Random.seed!(1234);
```

### Required Reading

- [Orthogonal projection](https://ubcmath.github.io/MATH307/orthogonality/projection.html) and [QR decomposition](https://ubcmath.github.io/MATH307/orthogonality/qr.html) of UBC Math 307 [lecture notes](https://ubcmath.github.io/MATH307)

\note{This page uses the `LinearAlgebra` package, which contains `norm` and `qr`.
```!
using LinearAlgebra
```
}

## Orthogonal and orthonormal vectors

For any two $n$-vectors $x$ and $y$, the cosine identity
\begin{equation*}
  x\T y = \|x\|\|y\|\cos(\theta)
\end{equation*}
relates the inner product of the vectors to the the angle between them. Thus, the vectors $x$ and $y$ are **orthogonal** when
\begin{equation*}
	x\T y = 0, \quad\text{or equivalently,}\quad \cos(\theta) = 0.
\end{equation*}
Orthogonal vectors are **orthonormal** when they're also normalized:
\begin{equation*}
	x\T y = 0, \quad \norm{x}=1, \quad \norm{y}=1.
\end{equation*}

These definitions extend readily to a collection of vectors: the set of $k$ vectors $\{x_1,\ldots,x_k\}$ in $\R^n$ is orthogonal (and orthonormal) whenever
\begin{equation*}
	x_i\T x_j = 0 \quad\text{for all}\quad i\ne j
\qquad (\text{and}\ \norm{x_i} = 1)
\end{equation*}

The 2-dimensional canonical unit vectors $e_1 = (1,0)$ and $e_2 = (0,1)$ are orthonormal:
```!
e₁ = [1, 0];  e₂ = [0, 1]
@show e₁'e₂
@show norm(e₁)
@show norm(e₂)
```

## Orthogonal matrices

A square $n$-by-$n$ matrix $Q$ is **orthogonal** when its columns are orthonormal:
\begin{equation*}
   Q = \begin{bmatrix} q_1 & q_2 & \cdots & q_n \end{bmatrix}
   \quad\text{and}\quad
   Q\T Q = QQ\T = I.
\end{equation*}
Because a matrix $B$ is the inverse of a matrix $A$ if and only if $BA =AB = I$, the inverse of an 
orthogonal matrix is its transpose:
\begin{equation*}
Q^{-1} = Q\T.
\end{equation*}

An orthogonal matrix rotates a vector, but doesn't change its length. Thus inner products are invariant under orthogonal transformations, i.e., for any vectors $x$ and $y$,
\begin{equation*}
(Qx)\T(Qy) = x\T Q\T Qy = x\T y,
\end{equation*}
and so
\begin{equation*}
  \| Q  x \|_2 = \sqrt{(Qx)\T (Qx)} = \sqrt{x\T x} = \| x \|_2.
\end{equation*}

Another property of orthogonal matrices is that their determinant is either $1$ or $-1$. This can be observed from the fact that for any compatible matrices $ A $ and $ B $, $\det( A  B ) = \det( A )\det( B )$ and $\det( A )= \det( A^T )$. Hence, 
\begin{equation*}
\det( Q \T  Q ) = \det( I ) \iff \det( Q )^2 = 1 \iff \det( Q ) = ±1.
\end{equation*}

## QR factorization

Every $m$-by-$n$ matrix $A$ has a **QR factorization**, which means it can be decomposed as the product
\begin{equation*}
  A  =  Q   R, 
\end{equation*}
where
\begin{equation*}
  Q = \begin{bmatrix} q_1 & \cdots & q_m \end{bmatrix}
  \textt{and}
  R = \begin{bmatrix} r_{11} & r_{12} & \cdots & r_{1n}
		   \\      0 & r_{22} & \cdots & r_{2n}
		   \\ \vdots &        & \ddots & \vdots
		   \\     0  &     0  & \cdots & r_{mn}
  \end{bmatrix},
\end{equation*}
respectively, are an $m$-by-$m$ orthogonal matrix and an $m$-by-$n$ upper-triangular matrix. This factorization isn't unique, unless we impose additional conditions on the factors $Q$ and $R$.

[Householder reflectors](https://en.wikipedia.org/wiki/Householder_transformation) and [Givens rotations](https://en.wikipedia.org/wiki/Givens_rotation) are the preferred methods for computing the QR factorization. These methods require $\Omicron(m^2n)$ flops to complete.

The columns of the orthogonal matrix $Q$ reveal an orthogonal basis for the range of $A,$ and so triangularity of $R$ gives a simple "recipe" for reconstructing the columns of $A=[a_1,\ldots,a_n]$ from those of $Q$:
\begin{align*}
   a_1 &= r_{11} q_1
\\ a_2 &= r_{12} q_1 + r_{22} q_2
\\ a_3 &= r_{13} q_1 + r_{23} q_2 + r_{33} q_3
\\     &\ \ \vdots
\\ a_n &= r_{1n} q_1 + r_{2n} q_2 + \cdots + r_{nn} q_n.
\end{align*}
From these formulas we deduce a key property of this factorization: for any index $k < n$, the span of the leading $k$ columns of $Q$ contain the span of the leading $k$ columns of $A$, i.e.,
\begin{equation*}
  \span\{a_1,\ldots,a_k\} \subseteq \span\{q_1,\ldots,q_k\}.
\end{equation*}
If we further assume that $A$ has full column rank, this inclusion tightens to an equality, which tells us that leading $k$ columns of $Q$ provide a tight basis for the corresponding columns.

\theorem{QR factorization for full column rank matrices}{If $A=QR$ is the QR factorization of an $m$-by-$n$ matrix $A$ with full column rank, then the first $n$ columns of $Q$ form a basis for the range of $A$ and the remaining $m-n$ columns form a basis for its orthogonal complement:
\begin{align*}
  \range(A) &= \range(Q_1) \quad \text{where}\quad Q_1 = \begin{bmatrix}q_1 & \cdots & q_n\end{bmatrix},
\\\range(A)^\perp=\Null(A) &= \range(Q_2) \quad \text{where}\quad Q_2 = \begin{bmatrix}q_{n+1} & \cdots & q_m\end{bmatrix},
\end{align*}
and for $k=1,\ldots,n$,
\begin{equation*}
   \span\{a_1,\ldots,a_k\} = \span\{q_1,\ldots,q_k\}
\end{equation*}
where the diagonal elements $r_{kk}\ne 0.$
}

This theorem tells us that a full column rank matrix can then be decomposed as
\begin{equation} \label{eq:qr-fact}
  A = Q R
    = \begin{bmatrix}Q_1 & Q_2\end{bmatrix}
      \begin{bmatrix}R_1 \\ 0 \end{bmatrix}
    = Q_1 R_1,
\end{equation}
where the triangular submatrix $R_1 = R[1{:}n,1{:}n]$ has no zeros on its diagonal. The factorization $A = Q_1 R_1$ is the _thin_ or _economy_ QR.


In Julia, we can compute the full QR decomposition of a matrix using via

```!
m, n = 4, 3
A = randn(m,n)
Q, R = qr(A)
```
and we can verify the orthogonality of `Q`:
```!
@show norm(Q'*Q - I)
@show norm(Q*Q' - I)
```

Julia [stores](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.QRCompactWY) the factor `Q` as a series of [Householder reflectors](https://en.wikipedia.org/wiki/Householder_transformation#QR_decomposition).

### Thin QR

To extract the "thin" QR factorization from `qr`, use the `Matrix` conversion function:

```!
F = qr(A)
Q₁ = Matrix(F.Q)
R₁ = F.R
@show size(Q₁)
@show size(R₁)
@show norm(A - Q₁*R₁)
```

## Solving least squares via QR

The QR factorization can be used to solve the least squares problem

\begin{equation} \label{eq:least-squares}
  \min_{x\in\R^n}\ \| A  x -b\|^2.
\end{equation}

Assume throughout this section that $A$ has full columns rank, and hence $m\geq n$. (The QR factorization also can be used to solve the more general problem, but we don't consider that case here.) Let $A = QR$ be the QR factorization shown in \eqref{eq:qr-fact}. Because the 2-norm is invariant under orthogonal rotation,

\begin{align} \label{eq:ls-derivation}
\|Ax -b\|^2
& = \| Q \T ( Ax -b)\|^2\\
&=\left\|\bmat{R_1\\0} x -\bmat{Q_1^T\\Q_2^T}b\right\|^2\\
& = \|R_1x  - Q_1^T b\|^2+\|Q_2^T b\|^2.
\end{align} 

Hence, minimizing $\|R_1 x - Q_1^T b\|$ also minimizes $\|Ax-b\|$. Because $ A $ is full rank,
$R_1$ is nonsingular, and the least-squares solution is obtained as the unique solution of the system 

\begin{equation} \label{eq:ls-via-qr}
  R_1 x  = Q_1^T b.
\end{equation}

The procedure for solving least-squares problems via the QR factorization can then be summarized as follows:

\algorithm{Least-squares via QR}{
\\ 
  - compute the thin QR factorization $A = Q_1 R_1$
  - backsolve the triangular linear system $R_1 x = Q_1^T b$

}

Here's a simple random example in Julia:

```!
m, n = 4, 3
A = randn(m, n)
b = randn(m)
F = qr(A)
Q₁, R₁ = Matrix(F.Q), F.R
x = R₁ \ Q₁'b   # <--- do NOT be tempted to use x = inv(R₁)*Q₁'b
; # hide
```
and we can verify that `x` is the least-squares solution by verifying that the residual $r=b-Ax$ is orthogonal to the columns of $A$: 
```!
r = b - A*x
@show norm(A'r)
```

The last line of the derivation shown in \eqref{eq:ls-derivation} asserts that the norm of the residual is equal to the norm of $Q_2^T b$. Let's check:

```!
Q₂ = F.Q[:,n+1:end]
norm(Q₂'b) ≈ norm(r)
```

## Accuracy of QR versus normal equations

A solution $x^*$ of a full-rank least-squares problem \eqref{eq:least-squares} can be obtained as the solution of the normal equations:
\begin{equation*}
  A\T A x = A\T b.
\end{equation*}
This approach, however, can in some cases suffer from numerical instability, which means that floating-point computations used by most computers yield solutions with unacceptable errors. The QR approach, on the other hand, often allows us to solve a wider-range of problems. The next example, borrowed from the course on [Fundamentals of Numerical Computation](https://fncbook.github.io/fnc/leastsq/demos/normaleqns-instab.html), illustrates this point.

The [Pythagorean identity](https://en.wikipedia.org/wiki/Pythagorean_trigonometric_identity) asserts that
\begin{equation*}
\sin^2(θ) + \cos^2(θ) = 1
\end{equation*}
for all values of $\theta$. This implies that the matrix
\begin{equation*}
  \bmat{
      \sin^2(\theta_1) & \cos^2(\theta_1) & 1
    \\ \vdots          & \vdots           & \vdots
    \\\sin^2(\theta_m) & \cos^2(\theta_m) & 1
},
\end{equation*}
for any values $\theta_1,\ldots,\theta_m$, doesn't have full column rank. Let's define a slightly perturbed version of this matrix that does have full column rank:
```!
θ = LinRange(0,3,400)
ε = 1e-7
A = @. [sin(θ)^2   cos(θ+ε)^2   θ^0]
; #hide
```
We can check that this matrix has full column rank:
```!
rank(A) == 3
```

Now create a right-hand side that corresponds to a known solution `xᵉ`:
```!
xᵉ = [1., 2., 1.]
b  = A*xᵉ 
; #hide
```

### Solution via the normal equations

Here's the solution obtained via the normal equations, and its corresponding relative error:
```!
xⁿ = A'A \ A'b
eⁿ = norm(xᵉ - xⁿ) / norm(xᵉ)
@show xⁿ
@show eⁿ
```
This solution has only about 1 digit of accuracy.

### Solution via QR factorization

Here's the solution obtained via the QR factorization, and its corresponding relative error:
```!
Q, R = qr(A); Q = Matrix(Q)
xʳ = R \ (Q'b)
eʳ = norm(xᵉ - xʳ) / norm(xᵉ)
@show xʳ
@show eʳ
```

Clearly, this solution is much more accurate than that obtained via the normal equations. The cost is essentially the same, so there's no real downside to using QR. In general, we can use QR implicitly via the backslash operator, which organizes its computations in a slightly different way, and obtains even better accuracy:
```!
xᵇ = A \ b
eᵇ = norm(xᵉ - xᵇ) / norm(xᵉ)
@show xᵇ
@show eᵇ
```


[^1]: See section 5.2 of Golub and Van Loan, _Matrix Computations_ (4th ed.), 2013.