# Column Orthogonalization

\blurb{The QR factorization of a matrix constructs an orthgonal basis for its columnspace. It's also one of the best computational tools for solving least-squares problems.}

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
```julia:qr
using LinearAlgebra
e₁ = [1, 0];  e₂ = [0, 1]
@show e₁'e₂
@show norm(e₁)
@show norm(e₂)
```
\show{qr}

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

```julia:./code/qr
using LinearAlgebra # gives `qr`
m, n = 4, 3
A = randn(m,n)
Q, R = qr(A)
```
\show{./code/qr}
and we can verify the orthogonality of `Q`:
```julia:./code/qr
@show norm(Q'*Q - I)
@show norm(Q*Q' - I)
```
\show{./code/qr}

Julia [stores](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.QRCompactWY) the factor `Q` as a series of [Householder reflectors](https://en.wikipedia.org/wiki/Householder_transformation#QR_decomposition).

### Thin QR

To extract the "thin" QR factorization from `qr`, use the `Matrix` conversion function:

```julia:./code/qr
F = qr(A)
Q₁ = Matrix(F.Q)
R₁ = F.R
@show size(Q₁)
@show size(R₁)
@show norm(A - Q₁*R₁)
```
\show{./code/qr}

## Solving least squares via QR

The QR factorization can be used to solve the least squares problem

\begin{equation*}
  \min_{x\in\R^n}\ \| A  x -b\|^2.
\end{equation*}

Assume throughout this section that $A$ has full columns rank, and hence $m\geq n$. (The QR factorization also can be used to solve the more general problem, but we don't consider that case here.) Let $A = QR$ be the QR factorization shown in \eqref{eq:qr-fact}. Because the 2-norm is invariant under orthogonal rotation,

\begin{align*}
\| A  x -b\|^2
& = \| Q \T ( Ax -b)\|^2\\
&=\left\|\bmat{R_1\\0} x -\bmat{Q_1^T\\Q_2^T}b\right\|^2\\
& = \|R_1x  - Q_1^T b\|^2+\|Q_2^T b\|^2.
\end{align*}

Hence, minimizing $\|R_1 x - Q_1^T b\|$ also minimizes $\|Ax-b\|$. Because $ A $ is full rank,
$R_1$ is nonsingular, and the least-squares solution is obtained as the unique solution of the system 

\begin{equation*}
  R_1 x  = Q_1^T b.
\end{equation*}

The procedure for solving least-squares problems via the QR factorization can then be summarized as follows:

\algorithm{Least-squares via QR}{
\\ 
  - compute the thin QR factorization $A = Q_1 R_1$
  - backsolve the triangular linear system $R_1 x = Q_1^T b$

}

Here's a simple random example in Julia:

```julia:./code/qr-ls
m, n = 4, 3
A = randn(m, n)
b = randn(m)
F = qr(A)
Q₁, R₁ = Matrix(F.Q), F.R
x = R₁ \ Q₁'b # do NOT use x = inv(R₁)*Q₁'b
```
and we can verify that `x` is the least-squares solution by verifying that the residual $r=b-Ax$ is orthogonal to the columns of $A$: 
```julia:./code/qr-ls
r = b - A*x
@show norm(A'r)
```
\show{./code/qr-ls}

```julia:./code/qr-ls
Q₂ = F.Q[:,n+1:end]
@show norm(Q₂'b)
@show norm(r)
```
\show{./code/qr-ls}

[^1]: See section 5.2 of Golub and Van Loan, _Matrix Computations_ (4th ed.), 2013.