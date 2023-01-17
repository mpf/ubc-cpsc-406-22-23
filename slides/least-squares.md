---
marp: true
headingDivider: 2
paginate: true
math: mathjax
footer: '[CPSC 406](https://friedlander.io/ubc-cpsc-406)'
---

# Linear Least Squares
<!-- _footer: '[CPSC 406 --- Computational Optimization](https://friedlander.io/ubc-cpsc-406), [University of British Columbia](https://ubc.ca)' -->

- formulate linear least-squares problems
- derive optimality conditions
- interpret as a projection
- solution via QR factorization

## Fitting a line

- $m$ pairs of data points $(t_1, b_1),(t_2,b_2),\ldots(t_m,b_m)$
- affine model of the data generation process:
  $$m_t(α,β):=α + tβ$$
- choose parameters $α$ (intercept) and $β$ (slope) that minimize least-squares errors ("residuals"):
  $$\min_{α,β}\ \sum_i^m (b_i - m_{t_i}(α,β))^2 = \|b-Ax\|^2$$
  with
  $$
        b = \begin{bmatrix}b_1\\\vdots\\b_m\end{bmatrix},
  \quad A = \begin{bmatrix}1 & t_1\\\vdots&\vdots\\1 & t_m\end{bmatrix},
  \quad x = \begin{bmatrix}\alpha \\ \beta\end{bmatrix} 
  $$

## Matrix formulation

- $m$ observations, $n$ model variables
- $A∈ℝ^{m×n}$ and $b∈ℝ^m$
- least-squares formulation
$$
\min_x\|Ax-b\|_2
\qquad\text{or}\qquad
\min_{x,r}\left\{\|r\|^2_2 \mid Ax +r = b\right\}
$$
- overdetermined ($m>n$)

## Orthogonal projection


The **orthogonal projection** of a vector $b$ onto a linear subspace $\mathcal{L}$:
$$
\mathbf{proj}_{\mathcal{L}}(b):=\underset{y}{\operatorname{argmin}} \{\|y-b\|_2 \mid y\in\mathcal{L}\}
$$

Let $(\bar x,\bar r)$ be a least-squares solution to
$$
\min_{x,r}\left\{\|r\|_2 \mid Ax +r = b\right\}
$$

- $\bar y:=A\bar x$ is the **orthogonal projection** of $b$ onto $\mathbf{range}(A)$
- $\bar r = b - A\bar x$ is the **orthogonal projection** of $b$ onto $\mathbf{range}(A)^\perp=\mathbf{null}(A^T)$

## Projection and normal equations



