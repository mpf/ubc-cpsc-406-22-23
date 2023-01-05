### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ c5054d49-06cb-4682-b070-71c3a397077f
md"""# Homework 3
!!! deadline
    Monday, 14 February, 11:59pm
"""

# ╔═╡ e0626ba4-87ae-45df-82c7-c7e0c767843d
md"Fill in your name and student ID:"

# ╔═╡ fe59b1e5-a4b2-45b1-91a8-c3e5fde92f5c
student = (name = "Harry Potter", ubcid = "123456789");

# ╔═╡ a70c0b8a-5d35-4c5c-bb79-262a83327717
md"""
**Student:** $(student.name) (ID# $(student.ubcid))
"""

# ╔═╡ f9527e44-3eed-42c4-8be5-47254955d498
md"""
## Q1. Implement Newton's method and variants (Beck, Exercise 5.2)
"""

# ╔═╡ 99df8118-7c73-11ec-201b-c75cb0f58794
md"""
## Q2. Cholesky factorization
Partition the $n$-by-$n$ positive definite matrix $A$ as
```math
A = \begin{bmatrix} a_{11} & w^T \\ w & K \end{bmatrix},
```
where, respectively, $a_{11}$ and $w$ are the top-left entry and the remaining first-column of $A$. The first iteration of the Cholesky process decomposes $A$ as
```math
\begin{align}
  A = \begin{bmatrix} \alpha & 0 \\ w/\alpha & I \end{bmatrix}
      \begin{bmatrix} 1 & 0 \\ 0 & K-ww^T/a_{11} \end{bmatrix}
      \begin{bmatrix} \alpha & w^T/\alpha \\ 0 & I \end{bmatrix} = R_1^TA_1R_1,
\end{align}
```
where $\alpha:=\sqrt{a_{11}}$, and $R_1$ and $A_1$ are defined by the matrices in the factorization.
1. This first step is valid because $a_{11}$ is positive. Why?
The next iterations of the factorization proceed inductively on the submatrix $K-ww^T/a_{11}$. But again, that's only possible if its first entry is positive.
2. Prove that this must be the case.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─a70c0b8a-5d35-4c5c-bb79-262a83327717
# ╟─c5054d49-06cb-4682-b070-71c3a397077f
# ╟─e0626ba4-87ae-45df-82c7-c7e0c767843d
# ╠═fe59b1e5-a4b2-45b1-91a8-c3e5fde92f5c
# ╟─f9527e44-3eed-42c4-8be5-47254955d498
# ╟─99df8118-7c73-11ec-201b-c75cb0f58794
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
