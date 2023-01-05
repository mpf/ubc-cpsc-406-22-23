### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ c5054d49-06cb-4682-b070-71c3a397077f
md"""# Homework 2
!!! deadline
    Friday, 4 February, 11:59pm
"""

# ╔═╡ fe59b1e5-a4b2-45b1-91a8-c3e5fde92f5c
# Fill in your name and student ID
student = (name = "Harry Potter", ubcid = "123456789");

# ╔═╡ a70c0b8a-5d35-4c5c-bb79-262a83327717
md"""
**Student:** $(student.name) (ID# $(student.ubcid))
"""

# ╔═╡ f9527e44-3eed-42c4-8be5-47254955d498
md"""
## Q1. Quadratic functions (Beck, Exercise 2.19)
"""

# ╔═╡ 99df8118-7c73-11ec-201b-c75cb0f58794
md"""
## Q2. Homogeneity of the directional derivative

Prove that if the directional derivative of a function $f:\Re^n\to\Re$ is positively homogeneous, that is, show that for all $\alpha\ge0$,
```math
f'(x;\alpha d) = \alpha f'(x;d).
```

"""

# ╔═╡ f295ace3-50cc-471a-b0b1-f86f265993dd
md"""
## Q3. Beck, Exercise 4.3

Julia doesn't have the command `hilb` mentioned here in Beck's exercise. But it's not hard to write yourself.
"""

# ╔═╡ e1b34810-7546-4e92-9f1f-318fbc1d49a4
md"""
## Q4. Non-linear least squares (Beck, Exercise 4.6)

1. Compute the gradient $\nabla f(x)$. (Note that this is a slightly different objective than defined in class.)

1. Part (ii) of this problem.
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
# ╠═fe59b1e5-a4b2-45b1-91a8-c3e5fde92f5c
# ╟─f9527e44-3eed-42c4-8be5-47254955d498
# ╟─99df8118-7c73-11ec-201b-c75cb0f58794
# ╟─f295ace3-50cc-471a-b0b1-f86f265993dd
# ╟─e1b34810-7546-4e92-9f1f-318fbc1d49a4
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
