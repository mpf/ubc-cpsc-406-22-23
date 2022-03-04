### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# ╔═╡ 01689882-1db2-4221-94ba-0b4a5d16e421
# Import Packages
begin
	import ForwardDiff
	using LinearAlgebra, Random
	using Printf
	Random.seed!(123)
end

# ╔═╡ 4afa77d1-ada2-47df-a7e6-542f5488a317
md"""
## Gradient method with backtracking
"""

# ╔═╡ 005dbe4d-77ff-48b9-9f0a-1aba91233a04
md"""
## Hybrid Gradient-Newton method
"""

# ╔═╡ 514cd782-71cf-4d81-b0d5-a80986f04ba6
md"""
## Damped Gauss-Newton
"""

# ╔═╡ 72937ebb-e880-424e-8017-36b59dee0b32
begin
	#s = 1
	#α = 0.5
	#β = 0.5
	ϵ = 1e-5
	maxits = 1e4
end

# ╔═╡ cc09d34b-28da-4fde-b918-3f551dd591ed
md"""
### Define $f$ and gradients
"""

# ╔═╡ 87c6db7d-955c-4a5c-98fd-b391dd31fe0f
begin
	f1(x) = -13+x[1]+((5-x[2])*x[2]-2)*x[2]
	f2(x) = -29+x[1]+((x[2]+1)*x[2]-14)*x[2]
	f(x) = f1(x)^2 + f2(x)^2
	∇f(x) = ForwardDiff.gradient(f,x)
	∇2f(x) = ForwardDiff.hessian(f,x)
	∇f1(x) = ForwardDiff.gradient(f1,x)
	∇f2(x) = ForwardDiff.gradient(f2,x)
	F(x) = [ f1(x) , f2(x) ]
	J(x) = ForwardDiff.jacobian(F,x)
end

# ╔═╡ a5972d4d-2f90-4ad9-83ca-fe6a7529bf1c
md"""
# Run algorithms
"""

# ╔═╡ 9c89fa69-c315-4fbb-a630-e07a4d214048
md"""
## Results
"""

# ╔═╡ b0cdbb43-afbd-4d32-9308-14157f34be50
md"""
### Gradient Method with backtracking
"""

# ╔═╡ ab5debc6-9abd-4ea1-a4b2-35e46435438d
md"""
Gradient Method converges to the global minimizer in roughly 2000 iterations for the starting points $(-50,7)$, $(20,7)$, and $(5, -10)$. For the starting point $(20,-18)$, Gradient Method converges to a local minimizer.
"""

# ╔═╡ 718de041-eeb2-42e3-827b-39d9b82edc46
md"""
### Hybrid Gradient-Newton method
"""

# ╔═╡ 536310cb-4bd0-4803-9013-03e563b4bf42
md"""
Hybrid Gradient-Newton converges to a global min in 9 iterations for the first two starting points, and then converges to a local min in <20 iterations for the other two starting points.
"""

# ╔═╡ 5179d962-3245-4d21-9bad-99c80a1b391e
md"""
### Damped Gauss-Newton
"""

# ╔═╡ 6a0ae461-5d3e-4939-869f-8d34d7fcea0c
md"""
Damped Gauss-Newton converges in 7 iterations to the global minimizer for the first two starting points, but doesn't converge for the other two starting points
"""

# ╔═╡ 00233933-c517-4b56-937f-2ac43917c75d
md""" 
## Helper Utilities
"""

# ╔═╡ 766c087e-b4e8-419e-813f-6479375e5b5e
begin
	struct State
		x::Vector
		f
		norm∇f
	end
	State(x::Vector, f::Real, ∇f::Vector) = State(x, f, norm(∇f))
end

# ╔═╡ a3edae1e-8e98-11ec-21b6-cd68c4acd72b
# Gradient method with backtracking
function grad_method_backtracking(f, g, x0; ϵ = 1e-5, s=1.0, α = 0.5, β=0.5,maxits = 1e4)
    x = copy(x0)
    fk = f(x); ∇fk = g(x)
    k = 0
    trace = [State(x, fk, ∇fk)]
    while norm(∇fk) > ϵ && k < maxits
        #α = 1.0
		tk = s
        while ((fk - f(x-tk*∇fk)) < -α*tk*dot(∇fk,-∇fk) )
            #α /= 2
			tk = tk*β
        end
        x = x - tk*∇fk
        fk = f(x); ∇fk = g(x)
        k += 1
		push!(trace, State(x, fk, ∇fk))
    end
    return x, trace
end

# ╔═╡ 24d814d0-cc8b-4fa8-901b-c5b7b1ec9343
function hybrid_gradient_newton(f, g, h, x0; ϵ = 1e-5, s=1.0, α = 0.5, β=0.5,maxits = 1e4)
    x = copy(x0)
    fk = f(x); ∇fk = g(x); Hk = h(x)
    k = 0
    trace = [State(x, fk, ∇fk)]
    while norm(∇fk) > ϵ && k < maxits
		if isposdef(Hk)
			#dk = (-∇fk \ Hk)'
			dk = - Hk \ ∇fk
		else
			dk = -∇fk
		end
		tk = s
        while ((fk - f(x+tk*dk)) < -α*tk*dot(∇fk,dk))
			tk = tk*β
        end
        x = x + tk*dk
        fk = f(x); ∇fk = g(x); Hk = h(x)
        k += 1
		push!(trace, State(x, fk, ∇fk))
    end
    return x, trace
end

# ╔═╡ 5c864ccd-b289-43aa-bc4e-b04033c61d57
function damped_gauss_newton(F,J,f,g,x0; ϵ = 1e-5, s=1.0, α = 0.5, β=0.5,maxits = 1e4)
    x = copy(x0)
	Fk = F(x); Jk = J(x);
	fk = f(x); ∇fk = g(x); 
    k = 0
    trace = [State(x, fk, ∇fk)]
    while norm(∇fk) > ϵ && k < maxits
		dk = Jk \ Fk
		tk = s
        while ((fk - f(x-tk*dk)) < -α*tk*dot(∇fk,dk))
			tk = tk*β
        end
        x = x - tk*dk
		Fk = F(x); Jk = J(x);
        fk = f(x); ∇fk = g(x);
        k += 1
		push!(trace, State(x, fk, ∇fk))
    end
    return x, trace
end

# ╔═╡ bddf3d4d-bd8e-44a1-b7ca-1f386853dcc0
begin
	struct Results
		p::Vector ## initial point
		x::Vector ## final point
		num_iter ## number of iterations ran
	end
	#Results(p,algo,f,∇f) = Results(p, algo(f,∇f,p)[1], length(algo(f,∇f,p)[2]))
	Results(trace::Vector{State}) = Results(first(trace).x, last(trace).x, length(trace))
end

# ╔═╡ 738fc79c-6f44-4c64-a743-efa57879ea30
begin
	x0 = [-50.0, 7.0]
	x1 = [20.0, 7.0]
	x2 = [20.0, -18.0]
	x3 = [5.0, -10.0]

	points = [x0, x1, x2, x3]
	
	gm_results = []
	hgn_results = []
	dgn_results = []
	
	for p in points
		push!(gm_results, Results(grad_method_backtracking(f,∇f,p)[2]))
		push!(hgn_results, Results(hybrid_gradient_newton(f, ∇f, ∇2f, p)[2] ))
		push!(dgn_results, Results(damped_gauss_newton(F,J,f,∇f,p)[2]))
	end
end

# ╔═╡ c0d1d15a-b341-45a5-ae6d-694d11162b4b
function log(states::Vector{State})
	io = IOBuffer()
	write(io, "| k | f | ∇f | \n")
	write(io, "|---|---|---| \n")
	for (k, s) in enumerate(states)
		f = @sprintf("%10.2e", s.f)
		nrmf = @sprintf("%10.2e", s.norm∇f)
		write(io, "| $(k-1) | $f | $nrmf |\n")
	end
	Markdown.parse(seekstart(io))
end

# ╔═╡ 86b5ba5f-a315-49cd-80a9-6ee93b51ae16
#function results_table(results::Vector{Results})
function results_table(results)
	io = IOBuffer()
	write(io, "|initial point | final point | number of iterations | \n")
	write(io, "|:---:|:---:|:---:| \n")
	for (i,r) in enumerate(results)
		p = r.p
		#[@sprintf("%10.2e",i) for i in r.p]
		#@sprintf("%10.2e", s.f)
		x = [ @sprintf("%.3f",i) for i in r.x ]
		num_iter = r.num_iter
		write(io, "| $p | $x | $num_iter | \n")
	end
	Markdown.parse(seekstart(io))
end

# ╔═╡ a4cfc6a7-777b-4a04-b7dd-d0a750c65cca
results_table(gm_results)

# ╔═╡ 1c6a5dc3-981b-48d8-b538-4d4af382e18b
results_table(hgn_results)

# ╔═╡ 6cef59b6-f04b-4f75-9f84-b4da2b9ecb89
results_table(dgn_results)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
ForwardDiff = "~0.10.25"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "f9982ef575e19b0e5c7a98c6e75ee496c0f73a93"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "84083a5136b6abf426174a58325ffd159dd6d94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "2cf929d64681236a2e074ffafb8d568733d2e6af"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "8d0c8e3d0ff211d9ff4a0c2307d876c99d10bdf1"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "95c6a5d0e8c69555842fc4a927fc485040ccc31c"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.5"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─01689882-1db2-4221-94ba-0b4a5d16e421
# ╟─4afa77d1-ada2-47df-a7e6-542f5488a317
# ╠═a3edae1e-8e98-11ec-21b6-cd68c4acd72b
# ╟─005dbe4d-77ff-48b9-9f0a-1aba91233a04
# ╠═24d814d0-cc8b-4fa8-901b-c5b7b1ec9343
# ╟─514cd782-71cf-4d81-b0d5-a80986f04ba6
# ╠═5c864ccd-b289-43aa-bc4e-b04033c61d57
# ╟─72937ebb-e880-424e-8017-36b59dee0b32
# ╟─cc09d34b-28da-4fde-b918-3f551dd591ed
# ╠═87c6db7d-955c-4a5c-98fd-b391dd31fe0f
# ╠═a5972d4d-2f90-4ad9-83ca-fe6a7529bf1c
# ╠═738fc79c-6f44-4c64-a743-efa57879ea30
# ╟─9c89fa69-c315-4fbb-a630-e07a4d214048
# ╟─b0cdbb43-afbd-4d32-9308-14157f34be50
# ╟─a4cfc6a7-777b-4a04-b7dd-d0a750c65cca
# ╟─ab5debc6-9abd-4ea1-a4b2-35e46435438d
# ╟─718de041-eeb2-42e3-827b-39d9b82edc46
# ╟─1c6a5dc3-981b-48d8-b538-4d4af382e18b
# ╟─536310cb-4bd0-4803-9013-03e563b4bf42
# ╟─5179d962-3245-4d21-9bad-99c80a1b391e
# ╠═6cef59b6-f04b-4f75-9f84-b4da2b9ecb89
# ╟─6a0ae461-5d3e-4939-869f-8d34d7fcea0c
# ╟─00233933-c517-4b56-937f-2ac43917c75d
# ╠═bddf3d4d-bd8e-44a1-b7ca-1f386853dcc0
# ╠═766c087e-b4e8-419e-813f-6479375e5b5e
# ╠═c0d1d15a-b341-45a5-ae6d-694d11162b4b
# ╠═86b5ba5f-a315-49cd-80a9-6ee93b51ae16
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
