# This file was generated, do not modify it. # hide
D = finiteDiff(n)
b̂ = [b; zeros(n-1)]
plot(t, b, w=1, leg =:topleft, label="noisy data")
for λ in LogRange(1e0, 1e4, 3) 
    Â = [ I; √λ*D ]
    xLS = Â \ b̂
    plot!(t, xLS, w=2, label="regularized solution: λ = $(λ)")
end
savefig(joinpath(@OUTPUT,"ls-reg-noisy3")) # hide