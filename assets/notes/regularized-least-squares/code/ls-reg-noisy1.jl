# This file was generated, do not modify it. # hide
n = 300
t = LinRange(0, 4, n)
x = @. sin(t) + t*cos(t)^2
w = 0.1*randn(n)
b = x + w
plot(b, leg=:topleft, label="noisy signal")
savefig(joinpath(@OUTPUT,"ls-reg-noisy1")) # hide