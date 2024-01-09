# This file was generated, do not modify it. # hide
using Plots; gr()

x = y = -0.25:0.001:1.5
z = Surface((x,y)->f([x,y]), x, y)
contourf(x,y,z, levels=collect([10^k for k ∈ -4.0:3.0]), colorbar=:none, leg=false)
scatter!([1.],[1.], c="gold", ms=10, shape=:star5)
savefig(joinpath(@OUTPUT,"surface")) # hide