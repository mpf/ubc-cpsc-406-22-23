# This file was generated, do not modify it. # hide
@df dfg scatter(:avgppb) 
p(t) = α + β*t
plot!(p, lw=3)
plot!(leg=false, xlab="Month index from 1999 to 2021", ylab="Monthly average N₂O ppb")
savefig(joinpath(@OUTPUT,"ls-n2o-2")) # hide