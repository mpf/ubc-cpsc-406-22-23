# This file was generated, do not modify it. # hide
@df dfg scatter(:avgppb - p.(t))
plot!(leg=false, xlab="Month index from 1999 to 2021", ylab="ppb deviation from mean") # hide
savefig(joinpath(@OUTPUT,"ls-n2o-3")) # hide