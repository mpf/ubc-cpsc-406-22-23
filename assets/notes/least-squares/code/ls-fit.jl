# This file was generated, do not modify it. # hide
using StatsPlots; pyplot()
x₁, x₄ = coef(ols)[[1,4]]
@df mtcars scatter(:WT, :MPG, xlabel="weight (1000s lbs)", ylabel="MPG", leg=false)
@df mtcars plot!(:WT, x₁ .+ x₄*:WT)
savefig(joinpath(@OUTPUT,"ls-fit")) # hide