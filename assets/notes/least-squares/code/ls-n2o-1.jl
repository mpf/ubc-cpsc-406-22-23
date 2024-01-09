# This file was generated, do not modify it. # hide
using CSV, DataFrames, DataFramesMeta, Statistics
df = CSV.read(fname, DataFrame, comment="#", normalizenames=true, delim=" ", ignorerepeated=true)
dfg = @chain df begin
    @select(:year = :N2OcatsMLOyr, :month = :N2OcatsMLOmon, :ppb = :N2OcatsMLOm)
    @subset((!isnan).(:ppb))
    @by([:year, :month], :avgppb = mean(:ppb))
end 
describe(dfg, :min, :max, cols=[:year, :avgppb])