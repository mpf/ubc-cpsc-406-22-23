# This file was generated, do not modify it. # hide
m = 13                 # number of beacons
b = 2rand(2, m) .- 1   # place beacons uniformly in the unit box 
x = zeros(2)           # true position 
x0 = 2rand(2) .- 1     # initial guess of the unknown position
ν = .5*rand(m)
δ = [ norm(x - b[:,i]) + ν[i] for i in 1:m]
; # hide
