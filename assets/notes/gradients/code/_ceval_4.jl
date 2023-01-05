# This file was generated, do not modify it. # hide
a = collect(1:5)
ForwardDiff.gradient(x->a'x, rand(5))
; #hide
