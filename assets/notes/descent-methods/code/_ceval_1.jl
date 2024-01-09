# This file was generated, do not modify it. # hide
using LinearAlgebra
import ForwardDiff:gradient,hessian

f(x) = (1 - x[1])^2 + 100.0(x[2] - x[1]^2)^2
∇f(x) = gradient(f, x)
∇²f(x) = hessian(f, x)
;#hide
