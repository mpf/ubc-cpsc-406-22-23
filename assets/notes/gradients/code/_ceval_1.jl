# This file was generated, do not modify it. # hide
a, b = 1, 100
f(x) = (a - x[1])^2 + b*(x[2] - x[1]^2)^2
∇f(x) = [-2(a - x[1]) - 4b*(x[2] - x[1]^2)*x[1] , 2b*(x[2] - x[1]^2) ]
; #hide
