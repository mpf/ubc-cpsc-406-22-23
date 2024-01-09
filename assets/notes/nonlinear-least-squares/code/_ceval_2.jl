# This file was generated, do not modify it. # hide
r(x, δ, b) = [ δ[i] - norm(x - b[:,i]) for i in 1:length(δ) ]
; #hide
