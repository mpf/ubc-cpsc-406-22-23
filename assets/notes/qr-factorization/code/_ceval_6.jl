# This file was generated, do not modify it. # hide
F = qr(A)
Q₁ = Matrix(F.Q)
R₁ = F.R
@show size(Q₁)
@show size(R₁)
@show norm(A - Q₁*R₁)
