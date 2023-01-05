# This file was generated, do not modify it. # hide
m, n = 4, 3
A = randn(m, n)
b = randn(m)
F = qr(A)
Q₁, R₁ = Matrix(F.Q), F.R
x = R₁ \ Q₁'b   # <--- do NOT be tempted to use x = inv(R₁)*Q₁'b
; # hide
