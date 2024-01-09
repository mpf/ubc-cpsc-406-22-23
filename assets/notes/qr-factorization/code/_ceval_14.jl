# This file was generated, do not modify it. # hide
Q, R = qr(A); Q = Matrix(Q)
xʳ = R \ (Q'b)
eʳ = norm(xᵉ - xʳ) / norm(xᵉ)
@show xʳ
@show eʳ
