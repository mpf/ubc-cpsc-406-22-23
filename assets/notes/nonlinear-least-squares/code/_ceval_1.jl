# This file was generated, do not modify it. # hide
function gauss_newton(r, J, x; ε=1e-4, maxits=30)
  err = []
  for i = 1:maxits
      rk, Jk = r(x), J(x)
      push!(err, norm(Jk'rk))
      err[end] ≤ ε && break
      x = x - Jk\rk 
    end
    return x, err
end
; # hide
