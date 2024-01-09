# This file was generated, do not modify it. # hide
J(x) = ForwardDiff.jacobian(x->r(x, δ, b), x)
xs, err = gauss_newton(x->r(x, δ, b), J, x0)
plot(err, yaxis=:log)
savefig("convergence"); # hide
