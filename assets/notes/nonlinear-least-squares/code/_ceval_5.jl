# This file was generated, do not modify it. # hide
function plotmap(b, x, x0)
scatter(xlim=(-1,+1), ylim=(-1,+1), leg=:outertopright,frame=:box, aspect_ratio=:equal)
scatter!(b[1,:], b[2,:], label="beacons", shape=:square, ms=7)
scatter!(x[1,:], x[2,:], label="true position", shape=:xcross, ms=7)
scatter!(x0[1,:], x0[2,:], label="initial guess", c="yellow", ms=7)
end
plotmap(b, x, x0)
savefig("map"); # hide
