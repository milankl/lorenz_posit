# FRACTAL DIMENSION
# regard the unit cube [0,1]
using PyPlot

# assume the attractor is a point
x0,y0,z0 = [0.5,0.5,0.5]

N = 1000000     # number of scattered monte carlo points
#ρ = Array{Float64}(N)   # shortest distance to a point on the attractor

ρ = sort(sqrt.((rand(N)-x0).^2 + (rand(N)-y0).^2 + (rand(N)-z0).^2))

loglog(ρ,1:N)
s = [1e-2,1]
loglog(s,1e6*s.^1,label="x^1")
loglog(s,1e6*s.^2,label="x^2")
loglog(s,1e6*s.^3,label="x^3")
legend()
