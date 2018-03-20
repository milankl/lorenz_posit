using SigmoidNumbers
using JLD

cd("/home/kloewer/julia/lorenz_posit/predict/")
include("lorenz_integration_posits.jl")

# create initial conditions
Nlong = 10000
Δt = 0.01

σ,ρ,β = 10.,28.,8./3.

error = zeros(3)

# for s = 10
XYZ0 = time_integration(Nlong,Float64,[.5,.5,15.],σ,ρ,β,1/10.,Δt)
xi,xj = rand(1:10000,10000),rand(1:10000,10000)
error[1] = mean(sqrt.(sum((XYZ0[:,xi]-XYZ0[:,xj]).^2,1)[:]))

# for s = 1
XYZ0 = time_integration(Nlong,Float64,[.5,.5,15.],σ,ρ,β,1.,Δt)
xi,xj = rand(1:10000,10000),rand(1:10000,10000)
error[2] = mean(sqrt.(sum((XYZ0[:,xi]-XYZ0[:,xj]).^2,1)[:]))

# for s = 1/10
XYZ0 = time_integration(Nlong,Float64,[.5,.5,15.],σ,ρ,β,10.,Δt)
xi,xj = rand(1:10000,10000),rand(1:10000,10000)
error[3] = mean(sqrt.(sum((XYZ0[:,xi]-XYZ0[:,xj]).^2,1)[:]))

println(error)
save("data/mean_error.jld","meanerror",error)
