using JLD

cd("/home/kloewer/julia/lorenz_posit/dec_accuracy/")
include("lorenz_integration_floats.jl")

N = 10000000
Δt = 0.005

# Lorenz 63 coefficients
σ,ρ,β = 10.,28.,8./3.

# start somewhere - last argument is the scale
s = 1.
xyz = time_integration(N,[.5,.5,15.],σ,ρ,β,Δt,s)
save("data/lorenz_hr.jld","xyz",xyz)
