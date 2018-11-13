using SigmoidNumbers
using JLD

cd("/home/kloewer/julia/lorenz_posit/frac_dim/")

include("frac_dim_boxcount.jl")
include("lorenz_integration_posits.jl")

frac_dims = Dict()
frac_dims["Nt"] = 10000000
frac_dims["P_nbits"] = [16,]#10:32
es = 1
frac_dims["P_ebits"] = ones(frac_dims["P_nbits"])*es
frac_dims["fracdim"] = Array{Float64}(length(frac_dims["P_nbits"]))

Nt = frac_dims["Nt"]

#for ip = 1:length(frac_dims["P_nbits"])
ip = 1

nbits = frac_dims["P_nbits"][ip]
ebits = frac_dims["P_ebits"][ip]

#Posit environment
P = Posit{nbits,ebits}

# create initial conditions
Nt0 = 1000      # length of spin-up
Δt = 0.001    # time step

σ = 10.      # Lorenz 63 coefficients
β = 8/3
ρ = 28.
s = 0.1

# start somewhere
xyz = [0.6,0.53,15.]

print("Posit($nbits,$ebits) with $Nt time steps.")

# do integration
xyzc = time_integration_opt(Nt0,Float64,xyz,σ,ρ,β,s,Δt)  # spin-up in full precision
xyzc = time_integration_opt(Nt,P,xyzc[:,end],σ,ρ,β,s,Δt)

# rescale - use boundary to avoid a value that is exactly zero
rescale(x,boundary=1e-10) = (x - (minimum(x)-boundary))/(maximum(x) - (minimum(x)-boundary))

xc = rescale(xyzc[1,:])
yc = rescale(xyzc[2,:])
zc = rescale(xyzc[3,:])

# BOX COUNTING ALGORITHM
frac_dims["fracdim"][ip] = box_count(xc,yc,zc)
println(" Frac Dim is $(frac_dims["fracdim"][ip])")

#end

# save
#save("data/frac_dims_posits_$(es)ebit.jld",frac_dims)
