using SigmoidNumbers
using JLD

cd("/home/kloewer/julia/lorenz_posit/frac_dim/")

include("frac_dim_boxcount.jl")
include("lorenz_integration_posits.jl")

frac_dims = Dict()
frac_dims["Nt"] = 1000000
frac_dims["P_nbits"] = 10:32
frac_dims["P_ebits"] = zeros(frac_dims["P_nbits"])
frac_dims["fracdim"] = Array{Float64}(length(frac_dims["P_nbits"]))

Nt = frac_dims["Nt"]

for ip = 1:length(frac_dims["P_nbits"])

    nbits = frac_dims["P_nbits"][ip]
    ebits = frac_dims["P_ebits"][ip]

    #Posit environment
    P = Posit{nbits,ebits}

    # create initial conditions
    Nt0 = 1000      # length of spin-up
    global Δt = P(0.01)    # time step

    global σ = P(10.)      # Lorenz 63 coefficients
    global β = P(8./3.)
    global ρ = P(28.)

    # start somewhere
    xyz = P.([0.6,0.5,15.])

    print("Posit($nbits,$ebits) with $Nt time steps.")

    # do integration
    xyzc = time_integration(Nt0,P,xyz)  # spin-up
    xyzc = time_integration(Nt,P,P.(xyzc[:,end]))

    # rescale - use boundary to avoid a value that is exactly zero
    rescale(x,boundary=1e-10) = (x - (minimum(x)-boundary))/(maximum(x) - (minimum(x)-boundary))

    xc = rescale(xyzc[1,:])
    yc = rescale(xyzc[2,:])
    zc = rescale(xyzc[3,:])

    # BOX COUNTING ALGORITHM
    frac_dims["fracdim"][ip] = box_count(xc,yc,zc)
    println(" Frac Dim is $(frac_dims["fracdim"][ip])")

end

# save
save("data/frac_dims_posits_0ebit.jld",frac_dims)
