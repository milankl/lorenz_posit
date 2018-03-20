using JLD

cd("/home/kloewer/julia/lorenz_posit/frac_dim/")

include("frac_dim_boxcount.jl")
include("lorenz_integration_posits.jl")

frac_dims = Dict()
frac_dims["Nt"] = 1000000
frac_dims["sbits"] = 4:26
frac_dims["ebits"] = 11*ones(frac_dims["sbits"])
frac_dims["nbits"] = frac_dims["sbits"] + frac_dims["ebits"] + 1
frac_dims["fracdim"] = Array{Float64}(length(frac_dims["nbits"]))

for ip = 1:length(frac_dims["nbits"])

    nbits = frac_dims["nbits"][ip]
    ebits = frac_dims["ebits"][ip]

    # change global precision for bigfloats
    setprecision(frac_dims["sbits"][ip])

    # create initial conditions
    Nt0 = 1000      # length of spin-up
    global Δt = BigFloat("0.01")    # time step

    global σ = BigFloat("10.")      # Lorenz 63 coefficients
    global β = BigFloat("8")/BigFloat("3")
    global ρ = BigFloat("28.")

    # start somewhere
    xyz = BigFloat.([0.6,0.5,15.])

    for (iNt,Nt) in enumerate(frac_dims["Nt"])

        print("Float($nbits,$ebits) with $Nt time steps.")

        # do integration
        xyzc = time_integration(Nt0,BigFloat,xyz)  # spin-up
        xyzc = time_integration(Nt,BigFloat,BigFloat.(xyzc[:,end]))

        # rescale - use boundary to avoid a value that is exactly zero
        rescale(x,boundary=1e-10) = (x - (minimum(x)-boundary))/(maximum(x) - (minimum(x)-boundary))

        xc = rescale(xyzc[1,:])
        yc = rescale(xyzc[2,:])
        zc = rescale(xyzc[3,:])

        # BOX COUNTING ALGORITHM
        frac_dims["fracdim"][ip,iNt] = box_count(xc,yc,zc)
        println(" Frac Dim is $(frac_dims["fracdim"][ip,iNt])")
    end
end

# save
save("data/frac_dims_bigfloats.jld",frac_dims)
