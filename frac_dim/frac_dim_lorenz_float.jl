using JLD

cd("/home/kloewer/julia/lorenz_posit/frac_dim/")

include("frac_dim_boxcount.jl")
include("lorenz_integration_posits.jl")

frac_dims = Dict()
frac_dims["Nt"] = [1000000,2000000,3000000]
frac_dims["nbits"] = [16,32,64]
frac_dims["ebits"] = [5,8,11]
frac_dims["fracdim"] = Array{Float64}(length(frac_dims["nbits"]),length(frac_dims["Nt"]))

for ip = 1:length(frac_dims["nbits"])

    nbits = frac_dims["nbits"][ip]
    ebits = frac_dims["ebits"][ip]

    #Posit environment
    if ip == 1
        P = Float16
    elseif ip == 2
        P = Float32
    elseif ip == 3
        P = Float64
    end

    # create initial conditions
    Nt0 = 1000      # length of spin-up
    Δt = 0.01    # time step

    σ = 10.      # Lorenz 63 coefficients
    β = 8./3.
    ρ = 28.
    s = 0.1

    # start somewhere
    xyz = [0.6,0.5,15.]

    for (iNt,Nt) in enumerate(frac_dims["Nt"])

        print("Float($nbits,$ebits) with $Nt time steps.")

        # do integration
        xyzc = time_integration_opt(Nt0,Float64,xyz,σ,ρ,β,s,Δt)  # spin-up
        xyzc = time_integration_opt(Nt,P,xyzc[:,end],σ,ρ,β,s,Δt)

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
save("data/frac_dims_floats.jld",frac_dims)
