function rhs(du,x,y,z,σ,ρ,β,s)
    du[1] = σ*(y-x)/s
    du[2] = (x*(ρ-z))/s - y
    du[3] = (x*y - β*z)/s
    return du
end

# time integration
function time_integration(Nt,P,xyz,σ,ρ,β,s,Δt)
    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    σ,ρ,β = s*σ,s*ρ,s*β
    Δt = s*Δt
    s_float = s

    # preallocate for storing results - store without scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)/s

    # Runge Kutta 4th order coefficients
    RKα = P.([s/6.,s/3.,s/3.,s/6.])
    RKβ = P.([s/2.,s/2.,s])

    # convert everything to the desired number system determined by P
    xyz = P.(xyz)
    σ,ρ,β,s,s²,Δt = P.([σ,ρ,β,s,s^2,Δt])

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    for i = 1:Nt
        xyz1 = deepcopy(xyz)

        for rki = 1:4
            dxyz = rhs(dxyz,xyz1[1],xyz1[2],xyz1[3],σ,ρ,β,s)

            if rki < 4
                xyz1[:] = xyz + (RKβ[rki]*Δt*dxyz)/s²
            end

            # sum the RK steps on the go
            xyz0 += (RKα[rki]*Δt*dxyz)/s²
        end

        xyz = deepcopy(xyz0)

        # store as 64bit undo scaling
        Xoutput[:,i+1] = Float64.(xyz)/s_float
    end

    return Xoutput
end

# SAME FOR INTS Integer(round(x)) to have desired rounding mode

function rhs_int(du,x,y,z,σ,ρ,β,s)
    du[1] = σ*(y-x)/s    # this will convert IntXX to Floats
    du[2] = x*(ρ-z)/s - y
    du[3] = (x*y - β*z)/s
    return du
end

# time integration
function time_integration_int(Nt,I,xyz,σ,ρ,β,s,Δt)
    # create conversion from float to Integer system with specified rounding mode
    Iround(x) = I(round(x))

    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    σ,ρ,β = s*σ,s*ρ,s*β
    Δt = s*Δt
    s_float = s
    s² = s^2

    # preallocate for storing results - store without scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)/s

    # Runge Kutta 4th order coefficients
    #RKα = Iround.([s/6.,s/3.,s/3.,s/6.])
    #RKβ = Iround.([s/2.,s/2.,s])
    RKα = [s/6.,s/3.,s/3.,s/6.]
    RKβ = [s/2.,s/2.,s]

    # convert everything to the desired number system determined by P
    xyz = Iround.(xyz)
    #σ,ρ,β,s,s²,Δt = Iround.([σ,ρ,β,s,s^2,Δt])

    # preallocate memory for intermediate results
    xyz0 = Float64.(deepcopy(xyz))
    xyz1 = Float64.(deepcopy(xyz))
    dxyz = Float64.(zeros(xyz))   # rhs

    for i = 1:Nt
        xyz1 = Float64.(deepcopy(xyz))

        for rki = 1:4
            dxyz = rhs_int(dxyz,Float64(xyz1[1]),Float64(xyz1[2]),Float64(xyz1[3]),σ,ρ,β,s)

            if rki < 4
                xyz1[:] = Float64.(xyz) + (RKβ[rki]*Δt*dxyz)/s²
            end

            # sum the RK steps on the go
            xyz0 += (RKα[rki]*Δt*dxyz)/s²
        end

        xyz = Iround.(deepcopy(xyz0))

        # store as 64bit undo scaling
        Xoutput[:,i+1] = Float64.(xyz)/s_float
    end

    return Xoutput
end
