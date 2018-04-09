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

# SAME FOR INTS REPLACE "/" by div operator to have floor rounding

function rhs_int(du,x,y,z,σ,ρ,β,s)
    du[1] = div(σ*(y-x),s)
    du[2] = div(x*(ρ-z),s) - y
    du[3] = div(x*y - β*z,s)
    return du
end

# time integration
function time_integration_int(Nt,P,xyz,σ,ρ,β,s,Δt)
    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    σ,ρ,β = s*σ,s*ρ,s*β
    Δt = s*Δt
    s_float = s

    # preallocate for storing results - store without scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)/s

    # Runge Kutta 4th order coefficients
    RKα = P.(round.([s/6.,s/3.,s/3.,s/6.]))
    RKβ = P.(round.([s/2.,s/2.,s]))

    # convert everything to the desired number system determined by P
    xyz = P.(round.(xyz))
    σ,ρ,β,s,s²,Δt = P.(round.([σ,ρ,β,s,s^2,Δt]))

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    for i = 1:Nt
        xyz1 = deepcopy(xyz)

        for rki = 1:4
            dxyz = rhs_int(dxyz,xyz1[1],xyz1[2],xyz1[3],σ,ρ,β,s)

            if rki < 4
                xyz1[:] = xyz + div.(RKβ[rki]*Δt*dxyz,s²)
            end

            # sum the RK steps on the go
            xyz0 += div.(RKα[rki]*Δt*dxyz,s²)
        end

        xyz = deepcopy(xyz0)

        # store as 64bit undo scaling
        Xoutput[:,i+1] = Float64.(xyz)/s_float
    end

    return Xoutput
end
