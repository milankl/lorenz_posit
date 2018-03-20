function rhs(du,x,y,z,σ,ρ,β,s)
    du[1] = σ*(y-x)/s
    du[2] = (x*(ρ-z))/s - y
    du[3] = (x*y - β*z)/s
    return du
end

# time integration
function time_integration(Nt,xyz,σ,ρ,β,Δt,s)
    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    σ,ρ,β = s*σ,s*ρ,s*β
    Δt = s*Δt

    # preallocate for storing results - store with scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = xyz

    # Runge Kutta 4th order coefficients
    RKα = [s/6.,s/3.,s/3.,s/6.]
    RKβ = [s/2.,s/2.,s]

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    for i = 1:Nt
        xyz1 = deepcopy(xyz)

        for rki = 1:4
            dxyz = rhs(dxyz,xyz1[1],xyz1[2],xyz1[3],σ,ρ,β,s)

            if rki < 4
                xyz1[:] = xyz + (RKβ[rki]*Δt*dxyz)/s^2
            end

            # sum the RK steps on the go
            xyz0 += (RKα[rki]*Δt*dxyz)/s^2
        end

        xyz = deepcopy(xyz0)

        # store as 64bit - remove the scaling
        Xoutput[:,i+1] = xyz
    end

    return Xoutput
end
