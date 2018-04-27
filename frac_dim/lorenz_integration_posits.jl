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

function rhs_opt(du,x,y,z,ρ,β,s)
    du[1] = y-x
    du[2] = x*(ρ-z/s) - y
    du[3] = x*(y/s) - β*z
    return du
end

# time integration
function time_integration_opt(Nt,P,xyz,σ,ρ,β,s,Δt)
    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    s_float = s

    # preallocate for storing results - store without scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)/s

    # Runge Kutta 4th order coefficients including time step and sigma for x
    RKα = zeros(3,4)
    RKα[2,:] = 1./([1/6.,1/3.,1/3.,1/6.]*Δt)
    RKα[3,:] = RKα[2,:]
    RKα[1,:] = RKα[2,:]/σ
    RKαz = RKα[2,:]/β

    RKβ = zeros(3,3)
    RKβ[2,:] = 1./([1/2.,1/2.,1.]*Δt)
    RKβ[3,:] = RKβ[2,:]
    RKβ[1,:] = RKβ[2,:]/σ
    RKβz = RKβ[2,:]/β

    # convert everything to the desired number system determined by P
    xyz = P.(xyz)
    ρ,β,s = P.([ρ,β,s])
    RKα = P.(RKα)
    RKβ = P.(RKβ)

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    for i = 1:Nt
        xyz1 = deepcopy(xyz)

        for rki = 1:4
            dxyz = rhs_opt(dxyz,xyz1[1],xyz1[2],xyz1[3],ρ,β,s)

            if rki < 4
                xyz1[:] = xyz + dxyz./RKβ[:,rki]
            end

            # sum the RK steps on the go
            xyz0 += dxyz./RKα[:,rki]
        end

        xyz = deepcopy(xyz0)

        # store as 64bit undo scaling
        Xoutput[:,i+1] = Float64.(xyz)/s_float
    end

    return Xoutput
end

# SAME FOR INTS Integer(round(x)) to have desired rounding mode

function rhs_int(du,Iround,x,y,z,ρ,β1,β2,s)
    du[1] = y-x
    du[2] = x*(ρ-Iround(z/s)) - y
    du[3] = x*(Iround(y/s)) - Iround(z/β2)*β1
    return du
end

# time integration
function time_integration_int(Nt,I,xyz,σ,ρ,β1,β2,s,Δt)
    # create conversion from float to Integer system with specified rounding mode
    Iround(x) = I(round(x))

    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    s_float = s

    # preallocate for storing results - store without scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)/s

    # Runge Kutta 4th order coefficients including time step and sigma for x
    RKα = zeros(3,4)
    RKα[2,:] = 1./([1/6.,1/3.,1/3.,1/6.]*Δt)
    RKα[3,:] = RKα[2,:]
    RKα[1,:] = RKα[2,:]/σ
    RKαz = RKα[2,:]/(β1/β2)

    RKβ = zeros(3,3)
    RKβ[2,:] = 1./([1/2.,1/2.,1.]*Δt)
    RKβ[3,:] = RKβ[2,:]
    RKβ[1,:] = RKβ[2,:]/σ
    RKβz = RKβ[2,:]/(β1/β2)

    # convert everything to the desired number system determined by P
    xyz = Iround.(xyz)
    ρ,β1,β2,s = Iround.([ρ,β1,β2,s])
    RKα = Iround.(RKα)
    RKβ = Iround.(RKβ)
    RKαz = Iround.(RKαz)
    RKβz = Iround.(RKβz)

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    zeroo = Iround(0)

    for i = 1:Nt
        xyz1 = deepcopy(xyz)

        for rki = 1:4
            #z = xyz1[3]
            dxyz = rhs_int(dxyz,Iround,xyz1[1],xyz1[2],xyz1[3],ρ,β1,β2,s)

            # sum the RK steps on the go
            xyz0 += Iround.(dxyz./RKα[:,rki])# - Iround.([zeroo,zeroo,z]./RKαz[rki])

            if rki < 4
                xyz1[:] = xyz + Iround.(dxyz./RKβ[:,rki])# - Iround.([zeroo,zeroo,z]./RKβz[rki])
            end
        end

        xyz = deepcopy(xyz0)

        # store as 64bit undo scaling
        Xoutput[:,i+1] = Float64.(xyz)/s
    end

    return Xoutput
end

# time integration
function time_integration_fullrhs(Nt,I,xyz,σ,ρ,β,s,Δt)
    # create conversion from float to Integer system with specified rounding mode
    #Iround(x) = I(round(x))

    # apply the scale to the initial conditions and L63 parameters
    xyz = s*xyz
    s_float = s

    # preallocate for storing results - store without scaling
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)/s

    # Runge Kutta 4th order coefficients including time step and sigma for x
    RKα = zeros(3,4)
    RKα[2,:] = 1./([1/6.,1/3.,1/3.,1/6.]*Δt)
    RKα[3,:] = RKα[2,:]
    RKα[1,:] = RKα[2,:]/σ
    RKαz = RKα[2,:]/β

    RKβ = zeros(3,3)
    RKβ[2,:] = 1./([1/2.,1/2.,1.]*Δt)
    RKβ[3,:] = RKβ[2,:]
    RKβ[1,:] = RKβ[2,:]/σ
    RKβz = RKβ[2,:]/β

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    for i = 1:Nt
        xyz1 = deepcopy(xyz)

        for rki = 1:4
            dxyz = rhs_opt(dxyz,xyz1[1],xyz1[2],xyz1[3],ρ,β,s)

            # sum the RK steps on the go
            xyz0 += dxyz./RKα[:,rki]

            if rki < 4
                xyz1[:] = xyz + dxyz./RKβ[:,rki]
            end
        end

        xyz = Float64.(I.(deepcopy(xyz0)))

        # store as 64bit undo scaling
        Xoutput[:,i+1] = Float64.(xyz)/s
    end

    return Xoutput
end
