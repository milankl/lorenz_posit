function rhs(du,x,y,z,σ,ρ,β,s)
    du[1] = σ*(y-x)
    du[2] = x*(ρ-s*z) - y
    du[3] = s*x*y - β*z
    return du
end

# time integration
function time_integration(Nt,P,xyz,σ,ρ,β,s,Δt)
    # preallocate for storing results
    Xoutput = Array{Float64}(3,Nt+1)
    Xoutput[:,1] = Float64.(xyz)

    # Runge Kutta 4th order coefficients
    RKα = P.([1/6.,1/3.,1/3.,1/6.])
    RKβ = P.([0.5,0.5,1.])

    # convert everything to the desired number system determined by P
    xyz = P.(xyz)
    σ,ρ,β,s,Δt = P.([σ,ρ,β,s,Δt])

    # preallocate memory for intermediate results
    xyz0 = deepcopy(xyz)
    xyz1 = deepcopy(xyz)
    dxyz = zeros(xyz)   # rhs

    for i = 1:Nt
        xyz1 = copy(xyz)    # deep copy apparently unnecessary for vectors

        for rki = 1:4
            dxyz = rhs(dxyz,xyz1[1],xyz1[2],xyz1[3],σ,ρ,β,s)

            if rki < 4
                xyz1[:] = xyz + RKβ[rki]*Δt*dxyz
            end

            # sum the RK steps on the go
            xyz0 += RKα[rki]*Δt*dxyz
        end

        xyz = copy(xyz0)

        # store as 64bit
        Xoutput[:,i+1] = Float64.(xyz)
    end

    return Xoutput
end
