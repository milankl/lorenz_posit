function lorenz_hist(xyz,s,σ,ρ,β,Δt,bins)

    # for readability
    x = xyz[1,:]
    y = xyz[2,:]
    z = xyz[3,:]

    # preallocate
    H = zeros(length(bins)-1)

    # histogram of prognostic variables
    H = hist_add(H,x)
    H = hist_add(H,y)
    H = hist_add(H,z)

    #= Lorenz 63 rescaled is

        x_dot = σ(y-x)/s
        y_dot = x(ρ-z)/s - y
        z_dot = (xy - βz)/s

    and for the time stepping

        x_i+1 = x_i + (x*Δt)/s

    where also Δt is rescaled with s. Similar for RK4 steps, but divide by s^2 as
    the RK4 coefficients are also scaled by s. (But ignore here...)
    =#

    # rescale all coefficients
    σ = s*σ
    ρ = s*ρ
    β = s*β
    Δt = s*Δt

    # precompute RHS
    dx = σ*(y-x)/s
    dy = x.*(ρ-z)/s - y
    dz = (x.*y-β*z)/s

    # add to histogram all substeps of the RHS
    H = hist_add(H,y-x)             # x_dot
    H = hist_add(H,σ*(y-x))
    H = hist_add(H,dx)

    H = hist_add(H,ρ-z)             # y_dot
    H = hist_add(H,x.*(ρ-z))
    H = hist_add(H,x.*(ρ-z)/s)
    H = hist_add(H,dy)

    H = hist_add(H,x.*y)             # z_dot
    H = hist_add(H,β*z)
    H = hist_add(H,x.*y-β*z)
    H = hist_add(H,dz)

    # Euler time stepping procedure
    # ignore the addition to update the state vector as this is x,y,z
    H = hist_add(H,Δt*dx)
    H = hist_add(H,Δt*dy)
    H = hist_add(H,Δt*dz)

    H = hist_add(H,Δt*dx/s)
    H = hist_add(H,Δt*dy/s)
    H = hist_add(H,Δt*dz/s)

    return H
end

function lorenz_hist_opt(xyz,s,σ,ρ,β,Δt,bins)

    # for readability
    x = xyz[1,:]
    y = xyz[2,:]
    z = xyz[3,:]

    # preallocate
    H = zeros(length(bins)-1)

    # histogram of prognostic variables
    H = hist_add(H,x)
    H = hist_add(H,y)
    H = hist_add(H,z)

    #= Lorenz 63 rescaled is

        x_dot = σ(y-x)/s
        y_dot = x(ρ-z)/s - y
        z_dot = (xy - βz)/s

    and for the time stepping

        x_i+1 = x_i + (x*Δt)/s

    where also Δt is rescaled with s. Similar for RK4 steps, but divide by s^2 as
    the RK4 coefficients are also scaled by s. (But ignore here...)
    =#

    # precompute RHS
    dx = y-x
    dy = x.*(ρ-z/s) - y
    dz = x.*(y/s) - β*z

    # add to histogram all substeps of the RHS
    H = hist_add(H,dx)              # x_dot

    H = hist_add(H,z/s)             # y_dot
    H = hist_add(H,ρ-z/s)
    H = hist_add(H,x.*(ρ-z/s))
    H = hist_add(H,dy)

    H = hist_add(H,y/s)             # z_dot
    H = hist_add(H,dz)

    # Euler time stepping procedure
    # ignore the addition to update the state vector as this is x,y,z
    H = hist_add(H,σ*Δt*dx)
    H = hist_add(H,Δt*dy)
    H = hist_add(H,β*Δt*dz)

    return H
end

function hist_add(H0,vec)
    Htemp,bin_edges = np.histogram(abs.(vec),bins)
    H0 += Htemp
    return H0
end
