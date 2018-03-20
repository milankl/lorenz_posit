using PyPlot

include("frac_dim_boxcount.jl")

function rhs(x,y,z)
    dxdt = σ*(y-x)
    dydt = x*(ρ-z) - y
    dzdt = x*y - β*z
    return dxdt,dydt,dzdt
end

function RK4(x,y,z)
    k1 = rhs(x,y,z)
    k2 = rhs(x+0.5*Δt*k1[1],y+0.5*Δt*k1[2],z+0.5*Δt*k1[3])
    k3 = rhs(x+0.5*Δt*k2[1],y+0.5*Δt*k2[2],z+0.5*Δt*k2[3])
    k4 = rhs(x+Δt*k3[1],y+Δt*k3[2],z+Δt*k3[3])
    dx,dy,dz = (k1 .+ 2.*k2 .+ 2.*k3 .+ k4)./6.
end

# time integration
function time_integration(Nt,x,y,z)
    # preallocate for storing results
    X = Array{Float64}(Nt+1)
    Y = Array{Float64}(Nt+1)
    Z = Array{Float64}(Nt+1)

    X[1] = Float64(x)
    Y[1] = Float64(y)
    Z[1] = Float64(z)

    for i = 1:Nt
        dx,dy,dz = RK4(x,y,z)
        x += Δt*dx
        y += Δt*dy
        z += Δt*dz

        # store as 64bit
        X[i+1] = Float64(x)
        Y[i+1] = Float64(y)
        Z[i+1] = Float64(z)
    end

    return X,Y,Z
end

## create initial conditions
Nt0 = 100000
Nt = 10000000
Δt = 0.0001

σ = 10.
β = 8./3.
ρ = 99.65

# start somewhere
x,y,z = 0.6,0.5,15.

xc,yc,zc = time_integration(Nt0,x,y,z)  # spin-up
xc,yc,zc = time_integration(Nt,xc[end],yc[end],zc[end])

##

# rescale use boundary to avoid a value that is exactly zero
rescale(x,boundary=1e-10) = (x - (minimum(x)-boundary))/(maximum(x) - (minimum(x)-boundary))

xrestart,yrestart,zrestart = xc[end],yc[end],zc[end]

xc = rescale(xc)
yc = rescale(yc)
zc = rescale(zc)

# largest distance between consequtive points
Δx = maximum(sqrt.((xc[1:end-1]-xc[2:end]).^2 + (yc[1:end-1]-yc[2:end]).^2 + (zc[1:end-1]-zc[2:end]).^2))

slope = box_count(xc,yc,zc)

print(slope)
