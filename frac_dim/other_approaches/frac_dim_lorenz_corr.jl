using PyPlot

function rhs(x,y,z)
    dx = σ*(y-x)
    dy = x*(ρ-z) - y
    dz = x*y - β*z
    return dx,dy,dz
end

function RK4(x,y,z)
    k1 = rhs(x,y,z)
    k2 = rhs(x+0.5*Δt*k1[1],y+0.5*Δt*k1[2],z+0.5*Δt*k1[3])
    k3 = rhs(x+0.5*Δt*k2[1],y+0.5*Δt*k2[2],z+0.5*Δt*k2[3])
    k4 = rhs(x+Δt*k3[1],y+Δt*k3[2],z+Δt*k3[3])
    dx,dy,dz = (k1 .+ 2*k2 .+ 2*k3 .+ k4)/6.
end



# time integration
function time_integration(x,y,z)
    # preallocate for storing results
    X = Array{Float64}(N+1)
    Y = Array{Float64}(N+1)
    Z = Array{Float64}(N+1)

    X[1] = Float64(x)
    Y[1] = Float64(y)
    Z[1] = Float64(z)


    for i = 1:N
        dx,dy,dz = RK4(x,y,z)
        x,y,z .+= Δt*(dx,y,dz)

        # store as 64bit
        X[i+1] = Float64(x)
        Y[i+1] = Float64(y)
        Z[i+1] = Float64(z)
    end

    return X,Y,Z
end

# create initial conditions
Nt = 15000
Δt = 0.01

σ = 10.
β = 8./3.
ρ = 28.

# start somewhere
xs,ys,zs = 0.6,0.5,15.

x,y,z = time_integration(Nt,xs,ys,zs,σ,ρ,β,Δt)

# # rescale use boundary to avoid a value that is exactly zero
# rescale(x,boundary=1e-10) = (x - (minimum(x)-boundary))/(maximum(x) - (minimum(x)-boundary))
#
# xc = rescale(xc)
# yc = rescale(yc)
# zc = rescale(zc)

# largest distance between consequtive points
Δx_max = sqrt(maximum((x[1:end-1]-x[2:end]).^2 + (y[1:end-1]-y[2:end]).^2 + (z[1:end-1]-z[2:end]).^2))
Δx_min = sqrt(minimum((x[1:end-1]-x[2:end]).^2 + (y[1:end-1]-y[2:end]).^2 + (z[1:end-1]-z[2:end]).^2))

## CORRELATION DIMENSION

lvec = (2.0.^(1:0.5:10))
lvec_sq = lvec.^2
Nl = zeros(Int,size(lvec))

for i = 1:Nt-1
    d = (x[i]-x[i+1:end]).^2 + (y[i]-y[i+1:end]).^2 + (z[i]-z[i+1:end]).^2
    for (il,l) in enumerate(lvec)
        Nl[il] += sum(d .< l)
    end
end

τ,μ = linreg(log10.(lvec),log10.(Nl/Nt^2))

## PLOTTING

fig,ax1 = subplots(1,1)

ax1[:plot](lvec,Nl/Nt^2,".-")
ax1[:plot]([Δx_max,Δx_max],[minimum(Nl),maximum(Nl)]./Nt^2,"k")
ax1[:plot]([Δx_min,Δx_min],[minimum(Nl),maximum(Nl)]./Nt^2,"k--")

ax1[:plot](lvec,10^τ*lvec.^μ,"C2--")
ax1[:plot](lvec,10^τ*lvec.^(1/2.),"C3--")

ax1[:set_xscale]("log",basex=2)
ax1[:set_yscale]("log",basex=2)

title("$(1/μ)")

xlabel("l")
ylabel(L"C(l) = n/N^2")
