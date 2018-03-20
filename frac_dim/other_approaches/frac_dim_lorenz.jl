using PyPlot

function rhs(x,y,z,σ,ρ,β)
    dx = σ*(y-x)
    dy = x*(ρ-z) - y
    dz = x*y - β*z
    return dx,dy,dz
end

# time integration
function time_integration(N,x,y,z,σ,ρ,β,Δt)
    # preallocate for storing results
    X = Array{Float64}(N+1)
    Y = Array{Float64}(N+1)
    Z = Array{Float64}(N+1)

    X[1] = Float64(x)
    Y[1] = Float64(y)
    Z[1] = Float64(z)

    for i = 1:N
        dx,dy,dz = rhs(x,y,z,σ,ρ,β)
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

# create initial conditions
Nt = 100000
Δt = 0.002

σ = 10.
β = 8./3.
ρ = 28.

# start somewhere
x,y,z = 0.5,0.5,15.

xc,yc,zc = time_integration(Nt,x,y,z,σ,ρ,β,Δt)

# rescale
#xc = (xc - minimum(xc))/(maximum(xc) - minimum(xc))
#yc = (yc - minimum(yc))/(maximum(yc) - minimum(yc))
#zc = (zc - minimum(zc))/(maximum(zc) - minimum(zc))

# largest distance between consequtive points
Δx = maximum(sqrt.((xc[1:end-1]-xc[2:end]).^2 + (yc[1:end-1]-yc[2:end]).^2 + (zc[1:end-1]-zc[2:end]).^2))

boxx = [minimum(xc),maximum(xc)]
boxy = [minimum(yc),maximum(yc)]
boxz = [minimum(zc),maximum(zc)]

L = maximum([boxx[2]-boxx[1],boxy[2]-boxy[1],boxz[2]-boxz[1]])
llc = [boxx[2]+boxx[1],boxy[2]+boxy[1],boxz[2]+boxz[1]]/2. - L/2.

##

o = 10
N = 10000     # number of scattered monte carlo points

slopes = Array{Float64}(o)
intercepts = Array{Float64}(o)

p = [0.05,0.5]

for oi=1:o
    ρ = Array{Float64}(N)   # shortest distance to a point on the attractor

    # monte carlo seeding
    for i = 1:N
        x,y,z = L*rand(3) + llc
        ρ[i] = sqrt(minimum((x-xc).^2 + (y-yc).^2 + (z-zc).^2))
    end

    sort!(ρ)

    # slope estimation
    p_vec = Int64(round(p[1]*N)):Int64(round(p[2]*N))
    ρ_at_p = ρ[p_vec]

    intercepts[oi],slopes[oi] = linreg(log10.(ρ_at_p),log10.(p_vec))
end

μ = mean(slopes)
μ_std = 2*std(slopes)
β = mean(intercepts)

println("$(3-μ) ± $μ_std")

#

fig,ax1 = subplots(1,1)

ax1[:loglog](ρ,1:N,lw=3,label=L"x^2")
ax1[:loglog](ρ[Int64.(round.(p*N))],p*N,"C1v")
ax1[:loglog](ρ,10^β*ρ.^μ,"C2--")

ax1[:loglog]([Δx,Δx],[1,N],"k")
