# FRACTAL DIMENSION
# regard the unit cube [0,1]
using PyPlot

# assume the attractor is a circle of radius r in the x,y plane represented by m points
m = 100
r = 0.5
Δx = 2π*r/m     # distance of points on the attractor
c = r*e.^(2*π*1im*(1:m)/m)
xc,yc,zc = real(c),imag(c),0.5*ones(m)

o = 10
N = 100000     # number of scattered monte carlo points

slopes = Array{Float64}(o)
intercepts = Array{Float64}(o)

for oi=1:o
    ρ = Array{Float64}(N)   # shortest distance to a point on the attractor

    # monte carlo seeding
    for i = 1:N
        x,y,z = rand(3)
        ρ[i] = sqrt(minimum((x-xc).^2 + (y-yc).^2 + (z-zc).^2))
    end

    sort!(ρ)

    # slope estimation
    p = [0.01,0.5]
    p_vec = Int64(round(p[1]*N)):Int64(round(p[2]*N))
    ρ_at_p = ρ[p_vec]

    intercepts[oi],slopes[oi] = linreg(log10.(ρ_at_p),log10.(p_vec))
end

μ = mean(slopes)
μ_std = 2*std(slopes)
β = mean(intercepts)

println("$μ ± $μ_std")

#=
fig,ax1 = subplots(1,1)

ax1[:loglog](ρ,1:N,lw=3,label=L"x^2")
ax1[:loglog](ρ[Int64.(round.(p*N))],p*N,"C1v")
ax1[:loglog](ρ,10^β*ρ.^μ,"C2--")

ax1[:loglog]([Δx,Δx],[1,N],"k")
=#
