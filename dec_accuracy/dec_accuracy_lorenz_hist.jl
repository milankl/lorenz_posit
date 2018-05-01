using SigmoidNumbers
using JLD
using PyPlot
using PyCall
@pyimport numpy as np

cd("/home/kloewer/julia/lorenz_posit/dec_accuracy/")
include("repr_numbers.jl")
include("lorenz_hist_functions.jl")

nbits = 16
febits = 5
pebits = 1

flist = representable_floats(nbits,febits)
plist = representable_posits(nbits,pebits)
#ilist = 1:(2^(nbits-1)-1)
ilist = [1,2^(nbits-1)-2,2^(nbits-1)-1]

# arithmetic mean between all representable numbers
f_am = (flist[1:end-1]+flist[2:end])/2.
p_am = (plist[1:end-1]+plist[2:end])/2.
i_am = ilist[1:end-1]+0.5

# worst case decimal accuracy
f_wda = -log10.(abs.(log10.(f_am./flist[1:end-1])))
p_wda = -log10.(abs.(log10.(p_am./plist[1:end-1])))
i_wda = -log10.(abs.(log10.(i_am./ilist[2:end])))

# extend first and last point, taking no overflow/underflow into account
p0 = plist[1]/16    # something much smaller than minpos
pinf = plist[end]*16    # something much bigger than maxpos

p_wda_0 = -log10.(abs.(log10.(p0/plist[2]))) # worst-case decimal accuracy for these extreme values
p_wda_inf = -log10.(abs.(log10.(pinf/plist[end])))

# worst-case decimal accuracy of interpolated on minpos/maxpos
p_wda_minpos = p_wda_0 + log10(plist[1]/p0)/log10(p0/p_am[1])*(p_wda_0-p_wda[1])
p_wda_maxpos = p_wda_inf + log10(plist[end]/pinf)/log10(pinf/p_am[end])*(p_wda_inf-p_wda[end])

p_wda = vcat(p_wda_0,p_wda_minpos,p_wda,p_wda_maxpos,p_wda_inf)
p_am = vcat(p0,plist[1],p_am,plist[end],pinf)

# extend with zeros due to overflow
f_wda = vcat(0.55,f_wda)        # extrapolate somehow
f_am = vcat(flist[1],f_am)
i_wda = vcat(i_wda)
i_am = vcat(i_am)

# calculate histogram
σ,β,ρ = 10.,8./3.,28.
Δt = 0.01
bins = 10.0.^(-8:0.1:9)

H1 = lorenz_hist_opt(load("data/lorenz_scale1.jld")["xyz"],1.,σ,β,ρ,Δt,bins)
H10 = lorenz_hist_opt(load("data/lorenz_scale10.jld")["xyz"],10.,σ,β,ρ,Δt,bins)
H100 = lorenz_hist_opt(load("data/lorenz_scale100.jld")["xyz"],100.,σ,β,ρ,Δt,bins)
H_10 = lorenz_hist_opt(load("data/lorenz_scale-10.jld")["xyz"],1/10.,σ,β,ρ,Δt,bins)
H_100 = lorenz_hist_opt(load("data/lorenz_scale-100.jld")["xyz"],1/100.,σ,β,ρ,Δt,bins)

# PLOTTING
fig,(ax1,ax2) = subplots(2,1,figsize=(8,6),sharex=true)

ax1[:plot](f_am,f_wda,"k",label="Float($nbits,$febits)")
ax1[:plot](p_am,p_wda,"C1",label="Posit($nbits,$pebits)")
ax1[:plot](i_am,i_wda,"C2",label="Int$nbits")

ax1[:fill_between](f_am,0.,f_wda,facecolor="k",alpha=0.3)
ax1[:fill_between](i_am,0.,i_wda,facecolor="C2",alpha=0.3)
ax1[:fill_between](p_am,0.,p_wda,where=((p_am .>= plist[1]).*(p_am .<= plist[end])),facecolor="C1",alpha=0.3)
ax1[:fill_between](p_am,0.,p_wda,where=((p_am .<= plist[1]).|(p_am .>= plist[end])),facecolor="C1",alpha=0.1)

ax1[:legend](loc=1)

ax1[:set_xscale]("log",basex=10)
ax1[:set_ylim](0,6)

ax2[:set_xlabel]("x")
ax1[:set_ylabel]("Worst-case decimal precision")
ax2[:set_ylabel]("N")

# ax11 = ax1[:twiny]()
# ax11[:set_xscale]("log",basex=10)
# ax11[:set_xlim](ax1[:get_xlim]())

xtik = 10.0.^(-10:2:10)
ax11[:set_xticks](xtik)

ax2[:set_xscale]("log",basex=10)
ax2[:set_xlim](ax1[:get_xlim]())
ax2[:set_xticks](xtik)

# ax2[:plot](bins[1:end-1],H100,label=L"s=10^2",drawstyle="steps-post")
ax2[:plot](bins[1:end-1],H10,"C0",label=L"s=10^1",drawstyle="steps-post")
ax2[:plot](bins[1:end-1],H1,"C3",label=L"s=10^0",drawstyle="steps-post")
ax2[:plot](bins[1:end-1],H_10,"C4",label=L"s=10^{-1}",drawstyle="steps-post")
# ax2[:plot](bins[1:end-1],H_100,label=L"s=10^{-2}",drawstyle="steps-post")


ax2[:set_title]("Numbers subject to rounding errors in rescaled Lorenz 63")

ax2[:legend](loc=1)

tight_layout()
savefig("figs/dec_acc_hist_opt.pdf")
close(fig)
