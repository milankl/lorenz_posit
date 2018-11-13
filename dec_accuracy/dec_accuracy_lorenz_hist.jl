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
ilist = [1,2^(nbits-1)-2,2^(nbits-1)-1]

# arithmetic mean between all representable numbers
f_am = (flist[1:end-1]+flist[2:end])/2.
i_am = ilist[1:end-1]+0.5

# worst case decimal accuracy
f_wda = -log10.(abs.(log10.(f_am./flist[1:end-1])))
i_wda = -log10.(abs.(log10.(i_am./ilist[2:end])))

# extend with zeros due to overflow
f_wda = vcat(0.55,f_wda)        # extrapolate somehow
f_am = vcat(flist[1],f_am)
i_wda = vcat(i_wda)
i_am = vcat(i_am)

function wcdp_posit(nbits,ebits)
    plist = representable_posits(nbits,ebits)
    p_am = (plist[1:end-1]+plist[2:end])/2.
    p_wda = -log10.(abs.(log10.(p_am./plist[1:end-1])))

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

    return p_am,p_wda,plist
end

p0_am, p0_wda, p0list = wcdp_posit(16,0)
p1_am, p1_wda, p1list = wcdp_posit(16,1)
p2_am, p2_wda, p2list = wcdp_posit(16,2)

# calculate histogram\mathnormal{
σ,β,ρ = 10.,8/3,28.
Δt = 0.01
bins = 10.0.^(-8:0.1:9)

H1 = lorenz_hist_opt(load("data/lorenz_scale1.jld")["xyz"],1.,σ,β,ρ,Δt,bins)
H10 = lorenz_hist_opt(load("data/lorenz_scale10.jld")["xyz"],10.,σ,β,ρ,Δt,bins)
H100 = lorenz_hist_opt(load("data/lorenz_scale100.jld")["xyz"],100.,σ,β,ρ,Δt,bins)
H_10 = lorenz_hist_opt(load("data/lorenz_scale-10.jld")["xyz"],1/10.,σ,β,ρ,Δt,bins)
H_100 = lorenz_hist_opt(load("data/lorenz_scale-100.jld")["xyz"],1/100.,σ,β,ρ,Δt,bins)

## PLOTTING
fig,(ax1,ax2) = subplots(2,1,figsize=(6,6),sharex=false)

ax1[:plot](f_am,f_wda,"k",label="Float$nbits",lw=2)
ax1[:plot](p0_am,p0_wda,"C1",label="Posit($nbits,0)",lw=1.4)
ax1[:plot](p1_am,p1_wda,"C2",label="Posit($nbits,1)",lw=1.2)
ax1[:plot](p2_am,p2_wda,"C3",label="Posit($nbits,2)",lw=0.8)
ax1[:plot](i_am,i_wda,"C0",label="Int$nbits",lw=2)

ax1[:fill_between](f_am,-0.1,f_wda,edgecolor="k",facecolor="none",linestyle="--")
ax1[:fill_between](i_am,-0.1,i_wda,edgecolor="C0",facecolor="none",linestyle="--")
ax1[:fill_between](p0_am,-0.1,p0_wda,where=((p0_am .>= p0list[1]).*(p0_am .<= p0list[end])),edgecolor="C1",facecolor="none",linestyle="--")
ax1[:fill_between](p1_am,-0.1,p1_wda,where=((p1_am .>= p1list[1]).*(p1_am .<= p1list[end])),edgecolor="C2",facecolor="none",linestyle="--")
ax1[:fill_between](p2_am,-0.1,p2_wda,where=((p2_am .>= p2list[1]).*(p2_am .<= p2list[end])),edgecolor="C3",facecolor="none",linestyle="--")
#ax1[:fill_between](p_am,0.,p_wda,where=((p_am .<= plist[1]).|(p_am .>= plist[end])),facecolor="C1",alpha=0.1)

ax1[:legend](loc=2,fontsize=9)

ax1[:set_xlim](1e-16,1e16)
ax1[:set_xscale]("log",basex=10)
ax1[:set_ylim](0,6)


ax2[:set_xlabel]("value")
ax1[:set_ylabel]("decimal places")
ax2[:set_ylabel](L"$N$ [$10^4$]")


# ax11 = ax1[:twiny]()
# ax11[:set_xscale]("log",basex=10)
# ax11[:set_xlim](ax1[:get_xlim]())

xtik = 10.0.^(-16:4:16)
#ax1[:set_xticks](xtik)

ax2[:set_xscale]("log",basex=10)
ax2[:set_xlim](1e-16,1e16)
ax2[:set_xticks](xtik)
ax1[:set_xticks](xtik)

bins[1] = 1e-16
bins[end-1] = 1e16

ax2[:plot](bins[1:end-1],H100/1e4,"C8",label=L"s=10^2",drawstyle="steps-post")
#ax2[:plot](bins[1:end-1],H10/1e4,"C7",label=L"s=10^1",drawstyle="steps-post")
ax2[:plot](bins[1:end-1],H1/1e4,"C6",label=L"s=10^0",drawstyle="steps-post")
#ax2[:plot](bins[1:end-1],H_10/1e4,"C4",label=L"s=10^{-1}",drawstyle="steps-post")
ax2[:plot](bins[1:end-1],H_100/1e4,"C5",label=L"s=10^{-2}",drawstyle="steps-post")

ax1[:set_title]("Decimal precision",loc="left")
ax2[:set_title]("Numbers subject to rounding errors in rescaled L63",loc="left")

ax1[:set_title]("a",loc="right",fontweight="bold")
ax2[:set_title]("b",loc="right",fontweight="bold")


ax2[:legend](loc=1,fontsize=9)

tight_layout()
savefig("/home/kloewer/julia/lorenz_posit/dec_accuracy/figs/dec_acc_hist.png",dpi=300)
close(fig)
