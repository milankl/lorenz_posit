using SigmoidNumbers
using PyPlot

include("repr_numbers.jl")

nbits = 16
febits = 5
pebits = 1

flist = representable_floats(nbits,febits)
plist = representable_posits(nbits,pebits)
ilist = (1:(2^(nbits-1)-1))

# arithmetic mean between all representable numbers
f_am = (flist[1:end-1]+flist[2:end])/2.
p_am = (plist[1:end-1]+plist[2:end])/2.
i_am = ilist[2:end]

# worst case decimal accuracy
f_wda = -log10.(abs.(log10.(f_am./flist[2:end])))
p_wda = -log10.(abs.(log10.(p_am./plist[2:end])))
i_wda = -log10.(abs.(log10.(i_am./ilist[1:end-1])))

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

i_wda = vcat(0.13,i_wda)    # somehow interpolate to 1...
i_am = vcat(1.,i_am)

# extend with zeros due to under/overflow
f_wda = vcat(0.73,f_wda)    # somehow interpolate to 1...
f_am = vcat(flist[1],f_am)

# PLOTTING
fig,ax1 = subplots(1,1,sharex=true,figsize=(9,4))

ax1[:plot](f_am,f_wda,"k",label="Float($nbits,$febits)")
ax1[:plot](p_am,p_wda,"C1",label="Posit($nbits,$pebits)")
ax1[:plot](i_am,i_wda,"C2",label="Int$nbits")


ax1[:fill_between](f_am,0.,f_wda,facecolor="k",alpha=0.3)
ax1[:fill_between](p_am,0.,p_wda,where=((p_am .>= plist[1]).*(p_am .<= plist[end])),facecolor="C1",alpha=0.3)
ax1[:fill_between](p_am,0.,p_wda,where=((p_am .<= plist[1]).|(p_am .>= plist[end])),facecolor="C1",alpha=0.1)
ax1[:fill_between](i_am,0.,i_wda,facecolor="C2",alpha=0.3)

legend()

ax1[:set_xscale]("log",basex=2)
ax1[:set_ylim](0,5.5)

ax1[:set_xlabel]("x")
ax1[:set_ylabel]("Worst-case decimal precision")

ax2 = ax1[:twiny]()
ax2[:set_xscale]("log",basex=10)
ax2[:set_xlim](ax1[:get_xlim]())
ax2[:set_xticks](10.0.^(-8:2:8))

tight_layout()
savefig("figs/dec_acc_16bit.pdf")
close(fig)
