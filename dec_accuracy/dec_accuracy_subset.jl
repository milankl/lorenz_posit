using SigmoidNumbers
using PyPlot

include("repr_numbers.jl")

nbits = 32

f_am5,f_wda5 = wc_dec_acc_float(nbits,5)
f_am8,f_wda8 = wc_dec_acc_float(nbits,8)

p_am1,p_wda1 = wc_dec_acc_posit(nbits,1)
p_am2,p_wda2 = wc_dec_acc_posit(nbits,2)

i_am,i_wda = wc_dec_acc_int(nbits)

# PLOTTING
fig,ax1 = subplots(1,1,sharex=true)

#ax1[:plot](f_am5,f_wda5,"grey",label="Floats($nbits,5)")
ax1[:plot](f_am8,f_wda8,"C0",label="Floats($nbits,8)")

#ax1[:plot](p_am1,p_wda1,"C1",label="Posits($nbits,1)")
ax1[:plot](p_am2,p_wda2,"C1",label="Posits($nbits,2)")
ax1[:plot](i_am,i_wda,"C2",label="Int$nbits")

ax1[:fill_between](f_am8,0.,f_wda8,color="C0",alpha=.3)
ax1[:fill_between](p_am2,0.,p_wda2,color="C1",alpha=.3)
ax1[:fill_between](i_am,0.,i_wda,color="C2",alpha=.3)


legend(loc=1)

ax1[:set_xlim](10e-50,10e50)
ax1[:set_xscale]("log",basex=2)
ax1[:set_ylim](2,10)

ax1[:set_xlabel]("x")
ax1[:set_ylabel]("Worst-case decimal accuracy")

ax2 = ax1[:twiny]()
ax2[:set_xscale]("log",basex=10)
ax2[:set_xlim](ax1[:get_xlim]())

tight_layout()
