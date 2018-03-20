using SigmoidNumbers
using PyCall
@pyimport matplotlib.pyplot as plt

include("repr_numbers.jl")

nbits = 32

f_am5,f_wda5 = wc_dec_acc_float(nbits,5)
f_am8,f_wda8 = wc_dec_acc_float(nbits,8)

p_am1,p_wda1 = wc_dec_acc_posit(nbits,1)
p_am2,p_wda2 = wc_dec_acc_posit(nbits,2)

# PLOTTING
fig,ax1 = plt.subplots(1,1,sharex=true)

ax1[:plot](f_am5,f_wda5,"grey",label="Floats($nbits,5)")
ax1[:plot](f_am8,f_wda8,"k",label="Floats($nbits,8)")

ax1[:plot](p_am1,p_wda1,"C1",label="Posits($nbits,1)")
ax1[:plot](p_am2,p_wda2,"C2",label="Posits($nbits,2)")

plt.legend(loc=1)

ax1[:set_xlim](10e-80,10e80)
ax1[:set_xscale]("log",basex=2)
ax1[:set_ylim](2,10)


ax1[:set_xlabel]("x")
ax1[:set_ylabel]("Worst-case decimal accuracy")

ax2 = ax1[:twiny]()
ax2[:set_xscale]("log",basex=10)
ax2[:set_xlim](ax1[:get_xlim]())

plt.tight_layout()
plt.show()
