using JLD
using PyPlot

cd("/home/kloewer/julia/lorenz_posit/frac_dim/")

# loading data
Fdims = load("data/frac_dims_floats.jld")
Bdims = load("data/frac_dims_bigfloats.jld")

Pdims1 = load("data/frac_dims_posits_1ebit.jld")
Pdims2 = load("data/frac_dims_posits_2ebit.jld")
Pdims3 = load("data/frac_dims_posits_3ebit.jld")


# adjust nbits of Bdims
Bnbits = 10:32


fig,ax = subplots(1,1,figsize=(8.5,4))
ms = 8  # markersize

ax[:plot]([0,66],[2.06,2.06],"k",lw=2,label="D = 2.06 ± 0.01")

ax[:plot](Pdims1["P_nbits"],Pdims1["fracdim"],"C2.-",lw=.5,ms=ms,alpha=.6,label="Posits 1ebit")
ax[:plot](Pdims2["P_nbits"],Pdims2["fracdim"],"C1.-",lw=.5,ms=ms,alpha=.6,label="Posits 2ebits")
ax[:plot](Pdims3["P_nbits"],Pdims3["fracdim"],"C3.-",lw=.5,ms=ms,alpha=.6,label="Posits 3ebits")

ax[:plot](Bnbits,Bdims["fracdim"],"k.-",ms=ms,lw=.5,alpha=.6,label="Floats (BigFloat)")
ax[:plot](Fdims["nbits"],Fdims["fracdim"][:,1],"k<",ms=ms,alpha=.6,label="Floats (native)")

ax[:set_yticks]([0,.5,1,1.5,2])
ax[:set_xticks](10:32)
ax[:set_xticklabels](10:32,rotation=90)

ax[:text](16,0,"HALF",rotation=90,ha="center",va="bottom")
ax[:text](32,0,"SINGLE",rotation=90,ha="center",va="bottom")

ax[:set_title]("Fractal Dimension Lorenz 63")
ax[:set_ylabel]("fractal dimension D")
ax[:set_xlabel]("number of bits")

ax[:set_xlim](9.5,33)
ax[:set_ylim](-.1,2.1)

ax[:legend](loc=8)

tight_layout()
savefig("figs/fractal_dim.pdf")
close(fig)
