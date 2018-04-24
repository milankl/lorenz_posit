using JLD
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport matplotlib.patches as pat

cd("/home/kloewer/julia/lorenz_posit/predict/")
er = load("data/mean_error.jld")["meanerror"]

fig,ax1 = subplots(1,1,figsize=(9,5))

# 32BIT VERSION

DAT5 = load("data/RMSE_16bit_opt_all_s0.1.jld")
DAT4 = load("data/RMSE_15bit_opt_s0.1.jld")
DAT3 = load("data/RMSE_14bit_opt_s0.1.jld")
DAT2 = load("data/RMSE_13bit_opt_s0.1.jld")
DAT1 = load("data/RMSE_12bit_opt_s0.1.jld")

d = 500

FF32 = zeros(5) # lazy: do not rename
PP32 = zeros(5,3)
δ32 = 0.05

for (iD,D) in enumerate([DAT1,DAT2,DAT3,DAT4,DAT5])

    F = zeros(6)
    P = zeros(3,6)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d]/er[2],p)
        P[:,ip] = np.percentile(D["RMSE_P"][2:4,:,d]/er[2],p,axis=1)
    end

    F[6] = mean(D["RMSE_F"][:,d])
    P[:,6] = mean(D["RMSE_P"][2:4,:,d],2)

    # store for accessing outside looponstruction of finite element matr
    FF32[iD],PP32[iD,:] = F[3],P[:,3]

    x = iD-3-δ32
    if iD == 5
        ax1[:plot](x,F[3],"ko")
        ax1[:plot](ones(2)*x,[F[2],F[4]],"k",alpha=.3,lw=2)
    end

    for i = 1:3
        x = iD+0.025*i-3-δ32
        ax1[:plot](x,P[i,3],"C"*"$i"*"o")
        ax1[:plot](ones(2)*x,[P[i,2],P[i,4]],"C"*"$i",alpha=.3,lw=2)
    end
end

#plotting the line connection
ax1[:plot]((-2:length(FF32)-3)-δ32,ones(FF32)*FF32[end],"k",zorder=1)
for i in [1,2,3]
    ax1[:plot]((-2:length(FF32)-3)-δ32+0.025*i,PP32[:,i],"C"*"$i",zorder=1)
end

ax1[:add_patch](pat.Rectangle([-4,1],15,2,color="k",alpha=.2))
ax1[:set_title]("Forecast Error in Lorenz 63, t=$(d/100)",loc="right")

ax1[:set_ylabel]("Error (normalized)")
ax1[:set_xlim](-2.5,3)

ax1[:set_xticks](-2:3)
xtiks = ["12","13","14","15","16"]
ax1[:set_xticklabels](xtiks)
ax1[:set_xlabel]("number of bits")

# LEGEND


ax1[:plot](0,-1,"k",label="Float(16,5)")
ax1[:plot](0,-1,"C1",label="Posit(16,1)")
ax1[:plot](0,-1,"C2",label="Posit(16,2)")
ax1[:plot](0,-1,"C3",label="Posit(16,3)")
ax1[:set_ylim](0,1.1)
ax1[:legend](loc=1)

tight_layout()
savefig("figs/forecast_error_16bits.pdf")
close(fig)
