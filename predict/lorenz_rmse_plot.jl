using JLD
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport matplotlib.patches as pat

cd("/home/kloewer/julia/lorenz_posit/predict/")
er = load("data/mean_error.jld")["meanerror"]

fig,(ax1,ax2) = subplots(2,1,figsize=(9,9),sharex=true)

# 16BIT VERSION

DAT1 = load("data/RMSE_16bit_s-10.jld")
DAT2 = load("data/RMSE_16bit_s1.jld")
DAT3 = load("data/RMSE_16bit_s10.jld")

d = 400

for (iD,D) in enumerate([DAT1,DAT2,DAT3])

    F = zeros(6)
    P = zeros(3,6)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d]/er[iD],p)
        P[:,ip] = np.percentile(D["RMSE_P"][:,:,d]/er[iD],p,axis=1)
    end

    F[6] = mean(D["RMSE_F"][:,d])
    P[:,6] = mean(D["RMSE_P"][:,:,d],2)

    ax1[:errorbar](F[3],13-iD,xerr=reshape(abs.(F[[1,5]]-F[3]),2,1),color="k",fmt="--",lw=.5)
    ax1[:errorbar](F[3],13-iD,xerr=reshape(abs.(F[[2,4]]-F[3]),2,1),color="k",fmt="s--",capsize=5,lw=2)

    for i = 1:3
        ax1[:errorbar](P[i,3],13-3i-iD,xerr=reshape(abs.(P[i,[1,5]]-P[i,3]),2,1),color="C"*"$i",fmt="--",lw=.5)
        ax1[:errorbar](P[i,3],13-3i-iD,xerr=reshape(abs.(P[i,[2,4]]-P[i,3]),2,1),color="C"*"$i",fmt="s--",capsize=5,lw=2)
    end
end

# 32BIT VERSION

DAT1 = load("data/RMSE_32bit_s-10.jld")
DAT2 = load("data/RMSE_32bit_s1.jld")
DAT3 = load("data/RMSE_32bit_s10.jld")

d = 1400

for (iD,D) in enumerate([DAT1,DAT2,DAT3])

    F = zeros(6)
    P = zeros(3,6)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d]/er[iD],p)
        P[:,ip] = np.percentile(D["RMSE_P"][:,:,d]/er[iD],p,axis=1)
    end

    F[6] = mean(D["RMSE_F"][:,d])
    P[:,6] = mean(D["RMSE_P"][:,:,d],2)

    ax2[:errorbar](F[3],13-iD,xerr=reshape(abs.(F[[1,5]]-F[3]),2,1),color="k",fmt="--",lw=.5)
    ax2[:errorbar](F[3],13-iD,xerr=reshape(abs.(F[[2,4]]-F[3]),2,1),color="k",fmt="s--",capsize=5,lw=2)

    for i = 1:3
        ax2[:errorbar](P[i,3],13-3*i-iD,xerr=reshape(abs.(P[i,[1,5]]-P[i,3]),2,1),color="C"*"$i",fmt="--",lw=.5)
        ax2[:errorbar](P[i,3],13-3*i-iD,xerr=reshape(abs.(P[i,[2,4]]-P[i,3]),2,1),color="C"*"$i",fmt="s--",capsize=5,lw=2)
    end
end

ax1[:add_patch](pat.Rectangle([1,0],2,15,color="k",alpha=.2))
ax2[:add_patch](pat.Rectangle([1,0],2,15,color="k",alpha=.2))

ax1[:set_title]("Forecast Error in Lorenz 63, at 16bit",loc="right")
ax2[:set_title]("Forecast Error in Lorenz 63, at 32bit",loc="right")

xlabel("Error (normalized)")

ax1[:set_yticks]([2,5,8,11])
ax1[:set_yticklabels](["Posit(16,2)","Posit(16,1)","Posit(16,0)","Float(16,5)"])

ax2[:set_yticks]([2,5,8,11])
ax2[:set_yticklabels](["Posit(32,3)","Posit(32,2)","Posit(32,1)","Float(32,8)"])

ax1[:set_xlim](-.1,1.6)
ax2[:set_xticks](0:0.25:1.5)

tight_layout()
savefig("figs/forecast_error.pdf")
close(fig)
