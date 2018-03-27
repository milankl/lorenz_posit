using JLD
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport matplotlib.patches as pat

cd("/home/kloewer/julia/lorenz_posit/predict/")
er = load("data/mean_error.jld")["meanerror"]

fig,(ax1,ax2) = subplots(2,1,figsize=(9,9),sharex=true,sharey=true)

# 16BIT VERSION

DAT1 = load("data/RMSE_16bit_s1000.0.jld")
DAT2 = load("data/RMSE_16bit_s100.0.jld")
DAT3 = load("data/RMSE_16bit_s10.0.jld")
DAT4 = load("data/RMSE_16bit_s1.0.jld")
DAT5 = load("data/RMSE_16bit_s0.1.jld")
DAT6 = load("data/RMSE_16bit_s0.01.jld")
DAT7 = load("data/RMSE_16bit_s0.001.jld")

d = 450
δ16 = 0.05     #shift

FF16 = zeros(7)
PP16 = zeros(7,3)

for (iD,D) in enumerate([DAT1,DAT2,DAT3,DAT4,DAT5,DAT6,DAT7])

    F = zeros(6)
    P = zeros(3,6)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d]/er[2],p)
        P[:,ip] = np.percentile(D["RMSE_P"][:,:,d]/er[2],p,axis=1)
    end

    F[6] = mean(D["RMSE_F"][:,d])
    P[:,6] = mean(D["RMSE_P"][:,:,d],2)

    # store for accessing outside loop
    FF16[iD],PP16[iD,:] = F[3],P[:,3]

    x = iD-δ16-2
    ax1[:plot](x,F[3],"ko")
    ax1[:plot](ones(2)*x,[F[2],F[4]],"k",alpha=.3,lw=2)

    for i = 1:3
        x = iD+.05*(i-1)-δ16-2
        ax1[:plot](x,P[i,3],"C"*"$i"*"o")
        ax1[:plot](ones(2)*x,[P[i,2],P[i,4]],"C"*"$i",alpha=.3,lw=2)
    end
end

ax1[:set_title]("Forecast Error in Lorenz 63, at 16bit, t=$(d/100)",loc="right")

# 32BIT VERSION

DAT1 = load("data/RMSE_32bit_s10000.0.jld")
DAT2 = load("data/RMSE_32bit_s1000.0.jld")
DAT3 = load("data/RMSE_32bit_s100.0.jld")
DAT4 = load("data/RMSE_32bit_s10.0.jld")
DAT5 = load("data/RMSE_32bit_s1.0.jld")
DAT6 = load("data/RMSE_32bit_s0.1.jld")
DAT7 = load("data/RMSE_32bit_s0.01.jld")
DAT8 = load("data/RMSE_32bit_s0.001.jld")
DAT9 = load("data/RMSE_32bit_s0.0001.jld")

d = 1550
δ32 = 0.1     #shift

FF32 = zeros(9)
PP32 = zeros(9,3)

for (iD,D) in enumerate([DAT1,DAT2,DAT3,DAT4,DAT5,DAT6,DAT7,DAT8,DAT9])

    F = zeros(6)
    P = zeros(3,6)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d]/er[2],p)
        P[:,ip] = np.percentile(D["RMSE_P"][:,:,d]/er[2],p,axis=1)
    end

    F[6] = mean(D["RMSE_F"][:,d])
    P[:,6] = mean(D["RMSE_P"][:,:,d],2)

    # store for accessing outside loop
    FF32[iD],PP32[iD,:] = F[3],P[:,3]

    x = iD-3-δ32
    ax2[:plot](x,F[3],"ko")
    ax2[:plot](ones(2)*x,[F[2],F[4]],"k",alpha=.3,lw=2)

    for i = 1:3
        x = iD+0.05*i-3-δ32
        ax2[:plot](x,P[i,3],"C"*"$i"*"o")
        ax2[:plot](ones(2)*x,[P[i,2],P[i,4]],"C"*"$i",alpha=.3,lw=2)
    end
end

#plotting the line connection
ax1[:plot]((-1:length(FF16)-2)-δ16,FF16,"k")
for i in [1,2,3]
    ax1[:plot]((-1:length(FF16)-2)-δ16+0.05*(i-1),PP16[:,i],"C"*"$i")
end

ax2[:plot]((-2:length(FF32)-3)-δ32,FF32,"k")
for i in [1,2,3]
    ax2[:plot]((-2:length(FF32)-3)-δ32+0.05*i,PP32[:,i],"C"*"$i")
end


ax1[:add_patch](pat.Rectangle([-4,1],15,2,color="k",alpha=.2))
ax2[:add_patch](pat.Rectangle([-4,1],15,2,color="k",alpha=.2))

ax2[:set_title]("Forecast Error in Lorenz 63, at 32bit, t=$(d/100)",loc="right")

ax1[:set_ylabel]("Error (normalized)")
ax2[:set_ylabel]("Error (normalized)")

ax1[:set_xlim](-2.5,6.5)
ax2[:set_xlim](-2.5,6.5)

ax1[:set_xticks](-2:6)
ax2[:set_xticks](-2:6)
ax2[:set_xticklabels]([L"×10^4",L"×10^3",L"×10^2",L"×10^1",L"×1",L"×10^{-1}",L"×10^{-2}",L"×10^{-3}",L"×10^{-4}"])
ax2[:set_xlabel]("Scaling")

# LEGEND
ax1[:plot](0,-1,"k",label="Float(16,5)")
ax1[:plot](0,-1,"C1",label="Posit(16,1)")
ax1[:plot](0,-1,"C2",label="Posit(16,2)")
ax1[:plot](0,-1,"C3",label="Posit(16,3)")

ax1[:legend](loc=3)

ax2[:plot](0,-1,"k",label="Float(32,8)")
ax2[:plot](0,-1,"C1",label="Posit(32,1)")
ax2[:plot](0,-1,"C2",label="Posit(32,2)")
ax2[:plot](0,-1,"C3",label="Posit(32,3)")

ax2[:legend](loc=9)

ax1[:set_ylim](0,1.1)

tight_layout()
savefig("figs/forecast_error.pdf")
close(fig)
