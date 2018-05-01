using JLD
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport matplotlib.patches as pat

cd("/home/kloewer/julia/lorenz_posit/predict/")
er = load("data/mean_error.jld")["meanerror"]

fig,(ax1,ax2) = subplots(2,1,figsize=(9,9))

# 16BIT

DAT12 = load("data/fullrhs/RMSE_16bit_opt_s1.0e7.jld")
DAT11 = load("data/fullrhs/RMSE_16bit_opt_s1.0e6.jld")
DAT10 = load("data/fullrhs/RMSE_16bit_opt_s100000.0.jld")
DAT9 = load("data/fullrhs/RMSE_16bit_opt_s10000.0.jld")
DAT8 = load("data/fullrhs/RMSE_16bit_opt_s600.0.jld")
DAT7 = load("data/fullrhs/RMSE_16bit_opt_s100.0.jld")
DAT6 = load("data/fullrhs/RMSE_16bit_opt_s10.0.jld")
DAT5 = load("data/fullrhs/RMSE_16bit_opt_s1.0.jld")
DAT4 = load("data/fullrhs/RMSE_16bit_opt_s0.1.jld")
DAT3 = load("data/fullrhs/RMSE_16bit_opt_s0.01.jld")
DAT2 = load("data/fullrhs/RMSE_16bit_opt_s0.001.jld")
DAT1 = load("data/fullrhs/RMSE_16bit_opt_s0.0001.jld")

d16 = 800

FF16 = zeros(12)
PP16 = zeros(12,3)
I16 = zeros(12)
δ32 = 0.05

for (iD,D) in enumerate([DAT1,DAT2,DAT3,DAT4,DAT5,DAT6,DAT7,DAT8,DAT9,DAT10,DAT11,DAT12])

    F = zeros(5)
    P = zeros(3,5)
    I = zeros(5)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d16]/er[2],p)
        I[ip] = np.percentile(D["RMSE_I"][:,d16]/er[2],p)
        P[:,ip] = np.percentile(D["RMSE_P"][:,:,d16]/er[2],p,axis=1)
    end

    # store for accessing outside looponstruction of finite element matr
    FF16[iD],PP16[iD,:],I16[iD,:] = F[3],P[:,3],I[3]

    x = iD-3-δ32
    ax1[:plot](x,F[3],"ko")
    ax1[:plot](ones(2)*x,[F[2],F[4]],"k",alpha=.3,lw=2)
    if iD < 9
        ax1[:plot](x-δ32,I[3],"C0o")
        ax1[:plot](ones(2)*(x-δ32),[I[2],I[4]],"C0",alpha=.3,lw=2)
    end

    for i = 1:3
        x = iD+0.05*i-3-δ32
        ax1[:plot](x,P[i,3],"C"*"$i"*"o")
        ax1[:plot](ones(2)*x,[P[i,2],P[i,4]],"C"*"$i",alpha=.3,lw=2)
    end
end

# 32BIT VERSION

DAT13 = load("data/fullrhs/RMSE_32bit_opt_s4.0e7.jld")
DAT12 = load("data/fullrhs/RMSE_32bit_opt_s1.0e7.jld")
DAT11 = load("data/fullrhs/RMSE_32bit_opt_s1.0e6.jld")
DAT10 = load("data/fullrhs/RMSE_32bit_opt_s100000.0.jld")
DAT9 = load("data/fullrhs/RMSE_32bit_opt_s10000.0.jld")
DAT8 = load("data/fullrhs/RMSE_32bit_opt_s1000.0.jld")
DAT7 = load("data/fullrhs/RMSE_32bit_opt_s100.0.jld")
DAT6 = load("data/fullrhs/RMSE_32bit_opt_s10.0.jld")
DAT5 = load("data/fullrhs/RMSE_32bit_opt_s1.0.jld")
DAT4 = load("data/fullrhs/RMSE_32bit_opt_s0.1.jld")
DAT3 = load("data/fullrhs/RMSE_32bit_opt_s0.01.jld")
DAT2 = load("data/fullrhs/RMSE_32bit_opt_s0.001.jld")
DAT1 = load("data/fullrhs/RMSE_32bit_opt_s0.0001.jld")

d32 = 2000

FF32 = zeros(13)
PP32 = zeros(13,3)
I32 = zeros(13)
δ32 = 0.05

for (iD,D) in enumerate([DAT1,DAT2,DAT3,DAT4,DAT5,DAT6,DAT7,DAT8,DAT9,DAT10,DAT11,DAT12,DAT13])

    F = zeros(5)
    P = zeros(3,5)
    I = zeros(5)

    perct = [10,25,50,75,90]

    for (ip,p) in enumerate(perct)
        F[ip] = np.percentile(D["RMSE_F"][:,d32]/er[2],p)
        I[ip] = np.percentile(D["RMSE_I"][:,d32]/er[2],p)
        P[:,ip] = np.percentile(D["RMSE_P"][:,:,d32]/er[2],p,axis=1)
    end

    # store for accessing outside looponstruction of finite element matr
    FF32[iD],PP32[iD,:],I32[iD,:] = F[3],P[:,3],I[3]

    x = iD-3-δ32
    ax2[:plot](x,F[3],"ko")
    ax2[:plot](ones(2)*x,[F[2],F[4]],"k",alpha=.3,lw=2)
    ax2[:plot](x-δ32,I[3],"C0o")
    ax2[:plot](ones(2)*(x-δ32),[I[2],I[4]],"C0",alpha=.3,lw=2)

    for i = 1:3
        x = iD+0.05*i-3-δ32
        ax2[:plot](x,P[i,3],"C"*"$i"*"o")
        ax2[:plot](ones(2)*x,[P[i,2],P[i,4]],"C"*"$i",alpha=.3,lw=2)
    end
end

#plotting the line connection
ax1[:plot]((-2:length(FF16)-3)-δ32,FF16,"k",zorder=1)
for i in [1,2,3]
    ax1[:plot]((-2:length(FF16)-3)-δ32+0.05*i,PP16[:,i],"C"*"$i",zorder=1)
end
x = Array(-2:length(I16)-3)-2*δ32
ax1[:plot](x[1:end-4],I16[1:end-4],"C0")

#plotting the line connection
ax2[:plot]((-2:length(FF32)-3)-δ32,FF32,"k",zorder=1)
for i in [1,2,3]
    ax2[:plot]((-2:length(FF32)-3)-δ32+0.05*i,PP32[:,i],"C"*"$i",zorder=1)
end
x = Array(-2:length(I32)-3)-2*δ32
ax2[:plot](x,I32,"C0")

ax1[:add_patch](pat.Rectangle([-4,1],15,2,color="k",alpha=.2))
ax1[:set_title]("Forecast Error in Lorenz 63, 16bit, RHS at Float64, t=$(d16/100)",loc="right")

ax2[:add_patch](pat.Rectangle([-4,1],15,2,color="k",alpha=.2))
ax2[:set_title]("Forecast Error in Lorenz 63, 32bit, RHS at Float64, t=$(d32/100)",loc="right")

ax1[:set_ylabel]("Error (normalized)")
ax2[:set_ylabel]("Error (normalized)")
ax1[:set_xlim](-2.5,10.5)
ax2[:set_xlim](-2.5,10.5)

ax1[:set_xticks](-2:9)
ax2[:set_xticks](-2:10)
xtiks1 = [L"×10^7",L"×10^6",L"×10^5",L"×10^4",L"×600",L"×10^2",L"×10^1",L"×1",L"×10^{-1}",L"×10^{-2}",L"×10^{-3}",L"×10^{-4}"][end:-1:1]
xtiks2 = [L"×4\cdot 10^7",L"×10^7",L"×10^6",L"×10^5",L"×10^4",L"×10^3",L"×10^2",L"×10^1",L"×1",L"×10^{-1}",L"×10^{-2}",L"×10^{-3}",L"×10^{-4}"][end:-1:1]
ax1[:set_xticklabels](xtiks1)
ax2[:set_xticklabels](xtiks2)
ax2[:set_xlabel]("Scaling")

# LEGEND
ax1[:plot](0,-1,"k",label="Float(16,5)")
ax1[:plot](0,-1,"C1",label="Posit(16,1)")
ax1[:plot](0,-1,"C2",label="Posit(16,2)")
ax1[:plot](0,-1,"C3",label="Posit(16,3)")
ax1[:plot](0,-1,"C0",label="Int16")
ax1[:set_ylim](0,1.1)
ax1[:legend](loc=1)

# LEGEND
ax2[:plot](0,-1,"k",label="Float(32,8)")
ax2[:plot](0,-1,"C1",label="Posit(32,1)")
ax2[:plot](0,-1,"C2",label="Posit(32,2)")
ax2[:plot](0,-1,"C3",label="Posit(32,3)")
ax2[:plot](0,-1,"C0",label="Int32")
ax2[:set_ylim](0,1.1)
ax2[:legend](loc=1)

tight_layout()
savefig("figs/forecast_error_fullrhs.pdf")
close(fig)
