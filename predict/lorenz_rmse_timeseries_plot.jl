using JLD
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport matplotlib.patches as pat

cd("/home/kloewer/julia/lorenz_posit/predict/")
er = load("data/mean_error.jld")["meanerror"][2]

fig,axs = subplots(3,1,figsize=(9,9),sharex=true,sharey=true)

DAT1 = load("data/RMSE_32bit_s1000.0.jld")
DAT2 = load("data/RMSE_32bit_s1.0.jld")
DAT3 = load("data/RMSE_32bit_s0.001.jld")

N = size(DAT1["RMSE_F"])[2]
Δt = 0.01
time = 0:Δt:(N-1)*Δt
tmax=25

for (iD,D) in enumerate([DAT1,DAT2,DAT3])

    F = zeros(N,3)  # N x 3 different percentiles
    P = zeros(N,3,3) # N x 3 exponent bits x 3 different percentiles

    perct = [25,50,75]

    for (ip,p) in enumerate(perct)
        F[:,ip] = np.percentile(D["RMSE_F"]/er,p,axis=0)
        P[:,:,ip] = transpose(np.percentile(D["RMSE_P"]/er,p,axis=1))
    end

    axs[iD][:plot](time,F[:,2],"k")
    axs[iD][:fill_between](time,F[:,1],F[:,3],color="k",alpha=.1)

    for i = 1:3
        axs[iD][:plot](time,P[:,i,2],"C"*"$i")
        axs[iD][:fill_between](time,P[:,i,1],P[:,i,3],color="C"*"$i",alpha=.1)
    end
end

axs[1][:set_xlim](0,tmax)
axs[1][:set_ylim](0,1.5)

axs[1][:set_title]("Forecast Error in Lorenz 63, at 32bit",loc="right")
axs[3][:set_xlabel]("time")

# ax2[:set_title]("Forecast Error in Lorenz 63, at 32bit",loc="right")
for i=1:3
    axs[i][:set_ylabel]("Error (normalized)")
    axs[i][:plot]([0,tmax],[1,1],"k--",lw=2)
end

axs[1][:text](20,.05,"UPSCALED ×10")
axs[2][:text](20,.05,"UNSCALED ×1")
axs[3][:text](20,.05,"DOWNSCALED ×1/10")

axs[2][:plot](-1,-1,"k",label="Float(32,8)")
axs[2][:plot](-1,-1,"C1",label="Posit(32,1)")
axs[2][:plot](-1,-1,"C2",label="Posit(32,2)")
axs[2][:plot](-1,-1,"C3",label="Posit(32,3)")

axs[2][:legend](loc=2,ncol=2)

tight_layout()
savefig("figs/forecast_error_32bit_timeseries.pdf")
close(fig)
