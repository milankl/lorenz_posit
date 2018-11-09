using SigmoidNumbers
using PyPlot

include("lorenz_integration_posits.jl")

Nt0 = 10000      # length of spin-up
Nt = 100000
Δt = 0.01    # time step

σ = 10.      # Lorenz 63 coefficients
β1 = 8
β2 = 3
β = 8./3.
ρ = 28.
s = 0.1

# start somewhere
xyz = [0.6,0.5,15.]

P = Posit{16,1}

# do integration
xyz0 = time_integration_opt(Nt0,Float64,xyz,σ,ρ,β,s,Δt)  # spin-up in full precision

xyzT = time_integration_opt(Nt,Float64,xyz0[:,end],σ,ρ,β,1,Δt)
xyzP = time_integration_opt(Nt,P,xyz0[:,end],σ,ρ,β,s,Δt)
xyzF = time_integration_opt(Nt,Float16,xyz0[:,end],σ,ρ,β,1,Δt)
xyzI = time_integration_int(Nt,Int16,xyz0[:,end],σ,ρ,β1,β2,100.,Δt)

## PLOTTING
fig,axs = subplots(2,2,figsize=(7,6),sharex=true,sharey=true)

axs[1,1][:plot](xyzT[1,:],xyzT[3,:],"k",lw=0.2)
axs[1,2][:plot](xyzF[1,:],xyzF[3,:],"k",lw=0.2)
axs[2,1][:plot](xyzP[1,:],xyzP[3,:],"k",lw=0.2)
axs[2,2][:plot](xyzI[1,:],xyzI[3,:],"k",lw=0.2)

axs[2,1][:set_xlim](-20,20)
axs[2,1][:set_ylim](0,50)

axs[2,1][:set_xlabel](L"x")
axs[2,2][:set_xlabel](L"x")
axs[2,1][:set_ylabel](L"z")
axs[1,1][:set_ylabel](L"z")

axs[1,1][:set_title]("a",loc="right",fontweight="bold")
axs[1,2][:set_title]("b",loc="right",fontweight="bold")
axs[2,1][:set_title]("c",loc="right",fontweight="bold")
axs[2,2][:set_title]("d",loc="right",fontweight="bold")


axs[1,1][:set_title]("Float64",loc="left")
axs[1,2][:set_title]("Float16",loc="left")
axs[2,1][:set_title]("Posit(16,1)",loc="left")
axs[2,2][:set_title]("Int16",loc="left")



tight_layout()
savefig("/home/kloewer/julia/lorenz_posit/trajectory/figs/lorenz_attractor.png",dpi=300)
close(fig)
