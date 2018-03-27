using SigmoidNumbers
using JLD

cd("/home/kloewer/julia/lorenz_posit/predict/")
include("lorenz_integration_posits.jl")

# create initial conditions
Nlong = 10000
Δt = 0.01

σ,ρ,β = 10.,28.,8./3.
s = 1000.
nbit = 16

# start somewhere
XYZ0 = time_integration(Nlong,Float64,[.5,.5,15.],σ,ρ,β,1.,Δt)

#
N = 1000
time = 0:Δt:(N*Δt)
M = 100    # number of independent forecasts, one from M different start dates

es = [1,2,3]    # number of exponent bits for posits
nes = length(es)

# preallocate
RMSE_float = Array{Float64}(M,N+1)
RMSE_posit = Array{Float64}(nes,M,N+1)

if nbit == 16
    rpf = Float16   # float environment (rpf: reduced precision float)
elseif nbit == 32
    rpf = Float32
end

println("Computing RMSE for s=$s")
for j = 1:M

    if ((j+1)/M*100 % 10) < (j/M*100 % 10)
        progress = Int(round((j+1)/M*100))
        print("$progress%,")
    end

    # pick a random start
    randi = rand(Int(Nlong/2):Nlong)
    xyz0 = XYZ0[:,randi]

    # truth without scaling
    xyz_true = time_integration(N,Float64,xyz0,σ,ρ,β,1.,Δt)

    # single/half precision
    xyz_rpf = time_integration(N,rpf,xyz0,σ,ρ,β,s,Δt)
    RMSE_float[j,:] = sqrt.(sum((xyz_rpf-xyz_true).^2,1))

    for ip = 1:nes
        # posit
        P = Posit{nbit,es[ip]}
        xyz_p = time_integration(N,P,xyz0,σ,ρ,β,s,Δt)

        # calculate RMSE
        RMSE_posit[ip,j,:] = sqrt.(sum((xyz_p-xyz_true).^2,1))
    end
end

save("data/RMSE_$(nbit)bit_s$s.jld","RMSE_F",RMSE_float,"RMSE_P",RMSE_posit)
println("Data stored for s=$s,n=$nbit")
