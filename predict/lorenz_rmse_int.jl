using SigmoidNumbers
using JLD

cd("/home/kloewer/julia/lorenz_posit/predict/")
include("lorenz_integration_posits.jl")

# create initial conditions
Nlong = 10000
Δt = 0.01

σ,ρ,β = 10.,28.,8./3.
s = 10000.
nbit = 32

# start somewhere
XYZ0 = time_integration(Nlong,Float64,[.5,.5,15.],σ,ρ,β,1.,Δt)

#
N = 3000
time = 0:Δt:(N*Δt)
M = 1000    # number of independent forecasts, one from M different start dates

# preallocate
RMSE_int = zeros(M,N+1)
xyz_int = zeros(3,N+1)

if nbit == 16
    rint = Int16
elseif nbit == 32
    rint = Int32
elseif nbit == 64
    rint = Int64
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

    # single/half precision integers
    xyz_int[:] = time_integration_int(N,rint,xyz0,σ,ρ,β,s,Δt)
    RMSE_int[j,:] = sqrt.(sum((xyz_int-xyz_true).^2,1))

end

save("data/RMSE_$(nbit)bit_int_s$s.jld","RMSE_I",RMSE_int)
println("Data stored for s=$s,n=$nbit")
