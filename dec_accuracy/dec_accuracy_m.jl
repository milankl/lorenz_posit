using SigmoidNumbers
using PyPlot

# next/prev for posits is nextfloat/prevfloat for Float16,Float32,Float64

nbits = 16
N = 2^(nbits-1)-1

allposits = Array{Float64,1}(N)

fig,(ax1,ax2,ax3) = subplots(3,1,sharex=true)

P = Posit{nbits,0}
positzero = P(0)
p = positzero    # start from the smallest posit

for i=1:N
    p = next(p)
    allposits[i] = Float64(p)
end

acc = 2./(allposits[2:end]-allposits[1:end-1])./abs.(allposits[1:end-1])
ax1[:loglog](allposits[1:end-1],acc)
for j=1:N
    ax2[:axvline](allposits[j],lw=0.2)
end

P = Posit{nbits,1}
positzero = P(0)
p = positzero    # start from the smallest posit

for i=1:N
    p = next(p)
    allposits[i] = Float64(p)
end

acc = 2./(allposits[2:end]-allposits[1:end-1])./abs.(allposits[1:end-1])
ax1[:loglog](allposits[1:end-1],acc)
for j=1:N
    ax3[:axvline](allposits[j],lw=0.2,color="C1")
end

tight_layout()
