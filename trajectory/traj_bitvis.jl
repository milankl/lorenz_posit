using SigmoidNumbers
using PyPlot
using PyCall
@pyimport matplotlib.patches as pat

include("lorenz_integration_posits.jl")

Nt0 = 10000      # length of spin-up
Δt = 0.01    # time step

σ = 10.      # Lorenz 63 coefficients
β = 8./3.
ρ = 28.
s = 0.1

# start somewhere
xyz = [0.6,0.5,15.]
nbits = 32
npbits = 1
P = Posit{nbits,npbits}

# do integration
xyz0 = time_integration_opt(Nt0,Float64,xyz,σ,ρ,β,s,Δt)  # spin-up in full precision

istart,iend = 300,420
Nn = iend-istart+1
nbits = 32

#convert to bits
x = (xyz0[1,istart:iend]+1)*2.1

# floats
Bf = zeros(Int,nbits,Nn)
for i = 1:Nn
    bitpattern = bits(Float32(x[i]))
    ic = 0
    for c in bitpattern
        ic += 1
        if c == '1'
            Bf[ic,i] = 1
        end
    end
end

#posits
Bp = zeros(Bf)
for i = 1:Nn
    bitpattern = bits(P(x[i]))
    ic = 0
    for c in bitpattern
        ic += 1
        if c == '1'
            Bp[ic,i] = 1
        end
    end
end

# PLOTTING
fig,(ax1,ax2,ax3) = subplots(3,1,figsize=(8,7))

ax1[:plot](x)
ax2[:spy](Bf)
ax3[:spy](Bp)

ax1[:grid]()

ax1[:set_ylabel]("x")
ax1[:set_xlim](0,Nn)
ax1[:set_xticks]([])
ax3[:set_xlabel]("time")
ax1[:set_yticks]([-32,-16,-8,-4,0,4,8,16,32])

ax2[:set_xticks]([])
ax2[:set_yticks]([0.5,8.5])
ax2[:set_yticklabels]([])
ax2[:set_ylabel]("Float32")
ax2[:add_patch](pat.Rectangle([-.5,.5],200,8,color="C0",alpha=.5))

ax3[:set_xticks]([])
ax3[:set_yticks]([0.5])
ax3[:set_yticklabels]([])
ax3[:set_ylabel]("Posit(32,1)")

# add regime bit information
switch01(x::Int) = mod(x+1,2)
rlength(x::Array) = findfirst(x[2:end],switch01(x[2]))

for i = 1:Nn
    rl = rlength(Bp[:,i])
    ax3[:add_patch](pat.Rectangle([i-1.5,.5],1,rl,color="C1",alpha=.5,ls="None"))
    ax3[:add_patch](pat.Rectangle([i-1.5,.5+rl],1,1,color="C0",alpha=.5))
end

tight_layout()
savefig("figs/traj_bitvis.pdf")
close(fig)
