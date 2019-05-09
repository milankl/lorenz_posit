using SoftPosit
using PyPlot
using StatsBase

N = 100000
x = 3*randn(Float64,N)
cfloat = fill(0,N)
cposit1 = fill(0,N)
cposit2 = fill(0,N)
l = 0:2^16-1
#
# function posit16_2int(x::Real,P16)
#     return Int(parse(UInt16,bitstring(P16(x)),base=2))
# end
#
# function float16_2int(x::Real)
#     return Int(reinterpret(UInt16,Float16(x)))
# end

function real2int(x::Real,T::DataType)
    return Int(reinterpret(UInt16,T(x)))
end

for i ∈ 1:N
    cfloat[i] = real2int(x[i],Float16)
    cposit1[i] = real2int(x[i],Posit16)
    #cposit2[i] = real2int(x[i],Posit16_2)
end

hfloat = fit(Histogram,cfloat,0:2^16)
hposit1 = fit(Histogram,cposit1,0:2^16)
# hposit2 = fit(Histogram,cposit2,0:2^16)

wfloat = hfloat.weights
wposit1 = hposit1.weights
# wposit2 = hposit2.weights

## orientation ticks
tik = [0,1,-1,4,-4,1/4,-1/4,Inf]

otickf = [real2int(t,Float16) for t in tik]
otick1 = [real2int(t,Posit16) for t in tik]
# otick2 = [real2int(t,Posit16_2) for t in tik]

## visualising nans
nan_pos0 = Int(parse(UInt16,"0111110000000001",base=2))
nan_pos1 = Int(parse(UInt16,"0111111111111111",base=2))

nan_neg0 = Int(parse(UInt16,"1111110000000001",base=2))
nan_neg1 = Int(parse(UInt16,"1111111111111111",base=2))

##
fig,(ax1,ax2) = subplots(2,1,sharey=true)

ax1[:plot](l,wfloat,"k")
#ax2[:plot](l,wposit0,label="Posit(16,0)")
ax2[:plot](l,wposit1,label="Posit(16,1)")
#ax2[:plot](l,wposit2,label="Posit(16,2)")

ax1[:fill_between]([nan_pos0,nan_pos1],0,40,alpha=0.2,color="red")
ax1[:fill_between]([nan_neg0,nan_neg1],0,40,alpha=0.2,color="red")

ax1[:set_title]("Float16",loc="left")
ax2[:set_title]("Posit16",loc="left")

ax1[:set_ylabel]("N")
ax2[:set_ylabel]("N")

ax1[:set_ylim](0,40)
ax1[:set_xlim](0,2^16-1)
ax2[:set_xlim](0,2^16-1)

ax1[:set_xticks](otickf)
ax1[:set_xticklabels](string.(tik))

ax2[:set_xticks](otick1)
ax2[:set_xticklabels](string.(tik))
