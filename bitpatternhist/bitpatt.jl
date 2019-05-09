using SoftPosit
using PyPlot
using StatsBase
using Printf
using BFloat16s

N = 10000000
x = 3*randn(Float64,N)
cfloat = fill(0,N)
cbfloat = fill(0,N)
cint = fill(0,N)
cposit1 = fill(0,N)
cposit2 = fill(0,N)
l = 0:2^16-1

function real2int(x::Real,T)
    return Int(reinterpret(UInt16,T(x)))
end

function real2int_positX2(x::Real,T)
    return Int(reinterpret(UInt32,T(x)) >> 16)
end

for i ∈ 1:N
    cfloat[i] = real2int(x[i],Float16)
    cbfloat[i] = real2int(Float32(x[i]),BFloat16)
    cposit1[i] = real2int(x[i],Posit16)
    cposit2[i] = real2int_positX2(x[i],Posit16_2)
    cint[i] = real2int(x[i],Int16∘round)
end

hfloat = fit(Histogram,cfloat,0:2^16)
hbfloat = fit(Histogram,cbfloat,0:2^16)
hposit1 = fit(Histogram,cposit1,0:2^16)
hposit2 = fit(Histogram,cposit2,0:2^16)
hint = fit(Histogram,cint,0:2^16)

wfloat = hfloat.weights
wbfloat = hbfloat.weights
wposit1 = hposit1.weights
wposit2 = hposit2.weights
wint = hint.weights

shannon_float = entropy(wfloat/N,2)
shannon_bfloat = entropy(wbfloat/N,2)
shannon_posit1 = entropy(wposit1/N,2)
shannon_posit2 = entropy(wposit2/N,2)
shannon_int = entropy(wint/N,2)

## orientation ticks
tik = [0,1,-1,4,-4,1/4,-1/4,Inf]
stik = ['0','1',"-1",'4',"-4",'¼',"-¼","∞"]

otickf = [real2int(t,Float16) for t in tik]
otickb = [real2int(Float32(t),BFloat16) for t in tik]
otickp1 = [real2int(t,Posit16) for t in tik]
otickp2 = [real2int_positX2(t,Posit16_2) for t in tik]

## visualising nans
nan_pos0 = Int(parse(UInt16,"0111110000000001",base=2))
nan_pos1 = Int(parse(UInt16,"0111111111111111",base=2))

nan_neg0 = Int(parse(UInt16,"1111110000000001",base=2))
nan_neg1 = Int(parse(UInt16,"1111111111111111",base=2))

nanb_pos0 = Int(parse(UInt16,"0111111110000001",base=2))
nanb_pos1 = Int(parse(UInt16,"0111111111111111",base=2))

nanb_neg0 = Int(parse(UInt16,"1111111110000001",base=2))
nanb_neg1 = Int(parse(UInt16,"1111111111111111",base=2))


##
fig,(ax1,ax2,ax3,ax4) = subplots(4,1)

ax1.plot(l,100*wfloat/N,"k")
ax2.plot(l,100*wbfloat/N,"k")
ax3.plot(l,100*wposit1/N,label="Posit(16,1)")
ax4.plot(l,100*wposit2/N,label="Posit(16,1)")

ax1.fill_between([nan_pos0,nan_pos1],0,ax1.get_ylim()[2],alpha=0.2,color="red")
ax1.fill_between([nan_neg0,nan_neg1],0,ax1.get_ylim()[2],alpha=0.2,color="red")

ax2.fill_between([nanb_pos0,nanb_pos1],0,ax1.get_ylim()[2],alpha=0.2,color="red")
ax2.fill_between([nanb_neg0,nanb_neg1],0,ax1.get_ylim()[2],alpha=0.2,color="red")

ax1.set_title("Float16, H="*@sprintf("%.2f",shannon_float),loc="left")
ax2.set_title("BFloat16, H="*@sprintf("%.2f",shannon_bfloat),loc="left")
ax3.set_title("Posit16,1 H="*@sprintf("%.2f",shannon_posit1),loc="left")
ax4.set_title("Posit16,2 H="*@sprintf("%.2f",shannon_posit2),loc="left")

ax1.set_ylabel("%")
ax2.set_ylabel("%")
ax3.set_ylabel("%")
ax4.set_ylabel("%")

ax1.set_xlim(0,2^16-1)
ax2.set_xlim(0,2^16-1)
ax3.set_xlim(0,2^16-1)
ax4.set_xlim(0,2^16-1)

ax2.set_ylim(ax1.get_ylim())
ax3.set_ylim(ax1.get_ylim())
ax4.set_ylim(ax1.get_ylim())
ax1.set_xticks(otickf)
ax1.set_xticklabels(string.(stik))

ax2.set_xticks(otickb)
ax2.set_xticklabels(string.(stik))

ax3.set_xticks(otickp1)
ax3.set_xticklabels(string.(stik))

ax4.set_xticks(otickp2)
ax4.set_xticklabels(string.(stik))
