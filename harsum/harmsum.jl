using SigmoidNumbers
using PyPlot

function harmsum(Numtype::DataType,imax::Int)
    s = Numtype(1.0)
    for i ∈ 2:imax
        s_old = s
        s += Numtype(1/i)
        if s == s_old
            println((Float64(s),i))
            break
        end
    end
end


N = 100000
x = randn(Float64,N)
cfloat = fill(0,N)
cposit = fill(0,N)
l = 0:2^16-1

for i ∈ 1:N
    cfloat[i] = Int(reinterpret(UInt16,Float16(x[i])))
    cposit[i] = Int(parse(UInt16,"0b"*bits(Posit{16,1}(x[i]))))
end

hfloat = fit(Histogram,cfloat,0:2^16)
hposit = fit(Histogram,cposit,0:2^16)
wfloat = hfloat.weights
wposit = hposit.weights
