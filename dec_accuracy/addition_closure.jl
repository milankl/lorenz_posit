using SigmoidNumbers
using PyPlot

nbits = 8
N = 2^nbits
P = Posit{nbits,1}
neg_maxreal = next(P(Inf))

# preallocate
C = Array{Float64}(N,N)

xi = neg_maxreal


for i=1:N
    xj = neg_maxreal
    xir = Float64(xi)
    for j=1:N

        xjr = Float64(xj)

        try
            C[i,j] = abs(log10(Float64(xi*xj)/(xir*xjr)))
        catch
            C[i,j] = NaN
        end

        xj = next(xj)
    end
    xi = next(xi)
end

pcolor(C,vmin=0,vmax=.1)
