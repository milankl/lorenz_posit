using FastSigmoid
using PyCall

setprecision(8)

@pyimport numpy.random as random

N = Int(1e5)
r = random.randn(N)
#Posits = Posit{16,1}

#r = convert(Array{Posits,1},r)
#r = convert(Array{BigFloat,1},r)

tic()
for n=1:100
    r + r
end
toc()
