using SigmoidNumbers

nbits = 8
N = 2^nbits
P = Posit{nbits,0}
neg_maxreal = next(P(Inf))

# preallocate
C = Array{Float64}(N,N)

xi = neg_maxreal
xj = neg_maxreal

for i=0:N-1
   for j=0:N-1
      xir = Float64(xi)
      xjr = Float64(xj)

      C[i,j] = 
