using TimeIt
using SigmoidNumbers

i = [1,2,3,4,5]
j = [1,2,3,4,5]
v1 = [1,2,4,8,16]

P = Posit{8,0}

v1 = P.(v1)
x = zeros(v1)
#v1 = convert(Array{BigFloat},v1)
S = sparse(i,j,v1)

# compile
A_mul_B!(x,S,v1)

#=
# measure time
function sparsemat(S,v1)
    #tic()
    for n = 1:100000
        x = S*v1
    end
    #toc()
end

sparsemat(S,v1)
=#
