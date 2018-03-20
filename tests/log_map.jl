using PyPlot

# test COMMENT

N = 100
x = zeros(N+1)
x2 = Array{UInt8}(N+1)

s = 10

r = 3.57
ρ = UInt8(round(r*s))

x[1] = .12
x2[1] = UInt8(round(x[1]*s))

rhs(x) = r*x*(1-x)
rhs2(x) = div(ρ*x*(s-x),s^2)

for i=1:N
    x[i+1] = rhs(x[i])
    x2[i+1] = rhs2(x2[i])
end

plot(x)
plot(x2)
