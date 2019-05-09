using SoftPosit

function harmsum(Numtype::DataType,imax::Int = 2000)
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

harmsum(Posit16,2000)
harmsum(Float16,2000)
