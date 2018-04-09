function representable_ints(nbits)
    return 1:(2^(nbits-1)-1)
end

function representable_floats(nbits,ebits)
    #= returns an array of all representable positive floats exclusive 0 and ComplexInfinity.
    nbits is the total number of bits (sign,exponent,significand)
    ebits is the number of exponent bits
    =#

    # derived constants
    bias = 2^(ebits-1) - 1
    sbits = nbits-ebits-1     # number of significand bits (minus one for sign bit)

    #=
    setprecision(256)
    repr_floats = BigFloat[0.]           # represented as 000...0000
    bigzero = big"0.0"
    one = big"1.0"
    two = big"2.0"
    =#

    repr_floats = []           # represented as 000...0000

    # subnormal numbers (exponent is 0)
    for f = 1:2^sbits-1
        append!(repr_floats,2.0^(1-bias)*(f/2^sbits))
    end

    # all other numbers
    for e = 1:2^ebits-2     # the exponent cannot be 11...111 this represents NaN
        for f = 0:2^sbits-1
            append!(repr_floats,2.0^(e-bias)*(1.+f/2^sbits))
        end
    end

    return Float64.(repr_floats)
end

function representable_flaots_expo(nbits,ebits)
    bias = 2^(ebits-1) - 1
    sbits = sbits = nbits-ebits-1
    smallsubnormal = 2.0^(1-bias)/2^sbits
    return vcat(smallsubnormal,2.0.^((1:2^ebits-2)-bias))
end

function representable_posits(nbits,ebits)
    #= returns an array of all representable positive posits exclusive 0 and ComplexInfinity.
    nbits is the total number of bits (sign,regime,exponent,fraction)
    ebits is the number of exponent bits
    =#

    N = 2^(nbits-1)-2   # number of posits in the exclusive range (0,ComplexInfinity)

    P = Posit{nbits,ebits}

    # minpos, smallest representable positive number
    p0 = next(P(0.))
    posit_list = [p0]

    for i = 1:N
        append!(posit_list,next(posit_list[end]))
    end

    return Float64.(posit_list)
end

function representable_posits_regime(nbits,ebits)
    useed = 2^(2^ebits)
    maxexp = nbits - 2
    return Float64(useed).^(-maxexp:maxexp)
end

function representable_posits_regimeexpo(nbits,ebits)
    eseed = 2^ebits
    maxexp = nbits-2
    if ebits == 0
        return 2.0.^(-maxexp:maxexp)
    else
        exprange = -maxexp:maxexp
        erange = eseed*exprange[1:ebits]
        exps = vcat(erange,eseed*exprange[ebits+1]:eseed*exprange[end-ebits],eseed*exprange[end-(ebits-1):end])
        return 2.0.^exps
    end
end

#=
function representable_posits_regimeexpo(nbits,ebits)
    N = 2^(nbits-1)-2   # number of posits in the exclusive range (0,ComplexInfinity)

    P = Posit{nbits,ebits}

    # minpos, smallest representable positive number
    pnext = next(P(0.))
    posit_list = [pnext]

    for i = 1:N
        pnext = next(pnext)
        pbits = split(bits(pnext," "))
        if length(pbits) < 4
            append!(posit_list,pnext)
        elseif !contains(pbits[4],"1")  # fraction bits are all zero
            append!(posit_list,pnext)
        end
    end
    return Float64.(posit_list)
end
=#

function wc_dec_acc_posit(nbits,ebits)
    # WORST CASE DECIMAL ACCURACY FOR FLOATS (ONLY REGIME BIT SUBSET)
    P = Posit{nbits,ebits}

    # list of posits (subset of all posits: only regime bits), and the next representable posit
    plist = representable_posits_regime(nbits,ebits)
    plist_next = Float64.(next.(P.(plist[1:end-1])))

    # arithmetic mean between representable numbers
    p_am = (plist[1:end-1]+plist_next)/2.

    # worst case decimal accuracy
    p_wda = -log10.(abs.(log10.(p_am./plist_next)))

    # extend first and last point, taking no overflow/underflow into account
    p0 = plist[1]/16    # something much smaller than minpos
    pinf = plist[end]*16    # something much bigger than maxpos

    # worst-case decimal accuracy for these extreme values
    p_wda_0 = -log10.(abs.(log10.(p0/plist[2])))
    p_wda_inf = -log10.(abs.(log10.(pinf/plist[end])))

    # worst-case decimal accuracy of interpolated on minpos/maxpos
    p_wda_minpos = p_wda_0 + log10(plist[1]/p0)/log10(p0/p_am[1])*(p_wda_0-p_wda[1])
    p_wda_maxpos = p_wda_inf + log10(plist[end]/pinf)/log10(pinf/p_am[end])*(p_wda_inf-p_wda[end])

    p_wda_r = vcat(p_wda_0,p_wda_minpos,p_wda,p_wda_maxpos,p_wda_inf)
    p_amr = vcat(p0,plist[1],p_am,plist[end],pinf)

    return p_am,p_wda
end

function wc_dec_acc_float(nbits,ebits)

    # list of representable floats (subset of all floats, only exponent bits flip over)
    # flist_next is the next bigger float (significand bits flip one over)
    flist = representable_flaots_expo(nbits,ebits)
    flist_next = flist[1:end-1]*(1.+1./2^(nbits-ebits-1))
    flist_next[1] = flist[1]*2    # nextfloat of smallsubnormal is *2

    # arithmetic mean between representable numbers
    f_am = (flist[1:end-1]+flist_next)/2.

    # worst case decimal accuracy
    f_wda = -log10.(abs.(log10.(f_am./flist_next)))

    # extend with zeros due to overflow
    f_wda = vcat(0,f_wda,0)
    f_am = vcat(flist[1],f_am,flist[end])

    return f_am,f_wda
end

function wc_dec_acc_int(nbits)

    i_am = [1.,2.,2^(nbits-1)-1,2^(nbits-1)-1]
    i_wda = [0.,-log10(abs(log10(2))),-log10(abs(log10((2^(nbits-1)-1)/(2^(nbits-1)-2)))),0]

    return i_am,i_wda
end
