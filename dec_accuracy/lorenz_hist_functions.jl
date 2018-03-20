function lorenz_hist(xyz,s,bins)

    # for readability
    x = xyz[1,:]
    y = xyz[2,:]
    z = xyz[3,:]

    # preallocate
    H = zeros(length(bins)-1)

    # histogram of prognostic variables
    H = hist_add(H,x)
    H = hist_add(H,y)
    H = hist_add(H,z)

    # add to histogram all substeps of the RHS
    H = hist_add(H,y-x)             # first line
    H = hist_add(H,σ*(y-x))
    H = hist_add(H,s*z)             # second line
    H = hist_add(H,ρ-s*z)             # second line
    H = hist_add(H,x.*(ρ-s*z))
    H = hist_add(H,x.*(ρ-s*z)-y)
    H = hist_add(H,x.*y)             # third line
    H = hist_add(H,s*x.*y)
    H = hist_add(H,β*z)
    H = hist_add(H,s*x.*y-β*z)

    return H
end

function hist_add(H0,vec)
    Htemp,bin_edges = np.histogram(abs.(vec),bins)
    H0 += Htemp
    return H0
end
