function write_to_XSF(
    isosurf::Array{ComplexF64, 3},
    xtal::PeriodicAtomList{3},
    filename::String
    )
    xcrystal = Crystal(xtal, 1, SVector{3, Float64}(0,0,0)) # space group number is 1 for now; origin is [0,0,0]
    xdata = RealDataGrid(real(isosurf), xtal.basis)
    xcrystalwithdatasets = CrystalWithDatasets(xcrystal, Dict("DENSITY"=>xdata))
    f = open(filename,"w")
    writeXSF(f, xcrystalwithdatasets, periodic=true)
    close(f)
end

"""
    psi_to_isosurf(
        print_limit::Int,
        xtal::PeriodicAtomList{3},
        kpt::Vector{Int},
        occ_states,
        psi_up::Matrix{ComplexF64}
    )
Transforms a raMO function (psi_up) to an electron density grid.
    - psi_up is a vector of complex coefficients with length = number of occupied states
    - xtal corresponds to the supercell
    - kpt corresponds to the supercell size
    - occ_states contains information about occupied states.
        occ_states.coeff = complex coefficients (dims = num_pw * num_states)
"""
function psi_to_isosurf(
    print_limit::Int,
    xtal::PeriodicAtomList{3},
    kpt::Vector{Int},
    occ_states,
    psi_up::Matrix{ComplexF64},
    )
    ratio = Vector{Float64}(undef,3)
    ratio[1] = norm(xtal.basis[1:3]); ratio[2] = norm(xtal.basis[4:6]); ratio[3] = norm(xtal.basis[7:9])
    ratio = ratio/maximum(ratio)
    ratiofactor = (print_limit/(ratio[1]*ratio[2]*ratio[3]*size(occ_states.coeff)[1]))^(1/3) #size(coeff)[1] = num plane waves
    ngfftsize = Int.(floor.(ratio*ratiofactor))
    NX = zeros(ngfftsize[1],ngfftsize[2],ngfftsize[3])
    NY = zeros(ngfftsize[1],ngfftsize[2],ngfftsize[3])
    NZ = zeros(ngfftsize[1],ngfftsize[2],ngfftsize[3])
    for i in 1:ngfftsize[1]
        for j in 1:ngfftsize[2]
            for k in 1:ngfftsize[3]
                NX[i,j,k] = kpt[1]*(i-1)/(ngfftsize[1]-1)
                NY[i,j,k] = kpt[2]*(j-1)/(ngfftsize[2]-1)
                NZ[i,j,k] = kpt[3]*(k-1)/(ngfftsize[3]-1)
            end
        end
    end
    NX = vec(NX); NY = vec(NY); NZ = vec(NZ)
    # Generate pw_raMO: express psi_up in terms of the occupied coefficients.
    # States at the same kpt are combined such that pw_raMO has (dims = num_pw * num_kpts)
    reduced_kpt = unique(occ_states.kpt)
    pw_raMO = Array{ComplexF32}(undef,size(occ_states.coeff)[1], length(reduced_kpt))
    for k in eachindex(reduced_kpt)
        m = findfirst(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        n = findlast(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        pw_raMO[:,k] = occ_states.coeff[:,m:n]*psi_up[m:n]
    end
    # Generate isosurface
    # Ψk(r) = exp(ikr)*sum(ck*exp(iGr))
    isosurf = zeros(ngfftsize[1]*ngfftsize[2]*ngfftsize[3])
    for i in eachindex(reduced_kpt)
        pw = Matrix{ComplexF32}(undef,size(occ_states.coeff)[1],length(NX))
        for j in eachindex(occ_states.G[:,1])
            # (kx+Gx)*NX + (ky+Gy)*NY + (kz+Gz)*NZ
            GX = (occ_states.G[j,1][1].+reduced_kpt[i][1])*NX
            GY = (occ_states.G[j,1][2].+reduced_kpt[i][2])*NY
            GZ = (occ_states.G[j,1][3].+reduced_kpt[i][3])*NZ
            G = GX+GY+GZ
            pw[j,:] = @. exp(pi*2*im*G)
        end
        isosurf += pw'*pw_raMO[:,i]
        #@info string("K-point number ", i, " printed.")
    end
    # This repeats the MATLAB code. But basically we create the isosurface by |Ψ|^2
    realimag = hcat(real(isosurf), imag(isosurf))
    (evalue, evec) = eigen(realimag'*realimag)
    isosurf *= complex(evec[2], evec[1])
    # Reverse for right symmetry for printing
    isosurf = reverse(reshape(isosurf,ngfftsize[1],ngfftsize[2],ngfftsize[3]))
    return isosurf
end


"""
    psi_to_isosurf2(occ_states, psi, kpt, grange)

Converts the raMO (psi) into an electron density grid using FFT.
"""
function psi_to_isosurf2(
    occ_states,
    psi,
    kpt,
    grange,
    )
    reduced_kpt = unique(occ_states.kpt)
    pw_raMO = Array{ComplexF32}(undef,size(occ_states.coeff)[1],length(reduced_kpt))
    for k in eachindex(reduced_kpt)
        m = findfirst(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        n = findlast(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        pw_raMO[:,k] = occ_states.coeff[:,m:n]*psi[m:n]
    end
    gridsize = grange
    real_gridsize = collect(gridsize).*kpt
    rgrid = [SVector(
        kpt[1]*(i-1)/(real_gridsize[1]),
        kpt[2]*(j-1)/(real_gridsize[2]),
        kpt[3]*(k-1)/(real_gridsize[3])) for i in 1:real_gridsize[1], j in 1:real_gridsize[2], k in 1:real_gridsize[3]]
    bin = FFTBins(gridsize);
    isosurf = zeros(ComplexF64, Tuple(real_gridsize));
    index = [CartesianIndex(Tuple(occ_states.G[i,1])) for i in eachindex(occ_states.G[:,1])]
    index2 = [findfirst(x->x == i, bin) for i in index] # cartesian indexing for states
    recip = zeros(ComplexF64, gridsize)
    for k in eachindex(reduced_kpt)
        for i in eachindex(pw_raMO[:,k])
            recip[index2[i]] = pw_raMO[i,k]
        end
        u = ifft(recip)
        u = repeat(u, outer = Tuple(kpt))
        kr = [dot(reduced_kpt[k],r) for r in rgrid]
        phase = [exp(2*pi*im*r) for r in kr]
        isosurf .+= u.*phase
    end
    isosurf = vec(isosurf)
    realimag = hcat(real(isosurf), imag(isosurf))
    (evalue, evec) = eigen(realimag'*realimag)
    isosurf *= complex(evec[2], evec[1]) 
    isosurf = reshape(ComplexF64.(isosurf), Tuple(real_gridsize))
    return isosurf
end

"""
output_files(
    run_name,
    num_electrons_left,
    num_raMO,
    super,
    isosurf,
    psi_previous,
    psi_up
    )
Saves files to disk.
    - raMO xsf
    - raMO coefficients
    - remainder coefficients
"""
function output_files(
    run_name,
    num_electrons_left,
    num_raMO,
    super,
    isosurf,
    psi_previous,
    psi_up
    )
    write_to_XSF(isosurf, super.atomlist, string(run_name, "_", num_raMO, "_", num_electrons_left, ".xsf"))
    # for now, write only one spin as save state
    writedlm(string(run_name, "_psi_prev_", num_raMO, "_", num_electrons_left, ".txt"), psi_previous[:,:,1])
    # same with the raMO function itself
    writedlm(string(run_name, "_psi_", num_raMO, "_", num_electrons_left, ".txt"), psi_up)
end