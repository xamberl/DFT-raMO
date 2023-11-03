"""
    function write_to_XSF(
        isosurf::AbstractArray{<:Number, 3},
        xtal::AbstractAtomList{3},
        filename::AbstractString
    )

Writes an XSF using the real components of `isosurf`.
"""
function write_to_XSF(
    isosurf::AbstractArray{<:Number, 3},
    xtal::AbstractAtomList{3},
    filename::AbstractString
)
    xcrystal = Crystal(xtal, 1, SVector{3, Float64}(0,0,0)) # space group number is 1 for now; origin is [0,0,0]
    xdata = RealDataGrid(real(isosurf), basis(xtal))
    xcrystalwithdatasets = CrystalWithDatasets(xcrystal, Dict("DENSITY"=>xdata))
    open(io -> writeXSF(io, xcrystalwithdatasets, periodic=true), filename, write=true)
end

"""
    raMO_to_density(occ_states, psi, kpt, grange)

Converts the raMO (psi) into an electron density grid using FFT.
"""
function raMO_to_density(
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
    #return abs2.(isosurf)
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

Saves files to disk:
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
    open(string(run_name, "_", num_raMO, "_", num_electrons_left, ".chkpt"), "w") do io
        sz = Int64.(size(psi_previous))
        write(io, sz[1], sz[2], sz[3])
        write(io, psi_previous)
    end
    # same with the raMO function itself
    open(string(run_name, "_", num_raMO, "_", num_electrons_left, ".raMO"), "w") do io
        write(io, psi_up)
    end
end
