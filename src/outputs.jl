function write_to_XSF(
    print_limit::Int,
    xtal::PeriodicAtomList{3},
    kpt::Vector{Int},
    occ_states,
    psi_up::Matrix{ComplexF64},
    filename::String
    )
    print_limit = 2000000 #180000000
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
    reduced_kpt = unique(occ_states.kpt)
    target_overlap = Array{ComplexF32}(undef,size(occ_states.coeff)[1], length(reduced_kpt))
    for k in eachindex(reduced_kpt)
        m = findfirst(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        n = findlast(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        target_overlap[:,k] = occ_states.coeff[:,m:n]*psi_up[m:n]
    end
    @info "Finished generating target_overlap"
    isosurf = zeros(ngfftsize[1]*ngfftsize[2]*ngfftsize[3])
    for i in eachindex(reduced_kpt)
        #Gx = occ_states.G[:,1]
        pw = Matrix{ComplexF32}(undef,size(occ_states.coeff)[1],length(NX))
        for j in eachindex(occ_states.G[:,1])
            # (kx+Gx)*NX + (ky+Gy)*NY + (kz+Gz)*NZ
            GX = (occ_states.G[j,1][1].+reduced_kpt[i][1])*NX
            GY = (occ_states.G[j,1][2].+reduced_kpt[i][2])*NY
            GZ = (occ_states.G[j,1][3].+reduced_kpt[i][3])*NZ
            G = GX+GY+GZ
            pw[j,:] = @. exp(pi*2*im*G)
        end
        isosurf += pw'*target_overlap[:,i]
        @info string("K-point number ", i, " printed.")
    end
    realimag = hcat(real(isosurf), imag(isosurf))
    (evalue, evec) = eigen(realimag'*realimag)
    isosurf *= complex(evec[2], evec[1])
    isosurf = reverse(reshape(isosurf,ngfftsize[1],ngfftsize[2],ngfftsize[3]), dims=3)
    xcrystal = Crystal(xtal, 1, SVector{3, Float64}(0,0,0)) # space group number is 1 for now; origin is [0,0,0]
    xdata = RealDataGrid(real(isosurf), xtal.basis)
    xcrystalwithdatasets = CrystalWithDatasets(xcrystal, Dict("DENSITY"=>xdata))
    f = open(filename,"w")
    writeXSF(f, xcrystalwithdatasets, periodic=false)
    close(f)
end