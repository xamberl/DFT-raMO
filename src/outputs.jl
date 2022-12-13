function write_to_XSF(
    super::AtomList{3},
    coeff::Matrix{ComplexF64},
    kpoint_repeating::Vector{Bool},
    kptlist::Vector{Vector{Float64}},
    psi_up::Matrix{ComplexF64},
    G::Vector{Vector{Int64}},
    energy::Vector{Float64}
    )
    printing_size_limit = 180000000
    ratios = lengths(super.basis)/maximum(lengths(super.basis))
    ratiofactor = (printing_size_limit/(ratios[1]*ratios[2]*ratios[3]*size(coeff)[1]))^(1/3) #size(coeff)[1] = num plane waves
    ngfftsize = Int.(floor.(ratios*ratiofactor))
    NX = zeros(ngfftsize[1],ngfftsize[2],ngfftsize[3])
    NY = zeros(ngfftsize[1],ngfftsize[2],ngfftsize[3])
    NZ = zeros(ngfftsize[1],ngfftsize[2],ngfftsize[3])
    for i in 1:ngfftsize[1]
        for j in 1:ngfftsize[2]
            for k in 1:ngfftsize[3]
                NX[i,j,k] = lengths(super.basis)[1]*(i-1)/(ngfftsize[1]-1)
                NY[i,j,k] = lengths(super.basis)[2]*(j-1)/(ngfftsize[2]-1)
                NZ[i,j,k] = lengths(super.basis)[3]*(k-1)/(ngfftsize[3]-1)
            end
        end
    end
    NX = vec(NX); NY = vec(NY); NZ = vec(NZ)
    reduced_kpt = unique(kptlist)
    kpoint_zero = vcat(findall(!iszero, kpoint_repeating),length(kpoint_repeating)+1)
    target_overlap = []
    for k in 1:length(kpoint_zero)-1
        # Finds overlap of plane waves at unique kpoints
        range = vcat(kpoint_zero[k]:1:(kpoint_zero[k+1])-1)
        push!(target_overlap, coeff[:,range]*psi_up[range,:])
    end
    @info "Finished generating target_overlap"
    isosurf = zeros(ngfftsize[1]*ngfftsize[2]*ngfftsize[3])
    for i in eachindex(reduced_kpt)
    #i = 1
    pw = Matrix{ComplexF16}(undef,length(G),length(NX))
    for j in eachindex(G)
        # (kx+Gx)*NX + (ky+Gy)*NY + (kz+Gz)*NZ
        pw[j,:] = sum((G[j].+reduced_kpt[i]).*[NX, NY, NZ])
        pw[j,:] = @. exp(pi*2*im*pw[j,:])
        @info j
    end
    isosurf += pw'*target_overlap[i]
    @info string("K-point number ", i, " printed.")
    #return isosurf
    end
    isosurf1 = isosurf
    
    isosurf = reshape(isosurf,ngfftsize[1],ngfftsize[2],ngfftsize[3]) 
    
    # Wave functions read in are complex. We rotate them back to the real axis.
    # We also deal with degenerate states
    degen_tol = 0.001
    energy = vcat(energy,10000)
    num_print = 1 # This variable is in case we have more than one to print. This will need to be an input argument. For now, set to one.
    num_degen = 0
    temp = 0
    for i in 1:num_print
        while temp == num_degen
            if abs(energy[i+temp+1]-energy[i] < degen_tol*abs(energy[i]))
                num_degen += 1
                temp += 1
            else
                num_degen += 1
            end
        end
        @info num_degen
        temp_matrix = reshape(isosurf[:,:,:,i:i+num_degen-1],ngfftsize[1]*ngfftsize[2]*ngfftsize[3],num_degen)
        # Pick out 1000 random vectors from the matrix
        randos = randperm(ngfftsize[1]*ngfftsize[2]*ngfftsize[3])[1:1000]
        @info randos
        # Separate into real and imaginary components
        realimag = hcat(real(temp_matrix[randos,:]), imag(temp_matrix[randos,:]))
        (eig_energy, eig_vect) = eigen(realimag'*realimag)
        vect = eig_vect[:,1:num_degen]
        @info vect
        vect = complex.(vect[num_degen+1:2*num_degen,:],vect[1:num_degen,:])
        @info vect
        temp_matrix = temp_matrix*vect
        #isosurf2 = Array{ComplexF64}(undef,)
        global isosurf[:,:,:,i:i+num_degen-1] = reshape(temp_matrix, ngfftsize[1],ngfftsize[2],ngfftsize[3], num_degen)
    end
    return isosurf1, isosurf
end  

#reshape(isosurf,ngfftsize[1],ngfftsize[2],ngfftsize[3]) 

#==
# Wave functions read in are complex. We rotate them back to the real axis.
# We also deal with degenerate states
degen_tol = 0.001
vcat(energy,10000)
num_print = 1 # This variable is in case we have more than one to print. This will need to be an input argument. For now, set to one.
num_degen = 0
temp = 0
for i in 1:num_print
    while temp == num_degen
        if abs(energy[i+temp+1]-energy[i] < degen_tol*abs(energy[i]))
            num_degen += 1
            temp += 1
        else
            num_degen += 1
        end
    end
end
temp_matrix = reshape(isosurf[:,:,j:j+num_degen-1],ngfftsize[1]*ngfftsize[2]*ngfftsize[3],num_degen)
# Pick out 1000 random vectors from the matrix
randos = randperm(ngfftsize[1]*ngfftsize[2]*ngfftsize[3])[1:1000]
# Separate into real and imaginary components
realimag = hcat(real(temp_matrix[randos,:]), imag(temp_matrix[randos,:]))
(eig_energy, eig_vect) = eigen(realimag'*realimag)
vect = eig_vect[1:num_degen]
vect = complex.(vect[num_degen+1:2*num_degen,:],vect[1:num_degen,:])
temp_matrix = temp_matrix*vect
isosurf[:,:,:,j:j+num_degen-1] = reshape(temp_matrix, ngfftsize[1],ngfftsize[2],ngfftsize[3], num_degen)
return isosurf==#
