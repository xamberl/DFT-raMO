function write_to_XSF(
    isosurf::Array{ComplexF64, 3},
    xtal::PeriodicAtomList{3},
    filename::String
    )
    xcrystal = Crystal(xtal, 1, SVector{3, Float64}(0,0,0)) # space group number is 1 for now; origin is [0,0,0]
    xdata = RealDataGrid(real(isosurf), xtal.basis)
    xcrystalwithdatasets = CrystalWithDatasets(xcrystal, Dict("DENSITY"=>xdata))
    f = open(filename,"w")
    writeXSF(f, xcrystalwithdatasets, periodic=false)
    close(f)
end

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
    reduced_kpt = unique(occ_states.kpt)
    target_overlap = Array{ComplexF32}(undef,size(occ_states.coeff)[1], length(reduced_kpt))
    for k in eachindex(reduced_kpt)
        m = findfirst(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        n = findlast(x->x == reduced_kpt[k], occ_states.kpt[1,:])
        target_overlap[:,k] = occ_states.coeff[:,m:n]*psi_up[m:n]
    end
    #@info "Finished generating target_overlap"
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
        #@info string("K-point number ", i, " printed.")
    end
    realimag = hcat(real(isosurf), imag(isosurf))
    (evalue, evec) = eigen(realimag'*realimag)
    isosurf *= complex(evec[2], evec[1])
    isosurf = reverse(reshape(isosurf,ngfftsize[1],ngfftsize[2],ngfftsize[3]))
    return isosurf
end

function Psphere(
    datagrid::RealDataGrid,
    origin::Vector{Float64},
    rsphere::Float64
    )

    # Volume of cell in BOHR
    vol = volume(datagrid)

    # Volume of voxel, length of voxel dimensions (angle not stored)
    vox_vol = vol/length(datagrid)

    # Find unnormalized electron count
    total_eden = sum(datagrid.data.^2)*vox_vol

    # Calculate unnormalized electron count in defined radius
    # Calculate bounds of box of sphere
    lowerbound = (origin.-rsphere)*Electrum.ANG2BOHR
    upperbound = (origin.+rsphere)*Electrum.ANG2BOHR
    # Calculate indices of box in the voxel grid
    lower_index = round.(Int, (lowerbound'/datagrid.basis)'.*size(datagrid))
    upper_index = round.(Int, (upperbound'/datagrid.basis)'.*size(datagrid))
    sphere_sum = 0
    count = 0
    for x in lower_index[1]:upper_index[1]
        for y in lower_index[2]:upper_index[2]
            for z in lower_index[3]:upper_index[3]
                # Transform grid index into Cartesian (Angstrom)
                grid_cart = (([x, y, z]./size(datagrid))'*datagrid.basis)*Electrum.BOHR2ANG
                # check if grid_cart is within sphere
                if ((grid_cart[1]-origin[1])^2+(grid_cart[2]-origin[2])^2+(grid_cart[3]-origin[3])^2)^.5 <= rsphere
                    # Check if point is out of bounds. If so, use adjacent cell's value.
                    (i1, i2, i3) = (x, y, z)
                    if x < 1
                        i1 = size(datagrid)[1]+x
                    elseif x > size(datagrid)[1]
                        i1 = x-size(datagrid)[1]
                    end
                    if y < 1
                        i2 = size(datagrid)[2]+y
                    elseif y > size(datagrid)[2]
                        i2 = y-size(datagrid)[2]
                    end
                    if z < 1
                        i3 = size(datagrid)[3]+z
                    elseif z > size(datagrid)[3]
                        i3 = z-size(datagrid)[3]
                    end
                    sphere_sum += datagrid.data[i1,i2,i3]^2
                    count += 1
                end
            end
        end
    end
    # Calculate partial electron density
    sphere_sum = sphere_sum*vox_vol
    psphere = sphere_sum/total_eden
    return(sphere_sum, total_eden, psphere)
end