using LinearAlgebra

# Calculates Psphere from an xsf file, a list of coordinates, 
# and a specified radii from the coordinates
# TODO: I think you can call Xtal.readXSF() to handle most of this. -BF
function DFTraMO_overlap(
    filename::AbstractString,
    originsite::AbstractVector{<:Real},
    radius::Real
)
    # Opens xsf file
    f = open(filename,"r")
    itr = eachline(f)
    
    # Grab cell vectors
    readuntil(f,"PRIMVEC")
    iterate(itr)[1]
    primvec = Array{Float64}(undef,3,3)
    for n in 1:3
        primvec[n,:] .= parse.(Float64,split(iterate(itr)[1]))
    end
    
    # Calculate length of cell vectors in angstroms
    r_a = norm(primvec[1,:]'*primvec[1,:])
    r_b = norm(primvec[2,:]'*primvec[2,:])
    r_c = norm(primvec[3,:]'*primvec[3,:])
    
    # Calculate volume
    vol = dot(cross(primvec[1,:],primvec[2,:]),primvec[3,:])
    
    # Grab grid sizes
    readuntil(f,"DATAGRID_3D_DENSITY");
    readline(f)
    grid = parse.(Int,split(iterate(itr)[1]))
    gridsize = grid[1]*grid[2]*grid[3]
    for n in 1:4
        local ln = iterate(itr)
    end
    
    # Grab 3D data
    s = iterate(itr)[1]
    data = zeros(0)
    while !(s == "END_DATAGRID_3D")
        append!(data, parse(Float64, s))
        s = iterate(itr)[1]
    end
    data = reshape(data,(grid[1],grid[2],grid[3]))
    close(f)
    
    # Volume of voxel, length of voxel dimensions (angle not stored)
    vox_vol = vol/gridsize
    vox_a = r_a/grid[1]
    vox_b = r_b/grid[2]
    vox_c = r_c/grid[3]
    
    # Find unnormalized electron count
    total = sum(data.^2*vox_vol)
    
    # Calculate unnormalized electron count in defined radius
    origin = originsite # Fill this in from dummy atoms. Cartesian to fractional.
    r =  radius # in angstroms
    # Calculate cartesian bounds of box of sphere
    lowerbound = origin.-radius
    upperbound = origin.+radius
    # Calculate indices of box in the voxel grid
    lower_index = round.(Int, (lowerbound'/primvec)'.*grid)
    upper_index = round.(Int, (upperbound'/primvec)'.*grid)
    sphere_sum = 0
    count = 0
    for x in lower_index[1]:upper_index[1]
        for y in lower_index[2]:upper_index[2]
            for z in lower_index[3]:upper_index[3]
                # println([x, y, z])
                # Transform grid index into Cartesian
                grid_cart = (([x, y, z]./grid)'*primvec)[1,:]
                # check if grid_cart is within sphere
                if norm(grid_cart - origin) <= radius
                    # Check if point is out of bounds. If so, use adjacent cell's value.
                    (i1, i2, i3) = (x, y, z)
                    if x < 1
                        i1 = grid[1]+x
                    elseif x > grid[1]
                        i1 = x-grid[1]
                    end
                    if y < 1
                        i2 = grid[2]+y
                    elseif y > grid[2]
                        i2 = y-grid[2]
                    end
                    if z < 1
                        i3 = grid[3]+z
                    elseif z > grid[3]
                        i3 = z-grid[3]
                    end
                    sphere_sum += data[i1,i2,i3]^2
                    count += 1
                end
            end
        end
    end
    
    # Calculate partial electron density
    total_sphere_sum = sum(sphere_sum*vox_vol)

    # Returns unnormalized partial electron density, unnormalized total electron density, Psphere
    return (total_sphere_sum, total, total_sphere_sum/total)
end

# Reads a void list (format: string float float float) and exports their xyz coordinates
# as a matrix for DFTraMO_overlap()
function readvoidlist(filename::AbstractString)
    open(filename) do f
        num_lines = countlines(f)
        seekstart(f)
        voidlist = zeros(Float64, num_lines, 3)
        for i in 1:num_lines
            voidlist[i,:] = parse.(Float64, split(readline(f))[2:4])
        end
        return voidlist
    end
end

"""
    writePsphere(
        start_xsf::Integer,
        end_xsf::Integer,
        xsf_prefix::AbstractString,
        xsf_suffix::AbstractString,
        voidlistfile::AbstractString,
        radius::Real,
        outputfilename::AbstractString
    ) -> Nothing

Computes and writes Psphere values for a set of DFT-raMO xsfs to a text file.

# Arguments
- `start_xsf`: Highest (first) electron count in DFT-raMO run.
- `end_xsf`: Lowest (last) electron count in DFT-raMO run.
- `xsf_prefix`: Prefix before electron count in xsf filenames.
- `xsf_suffix`: Suffix after electron count in xsf filenames.
- `voidlistfile`: Filename for coordinates for origins of Pspheres.
- `radius`: in angstroms.
- `outputfilename`: Filename for Psphere outputs.
"""
function writePsphere(
    start_xsf::Integer,
    end_xsf::Integer,
    xsf_prefix::AbstractString,
    xsf_suffix::AbstractString,
    voidlistfile::AbstractString,
    radius::Real,
    outputfilename::AbstractString
)
    voidlist = readvoidlist(voidlistfile)
    open(outputfilename, write=true) do f
        println(f, "Spherical DFTraMO overlaps with radius ", radius)
        println(f, "integrated_sphere \tintegrated_total_XSF \tsphere/total \torigin_coord")
        count = 1
        for i in start_xsf:-2:end_xsf
            (integrated_sphere, integrated_total, sph_total) = DFTraMO_overlap(
                string(xsf_prefix, i, xsf_suffix),
                voidlist[count,:],
                radius
            )
            println(
                f, integrated_sphere, "\t", integrated_total, "\t",
                sph_total, "\t", voidlist[count,:]
            )
            count += 1
        end
    end
end
