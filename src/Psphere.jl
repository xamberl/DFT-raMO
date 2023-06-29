"""
    Psphere(
        datagrid::RealDataGrid,
        origin::Vector{Float64},
        rsphere::Float64
        )
Finds the fraction of electron density in a sphere of radius `rsphere`` (Å) at point `origin` (fractional).
"""
function Psphere(
    datagrid::RealDataGrid,
    origin::Vector{Float64},
    rsphere::Float64
    )
    # Calculate volume of single voxel
    vox_vol = volume(datagrid)/length(datagrid)
    # Find unnormalized electron count
    total_eden = sum(datagrid.data.^2)*vox_vol
    # Calculate bounds of box that inscribes the sphere
    lowerbound = (origin.-rsphere)*Electrum.ANG2BOHR
    upperbound = (origin.+rsphere)*Electrum.ANG2BOHR
    # Calculate indices of box in the voxel grid
    lower_i = round.(Int, (lowerbound'/datagrid.basis)'.*size(datagrid))
    upper_i = round.(Int, (upperbound'/datagrid.basis)'.*size(datagrid))
    sphere_sum = 0
    for x in lower_i[1]:upper_i[1], y in lower_i[2]:upper_i[2], z in lower_i[3]:upper_i[3]
        # Transform grid index into Cartesian (Å)
        grid_cart = (([x, y, z]./size(datagrid))'*datagrid.basis)*Electrum.BOHR2ANG
        # check if grid_cart is within sphere
        if norm(grid_cart'.-origin) <= rsphere#((grid_cart[1]-origin[1])^2+(grid_cart[2]-origin[2])^2+(grid_cart[3]-origin[3])^2)^.5 <= rsphere
            # Check if point is out of bounds. If so, use adjacent cell's value.
            coord = [x,y,z]
            for n in eachindex(coord)
                if coord[n] < 1
                    coord[n] = size(datagrid)[n]+coord[n]
                elseif coord[n] > size(datagrid)[n]
                    coord[n] = coord[n]-size(datagrid)[n]
                end
            end
            sphere_sum += datagrid.data[CartesianIndex(Tuple(coord))]^2
        end
    end
    # Calculate partial electron density
    return(sphere_sum*vox_vol, total_eden, sphere_sum*vox_vol/total_eden)
end

"""
    psphere_eval(psphere::Vector{Float64}, super::Supercell, site_list::Vector{Int})
    
Prints Psphere analysis to the terminal.
"""
function psphere_eval(psphere::Vector{Float64}, super::Supercell, site_list)
    m = maximum(psphere)
    a = findall(x->x<0.15*m, psphere)
    if !isempty(a)
        println("The following atoms have Pspheres <15% of the maximum (", m, "):")
        for n in a
            site = Vector(super.atomlist.basis*Electrum.BOHR2ANG*super.atomlist[site_list[n]].pos)
            println("Atom ", site_list[n], ": ", @sprintf("Psphere: %.3f", psphere[n]), @sprintf(" at site [%.3f, %.3f, %.3f]", site[1], site[2], site[3]))
        end
    end
    return a
end