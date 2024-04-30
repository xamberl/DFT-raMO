"""
    Psphere(
        datagrid::RealDataGrid,
        origin::AbstractVector{<:Real},
        rsphere::Real
    )

Finds the fraction of electron density in a sphere of radius `rsphere`` (Å) at point `origin` (fractional).
"""
function Psphere(
    datagrid::RealDataGrid,
    origin::AbstractVector{<:Real},
    rsphere::Real
)
    # Calculate volume of single voxel
    vox_vol = volume(datagrid)/length(datagrid)
    # Find unnormalized electron count
    total_eden = sum(datagrid.data.^2)*vox_vol
    # Calculate bounds of box that inscribes the sphere
    lowerbound = (origin.-rsphere)*Electrum.ANG2BOHR
    upperbound = (origin.+rsphere)*Electrum.ANG2BOHR
    # Calculate indices of box in the voxel grid
    i1 = round.(Int, (lowerbound'/datagrid.basis)'.*size(datagrid))
    i2 = round.(Int, (upperbound'/datagrid.basis)'.*size(datagrid))
    # check whether i1 or i2 is the lower bound
    if i1[1] < i2[1]
        lower_i = i1
        upper_i = i2
    else
        lower_i = i2
        upper_i = i1
    end
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
    psphere_eval(psphere::AbstractVector{<:Real}, site::AbstractVector{SVector{3, Float64}})
    
Prints Psphere analysis to the terminal.
"""

function psphere_eval(psphere::AbstractVector{<:Real}, site::AbstractVector{SVector{3, Float64}})
    m = maximum(psphere)
    a = findall(x->x<0.15*m, psphere)
    if !isempty(a)
        println("The following sites have Pspheres <15% of the maximum (", m, "):")
        for n in a
            println("Site ", n, ": ", @sprintf("Psphere: %.3f", psphere[n]), @sprintf(" at site [%.3f, %.3f, %.3f]", site[n][1], site[n][2], site[n][3]))
        end
    end
    return a
end

"""
    print_psphere_terminal(iter, psphere, site)

Prints Psphere values to the terminal. The color of Psphere is arbitrarily color-coded.
    Green:        Psphere ≥ 0.5
    Yellow: 0.1 ≤ Psphere < 0.5
    Red:          Psphere < 0.1
"""
function print_psphere_terminal(iter, num_raMO, psphere, site, io)
    col = Crayon(foreground = :green)
    psphere < 0.5 ? col = Crayon(foreground = :light_yellow) : nothing
    psphere < 0.1 ? col = Crayon(foreground = :light_red) : nothing
    println(
        iter,
        num_raMO,
        ", Psphere: ",
        col, @sprintf("%.3f", psphere),
        Crayon(foreground = :default),
        @sprintf(" at site [%.3f, %.3f, %.3f]", site[1], site[2], site[3]))
    println(io, num_raMO, "\t", psphere, "\t", @sprintf(" at site [%.3f, %.3f, %.3f]", site[1], site[2], site[3]))
    flush(stdout)
end
    
function psphere_graph(psphere::AbstractVector{<:Real}, num_raMO::Integer, rsphere::Real)
    @info "UnicodePlots not loaded: consider loading it for Psphere logging."
end

"""
    mp_lcao(sites::Vector{Int}, xtal::PeriodicAtomList)

Calculates the geometric midpoint for LCAOs. See
https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
"""
function mp_lcao(sites::Vector{Int}, xtal::PeriodicAtomList)
    θ = reshape(reduce(vcat, [(xtal[s].pos*2*pi) for s in sites]), 3, length(sites))
    ξ = cos.(θ); ξ = [mean(ξ[x, :]) for x in 1:3]
    ζ = sin.(θ); ζ = [mean(ζ[x, :]) for x in 1:3]
    return [atan(-ζ[x], -ξ[x]) + pi for x in 1:3]./(2*pi)
end