"""
    make_target_AO(make_target_AO(atom_site::Int, orbital_to_use::Int, total_num_orbitals::Int) -> Vector

Returns a vector of length total_num_orbitals with "1" in the corresponding atomic orbital
"""
function make_target_AO(atom_site::Integer, target_orbital::Integer, super::Supercell)
    # psi_target has a length of total number of orbitals in the supercell
    psi_target = zeros(sum(super.orbitals))
    psi_target[sum(super.orbitals[1:atom_site])-super.orbitals[atom_site]+target_orbital] = 1
    return psi_target
end

"""
    make_target_cluster_sp(site_list::Vector{Vector{Real}}, radius::Real, site_num::Int, super::Supercell)

Returns a vector of length total_num_orbitals with "1" in the corresponding s & p atomic orbitals
if they are within specified radius to the void.
"""
function make_target_cluster_sp(site_list::Vector{Vector{Float64}}, radius::Real, site_num::Int, super::Supercell)
    # psi_target has a length of total number of orbitals in the supercell
    psi_target = zeros(sum(super.orbitals))
    # Loops through every atom and checks to see if it's within the radius to the void site
    for n in 1:length(super.atomlist)
        for j in -1:2
            for k in -1:2
                for l in -1:2
                    # Translation of atom in cartesian coordinates
                    new_pos = (super.atomlist.basis*super.atomlist[n].pos .+ super.atomlist.basis*[j,k,l])*Electrum.BOHR2ANG
                    check_distance = norm(new_pos-site_list[site_num])
                    if check_distance <= radius
                        # Δr is used to weigh the p orbitals in each direction
                        # Is there a reason why this is negative?
                        Δr = -(new_pos-site_list[site_num])/check_distance
                        # Fill in corresponding s orbital
                        psi_target[sum(super.orbitals[1:n])-super.orbitals[n]+1] = 0.5^0.5
                        # Fill in corresponding p orbitals, if any
                        if (super.orbitals[n] > 1)
                            for i in 1:3
                                psi_target[sum(super.orbitals[1:n])-super.orbitals[n]+1+i] = 0.5^0.5*Δr[i]
                            end
                        end
                    end
                end
            end
        end
    end
    iszero(norm(psi_target)) ? error("Target is empty. Try increasing radii.") : nothing
    return psi_target/norm(psi_target)
end

"""
    make_target_hybrid()

Makes a custom hybrid with manually specified coefficients from a cluster file that corresponds to an atom
"""
#==
Currently echoes MATLAB version of DFT-raMO. Cluster file should have one Float per line,
in multiples of 9.
==#
function make_target_hybrid(cluster_list::Vector{Vector{Float64}}, atom_num::Int, super::Supercell)
    num_targets = length(cluster_list)/9
    psi_target = zeros(sum(super.orbitals), num_targets)
    AO_number = sum(super.orbitals[1:atom_num])-super.orbitals[atom_num]
    # Loops through targets
    for n in 1:num_targets
        # Assigns s coefficient
        psi_target[AO_number+1,n] = cluster_list[(n-1)*9+1]
        # Assigns p coefficients
        if super.orbitals[atom_num] > 1
            for p in 1:3
                psi_target[AO_number+1+p,n] = cluster_list[(n-1)*9+1+p]
            end
        end
        # Assigns d coefficients
        if super.orbitals[atom_num] > 4
            for d in 1:3
                psi_target[AO_number+4+d,n] = cluster_list[(n-1)*9+4+d]
            end
        end
    end
    return psi_target
end