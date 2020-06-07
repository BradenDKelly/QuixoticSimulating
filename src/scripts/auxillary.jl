#####
#
#            Hodgepodge of functions  - should be NO STRUCTS HERE BRADEN.
#                    BRADEN - THERE SHOULD BE ZERO STRUCTS HERE.
#                     BRADEN YOU BETTER DAMN WELL BE LISTENING!
#
##########
using StaticArrays

include("structs.jl")

abstract type ForceField end
abstract type Gromacs <: ForceField end
abstract type GAFF <: ForceField end
abstract type AngleType <: ForceField end
abstract type DihedralType <: ForceField end

function potential_lrc(ρ, r_cut)
    """Calculates long-range correction for Lennard-Jones potential per atom."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * ((8.0 / 9.0) * sr3^3 - (8.0 / 3.0) * sr3) * ρ
end

function random_translate_vector(dr_max::Float64, old::SVector, box::Float64)

    # In: dr_max, old
    # Out: SVector{3,Float64}

    zeta = rand(Float64, 3)    # Three uniform random numbers in range (0,1)
    zeta = zeta .- 0.5         #! Now in range (-1,+1)
    return PBC(old + zeta .* dr_max, box) #! Move to new position

end #random_translate_vector

"""Evaluates if we keep a translation move or not"""
function Metropolis(delta)
    if delta < 0.0
        return true
    elseif exp(-delta) > rand()
        return true
    else
        return false
    end
end

#test_LJ()

""" returns max and min values in an array of SVectors"""
function maxmin(array::Vector)
    maxVal = maximum(array[1])
    minVal = minimum(array[1])

    for svector in array
        if maximum(svector) > maxVal
            maxVal = maximum(svector)
        end
        if minimum(svector) < minVal
            minVal = minimum(svector)
        end
    end

    return minVal, maxVal
end

"""Returns the center of mass of an array of atoms and array of masses"""
function COM(atoms::Vector, masses)
    totalMass = sum(masses)
    numerator = sum(atoms .* masses)
    #println("From inside COM: ", numerator ./ totalMass)
    return numerator ./ totalMass
end
#COM([[1,2,3],[2,3,4],[0,1,2]],[1,1,100])

"""Returns the Fortran equivalent of the same name"""
function MATMUL(ai, db)
    # In: ai:3x3 SMatrix
    #     db: 3x1 vector
    # Out: SVector{3}
    return SVector(dot(db, ai[:, 1]), dot(db, ai[:, 2]), dot(db, ai[:, 3]))
end

# TODO REmove duplicate COM and shift

"Calculates the Center of Mass of a molecule"
function Center_of_Mass( atom_coords, mol_mass) #result (COM)
#!============================================================================================================

 denominator = 0.0       #! sum of atom masses
 numerator = zeros(3)         #! weighted sum of atom masses and positions

for i = 1:length(mol_mass)
	numerator .+=  atom_coords[i] * mol_mass[i]
	denominator += mol_mass[i]

end

return SVector(numerator / denominator...)

end #Center_of_Mass
#!=============================

#!============================================================================================================
function Shift_COM_to_Zero!( atm_coords,COM )

#!============================================================================================================

# real, dimension(3), intent(in)                         			:: COM
# integer, intent(in)                                    			:: Atoms_in_molecule
# real, intent(inout), dimension(3,atoms_in_molecule)    			:: atm_coords							!atm = atom, qq = charge
# integer                                                			:: i
#!-----------------------------------------------------------------------------------------------------------

#! This is meant to take the body fixed coordinates of the atoms for a molecule and shift them so that the center of mass is at zero. This is necessary for proper rotation of the molecule.
#! This should only be called at the start of the run after getting *.pdb file coordinates.

for i = 1:length(atm_coords)

    atm_coords[i] = SVector(atm_coords[i] .- COM[:])

end

return nothing

end #Shift_COM_to_Zero
