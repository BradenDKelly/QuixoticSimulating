


""" Calculates total LJ && Wolf Summation potential energy of the system.
 Double Counts. """
function potential(
    system,
    tot,
    ewald::EWALD,
    qq_q::Vector,
    qq_r::Vector,
    string,
    string2 #triggers wolf summations using double strings
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = system.r_cut
    #@assert r_cut == 10.0
    box = system.box
    thisMol_thisAtom = system.thisMol_theseAtoms
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    totReal = 0.0
    for i = 1:length(rm)
        ener, overlap =
            EwaldReal(qq_r, qq_q,ewald.kappa, box,thisMol_thisAtom, i, system)
        totReal += ener
        if overlap
            exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    reall = totReal
    println("Coulomb contribution is: ", reall)

    prefactor = 0.0
    for i=1:length(qq_q)
        for j=1:length(qq_q)
            prefactor += qq_q[i]*qq_q[j] * erfc(ewald.kappa*r_cut) / r_cut
        end
    end
    prefactor *= -1
    prefactor2 = (erfc(ewald.kappa*r_cut) / 2 / r_cut +ewald.kappa / sqrt(π) ) *
                    dot(qq_q,qq_q)
    tot.energy += (prefactor - prefactor2)*ewald.factor
    tot.coulomb += (prefactor - prefactor2)*ewald.factor

    println("total energy is: ", tot.energy)
    println("total coulomb is: ", tot.coulomb)
    println("total LJ energy is: ", tot.energy - tot.coulomb)


    return tot

end

""" Calculates total LJ && Bare Coulomb potential energy of the system.
 Double Counts. """
function potential(
    system,
    tot,
    ewald::EWALD,
    qq_q::Vector,
    qq_r::Vector,
    string,
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = system.r_cut
    #@assert r_cut == 10.0
    box = system.box
    thisMol_thisAtom = system.thisMol_theseAtoms
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real

    totReal = 0.0
    for i = 1:length(rm)
        ener, overlap =
            CoulombReal(qq_r, qq_q, box, i, system)
        totReal += ener
        if overlap
            exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    reall = totReal
    println("Coulomb contribution is: ", reall)


    return tot, LJ, reall

end

"Real Ewald contribution with molecular cutoff"
#function EwaldReal(diff,coord1, coord2,q1,q2, L, rcut2, kappa)
function CoulombReal(
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    box::Float64,
    chosenOne::Int64,
    system::Requirements
)
    ####
    #
    #    Some prep stuff
    #    --------------------------------
    #    system: a struct with arrays for molecular and atomic coords and other stuff
    #    chosenOne: molecule that was moved
    #    qq_q: array of atomic partial charges
    #    qq_r: array of atomic coordinates. Each index is a Vector with x, y, z
    #
    #############
    ri = system.rm[chosenOne]                       # x, y, z coordinates for COM of molecule "chosenOne"
    thisMol_thisAtom = system.thisMol_theseAtoms    # length number of molecules
    start_a = thisMol_thisAtom[chosenOne][1]        # this atom starts molecule "chosenOne"
    end_a = thisMol_thisAtom[chosenOne][2]          # this atom ends molecule "chosenOne"
    rMol = system.rm                                # array of molecule COM coords
    r_cut = system.r_cut                            # site-site cutoff

    diameter = r_cut * 0.25  + 5.0        # made up diameter of molecule
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff
    rm_cut_box_sq = rm_cut_box^2          # Molecular cutoff squared
    r_cut_sq = r_cut^2                    # atomic cutoff squared

    # quick double check on atomic cutoffs
    @assert r_cut == 10.0
    #@assert r_cut_sq == 100.0

    pot = 0.0
    # Allocate for speed
    rab = SVector{3, Float64}(0.0, 0.0, 0.0)
    rij = SVector{3, Float64}(0.0, 0.0, 0.0)
    rj = SVector{3, Float64}(0.0, 0.0, 0.0)
    ra = SVector{3, Float64}(0.0, 0.0, 0.0)
    rb = SVector{3, Float64}(0.0, 0.0, 0.0)

    overlap = false
    ovr = 1.0

    for (j, rj) in enumerate(rMol) # cycle through all molecules
        if j == chosenOne # skip self interactions
            continue
        end

        #Molecular mirror image separation
        for k = 1:3
            rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij2 = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

        #Molecular cutoff check
        if rij2 < rm_cut_box_sq

            #Loop over all atoms in molecule A
            for a = start_a:end_a
                ra = qq_r[a]

                start_b = thisMol_thisAtom[j][1]  # first atom in molecule j
                end_b = thisMol_thisAtom[j][2]    # last atom in molecule j

                # Loop over all atoms in molecule B
                for b=start_b:end_b
                    rb = qq_r[b]       # coordinates of charge b

                    # Atomic mirror image separation
                    for k = 1:3
                        rab = @set rab[k] = vector1D(ra[k], rb[k], box)
                    end
                    rab2 = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    # Overlap Check
                    if (rab2 < ovr) && (qq_q[a] * qq_q[b] < 0)
                        println("got one")
                        return 0.0, true
                    # Atomic cutoff check
                    elseif rab2 < r_cut_sq
                        #println(a, "   ", b)
                        #rab_mag =
                        pot += qq_q[a] * qq_q[b] / sqrt(rab2)
                    else
                        pot += 0.0
                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end # if statement for molecular cutoff
    end # loop over molecules
    return pot , overlap
end

""" Calculates total LJ && Electrostatic potential energy of the system.
 Double Counts. """
function potential(
    system,
    tot,
    ewalds::EWALD,
    qq_q::Vector,
    qq_r::Vector,
    #kfacs::Vector,
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = system.r_cut
    box = system.box
    thisMol_thisAtom = system.thisMol_theseAtoms
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real, recip = 0.0, 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)

    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    #ener = EwaldReal(qq_r, qq_q, ewald.kappa, box, thisMol_thisAtom, 2, system)
    #=
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    kappa::Real,
    box::Float64,
    thisMol_thisAtom::Vector,
    chosenOne::Int64,
    system::Requirements
    =#
    totReal = 0.0
    for i = 1:length(rm)
        ener, overlap =
            EwaldReal(qq_r, qq_q, ewalds.kappa, box, thisMol_thisAtom, i, system)
        totReal += ener
        if overlap
            exit()
        end
    end
    tot.energy += totReal / 2 * ewalds.factor # divide by 2 to account for double counting
    tot.coulomb += totReal / 2 * ewalds.factor
    reall = totReal / 2 * ewalds.factor

    println("Total real Ewald is: ", totReal / 2 * ewalds.factor)

    # Recipricol
    #@time recipEnergy = RecipLong(system, ewald, qq_r, qq_q) * ewald.factor #RecipLong(system, ewald, qq_r, qq_q, kfacs)
    recipEnergy, ewalds = RecipLong(system, ewalds, qq_r, qq_q)  #RecipLong(system, ewald, qq_r, qq_q, kfacs)
    recipEnergy *= ewalds.factor
    #@btime RecipLong(system, ewald, qq_r, qq_q) * ewald.factor
    #println("first: ", recipEnergy)

    println("Total recipricol Ewald is: ", recipEnergy)

    tot.energy += recipEnergy #* ewald.factor
    tot.coulomb += recipEnergy #* ewald.factor
    tot.recipOld = recipEnergy #* ewald.factor
    tot.recip = recipEnergy #* ewald.factor
    # Self

    println("Self energy: ", EwaldSelf(ewalds, qq_q))
    tot.energy += EwaldSelf(ewalds, qq_q)
    # tinfoil boundary
    #@time tinfoil = TinfoilBoundary(system, ewald, qq_q, qq_r)
    #println("tinfoil boundaries: ", tinfoil * ewald.factor)
    #tot.energy += tinfoil * ewald.factor
    #tot.coulomb += tinfoil * ewald.factor
    return tot, LJ, reall, recipEnergy, ewalds

end


""" Calculates total LJ potential energy of the system. Double Counts. """
function potential(system, tot)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0

    #########
    #
    #     Calculate LJ
    #
    ##################

    for i = 1:length(system.rm)
        ener, vir = LJ_poly_ΔU(i, system)
        tot.energy += ener
        tot.virial += vir
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2

    return tot

end



""" Calculates LJ potential between particle 'i' and the other N-1 particles"""
function LJ_ΔU(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    r = system.r
    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    diff = SVector{3}(0.0, 0.0, 0.0)
    rcut_sq = system.r_cut^2
    pot, vir = 0.0, 0.0

    for (j, atom) in enumerate(r)
        if j == i
            continue
        end

        @inbounds for k = 1:3
            diff = @set diff[k] = vector1D(r[i][k], atom[k], box)
        end

        rij_sq = diff[1] * diff[1] + diff[2] * diff[2] + diff[3] * diff[3]

        if rij_sq > rcut_sq
            pot += 0.0
            vir += 0.0
        else
            # possibly needed in the case of random starting config to avoid
            # numerical overflow from overlapped particles. Not needed
            # if starting from lattice.
            #if rij_sq < 0.6^2
            #    rij_sq = 0.6^2
            #end

            sr2 = σ[j]^2 / rij_sq
            sr6 = sr2^3
            sr12 = sr6^2
            pot += ϵ[j] * (sr12 - sr6)  # LJ pair potentials (cut but not shifted)
            vir += ϵ[j] * (2 * sr12 - sr6)  # LJ pair virials
        end

    end

    return pot * 4.0, vir * 24.0 / 3.0
end

""" LJ energy of particle i with all other particles. Real units of σ and ϵ """
function LJ_poly_ΔU(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    ri = system.rm[i]
    rMol = system.rm

    startAtom = system.thisMol_theseAtoms[i][1]
    endAtom = system.thisMol_theseAtoms[i][2]
    list_type = system.atomTypes
    moli_type = system.atomTypes[startAtom:endAtom]

    sF = system.thisMol_theseAtoms

    ra = system.ra[startAtom:endAtom]
    rb = system.ra

    table = system.table
    box = system.box
    r_cut = system.r_cut
    diameter = 0 #r_cut * 0.25 + 5 #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2                   # Potential cutoff squared in sigma=1 units

    rᵢⱼ = rab = SVector{3}(0.0, 0.0, 0.0)

    pot, vir = 0.0, 0.0

    # cycle through all molecules
    @inbounds for (j, rj) in enumerate(rMol)
        if j == i
            continue # skip if molecule j is same as i
        end

        """Molecular mirror image separation"""
        @inbounds for k = 1:3
            rᵢⱼ = @set rᵢⱼ[k] = vector1D(ri[k], rj[k], box)
        end

        rᵢⱼ² = rᵢⱼ[1] * rᵢⱼ[1] + rᵢⱼ[2] * rᵢⱼ[2] + rᵢⱼ[3] * rᵢⱼ[3]

        if rᵢⱼ² < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            @inbounds for a = 1:3

                """Loop over all atoms in molecule B"""
                @inbounds for b = sF[j][1]:sF[j][2]

                    """Atomic mirror image separation"""
                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab² = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    ϵᵢⱼ = table.ϵᵢⱼ[moli_type[a], list_type[b]]
                    if rab² < (r_cut_sq + 100) && ϵᵢⱼ > 0.001
                        # this uses only molecular cutoff

                        σᵢⱼ = table.σᵢⱼ[moli_type[a], list_type[b]]

                        σ² = σᵢⱼ^2 / rab²             # (sigma/rab)**2
                        σ⁶ = σ²^3
                        σ¹² = σ⁶^2
                        pot += ϵᵢⱼ * (σ¹² - σ⁶)
                        virab = ϵᵢⱼ * (2.0 * σ¹² - σ⁶)
                        fab = rab * virab * σ²
                        vir += dot(rᵢⱼ, fab)

                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end
    end

    return pot * 4, vir * 24 / 3.0
end # lj_poly_ΔU

""" Calculates LJ potential between particle 'i' and the other N-1 particles.
 Cut & Shifted Potential"""
function LJ_poly_ΔU_shifted(i::Int, system::Requirements)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # Uses σ = ϵ = 1 units
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    # Cutoff distance and force-shift parameters (all private) chosen as per the reference:
    # S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, Phys Rev E, 65, 041205 (2002)
    r_cut = 2.612 # in sigma=1 units, where r_cut = 1.2616 nm, sigma = 0.483 nm
    sr_cut = 1.0 / r_cut
    sr_cut6 = sr_cut^6
    sr_cut12 = sr_cut6^2
    lambda1 = 4.0 * (7.0 * sr_cut6 - 13.0 * sr_cut12)
    lambda2 = -24.0 * (sr_cut6 - 2.0 * sr_cut12) * sr_cut
    sr2_ovr = 1.77

    diameter = 1.327441  #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2              # squared
    r_cut_sq = r_cut^2                   # Potential cutoff squared in sigma=1 units

    """
     ri   = deepcopy(system.rm[i])
     rMol = deepcopy(system.rm)
     rMol = deleteat!(rMol,i)

     startAtom = deepcopy(system.thisMol_theseAtoms[i][1])
     endAtom   = deepcopy(system.thisMol_theseAtoms[i][2])

     sF = deepcopy(system.thisMol_theseAtoms)
     sF = deleteat!(sF,i)

     ra = deepcopy(system.ra[startAtom:endAtom])
     rb = deepcopy(system.ra)
     rb = deleteat!(rb,startAtom:endAtom)
     #println(length(system.rm))
    """

    ri = system.rm[i]
    rMol = system.rm

    startAtom = system.thisMol_theseAtoms[i][1]
    endAtom = system.thisMol_theseAtoms[i][2]

    sF = system.thisMol_theseAtoms

    ra = system.ra[startAtom:endAtom]
    rb = system.ra

    ϵ = system.ϵ
    σ = system.σ
    box = system.box
    rcut_sq = system.r_cut^2

    rij = rab = SVector{3}(0.0, 0.0, 0.0)

    pot, vir = 0.0, 0.0

    for (j, rj) in enumerate(rMol)
        if j == i
            continue
        end

        @inbounds for k = 1:3
            rij = @set rij[k] = vector1D(ri[k], rj[k], box)
        end

        rij_sq = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

        if rij_sq < rm_cut_box_sq

            """Loop over all atoms in molecule A"""
            for a = 1:3
                """Loop over all atoms in molecule B"""
                for b = sF[j][1]:sF[j][2]
                    if (sF[j][2] - sF[j][1]) != 2
                        println("INdex wrong for molecule B")
                    end

                    @inbounds for k = 1:3
                        rab = @set rab[k] = vector1D(ra[a][k], rb[b][k], box)
                    end

                    rab_sq = rab[1] * rab[1] + rab[2] * rab[2] + rab[3] * rab[3]

                    sr2 = 1.0 / rab_sq             # (sigma/rab)**2
                    ovr = sr2 > sr2_ovr
                    if ovr
                        println("OVERLAP, OVERLAP!!!!!")
                    end
                    if rab_sq < r_cut_sq
                        #println(i, " ", j)
                        sr2 = 1.0 / rab_sq
                        rmag = sqrt(rab_sq)
                        sr6 = sr2^3
                        sr12 = sr6^2
                        pot += 4 * (sr12 - sr6) + lambda1 + lambda2 * rmag
                        virab = 24.0 * (2.0 * sr12 - sr6) - lambda2 * rmag
                        fab = rab * virab * sr2
                        vir += dot(rij, fab)

                    end # potential cutoff
                end # loop over atom in molecule b
            end # loop over atoms in molecule a
        end
    end

    return pot, vir / 3.0
end # lj_poly_ΔU


if 2==2

elseif occursin(lowercase(initialConfiguration), "cnf") #"cnf"  lowercase(initialConfiguration)
    rm, quat, box = ReadCNF("cnf_input.inp")
    nMol = length(rm)
    nAtoms = nMol * 3
    ρ = nMol / (box^3)
    println(" Initial configuration is from file: cnf_input.inp")
    println("boxsize is: ", box)
    println("density is: ", ρ)

    """
    A&T use a box centered at [0,0,0] whereas we use a box from 0->box in all
    dimensions. This next part shift all coordinates by the magnitude of the
    smallest coordinate so all coordinates are not between 0->box
    """
    xl = yl = zl = 0.0
    for (i, mol) in enumerate(rm)
        global xl, yl, zl
        if i == 1
            xl = mol[1]
            yl = mol[2]
            zl = mol[3]
        else
            if mol[1] < xl
                xl = mol[1]
            end
            if mol[2] < yl
                yl = mol[2]
            end
            if mol[3] < zl
                zl = mol[3]
            end
        end
    end
    xl = abs(xl)
    yl = abs(yl)
    zl = abs(zl)
    for (i, mol) in enumerate(rm)
        rm[i] = mol .+ [xl, yl, zl]
    end

elseif occursin(lowercase(initialConfiguration), "nist")
    # THis file stores a single configuration of 100-750 SPC/E water molecules.
    # This is for testing only. NIST has results for energy contributions.
    # load the configuration, get the atom types, charges, coordinates.
    # Return to here and test.
    filename = "spce_sample_config_periodic1.txt"
    filename = joinpath( pwd(),  "src", "topology_files","coord750.txt")
    println(filename)
    qq_r, qq_q, rm, ra, atomTracker, box, atomName, atomType =
        ReadNIST(filename)
    # make LJ table of values.
    σ_O = 0.316555789 * 10.0 # convert nm to Å
    σ_H = 0.0 # nm
    ϵ_O = 78.1974311 # K   (ϵ/kᵦ)
    ϵ_H = 0.0 # K   (ϵ/kᵦ)

    xl = yl = zl = 0.0
    for (i, mol) in enumerate(rm)
        global xl, yl, zl
        if i == 1
            xl = mol[1]
            yl = mol[2]
            zl = mol[3]
        else
            if mol[1] < xl
                xl = mol[1]
            end
            if mol[2] < yl
                yl = mol[2]
            end
            if mol[3] < zl
                zl = mol[3]
            end
        end
    end
    xl = abs(xl)
    yl = abs(yl)
    zl = abs(zl)
    for (i, mol) in enumerate(rm)
        rm[i] = mol .+ [xl, yl, zl]
    end
    for (i, part) in enumerate(ra)
        ra[i] = part .+ [xl, yl, zl]
        qq_r[i] = part .+ [xl, yl, zl]
    end
else
    rm = [SVector{3,Float64}(rand(), rand(), rand()) .* box for i = 1:nMol]

end


if occursin(lowercase(initialConfiguration), "cnf")
    ϵ = σ = ones(1, 1)
else
    #=
    ϵ = [ϵ_O, ϵ_H]
    σ = [σ_O, σ_H]
    molNames = ["Wat" for i = 1:length(rm)]
    molTypes = [1 for i = 1:length(rm)]
    =#
end

"""Table with LJ FF parameters, not in use atm"""
struct Tables2 <: ForceField
    ϵij::Array{Float64,2}
    σij::Array{Float64,2}
    αij::Array{Float64,2}
    Tables2(ϵij=zeros(2,2), σij=zeros(2,2),αij=zeros(2,2)) = new(ϵij, σij,αij)
end


#=
""" Recipricol space Ewald energy Move older style, now used"""
function RecipMove(
    system::Requirements,
    ewalds::EWALD,
    r_old::Vector,
    r_new::Vector,
    qq_q::Vector,
)
    #=
    Intended to calculate the recipricol for the molecule that was moved, not then
    entire system.

    This is twice as long as it should be.
    It calculates the fourier contibution for both old and new positions.
    TODO store old contributions
    =#
    #cfac = ewald.cfac
    kxyz = ewalds.kxyz[:]
    energy = float(0.0)
    L = system.box
    twopi = 2.0 * π
    nk = ewalds.nk
    n = length(r_old)
    @assert n == 3
    k_sq_max = ewalds.k_sq_max
    @assert k_sq_max == 27
    @assert nk == 5

    eikx_n = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky_n = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)

    eikz_n = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #similar(eiky_n) #  OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)
    eikx_o = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #similar(eikx_n) #  OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky_o = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #similar(eiky_n) # OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz_o = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #similar(eiky_n) # OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    for j = 1:n
        # calculate kx,ky,kz =1,2 explicitly
        eikx_n[j, 0] = 1.0 + 0.0 * im
        eiky_n[j, 0] = 1.0 + 0.0 * im
        eikz_n[j, 0] = 1.0 + 0.0 * im
        eikx_o[j, 0] = 1.0 + 0.0 * im
        eiky_o[j, 0] = 1.0 + 0.0 * im
        eikz_o[j, 0] = 1.0 + 0.0 * im

        eikx_n[j, 1] =
            cos(twopi * (r_new[j][1]) / L) + sin(twopi * (r_new[j][1]) / L) * im
        eiky_n[j, 1] =
            cos(twopi * (r_new[j][2]) / L) + sin(twopi * (r_new[j][2]) / L) * im
        eikz_n[j, 1] =
            cos(twopi * (r_new[j][3]) / L) + sin(twopi * (r_new[j][3]) / L) * im
        eikx_o[j, 1] =
            cos(twopi * (r_old[j][1]) / L) + sin(twopi * (r_old[j][1]) / L) * im
        eiky_o[j, 1] =
            cos(twopi * (r_old[j][2]) / L) + sin(twopi * (r_old[j][2]) / L) * im
        eikz_o[j, 1] =
            cos(twopi * (r_old[j][3]) / L) + sin(twopi * (r_old[j][3]) / L) * im

        eiky_n[j, -1] = conj(eiky_n[j, 1])
        eikz_n[j, -1] = conj(eikz_n[j, 1])
        eiky_o[j, -1] = conj(eiky_o[j, 1])
        eikz_o[j, -1] = conj(eikz_o[j, 1])
    end
    # calculate remaining positive kx, ky, kz by recurrence
    for k = 2:(nk)
        for j = 1:n
            eikx_n[j, k] = eikx_n[j, k-1] * eikx_n[j, 1]
            eiky_n[j, k] = eiky_n[j, k-1] * eiky_n[j, 1]
            eikz_n[j, k] = eikz_n[j, k-1] * eikz_n[j, 1]
            eikx_o[j, k] = eikx_o[j, k-1] * eikx_o[j, 1]
            eiky_o[j, k] = eiky_o[j, k-1] * eiky_o[j, 1]
            eikz_o[j, k] = eikz_o[j, k-1] * eikz_o[j, 1]

            #eikx[j,-k] = conj( eikx[j,k] )
            eiky_n[j, -k] = conj(eiky_n[j, k])
            eikz_n[j, -k] = conj(eikz_n[j, k])
            eiky_o[j, -k] = conj(eiky_o[j, k])
            eikz_o[j, -k] = conj(eikz_o[j, k])
        end
    end
    # TODO This can be sped up. It is generating too many temp Arrays
    # It should be faster if devectorized, but that makes even more temp arrays
    # TODO can probably vectorize better if I swap i and l loops.

    # Calculate difference between the twopi

    for i = 1:length(ewalds.cfac)
        for l = 1:n
            ewalds.sumQExpNew[i] +=
                qq_q[l] * (
                    eikx_n[l, kxyz[i][1]] *
                    eiky_n[l, kxyz[i][2]] *
                    eikz_n[l, kxyz[i][3]] -
                    eikx_o[l, kxyz[i][1]] *
                    eiky_o[l, kxyz[i][2]] *
                    eikz_o[l, kxyz[i][3]]
                )
        end


        energy +=
            ewalds.cfac[i] * (
                real(conj(ewalds.sumQExpNew[i]) * ewalds.sumQExpNew[i]) -
                real(conj(ewalds.sumQExpOld[i]) * ewalds.sumQExpOld[i])
            )
    end

    #
    return energy * ewalds.factor, ewalds
end # function Recipricol Move
=#



#function RecipLong(boxSize, n, kappa, kfac, nk, k_sq_max, r, qq_q)
#=
""" Recipricol space Ewald energy NOT USED"""
function RecipLong(
    system::Requirements,
    ewald::EWALD,
    r::Vector,
    qq_q::Vector{Float64},
    kfac::Vector,
)
    energy = 0.0
    L = system.box
    twopi = 2.0 * π
    nk = ewald.nk
    n = length(r)
    k_sq_max = ewald.k_sq_max

    eikx = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    @inbounds @simd for j = 1:n
        # calculate kx,ky,kz =1,2 explicitly
        eikx[j, 0] = 1.0 + 0.0 * im
        eiky[j, 0] = 1.0 + 0.0 * im
        eikz[j, 0] = 1.0 + 0.0 * im

        eikx[j, 1] =
            cos(twopi * (r[j][1]) / L) + sin(twopi * (r[j][1]) / L) * im
        eiky[j, 1] =
            cos(twopi * (r[j][2]) / L) + sin(twopi * (r[j][2]) / L) * im
        eikz[j, 1] =
            cos(twopi * (r[j][3]) / L) + sin(twopi * (r[j][3]) / L) * im

        #eikx[j,-1] = conj(eikx[j,1])
        eiky[j, -1] = conj(eiky[j, 1])
        eikz[j, -1] = conj(eikz[j, 1])
    end
    # calculate remaining positive kx, ky, kz by recurrence
    for k = 2:(nk)
        @inbounds @simd for j = 1:n
            eikx[j, k] = eikx[j, k-1] * eikx[j, 1]
            eiky[j, k] = eiky[j, k-1] * eiky[j, 1]
            eikz[j, k] = eikz[j, k-1] * eikz[j, 1]

            #eikx[j,-k] = conj( eikx[j,k] )
            eiky[j, -k] = conj(eiky[j, k])
            eikz[j, -k] = conj(eikz[j, k])
        end
    end

    term = 0.0 + 0.0 * im
    #fill!(term,0.0+0.0*im) #term = zeros(ComplexF64,n)
    # fun party time - lets bring out the loops
    # TODO use other code which is single loop and @simd it
    # TODO only evaluate moved charges - huge speedup will occur from this
    # ~ 50x increase in speed for a 300 spce water system. Bigger speedups within
    # larger systems.
    exit()
    @inbounds for kx = 0:nk
        if kx == 0
            factor = 1.0
        else
            factor = 2.0
        end

        for ky = -nk:nk
            for kz = -nk:nk
                k_sq = kx * kx + ky * ky + kz * kz

                if k_sq <= k_sq_max && k_sq > 0
                    term = 0.0 + 0.0 * im
                    @inbounds @simd for l = 1:n
                        term +=
                            qq_q[l] * eikx[l, kx] * eiky[l, ky] * eikz[l, kz]
                    end
                    energy += factor * kfac[k_sq] * real(conj(term) * term)  # / L
                end # if
            end # kz for loop
        end # ky for loop
    end # kx for loop
    #energy *= qq #/ L
    return energy
end # function
=#


#=
"Real Ewald contribution - without molecular cutoff"
#function EwaldReal(diff,coord1, coord2,q1,q2, L, rcut2, kappa)
@fastmath function EwaldReal(
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    kappa::Real,
    box::Real,
    thisMol_thisAtom::Vector,
    chosenOne::Int,
    r_cut::Real,
)
    ####
    #
    #    Some prep stuff
    #
    #############
    start_i = thisMol_thisAtom[chosenOne][1]
    end_i = thisMol_thisAtom[chosenOne][2]
    #r_cut = 10.0
    r_cut_sq = r_cut^2
    pot = 0.0
    rij = SVector{3}(0.0, 0.0, 0.0)
    rᵢ = similar(rij)
    rⱼ = similar(rij)



    """Loop over all atoms in molecule A"""
    @inbounds for i = start_i:end_i
        rᵢ = qq_r[i]
        """Loop over all atoms in molecule B"""
        @inbounds for (j, rⱼ) in enumerate(qq_r)
            if j in start_i:end_i
                continue
            end # keep mol i from self interacting
            @inbounds for k = 1:3

                rij = @set rij[k] = vector1D(rᵢ[k], rⱼ[k], box)
            end
            rᵢⱼ² = rij[1] * rij[1] + rij[2] * rij[2] + rij[3] * rij[3]

            if rᵢⱼ² < r_cut_sq
                rᵢⱼ = sqrt(rᵢⱼ²)
                pot += qq_q[i] * qq_q[j] * erfc(kappa * rᵢⱼ) / rᵢⱼ

            else
                pot += 0.0
            end # potential cutoff
        end # loop over atom in molecule b
    end # loop over atoms in molecule a
    return pot
end
=#

#=
PrintOutput(
    system,
    totProps,
    atomType,
    atomName,
    qq_r,
    qq_q,
    box,
    1,
    "xyz_quat_final",
)
=#

"""Calculates real, self and tinfoil contributions to Coulombic energy for
a single molecule interacting with the rest of the system"""
function EwaldShort(
    i::Int64,
    system::Requirements,
    ewald::EWALD,
    box::Float64,
    qq_r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
    tinfoil = false,
)
    partial_e = 0.0
    partial_v = 0.0
    overlap = false
    # Calculate new ewald REAL energy
    realEwald, overlap = EwaldReal(
        qq_r,
        qq_q,
        ewald.kappa,
        box,
        system.thisMol_theseAtoms,
        i,
        system,
    )
    realEwald *= ewald.factor
    partial_e += realEwald
    partial_v += (realEwald / 3)

    # Calculate Self interaction energy
    #selfEnergy = EwaldSelf(ewald, qq_q)
    #partial_e += selfEnergy
    #partial_v += (selfEnergy / 3)

    #Calculate tinfoil Boundaries
    #=
    if tinfoil
        tfb = TinfoilBoundary(system, ewald, qq_q, qq_r) * ewald.factor
        partial_e += tfb
        partial_e += (tfb / 3)
    end
    =#
    return partial_e, partial_v, overlap
end


""" Recipricol space Ewald energy parallel and molecular version"""
function RecipLong(
    system::Requirements,
    ewald::EWALD,
    r::Vector{SVector{3,Float64}},
    qq_q::Vector{Float64},
)
    #=
    Intended to calculate the recipricol for the molecule that was moved, not then
    entire system.
    =#
    #cfac = ewald.cfac
    kxyz = ewald.kxyz[:]
    energy = 0.0
    L = system.box
    twopi = 2.0 * π
    nk = ewald.nk
    n = length(r)
    k_sq_max = ewald.k_sq_max

    eikx = OffsetArray{Complex{Float64}}(undef, 1:n, 0:nk) #SArray{n,nk+2}([0.0 + 0.0*im for i=1:n for j=1:(nk+2) ]...)   #undef,n,nk+2)
    eiky = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk) #zeros(ComplexF64,n,2*nk+2) #SArray{n,2*nk+2}([0.0 + 0.0*im for i=1:n for j=1:(2*nk+2) ]...)  #(undef,n,2*nk+2)
    eikz = OffsetArray{Complex{Float64}}(undef, 1:n, -nk:nk)

    for j = 1:n
        # calculate kx,ky,kz =1,2 explicitly
        eikx[j, 0] = 1.0 + 0.0 * im
        eiky[j, 0] = 1.0 + 0.0 * im
        eikz[j, 0] = 1.0 + 0.0 * im

        eikx[j, 1] =
            cos(twopi * (r[j][1]) / L) + sin(twopi * (r[j][1]) / L) * im
        eiky[j, 1] =
            cos(twopi * (r[j][2]) / L) + sin(twopi * (r[j][2]) / L) * im
        eikz[j, 1] =
            cos(twopi * (r[j][3]) / L) + sin(twopi * (r[j][3]) / L) * im

        eiky[j, -1] = conj(eiky[j, 1])
        eikz[j, -1] = conj(eikz[j, 1])
    end
    # calculate remaining positive kx, ky, kz by recurrence
    for k = 2:(nk)
        for j = 1:n
            eikx[j, k] = eikx[j, k-1] * eikx[j, 1]
            eiky[j, k] = eiky[j, k-1] * eiky[j, 1]
            eikz[j, k] = eikz[j, k-1] * eikz[j, 1]

            #eikx[j,-k] = conj( eikx[j,k] )
            eiky[j, -k] = conj(eiky[j, k])
            eikz[j, -k] = conj(eikz[j, k])
        end
    end
    # TODO This can be sped up. It is generating too many temp Arrays
    # It should be faster if devectorized, but that makes even more temp arrays
    term = 0.0 + 0.0 * im
    for i = 1:length(ewald.cfac)
        term = 0.0 + 0.0 * im
        for l = 1:n
            term +=
                qq_q[l] *
                eikx[l, kxyz[i][1]] *
                eiky[l, kxyz[i][2]] *
                eikz[l, kxyz[i][3]]#* eiky[l,kxyz[i][2]] * eikz[l,kxyz[i][3]]
        end
        #term = sum( qq_q[:] .* eikx[:,kxyz[i][1]] .* eiky[:,kxyz[i][2]] .* eikz[:,kxyz[i][3]] )
        energy += ewald.cfac[i] * real(conj(term) * term)
        ewald.sumQExpNew[i] = term
        ewald.sumQExpOld[i] = term
    end
    return energy, ewald
end # function


#=
println(
    blk,
    " Energy: ",
    averages.energy / totProps.totalStepsTaken / nMol,
    " Ratio: ",
    totProps.numTranAccepted / totProps.totalStepsTaken,
    " Pressure: ",
    ρ * temperature + averages.virial / box^3 / totProps.totalStepsTaken, #  Pressure(total, ρ, temperature, box^3),
    " Pcut: ",
    pressure_lrc(ρ, system.r_cut),
    " Ecorr: ",
    potential_lrc(ρ, r_cut),
    "p delta: ",
    pressure_delta(ρ, system.r_cut),
    " A & T p_c: ",
    ρ * temperature +
    averages.virial / box^3 / totProps.totalStepsTaken +
    pressure_delta(ρ, system.r_cut),
    " instant energy: ",
    total.energy / nMol,
    " instant pressure: ",
    Pressure(total, ρ, temperature, box^3) +
    pressure_delta(ρ, system.r_cut),
    " overlap count: ", ovr_count
)
=#
