##########
#
#                         Energy Routines
#
###############
using Setfield
using Distributions
using BenchmarkTools

include("auxillary.jl")
include("ewalds.jl")

################################################################################
"""polyatomic LJ with mos and soa"""
function LJ_poly_ΔU(i, moa::StructArray, soa::StructArray,
                            vdwTable, r_cut, box)
    # Calculates Lennard-Jones energy of 1 particle, "i", interacting with the
    # other N-1 particles in the system
    # input:  i, Requirements (A struct with ϵ, σ, r_cut,, box and r)    #
    # output  energy, virial both scalars

    ri = moa.COM[i]   # molecule i
    rMol = moa.COM    # all molecules

    startAtom = moa.firstAtom[i] # first atom in molecule i
    endAtom = moa.lastAtom[i]    # lastatom in molecule i
    list_type = soa.atype
    moli_type = soa.atype[startAtom:endAtom]

    sFs, sFe = moa.firstAtom, moa.lastAtom

    ra = soa.coords[startAtom:endAtom]
    rb = soa.coords

    table = vdwTable
    box = box
    r_cut = r_cut
    diameter = 0 #r_cut * 0.25 + 5 #2.0 * sqrt( maximum( sum(db^2,dim=1) ) )
    rm_cut_box = (r_cut + diameter)       # Molecular cutoff in box=1 units
    rm_cut_box_sq = rm_cut_box^2          # squared
    r_cut_sq = r_cut^2                    # Potential cutoff squared in sigma=1 units

    rᵢⱼ = rab = SVector{3}(0.0, 0.0, 0.0)

    pot, vir = 0.0, 0.0

    # cycle through all molecules
     for (j, rj) in enumerate(rMol)
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
             for a = 1:(endAtom-startAtom+1)

                """Loop over all atoms in molecule B"""
                for b = sFs[j]:sFe[j]

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

"""Rotate each molecule X times, save the lowest energy rotation. Loop Through
All Molecules. Can do this loop Y times. """
function EnergyMinimize(thisSystem::Requirements, db, quatVec::Vector, qq_q, qq_r, ewald::EWALD)
    sys = deepcopy(thisSystem)
    qV = deepcopy(quatVec)
    nLoops = 5
    nMols = length(sys.rm)
    nTrials = 3
    ra = [SVector(0.0, 0.0, 0.0) for i = 1:3]
    println("EnergyMinimize in energy.jl, line 360, has some hardcoding")
    @inbounds for i = 1:nLoops
        for j = 1:nMols
            com = sys.rm[j]
            lowE, lowV = LJ_poly_ΔU(j, sys)
            partial_ewald_e, partial_ewald_v = EwaldShort(
                j, system,ewald, temperature, box, qq_r, qq_q, false )
            lowE += partial_ewald_e
            lowV += partial_ewald_v
            savedQuat = qV[j]
            for k = 1:nTrials
                ei = random_rotate_quaternion(0.05, savedQuat)
                ai = q_to_a(ei)
                for a = 1:at_per_mol # Loop over all atoms
                    ra[a] = com + MATMUL(ai, db[:, a])
                end # End loop over all atoms
                sys.ra[sys.thisMol_theseAtoms[j][1]:sys.thisMol_theseAtoms[j][2]] =
                    ra
                qq_r[sys.thisMol_theseAtoms[j][1]:sys.thisMol_theseAtoms[j][2]] =
                        ra
                newE, newV = LJ_poly_ΔU(j, sys)
                partial_ewald_e, partial_ewald_v = EwaldShort(
                    j, system,ewald, temperature, box, qq_r, qq_q, false )
                newE += partial_ewald_e
                newV += partial_ewald_v
                if newE < lowE
                    savedQuat = ei
                end
            end
            qV[j] = savedQuat
        end
        println("loop: ", i)
    end

    return deepcopy(qV)

end

function LennardJones(rij)
    return 4 * 1 * ((1 / rij)^12 - (1 / rij)^6)
end

# ================================================================================================================================
function press_corr(system::Requirements, num_atom_types = 2, b = [100, 200])
    # ================================================================================================================================

    # Adds lennard-jones or Buckingham (EXP-6) tail correction to pressure
    # ======================================================================================================
    #                 Correction to Each atom type
    # ======================================================================================================

    corp = 0.0
    corpi = 0.0
    Rc = system.r_cut
    vol = system.box^3
    for i = 1:num_atom_types # cycle through all of the atoms, this counts as being atom 1, cycle "k" is atom 2

        for j = 1:num_atom_types #for each atom, we must cycle through all of the molecules it interacts with

            #    			write(*,*) " coru pre-error error : ", i,coru, cori(i)
            epsij = system.table.ϵᵢⱼ[i, j]  #EPS( i, j )
            sigij = system.table.σᵢⱼ[i, j]  #SIG( i, j )

            # ============================================================================
            #					Lennard-Jones Force-Field
            # ============================================================================

            sig3 = (sigij * sigij * sigij)
            sigor3 = sig3 / (Rc * Rc * Rc)
            sigor9 = sigor3 * sigor3 * sigor3

            corp =
                corp +
                b[i] * b[j] * epsij * sig3 * ((2.0 / 3.0) * sigor9 - sigor3)
            #corpi(i) = corpi(i) + b(i)*b(j) * epsij * sig3 * ( (2.0/3.0) * sigor9 - sigor3 ) # individual molecule corrections (Not used, it is here to trick you)

        end

    end


    # ============================================================================
    #       Add on constants kept out of loop for speed purposes
    # ============================================================================

    #if(FF_Flag .eq. 1) then # use Lennard-Jones force-field

    corp = 16 * π / (3.0 * vol * vol) * corp
    #corpi = 16*π/(3.0*vol*vol) * corpi
    return corp
    #elseif(FF_Flag .eq. 2) then # use EXP-6 force-field

    #    corp  = 4.0 * 3.141592654 / (3.0 * vol * vol) * corp
    #    corpi = 4.0 * 3.141592654 / (3.0 * vol * vol) * corpi

    #end
end #Subroutine press_corr

# ================================================================================================================================
function ener_corr(system::Requirements, num_atom_types = 2, b = [100, 200])
    # ================================================================================================================================

    # -Adds lennard-jones or EXP-6 tail correction to potential energy.
    # -Works for polyatomic molecules. Calculates total potential for a pure
    # monatomic molecule as well as individual species tail corrections for
    # mixtures.

    # ============================================================
    #         Written by Braden Kelly a.k.a. Zarathustra sometime
    #                  prior to your reading this.
    # ============================================================

    coru = 0.0
    Rc = system.r_cut
    vol = system.box^3
    #println(system.table.ϵᵢⱼ)
    for i = 1:num_atom_types # cycle through all of the atoms, this counts as being atom 1, cycle "k" is atom 2

        for j = 1:num_atom_types #for each atom, we must cycle through all of the molecules it interacts with

            #    			write(*,*) " coru pre-error error : ", i,coru, cori(i)
            epsij = system.table.ϵᵢⱼ[i, j]  #EPS( i, j )
            sigij = system.table.σᵢⱼ[i, j]  #SIG( i, j )

            # ============================================================================
            #					Lennard-Jones Force-Field
            # ============================================================================


            sig3 = (sigij * sigij * sigij)
            sigor3 = sig3 / (Rc * Rc * Rc)
            sigor9 = sigor3 * sigor3 * sigor3

            coru += b[i] * b[j] * epsij * sig3 * ((1.0 / 3.0) * sigor9 - sigor3)

        end

    end

    coru = 8.0 * π / (3.0 * vol) * coru
    print("ener_corr value: ", coru, "vol: ", vol)
    return coru

end #Subroutine ener_corr

"""Calculate total potential energy using moa and soa for WolfSummation"""
function potential(
    moa::StructArray,
    soa::StructArray,
    tot::Properties,
    ewald::EWALD,
    vdwTable::Tables,
    sim_props::Properties2#triggers wolf summations using double strings
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = sim_props.LJ_rcut
    #@assert r_cut == 10.0
    box = sim_props.box
    total.energy = 0.0
    total.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0

    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing

    for i = 1:length(moa.COM)
        ener, vir = LJ_poly_ΔU(i, moa, soa, vdwTable, r_cut, box)
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
    for i = 1:length(moa.COM)
        ener, overlap =
            #EwaldReal(qq_r, qq_q,ewald.kappa, box,thisMol_thisAtom, i, system)
             EwaldReal(i, moa, soa, ewald,totProps.qq_rcut, box)
        totReal += ener
        if overlap
            println("overlap after EwaldReal")
            #exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb += totReal
    reall = totReal
    println("Coulomb contribution is: ", reall)

    prefactor = 0.0
    for i=1:length(soa.charge)
        for j=1:length(soa.charge)
            prefactor += soa.charge[i]*soa.charge[j] * erfc(ewald.kappa*r_cut) / r_cut
        end
    end
    prefactor *= -1
    prefactor2 = (erfc(ewald.kappa*r_cut) / 2 / r_cut +ewald.kappa / sqrt(π) ) *
                    dot(soa.charge,soa.charge)
    tot.energy += (prefactor - prefactor2)*ewald.factor
    tot.coulomb += (prefactor - prefactor2)*ewald.factor

    println("total energy is: ", tot.energy)
    println("total coulomb is: ", tot.coulomb)
    println("total LJ energy is: ", tot.energy - tot.coulomb)


    return tot

end

"""Calculate total potential energy using moa and soa for Ewald Summation"""
function potential(
    moa::StructArray,
    soa::StructArray,
    tot::Properties,
    ewalds::EWALD,
    vdwTable::Tables,
    sim_props::Properties2,
    coulomb_style::String #triggers wolf summations using double strings
)

    #tot = Properties(0.0,0.0)
    ener, vir = 0.0, 0.0
    r_cut = sim_props.LJ_rcut
    #@assert r_cut == 10.0
    box = sim_props.box
    tot.energy = 0.0
    tot.virial = 0.0
    #########
    #
    #     Calculate LJ
    #
    ##################
    LJ, real = 0.0, 0.0
    #ener, vir = LJ_poly_ΔU(2, system) # turn on for timing
    println("tot_vir is $(tot.virial) remainder is $(tot.lj_virial)")
    for i = 1:length(moa.COM)
        ener, vir = LJ_poly_ΔU(i, moa, soa, vdwTable, r_cut, box)
        tot.energy += ener
        tot.virial += vir
        LJ += ener
    end
    tot.energy = tot.energy / 2
    tot.virial = tot.virial / 2
    tot.lj_virial = tot.virial
    LJ = LJ / 2
    println("Total LJ energy is: ", tot.energy)
    println("LJ tot_vir is $(tot.virial) remainder is $(tot.lj_virial)")
    #########
    #
    #     Calculate EWALD
    #
    ##################

    # Real
    totReal = 0.0
    for i = 1:length(moa.COM)
        ener, overlap =
            #EwaldReal(qq_r, qq_q,ewald.kappa, box,thisMol_thisAtom, i, system)
             EwaldReal(i, moa, soa, ewald,totProps.qq_rcut, box)
        totReal += ener
        if overlap
            println("overlap after EwaldReal")
            #exit()
        end
    end
    totReal *= ewald.factor / 2
    tot.energy += totReal   # divide by 2 to account for double counting
    tot.coulomb = totReal
    tot.virial += totReal / 3.0
    tot.real_virial = totReal / 3.0
    reall = totReal
    println("Coulomb contribution is: ", reall)
    println("real tot_vir is $(tot.virial) remainder is $(tot.lj_virial + tot.real_virial)")
    recipEnergy, ewalds = RecipLong(ewalds, soa.coords, soa.charge, box)  #RecipLong(system, ewald, qq_r, qq_q, kfacs)
    recipEnergy *= ewalds.factor


    println("Total recipricol Ewald is: ", recipEnergy)
    tot.energy += recipEnergy   # divide by 2 to account for double counting
    tot.coulomb += recipEnergy
    tot.virial += recipEnergy / 3.0
    tot.recip_virial = recipEnergy / 3.0
    println("recip tot_vir is $(tot.virial) remainder is $(tot.lj_virial + tot.real_virial+tot.recip_virial)")
    selfEnergy = EwaldSelf(ewalds, soa.charge)
    println("Self energy: ", selfEnergy)
    tot.energy += selfEnergy
    tot.coulomb += selfEnergy
    tot.virial += selfEnergy / 3.0
    tot.self_virial = selfEnergy / 3.0

    intraEnergy = EwaldIntra(ewald, soa, moa)
    println("Intra energy: ", intraEnergy)
    tot.energy += intraEnergy
    tot.coulomb += intraEnergy
    tot.virial += intraEnergy / 3.0
    tot.intra_virial = intraEnergy / 3.0



    println("total energy is: ", tot.energy)
    println("total coulomb is: ", tot.coulomb)
    println("total LJ energy is: ", tot.energy - tot.coulomb)


    return tot

end
