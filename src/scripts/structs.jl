macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

@def aij begin
  ai::Int64
  aj::Int64

end

@def aijk begin
    @aij
    ak::Int64
end

@def aijkl begin
    @aijk
    al::Int64
end

@def xyz begin
  x::Float64
  y::Float64
  z::Float64
end

"""Struct for FF parameters (currently LJ and EXP6) """
mutable struct Tables #{T<:Vector} #<: ForceField
    # passed two 1D arrays, convert them both to 2D matrices and
    # apply geometric and arithmetic mixing rules
    ϵᵢⱼ::Array{Float64,2}
    σᵢⱼ::Array{Float64,2}
    function Tables(a::Vector{T}, b::Vector{T}) where {T}
        e = [sqrt(a[i] * a[j]) for i = 1:length(a), j = 1:length(a)]
        s = [(b[i] + b[j]) / 2 for i = 1:length(b), j = 1:length(b)]
        new(e, s)
    end
end

"""Struct with energies, old and new, short and long"""
mutable struct Properties
    energy::Float64
    virial::Float64
    lj_virial::Float64
    real_virial::Float64
    recip_virial::Float64
    self_virial::Float64
    intra_virial::Float64
    coulomb::Float64
    recip::Float64
    recipOld::Float64
    old_e::Float64
    old_v::Float64
end

"""Struct for tracking and optimizing translation moves"""
mutable struct Moves
    naccepp::Int           # previous number of attempts
    naccept::Int           # current number of attempts
    attempp::Int
    attempt::Int
    set_value::Float64     # desired success rate
    d_max::Float64
end

"""structure for energy calculation properties """

mutable struct Requirements
    rm::Vector{SVector{3,Float64}}
    ra::Vector{SVector{3,Float64}}
    nMols::Int
    nAtoms::Int
    nCharges::Int
    thisMol_theseAtoms::Vector{SVector{2,Int64}}
    molNames::Vector     # strings
    molTypes::Vector     # integers
    atomNames::Vector    # strings
    atomTypes::Vector    # integers
    #ϵ::Array{Float64,2}
    #σ::Array{Float64,2}
    table::Tables
    box::Float64
    r_cut::Float64
end

""" structure for generic simulation properties"""
mutable struct Properties2
    temperature::Float64
    ρ::Float64
    pressure::Float64
    dr_max::Float64
    dϕ_max::Float64
    move_accept::Float64
    numTranAccepted::Int
    totalStepsTaken::Int
    quat::Vector{SVector{4,Float64}}
    LJ_rcut::Float64
    qq_rcut::Float64
    box::Float64
end

#include("constraints.jl")
abstract type ForceField end
abstract type Gromacs <: ForceField end

abstract type GAFF <: ForceField end


abstract type AngleType <: ForceField end
abstract type DihedralType <: ForceField end

abstract type AtomInfo end


abstract type TypeStructArray end
#analyze_malloc(raw"C:\Users\Zarathustra\Documents\JuliaScripts")

################################################################################
#
#                Topology Data Structures (Objects)
#
################################################################################
#=; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
      1             2               yes             0.5       0.8333 =#
struct Defaults <: Gromacs
    nbfunc::Int64
    combRule::Int64
    genPairs::String
    fudgeLJ::Float64
    fudgeQQ::Float64
end
# Constructor for Defaults
Defaults( nbfunc::Int64,combRule::Int64,genPairs::String,fudgeLJ::Float64,
    fudgeQQ::Float64) = Defaults(nbfunc,combRule,genPairs,fudgeLJ,fudgeQQ)

#=
Gromacs Style
; name   atomicnr    mass      charge     ptype  sigma(nm)     epsilon(kJ/mol)
       O1        O1    15.99940   -0.8340      A   0.315061     0.6364000
=#
struct Atomtypes <: Gromacs
    name::String # = line[0]
    atomicnr::Any # Can be a number or a string
    mass::Float64 # = line[2]   # over rode by Atoms
    charge::Float64 # = line[3] # over rode by Atoms/ Not worth including
    ptype::String # = line[4]
    σ::Float64 # = line[5]
    ϵ::Float64 # = line[6]

end

#= ;   nr   type  resnr residue  atom   cgnr     charge       mass
          1  O1      1    SOL      O1      1      -0.8340        15.99940
=#
mutable struct Atoms <: Gromacs # mutable since charge may get updated
    nr::Int64 #  = int( line[0] )
    type::String # = line[1]
    resnr::Int64 # =int( line[2] )
    resnm::String # = line[3]
    atomnm::String # = line[4]
    cgnr::Int64 # = int( line[5] )
    charge::Float64 # = float( line[6] )    # overwritten if also in itp file
    mass::Float64 # = float( line[7] )

end
################################################################################
#
#                                BONDS
#
################################################################################
#=;   ai     aj funct   r             k
      1      2   1    1.0860e-01    2.8937e+05 ;     C1 - H1 =#
mutable struct Bonds <: GAFF
    @aij
    funct::Int64 # = int( line[2] )
    bondLength::Float64 # = float( line[3] )
    kparam::Float64 # = float( line[4] )
end

################################################################################
#
#                                Pairs
#
################################################################################
#=;   ai     aj    funct
        1      6      1  =#
struct Pairs <: Gromacs
    @aij
    funct::Int64
end
################################################################################
#
#                                ANGLES
#
################################################################################
#=;   ai     aj     ak    funct   theta         cth
 1      3      4      1    1.1988e+02    4.0317e+02  =#
struct GAFFAngles <: AngleType
    funct::Int64 # = int(line[3])
    theta::Float64 # = float(line[4])
    kparam::Float64 # = float(line[5])
end
struct OPLSAngles <: AngleType
    #= Insert demo =#
end

struct CHARMMAngles <: AngleType
    #= Insert demo =#
end

struct TrAPPEAngles <: AngleType
    #= Insert demo =#
end

mutable struct Angles{T<:AngleType}
    # parametric type storing generic atom indices and specific FF parameters
    @aijk
    params::T
end
# constructor for OPLS angles (incomplete)
Angles(ai::Int64,aj::Int64,ak::Int64)=Angles( ai,aj,ak,OPLSAngles() )
# constructor for GAFF angles
Angles(ai::Int64,aj::Int64,ak::Int64,funct::Int64,theta::Float64,
    kparam::Float64)=Angles( ai,aj,ak,GAFFAngles(funct,theta,kparam) )
# constructor for CHARMM angles (incomplete)
Angles(ai::Int64,aj::Int64,ak::Int64)=Angles( ai,aj,ak,CHARMMAngles() )
# constructor for TrAPPE angles (incomplete)
Angles(ai::Int64,aj::Int64,ak::Int64)=Angles( ai,aj,ak,TrAPPEAngles() )
# constructor for AMBER angles (incomplete)
Angles(ai::Int64,aj::Int64,ak::Int64)=Angles( ai,aj,ak,AMBERAngles() )
################################################################################
#
#                                DIHEDRALS
#
################################################################################
#= ; treated as RBs in GROMACS to use combine multiple AMBER torsions per quartet
;    i j k l func   C0   C1     C2     C3   C4    C5
     1 3 5 6  3  30.334 0.00 -30.334  0.00 0.00 0.000 ;=#
struct Improper <: DihedralType
    kparams::Vector{Float64}
end
#=; treated as propers in GROMACS to use correct AMBER analytical function
   ;    i  j  k  l  func phase     kd   pn
        1  5  3  4   1   180.0  4.6024   2  =#
struct Proper <: DihedralType
    # parametric type, only include here things specific to proper dihedrals
    phase::Float64
    kd::Float64
    pn::Int64
end

mutable struct Dihedrals{T<:DihedralType}
    # parametric abstract type, includes all general info plus specific{T}
    @aijkl
    funct::Int64
    params::T
end
# Improper Dihedral constructor (GAFF)
Dihedrals(ai::Int64,aj::Int64,ak::Int64,al::Int64,funct::Int64,
        kparam::Vector{Float64}) = Dihedrals(ai,aj,ak,al,funct,Improper(kparam))
# Proper Dihedral constructor (GAFF)
Dihedrals(ai::Int64,aj::Int64,ak::Int64,al::Int64,funct::Int64, phase::Float64,
        kd::Float64,pn::Int64) = Dihedrals(ai,aj,ak,al,funct,Proper(phase,kd,pn))

"""Struct holding molecular parameters for atoms/bonds/angles/dihedrals"""
struct MolParam <: ForceField
    name::AbstractString
    nrexcl::Int64
    atoms::Vector{Atoms}
    bonds::Vector{Bonds}
    angles::Vector{Angles}
    dihedrals::Vector{Dihedrals}
    #MolParam() = new()           # avoids needing to instantiate all variables initially
end

"""Struct holding all force-field parameters such as atom types, FF etc."""
struct FFParameters <: Gromacs
    defaults::Defaults # contains mixing rule, fudge factors
    atomTypes::Vector{Atomtypes} # one for each atom type
    molParams::Vector{MolParam}  # one object for each molecule
    system::AbstractString
    molecules::Dict{String,Int64} # name of molecule and number of molecules in system
    #molName2Num::Dict{String,Int64}
    #FFParameters()=new() # dummy initializer, necessary so I can define topology = FFParameters()
end
mutable struct XYZ <: FieldVector{3, Float64}
    @xyz
end

"""Struct holding molecule type specific information for rotations"""
struct BodyFixed{T}
    r::Vector{SVector{3,T}}
    mass::Vector{T}
    atype::Vector{Int64}
end

"""Function for generating body_fixed coordinates of each molecule"""
function BodyFixed(m,systemTop)
    # input: struct of type Topology (moleculeList)
    # output: struct of type BodyFixed

    # ASSUMPTIONS: Atom Centered Charges
    natoms = length(m.r)             # number of atoms in molecule
    mass = []
    atype= []
    for i=1:natoms
        mol = FindMolType(String(m.resnm[i]),moleculeList )  # sort of rhetorical
        atom = FindAtomInMol( String(m.atomnm[i]),mol )
        atomTypeNumber = FindNumericAtomType(mol,atom,systemTop)
        push!(mass, systemTop.molParams[mol].atoms[atom].mass)
        push!(atype,atomTypeNumber)
    end



    return BodyFixed( [m.r[i] for i=1:natoms],
                      [mass[i] for i=1:natoms ],
                      [atype[i] for i=1:natoms] )

end # BodyFixed

"pdb information, could be an entire system or one molecule"
struct Topology
    name::AbstractString
    box::Vector{Float64}
    r::Vector{SVector{3,Float64}}
    atomnm::Vector{AbstractString}
    resnm::Vector{AbstractString}
    resnr::Vector{Int64}
    elem::Vector{AbstractString}
    #conect::Vector{Tuples}
end

"""Struct for MD with information for 1 atom"""
struct PerAtomStruct <: AtomInfo
    atomID::String
    atomType::Int64
    molID::String
    molNumber::Int64
    r::XYZ                # coordinates
    v::XYZ                # velocities
    f::XYZ                # forces
    mass::Float64         # mass
    qq::Float64           # charge
    #Neighborlist::Vector{Int64}   # indices of neighbors
end

"""Molecular dynamics struct for positions, velocity, forces etc."""
struct ParticleAtom{T,P} <: TypeStructArray
    r::SVector{3,T}
    v::SVector{3,T}
    f::SVector{3,T}
    at::P
    mt::P
    mass::T
    qq::T
end

"""Monte Carlo Atomic properties like number, mass, coords, charge, molecule_parent"""
struct ParticleAtomKMC <: TypeStructArray
    molNum::Int64
    molType::Int64
    atype::Int64
    mass::Float64
    coords::SVector{3,Float64}
    charge::Float64
end

"""Struct with intramolecular parameters for bonds/angles/dihedrals"""
struct IntraForceField <: ForceField
    bonds::Vector{Bonds}          # Int64 Int64 bondtype
    angles::Vector{Angles}        # Int64 Int64 Int64 angletype
    dihedrals::Vector{Dihedrals}  # Int64 Int64 Int64 Int64 dihedraltype
end


""" Struct with PVT properties for the total system"""
mutable struct Properties22{F}
    volume::F
    temperature::F
    density::F
    pressure::F
    KE::F
    PE::F
end

"""Struct with atomic and orientational parameters for molecules"""
struct Molecule
    firstAtom::Int64
    lastAtom::Int64
    COM::SVector{3,Float64}
    quat::SVector{4,Float64}
    weight::Float64
    molType::Int64
end

"""Struct for KMC for numerical overflow limits"""
struct OverFlowTable{T}
    vdw::Array{T,2}
    qq::Array{T,2}
end

"""Struct for holding integer total quantities of atoms/molecules/charges etc."""
struct Numbers{I}
    atoms::I
    molecules::I
    atomTypes::I
    molTypes::I
    charges::I
end

"""
function Tables(a::Vector{T}, b::Vector{T}) where T
    e = [sqrt(a[i]*a[j]) for i=1:length(a), j=1:length(a)]
    s = [(b[i] + b[j]) / 2 for i=1:length(b), j=1:length(b)]
    return Tables(e,s)
end
"""
