#=
ERROR: LoadError: MethodError: no method matching FindNumericAtomType(::Int64, ::Nothing, ::FFParameters)
Closest candidates are:
  FindNumericAtomType(::Int64, ::Int64, ::FFParameters) at C:\Users\Zarathustra\github\QuixoticSimulating\src\scripts\setup.jl:21

This typically occurs when the atom names are different in the topology file and
pdb for individual molecules.

i.e., in spce.pdb oxygen is given the name OW, but in topol.top
oxygen in water is called O or O1 or something other than OW.

Another common mistake I make is calling water SOL in one file around
WAT in another :S
=#
