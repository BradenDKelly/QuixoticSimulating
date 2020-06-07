#=
Pressure gets its own folder because of how nuanced pressure calculations can
be. Currently it is quite simple, but in the future there should be both
virial and thermodynamic routes, of which there are at least 2 different
thermodynamic routes. Also, pressure tensors would be great, particularly for
boxes that are not square.
end
=#
include("constants.jl")

""" Calculates pressure including the tail correction for LJ"""
function Pressure(vir, ρ, T, vol, r_cut)
    return ρ * T + vir.virial / vol + pressure_lrc(ρ, r_cut)
end

"""Calculates pressure without tail correction for LJ"""
function Pressure(vir, ρ, T, vol)
    return ρ * T + vir.virial / vol
end

function pressure_lrc(ρ, r_cut)
    """Calculates long-range correction for Lennard-Jones pressure."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * ((32.0 / 9.0) * sr3^3 - (16.0 / 3.0) * sr3) * ρ^2
end

function pressure_delta(ρ, r_cut)
    """Calculates correction for Lennard-Jones pressure due to discontinuity in the potential at r_cut."""
    # density, r_cut, and the results, are in LJ units where sigma = 1, epsilon = 1
    sr3 = 1.0 / r_cut^3
    return π * (8.0 / 3.0) * (sr3^3 - sr3) * ρ^2
end

"""Calculates pressure for LJ and electrostatic"""
function Pressure(lj_vir, qq_ener, num_mol::Int, vol::Float64, T::Float64, r_cut::Float64)
    # ideal gas contribution to pressure_lrc
    ig_contrib = (num_mol * kb * T) / vol   # kJ/Å³
    ig_contrib *= 1e30                      # kJ/m³ = kPa
    ig_contrib /= 100.0                     # bar
    # LJ contribution to pressure, included LRC
    lj_contrib = vir.virial / vol           # K/m³          #+ pressure_lrc(ρ, r_cut)
    lj_contrib /= R                         # kJ/m³
    lj_contrib /= 100.00                    # bar
    # electrostatic contribution to pressure
    qq_contib = qq_ener / 3 / vol           # K/m³
    qq_contib /= R                          # kJ/m³
    qq_contib /= 100.00                     # bar

    println("Pressure results: ig is $ig_contrib, lj is $lj_contrib, qq is $qq_contrib")
    return ig_contrib + lj_contrib + qq_contrib
end

"""Calculates pressure for LJ and electrostatic"""
function Pressure(total::Properties, num_mol::Int, vol::Float64, T::Float64)
    # ideal gas contribution to pressure_lrc
    ig_contrib = (num_mol * kb * T) / vol   # kJ/Å³
    ig_contrib *= 1e30                      # kJ/m³ = kPa
    ig_contrib /= 100.0                     # bar
    println("Pressure results: ig is $ig_contrib")
    # LJ contribution to pressure, included LRC
    lj_vir_contrib = total.lj_virial / vol           # K/Å³          #+ pressure_lrc(ρ, r_cut)
    lj_vir_contrib *= 1e30                      # K/m³ = kPa
    lj_vir_contrib *= R                         # kJ⋅mol^-1/m³ = kPa/mol
    lj_vir_contrib /= 100.00                    # bar/mol
    lj_vir_contrib /= na                        # bar
    println("lj contribution is $lj_vir_contrib")
    real_vir_contrib = total.real_virial / vol    # K/Å³          #+ pressure_lrc(ρ, r_cut)
    real_vir_contrib *= 1e30                      # K/m³ = kPa
    real_vir_contrib *= R                         # kJ⋅mol^-1/m³ = kPa/mol
    real_vir_contrib /= 100.00                    # bar/mol
    real_vir_contrib /= na                        # bar
    println("real contribution is $real_vir_contrib")
    recip_vir_contrib = total.recip_virial / vol    # K/Å³          #+ pressure_lrc(ρ, r_cut)
    recip_vir_contrib *= 1e30                      # K/m³ = kPa
    recip_vir_contrib *= R                         # kJ⋅mol^-1/m³ = kPa/mol
    recip_vir_contrib /= 100.00                    # bar/mol
    recip_vir_contrib /= na                        # bar
    println("recipricol contribution is $recip_vir_contrib")
    self_vir_contrib = total.self_virial / vol    # K/Å³          #+ pressure_lrc(ρ, r_cut)
    self_vir_contrib *= 1e30                      # K/m³ = kPa
    self_vir_contrib *= R                         # kJ⋅mol^-1/m³ = kPa/mol
    self_vir_contrib /= 100.00                    # bar/mol
    self_vir_contrib /= na                        # bar
    println("self contribution is $self_vir_contrib")
    intra_vir_contrib = total.intra_virial / vol    # K/Å³          #+ pressure_lrc(ρ, r_cut)
    intra_vir_contrib *= 1e30                      # K/m³ = kPa
    intra_vir_contrib *= R                         # kJ⋅mol^-1/m³ = kPa/mol
    intra_vir_contrib /= 100.00                    # bar/mol
    intra_vir_contrib /= na                        # bar
    println("intra contribution is $intra_vir_contrib")

    tot = ig_contrib + lj_vir_contrib + real_vir_contrib +
            recip_vir_contrib + self_vir_contrib + intra_vir_contrib
    tot_e = lj_vir_contrib + real_vir_contrib +
            recip_vir_contrib + self_vir_contrib + intra_vir_contrib
    println("Comparison total.vir $(total.virial) sum is $(tot_e*na*100/R/1e30*vol) ")
    println("Total pressure is $tot")
    return tot
end
