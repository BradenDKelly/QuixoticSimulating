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

    tot = ig_contrib + total.virial *1e30 * R / 100 / na / vol
    println("Total pressure is $tot")
    return tot
end


##########################################################
##########################################################
#
#             Pressure Tensor
#
##########################################################
##########################################################

SUBROUTINE Compute_System_Total_Force(this_box)

   !****************************************************************************
   ! The subroutine calculates the total forces of a given box.
   ! The identity of the box is passed to the routine.
   ! The forces are then used to compute the pressure tensor.
   !
   ! CALLS
   !
   ! CALLED BY
   !
   ! Volume_Change
   ! Main
   ! Write_Properties_Buffer
   !
   !****************************************************************************

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: this_box

   !----------------------------------------------------------------------------

   INTEGER ::  is, im_1, im_2, is_1, is_2, this_im_1, this_im_2

   REAL(DP) :: rcom, rx, ry, rz, w_lrc

   REAL(DP),DIMENSION(3,3) :: tv_pair, tc_pair, w_inter_vdw, w_inter_charge

   LOGICAL :: get_interaction

   W_tensor_vdw(:,:,this_box) = 0.0_DP
   W_tensor_charge(:,:,this_box) = 0.0_DP
   W_tensor_recip(:,:,this_box) = 0.0_DP
   W_tensor_elec(:,:,this_box) =  0.0_DP

   DO is = 1, nspecies
      imLOOP1: DO im_1 = 1, nmols(is,this_box)
         this_im_1 = locate(im_1,is,this_box)
         IF (.NOT. molecule_list(this_im_1,is)%live) CYCLE imLOOP1

         w_inter_vdw(:,:) = 0.0_DP
         w_inter_charge(:,:) = 0.0_DP

         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP SCHEDULE(DYNAMIC) &
         !$OMP PRIVATE(im_2, this_im_2, get_interaction) &
         !$OMP PRIVATE(rcom, rx, ry, rz, tv_pair, tc_pair) &
         !$OMP REDUCTION(+:w_inter_vdw, w_inter_charge)
         imLOOP2: DO im_2 = im_1 + 1, nmols(is,this_box)
            this_im_2 = locate(im_2,is,this_box)
            IF (.NOT. molecule_list(this_im_2,is)%live) CYCLE imLOOP2

            CALL Check_MoleculePair_Cutoff(this_im_1,is,this_im_2,is,get_interaction, &
                                   rcom,rx,ry,rz)

            IF (.NOT. Get_Interaction) CYCLE imLOOP2

            CALL Compute_MoleculePair_Force(this_im_1,is,this_im_2,is, &
                   this_box,tv_pair,tc_pair,rx,ry,rz)

            w_inter_vdw(:,:) = w_inter_vdw(:,:) + tv_pair(:,:)
            w_inter_charge(:,:) = w_inter_charge(:,:) + tc_pair(:,:)

         END DO imLOOP2
         !$OMP END PARALLEL DO

         W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + w_inter_vdw(:,:)
         W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + w_inter_charge(:,:)

      END DO imLOOP1
   END DO

   DO is_1 = 1, nspecies
      imLOOP3: DO im_1 = 1, nmols(is_1,this_box)
         this_im_1 = locate(im_1,is_1,this_box)
         IF (.NOT. molecule_list(this_im_1,is_1)%live) CYCLE imLOOP3

         DO is_2 = is_1 + 1, nspecies

            w_inter_vdw(:,:) = 0.0_DP
            w_inter_charge(:,:) = 0.0_DP

            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP SCHEDULE(DYNAMIC) &
            !$OMP PRIVATE(im_2, this_im_2, get_interaction) &
            !$OMP PRIVATE(rcom, rx, ry, rz, tv_pair, tc_pair) &
            !$OMP REDUCTION(+:w_inter_vdw,w_inter_charge)
            imLOOP4: DO im_2 = 1, nmols(is_2,this_box)
               this_im_2 = locate(im_2,is_2,this_box)
               IF (.NOT. molecule_list(this_im_2,is_2)%live) CYCLE imLOOP4

               ! Check to see if the interaction needs to be computed between the molecules
               CALL Check_MoleculePair_Cutoff(this_im_1,is_1,this_im_2,is_2,get_interaction,rcom,rx,ry,rz)

               IF (.NOT. get_interaction ) CYCLE imLOOP4

               CALL Compute_MoleculePair_Force(this_im_1,is_1,this_im_2,is_2,this_box,tv_pair,tc_pair,rx,ry,rz)

               !                W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + tv_pair(:,:)
!                W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + tc_pair(:,:)

               w_inter_vdw(:,:) = w_inter_vdw(:,:) + tv_pair(:,:)
               w_inter_charge(:,:) = w_inter_charge(:,:) + tc_pair(:,:)

            END DO imLOOP4

            W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + w_inter_vdw(:,:)
            W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + w_inter_charge(:,:)
         END DO

      END DO imLOOP3
    END DO

    if (int_charge_sum_style(this_box) == charge_ewald)

       CALL Compute_System_Ewald_Reciprocal_Force(this_box)
       W_tensor_elec(:,:,this_box) =  W_tensor_recip(:,:,this_box)
   end

    if (int_vdw_sum_style(this_box) == vdw_cut_tail)

       CALL Compute_LR_Force(this_box,w_lrc)
       virial(this_box)%lrc = w_lrc

   end

    W_tensor_elec(:,:,this_box) = (W_tensor_elec(:,:,this_box) + W_tensor_charge(:,:,this_box)) * charge_factor
    W_tensor_total(:,:,this_box) = W_tensor_vdw(:,:,this_box) + W_tensor_elec(:,:,this_box)

end Compute_System_Total_Force
#-----------------------------------------------------------------------------



function Compute_MoleculePair_Force(im,is,jm,js,this_box,tens_vdw,tens_charge,rabx,raby,rabz)
#=
        Borrowed from Cassandra (they use globals, take care with ins/outs)

    INTEGER, INTENT(IN) :: im, is, jm, js, this_box
    !---------------------------------------------------------------------------

    INTEGER :: ia, ja

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rijsq, wij_vdw, wij_qq
    REAL(DP) :: rabx, raby, rabz
    REAL(DP),DIMENSION(3,3) :: tens_vdw, tens_charge

    REAL(DP) :: ffc, wxy, wxz, wyz

    LOGICAL :: get_vdw, get_qq
=#
    #tens_vdw(:,:) = 0.0_DP
    #tens_charge(:,:) = 0.0_DP
    tens_vdw = zeros(3,3)
    tens_charge = zeros(3,3)

    for ia = 1: natoms(is)
       for ja = 1: natoms(js)

          # Obtain the minimum image separation
          rxijp = atom_list(ia,im,is)%rxp - atom_list(ja,jm,js)%rxp
          ryijp = atom_list(ia,im,is)%ryp - atom_list(ja,jm,js)%ryp
          rzijp = atom_list(ia,im,is)%rzp - atom_list(ja,jm,js)%rzp

          # Now get the minimum image separation
          CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp, &
                  rxij,ryij,rzij)

          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          # Now figure out what needs to be computed, then call pair_energy
          CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)

          # Compute vdw and q-q energy using if required
          if (get_vdw || get_qq)
             # calculate force between two atoms. return virial for lj and real
             Compute_AtomPair_Force(rijsq,is,im,ia,js,jm,ja,
                  get_vdw,get_qq,Wij_vdw,Wij_qq)

             ffc = Wij_vdw/rijsq

             wxy = ffc*(0.5*(rxij*raby+ryij*rabx))
             wxz = ffc*(0.5*(rxij*rabz+rzij*rabx))
             wyz = ffc*(0.5*(ryij*rabz+rzij*raby))

             tens_vdw[1,1] += ffc*rxij*rabx
             tens_vdw[1,2] += wxy
             tens_vdw[1,3] += wxz
             tens_vdw[2,1] += wxy
             tens_vdw[2,2] += ffc*ryij*raby
             tens_vdw[2,3] += wyz
             tens_vdw[3,1] += wxz
             tens_vdw[3,2] += wyz
             tens_vdw[3,3] += ffc*rzij*rabz

             ffc = Wij_qq/rijsq

             wxy = ffc*(0.5*(rxij*raby+ryij*rabx))
             wxz = ffc*(0.5*(rxij*rabz+rzij*rabx))
             wyz = ffc*(0.5*(ryij*rabz+rzij*raby))

             tens_charge[1,1] += ffc*rxij*rabx
             tens_charge[1,2] += wxy
             tens_charge[1,3] += wxz
             tens_charge[2,1] += wxy
             tens_charge[2,2] += ffc*ryij*raby
             tens_charge[2,3] += wyz
             tens_charge[3,1] += wxz
             tens_charge[3,2] += wyz
             tens_charge[3,3] += ffc*rzij*rabz

         end # if

     end

 end

end #Compute_MoleculePair_Force

    #!-----------------------------------------------------------------------------
function Compute_AtomPair_Force
       (rijsq,is,im,ia,js,jm,ja,get_vdw,get_qq,Wij_vdw,Wij_qq)
#=
    ! LJ potential:  Wij = -rij/3 * d Eij / d rij.
    ! Use the virial in: P = NkT + < W >

    ! Computes the vdw and q-q pair force between atoms ia and ja of molecules
    ! im and jm and species is and js, given their separation rijsq. I have
    ! passed each component of separation but right now this is unnecessary.
    ! It also computes the real space part of the Ewald sum if necessary.

    ! Called by: Compute_System_Total_Force
    !-----------------------------------------------------------------------------
    ! Passed to
    REAL(DP) :: rxij,ryij,rzij,rijsq
    INTEGER :: is,im,ia,js,jm,ja
    LOGICAL :: get_vdw,get_qq

    ! Returned
    REAL(DP) :: Wij_vdw,Wij_qq

    ! Local
    INTEGER :: ibox
    ! LJ potential
    INTEGER :: itype, jtype
    REAL(DP) :: eps, sig, Eij_vdw
    REAL(DP) :: SigByR2,SigByR6,SigByR12
    REAL(DP) :: SigByR2_shift,SigByR6_shift,SigByR12_shift
    REAL(DP) :: roffsq_rijsq, roffsq_rijsq_sq, factor2, fscale
    ! Mie potential
    REAL(DP) :: SigByR, SigByRn, SigByRm, mie_coeff, mie_n, mie_m
    ! Coulomb potential
    REAL(DP) :: qi, qj, erfc_val, prefactor
    REAL(DP) :: rij, ewald_constant, exp_const

    Wij_vdw = 0.0_DP
    Wij_qq = 0.0_DP
    !-----------------------------------------------------------------------------
=#
    ibox = molecule_list(im,is)%which_box

    # If either atom is not yet present, then don't try to compute an energy
    #ExistCheck: &
    if (atom_list(ia,im,is)%exist && atom_list(ja,jm,js)%exist)

       # Determine atom type indices
       itype = nonbond_list(ia,is)%atom_type_number
       jtype = nonbond_list(ja,js)%atom_type_number

       #VDW_calc: &
       if (get_vdw .AND. itype /= 0 .AND. jtype /=0)

         if (int_vdw_style(ibox) == vdw_lj)
           # For now, assume all interactions are the same.
           # Use the lookup table created in Compute_Nonbond_Table
           eps = vdw_param1_table(itype,jtype)
           sig = vdw_param2_table(itype,jtype)

           SigByR2 = (sig^2.0) / rijsq
           SigByR6  = SigByR2 * SigByR2 * SigByR2
           SigByR12 = SigByR6 * SigByR6

           # Default potential for vdw_cut, vdw_cut_tail, vdw_cut_shift
           Wij_vdw = (24.0 * eps) * (2.0*SigByR12 - SigByR6)

           if (int_vdw_sum_style(ibox) == vdw_cut_switch)
             if (rijsq > ron_switch_sq(ibox) &&
                 rijsq <= roff_switch_sq(ibox))
               roffsq_rijsq = roff_switch_sq(ibox) - rijsq
               roffsq_rijsq_sq = roffsq_rijsq * roffsq_rijsq
               factor2 = switch_factor2(ibox) + 2.0_DP * rijsq
               fscale = roffsq_rijsq_sq * factor2 * switch_factor1(ibox)
               Eij_vdw = 4.0_DP * eps * (SigByR12 - SigByR6)
               Eij_vdw = fscale * Eij_vdw
               Wij_vdw = fscale / 3.0_DP * Wij_vdw
               Wij_vdw = Wij_vdw + 8.0_DP * rijsq * rijsq * roffsq_rijsq &
                       * Eij_vdw * switch_factor1(ibox) / 3.0_DP
             elseif (rijsq > roff_switch_sq(ibox))
               Wij_vdw = 0.0_DP
           end # if
        elseif (int_vdw_sum_style(ibox) == vdw_charmm)
             # Use the CHARMM LJ potential
             Wij_vdw = (12.0_DP * eps) * (SigByR12 - SigByR6)
         end # if
        elseif (int_vdw_style(ibox) == vdw_mie)
           eps = vdw_param1_table(itype,jtype)
           sig = vdw_param2_table(itype,jtype)
           mie_n = vdw_param3_table(itype,jtype) # repulsive exponent
           mie_m = vdw_param4_table(itype,jtype) # dispersive exponent
           rij = sqrt(rijsq)

           mie_coeff = mie_n/(mie_n-mie_m)*(mie_n/mie_m)**(mie_m/(mie_n-mie_m))
           SigByR = sig/rij
           SigByRn = SigByR ** mie_n
           SigByRm = SigByR ** mie_m
           Wij_vdw = (mie_coeff * eps) *(mie_n * SigByRn - mie_m * SigByRm)

         # Add other potential types here
     end

    end

       #qq_calc:
       if (get_qq)

         qi = nonbond_list(ia,is)%charge
         qj = nonbond_list(ja,js)%charge

         rij = sqrt(rijsq)
         prefactor = qi * qj / rij
         if (int_charge_sum_style(ibox) == charge_ewald)
           ewald_constant = 2.0 * alpha_ewald(ibox) / rootPI
           exp_const = DEXP(-alpha_ewald(ibox)*alpha_ewald(ibox)*rijsq)
           # May need to protect against very small rij
           erfc_val = erfc(alpha_ewald(ibox) * rij)
           Wij_qq = ( prefactor * erfc_val &
                  + qi * qj * ewald_constant * exp_const )

         else (int_charge_sum_style(ibox) == charge_dsf)

           Wij_qq = erfc(alpha_dsf(ibox)*rij)/(rijsq) + &
                    2.0 * alpha_dsf(ibox)/rootPI * &
                    DEXP(-alpha_dsf(ibox)*alpha_dsf(ibox) * rijsq) / rij - &
                    dsf_factor2(ibox)
           Wij_qq = qi*qj*Wij_qq*rij


       elseif (int_charge_sum_style(ibox) == charge_cut)
           Wij_qq = prefactor * charge_factor
       end

     end # qq_calc

    end # ExistCheck
end

function Compute_System_Ewald_Reciprocal_Force(this_box)
#=
  !***************************************************************************
  ! This subroutine computes the long range forces due to electrostatics
  !
  ! Based on APSS code reciprocal_ewald.f90
  !
  ! Added by Tom Rosch on 06/11/09
  ! (See Wheeler, Mol. Phys. 1997 Vol. 92 pg. 55)
  !
  !***************************************************************************

  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

!    !$ include 'omp_lib.h'

  INTEGER :: i, is, im, ia, this_locate, this_box

  REAL(DP) :: charge
  REAL(DP) :: qw(9), qwxy, qwxz, qwyz, un, const_val
  REAL(DP) :: xcmi, ycmi, zcmi, piix, piiy, piiz, arg, factor
  REAL(DP) :: recip_11, recip_21, recip_31, recip_22, recip_23, recip_33
 =#

  const_val = 1.0_DP/(2.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))
  qw(:) = 0.0_DP
  W_tensor_recip(:,:,this_box) = 0.0_DP

  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP SCHEDULE(STATIC) &
  !$OMP PRIVATE(i, un, qwxy, qwxz, qwyz) &
  !$OMP REDUCTION(+:qw)
  DO i = 1, nvecs(this_box)

     un = Cn(i,this_box) * (cos_sum(i,this_box) * cos_sum(i,this_box) + sin_sum(i,this_box) * sin_sum(i,this_box))

     qwxy =  un * ( -2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
            *hx(i,this_box)*hy(i,this_box) )
     qwxz =  un * ( -2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
            *hx(i,this_box)*hz(i,this_box) )
     qwyz =  un * ( -2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
            *hy(i,this_box)*hz(i,this_box) )

     qw(1) = qw(1) + &
             ( un * ( 1.0_DP - 2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
             *hx(i,this_box)*hx(i,this_box)))
     qw(2) = qw(2) + qwxy
     qw(3) = qw(3) + qwxz
     qw(5) = qw(5) + &
             ( un * ( 1.0_DP - 2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
             *hy(i,this_box)*hy(i,this_box)))
     qw(6) = qw(6) + qwyz
     qw(9) = qw(9) + &
             ( un * ( 1.0_DP - 2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
             *hz(i,this_box)*hz(i,this_box)))

  END DO
  !$OMP END PARALLEL DO

  W_tensor_recip(1,1,this_box) = qw(1)
  W_tensor_recip(2,1,this_box) = qw(2)
  W_tensor_recip(3,1,this_box) = qw(3)
  W_tensor_recip(2,2,this_box) = qw(5)
  W_tensor_recip(3,2,this_box) = qw(6)
  W_tensor_recip(3,3,this_box) = qw(9)

  DO is = 1, nspecies

     DO im = 1, nmols(is,this_box)

        this_locate = locate(im,is,this_box)
        IF( .NOT. molecule_list(this_locate,is)%live) CYCLE

        xcmi = molecule_list(this_locate,is)%xcom
        ycmi = molecule_list(this_locate,is)%ycom
        zcmi = molecule_list(this_locate,is)%zcom

        DO ia = 1, natoms(is)

           piix = atom_list(ia,this_locate,is)%rxp - xcmi
           piiy = atom_list(ia,this_locate,is)%ryp - ycmi
           piiz = atom_list(ia,this_locate,is)%rzp - zcmi
           charge = nonbond_list(ia,is)%charge

           recip_11 = 0.0_DP
           recip_21 = 0.0_DP
           recip_31 = 0.0_DP
           recip_22 = 0.0_DP
           recip_23 = 0.0_DP
           recip_33 = 0.0_DP

           !$OMP PARALLEL DO DEFAULT(SHARED) &
           !$OMP SCHEDULE(STATIC) &
           !$OMP PRIVATE(i,arg,factor) &
           !$OMP REDUCTION(+:recip_11, recip_21, recip_31) &
           !$OMP REDUCTION(+:recip_22, recip_23, recip_33)
           DO i = 1, nvecs(this_box)

              arg = hx(i,this_box)*atom_list(ia,this_locate,is)%rxp + &
                    hy(i,this_box)*atom_list(ia,this_locate,is)%ryp + &
                    hz(i,this_box)*atom_list(ia,this_locate,is)%rzp

              factor = Cn(i,this_box)*2.0_DP*(-cos_sum(i,this_box)*DSIN(arg) + &
                       sin_sum(i,this_box)*DCOS(arg))*charge

              recip_11 = recip_11 + factor*hx(i,this_box)*piix
              recip_21 = recip_21 + factor* 0.5_DP*(hx(i,this_box)*piiy+hy(i,this_box)*piix)
              recip_31 = recip_31 + factor* 0.5_DP*(hx(i,this_box)*piiz+hz(i,this_box)*piix)
              recip_22 = recip_22 + factor*hy(i,this_box)*piiy
              recip_23 = recip_23 + factor* 0.5_DP*(hy(i,this_box)*piiz+hz(i,this_box)*piiy)
              recip_33 = recip_33 + factor*hz(i,this_box)*piiz

           END DO
           !$OMP END PARALLEL DO

           W_tensor_recip(1,1,this_box) = W_tensor_recip(1,1,this_box) + recip_11
           W_tensor_recip(2,1,this_box) = W_tensor_recip(2,1,this_box) + recip_21
           W_tensor_recip(3,1,this_box) = W_tensor_recip(3,1,this_box) + recip_31
           W_tensor_recip(2,2,this_box) = W_tensor_recip(2,2,this_box) + recip_22
           W_tensor_recip(2,3,this_box) = W_tensor_recip(2,3,this_box) + recip_23
           W_tensor_recip(3,3,this_box) = W_tensor_recip(3,3,this_box) + recip_33


       end

    end

 end

W_tensor_recip(1,2,this_box) = W_tensor_recip(2,1,this_box)
W_tensor_recip(1,3,this_box) = W_tensor_recip(3,1,this_box)
W_tensor_recip(3,2,this_box) = W_tensor_recip(2,3,this_box)

end Compute_System_Ewald_Reciprocal_Force

!-----------------------------------------------------------------------------
