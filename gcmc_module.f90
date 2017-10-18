! Module GCMC. Contains the routines needed for Grand Cananonical Monte Carlo. 
! Implemmented from Frenkel p. 126- 134
    module gcmc_module
#   include 'control_simulation.h'
    use commons
    implicit none
    real(kind=8) :: beta,inv_vol,zz,inv_zz,mu,big_lambda,vol

        contains

            subroutine init_gcmc()
                    implicit none
                    print *,"  * Initializing Grand Canonical Monte Carlo (GCMC) parameters   "
                    mu = 0.0001 ! Chemical potential
                    print *," Chemical Potential: mu= ",mu
                    beta = 1./temp 
!                    pi = 4.0*atan(1.0)
                    inv_vol = inv_boundary(1)*inv_boundary(2)*inv_boundary(3)
                    vol = 1./inv_vol
                     big_lambda = 1./sqrt(2*pi*temp) ; print *, 'Thermal wavelength:i big_lambda= ', big_lambda !thermal wave length, without mass yet, planck constant h = 1
                    zz = exp(beta*mu)/big_lambda ; print *,'zz=',zz  
                    inv_zz = 1./zz

            end subroutine init_gcmc

! Grand Canonical MC particle exchange 

            subroutine mc_exc()
                    use commons
                    use ziggurat
                    implicit none
                    integer :: o_part
                    real(kind=8) :: arg,eno,enn,rn(3)

                    if (uni() < 0.5) then ! Decide if removing or adding particles
                        if (n_liq == 0 ) then
                            print *,"   * Warning: zero particles !!"
                            return
                        end if
! Choose a random liquid particle                         
                        o_part = int(dble(n_liq)*uni()) + part_init_d + 1

                        call ener(1,r0(:,o_part),o_part,eno)
                        print *,'eno mc_exc ',eno
                        arg = dble(n_liq)*exp(beta*eno)*inv_vol*inv_zz
                        print*,'n_liq,eno,arg',n_liq,eno,arg
                        if (uni() < arg) then
                            r0(:,o_part)=r0(:,n_liq) ! put n_liq in the index of the removed particle. 
                            ! Effectively removing o_part
                            n_liq = n_liq - 1 ! update  current number of particles. 
                            print *, "Removing a particle: N=",n_liq
                        end if

                    else ! Add particles 

                        rn(1) = uni()*boundary(1) 
                        rn(2) = uni()*boundary(2) 
                        rn(3) = uni()*boundary(3) 

                        call ener(2,rn(:),o_part,enn) ! o_part not used here 
                        arg = zz*vol*exp(-beta*enn)/(dble(n_liq)+1)
!                        print *,'arg=',zz,beta,enn,n_liq,arg ; stop ! debug
                        if (uni() < arg) then
                            n_liq = n_liq + 1
                            r0(:,n_liq) = rn(:)
                            a_type(n_liq) = 3 ! we are doing GCMC only with liquid particles
                            print *, "Adding a particle: N= ",n_liq
                        end if


                    end if ! decide if adding or removing particles


            end subroutine mc_exc

            subroutine ener(mode,r,oo_part,eno) ! calculates the energy of removing the chosen particle at position r
                    use commons
                    implicit none
!                    integer, intent (in) :: mode
                    real(kind=8), intent(in) :: r(3)
                    integer, intent(in)      :: mode,oo_part
                    real(kind=8), intent(out) :: eno
                    real (kind=8) :: delta_v(3),r_versor(3),g_rand,rrc,f_ipart(3) ! Needed for DPD_EMBEDDED only
                    real (kind=8) :: l_eps,r_cut2,r_61
                    real (kind=8) :: f_cou_real(3)
                    real(kind=8)  :: inv_r_2,inv_sqrt_r_2
                     logical, parameter :: debug=.false.
                     
!
! -------   LJ Energy calculation
!
!mode =1: compute energy for a chosen EXISTING particle
! mode=2: compute energy for NEW particle 

          eno = 0.

          i_type = a_type(part_init_d+1)  ! chosen particle always liquid

        !  print *, ff_list(:,1) ; stop
        select case (mode)

        case(1) ! energy for an existing particle
            

!          do j_part = part_init_d+1, part_init_d+n_liq
              do j_part = part_init_d+1, part_init_d + oo_part-1  ! energy for i_part < oo_part

! Get interaction range acording to pair                  
              !              j_part = ff_list(i_neigh,i_part)
              j_type = a_type(j_part)
!              print *,'j_type ',j_type
!              print *,'j_part ',size(a_type),i_type,j_type
               

              r_cut2 = range_2(i_type,j_type)

              delta_r(1) = r(1) - r0(1,j_part)
              delta_r(2) = r(2) - r0(2,j_part)
              delta_r(3) = r(3) - r0(3,j_part)

              !----- PBC ----

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  Check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                  l_eps = epsil(i_type,j_type)


                  r_61= sigma_2(i_type,j_type)*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  eno = eno + l_eps*pot_loc

              end if

          end do
! Note: splitting loops oo_part is explictedly avoided

          do j_part = part_init_d + oo_part+1, part_init_d+n_liq  ! Energy count of i_part > oo_part

              j_type = a_type(j_part)

              delta_r(1) = r(1) - r0(1,j_part)
              delta_r(2) = r(2) - r0(2,j_part)
              delta_r(3) = r(3) - r0(3,j_part)

              !----- PBC ----

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  Check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                  l_eps = epsil(i_type,j_type)


                  r_61= sigma_2(i_type,j_type)*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  eno = eno + l_eps*pot_loc

              end if
          end do



        case(2) ! Energy for a NEW randomly located particle


          do j_part = part_init_d + 1, part_init_d+n_liq  ! Energy count of i_part > oo_part

              r_cut2 = range_2(i_type,j_type)


              delta_r(1) = r(1) - r0(1,j_part)
              delta_r(2) = r(2) - r0(2,j_part)
              delta_r(3) = r(3) - r0(3,j_part)

              !----- PBC ----

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  Check whether interaction takes place

              
             
              if( r_2 .lt. r_cut2 ) then
                  if(r_2 < 0.49 ) then ! if new particle is too close to an existing one. Return with huge energy
                      eno = log(huge(r_2))
                      return 
                  end if
                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                  l_eps = epsil(i_type,j_type)


                  r_61= sigma_2(i_type,j_type)*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  eno = eno + l_eps*pot_loc

              end if
          end do
        end select

!          print *,'eno',eno !; stop
!          return eno

          ! Debugging
          if(debug) then
          end if

  end subroutine ener

    end module gcmc_module
