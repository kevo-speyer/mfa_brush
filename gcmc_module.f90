! Module GCMC. Contains the routines needed for Grand Cananonical Monte Carlo. 
! Implemmented from Frenkel p. 126- 134
    module gcmc_module
#   include 'control_simulation.h'
    implicit none
    real(kind=8) :: beta,inv_vol,zz,inv_zz,mu,big_lambda,vol

        contains

            subroutine init_gcmc()
                    use commons
                    implicit none
                    print *,"  * Initializing Grand Canonical Monte Carlo (GCMC) parameters   "
                    mu = 12000.0 ! Chemical potential
                    print *," Chemical Potential: mu=",mu
                    beta = 1./temp 
!                    pi = 4.0*atan(1.0)
                    inv_vol = inv_boundary(1)*inv_boundary(2)*inv_boundary(3)
                    vol = 1./inv_vol
                    big_lambda = 1./sqrt(2*pi*temp) !thermal wave length, without mass yet, planck constant h = 1
                    zz = exp(beta*mu)/big_lambda
                    inv_zz = 1./zz

            end subroutine init_gcmc

! Grand Canonical MC particle exchange 

            subroutine mc_exc()
                    use commons
                    use ziggurat
                    integer :: o_part
                    real(kind=8) :: arg,eno,enn,rn(3)


                    if (uni() < 0.5) then ! Decide if removing or adding particles
                        if (n_liq == 0 ) then
                            print *,"   * Warning: zero particles !!"
                            return
                        end if
! Choose a random liquid particle                         
                        o_part = int(dble(n_liq)*uni()) + part_init_d + 1
                        call ener(r0(:,o_part),eno)
                        arg = dble(n_liq)*exp(beta*eno)*inv_vol*inv_zz

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
                        call ener(rn(:),enn)
                        arg = zz*vol*exp(-beta*enn)/(dble(n_liq)+1)
                        if (uni() < arg) then
                            r0(:,n_liq+1) = rn(:)
                            n_liq = n_liq + 1
                            print *, "Adding a particle: N=",n_liq
                        end if


                    end if ! decide if adding or removing particles


            end subroutine mc_exc

            subroutine ener(r,eno) ! calculates the energy of removing the chosen particle at position r
                    use commons
                    implicit none
!                    integer, intent (in) :: mode
                    real(kind=8), intent(in) :: r(3)
                    real(kind=8), intent(out) :: eno
                    real (kind=8) :: delta_v(3),r_versor(3),g_rand,rrc,f_ipart(3) ! Needed for DPD_EMBEDDED only
                    real (kind=8) :: l_eps,r_cut2,r_61
                    real (kind=8) :: f_cou_real(3)
                    real(kind=8)  :: inv_r_2,inv_sqrt_r_2
                     logical, parameter :: debug=.false.
                     
!
! -------   LJ Energy calculation
!


          eno = 0.

          i_type = a_type(part_init_d+1)  ! chosen particle always liquid

        !  print *, ff_list(:,1) ; stop

          do j_part = part_init_d+1, part_init_d+n_liq

              !              j_part = ff_list(i_neigh,i_part)
              j_type = a_type(j_part)

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

!          return eno

          ! Debugging
          if(debug) then
          end if

  end subroutine ener

    end module gcmc_module
