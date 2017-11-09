! Module GCMC. Contains the routines needed for Grand Cananonical Monte Carlo. 
! Implemmented from Frenkel p. 126- 134
    module gcmc_module
#   include 'control_simulation.h'
    use commons
    implicit none
    real(kind=8) :: beta,zz,inv_zz,big_lambda,vol

        contains

            subroutine init_gcmc()
                    implicit none
                    print *,"  * Initializing Grand Canonical Monte Carlo (GCMC) parameters   "
                    !mu = -0.9 ! Chemical potential ! now in mfa_input
                    print *,"   Chemical Potential: mu= ",mu
                    beta = 1./temp 
                    inv_vol = inv_boundary(1)*inv_boundary(2)*inv_boundary(3)
                    vol = 1./inv_vol
                    big_lambda = 1./sqrt(2*pi*temp) ; print *, '   Thermal wavelength: big_lambda= ', big_lambda !thermal wave length, without mass yet, planck constant h = 1
                    zz = exp(beta*mu)/big_lambda**3 ; print *,'   zz=',zz  
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
                            print *,"   * Warning: zero liquid particles !!"
                            return
                        end if
! Choose a random liquid particle                         
                        o_part = int(dble(n_liq)*uni()) + part_init_d + 1

                        call ener(1,r0(:,o_part),o_part,eno)
!deb                        print *,'eno mc_exc ',eno
                        print *,'energy',eno
                        arg = dble(n_liq)*exp(beta*eno)*inv_vol*inv_zz
!deb                        print*,'n_liq,eno,arg',n_liq,eno,arg
                        if (uni() < arg) then
!                            print *,'o_part',o_part
                            r0(:,o_part)=r0(:,n_liq) ! put n_liq in the index of the removed particle. 
!        Effectively removing o_part
                            n_liq = n_liq - 1 ! update  current number of particles. 
                            n_chain_d = n_chain_d - 1 
                            n_part = n_part - 1
                            n_mon_tot = n_mon_tot - 1
                            print *, "Removing a particle: N=",n_liq
                        end if

                    else ! Add particles 

                        rn(1) = uni()*boundary(1) 
                        rn(2) = uni()*boundary(2) 
                        rn(3) = uni()*boundary(3) 
                        call ener(2,rn(:),o_part,enn) ! o_part not used here 
                        print *,'energy',enn 
                        arg = zz*vol*exp(-beta*enn)/(dble(n_liq+1))
!                        print *,'arg=',zz,beta,enn,n_liq,arg ; stop ! debug
                        if (uni() < arg) then
! WARN: analyze if n_liq and n_chain_d are not exactly the same all over the program
                            n_liq = n_liq + 1
                            n_chain_d = n_chain_d + 1 
                            n_part = n_part +1 
                            n_mon_tot = n_mon_tot +1
                            r0(:,n_part) = rn(:)
                            a_type(n_part) = 3 ! we are doing GCMC only with liquid particles

! Give the inserted particle a random velocity from a Maxwell distribution 

                             v(1,n_part)  =  sqrt(temp/mass_type(a_type(n_part)))*rnor()
                             v(2,n_part)  =  sqrt(temp/mass_type(a_type(n_part)))*rnor()
                             v(3,n_part)  =  sqrt(temp/mass_type(a_type(n_part)))*rnor()

                             print '(a,3f15.5)','velocity',v(:,n_part)



                            print '(a,i6,3f17.6)', "Adding a particle: N_liq= ",n_liq,rn(:)
                        end if


                    end if ! decide if adding or removing particles

                    print *,"  * Liquid particles = ",n_liq

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
                    real(kind=8)  :: inv_r_2,inv_sqrt_r_2,v_wall
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
#             if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#             endif

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

!  Add energy from wall 
              call ener_wall(r(3),eno)

            print *, '  *  Energy of removed particle ',eno



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
                  if(r_2 < 0.84 ) then ! if new particle is too close to an existing one. Return with huge energy
                      eno = 1.E10  !log(huge(r_2))
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

!  Add energy from wall 
              call ener_wall(r(3),eno)


            print *, ' Energy of added particle = ',eno
        end select

!          print *,'eno',eno !; stop
!          return eno

          ! Debugging
          if(debug) then
          end if

  end subroutine ener

    subroutine ener_wall(zp,eno) 
            real(kind=8) ,intent(in) :: zp
            real(kind=8) ,intent(inout) :: eno
            real(kind=8)  :: v_wall,inv_z

!       Compute energy from wall                  

! WARN: assuming i_type - 3 for liquid particles 

!           Energy  from bottom wall  (copied from v_fluid_wall)

              r_dummy = sigma_wall(3)/zp
              v_wall =  abs(a_wall(3))*r_dummy**9 - a_w*r_dummy**3    !int with bottom wall

                 eno = eno + v_wall 

!           Energy from top wall 
              inv_z = 1./(z_space_wall-zp)
              r_dummy = sigma_wall(3)*inv_z

              v_wall = v_wall + abs(a_wall(3))*r_dummy**9 - a_w*r_dummy**3        


                 eno = eno + v_wall 


    end subroutine ener_wall

    subroutine write_out_system_input()
! writes out a version of system_input with the current number of particles
! This should be used as system_input for a new run to account for the change in number of particles in GCMC

        open(unit=174,file='system_input.new',status='unknown')
    
        write(174,*) "# This file sets the system characteristics for mfa_prog"
        write(174,*) "# It is in free format but the number of files must be held fixed"
        write(174,*) "4"
        write(174,'(a)') "#Grafted polymers: n_mon n_chain wall:[number of cells in x] [number of cells in y]" 
        write(174,'(4i8)') n_mon,n_chain, n_cell_w_x,n_cell_w_y
        write(174,'(a)') "# The drop or melt:  [n_mon_d] [n_chain_d] [drop_cell_x] [drop_cell_y] [z_skin]"
        write(174,'(2i8,3f16.5)') n_mon_d,n_chain_d,drop_cell_x,drop_cell_y,z_skin
        write(174,*) "# non additive sigmas [sig23] [sig24] [sig34]"
        write(174,'(3f16.5)') delta_sig(1:3)
        write(174,*) "# Parameters A and s_w for wall-fluid interaction with: V=A (s_w/z) -(s_w/z)^6"
        write(174,'(4f16.5)')   a_w,sigma_w,a_w4,sigma_w4  
        write(174,*) "# External constant force for Poseuille flow.  [Fx]   [Fy]   [Fz]"
        write(174,'(5f16.5)')  const_force(:)
    
        close(unit=174)


    end subroutine write_out_system_input

    end module gcmc_module
