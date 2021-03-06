      subroutine fluid_fluid()
#include 'control_simulation.h'
      use commons
#ifdef _OPENMP
      use omp_lib
      use Par_Zig_mod, only: par_rnor,par_uni
#else
      use ziggurat, only: rnor,uni
#endif

     implicit none
     real (kind=8) :: delta_v(3),r_versor(3),g_rand,rrc,f_ipart(3) ! Needed for DPD_EMBEDDED only
     real (kind=8) :: l_eps,r_cut2,r_61
     real (kind=8) :: f_cou_real(3)
     real(kind=8)  :: inv_r_2,inv_sqrt_r_2
! cache blocking 
     integer :: ii_part,ii_neigh,i,j
      logical, parameter :: debug=.false.
      
#       if SYSTEM == 2 || SYSTEM == 3
            v_coul = 0.
#       endif

!!#ifdef DPD_EMBEDDED
! NOTE: this is at the beggining of dpd_forces_ll when the DPD forces are not
! embededd here. This must be added if dpd_forces is never called in the program
! flow. 

! Force switch-on update  

! 
!---  Define effective minimium distance used for force "switch on" at the beginning 
!     This is done for each time step

#ifdef FORCE_SWITCH_ON

      if(r_time.lt.r_2_min_time) then
       r_2_min = (1.-r_time/r_2_min_time)*r_2_min_init !*range_2(:,:)!ORI*r_2_min_init
      else
       r_2_min = 0.
      end if
#endif      
      v_fluid_fluid = 0.
!
! LJ fluid-fluid force and V  calculation 
!

!DEBUG

!BEGIN PARALLEL ZONE

!Warning: Paralelization not adapted for SYMMETRY=1 
# if BIN_TYPE == 0 || BIN_TYPE == 1
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(r_2_min,force,v_fluid_fluid,n_part,a_type,ff_list, range_2, r0,inv_boundary, boundary, epsil, sigma_2, e_shift,sig,mass,friction,v, inv_range_2, r_cut_dpd_2)
#ifdef _OPENMP
    ith=omp_get_thread_num()
#endif


!$OMP DO SCHEDULE(STATIC,10) REDUCTION(+:force,v_fluid_fluid)     
    do i_part = 1,n_part  !n_mon_tot= brushes + droplet/melt
          i_dummy = ff_list(0,i_part)

# elif BIN_TYPE == 2  /* cell_list.f90 */
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(r_2_min,force,v_fluid_fluid,n_part,a_type,l_cell,n_cells_tot,n_nei_cells,part_in_cell,n_cells,r_nei,l_nei,cell_neigh_ls, range_2, r0,inv_boundary, boundary, epsil, sigma_2, e_shift,sig,mass,friction,v, inv_range_2, r_cut_dpd_2)
#ifdef _OPENMP
    ith=omp_get_thread_num()
#endif

!$OMP DO SCHEDULE(STATIC,10) REDUCTION(+:force,v_fluid_fluid)
    do i_cell = 1, n_cells_tot ! loop over all cells
        i_part = part_in_cell(i_cell) ! get first particle in i_cell
        
        do while(i_part .ne. 0)  ! loop over particles in i_cell
            
# endif

          i_type = a_type(i_part)
          f_ipart(:) =0.
          q_part=i_part !dummy variable for OpenMP paralelization


# if BIN_TYPE == 0 || BIN_TYPE == 1
             do i_neigh = 1, i_dummy
                 j_part = ff_list(i_neigh,i_part)
 

# elif BIN_TYPE == 2  /* cell_list.f90 */
           do j_dummy = 1, n_nei_cells ! loop neighbor cells
                j_cell = cell_neigh_ls(i_cell,j_dummy) ! j_cell = neighbour cell
        
                j_part = part_in_cell(j_cell)
                do while(j_part .ne. 0) ! loop over neighbour particles in j_cell
                                    
# endif
            
             j_type = a_type(j_part)

                  
              r_cut2 = range_2(i_type,j_type)

              ! HPC

              delta_r(1) = r0(1,i_part) - r0(1,j_part)
              delta_r(2) = r0(2,i_part) - r0(2,j_part)
              delta_r(3) = r0(3,i_part) - r0(3,j_part)

              !----- PBC ----
              ! HPC

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then
#           ifdef FORCE_SWITCH_ON 
                  r_2 = max(r_2,r_2_min) 
#           endif                 

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                  l_eps = epsil(i_type,j_type)


                  r_61= sigma_2(i_type,j_type)*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  v_fluid_fluid = v_fluid_fluid + l_eps*pot_loc

                  r_dummy = l_eps*(-12*r_12+6*r_6)*inv_r_2

                  force_loc(1) = r_dummy*delta_r(1)
                  force_loc(2) = r_dummy*delta_r(2)
                  force_loc(3) = r_dummy*delta_r(3)
                  ! HPC

#       if BIN_TYPE == 0 ||  BIN_TYPE == 2/* binning.f90 */
                  f_ipart(1) = f_ipart(1) -  force_loc(1)
                  f_ipart(2) = f_ipart(2) -  force_loc(2)
                  f_ipart(3) = f_ipart(3) -  force_loc(3)
#       endif

                  force(1,j_part) = force(1,j_part) + force_loc(1)
                  force(2,j_part) = force(2,j_part) + force_loc(2)
                  force(3,j_part) = force(3,j_part) + force_loc(3)

#   if SYMMETRY == 1 /* bulk */                  
! Does not work with OMP
! Potential part of press tensor (see viscosity.f90) ~ virial contribution
         do i = 1,3
             do j = 1,3               
                 press_tensor(i,j) =  press_tensor(i,j) - force_loc(i)*delta_r(j)              
             end do 
         end do 
!WRONG?  end if
    
#   endif

!
!  --- DPD calculation if DPD has not its own cut-off  
!

#if THERMOSTAT == 0
#     ifndef DPD_CUT_OFF
#           if SYMMETRY == 0
#               if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3 /* if there are brushes with fix heads*/
!        Exclude heads from DPD calculation 

!#                  ifndef FREE_HEADS
!#                  endif

#                if WALL == 1 /* explicit wall*/
                  if(i_type.eq.n_type) cycle       ! Excluding thermostat for wall-fluid
                  if(j_type.eq.n_type) cycle       ! Excluding heads thermostat for wall-fluid
#                endif                    


#                endif

#                  if PINNED == 1
                  if(i_type.eq.3) cycle       !  Excluding  fixed particles
                  if(i_type.eq.4) cycle       ! 
                  if(j_type.eq.3) cycle       ! 
                  if(j_type.eq.4) cycle       ! 
#                  endif 
#         endif /* SYMMETRY */

                  call dpd_forces(inv_sqrt_r_2,force)

#       endif /* not DPD_CUT_OFF */
#endif /*THERMOSTAT = 0 */


#       if SYSTEM == 2  || SYSTEM  == 3 /* charged systems */

! Coulomb potential and interaction in REAL space 

                  call ewald_real()

#       endif

              end if ! if the particle is inside the interaction sphere 

!
! DPD has its own cutoff radius 
!

#ifdef DPD_CUT_OFF
#       if THERMOSTAT == 0
#                      if SYMMETRY == 0
#                              if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3
              !                               Exclude heads from DPD calculation 
#                                   ifndef FREE_HEADS
              if(i_type.eq.1) then
#                   if BIN_TYPE == 2
                  j_part = r_nei(j_part)
#                   endif

                  cycle       ! Excluding heads !!!
              end if

              if(j_type.eq.1) then
#                   if BIN_TYPE == 2
                  j_part = r_nei(j_part)
#                   endif
                  cycle       ! Excluding heads !!!
              end if
#                                   endif
#                                   if WALL == 1 /* explicit wall*/
              if(i_type.eq.n_type) cycle       ! Excluding thermostat for wall-fluid
              if(j_type.eq.n_type) cycle       ! Excluding heads thermostat for wall-fluid
#                                   endif                    
#                                  if PINNED == 1

              if(i_type.eq.3) cycle       !  Excluding  fixed particles
              if(i_type.eq.4) cycle       ! 
              if(j_type.eq.3) cycle       ! 
              if(j_type.eq.4) cycle       ! 
#                                  endif 

#                              endif
#                       endif 

              if(r_2 < r_cut_dpd_2 ) then 
                  
                  
                  call dpd_forces(inv_sqrt_r_2,force)

              end if

#       endif /* THEMOSTAT = 0 */
#endif /* DPD_CUT_OFF */

#if BIN_TYPE == 0 || BIN_TYPE == 1
          end do ! loop over particles neighbors (j_part)
#elif BIN_TYPE == 2
                   
                    j_part = r_nei(j_part)

                end do ! loop over particles in j_cell
            end do ! loop over neighbor cells 

           j_part =  part_in_cell(i_cell)
            do while(j_part .ne. 0)  ! loop over particles in the same cell
                if (i_part .lt. j_part) then !Count interactions just once

             j_type = a_type(j_part)

              r_cut2 = range_2(i_type,j_type)

              ! HPC

              delta_r(1) = r0(1,i_part) - r0(1,j_part)
              delta_r(2) = r0(2,i_part) - r0(2,j_part)
              delta_r(3) = r0(3,i_part) - r0(3,j_part)

              !----- PBC ----
              ! HPC

              delta_r(1) = delta_r(1) - boundary(1)*int(2.*delta_r(1)*inv_boundary(1))
              delta_r(2) = delta_r(2) - boundary(2)*int(2.*delta_r(2)*inv_boundary(2))
#       if SYMMETRY == 1
              delta_r(3) = delta_r(3) - boundary(3)*int(2.*delta_r(3)*inv_boundary(3))
#       endif

              r_2 =  delta_r(1)*delta_r(1)  + delta_r(2)*delta_r(2) +  delta_r(3)*delta_r(3)

              !-----  check whether interaction takes place

              if( r_2 .lt. r_cut2 ) then
#           ifdef FORCE_SWITCH_ON 
                  r_2 = max(r_2,r_2_min) 
#           endif                 

                  inv_r_2 = 1./r_2
                  inv_sqrt_r_2 = sqrt(1./r_2)
                  l_eps = epsil(i_type,j_type)


                  r_61= sigma_2(i_type,j_type)*inv_r_2 
                  r_6 = r_61*r_61*r_61

                  r_12 = r_6*r_6
                  pot_loc = (r_12-r_6) - e_shift(i_type,j_type)

                  v_fluid_fluid = v_fluid_fluid + l_eps*pot_loc

                  r_dummy = l_eps*(-12*r_12+6*r_6)*inv_r_2

                  force_loc(1) = r_dummy*delta_r(1)
                  force_loc(2) = r_dummy*delta_r(2)
                  force_loc(3) = r_dummy*delta_r(3)
                  ! HPC

                 f_ipart(1) = f_ipart(1) -  force_loc(1)
                  f_ipart(2) = f_ipart(2) -  force_loc(2)
                  f_ipart(3) = f_ipart(3) -  force_loc(3)

                  force(1,j_part) = force(1,j_part) + force_loc(1)
                  force(2,j_part) = force(2,j_part) + force_loc(2)
                  force(3,j_part) = force(3,j_part) + force_loc(3)

#   if SYMMETRY == 1 /* bulk */                  
! Does not work with OMP
! Potential part of press tensor (see viscosity.f90) ~ virial contribution
         do i = 1,3
             do j = 1,3               
                 press_tensor(i,j) =  press_tensor(i,j) - force_loc(i)*delta_r(j)              
             end do 
         end do 
!WRONG?  end if
    
#   endif
 
!
!  --- DPD calculation if DPD has not its own cut-off  
!

#if THERMOSTAT == 0
#     ifndef DPD_CUT_OFF
#           if SYMMETRY == 0
#               if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3 /* if there are brushes with fix heads*/
!        Exclude heads from DPD calculation 

!#                  ifndef FREE_HEADS
!#                  endif

#                if WALL == 1 /* explicit wall*/
                  if(i_type.eq.n_type) cycle       ! Excluding thermostat for wall-fluid
                  if(j_type.eq.n_type) cycle       ! Excluding heads thermostat for wall-fluid
#                endif                    


#                endif

#                  if PINNED == 1
                  if(i_type.eq.3) cycle       !  Excluding  fixed particles
                  if(i_type.eq.4) cycle       ! 
                  if(j_type.eq.3) cycle       ! 
                  if(j_type.eq.4) cycle       ! 
#                  endif 
#         endif /* SYMMETRY */

                  call dpd_forces(inv_sqrt_r_2,force)

#       endif /* not DPD_CUT_OFF */
#endif /*THERMOSTAT = 0 */


#       if SYSTEM == 2  || SYSTEM  == 3 /* charged systems */

! Coulomb potential and interaction in REAL space 

                  call ewald_real()

#       endif

              end if ! if the particle is inside the interaction sphere 

!
! DPD has its own cutoff radius 
!

#ifdef DPD_CUT_OFF
#       if THERMOSTAT == 0
#                      if SYMMETRY == 0
#                              if SYSTEM == 0 || SYSTEM == 1 || SYSTEM == 3
              !                               Exclude heads from DPD calculation 
#                                   ifndef FREE_HEADS
              if(i_type.eq.1) then        ! Excluding heads !!!
#                   if BIN_TYPE == 2
                                    j_part = r_nei(j_part)
#                   endif
                  cycle
              end if

              if(j_type.eq.1) then    ! Excluding heads !!!
#                   if BIN_TYPE == 2
                                    j_part = r_nei(j_part)
#                   endif
                    cycle
              end if
#                                   endif
#                                   if WALL == 1 /* explicit wall*/
              if(i_type.eq.n_type) cycle       ! Excluding thermostat for wall-fluid
              if(j_type.eq.n_type) cycle       ! Excluding heads thermostat for wall-fluid
#                                   endif                    
#                                  if PINNED == 1

              if(i_type.eq.3) cycle       !  Excluding  fixed particles
              if(i_type.eq.4) cycle       ! 
              if(j_type.eq.3) cycle       ! 
              if(j_type.eq.4) cycle       ! 
#                                  endif 

#                              endif
#                       endif 

              if(r_2 < r_cut_dpd_2 ) then 
                  
                  
                  call dpd_forces(inv_sqrt_r_2,force)

              end if

#       endif /* THEMOSTAT = 0 */
#endif /* DPD_CUT_OFF */

           
                
                end if !  i_part < j_part  
            j_part = r_nei(j_part) !gets next particle in this cell     
            end do  ! loop over neighbor particles in i_cell, j_part
 
#endif


#       if BIN_TYPE == 0 || BIN_TYPE == 2/* binning.f90 */
 
          force(1,i_part) = force(1,i_part)  + f_ipart(1)
          force(2,i_part) = force(2,i_part)  + f_ipart(2)
          force(3,i_part) = force(3,i_part)  + f_ipart(3)

#       endif

#   if THERMOSTAT == 1 /* Langevin */
! COMMENT to check NVE ensemble          
! call lgv_forces(force)
#   endif /* Langevin */

#if BIN_TYPE == 0 || BIN_TYPE == 1
end do ! loop over particles

!$OMP END DO 
!$OMP END PARALLEL 
!END PARALLEL ZONE

#elif BIN_TYPE == 2

        i_part = r_nei(i_part)
 
             end do    ! loop over particles in i_cell

      end do ! loop over all cells

!$OMP END DO 
!$OMP END PARALLEL 
!END PARALLEL ZONE


#endif


#ifdef FLUKT 
    if((n_mon_d.or.n_chain_d.or.n_mon_e.or.n_chain_e).gt.2) then
    write(*,*) "the file fort.555 gonna be too big! undefine flag FLUKT or decrease number of particles4 or -3"
    stop
    endif
    if (mod(i_time, 100).eq.0) then
    write(555,'(1i,18f15.4)') i_time,  r0(part_init_d+1:n_mon_tot, :), & 
                              v(:,part_init_d+1:n_mon_tot), force(part_init_d+1:n_mon_tot, :)
    end if
#endif
#                 if BIN_TYPE == 1                  
                  v_fluid_fluid = 0.5*v_fluid_fluid
#                 endif 
#if SYSTEM == 2 || SYSTEM == 3
        v_coul = v_coul + v_coul_self
#endif

! Debugging
if(debug) then
!print '(3f15.5)', force(:,1:n_chain*n_mon:n_mon)
          print *,"v_fluid_fluid=",i_time,v_fluid_fluid/dble(n_mon_tot)
r_dummy =     sum( sum(force(:,:)**2,dim=1) ,dim =1 )
print *,"[fluid_fluid]quad_m_force", i_time,sqrt(r_dummy)/dble(n_mon_tot)
print *, "[fluid_fluid]V=",v_fluid_fluid
!do i_part=1,n_mon_tot
!write(88,'(3i,3f15.4)') i_time, i_part,a_type(i_part),force(:,i_part)
!end do
end if
          !
end subroutine fluid_fluid
