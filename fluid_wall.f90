
      subroutine fluid_wall(inter_type)
        
! Computes fluid-wall interactions:
!
!  inter_type= 1 : Lennard Jones interaction between wall and fluid particles 
!  inter_type= 2 : Aw [ (s_w/z)^9-(s_w/z)^3 ] 
!  inter_type= 3 : Aw [ (s_w/z)^9-(s_w/z)^3 ] in bottom wall and top wall=hard wall
!  inter_type= 4 : Hard walls in top and bottom walls
#include 'control_simulation.h'
      use commons ; implicit none
      real (kind=8) :: inv_z
      integer, intent(in) :: inter_type
      logical, parameter :: debug_fw=.false.
!      real (kind=8) :: t1,t2 ! debug

!prof      call cpu_time(t1)
!
      v_fluid_wall = 0.

      select case(inter_type)

      case(1)  ! Lennard Jones interaction between fluid and wall atoms 
#       ifdef STARS
          ftw(:)=0.
          fbw(:)=0.
#       endif
#       ifdef FORCE_SWITCH_ON
          if(r_time.lt.r_2_min_time) then
              r_2_min = (1.-r_time/r_2_min_time)*r_2_min_init!*range_2(:,:)!ORI*r_2_min_init
          else
              r_2_min = 0.
          end if
#       endif


      case(2)     !  (1/z)^9-(1/z)^3 interaction 

#   if WALL==2    
#       ifndef ASYM_WALLS /* Asymmetric walls */         
!DEC$ IVDEP
          do i_part = 1,n_mon_tot

              ! *** Bottom wall interaction

#       ifndef NO_WARNS
              if (r0(3,i_part) < 0. ) then
                  print*,"[fluid_wall] Warn:",i_part," particle inside the wall ! z=",r0(3,i_part),"t=",i_time
                  print '(a/)',"[fluid_wall] Forcing a move in z "
                  r0(3,i_part) = 0.1 
              end if
#       endif

              i_type = a_type(i_part)

              ! HPC 
              inv_z = 1./r0(3,i_part) 
              r_dummy = sigma_wall(i_type)*inv_z
              v_fluid_wall = v_fluid_wall + abs(a_wall(i_type))*r_dummy**9 - a_w*r_dummy**3    !int with bottom wall
              ! HPC
              force(3,i_part) = force(3,i_part) +   9.*abs(a_wall(i_type))*(sigma_wall(i_type))**9*(inv_z)**10 
              force(3,i_part) = force(3,i_part)   - 3.*a_wall(i_type)*(sigma_wall(i_type))**3*(inv_z)**4                              

              ! ***  Top wall interaction. Warn: the force sign must be in opposite direction    

#       ifndef NO_WARNS
              if (r0(3,i_part) > z_space_wall ) then
                  print*,"[fluid_wall] Warn: ",i_part,"particle inside the wall ! z=",r0(3,i_part),"t=",i_time, "type= ",a_type(i_part)
                  print '(a/)',"[fluid_wall] Forcing a movement in z "
                  r0(3,i_part) = z_space_wall - 0.1 
              end if
#      endif
              ! HPC      
              inv_z = 1./(z_space_wall-r0(3,i_part))
              r_dummy = sigma_wall(i_type)*inv_z

              v_fluid_wall = v_fluid_wall + abs(a_wall(i_type))*r_dummy**9 - a_w*r_dummy**3        

              !        note: different sign in top wall 

              force(3,i_part) = force(3,i_part)    - 9.*abs(a_wall(i_type))*(sigma_wall(i_type))**9*(inv_z)**10  
              force(3,i_part) = force(3,i_part)    + 3*a_wall(i_type)*(sigma_wall(i_type))**3*(inv_z)**4

          end do
#       endif /* not def ASYM_WALLS*/
#       if ASYM_WALLS == 1

          do i_part = 1,n_mon_tot

              ! *** Bottom wall interaction

#       ifndef NO_WARNS
              if (r0(3,i_part) < 0. ) then
                  print*,"[fluid_wall] Warn:",i_part," particle inside the wall ! z=",r0(3,i_part),"t=",i_time
                  print '(a/)',"[fluid_wall] Forcing a move in z "
                  r0(3,i_part) = 0.1 
              end if
#       endif

              i_type = a_type(i_part)

              ! HPC 
              inv_z = 1./r0(3,i_part) 
              r_dummy = sigma_wall(2,i_type)*inv_z
              v_fluid_wall = v_fluid_wall + abs(a_wall(2,i_type))*r_dummy**9 - a_wall(2,i_type)*r_dummy**3    !int with bottom wall
              ! HPC
              force(3,i_part) = force(3,i_part) +   9.*abs(a_wall(2,i_type))*(sigma_wall(2,i_type))**9*(inv_z)**10 
              force(3,i_part) = force(3,i_part)   - 3.*a_wall(2,i_type)*(sigma_wall(2,i_type))**3*(inv_z)**4                              

              ! ***  Top wall interaction (possibly DIFFERENT from bottom wall)

#       ifndef NO_WARNS
              if (r0(3,i_part) > z_space_wall ) then
                  print*,"[fluid_wall] Warn: ",i_part,"particle inside the wall ! z=",r0(3,i_part), \ 
                  "t=",i_time, "type= ",a_type(i_part)
                  print '(a/)',"[fluid_wall] Forcing a movement in z "
                  r0(3,i_part) = z_space_wall - 0.1 
              end if
#      endif
              ! HPC      
              inv_z = 1./(z_space_wall-r0(3,i_part))
              r_dummy = sigma_wall(1,i_type)*inv_z

              v_fluid_wall = v_fluid_wall + abs(a_wall(1,i_type))*r_dummy**9 - a_wall(1,i_type)*r_dummy**3        

              !        note: different sign in top wall 

              force(3,i_part) = force(3,i_part)    - 9.*abs(a_wall(1,i_type))*(sigma_wall(1,i_type))**9*(inv_z)**10  
              force(3,i_part) = force(3,i_part)    + 3*a_wall(1,i_type)*(sigma_wall(1,i_type))**3*(inv_z)**4

          end do
#endif          /* Asymmetric walls == 1 */          
#       if ASYM_WALLS == 2

          do i_part = 1,n_mon_tot

              ! *** Bottom wall interaction

#       ifndef NO_WARNS
              if (r0(3,i_part) < 0. ) then
                  print*,"[fluid_wall] Warn:",i_part," particle inside the wall ! z=",r0(3,i_part),"t=",i_time
                  print '(a/)',"[fluid_wall] Forcing a move in z "
                  r0(3,i_part) = 0.1 
              end if
#       endif

              i_type = a_type(i_part)

              ! HPC 
              inv_z = 1./r0(3,i_part) 
              r_dummy = sigma_wall(2,i_type)*inv_z
              v_fluid_wall = v_fluid_wall + abs(a_wall(2,i_type))*r_dummy**9 - a_wall(2,i_type)*r_dummy**3    !int with bottom wall
              ! HPC
              force(3,i_part) = force(3,i_part) +   9.*abs(a_wall(2,i_type))*(sigma_wall(2,i_type))**9*(inv_z)**10 
              force(3,i_part) = force(3,i_part)   - 3.*a_wall(2,i_type)*(sigma_wall(2,i_type))**3*(inv_z)**4                              

              ! ***  NO Top wall interaction: adiabatic wall in thermal_walls.f90

          end do
#endif          /* Asymmetric walls == 2 */          

#endif /* WALL = 2*/













      case(3) !-----------  ! Only bottom wall with 9-3 potential and top wall pure repulsive (hard wall)

#if WALL==3
          ! Typical for DROPS

          !------- Bottom Wall interaction 

          do i_part = 1,n_mon_tot

#ifndef NO_WARNS
              if (r0(3,i_part) <0. ) then
                  print*,"Warn:",i_part," particle inside the wall ! z=",r0(3,i_part),"t=",i_time
                  r0(3,i_part) = 0.1 
              end if
#endif

              inv_z = 1./r0(3,i_part) 

              r_dummy = sigma_w*inv_z


              v_fluid_wall = v_fluid_wall + abs(a_w)*r_dummy**9 - a_w*r_dummy**3    !int with bottom wall

             

              ! HPC
              force(3,i_part) = force(3,i_part) +   9.*abs(a_w)*(sigma_w)**9*(inv_z)**10 
              force(3,i_part) = force(3,i_part)   - 3.*a_w*(sigma_w)**3*(inv_z)**4                              

              ! 
              !------- Top  Wall interaction 

              if (r0(3,i_part) > (z_space_wall-z_head) .and. v(3,i_part) > 0. ) then
                  v(3,i_part) = -v(3,i_part)
              end if

          end do
#endif
      case(4) !-----------  Hard Walls 

          ! The hard walls are located in z=0 and z=D


          do i_part = 1,n_mon_tot

              ! Bottom wall
              if ( r0(3,i_part)   < 0.01 .and. v(3,i_part) < 0.  ) then
                  v(3,i_part) = -v(3,i_part)
              end if
              ! Top wall 
              if ( r0(3,i_part) > (z_space_wall-0.01) .and. v(3,i_part) > 0. ) then
                  v(3,i_part) = -v(3,i_part)
              end if
          end do

      case   default
          print *," sub FLUID_WALL: inter_type value must be 1: LJ or 2: int -1/z^3 + 1/z^9 3: drop 4: hard top and lower walls "
          print*, "Change it! Stopping here "
          stop
      end select

      !profile
      !        call cpu_time(t2)
      !        print*,"fw time=", t2-t1

      if(debug_fw) then
          do i_part=1,n_mon_tot
              if(r0(3,i_part)>(z_space_wall-1.).or.r0(3,i_part)<1.)               & 
              print*,"[fluid_wall]z_too_close",i_time,a_type(i_part),r0(3,i_part)
          end do
          r_dummy =     sum( sum(force(:,:)**2,dim=1) ,dim =1 )
          print *,"[fluid_wall]quad_m_force", i_time,sqrt(r_dummy)/real(n_mon_tot)
          print*,"v_fluid_wall=",i_time,v_fluid_wall   
      end if 



         end subroutine fluid_wall
