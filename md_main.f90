      program md_pb 
#include 'control_simulation.h'
      use commons 
      use util  
#   ifdef GCMC
      use gcmc_module
      use ziggurat
#   endif
      
      implicit none
      integer  :: tot_time,i
      real (kind=8) :: t0,t1
      integer :: i_total,npav,mc_rand ! aux 
      logical,parameter :: debug =.false.

! [12/2015] Kevin Speyer: Parallelization with OpenMP
! [11/2013]: Droplet compatibility fixed
! [7/2013] : Bending Force for the brush added       
! March 2013: Program passes to Kevin Speyer
! [4/2010] Heavy changes in data structures and improvement of some routines 
! [5/2009]: Additions by Leonid Spirin: pinning, stars, LJ wall and shear
! protocols 
! [10/07/08]: Extended DPD cutoff, for better Temperature control.WARN: It changes
! a lot the friction 
! [20/11/07] Coulomb interactions added. It performs Ewald sum with dipole
! corrections to take into account the walls in Z. SYSTEM = 2 added.
! [18/11/07] Now with SYMMETRY = 1 the program has PBC in three directions for
! BULK simulations 
! [10/10/07] Particle 4 interaction added to the program. It can be also a  chain.
! [09/05/07] Fix an old error when writing confs. f_film discarded
! [2006] Incorporation of C-preprocessor directives
! [18/6/04] Pure velocity verlet  for integration and DPD thermostat implementation of
!           Random and Dissipative forces
! [10/6/04] Shear: constant velocity for the heads of brushes. Force on heads re written.    
! [6/6/04]  Added a 9-3 flat wall potential (instead of explicit wall atoms) 
! [5/5/04]  Dynamic allocation of forces and position. New input file: system_input
! [14/4/04] claudio: using now ifort as the main compiler
! for development
! [then] Torsten Kreer : Brushes and melt with shear
! [originally from] Martin Muesser

         call get_walltime(t0)
         print '(/a/)'," *** MFA_PROG *** MD, DPD, and more. [Compilation Date: "//__DATE__//"]"


! ---- Write out current program compilation settings ----------     
!
           call messages()


! --------------- ENDS writing  out of compilation settings ---------------------------     
      
! ----  Inizialization  routines




           call init_system() ! physical system. Reads system_input
           call init_params()   ! md simulation parameters, box dimensions. Read mfa_input
           call init_config()   ! initial conf or read an old one
           call init_obser()    ! Measurable variables init.
           
! GCMC + MD parameters
#       ifdef GCMC
!now in init_params            n_cycl = 500000  ! cycles of MD+GCMC
!now in init_params            nexc= 100   ! MC steps per cycle
          
            npav =  n_part+n_ghost ! First value of npav. It is the expected number of particles. 
            call init_gcmc()

#       endif
#ifdef BENDING     
          call bending(0) ! writes brush bending constants  to log
#endif

#ifdef BENDING_MELT     
          call bending_melt(0) ! writes mel  bending constants  to log
#endif


#ifdef ORIENTATION     
          call orientation(0) ! writes brush bending orientation to log
#endif

          tot_time = n_relax + n_obser 
!          tot_time = 100 ! MD steps per cycle
            
#         ifdef GCMC
          do i_total = 1,n_cycl                       ! total loop: combined MD + GCMC

               print*," * MD+GCMC cycle number: i_total=",i_total

!deb               mc_rand = int(uni()*(npav+nexc)) + 1

!deb               print *,"mc_rand=",mc_rand,n_part ! debug



           !    if(mc_rand <= n_part) then ! Decide if MD or GCMC
               if(uni() <= 0.5) then ! Decide if MD or GCMC

#          endif

           do i_time = 1 , tot_time !  MAIN TIME LOOP 

           r_time = dble(s_time+i_time-1)*dt


! ----  Propagate coordinates 


           call verlet_positions()



           call check_skin  !calculates if it is necessary to update the verlet list. If it is,
#          if BIN_TYPE == 0
 
           !if (f_skin.eq.1) call binning
!           open(unit=232)
!           do i_part = 1,n_part
!               write(232,'(3i6,3f17.5)') i_total,i_time,i_part,r0(:,i_part)
!           end do
               close(232)
            call binning() ! to debug GCMC
#           elif BIN_TYPE == 1 

           if (f_skin.eq.1) call my_binning
#           endif


! ----- Forces to zero

           force(:,:) = 0.0  

! Note: After r_2_min_time, the force is completely switched on 
! This r_2_min is used for calculating the v_fluid_fluid and v_fluid_wall
!
! Stuff for viscosity (it is spread in all routines that calculate forces and in
! viscosity.f90 )     ![VISC] 
!
!! WARN: We are assuming here that THIS routine is the first one in calculating
!forces. That's why
! the pressure tensor is put to zero in each time step. IF this is not the first
! routine computing forces. 
! this should not be here.
#   if SYMMETRY == 1
           press_tensor(:,:) = 0.
#   endif

 
! NOTE: fluid-fluid calculates LJ always, DPD and LGV forces and calls ewald in
! real space
           call fluid_fluid() ! Calculates various forces, including thermostat and LJ


#   if SYMMETRY == 0
#       if WALL != 1 
           call fluid_wall(wall_flag) ! 1= wall atoms, 2= 9-3 potenti , 3 and 4 also valid
#       endif

                                call wall_wall(wall_flag)  ! 1= wall atoms, 2= 9-3 potential
#       if WALL == 1
                                call intra_wall
#       endif
#   endif
           call intra_molec
#ifdef BENDING        
           call bending(1)  ! adds brush bending forces and bending energy
#endif


#ifdef BENDING_MELT        
           call bending_melt(1) ! adds melt bending forces and bending energy
#endif



#ifdef ORIENTATION        
           call orientation(1) ! adds brush orientation bending forces and bending energy
#endif


#   if SYSTEM == 2 || SYSTEM == 3
           call ewald_k(1)  ! coulomb force calculation in K-space for Ewald sum
#         if SYMMETRY == 0
           call dipolar_correction()
#         endif
#   endif

#           ifdef POISEUILLE
           call constant_force() ! Poiseuille flow generation
#           endif


           ! -----  Update  velocities

           call verlet_velocities()

#ifdef DPD_VV                     

!Note: this recalculates Fd with the new velocities and updates F for the begining og the next cycle 
!       with this new value.  

           call new_dpd_fd()  
#endif
#       ifndef GCMC
!----  Observe system after equilibration

           if(i_time.gt.n_relax) call observation 

!----  Make safety copies  to recover from crashes and write out of configurations

           if(mod(i_time,n_safe).eq.0) then
               call store_config(2)  ! writes out conf_xmol and conf_new
#       if STORE == 0
               call store_config(3)  ! Writes out film_xmol and vel.dat
#       elif STORE == 1
               call store_config(4)  ! Writes out film_xmol and vel.dat UNFOLDED
#       endif
           end if

#       endif /* if not defined GCMC */
       end do   ! --------------  ENDS TIME LOOP ----------------

#   ifdef GCMC 

          else ! mcrand: Do MC

              ! Loop of grand canonical MC
              ! Note: we develop GCMC algorithm in serial version first

              !              print *,"hola MC"
              do i = 1,nexc
                  print*,"GCMC step: i_time=",n_cycl*(nexc -1)+ i !debug
                  call mc_exc()
              end do

          end if ! mc_rand 

!----  Observe the system 

!       if(mod(i_total,n_safe)==0) then  ! warn: check if this is correct  
            call observation 
!       end if

!----  Make safety copies  to recover from crashes and write out of configurations

! Note: after each cycle, write

       !ori            if(mod(i_time,n_safe) == 0) then
!            if(mod(i_total,n_safe) == 0) then
                    call store_config(2)  ! writes out conf_xmol and conf_new
#       if STORE == 0
                    call store_config(3) ! Writes out ! film_xmol ! and ! vel.dat
#       elif STORE == 1
                    call store_config(4) ! Writes ! out ! film_xmol ! and ! vel.dat ! UNFOLDED
#       endif
!            end if

               end do ! n_cycl: ! Ends ! loop ! of ! cycles ! MD ! + ! GCMC ! moves 

! Write out system_input.new with updated number of particles                
              
            call write_out_system_input()

#       endif /* ifdef GCMC */


               call obser_out()

       close(20) ! closing mfa_output
!
       call get_walltime(t1)

       print *,'    WALL TIME (s)= ',t1-t0

  end program md_pb
