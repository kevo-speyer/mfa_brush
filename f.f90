subroutine verlet_velocities()
! * Updates velocities, forces and accels using Velocity verlet
! * Updates wall positions
! in the DPD scheme
use commons
 

 




 



 
 



 

 



 

 













 




                        




 


                     





                     








                  
                  
                  
 

    


  implicit none
  real (kind=8) ::  T_inst, T_dummy
  real (kind=8) ::  T_inst_mean = 0.
  integer :: i_head
  logical, parameter :: debug = .false.





!   ---  Zero all for the heads
                f_on_heads(:,:) =  force(:,1:part_init_d:n_mon)




! Update accelerations with the new force values

        do i_part = 1,n_mon_tot
            a(1,i_part) = force(1,i_part)*inv_mass(i_part)
            a(2,i_part) = force(2,i_part)*inv_mass(i_part)
            a(3,i_part) = force(3,i_part)*inv_mass(i_part)
!        end do

! Update velocities

!        do i_part = 1 , n_mon_tot
!            v(1,i_part) = v_half(1,i_part) + 0.5*dt*a(1,i_part)       ! old force(1,i_part)*inv_mass(i_part)
!            v(2,i_part) = v_half(2,i_part) + 0.5*dt*a(2,i_part)       ! old force(2,i_part)*inv_mass(i_part)
!            v(3,i_part) = v_half(3,i_part) + 0.5*dt*a(3,i_part)       ! old force(3,i_part)*inv_mass(i_part)
            v(1,i_part) = v(1,i_part) + 0.5*dt*a(1,i_part)    
            v(2,i_part) = v(2,i_part) + 0.5*dt*a(2,i_part)    
            v(3,i_part) = v(3,i_part) + 0.5*dt*a(3,i_part)    
        end do

! NOTE: if DPD_VV is defined the energies and temp are calculated afterwards in new_dpd_fd.f90










!!KEVIN 10/2013 Droplet compatibility added (2 ==3 and 0 == 1)
!!Droplet not compatible with moving wall



!      Fix the head beads at the original position
        do i_chain = 1 ,n_chain
                i_head = 1+n_mon*(i_chain-1)
                v(:,i_head) = 0. ! 0 the velocity of the heads
                force(:,i_head) = 0.! zero force
                a(:,i_head) = 0. ! zero accels
        end do


!
!   ---- Measurement of some quantities
!

!---  Kinetic energies and temperatures
!
!
! Instantaneous temperature

      T_inst = 0.
      do i_part = 1,n_mon_tot
          T_dummy =   mass(i_part)*( v(1,i_part)**2+ v(2,i_part)**2 + v(3,i_part)**2 )
          T_inst = T_inst + T_dummy
      end do
      t_fluid = T_inst
  



! *** Write out Instantanous T and energies
!

      t_fluid = 0.5*t_fluid ! *inv_dt_2
!
      v_total = v_wall_wall + v_fluid_wall + v_fluid_fluid + v_intra_molec






!        print *,"v_fluid_fluid=",v_fluid_fluid/dble(n_part) !; stop

!
! Update wall velocity

      if(f_twall(n_dim).eq.2) then
          v_total=v_total-r0_twall(n_dim)*va_spring_twall(n_dim)*0.5
      end if
!
      t_total = t_wall + t_fluid
      e_total = v_total + t_total
      
!
!!  Temperature and energy calculation and writing out
!

   time_ave_count_2 = time_ave_count_2 + 1

   if(time_ave_count_2.eq.10) then
       time_ave_count_2 = 0




        write(61,'(i7,3g17.5)') i_time,v_fluid_fluid,v_intra_molec,v_fluid_wall
        write(60,'(i7,3g17.5)') i_time,e_total*inv_N,v_total*inv_N, t_total*inv_N


    end if



!
!   *** DEBUGGING STUFF ***
!
if(debug) then

print '(a,4e15.5)',"[Venergies]" , v_wall_wall , v_fluid_wall , v_fluid_fluid , v_intra_molec
print*,"before_vel^2= ",i_time,sum(v(:,1:n_mon_tot)**2, dim=1)*inv_N
print*,"V_tot_Ekin",i_time,v_total*inv_N,t_total*inv_N
print*,"heads_vel= ",sum(v(:,1:n_mon*n_chain:n_mon), dim=1)
print*,"accel^2= ",i_time,sum(a(:,1:n_mon_tot)**2, dim=1)/real(n_mon_tot)
        write(77,'(i7,3f17.5)') i_time,v_total*inv_N, e_total*inv_N,t_total*inv_N
end if

end subroutine verlet_velocities
