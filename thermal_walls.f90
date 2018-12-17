subroutine thermal_walls(mode)
#include 'control_simulation.h'
use commons
use ziggurat

implicit none
integer, intent(in) :: mode
real(kind=8) :: v_new(3)
real(kind=8) :: fac


    select case(mode)
    case(0) ! initialization
        print *, "  * Initializing thermal walls" 

! variables for density of particles getting in and out

        rho_top_in  = 0.0
        rho_top_out = 0.0
        rho_bot_in  = 0.0
        rho_bot_out = 0.0

    case(2)

    do i_part = 1 , n_mon_tot
 
! Top wall 

!        If (r0(3,i_part)>(z_space_wall-thermal_skin) .and. (v(3,i_part)>0.0)) then !! change velocity if the particle is above interwall spacing -thermal skin  commons
        If (r0(3,i_part)>(z_space_wall-thermal_skin) ) then
            if ( (v(3,i_part)>0.0)) then !! change velocity if the particle is above interwall spacing -thermal skin  commons

            fac= sqrt(kb*top_thermal_wall*inv_mass(i_part))
            v_new(:) = fac*(/rnor(),rnor(),-sqrt(2.)*sqrt(-log( uni() ) )/)

! Before changing particle velocity compute force=dP/dt= m*(v_new-v_old)/dt. Contribution to the pressure
       
!            z-direction
            mean_force_top_wall = mean_force_top_wall+ mass(i_part)*inv_dt*(v_new(3) -v(3,i_part))
            v(:,i_part) = v_new(:)
            end if
       end if

!       Bottom wall       

!ori       If ( (r0(3,i_part)<thermal_skin) .and. (v(3,i_part)<0.0)) then
       If ( (r0(3,i_part) < thermal_skin) ) then
           if (v(3,i_part) < 0.0) then ! particles getting in the wall zone
            fac= sqrt(kb*bottom_thermal_wall*inv_mass(i_part))
            v_new(:) = fac*(/rnor(),rnor(),sqrt(2.)*sqrt(-log( uni() ) )/)

! Before changing particle velocity compute force=dP/dt= m*(v_new-v_old)/dt. Contribution to the pressure
       
! z-direction
            mean_force_bot_wall = mean_force_bot_wall+ mass(i_part)*inv_dt*(v_new(3) -v(3,i_part))

            v(:,i_part) = v_new(:)

!       print '(a,3g17.4)','bottom wall: vx,vy,vz',(v_new(1)**2+v_new(2)**2+v_new(3)**2)/3.
!deb       print '(a,3g17.4)','vz wall: vx,vy,vz',v_new(3)
            end if
       end if
  end do

  end select 
end subroutine thermal_walls
