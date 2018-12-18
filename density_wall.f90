subroutine density_wall(mode)
#include 'control_simulation.h'
use commons
use ziggurat

! Computes wall quantities for the verification of contact theorem. (With Nacho)
implicit none
integer, intent(in) :: mode
real(kind=8) :: v_new(3)
real(kind=8) :: fac


    select case(mode)
    case(0) ! initialization
        print *, "  * Initializing density wall: contact theorem" 

! variables for density of particles getting in and out

        rho_top_in  = 0.0
        rho_top_out = 0.0
        rho_bot_in  = 0.0
        rho_bot_out = 0.0

        temp_top_in(:)  = 0.0
        temp_top_out(:) = 0.0
        temp_bot_in (:) = 0.0
        temp_bot_out(:) = 0.0

    case(2)

! Criterium: density bin starts just to the side of thermal skin region        
    do i_part = 1 , n_mon_tot
 
!    Top wall 

        If ( (r0(3,i_part)>(z_space_wall-dens_wall_skin-thermal_skin)) .and. (r0(3,i_part)<=z_space_wall-thermal_skin )) then
            if ( (v(3,i_part)>0.0)) then !! change velocity if the particle is above interwall spacing -thermal skin  commons
            rho_top_in = rho_top_in + 1.0
                   temp_top_in(1) = temp_top_in(1) + mass(i_part)*v(2,i_part)**2
                   temp_top_in(2) = temp_top_in(2) + mass(i_part)*v(3,i_part)**2
            else ! particles is going out of the wall region 
                   rho_top_out = rho_top_out + 1.0 
                   temp_top_out(1) = temp_top_out(1) + mass(i_part)*v(2,i_part)**2
                   temp_top_out(2) = temp_top_out(2) + mass(i_part)*v(3,i_part)**2
            end if
       end if

!       Bottom wall       

       If ( (r0(3,i_part) < thermal_skin+ dens_wall_skin) .and.( r0(3,i_part) >= thermal_skin)) then
           if (v(3,i_part) < 0.0) then ! particles getting in the wall zone
           rho_bot_in = rho_bot_in + 1.0
                temp_bot_in(1) = temp_bot_in(1) + mass(i_part)*v(2,i_part)**2 ! Ty
                temp_bot_in(2) = temp_bot_in(2) + mass(i_part)*v(3,i_part)**2 !Tz
            else ! particles going out of the walls
                rho_bot_out = rho_bot_out + 1.0
                temp_bot_out(1) = temp_bot_out(1) + mass(i_part)*v(2,i_part)**2 ! Ty
                temp_bot_out(2) = temp_bot_out(2) + mass(i_part)*v(3,i_part)**2 !Tz

            end if
       end if
  end do

  end select 
end subroutine density_wall
