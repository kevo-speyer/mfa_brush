subroutine density_wall(mode)
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
 
!    Top wall 

        If (r0(3,i_part)>(z_space_wall-dens_wall_skin) ) then
            if ( (v(3,i_part)>0.0)) then !! change velocity if the particle is above interwall spacing -thermal skin  commons
            rho_top_in = rho_top_in + 1.0
            else ! particles is going out of the wall region 
                   rho_top_out = rho_top_out + 1.0 
            end if
       end if

!       Bottom wall       

       If ( (r0(3,i_part) < dens_wall_skin) ) then
           if (v(3,i_part) < 0.0) then ! particles getting in the wall zone
           rho_bot_in = rho_bot_in + 1.0
            else ! particles going out of the walls
           rho_bot_out = rho_bot_out + 1.0

            end if
       end if
  end do

  end select 
end subroutine density_wall
