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
!        fac_top = 
!        fac_bot = 

    case(2)
!     print *,kb,thermal_skin,top_thermal_wall,bottom_thermal_wall !deb

   !     stop ! deb
    do i_part = 1 , n_part
 

        If (r0(3,i_part)>z_space_wall-thermal_skin) then !! change velocity if the particle is above interwall spacing -thermal skin  commons
       fac= sqrt(kb*top_thermal_wall*inv_mass(i_part))
       v_new(:) = fac*(/rnor(),rnor(),-sqrt(2.)*sqrt(-log( uni() ) )/)
       v(:,i_part) = v_new(:)
       endif


       If (r0(3,i_part)<thermal_skin) then
       fac= sqrt(kb*bottom_thermal_wall*inv_mass(i_part))
       v_new(:) = fac*(/rnor(),rnor(),sqrt(2.)*sqrt(-log( uni() ) )/)
       v(:,i_part) = v_new(:)
!       print*,'bottom wall',v_new(:)
       end if
  end do

  end select 
end subroutine thermal_walls
