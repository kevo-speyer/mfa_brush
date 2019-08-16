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

    dp_wall_top  = 0.0
    dp_wall_bot  = 0.0

    do i_part = 1 , n_mon_tot
 
! Top wall 
        if (r0(3,i_part)>(z_space_wall-thermal_skin) .and. (v(3,i_part)>0.0)) then !! change velocity if the particle is above interwall spacing -thermal skin  commons
       fac= sqrt(kb*top_thermal_wall*inv_mass(i_part))
       v_new(:) = fac*(/rnor(),rnor(),-sqrt(2.)*sqrt(-log( uni() ) )/)
!       v_new(:) = fac*(/rnor(),rnor(),-1.0*sqrt(-log( uni() ) )/)

       dp_wall_top =  dp_wall_top + mass(i_part)*(v(3,i_part)-v_new(3))

       v(:,i_part) = v_new(:)
!deb       print '(a,3g17.4)','vz wall: vx,vy,vz',v_new(3)
       endif

! Bottom wall 

       if ( (r0(3,i_part)<thermal_skin) .and. (v(3,i_part)<0.0)) then
       fac= sqrt(kb*bottom_thermal_wall*inv_mass(i_part))
       v_new(:) = fac*(/rnor(),rnor(),sqrt(2.)*sqrt(-log( uni() ) )/)
!       v_new(:) = fac*(/rnor(),rnor(),sqrt(-log( uni() ) )/)
       v(:,i_part) = v_new(:)
!       print '(a,3g17.4)','bottom wall: vx,vy,vz',(v_new(1)**2+v_new(2)**2+v_new(3)**2)/3.
!deb       print '(a,3g17.4)','vz wall: vx,vy,vz',v_new(3)
       dp_wall_bot =  dp_wall_bot - mass(i_part)*(v(3,i_part)-v_new(3))
       end if
  end do

  end select 
end subroutine thermal_walls
