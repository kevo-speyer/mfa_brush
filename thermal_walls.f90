subroutine thermal_walls()
#include 'control_simulation.h'
use commons
use ziggurat

implicit none




    do i_part = 1 , n_part
 r0(3,i_part )


! Si está en el skin

! If top wall
    fac= sqrt(kb*temp_top/m)
       v_new(:) = fac*(/rnor(),rnor(),sqrt(2.)*sqrt(-log( uni() ) )/)

       v(:,i_part) = v_new(:)
    end do












end subroutine thermal_walls
