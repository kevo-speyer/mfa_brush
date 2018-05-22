subroutine thermal_walls()
#include 'control_simulation.h'
use commons
use ziggurat

implicit none
real(kind=8) :: v_new(3)
real(kind=8) :: fac


    do i_part = 1 , n_part
 


! Top wall
        If (r0(3,i_part)>z_space_wall-thermal_skin) then !! change velocity if the particle is above interwall spacing -thermal skin    CAMBIAR 40
            !! INTERWALL SPACING

            fac= sqrt(kb*top_thermal_wall*inv_mass(i_part))
            v_new(:) = fac*(/rnor(),rnor(),-sqrt(2.)*sqrt(-log( uni() ) )/)
            v(:,i_part) = v_new(:)
!            print '(a,3f10.5)',"top wall",v_new(:)
        endif
!   Bottom wall 
! Claudio: va menor
        If (r0(3,i_part)<thermal_skin) then
            fac= sqrt(kb*bottom_thermal_wall/inv_mass(i_part))
            v_new(:) = fac*(/rnor(),rnor(),sqrt(2.)*sqrt(-log( uni() ) )/)
            v(:,i_part) = v_new(:)
!            print '(a,3f10.5)',"bottom wall",v_new(:)
        end if
    end do












end subroutine thermal_walls
