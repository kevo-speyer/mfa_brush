subroutine constant_force(case_pois)
#include 'control_simulation.h' 
! This routine adds an external constant force to the  monomers.
! Gravity or force in x. If the masses are different, the force should be
! thought as m*g force, where g is read from the force parameter in
! system_input
! The force magnitud is constant and read form system_input
! * the force is added to the particles with a_type = 3 
    use commons
!use ziggurat, only: rnor,uni
!use util ! debug
    implicit none
    logical, parameter :: debug=.false.
    integer, intent(in) :: case_pois    

select case (case_pois)
case(1)
    !Apply constant force to all particles 
    do i_part = 1, n_mon_tot
          force(:,i_part) = force(:,i_part) + ext_force(:,i_part)
    end do

case(2)
!Concentrate melt in the middle of the simulation box in x
    
    do i_part = 1, n_mon_tot
          force(1,i_part) = force(1,i_part) - ext_force(1,i_part)
    end do


    do i_part = part_init_d+1, n_mon_tot
       force(1,i_part) = force(1,i_part) - ext_force(1,i_part) * ( 2 * int( r0(1,i_part) * 2 / boundary(1) ) - 1.)
    end do

end select


end subroutine constant_force
