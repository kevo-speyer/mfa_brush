

subroutine magnetic_force(mode)
#include 'control_simulation.h'
      use commons
      implicit none 
      integer, intent(in) :: mode
!      real (kind=8), save :: m_sat,V_part,r_loop_min,r_loop_max,i_loop,cte,mag_force,mag_m_z
!      real (kind=8), parameter :: mu_0 = 1.0
!      integer :: n_loop
     
        select case(mode)
        case(0)   ! Initialization 

            m_sat =1.0  ! Saturation magnetization
            V_part = 4.*pi*(sigma(4,4)/2.)**3/3. ! Volume of the particle
            r_loop_min =1.0 ! Minimum radius of the loop
            r_loop_max=2.0  ! maximun radius of the loop
            i_loop=1.0 ! Current of the loop
            n_loop=1   ! Number of turns of the loop
            cte=n_loop*mu_0/(2*pi*2*(r_loop_max-r_loop_min))
            mag_m_z=m_sat*V_part !magnetic moment of the particle

            print *,"Magnetic momemtum of particles =",mag_m_z

        case(1)  ! Compute magnetic force on particles 
!            do i_part = part_init_e + 1, n_part
                i_part= part_init_e + 1
                r0(3,i_part) = 3.0  ! debug
                mag_force=cte*mag_m_z*r0(3,i_part)* &
                    (-r_loop_min/((r_loop_min**2+r0(3,i_part)**2)**(3/2))  &
                     -1/(r_loop_min*sqrt(r_loop_min**2+r0(3,i_part)**2)+r_loop_min**2+r0(3,i_part)**2) & 
                     +r_loop_max/((r_loop_max**2+r0(3,i_part)**2)**(3/2))     & 
                     +1/(sqrt(r_loop_max**2+r0(3,i_part)**2)*(sqrt(r_loop_max**2+r0(3,i_part)**2)+r_loop_max)))
               
                  print *,'mag_force=',mag_force ! debug

                  stop !debug

                force(3,i_part) = force(3,i_part) + mag_force

!            end do
        end select
                

end subroutine magnetic_force
