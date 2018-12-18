

subroutine magnetic_force(mode)
#include 'control_simulation.h'
#ifdef MAGNETIC
      use commons
      implicit none 
      integer, intent(in) :: mode
!      real (kind=8), save :: m_sat,V_part,r_loop_min,r_loop_max,i_loop,cte,mag_force,mag_m_z
!      real (kind=8), parameter :: mu_0 = 1.0
!      integer :: n_loop
     
        select case(mode)
        case(0)   ! Initialization 

         !   m_sat =1.0  ! Saturation magnetization
            V_part = 4.*pi*(sigma(4,4)/2.)**3/3. ! Volume of the particle
         !   r_loop_min =1.0 ! Minimum radius of the loop
         !   r_loop_max=2.0  ! maximun radius of the loop
         !   i_loop=1.0 ! Current of the loop
         !   n_loop=1   ! Number of turns of the loop
            cte=dble(n_loop)*mu_0*i_loop/(2*pi*2*(r_loop_max-r_loop_min))
            mag_m_z=m_sat*V_part !magnetic moment of the particle

            print *," Magnetic momemtum of particles =",mag_m_z
            print *,"Volume=",V_part
            print *,"cte=",cte
            print *,"m_sat=",m_sat
            print *,"r_loop_min=",r_loop_min
            print *,"r_loop_max=",r_loop_max
            print*,"i_loop=",i_loop
            print*,"n_loop=",n_loop
            print*,"sigma=",sigma(4,4)
            print*,"a_eff=",a_eff
            print*,"z_sk=",z_sk

        case(1)  ! Compute magnetic force on particles. Calculate the magnetic field of a circular loop on the z axis using Biot Savart law
           ! and aproximate tha field of the spiral as the radial integral of the magnetic field of a single loop (asuming that the
           ! distance bwteen loops is small compare to the radius of the spiral ). Then calculate the force by taking - gradient(field)  
               do i_part = part_init_e + 1, n_part
!               i_part= part_init_e + 1 !debug
 !               r0(3,i_part) = 5  ! debug
                mag_force=cte*mag_m_z*r0(3,i_part)* &
                    (-r_loop_min/((r_loop_min**2+r0(3,i_part)**2)**(1.5))  &
                    -1/(r_loop_min*sqrt(r_loop_min**2+r0(3,i_part)**2)+r_loop_min**2+r0(3,i_part)**2) & 
                     +r_loop_max/((r_loop_max**2+r0(3,i_part)**2)**(1.5))     & 
                     +1/(sqrt(r_loop_max**2+r0(3,i_part)**2)*(sqrt(r_loop_max**2+r0(3,i_part)**2)+r_loop_max)))
               
! v_mag = v_mag + fljsjfls                
  !                print *,'mag_force=',mag_force ! debug
!               stop !debug

                force(3,i_part) = force(3,i_part) + mag_force

              end do
       case(2) ! simple magnetic force model.Calculate the magnetic force using a simple model of a potential that decays as 1/(zdistance)**2
!  mag_energy = 0.0           
            do i_part=part_init_e+1,n_part
                mag_force=a_eff/(r0(3,i_part)+z_sk)**2
                force(3,i_part) = force(3,i_part) + mag_force
! v_mag = v_mag + fljsjfls                
              end do 

        end select
                
#endif 
end subroutine magnetic_force
