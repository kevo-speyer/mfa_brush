! Routine that allows to add rigidity to the grafted chains
subroutine bending(mode_bend)
#include "control_simulation.h"
use commons
implicit none
integer, intent(in) :: mode_bend
integer ::i,j,l
real(kind=8) ::r_prev(3),r_next(3),dot_pr=0,r_prev_2=0,r_next_2=0, cos_alpha,dir_next(3)=0,dir_prev(3)=0,dir_prev_2=0, dir_next_2=0, F_bend(6),F_mod,delta_alpha
#ifdef BENDING
select case (mode_bend)

case(0)  ! Init  variables 
print*,""
print*," * Simulation with bending stiffnes"
print*," * Bending elastic constant ",k_bend
print*," * Equilubrium angle:",alpha_eq
print*,""

case(1)
v_bend=0  !Reset Bending Potential Energy
Do l=0,n_chain-1 !n_chain , loop over chains
    r_next=r0(:,l*n_mon+2)-r0(:,l*n_mon+1)
    !This lines below correct for perioduc boundry conditions
    r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
    r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
   ! #if SYMMETRY == 1
   ! r_next(3) = r_next(3) - boundary(3) * int(2.*r_next(3)*inv_boundary(3))
   ! #endif
    !End correction fo periodic boundry conditions
    Do i=2,n_mon-1 !n_mon, loop over particles in chain
        !Reset dummy variables
        dot_pr=0
        r_prev_2=0
        r_next_2=0
        dir_prev_2=0
        dir_next_2=0
        !Calculation of bending forces for the neighbors of i
        r_prev=r_next
        r_next=r0(:,l*n_mon+i+1)-r0(:,l*n_mon+i)
        !This lines below correct for perioduc boundry conditions
        r_next(1) = r_next(1) - boundary(1) * int(2.*r_next(1)*inv_boundary(1))
        r_next(2) = r_next(2) - boundary(2) * int(2.*r_next(2)*inv_boundary(2))
        !#if SYMMETRY == 1
        !r_next(3) = r_next(3) - boundary(3) * int(2.*r_next(3)*inv_boundary(3))
        !#endif
        !End correction fo periodic boundry conditions
        Do j=1,3
            dot_pr=dot_pr+r_next(j)*r_prev(j)
            r_prev_2=r_prev_2+r_prev(j)*r_prev(j)
            r_next_2=r_next_2+r_next(j)*r_next(j)
        End do
        cos_alpha=dot_pr/sqrt(r_prev_2*r_next_2)
        dir_next=r_prev*r_next_2-r_next*dot_pr
        dir_prev=-r_next*r_prev_2+r_prev*dot_pr
        Do j=1,3
            dir_prev_2=dir_prev_2+dir_prev(j)*dir_prev(j)
            dir_next_2=dir_next_2+dir_next(j)*dir_next(j)
        End do
        if(dir_next_2.lt.0.000001) then !if alpha is small, use approximation
            ! old ! 
            !dir_next=dir_next*0
            ! DEBUG !
            ! print*, "SMALL_ANGLE_APPROXIMATION_DONE"
            ! new
            dir_next = r_prev / sqrt(r_prev_2) - r_next / sqrt(r_next_2)
            dir_next_2=norm2(dir_next) 
            if(dir_next_2.ne.0.0) then
                dir_next = dir_next / dir_next_2 / sqrt(r_next_2)
            else
                dir_next=dir_next*0    
            end if
        else
        !Below I divide by |r_next|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_next|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_next=dir_next/sqrt(dir_next_2*r_next_2)
        end if
        if(dir_prev_2.lt.0.000001) then !if alpha is small make force 0
            ! old ! 
            !dir_prev=dir_prev*0
            !new
            dir_prev = r_prev / sqrt(r_prev_2) - r_next / sqrt(r_next_2)
            dir_prev_2=norm2(dir_prev) 
            if(dir_prev_2.ne.0.0) then
                dir_prev = dir_prev / dir_prev_2 / sqrt(r_prev_2)
            else
                dir_prev=dir_prev*0    
            end if
       else
        !Below I divide by |r_prev|, because the bending force is proportional to 
        !k_bend*delta_alpha/|r_prev|. This comes from V_vend=1/2*k_bend*delta_alpha_2
            dir_prev=dir_prev/sqrt(dir_prev_2*r_prev_2)
        end if
        if(cos_alpha.ge.1.0) then
            print*, "Error cos_alpha =",cos_alpha, " >= 0"
            delta_alpha=0
        else
            delta_alpha=acos(cos_alpha)-alpha_eq
        end if
        F_mod=k_bend*delta_alpha
        v_bend=v_bend+.5*F_mod*delta_alpha
        Do j=1,3
            F_bend(j)=F_mod*dir_prev(j)
            F_bend(3+j)=F_mod*dir_next(j)
        End do
        Do j=1,3
            force(j,l*n_mon+i-1) = force(j,l*n_mon+i-1) + F_bend(j)
            force(j,l*n_mon+i)   = force(j,l*n_mon+i) - F_bend(j) - F_bend(3+j)
            force(j,l*n_mon+i+1) = force(j,l*n_mon+i+1) + F_bend(3+j)
        End do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  TEST TO CHECK IF THE ROUTINE IS CALCULATING FORCES AND ENERGY
        !  CORRECTLY. 
        !print*, r0, "positions"
        !print*, k_bend, "bending constant"
        !print*, F_bend, "bending FORCE"
        !print*, v_bend, "bending ENERGY"
        !print*,delta_alpha,F_mod/sqrt(r_next_2),dot_product(F_bend(4:6),dir_next*sqrt(r_next_2))
        ! ERASE AFTER CHECKING
        !!!!!!!!!!!!!!!!!!!!!!!!!!!
    End Do
End do
end select

#endif /*close BENDING*/
end subroutine bending
