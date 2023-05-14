module utils_velocity
  implicit none
  public :: determine_velocity
   private
 contains  
subroutine determine_velocity(eng_tot,pos_wrt_bh)
 
     real, intent(in):: eng_tot,pos_wrt_bh
     real :: vel_star,pe,ke,tot_e
     tot_e = eng_tot
     do while(tot_e/=0.)
       !ideally use the potential energy from GR calculations
       pe = -1/pos_wrt_bh
       ke = 0.5*vel_star**2
       tot_e = pe+ke
       if (tot_e /= 0.) then
          vel_star = vel_star + 0.01 
       endif
       print*,vel_star,"vel star",tot_e,"tot e"
     enddo
  end subroutine  determine_velocity

end module utils_velocity
