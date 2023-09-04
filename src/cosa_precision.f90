!-----------------------------------------------------------------------
!
! Module to hold numerical precision configuration for COSA
!
! Adrian Jackson, EPCC, The University of Edinburgh
! June 2023
!
!-----------------------------------------------------------------------


module cosa_precision

  implicit none
       
     integer, public, parameter :: cosa_int = selected_int_kind(9)
     integer, public, parameter :: cosa_long = selected_int_kind(15)
     integer, public, parameter :: cosa_real = selected_real_kind(15, 307)
     integer, public, parameter :: cosa_single = selected_real_kind(6, 37)
     
   contains
  
end module cosa_precision
