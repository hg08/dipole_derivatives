MODULE wannier_center_module
!
! Purpose:
!   To define the derived data type for a wannier center
IMPLICIT NONE
TYPE :: wannier_center
  CHARACTER(LEN=2) :: wannier_center_name
  INTEGER :: molecular_id
  REAL :: charge
  REAL, DIMENSION(3) :: coord 
  REAL, DIMENSION(3) :: image_coord 
END TYPE wannier_center

! The array atom_info can be shared by subroutines  
TYPE(wannier_center), ALLOCATABLE, DIMENSION(:,:) :: wannier_center_info
END MODULE wannier_center_module
