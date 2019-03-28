MODULE parameter_shared
!
! Purpose:
!   To declare data to share between routines.
IMPLICIT NONE
SAVE 
!character(LEN=50) :: filename, pos_filename
character(LEN=50) :: sampled_pos_filename
!INTEGER :: nat ! number of atoms
INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
REAL, ALLOCATABLE, DIMENSION(:) :: sampled_time, sampled_energy
!CHARACTER(LEN=2) :: str_of_center_atoms
!INTEGER :: num_of_kind_center_atoms 
END MODULE parameter_shared
