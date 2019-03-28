SUBROUTINE label(nat,n_samples)
  
  ! 0a) Purpose: 
  ! The subroutine label to each atom. The total number of the molecular id is T*N, where T is the 
  ! total steps of the reduced trajectory obtained from sampling, and N is the total number of atoms.
  !
  ! 0b) Record of revisions:
  !    Date             Programmer                  description of Change
  !    ====             ==========                  ====================
  ! 2019.03.27          Gang Huang                  Original code 

  !==============
  !1) Declaration
  !==============
  USE wannier_center_module
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  INTEGER ::i, iatom
  INTEGER,INTENT(IN) :: nat ! number of atoms
  INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
  
  !===============
  ! Initialization
  !===============
  i =0
  iatom = 0

  !=======================================
  ! Labeling each atom with a molecular id
  !=======================================
  do i =1,n_samples
    DO iatom = 1, nat
       wannier_center_info(iatom, i)%molecular_id = (i-1)*nat + iatom
    ENDDO
  enddo

  !===========
  !For testing
  !===========
  WRITE(*,*) 'molecular id:', wannier_center_info(1, 1)%molecular_id
  WRITE(*,*) 'molecular id:', wannier_center_info(1, 2)%molecular_id

END SUBROUTINE label
