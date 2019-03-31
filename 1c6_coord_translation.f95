SUBROUTINE coord_translation(guest_atom, host_atom, nat, n_samples, a, b, c)
  
  ! 0a) Purpose: 
  ! The subroutine  coord_transform implement the coordinate translation for the atoms and wannier centers, considering the periodic
  ! boundary condition (PBC) and CLUSTERING.  
  ! 0b) The subroutine 1c6_coord_translation is different from the subroutine 1c5_coord_translation and before in the basic idea to
  ! consider the PBC
  ! 0c) 1c6_coord_translation is much accurater than 1c4_coord_translation. I tested for 20 steps on 1c6_coord_translation. No mistake.
  !     But 1c4_coord_translation gives 4 wrong clustering.
  !
  ! 0e) Record of revisions:
  !    Date             Programmer                  description of Change
  !    ====             ==========                  ====================
  ! 2019.03.29          Gang Huang                  coord_translation (version 1.4.2). 

  !==============
  !1) Declaration
  !==============
  USE wannier_center_module
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  REAL(kind=4),PARAMETER :: R_OD_c = 1.32   ! MAX(1.32,0.6,0.7)  in Angstrom
  REAL(kind=4) :: R_cutoff                  ! The cutoff of radium which is used to determine whether a guest atom is belong to a
                                            ! molecule defined by the host atom. (in Angstrom)
  REAL(kind=4) :: R_cutoff_square
  REAL(kind=4) :: tmp 
  REAL(kind=4) :: min_dist_image 

  INTEGER ::i, iatom, jatom
  INTEGER,INTENT(IN) :: nat                         ! number of atoms
  INTEGER, INTENT(IN) :: n_samples                  ! n_samples = INT(nmo/ns)
  REAL(kind=4), INTENT(IN) :: a, b, c               ! Parameters of the simulation box 
  CHARACTER(2), INTENT(IN) :: host_atom, guest_atom ! to define the types of host and guest atoms
  INTEGER :: u, ierror
  !===============
  ! Initialization
  !===============
  i =0
  iatom = 0 
  jatom = 0
  u = 8        ! UNIT for WRITING
  ierror = 2   
  R_cutoff = 0.0
  R_cutoff_square = 0.0
  tmp = 0.0

  !===========================================================================
  ! Translate the coordinates of gust atoms
  ! Purpose: to consider the PBC when calculate the distance between two atoms
  !===========================================================================
  
  ! To determine the Radium cutoff for the host-guest pair
  if (TRIM(host_atom)=="O" .AND. TRIM(guest_atom)=="H") then
    R_cutoff = 1.32                 ! in Angstrom 
    R_cutoff_square = R_cutoff**2
  else if (TRIM(host_atom)=="H" .AND. TRIM(guest_atom)=="X") then
    R_cutoff = 0.6                  ! in Angstrom 
    R_cutoff_square = R_cutoff**2
  else if (TRIM(host_atom)=="O" .AND. TRIM(guest_atom)=="X") then
    R_cutoff = 0.7                  ! in Angstrom 
    R_cutoff_square = R_cutoff**2
  endif
  
  ! coord.translation AND clustering
  outer: do i =1,n_samples
  !outer: do i =1, 2
    inner: DO iatom = 1, nat
      oxygen: if (TRIM(wannier_center_info(iatom, i)%wannier_center_name) == TRIM(host_atom)) then
        inner2: DO jatom = 1, nat
          hydrogen: if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
            !=====================================================================================
             tmp = (wannier_center_info(jatom, i)%coord(1) - wannier_center_info(iatom, i)%coord(1))**2  &
                  +(wannier_center_info(jatom, i)%coord(2) - wannier_center_info(iatom, i)%coord(2))**2  &
                  +(wannier_center_info(jatom, i)%coord(3) - wannier_center_info(iatom, i)%coord(3))**2 
             if (tmp .LE. R_cutoff_square) then                
               wannier_center_info(jatom, i)%molecular_id = wannier_center_info(iatom, i)%molecular_id
             !elseif ((tmp .GE. 103.75) .OR. (tmp .LE. 205.76)) then
             elseif (  &
               min_dist_image(wannier_center_info(iatom, i),wannier_center_info(jatom, i),a,b,c) .LE. R_cutoff_square  &
                    ) then 
               wannier_center_info(jatom, i)%molecular_id = wannier_center_info(iatom, i)%molecular_id
             endif
          endif hydrogen
        ENDDO inner2 
      endif oxygen
    ENDDO inner
  enddo outer

  !===========
  !For testing
  !===========
  OPEN (UNIT=u,FILE='test_output_'//TRIM(guest_atom)//'_molecular_id.dat',  &
        STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)
    WRITE(u, *) ' Guest_atom ', ' molecular ID '
  do i = 1, n_samples 
  !do i = 1, 2 
    guest2: do jatom = 1, nat 
      if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
        WRITE(u,*) wannier_center_info(jatom, i)%wannier_center_name, &
                   wannier_center_info(jatom, i)%molecular_id
      end if
    end do guest2
  end do
  CLOSE(UNIT=u)

END SUBROUTINE coord_translation 

! function
REAL(4) FUNCTION min_dist_image(p1,p2,a,b,c)
!
! Purpose:
!   To return the minimum of the distance between the image of a host-gust pair
!
USE wannier_center_module
IMPLICIT NONE

!================
! Data Dictionary
REAL(kind=4), INTENT(IN) :: a, b, c               ! Parameters of the simulation box 
INTEGER :: l,m,n
TYPE(wannier_center), INTENT(IN) :: p1,p2         ! p denotes one atom or wannier center
REAL(4), DIMENSION(3,3,3) :: arr
!Initialization
l = 0
m = 0
n = 0
arr(:,:,:) = 0
!==========
!Operations
do l = -1,1
  do m = -1,1
    do n = -1,1
      arr(l+2,m+2,n+2) =           &
        (p1%coord(1) + l*a - p2%coord(1))**2  &
        +(p1%coord(2) + m*b - p2%coord(2))**2 &
        +(p1%coord(3) + n*c - p2%coord(3))**2 
    enddo
  enddo
enddo
! Return the result
min_dist_image = MINVAL(arr)

END FUNCTION min_dist_image

