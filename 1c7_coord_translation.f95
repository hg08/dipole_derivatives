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
  ! 2019.03.29          Gang Huang                  coord_translation (version 1.7). 

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

  INTEGER ::i, iatom, jatom
  INTEGER,INTENT(IN) :: nat                         ! number of atoms
  INTEGER, INTENT(IN) :: n_samples                  ! n_samples = INT(nmo/ns)
  REAL(kind=4), INTENT(IN) :: a, b, c               ! Parameters of the simulation box 
  CHARACTER(2), INTENT(IN) :: host_atom, guest_atom ! to define the types of host and guest atoms
  INTEGER :: u, ierror
  INTEGER :: ll,mm,nn
  REAL(4) :: min_dist_image_value
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
  ll = 0
  mm = 0
  nn = 0
  min_dist_image_value = 0
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
    ! Loop for each host atom
    inner: DO iatom = 1, nat
      oxygen: if (TRIM(wannier_center_info(iatom, i)%wannier_center_name) == TRIM(host_atom)) then
       ! Loop for each guest atom
        inner2: DO jatom = 1, nat
          hydrogen: if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
            !=============================================================================================
            !NOTICE: we should use wannier_center_info(:, :)%image_coord(:) to calculate the distance, NOT
            !wannier_center_info(:,:)%coord(:)
             tmp = (wannier_center_info(jatom, i)%image_coord(1) - wannier_center_info(iatom, i)%image_coord(1))**2  &
                  +(wannier_center_info(jatom, i)%image_coord(2) - wannier_center_info(iatom, i)%image_coord(2))**2  &
                  +(wannier_center_info(jatom, i)%image_coord(3) - wannier_center_info(iatom, i)%image_coord(3))**2 
             if (tmp .LE. R_cutoff_square) then                
               wannier_center_info(jatom, i)%molecular_id = wannier_center_info(iatom, i)%molecular_id
               !wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
               !wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(1) 
               !wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(1) 
             !elseif ((tmp .GE. 103.75) .OR. (tmp .LE. 205.76)) then
             else
               CALL min_dist_image(wannier_center_info(iatom, i),wannier_center_info(jatom, i),a,b,c,min_dist_image_value,ll,mm,nn) 
               if (min_dist_image_value .LE. R_cutoff_square) then
                 !WRITE(*,*) ll
                 wannier_center_info(jatom, i)%molecular_id = wannier_center_info(iatom, i)%molecular_id
                 wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - ll*a
                 wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - mm*b
                 wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - nn*c
               end if
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
SUBROUTINE min_dist_image(p1,p2,a,b,c,min_value,ll,mm,nn)
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
INTEGER, INTENT(INOUT) :: ll,mm,nn
INTEGER, DIMENSION(3) :: min_location
TYPE(wannier_center), INTENT(IN) :: p1,p2         ! p denotes one atom or wannier center
REAL(4), DIMENSION(3,3,3) :: arr
REAL(4), INTENT(INOUT) :: min_value
!Initialization
l = 0
m = 0
n = 0
arr(:,:,:) = 0
!======================================================================================================
!Operations: to calculate all the image's coordinates.
! There are totally 27 images, therefore, there are 27 distances to be calculated, ie., the array arr
! has 27 elements.
! Be carefule, if I change the image_coord(:) into coord(:), the calculated dipole moment will be wrong!
! In that case, because of PBC, I will assign wrong coordinate for some D and X atoms. This mistake can 
! be avoided if I use image_coord(:) to calculate distance between atoms (or wannier centers).
!======================================================================================================
do l = -1,1
  do m = -1,1
    do n = -1,1
      arr(l+2,m+2,n+2) =           &
        (p1%image_coord(1) + l*a - p2%image_coord(1))**2  &
        +(p1%image_coord(2) + m*b - p2%image_coord(2))**2 &
        +(p1%image_coord(3) + n*c - p2%image_coord(3))**2 
    enddo
  enddo
enddo
!==================
! Return the result
!==================
! The smallest distance between host atom and the guest's image
min_value = MINVAL(arr)
! The index of the element of the array, which deterimine the very image that is closest to the host atom 
min_location = MINLOC(arr)
! the interger ll,mm, nn can help us to find out the location of the guest atom's image
ll = min_location(1)-2
mm = min_location(2)-2
nn = min_location(3)-2

END SUBROUTINE min_dist_image
