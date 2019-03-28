SUBROUTINE cluster(guest_atom, host_atom, nat, n_samples)
  
  ! 0a) Purpose: 
  ! The subroutine (cluster) divides the atoms and wannier centers into different groups. Each group is a molecule.
  ! 0b) Record of revisions:
  !    Date             Programmer                  description of Change
  !    ====             ==========                  ====================
  ! 2019.03.28          Gang Huang                  Original code 

  !==============
  !1) Declaration
  !==============
  USE wannier_center_module
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  REAL(kind=4) :: R_cutoff             ! The cutoff of radium which is used to determine whether a guest atom is belong to a
                                       ! molecule defined by the host atom. (in Angstrom)
  REAL(kind=4) :: R_cutoff_square
  INTEGER ::i, iatom, jatom
  INTEGER,INTENT(IN) :: nat                         ! number of atoms
  INTEGER, INTENT(IN) :: n_samples                  ! n_samples = INT(nmo/ns)
  CHARACTER(2), INTENT(IN) :: host_atom, guest_atom ! to define the types of host and guest atoms
  INTEGER :: u, ierror
  !===============
  ! Initialization
  !===============
  i =0
  iatom = 0
  jatom = 0
  u = 9        ! UNIT for WRITING
  ierror = 2   

  !=============================================================================
  ! To find the right group for each D atom
  ! Purpose: Assign a molecular ID to each D atom, based on the distance between
  ! atoms' images.
  !=============================================================================
  ! To determine the Radium cutoff for the host-guest pair
  if (TRIM(host_atom)=="O" .AND. TRIM(guest_atom)=="H") then
    R_cutoff = 1.02                 ! in Angstrom 
    R_cutoff_square = R_cutoff**2
  else if (TRIM(host_atom)=="H" .AND. TRIM(guest_atom)=="X") then
    R_cutoff = 0.6                  ! in Angstrom 
    R_cutoff_square = R_cutoff**2
  endif
  
  ! clustering
  outer: do i =1,n_samples
    inner: DO iatom = 1, nat
      hostatom: if (TRIM(wannier_center_info(iatom, i)%wannier_center_name) == TRIM(host_atom)) then
        inner2: DO jatom = 1, nat
          guestatom: if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
            ! define logical varibles
            if (                                 &                                                    
              (wannier_center_info(jatom, i)%image_coord(1) - wannier_center_info(iatom, i)%coord(1))**2  &
              +(wannier_center_info(jatom, i)%image_coord(2) - wannier_center_info(iatom, i)%coord(2))**2 &
              +(wannier_center_info(jatom, i)%image_coord(3) - wannier_center_info(iatom, i)%coord(3))**2 &
              < R_cutoff_square               &
               ) then
              wannier_center_info(jatom, i)%molecular_id = wannier_center_info(iatom, i)%molecular_id 
              !===========
              !For testing
              !===========
              !WRITE(*,*) wannier_center_info(jatom, i)%wannier_center_name, wannier_center_info(jatom, i)%molecular_id
            endif
          endif guestatom
        ENDDO inner2 
      endif hostatom 
    ENDDO inner
  enddo outer

  !===========
  !For testing
  !===========
  OPEN (UNIT=u,FILE='test_output_'//TRIM(guest_atom)//'_molecular_id.dat',  &
        STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)
    WRITE(u, *) ' Guest_atom ', ' molecular ID '
  do i = 1, 1 
    guest: do jatom = 1, nat 
      if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
        WRITE(u,*) wannier_center_info(jatom, i)%wannier_center_name, &
                   wannier_center_info(jatom, i)%molecular_id
      end if
    end do guest
  end do

END SUBROUTINE cluster 
