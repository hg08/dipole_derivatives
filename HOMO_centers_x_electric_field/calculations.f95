SUBROUTINE sample(pos_filename,nmo_start,nmo_end,nat,ns,str_of_center_atoms,n_samples)
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) The subroutine sample.f95 reduce the size of the trajectory. 
  ! 0b) Date     Programmer     version
  !  20181225    Gang Huang     3c
  !
  ! Declaration
  USE parameter_shared
  USE atom_module
  USE wannier_center_module
  USE traj
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  integer,parameter :: rk=4              
  INTEGER :: iatom, imovie
  INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
  INTEGER,INTENT(IN) :: nmo_start  
  INTEGER,INTENT(IN) :: nmo_end  
  INTEGER,INTENT(IN) :: nat ! number of atoms
  INTEGER, INTENT(INOUT) :: n_samples  !n_samples = INT(nmo/ns)
  character(LEN=*), INTENT(IN) :: pos_filename
  !INTEGER,INTENT(IN) :: num_of_kind_center_atoms
  CHARACTER(LEN=*),INTENT(IN) :: str_of_center_atoms
  
  ! Initialization
  iatom = 0
  imovie =0
  i =0

  !=======================
  !read in trajectory file 
  !=======================
  open(10,file=trim(pos_filename))
  ! READ the first two steps to get information of the system
  ! Now starting read data
  n_samples = sampling_number(nmo_start,nmo_end,ns)
  allocate(wannier_center_info(nat,n_samples))
  allocate(sampled_movie(n_samples))
  allocate(sampled_time(n_samples))
  allocate(sampled_energy(n_samples))
  to_skip: IF (nmo_start == 0) THEN
    CALL read_traj_wannier_center(10,nmo_start,nmo_end,ns,nat,n_samples)
  ELSE IF (nmo_start > 0 .AND. nmo_start < nmo_end) THEN
    CALL skip_lines(10, (nat+2)*(nmo_start-1)) ! use CALL to rerun a SUBROUTINE
    CALL read_traj_wannier_center(10,nmo_start,nmo_end,ns,nat,n_samples)
  ELSE
    WRITE(*,*) 'Please enter a right value of starting step'
  END IF to_skip

  close(10)
  write(6,*) 'End of trajectory reading.'
  !=============
  !write in file
  !=============
  sampled_pos_filename = 'traj_pos_sampled.dat'
  open(10,file=sampled_pos_filename)
  do i =1,n_samples
    write (10,'(I8)') nat
    WRITE(10,100) 'i =',i-1,', time =',sampled_time(i),', E =',sampled_energy(i)
    100 FORMAT (1X,A3,I10,A8,F10.3,A5,F20.10)
    DO iatom = 1, nat
       WRITE(10,*) TRIM(wannier_center_info(iatom, i)%wannier_center_name), &
         wannier_center_info(iatom,i)%coord(1), &
         wannier_center_info(iatom,i)%coord(2), &
         wannier_center_info(iatom,i)%coord(3)
    ENDDO
  enddo
  write(6,*)'Sampled trajectory is written: ',sampled_pos_filename
  close(10)

  deallocate(sampled_movie, sampled_time,sampled_energy)
  
  WRITE(*,*) 'list_of_center_atoms: ', str_of_center_atoms
  WRITE(*,*) 'Len of list_of_center_atoms: ', LEN(str_of_center_atoms)
END SUBROUTINE sample

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
  ! Purpose:
  !   To return the minimum of the distance between the image of a host-gust pair
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

SUBROUTINE calculate_dipole_moment(host_atom, nat, n_samples)  
  ! 0a) Purpose: 
  ! The subroutine calculate_dipole_moment is used to calculate dipole moment for each water molecule (D2O)  
  !
  ! 0b) Record of revisions:
  !    Date             Programmer                  description of Change
  !    ====             ==========                  ====================
  ! 2019.03.29          Gang Huang                  Original code 

  !==============
  !1) Declaration
  !==============
  USE wannier_center_module
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  REAL(kind=4), ALLOCATABLE, DIMENSION(:,:,:) :: mu             ! type: 3 * n_molecule * n_samples

  INTEGER ::i, iatom, jatom
  INTEGER :: indx                                   ! index of molecules
  INTEGER :: istat
  INTEGER,INTENT(IN) :: nat                         ! number of atoms
  INTEGER, INTENT(IN) :: n_samples                  ! n_samples = INT(nmo/ns)
  CHARACTER(2), INTENT(IN) :: host_atom             ! to define the types of host and guest atoms
  INTEGER :: u, ierror
  REAL(kind=4) :: total_charge
  !===============
  ! Initialization
  !===============
  i =0
  iatom = 0 
  jatom = 0
  u = 8        ! UNIT for WRITING
  ierror = 2   
  istat = 0
  total_charge = 0.0
  !============================
  ! Calculate the dipole moment
  !============================
  
  indx = 0
  ALLOCATE(mu(n_samples, 32, 3),STAT=istat)
  mu(:,:,:) = 0.0
 

  !==================================================================
  ! For testing
  ! Open a file to write down the molecular charge for each molecule.
  OPEN (UNIT=u,FILE='test_output_molecular_charge.dat',STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)
  !Loop the trajectory steps 
  do i = 1, n_samples 
    indx = 0
    Oxygen: do iatom = 1, nat 
      O: if (TRIM(wannier_center_info(iatom, i)%wannier_center_name) == TRIM(host_atom)) then
        indx = indx + 1
        mu(i,indx,1) = wannier_center_info(iatom, i)%coord(1) * wannier_center_info(iatom, i)%charge  
        mu(i,indx,2) = wannier_center_info(iatom, i)%coord(2) * wannier_center_info(iatom, i)%charge
        mu(i,indx,3) = wannier_center_info(iatom, i)%coord(3) * wannier_center_info(iatom, i)%charge
        total_charge = wannier_center_info(iatom, i)%charge
        guest3: do jatom = 1, nat
          if ( TRIM(wannier_center_info(jatom, i)%wannier_center_name) .NE. TRIM(host_atom) .AND.      &
            wannier_center_info(jatom, i)%molecular_id == wannier_center_info(iatom, i)%molecular_id & 
             ) then
            mu(i,indx,1) = mu(i,indx,1) + wannier_center_info(jatom, i)%image_coord(1) * wannier_center_info(jatom, i)%charge  
            mu(i,indx,2) = mu(i,indx,2) + wannier_center_info(jatom, i)%image_coord(2) * wannier_center_info(jatom, i)%charge
            mu(i,indx,3) = mu(i,indx,3) + wannier_center_info(jatom, i)%image_coord(3) * wannier_center_info(jatom, i)%charge
            total_charge = total_charge + wannier_center_info(jatom, i)%charge
          endif
          ! For testing 
          ! The next line can be deleted
          WRITE(u,*) "total charge:",iatom, jatom, total_charge
        enddo guest3
      end if O
    end do Oxygen
  end do
  CLOSE(UNIT=u)

  !===========
  !For testing
  !===========
  OPEN (UNIT=u,FILE='test_output_molecular_dipole_moments.dat',  &
        STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)
    WRITE(u, *) ' Molecule ID ', ' x ', ' y ', ' z ', 'Dipole moment'
  do i = 1, n_samples
    do indx = 1, 32 
      WRITE(u,*) indx,  mu(i,indx,1), mu(i,indx,2), mu(i,indx,3), SQRT(mu(i,indx,1)**2 + mu(i,indx,2)**2 +  mu(i,indx,3)**2)
    end do
  end do 
  CLOSE(u)

  DEALLOCATE(mu)
END SUBROUTINE calculate_dipole_moment
