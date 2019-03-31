MODULE traj
!
! Purpose: 
! To declear data related to the simulation and traj.
IMPLICIT NONE
!INTEGER :: n_samples  !n_samples = INT(nmo/ns)
!INTEGER :: nmo_start, nmo_end  ! To get the total number of moves
INTEGER :: i_sample, i_input,i ! dummy index
CONTAINS
  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
    !
    ! Purpose:
    !  To calculate the total numbers of samples one want to include in their analysis.
    ! Data dictionary
    INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves

    write(*,*) 'In Function sampling_number, mo_end:', nmo_end
    ! no. of samples = INT({no. of moves}/ns)
    positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
      write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
    ELSE IF (nmo_end < nmo_start) THEN
      write(*,*) 'Please note that starting step shoud not larger than  ending step.'
    ELSE IF (ns ==0) THEN
      sampling_number = nmo_end-(nmo_start-1)
    ELSE IF (nmo_end-(nmo_start-1) <= ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns + 1)
    ELSE IF (nmo_end-(nmo_start-1) > ns) THEN
      sampling_number = INT((nmo_end-(nmo_start-1))/ns)
    END IF positive
  END FUNCTION sampling_number

  SUBROUTINE read_traj(indx,nmo_start,nmo_end,ns,nat,n_samples)
    !
    ! Purpose:
    ! To read info from the trajectory file (format: ***.xyz)
    USE atom_module
    USE parameter_shared

    INTEGER :: iatom, imovie, i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    !TYPE(atom), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: atom_info
    i_sample = 1
    outer: DO imovie=1,nmo_end-nmo_start
      read(indx,*)!Neglect data of this line
      read(indx,120) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
      120 FORMAT (5X,I8,9X,F12.3,6X,F20.10)
      write(*,*) 'the step:', imovie
      inner: do iatom= 1,nat
        read (indx,*) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
          atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
        !WRITE (*,*) & 
        !atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
        !atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
      enddo inner
      ! ns should be non-negative
      IF (ns < 0) THEN
        WRITE(*,*) 'Please note that ns should be non-negative.'
      ELSE IF (ns == 0) THEN
        CYCLE
      ELSE
        CALL skip_lines(indx, (nat+2)*(ns-1)) ! use CALL to run a SUBROUTINE
      ENDIF 

      i_sample = i_sample + 1

      ! To check if the sampling is finished
      check: IF (i_sample > n_samples) THEN 
        WRITE (*,*)'The total number of sample points are: ', n_samples
        EXIT
      END IF check
    ENDDO outer

  END SUBROUTINE read_traj

  !similar to read_traj, but not extracting sampled_movie(i_sample), sampled_time(i_sample) or sampled_energy(i_sample)
  SUBROUTINE read_traj_wannier_center(indx,nmo_start,nmo_end,ns,nat,n_samples)
    !
    ! Purpose:
    ! To read info from the ion+center file (format: ion+center.xyz)
    USE wannier_center_module
    USE parameter_shared

    INTEGER :: iatom, imovie, i_sample
    INTEGER, INTENT(IN) :: nat
    INTEGER, INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)
    INTEGER, INTENT(IN) :: indx
    INTEGER, INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
    INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
    i_sample = 1
    outer: DO imovie=1,nmo_end-nmo_start
      read(indx,*)!Neglect data of this line
      read(indx,*)!Neglect data of this line 
      write(*,*) 'the step:', imovie
      inner: do iatom= 1,nat
        read (indx,*) wannier_center_info(iatom, i_sample)%wannier_center_name, wannier_center_info(iatom,i_sample)%coord(1), & 
          wannier_center_info(iatom,i_sample)%coord(2), wannier_center_info(iatom,i_sample)%coord(3)
        ! Initialize the image_coord attribute
        wannier_center_info(iatom,i_sample)%image_coord(1) = wannier_center_info(iatom,i_sample)%coord(1) 
        wannier_center_info(iatom,i_sample)%image_coord(2) = wannier_center_info(iatom,i_sample)%coord(2)
        wannier_center_info(iatom,i_sample)%image_coord(3) = wannier_center_info(iatom,i_sample)%coord(3)
        ! Initialize the charge
        charge: if(TRIM(wannier_center_info(iatom, i_sample)%wannier_center_name) == "O") then 
          wannier_center_info(iatom,i_sample)%charge = 6 
        else if(TRIM(wannier_center_info(iatom, i_sample)%wannier_center_name) == "H") then 
          wannier_center_info(iatom,i_sample)%charge = +1 
        else
          wannier_center_info(iatom,i_sample)%charge = -2
        end if charge
      enddo inner
      ! ns should be non-negative
      IF (ns < 0) THEN
        WRITE(*,*) 'Please note that ns should be non-negative.'
      ELSE IF (ns == 0) THEN
        CYCLE
      ELSE
        CALL skip_lines(indx, (nat+2)*(ns-1)) ! use CALL to run a SUBROUTINE
      ENDIF 

      i_sample = i_sample + 1

      ! To check if the sampling is finished
      check: IF (i_sample > n_samples) THEN 
        WRITE (*,*)'The total number of sample points are: ', n_samples
        EXIT
      END IF check
    ENDDO outer

  END SUBROUTINE read_traj_wannier_center

  SUBROUTINE skip_lines(indx, i_input)
    !
    ! Purpose: 
    ! To skip lines when read data from the input
    IMPLICIT NONE
    INTEGER :: i
    INTEGER,INTENT(IN) :: i_input,indx
    do i=1,i_input
       read(indx,*) !Neglect (nat+2)*(ns-1) lines
    enddo    
  END SUBROUTINE skip_lines

END MODULE traj
