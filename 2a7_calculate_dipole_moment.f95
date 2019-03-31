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
  !===============
  ! Initialization
  !===============
  i =0
  iatom = 0 
  jatom = 0
  u = 8        ! UNIT for WRITING
  ierror = 2   
  istat = 0
  !============================
  ! Calculate the dipole moment
  !============================
  
  indx = 0
  ALLOCATE(mu(n_samples, 32, 3),STAT=istat)
  mu(:,:,:) = 0.0
  
  do i = 1, n_samples 
    indx = 0
    Oxygen: do iatom = 1, nat 
      O: if (TRIM(wannier_center_info(iatom, i)%wannier_center_name) == TRIM(host_atom)) then
        indx = indx + 1
        mu(i,indx,1) = wannier_center_info(iatom, i)%coord(1) * wannier_center_info(iatom, i)%charge  
        mu(i,indx,2) = wannier_center_info(iatom, i)%coord(2) * wannier_center_info(iatom, i)%charge
        mu(i,indx,3) = wannier_center_info(iatom, i)%coord(3) * wannier_center_info(iatom, i)%charge
        guest3: do jatom = 1, nat
          if (wannier_center_info(jatom, i)%molecular_id == wannier_center_info(iatom, i)%molecular_id) then
            mu(i,indx,1) = mu(i,indx,1) + wannier_center_info(jatom, i)%image_coord(1) * wannier_center_info(jatom, i)%charge  
            mu(i,indx,2) = mu(i,indx,2) + wannier_center_info(jatom, i)%image_coord(2) * wannier_center_info(jatom, i)%charge
            mu(i,indx,3) = mu(i,indx,3) + wannier_center_info(jatom, i)%image_coord(3) * wannier_center_info(jatom, i)%charge
          endif
        enddo guest3
      end if O
    end do Oxygen
  end do

  !===========
  !For testing
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
