PROGRAM main 
  ! 0) Purpose:     
  ! To calculate the dipole derivatives.
  ! 0a) Date     Programmer     version
  !    =====     ==========     =======
  !  2019/03/27    Gang Huang     Original code
  ! EXPLAIN: host atom: For a water molecule D2O, the index of each O atom can characterize the water molecule, because there is
  ! only one O atom in each water molecule. In this case, we call O atom as the 'host atom' of the water molecule. 
  !
  !============
  ! Declaration
  !============
  USE parameter_shared
  USE atom_module
  USE wannier_center_module
  USE count_time
  USE traj  ! where we have difined n_samples
  USE types
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  integer,parameter :: rk=4              
  INTEGER :: iatom, imovie
  INTEGER :: ns  ! Get one sample from the trajectory every ns step.
  CHARACTER(2) :: atom_type
  INTEGER :: nmo_start
  INTEGER :: nmo_end
  INTEGER :: nat ! number of atoms
  INTEGER :: n_samples   ! n_samples = INT(nmo/ns)
  REAL(rk) :: delta_t0   ! For reading data
  REAL(rk) :: a,b,c      ! The parameters of simulation box
  character(LEN=50) :: filename, pos_filename
  CHARACTER(LEN=50) :: str_of_host_atoms
  CHARACTER(LEN=2) :: guest_atom, host_atom

  !===============
  ! Initialization
  !===============
  iatom = 0
  imovie =0
  i =0
  atom_type = ''
  a = 0.0 ! Angstrom
  b = 0.0 ! Angstrom
  c = 0.0 ! Angstrom
  guest_atom = "H"
  host_atom = "O"
  n_samples = 0

  call system_clock(begin_time,rat)
  
  !==================================================
  ! To read data: the required controlling parameters
  !==================================================
  CALL read_input(delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,str_of_host_atoms)
  
  !========================
  ! Sampling the trajectory
  !========================
  CALL sample(pos_filename,nmo_start,nmo_end,nat,ns,str_of_host_atoms,n_samples)
  ! After running the sample() subroutine, now we have the new sampled traj. file (stored in wannier_center_info).

  !====================================
  ! label each atom with a molecular id
  !====================================
  CALL label(nat,n_samples)
 
  !====================================================
  ! to finish coordinate translation of the guest atoms
  ! ie. D and X and to finich clustering.
  !====================================================
  ! Parameters for the simulation box
  a = 10.1 ! Angstrom
  b = 10.1 ! Angstrom
  c = 10.1 ! Angstrom
  ! coordinate translation and clustering for Deuterium atoms (D)
  host_atom = "O"
  guest_atom = "H"
  CALL coord_translation(guest_atom, host_atom, nat, n_samples, a, b, c)
  ! coordinate translation and clustering for Wannier centers (X)
  guest_atom = "X"
  CALL coord_translation(guest_atom, host_atom, nat, n_samples, a, b, c)

  ! Calculate dipole moments for molecules
  CALL calculate_oh_dipole_moment(host_atom, nat, n_samples)
  
  !================================================
  ! to write out the 3-D array ---- the linked list
  !================================================
  !ALLOCATE(host_atoms_array,STAT=istat)
  !WRITE(*,*)"istat =", istat
  !CALL position_of_atom(str_of_host_atoms, nat, n_samples, host_atoms_array)
  !WRITE(*,*) host_atoms_array%coord_array(1)
  !WRITE(*,*) host_atoms_array%coord_array(2)
  !WRITE(*,*) host_atoms_array%coord_array(3)
  !DEALLOCATE(host_atoms_array,STAT=istat)

  call system_clock(end_time,rat)
  write(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM main 
