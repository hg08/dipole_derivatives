SUBROUTINE read_input(delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,str_of_host_atoms) 
  !2018/12/27
  ! 0) Purpose:     
  ! 0a) Date     Programmer     version
  !  20181225    Gang Huang     3c
  ! 1a) The subroutine sample.f95 reduce the size of the trajectory. 
  !
  ! Declaration
  USE parameter_shared
  !USE atom_module
  USE traj
  IMPLICIT NONE

!==========
!parameters
!==========
  integer,parameter :: rk=8              
  INTEGER :: iatom, imovie
  INTEGER, INTENT(OUT) :: nat ! number of atoms
  !INTEGER, INTENT(OUT) :: n_samples  !n_samples = INT(nmo/ns)
  !TYPE(atom), INTENT(OUT), ALLOCATABLE, DIMENSION(:,:) :: atom_info
  !real,allocatable,dimension (:,:)   :: x,y,z
  REAL(rk), INTENT(OUT)  :: delta_t0  ! For reading data
  character(LEN=50), INTENT(OUT) :: filename
  character(LEN=50), INTENT(OUT) :: pos_filename
  CHARACTER(LEN=2),INTENT(OUT) :: str_of_host_atoms
  INTEGER, INTENT(OUT) :: nmo_start 
  INTEGER, INTENT(OUT) :: nmo_end
  INTEGER, INTENT(OUT) :: ns
  !REAL(rk), INTENT(OUT)  :: radium_OD_c  ! This value is defined from the O-D bond length in a D2O molecule
  !CHARACTER(LEN=1),INTENT(OUT) :: axis ! For define the normal direction
  

  ! Initialization
  iatom = 0
  imovie =0
  i =0

!==================
!read data in input
!==================
  write(6,*)'What is the timestep (ps):'
  read(5,*)delta_t0
  write(6,*)'What is the name of the system:'
  read(5,*)filename
  write(6,*)'What is the name of the trajecotry file:'
  read(5,*)pos_filename     
  write(6,*)'What is the initial step of the trajecotry:'
  read(5,*)nmo_start !value of the first movie steps
  write(6,*)'What is the end step of the trajecotry:'
  read(5,*)nmo_end !value of the last movie steps
  write(6,*)'What is the total number of atoms in the system:'
  read(5,*)nat!number of the total atoms.
  write(6,*)'What is the frequency for sampling? (Your sample the trajectory every ns step. What is the ns you want it to be?):'
  read(5,*)ns! [ns*0.0005] ps is the new time step for calculation,if delta_t0=0.0005 ps.
  WRITE (6,*)'What is the name of the host atom?(Eg. "I", "N", or "Na")'
  read(5,*) str_of_host_atoms
  !WRITE (6,*)'What is the radium cutoff of O-D bond in a D2O water molecule?'
  !read(5,*) radium_OD_c
  !WRITE (6,*)'What is the normal direction?(Eg. "x", "y", or "z")'
  !read(5,*)axis

END SUBROUTINE read_input 
