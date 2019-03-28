!=========================================================================
! The following algorithm is written before I starting writing a real code
! It is important to divide the problem into several piece of smaller ones.
!   Date           Programer                  Description
!  ======         ==========                  ===========
! 20190327         Gang Huang                 Steps to do
!=========================================================================
! 0) Clearly describe the problem
! 1) Define the input and output
! input: the trajectory file ion+center.xyz.
! input file : ion+center.xyz
! output: dipole moments of all the water molecules for all the sampled time points. It is an array with the shape (STEPS * N * 3),
! where STEPS is the number of sampled time points, N is the number of molecules, 3 means there are 3 coordinates in 3D space.
! output file: dipole_moment.dat
!
! 2) Design the algorithm which should be implemented in this program.
! 2a) For any atom, find out which molecule is it belong to 
!  2a1) host molecule's symbol: symb_host (indics of molecules must be defined), define a molecular ID number
!  2a2) For any atom, find out its host (group), we have to define a attribute for atom 
!  2a2a) PBC is considered before define the distance between the host atom and the guest atom
!       if (x(host) < (-a/2)+R_OD_c .and. x(host) > (-a/2)    .AND.      x(guest) < a/2 .and. x(host) > a/2-R_OD_c    .... A1
!　　.AND.  y(host) < (-b/2)+R_OD_c .and. y(host) > (-b/2)    .AND.      y(guest) < b/2 .and. y(host) > b/2-R_OD_c    .... B1
!　　.AND.  z(host) < (-c/2)+R_OD_c .and. z(host) > (-c/2)    .AND.      z(guest) < c/2 .and. z(host) > c/2-R_OD_c    .... C1 ) then
!    x_D' = x_D - a
!    y_D' = y_D - b
!    z_D' = z_D - c
!        elseif (x(guest) < (-a/2)+R_OD_c .and. x(host) > (-a/2) .AND.    x(host) < a/2 .and. x(host) > a/2-R_OD_c    .... A1
!　　.AND.  y(guett) < (-b/2)+R_OD_c .and. y(host) > (-b/2)    .AND.      y(host) < b/2 .and. y(host) > b/2-R_OD_c    .... B1
!　　.AND.  z(guest) < (-c/2)+R_OD_c .and. z(host) > (-c/2)    .AND.      z(host) < c/2 .and. z(host) > c/2-R_OD_c    .... C1 ) then
!    x_D' = x_D + a
!    y_D' = y_D + b
!    z_D' = z_D + c
!        elseif (x(guest) < (-a/2)+R_OD_c .and. x(host) > (-a/2) .AND.    x(host) < a/2 .and. x(host) > a/2-R_OD_c    .... A1
!　　.AND.  y(guett) < (-b/2)+R_OD_c .and. y(host) > (-b/2)    .AND.      y(host) < b/2 .and. y(host) > b/2-R_OD_c    .... B1
!　　.AND.  z(guest) < (-c/2)+R_OD_c .and. z(host) > (-c/2)    .AND.      z(host) < c/2 .and. z(host) > c/2-R_OD_c    .... C1 ) then
!    x_D' = x_D + a
!    y_D' = y_D + b
!    z_D' = z_D + c
!        elseif ( (.NOT. A1) .AND. (.NOT. A2) ... A0
!        .AND.    (.NOT. B1) .AND. (.NOT. B2) ... B0
!        .AND.    (.NOT. C1) .AND. (.NOT. C2) ... C0 ) then 
!    x_D' = x_D 
!    y_D' = y_D 
!    z_D' = z_D 
!         ...
!      endif
      
! 2b) For and wannier center, find out which atom is it belong to and therefore, find out the host of the wannier center, i.e., to
! find out which molecule is this wannier center belong to
!  2b1) PBC is also considered before define the distance between the host atom (eg, D atom) and the guest atom (X, ie., wannier
!  center). This step is the same to the 2a2a-th step.
! 3) Calculate the dipole moment for each water molecule
!     $ D = \sum_I R_I * Q_I $
! 5) Translate the algorithm to Fortran codes
! 6) Test the codes.
! 7) Finish up.
