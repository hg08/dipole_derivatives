SUBROUTINE coord_translation(guest_atom, host_atom, nat, n_samples, a, b, c)
  
  ! 0a) Purpose: 
  ! The subroutine  coord_transform implement the coordinate translation for the atoms and wannier centers, considering the periodic
  ! boundary condition (PBC).  
  ! The subroutine 1c2_coord_translation is different from the subroutine 1c_coord_translation in the definition of D1,D2, M1,M2, N1,
  ! N2. The new definition is two times larger that the original ones defined in the subroutine 1c_coord_translation.
  ! 0b) Record of revisions:
  !    Date             Programmer                  description of Change
  !    ====             ==========                  ====================
  ! 2019.03.28          Gang Huang                  version 1.2 for coord_translation

  !==============
  !1) Declaration
  !==============
  USE wannier_center_module
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  REAL(kind=4),PARAMETER :: R_OD_c = 1.17           ! in Angstrom
  INTEGER ::i, iatom, jatom
  INTEGER,INTENT(IN) :: nat                         ! number of atoms
  INTEGER, INTENT(IN) :: n_samples                  ! n_samples = INT(nmo/ns)
  REAL(kind=4), INTENT(IN) :: a, b, c               ! Parameters of the simulation box 
  CHARACTER(2), INTENT(IN) :: host_atom, guest_atom ! to define the types of host and guest atoms
  LOGICAL :: A1, A2, B1, B2, C1, C2, A0, B0, C0
  INTEGER :: u, ierror
  !===============
  ! Initialization
  !===============
  i =0
  iatom = 0
  jatom = 0
  u = 8        ! UNIT for WRITING
  ierror = 2   
  !===========================================================================
  ! Translate the coordinates of gust atoms
  ! Purpose: to consider the PBC when calculate the distance between two atoms
  !===========================================================================
  outer: do i =1,n_samples
    inner: DO iatom = 1, nat
      oxygen: if (TRIM(wannier_center_info(iatom, i)%wannier_center_name) == TRIM(host_atom)) then
        inner2: DO jatom = 1, nat
          hydrogen: if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
              ! define logical varibles
              A1 =  wannier_center_info(iatom, i)%coord(1) < (-a/2)+R_OD_c  &
              .AND. wannier_center_info(iatom, i)%coord(1) > (-a/2)-R_OD_c  &
              .AND. wannier_center_info(jatom, i)%coord(1) < (a/2)+R_OD_c   &
              .AND. wannier_center_info(jatom, i)%coord(1) > (a/2)-R_OD_c  
              B1 =  wannier_center_info(iatom, i)%coord(2) < (-b/2)+R_OD_c  &
              .AND. wannier_center_info(iatom, i)%coord(2) > (-b/2)-R_OD_c  &
              .AND. wannier_center_info(jatom, i)%coord(2) < (b/2)+R_OD_c   &
              .AND. wannier_center_info(jatom, i)%coord(2) > (b/2)-R_OD_c
              C1 =  wannier_center_info(iatom, i)%coord(3) < (-c/2)+R_OD_c  &
              .AND. wannier_center_info(iatom, i)%coord(3) > (-c/2)-R_OD_c  &
              .AND. wannier_center_info(jatom, i)%coord(3) < (c/2)+R_OD_c   &
              .AND. wannier_center_info(jatom, i)%coord(3) > (c/2)-R_OD_c 
              A2 =  wannier_center_info(jatom, i)%coord(1) < (-a/2)+R_OD_c  &
              .AND. wannier_center_info(jatom, i)%coord(1) > (-a/2)-R_OD_c  &
              .AND. wannier_center_info(iatom, i)%coord(1) < (a/2)+R_OD_c   &
              .AND. wannier_center_info(iatom, i)%coord(1) > (a/2)-R_OD_c  
              B2 =  wannier_center_info(jatom, i)%coord(2) < (-b/2)+R_OD_c  &
              .AND. wannier_center_info(jatom, i)%coord(2) > (-b/2)-R_OD_c  &
              .AND. wannier_center_info(iatom, i)%coord(2) < (b/2)+R_OD_c   &
              .AND. wannier_center_info(iatom, i)%coord(2) > (b/2)-R_OD_c   
              C2 =  wannier_center_info(jatom, i)%coord(3) < (-c/2)+R_OD_c  &
              .AND. wannier_center_info(jatom, i)%coord(3) > (-c/2)-R_OD_c  &
              .AND. wannier_center_info(iatom, i)%coord(3) < (c/2)+R_OD_c   &
              .AND. wannier_center_info(iatom, i)%coord(3) > (c/2)-R_OD_c 
              A0 = (.NOT. A1) .AND. (.NOT. A2)
              B0 = (.NOT. B1) .AND. (.NOT. B2)
              C0 = (.NOT. C1) .AND. (.NOT. C2)
              !==============================
              !For testing
              !WRITE(*,*) "A0B0C0:", A0,B0,C0
              !==============================
            !=====================================================================================
            ! There are 27 cases for the configuration of the host and guest coordinates. They are:
            ! Class 1: A0B0C0 (1 case)
            ! Class 2: A0B0C1, A0B0C2, A0B1C0, A0B2C0, A1B0C0,A2B0C0 (6 cases)
            ! Class 3: A0B1C1,A0B1C2,A0B2C1,A0B2C2,
            !          A1B0C1,A1B0C2,A2B0C1,A2B0C2,
            !          A1B1C0,A1B2C0,A2B1C0,A2B2C0 (12 cases)
            ! Class 4: A1B1C1,A1B1C2,A1B2C1,A1B2C2,
            !          A2B1C1,A2B1C2,A2B2C1,A2B2C2 (8 cases)
            !=====================================================================================
            ! Class 1 (1/1): the case of A0B0C0 (A0B0C0 is the most important part among all the 27 cases)  
            if (                                 &                               
              A0 .AND. B0 .AND. C0               &
               ) then
               ! In the case of A0B0C0, the image coordinates of the guest atom will be the same as the original coordinates, ie. the coordinate will 
               ! directly be used to calculate the distance between the host and guest atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) 
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3)  
            ! Class 2 (1/6): the case of A0B0C1 
            elseif (                            & 
              A0 .AND. B0 .AND. C1              & 
               ) then
               ! In the case of A0B0C1, the image coordinates of the guest atom will be the same as the original coordinates for x
               ! and y, ie. the coordinate x,y will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinate z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) 
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c  
              !===========
              !For testing
              !===========
              WRITE(*,*) "5 image x: ", wannier_center_info(5, 1)%image_coord(1), wannier_center_info(5, 1)%coord(1)  
              WRITE(*,*) "image y: ", wannier_center_info(5, 1)%image_coord(2), wannier_center_info(5, 1)%coord(2)
              WRITE(*,*) "image z: ", wannier_center_info(5, 1)%image_coord(3), wannier_center_info(5, 1)%coord(3)

            ! Class 2 (2/6): the case of A0B0C2 
            elseif (                             &
              A0 .AND. B0 .AND. C2               &
               ) then
               ! In the case of A0B0C2, the image coordinates of the guest atom will be the same as the original coordinates for x
               ! and y, ie. the coordinate x,y will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinate z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) 
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c  
            ! Class 2 (3/6): the case of A0B1C0
            elseif (                             &
              A0 .AND. B1 .AND. C0               &
               ) then
               ! In the case of A0B1C0, the image coordinates of the guest atom will be the same as the original coordinates for x
               ! and z, ie. the coordinate x and z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinate y will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3)   
            ! Class 2 (4/6): the case of A0B2C0
            elseif (                             &
              A0 .AND. B2 .AND. C0               &
               ) then
               ! In the case of A0B2C0, the image coordinates of the guest atom will be the same as the original coordinates for x
               ! and z, ie. the coordinate x and z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinate y will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3)   
            ! Class 2 (5/6): the case of A1B0C0
            elseif (                             &
              A1 .AND. B0 .AND. C0               &
               ) then
               ! In the case of A1B0C0, the image coordinates of the guest atom will be the same as the original coordinates for y
               ! and z, ie. the coordinate y and z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinate x will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) 
            ! Class 2 (6/6): the case of A2B0C0
            elseif (                             &
              A2 .AND. B0 .AND. C0               &
               ) then
               ! In the case of A2B0C0, the image coordinates of the guest atom will be the same as the original coordinates for y
               ! and z, ie. the coordinate y and z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinate x will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) 
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3)   
            ! Class 3 (1/12): the case of A0B1C1
            elseif (                             &
              A0 .AND. B1 .AND. C1               &
               ) then
               ! In the case of A0B1C1, the image coordinates of the guest atom will be the same as the original coordinates for x,
               ! ie. the coordinate x will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 3 (2/12): the case of A0B1C2
            elseif (                             & 
              A0 .AND. B1 .AND. C2               &
               ) then
               ! In the case of A0B1C2, the image coordinates of the guest atom will be the same as the original coordinates for x,
               ! ie. the coordinate x will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 3 (3/12): the case of A0B2C1
            elseif (                            &
              A0 .AND. B2 .AND. C1              &
               ) then
               ! In the case of A0B2C1, the image coordinates of the guest atom will be the same as the original coordinates for x,
               ! ie. the coordinate x will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 3 (4/12): the case of A0B2C2
            elseif (                            &
              A0 .AND. B2 .AND. C2              &
               ) then
               ! In the case of A0B2C2, the image coordinates of the guest atom will be the same as the original coordinates for x,
               ! ie. the coordinate x will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) 
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 3 (5/12): the case of A1B0C1
            elseif (                            &
              A1 .AND. B0 .AND. C1              &
               ) then
               ! In the case of A1B0C1, the image coordinates of the guest atom will be the same as the original coordinates for y,
               ! ie. the coordinate y will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2)
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 3 (6/12): the case of A1B0C2
            elseif (                            &
              A1 .AND. B0 .AND. C2              &
               ) then
               ! In the case of A1B0C2, the image coordinates of the guest atom will be the same as the original coordinates for y,
               ! ie. the coordinate y will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2)
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 3 (7/12): the case of A2B0C1
            elseif (                            &
              A2 .AND. B0 .AND. C1              &
               ) then
               ! In the case of A2B0C1, the image coordinates of the guest atom will be the same as the original coordinates for y,
               ! ie. the coordinate y will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2)
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 3 (8/12): the case of A2B0C2
            elseif (                            &
              A2 .AND. B0 .AND. C2              &
               ) then
               ! In the case of A2B0C2, the image coordinates of the guest atom will be the same as the original coordinates for y,
               ! ie. the coordinate y will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2)
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 3 (9/12): the case of A1B1C0
            elseif (                            &
              A1 .AND. B1 .AND. C0              &
               ) then
               ! In the case of A1B1C0, the image coordinates of the guest atom will be the same as the original coordinates for z,
               ! ie. the coordinate z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and y will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) 
            ! Class 3 (10/12): the case of A1B2C0
            elseif (                            &
              A1 .AND. B2 .AND. C0              &
               ) then
               ! In the case of A1B2C0, the image coordinates of the guest atom will be the same as the original coordinates for z,
               ! ie. the coordinate z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and y will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) 
            ! Class 3 (11/12): the case of A2B1C0
            elseif (                            &
              A2 .AND. B1 .AND. C0              &
               ) then                            
               ! In the case of A2B1C0, the imagecoordinates of the guest atom will be the same as the original coordinates for z,
               ! ie. the coordinate z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and y will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) 
            ! Class 3 (12/12): the case of A2B2C0
            elseif (                            &
              A2 .AND. B2 .AND. C0              &
               ) then
               ! In the case of A2B2C0, the image coordinates of the guest atom will be the same as the original coordinates for z,
               ! ie. the coordinate z will directly be used to calculate the distance between the host and the guest atom,
               ! while the coordinates x and y will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) 
            ! Class 4 (1/8): the case of A1B1C1
            !===========================================================================================================================
            ! My code is based on this pseudo-code, I still keep one example for A1B1C1 case here.
            ! pseudo-code for the case of A1B1C1: the "ELSEIF(...)THEN" part
            ! elseif (
            !       x(host) < (-a/2)+R_OD_c .and. x(host) > (-a/2)-R_OD_c .AND. x(guest) < a/2+R_OD_c .and. x(host) > a/2-R_OD_c .... A1
            ! .AND. y(host) < (-b/2)+R_OD_c .and. y(host) > (-b/2)-R_OD_c .AND. y(guest) < b/2+R_OD_c .and. y(host) > b/2-R_OD_c .... B1
            ! .AND. z(host) < (-c/2)+R_OD_c .and. z(host) > (-c/2)-R_OD_c .AND. z(guest) < c/2+R_OD_c .and. z(host) > c/2-R_OD_c .... C1 
            !        ) then
            !===========================================================================================================================
            elseif (                             &
              A1 .AND. B1 .AND. C1               &
              ) then
              ! In the case of A1B1C2, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
              ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c 
              !======================================
              ! The transformation written as formula
              ! x_D' = x_D - a
              ! y_D' = y_D - b
              ! z_D' = z_D - c
              !======================================
            ! Class 4 (2/8): the case of A1B1C2
            elseif (                             &
              A1 .AND. B1 .AND. C2               &
               ) then
               ! In the case of A1B1C2, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 4 (3/8): the case of A1B2C1
            elseif (                            &
              A1 .AND. B2 .AND. C1              &
               ) then
               ! In the case of A1B2C1, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 4 (4/8): the case of A1B2C2
            elseif (                            &
              A1 .AND. B2 .AND. C2              &
               ) then
               ! In the case of A1B2C2, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) - a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 4 (5/8): the case of A2B1C1
            elseif (                            &
              A2 .AND. B1 .AND. C1              &
               ) then
               ! In the case of A2B1C1, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 4 (6/8): the case of A2B1C2
            elseif (                            &
              A2 .AND. B1 .AND. C2              &
               ) then
               ! In the case of A2B1C2, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) - b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            ! Class 4 (7/8): the case of A2B2C1
            elseif (                            &
              A2 .AND. B2 .AND. C1              &
               ) then
               ! In the case of A2B1C2, the coordinates x, y and z will be transformed before it is used to calculate the distance between the host and the gust
               ! atoms.
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) - c
            ! Class 4 (8/8): the case of A2B2C2
            else
              wannier_center_info(jatom, i)%image_coord(1) = wannier_center_info(jatom, i)%coord(1) + a
              wannier_center_info(jatom, i)%image_coord(2) = wannier_center_info(jatom, i)%coord(2) + b
              wannier_center_info(jatom, i)%image_coord(3) = wannier_center_info(jatom, i)%coord(3) + c
            endif
          endif hydrogen
        ENDDO inner2 
      endif oxygen
    ENDDO inner
  enddo outer

  !===========
  !For testing
  !===========
  !WRITE(*,*) 'molecular id:', wannier_center_info(1, 1)%molecular_id
  !WRITE(*,*) 'transformed coord:', wannier_center_info(1, 1)%image_coord(1)
  !WRITE(*,*) 'molecular id:', wannier_center_info(1, 2)%molecular_id
  !WRITE(*,*) 'transformed coord:', wannier_center_info(1, 2)%image_coord(1)
  OPEN (UNIT=u,FILE='TEST_OUTPUT_'//TRIM(guest_atom)//'.dat',  &
    STATUS='REPLACE',ACTION='WRITE',IOSTAT=ierror)
    WRITE(u, *) 'Coord. of Guest_atom ', ' image_coord ', '  coord  ', '  difference '
  do i = 1, n_samples 
    guest: do jatom = 1, nat 
      if (TRIM(wannier_center_info(jatom, i)%wannier_center_name) == TRIM(guest_atom)) then
        WRITE(u,*) wannier_center_info(jatom, i)%wannier_center_name, ' x ', &
                  wannier_center_info(jatom, i)%image_coord(1), wannier_center_info(jatom, i)%coord(1), &
                  wannier_center_info(jatom, i)%image_coord(1) - wannier_center_info(jatom, i)%coord(1)
        WRITE(u,*) wannier_center_info(jatom, i)%wannier_center_name, ' y ', &
                  wannier_center_info(jatom, i)%image_coord(2), wannier_center_info(jatom, i)%coord(2), &
                  wannier_center_info(jatom, i)%image_coord(2) - wannier_center_info(jatom, i)%coord(2)
        WRITE(u,*) wannier_center_info(jatom, i)%wannier_center_name, ' z ', &
                  wannier_center_info(jatom, i)%image_coord(3), wannier_center_info(jatom, i)%coord(3), &
                  wannier_center_info(jatom, i)%image_coord(3) - wannier_center_info(jatom, i)%coord(3)
      end if
    end do guest
  end do

END SUBROUTINE coord_translation 
