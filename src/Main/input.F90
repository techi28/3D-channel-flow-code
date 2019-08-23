!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	     INPUT FILE                               !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE input(No_vp,No_cp,nbe,nbev,xv,yv,zbv)

!---------------------------------------------------------------------!
!                                                                     !
!    This program reads the data files:                               !
!                  -   input_data.dat                                 !
!    corresponding to the data of the physical domain.                !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  |             |           |                                     |  !
!  | <-- No_vp   |(N_CELL0,3)| Index number of the 3 cell vertices |  !
!  | <-- No_cp   |(N_CELL,3) | Index number of surrounding 3 cells |  !
!  | <-- nbe     | N_CELL    | Type of boundary cell (inside or bc)|  !
!  | <-- xv      | N_VERT    | x-coordinate of the vertex          |  !
!  | <-- yv      | N_VERT    | y-coordinate of the vertex          |  !
!  | <-- zbv     | N_VERT    | The bottom of the river at vertex   |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  | IDISPLAY    | Display the data at the windows                 |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!


!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!

!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            implicit none
#     endif
!     =============== END ================
!     ====================================
!     ____________________________________
!    |                                    |
!    |     Declaration of variables       |
!    |____________________________________|

      integer,dimension(:,:):: No_vp
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe
      integer,dimension(:)  :: nbev
      real*8 ,dimension(:)  :: xv
      real*8 ,dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      integer,parameter :: Idisplay = 0
      character*80 title

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: input'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Read & write:  data.txt                         !
!                                                                     !
!*********************************************************************!

!---------------------------------------------------------------------!
!     The file "data.txt" has the following format:                   !
!                                                                     !
!         --------NUMBERING OF CELL VERTICES----                      !
!              1       2     323                                      !
!              2       3     323                                      !
!              3       4     324                                      !
!              ....                                                   !
!         --------COORDINATES OF CELL VERTICES----                    !
!              0.00000      -20.00000                                 !
!              4.68750      -20.00000                                 !
!              9.37500      -20.00000                                 !
!              ....                                                   !
!         --------NUMBERING OF CELL CENTERS----                       !
!              0       2    6401                                      !
!              0    6402       1                                      !
!              0    6403    6402                                      !
!              ....                                                   !
!         --------BOTTOM LEVELS ---------------                       !
!              0.000   0.000   0.000   0.000   0.000   ...            !
!              0.000   0.080   0.310   0.550   0.780   ...            !
!              2.190   2.420   2.660   2.890   3.130   ...            !
!              ....                                                   !
!         --------TYPE OF BOUNDARY CELLS ---------------              !
!              1      1      1      1      1      1    ...            !
!              1      1      1      1      1      1    ...            !
!              2      0      0      0      0      0    ...            !
!              ....                                                   !
!         --------TYPE OF BOUNDARY VERTEX-----------                  !
!                                                                     !
!---------------------------------------------------------------------!

71    format(A62)
72    format(A25)
73    format(A46,I10)
74    format(A46)
76    format(A46,F15.10)
77    format(A46,E15.5)

!      ________________________________________________________
!     |                                                        |
!     |                    Read data file                      |
!     |________________________________________________________|

      open(22,file='data.txt',status='old')

      read(22,74) title
      read(22,74) title
      read(22,74) title
!     ________________________________________________________
!     Read: Numbering of the cell vertices for each element

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,88) ((No_vp(i,j),j=1,3),i=1,N_CELL0)
88    format(3I8)
!     ________________________________________________________
!     Read: Coordinates of cell vertices

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,87) (xv(nv),yv(nv),nv=1,N_VERT)
87    format(2F15.5)
!     ________________________________________________________
!     Read: Numbering of the surrounding three cell centers

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,88) ((No_cp(i,j),j=1,3),i=1,N_CELL0)
!     ________________________________________________________
!     Read: Bottom levels

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,86) (zbv(i),i=1,N_VERT)
86    format(10F8.3)
!     ________________________________________________________
!     Read: Type of cell (inside or boundary)

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,85) (nbe(i),i=1,N_CELL0)
85    format(20I7)
!     ________________________________________________________
!     Read: Inputting for wall boundary vertices

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,85)  (nbev(i),i=1,N_VERT)

      close(22)

!      ________________________________________________________
!     |                                                        |
!     |              Displaying all the results                |
!     |________________________________________________________|
       N_BC = 0
       do nv=1,N_VERT
          if(nbev(nv) .gt. 0) N_BC = N_BC + 1
       enddo

!   --------------------------------------------------------
!     modify nbev value according to its position
       if(ChooseBoundary .eq. 0) then
	       do i=1,N_VERT
		    if(nbev(i) .eq. 1) then
                       nbev(i) = 1 ! top
		    elseif(nbev(i) .eq. 2) then
                       nbev(i) = 1 ! bottom
		    elseif(nbev(i) .eq. 3) then
                      nbev(i) = 4 ! right
		    elseif(nbev(i) .eq. 4) then
                     nbev(i) = 3 ! left
            endif
	       enddo
      endif


      if (Idisplay.eq.1) then
          write(*,*)'                                              '
          write(*,*)'  Number of inside cell centers      N_CELL0 = ',N_CELL0
          write(*,*)'  Number of vertices                 N_VERT  = ',N_VERT
          write(*,*)'  Number of Boundary Vertices        N_BC    = ',N_BC
          write(*,*)'                                             '
      endif


!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

      if (Idisplay.ne.1) then
         write(*,'(t12,20a)') ' --- data.txt has been successfully read ---'
      endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: input'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	RETURN
	END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF INPUT                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

