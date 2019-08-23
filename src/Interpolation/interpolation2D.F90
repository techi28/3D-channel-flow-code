!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE interpolation2D(funv,xv,yv,No_vp,nbev,&
                                 fun,xc,yc,No_cp,nbe)
!---------------------------------------------------------------------!
!                                                                     !
!    This program interpolate the values phi from the center of the   !
!    cell to the vertex of the triangles. It is a 2D interpolation    !
!    in the horizontal direction. It is important to remark that we   !
!    need the exact values in the boundaries of our domain.           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | <-- funv    |(N_VERT)   | Function phi at the vertex          |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | --> fun     |(N_CELL)   | Function phi at the cell center     |  !
!  | --> No_vp   |(N_CELL0,3)| Numbering of cell vertices          |  !
!  | --> nbev    |(N_VERT)   | Type of tag about the kind of vertex|  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !
!  |______________|________________________________________________|  !
!  |--- N_CELL0   | Number of the cell centers inside the domain   |  !
!  |--- N_CELL    | Total number of cell centers                   |  !
!  |--- N_VERT    | Number of the computing vertices               |  !
!  |--- NZ        | Number of points in the vertical direction     |  !
!  |______________|________________________________________________|  !
!  |    dlCV      |(CELL,3)  | Distance from the center to vertex  |  !
!  |    dlVsum    | N_VERT   | Total sum of dlCV distances at vert.|  !
!  |    areaCELL  |(N_CELL,3)| Area of each cell                   |  !
!  |    areaVsum  | N_VERT   | Total area related to each vertex.  |  !
!  |______________|__________|_____________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !
!  |______________|________________________________________________|  !
!  | * nv,i,k     | Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
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

!      ________________________________________________________
!     |                                                        |
!     |   Keys and common parameters                           |
!     |________________________________________________________|

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
!      ________________________________________________________
!     |                                                        |
!     |    Declaration of variables                            |
!     |________________________________________________________|

      real*8 ,dimension(:)  :: funv(N_VERT)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)

      real*8 ,dimension(:)  :: fun(N_CELL)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: interpolation2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

      do nv=1,N_VERT
         funv(nv) = 0.0d0
      enddo

!      ________________________________________________________
!     |                                                        |
!     |  USING NEW TECHNIQUE OF DISTANCE WEIGHTING COLLECTING  |
!     |________________________________________________________|

         do nv=1,N_VERT
	        do j=1,Dimsurrounding(nv)
               nc = surrounding(nv,j)
               funv(nv)= funv(nv) + weight(nv,j)*fun(nc)
            enddo
            funv(nv) = funv(nv)/dlVsum(nv)
         enddo

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: interpolation2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
