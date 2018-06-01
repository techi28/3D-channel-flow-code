!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE interpolationEta(etav,eta,&
                                  No_vp)       
!---------------------------------------------------------------------!
!                                                                     !
!    This program interpolate the values phi from the center of the   !
!    cell to the vertex of the triangles including all the possible   !
!    vertices.                                                        !
!                                                                     !
!    We have two option here: we can choose the area of the triangle  !
!    or the distance the cell center-vertex as weighting value. This  !
!    option are declared in "cppdefs.h" as:                           !
!             -     #define KeyInterpoDist                            !
!             -     #define KeyInterpoArea                            !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <-- phiv    | (N_VERT)  | Function phi at the vertex          |  !  
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> phi     |(N_CELL)   | Function phi at the cell center     |  !
!  | --> No_vp   |(N_CELL0,3)| Numbering of cell vertices          |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  ! 
!  |--- N_CELL0   | Number of the cell centers inside the domain   |  !
!  |--- N_VERT    | Number of the computing vertices               |  !
!  |______________|________________________________________________|  !  
!  |   dlCV       |(CELL,3)  | Distance from the center to vertex  |  !
!  |   dlVsum     |(N_VERT)  | Total sum of dlCV distances at vert.|  !
!  |   areaCELL   |(N_CELL,3)| Area of each cell                   |  !
!  |   areaVsum   |(N_VERT)  | Total area related to each vertex.  |  !  
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

      real*8 ,dimension(:)  :: etav
      real*8 ,dimension(:)  :: eta
      integer,dimension(:,:):: No_vp 

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: interpolationETA'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

      do nv=1,N_VERT
            etav(nv) = 0.0d0
      enddo

!      ________________________________________________________
!     |                                                        |
!     |  USING CELL-VERTEX DISTANCE WEIGHTING TO INTERPOLATE   |
!     |________________________________________________________|
      
#     ifdef KeyInterpoDist
            do i=1,N_CELL0	
	       do j=1,3
	          nv = No_vp(i,j)
                  etav(nv)= dlCV(i,j)*eta(i) + etav(nv)
               enddo
            enddo

            do nv=1,N_VERT
               etav(nv) = etav(nv)/dlVsum(nv)
            enddo
#     endif

!      ________________________________________________________
!     |                                                        |
!     |         USING AREA WEIGHTING TO INTERPOLATE            |
!     |________________________________________________________|
    
#     ifdef KeyInterpoArea
            do i=1,N_CELL0	
	       do j=1,3
	          nv = No_vp(i,j)
                  etav(nv)= areaCELL(i)*eta(i) + etav(nv)
               enddo
            enddo

            do nv=1,N_VERT
               etav(nv) = etav(nv)/areaVsum(nv)
            enddo
#     endif

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: interpolationETA'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
