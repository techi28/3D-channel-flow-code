!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                 phi,xc,yc,sig,dsig,No_cp,nbe)
       
!---------------------------------------------------------------------!
!                                                                     !
!    This program interpolate the values phi from the center of the   !
!    elements to the vertex of the prism. It is important to mention  !
!    that we DO NOT update the vertex values at the boundary, then    !
!    the program is used to interpolate main variables as u,v,w,p     !
!    where the boundary conditions are already defined.               !  
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
!  | <-- phiv    |(N_VERT,NZ)| Function phi at the vertices        |  !  
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size    | Description                        |  !  
!  |_____________|____________|____________________________________|  !  
!  | --> phi     |(N_CELL,NZ) |Function phi at the cell center     |  !
!  | --> sigv    |(NZ-1)      |sigma coordinate vetex points       |  !
!  | --> dsigv   |(NZ-1)      |Increment of the sigma vertex points|  !  
!  | --> No_vp   |(N_CELL0,3) |Numbering of cell vertices          |  !
!  | --> nbev    |(N_VERT)    |Type of tag about the kind of vertex|  !
!  |_____________|____________|____________________________________|  !
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
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,k       | Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |            -   interpolation2D                                |  !
!  |_______________________________________________________________|  !
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

      real*8, dimension(:,:):: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  

      real*8, dimension(:,:):: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 

!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________|

      real*8,dimension(:),  allocatable :: phiB
      real*8,dimension(:),  allocatable :: phivB
      real*8,dimension(:,:),allocatable :: phivdz0
      real*8 :: dzT,dzB

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: interpolation3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(phiB(N_CELL),phivB(N_VERT),phivdz0(N_VERT,NZ))

!*********************************************************************!
!                                                                     !
!                           Interpolation                             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                Interpolation inside domain             |  
!     |________________________________________________________|

!     -------------------------------------------------
!     Interpolation vertices at position sig(k)

      DO k=1,NZ
         !do i=1,N_CELL
         !    phiB(i)= phi(i,k)
         !enddo
         call interpolation2D(phivB,xv,yv,No_vp,nbev, &
                              phi(1:N_CELL,k),xc,yc,No_cp,nbe)
         do nv=1,N_VERT
            phivdz0(nv,k) = phivB(nv)
         enddo
      ENDDO             

!     -------------------------------------------------
!     Interpolation vertices at position sigv(k=2:NZ)
      DO k=1,NZ-1        
         dzT = abs(sigv(k)-sig(k+1))
         dzB = abs(sigv(k)-sig(k))
         do nv=1,N_VERT
            phiv(nv,k)=(dzB*phivdz0(nv,k+1)+dzT*phivdz0(nv,k))/(dzT+dzB)
         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(phiB,phivB,phivdz0)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: interpolation3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
