!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                MODULE WITH COMMON GEOMETRIC VRAIBLES                !
!                      Miguel Angel Uh Zapata                         !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!   
!                                                                     !
!      SUBROUTINES:     -  alloc_geometry                             !
!                       -  dealloc_geometry                           !
!                                                                     !
!---------------------------------------------------------------------!
!      ________________________________________________________       !
!      sigma-transform variables                                      !
!                                                                     !
!      sigmat    :  transformation sigma in the time                  !
!      sigmax    :  transformation sigma in the x direction           !
!      sigmay    :  transformation sigma in the y direction           !
!      sigmaz    :  transformation sigma in the z direction           !         
!      ________________________________________________________       !
!      Horizontal geometrical vector variables                        !
!                                                                     !
!      dxVV      :  dx length of each side of the triangle            !
!      dyVV      :  dy length of each side of the triangle            !
!      dlVV      :  Distance from vertex to vertex (anticlockwise)    !
!      dxVV2     :  = dxVV*dxVV                                       !
!      dyVV2     :  = dyVV*dyVV                                       !
!                                                                     !
!      dxCC      :  Length dx from cell center to cell neigbor        !
!      dyCC      :  Length dy from cell center to cell neigbor        !
!      dlCC      :  Distance from the cell center to the cell neigbor !
!                                                                     !
!      xe,ye     :  (x,y) coordinates of the edge point               !
!      dxCE      :  Length dx from cell center to the edge point      !
!      dyCE      :  Length dy from cell center to the edge point      !
!      dlCE      :  Distance from the cell center to the edge point   !
!                                                                     !
!      dlCV      :  Distance from the cell center to each vertex      !
!      dlVsum    :  Total sum of dlCV distances at each vertex        !
!      areaCell  :  (N_CELL) Area of each cell                        !
!      areaVsum  :  Total area related to each vertex.                !
!      sum_xc2   :  (N_CELL0) Sum(xc(i)-xc(neighbor))^2               ! 
!      sum_yc2   :  (N_CELL0) Sum(yc(i)-yc(neighbor))^2               !
!      sum_xcyc  :  (N_CELL0) Sum(xc(i)-xc(neig))*(yc(i)-yc(neig))    !
!                                                                     !
!      VolPrism  :  (N_CELL,NZ) Volume of each prism element          !
!      dsigPrism :  (N_CELL)    Height of each prism element          !
!      ________________________________________________________       !
!                                                                     !
!      normxc    : (N\_CELL0,3) normal at cell-center point           !
!      normyc    : (N\_CELL0,3) normal at cell-center point           !
!      normxv    : (N\_VERT)    normal at vertex point                !
!      normyv    : (N\_VERT)    normal at vertex point                !
!      dn        : (N\_VERT)  distance of normal points               !
!      xn1,yn1   : (N\_VERT)  point in the normal direction:  dn      !
!      xn2,yn2   : (N\_VERT)  point in the normal direction:  2*dn    !
!         icxn1  : (N\_VERT)  Index of normal points boundary related !
!         icxn2  : (N\_VERT)  Index of normal points boundary related !
!      ________________________________________________________       !
!                                                                     !
!    surrounding    : (N\_VERT,10)  Indexes of cells surrounding      !
!    Dimsurrounding : (N\_VERT)     Number of cells surrounding       !
!            weight : (N\_VERT,10)) Weight of each surrounding cells  !
!      ________________________________________________________       !
!      Diffusion geometry                                             !
!                                                                     !
!      dhCE     :  (N_CELL0,3) Perpendicular distance from the        !
!                              cell center to the edge point          !
!      xme,yme  :  (N_CELL0,3) Edge points of a triangle              !
!                                                                     !
!---------------------------------------------------------------------!

      MODULE geometry

#     include "cppdefs.h"
      implicit none
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
#        include "common.mpf"
#     endif
!     =============== END ================    
!     ====================================

!     ____________________________________
!    |                                    |
!    |               Variables            |
!    |____________________________________|

      real*8, dimension(:,:), allocatable :: sigmat,sigmax,sigmay,sigmaz
      real*8, dimension(:,:), allocatable :: VolPrism       
      real*8, dimension(:,:), allocatable :: dxVV,dyVV,dlVV
      real*8, dimension(:,:), allocatable :: dxCC,dyCC,dlCC
      real*8, dimension(:,:), allocatable :: dxCE,dyCE,dlCE,dlCV 
      real*8, dimension(:,:), allocatable :: xe,ye,xme,yme,xmo,ymo,doe
      real*8, dimension(:),   allocatable :: areaCell
      real*8, dimension(:),   allocatable :: sum_xc2,sum_yc2,sum_xcyc  
      real*8, dimension(:),   allocatable :: dlVsum,areaVsum
      real*8, dimension(:,:), allocatable :: normxc,normyc
      real*8, dimension(:),   allocatable :: normxv,normyv
      real*8, dimension(:),   allocatable :: dn
      real*8, dimension(:),   allocatable :: xn1,yn1,xn2,yn2
      real*8, dimension(:,:), allocatable :: dhCE
      integer,dimension(:),   allocatable :: icxn1,icxn2
      integer,dimension(:,:), allocatable :: surrounding
      integer,dimension(:),   allocatable :: Dimsurrounding
      real*8, dimension(:,:), allocatable :: weight      
      integer,dimension(:,:), allocatable :: surroundingLSM
      integer,dimension(:),   allocatable :: DimsurroundingLSM           
      real*8, dimension(:,:), allocatable :: weightLSM
      real*8, dimension(:),   allocatable :: dlVsumLSM            
      integer,dimension(:),   allocatable :: ColorCell
      real*8, dimension(:,:), allocatable :: aGx,aGy,LSM          
      real*8, dimension(:,:), allocatable :: U1FACE,U2FACE,U3FACE
      real*8, dimension(:,:), allocatable :: UTFACE,UBFACE
      real*8, dimension(:),   allocatable :: wfT,wfvT
      real*8, dimension(:),   allocatable :: wfB,wfvB
      real*8, dimension(:),   allocatable :: dw1dt
!     ------------------------
      real*8, dimension(:,:), allocatable :: pf_Last
      real*8, dimension(:),   allocatable :: Shu_Last,Shv_Last
      real*8, dimension(:,:), allocatable :: walld,eddy
!     ------------------------
      real*8, dimension(:,:), allocatable :: Vortx,Vortxv
      real*8, dimension(:,:), allocatable :: Vorty,Vortyv
      real*8, dimension(:,:), allocatable :: Vortz,Vortzv
      real*8, dimension(:,:), allocatable :: Vort,Vortv
!     ------------------------
      integer,dimension(:,:), allocatable :: TagPeriodicBC
                  
      CONTAINS

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                    ALLOCATE GEOMETRIC VARIABLES                     !
!                             Jul 2014                                !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE alloc_geometry

#     include "common.mpf"

      allocate(sigmat(N_CELL,NZglobal),                             &
               sigmax(N_CELL,NZglobal),                             &
               sigmay(N_CELL,NZglobal),                             &  
               sigmaz(N_CELL,NZglobal),                             &
               VolPrism(N_CELL,NZglobal),                           &
               dxVV(N_CELL0,3),dyVV(N_CELL0,3),dlVV(N_CELL0,3),     &
               dxCC(N_CELL0,3),dyCC(N_CELL0,3),dlCC(N_CELL0,3),     &
               dxCE(N_CELL0,3),dyCE(N_CELL0,3),dlCE(N_CELL0,3),     &
               dlCV(N_CELL,3),                                      &
               xe(N_CELL0,3),ye(N_CELL0,3),                         &
               xme(N_CELL0,3),yme(N_CELL0,3),                       &
               areaCell(N_CELL),                                    &
               xmo(N_CELL0,3),ymo(N_CELL0,3),doe(N_CELL0,3),        &
               sum_xc2(N_CELL0),sum_yc2(N_CELL0),sum_xcyc(N_CELL0), &
               dlVsum(N_VERT),areaVsum(N_VERT),                     &
               normxc(N_CELL0,3),normyc(N_CELL0,3),                 &
               normxv(N_VERT),normyv(N_VERT),                       &
               dn(N_VERT),                                          &
               xn1(N_VERT),yn1(N_VERT),xn2(N_VERT),yn2(N_VERT),     &
               dhCE(N_CELL0,3),                                     &
               icxn1(N_VERT),icxn2(N_VERT))
               
      allocate(surrounding(N_VERT,10),                              &
               Dimsurrounding(N_VERT),                              &
               weight(N_VERT,10),                                   &
               surroundingLSM(N_VERT,40),                           &
               DimsurroundingLSM(N_VERT),                           &
               weightLSM(N_VERT,40),                                &
               dlVsumLSM(N_VERT),                                   &                                              
               ColorCell(N_CELL0),                                  &
               aGx(N_CELL0,0:3),aGy(N_CELL0,0:3),LSM(N_CELL0,0:3))
                              
      allocate(U1FACE(N_CELL0,NZglobal),                            &
               U2FACE(N_CELL0,NZglobal),                            &
               U3FACE(N_CELL0,NZglobal),                            &
               UTFACE(N_CELL0,NZglobal),                            &
               UBFACE(N_CELL0,NZglobal))

      allocate(eddy(N_CELL,NZglobal))
      allocate(walld(N_CELL,NZglobal))
      allocate(wfT(N_CELL),wfvT(N_VERT))
      allocate(wfB(N_CELL),wfvB(N_VERT))
      allocate(dw1dt(N_CELL0))
      allocate(pf_Last(N_CELL,NZglobal))
      allocate(Shu_Last(N_CELL))
      allocate(Shv_Last(N_CELL))

      allocate(Vortx(N_CELL,NZglobal))  
      allocate(Vorty(N_CELL,NZglobal))  
      allocate(Vortz(N_CELL,NZglobal))  
      allocate(Vort(N_CELL,NZglobal))  
      allocate(Vortxv(N_VERT,NZglobal-1))      
      allocate(Vortyv(N_VERT,NZglobal-1))
      allocate(Vortzv(N_VERT,NZglobal-1))
      allocate(Vortv(N_VERT,NZglobal-1))
      
      allocate(TagPeriodicBC(N_VERT,2))

      END SUBROUTINE alloc_geometry


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                    DEALLOCATE GEOMETRIC VARIABLES                   !
!                             Jul 2014                                !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE dealloc_geometry

      deallocate(sigmat,sigmax,sigmay,sigmaz, &
               VolPrism,                      &
               dxVV,dyVV,dlVV,                &
               dxCC,dyCC,dlCC,                &
               dxCE,dyCE,dlCE,dlCV,           &
               xe,ye,xme,yme,areaCell,        &
               xmo,ymo,doe,                   &
               sum_xc2,sum_yc2,sum_xcyc,      &
               dlVsum,areaVsum,               &
               normxc,normyc,normxv,normyv,dn,&
               xn1,yn1,xn2,yn2,               &
               dhCE,                          &
               icxn1,icxn2)
               
      deallocate(surrounding,Dimsurrounding,  &
               weight,                        &
               surroundingLSM,                &
               DimsurroundingLSM,             &
               weightLSM,                     &
               dlVsumLSM,                     &                              
               ColorCell,                     &
               aGx,aGy,LSM)
               
      deallocate(U1FACE,U2FACE,U3FACE,        &
               UTFACE,UBFACE)

      deallocate(walld,eddy)
      deallocate(wfT,wfvT)
      deallocate(wfB,wfvB)
      deallocate(dw1dt)
      deallocate(pf_Last)
      deallocate(Shu_Last)
      deallocate(Shv_Last)

      deallocate(Vortx,Vortxv)
      deallocate(Vorty,Vortyv)
      deallocate(Vortz,Vortzv)      
      deallocate(Vort,Vortv)

      deallocate(TagPeriodicBC)

      END SUBROUTINE dealloc_geometry

      END MODULE geometry

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END MODULE GEOMETRY                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
