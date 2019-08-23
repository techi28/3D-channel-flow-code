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
      real*8, dimension(:,:), allocatable :: dxVV,dyVV,dlVV
      real*8, dimension(:,:), allocatable :: dxCC,dyCC,dlCC
      real*8, dimension(:,:), allocatable :: dxCVi,dyCVi
      real*8, dimension(:,:), allocatable :: dxCE,dyCE,dlCE
      real*8, dimension(:,:), allocatable :: xe,ye,xme,yme,xmo,ymo,doe
      real*8, dimension(:),   allocatable :: areaCell
      real*8, dimension(:),   allocatable :: sum_xc2,sum_yc2,sum_xcyc
      real*8, dimension(:),   allocatable :: sum_xv2,sum_yv2,sum_xvyv
      real*8, dimension(:),   allocatable :: dlVsum,areaVsum
      real*8, dimension(:,:), allocatable :: normxc,normyc
      real*8, dimension(:),   allocatable :: normxv,normyv
      real*8, dimension(:),   allocatable :: dn
      real*8, dimension(:),   allocatable :: xn1,yn1,xn2,yn2
      real*8, dimension(:,:), allocatable :: dhCE
      integer,dimension(:),   allocatable :: icxn1,icxn2
      integer,dimension(:,:), allocatable :: surrounding
      integer,dimension(:),   allocatable :: Dimsurrounding
      integer,dimension(:,:), allocatable :: Idxsurrounding
      integer,dimension(:),   allocatable :: DimIdxsurrounding
      real*8, dimension(:,:), allocatable :: weight
      integer,dimension(:),   allocatable :: ColorCell
      integer, dimension(:), allocatable :: vertex_pair
      integer, dimension(:), allocatable :: vertex_edge
      integer, dimension(:), allocatable :: tagXV
      integer, dimension(:), allocatable :: tagYV
      integer, dimension(:), allocatable :: tagXC
      integer, dimension(:), allocatable :: tagYC
      real*8, dimension(:,:), allocatable :: xnn,ynn,xe1,ye1,xe2,ye2
      real*8, dimension(:,:), allocatable :: dle1,dle2,dle3
      real*8, dimension(:,:), allocatable :: U1FACE,U2FACE,U3FACE
      real*8, dimension(:,:), allocatable :: UTFACE,UBFACE
      real*8, dimension(:,:), allocatable :: UOFLOW,VOFLOW,WOFLOW
      real*8, dimension(:), allocatable :: tago
      integer,dimension(:), allocatable :: tagsp
      integer,dimension(:),allocatable :: bcuc,bcvc,bcwc,bcpc
      integer,dimension(:),allocatable :: bcuv,bcvv,bcwv,bcpv
      real*8, dimension(:,:),allocatable :: rhsuf,rhsvf,rhswf
      integer :: bctop,bcbot
      real*8 :: TotalVolume

#     ifdef KeyPPECenter
      real*8,  dimension(:,:),allocatable   :: AmC0
      real*8,  dimension(:,:),allocatable   :: AmC1
      real*8,  dimension(:,:),allocatable   :: AmC2
      real*8,  dimension(:,:),allocatable   :: AmC3
      real*8,  dimension(:,:),allocatable   :: AmCT
      real*8,  dimension(:,:),allocatable   :: AmCB
#     endif

      CONTAINS

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                    ALLOCATE GEOMETRIC VARIABLES                     !
!                             Jul 2014                                !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE alloc_geometry

#     include "common.mpf"

      allocate(sigmat(N_CELL,NZ),                                   &
               sigmax(N_CELL,NZ),                                   &
               sigmay(N_CELL,NZ),                                   &
               sigmaz(N_CELL,NZ),                                   &
               dxVV(N_CELL0,3),dyVV(N_CELL0,3),dlVV(N_CELL0,3),     &
               dxCC(N_CELL0,3),dyCC(N_CELL0,3),dlCC(N_CELL0,3),     &
               dxCVi(N_CELL0,3),dyCVi(N_CELL0,3),                   &
               dxCE(N_CELL0,3),dyCE(N_CELL0,3),dlCE(N_CELL0,3),     &
               xe(N_CELL0,3),ye(N_CELL0,3),                         &
               xme(N_CELL0,3),yme(N_CELL0,3),                       &
               xmo(N_CELL0,3),ymo(N_CELL0,3),doe(N_CELL0,3),        &
               areaCell(N_CELL),                                    &
               sum_xc2(N_CELL0),sum_yc2(N_CELL0),sum_xcyc(N_CELL0), &
               sum_xv2(N_CELL0),sum_yv2(N_CELL0),sum_xvyv(N_CELL0), &
               dlVsum(N_VERT),areaVsum(N_VERT),                     &
               normxc(N_CELL0,3),normyc(N_CELL0,3),                 &
               normxv(N_VERT),normyv(N_VERT),                       &
               dn(N_VERT),                                          &
               xn1(N_VERT),yn1(N_VERT),xn2(N_VERT),yn2(N_VERT),     &
               dhCE(N_CELL0,3),                                     &
               icxn1(N_VERT),icxn2(N_VERT),                         &
               surrounding(N_VERT,10),Dimsurrounding(N_VERT),       &
               Idxsurrounding(N_VERT,2),DimIdxsurrounding(N_VERT),  &
               weight(N_VERT,10),                                   &
               ColorCell(N_CELL0), tagsp(N_CELL),                   &
               tagXV(N_VERT),tagYV(N_VERT),                         &
               tagXC(N_CELL),tagYC(N_CELL),                         &
               xnn(N_CELL0,3),ynn(N_CELL0,3),                       &
               xe1(N_CELL0,3),ye1(N_CELL0,3),                       &
               xe2(N_CELL0,3),ye2(N_CELL0,3),                       &
               dle1(N_CELL0,3),dle2(N_CELL0,3),dle3(N_CELL0,3),     &
               U1FACE(N_CELL0,NZ),U2FACE(N_CELL0,NZ),U3FACE(N_CELL0,NZ),&
               UTFACE(N_CELL0,NZ),UBFACE(N_CELL0,NZ))

#    ifdef KeyPPECenter
      allocate(AmC0(N_CELL0,NZ),AmC1(N_CELL0,NZ),AmC2(N_CELL0,NZ), &
               AmC3(N_CELL0,NZ),AmCT(N_CELL0,NZ),AmCB(N_CELL0,NZ))
#    endif

      END SUBROUTINE alloc_geometry


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                    DEALLOCATE GEOMETRIC VARIABLES                   !
!                             Jul 2014                                !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE dealloc_geometry

      deallocate(sigmat,sigmax,sigmay,sigmaz, &
               dxVV,dyVV,dlVV,                &
               dxCC,dyCC,dlCC,                &
               dxCVi,dyCVi,                   &
               dxCE,dyCE,dlCE,doe,            &
               xe,ye,xme,yme,xmo,ymo,areaCell,&
               sum_xc2,sum_yc2,sum_xcyc,      &
               sum_xv2,sum_yv2,sum_xvyv,      &
               dlVsum,areaVsum,               &
               normxc,normyc,normxv,normyv,dn,&
               xn1,yn1,xn2,yn2,               &
               dhCE,                          &
               icxn1,icxn2,                   &
               surrounding,Dimsurrounding,    &
               Idxsurrounding,                &
               DimIdxsurrounding, weight,     &
               ColorCell,tagsp,               &
               tagXV,tagYV,tagXC,tagYC,       &
               xnn,ynn,xe1,ye1,xe2,ye2,       &
               dle1,dle2,dle3,                &
               U1FACE,U2FACE,U3FACE,UTFACE,UBFACE)
!   ----------------------------------------------
!    deallocate the boundary bc info
     deallocate(bcuc,bcvc,bcwc,bcpc, &
                bcuv,bcvv,bcwv,bcpv)
!   ----------------------------------------------
!    deallocate variables only exist for certain option
#   ifndef KeyParallel
#   ifdef KeyTestpBC
       deallocate(vertex_pair,vertex_edge)
#   endif
#   endif
!   ----------------------------------------
!   deallocate history RHS value
       deallocate(rhsuf,rhsvf,rhswf)

#     ifdef KeyPPECenter
      deallocate(AmC0,AmC1,AmC2,AmC3,AmCT,AmCB)
#     endif

      END SUBROUTINE dealloc_geometry

      END MODULE geometry

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END MODULE GEOMETRY                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
