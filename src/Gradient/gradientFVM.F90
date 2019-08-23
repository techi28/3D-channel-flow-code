!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    HORITZONTAL GRADIENT USING FVM                   !
!                              Oct 2016                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE grandientFVM(dphidx,dphidy,dphids,phi,phiv,   &
                              xc,yc,sig,dsig,No_cp,nbe,        &
                              xv,yv,sigv,dsigv,No_vp,nbev,Option)
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
!    |      Declaration of variables      |
!    |____________________________________|
!     --------------------------------------
      real*8, dimension(:,:) :: dphidx(N_CELL,NZ)
      real*8, dimension(:,:) :: dphidy(N_CELL,NZ)
      real*8, dimension(:,:) :: dphids(N_CELL,NZ)
      real*8, dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:,:) :: phiv(N_VERT,NZ-1)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      real*8, dimension(:)   :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0)
!     --------------------------------------
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)
      integer :: Option
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      real*8, dimension(:,:),allocatable :: funf1,funf2,funf3,funfT,funfB
      real*8, dimension(:) :: phif(1:3)
      real*8 :: phifT,phifB,Vol,residx,residy
      integer :: flg
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: grandientFVM'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      Option = 1  ! calculate direct gradient using face value (cell center)
!      Option = 2  ! Two step gradient calculation (smooth gradient for face)
       flg = Option
!      NOTE: THE OPTION WAS INTENTIALLY SET AS 1, AS THERE SEEMS TO BE
!      STABILITY PROBLEMS IF THE VALUE IS 2. (EVEN WITH IMPLICIT METHOD)
!      WILL INVESTIGATE LATER.
!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
     do k=1,NZ
        do i=1,N_CELL
            dphidx(i,k) = 0.
            dphidy(i,k) = 0.
            dphids(i,k) = 0.
        enddo
     enddo

      allocate(funf1(N_CELL0,NZ),         &
               funf2(N_CELL0,NZ),         &
               funf3(N_CELL0,NZ),         &
               funfT(N_CELL0,NZ),         &
               funfB(N_CELL0,NZ))

     call calcul_face(funf1,funf2,funf3,funfT,funfB, &
                      phi,phiv,                      &
                      xc,yc,sig,dsig,No_cp,nbe,      &
                      xv,yv,sigv,dsigv,No_vp,nbev,flg)
!    -------------------------
       DO k=2,NZ-1
          do i=1,N_CELL0
             Vol = areaCell(i)
             residx = 0.0d0
             residy = 0.0d0
             phif(1) = funf1(i,k)
             phif(2) = funf2(i,k)
             phif(3) = funf3(i,k)
!           ___________________________________________________
!           Horizontal neighbors
              do j=1,3
                    residx = residx + phif(j)*dlVV(i,j)*xnn(i,j)
                    residy = residy + phif(j)*dlVV(i,j)*ynn(i,j)
	          enddo
!           -------------------------------------------------
!           Final source term
              dphidx(i,k) = residx/Vol
              dphidy(i,k) = residy/Vol
              dphids(i,k) = (funfT(i,k)-funfB(i,k))/dsigv(k-1)
        enddo
       ENDDO

       call BCpcenter(dphidx,xc,yc,sig,dsig,No_cp,nbe)
       call BCpcenter(dphidy,xc,yc,sig,dsig,No_cp,nbe)
       call BCpcenter(dphids,xc,yc,sig,dsig,No_cp,nbe)

       deallocate(funf1,funf2,funf3,funfT,funfB)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: grandientFVM'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
