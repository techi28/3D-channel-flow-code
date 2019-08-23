!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                  INTERPOLATION FACE NORMAL VELOCITY                 !
!                              MAR 2016                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      SUBROUTINE interpolateU(un1,un2,un3,unT,unB,        &
                              phiu,phiv,phiw,             &
                              xc,yc,sig,dsig,No_cp,nbe,   &
                              xv,yv,sigv,dsigv,No_vp,nbev)

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
      real*8, dimension(:,:) :: un1(N_CELL0,NZ)
      real*8, dimension(:,:) :: un2(N_CELL0,NZ)
      real*8, dimension(:,:) :: un3(N_CELL0,NZ)
      real*8, dimension(:,:) :: unT(N_CELL0,NZ)
      real*8, dimension(:,:) :: unB(N_CELL0,NZ)
      real*8, dimension(:,:) :: phiu(N_CELL,NZ)
      real*8, dimension(:,:) :: phiv(N_CELL,NZ)
      real*8, dimension(:,:) :: phiw(N_CELL,NZ)
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      real*8, dimension(:,:),allocatable :: Ui1,Ui2,Ui3,UiT,UiB
      real*8, dimension(:,:),allocatable :: Vi1,Vi2,Vi3,ViT,ViB
      real*8, dimension(:,:),allocatable :: Wi1,Wi2,Wi3,WiT,WiB
      real*8, dimension(:,:),allocatable :: phiuv,phivv,phiwv
      integer :: jc1,jc2,jc3
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: interpolateU'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|
      allocate(Ui1(N_CELL0,NZ),Ui2(N_CELL0,NZ),Ui3(N_CELL0,NZ),&
               Vi1(N_CELL0,NZ),Vi2(N_CELL0,NZ),Vi3(N_CELL0,NZ),&
               Wi1(N_CELL0,NZ),Wi2(N_CELL0,NZ),Wi3(N_CELL0,NZ),&
               UiT(N_CELL0,NZ),UiB(N_CELL0,NZ),                &
               ViT(N_CELL0,NZ),ViB(N_CELL0,NZ),                &
               WiT(N_CELL0,NZ),WiB(N_CELL0,NZ))
!    --------------------------------------------------
      allocate(phiuv(N_VERT,NZ-1),         &
               phivv(N_VERT,NZ-1),         &
               phiwv(N_VERT,NZ-1))
!    --------------------------------------------------
      call BCglobalVC(phiu,phiv,phiw,                &
                      phiuv,phivv,phiwv,             &
                      xc,yc,sig,dsig,No_cp,nbe,      &
                      xv,yv,sigv,dsigv,No_vp,nbev,2)
!*********************************************************************!
!                                                                     !
!                          face value calculation                     !
!                                                                     !
!*********************************************************************!

      do k=1,NZ
         do i=1,N_CELL0
            Ui1(i,k) = 0.0d0
            Ui2(i,k) = 0.0d0
            Ui3(i,k) = 0.0d0
            Vi1(i,k) = 0.0d0
            Vi2(i,k) = 0.0d0
            Vi3(i,k) = 0.0d0
            WiT(i,k) = 0.0d0
            WiB(i,k) = 0.0d0
         enddo
      enddo

      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            Ui1(i,k) = 0.5d0*(phiu(i,k)+phiu(jc1,k))
            Ui2(i,k) = 0.5d0*(phiu(i,k)+phiu(jc2,k))
            Ui3(i,k) = 0.5d0*(phiu(i,k)+phiu(jc3,k))
            Vi1(i,k) = 0.5d0*(phiv(i,k)+phiv(jc1,k))
            Vi2(i,k) = 0.5d0*(phiv(i,k)+phiv(jc2,k))
            Vi3(i,k) = 0.5d0*(phiv(i,k)+phiv(jc3,k))
            WiT(i,k) = 0.5d0*(phiw(i,k)+phiw(i,k+1))
            WiB(i,k) = 0.5d0*(phiw(i,k)+phiw(i,k-1))
         enddo
      enddo

!       call calcul_face(Ui1,Ui2,Ui3,UiT,UiB,phiu,phiuv, &
!                        xc,yc,sig,dsig,No_cp,nbe,       &
!                        xv,yv,sigv,dsigv,No_vp,nbev,2)
!       call calcul_face(Vi1,Vi2,Vi3,ViT,ViB,phiv,phivv, &
!                        xc,yc,sig,dsig,No_cp,nbe,       &
!                        xv,yv,sigv,dsigv,No_vp,nbev,2)
!       call calcul_face(Wi1,Wi2,Wi3,WiT,WiB,phiw,phiwv, &
!                        xc,yc,sig,dsig,No_cp,nbe,       &
!                        xv,yv,sigv,dsigv,No_vp,nbev,2)

      DO k=2,NZ-1
         do i=1,N_CELL0
!           ___________________________________________________
!           Horizontal neighbors
            un1(i,k) = Ui1(i,k)*xnn(i,1)+Vi1(i,k)*ynn(i,1)
            un2(i,k) = Ui2(i,k)*xnn(i,2)+Vi2(i,k)*ynn(i,2)
            un3(i,k) = Ui3(i,k)*xnn(i,3)+Vi3(i,k)*ynn(i,3)
!          ___________________________________________________
!           Top and Bottom
            unT(i,k) = WiT(i,k)
            unB(i,k) = -1.0d0*WiB(i,k)
        enddo
      ENDDO


!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

         deallocate(Ui1,Ui2,Ui3, &
                    Vi1,Vi2,Vi3, &
                    Wi1,Wi2,Wi3, &
                    UiT,UiB,     &
                    ViT,ViB,     &
                    WiT,WiB)

        deallocate(phiuv,phivv,phiwv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: interpolateU'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
