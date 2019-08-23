!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   CALCULATION OF THE ADVECTION TERM                 !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE advection3DNEW(rhsu,rhsv,rhsw,              &
                                phiu,phiv,phiw,              &
                                phiuv,phivv,phiwv,           &
                                xc,yc,sig,dsig,No_cp,nbe,    &
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
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|
      real*8,dimension(:,:) :: phiu(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_CELL,NZ)
      real*8,dimension(:,:) :: phiw(N_CELL,NZ)
!     ----------------------------------------
      real*8,dimension(:,:) :: phiuv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: phivv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: phiwv(N_VERT,NZ-1)
!     ----------------------------------------
      real*8,dimension(:,:) :: rhsu(N_CELL,NZ)
      real*8,dimension(:,:) :: rhsv(N_CELL,NZ)
      real*8,dimension(:,:) :: rhsw(N_CELL,NZ)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ----------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
      integer :: Option
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      real*8, dimension(:,:),allocatable :: UFi1,UFi2,UFi3,UFiT,UFiB
      real*8, dimension(:,:),allocatable :: Ui1,Ui2,Ui3
      real*8, dimension(:,:),allocatable :: Vi1,Vi2,Vi3
      real*8, dimension(:,:),allocatable :: Wi1,Wi2,Wi3
      real*8, dimension(:,:),allocatable :: UiT,UiB
      real*8, dimension(:,:),allocatable :: ViT,ViB
      real*8, dimension(:,:),allocatable :: WiT,WiB
      real*8, dimension(:) :: coef1(1:5)
      real*8, dimension(:) :: coef2(1:5)
      real*8 :: uflux,vflux,wflux,Vol
      real*8 :: ucor,vcor,wcor
      integer :: jc1,jc2,jc3,kk

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Advection3DNEW'
#     endif
!*********************************************************************!
!                                                                     !
!                      FACE NORMAL VALUE calculation                  !
!                                                                     !
!*********************************************************************!

      allocate(UFi1(N_CELL0,NZ),         &
               UFi2(N_CELL0,NZ),         &
               UFi3(N_CELL0,NZ),         &
               UFiT(N_CELL0,NZ),         &
               UFiB(N_CELL0,NZ))

      do k=1,NZ
         do i=1,N_CELL0
            UFi1(i,k) = 0.0d0
            UFi2(i,k) = 0.0d0
            UFi3(i,k) = 0.0d0
            UFiT(i,k) = 0.0d0
            UFiB(i,k) = 0.0d0
         enddo
      enddo

      DO k=2,NZ-1
         do i=1,N_CELL0
            UFi1(i,k) = U1FACE(i,k)*dlVV(i,1)*dsigv(k-1)
            UFi2(i,k) = U2FACE(i,k)*dlVV(i,2)*dsigv(k-1)
            UFi3(i,k) = U3FACE(i,k)*dlVV(i,3)*dsigv(k-1)
            UFiT(i,k) = UTFACE(i,k)*areaCell(i)
            UFiB(i,k) = UBFACE(i,k)*areaCell(i)
        enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                    centre velocity interpolation                    !
!                                                                     !
!*********************************************************************!
      allocate(Ui1(N_CELL0,NZ),         &
               Ui2(N_CELL0,NZ),         &
               Ui3(N_CELL0,NZ),         &
               Vi1(N_CELL0,NZ),         &
               Vi2(N_CELL0,NZ),         &
               Vi3(N_CELL0,NZ),         &
               Wi1(N_CELL0,NZ),         &
               Wi2(N_CELL0,NZ),         &
               Wi3(N_CELL0,NZ),         &
               UiT(N_CELL0,NZ),         &
               UiB(N_CELL0,NZ),         &
               ViT(N_CELL0,NZ),         &
               ViB(N_CELL0,NZ),         &
               WiT(N_CELL0,NZ),         &
               WiB(N_CELL0,NZ))

      do k=1,NZ
         do i=1,N_CELL0
            Ui1(i,k) = 0.0d0
            Ui2(i,k) = 0.0d0
            Ui3(i,k) = 0.0d0
            Vi1(i,k) = 0.0d0
            Vi2(i,k) = 0.0d0
            Vi3(i,k) = 0.0d0
            Wi1(i,k) = 0.0d0
            Wi2(i,k) = 0.0d0
            Wi3(i,k) = 0.0d0
            UiT(i,k) = 0.0d0
            UiB(i,k) = 0.0d0
            ViT(i,k) = 0.0d0
            ViB(i,k) = 0.0d0
            WiT(i,k) = 0.0d0
            WiB(i,k) = 0.0d0
         enddo
      enddo

      call calcul_face(Ui1,Ui2,Ui3,UiT,UiB,phiu,phiuv,    &
                       xc,yc,sig,dsig,No_cp,nbe,          &
                       xv,yv,sigv,dsigv,No_vp,nbev,2)
!    -----------------------------------------------------------
      call calcul_face(Vi1,Vi2,Vi3,ViT,ViB,phiv,phivv,    &
                       xc,yc,sig,dsig,No_cp,nbe,          &
                       xv,yv,sigv,dsigv,No_vp,nbev,2)
!    -----------------------------------------------------------
      call calcul_face(Wi1,Wi2,Wi3,WiT,WiB,phiw,phiwv,    &
                       xc,yc,sig,dsig,No_cp,nbe,          &
                       xv,yv,sigv,dsigv,No_vp,nbev,2)
!*********************************************************************!
!                                                                     !
!                   advection contribution                            !
!                                                                     !
!*********************************************************************!
      do k=2, NZ-1
        do i=1,N_CELL0
!         ----------------------------
!         For u
          uflux = Ui1(i,k)*UFi1(i,k) &
                + Ui2(i,k)*UFi2(i,k) &
                + Ui3(i,k)*UFi3(i,k) &
                + UiT(i,k)*UFiT(i,k) &
                + UiB(i,k)*UFiB(i,k)
!         ----------------------------
!         For v
          vflux = Vi1(i,k)*UFi1(i,k) &
                + Vi2(i,k)*UFi2(i,k) &
                + Vi3(i,k)*UFi3(i,k) &
                + ViT(i,k)*UFiT(i,k) &
                + ViB(i,k)*UFiB(i,k)
!         ----------------------------
!         For w
          wflux = Wi1(i,k)*UFi1(i,k) &
                + Wi2(i,k)*UFi2(i,k) &
                + Wi3(i,k)*UFi3(i,k) &
                + WiT(i,k)*UFiT(i,k) &
                + WiB(i,k)*UFiB(i,k)

!        -----------------------------
!         correction for implicit method
#       ifdef KeyImplicit
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)

#       ifdef KeyAdvUpwind
            if(UFi1(i,k) .gt. 0.0d0) then
                coef1(1) = 1.0d0
                coef2(1) = 0.0d0
            else
                coef1(1) = 0.0d0
                coef2(1) = 1.0d0
            endif

            if(UFi2(i,k) .gt. 0.0d0) then
                coef1(2) = 1.0d0
                coef2(2) = 0.0d0
            else
                coef1(2) = 0.0d0
                coef2(2) = 1.0d0
            endif

            if(UFi3(i,k) .gt. 0.0d0) then
                coef1(3) = 1.0d0
                coef2(3) = 0.0d0
            else
                coef1(3) = 0.0d0
                coef2(3) = 1.0d0
            endif

            if(UFiT(i,k) .gt. 0.0d0) then
                coef1(4) = 1.0d0
                coef2(4) = 0.0d0
            else
                coef1(4) = 0.0d0
                coef2(4) = 1.0d0
            endif

            if(UFiB(i,k) .gt. 0.0d0) then
                coef1(5) = 1.0d0
                coef2(5) = 0.0d0
            else
                coef1(5) = 0.0d0
                coef2(5) = 1.0d0
            endif
#       else
            do kk=1,5
                coef1(kk) = 0.5d0
                coef2(kk) = 0.5d0
            enddo
#       endif

            ucor =     (coef1(1)*phiu(i,k)+coef2(1)*phiu(jc1,k))*UFi1(i,k)     &
                    +  (coef1(2)*phiu(i,k)+coef2(2)*phiu(jc2,k))*UFi2(i,k)     &
                    +  (coef1(3)*phiu(i,k)+coef2(3)*phiu(jc3,k))*UFi3(i,k)     &
                    +  (coef1(4)*phiu(i,k)+coef2(4)*phiu(i,k+1))*UFiT(i,k)     &
                    +  (coef1(5)*phiu(i,k)+coef2(5)*phiu(i,k-1))*UFiB(i,k)

            vcor =     (coef1(1)*phiv(i,k)+coef2(1)*phiv(jc1,k))*UFi1(i,k)     &
                    +  (coef1(2)*phiv(i,k)+coef2(2)*phiv(jc2,k))*UFi2(i,k)     &
                    +  (coef1(3)*phiv(i,k)+coef2(3)*phiv(jc3,k))*UFi3(i,k)     &
                    +  (coef1(4)*phiv(i,k)+coef2(4)*phiv(i,k+1))*UFiT(i,k)     &
                    +  (coef1(5)*phiv(i,k)+coef2(5)*phiv(i,k-1))*UFiB(i,k)

            wcor =     (coef1(1)*phiw(i,k)+coef2(1)*phiw(jc1,k))*UFi1(i,k)     &
                    +  (coef1(2)*phiw(i,k)+coef2(2)*phiw(jc2,k))*UFi2(i,k)     &
                    +  (coef1(3)*phiw(i,k)+coef2(3)*phiw(jc3,k))*UFi3(i,k)     &
                    +  (coef1(4)*phiw(i,k)+coef2(4)*phiw(i,k+1))*UFiT(i,k)     &
                    +  (coef1(5)*phiw(i,k)+coef2(5)*phiw(i,k-1))*UFiB(i,k)

            uflux = (uflux - ucor)  + 0.5d0*ucor
            vflux = (vflux - vcor)  + 0.5d0*vcor
            wflux = (wflux - wcor)  + 0.5d0*wcor
#        endif
!        -----------------------------
!        rhs
            Vol = areaCell(i)*dsigv(k-1)
            rhsu(i,k) = rhsu(i,k) - uflux/Vol
            rhsv(i,k) = rhsv(i,k) - vflux/Vol
            rhsw(i,k) = rhsw(i,k) - wflux/Vol
         enddo
      enddo
!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(UFi1,UFi2,UFi3,UFiT,UFiB)

      deallocate(Ui1,Ui2,Ui3, &
                 Vi1,Vi2,Vi3, &
                 Wi1,Wi2,Wi3, &
                 WiT,WiB,ViT,ViB,UiT,UiB)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Advection3DNEW'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
