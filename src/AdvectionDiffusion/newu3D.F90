!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!              CALCULATION OF NEW FACE NORMAL VELOCITY                !
!                           Nov 2016                                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE newu3D(fun,funv,flag,            &
                        dpdx,dpdy,dpds,           &
                        xc,yc,sig,dsig,No_cp,nbe, &
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
      real*8, dimension(:,:) :: fun(N_CELL,NZ)
      real*8, dimension(:,:) :: funv(N_VERT,NZ-1)
!     --------------------------------------
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ)
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
      integer :: flag
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      real*8, dimension(:,:),allocatable :: dpdn1,dpdn2,dpdn3
      real*8, dimension(:,:),allocatable :: dpdnT,dpdnB
      real*8, dimension(:) :: dfundn(3)
      real*8, dimension(:) :: phi(3)
      real*8 :: dfundT,dfundB,phiT,phiB
      real*8 :: funA,funB,delta1,delta2
      real*8 :: diff1,diff2,diff3,co1,co2,co3
      real*8 :: c1
      integer :: jc,jv1,jv2,jj
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: newU3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
!      For Adams-Bashforth time stepping
       c1 = dt

#      ifdef KeyOldPressure
       if(flag .eq. 1)  c1 = 1.55d0*dt
#      endif
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|
      allocate(dpdn1(N_CELL0,NZ),         &
               dpdn2(N_CELL0,NZ),         &
               dpdn3(N_CELL0,NZ),         &
               dpdnT(N_CELL0,NZ),         &
               dpdnB(N_CELL0,NZ))
!*********************************************************************!
!                                                                     !
!                pressure gradient at wall normal                     !
!                                                                     !
!*********************************************************************!

      do k=1,NZ
         do i=1,N_CELL0
            dpdn1(i,k) = 0.0d0
            dpdn2(i,k) = 0.0d0
            dpdn3(i,k) = 0.0d0
            dpdnT(i,k) = 0.0d0
            dpdnB(i,k) = 0.0d0
         enddo
      enddo

#    ifdef KeyKimChoi

     call BCglobalP(fun,funv,xc,yc,sig,dsig,No_cp,nbe,xv,yv,sigv,dsigv,No_vp,nbev,2)

     do k=2,NZ-1
        do i=1,N_CELL0
!         __________________________________________
!         For XY Plane
          do j =1,3
                jc = No_cp(i,j)
                jj = j+1
                if(jj .gt. 3) jj = jj -3
                jv1 = No_vp(i,j)
                jv2 = No_vp(i,jj)
                co1 = dle1(i,j)
                co3 = dle3(i,j)
                if(bcpc(jc) .ne. 1) then
                    co2 = dle2(i,j)
                else
                    dfundn(j) = 0.
                    goto 200
                endif
!             --------------
!              principle diffusion
                diff1 = (fun(jc,k) -fun(i,k))/dlCC(i,j)
!             --------------
!              cross diffusion
                funA = funv(jv1,k)
                funB = funv(jv2,k)
                diff2 = (funB-funA)/dlVV(i,j)
                diff3 = (xmo(i,j)-xme(i,j))*(dpdx(jc,k)-dpdx(i,k)) &
                       +(ymo(i,j)-yme(i,j))*(dpdy(jc,k)-dpdy(i,k))
!             --------------
!             both principle and cross term
                dfundn(j) = diff1*co1 + diff2*co2 + diff3*co3
!             --------------
200             continue
         enddo

             dpdn1(i,k) = dfundn(1)
             dpdn2(i,k) = dfundn(2)
             dpdn3(i,k) = dfundn(3)
!         ___________________________________________
!         For top and bottom
             dpdnT(i,k) = (fun(i,k+1)-fun(i,k))/dsig(k)
             dpdnB(i,k) = (fun(i,k-1)-fun(i,k))/dsig(k)
        enddo
     enddo

#    else

     do k=2,NZ-1
        do i=1,N_CELL0
!         __________________________________________
!         For XY Plane
        do j =1,3
                jc = No_cp(i,j)
                co1 = 1.0d0/(xe1(i,j)*xnn(i,j)+ye1(i,j)*ynn(i,j))
                if(bcpc(jc) .ne. 1) then
                    co2 = xe2(i,j)
                    co3 = ye2(i,j)
                    delta1 = 0.5d0
                    delta2 = 0.5d0
                else
                    dfundn(j) = 0.
                    goto 100
                endif
!             --------------
!              principle diffusion
                diff1 = fun(jc,k) -fun(i,k)
!             --------------
!              cross diffusion
                diff2 = (fun(jc,k) - fun(i,k))*xe1(i,j)/dlCC(i,j)
                diff3 = (fun(jc,k) - fun(i,k))*ye1(i,j)/dlCC(i,j)
!             --------------
!             both principle and cross term
                dfundn(j) = diff1*co1 + diff2*co2 + diff3*co3
!             --------------
100             continue
         enddo

             dpdn1(i,k) = dfundn(1)
             dpdn2(i,k) = dfundn(2)
             dpdn3(i,k) = dfundn(3)
!         ___________________________________________
!         For top and bottom
             dpdnT(i,k) = (fun(i,k+1)-fun(i,k))/dsig(k)
             dpdnB(i,k) = (fun(i,k-1)-fun(i,k))/dsig(k)
        enddo
     enddo

#    endif
!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
       DO k=2, NZ-1
         do i=1, N_CELL0
           phi(1) = U1FACE(i,k)
           phi(2) = U2FACE(i,k)
           phi(3) = U3FACE(i,k)
           phiT = UTFACE(i,k)
           phiB = UBFACE(i,k)
           dfundn(1) = dpdn1(i,k)
           dfundn(2) = dpdn2(i,k)
           dfundn(3) = dpdn3(i,k)
           dfundT = dpdnT(i,k)
           dfundB = dpdnB(i,k)
           do j=1,3
             phi(j) = phi(j) - dfundn(j)*c1
           enddo
            phiT = phiT-dfundT*c1
            phiB = phiB-dfundB*c1
            U1FACE(i,k) = phi(1)
            U2FACE(i,k) = phi(2)
            U3FACE(i,k) = phi(3)
            UTFACE(i,k) = phiT
            UBFACE(i,k) = phiB
         enddo
       ENDDO

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(dpdn1,dpdn2,dpdn3, &
                 dpdnT,dpdnB)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: newU3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!===================================================================

      SUBROUTINE newu3D_test(funu,funv,funw,funp,         &
                             dpdx,dpdy,dpds,              &
                             xc,yc,sig,dsig,No_cp,nbe,    &
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
      real*8, dimension(:,:) :: funu(N_CELL,NZ)
      real*8, dimension(:,:) :: funv(N_CELL,NZ)
      real*8, dimension(:,:) :: funw(N_CELL,NZ)
      real*8, dimension(:,:) :: funp(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ)
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
      real*8, dimension(:,:),allocatable :: dpdn1,dpdn2,dpdn3
      real*8, dimension(:,:),allocatable :: dpdnT,dpdnB
      real*8, dimension(:) :: dfundn(3)
      real*8, dimension(:) :: phi(3)
      real*8 :: dfundT,dfundB,phiT,phiB
      real*8 :: funA,funB,delta1,delta2
      real*8 :: diff1,diff2,diff3
      real*8 :: c1,co1,co2,co3
      integer :: jc,jc1,jc2,jc3,Option
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: newU3D_test'
#     endif
!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
      Option = 1 ! using pressure gradient difference between centre and face
!     Option = 2 ! using intermediate velocity and pressure normal gradient
!     Option = 3 ! using latest velocity field (geometry average)
      c1 = dt

      if(Option .eq. 1) then
!      ________________________________________________________
!     |                                                        |
!     |           Face Normal Velocity interpolation           |
!     |________________________________________________________|
      call interpolateU(U1FACE,U2FACE,U3FACE,UTFACE,UBFACE, &
                        funu,funv,funw,                     &
                        xc,yc,sig,dsig,No_cp,nbe,           &
                        xv,yv,sigv,dsigv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!                correction using pressure gradient                   !
!                                                                     !
!*********************************************************************!
      allocate(dpdn1(N_CELL0,NZ),         &
               dpdn2(N_CELL0,NZ),         &
               dpdn3(N_CELL0,NZ),         &
               dpdnT(N_CELL0,NZ),         &
               dpdnB(N_CELL0,NZ))

      do k=1,NZ
         do i=1,N_CELL0
            dpdn1(i,k) = 0.0d0
            dpdn2(i,k) = 0.0d0
            dpdn3(i,k) = 0.0d0
            dpdnT(i,k) = 0.0d0
            dpdnB(i,k) = 0.0d0
         enddo
      enddo

     do k=2,NZ-1
        do i=1,N_CELL0
!         __________________________________________
!         For XY Plane
        do j =1,3
                jc = No_cp(i,j)
                    co1 = dsqrt(xe1(i,j)*xe1(i,j)+ye1(i,j)*ye1(i,j))/dlCC(i,j)
                if(bcpc(jc) .ne. 1) then
                    co2 = -xe1(i,j)
                    co3 = -ye1(i,j)
                    delta1 = 0.5d0
                    delta2 = 0.5d0
                else
                    dfundn(j) = 0.
                    goto 100
                endif
!             --------------
!              principle diffusion
                diff1 = funp(jc,k) - funp(i,k)
!             --------------
!              cross diffusion
                funA = dpdx(i,k)
                funB = dpdx(jc,k)
                diff2 = delta1*funB+delta2*funA
                funA = dpdy(i,k)
                funB = dpdy(jc,k)
                diff3 = delta1*funB+delta2*funA
!             --------------
!             both principle and cross term
                dfundn(j) = co1*diff1 + diff2*co2 + diff3*co3
!             --------------
100             continue
         enddo

             dpdn1(i,k) = dfundn(1)
             dpdn2(i,k) = dfundn(2)
             dpdn3(i,k) = dfundn(3)
!         ___________________________________________
!         For top and bottom
             dpdnT(i,k) = (funp(i,k+1)-funp(i,k))/dsig(k)
             dpdnB(i,k) = (funp(i,k-1)-funp(i,k))/dsig(k)
             dpdnT(i,k) = dpdnT(i,k)-0.5d0*(dpds(i,k+1) + dpds(i,k))
             dpdnB(i,k) = dpdnB(i,k)+0.5d0*(dpds(i,k-1) + dpds(i,k))
        enddo
     enddo
!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
       DO k=2, NZ-1
         do i=1, N_CELL0
           phi(1) = U1FACE(i,k)
           phi(2) = U2FACE(i,k)
           phi(3) = U3FACE(i,k)
           phiT = UTFACE(i,k)
           phiB = UBFACE(i,k)
           dfundn(1) = dpdn1(i,k)
           dfundn(2) = dpdn2(i,k)
           dfundn(3) = dpdn3(i,k)
           dfundT = dpdnT(i,k)
           dfundB = dpdnB(i,k)
           do j=1,3
             phi(j) = phi(j) - dfundn(j)*c1
           enddo
            phiT = phiT -dpdnT(i,k)*c1
            phiB = phiB -dpdnB(i,k)*c1
            U1FACE(i,k) = phi(1)
            U2FACE(i,k) = phi(2)
            U3FACE(i,k) = phi(3)
            UTFACE(i,k) = phiT
            UBFACE(i,k) = phiB
         enddo
       ENDDO

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(dpdn1,dpdn2,dpdn3, &
                 dpdnT,dpdnB)
      endif

       if(Option .eq. 2) then
            DO k=2, NZ-1
                do i=1, N_CELL0
                    jc1 = No_cp(i,1)
                    jc2 = No_cp(i,2)
                    jc3 = No_cp(i,3)
                    phi(1) = U1FACE(i,k)
                    phi(2) = U2FACE(i,k)
                    phi(3) = U3FACE(i,k)
                    co1 = 1.0d0/(xe1(i,1)*xnn(i,1)+ye1(i,1)*ynn(i,1))
                    co2 = 1.0d0/(xe1(i,2)*xnn(i,2)+ye1(i,2)*ynn(i,2))
                    co3 = 1.0d0/(xe1(i,3)*xnn(i,3)+ye1(i,3)*ynn(i,3))
                    dfundn(1) = co1*(funp(jc1,k)-funp(i,k))
                    dfundn(2) = co2*(funp(jc2,k)-funp(i,k))
                    dfundn(3) = co3*(funp(jc3,k)-funp(i,k))

                    do j=1,3
                        phi(j) = phi(j) - dfundn(j)*c1
                    enddo

                    U1FACE(i,k) = phi(1)
                    U2FACE(i,k) = phi(2)
                    U3FACE(i,k) = phi(3)
                    UTFACE(i,k) = 0.
                    UBFACE(i,k) = 0.

                enddo
            ENDDO
       endif

      if(Option .eq. 3) then
!      ________________________________________________________
!     |                                                        |
!     |           Face Normal Velocity interpolation           |
!     |________________________________________________________|
      call interpolateU(U1FACE,U2FACE,U3FACE,UTFACE,UBFACE, &
                        funu,funv,funw,                     &
                        xc,yc,sig,dsig,No_cp,nbe,           &
                        xv,yv,sigv,dsigv,No_vp,nbev)
      endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: newU3D_test'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
