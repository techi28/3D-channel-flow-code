!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!     SOLUTION OF THE ADVECTION-DIFFUSION PROBLEM WITH FREE SURFACE   !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AdvDiffSolve(funu,funv,funw,                   &
                              funuv,funvv,funwv,                &
                              rhsu,rhsv,rhsw,                   &
                              Gamx,Gamy,Gamz,                   &
                              xc,yc,sig,dsig,No_cp,nbe,         &
                              xv,yv,sigv,dsigv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

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
      real*8,dimension(:,:) :: funu(N_CELL,NZ)
      real*8,dimension(:,:) :: funv(N_CELL,NZ)
      real*8,dimension(:,:) :: funw(N_CELL,NZ)
      real*8,dimension(:,:) :: funuv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: funvv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: funwv(N_VERT,NZ-1)
!     ----------------------------------------
      real*8,dimension(:,:) :: rhsu(N_CELL,NZ)
      real*8,dimension(:,:) :: rhsv(N_CELL,NZ)
      real*8,dimension(:,:) :: rhsw(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamz(N_CELL,NZ)
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
!     ----------------------------------------
!     Local variables
      real*8, dimension(:,:),allocatable :: Au0,Au1,Au2,Au3,AuT,AuB
      real*8, dimension(:,:),allocatable :: UFi1,UFi2,UFi3,UFiT,UFiB
      real*8, dimension(:) :: AA(1:3)
      real*8 :: AA0,AAB,AAT
      real*8 :: co1,Vol
      integer :: jc,IDisplay
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'Begin subroutine: AdvDiffSOLVE'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IDisplay = 1
#     ifdef KeyParallel
        if(rang_topo .ne. 0) IDisplay = 0
#     endif
!*********************************************************************!
!                                                                     !
!                         allocate matrix                             !
!                                                                     !
!*********************************************************************!
       allocate(Au0(N_CELL0,NZ),Au1(N_CELL0,NZ),Au2(N_CELL0,NZ), &
                Au3(N_CELL0,NZ),AuT(N_CELL0,NZ),AuB(N_CELL0,NZ))
!*********************************************************************!
!                                                                     !
!                           construct A                               !
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
!                   advection contribution                            !
!                                                                     !
!*********************************************************************!
      do k=2, NZ-1
        do i=1,N_CELL0
!        -----------------------------
!        rhs
          Vol = areaCell(i)*dsigv(k-1)
!         ----------------------------
!         For u
#         ifdef KeyAdvUpwind

          Au1(i,k) = 0.5d0*dt*UFi1(i,k)/Vol
          Au2(i,k) = 0.5d0*dt*UFi2(i,k)/Vol
          Au3(i,k) = 0.5d0*dt*UFi3(i,k)/Vol
          AuT(i,k) = 0.5d0*dt*UFiT(i,k)/Vol
          AuB(i,k) = 0.5d0*dt*UFiB(i,k)/Vol

!          Au1(i,k) = dt*UFi1(i,k)/Vol
!          Au2(i,k) = dt*UFi2(i,k)/Vol
!          Au3(i,k) = dt*UFi3(i,k)/Vol
!          AuT(i,k) = dt*UFiT(i,k)/Vol
!          AuB(i,k) = dt*UFiB(i,k)/Vol

          Au0(i,k) = 1.0d0

          if(UFi1(i,k) .gt. 0.0d0) then
            Au0(i,k) = Au0(i,k) + Au1(i,k)
            Au1(i,k) = 0.0d0
          endif

          if(UFi2(i,k) .gt. 0.0d0) then
            Au0(i,k) = Au0(i,k) + Au2(i,k)
            Au2(i,k) = 0.0d0
          endif

          if(UFi3(i,k) .gt. 0.0d0) then
            Au0(i,k) = Au0(i,k) + Au3(i,k)
            Au3(i,k) = 0.0d0
          endif

          if(UFiT(i,k) .gt. 0.0d0) then
            Au0(i,k) = Au0(i,k) + AuT(i,k)
            AuT(i,k) = 0.0d0
          endif

          if(UFiB(i,k) .gt. 0.0d0) then
            Au0(i,k) = Au0(i,k) + AuB(i,k)
            AuB(i,k) = 0.0d0
          endif

#         else
          Au1(i,k) = 0.25d0*dt*UFi1(i,k)/Vol
          Au2(i,k) = 0.25d0*dt*UFi2(i,k)/Vol
          Au3(i,k) = 0.25d0*dt*UFi3(i,k)/Vol
          AuT(i,k) = 0.25d0*dt*UFiT(i,k)/Vol
          AuB(i,k) = 0.25d0*dt*UFiB(i,k)/Vol

          Au0(i,k) =  Au1(i,k) + Au2(i,k) + Au3(i,k)  &
                    + AuT(i,k) + AuB(i,k) + 1.0d0
#        endif
         enddo
      enddo
!*********************************************************************!
!                                                                     !
!                   diffusion contribution                            !
!                                                                     !
!*********************************************************************!
    do k=2,NZ-1
        do i=1,N_CELL0
!         ---------------------------
          Vol = areaCell(i)*dsigv(k-1)
!         __________________________________________
!         For XY Plane
           do j =1, 3
                jc = No_cp(i,j)
                co1 = dsqrt(xe1(i,j)*xe1(i,j)+ye1(i,j)*ye1(i,j))/dlCC(i,j)
!             --------------
!              principle diffusion
                AA(j) =  co1*dlVV(i,j)*dsigv(k-1)
           enddo
!         __________________________________________
!         Bottom
                 AAB = areaCell(i)/dsig(k-1)
!         __________________________________________
!         Top
                 AAT = areaCell(i)/dsig(k)
!         __________________________________________
!         coef for center
            AA0 = -1.0d0*(AAT+AAB+AA(1)+AA(2)+AA(3))

!          |----------------------------------------------------------|
!          |              Final coefficients contributions            |
!          |----------------------------------------------------------|
	        Au0(i,k) = Au0(i,k) - 0.5d0*dt*Gamx(i,k)*AA0/Vol
            Au1(i,k) = Au1(i,k) - 0.5d0*dt*Gamx(i,k)*AA(1)/Vol
            Au2(i,k) = Au2(i,k) - 0.5d0*dt*Gamx(i,k)*AA(2)/Vol
            Au3(i,k) = Au3(i,k) - 0.5d0*dt*Gamx(i,k)*AA(3)/Vol
            AuT(i,k) = AuT(i,k) - 0.5d0*dt*Gamx(i,k)*AAT/Vol
            AuB(i,k) = AuB(i,k) - 0.5d0*dt*Gamx(i,k)*AAB/Vol

         enddo
      ENDDO
!*********************************************************************!
!                                                                     !
!                         solve Am = b                                !
!                                                                     !
!*********************************************************************!
      if(IDisplay .eq. 1) print*, '       SOLVE FOR U FIELD'
       call BICGSTAB_SOLVE(Au0,Au1,Au2,Au3,AuT,AuB,           &
                           rhsu,funu,                         &
                           xc,yc,sig,dsig,No_cp,nbe,tagBCu)
      if(IDisplay .eq. 1) print*, '       SOLVE FOR V FIELD'
       call BICGSTAB_SOLVE(Au0,Au1,Au2,Au3,AuT,AuB,           &
                           rhsv,funv,                         &
                           xc,yc,sig,dsig,No_cp,nbe,tagBCv)
      IF(I3D .EQ. 1) THEN
        if(IDisplay .eq. 1) print*, '       SOLVE FOR W FIELD'
        call BICGSTAB_SOLVE(Au0,Au1,Au2,Au3,AuT,AuB,       &
                            rhsw,funw,                     &
                            xc,yc,sig,dsig,No_cp,nbe,tagBCw)
      ENDIF
!*********************************************************************!
!                                                                     !
!                         finalize routine                            !
!                                                                     !
!*********************************************************************!
        deallocate(Au0,Au1,Au2,Au3,AuT,AuB)
        deallocate(UFi1,UFi2,UFi3,UFiT,UFiB)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), ' End   subroutine AdvDiffSOLVE'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF Advection-Diffusion 3D                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!   ===================================================================
      SUBROUTINE update_rhsf(rhsu,rhsv,rhsw,                    &
                             phiu,phiv,phiw,                    &
                             phiuv,phivv,phiwv,                 &
                             dpdx,dpdy,dpds,                    &
                             Gamx,Gamy,Gamz,                    &
                             xc,yc,sig,dsig,No_cp,nbe,          &
                             xv,yv,sigv,dsigv,No_vp,nbev)
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the advection contribution to the   !
!    general linear system. We have two options to the final values   !
!    of Am and AmG, depending if we choose between an implicit or     !
!    explicit scheme (cppdefs.h).                                     !
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
!    |      Declaration of variables      |
!    |____________________________________|
      real*8, dimension(:,:) :: rhsu(N_CELL,NZ)
      real*8, dimension(:,:) :: rhsv(N_CELL,NZ)
      real*8, dimension(:,:) :: rhsw(N_CELL,NZ)
      real*8, dimension(:,:) :: phiu(N_CELL,NZ)
      real*8, dimension(:,:) :: phiv(N_CELL,NZ)
      real*8, dimension(:,:) :: phiw(N_CELL,NZ)
      real*8, dimension(:,:) :: phiuv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: phivv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: phiwv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamz(N_CELL,NZ)
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
      real*8, dimension(:,:), allocatable :: dudx,dudy,duds
      real*8, dimension(:,:), allocatable :: dvdx,dvdy,dvds
      real*8, dimension(:,:), allocatable :: dwdx,dwdy,dwds
      real*8, dimension(:,:), allocatable :: source
      integer :: jc1,jc2,jc3
      real*8  :: Vol,residu,residv,residw
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: update_rhsf'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#       ifdef KeyOldPressure
         do k=1,NZ
            do i=1,N_CELL
                  rhsu(i,k) =  bdx - dpdx(i,k)
                  rhsv(i,k) =  bdy - dpdy(i,k)
                  rhsw(i,k) =  bdz - dpds(i,k)
            enddo ! end i loop
         enddo ! end k loop
#       else
         do k=1,NZ
            do i=1,N_CELL
                  rhsu(i,k) =  bdx
                  rhsv(i,k) =  bdy
                  rhsw(i,k) =  bdz
            enddo ! end i loop
         enddo ! end k loop
#       endif

!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
#       ifdef KeyFullStress
        allocate(dudx(N_CELL,NZ),dudy(N_CELL,NZ),duds(N_CELL,NZ), &
                 dvdx(N_CELL,NZ),dvdy(N_CELL,NZ),dvds(N_CELL,NZ), &
                 dwdx(N_CELL,NZ),dwdy(N_CELL,NZ),dwds(N_CELL,NZ))

        call grandientLSM(dudx,dudy,duds,phiu,xc,yc,sig,dsig,No_cp,nbe)
        call grandientLSM(dvdx,dvdy,dvds,phiv,xc,yc,sig,dsig,No_cp,nbe)
        call grandientLSM(dwdx,dwdy,dwds,phiw,xc,yc,sig,dsig,No_cp,nbe)

        if(ZPB .ne. 1) then
            do i=1,N_CELL0
                if(bctop .eq. 1) then
                    dwdx(i,NZ) = -1.0d0*dwdx(i,NZ-1)
                    dwdy(i,NZ) = -1.0d0*dwdy(i,NZ-1)
                endif

                if(bcbot .eq. 1) then
                    dwdx(i,1) = -1.0d0*dwdx(i,2)
                    dwdy(i,1) = -1.0d0*dwdy(i,2)
                endif
            enddo
        endif
!      ________________________________________________________
!     |                                                        |
!     |                     Main routine                       |
!     |________________________________________________________|

        do k=2,NZ-1
            do i=1,N_CELL0
               Vol = areaCell(i)*dsigv(k-1)
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
           !---------------
           ! for u
             residu =  0.5d0*(dudx(i,k)+dudx(jc1,k))*dlVV(i,1)*dsigv(k-1)*xnn(i,1) &
                     + 0.5d0*(dvdx(i,k)+dvdx(jc1,k))*dlVV(i,1)*dsigv(k-1)*ynn(i,1) &
                     + 0.5d0*(dudx(i,k)+dudx(jc2,k))*dlVV(i,2)*dsigv(k-1)*xnn(i,2) &
                     + 0.5d0*(dvdx(i,k)+dvdx(jc2,k))*dlVV(i,2)*dsigv(k-1)*ynn(i,2) &
                     + 0.5d0*(dudx(i,k)+dudx(jc3,k))*dlVV(i,3)*dsigv(k-1)*xnn(i,3) &
                     + 0.5d0*(dvdx(i,k)+dvdx(jc3,k))*dlVV(i,3)*dsigv(k-1)*ynn(i,3) &
                     + 0.5d0*(dwdx(i,k)+dwdx(i,k+1))*areaCell(i) &
                     - 0.5d0*(dwdx(i,k)+dwdx(i,k-1))*areaCell(i)
           !---------------
           ! for v
             residv =  0.5d0*(dudy(i,k)+dudy(jc1,k))*dlVV(i,1)*dsigv(k-1)*xnn(i,1) &
                     + 0.5d0*(dvdy(i,k)+dvdy(jc1,k))*dlVV(i,1)*dsigv(k-1)*ynn(i,1) &
                     + 0.5d0*(dudy(i,k)+dudy(jc2,k))*dlVV(i,2)*dsigv(k-1)*xnn(i,2) &
                     + 0.5d0*(dvdy(i,k)+dvdy(jc2,k))*dlVV(i,2)*dsigv(k-1)*ynn(i,2) &
                     + 0.5d0*(dudy(i,k)+dudy(jc3,k))*dlVV(i,3)*dsigv(k-1)*xnn(i,3) &
                     + 0.5d0*(dvdy(i,k)+dvdy(jc3,k))*dlVV(i,3)*dsigv(k-1)*ynn(i,3) &
                     + 0.5d0*(dwdy(i,k)+dwdy(i,k+1))*areaCell(i) &
                     - 0.5d0*(dwdy(i,k)+dwdy(i,k-1))*areaCell(i)
           !---------------
           ! for w
             residw =  0.5d0*(duds(i,k)+duds(jc1,k))*dlVV(i,1)*dsigv(k-1)*xnn(i,1) &
                     + 0.5d0*(dvds(i,k)+dvds(jc1,k))*dlVV(i,1)*dsigv(k-1)*ynn(i,1) &
                     + 0.5d0*(duds(i,k)+duds(jc2,k))*dlVV(i,2)*dsigv(k-1)*xnn(i,2) &
                     + 0.5d0*(dvds(i,k)+dvds(jc2,k))*dlVV(i,2)*dsigv(k-1)*ynn(i,2) &
                     + 0.5d0*(duds(i,k)+duds(jc3,k))*dlVV(i,3)*dsigv(k-1)*xnn(i,3) &
                     + 0.5d0*(dvds(i,k)+dvds(jc3,k))*dlVV(i,3)*dsigv(k-1)*ynn(i,3) &
                     + 0.5d0*(dwds(i,k)+dwds(i,k+1))*areaCell(i) &
                     - 0.5d0*(dwds(i,k)+dwds(i,k-1))*areaCell(i)
           !---------------
           ! final contribution
             rhsu(i,k) = rhsu(i,k) + Gamx(i,k)*residu/Vol
             rhsv(i,k) = rhsv(i,k) + Gamy(i,k)*residv/Vol
             rhsw(i,k) = rhsw(i,k) + Gamz(i,k)*residw/Vol
            enddo
        enddo

        deallocate(dudx,dudy,duds, &
                   dvdx,dvdy,dvds, &
                   dwdx,dwdy,dwds)

!        allocate(source(N_CELL,NZ))
!        do k=2,NZ-1
!            do i=1,N_CELL0
!                source(i,k) = U1FACE(i,k)*dlVV(i,1)*dsigv(k-1)  &
!                            + U2FACE(i,k)*dlVV(i,2)*dsigv(k-1)  &
!                            + U3FACE(i,k)*dlVV(i,3)*dsigv(k-1)  &
!                            + UTFACE(i,k)*areaCell(i)           &
!                            + UBFACE(i,k)*areaCell(i)
!                source(i,k) = source(i,k)/(areaCell(i)*dsigv(k-1))
!            enddo
!        enddo
!
!        call BCpcenter(source,xc,yc,sig,dsig,No_cp,nbe)
!
!
!        if(ZPB .ne.1) then
!            if(bctop .eq. 1) then
!                if(TopBC .eq. 0) then
!                    do i=1,N_CELL0
!                        source(i,NZ) = -U1FACE(i,NZ-1)*dlVV(i,1)*dsigv(NZ-1)  &
!                                      - U2FACE(i,NZ-1)*dlVV(i,2)*dsigv(NZ-1)  &
!                                      - U3FACE(i,NZ-1)*dlVV(i,3)*dsigv(NZ-1)  &
!                                      + UTFACE(i,NZ-1)*areaCell(i)            &
!                                      + UBFACE(i,NZ-1)*areaCell(i)
!                        source(i,NZ) = source(i,NZ)/(areaCell(i)*dsigv(NZ-1))
!                    enddo
!                endif
!            endif
!
!            if(bcbot .eq. 1) then
!                if(BotBC .eq. 0) then
!                    do i=1,N_CELL0
!                        source(i,1) =  -U1FACE(i,2)*dlVV(i,1)*dsigv(1)  &
!                                      - U2FACE(i,2)*dlVV(i,2)*dsigv(1)  &
!                                      - U3FACE(i,2)*dlVV(i,3)*dsigv(1)  &
!                                      + UTFACE(i,2)*areaCell(i)         &
!                                      + UBFACE(i,2)*areaCell(i)
!                        source(i,1) = source(i,1)/(areaCell(i)*dsigv(1))
!                    enddo
!                endif
!            endif
!        endif
!
!        do k=2,NZ-1
!            do i=1,N_CELL0
!               Vol = 3.0d0*areaCell(i)*dsigv(k-1)
!               jc1 = No_cp(i,1)
!               jc2 = No_cp(i,2)
!               jc3 = No_cp(i,3)
!           !---------------
!           ! for u
!             residu =  0.5d0*(source(i,k)+source(jc1,k))*dlVV(i,1)*dsigv(k-1)*xnn(i,1) &
!                     + 0.5d0*(source(i,k)+source(jc2,k))*dlVV(i,2)*dsigv(k-1)*xnn(i,2) &
!                     + 0.5d0*(source(i,k)+source(jc3,k))*dlVV(i,3)*dsigv(k-1)*xnn(i,3)
!           !---------------
!           ! for v
!             residv =  0.5d0*(source(i,k)+source(jc1,k))*dlVV(i,1)*dsigv(k-1)*ynn(i,1) &
!                     + 0.5d0*(source(i,k)+source(jc2,k))*dlVV(i,2)*dsigv(k-1)*ynn(i,2) &
!                     + 0.5d0*(source(i,k)+source(jc3,k))*dlVV(i,3)*dsigv(k-1)*ynn(i,3)
!           !---------------
!           ! for w
!             residw = 0.5d0*(source(i,k)+source(i,k+1))*areaCell(i) &
!                     - 0.5d0*(source(i,k)+source(i,k-1))*areaCell(i)
!           !---------------
!           ! final contribution
!             rhsu(i,k) = rhsu(i,k) - 2.0d0*Gamx(i,k)*residu/Vol
!             rhsv(i,k) = rhsv(i,k) - 2.0d0*Gamy(i,k)*residv/Vol
!             rhsw(i,k) = rhsw(i,k) - 2.0d0*Gamz(i,k)*residw/Vol
!            enddo
!        enddo
!       deallocate(source)


#     endif

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: update_rhsf'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
