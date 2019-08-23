!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   Update the RHS for PPE Solver                     !
!                              Oct 2016                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Source3D(source,phiu,phiv,phiw,      &
                          xc,yc,sig,dsig,No_cp,nbe,   &
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
!     --------------------------------------
      real*8, dimension(:,:) :: source(N_CELL,NZ)
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
!     --------------------------------------
      real*8, dimension(:,:), allocatable :: dudx,dudy,duds
      real*8, dimension(:,:), allocatable :: dvdx,dvdy,dvds
      real*8, dimension(:,:), allocatable :: dwdx,dwdy,dwds
      real*8, dimension(:) :: phi(1:3)
      real*8 :: phiT,phiB,Vol,resid
      integer :: elem,Option
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Source3D'
#     endif
!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
      Option = 1 ! Use face normal velocity
!      Option = 2 ! Use LSM method to get gradient

     if(Option .eq. 1) then

     call  interpolateU(U1FACE,U2FACE,U3FACE,UTFACE,UBFACE,   &
                         phiu,phiv,phiw,                      &
                         xc,yc,sig,dsig,No_cp,nbe,            &
                         xv,yv,sigv,dsigv,No_vp,nbev)
!   ---------------------------------------------------------
     do k=1,NZ
       do i=1,N_CELL
          source(i,k) = 0.
       enddo
     enddo
!    -------------------------
      do k=2,NZ-1
         do i=1,N_CELL0
             resid = 0.0d0
             phi(1) = U1FACE(i,k)
             phi(2) = U2FACE(i,k)
             phi(3) = U3FACE(i,k)
             phiT = UTFACE(i,k)
             phiB = UBFACE(i,k)
!           ___________________________________________________
!           Horizontal neighbors
              do j=1,3
                    resid = resid + phi(j)*dlVV(i,j)*dsigv(k-1)
	          enddo
!           ___________________________________________________
!           Top & bottom neighbors
            IF(I3D .EQ. 1) THEN
               resid = resid + phiT*areaCell(i)
               resid = resid + phiB*areaCell(i)
            ENDIF

              source(i,k) = resid
        enddo
      enddo

      endif


    if(Option .eq. 2) then

        allocate(dudx(N_CELL,NZ),dudy(N_CELL,NZ),duds(N_CELL,NZ), &
                 dvdx(N_CELL,NZ),dvdy(N_CELL,NZ),dvds(N_CELL,NZ), &
                 dwdx(N_CELL,NZ),dwdy(N_CELL,NZ),dwds(N_CELL,NZ))

        call grandientLSM(dudx,dudy,duds,phiu,xc,yc,sig,dsig,No_cp,nbe)
        call grandientLSM(dvdx,dvdy,dvds,phiv,xc,yc,sig,dsig,No_cp,nbe)
        call grandientLSM(dwdx,dwdy,dwds,phiw,xc,yc,sig,dsig,No_cp,nbe)

        do k=1,NZ
            do i=1,N_CELL
                source(i,k) = 0.
            enddo
        enddo

        DO k=2,NZ-1
            do i=1,N_CELL0

                Vol = areaCell(i)*dsigv(k-1)
                resid = dudx(i,k) + dvdy(i,k) + dwds(i,k)
                source(i,k) = resid*Vol
            enddo
        ENDDO

        deallocate(dudx,dudy,duds, &
                   dvdx,dvdy,dvds, &
                   dwdx,dwdy,dwds)

    endif
!*********************************************************************!
!                                                                     !
!                    Finalization of the subroutine                   !
!                                                                     !
!*********************************************************************!
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Source3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!   ===================================================================
      SUBROUTINE update_rhsp(rhs,phi,phiv,dphidx,dphidy,dphids,  &
                             xc,yc,sig,dsig,No_cp,nbe,           &
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
!     --------------------------------------
      real*8, dimension(:,:) :: rhs(N_CELL,NZ)
      real*8, dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:,:) :: dphidx(N_CELL,NZ)
      real*8, dimension(:,:) :: dphidy(N_CELL,NZ)
      real*8, dimension(:,:) :: dphids(N_CELL,NZ)
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      real*8, dimension(:) :: dphidn(1:3)
      real*8 :: funA,funB,co2,co3,diff2,diff3,source
      integer ::jc,jj,jv1,jv2
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: update_rhsp'
#     endif
!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

#     ifdef KeyKimChoi

        call BCglobalP(phi,phiv,xc,yc,sig,dsig,No_cp,nbe, &
                       xv,yv,sigv,dsigv,No_vp,nbev,2)

         do k=2, NZ-1
            do i=1, N_CELL0
                do j=1,3
                    jc = No_cp(i,j)
                    co2 = dle2(i,j)
                    co3 = dle3(i,j)
                    jj = j+1
                    if(jj .gt. 3) jj = jj-3
                    jv1 = No_vp(i,j)
                    jv2 = No_vp(i,jj)
!              -----------------
!              cross diffusion term for dphidn
                    funA = phiv(jv1,k)
                    funB = phiv(jv2,k)
                    diff2 = (funB-funA)/dlVV(i,j)
                    diff3 = (xmo(i,j)-xme(i,j))*(dphidx(jc,k)-dphidx(i,k)) &
                           +(ymo(i,j)-yme(i,j))*(dphidy(jc,k)-dphidy(i,k))
                    dphidn(j) =  diff2*co2 + diff3*co3
                enddo
!          -------------------
!          calculate extra source term by mass flux
               source =  dphidn(1)*dlVV(i,1)*dsigv(k-1) + &
                         dphidn(2)*dlVV(i,2)*dsigv(k-1) + &
                         dphidn(3)*dlVV(i,3)*dsigv(k-1)
!          -----------------------------------------------
!          update the right hand side term
               rhs(i,k) = rhs(i,k) - source
             enddo
            enddo
#     else
!      ________________________________________________________
!     |                                                        |
!     |                     Main routine                       |
!     |________________________________________________________|

     do k=2, NZ-1
        do i=1, N_CELL0
            do j=1,3
                jc = No_cp(i,j)
                co2 = xe2(i,j)
                co3 = ye2(i,j)
!              -----------------
!              cross diffusion term for dphidn
                funA = dphidx(i,k)
                funB = dphidx(jc,k)
                diff2 = 0.5d0*(funB+funA)
                funA = dphidy(i,k)
                funB = dphidy(jc,k)
                diff3 = 0.5d0*(funB+funA)
                dphidn(j) =  diff2*co2 + diff3*co3
            enddo
!          -------------------
!          calculate extra source term by mass flux
           source =  dphidn(1)*dlVV(i,1)*dsigv(k-1) + &
                     dphidn(2)*dlVV(i,2)*dsigv(k-1) + &
                     dphidn(3)*dlVV(i,3)*dsigv(k-1)
!          -----------------------------------------------
!          update the right hand side term
           rhs(i,k) = rhs(i,k) - source
         enddo
        enddo
#     endif

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: update_rhsp'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
