!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                  INTERPOLATION THE VALUE AT FACE CENTER             !
!                           (Geometry Average)                        !
!                              MAR 2016                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE calcul_face(funf1,funf2,funf3,funfT,funfB, &
                             fun,funv,                      &
                             xc,yc,sig,dsig,No_cp,nbe,      &
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
!     -------------------------------------
      real*8, dimension(:,:) :: funf1(N_CELL0,NZ)
      real*8, dimension(:,:) :: funf2(N_CELL0,NZ)
      real*8, dimension(:,:) :: funf3(N_CELL0,NZ)
      real*8, dimension(:,:) :: funfT(N_CELL0,NZ)
      real*8, dimension(:,:) :: funfB(N_CELL0,NZ)
      real*8, dimension(:,:) :: fun(N_CELL,NZ)
      real*8, dimension(:,:) :: funv(N_VERT,NZ-1)
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
      real*8, dimension(:,:),allocatable :: dphidx,dphidy,dphids
      real*8, dimension(:) :: phi(1:3)
      real*8, dimension(:) :: phif(1:3)
      real*8 :: delta1,delta2,phi1,phi2
      integer :: jc,jc1,jc2,jc3
      integer :: jj,nv1,nv2
      real*8  :: extra,phifT,phifB
      real*8 :: Vol, residx,residy

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: calcul_face'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#      ifdef KeyKimChoi
            DO k=2,NZ-1
                do i=1,N_CELL0
                    do j=1,3
!               ---------------
!               principle value (central at intersection point)
                    jc = No_cp(i,j)
                    jj = j+1
                    if(jj .gt. 3) jj = jj- 3
                    nv1 = No_vp(i,j)
                    nv2 = No_vp(i,jj)
                    phi1 = fun(i,k)
                    phi2 = fun(jc,k)
                    delta1 = 0.5d0
                    delta2 = 0.5d0
                    extra = doe(i,j)*(funv(nv2,k)-funv(nv1,k))/dlVV(i,j)
                    phi(j) = phi2*delta2+phi1*delta1+extra
                    enddo
!             ----------------
!              Final value
                    funf1(i,k) = phi(1)
                    funf2(i,k) = phi(2)
                    funf3(i,k) = phi(3)
                    funfT(i,k) = 0.5d0*(fun(i,k)+fun(i,k+1))
                    funfB(i,k) = 0.5d0*(fun(i,k)+fun(i,k-1))
                enddo
            ENDDO


#      else
!       Option = 1  ! geometry average
!       Option = 2  ! use gradient contribution
!*********************************************************************!
!                                                                     !
!                          Mass flux calculation                      !
!                                                                     !
!*********************************************************************!
        if(Option .eq. 1) then

            DO k=2,NZ-1
                do i=1,N_CELL0
                    do j=1,3
!               ---------------
!               principle value (central at intersection point)
                    jc = No_cp(i,j)
                    phi1 = fun(i,k)
                    phi2 = fun(jc,k)
                    delta1 = 0.5d0
                    delta2 = 0.5d0
                    phi(j) = (phi2*delta2+phi1*delta1)
                    enddo
!             ----------------
!              Final value
                    funf1(i,k) = phi(1)
                    funf2(i,k) = phi(2)
                    funf3(i,k) = phi(3)
                    funfT(i,k) = 0.5d0*(fun(i,k)+fun(i,k+1))
                    funfB(i,k) = 0.5d0*(fun(i,k)+fun(i,k-1))
                enddo
            ENDDO
        endif
!*********************************************************************!
!                                                                     !
!                    get gradient for face center                     !
!                                                                     !
!*********************************************************************!

        if(Option .eq. 2) then

         allocate(dphidx(N_CELL,NZ),dphidy(N_CELL,NZ),dphids(N_CELL,NZ))

         call grandientLSM(dphidx,dphidy,dphids,fun,xc,yc,sig,dsig,No_cp,nbe)

            DO k=2,NZ-1
                do i=1,N_CELL0
#                   ifdef KeyAdvUpwind
                        phif(1) = U1FACE(i,k)
                        phif(2) = U2FACE(i,k)
                        phif(3) = U3FACE(i,k)
                        phifT = UTFACE(i,k)
                        phifB = UBFACE(i,k)
#                   endif
                    do j=1,3
!               ---------------
!               principle value (central at intersection point)
                        jc = No_cp(i,j)
                        phi1 = fun(i,k)                          &
!                            + dphidx(i,k)*(xme(i,j) - xmo(i,j))  &
!                            + dphidy(i,k)*(yme(i,j) - ymo(i,j))

                            + dphidx(i,k)*(xme(i,j)-xc(i))   &
                            + dphidy(i,k)*(yme(i,j)-yc(i))
                        phi2 = fun(jc,k)                           &
!                            + dphidx(jc,k)*(xme(i,j) - xmo(i,j))   &
!                            + dphidy(jc,k)*(yme(i,j) - ymo(i,j))

                            + dphidx(jc,k)*(xme(i,j)-xc(jc))  &
                            + dphidy(jc,k)*(yme(i,j)-yc(jc))
!                       -------------------------
#                   ifdef KeyAdvUpwind
                        if(phif(j) .gt. 0.0d0) then
                            phi(j) = phi1
                        else
                            phi(j) = phi2
                        endif
#                   else
                        if(tagsp(jc) .eq. 1) then
                            phi(j) = 0.5d0*(fun(i,k)+fun(jc,k))
                        else
                            phi(j) = 0.5d0*(phi2+phi1)
                        endif
#                   endif
                    enddo
!                   ----------------
!                   Final value
                    funf1(i,k) = phi(1)
                    funf2(i,k) = phi(2)
                    funf3(i,k) = phi(3)
#                  ifdef KeyAdvUpwind
                     if(phifT .gt. 0.0d0) then
                         funfT(i,k) = fun(i,k) + dphids(i,k)*(sigv(k)-sig(k))
                     else
                         funfT(i,k) = fun(i,k+1) + dphids(i,k-1)*(sigv(k)-sig(k+1))
                     endif

                     if(phifB .gt. 0.0d0) then
                         funfB(i,k) = fun(i,k) + dphids(i,k)*(sigv(k-1)-sig(k))
                     else
                         funfB(i,k) = fun(i,k-1) + dphids(i,k-1)*(sigv(k-1)-sig(k-1))
                     endif
#                  else
                    funfT(i,k) = 0.5d0*(fun(i,k)+fun(i,k+1))
                    funfB(i,k) = 0.5d0*(fun(i,k)+fun(i,k-1))
#                  endif
                enddo
            ENDDO

            deallocate(dphidx,dphidy,dphids)
      endif

#     endif
!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: calcul_face'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
