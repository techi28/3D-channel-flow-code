!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    GRADIENT OF PRESSURE USING FVM                   !
!                              Mar 2016                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE gradientP(fun,funv,                   &
                           dfundx,dfundy,dfundsig,     &
                           xc,yc,sig,dsig,No_cp,nbe,   &
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
      real*8 ,dimension(:,:) :: fun(N_CELL,NZ)
      real*8 ,dimension(:,:) :: funv(N_VERT,NZ-1)
      real*8 ,dimension(:,:) :: dfundx(N_CELL,NZ)
      real*8 ,dimension(:,:) :: dfundy(N_CELL,NZ)
      real*8 ,dimension(:,:) :: dfundsig(N_CELL,NZ)
!     -----------------------------------------
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
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     local variables
      real*8 :: sumfx,sumfy,co,deter
      real*8 :: sum_xcaf2,sum_ycaf2,sum_xcycaf
      real*8 :: cox,coy
      integer :: jc,Option
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: GradientP '
#     endif
!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
!      Option = 1  ! FVM approximation
#      ifndef KeyMahesh
       Option = 2  ! Least square approximation
#      else
       Option = 3  ! LSM approach used by Malesh
#      endif

      do k=1,NZ
         do i=1,N_CELL
            dfundx(i,k)  = 0.0d0
            dfundy(i,k)  = 0.0d0
            dfundsig(i,k)= 0.0d0
         enddo
      enddo
!*********************************************************************!
!                                                                     !
!                           Main subroutine                           !
!                                                                     !
!*********************************************************************!
      if(Option .eq. 1) then

         call grandientFVM(dfundx,dfundy,dfundsig,fun,funv,   &
                           xc,yc,sig,dsig,No_cp,nbe,          &
                           xv,yv,sigv,dsigv,No_vp,nbev,2)

      elseif(Option .eq. 2) then

         call grandientLSM(dfundx,dfundy,dfundsig,fun,xc,yc,sig,dsig,No_cp,nbe)
!          call grandientFVM(dfundx,dfundy,dfundsig,fun,funv,  &
!                            xc,yc,sig,dsig,No_cp,nbe,         &
!                            xv,yv,sigv,dsigv,No_vp,nbev,2)

      else

        DO k=2,NZ-1
            do i=1,N_CELL0
                sumfx = 0.0d0
                sumfy = 0.0d0
                sum_xcaf2  = 0.
                sum_ycaf2  = 0.
                sum_xcycaf = 0.
                do j=1,3
                    jc = No_cp(i,j)
                    sumfx = sumfx + xnn(i,j)*dle1(i,j)*(fun(jc,k)-fun(i,k))/dlCC(i,j)
                    sumfy = sumfy + ynn(i,j)*dle1(i,j)*(fun(jc,k)-fun(i,k))/dlCC(i,j)
                    sum_xcaf2 = sum_xcaf2 + xnn(i,j)*xnn(i,j)
                    sum_ycaf2 = sum_ycaf2 + ynn(i,j)*ynn(i,j)
                    sum_xcycaf = sum_xcycaf + xnn(i,j)*ynn(i,j)
                enddo
                deter = sum_xcaf2*sum_ycaf2-sum_xcycaf*sum_xcycaf
                dfundx(i,k)=(sum_ycaf2*sumfx-sum_xcycaf*sumfy)/deter
                dfundy(i,k)=(sum_xcaf2*sumfy-sum_xcycaf*sumfx)/deter
                dfundsig(i,k) = (fun(i,k+1)-fun(i,k-1))/(dsigv(k-1)+dsigv(k))
            enddo
        ENDDO

        call BCpcenter(dfundx,xc,yc,sig,dsig,No_cp,nbe)
        call BCpcenter(dfundy,xc,yc,sig,dsig,No_cp,nbe)
        call BCpcenter(dfundsig,xc,yc,sig,dsig,No_cp,nbe)
      endif
!*********************************************************************!
!                                                                     !
!                  Finalizastion of the subroutine                    !
!                                                                     !
!*********************************************************************!
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: GradientP'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
