!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!     SOLUTION OF THE ADVECTION-DIFFUSION PROBLEM WITH FREE SURFACE   !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AdvDiffNS(rhsu,rhsv,rhsw,                   &
                           funu,funv,funw,                   &
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
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real*8, dimension(:,:), allocatable :: funuv,funvv,funwv
      integer :: DoThis

#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'Begin subroutine: AdvDiffNS'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      allocate(funuv(N_VERT,NZ-1),funvv(N_VERT,NZ-1),funwv(N_VERT,NZ-1))

      DO i=1,N_VERT
        do k=1,NZ-1
           funuv(i,k) = 0.
           funvv(i,k) = 0.
           funwv(i,k) = 0.
        enddo
      ENDDO
!     ---------------------------------------------------------------
!     get 2D vertex value for face value interpolation
      call BCglobalVC(funu,funv,funw,             &
                      funuv,funvv,funwv,          &
                      xc,yc,sig,dsig,No_cp,nbe,   &
                      xv,yv,sigv,dsigv,No_vp,nbev,2)
!*********************************************************************!
!                                                                     !
!                       Advection-Diffusion  3D                       !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                 Diffusion contribution                 |
!     |________________________________________________________|
!    -----------------------------------------------------------
!     switch for different combinations of adv and diff
!     0 (no adv and diff); 1 (only diff); 2 (only adv); 3 (all)
      DoThis = 3

!     ------------------------------------------------------------
      if ((DoThis .eq. 1) .or. (DoThis .eq. 3))then


      call diffusion3DNEW(rhsu,rhsv,rhsw,                   &
                          funu,funv,funw,                   &
                          funuv,funvv,funwv,                &
                          Gamx,Gamy,Gamz,                   &
                          xc,yc,sig,dsig,No_cp,nbe,         &
                          xv,yv,sigv,dsigv,No_vp,nbev,1)
      endif

!      ________________________________________________________
!     |                                                        |
!     |                 Advection contribution                 |
!     |________________________________________________________|
      if((DoThis .eq. 2) .or. (DoThis .eq. 3)) then

       call advection3DNEW(rhsu,rhsv,rhsw,                   &
                           funu,funv,funw,                   &
                           funuv,funvv,funwv,                &
                           xc,yc,sig,dsig,No_cp,nbe,         &
                           xv,yv,sigv,dsigv,No_vp,nbev,1)

      endif

      deallocate(funuv,funvv,funwv)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), ' End   subroutine AdvDiffNS'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF Advection-Diffusion 3D                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
