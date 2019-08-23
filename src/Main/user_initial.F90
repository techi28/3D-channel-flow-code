!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               INITIALIZATION OF FLOW AND PRESSURE FIELD             !
!                             December 2015                           !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE user_initial(ufn,vfn,wfn,pfn,            &
                              ufv,vfv,wfv,pfv,             &
                              xc,yc,sig,dsig,No_cp,nbe,    &
                              xv,yv,sigv,dsigv,No_vp,nbev, &
                              Hpr,h,etan,                  &
                              Hprv,hv,etav,                &
                              xct,yct,zct,                 &
                              xvt,yvt,zvt)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the initial condition of the time test   !
!    problems with free surface.                                      !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |           Name       |    Size     | Description              |  !
!  |______________________|_____________|__________________________|  !
!  | <--- ufn,vfn,wfn,pfn |(N_CELL,NZ)  | Initial solution center  |  !
!  | <--- ufv,vfv,wfv,pfv |(N_VERT,NZ-1)| Initial solution vertex  |  !
!  |______________________|_____________|__________________________|  !
!  | <--- Hpr,etan        |(N_CELL)     | Initial solution center  |  !
!  | <--- Hprv,etavv      |(N_VERT)     | Initial solution vertex  |  !
!  |______________________|_____________|__________________________|  !
!  | <--- xct,yct,zct     |(N_CELL,NZ)  | Coordinates center       |  !
!  | <--- xvt,yvt,zvt     |(N_VERT,NZ-1)| Coordinates vertex       |  !
!  |______________________|_____________|__________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc |(N_CELL)    | Coordinates of the cell centers     |  !
!  | ---> sig   |(NZ)        | Sigma value at the cell centers     |  !
!  | ---> dsig  |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | ---> No_cp |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | ---> nbe   |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xv,yv |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | ---> sigv  |(NZ-1)      | sigma of the vertex points          |  !
!  | ---> dsigv |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !
!  | ---> No_vp |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | ---> nbev  |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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

      real*8,dimension(:,:) :: ufn(N_CELL,NZ)
      real*8,dimension(:,:) :: vfn(N_CELL,NZ)
      real*8,dimension(:,:) :: wfn(N_CELL,NZ)
      real*8,dimension(:,:) :: pfn(N_CELL,NZ)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      real*8,dimension(:,:) :: xct(N_CELL,NZ)
      real*8,dimension(:,:) :: yct(N_CELL,NZ)
      real*8,dimension(:,:) :: zct(N_CELL,NZ)
      real*8,dimension(:,:) :: xvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: yvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: zvt(N_VERT,NZ-1)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z
      integer :: IDISPLAY
      real*8 :: funu,funv,funw,funp
      real*8 :: funu2,funv2,funw2,funp2
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: USER_INITIAL'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#     ifdef KeyParallel
       if(rang_topo .eq. 0) then
         IDISPLAY = 1
       else
         IDISPLAY = 0
       endif
#     else
         IDISPLAY = 1
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                       Free surface                     |
!     |________________________________________________________|

!     -----------------------------------
!     Cell-center
      do i=1,N_CELL
         x = xc(i)
         y = yc(i)
         h(i) = 0.0d0
         Hpr(i)  = 1.0d0
      enddo
!     -----------------------------------
!     Vertex
      do nv=1,N_VERT
         x = xv(nv)
         y = yv(nv)
         hv(nv) = 0.0d0
         Hprv(nv) = 1.0d0
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                 Velocity and pressure                  |
!     |________________________________________________________|

#    ifdef KeyTESTChannel
!     -----------------------------------
!     Cell-center
      do i=1,N_CELL
         do k=1,NZ
            x = xc(i)
            y = yc(i)
            z = sig(k)
            ufn(i,k) = funu(x,y,z)
            vfn(i,k) = funv(x,y,z)
            wfn(i,k) = funw(x,y,z)
            pfn(i,k) = funp(x,y,z)
         enddo
      enddo
!     -----------------------------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            x = xv(nv)
            y = yv(nv)
            z = sigv(k)
    	    ufv(nv,k) = funu(x,y,z)
            vfv(nv,k) = funv(x,y,z)
            wfv(nv,k) = funw(x,y,z)
            pfv(nv,k) = funp(x,y,z)
         enddo
      enddo
#    endif

!      ________________________________________________________
!     |                                                        |
!     |                 Coordinates (xt,yt,zt)                 |
!     |________________________________________________________|

!     -----------------------------------
!     Cell-center
      do i=1,N_CELL
         do k=1,NZ
            xct(i,k) = xc(i)
            yct(i,k) = yc(i)
            zct(i,k) = sig(k)*Hpr(i)-h(i)
         enddo
      enddo
!     -----------------------------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            xvt(nv,k) = xv(nv)
            yvt(nv,k) = yv(nv)
            zvt(nv,k) = sigv(k)*Hprv(nv)-hv(nv)
         enddo
      enddo

!    ------------------------------------------
!    Set up the global boundary condition
      call BCglobalVC(ufn,vfn,wfn,              &
                      ufv,vfv,wfv,              &
                      xc,yc,sig,dsig,No_cp,nbe, &
                      xv,yv,sigv,dsigv,No_vp,nbev,3)
!    ------------------------------------------
!    Set up the global pressure condition
      call BCglobalP(pfn,pfv,                  &
                    xc,yc,sig,dsig,No_cp,nbe, &
                    xv,yv,sigv,dsigv,No_vp,nbev,3)
      call interpolateU(U1FACE,U2FACE,U3FACE,UTFACE,UBFACE, &
                        ufn,vfn,wfn,                        &
                        xc,yc,sig,dsig,No_cp,nbe,           &
                        xv,yv,sigv,dsigv,No_vp,nbev)
!    ---------------------------------------------
!    allocate the history RHS value
     allocate(rhsuf(N_CELL,NZ),rhsvf(N_CELL,NZ),rhswf(N_CELL,NZ))
!       ----------
!       initial to zero
        do k=1,NZ
          do i=1,N_CELL
            rhsuf(i,k) = 0.0d0
            rhsvf(i,k) = 0.0d0
            rhswf(i,k) = 0.0d0
          enddo
        enddo
!    -------------------------------------------------------
!    Display the options
      if(IDISPLAY .EQ. 1) then
        print*, '      ==========================='

          if(I3D .eq. 1) then
               print*, '            DIM  OPT : 3D'
           else
               print*, '            DIM  OPT : 2D'
          endif

#         ifdef KeyImplicit
               print*, '            TIME OPT : Semi-implicit'
#         else
               print*, '            TIME OPT : Explicit'
#         endif

#         ifdef KeyPPECenter
               print*, '            PPE  OPT : Center Matrix'
#         endif

        print*, '      ==========================='

      endif

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
      call glob_cfl(ufn,vfn,wfn,pfn,               &
                    xc,yc,sig,dsig,No_cp,nbe)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: USER_INITIAL'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                     END OF User Defined Initial                     !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
