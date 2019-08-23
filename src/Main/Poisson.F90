!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!      SOLUTION OF THE POISSON EQUATION FOR FREE SURFACE PROBLEMS     !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Poisson(phi,phiv,                      &
                         rhs,Gamx,Gamy,Gamz,            &
                         xc,yc,sig,dsig,No_cp,nbe,      &
                         xv,yv,sigv,dsigv,No_vp,nbev,   &
                         Hpr,h,etan,                    &
                         Hprv,hv,etav)

!---------------------------------------------------------------------!
!                                                                     !
!     This subroutine solves the poisson equation of the pressure     !
!     for free surface problems.                                      !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |   Name     |   Size     | Description                         |  !
!  |____________|____________|_____________________________________|  !
!  | <-> phi    |(N_CELL,NZ) | Cell-center solution at t(n+1)      |  !
!  | <-> phiv   |(N_VERT,NZ) | Cell-vertex solution at t(n+1)      |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !
!  |____________|____________|_____________________________________|  !
!  | --> Gam    |(N_CELL,NZ) | Diffusive coefficients              |  !
!  | --> rhs    |(N_CELL,NZ) | right-hand side                     |  !
!  |____________|____________|_____________________________________|  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbev   |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!  | --> Hpr    |(N_CELL)    | Total depth = h + etan (cell-center)|  !
!  | --> h      |(N_CELL)    | Still depth            (cell-center)|  !
!  | --> etan   |(N_CELL)    | Free surface           (cell-center)|  !
!  |____________|____________|_____________________________________|  !
!  | --> Hprv   |(N_VERT)    | Total depth = h + etan (vertices)   |  !
!  | --> hv     |(N_VERT)    | Still depth            (vertices)   |  !
!  | --> etav   |(N_VERT)    | Free surface           (vertices)   |  !
!  |____________|____________|_____________________________________|  !
!  | --> rhof   |(N_CELL,NZ) | Density at the cell-centers         |  !
!  | --> dwfdtB |(N_CELL)    | d(w)/dt at the bottom boundary      |  !
!  |____________|____________|_____________________________________|  !
!  | --> rhofv  |(N_VERT,NZ) | Density at the vertices             |  !
!  | --> dwfvdtB|(N_VERT)    | d(w)/dt at the bottom boundary      |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |   - FS_BCpc                    ( FS_BCpressure.F90   )        |  !
!  |   - FS_BCpv                    ( FS_BCpressure.F90   )        |  !
!  |   - diffusion3D                ( diffusion3D.F90     )        |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input/Output variables                                       !
!---------------------------------------------------------------------!

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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     ----------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
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
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)

!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: Am0,Am1,Am2,Am3,AmT,AmB
      real*8,dimension(:,:),allocatable :: Bmv1T,Bmv2T,Bmv3T
      real*8,dimension(:,:),allocatable :: Bmv1B,Bmv2B,Bmv3B
!     --------------------------------------
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------

      real*8  :: errorsys,residu,som,Vol
      real*8  :: errorNeum,errorsysOld,bnrm
#     ifdef KeyParallel
      real*8  :: Perrorsys
#     endif
      integer :: iter,DisplayThis
!     --------------------------------------
      real*8,dimension(:,:) :: phiNew(N_CELL,NZ)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: sgm_Poisson'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(Am0(N_CELL0,NZ),Am1(N_CELL0,NZ),Am2(N_CELL0,NZ),       &
               Am3(N_CELL0,NZ),AmT(N_CELL0,NZ),AmB(N_CELL0,NZ),       &
               Bmv1T(N_CELL0,NZ),Bmv2T(N_CELL0,NZ),Bmv3T(N_CELL0,NZ), &
               Bmv1B(N_CELL0,NZ),Bmv2B(N_CELL0,NZ),Bmv3B(N_CELL0,NZ))

!      ________________________________________________________
!     |                                                        |
!     |                    Correct rhs values                  |
!     |________________________________________________________|

      do i=1,N_CELL
         Vol = AreaCell(i)*dsigv(1)
         rhs(i,1) = Vol*rhs(i,1)
         do k=2,NZ
            Vol = AreaCell(i)*dsigv(k-1)
            rhs(i,k) = Vol*rhs(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                     Initial guess                      |
!     |________________________________________________________|

!     ----------------------------------------------
!     Boundary Condition of the cell-center points
      call BCglobalP(phi,phiv,                      &
                     xc,yc,sig,dsig,No_cp,nbe,      &
                     xv,yv,sigv,dsigv,No_vp,nbev,3)
!     ________________________________________________________
!    |                                                        |
!    |            Matrix Am & Bm of the diffusion term        |
!    |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0
            Am0(i,k)  = 0.0d0
            Am1(i,k)  = 0.0d0
            Am2(i,k)  = 0.0d0
            Am3(i,k)  = 0.0d0
            AmT(i,k)  = 0.0d0
            AmB(i,k)  = 0.0d0
            Bmv1T(i,k)= 0.0d0
            Bmv2T(i,k)= 0.0d0
            Bmv3T(i,k)= 0.0d0
            Bmv1B(i,k)= 0.0d0
            Bmv2B(i,k)= 0.0d0
            Bmv3B(i,k)= 0.0d0
         enddo
      enddo

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                       Bmv1T,Bmv2T,Bmv3T,                   &
                       Bmv1B,Bmv2B,Bmv3B,                   &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!                 Solution of the system (S.0.R.) 3D                  !
!                                                                     !
!*********************************************************************!

!     ________________________________________________________
!    |                                                        |
!    |                    Update coefficients                 |
!    |________________________________________________________|
      do k=2,NZ-1
         do i=1,N_CELL0
            if (Am0(i,k)==0) then
               print*,'Error!!!: Division over zero Am0 i=',i,' k=',k
               stop
            endif
            Am1(i,k)  = Am1(i,k)/Am0(i,k)
            Am2(i,k)  = Am2(i,k)/Am0(i,k)
            Am3(i,k)  = Am3(i,k)/Am0(i,k)
            AmT(i,k)  = AmT(i,k)/Am0(i,k)
            AmB(i,k)  = AmB(i,k)/Am0(i,k)
            Bmv1T(i,k)= Bmv1T(i,k)/Am0(i,k)
            Bmv2T(i,k)= Bmv2T(i,k)/Am0(i,k)
            Bmv3T(i,k)= Bmv3T(i,k)/Am0(i,k)
            Bmv1B(i,k)= Bmv1B(i,k)/Am0(i,k)
            Bmv2B(i,k)= Bmv2B(i,k)/Am0(i,k)
            Bmv3B(i,k)= Bmv3B(i,k)/Am0(i,k)
            rhs(i,k)  = rhs(i,k)/Am0(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|
      do k=1,NZ
         do i=1,N_CELL
	    phiNew(i,k) = 0.0d0
         enddo
      enddo

      iter=0
      bnrm = 0.0d0
      do k=2,NZ-1
       do i=0,N_CELL0
         bnrm = bnrm + rhs(i,k)*rhs(i,k)
        enddo
      enddo
#     ifdef KeyParallel
       call SUM_parallel(bnrm,errorsys)
       bnrm = dsqrt(errorsys)
#     else
       bnrm = dsqrt(bnrm)
#     endif
111   continue
      iter=iter+1
!      ________________________________________________________
!     |                                                        |
!     |           New right-hand side: rhs - Bm(vertex)        |
!     |________________________________________________________|

#     ifdef Key_PoissonJacobi
      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
!           _________________________________________________
!           Bm(vertex)
            residu = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                               +Bmv2T(i,k)*phiv(jv2,k)   &
                               +Bmv3T(i,k)*phiv(jv3,k)   &
                               +Bmv1B(i,k)*phiv(jv1,k-1) &
                               +Bmv2B(i,k)*phiv(jv2,k-1) &
                               +Bmv3B(i,k)*phiv(jv3,k-1) &
!           _________________________________________________
!           Am(center)
                               +Am1(i,k)*phi(jc1,k)      &
                               +Am2(i,k)*phi(jc2,k)      &
                               +Am3(i,k)*phi(jc3,k)      &
                               +AmT(i,k)*phi(i,k+1)      &
                               +AmB(i,k)*phi(i,k-1)      &
                               +phi(i,k))
	    errorsys = errorsys + abs(residu)*abs(residu)
            if((k .eq. 2) .and. (i .eq.1)) then
              phiNew(i,k) = 0.
            else
	      phiNew(i,k) = phi(i,k) + relaxSOR*residu
            endif
         enddo
      enddo

      do k=2,NZ-1
         do i=1,N_CELL0
	    phi(i,k) = phiNew(i,k)
         enddo
      enddo
#     endif
!      ________________________________________________________
!     |                                                        |
!     |           New right-hand side: rhs - Bm(vertex)        |
!     |________________________________________________________|
#     ifdef Key_PoissonJSOR
        errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
!           _________________________________________________
!           Bm(vertex)
            residu = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                               +Bmv2T(i,k)*phiv(jv2,k)   &
                               +Bmv3T(i,k)*phiv(jv3,k)   &
                               +Bmv1B(i,k)*phiv(jv1,k-1) &
                               +Bmv2B(i,k)*phiv(jv2,k-1) &
                               +Bmv3B(i,k)*phiv(jv3,k-1) &
!           _________________________________________________
!           Am(center)
                               +Am1(i,k)*phi(jc1,k)      &
                               +Am2(i,k)*phi(jc2,k)      &
                               +Am3(i,k)*phi(jc3,k)      &
                               +AmT(i,k)*phi(i,k+1)      &
                               +AmB(i,k)*phi(i,k-1)      &
                               +phi(i,k))
	    errorsys = errorsys + abs(residu)*abs(residu)
	    phi(i,k) = phi(i,k) + relaxSOR*residu
         enddo
      enddo
#     endif
!     ________________________________________________________
!     Boundary conditions cell-centers
      call BCglobalP(phi,phiv,                      &
                     xc,yc,sig,dsig,No_cp,nbe,      &
                     xv,yv,sigv,dsigv,No_vp,nbev,3)
!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |
!     |________________________________________________________|
#     ifdef KeyParallel
      call SUM_parallel(errorsys,Perrorsys)
      if(rang_topo .eq. 0) then
        DisplayThis = 1
      else
        DisplayThis = 0
      endif
      errorsys = dsqrt(Perrorsys)      
#     else
      DisplayThis = 1
      errorsys = dsqrt(errorsys)
#     endif
!     ------------------------------------------------------
      if(iter .eq. 1) then
        bnrm = errorsys
        errorsys = 1.
      else
       errorsys = errorsys/bnrm
!     ------------------------------------------------------
      endif
      if (errorsys.lt.eps) then
        if(DisplayThis .eq. 1) then
         write(*,*) ' '
#     ifdef KeyFS_PoissonJSOR
         write(*,7) 'Solution J.S.0.R. 3D  : iters =',iter,&
                    ', error =',errorsys
#     else
         write(*,7) 'Solution Jacobi 3D  : iters =',iter,&
                    ', error =',errorsys
#     endif
         write(*,*) ' '
         endif
      elseif (errorsys.gt.1.0d12) then
        if(DisplayThis .eq. 1) then
         write(*,*) ' '
         write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                    ', error =',errorsys
         write(*,*) ' '
         endif
         stop
      elseif(iter.gt.MaxIters) then
        if(DisplayThis .eq. 1) then
         write(*,*) ' '
         write(*,7) 'Non-convergence: iters =',iter,&
                    ', error =',errorsys
         write(*,*) ' '
         endif
      else
         if ((0.eq.mod(iter,1000)) .or. (iter .eq. 1)) then
            if(DisplayThis .eq. 1) print*, 'iter=',iter,'Error=',errorsys
         endif
         goto 111
      endif

119   continue

      7 format(t10,a24,i5,a9,e10.3)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         De-allocate                    |
!     |________________________________________________________|

      deallocate(Am0,Am1,Am2,       &
                 Am3,AmT,AmB,       &
                 Bmv1T,Bmv2T,Bmv3T, &
                 Bmv1B,Bmv2B,Bmv3B)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- End   subroutine: Poisson'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF Poisson for Free surface                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
