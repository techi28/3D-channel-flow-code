!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      Navier Stokes Equation Solver                  !
!                            Dec  2015                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE  hydro(ufnp,ufn,uf,ufv,                   &
                        vfnp,vfn,vf,vfv,                   &
                        wfnp,wfn,wf,wfv,                   &
                        pfnp,pfn,pf,pfv,                   &
!                       ------------------------------------
                        rhof,viscof,                       &
                        rhofv,viscofv,                     &
!                       ------------------------------------
                        xc,yc,sig,dsig,No_cp,nbe,          &
                        xv,yv,sigv,dsigv,No_vp,nbev,       &
!                       ------------------------------------
                        Hpr,h,eta,etan,                    &
                        Hprv,hv,etav,                      &
!                       ------------------------------------
                        flag_ab)
!  ------------------------------------------------------------------ !
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- ufnp   |(N_CELL,NZ)  | u Approximate solution cell-center|  !
!  | <--> ufv    |(N_VERT,NZ-1)| u Approximate solution vertex     |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- vfnp   |(N_CELL,NZ)  | v Approximate solution cell-center|  !
!  | <--> vfv    |(N_VERT,NZ-1)| v Approximate solution vertex     |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- wfnp   |(N_CELL,NZ)  | w Approximate solution cell-center|  !
!  | <--> wfv    |(N_VERT,NZ-1)| w Approximate solution vertex     |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- pfnp   |(N_CELL,NZ)  | p Approximate solution cell-center|  !
!  | <--> pfv    |(N_VERT,NZ-1)| p Approximate solution vertex     |  !
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !
!  |____________|_____________|____________________________________|  !
!  | ---> uf    |(N_CELL,NZ)  | u solution at the current RK step  |  !
!  | ---> vf    |(N_CELL,NZ)  | v solution at the current RK step  |  !
!  | ---> wf    |(N_CELL,NZ)  | w solution at the current RK step  |  !
!  | ---> pf    |(N_CELL,NZ)  | p solution at the current RK step  |  !
!  |____________|_____________|____________________________________|  !
!  | ---> ufn   |(N_CELL,NZ)  | u solution at time n               |  !
!  | ---> vfn   |(N_CELL,NZ)  | v solution at time n               |  !
!  | ---> wfn   |(N_CELL,NZ)  | w solution at time n               |  !
!  | ---> pfn   |(N_CELL,NZ)  | p solution at time n               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc |(N_CELL)     | Coordinates of the cell centers    |  !
!  | ---> sig   |(NZ)         | Sigma value at the cell centers    |  !
!  | ---> dsig  |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | ---> No_cp |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | ---> nbe   |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | ---> sigv  |(NZ-1)       | sigma of the vertex points         |  !
!  | ---> dsigv |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !
!  | ---> No_vp |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | ---> nbev  |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
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
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
#     endif
!     =============== END ================
!     ====================================
#     ifdef KeyIBM
            USE ibm
#     endif
!     -------------
#     ifdef KeyLES
            USE les
#     endif
            implicit none
!     ____________________________________
!    |                                    |
!    |      Declaration of variables      |
!    |____________________________________|

      real*8, dimension(:,:):: ufnp(N_CELL,NZ)
      real*8, dimension(:,:):: ufn(N_CELL,NZ)
      real*8, dimension(:,:):: uf(N_CELL,NZ)
      real*8, dimension(:,:):: ufv(N_VERT,NZ-1)
!     -------------------------------------
      real*8, dimension(:,:):: vfnp(N_CELL,NZ)
      real*8, dimension(:,:):: vfn(N_CELL,NZ)
      real*8, dimension(:,:):: vf(N_CELL,NZ)
      real*8, dimension(:,:):: vfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8, dimension(:,:):: wfnp(N_CELL,NZ)
      real*8, dimension(:,:):: wfn(N_CELL,NZ)
      real*8, dimension(:,:):: wf(N_CELL,NZ)
      real*8, dimension(:,:):: wfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8, dimension(:,:):: pfnp(N_CELL,NZ)
      real*8, dimension(:,:):: pfn(N_CELL,NZ)
      real*8, dimension(:,:):: pf(N_CELL,NZ)
      real*8, dimension(:,:):: pfv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: rhof(N_CELL,NZ)
      real*8,dimension(:,:) :: viscof(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:,:) :: rhofv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: viscofv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      integer :: flag_ab

!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

!     --------------------------------------
      real*8,dimension(:,:),allocatable :: ufstar,vfstar,wfstar
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: Gamu,Gamv,Gamw,Gamp
      real*8,dimension(:,:),allocatable :: dpdx,dpdy,dpds
      real*8,dimension(:,:),allocatable :: rhsu,rhsv,rhsw,rhsp
!     --------------------------------------
      real*8 :: c1,c2
      integer :: Display,jc1,jc2,jc3,iter
!     ----------------
#     ifdef KeyParallel
      real*8 :: t1,t2,tmax
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)') '>>>>> Begin subroutine: hydro'
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

      allocate(ufstar(N_CELL,NZ),vfstar(N_CELL,NZ),wfstar(N_CELL,NZ), &
               dpdx(N_CELL,NZ),dpdy(N_CELL,NZ),dpds(N_CELL,NZ),       &
               Gamu(N_CELL,NZ),Gamv(N_CELL,NZ),Gamw(N_CELL,NZ),       &
               Gamp(N_CELL,NZ),rhsp(N_CELL,NZ),                       &
               rhsu(N_CELL,NZ),rhsv(N_CELL,NZ),rhsw(N_CELL,NZ))
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |           ---- 1) UPDATE FREE SURFACE    ----          |
!     |                     °  (eta)^(new)                     |
!     |                     °  (Hpr)^(new)                     |
!     |________________________________________________________|
!      ********************************************************
#     ifdef KeyParallel
           if(rang_topo .eq. 0) then
            Display = 1
            else
            Display = 0
           endif
#     else
            Display = 1
#     endif
!       _______________________________________________
!     |                                               |
!     |  1.2)  Update the time advance coeff          |
!     |_______________________________________________|
!     -------------------------
!     USE AB2 time stepping
         if(flag_ab .eq. 1) then
            c1 = 1.55d0
            c2 = -0.55d0
         else
            c1 = 1.0d0
            c2 = 0.0d0
         endif
!       _______________________________________________
!     |                                               |
!     |  1.2)  Sigma transformation                   |
!     |_______________________________________________|
      do k=1,NZ
         do i=1,N_CELL
            sigmat(i,k) = 0.0d0
            sigmax(i,k) = 0.0d0
            sigmay(i,k) = 0.0d0
            sigmaz(i,k) = 1.0d0
         enddo
      enddo
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |       ----  2) VELOCITY: PREDICTION STEP  ----         |
!     |                    °  uf^(*)                           |
!     |                    °  vf^(*)                           |
!     |                    °  wf^(*)                           |
!     |________________________________________________________|
!      ********************************************************
!      _______________________________________________
!     |                                               |
!     |  2.1)   eddy viscosity for turbulent flow     |
!     |_______________________________________________|

      if (Iturbu .eq. 1) then
!     ------------
#     ifdef KeyLES
            call sgles(uf,vf,wf,                      &
                       xc,yc,sig,dsig,No_cp,nbe,      &
                       xv,yv,sigv,dsigv,No_vp,nbev,   &
                       Hpr,h,eta,etan,                &
                       Hprv,hv,etav)
#      endif

      endif
!     =============== END ================
!     ====================================
!      _______________________________________________
!     |                                               |
!     |  2.2)  Right-hand side values                 |
!     |_______________________________________________|

!     ________________________________________________
!     Dynamic case: p_hydro + gravity + p_dynamic(old)
!        -----------------
!        critical to loop over all cells, otherwise need to communicate
         do k=1,NZ
            do i=1,N_CELL
               rhof(i,k) = 1.0d0
#        ifdef KeyLES
!        --------------------------------------
!         LES turbulence model
                  Gamu(i,k) = eddy(i,k)+1./Re
                  Gamv(i,k) = eddy(i,k)+1./Re
                  Gamw(i,k) = eddy(i,k)+1./Re
#        else
!        -------------------------------------------
!          No turbulence Model
                  Gamu(i,k) = 1./Re
                  Gamv(i,k) = 1./Re
                  Gamw(i,k) = 1./Re
#        endif
!        -------------------------------------------
                  rhsu(i,k) =  bdx
                  rhsv(i,k) =  bdy
                  rhsw(i,k) =  bdz
            enddo ! end i loop
         enddo ! end k loop
!      _______________________________________________
!     |                                               |
!     |  2.5) Update velocity components: uf,vf,wf^(*)|
!     |_______________________________________________|
      do k=1,NZ
        do i=1,N_CELL
          ufstar(i,k) = uf(i,k)
          vfstar(i,k) = vf(i,k)
          wfstar(i,k) = wf(i,k)
        enddo
      enddo
!     _________________________________________________
!     Update advection-diffusion equation
        call AdvDiffNS(rhsu,rhsv,rhsw,                 &
                     uf,vf,wf,                         &
                     Gamu,Gamv,Gamw,                   &
                     xc,yc,sig,dsig,No_cp,nbe,         &
                     xv,yv,sigv,dsigv,No_vp,nbev)
!    ---------------------------------------------------------------
!    Get intermediate velocity
      do k=2, NZ-1
        do i=1, N_CELL0
          ufstar(i,k) = ufstar(i,k) + dt*(c1*rhsu(i,k)+c2*rhsuf(i,k))
          vfstar(i,k) = vfstar(i,k) + dt*(c1*rhsv(i,k)+c2*rhsvf(i,k))
          rhsuf(i,k) = rhsu(i,k)
          rhsvf(i,k) = rhsv(i,k)
          IF(I3D .EQ. 1) THEN
          wfstar(i,k) = wfstar(i,k) + dt*(c1*rhsw(i,k)+c2*rhswf(i,k))
          rhswf(i,k) = rhsw(i,k)
          ENDIF
        enddo
      enddo
!      --------------------------------------------------------
!      Add Immersed Boundary Feedback Force
#      ifdef KeyIBM
       call ibp_interpolate(ufstar,vfstar,wfstar,pf,  &
                            xc,yc,sig,dsig,No_cp,nbe, &
                            xv,yv,sigv,dsigv,No_vp,nbev)
#     endif
!      --------------------------------------------------------
       call BCglobalVC(ufstar,vfstar,wfstar,              &
                       ufv,vfv,wfv,                       &
                       xc,yc,sig,dsig,No_cp,nbe,          &
                       xv,yv,sigv,dsigv,No_vp,nbev,3)
!      _______________________________________________
!     |                                               |
!     |  4.1)  Coefficients: diffusive coeff = H*Gamp |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL
            Gamp(i,k) = 1.0d0
         enddo
      enddo
!      _______________________________________________
!     |                                               |
!     |  4.2)  Pressure by projection method:         |
!     |_______________________________________________|

           call source3D(rhsp,ufstar,vfstar,wfstar,      &
                         xc,yc,sig,dsig,No_cp,nbe,       &
                         xv,yv,sigv,dsigv,No_vp,nbev)

            do k=2,NZ-1
               do i=1,N_CELL0
#                ifdef KeyOldPressure
                 rhsp(i,k) = (1.0d0/c1)*rhsp(i,k)/dt
#                else
                 rhsp(i,k) = rhsp(i,k)/dt
#                endif
               enddo
            enddo
!      _______________________________________________
!     |                                               |
!     |  4.3)               PDE SOLVER                |
!     |_______________________________________________|

!        __________________________________________________
!        Velocity correction
          call gradientP(pf,pfv,                     &
                         dpdx,dpdy,dpds,             &
                         xc,yc,sig,dsig,No_cp,nbe,   &
                         xv,yv,sigv,dsigv,No_vp,nbev)
!        DO iter = 1,10
!        __________________________________________________
!        Initial guess
         do k=1,NZ
            do i=1,N_CELL
                pfnp(i,k) = pf(i,k)
            enddo
         enddo
!        --------------------------------------------------
!        ADD cross diffusion contribution to rhsp
#        ifdef KeyPPERHS
          call update_rhsp(rhsp,pf,pfv,dpdx,dpdy,dpds,     &
                           xc,yc,sig,dsig,No_cp,nbe,       &
                           xv,yv,sigv,dsigv,No_vp,nbev)
#        endif
!        ---------------------------------
#        ifdef KeyParallel
            call MPI_Barrier(comm3D,code)
            t1 = MPI_Wtime()
#        endif
           call bicgstab(pfnp,pfv,rhsp,                 &
                         xc,yc,sig,dsig,No_cp,nbe,      &
                         xv,yv,sigv,dsigv,No_vp,nbev,   &
                         Hpr,h,etan,                    &
                         Hprv,hv,etav)

#       ifdef KeyParallel
         t2 = MPI_Wtime()
         tmax = t2 -t1
         call MPI_Barrier(comm3D,code)
         call MAX_parallel(tmax, tPoisson)
         if(rang_topo .eq. 0) print*, 'Time for Poisson =',tPoisson
#       endif
!        __________________________________________________
!        Velocity correction
          call gradientP(pfnp,pfv,                   &
                         dpdx,dpdy,dpds,             &
                         xc,yc,sig,dsig,No_cp,nbe,   &
                         xv,yv,sigv,dsigv,No_vp,nbev)
!        -----------------
!         update cell centre velocity
#         ifdef KeyOldPressure
          do k=2,NZ-1
            do i=1,N_CELL0
               ufnp(i,k) = ufstar(i,k)-c1*dt*dpdx(i,k)
               vfnp(i,k) = vfstar(i,k)-c1*dt*dpdy(i,k)
               rhsuf(i,k) = rhsuf(i,k) - dpdx(i,k)
               rhsvf(i,k) = rhsvf(i,k) - dpdy(i,k)
               IF(I3D .EQ. 1) THEN
               wfnp(i,k) = wfstar(i,k)-c1*dt*dpds(i,k)
               rhswf(i,k) = rhswf(i,k) - dpds(i,k)
               ENDIF
            enddo
          enddo
#         else
          do k=2,NZ-1
            do i=1,N_CELL0
               ufnp(i,k) = ufstar(i,k)-dt*dpdx(i,k)
               vfnp(i,k) = vfstar(i,k)-dt*dpdy(i,k)
               IF(I3D .EQ. 1) THEN
               wfnp(i,k) = wfstar(i,k)-dt*dpds(i,k)
               ENDIF
            enddo
          enddo
#         endif
!    ------------------------------------------
!    Set up the global boundary condition
        call BCglobalVC(ufnp,vfnp,wfnp,              &
                        ufv,vfv,wfv,                 &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev,3)
!    ---------------
!    update new velocity
!      call interpolateU(U1FACE,U2FACE,U3FACE,UTFACE,UBFACE, &
!                        ufnp,vfnp,wfnp,                     &
!                        xc,yc,sig,dsig,No_cp,nbe,           &
!                        xv,yv,sigv,dsigv,No_vp,nbev)

        call newu3D(pfnp,pfv,flag_ab,             &
                    dpdx,dpdy,dpds,               &
                    xc,yc,sig,dsig,No_cp,nbe,     &
                    xv,yv,sigv,dsigv,No_vp,nbev)
!        ---------------
!        reset the time stepping flag
            if(FinRK .eq. 1) then
                if(time_ab2 .eq. 1) then
                    if(flag_ab .ne. 1) then
                        flag_ab = 1
                        if(Display .eq. 1) then
                            print*, 'Switch from Euler to AB2'
                        endif
                    endif
                endif
            endif
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

      deallocate(ufstar,vfstar,wfstar,  &
                 dpdx,dpdy,dpds,        &
                 rhsu,rhsv,rhsw,rhsp,   &
                 Gamu,Gamv,Gamw,Gamp)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)') '<<<<< End   subroutine: hydro'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!          END OF FREE SURFACE TEST WITH NAVIER-STOKES EQUATION       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
