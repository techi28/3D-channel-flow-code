!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            WATER ELELEVATION BY THE FREE SURFACE EQUATION           !
!                              Nov 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE WaterLevel(Hpr_new,Hpr,Hprv,              &
                            eta_new,eta,etav,              &
!                           -------------------------------
                            u_new,v_new,                   &
                            ufn,vfn,                       &
!                           -------------------------------  
                            h,hv,                          &
!                           -------------------------------
                            xc,yc,sig,dsig,No_cp,nbe,      &
                            xv,yv,sigv,dsigv,No_vp,nbev,   &
!                           -------------------------------
                            gamma,                         &
!                           -------------------------------
                            No_sp)

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the free surface H = h + eta using the      !
!    equation:                                                        !
!                dH/dt +  d(H*qu)/dx + d(H*qv)/dy = 0                 ! 
!    where qu and qv are the averate velocity of u and v in the       !
!    sigma direction respectively.                                    !
!    We consider two options:                                         !
!    -------------------------------------------------------------    !
!    ChooseFreeMethod = 1:                                            !
!    The function H in the advection part is consider at the previous !
!    value and the derivatives are calculated by LSM                  !
!                   dHnew/dt + ADV(Hold,qu,qv) = 0                    !
!                                                                     !
!    and solve everything in a explicit fashion.                      !
!    -------------------------------------------------------------    !
!    ChooseFreeMethod = 2:                                            !
!    We solve the advection equation:                                 !
!                   dHnew/dt + ADV(Hnew,qu,qv) = 0                    !
!                                                                     !
!    We have four methods (chosen at file cppdefs.h) to solve the     !
!    problem :                                                        !
!    1) Explicit    : Only the diagonal elements of ADV are           !
!                     calculated at time(n+1).                        !
!    2) Implicit:    In this case the cell-center and neighbors       !
!                     are calculated at time(n+1), the gradient       !
!                     GF is calculated at (n). We use SOR to solve    !
!                     the linear system.                              !
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !
!  |_____________|_____________|___________________________________|  !
!  | <-- Hpr_new |  (N_CELL)   | Water depth solution cell-center  |  !
!  | <-- Hprv    |  (N_VERT)   | Water depth solution vertex       |  !
!  | <-- eta_new |  (N_CELL)   | Water elevation cell-center       |  !
!  | <-- etav    |  (N_VERT)   | Water elevation vertex            |  !
!  |_____________|_____________|___________________________________|  !
!  | <-- dHprdt  |  (N_CELL)   | d(Hpr)/dt for cell-centered values|  !
!  | <-- dHprdtv |  (N_VERT)   | d(Hpr)/dt for vertex values       |  !
!  | <-- detadt  |  (N_CELL)   | d(eta)/dt at cell-centered        |  !
!  | <-- detadtv |  (N_VERT)   | d(eta)/dt for vertex values       |  !
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !
!  |____________|_____________|____________________________________|  !
!  | --> h      |(N_CELL)     | Fixed water depth cell-center      |  !
!  | --> hv     |(N_VERT)     | Fixed water depth vertex           |  !
!  |____________|_____________|____________________________________|  !
!  | --> xc,yc  |(N_CELL)     | Coordinates of the cell centers    |  !
!  | --> sig    |(NZ)         | Sigma value at the cell centers    |  !
!  | --> dsig   |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | --> No_cp  |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | --> nbe    |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | --> xv,yv  |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | --> sigv   |(NZ-1)       | sigma of the vertex points         |  !
!  | --> dsigv  |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !
!  | --> No_vp  |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | --> nbev   |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!  SUBROUTINES:                                                       !
!         *  FS_ Advection2D       (Located below this subroutine)    !
!         *  Waterlevel_BC                                            !
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

      real*8, dimension(:)  :: Hpr_new(N_CELL)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
!     --------------------------------------
      real*8, dimension(:)  :: eta_new(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------
      real*8,dimension(:)   :: u_new(N_CELL,NZ)
      real*8,dimension(:)   :: v_new(N_CELL,NZ)
      real*8,dimension(:)   :: ufn(N_CELL,NZ)
      real*8,dimension(:)   :: vfn(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)  :: h(N_CELL)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     --------------------------------------
      real*8, dimension(:)  :: hv(N_VERT)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      real*8 :: gamma
!     --------------------------------------
      integer,dimension(:)  :: No_sp(N_SPmax)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8, dimension(:)  :: dHprdt_new(N_CELL)
      real*8, dimension(:)  :: dHprdx_new(N_CELL)
      real*8, dimension(:)  :: dHprdy_new(N_CELL)
      real*8, dimension(:)  :: dHprdtv(N_VERT)
      real*8, dimension(:)  :: dHprdxv(N_VERT)
      real*8, dimension(:)  :: dHprdyv(N_VERT)
      real*8, dimension(:)  :: detadt_new(N_CELL)
      real*8, dimension(:)  :: detadx_new(N_CELL)
      real*8, dimension(:)  :: detady_new(N_CELL)
      real*8, dimension(:)  :: detadtv(N_VERT)
      real*8, dimension(:)  :: detadxv(N_VERT)
      real*8, dimension(:)  :: detadyv(N_VERT)
!     --------------------------------------
      real*8,dimension(:)   :: Hprn(N_CELL)
      real*8,dimension(:)   :: etan(N_CELL)
!     --------------------------------------
      real*8,dimension(:)   :: Hnew(N_CELL)
      real*8,dimension(:)   :: Hvnew(N_VERT)
      real*8,dimension(:)   :: Hold(N_CELL)
      real*8,dimension(:)   :: Hvold(N_VERT)
!     --------------------------------------
      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0)
      real*8,dimension(:) :: AmG(N_CELL0)
      real*8,dimension(:) ::  bm(N_CELL0)
      real*8,dimension(:) ::  qu(N_CELL)
      real*8,dimension(:) ::  qv(N_CELL)
      real*8,dimension(:) :: rhs(N_CELL)
      real*8,dimension(:) :: Hprv_old(N_VERT)
      real*8,dimension(:) :: etav_old(N_VERT)
!     --------------------------------------
      real*8,dimension(:) :: Am0_new(N_CELL0),Am0n(N_CELL0)
      real*8,dimension(:) :: Am1_new(N_CELL0),Am1n(N_CELL0)
      real*8,dimension(:) :: Am2_new(N_CELL0),Am2n(N_CELL0)
      real*8,dimension(:) :: Am3_new(N_CELL0),Am3n(N_CELL0)
      real*8,dimension(:) :: AmG_new(N_CELL0),AmGn(N_CELL0)
      real*8,dimension(:) :: qu_new(N_CELL)
      real*8,dimension(:) :: qv_new(N_CELL)
      real*8,dimension(:) :: qun(N_CELL)
      real*8,dimension(:) :: qvn(N_CELL)
!     --------------------------------------
      real*8 :: dtoVol,const,Vol,som,somu,somv,c1,ct
      integer:: it,jc1,jc2,jc3,jc,nv1
!     --------------------------------------
      real*8 :: sumqux,sumquy,sumqvx,sumqvy,deter
      real*8 :: dqudx,dqvdy,som1,som2,som3,som4
      real*8 :: errorsys,residu,SUMerrorsys,FS_funeta
!     --------------------------------------
      integer :: Display
!     --------------------------------------
      real*8 :: xvR1,yvR1
      real*8 :: DR0,DR,etaMax
      integer:: n,pro
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      real*8, dimension(Nprocs) :: DR0V
      integer,dimension(Nprocs) :: nv1V
      integer,dimension(Nprocs) :: proV
#     endif
!     =============== END ================    
!     ====================================

!     ____________________________________
!    |                                    |
!    |           Choose method            |
!    |____________________________________|

      real*8, parameter :: WLrelaxSOR = 1.0d0
      integer,parameter :: ChooseFreeMethod = 2
      ! =1 LSM
      ! =2 FVM

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: FS_WaterLevel'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ====================================
!     =====    DISPLAY ITERATIONS  =======
      Display = 1
#     ifdef KeyParallel
         if (rang_topo.ne.0) Display = 0 
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!                       Components of the problem                     !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                       Initial values                   |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication2D(Hpr)
        call communication2D(eta)
        call communication3Dtype2(ufn)
        call communication3Dtype2(vfn)
        call communication3Dtype2(u_new)
        call communication3Dtype2(v_new)
#     endif
!     ====================================
!     ====================================

      Hprn = Hpr
      etan = eta
      do nv=1,N_VERT
         Hprv_old(nv) = Hprv(nv)
         etav_old(nv) = etav(nv)  
      enddo
!      ________________________________________________________
!     |                                                        |
!     |        Average velocity values (integrals) & rhs       |
!     |________________________________________________________|

      do i=1,N_CELL
!        ---------------------
         somu = 0.0d0
         somv = 0.0d0
         do k=1,NZ-1
            somu = somu + u_new(i,k)*dsig(k)
            somv = somv + v_new(i,k)*dsig(k)
         enddo
         qu_new(i) = somu 
         qv_new(i) = somv  
!        ---------------------
         somu = 0.0d0
         somv = 0.0d0
         do k=1,NZ-1
            somu = somu + ufn(i,k)*dsig(k)
            somv = somv + vfn(i,k)*dsig(k)
         enddo
         qun(i) = somu 
         qvn(i) = somv  
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication2D(qun)
        call communication2D(qvn)
        call communication2D(qu_new)
        call communication2D(qv_new)
#     endif
!     ====================================
!     ====================================

!*********************************************************************!
!                                                                     !
!                  Least Square Method solution (LSM)                 !
!                                                                     !
!*********************************************************************!

      IF (ChooseFreeMethod ==1) THEN
         !_______________________________________
         ! New values H*qu
         do i=1,N_CELL
            qu(i) = Hpr(i)*(gamma*qu_new(i)+(1.0d0-gamma)*qun(i))
            qv(i) = Hpr(i)*(gamma*qv_new(i)+(1.0d0-gamma)*qvn(i))
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
           call communication2D(qu)
           call communication2D(qv)
#        endif
!        ====================================
!        ====================================
         !_______________________________________
         ! Solution
         do i=1,N_CELL0
           !--------------------
           ! Derivatives by LSM
           sumqux = 0.0d0
           sumquy = 0.0d0
           sumqvx = 0.0d0
           sumqvy = 0.0d0
           do j=1,3
              jc = No_cp(i,j)
              sumqux = sumqux + dxCC(i,j)*(qu(jc)-qu(i))
              sumquy = sumquy + dyCC(i,j)*(qu(jc)-qu(i))
              sumqvx = sumqvx + dxCC(i,j)*(qv(jc)-qv(i))
              sumqvy = sumqvy + dyCC(i,j)*(qv(jc)-qv(i))
           enddo
           deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
           dqudx = (sum_yc2(i)*sumqux-sum_xcyc(i)*sumquy)/deter
           dqvdy = (sum_xc2(i)*sumqvy-sum_xcyc(i)*sumqvx)/deter
           !-------------------
           ! Update
           Hnew(i) = Hpr(i) - dt*(dqudx+dqvdy)
         enddo

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
           call communication2D(Hnew)
#        endif
!        ====================================
!        ====================================

!*********************************************************************!
!                                                                     !
!         Finite Volume Method: Explicit and implicit solution        !
!                                                                     !
!*********************************************************************!

      ELSEIF (ChooseFreeMethod ==2) THEN

!         ________________________________________________________
!        |                                                        |
!        |                 Full Explicit Solution                 |
!        |________________________________________________________|

#        if defined(KeyFreeFullExplicit)

            call FS_advection2D(Am0,Am1,Am2,Am3,AmG,&
                                qu_new,qv_new,Hpr,xc,yc,No_cp,nbe) 
!           ______________________
!           Update inside values
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
!              -------------------
!              New rhs
               som =  Am0(i)*Hpr(i)     &
                    + Am1(i)*Hpr(jc1)   &
                    + Am2(i)*Hpr(jc2)   &
                    + Am3(i)*Hpr(jc3)   &
                    + AmG(i) 
               bm(i) = -som/areaCell(i)
!              -------------------
!              Update
               Hnew(i) =  Hpr(i) + dt*bm(i)
            enddo

!         ________________________________________________________
!        |                                                        |
!        |                     Explicit Solution                  |
!        |________________________________________________________|

#        elif defined(KeyFreeExplicit)

            call FS_advection2D(Am0,Am1,Am2,Am3,AmG,&
                                qu_new,qv_new,Hpr,xc,yc,No_cp,nbe)
!           ______________________
!           Update inside values
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
!              -------------------
!              New diagonal values
               Am0(i)= 1.0d0 + dt/areaCell(i)*Am0(i)
!              -------------------
!              New rhs
               som =  Am1(i)*Hpr(jc1)   &
                    + Am2(i)*Hpr(jc2)   &
                    + Am3(i)*Hpr(jc3)   &
                    + AmG(i) 
               bm(i) = -som/areaCell(i)
!              -------------------
!              Update
               Hnew(i) = (Hpr(i)+dt*bm(i))/Am0(i)
            enddo
!         ________________________________________________________
!        |                                                        |
!        |                     Implicit Solution                  |
!        |________________________________________________________|

#        elif defined(KeyFreeImplicit)

            call FS_advection2D(Am0_new,Am1_new,Am2_new,Am3_new,AmG_new,&
                                qu_new,qv_new,Hpr,&
                                xc,yc,No_cp,nbe)
            call FS_advection2D(Am0n,Am1n,Am2n,Am3n,AmGn,&
                                qun,qvn,Hpr, &
                                xc,yc,No_cp,nbe)
            do i=1,N_CELL0
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
               ct = dt/areaCell(i)
               Am0(i) = gamma*ct*Am0_new(i) + 1.0d0
               Am1(i) = gamma*ct*Am1_new(i)
               Am2(i) = gamma*ct*Am2_new(i)
               Am3(i) = gamma*ct*Am3_new(i)
               bm(i)  = (-(1.0d0-gamma)*ct*Am0n(i) + 1.0d0)*Hpr(i)  &
                         -(1.0d0-gamma)*ct*Am1n(i)*Hpr(jc1)         &
                         -(1.0d0-gamma)*ct*Am2n(i)*Hpr(jc2)         &
                         -(1.0d0-gamma)*ct*Am3n(i)*Hpr(jc3)         &
                         -(1.0d0-gamma)*AmGn(i)                     &
                         - gamma*AmG_new(i)
            enddo
            
!           _________________________________________________________ 
!           Linear system solution  JSOR
#           ifdef KeyFS_PoissonJSOR
            Hnew = Hpr
            it=0
111         continue
            it=it+1
            errorsys = 0.0d0
            do i=1,N_CELL0
               !if (nbe(i).le.1) then
                  jc1 = No_cp(i,1)
                  jc2 = No_cp(i,2)
                  jc3 = No_cp(i,3)
                  som =  Am1(i)*Hnew(jc1)&
                       + Am2(i)*Hnew(jc2)&
                       + Am3(i)*Hnew(jc3)
                  residu   = (bm(i)-som)/Am0(i)-Hnew(i)
                  errorsys = errorsys + abs(residu)
                  Hnew(i)  = Hnew(i) + WLrelaxSOR*residu
               !endif
            enddo
!           ====================================
!           =====  START PARALLEL OPTION =======
#           ifdef KeyParallel
            call communication2D(Hnew)
#           endif
!           =============== END ================
!           ====================================
!           ----------------------------
            do nv=1,N_VERT
               !IF (nbev(nv).le.1) THEN
                  Hprv(nv) = 0.0d0
                  do j=1,Dimsurrounding(nv)
                     nc = surrounding(nv,j)
                     Hprv(nv)=Hprv(nv)+weight(nv,j)/dlVsum(nv)*Hnew(nc)
                  enddo
               !ENDIF
            enddo
!           ----------------------------
!           Free surface BC (reflecting)
            call Waterlevel_BC(Hnew,Hprv,eta,etav, &
                              h,xc,yc,No_cp,nbe,   &
                              hv,xv,yv,No_vp,nbev)
                              
!           ====================================
!           =====  START PARALLEL OPTION =======
#           ifdef KeyParallel
            call communication2D(Hnew)
            !--------------------------
            call SUM_parallel(errorsys,SUMerrorsys)
            errorsys = SUMerrorsys
#           endif
!           =============== END ================
!           ====================================

!           ----------------------------
!           Stop criteria
            if (errorsys.lt.eps) then
               IF (Display.eq.1) THEN
               write(*,6) 'Water level JS0R: iters =',it,&
                          ', error =',errorsys
               ENDIF
            elseif (errorsys.gt.1.0d5) then
               IF (Display.eq.1) THEN
               write(*,6) 'DIVERGENCE !!!!!!: iters =',it,&
                          ', error =',errorsys
               ENDIF
               stop
            elseif(it.gt.MaxIters) then
               IF (Display.eq.1) THEN
               write(*,6) 'Non-convergence: Maxiters =',it,&
                          ', error =',errorsys
               ENDIF
            else
               goto 111
            endif                   
6           format(t10,a24,i5,a9,e10.3)

#           endif

!           _________________________________________________________
!           Linear system solution  Jacobi
#           ifdef KeyFS_PoissonJacobi
            it=0
111         continue
            it=it+1
            errorsys = 0.0d0
            do i=1,N_CELL0
               !if (nbe(i).le.1) then
                  jc1 = No_cp(i,1)
                  jc2 = No_cp(i,2)
                  jc3 = No_cp(i,3)
                  som =  Am1(i)*Hpr(jc1)&
                       + Am2(i)*Hpr(jc2)&
                       + Am3(i)*Hpr(jc3)
                  residu   = (bm(i)-som)/Am0(i)-Hpr(i)
                  errorsys = errorsys + abs(residu)
                  Hnew(i)  = Hpr(i) + residu
               !endif
            enddo
!           ====================================
!           =====  START PARALLEL OPTION =======
#           ifdef KeyParallel
            call communication2D(Hnew)
#           endif	
!           =============== END ================
!           ====================================
!           ----------------------------
            do nv=1,N_VERT
               !IF (nbev(nv).le.1) THEN
                  Hprv(nv) = 0.0d0
                  do j=1,Dimsurrounding(nv)
                     nc = surrounding(nv,j)
                     Hprv(nv)=Hprv(nv)+weight(nv,j)/dlVsum(nv)*Hnew(nc)
                  enddo
               !ENDIF
            enddo
!           ----------------------------
!           Free surface BC (reflecting)
            call Waterlevel_BC(Hnew,Hprv,eta,etav, &
                              h,xc,yc,No_cp,nbe,   &
                              hv,xv,yv,No_vp,nbev)
!           ====================================
!           =====  START PARALLEL OPTION =======
#           ifdef KeyParallel 
            call SUM_parallel(errorsys,SUMerrorsys)
            errorsys = SUMerrorsys
#           endif	
!           =============== END ================
!           ====================================
!           ----------------------------
!           Update
            Hpr = Hnew
!           ----------------------------
!           Stop criteria
            if (errorsys.lt.eps) then
               IF (Display.eq.1) THEN
               write(*,6) 'Water level Jacobi: iters =',it,&
                          ', error =',errorsys
               ENDIF
            elseif (errorsys.gt.1.0d5) then
               IF (Display.eq.1) THEN
               write(*,6) 'DIVERGENCE !!!!!!: iters =',it,&
                          ', error =',errorsys
               ENDIF
               stop
            elseif(it.gt.MaxIters) then
               IF (Display.eq.1) THEN
               write(*,6) 'Non-convergence: Maxiters =',it,&
                          ', error =',errorsys
               ENDIF
            else
               goto 111
            endif                   
6           format(t10,a26,i5,a9,e10.3)
            
#           endif
#        endif

      ENDIF

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                      Final solution                    |
!     |________________________________________________________|

      do i=1,N_CELL0
         Hpr_new(i) = Hnew(i)
         eta_new(i) = Hnew(i) - h(i)  
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                    Boundary condition                  |
!     |________________________________________________________|
      
      call WaterLevel_BC(Hpr_new,Hprv,eta_new,etav, &
                         h,xc,yc,No_cp,nbe,         &
                         hv,xv,yv,No_vp,nbev)

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication2D(Hpr_new)
        call communication2D(eta_new)
#     endif
!     =============== END ================
!     ====================================
                         
!      ________________________________________________________
!     |                                                        |
!     |                 Inside vertex approximation            |
!     |________________________________________________________|

      do nv=1,N_VERT
         IF ((nbev(nv).le.1).or.(nbev(nv).eq.4)) THEN
         !IF (nbev(nv).le.1) THEN
            som1 = 0.0d0
            do j=1,Dimsurrounding(nv)
               nc = surrounding(nv,j)
               c1 = weight(nv,j)/dlVsum(nv)
               som1 = som1 + c1*Hpr_new(nc)
            enddo
            Hprv(nv) = som1
            etav(nv) = Hprv(nv)- hv(nv)
         ENDIF
      enddo

!      ________________________________________________________
!     |                                                        |
!     |         NUMERICAL TRICK: Maximum water level           |
!     |________________________________________________________|

!     For the cylinder the maximum water level will be given by
!     the point from the front of the cylinder. This point is
!     already calculated at the "FS_SaveReference" subroutine and
!     it is given by nv=No_sp(3); however, we can not applied 
!     directly because at each subdomain this value is different.
      
#     ifdef KeyMaximumEta
#        ifdef KeyStaticCylinder
         !IF (time.le.1.0d0) THEN
            xvR1 = -0.6d0
            yvR1 =  0.0d0
            nv1 = 1
            DR0 = sqrt((xvR1-xv(1))**2+(yvR1-yv(1))**2)      
            do nv=2,N_VERT
               DR = sqrt((xvR1-xv(nv))**2+(yvR1-yv(nv))**2)
               if (DR.lt.DR0) then
                  nv1 = nv
                  DR0 = DR
               endif 
            enddo
!           ====================================
!           =====  START PARALLEL OPTION =======           
#           ifdef KeyParallel
            pro = rang_topo
            call MPI_ALLGATHER(pro,1,MPI_INTEGER,proV,1,MPI_INTEGER,comm3D,code)
            call MPI_ALLGATHER(nv1,1,MPI_INTEGER,nv1V,1,MPI_INTEGER,comm3D,code)
            call MPI_ALLGATHER(DR0,1,MPI_DOUBLE_PRECISION,DR0V,1,MPI_DOUBLE_PRECISION,comm3D,code)
            pro = proV(1)
            nv1 = nv1V(1)
            DR0 = DR0V(1)
            do n=2,Nprocs
              if (DR0V(n).lt.DR0) then
                  nv1 = nv1V(n)
                  pro = proV(n)
                  DR0 = DR0V(n)
               endif 
            enddo
            etaMax = 1.25d0
            if (rang_topo.eq.pro) then
                etaMax = 1.25d0*abs(etav(nv1))
                print*,'        MAXIMUM W-L APPLIED: etaMax=',etaMax
            endif
            call MPI_Bcast(etaMax,1,MPI_DOUBLE_PRECISION,pro,comm3D,code)  
!           ====================================
!           ====================================
#           else
            etaMax = 1.2d0*abs(etav(nv1))
#           endif
!           =============== END ================
!           ====================================
            do nv=1,N_VERT
               etav(nv) = min(etav(nv),etaMax)
               Hprv(nv) = etav(nv) + hv(nv)
            enddo
            do i=1,N_CELL
               eta_new(i) = min(eta_new(i),etaMax)
               Hpr_new(i) = eta_new(i) + h(i)  
            enddo
         !ENDIF
#        endif
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                    Boundary condition                  |
!     |________________________________________________________|
      
      call WaterLevel_BC(Hpr_new,Hprv,eta_new,etav, &
                         h,xc,yc,No_cp,nbe,         &
                         hv,xv,yv,No_vp,nbev)

!      ________________________________________________________
!     |                                                        |
!     |                     Time derivatives                   |
!     |________________________________________________________|

      do i=1,N_CELL
         dHprdt_new(i) = (Hpr_new(i)-Hprn(i))/dt
         detadt_new(i) = (eta_new(i)-etan(i))/dt 
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication2D(dHprdt_new)
        call communication2D(detadt_new)
#     endif
!     ====================================
!     ====================================
      do nv=1,N_VERT
         som1 = 0.0d0
         som2 = 0.0d0
         do j=1,Dimsurrounding(nv)
            nc = surrounding(nv,j)
            c1 = weight(nv,j)/dlVsum(nv)
            som1 = som1 + c1*dHprdt_new(nc)
            som2 = som2 + c1*detadt_new(nc)
         enddo
         dHprdtv(nv) = som1
         detadtv(nv) = som2 
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                   Spatial derivatives                  |
!     |________________________________________________________|

!     ------------------------------------
!     Cell-centered
      call grandientLSM2D(detadx_new,detady_new,eta_new,No_cp,nbe)
      call grandientLSM2D(dHprdx_new,dHprdy_new,Hpr_new,No_cp,nbe)
!     ------------------------------------
!     Vertex 
      do nv=1,N_VERT
         som1 = 0.0d0
         som2 = 0.0d0
         som3 = 0.0d0
         som4 = 0.0d0
         do j=1,Dimsurrounding(nv)
            nc = surrounding(nv,j)
            c1 = weight(nv,j)/dlVsum(nv)
            som1 = som1 + c1*detadx_new(nc)
            som2 = som2 + c1*detady_new(nc)
            som3 = som3 + c1*dHprdx_new(nc)
            som4 = som4 + c1*dHprdy_new(nc)
         enddo
         detadxv(nv) = som1
         detadyv(nv) = som2
         dHprdxv(nv) = som3
         dHprdyv(nv) = som4
      enddo

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: FS_WaterLevel'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                CALCULATE THE ADVECTION COEFFICIENTS (2D)            !
!                             May 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_advection2D(Am0,Am1,Am2,Am3,AmG,& 
                                uu,vv,              &
                                phi,xc,yc,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the advection contribution to the   !
!    general linear system. We have two options to the final values   !
!    of Am and AmG, depending if we choose between an implicit or     !
!    explicit scheme (cppdefs.h).                                     !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> Am0   |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  | <--> Am1   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  | <--> Am2   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  | <--> Am3   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  | <--> AmG   |(N_CELL0,NZ)| Vector with Gradient terms          |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    | Description                         |  !  
!  |____________|____________|_____________________________________|  !
!  | --> uu     |(N_CELL,NZ) | Velocity component u                |  !
!  | --> vv     |(N_CELL,NZ) | Velocity component v                |  !
!  |____________|____________|_____________________________________|  !
!  | --> phi    |(N_CELL,NZ) | Function phi at the center element  |  !
!  |____________|____________|_____________________________________|  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Type of boundary cell (inside or bc)|  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name        |               Description                     |  !  
!  |_______________|_______________________________________________|  !  
!  | FluxLimiter   | (N_CELL0,NZ) Flux limiter (Vidovic 2006)      |  !
!  | FluxLimiterMin| Minimum value for the flux limiter            |  !
!  | FluxLimiterpos| Positive side of the flux limiter             |  !
!  | FluxLimiterneg| Negative side of the flux limiter             |  !
!  | FluxLimiterij | Pre-final value of the flux limiter           |  !
!  | phimin        | Maximum phi of all cell-centers related       |  !
!  | phimax        | Maximum phi of all cell-centers related       |  !
!  | d1,d2         | auxiliar diference values of phi              |  !
!  |_______________|_______________________________________________|  !
!  | dphidx        |(N_CELL,NZ) d(phi)/dx   = gradient component x |  !
!  | dphidy        |(N_CELL,NZ) d(phi)/dy   = gradient component y |  ! 
!  | dphidsig      |(N_CELL,NZ) d(phi)/dsig = gradient component z |  ! 
!  |_______________|_______________________________________________|  ! 
!  | C01           |(N_CELL0,NZ) Mass flux of horizontal neigh. 1  |  !
!  | C02           |(N_CELL0,NZ) Mass flux of horizontal neigh. 2  |  !
!  | C03           |(N_CELL0,NZ) Mass flux of horizontal neigh. 3  |  !
!  |_______________|_______________________________________________|  !
!  | aa0,aa(1:3)   | Auxiliar horizontal matrix coefficients       |  !
!  | C0j(1:3)      | Auxiliar mass flux values of the neighbors    |  !
!  | C0pos,C0neg   | Positive and negative values of the mass flux |  !
!  | GF            | Gradient (rhs contribution)                   |  ! 
!  | GFi,GFj       | Gradient at the cell center & neigborn        |  !
!  | GFpos,GFneg   | Positive and negative values of the gradient  |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!

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

      real*8, dimension(:) :: Am0(N_CELL0)
      real*8, dimension(:) :: Am1(N_CELL0)
      real*8, dimension(:) :: Am2(N_CELL0)
      real*8, dimension(:) :: Am3(N_CELL0)
      real*8, dimension(:) :: AmG(N_CELL0)
!     --------------------------------------
      real*8, dimension(:) :: uu(N_CELL)
      real*8, dimension(:) :: vv(N_CELL)
      real*8, dimension(:) :: phi(N_CELL)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8 ,dimension(:),allocatable :: dphidx
      real*8 ,dimension(:),allocatable :: dphidy
      real*8 :: sumfx,sumfy,deter
!     --------------------------------------
      real*8 ,dimension(:,:) :: Uij(N_CELL0,3)
      real*8 :: nxLength,nyLength
      real*8 :: u0j,v0j
!     --------------------------------------
      real*8 ,dimension(:),allocatable :: FluxLimiter
      real*8 :: FluxLimiterMin 
      real*8 :: FluxLimiterpos
      real*8 :: FluxLimiterneg
      real*8 :: FluxLimiterij
      real*8 :: phimin
      real*8 :: phimax
      real*8 :: d1,d2
!     --------------------------------------
      real*8 ,dimension(:) :: aa(1:3)
      real*8 :: aa0,aaB,aaT
      real*8 :: Ue,Uijneg,Uijpos
      real*8 :: GF,GFi,GFj
      real*8 :: GFpos,GFneg
      integer:: jc,elem
!     --------------------------------------
      real*8 :: x,y,dfdxExam2D,dfdyExam2D,funSolExam2D

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: Advection2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication2D(phi)
        call communication2D(uu)
        call communication2D(vv)
#     endif
!     ====================================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(dphidx(N_CELL),       &
               dphidy(N_CELL),       &
               FluxLimiter(N_CELL))


!*********************************************************************!
!                                                                     !
!                                Gradient                             !
!                                                                     !
!*********************************************************************!

      do i=1,N_CELL	
         dphidx(i)=0.0d0
         dphidy(i)=0.0d0
      enddo

      DO i=1,N_CELL0
!        __________________________________________________
!        Inside & boundary cell-centers
         sumfx = 0.0d0
         sumfy = 0.0d0
         do j=1,3
            jc = No_cp(i,j)
            sumfx = sumfx + dxCC(i,j)*(phi(jc)-phi(i))
            sumfy = sumfy + dyCC(i,j)*(phi(jc)-phi(i))
         enddo
         deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
         dphidx(i)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
         dphidy(i)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
!        __________________________________________________
!        Outside cell-centers  (ficticious)
!        ----------------------------------
         do j=1,3
           nc=No_cp(i,j)
!              ====================================
!              ==========  SEQUENTIAL =============
#              ifndef KeyParallel
               if (nc.lt.1.OR.nc.gt.N_CELL0) then
!              ====================================
!              =====  START PARALLEL OPTION =======
#              else
               elem = index_global(nc)
               if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!              =============== END ================
!              ====================================
              dphidx(nc) = -dphidx(i) 
              dphidy(nc) = -dphidy(i)
           endif
         enddo
      ENDDO

!        ====================================
!        ====================================
#        ifdef KeyParallel
            call communication2D(dphidx)
            call communication2D(dphidy)
#        endif
!        ====================================
!        ====================================

!*********************************************************************!
!                                                                     !
!                        Mass flux calculation                        !
!                                                                     !
!*********************************************************************!

      do i=1,N_CELL0
         do j=1,3
!           ---------------------------------
!           Index of the neighbor cell-center
            jc = No_cp(i,j)
!           ---------------------------------
!           Velocity at the faces 
            u0j = 0.5d0*(uu(i)+uu(jc))
            v0j = 0.5d0*(vv(i)+vv(jc))
!           ---------------------------------
!           Mass fluxes of the cell
            nxLength =  dyVV(i,j)
            nyLength = -dxVV(i,j)
            Uij(i,j) = u0j*nxLength+v0j*nyLength
         enddo
      enddo

!*********************************************************************!
!                                                                     !
!                       Flux limiter calculation                      !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                 Using the LED technique                |
!     |________________________________________________________|

      FluxLimiter = 1.0d0
#     ifdef KeyFluxLimiterWL
         do i=1,N_CELL0
!           -------------------------------
!           Maximum & minimum value of phi
            phimax = phi(i)
            phimin = phi(i)
            do j=1,3
               jc = No_cp(i,j)
               phimax = max(phimax,phi(jc))
               phimin = min(phimin,phi(jc))
            enddo
!           -------------------------------
!           Flux limiter calculation 
            FluxLimiterMin = 1.0d5
            do j=1,3 
               GFi   = dphidx(i)*dxCE(i,j)+dphidy(i)*dyCE(i,j)
               GFpos = 0.5d0*(GFi+abs(GFi))
               GFneg = 0.5d0*(GFi-abs(GFi))
               if(dabs(GFi).lt.1.0d-10) then
                  FluxLimiterij = 1.0d0
               else
                  d1 = (phimax-phi(i))/GFi
                  d2 = (phimin-phi(i))/GFi
                  FluxLimiterpos = min(1.0d0,d1)
                  FluxLimiterneg = min(1.0d0,d2)
                  if (GFi.gt.0.0d0) then
                      FluxLimiterij = FluxLimiterpos
                  else
                      FluxLimiterij = FluxLimiterneg
                  endif
               endif
               FluxLimiterMin = min(FluxLimiterMin,FluxLimiterij)
            enddo
            FluxLimiter(i) = FluxLimiterMin	  
!           -------------------------------
!           Apply to BC (We should not do it [2017])
            do j=1,3
               jc = No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (jc.lt.1.OR.jc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(jc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif	
!                 =============== END ================
!                 ====================================
                  FluxLimiter(jc) = FluxLimiter(i)
               endif
             enddo 
           enddo
#     endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(FluxLimiter)
#     endif 
!     ====================================
!     ====================================

!*********************************************************************!
!                                                                     !
!                       Advection contributions                       !
!                                                                     !
!*********************************************************************!

      do i=1,N_CELL0
!         __________________________________________________
!        |                                                  |
!        |                  Matrix & rhs                    |
!        |__________________________________________________|
           
         aa0 = 0.0d0
         GF  = 0.0d0
         do j=1,3
!           ---------------------------------
!           Number of the neighbor cell-center 
            jc = No_cp(i,j)
!           ---------------------------------
!           Mass fluxes of the cell
            Ue = Uij(i,j) 
            Uijpos =  0.5d0*(Ue+abs(Ue))
            Uijneg =  0.5d0*(Ue-abs(Ue))
!           ---------------------------------
!           Matrix coeffients
            aa0   =  Uijpos + aa0 
            aa(j) =  Uijneg          
!           ---------------------------------
!           Gradient at of the cell & the neighbor
            GFi =  dphidx(i)*(xme(i,j)-xc(i)) &
                 + dphidy(i)*(yme(i,j)-yc(i))

            GFj =  dphidx(jc)*(xme(i,j)-xc(jc)) &
                 + dphidy(jc)*(yme(i,j)-yc(jc))
!           ---------------------------------
!           Gradient coefficient (rhs)
            GF = GF + Uijpos*GFi*FluxLimiter(i)  &
                    + Uijneg*GFj*FluxLimiter(jc)
         enddo
!         __________________________________________________
!        |                                                  |
!        |                  Contributions                   |
!        |__________________________________________________|

!        --------------------------------------
!        Matrix
         Am0(i) =  aa0
         Am1(i) =  aa(1)
         Am2(i) =  aa(2)
         Am3(i) =  aa(3)
!        --------------------------------------
!        Known part to rhs
         AmG(i) =  GF
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

      deallocate(dphidx,dphidy, &
                 FluxLimiter)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: Advection2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
      
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF WATER ELEVATION                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
