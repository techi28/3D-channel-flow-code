!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    TEST NAVIER-STOKES EQUATION                      !
!                              May 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Hydro(ufnp,vfnp,wfnp,pfnp,                  & 
                       uf,vf,wf,pf,                          & 
                       ufv,vfv,wfv,pfv,                      &
!                      -------------------------------------  
                       Hprnp,etanp,                          &                             
                       Hpr,eta,                              &
                       Hprv,etav,                            &
!                      ------------------------------------- 
                       h,hv,                                 &                              
!                      -------------------------------------                               
                       xc,yc,sig,dsig,No_cp,nbe,             &
                       xv,yv,sigv,dsigv,No_vp,nbev,          &
!                      ------------------------------------- 
                       No_wb,No_qb,No_hb,No_sp,              &
!                      ------------------------------------- 
                       nRK)               

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the 3D Navier-Stokes equations with velo-   !
!    city variables u, v, w and pressure p. The program also calcu-   !
!    lates the free-surface elevation by introducing a new loop.      !
!---------------------------------------------------------------------!
!                                                                     !
!    There are three different methods to decouple the velocity and   !
!    the pressure:                                                    !
!    Key_ProjectionFirst:     Using the 1st projection meth.          !
!    Key_ProjectionSecond:    Using the 2nd projection meth.          !
!    Key_PredictorCorrector:  Using the predictor-corrector meth.     !
!                                                                     !
!---------------------------------------------------------------------!
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
!  | ---> uf    |(N_CELL,NZ)  | u Old solution of the equation     |  !
!  | ---> vf    |(N_CELL,NZ)  | v Old solution of the equation     |  !  
!  | ---> wf    |(N_CELL,NZ)  | w Old solution of the equation     |  !  
!  | ---> pf    |(N_CELL,NZ)  | p Old solution of the equation     |  !
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
!      ________________________________________________________       !
!      sigma-transform variables                                      !
!                                                                     !
!      sigmat    :  transformation sigma in the time                  !
!      sigmax    :  transformation sigma in the x direction           !
!      sigmay    :  transformation sigma in the y direction           !
!      sigmaz    :  transformation sigma in the z direction           !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - AdvDiffEqn3D                ( AdvDiffEqn3D.F90 )          |  !
!  |   - Poisson                     ( Poisson.F90 )               |  !
!  |   - grandientLSM                ( grandientLSM.F90 )          |  !
!  |   - interpolation3D             ( interpolation3D.F90 )       |  !
!  |   - BCvelcenter3D               ( BCvelocity.F90 )            |  !
!  |   - BCvelvertex3D               ( BCvelocity.F90 )            |  !
!  |   - BCprescenter3D              ( BCpressure.F90 )            |  !
!  |   - BCpresvertex3D              ( BCpressure.F90 )            |  !
!  |_______________________________________________________________|  !
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

      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)

      real*8,dimension(:,:) :: uf(N_CELL,NZ)
      real*8,dimension(:,:) :: vf(N_CELL,NZ)
      real*8,dimension(:,:) :: wf(N_CELL,NZ)
      real*8,dimension(:,:) :: pf(N_CELL,NZ)
            
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)   
!     --------------------------------------
      real*8, dimension(:)  :: Hprnp(N_CELL)
      real*8, dimension(:)  :: etanp(N_CELL)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     -------------------------------------      
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: hv(N_VERT)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 

      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     -------------------------------------
      integer,dimension(:)  :: No_wb(N_WB)
      integer,dimension(:)  :: No_qb(N_QB)
      integer,dimension(:)  :: No_hb(N_HB)
      integer,dimension(:)  :: No_sp(N_SPmax)
!     --------------------------------------      
      integer :: nRK
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8, dimension(:)  :: Hprn(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL) 
      real*8, dimension(:)  :: Hpr_new(N_CELL)
      real*8, dimension(:)  :: eta_new(N_CELL) 
      real*8, dimension(:)  :: Hpr_old(N_CELL)    
      real*8, dimension(:)  :: eta_old(N_CELL)
      real*8, dimension(:)  :: Hpr_gam(N_CELL)    
      real*8, dimension(:)  :: eta_gam(N_CELL)
!     -------- -----------------------------
      real*8, dimension(:)  :: detadx_gam(N_CELL)
      real*8, dimension(:)  :: detady_gam(N_CELL)           
!     -------- -----------------------------       
      real*8,dimension(:,:) :: Hufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Hvfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Hwfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Huf(N_CELL,NZ)
      real*8,dimension(:,:) :: Hvf(N_CELL,NZ)
      real*8,dimension(:,:) :: Hwf(N_CELL,NZ)      
      real*8,dimension(:,:) :: Hufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hvfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hwfv(N_VERT,NZ-1)
!     -------------------------------------  
      real*8,dimension(:,:) :: omef(N_CELL,NZ)
      real*8,dimension(:,:) :: Homef(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
!     --------------------------------------
      real*8,dimension(N_CELL,NZ) :: Sdudx,Sdudy,Sdudz
      real*8,dimension(N_CELL,NZ) :: Sdvdx,Sdvdy,Sdvdz
      real*8,dimension(N_CELL,NZ) :: Sdwdx,Sdwdy,Sdwdz
      real*8,dimension(N_CELL,NZ) :: Sdpdx,Sdpdy,Sdpdz                  
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: Hufstar,Hvfstar,Hwfstar      
      real*8,dimension(:,:),allocatable :: ufstar,vfstar,wfstar
      real*8,dimension(:,:),allocatable :: Dpfnp,Dpf,Dpfnpv,Dpfv
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: dudx,dudy,duds
      real*8,dimension(:,:),allocatable :: dvdx,dvdy,dvds
      real*8,dimension(:,:),allocatable :: dwdx,dwdy,dwds
      real*8,dimension(:,:),allocatable :: dpdx,dpdy,dpds
      real*8,dimension(:,:),allocatable :: dDpdx,dDpdy,dDpds
!     --------------------------------------
      real*8,dimension(:,:),allocatable :: Gamux,Gamuy,Gamuz
      real*8,dimension(:,:),allocatable :: Gamvx,Gamvy,Gamvz
      real*8,dimension(:,:),allocatable :: Gamwx,Gamwy,Gamwz            
      real*8,dimension(:,:),allocatable :: Gampx,Gampy,Gampz
      real*8,dimension(:,:),allocatable :: rhsu,rhsv,rhsw,rhsp
!     --------------------------------------
      real*8,dimension(:),allocatable :: Shu,Shv,Shw
      real*8,dimension(:),allocatable :: dHprdx,dHprdy,dHprdt 
      real*8,dimension(:),allocatable :: detadx,detady,detadt         
      real*8,dimension(:),allocatable :: dhdx,dhdy 
!     --------------------------------------
      real*8,dimension(:),allocatable :: dHprdxv,dHprdyv,dHprdtv 
      real*8,dimension(:),allocatable :: detadxv,detadyv,detadtv      
      real*8,dimension(:),allocatable :: dhdxv,dhdyv 
!     --------------------------------------      
      real*8,dimension(:),allocatable :: gradetadx,gradetady
      real*8,dimension(:),allocatable :: dgraxdx,dgraxdy
      real*8,dimension(:),allocatable :: dgraydx,dgraydy
!     --------------------------------------
      real*8, dimension(:,:) :: dpdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdsn(N_CELL,NZ)       
      real*8, dimension(:,:) :: Sdpdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdzn(N_CELL,NZ)
 !     --------------------------------------     
      real*8,dimension(:) :: w1vn(N_VERT)     
      real*8,dimension(:) :: dw1dtv(N_VERT)           
!     --------------------------------------     
      real*8 :: c1,nu,Vol
      integer:: tagBC
!     --------------------------------------
      real*8 :: funExamNSu,funExamNSv,funExamNSw,x,y,z
!     --------------------------------------          
      integer:: iconv
      integer:: jv1,jv2,jv3      
      real*8 :: epsp,SUMepsp
      real*8 :: uB,vB,uT,vT
!     ____________________________________
!    |                                    |
!    |       Free surface parameters      |
!    |____________________________________|

      real*8, parameter :: gamma   = 1.0d0
      real*8, parameter :: tolConv = 1d-8
      integer,parameter :: nconv   = 1         
            
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: The N-S Equations'
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

      allocate(Hufstar(N_CELL,NZ),Hvfstar(N_CELL,NZ),Hwfstar(N_CELL,NZ),&
               ufstar(N_CELL,NZ),vfstar(N_CELL,NZ),wfstar(N_CELL,NZ), &
               Dpfnp(N_CELL,NZ),Dpf(N_CELL,NZ),                       &
               Dpfnpv(N_VERT,NZ-1),Dpfv(N_VERT,NZ-1),                 &
!              --------------------------------------
               dudx(N_CELL,NZ),dudy(N_CELL,NZ),duds(N_CELL,NZ),       &
               dvdx(N_CELL,NZ),dvdy(N_CELL,NZ),dvds(N_CELL,NZ),       &
               dwdx(N_CELL,NZ),dwdy(N_CELL,NZ),dwds(N_CELL,NZ),       &
               dpdx(N_CELL,NZ),dpdy(N_CELL,NZ),dpds(N_CELL,NZ),       &
               dDpdx(N_CELL,NZ),dDpdy(N_CELL,NZ),dDpds(N_CELL,NZ),    &
!              --------------------------------------
               Gamux(N_CELL,NZ),Gamuy(N_CELL,NZ),Gamuz(N_CELL,NZ),    &
               Gamvx(N_CELL,NZ),Gamvy(N_CELL,NZ),Gamvz(N_CELL,NZ),    &
               Gamwx(N_CELL,NZ),Gamwy(N_CELL,NZ),Gamwz(N_CELL,NZ),    &
               Gampx(N_CELL,NZ),Gampy(N_CELL,NZ),Gampz(N_CELL,NZ),    &
!              --------------------------------------               
               rhsu(N_CELL,NZ),rhsv(N_CELL,NZ),rhsw(N_CELL,NZ),       &
               rhsp(N_CELL,NZ),                                       &
!              --------------------------------------
               Shu(N_CELL),Shv(N_CELL),Shw(N_CELL),                   &
               dHprdx(N_CELL),dHprdy(N_CELL),dHprdt(N_CELL),          &
               detadx(N_CELL),detady(N_CELL),detadt(N_CELL),          &
               dhdx(N_CELL),dhdy(N_CELL),                             &
!              --------------------------------------
               dHprdxv(N_VERT),dHprdyv(N_VERT),dHprdtv(N_VERT),       & 
               detadxv(N_VERT),detadyv(N_VERT),detadtv(N_VERT),       &               
               dhdxv(N_VERT),dhdyv(N_VERT),                           & 
!              --------------------------------------
               gradetadx(N_CELL),dgraxdx(N_CELL),dgraydx(N_CELL),     &
               gradetady(N_CELL),dgraxdy(N_CELL),dgraydy(N_CELL))
    
!      ________________________________________________________
!     |                                                        |
!     |              Initialization of variables               |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(uf)
         call communication3Dtype2(vf)
         call communication3Dtype2(wf)
         call communication3Dtype2(pf)
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________
!     |                                               |
!     |  I.1)  Initial free-surface & velocity        |
!     |_______________________________________________|

!     --------------
!     Free surface
      Hprn    = Hpr
      etan    = eta
      Hprnp   = Hpr
      etanp   = eta
      Hpr_new = Hpr
      eta_new = eta
      Hpr_gam = Hpr
      eta_gam = eta
      Hpr_old = Hpr
      eta_old = eta
!     --------------
!     Velocity
      ufnp = uf
      vfnp = vf
      wfnp = wf
!     --------------
!     Preessure
      pfnp = pf

!      _______________________________________________
!     |                                               |
!     |  I.2) Initial velocities (Hu,Hv,Hw)           |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL  
            Huf(i,k) = uf(i,k)*Hpr(i)
            Hvf(i,k) = vf(i,k)*Hpr(i)
            Hwf(i,k) = wf(i,k)*Hpr(i)
         enddo
      enddo
      do k=1,NZ-1
         do nv=1,N_VERT  
            Hufv(nv,k) = ufv(nv,k)*Hprv(nv)
            Hvfv(nv,k) = vfv(nv,k)*Hprv(nv)
            Hwfv(nv,k) = wfv(nv,k)*Hprv(nv)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(Huf)
         call communication3Dtype2(Hvf)
         call communication3Dtype2(Hwf)
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________
!     |                                               |
!     |  I.3) wfv^n at the bottom (for pressure BC)   |
!     |_______________________________________________|

      do nv=1,N_VERT 
         w1vn(nv) = wfv(nv,1)
      enddo

!      _______________________________________________
!     |                                               |
!     |  I.4) Spacial derivatives of h (2D)           |
!     |_______________________________________________|

!     --------------
!     Cell-centered
      call grandientLSM2D(dhdx,dhdy,h,No_cp,nbe)

!     --------------
!     Vertex values 
      call interpolation2D(dhdxv,xv,yv,No_vp,nbev, &
                           dhdx,xc,yc,No_cp,nbe)   
      call interpolation2D(dhdyv,xv,yv,No_vp,nbev, &
                           dhdy,xc,yc,No_cp,nbe)
            
!*********************************************************************!
!---------------------------------------------------------------------!
!             Opening: Iterative loop for free-surface update         !
!---------------------------------------------------------------------!
!*********************************************************************!
 
#     ifndef KeyFixedFreeSurface        

!      _______________________________________________
!     |                                               |
!     |  [FS.O.1]  Initialization of iterations       |
!     |_______________________________________________|

      iconv = 0
9999  continue
      iconv = iconv + 1

!      _______________________________________________
!     |                                               |
!     |  [FS.O.2]  Save current iteration             |
!     |_______________________________________________|

      eta_old = eta_new
      
!      _______________________________________________
!     |                                               |
!     |  [FS.O.3]  Update FS: dH/dt + Flux[H(u,v)] = 0|
!     |               °  (Hpr)^(new)                  |
!     |               °  (eta)^(new)                  | 
!     |_______________________________________________|
                   
      call WaterLevel(Hpr_new,Hprn,Hprv,              &
                      eta_new,etan,etav,              &
!                     ---------------------------------                                                                          
                      ufnp,vfnp,                      &
                      uf,vf,                          &
!                     ---------------------------------
                      h,hv,                           &
!                     --------------------------------- 
                      xc,yc,sig,dsig,No_cp,nbe,       &
                      xv,yv,sigv,dsigv,No_vp,nbev,    &
!                     ---------------------------------                      
                      gamma,                          &
!                     ---------------------------------
                      No_sp)
                    
!      _______________________________________________
!     |                                               |
!     |  [FS.O.4] Theta-method for FS variables       |
!     |_______________________________________________|
       
      Hpr_gam = gamma*Hpr_new+(1.0d0-gamma)*Hprn
      eta_gam = gamma*eta_new+(1.0d0-gamma)*etan

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication2D(Hpr_gam)
        call communication2D(eta_gam)
#     endif
!     ====================================
!     ====================================
       
#     endif 

!*********************************************************************!
!                                                                     !
!            Navier-Stokes solution with sigma-transformation         !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |             ---- 1) PRECALCULATIONS    ----            |
!     |                     °  dh/dx, dh/dy                    |
!     |                     °  d(eta)/dx^n, d(eta)/dy^n        |
!     |                     °  d(Hpr)/dx^n, d(Hpr)/dy^n        |
!     |                     °  sigma^(n)                       |
!     |                                                        |
!     |     If domain is fixed, the the values are the same    |
!     |     for each time step. Calculate once to efficiency.  |
!     |________________________________________________________|
!      ******************************************************** 

!      _______________________________________________
!     |                                               |
!     |  1.1) Time derivatives of eta & Hpr (2D)      |
!     |_______________________________________________|  
      
#     ifdef KeyFixedFreeSurface
         dHprdt  = 0.0d0  
         dHprdtv = 0.0d0               
         detadt  = 0.0d0   
         detadtv = 0.0d0
#     else
!        --------------
!        Cell-centered
         do i=1,N_CELL
            dHprdt(i) = (Hpr_new(i)-Hprn(i))/dt
            detadt(i) = (eta_new(i)-etan(i))/dt 
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication2D(dHprdt)
            call communication2D(detadt)
#        endif
!        =============== END ================
!        ====================================
!        --------------
!        Vertex values
         call interpolation2D(dHprdtv,xv,yv,No_vp,nbev, &
                              dHprdt,xc,yc,No_cp,nbe)   
         call interpolation2D(detadtv,xv,yv,No_vp,nbev, &
                              detadt,xc,yc,No_cp,nbe)
#     endif
      
!      _______________________________________________
!     |                                               |
!     |  1.2) Spacial derivatives of eta & Hpr (2D)   |
!     |_______________________________________________|
 
!     --------------
!     Cell-center
      call grandientLSM2D(detadx,detady,eta_new,No_cp,nbe)
      call grandientLSM2D(dHprdx,dHprdy,Hpr_new,No_cp,nbe)

!     --------------
!     Vertex
      call interpolation2D(detadxv,xv,yv,No_vp,nbev, &
                           detadx,xc,yc,No_cp,nbe)   
      call interpolation2D(detadyv,xv,yv,No_vp,nbev, &
                           detady,xc,yc,No_cp,nbe)
      
!      _______________________________________________
!     |                                               |
!     |  1.3)  Sigma transformation at sigma^(n)      |
!     |_______________________________________________|
       
      do k=1,NZ
         do i=1,N_CELL
            c1 = 1.0d0/Hpr_new(i)
            sigmat(i,k) = c1*(-sig(k)*detadt(i))
            sigmax(i,k) = c1*(dhdx(i)-sig(k)*dHprdx(i)) 
            sigmay(i,k) = c1*(dhdy(i)-sig(k)*dHprdy(i)) 
            sigmaz(i,k) = c1
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication3Dtype2(sigmat)
      call communication3Dtype2(sigmax)
      call communication3Dtype2(sigmay)
      call communication3Dtype2(sigmaz)
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________
!     |                                               |
!     |  1.4)  Special case: Test Only Poisson Eqn.   |
!     |_______________________________________________|

#     if defined(KeyTestOnlyPoisson)
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if (rang_topo.eq.0) then
         write(*,'(t10,a42)') 'WARNING!!! Only test Poisson executed!!!'
         write(*,'(t10,a42)') '      ---> Go directly to pressure step '
         write(*,*) ' '
         endif
!     ====================================
!     ====================================
#     else
         write(*,'(t10,a42)') 'WARNING!!! Only test Poisson executed!!!'
         write(*,'(t10,a42)') '      ---> Go directly to pressure step '
         write(*,*) ' '         
#     endif
      goto 11111
!     =============== END ================
!     ====================================
#     endif

!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |          ----  2) PREDICTION STEP: fluid  ----         |
!     |                                                        |
!     |                    °  Huf^(*), uf^(*)                  |
!     |                    °  Hvf^(*), vf^(*)                  |
!     |                    °  Hwf^(*), wf^(*)                  |
!     |________________________________________________________|
!      ********************************************************

!      _______________________________________________
!     |                                               |
!     |  2.1) New vertical velocity: omef & H*omef    |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL 
            omef(i,k) = sigmat(i,k) + uf(i,k)*sigmax(i,k) &
                                    + vf(i,k)*sigmay(i,k) &
                                    + wf(i,k)*sigmaz(i,k)
            Homef(i,k)=(sigmat(i,k) + uf(i,k)*sigmax(i,k) &
                                    + vf(i,k)*sigmay(i,k) &
                                    + wf(i,k)*sigmaz(i,k))*Hpr_new(i)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication3Dtype2(omef)
      call communication3Dtype2(Homef)
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________
!     |                                               |
!     |  2.2) Diffusive Coefficients: Gamu,Gamv,Gamw  |
!     |                NOT TURBULENCE                 |
!     |_______________________________________________|

      Gamux = 0.0
      Gamuy = 0.0
      Gamuz = 0.0
      Gamvx = 0.0
      Gamvy = 0.0
      Gamvz = 0.0
      Gamwx = 0.0
      Gamwy = 0.0
      Gamwz = 0.0
!     ________________________________________________
!     LES method 
#     ifdef KeyDiffusiveRe     
      do k=1,NZ
         do i=1,N_CELL
            nu = 1.0d0/Re
            Gamux(i,k) = nu
            Gamuy(i,k) = nu
            Gamuz(i,k) = nu
            Gamvx(i,k) = nu
            Gamvy(i,k) = nu
            Gamvz(i,k) = nu
            Gamwx(i,k) = nu
            Gamwy(i,k) = nu
            Gamwz(i,k) = nu
         enddo
      enddo
#     endif
!     ________________________________________________
!     LES method 
#     ifdef KeyDiffusiveLES
      call sgles(Gamux,Gamuy,Gamuz,              &
                 Gamvx,Gamvy,Gamvz,              &
                 Gamwx,Gamwy,Gamwz,              &
!                --------------------------------
                 ufnp,vfnp,wfnp,                 &
!                --------------------------------
                 Hpr_gam,Hpr_new,Hprn,Hprv,      &
!                --------------------------------
                 h,hv,                           &
!                --------------------------------
                 xc,yc,sig,dsig,No_cp,nbe,       &
                 xv,yv,sigv,dsigv,No_vp,nbev,    &
                 No_wb)  
#     endif
!     ________________________________________________
!     Sigma transform contribution
      do k=1,NZ
         do i=1,N_CELL
            Gamux(i,k) = Gamux(i,k)*Hpr_new(i)
            Gamuy(i,k) = Gamuy(i,k)*Hpr_new(i)
            Gamuz(i,k) = Gamuz(i,k)*Hpr_new(i)
            Gamvx(i,k) = Gamvx(i,k)*Hpr_new(i)
            Gamvy(i,k) = Gamvy(i,k)*Hpr_new(i)
            Gamvz(i,k) = Gamvz(i,k)*Hpr_new(i)
            Gamwx(i,k) = Gamwx(i,k)*Hpr_new(i)
            Gamwy(i,k) = Gamwy(i,k)*Hpr_new(i)
            Gamwz(i,k) = Gamwz(i,k)*Hpr_new(i)
         enddo
      enddo

!      _______________________________________________
!     |                                               |
!     |  2.3) Mass-flux (uu,vv,ww)                    |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL  
            uu(i,k) =   uf(i,k)*Hpr_new(i)
            vv(i,k) =   vf(i,k)*Hpr_new(i)
            ww(i,k) = omef(i,k)*Hpr_new(i)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication3Dtype2(uu)
      call communication3Dtype2(vv)
      call communication3Dtype2(ww)
#     endif
!     =============== END ================
!     ====================================
     
!      _______________________________________________
!     |                                               |
!     |  2.4)  Right-hand side: rhsu, rhsv, rhsw      |
!     |_______________________________________________|

      rhsu = 0.0d0
      rhsv = 0.0d0
      rhsw = 0.0d0 
!      ____________________________________
!     |                                    |
!     |   Dynamic pressure contribution    |
!     |____________________________________|

#     ifdef DynamicPressure
!        ----.----.----.----.----.----.----.----.----.
!        Decoupling Method 1:
#        ifdef Key_ProjectionFirst
            rhsu = 0.0d0
            rhsv = 0.0d0
            rhsw = 0.0d0                    
#        endif
!        ----.----.----.----.----.----.----.----.----.
!        Decoupling Method 2:
#        ifdef Key_ProjectionSecond
            call grandientLSM(dpdx,dpdy,dpds,pf,No_cp,nbe,sig)
            do k=1,NZ
               do i=1,N_CELL  
                  Sdpdx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
                  Sdpdy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
                  Sdpdz(i,k) = sigmaz(i,k)*dpds(i,k)
                  rhsu(i,k) = -(1.0d0/rho_f)*Sdpdx(i,k)*Hpr_new(i)
                  rhsv(i,k) = -(1.0d0/rho_f)*Sdpdy(i,k)*Hpr_new(i)
                  rhsw(i,k) = -(1.0d0/rho_f)*Sdpdz(i,k)*Hpr_new(i) 
               enddo
            enddo
#        endif
!        ----.----.----.----.----.----.----.----.----.
!        Decoupling Method 3:
#        ifdef Key_PredictorCorrector
            call grandientLSM(dpdx,dpdy,dpds,pf,No_cp,nbe,sig)
            do k=1,NZ
               do i=1,N_CELL  
                  Sdpdx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
                  Sdpdy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
                  Sdpdz(i,k) = sigmaz(i,k)*dpds(i,k)
                  rhsu(i,k) = -(1.0d0/rho_f)*Sdpdx(i,k)*Hpr_new(i)
                  rhsv(i,k) = -(1.0d0/rho_f)*Sdpdy(i,k)*Hpr_new(i)
                  rhsw(i,k) = -(1.0d0/rho_f)*Sdpdz(i,k)*Hpr_new(i) 
               enddo
            enddo
#        endif                
#     endif

!      ____________________________________
!     |                                    |
!     |  Hydrostatic pressure contribution |
!     |____________________________________|

#     ifndef KeyFixedFreeSurface
         call grandientLSM2D(detadx_gam,detady_gam,eta_gam,No_cp,nbe)
!        -------------------
!        Using Froude number
#        ifdef KeyFroudeEta
            do i=1,N_CELL
               Shu(i) = -1.0d0/(Fr*Fr)*detadx_gam(i)*Hpr_gam(i)
               Shv(i) = -1.0d0/(Fr*Fr)*detady_gam(i)*Hpr_gam(i)
               Shw(i) = 0.0d0
            enddo
!        -------------------
!        Using gravity (default)
#        else
            do i=1,N_CELL
               Shu(i) = -gra*detadx_gam(i)*Hpr_gam(i)
               Shv(i) = -gra*detady_gam(i)*Hpr_gam(i)
               Shw(i) = 0.0d0
            enddo
#        endif
         do i=1,N_CELL
            do k=1,NZ     
               rhsu(i,k) = rhsu(i,k) + Shu(i)
               rhsv(i,k) = rhsv(i,k) + Shv(i)
               rhsw(i,k) = rhsw(i,k) + Shw(i)  
            enddo
         enddo          
#     endif

!      ____________________________________
!     |                                    |
!     |     Body force (non-dimensional)   |
!     |____________________________________|

#     ifdef KeyBodyForce
         do i=1,N_CELL
            do k=1,NZ     
               rhsu(i,k) = rhsu(i,k) + forcex
               rhsv(i,k) = rhsv(i,k) + forcey
               rhsw(i,k) = rhsw(i,k) + forcez 
            enddo
         enddo          
#     endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(rhsu)
         call communication3Dtype2(rhsv)
         call communication3Dtype2(rhsw)
#     endif
!     =============== END ================
!     ==================================== 

!      _______________________________________________
!     |                                               |
!     |  2.5)  Advection-diffusion equation (*)       |
!     |_______________________________________________|

!      ____________________________________
!     |                                    |
!     |  Velocity comp. wf at F-S & Bottom |
!     |____________________________________|

!     --------------
!     Cell-centered
      do i=1,N_CELL
         uB = uf(i,1)
         vB = vf(i,1)
         wfB(i) = -uB*dhdx(i)-vB*dhdy(i)
         uT = uf(i,NZ)
         vT = vf(i,NZ)
         wfT(i) = detadt(i)+uT*detadx(i)+vT*detady(i)
      enddo
!     --------------
!     Vertex
      do nv=1,N_VERT
         uB = ufv(nv,1)
         vB = vfv(nv,1)  
         wfvB(nv) = -uB*dhdxv(nv)-vB*dhdyv(nv)
         uT = ufv(nv,NZ-1)
         vT = vfv(nv,NZ-1)
         wfvT(nv) = detadtv(nv)+uT*detadxv(nv)+vT*detadyv(nv)             
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(wfB)
         call communication2D(wfT)
#     endif
!     =============== END ================
!     ====================================

!      ____________________________________
!     |                                    |
!     |  Update: uf,vf,wf^(*)              |
!     |____________________________________|

!     ----.----.----.----.----.----.----.----.----.
!     Decoupling Method 1 & 2:
#     if defined(Key_ProjectionFirst) || defined(Key_ProjectionSecond) 
!       ------ Hustar,ustar >>>
        call AdvDiffVelocity_NS(Hufstar,Huf,Hufv,                    &
                                ufstar,uf,ufv,                       &
                                Hpr_new,eta_new,                     &
                                Hprv,etav,                           &
                                h,hv,                                &
                                rhsu,Gamux,Gamuy,Gamuz,uu,vv,ww,pf,  &
                                xc,yc,sig,dsig,No_cp,nbe,            &
                                xv,yv,sigv,dsigv,No_vp,nbev,         &
                                1)                            
!       ------ Hvstar,vstar >>> 
        call AdvDiffVelocity_NS(Hvfstar,Hvf,Hvfv,                    &
                                vfstar,vf,vfv,                       &
                                Hpr_new,eta_new,                     &
                                Hprv,etav,                           &
                                h,hv,                                &
                                rhsv,Gamvx,Gamvy,Gamvz,uu,vv,ww,pf,  & 
                                xc,yc,sig,dsig,No_cp,nbe,            &
                                xv,yv,sigv,dsigv,No_vp,nbev,         &
                                2)               
!       ------ Hwstar,wstar >>> 
        call AdvDiffVelocity_NS(Hwfstar,Hwf,Hwfv,                    &
                                wfstar,wf,wfv,                       &
                                Hpr_new,eta_new,                     &
                                Hprv,etav,                           &
                                h,hv,                                &
                                rhsw,Gamwx,Gamwy,Gamwz,uu,vv,ww,pf,  &
                                xc,yc,sig,dsig,No_cp,nbe,            &
                                xv,yv,sigv,dsigv,No_vp,nbev,         &
                                3)
#     endif
!     ----.----.----.----.----.----.----.----.----.
!     Decoupling Method 3:
#     ifdef Key_PredictorCorrector
          Hufstar = Huf
          Hvfstar = Hvf
          Hwfstar = Hwf
          ufstar  = uf
          vfstar  = vf
          wfstar  = wf
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Hufstar)
            call communication3Dtype2(Hvfstar)
            call communication3Dtype2(Hwfstar)
            call communication3Dtype2(ufstar)
            call communication3Dtype2(vfstar)
            call communication3Dtype2(wfstar)
#        endif	
!        =============== END ================    
!        ====================================
#     endif

!      _______________________________________________
!     |                                               |
!     |  2.6) dw/dt at the bottom (for pressure BC)   |
!     |_______________________________________________|

      do nv=1,N_VERT 
         dw1dtv(nv) = (wfv(nv,1)-w1vn(nv))/dt
      enddo
      do i=1,N_CELL0 
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)       
         dw1dt(i) = (dw1dtv(jv1)+dw1dtv(jv2)+dw1dtv(jv3))/3.0d0                                          
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(dw1dt)
#     endif
!     =============== END ================
!     ====================================
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |       ----  3.H) HYDROSTATIC PRESSURE STEP   ----      |
!     |              (Only hydrostatic pressure)               |
!     |                    °  pf^(new)                         |
!     |________________________________________________________|
!      ********************************************************

#     ifndef DynamicPressure
!        ________________________________________________
!        Velocity
         Hufnp = Hufstar
         ufnp  = ufstar 
!        -------         
         Hvfnp = Hvfstar
         vfnp  = vfstar
!        -------
         Hwfnp = Hwfstar
         wfnp  = wfstar

!        ________________________________________________
!        Hydrostatic pressure (cell-center)
         do k=1,NZ 
            do i=1,N_CELL
               pfnp(i,k) = pa + gra*rho_f*(1.d0-sig(k))*Hpr_gam(i)
            enddo
         enddo
         
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(pfnp)
#        endif	
!        =============== END ================    
!        ====================================

!        ________________________________________________
!        Vertex
         call interpolation3D(pfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              pfnp,xc,yc,sig,dsig,No_cp,nbe)
        
#     else                              
!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |             ----  3.D) PRESSURE STEP   ----            |
!     |             (hydrostatic + dynamic pressure)           |
!     |                    °  pf^(new)                         |
!     |________________________________________________________|
!      ********************************************************

11111  continue

!      _______________________________________________
!     |                                               |
!     |  3.1)  Coefficients: diffusive coeff = Gamp   |
!     |_______________________________________________|

      do k=1,NZ
         do i=1,N_CELL  
            Gampx(i,k) = 1.0d0*Hpr_new(i)
            Gampy(i,k) = 1.0d0*Hpr_new(i)
            Gampz(i,k) = 1.0d0*Hpr_new(i)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(Gampx)
         call communication3Dtype2(Gampy)
         call communication3Dtype2(Gampz)
#     endif	
!     =============== END ================
!     ====================================

!      _______________________________________________
!     |                                               |
!     |  3.2) Right-hand side  = rhsp                 |
!     |_______________________________________________|
      
!     ----.----.----.----.----.----.----.----.----.----
!     Decoupling Method 1:
#     ifdef Key_ProjectionFirst 
         call Calculate_rhsP1(rhsp,                         &
!                            -------------------------------
                             Hpr_new,Hpr,Hprv,              &
                             dHprdt,                        &
!                            -------------------------------
                             Hufstar,Hvfstar,Hwfstar,       &
                             ufstar,vfstar,wfstar,          &
!                            -------------------------------
                             Huf,Hvf,Hwf,                   &
                             uf,vf,wf,pf,                   &                             
!                            -------------------------------  
                             h,xc,yc,sig,dsig,No_cp,nbe,    &
                             hv,xv,yv,sigv,dsigv,No_vp,nbev,&
                             gamma)                       
#     endif 
!     ----.----.----.----.----.----.----.----.----.----
!     Decoupling Method 2:
#     ifdef Key_ProjectionSecond
         call Calculate_rhsP2(rhsp,                         &
!                            -------------------------------
                             Hpr_new,Hpr,Hprv,              &
                             dHprdt,                        &
!                            -------------------------------
                             Hufstar,Hvfstar,Hwfstar,       &
                             ufstar,vfstar,wfstar,          &
!                            -------------------------------
                             Huf,Hvf,Hwf,                   &
                             uf,vf,wf,pf,                   &                             
!                            -------------------------------  
                             h,xc,yc,sig,dsig,No_cp,nbe,    &
                             hv,xv,yv,sigv,dsigv,No_vp,nbev,&
                             gamma)                       
#     endif 
!     ----.----.----.----.----.----.----.----.----.----
!     Decoupling Method 3:
#     ifdef Key_PredictorCorrector
        call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig)
        call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig)
        call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig)
         do k=1,NZ
            do i=1,N_CELL
               Sdudx(i,k) = sigmax(i,k)*duds(i,k) + dudx(i,k)
               Sdvdx(i,k) = sigmax(i,k)*dvds(i,k) + dvdx(i,k)
               Sdwdx(i,k) = sigmax(i,k)*dwds(i,k) + dwdx(i,k)
               Sdudy(i,k) = sigmay(i,k)*duds(i,k) + dudy(i,k)
               Sdvdy(i,k) = sigmay(i,k)*dvds(i,k) + dvdy(i,k)
               Sdwdy(i,k) = sigmay(i,k)*dwds(i,k) + dwdy(i,k)
               Sdudz(i,k) = sigmaz(i,k)*duds(i,k)
               Sdvdz(i,k) = sigmaz(i,k)*dvds(i,k)
               Sdwdz(i,k) = sigmaz(i,k)*dwds(i,k)
               !------------------------------------
               rhsp(i,k) = -2.0d0*rho_f*Hpr_new(i)*VolPrism(i,k)*     &
                         (Sdudy(i,k)*Sdvdx(i,k)-Sdudx(i,k)*Sdvdy(i,k) &
                        + Sdudz(i,k)*Sdwdx(i,k)-Sdudx(i,k)*Sdwdz(i,k) &
                        + Sdvdz(i,k)*Sdwdy(i,k)-Sdvdy(i,k)*Sdwdz(i,k))
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(rhsp)
#        endif	
!        =============== END ================
!        ====================================
#     endif
      
!      _______________________________________________
!     |                                               |
!     |  3.3) Poisson solver: D(Gamma*Dp) = rhsp      |
!     |_______________________________________________|

      call Poisson(pfnp,pfv,                      &
                   rhsp,Gampx,Gampy,Gampz,        &
!                  --------------------------------
                   Hpr_new,eta_new,               &
                   Hprv,etav,                     &
!                  --------------------------------
                   h,hv,                          &
!                  -------------------------------- 
                   xc,yc,sig,dsig,No_cp,nbe,      &
                   xv,yv,sigv,dsigv,No_vp,nbev)

!      _______________________________________________
!     |                                               |
!     |  1.4)  Special case: Test Only Poisson Eqn.   |
!     |_______________________________________________|

#     if defined(KeyTestOnlyPoisson)
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if (rang_topo.eq.0) then
         write(*,'(t10,a40)') '      ---> Now, go to the end of Hydro'
         write(*,*) ' '      
         endif
!     ====================================
!     ====================================
#     else
         write(*,'(t10,a40)') '      ---> Now, go to the end of Hydro'
         write(*,*) ' '  
#     endif
!     =============== END ================
!     ====================================
       goto 11112
#     endif

!      ________________________________________________________
!     |********************************************************|
!     |                                                        |
!     |           ---- 4) CORRECTION STEP: fluid ----          |
!     |               (only with dynamic pressure)             | 
!     |                    °  uf^(new)                         | 
!     |                    °  vf^(new)                         |   
!     |                    °  wf^(new)                         |
!     |________________________________________________________|
!      ********************************************************

!      _______________________________________________
!     |                                               |
!     |  4.1.P)Update of the velocity components PROJ.|
!     |_______________________________________________|

!     ----.----.----.----.----.----.----.----.----.----
!     Decoupling Method 1:
#     ifdef Key_ProjectionFirst
           call Corrector_Flow1(Hufnp,ufnp,Hufv,ufv,           &
                                Hvfnp,vfnp,Hvfv,vfv,           &
                                Hwfnp,wfnp,Hwfv,wfv,           &
!                               --------------------------------
                                ufstar,vfstar,wfstar,          &
!                               --------------------------------
                                pfnp,pf,pfv,                   &
!                               --------------------------------
                                Hpr_new,Hpr,Hprv,              &
                                eta_new,eta,etav,              & 
!                               --------------------------------                                
                                dhdx,dhdy,                     &
                                dhdxv,dhdyv,                   &
!                               --------------------------------                                  
                                detadt,detadx,detady,          &                                                          
                                detadtv,detadxv,detadyv,       &                                                               
!                               --------------------------------  
                                h,xc,yc,sig,dsig,No_cp,nbe,    &
                                hv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                               --------------------------------                                
                                gamma)
#     endif 

!     ----.----.----.----.----.----.----.----.----.----
!     Decoupling Method 2:
#     ifdef Key_ProjectionSecond
           call Corrector_Flow2(Hufnp,ufnp,Hufv,ufv,           &
                                Hvfnp,vfnp,Hvfv,vfv,           &
                                Hwfnp,wfnp,Hwfv,wfv,           &
!                               --------------------------------
                                ufstar,vfstar,wfstar,          &
!                               --------------------------------
                                pfnp,pf,pfv,                   &
!                               --------------------------------
                                Hpr_new,Hpr,Hprv,              &
                                eta_new,eta,etav,              & 
!                               --------------------------------                                
                                dhdx,dhdy,                     &
                                dhdxv,dhdyv,                   &
!                               --------------------------------                                  
                                detadt,detadx,detady,          &                                                          
                                detadtv,detadxv,detadyv,       &                                                               
!                               --------------------------------  
                                h,xc,yc,sig,dsig,No_cp,nbe,    &
                                hv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                               --------------------------------                                
                                gamma) 
#     endif                                

!     ----.----.----.----.----.----.----.----.----.----
#     ifdef Key_PredictorCorrector
!        _______________________________________________
!       |                                               |
!       |  4.1.1)  Right-hand side: rhsu, rhsv, rhsw    |
!       |_______________________________________________|

         rhsu = 0.0d0
         rhsv = 0.0d0
         rhsw = 0.0d0 
!        ____________________________________
!       |                                    |
!       |   Dynamic pressure contribution    |
!       |____________________________________|

         call grandientLSM(dpdx,dpdy,dpds,pfnp,No_cp,nbe,sig) 
!        --------- 
         do k=1,NZ
            do i=1,N_CELL  
               Sdpdx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
               Sdpdy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
               Sdpdz(i,k) = sigmaz(i,k)*dpds(i,k)
               rhsu(i,k) = -(1.0d0/rho_f)*Sdpdx(i,k)*Hpr_new(i)
               rhsv(i,k) = -(1.0d0/rho_f)*Sdpdy(i,k)*Hpr_new(i)
               rhsw(i,k) = -(1.0d0/rho_f)*Sdpdz(i,k)*Hpr_new(i)
            enddo
         enddo
!        ____________________________________
!       |                                    |
!       |  Hydrostatic pressure contribution |
!       |____________________________________|

#        ifndef KeyFixedFreeSurface
            do i=1,N_CELL
               do k=1,NZ     
                  rhsu(i,k) = rhsu(i,k) + Shu(i)
                  rhsv(i,k) = rhsv(i,k) + Shv(i)
                  rhsw(i,k) = rhsw(i,k) + Shw(i)  
               enddo
            enddo          
#        endif
!        ____________________________________
!       |                                    |
!       |     Body force (non-dimensional)   |
!       |____________________________________|

#        ifdef KeyBodyForce
            do i=1,N_CELL
               do k=1,NZ     
                  rhsu(i,k) = rhsu(i,k) + forcex
                  rhsv(i,k) = rhsv(i,k) + forcey
                  rhsw(i,k) = rhsw(i,k) + forcez 
               enddo
            enddo          
#        endif

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(rhsu)
            call communication3Dtype2(rhsv)
            call communication3Dtype2(rhsw)
#        endif
!        =============== END ================
!        ==================================== 

!        _______________________________________________
!       |                                               |
!       |  4.1.2)  Advection-diffusion equation (n+1)   |
!       |_______________________________________________|

!        --------- Hufnp,Huf >>> 
         call AdvDiffVelocity_NS(Hufnp,Huf,Hufv,                      &
                                 ufnp,uf,ufv,                         &
                                 Hpr_new,eta_new,                     &
                                 Hprv,etav,                           &
                                 h,hv,                                &
                                 rhsu,Gamux,Gamuy,Gamuz,uu,vv,ww,pfnp,&
                                 xc,yc,sig,dsig,No_cp,nbe,            &
                                 xv,yv,sigv,dsigv,No_vp,nbev,         &
                                 1)                              
!        --------- Hvfnp,Hvf >>>
         call AdvDiffVelocity_NS(Hvfnp,Hvf,Hvfv,                      &
                                 vfnp,vf,vfv,                         &
                                 Hpr_new,eta_new,                     &
                                 Hprv,etav,                           &
                                 h,hv,                                &
                                 rhsv,Gamvx,Gamvy,Gamvz,uu,vv,ww,pfnp,&
                                 xc,yc,sig,dsig,No_cp,nbe,            &
                                 xv,yv,sigv,dsigv,No_vp,nbev,         &
                                 2)
!        --------- Hwfnp,Hwf >>>      
         call AdvDiffVelocity_NS(Hwfnp,Hwf,Hwfv,                      &
                                 wfnp,wf,wfv,                         &
                                 Hpr_new,eta_new,                     &
                                 Hprv,etav,                           &
                                 h,hv,                                &
                                 rhsw,Gamwx,Gamwy,Gamwz,uu,vv,ww,pfnp,&
                                 xc,yc,sig,dsig,No_cp,nbe,            &
                                 xv,yv,sigv,dsigv,No_vp,nbev,         &
                                 3)
#     endif
      
#     endif /* <<< End of dynamic pressure >>> */                                                            

!*********************************************************************!
!---------------------------------------------------------------------!
!             Closing: Iterative loop for free-surface update         !
!---------------------------------------------------------------------!
!*********************************************************************!

#     ifndef KeyFixedFreeSurface 

!      _______________________________________________
!     |                                               |
!     |  [FS.C.1]  Error (L2-norm)                    |
!     |_______________________________________________|

      epsp = 0.0d0
      do i=1,N_CELL0
         epsp = epsp + (dabs(eta_old(i)-eta_new(i))**2)
      enddo 

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call SUM_parallel(epsp,SUMepsp)
         epsp = SUMepsp
#     endif	
!     =============== END ================    
!     ==================================== 
     
      epsp = dsqrt(epsp)/float(N_CELL0global)
      
!      _______________________________________________
!     |                                               |
!     |  [FS.C.2]  Convergence criteria               |
!     |_______________________________________________|
      
      if ((epsp.gt.tolConv).and.(iconv.lt.nconv)) goto 9999

!      _______________________________________________
!     |                                               |
!     |  [FS.C.3]  Displaying                         |
!     |_______________________________________________|

112   format(t2,'Iterations=',I3,', Error=',e23.15)

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if(rang_topo.eq.0) then
         if (iconv.gt.1) write(*,112) iconv,epsp
         endif
#     else
         if (iconv.gt.1) write(*,112) iconv,epsp
#     endif
!     =============== END ================    
!     ====================================

!      _______________________________________________
!     |                                               |
!     |  [FS.C.4]  Final FS results                   |
!     |_______________________________________________|

      Hprnp = Hpr_new
      etanp = eta_new

#     endif

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |         Calculate Vorticity (vector & magnitud)        |
!     |________________________________________________________|

      call grandientLSM(dudx,dudy,duds,ufnp,No_cp,nbe,sig)
      call grandientLSM(dvdx,dvdy,dvds,vfnp,No_cp,nbe,sig)
      call grandientLSM(dwdx,dwdy,dwds,wfnp,No_cp,nbe,sig)
!     ----------------------------- 
      do k=1,NZ
         do i=1,N_CELL  
            Sdudy(i,k) = sigmay(i,k)*duds(i,k) + dudy(i,k)
            Sdudz(i,k) = sigmaz(i,k)*duds(i,k)
            Sdvdx(i,k) = sigmax(i,k)*dvds(i,k) + dvdx(i,k)
            Sdvdz(i,k) = sigmaz(i,k)*dvds(i,k)
            Sdwdx(i,k) = sigmax(i,k)*dwds(i,k) + dwdx(i,k)
            Sdwdy(i,k) = sigmay(i,k)*dwds(i,k) + dwdy(i,k)
            Vortx(i,k) = Sdwdy(i,k)-Sdvdz(i,k)
            Vorty(i,k) = Sdudz(i,k)-Sdwdx(i,k)
            Vortz(i,k) = Sdvdx(i,k)-Sdudy(i,k)
            Vort(i,k) = dsqrt((Sdwdy(i,k)-Sdvdz(i,k))**2  &
                             +(Sdudz(i,k)-Sdwdx(i,k))**2  &
                             +(Sdvdx(i,k)-Sdudy(i,k))**2)                                    
         enddo
      enddo
!     ____________________________
!     Vertex interpolation
      call interpolation3D(Vortxv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vortx,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Vortyv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vorty,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Vortzv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vortz,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Vortv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vort,xc,yc,sig,dsig,No_cp,nbe)                 
        
!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

11112 continue

      deallocate(Hufstar,Hvfstar,Hwfstar, &
                 ufstar,vfstar,wfstar,    &
                 Dpfnp,Dpf,Dpfnpv,Dpfv,   &
!                --------------------------
                 dudx,dudy,duds,          &
                 dvdx,dvdy,dvds,          &
                 dwdx,dwdy,dwds,          &
                 dpdx,dpdy,dpds,          &
                 dDpdx,dDpdy,dDpds,       &
!                --------------------------
                 Gamux,Gamuy,Gamuz,       &
                 Gamvx,Gamvy,Gamvz,       &
                 Gamwx,Gamwy,Gamwz,       &
                 Gampx,Gampy,Gampz,       &
!                --------------------------
                 rhsu,rhsv,rhsw,rhsp,     &
!                --------------------------
                 Shu,Shv,Shw,             &
                 dHprdx,dHprdy,dHprdt,    &
                 detadx,detady,detadt,    &
                 dhdx,dhdy,               &
!                --------------------------
                 dHprdxv,dHprdyv,dHprdtv, &
                 detadxv,detadyv,detadtv, &
                 dhdxv,dhdyv,             & 
!                --------------------------
                 gradetadx,dgraxdx,dgraydx,&
                 gradetady,dgraxdy,dgraydy)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: The N-S equations'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!---------------------------------------------------------------------!
!             AUXILIAR 1: RHS OF THE POISSON PRESSURE EQUATION        !
!                         (Projection 1st order)                      !
!                             April 2017                              !
!---------------------------------------------------------------------!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      SUBROUTINE Calculate_rhsP1(rhsp,                         &
!                               --------------------------------
                                Hpr_new,Hpr,Hprv,              &
                                dHprdt,                        &
!                               --------------------------------
                                Hufstar,Hvfstar,Hwfstar,       &
                                ufstar,vfstar,wfstar,          &
!                               --------------------------------
                                Hufn,Hvfn,Hwfn,                &
                                ufn,vfn,wfn,pfn,               &                                
!                               --------------------------------  
                                h,xc,yc,sig,dsig,No_cp,nbe,    &
                                hv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                               --------------------------------                               
                                gamma)   
                                           
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the right-hand side term corresponding   !
!    to the projection method.                                        !
!                                                                     !                                                              !
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

      real*8, dimension(:,:):: rhsp(N_CELL,NZ)
!     --------------------------------------      
      real*8, dimension(:)  :: Hpr_new(N_CELL)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT) 
      real*8, dimension(:)  :: dHprdt(N_CELL)
!     --------------------------------------
      real*8, dimension(:,:):: Hufstar(N_CELL,NZ)
      real*8, dimension(:,:):: Hvfstar(N_CELL,NZ)
      real*8, dimension(:,:):: Hwfstar(N_CELL,NZ) 
      real*8, dimension(:,:):: ufstar(N_CELL,NZ)
      real*8, dimension(:,:):: vfstar(N_CELL,NZ)
      real*8, dimension(:,:):: wfstar(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:):: Hufn(N_CELL,NZ)
      real*8, dimension(:,:):: Hvfn(N_CELL,NZ)
      real*8, dimension(:,:):: Hwfn(N_CELL,NZ) 
      real*8, dimension(:,:):: ufn(N_CELL,NZ)
      real*8, dimension(:,:):: vfn(N_CELL,NZ)
      real*8, dimension(:,:):: wfn(N_CELL,NZ)      
      real*8, dimension(:,:):: pfn(N_CELL,NZ)          
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8, dimension(:,:) :: omegaf(N_CELL,NZ)
      real*8, dimension(:,:) :: Homegaf(N_CELL,NZ)
!     ----------------------------------------      
      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0)
      real*8,dimension(:) :: AmG(N_CELL0)
!     ------------------------------------------
      real*8,dimension(:) :: dHprdx(N_CELL),dHprdy(N_CELL)
      real*8,dimension(:) :: detadx(N_CELL),detady(N_CELL)
      real*8,dimension(:) :: dHprdx_new(N_CELL),dHprdy_new(N_CELL)
      real*8,dimension(:) :: detadx_new(N_CELL),detady_new(N_CELL)
!     ------------------------------------------
      real*8, dimension(:,:) :: dudx(N_CELL,NZ)
      real*8, dimension(:,:) :: dudy(N_CELL,NZ)
      real*8, dimension(:,:) :: duds(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dvds(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dwds(N_CELL,NZ)
!     ------------------------------------------
      real*8, dimension(:,:) :: dudxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dudyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dudsn(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdsn(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdsn(N_CELL,NZ)      
!     ------------------------------------------      
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ)
!     ------------------------------------------             
      real*8, dimension(:,:) :: dHudx(N_CELL,NZ)
      real*8, dimension(:,:) :: dHvdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dHwds(N_CELL,NZ)      
!     ------------------------------------------      
      real*8 :: HH,HHi,HHj,dHHx,dHhy
      real*8 :: Ue,Uijpos,Uijneg,GFi,GFj,UiT,UiB
      real*8 :: x,y,z,nxArea,nyArea,u0j,v0j
      real*8 :: dtoVol,const,Vol,som,somu,somv,ct,Are
      integer:: it,jc1,jc2,jc3,jc,elem
!     ----------------------------------------
      real*8 :: sumqux,sumquy,sumqvx,sumqvy,deter
      real*8 :: dqudx,dqvdy,som1,som2,som3
      real*8 :: errorsys,residu,SUMerrorsys,FS_funeta
!     --------------------------------------
      real*8 :: Sdudx,Sdvdy,Sdwdz      
      real*8 :: Sdudxn,Sdvdyn,Sdwdzn      
      real*8, dimension(:,:) :: Tpx(N_CELL,NZ)
      real*8, dimension(:,:) :: Tpy(N_CELL,NZ)
      real*8, dimension(:,:) :: Tpz(N_CELL,NZ)
      real*8, dimension(:,:) :: dTpxdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dTpydy(N_CELL,NZ)
      real*8, dimension(:,:) :: dTpzds(N_CELL,NZ)
!     ------------------------------------------       
      real*8, dimension(N_CELL,NZ) :: Ui1,Ui2,Ui3,Ui4,Ui5
      real*8 ,dimension(1:3) :: Uij                 
!     ------------------------------------------
      real*8,dimension(N_CELL,NZ) :: aux 
      real*8 :: c1,c2,Du,Dun,Lpn
      integer :: ind1,ind2  
!     ------------------------------------------
!     ____________________________________
!    |                                    |
!    |   Declaration of parameters        |
!    |____________________________________|
             
#     ifdef Key_ProjectionFirst
      integer, parameter :: OptionRHSP = 4 ! <=4 FVM(1st)
#     else
      integer, parameter :: OptionRHSP = 6 ! Always=6
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: Calculate_rhsP1'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     _______________________________________________
!     Update (H*omega)^(*)

      do k=1,NZ
         do i=1,N_CELL 
            omegaf(i,k) = sigmat(i,k)                 & 
                         + ufstar(i,k)*sigmax(i,k)    &
                         + vfstar(i,k)*sigmay(i,k)    &
                         + wfstar(i,k)*sigmaz(i,k)

            Homegaf(i,k) = Hpr_new(i)*(sigmat(i,k)    & 
                         + ufstar(i,k)*sigmax(i,k)    &
                         + vfstar(i,k)*sigmay(i,k)    &
                         + wfstar(i,k)*sigmaz(i,k))
         enddo
      enddo 
 
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication3Dtype2(omegaf)
      call communication3Dtype2(Homegaf)
#     endif
!     =============== END ================
!     ====================================

      rhsp = 0.0d0
         
!      _______________________________________________
!     |                                               |
!     |  OptionRHSP = 1:                              |
!     |  By Least square technique: Ga et al. 2012    |
!     |_______________________________________________|

      IF (OptionRHSP.eq.1) THEN
        call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig)
        call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig)
        call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig)        
         do k=2,NZ-1
            do i=1,N_CELL0 
               Sdudx = sigmax(i,k)*duds(i,k) + dudx(i,k)
               Sdvdy = sigmay(i,k)*dvds(i,k) + dvdy(i,k)
               Sdwdz = sigmaz(i,k)*dwds(i,k)
               rhsp(i,k) = Hpr(i)*VolPrism(i,k)*rho_f/dt* &
                          (Sdudx+Sdvdy+Sdwdz)
            enddo
         enddo
!      _______________________________________________
!     |                                               |
!     |  OptionRHSP = 2:                              |
!     |  By Least square technique: Ponderate gamma   |
!     |_______________________________________________|
                            
      ELSEIF (OptionRHSP.eq.2) THEN
        !------------      
         c1 =  rho_f/(gamma*gamma*dt)
         call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig)
         call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig)
         call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig)
         call grandientLSM(dudxn,dudyn,dudsn,ufn,No_cp,nbe,sig)
         call grandientLSM(dvdxn,dvdyn,dvdsn,vfn,No_cp,nbe,sig)
         call grandientLSM(dwdxn,dwdyn,dwdsn,wfn,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0 
               Sdudx = sigmax(i,k)*duds(i,k) + dudx(i,k)
               Sdvdy = sigmay(i,k)*dvds(i,k) + dvdy(i,k)
               Sdwdz = sigmaz(i,k)*dwds(i,k)
               Sdudxn = sigmax(i,k)*dudsn(i,k) + dudxn(i,k)
               Sdvdyn = sigmay(i,k)*dvdsn(i,k) + dvdyn(i,k)
               Sdwdzn = sigmaz(i,k)*dwdsn(i,k)
               Du  = Sdudx+Sdvdy+Sdwdz
               Dun = Sdudxn+Sdvdyn+Sdwdzn
               rhsp(i,k) = c1*(gamma*Du + (1.0d0-gamma)*Dun) 
            enddo
         enddo
        !------------
         c2 = -rho_f*(1.0d0-gamma)/gamma                 
         call grandientLSM(dpdx,dpdy,dpds,pfn,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0                
               Tpx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
               Tpy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
               Tpz(i,k) = sigmax(i,k)*Tpx(i,k) &
                         +sigmay(i,k)*Tpy(i,k) &
                         +sigmaz(i,k)*(sigmaz(i,k)*dpds(i,k))
            enddo
         enddo         
         call grandientLSM(dTpxdx,aux,aux,Tpx,No_cp,nbe,sig)
         call grandientLSM(aux,dTpydy,aux,Tpy,No_cp,nbe,sig)
         call grandientLSM(aux,aux,dTpzds,Tpz,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0
               Lpn = dTpxdx(i,k)+dTpydy(i,k)+dTpzds(i,k)     
               rhsp(i,k) = rhsp(i,k) + c2*Lpn      
            enddo
         enddo
        !------------
         do k=2,NZ-1
            do i=1,N_CELL0                   
               rhsp(i,k) = Hpr(i)*VolPrism(i,k)*rhsp(i,k) 
            enddo
         enddo 
!      _______________________________________________
!     |                                               |
!     |  OptionRHSP = 3:                              |
!     |  By Least square technique: D*(Hu,Hv,Hw)      |
!     |_______________________________________________|                               

      ELSEIF (OptionRHSP.eq.3) THEN
         call grandientLSM(dHudx,dudy,duds,Hufstar,No_cp,nbe,sig)
         call grandientLSM(dvdx,dHvdy,dvds,Hvfstar,No_cp,nbe,sig)
         call grandientLSM(dwdx,dwdy,dHwds,Homegaf,No_cp,nbe,sig)        
         do k=2,NZ-1
            do i=1,N_CELL0
               rhsp(i,k) = VolPrism(i,k)*rho_f/dt*(dHprdt(i)   &
                                                  +dHudx(i,k)  &
                                                  +dHvdy(i,k)  &
                                                  +dHwds(i,k))         
            enddo
         enddo
!      _______________________________________________
!     |                                               |
!     |  OptionRHSP = 4:                              |
!     |  By first order Finite Volume Method (FVM)    |
!     |_______________________________________________|
    
      ELSEIF (OptionRHSP.eq.4) THEN
         Ui1 = 0.0d0    
         Ui2 = 0.0d0
         Ui3 = 0.0d0
         Ui4 = 0.0d0
         Ui5 = 0.0d0        
         DO k=2,NZ-1
            do i=1,N_CELL0
               do j=1,3
                  jc = No_cp(i,j)
                  u0j = 0.5d0*(Hufstar(i,k)+Hufstar(jc,k)) 
                  v0j = 0.5d0*(Hvfstar(i,k)+Hvfstar(jc,k))                
                  nxArea =  dyVV(i,j)*dsigv(k-1)
                  nyArea = -dxVV(i,j)*dsigv(k-1)
                  Uij(j) = u0j*nxArea+v0j*nyArea
               enddo
               Ui1(i,k) = Uij(1)    
               Ui2(i,k) = Uij(2)
               Ui3(i,k) = Uij(3)
               Ui4(i,k) = 0.5d0*(Homegaf(i,k)+Homegaf(i,k+1))*areaCell(i)
               Ui5(i,k) =-0.5d0*(Homegaf(i,k)+Homegaf(i,k-1))*areaCell(i)
            enddo
         ENDDO
         DO k=2,NZ-1  
            do i=1,N_CELL0         
               som = 0.0d0
               !-------------- 
               Uij(1) = Ui1(i,k)  
               Uij(2) = Ui2(i,k)
               Uij(3) = Ui3(i,k)              
               do j=1,3
                  jc = No_cp(i,j)
                  Uijpos =  0.5d0*(Uij(j)+abs(Uij(j)))
                  Uijneg = -0.5d0*(Uij(j)-abs(Uij(j)))               	
                  som    =  som + Uijpos - Uijneg
               enddo
               !--------------
               Uijpos =  0.5d0*(Ui4(i,k)+abs(Ui4(i,k)))
               Uijneg = -0.5d0*(Ui4(i,k)-abs(Ui4(i,k)))               	
               som    =  som + Uijpos - Uijneg 
               !--------------	           
               Uijpos =  0.5d0*(Ui5(i,k)+abs(Ui5(i,k)))
               Uijneg = -0.5d0*(Ui5(i,k)-abs(Ui5(i,k)))
               som    =  som + Uijpos - Uijneg
               !--------------
               rhsp(i,k) = (rho_f/dt)*(VolPrism(i,k)*dHprdt(i)+som)	           	                                                        
            enddo
         ENDDO
!      _______________________________________________
!     |                                               |
!     |  OptionRHSP = 5:                              |
!     |  By second order Finite Volume Method (FVM)   |
!     |_______________________________________________|                         
    
      ELSEIF (OptionRHSP.eq.5) THEN
         call grandientLSM2D(dHprdx,dHprdy,Hpr,No_cp,nbe)
         call grandientLSM2D(dHprdx_new,dHprdy_new,Hpr_new,No_cp,nbe)       
         DO k=2,NZ-1  
            do i=1,N_CELL0
                som = 0.0d0
               !--------------------------------
               ! Horizontal
               do j=1,3
                  jc = No_cp(i,j)
                  u0j = 0.5d0*(ufstar(i,k)+ufstar(jc,k)) 
                  v0j = 0.5d0*(vfstar(i,k)+vfstar(jc,k))                
                  nxArea =  dyVV(i,j)*dsigv(k-1)
                  nyArea = -dxVV(i,j)*dsigv(k-1)
                  Ue = u0j*nxArea+v0j*nyArea
                  Uijpos = 0.5d0*(Ue+abs(Ue))
                  Uijneg = 0.5d0*(Ue-abs(Ue))
                  !--------
                  dHHx = gamma*dHprdx_new(i)+(1.0d0-gamma)*dHprdx(i)
                  dHHy = gamma*dHprdy_new(i)+(1.0d0-gamma)*dHprdy(i)  
                  GFi =  dHHx*(xme(i,j)-xc(i))+dHHy*(yme(i,j)-yc(i))
                  !--------
                  dHHx = gamma*dHprdx_new(jc)+(1.0d0-gamma)*dHprdx(jc)
                  dHHy = gamma*dHprdy_new(jc)+(1.0d0-gamma)*dHprdy(jc) 
                  GFj =  dHHx*(xme(i,j)-xc(jc))+dHHy*(yme(i,j)-yc(jc))
                  !--------
                  HHi = gamma*Hpr_new(i) +(1.0d0-gamma)*Hpr(i)
                  HHj = gamma*Hpr_new(jc)+(1.0d0-gamma)*Hpr(jc)
                  som = som + Uijpos*(HHi+GFi) + Uijneg*(HHj+GFj)
               enddo
               !--------------------------------
               ! Top
               UiT = 0.5d0*(omegaf(i,k)+omegaf(i,k+1))*areaCell(i)
               HHi = gamma*Hpr_new(i)+(1.0d0-gamma)*Hpr(i)
               som = som + UiT*HHi
               !--------------------------------
               ! Bottom
               UiB = -0.5d0*(omegaf(i,k)+omegaf(i,k-1))*areaCell(i)
               HHi = gamma*Hpr_new(i)+(1.0d0-gamma)*Hpr(i)
               som = som + UiB*HHi
               !--------------
               rhsp(i,k) = Hpr(i)*(rho_f/dt)*(VolPrism(i,k)*dHprdt(i)+som)
            enddo
         ENDDO
      ENDIF

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(rhsp)
#     endif
!     =============== END ================
!     ==================================== 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: Calculate_rhsP1'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!---------------------------------------------------------------------!
!             AUXILIAR 2: RHS OF THE POISSON PRESSURE EQUATION        !
!                         (Projection 2nd order)                      !
!                             Nov 2017                                !
!---------------------------------------------------------------------!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      SUBROUTINE Calculate_rhsP2(rhsp,                         &
!                               --------------------------------
                                Hpr_new,Hpr,Hprv,              &
                                dHprdt,                        &
!                               --------------------------------
                                Hufstar,Hvfstar,Hwfstar,       &
                                ufstar,vfstar,wfstar,          &
!                               --------------------------------
                                Hufn,Hvfn,Hwfn,                &
                                ufn,vfn,wfn,pfn,               &                                
!                               --------------------------------  
                                h,xc,yc,sig,dsig,No_cp,nbe,    &
                                hv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                               --------------------------------                               
                                gamma)   
                                           
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the right-hand side term corresponding   !
!    to the projection method usinf dp=p^new-p^old.                   !
!                                                                     !                                                              !
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

      real*8, dimension(:,:):: rhsp(N_CELL,NZ)
!     --------------------------------------      
      real*8, dimension(:)  :: Hpr_new(N_CELL)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT) 
      real*8, dimension(:)  :: dHprdt(N_CELL)
!     --------------------------------------
      real*8, dimension(:,:):: Hufstar(N_CELL,NZ)
      real*8, dimension(:,:):: Hvfstar(N_CELL,NZ)
      real*8, dimension(:,:):: Hwfstar(N_CELL,NZ) 
      real*8, dimension(:,:):: ufstar(N_CELL,NZ)
      real*8, dimension(:,:):: vfstar(N_CELL,NZ)
      real*8, dimension(:,:):: wfstar(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:):: Hufn(N_CELL,NZ)
      real*8, dimension(:,:):: Hvfn(N_CELL,NZ)
      real*8, dimension(:,:):: Hwfn(N_CELL,NZ) 
      real*8, dimension(:,:):: ufn(N_CELL,NZ)
      real*8, dimension(:,:):: vfn(N_CELL,NZ)
      real*8, dimension(:,:):: wfn(N_CELL,NZ)      
      real*8, dimension(:,:):: pfn(N_CELL,NZ)          
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8, dimension(:,:) :: omegaf(N_CELL,NZ)
      real*8, dimension(:,:) :: Homegaf(N_CELL,NZ)
!     ----------------------------------------      
      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0)
      real*8,dimension(:) :: AmG(N_CELL0)
!     ------------------------------------------
      real*8,dimension(:) :: dHprdx(N_CELL),dHprdy(N_CELL)
      real*8,dimension(:) :: detadx(N_CELL),detady(N_CELL)
      real*8,dimension(:) :: dHprdx_new(N_CELL),dHprdy_new(N_CELL)
      real*8,dimension(:) :: detadx_new(N_CELL),detady_new(N_CELL)
!     ------------------------------------------
      real*8, dimension(:,:) :: rhsp1(N_CELL,NZ)
      real*8, dimension(:,:) :: rhsp2(N_CELL,NZ)
!     ------------------------------------------
      real*8, dimension(:,:) :: dudx(N_CELL,NZ)
      real*8, dimension(:,:) :: dudy(N_CELL,NZ)
      real*8, dimension(:,:) :: duds(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dvds(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dwds(N_CELL,NZ)
!     ------------------------------------------
      real*8, dimension(:,:) :: dudxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dudyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dudsn(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdsn(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdsn(N_CELL,NZ)      
!     ------------------------------------------      
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ)
!     ------------------------------------------             
      real*8, dimension(:,:) :: dHudx(N_CELL,NZ)
      real*8, dimension(:,:) :: dHvdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dHwds(N_CELL,NZ)      
!     ------------------------------------------      
      real*8 :: HH,HHi,HHj,dHHx,dHhy
      real*8 :: Ue,Uijpos,Uijneg,GFi,GFj,UiT,UiB
      real*8 :: x,y,z,nxArea,nyArea,u0j,v0j
      real*8 :: dtoVol,const,Vol,som,somu,somv,ct,Are
      integer:: it,jc1,jc2,jc3,jc,elem
!     ----------------------------------------
      real*8 :: sumqux,sumquy,sumqvx,sumqvy,deter
      real*8 :: dqudx,dqvdy,som1,som2,som3
      real*8 :: errorsys,residu,SUMerrorsys,FS_funeta
!     --------------------------------------
      real*8 :: Sdudx,Sdvdy,Sdwdz      
      real*8 :: Sdudxn,Sdvdyn,Sdwdzn      
      real*8, dimension(:,:) :: Tpx(N_CELL,NZ)
      real*8, dimension(:,:) :: Tpy(N_CELL,NZ)
      real*8, dimension(:,:) :: Tpz(N_CELL,NZ)
      real*8, dimension(:,:) :: dTpxdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dTpydy(N_CELL,NZ)
      real*8, dimension(:,:) :: dTpzds(N_CELL,NZ)
!     ------------------------------------------       
      real*8, dimension(N_CELL,NZ) :: Ui1,Ui2,Ui3,Ui4,Ui5
      real*8 ,dimension(1:3) :: Uij                 
!     ------------------------------------------
      real*8,dimension(N_CELL,NZ) :: aux 
      real*8 :: c1,c2,Du,Dun,Lpn
      integer :: ind1,ind2  
!     ------------------------------------------
!     ____________________________________
!    |                                    |
!    |   Declaration of parameters        |
!    |____________________________________|
             
      integer, parameter :: OptionRHSP2nd = 1 ! Always=2

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: Calculate_rhsP2'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     _______________________________________________
!     Update (H*omega)^(*)

      do k=1,NZ
         do i=1,N_CELL 
            omegaf(i,k) = sigmat(i,k)                 & 
                         + ufstar(i,k)*sigmax(i,k)    &
                         + vfstar(i,k)*sigmay(i,k)    &
                         + wfstar(i,k)*sigmaz(i,k)

            Homegaf(i,k) = Hpr_new(i)*omegaf(i,k)
         enddo
      enddo 
 
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication3Dtype2(omegaf)
      call communication3Dtype2(Homegaf)
#     endif
!     =============== END ================
!     ====================================

      rhsp  = 0.0d0
      rhsp1 = 0.0d0
      rhsp2 = 0.0d0

!      _______________________________________________
!     |                                               |
!     |  OptionRHSP = 1:                              |
!     |  Using dp=p^(n+1)-p^(n) & LSM                 |
!     |_______________________________________________|
                            
      IF (OptionRHSP2nd.eq.1) THEN

        !------------      
         c1 =  rho_f/dt
         call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig)
         call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig)
         call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig)
         do k=2,NZ-1
            do i=1,N_CELL0 
               Sdudx  = sigmax(i,k)*duds(i,k) + dudx(i,k)
               Sdvdy  = sigmay(i,k)*dvds(i,k) + dvdy(i,k)
               Sdwdz  = sigmaz(i,k)*dwds(i,k)
               Du  = Sdudx+Sdvdy+Sdwdz
               rhsp1(i,k) = c1*Du
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(rhsp1)
#        endif
!        =============== END ================
!        ==================================== 
        !------------
         c2 = - 1.0d0                 
         call grandientLSM(dpdx,dpdy,dpds,pfn,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0                
               Tpx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
               Tpy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
               Tpz(i,k) = sigmax(i,k)*Tpx(i,k) &
                         +sigmay(i,k)*Tpy(i,k) &
                         +sigmaz(i,k)*(sigmaz(i,k)*dpds(i,k))
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Tpx)
            call communication3Dtype2(Tpy)
            call communication3Dtype2(Tpz)
#        endif
!        =============== END ================
!        ====================================        
         call grandientLSM(dTpxdx,aux,aux,Tpx,No_cp,nbe,sig)
         call grandientLSM(aux,dTpydy,aux,Tpy,No_cp,nbe,sig)
         call grandientLSM(aux,aux,dTpzds,Tpz,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0
               Lpn = dTpxdx(i,k)+dTpydy(i,k)+dTpzds(i,k)     
               rhsp2(i,k) = c2*Lpn
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(rhsp2)
#        endif
!        =============== END ================
!        ====================================
        !------------
         do k=2,NZ-1
            do i=1,N_CELL0                     
               rhsp(i,k) = Hpr(i)*VolPrism(i,k)*(rhsp1(i,k)+rhsp2(i,k))
            enddo
         enddo
!      _______________________________________________
!     |                                               |
!     |  OptionRHSP2nd = 2:                           |
!     |  Using dp=p^(n+1)-p^(n) & LSM Ponderate gamma |
!     |_______________________________________________|
                            
      ELSEIF (OptionRHSP2nd.eq.2) THEN
        !------------      
         c1 =  rho_f/(gamma*gamma*dt)
         call grandientLSM(dudx,dudy,duds,ufstar,No_cp,nbe,sig)
         call grandientLSM(dvdx,dvdy,dvds,vfstar,No_cp,nbe,sig)
         call grandientLSM(dwdx,dwdy,dwds,wfstar,No_cp,nbe,sig)
         call grandientLSM(dudxn,dudyn,dudsn,ufn,No_cp,nbe,sig)
         call grandientLSM(dvdxn,dvdyn,dvdsn,vfn,No_cp,nbe,sig)
         call grandientLSM(dwdxn,dwdyn,dwdsn,wfn,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0 
               Sdudx  = sigmax(i,k)*duds(i,k) + dudx(i,k)
               Sdvdy  = sigmay(i,k)*dvds(i,k) + dvdy(i,k)
               Sdwdz  = sigmaz(i,k)*dwds(i,k)
               Sdudxn = sigmax(i,k)*dudsn(i,k) + dudxn(i,k)
               Sdvdyn = sigmay(i,k)*dvdsn(i,k) + dvdyn(i,k)
               Sdwdzn = sigmaz(i,k)*dwdsn(i,k)
               Du  = Sdudx+Sdvdy+Sdwdz
               Dun = Sdudxn+Sdvdyn+Sdwdzn
               rhsp1(i,k) = c1*(gamma*Du + (1.0d0-gamma)*Dun) 
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(rhsp1)
#        endif
!        =============== END ================
!        ====================================
        !------------
         c2 = -rho_f/gamma*(1.0d0-gamma) - 1.0d0/gamma                 
         call grandientLSM(dpdx,dpdy,dpds,pfn,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0                
               Tpx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
               Tpy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
               Tpz(i,k) = sigmax(i,k)*Tpx(i,k) &
                         +sigmay(i,k)*Tpy(i,k) &
                         +sigmaz(i,k)*(sigmaz(i,k)*dpds(i,k))
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Tpx)
            call communication3Dtype2(Tpy)
            call communication3Dtype2(Tpz)
#        endif
!        =============== END ================
!        ====================================        
         call grandientLSM(dTpxdx,aux,aux,Tpx,No_cp,nbe,sig)
         call grandientLSM(aux,dTpydy,aux,Tpy,No_cp,nbe,sig)
         call grandientLSM(aux,aux,dTpzds,Tpz,No_cp,nbe,sig)         
         do k=2,NZ-1
            do i=1,N_CELL0
               Lpn = dTpxdx(i,k)+dTpydy(i,k)+dTpzds(i,k)     
               rhsp2(i,k) = c2*Lpn
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(rhsp2)
#        endif
!        =============== END ================
!        ====================================
         do k=2,NZ-1
            do i=1,N_CELL0                     
               rhsp(i,k) = Hpr(i)*VolPrism(i,k)*(rhsp1(i,k)+rhsp2(i,k))
            enddo
         enddo 
         
      ENDIF

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(rhsp)
#     endif
!     =============== END ================
!     ==================================== 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: Calculate_rhsP2'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!---------------------------------------------------------------------!
!   AUXILIAR 2: CORRECTOR STEP OF THE FLUID FLOW FOR THE PROJECTION   !
!                             (First order)                           !
!                               May 2016                              !
!---------------------------------------------------------------------!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      SUBROUTINE Corrector_Flow1(Hufnp,ufnp,Hufv,ufv,          &
                                Hvfnp,vfnp,Hvfv,vfv,           &
                                Hwfnp,wfnp,Hwfv,wfv,           &
!                               --------------------------------                                
                                ufstar,vfstar,wfstar,          &
!                               --------------------------------
                                pfnp,pf,pfv,                   &
!                               --------------------------------
                                Hprnp,Hpr,Hprv,                &
                                etanp,eta,etav,                &
!                               --------------------------------                                
                                dhdx,dhdy,                     &
                                dhdxv,dhdyv,                   &
!                               --------------------------------                                  
                                detadt,detadx,detady,          &                                                          
                                detadtv,detadxv,detadyv,       &                                                                
!                               --------------------------------  
                                h,xc,yc,sig,dsig,No_cp,nbe,    &
                                hv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                               --------------------------------                                
                                gamma) 
                                           
!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the velocity components resulting from the  !
!    projection technique of first-order accurate. It uses the        !
!    gradient of p_new already updated by the Poisson equation.       !
!                                                                     !                                                              !
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

      real*8,dimension(:,:) :: Hufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Hvfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Hwfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)  
      real*8,dimension(:,:) :: Hufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hvfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hwfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)      
!     --------------------------------------  
      real*8,dimension(:,:) :: ufstar(N_CELL,NZ)
      real*8,dimension(:,:) :: vfstar(N_CELL,NZ)
      real*8,dimension(:,:) :: wfstar(N_CELL,NZ)      
!     --------------------------------------    
      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pf(N_CELL,NZ)      
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
!     --------------------------------------      
      real*8, dimension(:)  :: Hprnp(N_CELL)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT) 
      real*8, dimension(:)  :: etanp(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: etav(N_VERT)
!     ----------------------------------------      
      real*8,dimension(:)   :: dhdx(N_CELL)
      real*8,dimension(:)   :: dhdy(N_CELL)
      real*8,dimension(:)   :: dhdxv(N_VERT)
      real*8,dimension(:)   :: dhdyv(N_VERT)
!     ----------------------------------------
      real*8,dimension(:)   :: detadt(N_CELL)      
      real*8,dimension(:)   :: detadx(N_CELL)
      real*8,dimension(:)   :: detady(N_CELL)  
      real*8,dimension(:)   :: detadtv(N_VERT)      
      real*8,dimension(:)   :: detadxv(N_VERT)
      real*8,dimension(:)   :: detadyv(N_VERT)      
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
                 
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
     
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ) 
      real*8, dimension(:,:) :: Sdpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdz(N_CELL,NZ)        
      real*8, dimension(:,:) :: dpdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdsn(N_CELL,NZ) 
      real*8, dimension(:,:) :: Sdpdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdzn(N_CELL,NZ)       
!     --------------------------------------
      real*8 :: Sdpdx_gam,Sdpdy_gam,Sdpdz_gam
      real*8 :: uT,vT,uB,vB

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: Corrector_Flow1'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                 Update (Cell-center)                   |
!     |________________________________________________________|
 
         call grandientLSM(dpdxn,dpdyn,dpdsn,pf,No_cp,nbe,sig)
         call grandientLSM(dpdx,dpdy,dpds,pfnp,No_cp,nbe,sig)
         
         do k=1,NZ
            do i=1,N_CELL0
               Sdpdx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
               Sdpdy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
               Sdpdz(i,k) = sigmaz(i,k)*dpds(i,k)
               !----------
               Sdpdxn(i,k) = sigmax(i,k)*dpdsn(i,k) + dpdxn(i,k)
               Sdpdyn(i,k) = sigmay(i,k)*dpdsn(i,k) + dpdyn(i,k)
               Sdpdzn(i,k) = sigmaz(i,k)*dpdsn(i,k)
               !----------
               Sdpdx_gam = gamma*Sdpdx(i,k)+(1.0d0-gamma)*Sdpdxn(i,k)
               Sdpdy_gam = gamma*Sdpdy(i,k)+(1.0d0-gamma)*Sdpdyn(i,k)
               Sdpdz_gam = gamma*Sdpdz(i,k)+(1.0d0-gamma)*Sdpdzn(i,k)
               !----------              
               ufnp(i,k)  = ufstar(i,k) - (dt/rho_f)*Sdpdx_gam                                      
               vfnp(i,k)  = vfstar(i,k) - (dt/rho_f)*Sdpdy_gam
               wfnp(i,k)  = wfstar(i,k) - (dt/rho_f)*Sdpdz_gam
               !----------
               Hufnp(i,k) = ufnp(i,k)*Hpr(i)
               Hvfnp(i,k) = vfnp(i,k)*Hpr(i)
               Hwfnp(i,k) = wfnp(i,k)*Hpr(i)
            enddo
         enddo

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Hufnp)
            call communication3Dtype2(Hvfnp)
            call communication3Dtype2(Hwfnp)
            call communication3Dtype2(ufnp)            
            call communication3Dtype2(vfnp)            
            call communication3Dtype2(wfnp)
#        endif	
!        =============== END ================    
!        ====================================

!      ________________________________________________________
!     |                                                        |
!     |       Boundary conditions & vertex values: u           |
!     |________________________________________________________| 
           
!     --------------------------
!     BC: cell-center
      call BCvelocity3D(Hufnp,ufnp,                  &
                        Hufv,ufv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        1)
!     --------------------------
!     Vertex values
      call interpolation3D(Hufv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Hufnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(ufv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           ufnp,xc,yc,sig,dsig,No_cp,nbe)
!     --------------------------
!     BC: vertex
      call BCvelocity3D(Hufnp,ufnp,                  &
                        Hufv,ufv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -1)
!      ________________________________________________________
!     |                                                        |
!     |       Boundary conditions & vertex values: v           |
!     |________________________________________________________|

!     --------------------------
!     BC: cell-center
      call BCvelocity3D(Hvfnp,vfnp,                  &
                        Hvfv,vfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        2)
!     --------------------------
!     Vertex values              
      call interpolation3D(Hvfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                          Hvfnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(vfv,xv,yv,sigv,dsigv,No_vp,nbev, &
                          vfnp,xc,yc,sig,dsig,No_cp,nbe)
!     --------------------------
!     BC: vertex
      call BCvelocity3D(Hvfnp,vfnp,                  &
                        Hvfv,vfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -2)
!      ________________________________________________________
!     |                                                        |
!     |       Velocity componet w at the bottom & top          |
!     |________________________________________________________|
                                                         
!     ____________________________
!     Cell-center
      do i=1,N_CELL
         uB = ufnp(i,1)
         vB = vfnp(i,1)
         wfB(i) = -uB*dhdx(i)-vB*dhdy(i)
         uT = ufnp(i,NZ)
         vT = vfnp(i,NZ)
         wfT(i) = detadt(i)+uT*detadx(i)+vT*detady(i)
      enddo
!     ____________________________
!     Vertex
      do nv=1,N_VERT
         uB = ufv(nv,1)
         vB = vfv(nv,1)  
         wfvB(nv) = -uB*dhdxv(nv)-vB*dhdyv(nv)
         uT = ufv(nv,NZ-1)
         vT = vfv(nv,NZ-1)
         wfvT(nv) = detadtv(nv)+uT*detadxv(nv)+vT*detadyv(nv)              
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(wfB)
         call communication2D(wfT)
#     endif	
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |       Boundary conditions & vertex values: w           |
!     |________________________________________________________|

!     --------------------------
!     BC: cell-center
      call BCvelocity3D(Hwfnp,wfnp,                  &
                        Hwfv,wfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        3)
!     --------------------------
!     Vertex values
      call interpolation3D(Hwfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Hwfnp,xc,yc,sig,dsig,No_cp,nbe)                                   
      call interpolation3D(wfv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           wfnp,xc,yc,sig,dsig,No_cp,nbe)

!     --------------------------
!     BC: vertex
      call BCvelocity3D(Hwfnp,wfnp,                  &
                        Hwfv,wfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -3)                             
!      ________________________________________________________
!     |                                                        |
!     |                      Finalization                      |
!     |________________________________________________________|

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Hufnp)
            call communication3Dtype2(Hvfnp)
            call communication3Dtype2(Hwfnp)
            call communication3Dtype2(ufnp)            
            call communication3Dtype2(vfnp)            
            call communication3Dtype2(wfnp)
#        endif	
!        =============== END ================    
!        ====================================

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: Corrector_Flow1'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!
!---------------------------------------------------------------------!
!   AUXILIAR 3: CORRECTOR STEP OF THE FLUID FLOW FOR THE PROJECTION   !
!                             (Second order)                          !
!                               May 2016                              !
!---------------------------------------------------------------------!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss!

      SUBROUTINE Corrector_Flow2(Hufnp,ufnp,Hufv,ufv,          &
                                Hvfnp,vfnp,Hvfv,vfv,           &
                                Hwfnp,wfnp,Hwfv,wfv,           &
!                               --------------------------------                                
                                ufstar,vfstar,wfstar,          &
!                               --------------------------------
                                pfnp,pf,pfv,                   &
!                               --------------------------------
                                Hprnp,Hpr,Hprv,                &
                                etanp,eta,etav,                &
!                               --------------------------------                                
                                dhdx,dhdy,                     &
                                dhdxv,dhdyv,                   &
!                               --------------------------------                                  
                                detadt,detadx,detady,          &                                                          
                                detadtv,detadxv,detadyv,       &                                                                
!                               --------------------------------  
                                h,xc,yc,sig,dsig,No_cp,nbe,    &
                                hv,xv,yv,sigv,dsigv,No_vp,nbev,&
!                               --------------------------------                                
                                gamma) 
                                           
!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the velocity components resulting from the  !
!    projection technique of second-order accurate. It uses the       !
!    gradient of dp = p_new-p_old already updated by the Poisson eqn. !
!                                                                     !                                                              !
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

      real*8,dimension(:,:) :: Hufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Hvfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: Hwfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)  
      real*8,dimension(:,:) :: Hufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hvfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hwfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)      
!     --------------------------------------  
      real*8,dimension(:,:) :: ufstar(N_CELL,NZ)
      real*8,dimension(:,:) :: vfstar(N_CELL,NZ)
      real*8,dimension(:,:) :: wfstar(N_CELL,NZ)      
!     --------------------------------------    
      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pf(N_CELL,NZ)      
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
!     --------------------------------------      
      real*8, dimension(:)  :: Hprnp(N_CELL)
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT) 
      real*8, dimension(:)  :: etanp(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: etav(N_VERT)
!     ----------------------------------------      
      real*8,dimension(:)   :: dhdx(N_CELL)
      real*8,dimension(:)   :: dhdy(N_CELL)
      real*8,dimension(:)   :: dhdxv(N_VERT)
      real*8,dimension(:)   :: dhdyv(N_VERT)
!     ----------------------------------------
      real*8,dimension(:)   :: detadt(N_CELL)      
      real*8,dimension(:)   :: detadx(N_CELL)
      real*8,dimension(:)   :: detady(N_CELL)  
      real*8,dimension(:)   :: detadtv(N_VERT)      
      real*8,dimension(:)   :: detadxv(N_VERT)
      real*8,dimension(:)   :: detadyv(N_VERT)      
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
     
      real*8, dimension(:,:) :: dpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dpds(N_CELL,NZ) 
      real*8, dimension(:,:) :: Sdpdx(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdy(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdz(N_CELL,NZ)        
      real*8, dimension(:,:) :: dpdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: dpdsn(N_CELL,NZ) 
      real*8, dimension(:,:) :: Sdpdxn(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdyn(N_CELL,NZ)
      real*8, dimension(:,:) :: Sdpdzn(N_CELL,NZ)       
!     --------------------------------------
      real*8 :: Sdpdx_gam,Sdpdy_gam,Sdpdz_gam
      real*8 :: uT,vT,uB,vB

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: Corrector_Flow2'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |               Update (inside cell-centers)             |
!     |________________________________________________________|
 
!        ____________________________
!        Update the pressure
         do k=1,NZ
            do i=1,N_CELL0
               pfnp(i,k) = pfnp(i,k) + pf(i,k)
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel           
            call communication3Dtype2(pfnp)
#        endif	
!        =============== END ================
!        ====================================
!        ____________________________
!        Vertex interpolation of p                       
         call interpolation3D(pfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              pfnp,xc,yc,sig,dsig,No_cp,nbe)
!        ____________________________
!        Update the velocity                                        
         call grandientLSM(dpdx,dpdy,dpds,pfnp,No_cp,nbe,sig) 
         do k=1,NZ
            do i=1,N_CELL0
               Sdpdx(i,k) = sigmax(i,k)*dpds(i,k) + dpdx(i,k)
               Sdpdy(i,k) = sigmay(i,k)*dpds(i,k) + dpdy(i,k)
               Sdpdz(i,k) = sigmaz(i,k)*dpds(i,k)            
               ufnp(i,k)  = ufstar(i,k)-(dt/rho_f)*Sdpdx(i,k)
               vfnp(i,k)  = vfstar(i,k)-(dt/rho_f)*Sdpdy(i,k)
               wfnp(i,k)  = wfstar(i,k)-(dt/rho_f)*Sdpdz(i,k)
               Hufnp(i,k) = ufnp(i,k)*Hpr(i)
               Hvfnp(i,k) = vfnp(i,k)*Hpr(i)
               Hwfnp(i,k) = wfnp(i,k)*Hpr(i)
            enddo
         enddo

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Hufnp)
            call communication3Dtype2(Hvfnp)
            call communication3Dtype2(Hwfnp)
            call communication3Dtype2(ufnp)            
            call communication3Dtype2(vfnp)            
            call communication3Dtype2(wfnp)
#        endif	
!        =============== END ================    
!        ====================================

!      ________________________________________________________
!     |                                                        |
!     |       Boundary conditions & vertex values: u           |
!     |________________________________________________________| 
           
!     --------------------------
!     BC: cell-center
      call BCvelocity3D(Hufnp,ufnp,                  &
                        Hufv,ufv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        1)
!     --------------------------
!     Vertex values
      call interpolation3D(Hufv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Hufnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(ufv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           ufnp,xc,yc,sig,dsig,No_cp,nbe)
!     --------------------------
!     BC: vertex
      call BCvelocity3D(Hufnp,ufnp,                  &
                        Hufv,ufv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -1)
!      ________________________________________________________
!     |                                                        |
!     |       Boundary conditions & vertex values: v           |
!     |________________________________________________________|

!     --------------------------
!     BC: cell-center
      call BCvelocity3D(Hvfnp,vfnp,                  &
                        Hvfv,vfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        2)
!     --------------------------
!     Vertex values              
      call interpolation3D(Hvfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                          Hvfnp,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(vfv,xv,yv,sigv,dsigv,No_vp,nbev, &
                          vfnp,xc,yc,sig,dsig,No_cp,nbe)
!     --------------------------
!     BC: vertex
      call BCvelocity3D(Hvfnp,vfnp,                  &
                        Hvfv,vfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -2)
!      ________________________________________________________
!     |                                                        |
!     |       Velocity componet w at the bottom & top          |
!     |________________________________________________________|
                                                         
!     ____________________________
!     Cell-center
      do i=1,N_CELL
         uB = ufnp(i,1)
         vB = vfnp(i,1)
         wfB(i) = -uB*dhdx(i)-vB*dhdy(i)
         uT = ufnp(i,NZ)
         vT = vfnp(i,NZ)
         wfT(i) = detadt(i)+uT*detadx(i)+vT*detady(i)
      enddo
!     ____________________________
!     Vertex
      do nv=1,N_VERT
         uB = ufv(nv,1)
         vB = vfv(nv,1)  
         wfvB(nv) = -uB*dhdxv(nv)-vB*dhdyv(nv)
         uT = ufv(nv,NZ-1)
         vT = vfv(nv,NZ-1)
         wfvT(nv) = detadtv(nv)+uT*detadxv(nv)+vT*detadyv(nv)              
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(wfB)
         call communication2D(wfT)
#     endif	
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |       Boundary conditions & vertex values: w           |
!     |________________________________________________________|

!     --------------------------
!     BC: cell-center
      call BCvelocity3D(Hwfnp,wfnp,                  &
                        Hwfv,wfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        3)
!     --------------------------
!     Vertex values
      call interpolation3D(Hwfv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Hwfnp,xc,yc,sig,dsig,No_cp,nbe)                                   
      call interpolation3D(wfv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           wfnp,xc,yc,sig,dsig,No_cp,nbe)

!     --------------------------
!     BC: vertex
      call BCvelocity3D(Hwfnp,wfnp,                  &
                        Hwfv,wfv,                    &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -3)                             
                           
!      ________________________________________________________
!     |                                                        |
!     |                      Finalization                      |
!     |________________________________________________________|

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype2(Hufnp)
            call communication3Dtype2(Hvfnp)
            call communication3Dtype2(Hwfnp)
            call communication3Dtype2(ufnp)            
            call communication3Dtype2(vfnp)            
            call communication3Dtype2(wfnp)
#        endif	
!        =============== END ================    
!        ====================================

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: Corrector_Flow2'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    TEST NAVIER-STOKES EQUATION                      !
!                              May 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Vorticity(ufnp,vfnp,wfnp,            &                           
                           xc,yc,sig,dsig,No_cp,nbe,  &
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

      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)  
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(N_CELL,NZ) :: dudx,dudy,duds
      real*8,dimension(N_CELL,NZ) :: dvdx,dvdy,dvds
      real*8,dimension(N_CELL,NZ) :: dwdx,dwdy,dwds
!     --------------------------------------
      real*8 :: Sdudx,Sdudy,Sdudz
      real*8 :: Sdvdx,Sdvdy,Sdvdz
      real*8 :: Sdwdx,Sdwdy,Sdwdz
      
!      ________________________________________________________
!     |                                                        |
!     |         Calculate Vorticity (vector & magnitud)        |
!     |________________________________________________________|

      call grandientLSM(dudx,dudy,duds,ufnp,No_cp,nbe,sig)
      call grandientLSM(dvdx,dvdy,dvds,vfnp,No_cp,nbe,sig)
      call grandientLSM(dwdx,dwdy,dwds,wfnp,No_cp,nbe,sig)
!     ----------------------------- 
      do k=1,NZ
         do i=1,N_CELL  
            Sdudy = sigmay(i,k)*duds(i,k) + dudy(i,k)
            Sdudz = sigmaz(i,k)*duds(i,k)
            Sdvdx = sigmax(i,k)*dvds(i,k) + dvdx(i,k)
            Sdvdz = sigmaz(i,k)*dvds(i,k)
            Sdwdx = sigmax(i,k)*dwds(i,k) + dwdx(i,k)
            Sdwdy = sigmay(i,k)*dwds(i,k) + dwdy(i,k)
            Vortx(i,k) = Sdwdy-Sdvdz
            Vorty(i,k) = Sdudz-Sdwdx
            Vortz(i,k) = Sdvdx-Sdudy
            Vort(i,k)  = dsqrt((Sdwdy-Sdvdz)**2  &
                              +(Sdudz-Sdwdx)**2  &
                              +(Sdvdx-Sdudy)**2)                                    
         enddo
      enddo
!     ____________________________
!     Vertex interpolation
      call interpolation3D(Vortxv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vortx,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Vortyv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vorty,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Vortzv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vortz,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(Vortv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           Vort,xc,yc,sig,dsig,No_cp,nbe)

      RETURN
      END
                       
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF THE NAVIER-STOKES EQUATION                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
