!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          INITIAL VALUES                             !
!                             March 2017                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE input_initial(alphafn,ufn,vfn,wfn,pfn,viscof,rhof,  &
                         alphasn,usn,vsn,wsn,psn,viscos,rhos,        &
                         alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv,      &
                         alphasv,usv,vsv,wsv,psv,viscosv,rhosv,      &
                         etan,etav,                                  &
                         Hpr,Hprv,                                   &
                         h,hv,                                       &
                         xc,yc,sig,dsig,No_cp,nbe,                   &
                         xv,yv,sigv,dsigv,No_vp,nbev,                &
                         Heaviside,mask)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the initial condition of the main        !
!    variables of our problem. We also impose the initial volume      !
!    fraction and velocity ws in the input region (a rectangular      !
!    box or a cylinder)                                               !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !
!  |______________|___________|____________________________________|  !
!  | <-- alphafn  |(N_CELL,NZ)| Fluid control volume at (n)        |  !
!  | <-- ufn      |(N_CELL,NZ)| Velocity component u_f at (n)      |  !
!  | <-- vfn      |(N_CELL,NZ)| Velocity component v_f at (n)      |  !
!  | <-- wfn      |(N_CELL,NZ)| Velocity component w_f at (n)      |  !
!  | <-- pfn      |(N_CELL,NZ)| Pressure of the fluid at (n)       |  !
!  | <-- viscof   |(N_CELL,NZ)| Viscosity of the fluid             |  !
!  | <-- rof      |(N_CELL,NZ)| Density of the fluid               |  !
!  |______________|___________|____________________________________|  !
!  | <-- alphasn  |(N_CELL,NZ)| Solid control volume at (n)        |  !
!  | <-- usn      |(N_CELL,NZ)| Velocity component u_s at (n)      |  !
!  | <-- vsn      |(N_CELL,NZ)| Velocity component v_s at (n)      |  !
!  | <-- wsn      |(N_CELL,NZ)| Velocity component w_s at (n)      |  !
!  | <-- psn      |(N_CELL,NZ)| Pressure of the solid at (n)       |  !
!  | <-- viscos   |(N_CELL,NZ)| Viscosity of the solid             |  !
!  | <-- ros      |(N_CELL,NZ)| Density of the solid               |  !
!  |______________|___________|____________________________________|  !
!  | <-- etan     | N_CELL    | Free surface level eta at (n)      |  !
!  | <-- etav     | N_VERT    | Free surface eta at the vertex     |  !
!  | <-- Hpr      | N_CELL    | Total water depth: H = h + eta     |  !
!  | <-- Hprv     | N_VERT    | Total water depth at the vertex    |  !
!  |______________|___________|____________________________________|  !
!  | <-- Heaviside| N_CELL    | Step function=0 outside, =1 inisde |  !
!  | <-- mask     |(N_CELL,NZ)| Profile of the release input       |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   |         Description                 |  !
!  |_____________|___________|_____________________________________|  !
!  | --> h       | N_CELL    | Depth of the domain at each cell    |  !
!  | --> hv      | N_CELL    | Depth of the domain at the vertex   |  !
!  | --> xc      | N_CELL    | x-coordinate of the cell center     |  !
!  | --> yc      | N_CELL    | y-coordinate of the cell center     |  !
!  | --> xv      | N_VERT    | x-coordinate of the vertex          |  !
!  | --> yv      | N_VERT    | y-coordinate of the vertex          |  !
!  | --> sig     | NZ        | sigma of the cell centers           |  !
!  | --> sigv    | NZ-1      | sigma of the vertex                 |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !
!  |_____________|_________________________________________________|  !
!  |--- N_CELL0  | Number of the cell centers inside the domain    |  !
!  |--- N_CELL   | Total number of the cells                       |  !
!  |--- N_VERT   | Number of the computing vertices                |  !
!  |    NZ       | Number of vertical points                       |  !
!  |--- pi       | 3.14159....                                     |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name          |              Description                    |  !
!  |_________________|_____________________________________________|  !
!  | * gra           | gravity = 9.81                              |  !
!  | * pa            | atmospheric pressure                        |  !
!  | * Re            | Reynolds number                             |  !
!  | * rho_s         | Homogeneous density of the solid            |  !
!  | * rho_f         | Homogeneous density of the fluid            |  !
!  | * visco_f       | Homogeneous viscosity of the fluid          |  !
!  | * alphasmin     | Minimum alpha_s                             |  !
!  | * alphasmax     | Maximum alpha_s                             |  !
!  |_________________|_____________________________________________|  !
!  | * ISOLID        | Tag to calculate solid                      |  !
!  | * ICASE         | Tag about the release box                   |  !
!  | * IBED          | Tag about the granular bed option           |  !
!  | * ChooseREGION  | Choose mesh                                 |  !
!  | * ChooseCASE    | Choose the problem parameters               |  !
!  | * ChoosePROFILE | Choose injection profile                    |  !
!  | * xRegMin       | x coodinate west box boundary               |  !
!  | * xRegMax       | x coodinate east box boundary               |  !
!  | * yRegMin       | y coodinate south box boundary              |  !
!  | * yRegMax       | y coodinate north box boundary              |  !
!  | * kRegMin       | k index of the bottom box                   |  !
!  | * kRegMax       | k index of the top box                      |  !
!  | * xcReg         | xc coordinate of the circle center          |  !
!  | * ycReg         | yc coordinate of the circle center          |  !
!  | * rReg          | radius of the circle                        |  !
!  | * Aini          | initial volumen fraction                    |  !
!  | * Wini          | initial velocity ws                         |  !
!  | * rs            | Radius of the particules                    |  !
!  | * drejet        | Diameter of the release                     |  !
!  | * cini          | Concentration of the release                |  !
!  | * vrejet        | Volume of the release                       |  !
!  | * rmassrejet    | Total mass release in the simulation        |  !
!  |_________________|_____________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  | rr          | radius of a cell center to (xcReg,ycReg)        |  !
!  | kZBed       | index k until the sediment bed is               |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |                 -     interpolationEta.F90                    |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!        Common variables used                                        !
!    *   Common variables modified                                    !
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

      real*8,dimension(:,:) :: alphafn,ufn,vfn,wfn,pfn,viscof,rhof      
      real*8,dimension(:,:) :: alphasn,usn,vsn,wsn,psn,viscos,rhos   
      real*8,dimension(:,:) :: alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv
      real*8,dimension(:,:) :: alphasv,usv,vsv,wsv,psv,viscosv,rhosv
      real*8,dimension(:)   :: etan,etav   
      real*8,dimension(:)   :: Hpr,Hprv
      real*8,dimension(:)   :: h,hv
      real*8,dimension(:)   :: xc,yc,sig,dsig
      real*8,dimension(:)   :: xv,yv,sigv,dsigv
      integer,dimension(:,:):: No_cp,No_vp
      integer,dimension(:)  :: nbe,nbev  
      real*8,dimension(:)   :: Heaviside
      real*8,dimension(:,:) :: mask
      
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8  :: rr,aaReg,bbReg,tt
      integer :: kZBed 
!     ----------------------------------------
      real*8  :: x,y,z,c
!     ----------------------------------------
      real*8 :: funEta,funh
      real*8 :: funExamNSu,funExamNSv,funExamNSw,funExamNSp
      real*8 :: funExample1u,funExample1v,funExample1w
      real*8 :: funExample1p
      real*8 :: FS_funu,FS_funv,FS_funw,FS_funeta 
      real*8 :: Vref,Dref,nuref,Href
      real*8 :: c0
!     ----------------------------------------
      integer:: InitialRandom
      real*8 :: funu_ChannelFlow,funv_ChannelFlow
      real*8 :: funw_ChannelFlow,funp_ChannelFlow 
      integer :: fcG
!     ----------------------------------------
      real*8 :: funInflow

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: initial '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                Calculation of the initial conditions                !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                   General parameters                   |
!     |________________________________________________________|

      pa        = 0.0000d+05
      Re        = 7.5000d+01
      rho_f     = 1.0000d+00 
      rho_s     = 2.6500d+03
      visco_f   = 1.0476d-06
      Fr        = 0.2000d+00
!     ----------------------
      alphasmin = 1.0000d-08
      alphasmax = 0.6250d+00
!     ----------------------
      Uinflow   = 1.0d0

!*********************************************************************!
!                                                                     !
!              Initial conditions by different test cases             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     | ====================================================== |
!     |                  CASE: Estuary Gironde                 |
!     |________________________________________________________|

#     ifdef KeyEstuaryGironde
!        ______________________________________________________
!        Non-dimensional parameters
         H0 = 3.0d0
         Re = 1.46d4
         Fr = 0.21d0
!        ______________________________________________________
!        Free surface: eta & Total water depth: H
!        ------------
!        Cell center
         do i=1,N_CELL
            h(i)    = H0   !Miguel Nov2017
            etan(i) = 0.0d0
            Hpr(i)  = etan(i) + h(i)
         enddo
!        ------------
!        Vertex 
         do nv=1,N_VERT
            hv(nv)   = H0   !Miguel Nov2017
            etav(nv) = 0.0d0
            Hprv(nv) = etav(nv) + hv(nv)
         enddo
!        ______________________________________________________
!        Velocity BC inflow magnitude
         Uinflow = 1.0d0
!        ______________________________________________________
!        Velocity & pressure (Navier-Stokes Equations)
!        ------------
!        Cell center 
         ufn = 1.0d0
         vfn = 0.0d0
         wfn = 0.0d0
         pfn = 0.0d0
!        ------------
!        Vertex
         ufv = 1.0d0
         vfv = 0.0d0
         wfv = 0.0d0
         pfv = 0.0d0
!        _______________________________________________________
!        Velocity profile fixed bottom
#        ifdef KeyBCbottomNoSlip
!        ------------
!        Cell center
         do i=1,N_CELL
            ufn(i,1) = 0.0d0
            do k=2,NZ
               x = xc(i)
               y = yc(i)
               z = sig(k)*Hpr(i)-h(i)
               ufn(i,k) = funInflow(x,y,z,time)
            enddo
         enddo
!        ------------
!        Vertex
         do nv=1,N_VERT
            ufv(nv,1) = 0.0d0
            do k=2,NZ-1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               ufv(nv,k) = funInflow(x,y,z,time)
            enddo
         enddo
#        endif
#     endif
!      ________________________________________________________
!     | ====================================================== |
!     |                  CASE: Static Cylinder                 |
!     |________________________________________________________|

#     ifdef KeyStaticCylinder
         fcG = 1075 ! <<<--- Choose test case: Cylinder+FS
!        ______________________________________________________
!        Non-dimensional parameters
!        _________________________
!        General case
         if (fcG.eq.0) then
            Dref  = 0.06d0
            nuref = 1.1386d-6
            Href  = 0.180d0 !0.058d0
            Vref  = 0.15086d0 !0.43730554d0
            H0    = Href/Dref
            Re    = Vref*Dref/nuref
            Fr    = Vref/sqrt(gra*Href)
!        _______________________________
!        Experimental comparison: Johnson et al. 2013
!        -------- fcG1
         elseif (fcG.eq.1) then
            H0 = 1.0d0
            Re = 2.33d4
            Fr = 0.57d0
!        -------- fcG2
         elseif (fcG.eq.2) then
            H0 = 0.97d0
            Re = 0.79d4
            Fr = 0.20d0
!        -------- fcG3
         elseif (fcG.eq.3) then
            H0 = 3.0d0
            Re = 1.46d4
            Fr = 0.21d0
!        -------- fcG4
         elseif (fcG.eq.4) then
            H0 = 2.02d0
            Re = 1.67d4
            Fr = 0.29d0
!        -------- fcG5
         elseif (fcG.eq.5) then
            H0 = 2.06d0
            Re = 2.52d4
            Fr = 0.43d0
!        -------- fcG6
         elseif (fcG.eq.6) then
            H0 = 1.10d0
            Re = 2.98d4
            Fr = 0.69d0
!        -------- fcG7
         elseif (fcG.eq.7) then
            H0 = 1.03d0
            Re = 1.41d4
            Fr = 0.34d0
!        -------- fcG8
         elseif (fcG.eq.8) then
            H0 = 0.99d0
            Re = 3.29d4
            Fr = 0.81d0
!        -------- fcG9
         elseif (fcG.eq.9) then
            H0 = 0.96d0
            Re = 1.99d4
            Fr = 0.49d0
!        ======== fcH1
         elseif (fcG.eq.11) then
            H0 = 1.02d0
            Re = 2.73d4
            Fr = 0.65d0
!        ======== fcH2
         elseif (fcG.eq.12) then
            H0 = 1.04d0
            Re = 0.76d4
            Fr = 0.18d0
!        ======== fcH3
         elseif (fcG.eq.13) then
            H0 = 0.99d0
            Re = 2.30d4
            Fr = 0.56d0
!        ======== fcH4
         elseif (fcG.eq.14) then
            H0 = 1.06d0
            Re = 1.85d4
            Fr = 0.44d0
!        ======== fcH5
         elseif (fcG.eq.15) then
            H0 = 1.03d0
            Re = 1.50d4
            Fr = 0.36d0
!        ======== fcH6
         elseif (fcG.eq.16) then
            H0 = 1.03d0
            Re = 1.15d4
            Fr = 0.27d0
!        _______________________________
!        Formula comparison: Roulund et al. 2005
!        --------
         elseif (fcG.eq.1011) then
            H0 = 1.0d0
            Re = 800.0d0
            Fr = 0.011d0
!        --------
         elseif (fcG.eq.1025) then
            H0 = 1.0d0
            Re = 1010.72d0
            Fr = 0.025d0
!        --------
         elseif (fcG.eq.1050) then
            H0 = 1.0d0
            Re = 2021.44d0
            Fr = 0.050d0
!        --------
         elseif (fcG.eq.1075) then
            H0 = 1.0d0
            Re = 3032.16d0
            Fr = 0.075d0
!        --------
         elseif (fcG.eq.1100) then
            H0 = 1.0d0
            Re = 4042.87d0
            Fr = 0.100d0
!        _______________________________
!        Using gravity
         else
            H0 = 3.0d0
            Re = 7.5000d+01
            Fr = 1.0d0/sqrt(gra)
         endif    
!        ______________________________________________________
!        Free surface: eta & Total water depth: H
!        ------------
!        Cell center
         do i=1,N_CELL
            h(i)    = H0   !Miguel Nov2017
            etan(i) = 0.0d0
            Hpr(i)  = etan(i) + h(i)
         enddo
!        ------------
!        Vertex 
         do nv=1,N_VERT
            hv(nv)   = H0   !Miguel Nov2017
            etav(nv) = 0.0d0
            Hprv(nv) = etav(nv) + hv(nv)
         enddo
!        ______________________________________________________
!        Velocity BC inflow magnitude
         Uinflow = 1.0d0 
!        ______________________________________________________
!        Velocity & pressure (Non-dimensional)
!        ------------
!        Cell center 
         ufn = Uinflow
         vfn = 0.0d0
         wfn = 0.0d0
         pfn = 0.0d0
!        ------------
!        Vertex
         ufv = Uinflow
         vfv = 0.0d0
         wfv = 0.0d0
         pfv = 0.0d0
!        _______________________________________________________
!        Velocity profile fixed bottom
#        ifdef KeyBCbottomNoSlip
!        ------------
!        Cell center
         do i=1,N_CELL
            ufn(i,1) = 0.0d0
            do k=2,NZ
               x = xc(i)
               y = yc(i)
               z = sig(k)*Hpr(i)-h(i)
               ufn(i,k) = funInflow(x,y,z,time)
            enddo
         enddo
!        ------------
!        Vertex
         do nv=1,N_VERT
            ufv(nv,1) = 0.0d0
            do k=2,NZ-1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               ufv(nv,k) = funInflow(x,y,z,time)
            enddo
         enddo      
#        endif 
#     endif
!      ________________________________________________________
!     | ====================================================== |
!     |                  CASE: Static Channel                  |
!     |________________________________________________________|

#     ifdef KeyStaticChannel
!        ______________________________________________________
!        Non-dimensional parameters: body force
         forcex = 1.0d0
         forcey = 0.0d0
         forcez = 0.0d0
!        ______________________________________________________
!        Free surface: eta & Total water depth: H
#        ifdef KeyBCtopNoSlip
         H0 = 2.0d0
#        else
         H0 = 1.0d0
#        endif
!        ------------
!        Cell center
         do i=1,N_CELL
            if (IDEPTH.eq.1) h(i) = H0
            etan(i) = 0.0d0
            Hpr(i)  = etan(i) + h(i)
         enddo
!        ------------
!        Vertex 
         do nv=1,N_VERT
            if (IDEPTH.eq.1) hv(nv) = H0
            etav(nv) = 0.0d0
            Hprv(nv) = etav(nv) + hv(nv)
         enddo
!        ______________________________________________________
!        Velocity BC inflow magnitude
         Uinflow = 1.0d0
!        ______________________________________________________
!        Velocity & pressure (Non-dimensional)
         InitialRandom = 1
!        -------------------------
!        TURBULENT FLOW
         IF (InitialRandom.eq.1) THEN
            Re = 180.0d0
            c0 = 1.0d0
            !Re = 3300.0d0
            !c0 = 1.0d0/18.2d0
!           ------------
!           Cell center
            do i=1,N_CELL
               do k=1,NZ
                  x = xc(i)
                  y = yc(i)
                  z = sig(k) ! The function definition needs z[0,1]
                  ufn(i,k) = c0*funu_ChannelFlow(x,y,z) !<--(at Functions.F90)
                  vfn(i,k) = c0*funv_ChannelFlow(x,y,z) !<--(at Functions.F90)
                  wfn(i,k) = c0*funw_ChannelFlow(x,y,z) !<--(at Functions.F90)
                  pfn(i,k) = c0*funp_ChannelFlow(x,y,z) !<--(at Functions.F90)
               enddo
            enddo
!           ------------
!           Vertex
            do nv=1,N_VERT
               do k=1,NZ-1
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k) ! The function definition needs z[0,1]
                  ufv(nv,k) = c0*funu_ChannelFlow(x,y,z) !<--(at Functions.F90) 
                  vfv(nv,k) = c0*funv_ChannelFlow(x,y,z) !<--(at Functions.F90) 
                  wfv(nv,k) = c0*funw_ChannelFlow(x,y,z) !<--(at Functions.F90)           
                  pfv(nv,k) = c0*funp_ChannelFlow(x,y,z) !<--(at Functions.F90) 
               enddo
            enddo
!        -------------------------
!        LAMINAR FLOW
         ELSE
            Re = 10.0d0
!           ------------
!           Cell center 
            ufn = 4.0d0 
            vfn = 0.0d0
            wfn = 0.0d0
            pfn = 0.0d0
!           ------------
!           Vertex
            ufv = 4.0d0
            vfv = 0.0d0
            wfv = 0.0d0
            pfv = 0.0d0
         ENDIF
#     endif
!      ________________________________________________________
!     | ====================================================== |
!     |                  CASE: Standing Wave                   |
!     |________________________________________________________|

#     ifdef KeyStandingWave
!        ______________________________________________________
!        Non-dimensional parameters
         Re = 1.0000d+00 !Not used
!        ______________________________________________________
!        Free surface: eta & Total water depth: H
         h0 = 10.0d0
!        ------------
!        Cell center   
         do i=1,N_CELL
            h(i)    = 10.0d0
            etan(i) = FS_funeta(xc(i),yc(i),time)   !<---(at Functions.F90)
            Hpr(i)  = etan(i) + h(i)
         enddo
!        ------------
!        Vertex 
         do nv=1,N_VERT
            hv(nv)   = 10.0d0 
            etav(nv) = FS_funeta(xv(nv),yv(nv),time) !<--(at Functions.F90)
            Hprv(nv) = etav(nv) + hv(nv)
         enddo   
!        ______________________________________________________
!        Velocity & pressure (Navier-Stokes Equations)
!        ------------
!        Cell center
         do i=1,N_CELL
            do k=1,NZ
               x = xc(i)
               y = yc(i)
               z = sig(k)*Hpr(i)-h(i)
               ufn(i,k) = FS_funu(x,y,z,time) !<----(at Functions.F90)
               vfn(i,k) = FS_funv(x,y,z,time) !<----(at Functions.F90)
               wfn(i,k) = FS_funw(x,y,z,time) !<----(at Functions.F90)
               pfn(i,k) = 0.0d0
            enddo
         enddo         
!        ------------
!        Vertex
         do nv=1,N_VERT
            do k=1,NZ-1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               ufv(nv,k) = FS_funu(x,y,z,time) !<----(at Functions.F90)
               vfv(nv,k) = FS_funv(x,y,z,time) !<----(at Functions.F90)
               wfv(nv,k) = FS_funw(x,y,z,time) !<----(at Functions.F90)
               pfv(nv,k) = 0.0d0 
            enddo
         enddo 
#     endif

!      ________________________________________________________
!     | ====================================================== |
!     |                 CASE: KeyTaylorVortex                  |
!     |________________________________________________________|

#     if  defined(KeyTaylorVortex) 
!        ______________________________________________________
!        Non-dimensional parameters
         Re = 1.0000d+00
!        ______________________________________________________
!        Free surface: eta & Total water depth: H
!        ------------
!        Cell center
         do i=1,N_CELL
            h(i)    = funh(xc(i),yc(i))      !<--- (at Functions.F90)
            etan(i) = funEta(xc(i),yc(i))    !<--- (at Functions.F90)
            Hpr(i)  = etan(i) + h(i)
         enddo
!        ------------
!        Vertex 
         do nv=1,N_VERT
            hv(nv)   = funh(xv(nv),yv(nv))   !<--- (at Functions.F90)
            etav(nv) = funEta(xv(nv),yv(nv)) !<--- (at Functions.F90)
            Hprv(nv) = etav(nv) + hv(nv)
         enddo 
!        ______________________________________________________
!        Velocity & pressure (Navier-Stokes Equations)
!        ------------
!        Cell center
         do i=1,N_CELL
            do k=1,NZ
               x = xc(i)
               y = yc(i)
               z = sig(k)*Hpr(i)-h(i)
               ufn(i,k) = funExamNSu(x,y,z,time) !<---(at Functions.F90)
               vfn(i,k) = funExamNSv(x,y,z,time) !<---(at Functions.F90)
               wfn(i,k) = funExamNSw(x,y,z,time) !<---(at Functions.F90)
               pfn(i,k) = funExamNSp(x,y,z,time) !<---(at Functions.F90)
            enddo
         enddo
!        ------------
!        Vertex
         do nv=1,N_VERT
            do k=1,NZ-1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               ufv(nv,k) = funExamNSu(x,y,z,time) !<--(at Functions.F90) 
               vfv(nv,k) = funExamNSv(x,y,z,time) !<--(at Functions.F90)
               wfv(nv,k) = funExamNSw(x,y,z,time) !<--(at Functions.F90)
               pfv(nv,k) = funExamNSp(x,y,z,time) !<--(at Functions.F90)
            enddo
         enddo

#     endif

!      ________________________________________________________
!     | ====================================================== |
!     |                 CASE: TestOnlyPoisson                  |
!     |________________________________________________________|

#     if  defined(KeyTestOnlyPoisson)
!        ______________________________________________________
!        Non-dimensional parameters
         Re = 1.0000d+00 !<--- Not used
!        ______________________________________________________
!        Free surface: eta & Total water depth: H
!        ------------
!        Cell center
         do i=1,N_CELL
            h(i)    = funh(xc(i),yc(i))      !<--- (at Functions.F90)
            etan(i) = funEta(xc(i),yc(i))    !<--- (at Functions.F90)
            Hpr(i)  = etan(i) + h(i)
         enddo
!        ------------
!        Vertex 
         do nv=1,N_VERT
            hv(nv)   = funh(xv(nv),yv(nv))   !<--- (at Functions.F90)
            etav(nv) = funEta(xv(nv),yv(nv)) !<--- (at Functions.F90)
            Hprv(nv) = etav(nv) + hv(nv)
         enddo
!        ______________________________________________________
!        Velocity & pressure (Navier-Stokes Equations)
!        ------------
!        Cell center
         ufn = 1.0d0 !<--- Not used
         vfn = 0.0d0 !<--- Not used
         wfn = 0.0d0 !<--- Not used
         do i=1,N_CELL
            do k=1,NZ
               x = xc(i)
               y = yc(i)
               z = sig(k)*Hpr(i)-h(i)
               pfn(i,k) = funExamNSp(x,y,z,time) !<---(at Functions.F90)
            enddo
         enddo
!        ------------
!        Vertex
         ufv = 1.0d0 !<--- Not used
         vfv = 0.0d0 !<--- Not used
         wfv = 0.0d0 !<--- Not used
         do nv=1,N_VERT
            do k=1,NZ-1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               pfv(nv,k) = funExamNSp(x,y,z,time) !<--(at Functions.F90)
            enddo
         enddo
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                     Communication                      |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(Hpr)
         call communication2D(etan)
         call communication2D(h)
!        -------------------------
         call communication3D(ufn)
         call communication3D(vfn)
         call communication3D(pfn)
         call communication3D(wfn)
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!                           Other variables                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                Velocity solid: (us,vs,ws)= 0           |
!     |________________________________________________________|

!     ------------
!     Cell center
      do i=1,N_CELL
         do k=1,NZ
            usn(i,k) = ufn(i,k)
            vsn(i,k) = vfn(i,k)
            wsn(i,k) = wfn(i,k)
         enddo
      enddo
!     ------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            usv(nv,k) = ufv(nv,k)
            vsv(nv,k) = vfv(nv,k)
            wsv(nv,k) = wfv(nv,k)
         enddo
      enddo 
!      ________________________________________________________
!     |                                                        |
!     |                     Pressure solid                     |
!     |________________________________________________________|

!     ------------
!     Cell center
      do i=1,N_CELL
         do k=1,NZ
            !pfn(i,k) = pa + gra*rhof(i,k)*(1.d0-sig(k))*Hpr(i)
            psn(i,k) = pfn(i,k)
         enddo
      enddo
!     ------------
!     Vertex
      do i=1,N_VERT
         do k=1,NZ-1
            !pfv(nv,k) = pa + gra*rhofv(nv,k)*(1.d0-sigv(k))*Hprv(i)
            psv(nv,k) = pfv(nv,k)
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                     Volume fraction                    |
!     |________________________________________________________|

!     ------------
!     Cell center
      do k=1,NZ
         do i=1,N_CELL
             alphafn(i,k) = 1.0d0
             alphasn(i,k) = 0.0d0
         enddo
      enddo
!     ------------
!     Vertex 
      do k=1,NZ-1
         do nv=1,N_VERT
            alphafv(i,k) = 1.0d0
            alphasv(i,k) = 0.0d0
         enddo
      enddo 

!      ________________________________________________________
!     |                                                        |
!     |                  Density  &  Viscosity                 |
!     |________________________________________________________|

!     ------------
!     Cell center
      do k=1,NZ
         do i=1,N_CELL
            rhof(i,k)   = rho_f
            rhos(i,k)   = rho_s
            viscof(i,k) = visco_f
            viscos(i,k) = visco_f
         enddo
      enddo    
!     ------------
!     Vertex 
      do k=1,NZ-1
         do nv=1,N_VERT
            rhofv(nv,k)   = rho_f
            rhosv(nv,k)   = rho_s
            viscofv(nv,k) = visco_f
            viscosv(nv,k) = visco_f
         enddo
      enddo 
      
!*********************************************************************!
!                                                                     !
!                Initial condition of the release box                 !
!                                                                     !
!*********************************************************************!

      ChooseREGION   = 1     ! =1 rectangular box, =2 cylinder
      ChooseCASE     = 1     ! =1,2,3 Choose the problem parameters 
      ChoosePROFILE  = 2     ! =1 constant, =2 Poiseville, =3 Tube

!      ________________________________________________________
!     |                                                        |
!     |                        Region                          |
!     |________________________________________________________|

!      _________________________________
!     |                                 |
!     |            Parameters           |
!     |_________________________________|

!     --------------------------------------------------
!     Index sigma direction
      kRegMin = 50            
      kRegMax = 55

!     --------------------------------------------------
!     Parameters Rectangular Box
      xRegMin = 2.25d0        
      xRegMax = 2.75d0        
      yRegMin = 0.0d0         
      yRegMax = 1.0d0         

!     --------------------------------------------------
!     Parameters Cylinder
      xcReg  = 0.5d0*(xRegMin+xRegMax)
      ycReg  = 0.5d0*(yRegMin+yRegMax)
      aaReg  = 0.5d0*abs(xRegMax-xRegMin)
      bbReg  = 0.5d0*abs(yRegMax-yRegMin)
      rReg   = min(aaReg,bbReg)  

!      _________________________________
!     |                                 |
!     |    Heaviside step function      |
!     |_________________________________|

!     --------------------------------------------------
!     Rectangle box 
      IF (ChooseREGION.eq.1) THEN
         do i=1,N_CELL  
            Heaviside(i) = 0.0d0                
            if (((xc(i).ge.xRegMin).and.(xc(i).le.xRegMax)).and.& 
                ((yc(i).ge.yRegMin).and.(yc(i).le.yRegMax)))then
                Heaviside(i) = 1.0d0
            endif 
         enddo

!     --------------------------------------------------
!     Cylinder          
      ELSEIF (ChooseREGION.eq.2) THEN
         do i=1,N_CELL  
            Heaviside(i) = 0.0d0
            rr = sqrt( (xc(i)-xcReg)**2+(yc(i)-ycReg)**2)
            if (rr.le.rReg) then
               Heaviside(i) = 1.0d0
            endif 
         enddo
      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |                       Parameters                       |
!     |________________________________________________________|

!     --------------------------------------------------
!     HOGG TEST CASE
      if (ChooseCASE.eq.1) then   
         rs         = 200.0d-6      ! Radius of the particules
         drejet     = 0.008d0       ! Release diameter
         cini       = 0.0d0*rho_s   ! Release concentration
         vrejet     = 10d0          ! Release volume
         rmassrejet = 1.0d3         ! Total release mass
!     --------------------------------------------------
!     TEST 2
      elseif (ChooseCASE.eq.2) then    
         rs         = 45d-6            
         drejet     = 0.1d0             
         cini       = 100.d0           
         vrejet     = 4.6d0            
         rmassrejet = 4.64d-3*CINI/0.1 
!     --------------------------------------------------
!     TEST 3 
      elseif (ChooseCASE.eq.3) then 
         rs         = 45d-6                      
         drejet     = 0.11d0                     
         cini       = 350.d0                     
         vrejet     = 45.0d0                     
         rmassrejet = cini*vrejet*1d-3/1.53D0
      endif

!      ________________________________________________________
!     |                                                        |
!     |                         Profile                        |
!     |________________________________________________________|

!      _________________________________
!     |                                 |
!     |         CONSTANT profile        |
!     |_________________________________|
     
      IF (ChoosePROFILE.eq.1) THEN
         do k=1,NZ
            do i=1,N_CELL
                mask(i,k) = 1.0d0
                mask(i,k) = mask(i,k)*Heaviside(i)
            enddo
         enddo
!      _________________________________
!     |                                 |
!     |        POISEUILLE profile       |
!     |_________________________________|

      ELSEIF (ChoosePROFILE.eq.2) THEN
         do k=1,NZ
            do i=1,N_CELL                  
               mask(i,k)=(1.d0-(2.0d0*(xc(i)-0.5d0*(xRegMax+xRegMin))/&
                              (xRegMax-xRegMin))**2)&
                        *(1.d0-(2.0d0*(yc(i)-0.5d0*(yRegMax+yRegMin))/&
                              (yRegMax-yRegMin))**2)
               if (ChooseREGION.eq.2) then
                   mask(i,j) = -(((xc(i)-xcReg)/aaReg)**2 + &
                                 ((yc(j)-ycReg)/bbReg)**2 - 1.0d0) ; 
               endif
               mask(i,k) = mask(i,k)*Heaviside(i)
            enddo
         enddo
!      _________________________________
!     |                                 |
!     |           TUBE profile          |
!     |_________________________________|

      ELSEIF (ChoosePROFILE.eq.3) THEN
         do k=1,NZ
            do i=1,N_CELL
               mask(i,k)=1.d0-(2.0d0*(xc(i)-0.5d0*(xRegMax+xRegMin))/&
                              (xRegMax-xRegMin))**2
               mask(i,k) = mask(i,k)*Heaviside(i)
            enddo
         enddo
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |    Volume fraction & velocity of the initial release   |
!     |________________________________________________________|

      if (ICASE.eq.1) then
         Aini = cini/rho_s
         Winj = 0.5d0
         do k=1,NZ
            if ((k.ge.kRegMin).and.(k.le.kRegMax)) then
                do i=1,N_CELL
                   wfn(i,k)     = -Winj*mask(i,k)
                   alphasn(i,k) =  Aini*mask(i,k)
                   alphafn(i,k) =  1.0d0 - alphasn(i,k)
                enddo
            endif 
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                      Granular bed & NO-solid option                 !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |         Input of the granular bed at the bottom        |
!     |________________________________________________________|

      if (IBED.eq.1) then
         kZbed = 51
         do k=1,kZbed	
            do i=1,N_CELL
                alphasn(i,k) = 0.55d0
                alphafn(i,k) = 1.0d0 - alphasn(i,k)
            enddo
         enddo
      endif
!      ________________________________________________________
!     |                                                        |
!     |                     No solid option                    |
!     |________________________________________________________|

      if (ISOLID.eq.0) then
         do k=1,NZ
            do i=1,N_CELL
               alphasn(i,k) = 0.0d0
               alphafn(i,k) = 1.0d0 
            enddo
         enddo
      endif

!*********************************************************************!
!                                                                     !
!                             Turbulance                              !
!                                                                     !
!*********************************************************************!

!      energmin  = 1.d-6
!      dissipmin = 1.d-8

!      do k=1,NZ 
!         do i=1,N_CELL
!            energftn(i,k)   = energmin
!            dissipftn(i,k)  = dissipmin
!            energcofsn(i,k) = energmin
!            energstn(i,k)   = energftn(i,k)
!            dissipstn(i,k)  = dissipftn(i,k)
!         enddo
!      enddo


!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!
  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: initial'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                            END INITIAL                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!  
