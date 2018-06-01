!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   CALCULATION OF THE ADVECTION TERM                 !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE advection3D(Am0,Am1,Am2,Am3,AmT,AmB,AmG,    &
                             uu,vv,ww,                       &
                             phi,xc,yc,sig,No_cp,nbe,        &
                             sigv,dsigv)     
 
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
!  | <--> AmT   |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  | <--> AmB   |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  ! 
!  | <--> AmG   |(N_CELL0,NZ)| Vector with Gradient terms          |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | --> uu      |(N_CELL,NZ)| Velocity component u                |  !
!  | --> vv      |(N_CELL,NZ)| Velocity component v                |  !
!  | --> ww      |(N_CELL,NZ)| Velocity component w                |  ! 
!  |_____________|___________|_____________________________________|  !
!  | --> phi     |(N_CELL,NZ)| Function phi at the center element  |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding three cells|  !
!  | --> nbe     |(N_CELL)   | Type of boundary cell (inside or bc)|  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters and variables used:                            !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |    dt      |  Time step                                       |  ! 
!  |   areaCell |  Area of the cell                                |  ! 
!  |____________|__________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,nc,i,j,k| Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
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
!  | Ui1           |(N_CELL0,NZ) Mass flux of horizontal neigh. 1  |  !
!  | Ui2           |(N_CELL0,NZ) Mass flux of horizontal neigh. 2  |  !
!  | Ui3           |(N_CELL0,NZ) Mass flux of horizontal neigh. 3  |  !
!  | UiT           |(N_CELL0,NZ) Mass flux of the vert. top neigh. |  !
!  | UiB           |(N_CELL0,NZ) Mass flux of the vert. bott. neigh|  !
!  |_______________|_______________________________________________|  !
!  | aa0,aa(1:3)   | Auxiliar horizontal matrix coefficients       |  !
!  | aaT,aaB       | Auxiliar vertical matrix coefficients         |  !
!  | Uij(1:3)      | Auxiliar mass flux values of the neighbors    |  !
!  | Uipos,Uineg   | Positive and negative values of the mass flux |  !
!  | GF            | Gradient (rhs contribution)                   |  ! 
!  | GFi,GFj       | Gradient at the cell center & neigborn        |  !
!  | GFT,GFB       | Gradient at the top and bottom                |  !
!  | GFpos,GFneg   | Positive and negative values of the gradient  |  !
!  |_______________|_______________________________________________|  !
!  | jc            | neighborn loop index                          |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |                    -  grandientLSM                            |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
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

      real*8, dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmG(N_CELL0,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: uu(N_CELL,NZ)
      real*8, dimension(:,:) :: vv(N_CELL,NZ)
      real*8, dimension(:,:) :: ww(N_CELL,NZ)
      real*8, dimension(:,:) :: phi(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     --------------------------------------
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8 ,dimension(:,:),allocatable :: dphidx
      real*8 ,dimension(:,:),allocatable :: dphidy
      real*8 ,dimension(:,:),allocatable :: dphidsig
!     --------------------------------------
      real*8, dimension(:,:),allocatable :: Ui1,Ui2,Ui3,UiT,UiB
      real*8 ,dimension(:) :: Uij(1:3)
      real*8 :: nxArea,nyArea,nzArea
      real*8 :: u0j,v0j,w0B,w0T
      real*8 :: Uineg,Uipos
!     --------------------------------------  
      real*8 ,dimension(:,:),allocatable :: FluxLimiter
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
      real*8 :: GF,GFi,GFj,GFT,GFB
      real*8 :: GFpos,GFneg
      integer:: jc

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Advection3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(dphidx(N_CELL,NZ),       &
               dphidy(N_CELL,NZ),       &
               dphidsig(N_CELL,NZ),     &
               Ui1(N_CELL0,NZ),         &
               Ui2(N_CELL0,NZ),         &
               Ui3(N_CELL0,NZ),         &
               UiT(N_CELL0,NZ),         &
               UiB(N_CELL0,NZ),         &
               FluxLimiter(N_CELL,NZ))

!*********************************************************************!
!                                                                     !
!                               Gradient                              !
!                                                                     !
!*********************************************************************!

      call grandientLSM(dphidx,dphidy,dphidsig,&
                        phi,No_cp,nbe,sig) 

!*********************************************************************!
!                                                                     !
!                          Mass flux calculation                      !
!                                                                     !
!*********************************************************************!

      do k=1,NZ  
         do i=1,N_CELL0	
            Ui1(i,k) = 0.0d0    
            Ui2(i,k) = 0.0d0
            Ui3(i,k) = 0.0d0
            UiT(i,k) = 0.0d0
            UiB(i,k) = 0.0d0
         enddo
      enddo

      DO k=2,NZ-1  
         do i=1,N_CELL0
!           ___________________________________________________
!           Horizontal neighbors
            do j=1,3
!              ---------------------------------
!              Index of the neighbor cell-center
               jc = No_cp(i,j)
!              ---------------------------------
!              Velocity at the faces 
               u0j = 0.5d0*(uu(i,k)+uu(jc,k))
               v0j = 0.5d0*(vv(i,k)+vv(jc,k))
!              ---------------------------------
!              Mass fluxes of the cell
               nxArea =  dyVV(i,j)*dsigv(k-1)
               nyArea = -dxVV(i,j)*dsigv(k-1)
               Uij(j) = u0j*nxArea+v0j*nyArea
            enddo
!           ------------------------------------
!           Assignation of mass fluxes C
            Ui1(i,k) = Uij(1)    
            Ui2(i,k) = Uij(2)
            Ui3(i,k) = Uij(3)
!           ___________________________________________________
!           Top & bottom neighbors  
!           ---------------------------------
!           Velocity at the faces               
            w0T = 0.5d0*(ww(i,k)+ww(i,k+1))
            w0B = 0.5d0*(ww(i,k)+ww(i,k-1))
!           ------------------------------------
!           Assignation of mass fluxes C
            nzArea   =  areaCell(i)
            UiT(i,k) =  w0T*nzArea
            nzArea   = -areaCell(i)
            UiB(i,k) =  w0B*nzArea
        enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                       Flux limiter calculation                      !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                 Using the LED technique                |
!     |________________________________________________________|

#     ifdef KeyFluxLimiter
         DO k=2,NZ-1
            do i=1,N_CELL0	              
!              -------------------------------
!              Maximum & minimum value of phi	
               phimax = phi(i,k)
               phimin = phi(i,k)
               do j=1,3
                  jc = No_cp(i,j)
                  phimax = max(phimax,phi(jc,k))
                  phimin = min(phimin,phi(jc,k))
               enddo
!              -------------------------------
!              Flux limiter calculation 
               FluxLimiterMin = 1.0d5
               do j=1,3 
                  GFi =  dphidx(i,k)*(xme(i,j)-xc(i)) &
                       + dphidy(i,k)*(yme(i,j)-yc(i))
                  GFpos = 0.5d0*(GFi+abs(GFi))
                  GFneg = 0.5d0*(GFi-abs(GFi))

	          if(dabs(GFi).lt.1.0d-10) then
	             FluxLimiterij = 1.0d0
	          else
	             d1 = (phimax-phi(i,k))/GFi
		     d2 = (phimin-phi(i,k))/GFi
		     FluxLimiterpos = min(1.0d0,d1)
		     FluxLimiterneg = min(1.0d0,d2)
		     FluxLimiterij  =  FluxLimiterpos*GFpos &
                                     + FluxLimiterneg*GFneg
	          endif
	          FluxLimiterMin = min(FluxLimiterMin,FluxLimiterij)
	       enddo	
	       FluxLimiter(i,k) = FluxLimiterMin	  
            enddo
         ENDDO
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                    No Flux limiter                     |
!     |________________________________________________________|

#     ifndef KeyFluxLimiter 
         do k=1,NZ
            do i=1,N_CELL	
               FluxLimiter(i,k) = 1.0d0
            enddo
         enddo
#     endif


!*********************************************************************!
!                                                                     !
!                       Advection contributions                       !
!                                                                     !
!*********************************************************************!

      DO k=2,NZ-1  
         do i=1,N_CELL0	
!            __________________________________________________
!           |                                                  |
!           |                  Matrix & rhs                    |
!           |__________________________________________________|

!           ___________________________________________________
!           Horizontal neighbors  
!           ------------------------------------
!           Assignation of mass fluxes C
            Uij(1) = Ui1(i,k)  
            Uij(2) = Ui2(i,k)
            Uij(3) = Ui3(i,k)
           
            aa0 = 0.0d0
            GF  = 0.0d0
            do j=1,3
!              ---------------------------------
!              Number of the neighbor cell-center 
	       jc = No_cp(i,j)
!              ---------------------------------
!              Mass fluxes of the cell
	       Uipos =  0.5d0*(Uij(j)+abs(Uij(j)))
	       Uineg = -0.5d0*(Uij(j)-abs(Uij(j)))
!              ---------------------------------
!              Matrix coeffients                	
	       aa0   =  Uipos + aa0 
	       aa(j) = -Uineg               
!              ---------------------------------
!              Gradient at of the cell & the neighbor
	       GFi =  dphidx(i,k)*(xme(i,j)-xc(i)) &
                    + dphidy(i,k)*(yme(i,j)-yc(i))
	 
	       GFj =  dphidx(jc,k)*(xme(i,j)-xc(jc)) &
                    + dphidy(jc,k)*(yme(i,j)-yc(jc))
!              ---------------------------------
!              Gradient coefficient (rhs)
               GF = GF + Uipos*GFi*FluxLimiter(i,k)  &
                       - Uineg*GFj*FluxLimiter(jc,k) 
	    enddo
!           ___________________________________________________
!           Top neighbor
!           ---------------------------------
!           Mass fluxes of the cell
	    Uipos =  0.5d0*(UiT(i,k)+abs(UiT(i,k)))
	    Uineg = -0.5d0*(UiT(i,k)-abs(UiT(i,k)))
!           ---------------------------------
!           Matrix coeffients                	
	    aa0  =  Uipos + aa0 
	    aaT  = -Uineg               
!           ---------------------------------
!           Gradient at of the cell & the neighbor
	    GFi = dphidsig(i,k)*(sigv(k)-sig(k))
            GFT = dphidsig(i,k+1)*(sigv(k)-sig(k+1))
!           ---------------------------------
!           Gradient coefficient (rhs)
            GF = GF + Uipos*GFi-Uineg*GFT 
!           ___________________________________________________
!           Bottom neighbor
!           ---------------------------------
!           Mass fluxes of the cell
	    Uipos =  0.5d0*(UiB(i,k)+abs(UiB(i,k)))
	    Uineg = -0.5d0*(UiB(i,k)-abs(UiB(i,k)))
!           ---------------------------------
!           Matrix coeffients                	
	    aa0   =  Uipos + aa0 
	    aaB   = -Uineg               
!           ---------------------------------
!           Gradient at of the cell & the neighbor
	    GFi = dphidsig(i,k)*(sigv(k-1)-sig(k))	 
	    GFB = dphidsig(i,k-1)*(sigv(k-1)-sig(k-1))
!           ---------------------------------
!           Gradient coefficient (rhs)
            GF = GF + Uipos*GFi-Uineg*GFB 
!            __________________________________________________
!           |                                                  |
!           |                  Contributions                   |
!           |__________________________________________________|

!           --------------------------------------
!           Matrix
	    Am0(i,k) = Am0(i,k) + aa0    
            Am1(i,k) = Am1(i,k) + aa(1)   
            Am2(i,k) = Am2(i,k) + aa(2)   
            Am3(i,k) = Am3(i,k) + aa(3)
            AmT(i,k) = AmT(i,k) + aaT             
            AmB(i,k) = AmB(i,k) + aaB
!           --------------------------------------
!           Known part to rhs
	    AmG(i,k) = AmG(i,k) + GF    
        enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(dphidx,dphidy,dphidsig,           &
                 Ui1,Ui2,Ui3,UiT,UiB,              &
                 FluxLimiter)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Advection3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            CALCULATION OF THE ADVECTION TERM FOR THE VELOCITY       !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE advectionVelocity(Am0,Am1,Am2,Am3,AmT,AmB,AmG,    &
                                   uu,vv,ww,pfn,                   &
                                   phi,xc,yc,sig,No_cp,nbe,        &
                                   phiv,xv,yv,sigv,dsigv,No_vp,nbev)     
 
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
!  | <--> AmT   |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  | <--> AmB   |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  ! 
!  | <--> AmG   |(N_CELL0,NZ)| Vector with Gradient terms          |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | --> uu      |(N_CELL,NZ)| Velocity component u                |  !
!  | --> vv      |(N_CELL,NZ)| Velocity component v                |  !
!  | --> ww      |(N_CELL,NZ)| Velocity component w                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> pfn     |(N_CELL,NZ)| Pressure variable at the old time   |  ! 
!  |_____________|___________|_____________________________________|  !
!  | --> phi     |(N_CELL,NZ)| Function phi at the center element  |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding three cells|  !
!  | --> nbe     |(N_CELL)   | Type of boundary cell (inside or bc)|  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters and variables used:                            !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |    dt      |  Time step                                       |  ! 
!  |   areaCell |  Area of the cell                                |  ! 
!  |____________|__________________________________________________|  !
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
!  | Ui1           |(N_CELL0,NZ) Mass flux of horizontal neigh. 1  |  !
!  | Ui2           |(N_CELL0,NZ) Mass flux of horizontal neigh. 2  |  !
!  | Ui3           |(N_CELL0,NZ) Mass flux of horizontal neigh. 3  |  !
!  | UiT           |(N_CELL0,NZ) Mass flux of the vert. top neigh. |  !
!  | UiB           |(N_CELL0,NZ) Mass flux of the vert. bott. neigh|  !
!  |_______________|_______________________________________________|  !
!  | aa0,aa(1:3)   | Auxiliar horizontal matrix coefficients       |  !
!  | aaT,aaB       | Auxiliar vertical matrix coefficients         |  !
!  | Uij(1:3)      | Auxiliar mass flux values of the neighbors    |  !
!  | Uipos,Uineg   | Positive and negative values of the mass flux |  !
!  | GF            | Gradient (rhs contribution)                   |  ! 
!  | GFi,GFj       | Gradient at the cell center & neigborn        |  !
!  | GFT,GFB       | Gradient at the top and bottom                |  !
!  | GFpos,GFneg   | Positive and negative values of the gradient  |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |                    -  grandientLSM                            |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
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

      real*8, dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmG(N_CELL0,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: uu(N_CELL,NZ)
      real*8, dimension(:,:) :: vv(N_CELL,NZ)
      real*8, dimension(:,:) :: ww(N_CELL,NZ)
      real*8, dimension(:,:) :: pfn(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     ----------------------------------------
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
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

      real*8 ,dimension(:,:),allocatable :: dphidx
      real*8 ,dimension(:,:),allocatable :: dphidy
      real*8 ,dimension(:,:),allocatable :: dphidsig
      real*8,dimension(:,:),allocatable  :: dpdx,dpdy,dpdz
!     --------------------------------------
      real*8, dimension(:,:),allocatable :: Ui1,Ui2,Ui3,UiT,UiB
      real*8 ,dimension(:) :: Uij(1:3)
      real*8 :: nxArea,nyArea,nzArea
      real*8 :: u0j,v0j,w0B,w0T
      real*8 :: Uineg,Uipos
      real*8 :: uRC,vRC,wRC,A0i,A0j,A0T,A0B
      real*8 :: Vol0i,Vol0j,Vol0T,Vol0B
!     --------------------------------------  
      real*8 ,dimension(:,:),allocatable :: FluxLimiter
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
      real*8 :: GF,GFi,GFj,GFT,GFB
      real*8 :: GFpos,GFneg
      integer:: jc,jj,nv1,nv2
!     --------------------------------------  
      integer, parameter :: ChooseMasswithP   = 0
      integer, parameter :: ChooseMassKimChoi = 0
      integer, parameter :: ChooseVeloKimChoi = 0                     
!     --------------------------------------  

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Advection3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(dphidx(N_CELL,NZ),       &
               dphidy(N_CELL,NZ),       &
               dphidsig(N_CELL,NZ),     &
               Ui1(N_CELL0,NZ),         &
               Ui2(N_CELL0,NZ),         &
               Ui3(N_CELL0,NZ),         &
               UiT(N_CELL0,NZ),         &
               UiB(N_CELL0,NZ),         &
               FluxLimiter(N_CELL,NZ))

      allocate(dpdx(N_CELL,NZ),dpdy(N_CELL,NZ),dpdz(N_CELL,NZ))

!*********************************************************************!
!                                                                     !
!                               Gradients                             !
!                                                                     !
!*********************************************************************!

!     _________________________________________________
!     Gradient of the velocity component
      call grandientLSM(dphidx,dphidy,dphidsig,&
                        phi,No_cp,nbe,sig)

!     _________________________________________________
!     Gradient of the pressure p^(n)
      call grandientLSM(dpdx,dpdy,dpdz,&
                        pfn,No_cp,nbe,sig)

!*********************************************************************!
!                                                                     !
!                          Mass flux calculation                      !
!                                                                     !
!*********************************************************************!

      do k=1,NZ  
         do i=1,N_CELL0	
            Ui1(i,k) = 0.0d0    
            Ui2(i,k) = 0.0d0
            Ui3(i,k) = 0.0d0
            UiT(i,k) = 0.0d0
            UiB(i,k) = 0.0d0
         enddo
      enddo

      DO k=2,NZ-1  
         do i=1,N_CELL0	
!           ___________________________________________________
!           Horizontal neighbors    
            do j=1,3
!              ---------------------------------
!              Index of the neighbor cell-center  
	           jc = No_cp(i,j)	           
!              ---------------------------------
!              Velocity at the faces 
               u0j = 0.5d0*(uu(i,k)+uu(jc,k)) 
               v0j = 0.5d0*(vv(i,k)+vv(jc,k)) 
               
!              ---------------------------------
!              Rhie & Chow interpolation 
               IF (ChooseMasswithP.eq.1) THEN
                  !A0i = Am0(i,k)*AreaCell(i)*dsigv(k-1)/dt
                  !A0j = Am0(jc,k)*AreaCell(jc)*dsigv(k-1)/dt
                  A0i = AreaCell(i)*dsigv(k-1)/dt
                  A0j = AreaCell(jc)*dsigv(k-1)/dt                  
                  if (jc.gt.N_CELL0) A0j = A0i               
!                 -------------
                  if ((A0i.eq.0).or.(A0j.eq.0)) then
                     print*,'Error: Am(i,k) = 0 at one side'
                     stop
                  endif
!                 -------------
                  uRC = 0.5d0*(dpdx(i,k)/A0i + dpdx(jc,k)/A0j) &
                       -2.0d0*dyVV(i,j)*(pfn(jc,k)-pfn(i,k))/(A0i+A0j)
                  vRC = 0.5d0*(dpdy(i,k)/A0i+dpdy(jc,k)/A0j)  &
                       +2.0d0*dxVV(i,j)*(pfn(jc,k)-pfn(i,k))/(A0i+A0j)
!                 -------------
                  u0j = u0j + uRC
                  v0j = v0j + vRC
               ENDIF
               
!              >>>>>>>>>>>>>>>>>>>>>>>
!              Mass average Kim Choi              
               IF (ChooseMassKimChoi.eq.100) THEN
                  Vol0i = AreaCell(i)*dsigv(k-1)
                  Vol0j = AreaCell(jc)*dsigv(k-1)                  
                  if (jc.gt.N_CELL0) Vol0j = Vol0i 
                  uRC = 0.5d0*doe(i,j)*(dphidx(i,k)+dphidx(jc,k))
                  vRC = 0.5d0*doe(i,j)*(dphidy(i,k)+dphidy(jc,k))
                  u0j = (uu(i,k)*Vol0j+uu(jc,k)*Vol0i)/(Vol0i+Vol0j)!+uRC
                  v0j = (vv(i,k)*Vol0j+vv(jc,k)*Vol0i)/(Vol0i+Vol0j)!+vRC
               ENDIF
!              >>>>>>>>>>>>>>>>>>>>>>>

!              Mass fluxes of the cell
               nxArea =  dyVV(i,j)*dsigv(k-1)
               nyArea = -dxVV(i,j)*dsigv(k-1)
               Uij(j) = u0j*nxArea+v0j*nyArea
            enddo
!           ------------------------------------
!           Assignation of mass fluxes C
            Ui1(i,k) = Uij(1)    
            Ui2(i,k) = Uij(2)
            Ui3(i,k) = Uij(3)
!           ___________________________________________________
!           Top neighbor  
!           ------------------------------------
!           Velocity at the faces               
            w0T = 0.5d0*(ww(i,k)+ww(i,k+1)) 
!           ------------------------------------
!           Rhie & Chow interpolation 
            IF (ChooseMasswithP.eq.1) THEN
               A0i = Am0(i,k)*AreaCell(i)*dsigv(k-1)/dt
               A0T = Am0(i,k+1)*AreaCell(i)*dsigv(k)/dt
               if (k.ge.NZ-1) A0T = A0i
!              -----------
               if (A0T.eq.0) then
                  print*,'Error Am(i,k) = 0 at top' 
                  stop
               endif
!              ------------
               wRC = 0.5d0*(dpdz(i,k)/A0i + dpdz(i,k+1)/A0T) &
                     -2.0d0/(A0i+A0T)*(pfn(i,k+1)-pfn(i,k))/ &
                     (sig(k+1)-sig(k))
!              -------------
               w0T = w0T + wRC
            ENDIF
            
!           >>>>>>>>>>>>>>>>>>>>>>>
!           Mass average Kim Choi                
            IF (ChooseMassKimChoi.eq.100) THEN               
                  Vol0i = AreaCell(i)*dsigv(k-1)
                  Vol0T = AreaCell(i)*dsigv(k)
                  if (k.ge.NZ-1) Vol0T = Vol0i
                  w0T = (ww(i,k)*Vol0T+ww(i,k+1)*Vol0i)/(Vol0i+Vol0T)                
            ENDIF            
!           >>>>>>>>>>>>>>>>>>>>>>>

!           Assignation of mass fluxes C
            nzArea   =  areaCell(i)
            UiT(i,k) =  w0T*nzArea
!           ___________________________________________________
!           Bottom neighbor
!           ------------------------------------
!           Velocity at the faces               
            w0B = 0.5d0*(ww(i,k)+ww(i,k-1)) 
!           ------------------------------------
!           Rhie & Chow interpolation 
            IF (ChooseMasswithP.eq.1) THEN
               A0i = Am0(i,k)*AreaCell(i)*dsigv(k-1)/dt
               A0B = Am0(i,k-1)*AreaCell(i)*dsigv(k-1)/dt !<<<<<Check k=1 case
               if (k.le.2) A0B = A0i
!              ------------
               if (A0B.eq.0) then
                  print*,'Error Am(i,k) = 0 at bottom' 
                  stop
               endif
!              ------------
               wRC = 0.5d0*(dpdz(i,k)/A0i + dpdz(i,k-1)/A0B) &
                     -2.0d0/(A0i+A0B)*(pfn(i,k)-pfn(i,k-1))/&
                     (sig(k)-sig(k-1))
!              ------------
               w0B = w0B + wRC
            ENDIF
            
!           >>>>>>>>>>>>>>>>>>>>>>>
!           Mass average Kim & Choi                
            IF (ChooseMassKimChoi.eq.1) THEN 
               Vol0i = AreaCell(i)*dsigv(k-1)
               if (k.eq.2) then              
                  Vol0B = AreaCell(i)*dsigv(1)
               else
                  Vol0B = AreaCell(i)*dsigv(k-2)
               endif
               if (k.le.2) Vol0B = Vol0i
               w0B = (ww(i,k)*Vol0B+ww(i,k-1)*Vol0i)/(Vol0i+Vol0B)                
            ENDIF              
!           >>>>>>>>>>>>>>>>>>>>>>>

!           ------------------------------------
!           Assignation of mass fluxes C
            nzArea   = -areaCell(i)
            UiB(i,k) =  w0B*nzArea
        enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                       Flux limiter calculation                      !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                 Using the LED technique                |
!     |________________________________________________________|

#     ifdef KeyFluxLimiter
         DO k=2,NZ-1
            do i=1,N_CELL0	              
!              -------------------------------
!              Maximum & minimum value of phi	
               phimax = phi(i,k)
               phimin = phi(i,k)
               do j=1,3
                  jc = No_cp(i,j)
                  phimax = max(phimax,phi(jc,k))
                  phimin = min(phimin,phi(jc,k))
               enddo
!              -------------------------------
!              Flux limiter calculation 
               FluxLimiterMin = 1.0d5	
               do j=1,3 
                  GFi =  dphidx(i,k)*(xme(i,j)-xc(i)) &
                       + dphidy(i,k)*(yme(i,j)-yc(i))
                  GFpos = 0.5d0*(GFi+abs(GFi))
                  GFneg = 0.5d0*(GFi-abs(GFi))

	           if(dabs(GFi).lt.1.0d-10) then
	              FluxLimiterij = 1.0d0
	           else
	              d1 = (phimax-phi(i,k))/GFi
	              d2 = (phimin-phi(i,k))/GFi
	              FluxLimiterpos = min(1.0d0,d1)
	              FluxLimiterneg = min(1.0d0,d2)
	              FluxLimiterij  =  FluxLimiterpos*GFpos &
                                  + FluxLimiterneg*GFneg
	           endif
	           FluxLimiterMin = min(FluxLimiterMin,FluxLimiterij)
	        enddo	
	        FluxLimiter(i,k) = FluxLimiterMin	  
            enddo
         ENDDO
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                    No Flux limiter                     |
!     |________________________________________________________|

#     ifndef KeyFluxLimiter 
         do k=1,NZ
            do i=1,N_CELL	
               FluxLimiter(i,k) = 1.0d0
            enddo
         enddo
#     endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(FluxLimiter)
#     endif	
!     =============== END ================    
!     ====================================

!*********************************************************************!
!                                                                     !
!                       Advection contributions                       !
!                                                                     !
!*********************************************************************!

      DO k=2,NZ-1  
         do i=1,N_CELL0	
!            __________________________________________________
!           |                                                  |
!           |                  Matrix & rhs                    |
!           |__________________________________________________|

!           ___________________________________________________
!           Horizontal neighbors  
!           ------------------------------------
!           Assignation of mass fluxes C
            Uij(1) = Ui1(i,k)  
            Uij(2) = Ui2(i,k)
            Uij(3) = Ui3(i,k)
           
            aa0 = 0.0d0
            GF  = 0.0d0
            do j=1,3
!              ---------------------------------
!              Number of the neighbor cell-center 
               jc = No_cp(i,j)
!              ---------------------------------
!              Positive & negative mass fluxes 
!              ---------------------------------
               Uipos =  0.5d0*(Uij(j)+abs(Uij(j)))
               Uineg = -0.5d0*(Uij(j)-abs(Uij(j)))
!              ---------------------------------
!              Matrix coeffients                	
               aa0   =  Uipos + aa0 
               aa(j) = -Uineg               
!              ---------------------------------
!              Gradient at of the cell & the neighbor
               GFi =  dphidx(i,k)*(xme(i,j)-xc(i)) &
                    + dphidy(i,k)*(yme(i,j)-yc(i))
	 
               GFj =  dphidx(jc,k)*(xme(i,j)-xc(jc)) &
                    + dphidy(jc,k)*(yme(i,j)-yc(jc))
!              ---------------------------------
!              Gradient coefficient (rhs)
               GF = GF + Uipos*GFi*FluxLimiter(i,k)  &
                       - Uineg*GFj*FluxLimiter(jc,k) 
	        enddo	        
!           ___________________________________________________
!           Top neighbor
!           ---------------------------------
!           Mass fluxes of the cell
            Uipos =  0.5d0*(UiT(i,k)+abs(UiT(i,k)))
            Uineg = -0.5d0*(UiT(i,k)-abs(UiT(i,k)))
!           ---------------------------------
!           Matrix coeffients                	
	        aa0  =  Uipos + aa0 
	        aaT  = -Uineg               
!           ---------------------------------
!           Gradient at of the cell & the neighbor
	        GFi = dphidsig(i,k)*(sigv(k)-sig(k))
            GFT = dphidsig(i,k+1)*(sigv(k)-sig(k+1))
!           ---------------------------------
!           Gradient coefficient (rhs)
            GF = GF + Uipos*GFi-Uineg*GFT 
!           ___________________________________________________
!           Bottom neighbor
!           ---------------------------------
!           Mass fluxes of the cell
	        Uipos =  0.5d0*(UiB(i,k)+abs(UiB(i,k)))
	        Uineg = -0.5d0*(UiB(i,k)-abs(UiB(i,k)))
!           ---------------------------------
!           Matrix coeffients                	
	        aa0   =  Uipos + aa0 
	        aaB   = -Uineg               
!           ---------------------------------
!           Gradient at of the cell & the neighbor
	        GFi = dphidsig(i,k)*(sigv(k-1)-sig(k))	 
	        GFB = dphidsig(i,k-1)*(sigv(k-1)-sig(k-1))
!           ---------------------------------
!           Gradient coefficient (rhs)
            GF = GF + Uipos*GFi-Uineg*GFB 
!            __________________________________________________
!           |                                                  |
!           |         Face normal velocity contributions       |
!           |                    Kim & Choi                    |
!           |__________________________________________________|

            IF (ChooseVeloKimChoi.eq.1) THEN
            
            aa0 = 0.0d0
            GF  = 0.0d0
!           ___________________________________________________
!           Vertical                         
            do j=1,3 
               jc = No_cp(i,j) 
               jj = j+1
               if (jj.gt.3) jj=jj-3
               nv1 = No_vp(i,j)
               nv2 = No_vp(i,jj)                            
               Vol0i = AreaCell(i)*dsigv(k-1)
               Vol0j = AreaCell(jc)*dsigv(k-1)
               if (jc.gt.N_CELL0) Vol0j = Vol0i 
!              -----------------
!              Matrix coeffients                
               aa0   = Uij(j)*Vol0j/(Vol0i+Vol0j) + aa0
               aa(j) = Uij(j)*Vol0i/(Vol0i+Vol0j) 
!              -----------------
!              Gradient contribution
               !GF = GF + Uij(j)*0.5d0*doe(i,j)*(dphidx(i,k)+dphidx(jc,k))
               GF = GF + Uij(j)*doe(i,j)*(phiv(nv2,k)-phiv(nv1,k))/dlVV(i,j)               
            enddo
!           ___________________________________________________
!           Top             
            Vol0i = AreaCell(i)*dsigv(k-1)
            Vol0T = AreaCell(i)*dsigv(k)
            if (k.ge.NZ-1) Vol0T = Vol0i
	        aa0  =  UiT(i,k)*Vol0T/(Vol0i+Vol0T) + aa0 
	        aaT  =  UiT(i,k)*Vol0i/(Vol0i+Vol0T)
!           ___________________________________________________
!           Bottom                	
            Vol0i = AreaCell(i)*dsigv(k-1)
            if (k.eq.2) then
               Vol0B = AreaCell(i)*dsigv(1)
            else
               Vol0B = AreaCell(i)*dsigv(k-2)
            endif
            if (k.le.2) Vol0B = Vol0i
	        aa0  =  UiB(i,k)*Vol0B/(Vol0i+Vol0B) + aa0 
	        aaB  =  UiB(i,k)*Vol0i/(Vol0i+Vol0B)

            ENDIF            
!            __________________________________________________
!           |                                                  |
!           |                  Contributions                   |
!           |__________________________________________________|

!           --------------------------------------
!           Matrix
            Am0(i,k) = Am0(i,k) + aa0    
            Am1(i,k) = Am1(i,k) + aa(1)   
            Am2(i,k) = Am2(i,k) + aa(2)   
            Am3(i,k) = Am3(i,k) + aa(3)
            AmT(i,k) = AmT(i,k) + aaT             
            AmB(i,k) = AmB(i,k) + aaB
!           --------------------------------------
!           Known part to rhs
            AmG(i,k) = AmG(i,k) + GF    
        enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(dphidx,dphidy,dphidsig,           &
                 dpdx,dpdy,dpdz,                   &
                 Ui1,Ui2,Ui3,UiT,UiB,              &
                 FluxLimiter)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Advection3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	     END OF ADVECTION                         !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
