!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   INITIAL VALUES                             !
!                             Dec 2015                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE initial(alphafn,ufn,vfn,wfn,pfn,viscof,rhof,        &
                         alphasn,usn,vsn,wsn,psn,viscos,rhos,        &
                         alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv,      &
                         alphasv,usv,vsv,wsv,psv,viscosv,rhosv,      &
                         etan,etav,Hpr,Hprv,                         &
                         xct,yct,zct,                                &
                         xvt,yvt,zvt,                                &
                         xc,yc,sig,dsig,No_cp,nbe,                   &
                         xv,yv,sigv,dsigv,No_vp,nbev,                &
                         h,hv)
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the initial condition of the main        !
!    variables of our problem.                                        !                                                                     !
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
!  | <-- walld    | N_CELL    | Wall distance for LES option       |  !
!  | <-- eddy     | N_CELL    | Eddy viscosity for LES option      |  !
!  |______________|___________|____________________________________|  !
!  | <-- xct      |(N_CELL,NZ)| coordinate y at element center     |  !
!  | <-- yct      |(N_CELL,NZ)| coordinate x at element center     |  !
!  | <-- zct      |(N_CELL,NZ)| coordinate z at element center     |  !
!  | <-- xvt      |(N_VERT,NZ-1)| coordinate y at the vertices     |  !
!  | <-- yvt      |(N_VERT,NZ-1)| coordinate x at the vertices     |  !
!  | <-- zvt      |(N_VERT,NZ-1)| coordinate z at the vertices     |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   |         Description                 |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xc      | N_CELL    | x-coordinate of the cell center     |  !
!  | --> yc      | N_CELL    | y-coordinate of the cell center     |  !
!  | --> xv      | N_VERT    | x-coordinate of the vertex          |  !
!  | --> yv      | N_VERT    | y-coordinate of the vertex          |  !
!  | --> sig     | NZ        | sigma of the cell centers           |  !
!  | --> sigv    | NZ-1      | sigma of the vertex                 |  !
!  | --> h       | N_CELL    | Depth of the domain at each cell    |  !
!  | --> hv      | N_CELL    | Depth of the domain at the vertex   |  !
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
!  | * ICASE         | Tag about the release box                   |  !
!  | * ITURBU        | Tag about the turbulence model option       |  !
!  |_________________|_____________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  | rr          | radius of a cell center to (xcReg,ycReg)        |  !
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
      real*8,dimension(:,:) :: xct,yct,zct
      real*8,dimension(:,:) :: xvt,yvt,zvt
      real*8,dimension(:)   :: xc,yc,sig,dsig
      real*8,dimension(:)   :: xv,yv,sigv,dsigv
      integer,dimension(:,:):: No_cp,No_vp
      integer,dimension(:)  :: nbe,nbev
      real*8,dimension(:)   :: h,hv
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer :: IDISPLAY
!     ----------------------------------------
      real*8  :: x,y,z,c
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
!     |                   General parameters                   |
!     |________________________________________________________|

      pa        = 0.00d00
!      Re        = 1.0000d+02
      alphasmin = 1.0000d-08
      alphasmax = 0.6250d+00
      rho_f     = 1.0000d+00
      rho_s     = 2.6500d+00
      visco_f   = 1./Re
!      ________________________________________________________
!     |                                                        |
!     |                  Density  &  Viscosity                 |
!     |________________________________________________________|

!     ------------------------------------
!     Cell center
      do k=1,NZ
         do i=1,N_CELL
            rhof(i,k)   = rho_f
            rhos(i,k)   = rho_s
            viscof(i,k) = visco_f
            viscos(i,k) = visco_f
         enddo
      enddo
!     ------------------------------------
!     Vertex
      do k=1,NZ-1
         do nv=1,N_VERT
            rhofv(nv,k) = rho_f
            rhosv(nv,k) = rho_s
            viscofv(nv,k) = visco_f
            viscosv(nv,k) = visco_f
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |         Free surface: eta & Total water depth: H       |
!     |________________________________________________________|

!     ------------------------------------
!     Cell center
      do i=1,N_CELL
         etan(i) = 0.0d0
         Hpr(i)  = etan(i) + h(i)
      enddo
!     ------------------------------------
!     Vertex
      do nv=1,N_VERT
         etav(nv) = 0.0d0
         Hprv(nv) = etav(nv) + hv(nv)
      enddo
!     ------------------------------------
!     Box case
      if (ChooseDomBox.eq.1) then
         do i=1,N_CELL
            etan(i) = 0.0d0
            h(i)    = 0.0d0
            Hpr(i)  = 1.0d0
         enddo
         do nv=1,N_VERT
            etav(nv) = 0.0d0
            hv(nv)   = 0.0d0
            Hprv(nv) = 1.0d0
         enddo
      endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(etan)
         call communication2D(Hpr)
#     endif
!     =============== END ================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                     Volume fraction                    |
!     |________________________________________________________|

!     ------------------------------------
!     Cell center
      do k=1,NZ
         do i=1,N_CELL
             alphafn(i,k) = 1.0d0
             alphasn(i,k) = 0.0d0
         enddo
      enddo
!     ------------------------------------
!     Vertex
      do k=1,NZ-1
         do nv=1,N_VERT
            alphafv(i,k) = 1.0d0
            alphasv(i,k) = 0.0d0
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                       Velocity                         |
!     |________________________________________________________|

!     ____________________________________
!     Velocity: (u,v,w)= 0
!     ------------------------------------
!     Cell center
      do i=1,N_CELL
         do k=1,NZ
            ufn(i,k) = 0.0d0
            vfn(i,k) = 0.0d0
            wfn(i,k) = 0.0d0
            usn(i,k) = ufn(i,k)
            vsn(i,k) = vfn(i,k)
            wsn(i,k) = wfn(i,k)
         enddo
      enddo
!     ------------------------------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            ufv(nv,k) = 0.0d0
            vfv(nv,k) = 0.0d0
            wfv(nv,k) = 0.0d0
            usv(nv,k) = ufv(nv,k)
            vsv(nv,k) = vfv(nv,k)
            wsv(nv,k) = wfv(nv,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                        Pressure                        |
!     |________________________________________________________|
!    Hydrostatic pressure
     IF(Ihydro .gt. 0) THEN
!     ------------------------------------
!     Cell center
      do i=1,N_CELL
         do k=1,NZ
            pfn(i,k) = pa + gra*rhof(i,k)*(1.d0-sig(k))*Hpr(i)
            psn(i,k) = pfn(i,k)
         enddo
      enddo
!     ------------------------------------
!     Vertex
      do i=1,N_VERT
         do k=1,NZ-1
            pfv(nv,k) = pa + gra*rhofv(nv,k)*(1.d0-sigv(k))*Hprv(i)
            psv(nv,k) = pfv(nv,k)
         enddo
      enddo
      ELSE
!     No hydrostatic pressure
!     ------------------------------------
!     Cell center
      do i=1,N_CELL
         do k=1,NZ
            pfn(i,k) = pa
            psn(i,k) = pfn(i,k)
         enddo
      enddo
!     ------------------------------------
!     Vertex
      do i=1,N_VERT
         do k=1,NZ-1
            pfv(nv,k) = pa
            psv(nv,k) = pfv(nv,k)
         enddo
      enddo
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |                 Coordinates (xt,yt,zt)                 |
!     |________________________________________________________|

!     --------------------------------------------------
!     Cell
      do i=1,N_CELL
         do k=1,NZ
            xct(i,k) = xc(i)
            yct(i,k) = yc(i)
            zct(i,k) = sig(k)*Hpr(i)-h(i)
         enddo
      enddo

!     --------------------------------------------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            xvt(nv,k) = xv(nv)
            yvt(nv,k) = yv(nv)
            zvt(nv,k) = sigv(k)*Hprv(nv)-hv(nv)
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                     No solid option                    |
!     |________________________________________________________|

        if(IDISPLAY .eq. 1) print*,'      SOLID OPTION = NO'
         do k=1,NZ
            do i=1,N_CELL
               alphasn(i,k) = 0.0d0
               alphafn(i,k) = 1.0d0
            enddo
         enddo

!*********************************************************************!
!                                                                     !
!                             Turbulence                              !
!                                                                     !
!*********************************************************************!
      IF (Iturbu .eq. 0) THEN
        if(IDISPLAY .eq. 1) print*,'      TURB MODEL =  NONE'
      ELSEIF (Iturbu .eq. 1) THEN
        if(IDISPLAY .eq. 1) print*,'      TURB MODEL =  LES'
      ELSEIF (Iturbu .eq. 2) THEN
        if(IDISPLAY .eq. 1) then
            print*,'      TURB MODEL =  RANS'
            print*,'      WARNING: THIS OPTION IS STILL UNDER DEVELOP'
        endif
        STOP
      ENDIF
!   ------------------------------------------------------------------
!    Pressure solver
      IF (Ipress .eq. 1) THEN
        if(IDISPLAY .eq. 1) print*,'      POISSON SOLVER  =  SOR'
      ELSEIF (Ipress .eq. 2) THEN
        if(IDISPLAY .eq. 1) print*,'      POISSON SOLVER  =  BICGSTAB'
      ELSE
        if(IDISPLAY .eq. 1) then
            print*,'      UNDEFINED POISSON SOLVER. EXIT.'
        endif
        STOP
      ENDIF
!   ------------------------------------------------------------------
!   hydrostatic pressure
      IF (Ihydro .eq. 0) THEN
        if(IDISPLAY .eq. 1) print*,'      PRESSURE  =  DYNAMIC ONLY'
      ELSEIF (Ihydro .eq. 1) THEN
        if(IDISPLAY .eq. 1) print*,'      PRESSURE  =  DYNAMIC + STATIC'
      ELSE
        if(IDISPLAY .eq. 1) then
            print*,'       PRESSURE = STATIC ONLY'
        endif
      ENDIF

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
!                      	     END INITIAL                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
