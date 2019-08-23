!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     DEFINITION OF THE VARIABLES                     !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

	MODULE variables

!---------------------------------------------------------------------!
!                                                                     !
!      SUBROUTINES:     -  alloc_variables                            !
!                       -  dealloc_varaibles                          !
!                                                                     !
!---------------------------------------------------------------------!

	implicit none

!      ________________________________________________________
!     |                                                        |
!     |               Point location at time                   |
!     |________________________________________________________|

      real*8,dimension(:,:),allocatable :: &
	     xct, &         ! cell center x-coordinate location at time t
	     yct, &         ! cell center y-coordinate location at time t
	     zct            ! cell center z-coordinate location at time t
!     ---------------------------------------------------------
!     At the vertex points
      real*8,dimension(:,:),allocatable :: &
	     xvt, &         ! vertex x-coordinate location at time t
	     yvt, &         ! vertex y-coordinate location at time t
	     zvt            ! vertex z-coordinate location at time t
!      ________________________________________________________
!     |                                                        |
!     |                  FLuid variables                       |
!     |________________________________________________________|

        real*8,dimension(:,:),allocatable :: &
             alphafnp, &     ! Fluid control volume at(n+1)
             alphafn,  &     ! Fluid control volume at (n)
             alphaf,   &     ! Fluid control volume at RK step
             ufnp,     &     ! Velocity component u_f at (n+1)
             ufn,      &     ! Velocity component u_f at (n)
             uf,       &     ! Velocity component u_f at RK step
             suf,      &     ! Velocity component u_f stats
             ufBC,     &     ! Velocity component u_f at (n+1) at Boundary
             vfnp,     &     ! Velocity component v_f at (n+1)
             vfn,      &     ! Velocity component v_f at (n)
             vf,       &     ! Velocity component v_f at RK step
             svf,      &     ! Velocity component v_f stats
             vfBC,     &     ! Velocity component v_f at (n+1) at Boundary
             wfnp,     &     ! Velocity component w_f at (n+1)
             wfn,      &     ! Velocity component w_f at (n)
             wf,       &     ! Velocity component w_f at RK step
             swf,      &     ! Velocity component w_f stats
             wfBC,     &     ! Velocity component w_f at (n+1) at Boundary
             pfnp,     &     ! Pressure of the fluid at (n+1)
             pfn,      &     ! Pressure of the fluid at (n)
             pf,       &     ! Pressure of the fluid at RK step
             spf,      &     ! Pressure of the fluid stats
             Dpfnp,    &     ! Difference of pressure of the fluid at (n+1)
             Dpfn,     &     ! Difference of pressure of the fluid at (n)
             Dpf,      &     ! Difference of pressure of the fluid at RK step
             omenpuf,  &     ! omega at (n+1) for u_f
             omenpvf,  &     ! omega at (n+1) for v_f
             omenpwf,  &     ! omega at (n+1) for w_f
             rhof,     &     ! Density of the fluid
             viscof          ! Viscosity of the fluid
!     ---------------------------------------------------------
!     At the vertex points
      real*8,dimension(:,:),allocatable :: &
             alphafv,  &     ! Fluid control volume at vertex
             ufv,      &     ! Velocity component u_f at vertex
             vfv,      &     ! Velocity component v_f at vertex
             wfv,      &     ! Velocity component w_f at vertex
             pfv,      &     ! Pressure of the fluid at vertex
             Dpfv,     &     ! Difference of pressure of the fluid at vertex
             rhofv,    &     ! Density of the fluid vertex
             viscofv         ! Viscosity of the fluid vertex
!      ________________________________________________________
!     |                                                        |
!     |                  Solid variables                       |
!     |________________________________________________________|

      real*8,dimension(:,:),allocatable :: &
             alphasnp, &     ! Solid control volume at(n+1)
             alphasn,  &     ! Solid control volume at (n)
             alphas,   &     ! Solid control volume at RK step
             usnp,     &     ! Velocity component u_s at (n+1)
             usn,      &     ! Velocity component u_s at (n)
             us,       &     ! Velocity component u_s at RK step
             vsnp,     &     ! Velocity component v_s at (n+1)
             vsn,      &     ! Velocity component v_s at (n)
             vs,       &     ! Velocity component v_s at RK step
             wsnp,     &     ! Velocity component w_s at (n+1)
             wsn,      &     ! Velocity component w_s at (n)
             ws,       &     ! Velocity component w_s at RK step
             psnp,     &     ! Pressure of the solid at (n+1)
             psn,      &     ! Pressure of the solid at (n)
             ps,       &     ! Pressure of the solid at RK step
             omenpus,  &     ! omega at (n+1) for u_s
             omenpvs,  &     ! omega at (n+1) for v_s
             omenpws,  &     ! omega at (n+1) for w_s
             rhos,     &     ! Density of the solid
             viscos          ! Viscosity of the solid
!     ---------------------------------------------------------
!     At the vertex points
      real*8,dimension(:,:),allocatable :: &
             alphasv,  &     ! Solid control volume at vertex
             usv,      &     ! Velocity component u_s at vertex
             vsv,      &     ! Velocity component v_s at vertex
             wsv,      &     ! Velocity component w_s at vertex
             psv,      &     ! Pressure of the solid at vertex
             rhosv,    &     ! Density of the solid vertex
             viscosv         ! Viscosity of the solid vertex
!      ________________________________________________________
!     |                                                        |
!     |             Free surface & bottom levels               |
!     |________________________________________________________|

      real*8,dimension(:),allocatable :: &
             Hpr,      &     ! (N_CELL) Total water depth: H = h + eta
             h,        &     ! (N_CELL) The bottom depth
             zbc,      &     ! (N_CELL) Bottom levels at each cell
             etanp,    &     ! (N_CELL) Free surface level eta at (n+1)
             etan,     &     ! (N_CELL) Free surface level eta at (n)
             HprBC,    &     ! (N_CELL) Free surface level eta at (n+1) at Boundary
             eta,      &     ! (N_CELL) Free surface level eta at RK step
             wfsurf,   &     ! (N_CELL) Velocity component w_f at the free surface
             wssurf          ! (N_CELL) Velocity component w_s at the free surface
!     ---------------------------------------------------------
!     At the vertex points
      real*8,dimension(:),allocatable :: &
             Hprv,     &     ! (N_VERT) Total water depth at vertex
             hv,       &     ! (N_VERT) Depth of the river at vertex
             zbv,      &     ! (N_VERT) Bottom levels at each vertex
             etav,     &     ! (N_VERT) Free surface level eta at the vertices
             wfsurfv,  &     ! (N_VERT) w_f at the free surface vertices
             wssurfv         ! (N_VERT) w_s at the free surface vertices

!      ________________________________________________________
!     |                                                        |
!     |              Geometry (unstructure basic)              |
!     |________________________________________________________|

      integer,dimension(:,:),allocatable :: &
              No_vp,   &     ! (N_CELL0,3) Numbering of the vertices for each cell
              No_cp          ! (N_CELL,3)  Numbering of the surrounding three cells
      integer,dimension(:),allocatable :: &
              nbe,     &     ! (N_CELL0) Assigned tag of the type of cell
                             !          = 0  without any boundary vertex
                             !          = 1  with a wall boundary vertex
                             !          = 2  with a water level boundary  vertex
                             !          = 3  with a discharge normal boundary v.
              nbev           ! (N_VERT)  Assigned tag of the type of vertex:
                             !          = 0  inside the domain (not boundary)
                             !          = 1  on the wall boundary
                             !          = 2  on the water level boundary
                             !          = 3  on the discharge normal boundary
      real*8,dimension(:),allocatable :: &
              xc,      &     ! (N_CELL)    x-coordinate for the center of each cell
              yc,      &     ! (N_CELL)    y-coordinate for the center of each cell
              sig,     &     ! (NZglobal)  sigma (new z direction coordinate system)
              dsig           ! (NZglobal)  Increment in the sigma-direction
      real*8,dimension(:),allocatable :: &
              xv,      &     ! (N_VERT)    x-coordinate for each vertex
              yv,      &     ! (N_VERT)    y-coordinate for each vertex
              sigv,    &     ! (NZ-1)      sigma at the vertex points
              dsigv          ! (NZ-1)      Increment in the sigma-direction (vertex)

!      ________________________________________________________
!     |                                                        |
!     |                       Others                           |
!     |________________________________________________________|
!     ----------------------------------------
!     (N_CELL) Dispersion coeff. at the cell center
      real*8, dimension(:),allocatable  ::  ah
!     ----------------------------------------
!     Auxiliar variables used to choose saving data
      real*8,dimension(:,:),allocatable :: aux1,aux2,aux3,aux4
      real*8,dimension(:,:),allocatable :: auxv1,auxv2,auxv3,auxv4
!     ----------------------------------------
!     Conexions to display elements in tecplot
      integer,dimension(:), allocatable :: ic1tec,ic2tec,ic3tec
!     ----------------------------------------
!     Adams-Bashforth 2nd order flag
      integer :: flag_ab

      CONTAINS


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  SUBROUTINE: ALLOCATE VARIABLES                     !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

	SUBROUTINE alloc_variables

#	include "common.mpf"

!      ________________________________________________________
!     |                                                        |
!     |               Point location at time                   |
!     |________________________________________________________|

	    allocate(xct(N_CELL,NZ),   &
                 yct(N_CELL,NZ),   &
                 zct(N_CELL,NZ))

        allocate(xvt(N_VERT,NZ-1), &
                 yvt(N_VERT,NZ-1), &
                 zvt(N_VERT,NZ-1))
!      ________________________________________________________
!     |                                                        |
!     |                  FLuid variables                       |
!     |________________________________________________________|

	allocate(alphafnp(N_CELL,NZ),  &
                  alphafn(N_CELL,NZ),  &
                   alphaf(N_CELL,NZ),  &
                     ufnp(N_CELL,NZ),  &
                      ufn(N_CELL,NZ),  &
                       uf(N_CELL,NZ),  &
                      suf(N_CELL,NZ),  &
                     ufBC(N_CELL,NZ),  &
                     vfnp(N_CELL,NZ),  &
                      vfn(N_CELL,NZ),  &
                       vf(N_CELL,NZ),  &
                      svf(N_CELL,NZ),  &
                     vfBC(N_CELL,NZ),  &
                     wfnp(N_CELL,NZ),  &
                      wfn(N_CELL,NZ),  &
                       wf(N_CELL,NZ),  &
                      swf(N_CELL,NZ),  &
                     wfBC(N_CELL,NZ),  &
                     pfnp(N_CELL,NZ),  &
                      pfn(N_CELL,NZ),  &
                       pf(N_CELL,NZ),  &
                      spf(N_CELL,NZ),  &
                    Dpfnp(N_CELL,NZ),  &
                     Dpfn(N_CELL,NZ),  &
                      Dpf(N_CELL,NZ),  &
                  omenpuf(N_CELL,NZ),  &
                  omenpvf(N_CELL,NZ),  &
                  omenpwf(N_CELL,NZ),  &
                     rhof(N_CELL,NZ),  &
                   viscof(N_CELL,NZ))

        allocate(alphafv(N_VERT,NZ-1), &
                     ufv(N_VERT,NZ-1), &
                     vfv(N_VERT,NZ-1), &
                     wfv(N_VERT,NZ-1), &
                     pfv(N_VERT,NZ-1), &
                    Dpfv(N_VERT,NZ-1), &
                   rhofv(N_VERT,NZ-1), &
                 viscofv(N_VERT,NZ-1))
!      ________________________________________________________
!     |                                                        |
!     |                  Solid variables                       |
!     |________________________________________________________|

	allocate(alphasnp(N_CELL,NZ),  &
                  alphasn(N_CELL,NZ),  &
                   alphas(N_CELL,NZ),  &
                     usnp(N_CELL,NZ),  &
                      usn(N_CELL,NZ),  &
                       us(N_CELL,NZ),  &
                     vsnp(N_CELL,NZ),  &
                      vsn(N_CELL,NZ),  &
                       vs(N_CELL,NZ),  &
                     wsnp(N_CELL,NZ),  &
                      wsn(N_CELL,NZ),  &
                       ws(N_CELL,NZ),  &
                     psnp(N_CELL,NZ),  &
                      psn(N_CELL,NZ),  &
                       ps(N_CELL,NZ),  &
                  omenpus(N_CELL,NZ),  &
                  omenpvs(N_CELL,NZ),  &
                  omenpws(N_CELL,NZ),  &
                     rhos(N_CELL,NZ),  &
                   viscos(N_CELL,NZ))

        allocate(alphasv(N_VERT,NZ-1), &
                     usv(N_VERT,NZ-1), &
                     vsv(N_VERT,NZ-1), &
                     wsv(N_VERT,NZ-1), &
                     psv(N_VERT,NZ-1), &
                   rhosv(N_VERT,NZ-1), &
                 viscosv(N_VERT,NZ-1))
!      ________________________________________________________
!     |                                                        |
!     |             Free surface & bottom levels               |
!     |________________________________________________________|

	    allocate( etanp(N_CELL), &
                   etan(N_CELL), &
                    eta(N_CELL), &
                  HprBC(N_CELL), &
                    Hpr(N_CELL), &
                      h(N_CELL), &
                    zbc(N_CELL), &
                 wfsurf(N_CELL), &
                 wssurf(N_CELL))

	    allocate(   etav(N_VERT),  &
                    Hprv(N_VERT),  &
                      hv(N_VERT),  &
                     zbv(N_VERT),  &
                 wfsurfv(N_VERT),  &
                 wssurfv(N_VERT))
!      ________________________________________________________
!     |                                                        |
!     |                 Geometry (unstructure)                 |
!     |________________________________________________________|

	   allocate(No_cp(N_CELL ,3),  &
                     nbe(N_CELL0),  &
                       xc(N_CELL),  &
                       yc(N_CELL),  &
                          sig(NZ),  &
                         dsig(NZ),  &
                 No_vp(N_CELL0,3),  &
                     nbev(N_VERT),  &
                       xv(N_VERT),  &
                       yv(N_VERT),  &
                       sigv(NZ-1),  &
                      dsigv(NZ-1))
!      ________________________________________________________
!     |                                                        |
!     |                       Others                           |
!     |________________________________________________________|

	  allocate(ah(N_CELL))

       allocate(aux1(N_CELL,NZ),aux2(N_CELL,NZ),      &
                aux3(N_CELL,NZ),aux4(N_CELL,NZ),      &
                auxv1(N_VERT,NZ-1),auxv2(N_VERT,NZ-1),&
                auxv3(N_VERT,NZ-1),auxv4(N_VERT,NZ-1))

       allocate(ic1tec(5*N_VERT),&
                ic2tec(5*N_VERT),&
                ic3tec(5*N_VERT))

	END SUBROUTINE alloc_variables


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  SUBROUTINE: DEALLOCATE VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


        SUBROUTINE dealloc_variables

!      ________________________________________________________
!     |                                                        |
!     |               Point location at time                   |
!     |________________________________________________________|

	  deallocate(xct,yct,zct,  &
                   xvt,yvt,zvt)
!      ________________________________________________________
!     |                                                        |
!     |                  FLuid variables                       |
!     |________________________________________________________|

	  deallocate(alphafnp,alphafn,alphaf,    &
                   ufnp,ufn,uf,ufBC,         &
                   vfnp,vfn,vf,vfBC,         &
                   wfnp,wfn,wf,wfBC,         &
                   suf,svf,swf,spf,          &
                   omenpuf,omenpvf,omenpwf,  &
                   pfn,pfnp,pf,              &
                   Dpfn,Dpfnp,Dpf,           &
                   rhof,                     &
                   viscof)

	   deallocate(alphafv,                    &
                   ufv,vfv,wfv,pfv,Dpfv,      &
                   rhofv,viscofv)

!      ________________________________________________________
!     |                                                        |
!     |                  Solid variables                       |
!     |________________________________________________________|


	    deallocate(alphasnp,alphasn,alphas,  &
                   usnp,usn,us,              &
                   vsnp,vsn,vs,              &
                   wsnp,wsn,ws,              &
                   omenpus,omenpvs,omenpws,  &
                   psn,psnp,ps,              &
                   rhos,viscos)

	    deallocate(alphasv,                  &
                   usv,vsv,wsv,psv,          &
                   rhosv,viscosv)
!      ________________________________________________________
!     |                                                        |
!     |             Free surface & bottom levels               |
!     |________________________________________________________|

	    deallocate(etanp,etan,eta,HprBC,     &
                   Hpr,h,zbc,                &
                   wfsurf,wssurf)

	     deallocate(Hprv,hv,zbv,             &
                   etav,                     &
                   wfsurfv,wssurfv)
!      ________________________________________________________
!     |                                                        |
!     |                 Geometry (unstructure)                 |
!     |________________________________________________________|

	    deallocate(No_vp,No_cp,  &
                   xc,yc,        &
                   sig,dsig,     &
                   xv,yv,        &
                   sigv,dsigv,   &
                   nbe,nbev)
!      ________________________________________________________
!     |                                                        |
!     |                       Others                           |
!     |________________________________________________________|

	    deallocate(ah)
        deallocate(aux1,aux2,aux3,aux4,auxv1,auxv2,auxv3,auxv4)
        deallocate(ic1tec,ic2tec,ic3tec)



	END SUBROUTINE dealloc_variables

      END MODULE variables

