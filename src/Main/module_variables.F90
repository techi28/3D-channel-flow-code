!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     DEFINITION OF THE VARIABLES                     !
!                             Dic 2017                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

	MODULE variables 

!---------------------------------------------------------------------!   
!                                                                     !
!      SUBROUTINES:     -  alloc_variables                            !
!                       -  dealloc_varaibles                          !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!  ____________________________________                               !
!  Index                                                              !
!  nbe,nbev   :   (N_VERT) Assigned tag of the type of vertex.        !
!  No_wb      :   (N_WB)   Vertex numbering of the wall bdy.          !
!  No_qb      :   (N_QB)   Vertex numbering of the dischar. normal bdy!
!  No_hb      :   (N_HB)   Vertex numbering of water surface bdy.     !
!  No_sp      :   (2)      Vertex numbering of the sample             !
!  ____________________________________                               !
!  Geometry                                                           !
!  xc         :   (N_CELL)    x-coordinate for the cell-center        !    
!  yc         :   (N_CELL)    y-coordinate for the cell-center        !
!  sig        :   (NZ)        sigma-coordinate (new z direction)      !
!  dsig       :   (NZ)        Increment in the sigma-direction        !
!  xv         :   (N_VERT)    x-coordinate for each vertex            !
!  yv         :   (N_VERT)    y-coordinate for each vertex            !
!  sigv       :   (NZ-1)      sigma-coordinate at the vertex points   !
!  dsigv      :   (NZ-1)      Increment in the sigma-direction(vertex)!
!  xct        :   (N_CELL,NZ) cell-center x-coordinate at time t      !
!  yct        :   (N_CELL,NZ) cell-center y-coordinate at time t      !    
!  zct        :   (N_CELL,NZ) cell-center z-coordinate at time t      !
!  xvt        :   (N_VERT,NZ) vertex x-coordinate at time t           !
!  yvt        :   (N_VERT,NZ) vertex y-coordinate at time t           !    
!  zvt        :   (N_VERT,NZ) vertex z-coordinate at time t           !
!  ____________________________________                               !
!  Free surface & bottom                                              ! 
!  Hprnp       :   (N_CELL) Total water depth: H=h+eta (n+1)          !
!  Hprn        :   (N_CELL) Total water depth: H=h+eta (n)            !
!  Hpr         :   (N_CELL) Total water depth: H=h+eta (RK)           !
!  Hprv,       :   (N_VERT) Total water depth at vertex               !
!  etav,       :   (N_VERT) Free surface level eta at the vertices    !
!  etanp       :   (N_CELL) Free surface level eta at (n+1)           !
!  etan        :   (N_CELL) Free surface level eta at (n)             !
!  eta         :   (N_CELL) Free surface level eta at (RK)            !
!  h           :   (N_CELL) The bottom depth                          !
!  hv,         :   (N_VERT) Depth of the river at vertex              !
!  zbc         :   (N_CELL) Bottom levels at each cell                !
!  zbv,        :   (N_VERT) Bottom levels at each vertex              !
!  ____________________________________                               !
!  Fluid                                                              !
!  alphafnp    :   Fluid control volume at(n+1)                       !
!  alphafn     :   Fluid control volume at (n)                        !
!  alphaf      :   Fluid control volume at RK step                    !
!  alphafv     :   Fluid control volume at vertex                     !
!  rhof        :   Density of the fluid                               !
!  viscof      :   Viscosity of the fluid                             !
!  rhofv       :   Density of the fluid vertex                        !
!  viscofv     :   Viscosity of the fluid vertex                      ! 
!  ufnp        :   Velocity component u_f at (n+1)                    !
!  ufn         :   Velocity component u_f at (n)                      !
!  uf          :   Velocity component u_f at RK step                  !
!  ufv         :   Velocity component u_f at vertex                   !
!  vfnp        :   Velocity component v_f at (n+1)                    !           
!  vfn         :   Velocity component v_f at (n)                      !
!  vf          :   Velocity component v_f at RK step                  !
!  vfv         :   Velocity component v_f at vertex                   !       
!  wfnp        :   Velocity component w_f at (n+1)                    !
!  wfn         :   Velocity component w_f at (n)                      !
!  wf          :   Velocity component w_f at RK step                  !
!  wfv         :   Velocity component w_f at vertex                   !
!  pfnp        :   Pressure of the fluid at (n+1)                     !
!  pfn         :   Pressure of the fluid at (n)                       !
!  pf          :   Pressure of the fluid at RK step                   !
!  pfv         :   Pressure of the fluid at vertex                    !
!  ____________________________________                               !
!  Solid                                                              !
!  alphasnp    :   Solid control volume at(n+1)                       !
!  alphasn     :   Solid control volume at (n)                        !
!  alphas      :   Solid control volume at RK step                    !
!  alphasv     :   Solid control volume at vertex                     !
!  rhos        :   Density of the solid                               !
!  viscos      :   Viscosity of the solid                             !
!  rhosv       :   Density of the solid vertex                        !
!  viscosv     :   Viscosity of the solid vertex                      !
!  usnp        :   Velocity component u_s at (n+1)                    !
!  usn         :   Velocity component u_s at (n)                      !
!  us          :   Velocity component u_s at RK step                  !
!  usv         :   Velocity component u_s at vertex                   !
!  vsnp        :   Velocity component v_s at (n+1)                    !         
!  vsn         :   Velocity component v_s at (n)                      !
!  vs          :   Velocity component v_s at RK step                  !
!  vsv         :   Velocity component v_s at vertex                   !
!  wsnp        :   Velocity component w_s at (n+1)                    !
!  wsn         :   Velocity component w_s at (n)                      !
!  ws          :   Velocity component w_s at RK step                  !
!  wsv         :   Velocity component w_s at vertex                   !
!  psnp        :   Pressure of the solid at (n+1)                     !
!  psn         :   Pressure of the solid at (n)                       !
!  ps          :   Pressure of the solid at RK step                   ! 
!  psv         :   Pressure of the solid at vertex                    !
!  ____________________________________                               !
!  Fluid methodology                                                  !
!  Dpfnp       :   Difference of pressure of the fluid at (n+1)       !
!  Dpfn        :   Difference of pressure of the fluid at (n)         ! 
!  Dpf         :   Difference of pressure of the fluid at RK step     !
!  Dpfv        :   Difference of pressure of the fluid at vertex      !
!  omenpuf     :   omega at (n+1) for u_f                             !
!  omenpvf     :   omega at (n+1) for v_f                             !
!  omenpwf     :   omega at (n+1) for w_f                             !
!  ____________________________________                               !
!  Analytical & Errors                                                !
!  ____________________________________                               !
!  Boundary                                                           !
!  qb,         :   (N_BC)   Value for discharges at the bdy. vertex   !
!  wlb         :   (N_BC)   Value for water levels at the bdy. vertex !
!  wfsurf      :   (N_CELL) Velocity component w_f at the free surface!
!  wfsurfv     :   (N_VERT) w_f at the free surface vertex            !
!  ____________________________________                               !
!  Others                                                             !
!  Heaviside   : (N_CELL)    Step function: =1 inside & =0 outside    !
!  mask        : (N_CELL,NZ) Profile description of the input release !
!                                                                     !
!---------------------------------------------------------------------!

      implicit none

!     ____________________________________
!     Index
      integer,dimension(:,:),allocatable :: No_cp
      integer,dimension(:,:),allocatable :: No_vp
      integer,dimension(:),  allocatable :: No_wb,No_qb,No_hb,No_sp
      integer,dimension(:),  allocatable :: nbe,nbev
!     ____________________________________
!     Geometry      
      real*8, dimension(:),  allocatable :: xc,yc,sig,dsig 
      real*8, dimension(:),  allocatable :: xv,yv,sigv,dsigv
      real*8, dimension(:,:),allocatable :: xct,yct,zct
      real*8, dimension(:,:),allocatable :: xvt,yvt,zvt      
!     ____________________________________
!     Free surface & bottom 
      real*8,dimension(:),allocatable :: Hprnp,Hprn,Hpr,Hprv
      real*8,dimension(:),allocatable :: etanp,etan,eta,etav
      real*8,dimension(:),allocatable :: h,hv
      real*8,dimension(:),allocatable :: zbc,zbv
!     ____________________________________
!     Fluid 
      real*8,dimension(:,:),allocatable :: alphafnp,alphafn,alphaf
      real*8,dimension(:,:),allocatable :: alphafv
!     ------
      real*8,dimension(:,:),allocatable :: rhof,viscof
      real*8,dimension(:,:),allocatable :: rhofv,viscofv
!     ------
      real*8,dimension(:,:),allocatable :: ufnp,ufn,uf,ufv
      real*8,dimension(:,:),allocatable :: vfnp,vfn,vf,vfv        
      real*8,dimension(:,:),allocatable :: wfnp,wfn,wf,wfv
      real*8,dimension(:,:),allocatable :: pfnp,pfn,pf,pfv
!     ____________________________________
!     Solid 
      real*8,dimension(:,:),allocatable :: alphasnp,alphasn,alphas
      real*8,dimension(:,:),allocatable :: alphasv
!     ------
      real*8,dimension(:,:),allocatable :: rhos,viscos
      real*8,dimension(:,:),allocatable :: rhosv,viscosv
!     ------           
      real*8,dimension(:,:),allocatable :: usnp,usn,us,usv
      real*8,dimension(:,:),allocatable :: vsnp,vsn,vs,vsv               
      real*8,dimension(:,:),allocatable :: wsnp,wsn,ws,wsv
      real*8,dimension(:,:),allocatable :: psnp,psn,ps,psv
!     ____________________________________
!     Fluid methodology 
      real*8,dimension(:,:),allocatable :: omenpuf,omenpvf,omenpwf
      real*8,dimension(:,:),allocatable :: Dpfnp,Dpfn,Dpf,Dpfv
!     ____________________________________
!     Analytical & Errors 
      real*8,dimension(:),  allocatable :: HprvA,etavA      
      real*8,dimension(:,:),allocatable :: ufvA,vfvA,wfvA,pfvA
      real*8,dimension(:,:),allocatable :: uErr,vErr,wErr,pErr      
      real*8,dimension(:,:),allocatable :: uErrv,vErrv,wErrv,pErrv
!     ____________________________________
!     Boundary
      real*8,dimension(:),allocatable :: qb
      real*8,dimension(:),allocatable :: wlb
      real*8,dimension(:),allocatable :: wfsurf,wfsurfv
!     ____________________________________
!     Others
      real*8, dimension(:,:),allocatable :: mask
      real*8, dimension(:),  allocatable :: Heaviside
      integer,dimension(:),  allocatable :: ic1tec,ic2tec,ic3tec
      
      CONTAINS


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  SUBROUTINE: ALLOCATE VARIABLES                     !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

	SUBROUTINE alloc_variables

#	include "common.mpf"

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: alloc_variables'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ____________________________________
!     Index
      allocate(No_cp(N_CELL,3),No_vp(N_CELL0,3))
      allocate(nbev(N_VERT),nbe(N_CELL0)) 
      allocate(No_wb(N_WBMAX),No_qb(N_QBMAX),No_hb(N_HBMAX), &
               No_sp(N_SPmax))             

!     ____________________________________
!     Geometry      
      allocate(xc(N_CELL),yc(N_CELL),sig(NZ),dsig(NZ),       &
               xv(N_VERT),yv(N_VERT),sigv(NZ-1),dsigv(NZ-1))     
      allocate(xct(N_CELL,NZ),yct(N_CELL,NZ),zct(N_CELL,NZ), &
               xvt(N_VERT,NZ-1),yvt(N_VERT,NZ-1),zvt(N_VERT,NZ-1))  
                  
!     ____________________________________
!     Free surface & bottom 
      allocate(etanp(N_CELL),etan(N_CELL),eta(N_CELL),etav(N_VERT), &
               Hprnp(N_CELL),Hprn(N_CELL),Hpr(N_CELL),Hprv(N_VERT))
      allocate(h(N_CELL),hv(N_VERT))
      allocate(zbc(N_CELL),zbv(N_VERT))
      
!     ____________________________________
!     Fluid 
      allocate(alphafnp(N_CELL,NZ),alphafn(N_CELL,NZ), &
               alphaf(N_CELL,NZ),alphafv(N_VERT,NZ-1))
      allocate(rhof(N_CELL,NZ),viscof(N_CELL,NZ),&
               rhofv(N_VERT,NZ-1),viscofv(N_VERT,NZ-1))  
      allocate(ufnp(N_CELL,NZ),ufn(N_CELL,NZ),uf(N_CELL,NZ), &
               vfnp(N_CELL,NZ),vfn(N_CELL,NZ),vf(N_CELL,NZ), &
               wfnp(N_CELL,NZ),wfn(N_CELL,NZ),wf(N_CELL,NZ), &
               pfnp(N_CELL,NZ),pfn(N_CELL,NZ),pf(N_CELL,NZ), &
               ufv(N_VERT,NZ-1),vfv(N_VERT,NZ-1), &
               wfv(N_VERT,NZ-1),pfv(N_VERT,NZ-1))   
!     ____________________________________
!     Solid 
      allocate(alphasnp(N_CELL,NZ),alphasn(N_CELL,NZ), &
               alphas(N_CELL,NZ),alphasv(N_VERT,NZ-1))
      allocate(rhos(N_CELL,NZ),viscos(N_CELL,NZ), &
               rhosv(N_VERT,NZ-1),viscosv(N_VERT,NZ-1))  
      allocate(usnp(N_CELL,NZ),usn(N_CELL,NZ),us(N_CELL,NZ),  &
               vsnp(N_CELL,NZ),vsn(N_CELL,NZ),vs(N_CELL,NZ),  &
               wsnp(N_CELL,NZ),wsn(N_CELL,NZ),ws(N_CELL,NZ),  &
               psnp(N_CELL,NZ),psn(N_CELL,NZ),ps(N_CELL,NZ),  &
               usv(N_VERT,NZ-1),vsv(N_VERT,NZ-1),             &
               wsv(N_VERT,NZ-1),psv(N_VERT,NZ-1))
!     ____________________________________
!     Fluid methodology 
      allocate(omenpuf(N_CELL,NZ),omenpvf(N_CELL,NZ),omenpwf(N_CELL,NZ))
      allocate(Dpfnp(N_CELL,NZ),Dpfn(N_CELL,NZ),Dpf(N_CELL,NZ)) 
      allocate(Dpfv(N_VERT,NZ-1)) 
       
!     ____________________________________
!     Analytical & Errors 
      allocate(HprvA(N_VERT),etavA(N_VERT))
      allocate(ufvA(N_VERT,NZ-1),vfvA(N_VERT,NZ-1),   &       
               wfvA(N_VERT,NZ-1),pfvA(N_VERT,NZ-1))
      allocate(uErr(N_CELL,NZ),vErr(N_CELL,NZ),       &
               wErr(N_CELL,NZ),pErr(N_CELL,NZ),       &
               uErrv(N_VERT,NZ-1),vErrv(N_VERT,NZ-1), &
               wErrv(N_VERT,NZ-1),pErrv(N_VERT,NZ-1))
               
!     ____________________________________
!     Boundary
      allocate(qb(N_QB),wlb(N_HB))
      allocate(wfsurf(N_CELL),wfsurfv(N_VERT))
      
!     ____________________________________
!     Others
      allocate(mask(N_CELL,NZ),Heaviside(N_CELL)) 
      allocate(ic1tec(5*N_VERT),ic2tec(5*N_VERT),ic3tec(5*N_VERT))
        

      END SUBROUTINE alloc_variables


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  SUBROUTINE: DEALLOCATE VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE dealloc_variables

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: dealloc_variables'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ____________________________________
!     Index
      deallocate(No_cp)
      deallocate(No_vp)
      deallocate(No_wb,No_qb,No_hb,No_sp)
!      deallocate(nbe,nbev)
      print*,'Not deallocate(nbe,nbev)' 
!     ____________________________________
!     Geometry      
      deallocate(xc,yc,sig,dsig) 
      deallocate(xv,yv,sigv,dsigv)
      deallocate(xct,yct,zct)
      deallocate(xvt,yvt,zvt)    
!     ____________________________________
!     Free surface & bottom 
      deallocate(Hprnp,Hprn,Hpr,Hprv)
      deallocate(etanp,etan,eta,etav)
      deallocate(h,hv)
      deallocate(zbc,zbv) 
!     ____________________________________
!     Fluid 
      deallocate(alphafnp,alphafn,alphaf)
      deallocate(alphafv)
!     ------
      deallocate(rhof,viscof)
      deallocate(rhofv,viscofv)
!     ------
      deallocate(ufnp,ufn,uf,ufv)
      deallocate(vfnp,vfn,vf,vfv)
      deallocate(wfnp,wfn,wf,wfv)
      deallocate(pfnp,pfn,pf,pfv) 
!     ____________________________________
!     Solid 
      deallocate(alphasnp,alphasn,alphas)
      deallocate(alphasv)
!     ------
      deallocate(rhos,viscos)
      deallocate(rhosv,viscosv)
!     ------           
      deallocate(usnp,usn,us,usv)
      deallocate(vsnp,vsn,vs,vsv)             
      deallocate(wsnp,wsn,ws,wsv)
      deallocate(psnp,psn,ps,psv) 
!     ____________________________________
!     Fluid methodology 
      deallocate(omenpuf,omenpvf,omenpwf)
      deallocate(Dpfnp,Dpfn,Dpf,Dpfv) 
!     ____________________________________
!     Analytical & Errors 
      deallocate(HprvA,etavA)
      deallocate(ufvA,vfvA,wfvA,pfvA)
      deallocate(uErr,vErr,wErr,pErr)     
      deallocate(uErrv,vErrv,wErrv,pErrv) 
!     ____________________________________
!     Boundary
      deallocate(qb)
      deallocate(wlb)
      deallocate(wfsurf,wfsurfv) 
!     ____________________________________
!     Others
      deallocate(mask)
      deallocate(Heaviside)
      deallocate(ic1tec,ic2tec,ic3tec)
      

      END SUBROUTINE dealloc_variables

      END MODULE variables

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           END OF VARIABLES                          !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
