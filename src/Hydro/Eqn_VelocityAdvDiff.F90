!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             SOLUTION OF THE ADVECTION-DIFFUSION PROBLEM             !
!                             May 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AdvDiffVelocity_NS(HphiNew,Hphi,Hphiv,               &
                                    phiNew,phi,phiv,                  &
!                                   ----------------------------------
                                    Hpr,eta,                          &
                                    Hprv,etav,                        &
!                                   ----------------------------------
                                    h,hv,                             &
!                                   ----------------------------------
                                    rhs,Gamx,Gamy,Gamz,uu,vv,ww,pfn,  &
!                                   ----------------------------------  
                                    xc,yc,sig,dsig,No_cp,nbe,         &
                                    xv,yv,sigv,dsigv,No_vp,nbev,      & 
!                                   ----------------------------------                                                                       
                                    tagBC)              
                               
!---------------------------------------------------------------------!   
!                                                                     !
!     This subroutine calculates the 3D advection-diffusion equation  !
!     given the velocity profile: (uu,vv,ww),the diffusive components:!
!     (Gamx,Gamy,Gamy) and a right-hand side: rhs as follows:         !
!                                                                     !
!     d(phi)/dt  + ADV(uu,vv,ww,phi) = Diff(Gamx,Gamy,Gamz,phi) + rhs !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |   Name     |   Size     | Description                         |  !  
!  |____________|____________|_____________________________________|  ! 
!  | <--phiNew  |(N_CELL,NZ) | Cell-center solution at t(n+1)      |  !
!  | <--phivNew |(N_VERT,NZ) | Cell-vertex solution at t(n+1)      |  ! 
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | --> uu     |(N_CELL,NZ) | x-velocity component                |  !
!  | --> vv     |(N_CELL,NZ) | y-velocity component                |  !
!  | --> ww     |(N_CELL,NZ) | z-velocity component                |  !
!  | --> Gamx   |(N_CELL,NZ) | Diffusive coefficient in x          |  !
!  | --> Gamy   |(N_CELL,NZ) | Diffusive coefficient in y          |  !
!  | --> Gamz   |(N_CELL,NZ) | Diffusive coefficient in z          |  !
!  | --> rhs    |(N_CELL,NZ) | right-hand side of the problem      |  !
!  |____________|____________|_____________________________________|  !
!  | --> phi    |(N_CELL,NZ) | Cell-center solution at t(n)        |  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ)        | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> phiv   |(N_VERT,NZ) | Cell-vertex solution at t(n)        |  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sigv(k+1)-sigv(k)       |  !  
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbe    |(N_CELL0)   | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!  | --> tagBC  | integer    | Tag = 1:vel.u, =2:vel.v, =3:vel.w   |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  |    Am1     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  |    AmT     |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  |    AmB     |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  !
!  |____________|____________|_____________________________________|  !
!  |    AmG     |(N_CELL0,NZ)| Gradient contribution  (advection)  |  !
!  |____________|____________|_____________________________________|  ! 
!  |    Bmv1T   |(N_CELL0,NZ)| matrix coeff. vertex 1 top    (diff)|  ! 
!  |    Bmv2T   |(N_CELL0,NZ)| matrix coeff. vertex 2 top    (diff)|  ! 
!  |    Bmv3T   |(N_CELL0,NZ)| matrix coeff. vertex 3 top    (diff)|  ! 
!  |    Bmv1B   |(N_CELL0,NZ)| matrix coeff. vertex 1 bottom (diff)|  ! 
!  |    Bmv2B   |(N_CELL0,NZ)| matrix coeff. vertex 2 bottom (diff)|  ! 
!  |    Bmv3B   |(N_CELL0,NZ)| matrix coeff. vertex 3 bottom (diff)|  !
!  |____________|____________|_____________________________________|  !  
!  |    bm      |(N_CELL0,NZ)| right hand side of the method       |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - advection3D                ( advection3D.F90 )            |  !
!  |   - diffusion3D                ( diffusion3D.F90 )            |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - BCvelcenter3D              ( BCvelocity.F90 )             |  !
!  |   - BCvelvertex3D              ( BCvelocity.F90 )             |  !
!  |   - solSOR3D                   ( NewSOR3D.F90 )               |  !
!  |   - solGMRES3D                 ( NewSOR3D.F90 )               |  !
!  |   - gmres3DAdvDiff             ( gmres3DAdvDiff.F90 )         |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!                                                                     !
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

      real*8,dimension(:,:) :: HphiNew(N_CELL,NZ)
      real*8,dimension(:,:) :: Hphi(N_CELL,NZ)      
      real*8,dimension(:,:) :: Hphiv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: phiNew(N_CELL,NZ)
      real*8,dimension(:,:) :: phi(N_CELL,NZ)      
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     ----------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     ----------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
      real*8,dimension(:,:) :: pfn(N_CELL,NZ)
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
!     ----------------------------------------
      integer:: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:) :: Flow(N_CELL,NZ)
      real*8,dimension(:,:) :: Flowv(N_VERT,NZ-1)      
      real*8,dimension(:,:) :: Am0C(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1C(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2C(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3C(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmTC(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmBC(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmG(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)      
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
      real*8,dimension(:,:) :: bm(N_CELL0,NZ)  
!     --------------------------------------
      real*8 :: Vol,som
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'Begin subroutine: AdvDiffVelocity_NSIrre'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                       Advection-Diffusion  3D                       !
!                                                                     !
!*********************************************************************!

      Flow  = phi
      Flowv = phiv
                      
      HphiNew = Hphi
      phiNew  = phi
!      ________________________________________________________
!     |                                                        |
!     |                 Diffusion contribution                 |
!     |________________________________________________________|
	
      Am0 = 0.0d0
      Am1 = 0.0d0
      Am2 = 0.0d0
      Am3 = 0.0d0
      AmT = 0.0d0
      AmB = 0.0d0
!     -------------
      Bmv1T = 0.0d0 
      Bmv2T = 0.0d0 
      Bmv3T = 0.0d0 
      Bmv1B = 0.0d0 
      Bmv2B = 0.0d0 
      Bmv3B = 0.0d0 
!     -------------
#     ifdef KeyDiffusion            
      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       -Gamx,-Gamy,-Gamz,                   &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)
#     endif  
                                             
!     ----------                  
!     * REMARK: The negative sign in the diffusion coefficients
!               is because we move it to the left side of the
!               transport model.
                                                    
!      ________________________________________________________
!     |                                                        |
!     |          Advection (Convection) contribution           |
!     |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0	
            Am0C(i,k)  = 0.0d0
            Am1C(i,k)  = 0.0d0
            Am2C(i,k)  = 0.0d0
            Am3C(i,k)  = 0.0d0
            AmTC(i,k)  = 0.0d0
            AmBC(i,k)  = 0.0d0
!           -----------------
            AmG(i,k)  = 0.0d0             
         enddo
      enddo

      call advectionVelocity(Am0C,Am1C,Am2C,Am3C,AmTC,AmBC,AmG, &
                             uu,vv,ww,pfn,                      &
                             Flow,xc,yc,sig,No_cp,nbe,          &  
                             Flowv,xv,yv,sigv,dsigv,No_vp,nbev)
!      ________________________________________________________
!     |                                                        |
!     |                        Solution                        |
!     |________________________________________________________|

!     _________________________________________________________
!     FullExplicit             

#     ifdef KeyAdvFullExplicit
         do k=2,NZ-1
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               Vol = areaCell(i)*dsigv(k-1)
!              ----------------
!              New rhs
               som =  Am0C(i,k)*Flow(i,k)     &
                    + Am1C(i,k)*Flow(jc1,k)   &
                    + Am2C(i,k)*Flow(jc2,k)   &
                    + Am3C(i,k)*Flow(jc3,k)   &
                    + AmTC(i,k)*Flow(i,k+1)   &       
                    + AmBC(i,k)*Flow(i,k-1)   &
                    + AmG(i,k)                &
!                   -------------------------- 
                    + Am0(i,k)*phi(i,k)       &
                    + Am1(i,k)*phi(jc1,k)     &
                    + Am2(i,k)*phi(jc2,k)     &
                    + Am3(i,k)*phi(jc3,k)     &
                    + AmT(i,k)*phi(i,k+1)     &       
                    + AmB(i,k)*phi(i,k-1)     &
                    + Bmv1T(i,k)*phiv(jv1,k)  &
                    + Bmv2T(i,k)*phiv(jv2,k)  &
                    + Bmv3T(i,k)*phiv(jv3,k)  &
                    + Bmv1B(i,k)*phiv(jv1,k-1)&
                    + Bmv2B(i,k)*phiv(jv2,k-1)&
                    + Bmv3B(i,k)*phiv(jv3,k-1)  
!               ----------------
!               Update
                HphiNew(i,k) = Hphi(i,k) + dt*(-som/Vol + rhs(i,k))
                phiNew(i,k)  = HphiNew(i,k)/Hpr(i)
            enddo
         enddo
#     endif

!     _________________________________________________________
!     Explicit  

#     ifdef KeyAdvExplicit
         do k=2,NZ-1
            do i=1,N_CELL0 	
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               Vol = areaCell(i)*dsigv(k-1)
!              ----------------
!              New diagonal values
               Am0(i,k)=  dt/Vol*(Am0C(i,k)+Am0(i,k)) + Hpr(i) 
!              ----------------
!              New rhs
               som =  Am1C(i,k)*Flow(jc1,k)   &
                    + Am2C(i,k)*Flow(jc2,k)   &
                    + Am3C(i,k)*Flow(jc3,k)   &
                    + AmTC(i,k)*Flow(i,k+1)   &       
                    + AmBC(i,k)*Flow(i,k-1)   &
                    + AmG(i,k)                &
!                   -------------------------- 
                    + Am1(i,k)*phi(jc1,k)     &
                    + Am2(i,k)*phi(jc2,k)     &
                    + Am3(i,k)*phi(jc3,k)     &
                    + AmT(i,k)*phi(i,k+1)     &       
                    + AmB(i,k)*phi(i,k-1)     &
                    + Bmv1T(i,k)*phiv(jv1,k)  &
                    + Bmv2T(i,k)*phiv(jv2,k)  &
                    + Bmv3T(i,k)*phiv(jv3,k)  &
                    + Bmv1B(i,k)*phiv(jv1,k-1)&
                    + Bmv2B(i,k)*phiv(jv2,k-1)&
                    + Bmv3B(i,k)*phiv(jv3,k-1)    
               bm(i,k) = -som/Vol + rhs(i,k)
!              ----------------
!              Update
               phiNew(i,k)  = (Hphi(i,k)+dt*bm(i,k))/Am0(i,k)
               HphiNew(i,k) = phiNew(i,k)*Hpr(i)
            enddo
         enddo
#     endif

!     _________________________________________________________
!     Semi-Implicit  
#     ifdef KeyAdvSemiImplicit
!        ______________________
!        New Matrix & rhs
         do k=1,NZ
            do i=1,N_CELL0
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               Vol = areaCell(i)*dsigv(k-1)
!              ----------------
!              New matrix coeff.
               Am0(i,k)=  dt/Vol*(Am0C(i,k)+Am0(i,k)) + Hpr(i) 
               Am1(i,k) = dt/Vol*(Am1C(i,k)+Am1(i,k))
               Am2(i,k) = dt/Vol*(Am2C(i,k)+Am2(i,k))
               Am3(i,k) = dt/Vol*(Am3C(i,k)+Am3(i,k))
               AmT(i,k) = dt/Vol*(AmTC(i,k)+AmT(i,k))
               AmB(i,k) = dt/Vol*(AmBC(i,k)+AmB(i,k))      
!              ----------------
!              New rhs
               som =   AmG(i,k)                &
                     + Bmv1T(i,k)*phiv(jv1,k)  &
                     + Bmv2T(i,k)*phiv(jv2,k)  &
                     + Bmv3T(i,k)*phiv(jv3,k)  &
                     + Bmv1B(i,k)*phiv(jv1,k-1)&
                     + Bmv2B(i,k)*phiv(jv2,k-1)&
                     + Bmv3B(i,k)*phiv(jv3,k-1) 
               bm(i,k) = Hphi(i,k) + dt*(-som/Vol + rhs(i,k))
            enddo
         enddo
!        ______________________
!        Call the linear solver
         HphiNew = Hphi
         phiNew  = phi
         call SOR3Dvelocity_Irre(HphiNew,phiNew,               &
                                 Hphiv,phiv,                   &
                                 Hpr,eta,                      &
                                 Hprv,etav,                    &
                                 h,hv,                         &
                                 Am0,Am1,Am2,Am3,AmT,AmB,bm,   & 
                                 xc,yc,sig,dsig,No_cp,nbe,     &
                                 xv,yv,sigv,dsigv,No_vp,nbev,  &
                                 tagBC)       
#     endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(HphiNew)
         call communication3Dtype2(phiNew)
#     endif
!     =============== END ================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                     BC (old code)                      |
!     |________________________________________________________|

!     ----------------------------------------
!     Boundary Conditions        
      !call BCvelcenter3D(HphiNew,phiNew,                  &
      !                   xc,yc,sig,dsig,No_cp,nbe,        &
      !                   Hpr,h,eta,                       &
      !                   tagBC)
                                                      
!     ----------------------------------------
!     Boundary Conditions 
      !call BCvelvertex3D(Hphiv,phiv,                  &
      !                   xv,yv,sigv,dsigv,No_vp,nbev, &
      !                   Hprv,hv,etav,                &
      !                   tagBC)

!      ________________________________________________________
!     |                                                        |
!     |         Boundary condition & Vertex values             |
!     |________________________________________________________|


!     ----------------------------------------
!     BC: Cell-center
      call BCvelocity3D(HphiNew,phiNew,              &
                        Hphiv,phiv,                  &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        tagBC) ! positive for cell-center

!     ----------------------------------------
!     Vertex by interpolation
      call interpolation3D(Hphiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           HphiNew,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           phiNew,xc,yc,sig,dsig,No_cp,nbe)

!     ----------------------------------------
!     BC: Vertex
      call BCvelocity3D(HphiNew,phiNew,              &
                        Hphiv,phiv,                  &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        -tagBC) ! negative for vertex

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(HphiNew)
         call communication3Dtype2(phiNew)
#     endif
!     =============== END ================
!     ====================================

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), ' End   subroutine AdvDiffvelocity_NSIrre'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF Advection-Diffusion 3D                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             S.O.R. FOR CELL-CENTER VARIABLES & VELOCITY BC          !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SOR3Dvelocity_Irre(Hphi,phi,                          &
                                    Hphiv,phiv,                        &
                                    Hpr,eta,                           &
                                    Hprv,etav,                         &
                                    h,hv,                              &
                                    MAm0,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm, &
                                    xc,yc,sig,dsig,No_cp,nbe,          &
                                    xv,yv,sigv,dsigv,No_vp,nbev,       &
                                    tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the S.O.R. techni-  !
!    que for different relaxion factors: relax. This system is good   !
!    for the linear systems with coefficients "MAm" and righ-hand     !
!    side "Vbm" that only depends on cell-center phi values.          !
!    This program is almost identical as the solSOR3D but now we      !
!    can decide between the velocity BC that we want to impose.       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi   |(N_CELL,NZ)  | Solution & initial guess           |  !
!  |____________|_____________|____________________________________|  !
!  | ---> MAm0  |(N_CELL0,NZ) | matrix coefficient of element i    |  !
!  | ---> MAm1  |(N_CELL0,NZ) | matrix coeff. horizontal neigborn 1|  ! 
!  | ---> MAm2  |(N_CELL0,NZ) | matrix coeff. horizontal neigborn 2|  ! 
!  | ---> MAm3  |(N_CELL0,NZ) | matrix coeff. horizontal neigborn 3|  ! 
!  | ---> MAmT  |(N_CELL0,NZ) | matrix coeff. vertical top         |  ! 
!  | ---> MAmB  |(N_CELL0,NZ) | matrix coeff. vertical bottom      |  !
!  | ---> Vbm   |(N_CELL0,NZ) | right hand side of the method      |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> sig   | (NZ)        | sigma value at the cell centers    |  !
!  | ---> dsig  | (NZ)        | = sig(k+1)-sig(k+1)                |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> tagBC | integer     | Tag related to the veclocity BC    |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - BCvelcenter3D                ( BCvelocity.F90 )           |  !
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

      real*8,dimension(:,:) :: Hphi(N_CELL,NZ)
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: Hphiv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:)   :: Hpr(N_CELL)
      real*8,dimension(:)   :: eta(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     -------------------------------------
      real*8,dimension(:)   :: h(N_CELL) 
      real*8, dimension(:)  :: hv(N_VERT)   
!     -------------------------------------
      real*8,dimension(:,:) :: MAm0(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Vbm(N_CELL0,NZ)
!     -------------------------------------
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
!     -------------------------------------
      integer :: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer:: jc1,jc2,jc3
!     --------------------------------------
      real*8  :: errorsys,residu,som,SUMerrorsys
      integer :: it,Display

!      ____________________________________
!     |                                    |
!     |             Parameter              |
!     |____________________________________|

      real*8, parameter :: JrelaxSOR = 1.0

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: SOR3Dvelocity_Irre'
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
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

      it=0
111   continue
      it=it+1 

!     ________________________________________________________
!     Solution of the system SOR (JSOR in parallel)

!     --------------------
!     Update: phi
      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            som =  MAm1(i,k)*phi(jc1,k) &
                 + MAm2(i,k)*phi(jc2,k) &
                 + MAm3(i,k)*phi(jc3,k) &
                 + MAmT(i,k)*phi(i,k+1) &       
                 + MAmB(i,k)*phi(i,k-1) 
            residu = (Vbm(i,k)-som)/MAm0(i,k)-phi(i,k)
            errorsys = errorsys + abs(residu)
            phi(i,k) = phi(i,k) + JrelaxSOR*residu
            Hphi(i,k)= phi(i,k)*Hpr(i)
         enddo
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(Hphi)
         call communication3Dtype2(phi)
#     endif
!     =============== END ================
!     ====================================

      call BCvelocity3D(Hphi,phi,                    &
                        Hphiv,phiv,                  &
                        Hpr,eta,                     &
                        Hprv,etav,                   &
                        h,hv,                        &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        tagBC)

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(Hphi)
         call communication3Dtype2(phi)
         !--------------------------	 
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
#     endif
!     =============== END ================
!     ====================================

!     ________________________________________________________
!     Convergence criteria
     
      7 format(t10,a24,i5,a9,e10.3)
      if (errorsys.lt.eps) then
         IF (Display.eq.1) THEN
#        ifdef KeyParallel
         write(*,7) 'Velocity JS0R 3D: iters =',it,', error =',errorsys
#        else
         write(*,7) 'Velocity  S0R 3D: iters =',it,', error =',errorsys
#        endif
         ENDIF
      elseif (errorsys.gt.1.5d6) then
         IF (Display.eq.1) THEN
         write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         ENDIF
         stop
      elseif(it.gt.MaxIters) then
         IF (Display.eq.1) THEN
         write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
         ENDIF
      else
         goto 111
      endif

119   continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine:SOR3Dvelocity_Irre'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      End of S.O.R. Methods 3D                       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
