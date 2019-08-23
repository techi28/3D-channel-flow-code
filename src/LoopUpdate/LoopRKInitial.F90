!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   Initialization of the RK method                   !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE LoopRKInitial(eta,                         &
                               alphaf,uf,vf,wf,pf,          &
                               alphas,us,vs,ws,ps,          &  
                               etan,                        &
                               alphafn,ufn,vfn,wfn,pfn,     &
                               alphasn,usn,vsn,wsn,psn,     &
                               xc,yc,sig,Hpr,h,             &
                               xv,yv,sigv,Hprv,hv,          &
                               No_cp,nbe)     

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the values of the main variables for        !
!    the initial RK iteration .                                       ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !  
!  |______________|___________|____________________________________|  !
!  | <-- etan     |(N_CELL)   | Free surface level eta      (known)|  !
!  |______________|___________|____________________________________|  !
!  | <-- alphaf   |(N_CELL,NZ)| Fluid control volume        (known)|  !
!  | <-- uf       |(N_CELL,NZ)| Velocity component u_f      (known)|  !
!  | <-- vf       |(N_CELL,NZ)| Velocity component v_f      (known)|  !      
!  | <-- wf       |(N_CELL,NZ)| Velocity component w_f      (known)|  !
!  | <-- pf       |(N_CELL,NZ)| Pressure of the fluid       (known)|  !
!  | <-- alphas   |(N_CELL,NZ)| Solid control volume        (known)|  !
!  | <-- us       |(N_CELL,NZ)| Velocity component u_s      (known)|  !
!  | <-- vs       |(N_CELL,NZ)| Velocity component v_s      (known)|  !     
!  | <-- ws       |(N_CELL,NZ)| Velocity component w_s      (known)|  !
!  | <-- ps       |(N_CELL,NZ)| Pressure of the solid       (known)|  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !      
!  |______________|___________|____________________________________|  !
!  | --> etanp    |(N_CELL)   | Free surface level eta       t(n+1)|  !
!  |______________|___________|____________________________________|  !
!  | --> alphafnp |(N_CELL,NZ)| Fluid control volume           t(n)|  !
!  | --> ufnp     |(N_CELL,NZ)| Velocity component u_f         t(n)|  !
!  | --> vfnp     |(N_CELL,NZ)| Velocity component v_f         t(n)|  !      
!  | --> wfnp     |(N_CELL,NZ)| Velocity component w_f         t(n)|  !
!  | --> pfnp     |(N_CELL,NZ)| Pressure of the fluid          t(n)|  !
!  | --> alphasnp |(N_CELL,NZ)| Solid control volume           t(n)|  !
!  | --> usnp     |(N_CELL,NZ)| Velocity component u_s         t(n)|  !
!  | --> vsnp     |(N_CELL,NZ)| Velocity component v_s         t(n)|  !     
!  | --> wsnp     |(N_CELL,NZ)| Velocity component w_s         t(n)|  !
!  | --> psnp     |(N_CELL,NZ)| Pressure of the solid          t(n)|  !
!  |______________|___________|____________________________________|  !
!  | --> xc       |(N_CELL)   | x-coordinate of the cell center    |  !
!  | --> yc       |(N_CELL)   | y-coordinate of the cell center    |  !
!  | --> sig      |(NZ)       | sigma                              |  !
!  | --> Hpr      |(N_CELL)   | Total water depth: H = h + eta     |  !
!  | --> h        |(N_CELL)   | Depth of the domain at each cell   |  !
!  | --> xv       |(N_VERT)   | x-coordinate of the cell vertex    |  !
!  | --> yv       |(N_VERT)   | y-coordinate of the cell vertex    |  !   
!  | --> sigv     |(NZ-1)     | sigma of the vertex points         |  !
!  | --> Hprv     |(N_VERT)   | Total water depth at the vertex    |  !
!  | --> hv       |(N_VERT)   | Depth of the domain at each vertex |  !
!  |______________|___________|____________________________________|  !
!  | --> No_cp    |(N_CELL,3) | Node No. of surrounding three cell |  !
!  | --> nbe      |(N_CELL)   | Type of boundary cell: 0,1,2,3     |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  ! 
!  | --- N_CELL  | Total number of the cells                       |  !
!  | --- N_CELL0 | Number of inside cells                          |  !
!  |     NZ      | Points in the sigma direction                   |  !
!  |_____________|_________________________________________________|  !
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

!      ________________________________________________________
!     |                                                        |
!     |            Keys, subroutines and parameters            |
!     |________________________________________________________|

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

!      ________________________________________________________
!     |                                                        |
!     |               Definition of variables                  |
!     |________________________________________________________|

      real*8,dimension(:)   :: eta(N_CELL)  
      real*8,dimension(:,:) :: alphaf(N_CELL,NZ)
      real*8,dimension(:,:) :: uf(N_CELL,NZ)
      real*8,dimension(:,:) :: vf(N_CELL,NZ)
      real*8,dimension(:,:) :: wf(N_CELL,NZ)
      real*8,dimension(:,:) :: pf(N_CELL,NZ)
      real*8,dimension(:,:) :: alphas(N_CELL,NZ)
      real*8,dimension(:,:) :: us(N_CELL,NZ)
      real*8,dimension(:,:) :: vs(N_CELL,NZ)
      real*8,dimension(:,:) :: ws(N_CELL,NZ)
      real*8,dimension(:,:) :: ps(N_CELL,NZ)

      real*8,dimension(:)   :: etan(N_CELL)    
      real*8,dimension(:,:) :: alphafn(N_CELL,NZ)
      real*8,dimension(:,:) :: ufn(N_CELL,NZ)
      real*8,dimension(:,:) :: vfn(N_CELL,NZ)
      real*8,dimension(:,:) :: wfn(N_CELL,NZ)
      real*8,dimension(:,:) :: pfn(N_CELL,NZ)
      real*8,dimension(:,:) :: alphasn(N_CELL,NZ)
      real*8,dimension(:,:) :: usn(N_CELL,NZ)
      real*8,dimension(:,:) :: vsn(N_CELL,NZ)
      real*8,dimension(:,:) :: wsn(N_CELL,NZ)
      real*8,dimension(:,:) :: psn(N_CELL,NZ)

      real*8,dimension(:)   :: xc(N_CELL),yc(N_CELL),sig(NZ)
      real*8,dimension(:)   :: Hpr(N_CELL),h(N_CELL)   
      real*8,dimension(:)   :: xv(N_VERT),yv(N_VERT),sigv(NZ-1)
      real*8,dimension(:)   :: Hprv(N_VERT),hv(N_VERT)   
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)     

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '----> Begin subroutine: LoopRKInitial'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Updating                                 !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |              Main variables: (known)=(n)               |
!     |________________________________________________________|

      DO i=1,N_CELL
!        ----------------------------
!        Free surface  
         eta(i)=etan(i)
         do k=1,NZ
!           ----------------------------
!           Volume Fraction  
            alphaf(i,k) = alphafn(i,k)
            alphas(i,k) = alphasn(i,k)
!           ----------------------------
!           Velocity
            uf(i,k) = ufn(i,k)
            vf(i,k) = vfn(i,k)
            wf(i,k) = wfn(i,k)
            us(i,k) = usn(i,k)
            vs(i,k) = vsn(i,k)
            ws(i,k) = wsn(i,k)
!           ----------------------------
!           Pressure
            pf(i,k) = pfn(i,k)
            ps(i,k) = psn(i,k)
!           ----------------------------
!           Turbulance
        enddo  
      ENDDO 

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '<---- End   subroutine: LoopRKInitial'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	          END                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
