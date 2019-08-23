!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        Update of the RK step                        !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE LoopRKUpdate(eta,                          &
                              alphaf,uf,vf,wf,pf,           &
                              alphas,us,vs,ws,ps,           &  
                              etanp,                        &
                              alphafnp,ufnp,vfnp,wfnp,pfnp, &
                              alphasnp,usnp,vsnp,wsnp,psnp, &
                              xc,yc,sig,Hpr,h,              &
                              xv,yv,sigv,Hprv,hv,           &
                              No_cp,nbe)     

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the values of the main variables for        !
!    the RK step.                                                     ! 
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
!  | --> alphafnp |(N_CELL,NZ)| Fluid control volume         t(n+1)|  !
!  | --> ufnp     |(N_CELL,NZ)| Velocity component u_f       t(n+1)|  !
!  | --> vfnp     |(N_CELL,NZ)| Velocity component v_f       t(n+1)|  !      
!  | --> wfnp     |(N_CELL,NZ)| Velocity component w_f       t(n+1)|  !
!  | --> pfnp     |(N_CELL,NZ)| Pressure of the fluid        t(n+1)|  !
!  | --> alphasnp |(N_CELL,NZ)| Solid control volume         t(n+1)|  !
!  | --> usnp     |(N_CELL,NZ)| Velocity component u_s       t(n+1)|  !
!  | --> vsnp     |(N_CELL,NZ)| Velocity component v_s       t(n+1)|  !     
!  | --> wsnp     |(N_CELL,NZ)| Velocity component w_s       t(n+1)|  !
!  | --> psnp     |(N_CELL,NZ)| Pressure of the solid        t(n+1)|  !
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

      real*8,dimension(:)   :: etanp(N_CELL)    
      real*8,dimension(:,:) :: alphafnp(N_CELL,NZ)
      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: alphasnp(N_CELL,NZ)
      real*8,dimension(:,:) :: usnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vsnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wsnp(N_CELL,NZ)
      real*8,dimension(:,:) :: psnp(N_CELL,NZ)

      real*8,dimension(:)   :: xc(N_CELL),yc(N_CELL),sig(NZ)
      real*8,dimension(:)   :: Hpr(N_CELL),h(N_CELL)   
      real*8,dimension(:)   :: xv(N_VERT),yv(N_VERT),sigv(NZ-1)
      real*8,dimension(:)   :: Hprv(N_VERT),hv(N_VERT)   
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)   

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '----> Begin subroutine: LoopRKUpdate'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Updating                                 !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |              Main variables: (known)=(n+1)             |
!     |________________________________________________________|

      DO i=1,N_CELL
!        ----------------------------
!        Free surface  
         eta(i)=etanp(i)
         do k=1,NZ
!           ----------------------------
!           Volume Fraction  
            alphaf(i,k) = alphafnp(i,k)
            alphas(i,k) = alphasnp(i,k)
!           ----------------------------
!           Velocity
            uf(i,k) = ufnp(i,k)
            vf(i,k) = vfnp(i,k)
            wf(i,k) = wfnp(i,k)
            us(i,k) = usnp(i,k)
            vs(i,k) = vsnp(i,k)
            ws(i,k) = wsnp(i,k)
!           ----------------------------
!           Pressure
            pf(i,k) = pfnp(i,k)
            ps(i,k) = psnp(i,k)
!           ----------------------------
!           Turbulance
        enddo  
      ENDDO 

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '<---- End   subroutine: LoopRKUpdate'
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
