!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        Finalization RK method                       !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE LoopRKFinal(etanp,                        &
                             alphafnp,ufnp,vfnp,wfnp,pfnp, &
                             alphasnp,usnp,vsnp,wsnp,psnp, &
                             etan,                         &
                             alphafn,ufn,vfn,wfn,pfn,      &
                             alphasn,usn,vsn,wsn,psn,      &
                             eta,                          &
                             alphaf,uf,vf,wf,pf,           &
                             alphas,us,vs,ws,ps,           &  
                             xc,yc,sig,Hpr,h,              &
                             xv,yv,sigv,Hprv,hv,           &
                             No_cp,nbe)     

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the final values of the main variables   !
!    of the Runge-Kutta method.                                       ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !  
!  |______________|___________|____________________________________|  !
!  | <-- etanp    |(N_CELL)   | Free surface level eta       t(n+1)|  !
!  |______________|___________|____________________________________|  !
!  | <-- alphafnp |(N_CELL,NZ)| Fluid control volume         t(n+1)|  !
!  | <-- ufnp     |(N_CELL,NZ)| Velocity component u_f       t(n+1)|  !
!  | <-- vfnp     |(N_CELL,NZ)| Velocity component v_f       t(n+1)|  !      
!  | <-- wfnp     |(N_CELL,NZ)| Velocity component w_f       t(n+1)|  !
!  | <-- pfnp     |(N_CELL,NZ)| Pressure of the fluid        t(n+1)|  !
!  | <-- alphasnp |(N_CELL,NZ)| Solid control volume         t(n+1)|  !
!  | <-- usnp     |(N_CELL,NZ)| Velocity component u_s       t(n+1)|  !
!  | <-- vsnp     |(N_CELL,NZ)| Velocity component v_s       t(n+1)|  !     
!  | <-- wsnp     |(N_CELL,NZ)| Velocity component w_s       t(n+1)|  !
!  | <-- psnp     |(N_CELL,NZ)| Pressure of the solid        t(n+1)|  !
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
!  | --> eta      |(N_CELL)   | Free surface level eta      (known)|  !
!  |______________|___________|____________________________________|  !
!  | --> alphaf   |(N_CELL,NZ)| Fluid control volume        (known)|  !
!  | --> uf       |(N_CELL,NZ)| Velocity component u_f      (known)|  !
!  | --> vf       |(N_CELL,NZ)| Velocity component v_f      (known)|  !      
!  | --> wf       |(N_CELL,NZ)| Velocity component w_f      (known)|  !
!  | --> pf       |(N_CELL,NZ)| Pressure of the fluid       (known)|  !
!  | --> alphas   |(N_CELL,NZ)| Solid control volume        (known)|  !
!  | --> us       |(N_CELL,NZ)| Velocity component u_s      (known)|  !
!  | --> vs       |(N_CELL,NZ)| Velocity component v_s      (known)|  !     
!  | --> ws       |(N_CELL,NZ)| Velocity component w_s      (known)|  !
!  | --> ps       |(N_CELL,NZ)| Pressure of the solid       (known)|  !
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

      real*8,dimension(:)   :: xc(N_CELL),yc(N_CELL),sig(NZ)
      real*8,dimension(:)   :: Hpr(N_CELL),h(N_CELL)   
      real*8,dimension(:)   :: xv(N_VERT),yv(N_VERT),sigv(NZ-1)
      real*8,dimension(:)   :: Hprv(N_VERT),hv(N_VERT)   
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)   

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '----> Begin subroutine: LoopRKFinal'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Updating                                 !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Main variables: (n+1)= (n)+(known)         |
!     |________________________________________________________|

      DO i=1,N_CELL
!        _______________________________
!       |                               |
!       |        Inside elements        |
!       |_______________________________|
!        ----------------------------
!        Free surface  
         etanp(i) = 0.5d0*(etan(i)+eta(i))
         do k=1,NZ
!           ----------------------------
!           Volume Fraction  
            alphafnp(i,k) = 0.5d0*(alphafn(i,k)+alphaf(i,k))
            alphasnp(i,k) = 0.5d0*(alphasn(i,k)+alphas(i,k))
!           ----------------------------
!           Velocity
            ufnp(i,k) = 0.5d0*(ufn(i,k)+uf(i,k)) 
            vfnp(i,k) = 0.5d0*(vfn(i,k)+vf(i,k))
            wfnp(i,k) = 0.5d0*(wfn(i,k)+wf(i,k))
            usnp(i,k) = 0.5d0*(usn(i,k)+us(i,k))
            vsnp(i,k) = 0.5d0*(vsn(i,k)+vs(i,k))
            wsnp(i,k) = 0.5d0*(wsn(i,k)+ws(i,k))
!           ----------------------------
!           Pressure
            pfnp(i,k) = 0.5d0*(pfn(i,k)+pf(i,k))
            psnp(i,k) = 0.5d0*(psn(i,k)+ps(i,k))
!           ----------------------------
            ufnp(i,k) = uf(i,k) 
            vfnp(i,k) = vf(i,k)
            wfnp(i,k) = wf(i,k)
            pfnp(i,k) = pf(i,k)
        enddo 
      ENDDO 

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '<---- End   subroutine: LoopRKFinal'
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
