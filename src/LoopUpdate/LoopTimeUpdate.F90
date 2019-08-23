!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      Finalization one time step                     !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE LoopTimeUpdate(etan,                         &
                                alphafn,ufn,vfn,wfn,pfn,      &
                                alphasn,usn,vsn,wsn,psn,      & 
                                xct,yct,zct,                  &  
                                xvt,yvt,zvt,                  &                                 
                                etanp,                        &
                                alphafnp,ufnp,vfnp,wfnp,pfnp, &
                                alphasnp,usnp,vsnp,wsnp,psnp, &
                                xc,yc,sig,Hpr,h,              &
                                xv,yv,sigv,Hprv,hv,           &
                                No_cp,nbe)     

!---------------------------------------------------------------------!
!                                                                     !
!    This program updates the values of the main variables for        !
!    new time step simulation or for the exit saved file. It          !
!    also calculates the criteria to stop the input box values        !
!    according to the amount of mass released.                        ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !  
!  |______________|___________|____________________________________|  !
!  | <-- etan     |(N_CELL)   | Free surface level eta at (n)      |  !
!  |______________|___________|____________________________________|  !
!  | <-- alphafn  |(N_CELL,NZ)| Fluid control volume at t(n)       |  !
!  | <-- ufn      |(N_CELL,NZ)| Velocity component u_f at t(n)     |  !
!  | <-- vfn      |(N_CELL,NZ)| Velocity component v_f at t(n)     |  !      
!  | <-- wfn      |(N_CELL,NZ)| Velocity component w_f at t(n)     |  !
!  | <-- pfn      |(N_CELL,NZ)| Pressure of the fluid at t(n)      |  !
!  | <-- alphasn  |(N_CELL,NZ)| Solid control volume at t(n)       |  !
!  | <-- usn      |(N_CELL,NZ)| Velocity component u_s at t(n)     |  !
!  | <-- vsn      |(N_CELL,NZ)| Velocity component v_s at t(n)     |  !     
!  | <-- wsn      |(N_CELL,NZ)| Velocity component w_s at t(n)     |  !
!  | <-- psn      |(N_CELL,NZ)| Pressure of the solid at t(n)      |  !
!  |______________|___________|____________________________________|  !
!  | <-- xct      |(N_CELL,NZ)| coordinate y at element center     |  !
!  | <-- yct      |(N_CELL,NZ)| coordinate x at element center     |  !
!  | <-- zct      |(N_CELL,NZ)| coordinate z at element center     |  !
!  | <-- xvt      |(N_VERT,NZ)| coordinate y at element vertex     |  !
!  | <-- yvt      |(N_VERT,NZ)| coordinate x at element vertex     |  !
!  | <-- zvt      |(N_VERT,NZ)| coordinate z at element vertex     |  !
!  |______________|___________|____________________________________|  !
!  | <-- mask     |(N_CELL,NZ)| Profile of the input box velocity  |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !      
!  |______________|___________|____________________________________|  !
!  | --> etanp    | N_CELL    | Free surface level eta at (n+1)    |  !
!  |______________|___________|____________________________________|  !
!  | --> alphafnp |(N_CELL,NZ)| Fluid control volume at t(n+1)     |  !
!  | --> ufnp     |(N_CELL,NZ)| Velocity component u_f at t(n+1)   |  !
!  | --> vfnp     |(N_CELL,NZ)| Velocity component v_f at t(n+1)   |  !      
!  | --> wfnp     |(N_CELL,NZ)| Velocity component w_f at t(n+1)   |  !
!  | --> pfnp     |(N_CELL,NZ)| Pressure of the fluid at t(n+1)    |  !
!  | --> alphasnp |(N_CELL,NZ)| Solid control volume at t(n+1)     |  !
!  | --> usnp     |(N_CELL,NZ)| Velocity component u_s at t(n+1)   |  !
!  | --> vsnp     |(N_CELL,NZ)| Velocity component v_s at t(n+1)   |  !     
!  | --> wsnp     |(N_CELL,NZ)| Velocity component w_s at t(n+1)   |  !
!  | --> psnp     |(N_CELL,NZ)| Pressure of the solid at t(n+1)    |  !
!  |______________|___________|____________________________________|  !
!  | --> xc       | N_CELL    | x-coordinate of the cell center    |  !
!  | --> yc       | N_CELL    | y-coordinate of the cell center    |  !
!  | --> sig      | NZ        | sigma                              |  !
!  | --> Hpr      | N_CELL    | Total water depth: H = h + eta     |  !
!  | --> h        | N_CELL    | Depth of the domain at each cell   |  !
!  | --> xv       | N_VERT    | x-coordinate of the cell vertex    |  !
!  | --> yv       | N_VERT    | y-coordinate of the cell vertex    |  !   
!  | --> sigv     | NZ-1      | sigma of the vertex points         |  !
!  | --> Hprv     | N_VERT    | Total water depth at the vertex    |  !
!  | --> hv       | N_VERT    | Depth of the domain at each vertex |  !
!  |______________|___________|____________________________________|  !
!  | --> No_cp    |(N_CELL,3) | Node No. of surrounding three cell |  !
!  | --> nbe      | N_CELL    | Type of boundary cell: 0,1,2,3     |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  ! 
!  | --- N_CELL  | Total number of the cells                       |  !
!  | --- N_CELL0 | Number of inside cells                          |  !
!  |     NZ      | Points in the sigma direction                   |  !
!  |  rmassrejet | Total mass release in the simulation            |  ! 
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  | * i,k       |  Loop counters                                  |  !    
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  |   rmass     | Total mass released at the current time         |  !  
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

      real*8,dimension(:)   :: etan  
      real*8,dimension(:,:) :: alphafn,ufn,vfn,wfn,pfn
      real*8,dimension(:,:) :: alphasn,usn,vsn,wsn,psn
      real*8,dimension(:,:) :: xct,yct,zct  
      real*8,dimension(:,:) :: xvt,yvt,zvt   
      real*8,dimension(:)   :: etanp    
      real*8,dimension(:,:) :: alphafnp,ufnp,vfnp,wfnp,pfnp
      real*8,dimension(:,:) :: alphasnp,usnp,vsnp,wsnp,psnp
      real*8,dimension(:)   :: xc,yc,sig,Hpr,h   
      real*8,dimension(:)   :: xv,yv,sigv,Hprv,hv   
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe    

!      ________________________________________________________
!     |                                                        |
!     |            Definition of local variables               |
!     |________________________________________________________|

      real*8 :: rmass

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '----> Begin subroutine: LoopTimeUpdate'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Updating                                 !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                 Depth H, h for the box case            |
!     |________________________________________________________|

      if (ChooseDomBox.eq.1) then
         do i=1,N_CELL 
            h(i)   = 0.0d0
            Hpr(i) = 1.0d0  
         enddo
         do nv=1,N_VERT 
            hv(nv)   = 0.0d0
            Hprv(nv) = 1.0d0
         enddo
      endif

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
!     |                Main variables: (n)=(n+1)               |
!     |________________________________________________________|

      DO i=1,N_CELL
!        ----------------------------
!        Free surface  
         etan(i)=etanp(i)
         do k=1,NZ
!           ----------------------------
!           Volume Fraction  
            alphafn(i,k) = alphafnp(i,k)
            alphasn(i,k) = alphasnp(i,k)
!           ----------------------------
!           Velocity
            ufn(i,k) = ufnp(i,k)
            vfn(i,k) = vfnp(i,k)
            wfn(i,k) = wfnp(i,k)            
            usn(i,k) = usnp(i,k)
            vsn(i,k) = vsnp(i,k)
            wsn(i,k) = wsnp(i,k)
!           ----------------------------
!           Pressure
            pfn(i,k) = pfnp(i,k)
            psn(i,k) = psnp(i,k)
        enddo  
      ENDDO 

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '<---- End   subroutine: LoopTimeUpdate'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	          END                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
