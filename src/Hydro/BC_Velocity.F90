!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             3D BOUNDARY CONDITION CELL-CENTERS & VERTEX             !
!                             Nov 2017                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCvelocity3D(Hphi,phi,                      &
                              Hphiv,phiv,                    &
!                             --------------------------------
                              Hpr,eta,                       &
                              Hprv,etav,                     &
                              h,hv,                          &
!                             --------------------------------
                              xc,yc,sig,dsig,No_cp,nbe,      &
                              xv,yv,sigv,dsigv,No_vp,nbev,   &
!                             --------------------------------
                              tagBC)
                                 
!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center & vertex bound.  !
!    condition of the velocity. The tag called "tagBC" is used to     !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <--- phi    |(N_CELL,NZ)| Function at the cell-center         |  !
!  | <--- phiv   |(N_VERT,NZ-1)| Function at the vertices          |  !  
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | ---> xc,yc  |(N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> sig    |(NZ)       | sigma value at the cell centers     |  !
!  | ---> dsig   |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | ---> No_cp  |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe    |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | ---> xv,yv  |(N_VERT)   | Coordinates of the vertices         |  !
!  | ---> sigv   |(NZ-1)     | sigma value at the vertices         |  !
!  | ---> dsigv  |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | ---> No_vp  |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | ---> nbev   |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
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

      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: hv(N_VERT)
            
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
          
      integer :: tagBC      
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB
      real*8 :: z1,z2,z3,f1,f2,f3,a1,a2,a3
      real*8 :: funExamNSu,funExamNSv,funExamNSw,funInflow
      real*8 :: nxc,nyc,taux,tauy,ntau
      integer:: elem,ii,jj,jv1,jv2,jv3,i0,nvp
!      ____________________________________
!     |                                    |
!     |   Declaration of parameters        |
!     |____________________________________|

      integer, parameter :: DoImpermeability = 1 !(Free-slip BC at wall)
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCvelocity3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!        Verify the proper conditions for the different cases         !
!                                                                     !
!*********************************************************************! 

      IF (time.eq.0.0d0) THEN
         print*,'   ---->> VERIFYING BOUNDARY TAGS nbe(i):     <<<---'
         DO i=1,N_CELL0
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL INFLOW (nbe=2)     |
!          |___________________________________|                       

            if (nbe(i).eq.2) then
#              if defined(KeyStandingWave) 
                   print*,'StandingWave: all BC should be wall nbe=1'
                   stop
#              elif defined(KeyTaylorVortex)
                   print*,'TaylorVortex: all BC should be wall nbe=1'
                   stop            
#              endif
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL OUTFLOW (nbe=3)    |
!          |___________________________________| 

            elseif (nbe(i).eq.3) then
#              if defined(KeyStandingWave) 
                   print*,'StandingWave: all BC should be wall nbe=1'
                   stop
#              elif defined(KeyTaylorVortex)
                   print*,'TaylorVortex: all BC should be wall nbe=1'
                   stop            
#              endif
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL STRUCTURE (nbe=4)  |
!          |___________________________________| 

            elseif (nbe(i).eq.4) then
#              if defined(KeyEstuaryGironde)
                   print*,'StandingWave: no structure nbe=/4'
                   stop 
#              elif defined(KeyStaticChannel)
                   print*,'StaticChannel: no structure nbe=/4'
                   stop
#              elif defined(KeyStandingWave) 
                   print*,'StandingWave: no structure nbe=/4'
                   stop
#              elif defined(KeyTaylorVortex)
                   print*,'TaylorVortex: no structure nbe=/4'
                   stop
#              endif
            endif
         ENDDO
         print*,'   ---->> VERIFYING BOUNDARY TAGS nbe(i): OK! <<<---'
      ENDIF

!*********************************************************************!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!                                                                     !
!               Velocity component  Hu & u (cell-center)              !
!                                                                     !
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!*********************************************************************!

      IF (tagBC.eq.1) THEN
      
         DO i=1,N_CELL
            x = xc(i)
            y = yc(i)
!           ___________________________________
!          |                                   |
!          |       VERTICAL BOTTOM (k=1)       |
!          |___________________________________|                     

!           ___________________________________                               
!           [a] SLIP BOTTOM: ["du/dn=0",dv/dn=0,w=0]
#           ifdef KeyBCbottomSlip                 
               phi(i,1)  = phi(i,2)   !<<<-- du/dn=0
               Hphi(i,1) = Hphi(i,2) 
#           endif
!           ___________________________________                               
!           [b] NO-SLIP BOTTOM: ["u=0",v=0,w=0]                
#           ifdef KeyBCbottomNoSlip                 
               phi(i,1)  = -phi(i,2)  !<<<-- u=0
               Hphi(i,1) = -Hphi(i,2)
#           endif
!           ___________________________________                               
!           [c.1] EXACT BOTTOM:["u=fB",v=fB,wf=fB] 1st order 
#           ifdef KeyBCbottomExact_1st
               z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)                         
               fB = funExamNSu(x,y,z,time)
               phi(i,1)  = 2.0d0*fB-phi(i,2) !<<<-- u=fB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif
!           ___________________________________                               
!           [c.2] EXACT BOTTOM:["u=fB",v=fB,wf=fB] 2nd order
#           ifdef KeyBCbottomExact_2nd 
               k  = 1
               z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)                         
               fB = funExamNSu(x,y,z,time) 
               z1 = sig(1)*Hpr(i)-h(i)
               z2 = sig(2)*Hpr(i)-h(i)
               z3 = sig(3)*Hpr(i)-h(i)
               f2 = phi(i,2)
               f3 = phi(i,3)
               a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
               a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
               a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
               f1 = (fB-a2*f2-a3*f3)/a1
               phi(i,1) = f1 !<<<-- u=fB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif                         

!           ___________________________________
!          |                                   |
!          |        VERTICAL TOP (k=NZ)        |
!          |___________________________________|              

!           ___________________________________                               
!           [a] SLIP TOP: ["du/dn=0",dv/dn=0,w=0]
#           if defined(KeyBCtopSlip)  
                phi(i,NZ)  = phi(i,NZ-1)  !<<<-- du/dn=0
                Hphi(i,NZ) = Hphi(i,NZ-1)
#           endif
!           ___________________________________                               
!           [b] NO-SLIP TOP: ["u=0",v=0,w=0]
#           if defined(KeyBCtopNoSlip)  
                phi(i,NZ)  = -phi(i,NZ-1)  !<<<-- u=0
                Hphi(i,NZ) = -Hphi(i,NZ-1)
#           endif
!           ___________________________________                               
!           [c.1] EXACT TOP:["u=fB",v=fB,wf=fB] 1st order 
#           ifdef KeyBCtopExact_1st
               z  = 0.5d0*(sig(NZ)+sig(NZ-1))*Hpr(i)-h(i)            
               fB = funExamNSu(x,y,z,time)
               phi(i,NZ)  = 2.0d0*fB-phi(i,NZ-1) !<<<-- u=fB
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i)
#           endif
!           ___________________________________                               
!           [c.2] EXACT TOP:["u=fB",v=fB,wf=fB] 2nd order
#           ifdef KeyBCtopExact_2nd 
               k  = NZ
               z  = 0.5d0*(sig(NZ)+sig(NZ-1))*Hpr(i)-h(i)            
               fB = funExamNSu(x,y,z,time)
               z1 = sig(NZ-2)*Hpr(i)-h(i)
               z2 = sig(NZ-1)*Hpr(i)-h(i)
               z3 = sig(NZ)*Hpr(i)-h(i)
               f1 = phi(i,NZ-2)
               f2 = phi(i,NZ-1)
               a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
               a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
               a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
               f3 = (fB-a1*f1-a2*f2)/a3
               phi(i,NZ) = f3 !<<<-- u=fB
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i)
#           endif            

         ENDDO

         DO i=1,N_CELL0
!           ___________________________________
!          |                                   |
!          |      HORIZONTAL WALL (nbe=1)      |
!          |___________________________________|                       

            if (nbe(i).eq.1) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ______________________________________________                               
!                    [a.1] SLIP WALL X-DIRECTION: ["du/dn=0",v=0,dw/dn=0]                                                     		
#                    ifdef KeyBCwallXSlip
                        do k=1,NZ 
                           phi(nc,k)  = phi(i,k)    !<--- du/dn=0
                           Hphi(nc,k) = Hphi(i,k)
                        enddo
#                    endif
!                    ______________________________________________                               
!                    [a.2] SLIP WALL Y-DIRECTION: ["u=0",dv/dn=0,dw/dn=0]                                                     		
#                    ifdef KeyBCwallYSlip
                        do k=1,NZ 
                           phi(nc,k)  = -phi(i,k)    !<--- u=0
                           Hphi(nc,k) = -Hphi(i,k)
                        enddo
#                    endif
!                    ______________________________________________                               
!                    [b] NO-SLIP WALL: ["u=0",v=0,w=0]                                       
#                    ifdef KeyBCwallNoSlip
                        do k=1,NZ
                           phi(nc,k)  = -phi(i,k)   !<--- u=0
                           Hphi(nc,k) = -Hphi(i,k)
                        enddo                      
#                    endif
!                    ______________________________________________                                
!                    [c] EXACT WALL: ["u=fB",v=fB,wf=fB]                      
#                    ifdef KeyBCwallExact
                        x = xe(i,j)
                        y = ye(i,j)
                        do k=1,NZ                            
                           z  = sig(k)*Hpr(i)-h(i)                       
                           fB = funExamNSu(x,y,z,time)
                           phi(nc,k) = 2.0d0*fB-phi(i,k)
                           Hphi(nc,k) = phi(nc,k)*Hpr(i)
                        enddo  
#                    endif
!                    ______________________________________________ 
!                    [d] IMPERMEABILITY WALL:[*,*,dw/dn=0]
#                    ifdef KeyBCwallImpermeability
                        jj=j+1
                        if (jj.gt.3) jj=jj-3
                        jv1 = No_vp(i,j)
                        jv2 = No_vp(i,jj)
                        taux = xv(jv2)-xv(jv1) !dxVV(i,j)!
                        tauy = yv(jv2)-yv(jv1) !dyVV(i,j)!
                        ntau = sqrt(taux*taux+tauy*tauy)
                        nxc =  tauy/ntau
                        nyc = -taux/ntau  
                        if (nyc.eq.0) then                  
                           do k=1,NZ
                              phi(nc,k)  = -phi(i,k)  !<--- u=0
                              Hphi(nc,k) = -Hphi(i,k) 
                           enddo 
                        else
                           do k=1,NZ
                              phi(nc,k)  = phi(i,k)   !<--- du/dtau=0
                              Hphi(nc,k) = Hphi(i,k)
                           enddo      
                        endif  
#                    endif
!                    ___________________________________
!                    [e] PERIODIC WALL X-DIRECTION
#                    ifdef KeyBCperiodicX
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif                                              
!                    ___________________________________
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL INFLOW (nbe=2)     |
!          |___________________________________|                

            elseif (nbe(i).eq.2) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    EXACT INFLOW: u=fB
#                    ifdef KeyBCinflowExact
                        do k=1,NZ 
                           x = xe(i,j)
                           y = ye(i,j)
                           z = sig(k)*Hpr(i)-h(i) 
!                          --------------------------
                           fB = funInflow(x,y,z,time)
!                          --------------------------
                           phi(nc,k) = 2.0d0*fB-phi(i,k) 
                           Hphi(nc,k) = phi(nc,k)*Hpr(i)
                        enddo
#                    endif
!                    ___________________________________
!                    PERIODIC INFLOW Y-DIRECTION
#                    ifdef KeyBCperiodicY
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif
!                    ___________________________________ 
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL OUTFLOW (nbe=3)    |
!          |___________________________________|                

            elseif (nbe(i).eq.3) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    DEFAULT NEUMANN OUTFLOW: ["du/dn",dv/dn=0,dw/dn=0]
#                    ifdef KeyBCoutflowNeumann
                        do k=1,NZ
                           phi(nc,k) = phi(i,k)  !<<<-- du/dn=0
                           Hphi(nc,k)= Hphi(i,k)
                        enddo
#                    endif
!                    ___________________________________
!                    PERIODIC OUTFLOW Y-DIRECTION
#                    ifdef KeyBCperiodicY
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif
!                    ___________________________________
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL STRUCTURE (nbe=4)  |
!          |___________________________________|                

            elseif (nbe(i).eq.4) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    NO-SLIP STRUCTURE: ["u=0",v=0,w=0]
                     do k=1,NZ                    
                        phi(nc,k)  = -phi(i,k) !<<<-- u=0
                        Hphi(nc,k) = -Hphi(i,k)
                     enddo
!                    ___________________________________
                  endif 
               enddo
            endif                                  
         ENDDO
      
!*********************************************************************!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                     !
!                    Velocity component u (vertex)                    !
!                                                                     !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!*********************************************************************! 

      ELSEIF (tagBC.eq.-1) THEN

!         ___________________________________
!        |                                   |
!        |        (Impermeability)           |
!        |        HORIZONTAL (nbe=1)         |
!        |___________________________________| 

!        ___________________________________
!        IMPERMEABILITY WALL
#        ifdef KeyBCwallImpermeability
            DO i=1,N_CELL0
               x = xc(i)
               y = yc(i)
               if (nbe(i).eq.1) then
                  do j=1,3
                     nc=No_cp(i,j)
!                    ====================================
!                    ==========  SEQUENTIAL =============
#                    ifndef KeyParallel
                     if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                    ====================================
!                    =====  START PARALLEL OPTION =======
#                    else
                     elem = index_global(nc)
                     if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                    endif
!                    =============== END ================
!                    ====================================
                        jj=j+1
                        if (jj.gt.3) jj=jj-3
                        jv1 = No_vp(i,j)
                        jv2 = No_vp(i,jj)
                        taux = xv(jv2)-xv(jv1) 
                        tauy = yv(jv2)-yv(jv1) 
                        ntau = sqrt(taux*taux+tauy*tauy)
                        nxc =  tauy/ntau
                        nyc = -taux/ntau  
                        if (nyc.eq.0) then
                            do k=1,NZ-1
                               phiv(jv1,k)  = 0.0d0
                               phiv(jv2,k)  = 0.0d0
                               Hphiv(jv1,k) = 0.0d0
                               Hphiv(jv2,k) = 0.0d0 
                            enddo
                        endif
                     endif 
                  enddo
               endif
            ENDDO
#        endif

!         ___________________________________
!        |                                   |
!        |    VERTEX VERTICAL BOTTOM (k=1)   |
!        |___________________________________|  

!        ___________________________________
!        NO-SLIP BOTTOM: ["u=0",v=0,w=0] 
#        ifdef KeyBCbottomNoSlip 
            DO nv=1,N_VERT                
               phiv(nv,1)  = 0.0d0  !<<<-- u=0
               Hphiv(nv,1) = 0.0d0
            ENDDO
#        endif
!        ___________________________________
!        EXACT BOTTOM:["u=fB",v=fB,wf=fB] 1st or 2nd order
#        if defined( KeyBCbottomExact_1st) || defined(KeyBCbottomExact_2nd)
            DO nv=1,N_VERT 
               k = 1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               phiv(nv,k)  = funExamNSu(x,y,z,time)
               Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
            ENDDO                
#        endif              
!         ___________________________________
!        |                                   |
!        |    VERTEX VERTICAL TOP (k=NZ-1)   |
!        |___________________________________|

!        ___________________________________
!        NO-SLIP TOP: ["u=0",v=0,w=0] 
#        ifdef KeyBCtopNoSlip 
            DO nv=1,N_VERT                
               phiv(nv,NZ-1)  = 0.0d0  !<<<-- u=0
               Hphiv(nv,NZ-1) = 0.0d0
            ENDDO
#        endif
!        ___________________________________
!        EXACT TOP:["u=fB",v=fB,wf=fB] 1st or 2nd order
#        if defined( KeyBCtopExact_1st) || defined(KeyBCtopExact_2nd)
            DO nv=1,N_VERT 
               k = NZ-1
               x = xv(nv)
               y = yv(nv)                
               z = sigv(k)*Hprv(nv)-hv(nv)  
               phiv(nv,k)  = funExamNSu(x,y,z,time)
               Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
            ENDDO
#        endif
                  
         DO nv=1,N_VERT
!            ___________________________________
!           |                                   |
!           |      HORIZONTAL WALL (nbev=1)     |
!           |___________________________________|

            if (nbev(nv).eq.1) then
!              ______________________________________________
!              SLIP WALL Y-DIRECTION: ["u=0",dv/dn=0,dw/dn=0]
#              ifdef KeyBCwallYSlip
                  do k=1,NZ-1
                     phiv(nv,k)  = 0.0d0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              NO-SLIP WALL: ["u=0",v=0,w=0]
#              ifdef KeyBCwallNoSlip
                  do k=1,NZ-1
                     phiv(nv,k)  = 0.0d0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              EXACT WALL: ["u=fB",v=fB,w=fB]
#              ifdef KeyBCwallExact
                  do k=1,NZ-1 
                     x = xv(nv)
                     y = yv(nv)
                     z = sigv(k)*Hprv(nv)-hv(nv)
                     phiv(nv,k)  = funExamNSu(x,y,z,time)
                     Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
                  enddo
#              endif
!              ___________________________________
!              PERIODIC WALL X-DIRECTION
#              ifdef KeyBCperiodicX
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
                  if (TagPeriodicBC(nv,1).eq.1) then
                     nvp = TagPeriodicBC(nv,2)
                     do k=1,NZ-1 
                        phiv(nv,k)  = phiv(nvp,k)
                        Hphiv(nv,k) = Hphiv(nvp,k)
                     enddo
                  endif
!#                 endif
!                 NOTE: In parallel only works if the 
!                       corresponding vertex is in the
!                       same sub-domain.
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL INFLOW (nbev=2)    |
!           |___________________________________|

            elseif (nbev(nv).eq.2) then
!              ___________________________________
!              EXACT INFLOW: ["u=fB",v=0,w=0]
#              ifdef KeyBCinflowExact
                  do k=1,NZ-1 
                     x = xv(nv)
                     y = yv(nv)
                     z = sigv(k)*Hprv(nv)-hv(nv)
                     phiv(nv,k) = funInflow(x,y,z,time)
                     Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
                 enddo
#             endif
!              ___________________________________
!              PERIODIC INFLOW Y-DIRECTION
#              ifdef KeyBCperiodicY
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
                  if (TagPeriodicBC(nv,1).eq.1) then
                     nvp = TagPeriodicBC(nv,2)
                     do k=1,NZ-1 
                        phiv(nv,k)  = phiv(nvp,k)
                        Hphiv(nv,k) = Hphiv(nvp,k)
                     enddo
                  endif
!#                 endif
!                 NOTE: In parallel only works if the 
!                       corresponding vertex is in the
!                       same sub-domain.
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL OUTFLOW (nbev=3)   |
!           |___________________________________|

            elseif (nbev(nv).eq.3) then
!              ___________________________________
!              PERIODIC OUTFLOW Y-DIRECTION
#              ifdef KeyBCperiodicY
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
!                  if (TagPeriodicBC(nv,1).eq.1) then
!                     nvp = TagPeriodicBC(nv,2)
!                     do k=1,NZ-1 
!                        phiv(nv,k)  = phiv(nvp,k)
!                        Hphiv(nv,k) = Hphiv(nvp,k)
!                     enddo
!                  endif
!#                 endif
!                 NOTE: Not necessary 
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL STRUCTURE (nbev=4) |
!           |___________________________________|

            elseif (nbev(nv).eq.4) then !<<<-- Structure: u=0
!              ___________________________________
!              NO-SLIP STRUCTURE: ["u=0",v=0,w=0]
               do k=1,NZ-1 
                  phiv(nv,k)  = 0.0d0
                  Hphiv(nv,k) = 0.0d0
               enddo
            endif

         ENDDO

!*********************************************************************!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!                                                                     !
!               Velocity component  Hv & v (cell-center)              !
!                                                                     !
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!*********************************************************************!

      ELSEIF (tagBC.eq.2) THEN
      
         DO i=1,N_CELL
            x = xc(i)
            y = yc(i)
!           ___________________________________
!          |                                   |
!          |       VERTICAL BOTTOM (k=1)       |
!          |___________________________________| 

!           ___________________________________                               
!           SLIP BOOTOM: [du/dz,"dv/dz=0",w=0]
#           ifdef KeyBCbottomSlip                
               phi(i,1)  = phi(i,2)  !<<<-- dv/dn=0
               Hphi(i,1) = Hphi(i,2)
#           endif
!           ___________________________________                               
!           NO-SLIP BOOTOM: [u=0,"v=0",w=0]                
#           ifdef KeyBCbottomNoSlip                   
               phi(i,1)  = -phi(i,2)  !<<<-- v=0
               Hphi(i,1) = -Hphi(i,2)
#           endif
!           ___________________________________                               
!           EXACT BOTTOM:[u=fB,"v=fB",wf=fB] 1st order 
#           ifdef KeyBCbottomExact_1st
               z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)            
               fB = funExamNSv(x,y,z,time) 
               phi(i,1)  = 2.0d0*fB-phi(i,2) !<<<-- v=fB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif
!           ___________________________________                               
!           EXACT BOTTOM:[u=fB,"v=fB",wf=fB] 2nd order
#           ifdef KeyBCbottomExact_2nd
               k  = 1
               z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)            
               fB = funExamNSv(x,y,z,time) 
               z1 = sig(1)*Hpr(i)-h(i)
               z2 = sig(2)*Hpr(i)-h(i)
               z3 = sig(3)*Hpr(i)-h(i)
               f2 = phi(i,2)
               f3 = phi(i,3)
               a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
               a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
               a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
               f1 = (fB-a2*f2-a3*f3)/a1
               phi(i,1) = f1 !<<<-- v=fB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif         
            
!           ___________________________________
!          |                                   |
!          |        VERTICAL TOP (k=NZ)        |
!          |___________________________________| 
  
!           _____________________________________
!           SLIP TOP: [du/dz=0,"dv/dz=0",w=0]
#           if defined(KeyBCtopSlip)
                phi(i,NZ)  = phi(i,NZ-1)  !<<<-- dv/dn=0
                Hphi(i,NZ) = Hphi(i,NZ-1)           
#           endif
!           ___________________________________                               
!           NO-SLIP TOP: [u=0,"v=0",w=0]
#           if defined(KeyBCtopNoSlip)
                phi(i,NZ)  = -phi(i,NZ-1)  !<<<-- v=0
                Hphi(i,NZ) = -Hphi(i,NZ-1)           
#           endif
!           ___________________________________                               
!           EXACT TOP:[u=fB,"v=fB",wf=fB] 1st order 
#           ifdef KeyBCtopExact_1st
               z = 0.5d0*(sig(NZ)+sig(NZ-1))*Hpr(i)-h(i)            
               fB = funExamNSv(x,y,z,time)
               phi(i,NZ)  = 2.0d0*fB-phi(i,NZ-1) !<<<-- v=fB
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i) 
#           endif
!           ___________________________________                               
!           EXACT TOP:[u=fB,"v=fB",wf=fB] 2nd order
#           ifdef KeyBCtopExact_2nd 
               k = NZ
               z = 0.5d0*(sig(NZ)+sig(NZ-1))*Hpr(i)-h(i)            
               fB = funExamNSv(x,y,z,time)
               z1 = sig(NZ-2)*Hpr(i)-h(i)
               z2 = sig(NZ-1)*Hpr(i)-h(i)
               z3 = sig(NZ)*Hpr(i)-h(i)
               f1 = phi(i,NZ-2)
               f2 = phi(i,NZ-1)
               a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
               a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
               a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
               f3 = (fB-a1*f1-a2*f2)/a3
               phi(i,NZ) = f3 !<<<-- v=fB
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i)
#           endif                   

         ENDDO

         DO i=1,N_CELL0
!           ___________________________________
!          |                                   |
!          |      HORIZONTAL WALL (nbe=1)      |
!          |___________________________________| 

            if (nbe(i).eq.1) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ______________________________________________
!                    SLIP WALL X-DIRECTION: [du/dn=0,"v=0",dw/dn=0]                                       
#                    ifdef KeyBCwallXSlip
                        do k=1,NZ
                           phi(nc,k)  = -phi(i,k) !<--- v=0
                           Hphi(nc,k) = -Hphi(i,k)
                        enddo
#                    endif
!                    ______________________________________________
!                    SLIP WALL Y-DIRECTION: [u=0,"dv/dn=0",dw/dn=0]                                       
#                    ifdef KeyBCwallYSlip
                        do k=1,NZ
                           phi(nc,k)  = phi(i,k) !<--- dv/dn=0
                           Hphi(nc,k) = Hphi(i,k)
                        enddo
#                    endif
!                    ______________________________________________                               
!                    NO-SLIP WALL: [u=0,"v=0",w=0]                                                        		
#                    ifdef KeyBCwallNoSlip                    
                        do k=1,NZ
                           phi(nc,k)  = -phi(i,k) !<--- v=0
                           Hphi(nc,k) = -Hphi(i,k)
                        enddo
#                    endif  
!                    ______________________________________________                      
!                    EXACT WALL: [u=fB,"v=fB",w=fB]                      
#                    ifdef KeyBCwallExact
                        do k=1,NZ 
                           x = xe(i,j)
                           y = ye(i,j)
                           z = sig(k)*Hpr(i)-h(i)                      
                           fB = funExamNSv(x,y,z,time)
                           phi(nc,k) = 2.0d0*fB-phi(i,k)
                           Hphi(nc,k) = phi(nc,k)*Hpr(i)
                        enddo 
#                    endif
!                    ______________________________________________
!                    IMPERMEABILITY WALL:[*,*,dw/dn=0]
#                    ifdef KeyBCwallImpermeability
                        jj=j+1
	                    if (jj.gt.3) jj=jj-3
	                    jv1 = No_vp(i,j)
	                    jv2 = No_vp(i,jj)
                        taux = xv(jv2)-xv(jv1) !dxVV(i,j)
	                    tauy = yv(jv2)-yv(jv1) !dyVV(i,j)	                        
	                    ntau = sqrt(taux*taux+tauy*tauy)
                        nxc =  tauy/ntau
                        nyc = -taux/ntau  
                        if (nxc.eq.0) then                  
                           do k=1,NZ
                              phi(nc,k)  = -phi(i,k)  !<--- v=0
                              Hphi(nc,k) = -Hphi(i,k)
                           enddo 
                        else
                           do k=1,NZ
                              phi(nc,k)  = phi(i,k)   !<--- dv/dn=0
                              Hphi(nc,k) = Hphi(i,k)
                           enddo      
                        endif  
#                    endif 
!                    ______________________________________________
!                    PERIODIC WALL X-DIRECTION
#                    ifdef KeyBCperiodicX
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif                                              
!                    ___________________________________ 
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL INFLOW (nbe=2)     |
!          |___________________________________| 

            elseif (nbe(i).eq.2) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    EXACT INFLOW: [u=fB,"v=0",w=0]
#                    ifdef KeyBCinflowExact
                        do k=1,NZ 
                           phi(nc,k) = -phi(i,k) !<<<-- v=0                       
                           Hphi(nc,k)= -Hphi(i,k)
                        enddo
#                    endif
!                    ___________________________________
!                    PERIODIC INFLOW Y-DIRECTION
#                    ifdef KeyBCperiodicY
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif                                              
!                    ___________________________________ 
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL OUTFLOW (nbe=3)    |
!          |___________________________________|

            elseif (nbe(i).eq.3) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    DEFAULT NEUMANN OUTFLOW: [du/dn,"dv/dn=0",dw/dn=0]
#                    ifdef KeyBCoutflowNeumann
                        do k=1,NZ
                           phi(nc,k) = phi(i,k)  !<<<-- dv/dn=0
                           Hphi(nc,k)= Hphi(i,k)
                        enddo
#                    endif
!                    ___________________________________
!                    PERIODIC OUTFLOW Y-DIRECTION
#                    ifdef KeyBCperiodicY
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif
!                    ___________________________________
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL STRUCTURE (nbe=4)  |
!          |___________________________________|

            elseif (nbe(i).eq.4) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    NO-SLIP STRUCTURE: [u=0,"v=0",w=0]
                     do k=1,NZ 
                        phi(nc,k) = -phi(i,k) !<<<-- v=0
                        Hphi(nc,k)= -Hphi(i,k)
                     enddo
!                    ___________________________________
                  endif 
               enddo                                             
            endif
         ENDDO

!*********************************************************************!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                     !
!                    Velocity component v (vertex)                    !
!                                                                     !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!*********************************************************************! 

      ELSEIF (tagBC.eq.-2) THEN

!         ___________________________________
!        |                                   |
!        |        (Impermeability)           |
!        |        HORIZONTAL (nbev=1)        |
!        |___________________________________| 

!        ___________________________________
!        IMPERMEABILITY WALL
#        ifdef KeyBCwallImpermeability
            DO i=1,N_CELL0
               x = xc(i)
               y = yc(i)
               if (nbe(i).eq.1) then
                  do j=1,3
                     nc=No_cp(i,j)
!                    ====================================
!                    ==========  SEQUENTIAL =============
#                    ifndef KeyParallel
                     if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                    ====================================
!                    =====  START PARALLEL OPTION =======
#                    else
                     elem = index_global(nc)
                     if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                    endif
!                    =============== END ================
!                    ====================================
                        jj=j+1
                        if (jj.gt.3) jj=jj-3
                        jv1 = No_vp(i,j)
                        jv2 = No_vp(i,jj)
                        taux = xv(jv2)-xv(jv1) 
                        tauy = yv(jv2)-yv(jv1) 
                        ntau = sqrt(taux*taux+tauy*tauy)
                        nxc =  tauy/ntau
                        nyc = -taux/ntau  
                        if (nxc.eq.0) then
                            do k=1,NZ-1
                               phiv(jv1,k)  = 0.0d0
                               phiv(jv2,k)  = 0.0d0
                               Hphiv(jv1,k) = 0.0d0
                               Hphiv(jv2,k) = 0.0d0 
                            enddo
                        endif
                     endif 
                  enddo
               endif
            ENDDO
#        endif

!         ___________________________________
!        |                                   |
!        |    VERTEX VERTICAL BOTTOM (k=1)   |
!        |___________________________________|  

!        ___________________________________
!        NO-SLIP BOTTOM: [u=0,"v=0",w=0] 
#        ifdef KeyBCbottomNoSlip 
            DO nv=1,N_VERT                
               phiv(nv,1)  = 0.0d0  !<<<-- v=0
               Hphiv(nv,1) = 0.0d0
            ENDDO
#        endif
!        ___________________________________
!        EXACT BOTTOM:["u=fB",v=fB,wf=fB] 1st or 2nd order
#        if defined( KeyBCbottomExact_1st) || defined(KeyBCbottomExact_2nd)
            DO nv=1,N_VERT 
               k = 1
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               phiv(nv,k)  = funExamNSv(x,y,z,time)
               Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
            ENDDO                
#        endif                   
!         ___________________________________
!        |                                   |
!        |    VERTEX VERTICAL TOP (k=NZ-1)   |
!        |___________________________________| 

!        ___________________________________
!        NO-SLIP TOP: [u=0,"v=0",w=0] 
#        ifdef KeyBCtopNoSlip 
            DO nv=1,N_VERT                
               phiv(nv,NZ-1)  = 0.0d0  !<<<-- v=0
               Hphiv(nv,NZ-1) = 0.0d0
            ENDDO
#        endif
!        ___________________________________
!        EXACT TOP:["u=fB",v=fB,wf=fB] 1st or 2nd order
#        if defined( KeyBCtopExact_1st) || defined(KeyBCtopExact_2nd)
            DO nv=1,N_VERT 
               k = NZ-1
               x = xv(nv)
               y = yv(nv)                
               z = sigv(k)*Hprv(nv)-hv(nv)  
               phiv(nv,k)  = funExamNSv(x,y,z,time)
               Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
            ENDDO
#        endif  

         DO nv=1,N_VERT
!            ___________________________________
!           |                                   |
!           |      HORIZONTAL WALL (nbev=1)     |
!           |___________________________________|

            if (nbev(nv).eq.1) then
!              ______________________________________________
!              SLIP WALL X-DIRECTION: [du/dn=0,"v=0",dw/dn=0]
#              ifdef KeyBCwallXSlip
                  do k=1,NZ-1
                     phiv(nv,k)  = 0.0d0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              NO-SLIP WALL: [u=0,"v=0",w=0]
#              ifdef KeyBCwallNoSlip
                  do k=1,NZ-1
                     phiv(nv,k)  = 0.0d0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              EXACT WALL: [u=fB,"v=fB",w=fB]
#              ifdef KeyBCwallExact
                  do k=1,NZ-1 
                     x = xv(nv)
                     y = yv(nv)
                     z = sigv(k)*Hprv(nv)-hv(nv)
                     phiv(nv,k)  = funExamNSv(x,y,z,time)
                     Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
                 enddo
#              endif
!              ___________________________________
!              PERIODIC WALL X-DIRECTION
#              ifdef KeyBCperiodicX
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
                  if (TagPeriodicBC(nv,1).eq.1) then
                     nvp = TagPeriodicBC(nv,2)
                     do k=1,NZ-1 
                        phiv(nv,k)  = phiv(nvp,k)
                        Hphiv(nv,k) = Hphiv(nvp,k)
                     enddo
                  endif
!#                 endif
!                 NOTE: In parallel only works if the 
!                       corresponding vertex is in the
!                       same sub-domain.
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL INFLOW (nbev=2)    |
!           |___________________________________|

            elseif (nbev(nv).eq.2) then
!              ___________________________________
!              EXACT INFLOW: [u=fB,"v=0",w=0]
#              ifdef KeyBCinflowExact
                  do k=1,NZ-1 
                     phiv(nv,k)  = 0.0d0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              PERIODIC INFLOW Y-DIRECTION
#              ifdef KeyBCperiodicY
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
                  if (TagPeriodicBC(nv,1).eq.1) then
                     nvp = TagPeriodicBC(nv,2)
                     do k=1,NZ-1 
                        phiv(nv,k)  = phiv(nvp,k)
                        Hphiv(nv,k) = Hphiv(nvp,k)
                     enddo
                  endif
!#                 endif
!                 NOTE: In parallel only works if the 
!                       corresponding vertex is in the
!                       same sub-domain.
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL OUTFLOW (nbev=3)   |
!           |___________________________________|

            elseif (nbev(nv).eq.3) then
!              ___________________________________
!              PERIODIC OUTFLOW Y-DIRECTION
#              ifdef KeyBCperiodicY
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
!                  if (TagPeriodicBC(nv,1).eq.1) then
!                     nvp = TagPeriodicBC(nv,2)
!                     do k=1,NZ-1 
!                        phiv(nv,k)  = phiv(nvp,k)
!                        Hphiv(nv,k) = Hphiv(nvp,k)
!                     enddo
!                  endif
!#                 endif
!                 NOTE: Not necessary
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL STRUCTURE (nbev=4) |
!           |___________________________________|

            elseif (nbev(nv).eq.4) then !<<<-- Structure: u=0
!              ___________________________________
!              NO-SLIP STRUCTURE: [u=0,"v=0",w=0]
               do k=1,NZ-1 
                  phiv(nv,k)  = 0.0d0
                  Hphiv(nv,k) = 0.0d0
               enddo
            endif   
         ENDDO

!*********************************************************************!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!                                                                     !
!               Velocity component  Hw & w (cell-center)              !
!                                                                     !
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!
!*********************************************************************!

      ELSEIF (tagBC.eq.3) THEN
                                 
         DO i=1,N_CELL
            x = xc(i)
            y = yc(i)
!           ___________________________________
!          |                                   |
!          |       VERTICAL BOTTOM (k=1)       |
!          |___________________________________| 

!           ___________________________________                               
!           SLIP or NO-SLIP BOTTOM: "w=0"                
#           if defined(KeyBCbottomSlip) || defined(KeyBCbottomNoSlip)                
               phi(i,1)  = -phi(i,2)  !<<<-- w=0
               Hphi(i,1) = -Hphi(i,2)
#           endif
!           ___________________________________
!           APPROX BOTTOM: [*,*,"w=wfB"]
#           ifdef KeyBCbottomApprox
               phi(i,1)  = wfB(i)    !<<<-- w=wfB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif
!           ___________________________________                               
!           EXACT BOTTOM: [u=fB,v=fB,"w=fB"] 1st order 
#           ifdef KeyBCbottomExact_1st
               z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)           
               fB = funExamNSw(x,y,z,time)
               phi(i,1)  = 2.0d0*fB-phi(i,2) !<<<-- w=fB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif
!           ___________________________________                               
!           EXACT BOTTOM: [u=fB,v=fB,"w=fB"] 2nd order
#           ifdef KeyBCbottomExact_2nd 
               k  = 1
               z  = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)           
               fB = funExamNSw(x,y,z,time) 
               z1 = sig(1)*Hpr(i)-h(i)
               z2 = sig(2)*Hpr(i)-h(i)
               z3 = sig(3)*Hpr(i)-h(i)
               f2 = phi(i,2)
               f3 = phi(i,3)
               a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
               a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
               a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
               f1 = (fB-a2*f2-a3*f3)/a1
               phi(i,1) = f1 !<<<-- w=fB
               Hphi(i,1) = phi(i,1)*Hpr(i)
#           endif            

!           ___________________________________
!          |                                   |
!          |        VERTICAL TOP (k=NZ)        |
!          |___________________________________| 

!           ______________________________________________                              
!           FIXED-SURFACE or NO-SLIP TOP: [*,*,"w=0"]
#           if defined(KeyFixedFreeSurface) || defined(KeyBCtopNoSlip)
               phi(i,NZ) = -phi(i,NZ-1)  !<<<-- w=0
               Hphi(i,NZ)= -Hphi(i,NZ-1)
!           ______________________________________________                              
!           FREE-SURFACE, APPROX TOP:[*,*,"w=wfT"]
#           else
               phi(i,NZ)  = wfT(i)        !<<<-- w=wfT
!              phi(i,NZ) = 2.0d0*wfT(i) - phi(i,NZ-1)
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i)
#           endif
!           ______________________________________________                               
!           EXACT TOP: [u=fB,v=fB,"w=fB"] 1st order 
#           ifdef KeyBCtopExact_1st
               z  = 0.5d0*(sig(NZ)+sig(NZ-1))*Hpr(i)-h(i)            
               fB = funExamNSw(x,y,z,time)
               phi(i,NZ)  = 2.0d0*fB-phi(i,NZ-1) !<<<-- w=fB
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i)
#           endif
!           ______________________________________________                               
!           EXACT TOP: [u=fB,v=fB,"w=fB"] 2nd order
#           ifdef KeyBCtopExact_2nd 
               k  = NZ
               z  = 0.5d0*(sig(NZ)+sig(NZ-1))*Hpr(i)-h(i)            
               fB = funExamNSw(x,y,z,time)
               z1 = sig(NZ-2)*Hpr(i)-h(i)
               z2 = sig(NZ-1)*Hpr(i)-h(i)
               z3 = sig(NZ)*Hpr(i)-h(i)
               f1 = phi(i,NZ-2)
               f2 = phi(i,NZ-1)
               a1 = (z-z2)*(z-z3)/((z1-z2)*(z1-z3))
               a2 = (z-z1)*(z-z3)/((z2-z1)*(z2-z3))
               a3 = (z-z1)*(z-z2)/((z3-z1)*(z3-z2))
               f3 = (fB-a1*f1-a2*f2)/a3
               phi(i,NZ) = f3 !<<<-- w=fB
               Hphi(i,NZ) = phi(i,NZ)*Hpr(i)
#           endif

         ENDDO

         DO i=1,N_CELL0            
!           ___________________________________
!          |                                   |
!          |      HORIZONTAL WALL (nbe=1)      |
!          |___________________________________| 

            if (nbe(i).eq.1) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________                               
!                    DEFAULT SLIP WALL: [*,*,dw/dn=0]
                     do k=2,NZ-1
                        phi(nc,k)  = phi(i,k) !<<-- dw/dn=0
                        Hphi(nc,k) = Hphi(i,k)
                     enddo
!                    ___________________________________                               
!                    NO-SLIP WALL: [u=0,v=0,"w=0"]                                        
#                    ifdef KeyBCwallNoSlip
                        do k=1,NZ
                              phi(nc,k)  = -phi(i,k) !<<-- w=0
                              Hphi(nc,k) = -Hphi(i,k)
                        enddo
#                    endif 
!                    ___________________________________                               
!                    EXACT WALL: [u=fB,v=fB,"w=fB"]                      
#                    ifdef KeyBCwallExact 
                        do k=1,NZ
                           x  = xe(i,j)
                           y  = ye(i,j)
                           z  = sig(k)*Hpr(i)-h(i)                      
                           fB = funExamNSw(x,y,z,time)
                           phi(nc,k)  = 2.0d0*fB-phi(i,k)
                           Hphi(nc,k) = phi(nc,k)*Hpr(i)
                        enddo                           
#                    endif                           
!                    ___________________________________
!                    PERIODIC WALL X-DIRECTION 
#                    ifdef KeyBCperiodicX
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif                                              
!                    ___________________________________
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL INFLOW (nbe=2)     |
!          |___________________________________| 

            elseif (nbe(i).eq.2) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    EXACT INFLOW: [u=fB,v=0,"w=0"]
#                    ifdef KeyBCinflowExact
                        do k=1,NZ 
                           phi(nc,k) = -phi(i,k) !<<<-- w=0                       
                           Hphi(nc,k)= -Hphi(i,k)
                        enddo
#                    endif
!                    ___________________________________
!                    PERIODIC INFLOW Y-DIRECTION
#                    ifdef KeyBCperiodicY
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif
!                    ___________________________________ 
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL OUTFLOW (nbe=3)    |
!          |___________________________________|

            elseif (nbe(i).eq.3) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    DEFAULT NEUMANN OUTFLOW: [du/dn,dv/dn=0,"dw/dn=0"]
#                    ifdef KeyBCoutflowNeumann
                        do k=1,NZ
                           phi(nc,k) = phi(i,k)  !<<<-- dw/dn=0
                           Hphi(nc,k)= Hphi(i,k)
                        enddo
#                    endif
!                    ___________________________________
!                    PERIODIC OUTFLOW Y-DIRECTION
#                    ifdef KeyBCperiodicY
                        i0 = No_cp(nc,3)
!                       ====================================
!                       ==========  SEQUENTIAL =============
#                       ifndef KeyParallel
                           if (i0.ge.1.AND.i0.le.N_CELL0) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
!                       ====================================
!                       =====  START PARALLEL OPTION =======
#                       else
                           elem = index_global(i0)
                           if (elem.ge.1.AND.elem.le.N_CELL0global) then
                              do k=1,NZ
                                 phi(nc,k)  = phi(i0,k)
                                 Hphi(nc,k) = Hphi(i0,k)
                              enddo
                           endif
#                       endif
!                       =============== END ================
!                       ====================================
#                    endif
!                    ___________________________________    
                  endif 
               enddo
!           ___________________________________
!          |                                   |
!          |     HORIZONTAL STRUCTURE (nbe=4)  |
!          |___________________________________|

            elseif (nbe(i).eq.4) then
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
!                    ___________________________________
!                    NO-SLIP STRUCTURE: [u=0,v=0,"w=0"]
                     do k=1,NZ 
                        phi(nc,k) = -phi(i,k) !<<<-- w=0
                        Hphi(nc,k)= -Hphi(i,k)
                     enddo
!                    ___________________________________
                  endif 
               enddo               
            endif
         ENDDO
      
!*********************************************************************!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                     !
!                    Velocity component w (vertex)                    !
!                                                                     !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!*********************************************************************! 

      ELSEIF (tagBC.eq.-3) THEN

!         ___________________________________
!        |                                   |
!        |    VERTEX VERTICAL BOTTOM (k=1)   |
!        |___________________________________| 

!        ______________________________________________
!        SLIP or NO-SLIP BOTTOM: [*,*,"w=0"] 
#        if defined(KeyBCbottomSlip) || defined(KeyBCbottomNoSlip)
            DO nv=1,N_VERT                
               phiv(nv,1)  = 0.0d0  !<<<-- u=0
               Hphiv(nv,1) = 0.0d0
            ENDDO
#        endif
!        ______________________________________________
!        APPROX BOTTOM: [*,*,"w=wfB"]
#        ifdef KeyBCbottomApprox
            DO nv=1,N_VERT                
               phiv(nv,1)  = wfvB(nv)  !<<<-- u=wfvB
               Hphiv(nv,1) = phiv(nv,1)*Hprv(nv)
            ENDDO
#        endif
!        ______________________________________________
!        EXACT BOTTOM:["u=fB",v=fB,wf=fB] 1st or 2nd order
#        if defined( KeyBCbottomExact_1st) || defined(KeyBCbottomExact_2nd)
            DO nv=1,N_VERT 
               k = 1
               z  = sigv(k)*Hprv(nv)-hv(nv)
               x = xv(nv)
               y = yv(nv)         
               phiv(nv,k)  = funExamNSw(x,y,z,time)
               Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
            ENDDO                
#        endif
!         ___________________________________
!        |                                   |
!        |    VERTEX VERTICAL TOP (k=NZ-1)   |
!        |___________________________________| 

!        ______________________________________________                              
!        FIXED-SURFACE or NO-SLIP TOP: [*,*,"w=0"]
#        if defined(KeyFixedFreeSurface) || defined(KeyBCtopNoSlip)
            DO nv=1,N_VERT
               phiv(nv,NZ-1)  = 0.0d0 !<<<-- w=0
               Hphiv(nv,NZ-1) = 0.0d0
            ENDDO
!        ______________________________________________                               
!        FREE-SURFACE, APPROX TOP:[*,*,"w=wfT"]
#        else
            DO nv=1,N_VERT
               phiv(nv,NZ-1)  = wfvT(nv) !<<<-- w=wfvT
               Hphiv(nv,NZ-1) = phiv(nv,NZ-1)*Hprv(nv)
            ENDDO
#        endif
!        ______________________________________________
!        EXACT TOP:["u=fB",v=fB,wf=fB] 1st or 2nd order
#        if defined( KeyBCtopExact_1st) || defined(KeyBCtopExact_2nd)
            DO nv=1,N_VERT 
               k = NZ-1
               z = sigv(k)*Hprv(nv)-hv(nv) 
               x = xv(nv)
               y = yv(nv)           
               phiv(nv,k)  = funExamNSw(x,y,z,time)
               Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
            ENDDO
#        endif

         DO nv=1,N_VERT
!            ___________________________________
!           |                                   |
!           |      HORIZONTAL WALL (nbev=1)     |
!           |___________________________________|

            if (nbev(nv).eq.1) then
!              ___________________________________
!              NO-SLIP WALL: [u=0,v=0,"w=0"]
#              ifdef KeyBCwallNoSlip
                  do k=1,NZ-1
                     phiv(nv,k)  = 0.0d0 !<<<-- w=0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              EXACT WALL: [u=fB,v=fB,"w=fB"]
#              ifdef KeyBCwallExact
                  do k=1,NZ-1 
                     x = xv(nv)
                     y = yv(nv)
                     z = sigv(k)*Hprv(nv)-hv(nv)
                     phiv(nv,k)  = funExamNSw(x,y,z,time)
                     Hphiv(nv,k) = phiv(nv,k)*Hprv(nv)
                  enddo
#              endif
!              ___________________________________
!              PERIODIC WALL X-DIRECTION
#              ifdef KeyBCperiodicX
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
                  if (TagPeriodicBC(nv,1).eq.1) then
                     nvp = TagPeriodicBC(nv,2)
                     do k=1,NZ-1 
                        phiv(nv,k)  = phiv(nvp,k)
                        Hphiv(nv,k) = Hphiv(nvp,k)
                     enddo
                  endif
!#                 endif
!                 NOTE: In parallel only works if the 
!                       corresponding vertex is in the
!                       same sub-domain.
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL INFLOW (nbev=2)    |
!           |___________________________________|

            elseif (nbev(nv).eq.2) then
!              ___________________________________
!              EXACT INFLOW: [u=fB,v=0,"w=0"]
#              ifdef KeyBCinflowExact
                  do k=1,NZ-1 
                     phiv(nv,k)  = 0.0d0 !<<<-- w=0
                     Hphiv(nv,k) = 0.0d0
                  enddo
#              endif
!              ___________________________________
!              PERIODIC INFLOW Y-DIRECTION
#              ifdef KeyBCperiodicY
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
                  if (TagPeriodicBC(nv,1).eq.1) then
                     nvp = TagPeriodicBC(nv,2)
                     do k=1,NZ-1 
                        phiv(nv,k)  = phiv(nvp,k)
                        Hphiv(nv,k) = Hphiv(nvp,k)
                     enddo
                  endif
!#                 endif
!                 NOTE: In parallel only works if the 
!                       corresponding vertex is in the
!                       same sub-domain.
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL OUTFLOW (nbev=3)   |
!           |___________________________________|

            elseif (nbev(nv).eq.3) then
!              ___________________________________
!              PERIODIC OUTFLOW Y-DIRECTION
#              ifdef KeyBCperiodicY
!                 ====================================
!                 ==========  SEQUENTIAL =============
!#                 ifndef KeyParallel
!                  if (TagPeriodicBC(nv,1).eq.1) then
!                     nvp = TagPeriodicBC(nv,2)
!                     do k=1,NZ-1 
!                        phiv(nv,k)  = phiv(nvp,k)
!                        Hphiv(nv,k) = Hphiv(nvp,k)
!                     enddo
!                  endif
!#                 endif
!                 NOTE: Not necessary
!                 =============== END ================
!                 ====================================
#              endif
!            ___________________________________
!           |                                   |
!           |     HORIZONTAL STRUCTURE (nbev=4) |
!           |___________________________________|

            elseif (nbev(nv).eq.4) then 
!              ___________________________________
!              NO-SLIP STRUCTURE: [u=0,v=0,"w=0"]
               do k=1,NZ-1 
                  phiv(nv,k)  = 0.0d0 !<<<-- w=0
                  Hphiv(nv,k) = 0.0d0
               enddo
            endif
         ENDDO

      ENDIF
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCvelocity3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                END OF VELOCITY BOUNDARY CONDITIONS                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
