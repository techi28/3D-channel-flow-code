!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           WATER LEVEL BOUNDARY CONDITION AT CELL-CENTERS            !
!                             Nov 2017                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!   This program calculates the boundary condition of cell-centers    !
!   and vertices of all variables. There are two section to consider  !
!   The vertical boundary and the horizontal boundary.                !
!                                                                     !
!   VERTICAL                                                          !
!   (1) Free surface:                                                 !
!       Wind effects absent, resulting with tangential stress = 0,    !
!                         du/ds = 0, dv/ds =0                         !
!       We ensure that kinematic boundary condition is satisfied:     !
!                  w = deta/dt +u*deta/dx +v*deta/dy                  !
!   (2) Bottom:                                                       !
!       Free-slip boundary conditions are applied:                    !
!                         du/ds = 0, dv/ds =0                         !
!       We ensure that kinematic boundary condition is satisfied:     !
!                  w = -dh/dt -u*dh/dx-v*dh/dy                        !
!                                                                     !
!   HORIZONTAL                                                        !
!   It depends on the type of boundary element. This work uses the    !
!   strategy used by UFVM-ECOMOD and it divides the cells into the    !
!   following types.                                                  !
!          nbe = 0   INSIDE cell (not update)                         !
!          nbe = 1   WALL (slip and no-slip)                          !
!          nbe = 2   DISCHRAGE   (Known values at discharges vertex)  !
!          nbe = 3   WATER LEVEL (Known values at water level vertex) !
!          nbe = 4   FREE BOUNDARY (Radiation method)                 !
!                                                                     !
!---------------------------------------------------------------------!

      SUBROUTINE WaterLevel_BC(Hpr,Hprv,eta,etav,   &
                               h,xc,yc,No_cp,nbe,   &
                               hv,xv,yv,No_vp,nbev)

!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   ______________________________________________________________    !
!  |     Name    |    Size   | Description                        |   !
!  |_____________|___________|____________________________________|   !
!  | <-- Hpr,eta |(N_CELL,NZ)| Pressure at the c-center & boundary|   !
!  |_____________|___________|____________________________________|   !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | ---> xc,yc  |(N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> No_cp  |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe    |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> phiv    |(N_VERT,NZ)| Cell-vertex solution at t(n)        |  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the cell vertices    |  !
!  | --> No_vp   |(N_CELL0,3)| Numbering of the cell vertices      |  !
!  | --> nbe     |(N_CELL0)  | Tag: Type of vertex (inside or bc)  |  !
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

      real*8,dimension(:)   :: Hpr(N_CELL)
      real*8,dimension(:)   :: eta(N_CELL)
      real*8,dimension(:)   :: Hprv(N_VERT)
      real*8,dimension(:)   :: etav(N_VERT)
!     -------------------------------------
      real*8,dimension(:)   :: h(N_CELL)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8,dimension(:)   :: hv(N_VERT)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: con,ddetadx,ddetady,deter,sumfx,sumfy,HprI
      real*8 :: TF,eta_k,term0,term1,term2,DxL,DyL,DISL
      real*8 :: cn,na,n_x,n_y,rr,WLx,WLy,WLnn,WLext
      integer:: jv1,jv2,jv3,jj,jc
      real*8 :: FS_funeta,etaE,WaterLevel,FS_waveBC
      real*8 :: nxc,nyc,nzc,taux,tauy,ntau
      integer :: elem

!      ____________________________________
!     |                                    |
!     |     Declaration of parameters      |
!     |____________________________________|

      integer, parameter :: OptionBCeta = 0
      !  = 0   Reflecting (Neumann)
      !  = 1   Reflecting using vertex with normal ny=0
      !  = 2   Values from vertex
      !  = 3   Exact solution vertex
      !  = 4   Exact solution cell-center

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
      write(*,'(t22,60a)'), '----> Begin subroutine: WaterLevel_BC'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                         CELL-CENTERED VALUES                        !
!                                                                     !
!*********************************************************************!

#     if defined(KeyStaticCylinder) || defined(KeyStaticChannel)
      DO i=1,N_CELL0
         IF (nbe(i).ne.0) THEN
            do j=1,3
               nc=No_cp(i,j)
!              ====================================
!              ==========  SEQUENTIAL =============
#              ifndef KeyParallel
               if (nc.lt.1.OR.nc.gt.N_CELL0) then
!              ====================================
!              =====  START PARALLEL OPTION =======
#              else
               elem = index_global(nc)
               if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!              =============== END ================
!              ====================================
                  eta(nc) = eta(i)
                  Hpr(nc) = eta(nc) + h(nc)
                endif
            enddo
          ENDIF
      ENDDO
      goto 1000
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                  Wall Boundary: nbe = 1                |
!     |________________________________________________________|
     
!     _________________________________________________________
!     CELL-CENTERED POINTS 
    
      DO i=1,N_CELL0
         IF (nbe(i).eq.1) THEN
            do j=1,3
               nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
               if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#              else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!                 =============== END ================
!                 ====================================
                  jj = j+1
                  if (jj.gt.3) jj=jj-3 
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
!                 ______________________________________________
!                 FS case: StandingWave
!                    ------------------------
!                    Reflecting Opt 1
                     if (OptionBCeta.eq.0) then
                        eta(nc) = eta(i)
                        Hpr(nc) = eta(nc) + h(nc)
!                    ------------------------
!                    Reflecting Opt 2
                     elseif (OptionBCeta.eq.1) then
                        taux = xv(jv2)-xv(jv1)
                        tauy = yv(jv2)-yv(jv1)
                        ntau = sqrt(taux*taux+tauy*tauy)
                        nxc =  tauy/ntau
                        nyc = -taux/ntau
                        if (nyc.eq.0) then
                           HprI    = 0.5d0*(Hprv(jv1)+Hprv(jv2))
                           Hpr(nc) = 2.0d0*HprI - Hpr(i)
                           eta(nc) = Hpr(nc) - h(nc)
                        else
                           eta(nc) = eta(i)
                           Hpr(nc) = eta(nc) + h(nc)
                        endif
!                    -------------------------
!                    Values from vertex
                     elseif (OptionBCeta.eq.2) then
                        HprI   = 0.5d0*(Hprv(jv1)+Hprv(jv2))
                        Hpr(nc) = 2.0d0*HprI - Hpr(i)
                        eta(nc) = Hpr(nc) - h(nc) 
!                    -------------------------
!                    Exact solution vertex
                     elseif (OptionBCeta.eq.3) then
                        etaE =   0.5d0*FS_funeta(xv(jv1),yv(jv1),time)&
                               + 0.5d0*FS_funeta(xv(jv2),yv(jv2),time)
                        eta(nc) = 2.0d0*etaE - eta(i)
                        Hpr(nc) = eta(nc) + h(nc)
                        !-------
                        etav(jv1) = FS_funeta(xv(jv1),yv(jv1),time)
                        etav(jv2) = FS_funeta(xv(jv2),yv(jv2),time)
                        Hprv(jv1) = etav(jv1) + hv(jv1)
                        Hprv(jv2) = etav(jv2) + hv(jv2)
!                    -------------------------
!                    Exact solution c-center
                     elseif (OptionBCeta.eq.4) then
                        eta(i)  = FS_funeta(xc(i),yc(i),time)
                        eta(nc) = eta(i)
                        Hpr(nc) = eta(nc) + h(nc)
                     endif    
                endif
            enddo
          ENDIF
      ENDDO
!      ________________________________________________________
!     |                                                        |
!     |               Discharge Boundary: nbe = 2              |
!     |         Minimize wave reflexion by radiation method    |
!     |________________________________________________________|

!1000  continue

      TF = 1.0d9
      DO i=1,N_CELL0
         IF (nbe(i).eq.2) THEN
                !--------------------------
                ! Gradient eta: (d(eta)/dx,d(eta)/dy)
                sumfx = 0.0d0
                sumfy = 0.0d0
                do j=1,3
                   jc = No_cp(i,j)
                   sumfx = sumfx + dxCC(i,j)*(eta(jc)-eta(i))
                   sumfy = sumfy + dyCC(i,j)*(eta(jc)-eta(i))
                enddo
                deter  = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
                ddetadx = (sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
                ddetady = (sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
                !--------------------------
                ! Radiation method
                eta_k=eta(i)
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
                  jj = j+1
                  if (jj.gt.3) jj=jj-3
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
                  DxL = xv(jv2)-xv(jv1)
                  DyL = yv(jv2)-yv(jv1)
                  DISL= dsqrt(DxL**2+DyL**2)
                  term0 = ddetadx*DyL/DISL-ddetady*DxL/DISL
                  term1 = -dt*dsqrt(gra*Hpr(i))*term0
                  term2 = eta(i)-Dt*eta_k/TF
                  eta(nc) = (term1+term2)/(1-Dt/TF)
                  Hpr(nc) = eta(i) + h(i)
                endif
            enddo
         ENDIF
      ENDDO 

!      ________________________________________________________
!     |                                                        |
!     |                Water level boundary: nbe = 3           |
!     |                   (Known values of eta)                |
!     |________________________________________________________|

      IF (N_HB.GE.1) THEN
#        ifdef EllipticalShoal
         WaterLevel = FS_waveBC(time)
#        endif
#        ifdef SolitaryWave
         WaterLevel = 0.0d0 
#        endif
!        _________________________________________________________
!        WATER LEVEL GIVEN AT VERTEX POINTS 
         do nv=1,N_VERT  
            if (nbev(nv).eq.300) then
                etav(nv) = WaterLevel
                Hprv(nv) = etav(nv) + hv(nv)
            endif
         enddo
!        _________________________________________________________
!        CELL-CENTERED POINTS 
         DO i=1,N_CELL0
            IF (nbe(i).eq.3) THEN
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
               if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#              else
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!                  =============== END ================
!                  ====================================
                   jj = j+1
                   if (jj.gt.3) jj=jj-3
                   jv1 = No_vp(i,j)
                   jv2 = No_vp(i,jj)
                   eta(nc) = 0.5d0*(etav(jv1)+etav(jv2))
                   Hpr(nc) = eta(nc) + h(nc)
                endif
            enddo
          ENDIF
          ENDDO  
      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |                 Free Boundary: nbe = 3                 |
!     |         Minimize wave reflexion by radiation method    |
!     |________________________________________________________|

!1000  continue
    
      TF = 1.0d12
      DO i=1,N_CELL0
         IF (nbe(i).eq.300) THEN
                !--------------------------
                ! d(eta)/dx, d(eta)/dy
                sumfx = 0.0d0
                sumfy = 0.0d0
                do j=1,3
                   jc = No_cp(i,j)
                   sumfx = sumfx + dxCC(i,j)*(eta(jc)-eta(i))
                   sumfy = sumfy + dyCC(i,j)*(eta(jc)-eta(i))
                enddo
                deter  = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
                ddetadx = (sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
                ddetady = (sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
                !--------------------------
                ! Radiation method
                eta_k=eta(i)
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
                  jj = j+1
                  if (jj.gt.3) jj=jj-3
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
!                 ---------------
                  DxL = xv(jv2)-xv(jv1)
                  DyL = yv(jv2)-yv(jv1)
                  cn  = dsqrt(gra*Hpr(i))
                  na  = dsqrt(DxL**2+DyL**2)
                  n_x =  DyL/na
                  n_y = -DxL/na
                  rr  = Dt/TF
                  WLx = ddetadx
                  WLy = ddetady
                  WLnn = WLx*n_x+WLy*n_y
                  WLext= h0 + eta(nc)
!                 ---------------              
                  eta(nc) = (eta(nc)-dt*cn*WLnn+rr*WLext)/(1-rr)
                  Hpr(nc) = eta(nc) + h(nc)
                endif
            enddo
         ENDIF
      ENDDO 

!      ________________________________________________________
!     |                                                        |
!     |        HORIZONTAL STRUCTURE (nbe=4)   (Cylinder)       |
!     |________________________________________________________|

!1000  continue

      DO i=1,N_CELL0
         IF (nbe(i).eq.4) THEN
            do j=1,3
               nc=No_cp(i,j)
!              ====================================
!              ==========  SEQUENTIAL =============
#              ifndef KeyParallel
               if (nc.lt.1.OR.nc.gt.N_CELL0) then
!              ====================================
!              =====  START PARALLEL OPTION =======
#              else
               elem = index_global(nc)
               if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#              endif
!              =============== END ================
!              ====================================
                  eta(nc) = -eta(i)
                  Hpr(nc) = eta(nc) + h(nc)
                endif
            enddo
          ENDIF
      ENDDO
      DO nv=1,N_VERT  
         IF (nbev(nv).eq.4) THEN
             etav(nv) = 0.0
             Hprv(nv) = etav(nv) + hv(nv)
         ENDIF
      ENDDO

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

1000  continue

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
      write(*,'(t22,60a)'), '<---- End   subroutine: WaterLevel_BC'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
      
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF WATER ELEVATION                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
