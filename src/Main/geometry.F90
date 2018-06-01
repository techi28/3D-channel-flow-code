!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        GEOMETRY VARIABLES                           !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE calcul_geometry(xc,yc,sig,dsig,No_cp,nbe,h, &
                                 xv,yv,sigv,dsigv,No_vp,nbev,&
                                 ic1tec,ic2tec,ic3tec)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates all the geometric variables correspon-   !
!    ding to the cell center and vertex of each element. More spe-    !
!    cific, we calculate:                                             !
!          1) The coordinate of the sigma direction,                  !
!          2) The cell centers from the vertices,                     !
!          3) The area of each cell,                                  !
!          4) The length of the triangle sides,                       !
!          5) Distance of surronding cell center at each vertex,      !
!          6) Area of surronding cells at each vertex,                !
!          7) The location of the outside (fictitious) cells,         !
!          8) The distance from the cell center to their neighbors.   !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <-- xc      | N_CELL    | x-coordinate of the cell center     |  !
!  | <-- yc      | N_CELL    | y-coordinate of the cell center     |  ! 
!  | <-- sig     | NZ        | sigma                               |  !
!  | <-- dsig    | NZ        | Increment in the sigma-direction    |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Output & input variables:                                        !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <--> h      | N_CELL    | Depth of the domain at cell center  |  !
!  | <--> No_cp  |(N_CELL,3) | Node No. of surrounding three cell  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | --> xv      |(N_VERT)   | x-coordinate of the vertex          |  !
!  | --> yv      |(N_VERT)   | y-coordinate of the vertex          |  !
!  | --> sigv    |(NZ-1)     | sigma of the vertex points          |  !
!  | --> dsigv   |(NZ-1)     | Increment sigma of the vertex points|  !  
!  | --> No_vp   |(N_CELL0,3)| Numbering of cell vertices          |  !
!  | --> nbe     |(N_CELL0)  | Type of tag about the kind of cell  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  !
!  |--N_CELL0    | Number of interior cell centers of the domain   |  !
!  |--N_CELL     | Total number of the cells                       |  !
!  |--N_VERT     | Number of the computing vertices                |  ! 
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !
!  |* N_CELLghost| Number of ghost cell-centers                    |  !
!  |* N_CELLexact| Number of interior + ghost cell-centers         |  !
!  |_____________|_________________________________________________|  !
!  | *  dxVV     |(N_CELL0,3)| dx length of each side of the trian.|  !
!  | *  dyVV     |(N_CELL0,3)| dy length of each side of the trian.|  !
!  | *  dxCC     |(N_CELL0,3)| = dx from cell center to cell neigb.|  !
!  | *  dyCC     |(N_CELL0,3)| = dy from cell center to cell neigb.|  !
!  | *  dlCC     |(N_CELL0,3)| Distance cell-center to the  neigb. |  !
!  |_____________|_________________________________________________|  !
!  | *  xe       |(N_CELL0,3)|  x coordinate of the edge point     |  !
!  | *  ye       |(N_CELL0,3)|  y coordinate of the edge point     |  !
!  | *  dxCE     |(N_CELL0,3)|  =dx from cell-center to the edge   |  !
!  | *  dyCE     |(N_CELL0,3)|  =dy from cell-center to the edge   |  !
!  | *  dlCE     |(N_CELL0,3)|  Distance cell-center to the edge   |  !
!  |_____________|_________________________________________________|  !
!  | *  dlCV     |(N_CELL0,3)| distance from the center to vertex  |  !
!  | *  dlVsum   |(N_VERT)   | Total sum of dlCV distances at vert.|  !
!  | *  areaCELL |(N_CELL)   | Area of each cell                   |  !
!  | *  areaVsum |(N_VERT)   | Total area related to each vertex.  |  !
!  | *  sum_xc2  |(N_CELL0)  | Sum(xc(i)-xc(neighborn))^2          |  !
!  | *  sum_yc2  |(N_CELL0)  | Sum(yc(i)-yc(neighborn))^2          |  !
!  | *  sum_xcyx |(N_CELL0)  | Sum(xc(i)-xc(neig))*(yc(i)-yc(neig))|  !
!  | *  VolPrism |(N_CELL,NZ)| Volume of each prism element        |  !
!  | * dsigPrism |(N_CELL)   | Height of each prism element        |  !
!  |_____________|_________________________________________________|  !
!  | * dhCE      |(N_CELL0,NZglobal,5)  | Perpendicular distance   |  !
!  |             |                      | from cell-center to edge |  !
!  | * areaTriaH |(N_CELL0,3,NZglobal,8)| area of the triangle     |  !
!  | * nxTriaH   |(N_CELL0,3,NZglobal,8)| normal of the triangle x |  !
!  | * nyTriaH   |(N_CELL0,3,NZglobal,8)| normal of the triangle y |  !
!  | * nzTriaH   |(N_CELL0,3,NZglobal,8)| normal of the triangle z |  !
!  | * areaTriaT |(N_CELL0,NZglobal,6)  | area of the top triangles|  !
!  | * nxTriaT   |(N_CELL0,NZglobal,6)  | normal x top triangle    |  !
!  | * nyTriaT   |(N_CELL0,NZglobal,6)  | normal y top triangle    |  !
!  | * nzTriaT   |(N_CELL0,NZglobal,6)  | normal z top triangle    |  !
!  | * areaTriaB |(N_CELL0,NZglobal,6)  | area of the bottom trian.|  !
!  | * nxTriaB   |(N_CELL0,NZglobal,6)  | normal x bottom triangle |  !
!  | * nyTriaB   |(N_CELL0,NZglobal,6)  | normal y bottom triangle |  !
!  | * nzTriaB   |(N_CELL0,NZglobal,6)  | normal z bottom triangle |  !
!  |_____________|_________________________________________________|  !
!  | * i,j,nv,k  |   Loop counters                                 |  !    
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Local variables:                                                 !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | jv1,jv2,jv3 |  Vertex index                                   |  !
!  | jc,ii       |  Cell index                                     |  !
!  | dxCV,dyCV   |  dx & dy distance from the center to the  vertex|  !
!  | XO,YO       |  Intersec. point of vertex and perpendic. line  |  !
!  | dxody       |  dxVV/dyVV  (ratio between trian. side lengths) |  !
!  | dxody2      |  = (dxVV/dyVV)^2                                |  !
!  | sumX,sumY   |  Real variables to add terms                    |  !
!  | sigIni      |  Initial sigma value                            |  !
!  | sigFin      |  Final sigma value                              |  !
!  |_____________|_________________________________________________|  !
!  | xvp1(1:3)   |  Coordinates of the triangle vertex point 1     |  !
!  | xvp2(1:3)   |  Coordinates of the triangle vertex point 2     |  !
!  | xvp3(1:3)   |  Coordinates of the triangle vertex point 3     |  !
!  | a1,a2,a3    |  Displacement vector = xvp2-xvp1                |  !
!  | b1,b2,b3    |  Displacement vector = xvp2-xvp1                |  !
!  | atimesb     |  Cross product = |a x b|                        |  !
!  | areaT3D     |  Area of the triangle in 3D                     |  !
!  | nn1,nn2,nn3 |  Normal coordinates of the triangle in 3D       |  !
!  | xee,yee     |  Perpendicular intersection of each edge        |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
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
!     |   Keys and common parameters                           |
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
!     |    Declaration of variables                            |
!     |________________________________________________________|

      real*8, dimension(:)  :: xc
      real*8, dimension(:)  :: yc  
      real*8, dimension(:)  :: sig
      real*8, dimension(:)  :: dsig 
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe   
      real*8, dimension(:)  :: h  
!     --------------------------------
      real*8, dimension(:)  :: xv
      real*8, dimension(:)  :: yv  
      real*8, dimension(:)  :: sigv
      real*8, dimension(:)  :: dsigv
      integer,dimension(:,:):: No_vp
      integer,dimension(:)  :: nbev   
!     --------------------------------
      integer,dimension(:)  :: ic1tec
      integer,dimension(:)  :: ic2tec
      integer,dimension(:)  :: ic3tec
!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________|

      integer :: jv,jv1,jv2,jv3,jj
      integer :: ic,jc,ii,s,L
      real*8 ::  dxCV,dyCV
      real*8 :: XOj,YOj,dxody,dxody2
      real*8 :: sumX,sumY,sumXY,deter
      real*8, dimension(:) :: xvp1(1:3)
      real*8, dimension(:) :: xvp2(1:3)
      real*8, dimension(:) :: xvp3(1:3)
      real*8 :: a1,a2,a3
      real*8 :: b1,b2,b3
      real*8 :: atimesb
      real*8 :: areaT3D
      real*8 :: nn1,nn2,nn3
      real*8 :: xee,yee
!     ---------------------------------
      real*8, dimension(:) :: SumFun(1:N_VERT)
      real*8 :: taux,tauy,ntau
      real*8 :: norx,nory,nnor
      real*8 :: x,y,z1,z2,z3
      real*8 :: xp1,yp1,xp2,yp2,suma
!     ---------------------------------
      integer:: EdgePointOption
      real*8 :: xm1,xm2,xm3,xm4
      real*8 :: ym1,ym2,ym3,ym4      
!     ---------------------------------
      integer :: initc,lastc,neigc,num,elem
      integer :: jvert,kvert
      integer :: i0,j0,j1,j2,j3,tag,tagNormal
!     ---------------------------------      
      integer :: OptionNormalV
      integer :: ChooseLSMDista,Mix
      real*8  :: val
      integer,dimension(:) :: NumCellperVertex(N_VERT)
      integer,dimension(:) :: tagBC(N_CELL0) 
      integer,dimension(:) :: tagBCv(N_VERT)      
!     ---------------
      integer :: n1,n2,n3,n4,nvp
!     ---------------------------------
      real*8  :: xmin,xmax,dxl,dxl0
      real*8  :: ymin,ymax,dyl,dyl0
      real*8  :: xvjv3,yvjv3
      integer :: vcount,ncount,ncount1,ncount2
      integer :: jvv,jv1p,jv2p
      integer :: jv10,jv20,jv30
      integer :: XPeriodicBC,YPeriodicBC 
      integer :: jvSW,jvNW,jvSE,jvNE     

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: geometry'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                        Vertical domain mesh                         !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |    Homogeneous division of sigma in [0,1] (vertices)   |
!     |________________________________________________________|

      dsigv(1) = (sigFin-sigIni)/float(NZ-2)   
      sigv(1)  = sigIni                 
      do k=2,NZ-2                 
         dsigv(k) = dsigv(1) 
         sigv(k)  = sigv(k-1)+dsigv(k-1) 
      enddo
      dsigv(NZ-1) = dsigv(1)      
      sigv(NZ-1)  = sigFin                 
!      ________________________________________________________
!     |                                                        |
!     |            sigma values at the cell-centers            | 
!     |________________________________________________________|

       dsig(1) = dsigv(1)
       sig(1)  = sigv(1)-0.5d0*dsigv(1)
       do k=2,NZ-1
          dsig(k) = 0.5d0*(dsigv(k)+dsigv(k-1))
          sig(k)  = sig(k-1)+dsig(k-1)
       enddo
       dsig(NZ) = dsigv(NZ-1)
       sig(NZ)  = sig(NZ-1)+dsig(NZ-1)

!*********************************************************************!
!                                                                     !
!           Geometric variables for each cell in the domain           !
!                                                                     !
!*********************************************************************!

      do i=1,N_CELL0
         jv1=No_vp(i,1)
         jv2=No_vp(i,2)
         jv3=No_vp(i,3)
!         ________________________________________________________
!        |                                                        |
!        |                  Cell-centers (xc,yc)                  |
!        |________________________________________________________|
              
          xc(i)=(xv(jv1)+xv(jv2)+xv(jv3))/3.0d0
          yc(i)=(yv(jv1)+yv(jv2)+yv(jv3))/3.0d0
!         ________________________________________________________
!        |                                                        |
!        |                 Area of the triangle                   |
!        |________________________________________________________|

         areaCell(i) = 0.5d0*abs( xv(jv1)*(yv(jv2)-yv(jv3))  &
                                + xv(jv2)*(yv(jv3)-yv(jv1))  &
                                + xv(jv3)*(yv(jv1)-yv(jv2)))
!         ________________________________________________________
!        |                                                        |
!        |              Length of the triangle sides:             |
!        |                    dxVV, dyVV, dlVV,                   |
!        |________________________________________________________|

!        --------------------------------------------------
!        LENGTH OF THE SIDES OF THE TRIANGLES
         dxVV(i,1)=xv(jv2)-xv(jv1)
         dxVV(i,2)=xv(jv3)-xv(jv2)
         dxVV(i,3)=xv(jv1)-xv(jv3)
         dyVV(i,1)=yv(jv2)-yv(jv1)
         dyVV(i,2)=yv(jv3)-yv(jv2)
         dyVV(i,3)=yv(jv1)-yv(jv3)
!        --------------------------------------------------
!        SQUARE LENGTH OF SIDES OF THE TRIANGLES 
         dlVV(i,1) = dsqrt(dxVV(i,1)*dxVV(i,1)+dyVV(i,1)*dyVV(i,1))
         dlVV(i,2) = dsqrt(dxVV(i,2)*dxVV(i,2)+dyVV(i,2)*dyVV(i,2))
         dlVV(i,3) = dsqrt(dxVV(i,3)*dxVV(i,3)+dyVV(i,3)*dyVV(i,3))
      enddo       

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(xc)
         call communication2D(yc)
         call communication2D(areaCell)
#     endif
!     =============== END ================    
!     ====================================

!*********************************************************************!
!                                                                     !
!                      Outside (fictitious) cells                     !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |  Only at boundary cells: wall,discharge & water level  |
!     |________________________________________________________|

      ii = N_CELL0
      DO i=1,N_CELL0
         IF (nbe(i).gt.0) THEN
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
!                 -------------------------------------------
                  ii=ii+1
!                 -------------------------------------------
!                 Index
                  No_cp(i,j)  = ii
!                 -------------------------------------------
!                 Information used in the new BC formulation
                  No_cp(ii,1) = i    !<--- Cell i coming from
                  No_cp(ii,2) = j    !<--- Cell neighborn j of i
!                 -------------------------------------------
!                 Depth
                  h(ii) = h(i)
!                 -------------------------------------------
!                 Area
                  areaCell(ii)= areaCell(i)
!                 -------------------------------------------
!                 Intersection on the triangle boundary edge
                  jj=j+1
                  if (jj.gt.3) jj=jj-3
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
                  if (dabs(dyVV(i,j)).lt.1.0E-7) then
                     yee = yv(jv2)
                     xee = xc(i)
                  else
                     dxody  = dxVV(i,j)/dyVV(i,j)
                     dxody2 = dxody*dxody
                     yee= (yc(i)-(xv(jv2)-xc(i))*dxody &
                               +yv(jv2)*dxody2)/(1.0d0+dxody2)
                     xee= xv(jv2)+(yee-yv(jv2))*dxody
                  endif
!                 -------------------------------------------
!                 Save edge point (xe,ye) at boundary
                  xe(i,j) = xee
                  ye(i,j) = yee
!                 -------------------------------------------
!                 OUTSIDE CELL CENTER (xc,yc)
                  xc(ii)=2*xe(i,j)-xc(i)
                  yc(ii)=2*ye(i,j)-yc(i)
!                 -------------------------------------------
!                 Distance between cell-center and edge point
                  dxCE(i,j) = xe(i,j)-xc(i)
                  dyCE(i,j) = ye(i,j)-yc(i)
!                 -------------------------------------------
!                 Distance between cell-center & the edge line
                  dlCE(i,j) = dsqrt(dxCE(i,j)**2+dyCE(i,j)**2)
!                 -------------------------------------------
               endif
            enddo
         ENDIF
      ENDDO

!      ________________________________________________________
!     |                                                        |
!     |             Exact number of ghost cells                |
!     |________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         N_CELLexact = ii
         N_CELLghost = ii - N_CELL0
#     endif
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                      PERIODIC BC                       |
!     |________________________________________________________|

#     if defined(KeyBCperiodicX)|| defined(KeyBCperiodicY)
!     ___________________________________________________
!     CORNERS
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
                if ((dabs(XDIni-xv(nv)).lt.1.0d-7).and. &
                    (dabs(yDIni-yv(nv)).lt.1.0d-7)) then
                    jvSW = nv
                    !print*,'corner SW:',jvSW,xv(nv),yv(nv)
                endif
                if ((dabs(XDIni-xv(nv)).lt.1.0d-7).and. &
                    (dabs(yDFin-yv(nv)).lt.1.0d-7)) then
                    jvNW = nv
                    !print*,'corner NW:',jvNW,xv(nv),yv(nv)
                endif
                if ((dabs(XDFin-xv(nv)).lt.1.0d-7).and. &
                    (dabs(yDIni-yv(nv)).lt.1.0d-7)) then
                    jvSE = nv
                    !print*,'corner SE:',jvSE,xv(nv),yv(nv)
                endif
                if ((dabs(xDFin-xv(nv)).lt.1.0d-7).and. &
                    (dabs(yDFin-yv(nv)).lt.1.0d-7)) then
                    jvNE = nv
                    !print*,'corner NE',jvNE,xv(nv),yv(nv)
                endif
             endif
         enddo
!     ___________________________________________________
!     TAGS

      TagPeriodicBC = 0
      ii = N_CELL0
      DO i=1,N_CELL0
         IF (nbe(i).gt.0) THEN
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
!                 -------------------------------------------
                  ii=ii+1 
                  No_cp(ii,3) = -1
!                 -------------------------------------------
                  jj=j+1
                  if (jj.gt.3) jj=jj-3
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
!                 -------------------------------------------
!                 Limits of the rectangular domain
                  xmin = XDIni
                  xmax = XDFin
                  ymin = YDIni
                  ymax = YDFin
                  dxl0 = xmax-xmin
                  dyl0 = ymax-ymin
!                 -------------------------------------------
!                 [1] Tags for the direction periodic ghost cells
!                 --------------------------
!                 UpWall-DownWall: X-direction
#                 if defined(KeyBCperiodicX)
                     if (yc(ii).gt.ymax) then
                        dxl =  0.0d0
                        dyl = -dyl0
                        !--------
                        TagPeriodicBC(jv1,1) = 1 !<<-- Only up
                        TagPeriodicBC(jv2,1) = 1 !<<-- Only up
                     elseif(yc(ii).lt.ymin) then
                        dxl = 0.0d0
                        dyl = dyl0
                        !--------
                        !TagPeriodicBC(jv1,1) = 1 !<<-- Only down
                        !TagPeriodicBC(jv2,1) = 1 !<<-- Only down
                     endif
#                 endif
!                 --------------------------
!                 Inflow-Outflow: Y-direction
#                 if defined(KeyBCperiodicY)
                     if (xc(ii).gt.xmax) then
                        dxl = -dxl0
                        dyl =  0.0d0
                     elseif (xc(ii).lt.xmin) then
                        dxl = dxl0
                        dyl = 0.0d0
                        !--------
                        if ((yv(jv1).gt.ymin).and. &
                            (yv(jv1).lt.ymax)) then
                           TagPeriodicBC(jv1,1) = 1 !<<-- Only inflow
                        endif
                        if ((yv(jv2).gt.ymin).and. &
                            (yv(jv2).lt.ymax)) then
                           TagPeriodicBC(jv2,1) = 1 !<<-- Only inflow
                        endif
                     endif
#                 endif
                  IF (dxl*dyl.lt.1.0d-07) THEN
!                 -------------------------------------------
!                 [2.1] Find the vertex 
!                 --------------------------
!                 Index jv1
                  jv1p = 0
                  ncount1 = 0
                  xee = xv(jv1) + dxl
                  yee = yv(jv1) + dyl
                  do nv=1,N_VERT
                     if (nbev(nv).ne.0) then
                         if ((dabs(xee-xv(nv)).lt.1.0d-7).and. &
                             (dabs(yee-yv(nv)).lt.1.0d-7)) then
                             jv1p = nv
                             ncount1 = ncount1 + 1
                         endif
                      endif
                  enddo
!                 --------------------------
!                 Index jv2
                  jv2p = 0
                  ncount2 = 0
                  xee = xv(jv2) + dxl
                  yee = yv(jv2) + dyl
                  do nv =1,N_VERT
                     if (nbev(nv).ne.0) then
                         if ((dabs(xee-xv(nv)).lt.1.0d-7).and. &
                             (dabs(yee-yv(nv)).lt.1.0d-7)) then
                             jv2p = nv
                             ncount2 = ncount2 + 1
                         endif
                      endif
                  enddo
!                 -------------------------------------------
!                 [2.2] Assing vertex tags 
!                 --------------------------
!                 Index jv1
                  if (ncount1.ne.1) then
                     print*,'MESH ERROR! VERTEX PAIR FOR JV1'
                     stop
                  else
                     if (TagPeriodicBC(jv1,2).eq.0) then
                        TagPeriodicBC(jv1,2) = jv1p
                     endif
                  endif
                  !--------
                  ! Corners
                  if (jv1.eq.jvSW) then
                     TagPeriodicBC(jv1,1) = 1
                     TagPeriodicBC(jv1,2) = jvSE
                  endif
                  if (jv1.eq.jvNW) then
                     TagPeriodicBC(jv1,1) = 1
                     TagPeriodicBC(jv1,2) = jvNE
                  endif
!                 --------------------------
!                 Index jv2
                  if (ncount2.ne.1) then
                     print*,'MESH ERROR! VERTEX PAIR FOR JV2'
                     stop
                  else
                     if (TagPeriodicBC(jv2,2).eq.0) then
                        TagPeriodicBC(jv2,2) = jv2p
                     endif
                  endif
                  !--------
                  ! Corners
                  if (jv2.eq.jvSW) then
                     TagPeriodicBC(jv2,1) = 1
                     TagPeriodicBC(jv2,2) = jvSE
                  endif
                  if (jv2.eq.jvNW) then
                     TagPeriodicBC(jv2,1) = 1
                     TagPeriodicBC(jv2,2) = jvNE
                  endif
!                 -------------------------------------------
!                 [3] Find the third index vertex jvv and  
!                 the cell-center index i0 corresponding to ii
                  jvv = 0
                  ncount = 0
                  DO i0=1,N_CELL0
                     IF (nbe(i0).ne.0) THEN
                        jv10 = No_vp(i0,1)
                        jv20 = No_vp(i0,2)
                        jv30 = No_vp(i0,3)
!                       -------------
!                       Case 1
                        if (jv10*jv20 .eq. jv1p*jv2p) then
                            if (jv10 .eq. jv1p) then
                               jvv = jv30
                               No_cp(ii,3) = i0
                               ncount = ncount+1
                            elseif(jv10 .eq. jv2p) then
                               jvv = jv30
                               No_cp(ii,3) = i0
                               ncount = ncount+1
                            endif
                        endif
!                       -------------
!                       Case 2
                        if (jv30*jv20 .eq. jv1p*jv2p) then
                            if (jv30 .eq. jv1p) then
                               jvv = jv10
                               No_cp(ii,3) = i0
                               ncount = ncount+1
                            elseif (jv30 .eq. jv2p) then
                               jvv = jv10
                               No_cp(ii,3) = i0
                               ncount = ncount+1
                            endif
                        endif
!                       -------------
!                       Case 3
                        if (jv10*jv30 .eq. jv1p*jv2p) then
                           if (jv10 .eq. jv1p) then
                               jvv = jv20
                               No_cp(ii,3) = i0
                               ncount = ncount+1
                           elseif (jv10 .eq. jv2p) then
                               jvv = jv20
                               No_cp(ii,3) = i0
                               ncount = ncount+1
                           endif
                        endif
                     ENDIF
                  ENDDO
                  if (ncount.ne.1) then
                      print*,'MESH ERROR! TO FIND PERIDOIC CELL'
                      stop
                  endif
!                 -------------------------------------------
!                 [4.1] New OUTSIDE cell-center values
                  xvjv3 = xv(jvv) - dxl
                  yvjv3 = yv(jvv) - dyl
                  xc(ii) = (xv(jv1)+xv(jv2)+xvjv3)/3.0d0
                  yc(ii) = (yv(jv1)+yv(jv2)+yvjv3)/3.0d0
                  !i0 = No_cp(ii,3)
                  !xc(ii) = xc(i0) -dxl
                  !yc(ii) = yc(i0) -dyl
!                 -------------------------------------------   
                  !if (dyl.eq.0.0d0) write(*,'(a5,i5,i10,f7.2,f7.2,i10,f7.2,f7.2)')'x :',i,ii, &
                  !xc(ii),yc(ii),No_cp(ii,3),xc(No_cp(ii,3)),yc(No_cp(ii,3))
                  !if (dxl.eq.0.0d0) write(*,'(a5,i5,i10,f7.2,f7.2,i10,f7.2,f7.2)')'y :',i,ii, &
                  !xc(ii),yc(ii),No_cp(ii,3),xc(No_cp(ii,3)),yc(No_cp(ii,3))
!                 ------------------------------------------- 
                  ENDIF
               endif
            enddo
         ENDIF
      ENDDO
	
      DO nv=1,N_VERT
         if (TagPeriodicBC(nv,1).eq.1) then
            !nvp = TagPeriodicBC(nv,2)
            !write(*,'(i10,f7.2,f7.2,i10,f7.2,f7.2)') &
            !nv,xv(nv),yv(nv),nvp,xv(nvp),yv(nvp)
         endif
      ENDDO

#     endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(xc)
         call communication2D(yc)
#     endif
!     =============== END ================    
!     ====================================
            

!*********************************************************************!
!                                                                     !
!                New tags for the Static Cylinder Case                !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Tags of the static cylinder case:          |
!     |                      nbe=1: Wall                       |
!     |                      nbe=2: inflow                     |
!     |                      nbe=3: outflow                    | 
!     |                      nbe=4: structure                  |
!     |________________________________________________________|

#     ifdef KeyStaticCylinder       
!        _____________________________________________
!        Cell-Centers
         tagBC = 0
         DO ii = N_CELL0+1,N_CELL
            i = No_cp(ii,1)
            !---------------
            ! Walls
            if (yc(ii).lt.YDIni) then
                nbe(i) = 1
                tagBC(i) = 1
            endif
            if (yc(ii).gt.YDFin) then
                nbe(i) = 1
                tagBC(i) = 1
            endif            
            !---------------
            ! Inflow                    
            !if ((xc(ii).lt.XDIni).and.(tagBC(i).eq.0)) then
            if (xc(ii).lt.XDIni) then            
                nbe(i) = 2
            endif
            !---------------
            ! Outflow                
            if ((xc(ii).gt.XDFIn).and.(tagBC(i).eq.0)) then
                nbe(i) = 3
            endif 
            !---------------
            ! Structure
            if ((yc(ii).gt.YDIni).AND.(yc(ii).lt.YDFin).AND. &
                (xc(ii).gt.XDIni).AND.(xc(ii).lt.XDFIn)) then 
                nbe(i) = 4
            endif                      
         ENDDO
!        _____________________________________________
!        Vertex
         n1 = 0
         n2 = 0
         n3 = 0
         n4 = 0
         tagBCv = 0
         DO nv=1,N_VERT
            if (nbev(nv).ne.0) then
                !---------------
                ! Walls
                if (dabs(yv(nv)-YDIni).lt.1.0d-7) then
                    n1 = n1 + 1
                    !No_wb(n1) = nv
                    nbev(nv) = 1
                    tagBCv(nv) = 1
                endif
                if (dabs(yv(nv)-YDFin).lt.1.0d-7) then
                    n1 = n1 + 1
                    !No_wb(n1) = nv
                    nbev(nv) = 1
                    tagBCv(nv) = 1
                endif            
                !---------------            
                ! Inflow             
                if (xv(nv).le.(XDIni+1.0d-7).and. &
                    (tagBCv(nv).eq.0)) then
                    n2 = n2 + 1
                    !No_qb(n2) = nv
                    nbev(nv) = 2
                endif
                !---------------
                ! Outflow
                if (xv(nv).ge.(XDFin-1.0d-7).and. &
                    (tagBCv(nv).eq.0)) then
                    n3 = n3 + 1
                    !No_hb(n3) = nv
                    nbev(nv) = 3   
                endif
                !---------------
                ! Structure
                if ((yv(nv).gt.YDIni).AND.(yv(nv).lt.YDFin).AND. &
                    (xv(nv).gt.XDIni).AND.(xv(nv).lt.XDFin)) then
                    n4 = n4 + 1
                    nbev(nv) = 4
                endif                
            endif
         ENDDO
         N_WB = n1
         N_QB = n2
         N_HB = n3
         N_SB = n4
#     endif

!*********************************************************************!
!                                                                     !
!             Variables for the least square technique                !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |        dxCC, dyCC, dlCC, sum_xc, sum_yc & sum_xcyc     |
!     |________________________________________________________|

      do i=1,N_CELL0
         sumX = 0.0d0
         sumY = 0.0d0
         sumXY= 0.0d0
         do j=1,3
            jc = No_cp(i,j)
            dxCC(i,j) = xc(jc)-xc(i) 
            dyCC(i,j) = yc(jc)-yc(i)
            dlCC(i,j) = sqrt(dxCC(i,j)**2+dyCC(i,j)**2)
            sumX = sumX + dxCC(i,j)*dxCC(i,j)
            sumY = sumY + dyCC(i,j)*dyCC(i,j)
            sumXY= sumXY+ dxCC(i,j)*dyCC(i,j)
         enddo
         sum_xc2(i) = sumX
         sum_yc2(i) = sumY
         sum_xcyc(i)= sumXY
       enddo

!      ________________________________________________________
!     |                                                        |
!     |      df/dx= aGx0*f(i) + aGx1*f1 + aGx2*f2 + aGx3*f3    |
!     |      df/dy= aGy0*f(i) + aGy1*f1 + aGy2*f2 + aGy3*f3    |
!     |________________________________________________________|

      do i=1,N_CELL0
         do j=1,3
            jc = No_cp(i,j)
            jv = No_vp(i,j)
            dxCV = xv(jv)-xc(i) 
            dyCV = yv(jv)-yc(i) 
            deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
            aGx(i,j) = (sum_yc2(i)*dxCC(i,j)-sum_xcyc(i)*dyCC(i,j))/deter
            aGy(i,j) = (sum_xc2(i)*dyCC(i,j)-sum_xcyc(i)*dxCC(i,j))/deter
         enddo
         aGx(i,0) = -(aGx(i,1)+aGx(i,2)+aGx(i,3))
         aGy(i,0) = -(aGy(i,1)+aGy(i,2)+aGy(i,3))
       enddo
      
             
!*********************************************************************!
!                                                                     !
!             Variables used for the vertex interpolation             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |    Distance of surronding cell center at each vertex   |
!     |             New way to collect weight values           |
!     |                 surrounding, weight                    |
!     |________________________________________________________|

#     ifdef KeyInterpoNew
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         DO nv=1,N_VERT
            j0 = 0
            dlVsum(nv) = 0.0d0
!           ___________________
!           Inside cells
            do i=1,N_CELL0
               do j=1,3
                  if (nv.eq.No_vp(i,j)) then
                     j0 = j0 + 1
                     surrounding(nv,j0) = i
                     dxCV = xc(i)-xv(nv)
                     dyCV = yc(i)-yv(nv)
                     weight(nv,j0) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
                     dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                  endif
               enddo
            enddo
!           ___________________
!           Ghost cells
#           ifdef KeyUseInterGhost
            if (nbev(nv).ne.0) then
                do ii=N_CELL0+1,N_CELL0+N_CELLghost
                   nc = No_cp(ii,1)
                   do j=1,3
                      if (nv.eq.No_vp(nc,j)) then
                          j0 = j0 + 1
                          surrounding(nv,j0) = ii
                          dxCV = xc(ii)-xv(nv)
                          dyCV = yc(ii)-yv(nv)
                          weight(nv,j0) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
                          dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                      endif
                   enddo 
               enddo
            endif
#           endif
!           ___________________
!           Number of surroun
            Dimsurrounding(nv) = j0
         ENDDO
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         DO nv=1,N_VERT
            j0 = 0
            dlVsum(nv) = 0.0d0
!           ___________________
!           Inside cells
            do i=1,N_CELL0
               do j=1,3
                  if (nv.eq.No_vp(i,j)) then
                     j0 = j0 + 1
                     surrounding(nv,j0) = i
                     dxCV = xc(i)-xv(nv)
                     dyCV = yc(i)-yv(nv)
                     weight(nv,j0) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
                     dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                  endif
               enddo
            enddo
!           ___________________
!           Ghost cells
#           ifdef KeyUseInterGhost
            if (nbev(nv).ne.0) then
                do ii=N_CELL0+1,N_CELL0+N_CELLghost
                   nc = No_cp(ii,1)
                   do jj=1,3
                      if (nv.eq.No_vp(nc,jj)) then
                          j0 = j0 + 1
                          surrounding(nv,j0) = ii
                          dxCV = xc(ii)-xv(nv)
                          dyCV = yc(ii)-yv(nv)
                          weight(nv,j0) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
                          dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                      endif
                   enddo 
               enddo
            endif
#           endif
!           __________________
!           Overlaping cells
            kvert = index_globalv(nv)  
            do i = N_CELL0+N_CELLghost+1,N_CELL
               elem = index_global(i)
               do j=1,3
                  jvert = No_vp_global(elem,j)
                  if (jvert.eq.kvert) then
                     j0 = j0 + 1
                     surrounding(nv,j0) = i
                     dxCV = xc(i)-xv(nv)
                     dyCV = yc(i)-yv(nv)
                     weight(nv,j0) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
                     dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                  endif
               enddo
            enddo
!           ___________________
!           Number of surroun
            Dimsurrounding(nv) = j0
         ENDDO
#        endif
!        =============== END ================    
!        ====================================
#     else
!      ________________________________________________________
!     |                                                        |
!     |    Distance of surronding cell center at each vertex   |
!     |                    dxCV, dyCV, dlCV                    |
!     |________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
!        __________________
!        Initial values 
         do i=1,N_CELL
            do j=1,3
               dlCV(i,j) = 0.0d0
            enddo
         enddo
         do nv=1,N_VERT
            dlVsum(nv)=0.0d0
         enddo
!        __________________
!        Inside cells 
         do i=1,N_CELL0
            do j=1,3
               nv = No_vp(i,j)
               dxCV = xc(i)-xv(nv)
               dyCV = yc(i)-yv(nv)
               dlCV(i,j) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
               dlVsum(nv)= dlCV(i,j) + dlVsum(nv)
            enddo
         enddo
!        __________________
!        Ghost cells
#        ifdef KeyUseInterGhost
         do i=1,N_CELL0
            if (nbe(i).ne.0) then
               do j=1,3
                  nc = No_cp(i,j)
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
                     j1 = j
                     j2 = j+1
                     if (j2.gt.3) j2=j2-3
!                    -----------------
                     jv1 = No_vp(i,j1)
                     jv2 = No_vp(i,j2)
                     dlVsum(jv1) = dlVsum(jv1) + dlCV(i,j1)
                     dlVsum(jv2) = dlVsum(jv2) + dlCV(i,j2)
                  endif
               enddo 
            endif
         enddo
#        endif
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
!        __________________
!        Initial values 
         do i=1,N_CELL
            do j=1,3
               dlCV(i,j) = 0.0d0
            enddo
         enddo
         do nv=1,N_VERT
            dlVsum(nv)=0.0d0
         enddo
!        __________________
!        Inside cells 
         do i=1,N_CELL0
            do j=1,3
               nv = No_vp(i,j)
               dxCV = xc(i)-xv(nv)
               dyCV = yc(i)-yv(nv)
               dlCV(i,j) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
               dlVsum(nv)=  dlVsum(nv) + dlCV(i,j)
            enddo
         enddo
!        __________________
!        Ghost cells
#        ifdef KeyUseInterGhost
         do i=1,N_CELL0
            if (nbe(i).ne.0) then	
               do j=1,3
                  nc = No_cp(i,j)
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
                     jj=j+1
                     if (jj.gt.3) jj=jj-3
                     jv1 = No_vp(i,j)
                     jv2 = No_vp(i,jj)
                     dlVsum(jv1) = dlVsum(jv1) + dlCV(i,j)
                     dlVsum(jv2) = dlVsum(jv2) + dlCV(i,jj)
                  endif
               enddo 
            endif
         enddo
#        endif
!        __________________
!        Overlaping cells
         do i = N_CELL0+N_CELLghost+1,N_CELL
            elem = index_global(i)
            do j=1,3
               jvert = No_vp_global(elem,j)
               tag = 0
               do nv=1,N_VERT
                  kvert = index_globalv(nv)  
                  if (jvert.eq.kvert) then
                     dxCV = xc(i)-xv(nv)
                     dyCV = yc(i)-yv(nv)
                     dlCV(i,j) = 1.0d0/sqrt(dxCV*dxCV+dyCV*dyCV)
                     dlVsum(nv)= dlVsum(nv) + dlCV(i,j) 
                     tag = 1
                  endif
               enddo
               if (tag.eq.0) then
                  dlCV(i,j) = 0.0d0
               endif
            enddo
         enddo
#     endif
!     =============== END ================    
!     ====================================
#     endif

!      ________________________________________________________
!     |                                                        |
!     |        Area of surronding cells at each vertex         |
!     |________________________________________________________|

      do nv=1,N_VERT
         areaVsum(nv)= 0.0d0
      enddo

      do i=1,N_CELL0
         do j=1,3
            nv = No_vp(i,j)
            areaVsum(nv) = areaVsum(nv) + areaCell(i) 
         enddo        
      enddo

!      ________________________________________________________
!     |                                                        |
!     |      LSM of surronding cell center at each vertex      |
!     |             New way to collect weight values           |
!     |               surroundingLSM, weightLSM                |
!     |________________________________________________________|
      
#     ifdef KeyInterpoNewPressureLSM

      ChooseLSMDista = 1
      Mix = 0 ! =2 for MixBC (eliminate corners)
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel

!         --------------------------------
!         Number of cells per each vertex
          NumCellperVertex = 0
          do i=1,N_CELL0
             do j=1,3
                nv = No_vp(i,j)
                NumCellperVertex(nv) =  NumCellperVertex(nv) + 1
             enddo
          enddo
          
         DO nv=1,N_VERT
            j0 = 0
            dlVsumLSM(nv) = 0.0d0
!           ___________________
!           Inside cells
            do i=1,N_CELL0
               do j=1,3
                  if (nv.eq.No_vp(i,j)) then
!                    --------------                   
                     dxCV = xv(nv)-xc(i)
                     dyCV = yv(nv)-yc(i)
                     dlCV(i,j) = 1.0d0/sqrt(dxCV**2+dyCV**2)
!                    --------------                                 
                     if (ChooseLSMDista.eq.1) then
                        val = dlCV(i,j)
                     else
                        val = 1.0d0
                     endif            
                     dlVsumLSM(nv)= dlVsumLSM(nv) + val                                  
!                    -------------- corners
                     if (NumCellperVertex(nv).le.Mix) then
                         surroundingLSM(nv,j0+1) = i
                         weightLSM(nv,j0+1)      = val
                         j0 = j0 + 1
!                    --------------                         
                     else
                        surroundingLSM(nv,j0+1) = i                     
                        surroundingLSM(nv,j0+2) = No_cp(i,1)
                        surroundingLSM(nv,j0+3) = No_cp(i,2)
                        surroundingLSM(nv,j0+4) = No_cp(i,3)                      
                        weightLSM(nv,j0+1)=val*(aGx(i,0)*dxCV+aGy(i,0)*dyCV+1.0d0)
                        weightLSM(nv,j0+2)=val*(aGx(i,1)*dxCV+aGy(i,1)*dyCV)
                        weightLSM(nv,j0+3)=val*(aGx(i,2)*dxCV+aGy(i,2)*dyCV)             
                        weightLSM(nv,j0+4)=val*(aGx(i,3)*dxCV+aGy(i,3)*dyCV)
                        j0 = j0 + 4
                     endif
!                    --------------                         
                  endif
               enddo
            enddo
!           ___________________
!           Number of surroun
            DimsurroundingLSM(nv) = j0            
         ENDDO
!        ___________________
!        Final weights            
         DO nv=1,N_VERT
            do j=1,DimsurroundingLSM(nv)
               weightLSM(nv,j) = weightLSM(nv,j)/dlVsumLSM(nv)
            enddo
         ENDDO 
                    
#     endif
!     =============== END ================    
!     ====================================
#     endif                  

!*********************************************************************!
!                                                                     !
!    outward normal unit vector at each cell-center and vertex point  !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                    Normal cell-center                  |
!     |________________________________________________________|

      do i=1,N_CELL0
         do j=1,3
            normxc(i,j) = 0.0d0
            normyc(i,j) = 0.0d0
         enddo
      enddo
      do i=1,N_CELL0
         if (nbe(i).ne.0) then
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
                  jj=j+1
                  if (jj.gt.3) jj=jj-3
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
                  taux = xv(jv2)-xv(jv1)
                  tauy = yv(jv2)-yv(jv1)
                  ntau = sqrt(taux*taux+tauy*tauy)
                  taux = taux/ntau
                  tauy = tauy/ntau
                  normxc(i,j) =  tauy
                  normyc(i,j) = -taux
               endif
            enddo
         endif
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                      Normal vertex                     |
!     |________________________________________________________|
 
      OptionNormalV = 1
      
      normxv = 0.0d0
      normyv = 0.0d0
      SumFun = 0.0d0
      do i=1,N_CELL0
         if (nbe(i).ne.0) then	
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
                  jj=j+1
                  if (jj.gt.3) jj=jj-3
                  jv1 = No_vp(i,j)
                  jv2 = No_vp(i,jj)
!                 -------------------------------------
!                 Option 1 : Vertex same as cell-center
                  if (OptionNormalV.eq.1) then
                     normxv(jv1) =  normxc(i,j)
                     normyv(jv1) =  normyc(i,j)
                     normxv(jv2) =  normxc(i,j)
                     normyv(jv2) =  normyc(i,j)  
                  endif                 
!                 -------------------------------------
!                 Option 2, Part 1: Vertex average of cell-centers
                  if (OptionNormalV.eq.2) then
                     normxv(jv1) =  normxv(jv1) + normxc(i,j)
                     normyv(jv1) =  normyv(jv1) + normyc(i,j)
                     SumFun(jv1) =  SumFun(jv1) + 1.0d0
                     normxv(jv2) =  normxv(jv2) + normxc(i,j)
                     normyv(jv2) =  normyv(jv2) + normyc(i,j)
                     SumFun(jv2) =  SumFun(jv2) + 1.0d0
                  endif
               endif
            enddo
         endif
      enddo
!     -------------------------------------
!     Option 2, Part2: Vertex average of cell-centers 
      if (OptionNormalV.eq.2) then     
         do nv=1,N_VERT
            if (SumFun(nv).eq.0) then
               normxv(nv) = 0.0d0
               normyv(nv) = 0.0d0            
            else
               norx = normxv(nv)/SumFun(nv)
               nory = normyv(nv)/SumFun(nv)
               nnor = dsqrt(norx*norx+nory*nory)
               normxv(nv) = norx/nnor
               normyv(nv) = nory/nnor
            endif
         enddo
      endif

!*********************************************************************!
!                                                                     !
!           Normal points xn1, xn2, and their element location        !
!                                                                     !
!*********************************************************************!

      tagNormal= 0
#     ifdef KeyNeumannVertex
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
!      ________________________________________________________
!     |                                                        |
!     |               Distance dn of normal points             |
!     |________________________________________________________|

      dn = 0.0d0
      do i=1,N_CELL0
         do j=1,3
            nv=No_vp(i,j)
            if (nbev(nv).ne.0) then  
               !dn(nv) = 0.5d0*dlCV(i,j) 
               dn(nv) = 0.45d0*min(dlVV(i,1),dlVV(i,2),dlVV(i,3)) 
            endif
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |               Points and triangle location             |
!     |________________________________________________________|

      do nv=1,N_VERT
         xn1(nv) = xv(nv)
         yn1(nv) = yv(nv)
         xn2(nv) = xv(nv)
         yn2(nv) = yv(nv)
         icxn1(nv) = 200000!nv (auxiliar to start)
         icxn2(nv) = 100000!nv (auxiliar to start)
      enddo
      do nv=1,N_VERT
         if (nbev(nv).ne.0) then  
!           ______________________________________________________
!           First point: xn1 in the normal direction
            xn1(nv) = xv(nv)+dn(nv)*(-normxv(nv))
            yn1(nv) = yv(nv)+dn(nv)*(-normyv(nv))
            xn2(nv) = xv(nv)+2.0d0*dn(nv)*(-normxv(nv))
            yn2(nv) = yv(nv)+2.0d0*dn(nv)*(-normyv(nv))
!           -------------------------------------       
!           Localizing the element for point xn1
            x = xn1(nv)
            y = yn1(nv)
            s=0
300         continue
            s=s+1
            jv1 = No_vp(s,1)
            jv2 = No_vp(s,2)
            jv3 = No_vp(s,3)
            z1 = dxVV(s,1)*(y-yv(jv1))-dyVV(s,1)*(x-xv(jv1))
            z2 = dxVV(s,2)*(y-yv(jv2))-dyVV(s,2)*(x-xv(jv2))
            z3 = dxVV(s,3)*(y-yv(jv3))-dyVV(s,3)*(x-xv(jv3))
!           -----------------------------------
!           Conditional to be inside
            if ((z1.ge.0).and.(z2.ge.0).and.(z3.ge.0)) then
               icxn1(nv) = s
               goto 400
            elseif (s.gt.N_CELL0) then
               print*, ' '
               print*,'Geometry Error!: No element for normal point 1'
               print*,'nv=', nv,'xv=',xv(nv),'yv(nv)',yv(nv)
            else
               goto 300
            endif
400         continue
!           ______________________________________________________
!           Second point: xn2 in the normal direction
            xn2(nv) = xv(nv)+2.0d0*dn(nv)*(-normxv(nv))
            yn2(nv) = yv(nv)+2.0d0*dn(nv)*(-normyv(nv))
!           -----------------------------------
!           Localizing the element for point xn2
            x = xn2(nv)
            y = yn2(nv)
            s=0
500         continue
            s=s+1
            jv1 = No_vp(s,1)
            jv2 = No_vp(s,2)
            jv3 = No_vp(s,3)
            z1 = dxVV(s,1)*(y-yv(jv1))-dyVV(s,1)*(x-xv(jv1))
            z2 = dxVV(s,2)*(y-yv(jv2))-dyVV(s,2)*(x-xv(jv2))
            z3 = dxVV(s,3)*(y-yv(jv3))-dyVV(s,3)*(x-xv(jv3))
!           -----------------------------------
!           Conditional to be inside
            if ((z1.ge.0).and.(z2.ge.0).and.(z3.ge.0)) then
               icxn2(nv) = s
               goto 600
            elseif (s.gt.N_CELL0) then
               print*,'Geometry Error!: No element for normal point 2'
               print*,'nv=', nv,'xv=',xv(nv),'yv(nv)',yv(nv)
            else
               goto 500
            endif
600         continue
         endif
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         tagNormal= 0
         do nv=1,N_VERT
            if (nbev(nv).eq.2) then
               tagNormal = 1
            endif
         enddo
         if (tagNormal.eq.1) then  
            print*,'   '
            print*,'    ----------------------------------------------------'
            print*,'     WARNING!!! :Normal cells to vertex points not in   '
            print*,'                 parallel yet. This is only used in     '
            print*,'                 the case of Neumman BC. Theese values  '
            print*,'                 are related to the normal.             '
            print*,'    ----------------------------------------------------'
            print*,'                                                        '
         endif
#     endif
!     =============== END ================    
!     ====================================
#     endif

!*********************************************************************!
!                                                                     !
!        3 Edge points (xme(i,j),yme(i,j)) of each element i          !
!        & distance between the cell-centers and this edge point      !
!                                                                     !
!*********************************************************************!

      EdgePointOption=1

!     ________________________________________________________________
!     Option 1: Middle point between vertices (Original 2017)
      IF (EdgePointOption.eq.1) THEN
         do i=1,N_CELL0
            do j=1,3
               jj=j+1
               if (jj.gt.3) jj=jj-3
               jv1 = No_vp(i,j)
               jv2 = No_vp(i,jj)
               xme(i,j) = 0.5d0*(xv(jv2)+xv(jv1)) 
               yme(i,j) = 0.5d0*(yv(jv2)+yv(jv1))
            enddo
         enddo      
!     ________________________________________________________________
!     Option 2:  Intersection point between the line connecting
!                cell centers and the edge 
      ELSEIF (EdgePointOption.eq.2) THEN
         do i=1,N_CELL0
            do j=1,3
               jj=j+1
               if (jj.gt.3) jj=jj-3
               jv1 = No_vp(i,j)
               jv2 = No_vp(i,jj)
               jc  = No_cp(i,j)
               xm1 = xc(i)
               xm2 = xc(jc)
               xm3 = xv(jv1)
               xm4 = xv(jv2)
               ym1 = yc(i)
               ym2 = yc(jc)
               ym3 = yv(jv1)
               ym4 = yv(jv2)
               xme(i,j) = ((xm1*ym2-ym1*xm2)*(xm3-xm4) &
                         -(xm1-xm2)*(xm3*ym4-ym3*xm4))/ &
                          ((xm1-xm2)*(ym3-ym4)-(ym1-ym2)* &
                          (xm3-xm4))
               yme(i,j) = ((xm1*ym2-ym1*xm2)*(ym3-ym4) &
                         -(ym1-ym2)*(xm3*ym4-ym3*xm4))/ &
                          ((xm1-xm2)*(ym3-ym4)-(ym1-ym2)* &
                          (xm3-xm4))
            enddo
         enddo       
!     ________________________________________________________________
!     Option 3: Middle point between cell-centers 
      ELSEIF (EdgePointOption.eq.3) THEN
         do i=1,N_CELL0
            do j=1,3  
               jc = No_cp(i,j)
!              -------------------------------------------
!              Location of the point
               xme(i,j) = 0.5d0*(xc(i)+xc(jc)) 
               yme(i,j) = 0.5d0*(yc(i)+yc(jc))
            enddo
         enddo
!     ________________________________________________________________
!     Option 4: Perpendicular point to the edge  
      ELSEIF (EdgePointOption.eq.4) THEN
         do i=1,N_CELL0
            do j=1,3
!              -------------------------------------------
!              Location of the point	
               jj=j+1
               if (jj.gt.3) jj=jj-3
               jv1 = No_vp(i,j)
               jv2 = No_vp(i,jj)
               if (dabs(dyVV(i,j)).lt.1.0E-7) then
                  yee = yv(jv2)
                  xee = xc(i)
               else
                  dxody  = dxVV(i,j)/dyVV(i,j)
                  dxody2 = dxody*dxody
                  yee= (yc(i)-(xv(jv2)-xc(i))*dxody &
                       +yv(jv2)*dxody2)/(1.0d0+dxody2)
                  xee= xv(jv2)+(yee-yv(jv2))*dxody
               endif
               xme(i,j) = xee
               yme(i,j) = yee
            enddo
         enddo
!     ________________________________________________________________
!     Option 5: Average of the two triangle centroids around the edge 
      ELSEIF (EdgePointOption.eq.5) THEN
         do i=1,N_CELL0
            do j=1,3
               nc=No_cp(i,j)
               jj=j+1
               if (jj.gt.3) jj=jj-3
               jv1 = No_vp(i,j)
               jv2 = No_vp(i,jj)
!              -------------------------------------------
!              Location of the point
               xp1 = (xv(jv2)+xv(jv1)+xc(i))/3.0d0
               yp1 = (yv(jv2)+yv(jv1)+yc(i))/3.0d0
               xp2 = (xv(jv2)+xv(jv1)+xc(nc))/3.0d0
               yp2 = (yv(jv2)+yv(jv1)+yc(nc))/3.0d0
!              -------------------------------------------
!              Location of the point
               xme(i,j) = 0.5d0*(xp1+xp2) 
               yme(i,j) = 0.5d0*(yp1+yp2)
            enddo
         enddo
      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |   Distance between face center and intersection point  |
!     |________________________________________________________|

         do i=1,N_CELL0
            do j=1,3
               jj=j+1
               if (jj.gt.3) jj=jj-3
               jv1 = No_vp(i,j)
               jv2 = No_vp(i,jj)
!              ---------------------------------------
!              Intersection point between the line connecting
!              cell centers and the edge
               jc = No_cp(i,j)
               xm1 = xc(i)
               xm2 = xc(jc)
               xm3 = xv(jv1)
               xm4 = xv(jv2)
               ym1 = yc(i)
               ym2 = yc(jc)
               ym3 = yv(jv1)
               ym4 = yv(jv2)
               xmo(i,j) = ((xm1*ym2-ym1*xm2)*(xm3-xm4) &
                         -(xm1-xm2)*(xm3*ym4-ym3*xm4))/ &
                          ((xm1-xm2)*(ym3-ym4)-(ym1-ym2)* &
                          (xm3-xm4))
               ymo(i,j) = ((xm1*ym2-ym1*xm2)*(ym3-ym4) &
                         -(ym1-ym2)*(xm3*ym4-ym3*xm4))/ &
                          ((xm1-xm2)*(ym3-ym4)-(ym1-ym2)* &
                          (xm3-xm4))
!              --------------------------------------
!              Get distance between face center and intersection point
               xm1 = xme(i,j) - xmo(i,j)
               ym1 = yme(i,j) - ymo(i,j)
               doe(i,j) = dsqrt(xm1*xm1+ym1*ym1)
!              ------------------------------------
!              Determine the direction with e2
               xm2 = xm1*dxVV(i,j)+ym1*dyVV(i,j)
               if (xm2 .lt. 0.0d0) then
                  doe(i,j) = -doe(i,j)
               endif
            enddo
         enddo
               
!      ________________________________________________________
!     |                                                        |
!     |            Spetial treatment of boundary points        |
!     |________________________________________________________|
      
      do i=1,N_CELL0
         if (nbe(i).gt.0) then
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
                  !xme(i,j) = xe(i,j)
                  !yme(i,j) = ye(i,j)
               endif
            enddo
         endif
      enddo

!*********************************************************************!
!                                                                     !
!           Geometric variables for each prism element                !
!                                                                     !
!*********************************************************************!

      VolPrism = 0.0d0
      do i=1,N_CELL 
         VolPrism(i,1) = areaCell(i)*dsigv(1)
         do k=2,NZ
            VolPrism(i,k) = areaCell(i)*dsigv(k-1) 
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3Dtype2(VolPrism)
#     endif	
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!           Height of the piramid for each face of the prism          !
!                                                                     !
!*********************************************************************!

         do i=1,N_CELL0
 	    do j=1,3
!              ----------------------------------------------------
!              Vertex index	
               jj=j+1
               if (jj.gt.3) jj=jj-3
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)
!              ----------------------------------------------------
!              Intersection on the triangle boundary edge  i
	       if (dabs(dyVV(i,j)).lt.1.0E-7) then
	          yee = yv(jv2)
	          xee = xc(i)				
               else
	          dxody  = dxVV(i,j)/dyVV(i,j)
	          dxody2 = dxody*dxody
	          yee= (yc(i)-(xv(jv2)-xc(i))*dxody &
                              +yv(jv2)*dxody2)/(1.0d0+dxody2)
	          xee= xv(jv2)+(yee-yv(jv2))*dxody
	       endif
!              ----------------------------------------------------
!              Perpendicular distance from cell-center to the edge
	       dhCE(i,j) = dsqrt((xee-xc(i))**2+(yee-yc(i))**2)
!              ----------------------------------------------------
!              Intersection on the triangle boundary edge  j
               jc = No_cp(i,j)
               if (dabs(dyVV(i,j)).lt.1.0E-7) then
	          yee = yv(jv2)
	          xee = xc(jc)				
               else
	          dxody  = dxVV(i,j)/dyVV(i,j)
	          dxody2 = dxody*dxody
	          yee = (yc(jc)-(xv(jv2)-xc(jc))*dxody &
                                +yv(jv2)*dxody2)/(1.0d0+dxody2)
	          xee = xv(jv2)+(yee-yv(jv2))*dxody
	       endif
!              ----------------------------------------------------
!              Perpendicular distance from cell-center to the edge
	       dhCE(i,j) = dsqrt((xee-xc(jc))**2+(yee-yc(jc))**2)+dhCE(i,j)
            enddo
         enddo              

!*********************************************************************!
!                                                                     !
!             Index of tecplot cell-centers extra elements            !
!                                                                     !
!*********************************************************************!

      num = 0
      DO nv=1,N_VERT
!        ---------------------------------------------------------
!        Find the initial cell-center and its neighborn
         IF (nbev(nv).ne.0) THEN
            do i=1,N_CELL0
               if (nbe(i).ne.0) then
                  do j=1,3
                     jc = No_cp(i,j)
                     jv = No_vp(i,j)
                     if (jc.gt.N_CELL0) then
                        jj=j+1
                        if (jj.gt.3) jj=jj-3
	                jv1 = No_vp(i,j)
	                jv2 = No_vp(i,jj)
                        if (jv1.eq.nv.or.jv2.eq.nv) then
                           initc = jc
                           lastc = jc  
                           neigc = i
                           goto 100
                        endif 
                     endif
                  enddo
               endif
            enddo
         ELSE
            do i=1,N_CELL0
 	       do j=1,3
                  jc = No_cp(i,j)
                  jv = No_vp(i,j)
                  if (jv.eq.nv) then
                     initc = i
                     lastc = i
                     neigc = jc
                     goto 100
                  endif
               enddo
            enddo
         ENDIF
        
!        ---------------------------------------------------------
!        The three index of the tecplot element	
         100 continue
         IF ((neigc.ge.1).and.(neigc.le.N_CELL0)) THEN               
!           ------------------------------
!           First 2 index 
            num = num + 1 
            ic1tec(num) = initc
            ic2tec(num) = neigc      
!           ------------------------------
!           Finding the third index
 	    do j=1,3
               jc = No_cp(neigc,j)
               if (jc.ne.lastc) then
                  jj=j+1
                  if (jj.gt.3) jj=jj-3
	          jv1 = No_vp(neigc,j)
	          jv2 = No_vp(neigc,jj)
                  if ((jv1.eq.nv).or.(jv2.eq.nv)) then
                     ic3tec(num) = jc    
                     lastc = neigc         
                     neigc = jc     
                     if (jc.eq.initc) then
                        num=num-1
                        goto 200
                     else
                        goto 100    
                     endif
                  endif
               endif
            enddo           
         ENDIF
         200 continue         
      ENDDO
      Numictec = num
 
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: geometry'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	      END OF GEOMETRY                         !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
