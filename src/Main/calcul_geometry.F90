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
      real*8 ::  dxCV,dyCV,dlCV
      real*8 :: XOj,YOj,dxody,dxody2
      real*8 :: sumX,sumY,sumXY
      real*8, dimension(:) :: xvp1(1:3)
      real*8, dimension(:) :: xvp2(1:3)
      real*8, dimension(:) :: xvp3(1:3)
      real*8 :: a1,a2,a3
      real*8 :: b1,b2,b3
      real*8 :: atimesb
      real*8 :: areaT3D
      real*8 :: nn1,nn2,nn3
      real*8 :: xee,yee
      real*8 :: nnx,nny,tnx,tny
!     ---------------------------------
      real*8 :: taux,tauy,ntau
      real*8 :: norx,nory,nnor
      real*8 :: x,y,z1,z2,z3
      real*8 :: xp1,yp1,xp2,yp2,suma
!     ---------------------------------
      integer :: EdgePointOption
!     ---------------------------------
      integer :: initc,lastc,neigc,num,elem
      integer :: jvert,kvert
      integer :: i0,j0,j1,j2,j3,tag,tagNormal
      real*8 :: som1,som2,som3,som4,SOMsom
!     ---------------------------------
      real*8 :: xm1,xm2,xm3,xm4
      real*8 :: ym1,ym2,ym3,ym4
!     ---------------------------------
      real*8 :: xmin, xmax, dxl, det1
      real*8 :: ymin, ymax, dyl, det2
      integer :: kk, cc,s0,vcount,ncount
      integer :: jv1p,jv2p,jv10,jv20,jv30
      integer :: ncount1,ncount2,nvg
      real*8  :: nfn1,Volsum
!     ---------------------------------
#     ifdef KeyParallel
      integer :: Iini,Ifin,Ielem,Jelem,Kelem,ghost_global,ccount
      integer :: cloop,localc,globalc,vcc,count1,count2
#     endif
      integer:: irec
      character*50 filen
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
     Volsum = 0.
!    -------------------
      do i=1,N_CELL0
        jv1=No_vp(i,1)
        jv2=No_vp(i,2)
        jv3=No_vp(i,3)
!      ________________________________________________________
!     |                                                        |
!     |                  Cell-centers (xc,yc)                  |
!     |________________________________________________________|

        xc(i)=(xv(jv1)+xv(jv2)+xv(jv3))/3.0d0
        yc(i)=(yv(jv1)+yv(jv2)+yv(jv3))/3.0d0
!      ________________________________________________________
!     |                                                        |
!     |                 Area of the triangle                   |
!     |________________________________________________________|

        areaCell(i) = 0.5d0*dabs( xv(jv1)*(yv(jv2)-yv(jv3))  &
            + xv(jv2)*(yv(jv3)-yv(jv1))  &
            + xv(jv3)*(yv(jv1)-yv(jv2)))
!      ________________________________________________________
!     |                                                        |
!     |              Length of the triangle sides:             |
!     |                    dxVV, dyVV, dlVV,                   |
!     |________________________________________________________|

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
!        --------------------------------------------------
!        outward normal for automatic mesh
        xnn(i,1) = dyVV(i,1)/dlVV(i,1)
        xnn(i,2) = dyVV(i,2)/dlVV(i,2)
        xnn(i,3) = dyVV(i,3)/dlVV(i,3)
        ynn(i,1) = -dxVV(i,1)/dlVV(i,1)
        ynn(i,2) = -dxVV(i,2)/dlVV(i,2)
        ynn(i,3) = -dxVV(i,3)/dlVV(i,3)
!       ----------------------------------------------------
!        SUM THE VOLUME OF ALL DOMAIN
        Volsum = Volsum + areaCell(i)
!        ---------------------------------------------------
!        outward normal for bluekenue mesh(seems not needed)
!        nn1 = xv(jv3) -xv(jv1)
!        nn2 = yv(jv3)- yv(jv1)
!        nn3 = nn1*dyVV(i,1) + nn2*(-dxVV(i,1))
!        if(nn3 .gt. 0.) then
!            xnn(i,1) = -1.0d0*xnn(i,1)
!            ynn(i,1) = -1.0d0*ynn(i,1)
!            print*,'normal reversed for cell',i,'edge',1
!        endif
!        nn1 = xv(jv1) -xv(jv2)
!        nn2 = yv(jv1)- yv(jv2)
!        nn3 = nn1*dyVV(i,2) + nn2*(-dxVV(i,2))
!        if(nn3 .gt. 0.) then
!            xnn(i,2) = -1.0d0*xnn(i,2)
!            ynn(i,2) = -1.0d0*ynn(i,2)
!            print*,'normal reversed for cell',i,'edge',2
!        endif
!        nn1 = xv(jv2) -xv(jv3)
!        nn2 = yv(jv2)- yv(jv3)
!        nn3 = nn1*dyVV(i,3) + nn2*(-dxVV(i,3))
!        if(nn3 .gt. 0.) then
!            xnn(i,3) = -1.0d0*xnn(i,3)
!            ynn(i,3) = -1.0d0*ynn(i,3)
!            print*,'normal reversed for cell',i,'edge',3
!        endif
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication2D(xc)
      call communication2D(yc)
      call communication2D(areaCell)
      call SUM_parallel(Volsum,TotalVolume)
      TotalVolume = TotalVolume/NZBlock
#     else
      TotalVolume = Volsum
#     endif
!     =============== END ================
!     ====================================
!      ________________________________________________________
!     |                                                        |
!     |            Check Periodic BC implementation            |
!     |________________________________________________________|
!#     ifdef KeyTESTpBC
!       irec=60
!       filen='../output/Parallel/pb_  .txt'
!       write(filen(23:24),'(i2.2)') rang_topo+1
!       open(irec,file=filen)
!#     endif
!*********************************************************************!
!                                                                     !
!                      Outside (fictitious) cells                     !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |  Only at boundary cells: wall,discharge & water level  |
!     |________________________________________________________|

#     ifndef KeyParallel
         ii = N_CELL0
         vcount = 0
      DO i=1,N_CELL0
         IF (nbe(i).gt.0) THEN
           do j=1,3
             nc=No_cp(i,j)
             if (nc.lt.1.OR.nc.gt.N_CELL0) then
                  ii=ii+1
!                 -------------------------------------------
!                 Index
                  No_cp(i,j)  = ii
!                 -------------------------------------------
!                 Information used in the new BC formulation
                  No_cp(ii,1) = i    !<--- Cell i coming from
                  No_cp(ii,2) = j    !<--- Cell neighborn j of i
!                 -------------------------------------------
!                 Area
                  areaCell(ii)= areaCell(i)
!                 -------------------------------------------
!                 Depth
                  h(ii) = h(i)
!                 -------------------------------------------
!                 Perpendicur point on the triangle boundary edge
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
                    xc(ii)=2.0d0*xee-xc(i)
                    yc(ii)=2.0d0*yee-yc(i)
                    No_cp(ii,3) = -1 ! not periodic point
!          -------------------------------------------
!                      Periodic BC Implementation
#         ifdef KeyTESTpBC
                       xmin = XDIni
                       xmax = XDFin
                       ymin = YDIni
                       ymax = YDFin
                       dxl = xmax - xmin
                       dyl = ymax - ymin
                       jv1p = 0
                       jv2p = 0
                       cc = 0
!           --------------------------------------------
                    IF(XPB .EQ. 1) THEN
                         if (xc(ii).gt. xmax) then
                                dxl = -1.0d0*dxl
                                dyl = 0.
                         elseif(xc(ii) .lt. xmin) then
                                dxl = dxl
                                dyl = 0.
                         endif
                    ENDIF
!          ----------------------------------------------
                    IF(YPB .EQ. 1) THEN
                         if (yc(ii).gt. ymax) then
                                dxl = 0.
                                dyl = -1.0d0*dyl
                         elseif(yc(ii) .lt. ymin) then
                                dxl = 0.
                                dyl = dyl
                         endif
                    ENDIF
!          ================================================
                    IF(dxl*dyl .gt. 1.0E-7) THEN
                        vcount = vcount + 1
                    ELSE
!           -----------------------------------------------
!            Pair for jv1
                        ncount = 0
                        xee = xv(jv1) + dxl
                        yee = yv(jv1) + dyl
                        do nv =1,N_VERT
                            if(nbev(nv) .ne. 0) then
                                if(dabs(xee-xv(nv)) .lt. 1.0E-7) then
                                    if(dabs(yee-yv(nv)) .lt. 1.0E-7) then
                                        jv1p = nv
                                        ncount = ncount +1
                                    endif
                                endif
                            endif
                        enddo

                        if(ncount .ne. 1) print*, 'ERROR!VERTEX PAIR FOR JV1'
!           -----------------------------------------------
!            Pair for jv2
                                ncount = 0
                                xee = xv(jv2) + dxl
                                yee = yv(jv2) + dyl
                                do nv =1,N_VERT
                                    if(nbev(nv) .ne. 0) then
                                        if(dabs(xee-xv(nv)) .lt. 1.0E-7) then
                                            if(dabs(yee-yv(nv)) .lt. 1.0E-7) then
                                                jv2p = nv
                                                ncount = ncount+1
                                            endif
                                        endif
                                    endif
                                enddo

                        if(ncount .ne. 1) print*, 'ERROR!VERTEX PAIR FOR JV2'

!              print*, jv1,jv1p,jv2,jv2p
!              print*, xv(jv1),yv(jv1),xv(jv1p),yv(jv1p)
!              print*, xv(jv2),yv(jv2),xv(jv2p),yv(jv2p)
!           -----------------------------------------------
                                ncount = 0
                                do i0 =1, N_CELL0
                                    if(nbe(i0) .ne. 0) then

                                        jv10 = No_vp(i0,1)
                                        jv20 = No_vp(i0,2)
                                        jv30 = No_vp(i0,3)

                                        if(jv10*jv20 .eq. jv1p*jv2p) then
                                            if(jv10 .eq. jv1p) then
                                                cc = jv30
                                                No_cp(ii,3) = i0
                                                ncount = ncount+1
                                            elseif(jv10 .eq. jv2p) then
                                                cc = jv30
                                                No_cp(ii,3) = i0
                                                ncount = ncount+1
                                            endif
                                        endif

                                        if(jv30*jv20 .eq. jv1p*jv2p) then
                                            if(jv30 .eq. jv1p) then
                                                cc = jv10
                                                No_cp(ii,3) = i0
                                                ncount = ncount+1
                                            elseif(jv30 .eq. jv2p) then
                                                cc = jv10
                                                No_cp(ii,3) = i0
                                                ncount = ncount+1
                                            endif
                                        endif

                                        if(jv10*jv30 .eq. jv1p*jv2p) then
                                            if(jv10 .eq. jv1p) then
                                                cc = jv20
                                                No_cp(ii,3) = i0
                                                ncount = ncount+1
                                            elseif(jv10 .eq. jv2p) then
                                                cc = jv20
                                                No_cp(ii,3) = i0
                                                ncount = ncount+1
                                            endif
                                        endif
                                    endif
                                enddo
!                   ------------------------------------------------------
                        if (ncount .ne. 1) print*, 'ERROR FIND PERIDOIC CELL'
                        xee = xv(cc) - dxl
                        yee = yv(cc) - dyl
                        xc(ii)=(xv(jv1)+xv(jv2)+xee)/3.0d0
                        yc(ii)=(yv(jv1)+yv(jv2)+yee)/3.0d0
                        det1 = dabs(xc(ii)-xc(No_cp(ii,3)))/2.0d0
                        det2 = dabs(yc(ii)-yc(No_cp(ii,3)))/2.0d0
           ENDIF
#          endif
                endif
            enddo
        ENDIF
        ENDDO
            print*, '      No. of ghost cells not periodic =', vcount
#    endif
!    ________________________________________________________
!    PARALLEL OPTION (WITH PERIODIC BC)
#     ifdef KeyParallel
!     -------------------
          s = rang_topo + 1

#     ifndef KeyTESTpBC
          ii = N_CELL0
#     endif
            count1 = 0
            count2 = 0
            s0 = 0
            DO i =1, N_CELL0
                Ielem = Index_global(i)
                if(nbe(i) .ne. 0) then
                    do j =1, 3
                        Jelem = No_cp(i,j) ! = 0 boundary cell
                        Kelem = No_cp_global(Ielem,j) ! global neighbour index
                    if((Jelem .lt.1) .OR. (Kelem .gt. N_CELL0global)) then
!            ----------------------------------------------------
!            Determine the local index of ghost point
#                 ifdef KeyTESTpBC
                    Iini = DimghostCCdom(s,2)
                    Ifin = DimghostCCdom(s,3)
                    i0 = 0
                    ii = 0
                    ghost_global = 0
                    do kk = Iini, Ifin
                        Jelem = ghostCCdom(kk,1)
                        Kelem = CCdom(Jelem,1)
                        cc = ghostCCdom(kk,2)
                        if(Ielem .eq. Kelem) then
                            if(cc .eq. j) then !why essential?
                                s0 = s0 +1
                                i0 = i0 +1
                                ii = ghost_local(kk) + N_CELL0 !----> Crucial
                                ghost_global = kk
                            endif
                        endif
                    enddo
!              =====================================================================
!                Check if local index found!
                 if (i0 .ne. 1) then
                    kk = Iini + s0
                    Jelem = ghostCCdom(ghost_global,1)
                    Kelem = CCdom(Jelem,1)
                    cc = ghostCCdom(ghost_global,2)
                   print*,'Fatal Error in Calculate Geometry',i0,s,s0,Iini,Ifin,Ielem,Kelem,cc,j
                   stop
                 endif
!                =================================================================
#                 else
                  ii = ii +1
#                 endif
!                 -------------------------------------------
!                 Index
                  No_cp(i,j)  = ii
!                 -------------------------------------------
!                 Information used in the new BC formulation
                  No_cp(ii,1) = i    !<--- Cell i coming from
                  No_cp(ii,2) = j    !<--- Cell neighborn j of i
!                 -------------------------------------------
!                 Area
                  areaCell(ii)= areaCell(i)
!                 -------------------------------------------
!                 Depth
                  h(ii) = h(i)
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
                    xc(ii)=2.0d0*xee-xc(i)
                    yc(ii)=2.0d0*yee-yc(i)
!                  -------------------------------------------
                    No_cp(ii,3) = -1 ! not periodic point
!                  -------------------------------------------
#           ifdef KeyTESTpBC
                     cc = ghostCCdom(ghost_global,3)

                     if(CCdom(cc,1) .ne. Ielem) then
                            xmin = XDIni
                            xmax = XDFin
                            ymin = YDIni
                            ymax = YDFin
                            dxl = xmax - xmin
                            dyl = ymax - ymin
                            jv1p = 0
                            jv2p = 0
!           --------------------------------------------
                     IF(XPB .EQ. 1) THEN
                         if (xc(ii).gt. xmax) then
                                dxl = -1.0d0*dxl
                                dyl = 0.
                         elseif(xc(ii) .lt. xmin) then
                                dxl = dxl
                                dyl = 0.
                         endif
                     ENDIF
!          ----------------------------------------------
                    IF(YPB .EQ. 1) THEN
                         if (yc(ii) .gt. ymax) then
                                dxl = 0.
                                dyl = -1.0d0*dyl
                         elseif(yc(ii) .lt. ymin) then
                                dxl = 0.
                                dyl = dyl
                         endif
                    ENDIF

                    IF(dxl*dyl .gt. 1.0E-7) THEN
                       print*, 'ERROR IN PARALLEL MODE. CALCULATE_GEOMETRY'
                       stop
                    ELSE
!           -----------------------------------------------
!            Pair for jv1
                        ncount = 0
                        xee = xv(jv1) + dxl
                        yee = yv(jv1) + dyl
                        do nv =1,N_VERTglobal
                            if(dabs(xee-xv_global(nv)) .lt. 1.0E-7) then
                                if(dabs(yee-yv_global(nv)) .lt. 1.0E-7) then
                                    jv1p = nv
                                    ncount = ncount +1
                                endif
                            endif
                        enddo

                        if(ncount .ne. 1) print*, 'ERROR!VERTEX PAIR FOR JV1'

!           -----------------------------------------------
!            Pair for jv2
                        ncount = 0
                        xee = xv(jv2) + dxl
                        yee = yv(jv2) + dyl
                        do nv =1,N_VERTglobal
                            if(dabs(xee-xv_global(nv)) .lt. 1.0E-7) then
                                if(dabs(yee-yv_global(nv)) .lt. 1.0E-7) then
                                    jv2p = nv
                                    ncount = ncount+1
                                endif
                            endif
                        enddo

                        if(ncount .ne. 1) print*, 'ERROR!VERTEX PAIR FOR JV2'
!           -----------------------------------------------
                        ncount = 0
                        do i0 =1, N_CELL0global
                            if(nbe_global(i0) .ne. 0) then
                                jv10 = No_vp_global(i0,1)
                                jv20 = No_vp_global(i0,2)
                                jv30 = No_vp_global(i0,3)

                                if(jv10*jv20 .eq. jv1p*jv2p) then
                                    if(jv10 .eq. jv1p) then
                                        vcc = jv30
                                        globalc = i0
                                        ncount = ncount+1
                                    elseif(jv10 .eq. jv2p) then
                                        vcc = jv30
                                        globalc = i0
                                        ncount = ncount+1
                                    endif
                                endif

                                if(jv30*jv20 .eq. jv1p*jv2p) then
                                    if(jv30 .eq. jv1p) then
                                        vcc = jv10
                                        globalc = i0
                                        ncount = ncount+1
                                    elseif(jv30 .eq. jv2p) then
                                        vcc = jv10
                                        globalc = i0
                                        ncount = ncount+1
                                    endif
                                endif

                                if(jv10*jv30 .eq. jv1p*jv2p) then
                                    if(jv10 .eq. jv1p) then
                                        vcc = jv20
                                        globalc = i0
                                        ncount = ncount+1
                                    elseif(jv10 .eq. jv2p) then
                                        vcc = jv20
                                        globalc = i0
                                        ncount = ncount+1
                                    endif
                                endif
                            endif
                        enddo
!                   ------------------------------------------------------
                       if (ncount .ne. 1) print*, 'ERROR FIND PERIDOIC CELL',ncount
!                   ------------------------------------------------------
                        xee = xv_global(vcc) - dxl
                        yee = yv_global(vcc) - dyl
                        xc(ii)=(xv(jv1)+xv(jv2)+xee)/3.0d0
                        yc(ii)=(yv(jv1)+yv(jv2)+yee)/3.0d0

                        areaCell(ii) = 0.5d0*dabs(xv(jv1)*(yv(jv2)-yee)  &
                                      + xv(jv2)*(yee-yv(jv1))  &
                                      + xee*(yv(jv1)-yv(jv2)))
                      ENDIF

                        ccount = ghostCCdom(ghost_global,5)

                        if(ccount .ne. CCdom(cc,2)) then
                            print*, ' assumption is wrong'
                            stop
                        endif

                        if (ccount .ne. s) then
                            No_cp(ii,3) = 0
                            count1 = count1 +1
!                            write(irec,*) 'remote pb',ii,No_cp(ii,3),xc(ii),yc(ii),globalc,xc_global(globalc),yc_global(globalc)
                        else
                            localc = 0
                            do cloop=1,N_CELL0
                                if(Index_global(cloop) .eq. globalc) then
                                    localc = cloop
                                    count2 = count2 +1
                                endif
                            enddo
                            if(localc .eq. 0) print*, '   assumption is wrong!'
                            No_cp(ii,3) = localc
!                            write(irec,*) 'local pb',ii,No_cp(ii,3),xc(ii),yc(ii),localc,xc(localc),yc(localc)
                        endif

                    endif
#           endif
!     =============== END Periodic ================
!     =============================================
                endif !-----> End ghost cell
            enddo !-----> End j loop
        endif
        ENDDO !  ----->  End N_CELL0 loop
           print*,'pb on proc#: ', rang_topo, 'remote: ', count1, 'local: ', count2
!            write(irec,*) 'total local pb',  count2
!            write(irec,*) 'total remote pb', count1
#     endif
!     =============== END PARALLEL ================
!     =============================================

!#     ifdef KeyDbgPBC
!       close(irec)
!#     endif

!      ________________________________________________________
!     |                                                        |
!     |             Exact number of ghost cells                |
!     |________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         N_CELLexact = ii
         N_CELLghost = ii - N_CELL0
      print*, '      N_CELLghost =',N_CELLghost
      print*, '      N_CELLexact =',N_CELLexact
#     endif
!     =============== END ================
!     ====================================
      if(ChooseBoundary .eq. 1) then
!*********************************************************************!
!                                                                     !
!                       Assign TagXC,TagYC info                       !
!                                                                     !
!*********************************************************************!
        DO i=1,N_CELL
            TagXC(i) = 0
            TagYC(i) = 0
            xee = xc(i)
            yee = yc(i)
            ncount = 0
            if(xee .lt. XDIni) then
                tagXC(i) = 1
                ncount = ncount +1
            endif
            if(xee .gt. XDFIn) then
                tagXC(i) = 2
                ncount = ncount +1
            endif
            if(yee .lt. YDIni) then
                tagYC(i) = 1
                ncount = ncount +1
            endif
            if(yee .gt. YDFin) then
                tagYC(i) = 2
                ncount = ncount +1
            endif
            if(ncount .gt. 1) then
                print*, '      ERROR IN TAGXC,TAGYC. EXIT'
                stop
            endif
        ENDDO
!      ------------------------
#      ifdef KeyDbg
            print*, '       TagXC and TagYC initialized for all cells!'
#      endif
!*********************************************************************!
!                                                                     !
!                       Assign TagXV,TagYV info                       !
!                                                                     !
!*********************************************************************!
!      --------------------------
!      check local vertex
        do nv =1, N_VERT
            tagXV(nv) = 0
            tagYV(nv) = 0
            if (nbev(nv).ne.0) then
                xee = xv(nv)
                yee = yv(nv)
                ncount = 0
                if(dabs(xee - XDIni) .lt. 1.0d-7) then
                    tagXV(nv) = 1
                endif

                if(dabs(xee - XDFin) .lt. 1.0d-7) then
                    tagXV(nv) = 2
                endif

                if(dabs(yee - YDIni) .lt. 1.0d-7) then
                    tagYV(nv) = 1
                endif

                if(dabs(yee - YDFin) .lt. 1.0d-7) then
                    tagYV(nv) = 2
                endif
            endif
        enddo

#        ifndef KeyParallel
#        ifdef KeyTESTpBC
!     ----------------------------------------------------------
!          allocate variable for periodic vertex in serial mode
            allocate (vertex_pair(N_VERT),vertex_edge(N_VERT))
            do nv=1,N_VERT
                vertex_pair(nv) = -1
                vertex_edge(nv) = -1
            enddo
!     ----------------------------------------------------------
!         check vertex and boundary points
          do nv=1,N_VERT
                xee = xv(nv)
                yee = yv(nv)

                if(TagXV(nv)*TagYV(nv) .gt. 0) then
                    ncount = 2
                elseif(TagXV(nv) .gt. 0) then
                    ncount = 1
                elseif(TagYV(nv) .gt. 0) then
                    ncount = 1
                else
                    ncount = 0
                endif
!           ____________________________________
!            Boundary Point
            if(ncount .eq. 1) then
                vcount = 0
                do i=1,N_VERT
                    if (nbev(i).ne.0) then
!                    _________________________________________________
!                       periodic pair in y direction
                     if (YPB .eq. 1) then
                        if(dabs(xee - xv(i)) .lt. 1.0d-7) then
                            if(dabs(yee - yv(i)) .gt. (YDFin-YDIni-1.0d-7)) then
                             vcount = vcount +1
                             vertex_pair(nv) = i
                             endif
                        endif
                     endif
!                   _____________________________________________________
!                       periodic pair in x direction
                     if (XPB .eq. 1) then
                         if(dabs(yee - yv(i)) .lt. 1.0d-7) then
                            if(dabs(xee - xv(i)) .gt. (XDFin-XDIni-1.0d-7)) then
                             vcount = vcount +1
                             vertex_pair(nv) = i
                             endif
                         endif
                      endif

                    endif
                enddo
!            ---------------------------------------
!              Check the vertex pair
                if (vcount .eq. 0) then
                   vertex_pair(nv) = -1
                elseif (vcount .ne. 1) then
                   print*, 'ERROR! MISSING PERIODIC NB'
                endif
            endif ! end ncount = 1
!           _______________________________________
!           Edge Point
            if(ncount .eq. 2) then
                if(XPB*YPB .eq. 1) then
                  vertex_edge(nv) = 1
                else
                  vcount = 0
                  vertex_edge(nv) = 0
!               ------------------------------------------------------
                  do i=1,N_VERT
                    if (nbev(i).ne.0) then
!                    _________________________________________________
!                       periodic pair in y direction
                     if (YPB .eq. 1) then
                        if(dabs(xee - xv(i)) .lt. 1.0d-7) then
                            if(abs(yee - yv(i)) .gt. (YDFin-YDIni-1.0d-7)) then
                             vcount = vcount +1
                             vertex_pair(nv) = i
                             endif
                        endif
                     endif
!                   _____________________________________________________
!                       periodic pair in x direction
                     if (XPB .eq. 1) then
                         if(dabs(yee - yv(i)) .lt. 1.0d-7) then
                            if(dabs(xee - xv(i)) .gt. (XDFin-XDIni-1.0d-7)) then
                             vcount = vcount +1
                             vertex_pair(nv) = i
                             endif
                         endif
                      endif
                    endif
                  enddo
!            ---------------------------------------
!              Check the vertex pair
                  if (vcount .ne. 1) then
                     print*, 'ERROR! MISSING PERIODIC NB FOR CORNER'
                  endif

               endif ! end XPB*YPB = 1

               print*, '       CORNER VERTEX',nv, 'CORNER:',vertex_edge(nv),'PAIR:',vertex_pair(nv)
            endif ! end ncount = 2

            enddo

#      ifdef KeyDbg
        ncount = 0
        vcount = 0
        do nv=1,N_VERT
            if(vertex_edge(nv) .gt. 0) ncount = ncount + 1
            if(vertex_pair(nv) .gt. 0) vcount = vcount + 1
        enddo
        print*, '       TagXV and TagYV initialized for all cells!'
        print*, '       CORNER VERTEX',ncount,'PERIODIC VERTEX',vcount
        stop
#      endif

#      endif
#      endif
!*********************************************************************!
!                                                                     !
!                 End TagXC,TagYC,TagXV,TagYV info                    !
!                                                                     !
!*********************************************************************!
       endif
!      ________________________________________________________
!     |                                                        |
!     |    Distance of surronding cell center at each vertex   |
!     |             New way to collect weight values           |
!     |                 surrounding, weight                    |
!     |________________________________________________________|
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
                            weight(nv,j0) = 1.0d0/dsqrt(dxCV*dxCV+dyCV*dyCV)
                            dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                        endif
                    enddo
                enddo
!           ___________________
!           Number of surroun
                Dimsurrounding(nv) = j0
            ENDDO
!     ====================================
!     =====  START PARALLEL OPTION =======
!     ====================================
#     else
!       ------------------------------------------------
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
                        weight(nv,j0) = 1.0d0/dsqrt(dxCV*dxCV+dyCV*dyCV)
                        dlVsum(nv)=  dlVsum(nv) + weight(nv,j0)
                    endif
                enddo
            enddo
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
                        weight(nv,j0) = 1.0d0/dsqrt(dxCV*dxCV+dyCV*dyCV)
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
                    ntau = dsqrt(taux*taux+tauy*tauy)
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

      do nv=1,N_VERT
         normxv(nv) = 0.0d0
         normyv(nv) = 0.0d0
      enddo
!     ---------------------------------------------
!     For problems with straight bc
      if(XPB*YPB .ne. 1) then
        do nv=1,N_VERT
          if(XPB .eq. 0) then
            if(tagXV(nv) .eq. 1) then
                normxv(nv) = -1.0d0
            elseif(tagXV(nv) .eq. 2) then
                normxv(nv) = 1.0d0
            endif
          endif
!       -----------------------------
          if(YPB .eq. 0) then
            if(tagYV(nv) .eq. 1) then
                normyv(nv) = -1.0d0
            elseif(tagYV(nv) .eq. 2) then
                normyv(nv) = 1.0d0
            endif
          endif
!       -----------------------------
!        Unit vector
           if(tagXV(nv)*tagYV(nv) .gt. 0) then
            norx = normxv(nv)
            nory = normyv(nv)
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
#    ifndef KeyParallel
!     ====================================
!     ==========  SEQUENTIAL =============
!      ________________________________________________________
!     |                                                        |
!     |               Distance dn of normal points             |
!     |________________________________________________________|

      do nv=1,N_VERT
         dn(nv)  = 0.0d0
      enddo
      do i=1,N_CELL0
        do j=1,3
        nv=No_vp(i,j)
            if (nbev(nv).ne.0) then
               !dn(nv) = 0.5d0*dlCV(i,j)
               dn(nv) = 0.51d0*min(dlVV(i,1),dlVV(i,2),dlVV(i,3))
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
         icxn1(nv) = 2000!nv
         icxn2(nv) = 1000!nv
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
!           Second point: xn1 in the normal direction
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
#     else
#     endif
!*********************************************************************!
!                                                                     !
!        3 Edge points (xme(i,j),yme(i,j)) of each element i          !
!        & distance between the cell-centers and this edge point      !
!                                                                     !
!*********************************************************************!

        EdgePointOption = 1
!     ________________________________________________________________
!     Option 1: Middle point between vertices
        IF (EdgePointOption.eq.1) THEN
            do i=1,N_CELL0
                do j=1,3
                    jj=j+1
                    if (jj.gt.3) jj=jj-3
                    jv1 = No_vp(i,j)
                    jv2 = No_vp(i,jj)
                    xme(i,j) = 0.5d0*(xv(jv2)+xv(jv1))
                    yme(i,j) = 0.5d0*(yv(jv2)+yv(jv1))
!     ---------------------------------------
!     Intersection point between the line connecting
!     cell centers and the edge
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
!   --------------------------------------
!     Get distance between face center and intersection point
                    xm1 = xme(i,j) - xmo(i,j)
                    ym1 = yme(i,j) - ymo(i,j)
                    doe(i,j) = dsqrt(xm1*xm1+ym1*ym1)
!   ------------------------------------
!     determine the direction with e2
                    xm2 = xm1*dxVV(i,j)+ym1*dyVV(i,j)
                    if(xm2 .lt. 0.0d0) then
                        doe(i,j) = -doe(i,j)
                    endif
                enddo
            enddo
!     ________________________________________________________________
!     Option 2: Middle point between cell-centers
        ELSEIF (EdgePointOption.eq.2) THEN
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
!     Option 3: Perpendicular point to the edge
        ELSEIF (EdgePointOption.eq.3) THEN
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
!     Option 4: Average of the two triangle centroids around the edge
        ELSEIF (EdgePointOption.eq.4) THEN
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
!*********************************************************************!
!                                                                     !
!             Variables for the least square technique                !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |        dxCC, dyCC, dlCC, sum_xc, sum_yc & sum_xcyc     |
!     !        xe1,ye1, unit vector to neighbour cell center   !
!     !        xe2,ye2, (xnn-xe1),(ynn-ye1)                    !
!     |________________________________________________________|

        do i=1,N_CELL0
            sumX = 0.0d0
            sumY = 0.0d0
            sumXY= 0.0d0
            do j=1,3
                jc = No_cp(i,j)
                dxCC(i,j) = xc(jc)-xc(i)
                dyCC(i,j) = yc(jc)-yc(i)
                dlCC(i,j) = dsqrt(dxCC(i,j)**2+dyCC(i,j)**2)
                sumX = sumX + dxCC(i,j)*dxCC(i,j)
                sumY = sumY + dyCC(i,j)*dyCC(i,j)
                sumXY= sumXY+ dxCC(i,j)*dyCC(i,j)
!               -----------------------------
!               get vector e1
                nfn1 = dxCC(i,j)*xnn(i,j)+dyCC(i,j)*ynn(i,j)
                xe1(i,j) =  dxCC(i,j)/nfn1
                ye1(i,j) =  dyCC(i,j)/nfn1
!               -----------------------------
!               get vector e2
                xe2(i,j) = xnn(i,j)-xe1(i,j)
                ye2(i,j) = ynn(i,j)-ye1(i,j)

#               ifdef KeyKimChoi
!               -----------------------------
!               get vector n1
                xee = dxCC(i,j)/dlCC(i,j)
                yee = dyCC(i,j)/dlCC(i,j)
!               -----------------------------
!               get dle1 = 1/cos(theta)
                nfn1 = xee*xnn(i,j)+yee*ynn(i,j)
                xee =  xee/nfn1
                yee =  yee/nfn1
                dle1(i,j) = 1.0d0/nfn1
                dle3(i,j) = 0.!-1.0d0*dle1(i,j)/dlCC(i,j)
!               -----------------------------
!               get dle2 = tan(theta) with sign
                xee = xnn(i,j)-xee
                yee = ynn(i,j)-yee
!               ------------------------------
                nfn1 = xee*dxVV(i,j)+yee*dyVV(i,j)
                if(nfn1 .lt. 0.) then
                   dle2(i,j) = -1.0d0*dsqrt(xee*xee+yee*yee)
                else
                   dle2(i,j) = dsqrt(xee*xee+yee*yee)
                endif

#               ifdef KeyMahesh
                   dle2(i,j) = 0.
                   dle3(i,j) = 0.
                   doe(i,j)  = 0.
#               endif

#               endif
            enddo
            sum_xc2(i) = sumX
            sum_yc2(i) = sumY
            sum_xcyc(i)= sumXY
        enddo
!      ________________________________________________________
!     |                                                        |
!     |        dxCV, dyCV, dlCV, sum_xv, sum_yv & sum_xvyv     |
!     |________________________________________________________|
!    -----------------------------------------------------------------
        do i=1,N_CELL0
            sumX = 0.0d0
            sumY = 0.0d0
            sumXY= 0.0d0
            do j=1,3
                jc = No_vp(i,j)
                dxCV = xv(jc)-xc(i)
                dyCV = yv(jc)-yc(i)
                dxCVi(i,j) = dxCV
                dyCVi(i,j) = dyCV
                dlCV = dsqrt(dxCV**2+dyCV**2)
                sumX = sumX + dxCV*dxCV
                sumY = sumY + dyCV*dyCV
                sumXY= sumXY+ dxCV*dyCV
            enddo
            sum_xv2(i) = sumX
            sum_yv2(i) = sumY
            sum_xvyv(i)= sumXY
        enddo
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
100      continue
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
200      continue
      ENDDO

      Numictec = num
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
