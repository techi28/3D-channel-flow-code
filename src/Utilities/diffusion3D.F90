!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                             DIFFUSION 3D                            !
!                              March 2017                             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                             Bmv1T,Bmv2T,Bmv3T,                   &
                             Bmv1B,Bmv2B,Bmv3B,                   &
!                            --------------------------------------                             
                             Gamx,Gamy,Gamz,                      &
!                            --------------------------------------                             
                             xc,yc,sig,dsig,No_cp,nbe,            &
                             xv,yv,sigv,dsigv,No_vp,nbev)    
 
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the diffudion contribution to the   !
!    general linear system.                                           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> Am0   |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  | <--> Am1   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  | <--> Am2   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  | <--> Am3   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  | <--> AmT   |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  | <--> AmB   |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  ! 
!  |____________|____________|_____________________________________|  !
!  | <--> Bmv1T |(N_CELL0,NZ)| matrix coeff. vertex 1 Top          |  !
!  | <--> Bmv2T |(N_CELL0,NZ)| matrix coeff. vertex 2 Top          |  ! 
!  | <--> Bmv3T |(N_CELL0,NZ)| matrix coeff. vertex 3 Top          |  ! 
!  | <--> Bmv1B |(N_CELL0,NZ)| matrix coeff. vertex 1 Bottom       |  ! 
!  | <--> Bmv2B |(N_CELL0,NZ)| matrix coeff. vertex 2 Bottom       |  ! 
!  | <--> Bmv3B |(N_CELL0,NZ)| matrix coeff. vertex 3 Bottom       |  ! 
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ+1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ+1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters and variables used:                            !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |    AreaCell|  Area of the cell                                |  ! 
!  |____________|__________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,nc,i,j,k| Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name        |               Description                     |  !  
!  |_______________|_______________________________________________|  !  
!  | AA0,AAj(1:3)  | Auxiliar horizontal matrix coefficients       |  !
!  | AAT,AAB       | Auxiliar vertical matrix coefficients         |  !
!  | dtoVol        | = dt/Volume(i,k)                              |  !
!  |_______________|_______________________________________________|  !
!  | jc,jj         | neighborn loop index                          |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |                    -  grandientLSM                            |  !
!  |                    -  massflux                                |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
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

      real*8, dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv3B(N_CELL0,NZ)
!     -------------------------------------- 
      real*8, dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamz(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      real*8, dimension(:)   :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     -------------------------------------- 
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)       
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:) :: AA(1:3)
      real*8 :: AA0,AAB,AAT
      real*8 :: Bv1T,Bv2T,Bv3T,Bv1B,Bv2B,Bv3B
!     --------------------------------------
      real*8 :: nxAreaFaceij,nyAreaFaceij,nzAreaFaceij
      real*8 :: sigxij,sigyij,sigzij
      real*8 :: Gamxij,Gamyij,Gamzij,Gamsig
!     --------------------------------------
      real*8 :: VoluReg
      real*8 :: nxAreaoVol,nyAreaoVol,nzAreaoVol
      real*8 :: a1,a2,a3
      real*8 :: b1,b2,b3
      real*8 :: cxx,cyy,czz
      real*8 :: sumaix,sumaiy,sumaiz
      real*8 :: sumajx,sumajy,sumajz
!     --------------------------------------
      real*8 :: fv1Tx,fv2Tx,fv3Tx
      real*8 :: fv1Ty,fv2Ty,fv3Ty
      real*8 :: fv1Tz,fv2Tz,fv3Tz
      real*8 :: fv1Bx,fv2Bx,fv3Bx
      real*8 :: fv1By,fv2By,fv3By
      real*8 :: fv1Bz,fv2Bz,fv3Bz
!     --------------------------------------
      integer:: jc,s,ss
      integer:: jj,jv1,jv2,jv3

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: diffusion3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
!*********************************************************************!
!                                                                     !
!                     Calculation of the gradient                     !
!                                                                     !
!*********************************************************************!

      DO k=2,NZ-1
         do i=1,N_CELL0

            AA0   = 0.0d0
            AA(1) = 0.0d0
            AA(2) = 0.0d0
            AA(3) = 0.0d0
            AAT   = 0.0d0
            AAB   = 0.0d0
            Bv1T  = 0.0d0
            Bv2T  = 0.0d0
            Bv3T  = 0.0d0
            Bv1B  = 0.0d0
            Bv2B  = 0.0d0
            Bv3B  = 0.0d0
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |               Horizontal neighbors (j=1,3)               |
!          |----------------------------------------------------------|
!          |__________________________________________________________|
           
            do j=1,3
	       sumaix = 0.0d0
	       sumaiy = 0.0d0    
	       sumaiz = 0.0d0 
            
               sumajx = 0.0d0
	       sumajy = 0.0d0
	       sumajz = 0.0d0

               fv1Tx  = 0.0d0
               fv1Ty  = 0.0d0
               fv1Tz  = 0.0d0

               fv2Tx  = 0.0d0
               fv2Ty  = 0.0d0
               fv2Tz  = 0.0d0

               fv3Tx  = 0.0d0
               fv3Ty  = 0.0d0
               fv3Tz  = 0.0d0

               fv1Bx  = 0.0d0
               fv1By  = 0.0d0
               fv1Bz  = 0.0d0

               fv2Bx  = 0.0d0
               fv2By  = 0.0d0
               fv2Bz  = 0.0d0

               fv3Bx  = 0.0d0
               fv3By  = 0.0d0
               fv3Bz  = 0.0d0
!              _______________________________
!             |                               |
!             |        Index (interface)      |
!             |_______________________________|

!              -------------------------------
!              Neighborn
	       jc  = No_cp(i,j)
!              -------------------------------
!              Vertices
               jj=j+1
	       if (jj.gt.3) jj=jj-3
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)
!              _______________________________
!             |                               |
!             |      Volume of the region     |
!             |_______________________________|

               VoluReg = dlVV(i,j)*dsigv(k-1)*(dhCE(i,j))/3.0d0 
!              _______________________________
!             |                               |
!             |          Region faces         |
!             |_______________________________|

!              _________________________________________________
!              Region face 1 (lateral)        
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(yc(jc)-yv(jv1))*dsigv(k-1)/VoluReg
               nyAreaoVol =-0.5d0*(xc(jc)-xv(jv1))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell neighbor j 
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Cell-centers assignation 
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
!              -------------------------------
!              Vertex assignation             
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy                
               endif
!              _________________________________________________
!              Region face 2 (lateral)
!              -------------------------------
!              normal*Area
	       nxAreaoVol = 0.5d0*(yv(jv2)-yc(jc))*dsigv(k-1)/VoluReg
	       nyAreaoVol =-0.5d0*(xv(jv2)-xc(jc))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell neighbor j 
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Cell-centers assignation
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
!              -------------------------------
!              Vertex assignation
               if (j.eq.1) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.2) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
               elseif (j.eq.3) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy                
               endif
!              _________________________________________________
!              Region face 3 (lateral)
!              ----------------------------
!              normal*Area
	       nxAreaoVol =  0.5d0*(yc(i)-yv(jv2))*dsigv(k-1)/VoluReg
	       nyAreaoVol = -0.5d0*(xc(i)-xv(jv2))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell-center i  
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Cell-center i 
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
!              -------------------------------
!              Known vertex values 
               if (j.eq.1) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.2) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
               elseif (j.eq.3) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy                
               endif
!              _________________________________________________
!              Region face 4 (lateral)        
!              -------------------------------
!              normal*Area
	       nxAreaoVol =  0.5d0*(yv(jv1)-yc(i))*dsigv(k-1)/VoluReg
	       nyAreaoVol = -0.5d0*(xv(jv1)-xc(i))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell-center i  
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Assignation to cell-center i 
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
!              -------------------------------
!              Assignation vertex values
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy                
               endif
!              _________________________________________________
!              Region face 5 (inferior out)
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(jc)
               a2 = yv(jv1)-yc(jc)
               a3 = sigv(k-1)-sig(k)
               b1 = xv(jv2)-xc(jc)
               b2 = yv(jv2)-yc(jc)
               b3 = sigv(k-1)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Constant cell neighborn j   
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              -------------------------------
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
               sumajz = sumajz + czz
!              -------------------------------
!              Known vertex values
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy  
	          fv3Bz = fv3Bz + czz  
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
               endif
!              _________________________________________________
!              Region face 6 (superior out)
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(jc)
               a2 = yv(jv2)-yc(jc)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv1)-xc(jc)
               b2 = yv(jv1)-yc(jc)
               b3 = sigv(k)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Constant cell neighborn j   
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              -------------------------------
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
               sumajz = sumajz + czz
!              -------------------------------
!              Known vertex values
               if (j.eq.1) then
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
               elseif (j.eq.2) then
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy  
	          fv3Tz = fv3Tz + czz  
               elseif (j.eq.3) then
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
               endif
!              _________________________________________________
!              Region face 7 (inferior in)
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k-1)-sig(k)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k-1)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Constant cell-center i   
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              -------------------------------
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
               sumaiz = sumaiz + czz
!              -------------------------------
!              Known vertex values face 7
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy  
	          fv3Bz = fv3Bz + czz  
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
               endif
!              _________________________________________________
!              Region face 8 (superior in)
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Constant cell-center i   
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              -------------------------------
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
               sumaiz = sumaiz + czz
!              -------------------------------
!              Known vertex values face 8
               if (j.eq.1) then
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
               elseif (j.eq.2) then
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy  
	          fv3Tz = fv3Tz + czz  
               elseif (j.eq.3) then
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
               endif
!              _______________________________
!             |                               |
!             |         Contributions         |
!             |_______________________________|

!              _______________________________
!              normal*Area (interface)
               nxAreaFaceij =  dyVV(i,j)*dsigv(k-1)
               nyAreaFaceij = -dxVV(i,j)*dsigv(k-1) 
!              _______________________________
!              Gamma & sigma (interface)  
               Gamxij = 0.5d0*(Gamx(i,k)+Gamx(jc,k))
               Gamyij = 0.5d0*(Gamy(i,k)+Gamy(jc,k))
!              _______________________________
!              sigma (interface)
               sigxij = 0.5d0*(sigmax(i,k)+sigmax(jc,k))
               sigyij = 0.5d0*(sigmay(i,k)+sigmay(jc,k))
!              _______________________________
!              Coefficients
!              ---------------------------------
!              Matrix coeff. cell-center i
               AA0  =  Gamxij*(sumaix+sigxij*sumaiz)*nxAreaFaceij &
                     + Gamyij*(sumaiy+sigyij*sumaiz)*nyAreaFaceij + AA0  
!              ---------------------------------
!              Matrix coeff. cell-neighborn j
               AA(j)=  Gamxij*(sumajx+sigxij*sumajz)*nxAreaFaceij &
                     + Gamyij*(sumajy+sigyij*sumajz)*nyAreaFaceij
!              ---------------------------------
!              Matrix coeff. vertices
               if (j.eq.1) then
	          Bv1B =  Gamxij*(fv1Bx+sigxij*fv1Bz)*nxAreaFaceij &
                        + Gamyij*(fv1By+sigyij*fv1Bz)*nyAreaFaceij + Bv1B
	          Bv1T =  Gamxij*(fv1Tx+sigxij*fv1Tz)*nxAreaFaceij &
                        + Gamyij*(fv1Ty+sigyij*fv1Tz)*nyAreaFaceij + Bv1T
	          Bv2B =  Gamxij*(fv2Bx+sigxij*fv2Bz)*nxAreaFaceij &
                        + Gamyij*(fv2By+sigyij*fv2Bz)*nyAreaFaceij + Bv2B
	          Bv2T =  Gamxij*(fv2Tx+sigxij*fv2Tz)*nxAreaFaceij &
                        + Gamyij*(fv2Ty+sigyij*fv2Tz)*nyAreaFaceij + Bv2T
               elseif (j.eq.2) then
	          Bv2B =  Gamxij*(fv2Bx+sigxij*fv2Bz)*nxAreaFaceij &
                        + Gamyij*(fv2By+sigyij*fv2Bz)*nyAreaFaceij + Bv2B
	          Bv2T =  Gamxij*(fv2Tx+sigxij*fv2Tz)*nxAreaFaceij &
                        + Gamyij*(fv2Ty+sigyij*fv2Tz)*nyAreaFaceij + Bv2T
	          Bv3B =  Gamxij*(fv3Bx+sigxij*fv3Bz)*nxAreaFaceij &
                        + Gamyij*(fv3By+sigyij*fv3Bz)*nyAreaFaceij + Bv3B
	          Bv3T =  Gamxij*(fv3Tx+sigxij*fv3Tz)*nxAreaFaceij &
                        + Gamyij*(fv3Ty+sigyij*fv3Tz)*nyAreaFaceij + Bv3T
               elseif (j.eq.3) then
	          Bv3B =  Gamxij*(fv3Bx+sigxij*fv3Bz)*nxAreaFaceij &
                        + Gamyij*(fv3By+sigyij*fv3Bz)*nyAreaFaceij + Bv3B
	          Bv3T =  Gamxij*(fv3Tx+sigxij*fv3Tz)*nxAreaFaceij &
                        + Gamyij*(fv3Ty+sigyij*fv3Tz)*nyAreaFaceij + Bv3T
	          Bv1B =  Gamxij*(fv1Bx+sigxij*fv1Bz)*nxAreaFaceij &
                        + Gamyij*(fv1By+sigyij*fv1Bz)*nyAreaFaceij + Bv1B
	          Bv1T =  Gamxij*(fv1Tx+sigxij*fv1Tz)*nxAreaFaceij &
                        + Gamyij*(fv1Ty+sigyij*fv1Tz)*nyAreaFaceij + Bv1T
               endif
	    enddo 
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |                Vertical TOP neighbor: j=4                |
!          |----------------------------------------------------------|
!          |__________________________________________________________|

	    sumaix = 0.0d0
	    sumaiy = 0.0d0
	    sumaiz = 0.0d0

	    sumajx = 0.0d0
	    sumajy = 0.0d0
	    sumajz = 0.0d0

            fv1Tx  = 0.0d0
            fv1Ty  = 0.0d0
            fv1Tz  = 0.0d0

            fv2Tx  = 0.0d0
            fv2Ty  = 0.0d0
            fv2Tz  = 0.0d0

            fv3Tx  = 0.0d0
            fv3Ty  = 0.0d0
            fv3Tz  = 0.0d0
!           _______________________________
!          |                               |
!          |      Volume of the region     |
!          |_______________________________|

            VoluReg = (1.0d0/3.0d0)*AreaCell(i)*(dsig(k))

!           _______________________________
!          |                               |
!          |         Region Faces          |
!          |_______________________________|

            do s=1,3
!              -------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss) 
!              _________________________________________________
!              Region face 1,2,3

!              -------------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k)-sig(k)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0 
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Cell-center i contribution
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
               sumaiz = sumaiz + czz
!              ------------------------------
!              Vertex contribution
               if (s.eq.1) then
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
               elseif (s.eq.2) then
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
               elseif (s.eq.3) then
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz               
               endif
!              _________________________________________________
!              Region faces 4,5,6
!              -------------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k)-sig(k+1)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k)-sig(k+1)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Neighborn center j contribution
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
               sumajz = sumajz + czz
!              ------------------------------
!              Vertex contribution
               if (s.eq.1) then
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
               elseif (s.eq.2) then
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz 
               elseif (s.eq.3) then
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz                
               endif
            enddo
!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|
!           _______________________________
!           normal*Area (interface)
            nzAreaFaceij = areaCell(i)
!           _______________________________
!           Gamma (interface)  
            Gamxij = 0.5d0*(Gamx(i,k)+Gamx(i,k+1))
            Gamyij = 0.5d0*(Gamy(i,k)+Gamy(i,k+1))
            Gamzij = 0.5d0*(Gamz(i,k)+Gamz(i,k+1))
!           _______________________________
!           sigma (interface)
            sigxij = 0.5d0*(sigmax(i,k)+sigmax(i,k+1))
            sigyij = 0.5d0*(sigmay(i,k)+sigmay(i,k+1))
            sigzij = 0.5d0*(sigmaz(i,k)+sigmaz(i,k+1))
            Gamsig = Gamxij*(sigxij**2) &
                    +Gamyij*(sigyij**2) &
                    +Gamzij*(sigzij**2)
!           _______________________________
!           Coefficients
!           -------------------------------
!           Matrix coeff. cell-center i
            AA0 = (Gamxij*sigxij*sumaix &
                  +Gamyij*sigyij*sumaiy &
                  +Gamsig*sumaiz)*nzAreaFaceij &
                  +AA0        	
!           -------------------------------
!           Matrix coeff. cell-center j=T  
            AAT = (Gamxij*sigxij*sumajx &
                  +Gamyij*sigyij*sumajy &
                  +Gamsig*sumajz)*nzAreaFaceij
!           -------------------------------
!           Matrix coeff. vertices
	    Bv1T = ( Gamxij*sigxij*fv1Tx &
                    +Gamyij*sigyij*fv1Ty &
                    +Gamsig*fv1Tz)*nzAreaFaceij + Bv1T
	    Bv2T = ( Gamxij*sigxij*fv2Tx &
                    +Gamyij*sigyij*fv2Ty &
                    +Gamsig*fv2Tz)*nzAreaFaceij + Bv2T
	    Bv3T = ( Gamxij*sigxij*fv3Tx &
                    +Gamyij*sigyij*fv3Ty &
                    +Gamsig*fv3Tz)*nzAreaFaceij + Bv3T
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |              Vertical BOTTOM neighbor: (j=5)             |
!          |----------------------------------------------------------|
!          |__________________________________________________________|

	    sumaix = 0.0d0
	    sumajy = 0.0d0
	    sumaiz = 0.0d0
	    sumajx = 0.0d0
	    sumaiy = 0.0d0
	    sumajz = 0.0d0
            fv1Bx  = 0.0d0
            fv1By  = 0.0d0
            fv1Bz  = 0.0d0
            fv2Bx  = 0.0d0
            fv2By  = 0.0d0
            fv2Bz  = 0.0d0
            fv3Bx  = 0.0d0
            fv3By  = 0.0d0
            fv3Bz  = 0.0d0
!           _______________________________
!          |                               |
!          |      Volume of the region     |
!          |_______________________________|
            
            VoluReg = (1.0d0/3.0d0)*AreaCell(i)*(dsig(k-1)) 

!           _______________________________
!          |                               |
!          |         Region faces          |
!          |_______________________________|

            do s=1,3
!              -------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss) 
!              _________________________________________________
!              Region face 1,2,3
!              -------------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k-1)-sig(k-1)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k-1)-sig(k-1)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Cell neighborn j contribution
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
               sumajz = sumajz + czz
!              ------------------------------
!              Vertex contribution   
               if (s.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
               elseif (s.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz  
               elseif (s.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz                
               endif    
!              _________________________________________________
!              Region faces 4,5,6
!              -------------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k-1)-sig(k)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k-1)-sig(k)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              ------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Cell-center i contribution
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
               sumaiz = sumaiz + czz
!              ------------------------------
!              Vertex contribution
               if (s.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
               elseif (s.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz 
               elseif (s.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz                
               endif
            enddo
!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|
!           _______________________________
!           normal*Area (interface)
            nzAreaFaceij = -areaCell(i)
!           _______________________________
!           Gamma (interface)  
            Gamxij = 0.5d0*(Gamx(i,k)+Gamx(i,k-1))
            Gamyij = 0.5d0*(Gamy(i,k)+Gamy(i,k-1))
            Gamzij = 0.5d0*(Gamz(i,k)+Gamz(i,k-1))
!           _______________________________
!           sigma (interface)
            sigxij = 0.5d0*(sigmax(i,k)+sigmax(i,k-1))
            sigyij = 0.5d0*(sigmay(i,k)+sigmay(i,k-1))
            sigzij = 0.5d0*(sigmaz(i,k)+sigmaz(i,k-1))
            Gamsig = Gamxij*(sigxij**2) &
                    +Gamyij*(sigyij**2) &
                    +Gamzij*(sigzij**2)
!           _______________________________
!           Coefficients
!           -------------------------------
!           Matrix coeff. cell-center i
            AA0 = (Gamxij*sigxij*sumaix &
                  +Gamyij*sigyij*sumaiy &
                  +Gamsig*sumaiz)*nzAreaFaceij + AA0      	
!           -------------------------------
!           Matrix coeff. cell-center j=B  
            AAB = (Gamxij*sigxij*sumajx &
                  +Gamyij*sigyij*sumajy &
                  +Gamsig*sumajz)*nzAreaFaceij
!           -------------------------------
!           Matrix coeff. vertices
	    Bv1B = ( Gamxij*sigxij*fv1Tx &
                    +Gamyij*sigyij*fv1Ty &
                    +Gamsig*fv1Tz)*nzAreaFaceij + Bv1B
	    Bv2B = ( Gamxij*sigxij*fv2Tx &
                    +Gamyij*sigyij*fv2Ty &
                    +Gamsig*fv2Tz)*nzAreaFaceij + Bv2B
	    Bv3B = ( Gamxij*sigxij*fv3Tx &
                    +Gamyij*sigyij*fv3Ty &
                    +Gamsig*fv3Tz)*nzAreaFaceij + Bv3B
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |              Final coefficients contributions            |
!          |----------------------------------------------------------|
!          |__________________________________________________________|

!           --------------------------------------
!           Matrix cell-center contribution
	    Am0(i,k) = Am0(i,k) + AA0    
            Am1(i,k) = Am1(i,k) + AA(1)   
            Am2(i,k) = Am2(i,k) + AA(2)   
            Am3(i,k) = Am3(i,k) + AA(3)
            AmT(i,k) = AmT(i,k) + AAT             
            AmB(i,k) = AmB(i,k) + AAB
!           --------------------------------------
!           Matrix vextex contribution
	    Bmv1T(i,k) = Bmv1T(i,k) + Bv1T    
            Bmv2T(i,k) = Bmv2T(i,k) + Bv2T   
            Bmv3T(i,k) = Bmv3T(i,k) + Bv3T   
            Bmv1B(i,k) = Bmv1B(i,k) + Bv1B
            Bmv2B(i,k) = Bmv2B(i,k) + Bv2B             
            Bmv3B(i,k) = Bmv3B(i,k) + Bv3B  
         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: Diffusion3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      DIFFUSION 3D FOR NEUMANN BC                    !
!                               Nov 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE diffusionNeumann3D(Am0,Am1,Am2,Am3,AmT,AmB,      &
                                    Bmv1T,Bmv2T,Bmv3T,            &
                                    Bmv1B,Bmv2B,Bmv3B,Nm,         &
                                    Gamx,Gamy,Gamz,               &
                                    xc,yc,sig,dsig,No_cp,nbe,     &
                                    xv,yv,sigv,dsigv,No_vp,nbev)    
 
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine approximates the diffudion contribution to the   !
!    general linear system.                                           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output & Input variables:                                        !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> Am0   |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  | <--> Am1   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  | <--> Am2   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  | <--> Am3   |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  | <--> AmT   |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  | <--> AmB   |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  ! 
!  |____________|____________|_____________________________________|  !
!  | <--> Bmv1T |(N_CELL0,NZ)| matrix coeff. vertex 1 Top          |  !
!  | <--> Bmv2T |(N_CELL0,NZ)| matrix coeff. vertex 2 Top          |  ! 
!  | <--> Bmv3T |(N_CELL0,NZ)| matrix coeff. vertex 3 Top          |  ! 
!  | <--> Bmv1B |(N_CELL0,NZ)| matrix coeff. vertex 1 Bottom       |  ! 
!  | <--> Bmv2B |(N_CELL0,NZ)| matrix coeff. vertex 2 Bottom       |  ! 
!  | <--> Bmv3B |(N_CELL0,NZ)| matrix coeff. vertex 3 Bottom       |  ! 
!  |____________|____________|_____________________________________|  !
!  | <--> Nm    |(N_CELL0,NZ)| Term regarding Neumann BC           |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ+1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ+1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters and variables used:                            !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |    NZ      |  Number of vertical points                       |  ! 
!  |    AreaCell|  Area of the cell                                |  ! 
!  |____________|__________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,nc,i,j,k| Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name        |               Description                     |  !  
!  |_______________|_______________________________________________|  !  
!  | AA0,AAj(1:3)  | Auxiliar horizontal matrix coefficients       |  !
!  | AAT,AAB       | Auxiliar vertical matrix coefficients         |  !
!  | dtoVol        | = dt/Volume(i,k)                              |  !
!  |_______________|_______________________________________________|  !
!  | jc,jj         | neighborn loop index                          |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |                    -  grandientLSM                            |  !
!  |                    -  massflux                                |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
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

      real*8, dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv3B(N_CELL0,NZ)
!     -------------------------------------- 
      real*8, dimension(:,:) :: Nm(N_CELL0,NZ)
!     -------------------------------------- 
      real*8, dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamz(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      real*8, dimension(:)   :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     -------------------------------------- 
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)       
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:) :: AA(1:3)
      real*8 :: AA0,AAB,AAT
      real*8 :: Bv1T,Bv2T,Bv3T,Bv1B,Bv2B,Bv3B
!     --------------------------------------
      real*8 :: nxAreaFaceij,nyAreaFaceij,nzAreaFaceij
      real*8 :: sigxij,sigyij,sigzij
      real*8 :: Gamxij,Gamyij,Gamzij,Gamsig
!     --------------------------------------
      real*8 :: VoluReg
      real*8 :: nxAreaoVol,nyAreaoVol,nzAreaoVol
      real*8 :: a1,a2,a3
      real*8 :: b1,b2,b3
      real*8 :: cxx,cyy,czz
      real*8 :: sumaix,sumaiy,sumaiz
      real*8 :: sumajx,sumajy,sumajz
!     --------------------------------------
      real*8 :: fv1Tx,fv2Tx,fv3Tx
      real*8 :: fv1Ty,fv2Ty,fv3Ty
      real*8 :: fv1Tz,fv2Tz,fv3Tz
      real*8 :: fv1Bx,fv2Bx,fv3Bx
      real*8 :: fv1By,fv2By,fv3By
      real*8 :: fv1Bz,fv2Bz,fv3Bz
!     --------------------------------------
      integer:: jc,s,ss
      integer:: jj,jv1,jv2,jv3
!     --------------------------------------
      real*8  :: nnx,nny,nnz,x,y,z
      real*8  :: dfBdn,Neumanndfdn3D

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: diffusion3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!               Calculation of the gradient: NEUMANN BC               !
!                                                                     !
!*********************************************************************!

      DO k=2,NZ-1
         do i=1,N_CELL0

            AA0   = 0.0d0
            AA(1) = 0.0d0
            AA(2) = 0.0d0
            AA(3) = 0.0d0
            AAT   = 0.0d0
            AAB   = 0.0d0
            Bv1T  = 0.0d0
            Bv2T  = 0.0d0
            Bv3T  = 0.0d0
            Bv1B  = 0.0d0
            Bv2B  = 0.0d0
            Bv3B  = 0.0d0
            Nm(i,k) = 0.0d0
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |               Horizontal neighbors (j=1,3)               |
!          |----------------------------------------------------------|
!          |__________________________________________________________|
           
            do j=1,3
               sumajx = 0.0d0
	       sumaix = 0.0d0
	       sumajy = 0.0d0
	       sumaiy = 0.0d0             
               fv1Tx  = 0.0d0
               fv1Ty  = 0.0d0
               fv2Tx  = 0.0d0
               fv2Ty  = 0.0d0
               fv3Tx  = 0.0d0
               fv3Ty  = 0.0d0
               fv1Bx  = 0.0d0
               fv1By  = 0.0d0
               fv2Bx  = 0.0d0
               fv2By  = 0.0d0
               fv3Bx  = 0.0d0
               fv3By  = 0.0d0
!              _______________________________
!             |                               |
!             |        Index (interface)      |
!             |_______________________________|

!              -------------------------------
!              Neighborn
	       jc  = No_cp(i,j)
!              -------------------------------
!              Vertices
               jj=j+1
	       if (jj.gt.3) jj=jj-3
	       jv1 = No_vp(i,j)
	       jv2 = No_vp(i,jj)

!              xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                                  BOUNDARY FACES
!              xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               IF ((jc.lt.1).or.(jc.gt.N_CELL0)) THEN
                  nnx = normxc(i,j)
                  nny = normyc(i,j)
                  nnz = 0.0d0
                  x = 0.5d0*(xv(jv1)+xv(jv2))
                  y = 0.5d0*(yv(jv1)+yv(jv2))
                  z = sig(k)
                  dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
                  Nm(i,k) = Nm(i,k) + dfBdn*Gamx(i,k)&
                                     *dlVV(i,j)*dsigv(k-1)

!              xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                                  INSIDE FACES
!              xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

               ELSE
!              _______________________________
!             |                               |
!             |      Volume of the region     |
!             |_______________________________|

               VoluReg = dlVV(i,j)*dsigv(k-1)*(dhCE(i,j))/3.0d0 
!              _______________________________
!             |                               |
!             |          Region faces         |
!             |_______________________________|

!              _________________________________________________
!              Region face 1         
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(yc(jc)-yv(jv1))*dsigv(k-1)/VoluReg
               nyAreaoVol =-0.5d0*(xc(jc)-xv(jv1))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell neighbor j 
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Cell-centers assignation 
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
!              -------------------------------
!              Vertex assignation             
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy                
               endif
!              _________________________________________________
!              Region face 2  
!              -------------------------------
!              normal*Area
	       nxAreaoVol = 0.5d0*(yv(jv2)-yc(jc))*dsigv(k-1)/VoluReg
	       nyAreaoVol =-0.5d0*(xv(jv2)-xc(jc))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell neighbor j 
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Cell-centers assignation
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
!              -------------------------------
!              Vertex assignation
               if (j.eq.1) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.2) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
               elseif (j.eq.3) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy                
               endif
!              _________________________________________________
!              Region face 3
!              ----------------------------
!              normal*Area
	       nxAreaoVol =  0.5d0*(yc(i)-yv(jv2))*dsigv(k-1)/VoluReg
	       nyAreaoVol = -0.5d0*(xc(i)-xv(jv2))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell-center i  
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Cell-center i 
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
!              -------------------------------
!              Known vertex values 
               if (j.eq.1) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.2) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
               elseif (j.eq.3) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy                
               endif
!              _________________________________________________
!              Region face 4         
!              -------------------------------
!              normal*Area
	       nxAreaoVol =  0.5d0*(yv(jv1)-yc(i))*dsigv(k-1)/VoluReg
	       nyAreaoVol = -0.5d0*(xv(jv1)-xc(i))*dsigv(k-1)/VoluReg
!              -------------------------------
!              Constant cell-center i  
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
!              Assignation to cell-center i 
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
!              -------------------------------
!              Assignation vertex values
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy                
               endif
!              _________________________________________________
!              Region face 5,6
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(jc)
               a2 = yv(jv2)-yc(jc)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv1)-xc(jc)
               b2 = yv(jv1)-yc(jc)
               b3 = sigv(k)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
!              -------------------------------
!              Constant cell neighborn j   
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
               sumajx = sumajx + 2.0d0*cxx
               sumajy = sumajy + 2.0d0*cyy
!              -------------------------------
!              Known vertex values
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
!                 --------------------
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
!                 --------------------
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy  
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
!                 --------------------
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy                
               endif
!              _________________________________________________
!              Region face 7,8
!              ----------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k)-sig(k)
!              ----------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
!              -------------------------------
!              Constant cell-center i   
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
!              -------------------------------
               sumaix = sumaix + 2.0d0*cxx
               sumaiy = sumaiy + 2.0d0*cyy
!              -------------------------------
!              Known vertex values face 7,8
               if (j.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
!                 --------------------
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
               elseif (j.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
!                 --------------------
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy  
               elseif (j.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
!                 --------------------
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy                
               endif
!              _______________________________
!             |                               |
!             |         Contributions         |
!             |_______________________________|

!              _______________________________
!              normal*Area (interface)
               nxAreaFaceij =  dyVV(i,j)*dsigv(k-1)
               nyAreaFaceij = -dxVV(i,j)*dsigv(k-1) 
!              _______________________________
!              Gamma & sigma (interface)  
               Gamxij  = 0.5d0*(Gamx(i,k)+Gamx(jc,k))
               Gamyij  = 0.5d0*(Gamy(i,k)+Gamy(jc,k))
!              _______________________________
!              Coefficients
!              ---------------------------------
!              Matrix coeff. cell-center i
               AA0  =  (Gamxij*sumaix)*nxAreaFaceij &
                     + (Gamyij*sumaiy)*nyAreaFaceij + AA0  
!              ---------------------------------
!              Matrix coeff. cell-neighborn j
               AA(j)=  (Gamxij*sumajx)*nxAreaFaceij &
                     + (Gamyij*sumajy)*nyAreaFaceij
!              ---------------------------------
!              Matrix coeff. vertices
               if (j.eq.1) then
	          Bv1B =  (Gamxij*fv1Bx)*nxAreaFaceij &
                        + (Gamyij*fv1By)*nyAreaFaceij + Bv1B
	          Bv1T =  (Gamxij*fv1Tx)*nxAreaFaceij &
                        + (Gamyij*fv1Ty)*nyAreaFaceij + Bv1T
	          Bv2B =  (Gamxij*fv2Bx)*nxAreaFaceij &
                        + (Gamyij*fv2By)*nyAreaFaceij + Bv2B
	          Bv2T =  (Gamxij*fv2Tx)*nxAreaFaceij &
                        + (Gamyij*fv2Ty)*nyAreaFaceij + Bv2T
               elseif (j.eq.2) then
	          Bv2B =  (Gamxij*fv2Bx)*nxAreaFaceij &
                        + (Gamyij*fv2By)*nyAreaFaceij + Bv2B
	          Bv2T =  (Gamxij*fv2Tx)*nxAreaFaceij &
                        + (Gamyij*fv2Ty)*nyAreaFaceij + Bv2T
	          Bv3B =  (Gamxij*fv3Bx)*nxAreaFaceij &
                        + (Gamyij*fv3By)*nyAreaFaceij + Bv3B
	          Bv3T =  (Gamxij*fv3Tx)*nxAreaFaceij &
                        + (Gamyij*fv3Ty)*nyAreaFaceij + Bv3T
               elseif (j.eq.3) then
	          Bv3B =  (Gamxij*fv3Bx)*nxAreaFaceij &
                        + (Gamyij*fv3By)*nyAreaFaceij + Bv3B
	          Bv3T =  (Gamxij*fv3Tx)*nxAreaFaceij &
                        + (Gamyij*fv3Ty)*nyAreaFaceij + Bv3T
	          Bv1B =  (Gamxij*fv1Bx)*nxAreaFaceij &
                        + (Gamyij*fv1By)*nyAreaFaceij + Bv1B
	          Bv1T =  (Gamxij*fv1Tx)*nxAreaFaceij &
                        + (Gamyij*fv1Ty)*nyAreaFaceij + Bv1T
               endif

!             xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                                END OF CONDITIONAL
!             xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

              ENDIF

	    enddo 
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |                Vertical TOP neighbor: j=4                |
!          |----------------------------------------------------------|
!          |__________________________________________________________|

!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                             BOUNDARY TOP FACES
!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            IF (k.eq.NZ-1) THEN                         
               nnx = 0.0d0
               nny = 0.0d0
               nnz = 1.0d0
               x = xc(i)
               y = yc(i)
               z = 0.5d0*(sig(NZ-1)+sig(NZ))
               dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
               Nm(i,k) = Nm(i,k)+ dfBdn*Gamx(i,k)*AreaCell(i)

!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                              INSIDE TOP FACES
!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            ELSE

	    sumaix = 0.0d0
	    sumaiy = 0.0d0
	    sumaiz = 0.0d0

	    sumajx = 0.0d0
	    sumajy = 0.0d0
	    sumajz = 0.0d0

            fv1Tx  = 0.0d0
            fv1Ty  = 0.0d0
            fv1Tz  = 0.0d0

            fv2Tx  = 0.0d0
            fv2Ty  = 0.0d0
            fv2Tz  = 0.0d0

            fv3Tx  = 0.0d0
            fv3Ty  = 0.0d0
            fv3Tz  = 0.0d0
!           _______________________________
!          |                               |
!          |      Volume of the region     |
!          |_______________________________|

            VoluReg = (1.0d0/3.0d0)*AreaCell(i)*(dsig(k))

!           _______________________________
!          |                               |
!          |         Region Faces          |
!          |_______________________________|

            do s=1,3
!              -------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss) 
!              _________________________________________________
!              Region face 1,2,3

!              -------------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k)-sig(k)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k)-sig(k)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0 
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Cell-center i contribution
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
               sumaiz = sumaiz + czz
!              ------------------------------
!              Vertex contribution
               if (s.eq.1) then
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
               elseif (s.eq.2) then
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
               elseif (s.eq.3) then
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz               
               endif
!              _________________________________________________
!              Region faces 4,5,6
!              -------------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k)-sig(k+1)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k)-sig(k+1)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Neighborn center j contribution
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
               sumajz = sumajz + czz
!              ------------------------------
!              Vertex contribution
               if (s.eq.1) then
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
               elseif (s.eq.2) then
	          fv2Tx = fv2Tx + cxx
	          fv2Ty = fv2Ty + cyy
	          fv2Tz = fv2Tz + czz
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz 
               elseif (s.eq.3) then
	          fv3Tx = fv3Tx + cxx
	          fv3Ty = fv3Ty + cyy
	          fv3Tz = fv3Tz + czz
	          fv1Tx = fv1Tx + cxx
	          fv1Ty = fv1Ty + cyy
	          fv1Tz = fv1Tz + czz                
               endif
            enddo
!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|
!           _______________________________
!           normal*Area (interface)
            nzAreaFaceij = areaCell(i)
!           _______________________________
!           Gamma (interface)  
            Gamxij = 0.5d0*(Gamx(i,k)+Gamx(i,k+1))
            Gamyij = 0.5d0*(Gamy(i,k)+Gamy(i,k+1))
            Gamzij = 0.5d0*(Gamz(i,k)+Gamz(i,k+1))
!           _______________________________
!           sigma (interface)
            sigxij = 0.5d0*(sigmax(i,k)+sigmax(i,k+1))
            sigyij = 0.5d0*(sigmay(i,k)+sigmay(i,k+1))
            sigzij = 0.5d0*(sigmaz(i,k)+sigmaz(i,k+1))
            Gamsig = Gamxij*sigxij**2 &
                    +Gamyij*sigyij**2 &
                    +Gamzij*sigzij**2
!           _______________________________
!           Coefficients
!           -------------------------------
!           Matrix coeff. cell-center i
            AA0 = (Gamxij*sigxij*sumaix &
                  +Gamyij*sigyij*sumaiy &
                  +Gamsig*sumaiz)*nzAreaFaceij &
                  +AA0        	
!           -------------------------------
!           Matrix coeff. cell-center j=T  
            AAT = (Gamxij*sigxij*sumajx &
                  +Gamyij*sigyij*sumajy &
                  +Gamsig*sumajz)*nzAreaFaceij
!           -------------------------------
!           Matrix coeff. vertices
	    Bv1T = ( Gamxij*sigxij*fv1Tx &
                    +Gamyij*sigyij*fv1Ty &
                    +Gamsig*fv1Tz)*nzAreaFaceij + Bv1T
	    Bv2T = ( Gamxij*sigxij*fv2Tx &
                    +Gamyij*sigyij*fv2Ty &
                    +Gamsig*fv2Tz)*nzAreaFaceij + Bv2T
	    Bv3T = ( Gamxij*sigxij*fv3Tx &
                    +Gamyij*sigyij*fv3Ty &
                    +Gamsig*fv3Tz)*nzAreaFaceij + Bv3T

!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                              END OF CONDITIONAL
!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            ENDIF 
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |              Vertical BOTTOM neighbor: (j=5)             |
!          |----------------------------------------------------------|
!          |__________________________________________________________|

!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                             BOUNDARY BOTTOM FACES
!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            IF (k.eq.2) THEN                         
               nnx =  0.0d0
               nny =  0.0d0
               nnz = -1.0d0
               x = xc(i)
               y = yc(i)
               z = 0.5d0*(sig(1)+sig(2))
               dfBdn = Neumanndfdn3D(x,y,z,nnx,nny,nnz)
               Nm(i,k) = Nm(i,k)+ dfBdn*Gamx(i,k)*AreaCell(i)

!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                              INSIDE TOP FACES
!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            ELSE

	    sumaix = 0.0d0
	    sumajy = 0.0d0
	    sumaiz = 0.0d0
	    sumajx = 0.0d0
	    sumaiy = 0.0d0
	    sumajz = 0.0d0
            fv1Bx  = 0.0d0
            fv1By  = 0.0d0
            fv1Bz  = 0.0d0
            fv2Bx  = 0.0d0
            fv2By  = 0.0d0
            fv2Bz  = 0.0d0
            fv3Bx  = 0.0d0
            fv3By  = 0.0d0
            fv3Bz  = 0.0d0
!           _______________________________
!          |                               |
!          |      Volume of the region     |
!          |_______________________________|
            
            VoluReg = (1.0d0/3.0d0)*AreaCell(i)*(dsig(k-1)) 

!           _______________________________
!          |                               |
!          |         Region faces          |
!          |_______________________________|

            do s=1,3
!              -------------------------------
!              Vertex index               
	       ss=s+1
	       if (ss.gt.3) ss=ss-3
	       jv1 = No_vp(i,s)
	       jv2 = No_vp(i,ss) 
!              _________________________________________________
!              Region face 1,2,3
!              -------------------------------
!              Displacement vectors
               a1 = xv(jv2)-xc(i)
               a2 = yv(jv2)-yc(i)
               a3 = sigv(k-1)-sig(k-1)
               b1 = xv(jv1)-xc(i)
               b2 = yv(jv1)-yc(i)
               b3 = sigv(k-1)-sig(k-1)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              -------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Cell neighborn j contribution
               sumajx = sumajx + cxx
               sumajy = sumajy + cyy
               sumajz = sumajz + czz
!              ------------------------------
!              Vertex contribution   
               if (s.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
               elseif (s.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz  
               elseif (s.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz                
               endif    
!              _________________________________________________
!              Region faces 4,5,6
!              -------------------------------
!              Displacement vectors
               a1 = xv(jv1)-xc(i)
               a2 = yv(jv1)-yc(i)
               a3 = sigv(k-1)-sig(k)
               b1 = xv(jv2)-xc(i)
               b2 = yv(jv2)-yc(i)
               b3 = sigv(k-1)-sig(k)
!              -------------------------------
!              normal*Area
               nxAreaoVol = 0.5d0*(a2*b3-a3*b2)/VoluReg
               nyAreaoVol = 0.5d0*(a3*b1-a1*b3)/VoluReg
               nzAreaoVol = 0.5d0*(a1*b2-a2*b1)/VoluReg
!              ------------------------------
!              Coefficient
               cxx = nxAreaoVol/3.0d0
               cyy = nyAreaoVol/3.0d0
               czz = nzAreaoVol/3.0d0
!              ------------------------------
!              Cell-center i contribution
               sumaix = sumaix + cxx
               sumaiy = sumaiy + cyy
               sumaiz = sumaiz + czz
!              ------------------------------
!              Vertex contribution
               if (s.eq.1) then
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
               elseif (s.eq.2) then
	          fv2Bx = fv2Bx + cxx
	          fv2By = fv2By + cyy
	          fv2Bz = fv2Bz + czz
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz 
               elseif (s.eq.3) then
	          fv3Bx = fv3Bx + cxx
	          fv3By = fv3By + cyy
	          fv3Bz = fv3Bz + czz
	          fv1Bx = fv1Bx + cxx
	          fv1By = fv1By + cyy
	          fv1Bz = fv1Bz + czz                
               endif
            enddo
!           _______________________________
!          |                               |
!          |         Contributions         |
!          |_______________________________|
!           _______________________________
!           normal*Area (interface)
            nzAreaFaceij = -areaCell(i)
!           _______________________________
!           Gamma (interface)  
            Gamxij = 0.5d0*(Gamx(i,k)+Gamx(i,k-1))
            Gamyij = 0.5d0*(Gamy(i,k)+Gamy(i,k-1))
            Gamzij = 0.5d0*(Gamz(i,k)+Gamz(i,k-1))
!           _______________________________
!           sigma (interface)
            sigxij = 0.5d0*(sigmax(i,k)+sigmax(i,k-1))
            sigyij = 0.5d0*(sigmay(i,k)+sigmay(i,k-1))
            sigzij = 0.5d0*(sigmaz(i,k)+sigmaz(i,k-1))
            Gamsig = Gamxij*sigxij**2 &
                    +Gamyij*sigyij**2 &
                    +Gamzij*sigzij**2
!           _______________________________
!           Coefficients
!           -------------------------------
!           Matrix coeff. cell-center i
            AA0 = (Gamxij*sigxij*sumaix &
                  +Gamyij*sigyij*sumaiy &
                  +Gamsig*sumaiz)*nzAreaFaceij + AA0        	
!           -------------------------------
!           Matrix coeff. cell-center j=B  
            AAB = (Gamxij*sigxij*sumajx &
                  +Gamyij*sigyij*sumajy &
                  +Gamsig*sumajz)*nzAreaFaceij
!           -------------------------------
!           Matrix coeff. vertices
	    Bv1B = ( Gamxij*sigxij*fv1Tx &
                    +Gamyij*sigyij*fv1Ty &
                    +Gamsig*fv1Tz)*nzAreaFaceij + Bv1B
	    Bv2B = ( Gamxij*sigxij*fv2Tx &
                    +Gamyij*sigyij*fv2Ty &
                    +Gamsig*fv2Tz)*nzAreaFaceij + Bv2B
	    Bv3B = ( Gamxij*sigxij*fv3Tx &
                    +Gamyij*sigyij*fv3Ty &
                    +Gamsig*fv3Tz)*nzAreaFaceij + Bv3B

!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                              END OF CONDITIONAL
!           xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            ENDIF 
!           __________________________________________________________
!          |                                                          |
!          |----------------------------------------------------------|
!          |              Final coefficients contributions            |
!          |----------------------------------------------------------|
!          |__________________________________________________________|

!           --------------------------------------
!           Matrix cell-center contribution
	    Am0(i,k) = Am0(i,k) + AA0    
            Am1(i,k) = Am1(i,k) + AA(1)   
            Am2(i,k) = Am2(i,k) + AA(2)   
            Am3(i,k) = Am3(i,k) + AA(3)
            AmT(i,k) = AmT(i,k) + AAT             
            AmB(i,k) = AmB(i,k) + AAB
!           --------------------------------------
!           Matrix vextex contribution
	    Bmv1T(i,k) = Bmv1T(i,k) + Bv1T    
            Bmv2T(i,k) = Bmv2T(i,k) + Bv2T   
            Bmv3T(i,k) = Bmv3T(i,k) + Bv3T   
            Bmv1B(i,k) = Bmv1B(i,k) + Bv1B
            Bmv2B(i,k) = Bmv2B(i,k) + Bv2B             
            Bmv3B(i,k) = Bmv3B(i,k) + Bv3B  

         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: Diffusion3D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF DIFFUSION 3D                        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
