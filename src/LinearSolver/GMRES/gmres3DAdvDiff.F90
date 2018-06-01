!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               --   SOLUTION OF SYSTEM BY GMRES  --                  !
!                            April 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE gmres3DAdvDiff(phi,phiv,                      &
                                rhs,uu,vv,ww,Gamx,Gamy,Gamz,   &
                                xc,yc,sig,dsig,No_cp,nbe,      &
                                xv,yv,sigv,dsigv,No_vp,nbev,   &
                                Hpr,h,etan,                    &
                                Hprv,hv,etav)   

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the GMRES techni-   !
!    que. This program is just a simple modification of the one in    !
!    nsmp program in 2D modified by Tam Nguyen.                       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi   |(N_CELL,NZ)  | Solution & initial guess           |  !
!  | <--> phiv  |(N_VERT,NZ-1)| Solution at the vertex values      |  !
!  |____________|_____________|____________________________________|  !
!  | ---> rhs   |(N_CELL,NZ)  | Right-hand side of the system      |  !
!  | ---> Gamx  |(N_CELL,NZ)  | Diffusive coefficient in x         |  ! 
!  | ---> Gamy  |(N_CELL,NZ)  | Diffusive coefficient in y         |  ! 
!  | ---> Gamz  |(N_CELL,NZ)  | Diffusive coefficient in z         |  ! 
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> sig   | (NZ)        | sigma value at the cell centers    |  !
!  | ---> dsig  | (NZ)        | = sig(k+1)-sig(k+1)                |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv | (N_VERT)    | Coordinates of the vertices        |  !
!  | ---> sigv  | (NZ-1)      | sigma value at the vertices        |  !
!  | ---> dsigv | (NZ-1)      | = sigv(k+1)-sigv(k)                |  !
!  | ---> No_vp | (N_CELL0,3) | Numbering of the 3 cell vertices   |  !
!  | ---> nbev  | (N_VERT)    | Tag type of cell vertex            |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!   _______________________________________________________________   !
!  |     Name   |               Description                        |  !  
!  |____________|__________________________________________________|  !
!  |  nkr       | Dimension of the Hessenberg matrix               |  !
!  |  Ax        | (N_CELL)  Matrix-vector multiplication           |  !
!  |  ind       | Tag related to the gmres method                  |  ! 
!  |____________|__________________________________________________|  !

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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

!     --------------------------------------
      real*8,dimension(:),  allocatable :: xg
      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: hh
      real*8,dimension(:),  allocatable :: r
      real*8,dimension(:),  allocatable :: aux,aux1,vecv
!     -----------------------------------------
      data zero/1.d-8/
      integer::ii,il,ind,ip,it,nk
      real*8::zero,rec,maxi
      real*8::tem,dem,res,prods
      real*8::gamma,sing,cosg,raux,hij,hii,hipj,hipi
!     ----------------------------------------
      integer :: ns 
      integer :: nkr
                 ns = N_CELL0*(NZ-2)
                 nkr= DimHess
!     -----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: GMRES 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

	allocate(v(ns,nkr),hh(nkr,nkr),r(nkr))
	allocate(xg(ns),aux(ns),aux1(ns),vecv(ns))

!      ________________________________________________________
!     |                                                        |
!     |                     Initial guess                      |
!     |________________________________________________________|

!     ______________________________________________________
!     Cell-centers BC   

!     ----------------------------------------------
!     Boundary Condition of the cell-center points  
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
!     -----------------------------------------------                      
!     xg = phi 
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            xg(ii) = phi(i,k)
         enddo
      enddo 

!*********************************************************************!
!                                                                     !
!           Residual for the initial guess (no preconditioner)        !
!                                                                     !
!*********************************************************************!

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg   
      call pmatAdvDiff3D(aux1,xg,ns,                   &
                         uu,vv,ww,Gamx,Gamy,Gamz,      &
                         xc,yc,sig,dsig,No_cp,nbe,     &
                         xv,yv,sigv,dsigv,No_vp,nbev,  &
                         Hpr,h,etan,                   &
                         Hprv,hv,etav)  

!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs      
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            aux(ii)= rhs(i,k)-aux1(ii)
         enddo
      enddo 
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)	
      res = prods(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
      rec = zero 
      if (res.lt.rec) then
         write(*,*) ' '
         write(*,7) '* No GMRES iterations * ',0,', error =',res
         write(*,*) ' '
         goto 4
      endif
      rec = max(rec,res*eps)

!*********************************************************************!
!                                                                     !
!                             GMRES iterations                        !
!                                                                     !
!*********************************************************************!

      it = 0   
1     continue
      it = it + 1 
!      ________________________________________________________
!     |                                                        |
!     |               Orthogonalisation d'arnoldi              |  
!     |________________________________________________________|

      nk=nkr

      do i=1,ns
         aux(i) = aux(i)/res
      enddo

      DO j=1,nkr
!        ------------------------------------------- 
!        Vector v(:,j)
         do il=1,ns 
            v(il,j) = aux(il) 
         enddo
!        ------------------------------------------- 
!        Product: aux1 = A*aux
         call pmatAdvDiff3D(aux1,aux,ns,                 &
                            uu,vv,ww,Gamx,Gamy,Gamz,     &
                            xc,yc,sig,dsig,No_cp,nbe,    &
                            xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,etan,                  &
                            Hprv,hv,etav) 
         do i=1,ns
            aux(i)=aux1(i)
         enddo
!        ------------------------------------------- 
!        H(1:j,j)	  
         do i=1,j
	    do il=1,ns
	       vecv(il)=v(il,i)
	    enddo	
            tem= prods(ns,aux1,vecv)
            hh(i,j) = tem
            do il=1,ns
                aux(il)=aux(il)-tem*v(il,i)
            enddo
         enddo
!        ------------------------------------------- 
!        H(j+1,j)
         dem= prods(ns,aux,aux)
         dem= dsqrt(dem)
         hh(j+1,j) = dem
!        ------------------------------------------- 
!        Criteria: dimension of H
         if (dem.lt.rec) then
             nk=j
             goto 5
         endif
!        ------------------------------------------- 
!        Vector v(:,j+1)
         do i=1,ns
            aux(i) = aux(i)/dem
         enddo
      ENDDO
!      ________________________________________________________
!     |                                                        |
!     |        Triangularisation et modif. du second membre    |  
!     |________________________________________________________|
  
 5    continue

      r(1) = res
      do i=2,nk
         r(i)=0.0d0
      enddo

      do i=1,nk
         ip=i+1
         hii =hh(i,i)
         hipi=hh(ip,i) 
         gamma= dsqrt(hii*hii+hipi*hipi)
         cosg= hii/gamma
         sing=-hipi/gamma
         do j=i,nk
            hij=hh(i,j)
            hipj=hh(ip,j)
            hh(i,j) = cosg*hij - sing*hipj
            hh(ip,j)= sing*hij + cosg*hipj
         enddo
         raux=r(i)
         r(i) =cosg*raux
         r(ip)=sing*raux
      enddo

!      ________________________________________________________
!     |                                                        |
!     |            Solution of the tridiagonal system          |  
!     |________________________________________________________|

      call soln(nk,ns,r,hh,v,xg)

!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |  
!     |________________________________________________________|

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg   
      call pmatAdvDiff3D(aux1,xg,ns,                   &
                         uu,vv,ww,Gamx,Gamy,Gamz,      &
                         xc,yc,sig,dsig,No_cp,nbe,     &
                         xv,yv,sigv,dsigv,No_vp,nbev,  &
                         Hpr,h,etan,                   &
                         Hprv,hv,etav) 
!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs  
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            aux(ii)= rhs(i,k)-aux1(ii)
         enddo
      enddo            
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)
      res = prods(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria

7     format(t10,a26,i5,a9,e10.3)
6     format(t10,a12,i4,a9,e12.6)

      if (res.lt.eps) then
         write(*,*) ' '
         write(*,7) 'Solution GMRES 3D: iters =',it,&
                    ', error =',res     
         write(*,*) ' '
         goto 4
      elseif (it.ge.MaxIters) then
         write(*,*) ' '
         write(*,7) 'Non-convergence: iters =',it,&
                    ', error =',res
         write(*,*) ' '
         goto 4
      else
         !write(*,6)'iter GMRES = ',it,', error= ',res  
         goto 1
      endif

 4    continue

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                     Final solution                     |
!     |________________________________________________________|

!     ________________________________________________________
!     Cell-center values 
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            phi(i,k) = xg(ii)
         enddo
      enddo 
!     ----------------------------------------------
!     Boundary Condition of the cell-center points  
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
!     ______________________________________________________
!     Vertex interpolation & BC
      
!     ------------------------------------------------
!     Interpolation of the inside vertex points 
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav, &
                      phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|
!
	deallocate(v,hh,r,xg,aux,aux1,vecv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: GMRES 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END 

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	      END OF GMRES                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



!                      !wwwwwwwwwwwwwwwwwwwwwwww!
!                      !------------------------!
!                      !   Auxiliar subroutine  !
!                      !------------------------!
!                      !wwwwwwwwwwwwwwwwwwwwwwww!



!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               AUXILIAR 1:  Construction of vector A*x               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE pmatAdvDiff3D(Axg,xg,ns,                           &
                               uu,vv,ww,Gamx,Gamy,Gamz,             &
                               xc,yc,sig,dsig,No_cp,nbe,            &
                               xv,yv,sigv,dsigv,No_vp,nbev,         &
                               Hpr,h,etan,                          &
                               Hprv,hv,etav)
!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine calculates the resultant vector of the multi-    !
!    plication of the matrix A with the vector x.                     !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <-- Ax     | (ns)       | Matrix-vector multiplication        |  !
!  | --> xg     | (ns)       | Guess function vector               |  !
!  | --> ns     | integer    | size of teh vector ns = N_CELL*NZ   |  !  
!  | --> rhs    | (N_CELL)   | right-hand side of the system       |  !
!  | --> ind    | integer    | Tag related to the gmres method     |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!     !**********************************************************!
!     !                    Definitions                           !
!     !**********************************************************!
!     ____________________________________
!     Keys and common parameters 

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
!     ____________________________________
!     Declaration of variables  

      integer :: ns
      real*8,dimension(:)   :: Axg(ns)
      real*8,dimension(:)   :: xg(ns)
!     -------------------------------------
      real*8,dimension(:,:) :: uu(N_CELL,NZ)
      real*8,dimension(:,:) :: vv(N_CELL,NZ)
      real*8,dimension(:,:) :: ww(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     ____________________________________
!     Declaration of local variables

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     -------------------------------------
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmG(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
!     -------------------------------------
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
      integer:: ii
      real*8 :: Vol

!     !**********************************************************!
!     !                     Initialization                       !
!     !**********************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                   Cell-center values                   |
!     |________________________________________________________|

      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            phi(i,k) = xg(ii)
         enddo
      enddo 
!     ----------------------------------------------
!     Boundary Condition of the cell-center points  
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
!      ________________________________________________________
!     |                                                        |
!     |                     Vertex values                      |
!     |________________________________________________________|
      
!     ------------------------------------------------
!     Interpolation of the inside vertex points 
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav, &
                      phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)

!     ________________________________________________________
!    |                                                        |
!    |                   Matrix contributions                 |
!    |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0	
            Am0(i,k)  = 0.0d0
            Am1(i,k)  = 0.0d0
            Am2(i,k)  = 0.0d0
            Am3(i,k)  = 0.0d0
            AmT(i,k)  = 0.0d0
            AmB(i,k)  = 0.0d0
!           -----------------
            AmG(i,k)  = 0.0d0 
!           -----------------
            Bmv1T(i,k)= 0.0d0 
            Bmv2T(i,k)= 0.0d0 
            Bmv3T(i,k)= 0.0d0 
            Bmv1B(i,k)= 0.0d0 
            Bmv2B(i,k)= 0.0d0 
            Bmv3B(i,k)= 0.0d0 
         enddo
      enddo
!     ________________________________________________________
!     Advection contribution                 

      call advection3D(Am0,Am1,Am2,Am3,AmT,AmB,AmG,   &
                       uu,vv,ww,                      &
                       phi,xc,yc,sig,No_cp,nbe,       &  
                       sigv,dsigv)

!     ________________________________________________________
!     Diffusion contribution with a negative sign

      do k=1,NZ
         do i=1,N_CELL	
            Gamx(i,k) = -Gamx(i,k)
            Gamy(i,k) = -Gamy(i,k)
            Gamz(i,k) = -Gamz(i,k)
         enddo
      enddo

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             &
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)  
!      ________________________________________________________
!     |                                                        |
!     |      Multiplication A*xg = phi + dt*(ADV - DIFF)       |
!     |    (The negative sign is already included in Bmv)      |
!     |________________________________________________________|

      do k=2,NZ-1
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            Vol = areaCell(i)*dsigv(k-1)
	    Axg(ii) = phi(i,k)+ dt/Vol*( Am0(i,k)*phi(i,k)        &
                                       + Am1(i,k)*phi(jc1,k)      &
                                       + Am2(i,k)*phi(jc2,k)      &
                                       + Am3(i,k)*phi(jc3,k)      &
                                       + AmT(i,k)*phi(i,k+1)      &       
                                       + AmB(i,k)*phi(i,k-1)      &
!                                      ----------------------------
                                       + AmG(i,k)                 &
!                                      ----------------------------
                                       + Bmv1T(i,k)*phiv(jv1,k)   &
                                       + Bmv2T(i,k)*phiv(jv2,k)   &
                                       + Bmv3T(i,k)*phiv(jv3,k)   &
                                       + Bmv1B(i,k)*phiv(jv1,k-1) &
                                       + Bmv2B(i,k)*phiv(jv2,k-1) &
                                       + Bmv3B(i,k)*phiv(jv3,k-1))  
         enddo
      enddo

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      END OF AUXILIAR SUBROUTINES                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
