!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               --   SOLUTION OF SYSTEM BY GMRES  --                  !
!                            April 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE gmres2DAdvDiff(phi,phiv,            &
                                rhs,uu,vv,Gamx,Gamy, & 
                                xc,yc,No_cp,nbe,     &
                                xv,yv,No_vp,nbev)   

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the GMRES techni-   !
!    que. This program is just a simple modification of the one in    !
!    nsmp program in 2D modified by Tam Nguyen.                       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Outpu & Input variables:                                         !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> phi   | (N_CELL)   | Solution & initial guess cell-center|  !
!  | <--> phiv  | (N_CELL)   | Solution & initial guess vertex     |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> phi   | (N_CELL)   | Solution & initial guess cell-center|  !
!  | <--> phiv  | (N_CELL)   | Solution & initial guess vertex     |  !
!  | ---> rhs   | (N_CELL)   | Right-hand side of the system       |  ! 
!  | ---> uu    | (N_CELL)   | Component of the velocity profile u |  !
!  | ---> vv    | (N_CELL)   | Component of the velocity profile v |  !
!  | ---> Gamx  | (N_CELL)   | Diffusive coefficient x             |  !
!  | ---> Gamy  | (N_CELL)   | Diffusive coefficient y             |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc | (N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> sig   | (NZ)       | sigma value at the cell centers     |  !
!  | ---> dsig  | (NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | ---> No_cp | (N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe   | (N_CELL0)  | Tag type cell-center                |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xv,yv | (N_VERT)   | Coordinates of the vertices         |  !
!  | ---> sigv  | (NZ-1)     | sigma value at the vertices         |  !
!  | ---> dsigv | (NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | ---> No_vp | (N_CELL0,3)| Numbering of the 3 cell vertices    |  !
!  | ---> nbev  | (N_VERT)   | Tag type of cell vertex             |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!   _______________________________________________________________   !
!  |     Name   |               Description                        |  !  
!  |____________|__________________________________________________|  !
!  |  nkr       | Dimension of the Hessenberg matrix               |  !
!  |  Ax        | (N_CELL)  Matrix-vector multiplication           |  !
!  |____________|__________________________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

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
!    |                                    |
!    |      Declaration of variables      |
!    |____________________________________|

      real*8,dimension(:)   :: phi(N_CELL)
      real*8,dimension(:)   :: phiv(N_VERT)
!     -------------------------------------
      real*8,dimension(:)   :: rhs(N_CELL)
      real*8,dimension(:)   :: uu(N_CELL)
      real*8,dimension(:)   :: vv(N_CELL)
      real*8,dimension(:)   :: Gamx(N_CELL)
      real*8,dimension(:)   :: Gamy(N_CELL)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: h
      real*8,dimension(:),  allocatable :: r
      real*8,dimension(:),  allocatable :: xg
      real*8,dimension(:),  allocatable :: aux,aux1,vecv
!     ------------------------------------
      data zero/1.d-8/
      integer::il,ip,it,nk
      real*8::zero,rec
      real*8::tem,dem,res,prods
      real*8::gamma,sing,cosg,raux,hij,hii,hipj,hipi
!     ------------------------------------
      integer :: ns 
      integer :: nkr
                 ns = N_CELL0
                 nkr= DimHess
!     ------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: GMRESAdvDiff2D'
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

      allocate(v(ns,nkr),h(nkr,nkr),r(nkr))
      allocate(xg(ns),aux(ns),aux1(ns),vecv(ns))

!      ________________________________________________________
!     |                                                        |
!     |                       Initial guess                    |
!     |________________________________________________________|

!     ____________________________________
!     Cell-center BC
      call BCcellcenter2D(phi,xc,yc,No_cp,nbe)
!     _____________________________________
!     Initial guess assignation xg = phi
      do i=1,N_CELL0	
         xg(i) = phi(i)
      enddo

!*********************************************************************!
!                                                                     !
!           Residual for the initial guess (no preconditioner)        !
!                                                                     !
!*********************************************************************!

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg      
      call pmatAdvDiff2D(aux1,xg,                   &
                         uu,vv,Gamx,Gamy,           &
                         xc,yc,No_cp,nbe,           &
                         xv,yv,No_vp,nbev)  
!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs      
      do i=1,ns
        aux(i) = rhs(i)- aux1(i) 
      enddo
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)	
      res = prods(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
      rec = zero
      if (res.lt.rec) then
         write(*,7) '* No GMRES iterations * ',0,', error =',res
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
         call pmatAdvDiff2D(aux1,aux,              &
                            uu,vv,Gamx,Gamy,       &
                            xc,yc,No_cp,nbe,       &
                            xv,yv,No_vp,nbev)  
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
            h(i,j) = tem
            do il=1,ns
                aux(il)=aux(il)-tem*v(il,i)
            enddo
         enddo
!        ------------------------------------------- 
!        H(j+1,j)
         dem= prods(ns,aux,aux)
         dem= dsqrt(dem)
         h(j+1,j) = dem
!        ------------------------------------------- 
!        Criteria: dimension of H
         if (dem .lt. rec) then
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
         hii =h(i,i)
         hipi=h(ip,i) 
         gamma= dsqrt(hii*hii+hipi*hipi)
         cosg= hii/gamma
         sing=-hipi/gamma
         do j=i,nk
            hij=h(i,j)
            hipj=h(ip,j)
            h(i,j) = cosg*hij - sing*hipj
            h(ip,j)= sing*hij + cosg*hipj
         enddo
         raux=r(i)
         r(i) =cosg*raux
         r(ip)=sing*raux
      enddo
!      ________________________________________________________
!     |                                                        |
!     |            Solution of the tridiagonal system          |  
!     |________________________________________________________|

      call soln(nk,ns,r,h,v,xg)

!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |  
!     |________________________________________________________|

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg      
      call pmatAdvDiff2D(aux1,xg,                   &
                         uu,vv,Gamx,Gamy,           &
                         xc,yc,No_cp,nbe,           &
                         xv,yv,No_vp,nbev) 
!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs        
      do i=1,ns
        aux(i) = rhs(i)- aux1(i) 
      enddo
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)
      res = prods(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
7     format(t10,a26,i5,a9,e10.3)
      if (res.lt.rec) then
         write(*,7) 'Sol GMRES Adv2D: iters =',it,', error =',res
         goto 4
      elseif (res.gt.1.0d5) then
         write(*,7) 'DIVERGENCE !!! : iters =',it,', error =',res
         stop
      elseif (it.ge.MaxIters) then
         write(*,7) 'Non-convergence: iters =',it,', error =',res
         goto 4
      else
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
      do i=1,N_CELL0	
         phi(i) = xg(i)
      enddo
!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi,xc,yc,No_cp,nbe)
!     ________________________________________________________
!     Vertex values                      

!     ------------------------------------------------
!     Vertex interpolation
      call interpolation2D(phiv,xv,yv,No_vp,nbev, &
                           phi,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex2D(phiv,xv,yv,No_vp,nbev, &
                      phi,xc,yc,No_cp,nbe)

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|
!
	deallocate(v,h,r,xg,aux,aux1,vecv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: GMRESAdvDiff2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END 

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	      END OF GMRES adv                        !
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

      SUBROUTINE pmatAdvDiff2D(Axg,xg,               &
                               uu,vv,Gamx,Gamy,      & 
                               xc,yc,No_cp,nbe,      &
                               xv,yv,No_vp,nbev)    

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
!  | <-- Ax     | (N_CELL)   | Matrix-vector multiplication        |  !
!  | --> xg     | (N_CELL)   | Guess function vector               |  !
!  | --> uu     | (N_CELL)   | Component of the velocity profile u |  !
!  | --> vv     | (N_CELL)   | Component of the velocity profile v |  !
!  | --> Gamx   | (N_CELL)   | Diffusive coefficient x             |  !
!  | --> Gamy   | (N_CELL)   | Diffusive coefficient y             |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!     !**********************************************************!
!     !                    Definitions                           !
!     !**********************************************************!

!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

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
!     |   Declaration of variables         |
!     |____________________________________|

      real*8,dimension(:)   :: Axg(N_CELL0)
!     ------------------------------------
      real*8,dimension(:)   :: xg(N_CELL0)
!     ------------------------------------
      real*8,dimension(:,:) :: uu(N_CELL)
      real*8,dimension(:,:) :: vv(N_CELL)
      real*8,dimension(:,:) :: Gamx(N_CELL)
      real*8,dimension(:,:) :: gamy(N_CELL)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 

!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8, dimension(:) :: phi(N_CELL)
      real*8, dimension(:) :: phiv(N_VERT)
!     --------------------------------------
      real*8, dimension(:) :: Am0(N_CELL0)
      real*8, dimension(:) :: Am1(N_CELL0)
      real*8, dimension(:) :: Am2(N_CELL0)
      real*8, dimension(:) :: Am3(N_CELL0)
      real*8, dimension(:) :: AmG(N_CELL0)
      real*8, dimension(:) :: Bmv1(N_CELL0)
      real*8, dimension(:) :: Bmv2(N_CELL0)
      real*8, dimension(:) :: Bmv3(N_CELL0)
      real*8 :: Vol
!     -------------------------------------
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3


!     !**********************************************************!
!     !             Construction of the vector A*x               !
!     !**********************************************************!

!     ________________________________________________________
!    |                                                        |
!    |                   Cell-center values                   |
!    |________________________________________________________|

!     ------------------------------------------------
!     Read the current inside values 
      do i=1,N_CELL0	
         phi(i) = xg(i)
      enddo
!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi,xc,yc,No_cp,nbe)

!     ________________________________________________________
!    |                                                        |
!    |                     Vertex values                      |
!    |________________________________________________________|

!     ------------------------------------------------
!     Vertex interpolation
      call interpolation2D(phiv,xv,yv,No_vp,nbev, &
                           phi,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex2D(phiv,xv,yv,No_vp,nbev, &
                      phi,xc,yc,No_cp,nbe)

!      ________________________________________________________
!     |                                                        |
!     |         Construction of the Matrix coefficients        |
!     |________________________________________________________|

      do i=1,N_CELL0	
         Am0(i)  = 0.0d0
         Am1(i)  = 0.0d0
         Am2(i)  = 0.0d0
         Am3(i)  = 0.0d0
!        ---------------
         AmG(i)  = 0.0d0  
!        ---------------
         Bmv1(i) = 0.0d0 
         Bmv2(i) = 0.0d0 
         Bmv3(i) = 0.0d0 
      enddo
!     ________________________________________________________
!     Advection contribution                 

      call advection2D(Am0,Am1,Am2,Am3,AmG, &
                       uu,vv,               &
                       phi,xc,yc,No_cp,nbe)  
!     ________________________________________________________
!     Diffusion contribution with a negative sign            

      do i=1,N_CELL	
         Gamx(i) = -Gamx(i)
         Gamy(i) = -Gamy(i)
      enddo

      call diffusion2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                       Gamx,Gamy,                      &
                       xc,yc,No_cp,nbe,                &
                       xv,yv,No_vp,nbev)
!      ________________________________________________________
!     |                                                        |
!     |      Multiplication A*xg = phi + dt*(ADV - DIFF)       |
!     |    (The negative sign is already included in Bmv)      |
!     |________________________________________________________|

      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
         Vol = areaCell(i)
	 Axg(i) = phi(i) + dt/Vol*( Am0(i)*phi(i)     &
                                  + Am1(i)*phi(jc1)   &
                                  + Am2(i)*phi(jc2)   &
                                  + Am3(i)*phi(jc3)   &
!                                 ---------------
                                  + AmG(i)            &
!                                 ---------------
                                  + Bmv1(i)*phiv(jv1) &
                                  + Bmv2(i)*phiv(jv2) &
                                  + Bmv3(i)*phiv(jv3)) 
      enddo

      return
      end

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      END OF AUXILIAR SUBROUTINES                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
