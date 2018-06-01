!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               --   SOLUTION OF SYSTEM BY GMRES  --                  !
!                            April 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE gmres2D(phi2D,phiv2D,        &
                         rhs2D,Gamx2D,Gamy2D, & 
                         xc,yc,No_cp,nbe,     &
                         xv,yv,No_vp,nbev)   

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the GMRES techni-   !
!    que. This program is just a modification of the one in nsmp 2D   !
!    program modified by Tam Nguyen. This system is adecuate for the  !
!    Poisson problem because it consider the cell-center "Am" and the !
!    vertex "Bmv" matrix coefficients.                                !
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
!   Input variables:                                                  !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> xg    | (N_CELL)   | Solution & initial guess            |  !
!  | ---> rhs   | (N_CELL)   | Right-hand side of the system       |  ! 
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
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     | (N_CELL0)  | matrix coefficient of element i     |  !
!  |    Am1     | (N_CELL0)  | matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     | (N_CELL0)  | matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     | (N_CELL0)  | matrix coeff. horizontal neigborn 3 |  ! 
!  |____________|____________|_____________________________________|  !
!  |    Bmv1    | (N_CELL0)  | matrix coeff. vertex 1 top          |  ! 
!  |    Bmv2    | (N_CELL0)  | matrix coeff. vertex 2 top          |  ! 
!  |    Bmv3    | (N_CELL0)  | matrix coeff. vertex 3 top          |  ! 
!  |____________|____________|_____________________________________|  !  
!  |    Newrhs  | (N_CELL0)  | right hand side of the method       |  !
!  |____________|____________|_____________________________________|  !  
!  |    nkr     | Dimension of the Hessenberg matrix               |  !
!  |    Ax      | (N_CELL)  Matrix-vector multiplication           |  !  
!  |____________|__________________________________________________|  !
!                                                                     !
!    Subroutines & functions used:                                    !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - diffusion2D                ( diffusion2D.F90     )        |  !
!  |   - interpolation2D            ( interpolation2D.F90 )        |  !
!  |   - BCcellcenter2D             ( BCcellcenter2D.F90  )        |  !
!  |   - BCvertex2D                 ( BCvertex2D.F90      )        |  !
!  |   - pmat2D                     ( gmres2D.F90         )        |  !
!  |   - soln2D                     ( gmres2D.F90         )        |  !
!  |   - prods2D                    ( gmres2D.F90         )        |  !
!  |_______________________________________________________________|  !
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

      real*8,dimension(:)   :: phi2D(N_CELL)
      real*8,dimension(:)   :: phiv2D(N_VERT)
!     -------------------------------------
      real*8,dimension(:)   :: rhs2D(N_CELL)
      real*8,dimension(:)   :: Gamx2D(N_CELL)
      real*8,dimension(:)   :: Gamy2D(N_CELL)
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

      real*8,dimension(:)   :: Am0(N_CELL0)
      real*8,dimension(:)   :: Am1(N_CELL0)
      real*8,dimension(:)   :: Am2(N_CELL0)
      real*8,dimension(:)   :: Am3(N_CELL0) 
      real*8,dimension(:)   :: Bmv1(N_CELL0)
      real*8,dimension(:)   :: Bmv2(N_CELL0)
      real*8,dimension(:)   :: Bmv3(N_CELL0)
!     ------------------------------------
      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: h
      real*8,dimension(:),  allocatable :: r
      real*8,dimension(:),  allocatable :: xg
      real*8,dimension(:),  allocatable :: aux,aux1,vecv
!     ------------------------------------
      data zero/1.d-8/
      integer::il,ip,it,nk
      real*8::zero,rec
      real*8::tem,dem,res,prods2D
      real*8::gamma,sing,cosg,raux,hij,hii,hipj,hipi
!     ------------------------------------
      integer :: ns 
      integer :: nkr
                 ns = N_CELL0
                 nkr= DimHess
!     ------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: GMRES 2D'
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
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!     _____________________________________
!     Vertex interpolation and vertex BC

      call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                           phi2D,xc,yc,No_cp,nbe)

      call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                      phi2D,xc,yc,No_cp,nbe)

      do i=1,N_CELL0	
         xg(i) = phi2D(i)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |            Matrix Am & Bm of the diffusion term        |
!     |________________________________________________________|
        
      do i=1,N_CELL0 	
         Am0(i)  = 0.0d0 
         Am1(i)  = 0.0d0
         Am2(i)  = 0.0d0
         Am3(i)  = 0.0d0
         Bmv1(i) = 0.0d0 
         Bmv2(i) = 0.0d0
         Bmv3(i) = 0.0d0
      enddo

      call diffusion2D(Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                       Gamx2D,Gamy2D,                  &
                       xc,yc,No_cp,nbe,                &
                       xv,yv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!           Residual for the initial guess (no preconditioner)        !
!                                                                     !
!*********************************************************************!

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg      
      call pmat2D(aux1,                              &
                  xg,Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                  xc,yc,No_cp,nbe,                   &
                  xv,yv,No_vp,nbev)  
!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs      
      do i=1,ns
        aux(i) = rhs2D(i)- aux1(i) 
      enddo
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)	
      res = prods2D(ns,aux,aux)
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
         call pmat2D(aux1,                               &
                     aux,Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                     xc,yc,No_cp,nbe,                    &
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
            tem= prods2D(ns,aux1,vecv)
            h(i,j) = tem
            do il=1,ns
                aux(il)=aux(il)-tem*v(il,i)
            enddo
         enddo
!        ------------------------------------------- 
!        H(j+1,j)
         dem= prods2D(ns,aux,aux)
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

      call soln2D(nk,ns,r,h,v,xg)

!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |  
!     |________________________________________________________|

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg      
      call pmat2D(aux1,                              &
                  xg,Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                  xc,yc,No_cp,nbe,                   &
                  xv,yv,No_vp,nbev) 

!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs        
      do i=1,ns
        aux(i) = rhs2D(i)- aux1(i) 
      enddo
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)
      res = prods2D(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
7     format(t10,a26,i5,a9,e10.3)
6     format(t10,a12,i4,a9,e10.3)

      if (res.lt.rec) then
         write(*,*) ' '
         write(*,7) 'Sol GMRES 2D  : iters =',it,', error =',res
         write(*,*) ' '
         goto 4
      elseif (res.gt.1.0d5) then
         write(*,7) 'DIVERGENCE !!! : iters =',it,', error =',res
         stop
      elseif (it.ge.MaxIters) then
         write(*,7) 'Non-convergence: iters =',it,', error =',res
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
      do i=1,N_CELL0	
         phi2D(i) = xg(i)
      enddo
!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)
!     ________________________________________________________
!     Vertex values                      

!     ------------------------------------------------
!     Vertex interpolation
      call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                           phi2D,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                      phi2D,xc,yc,No_cp,nbe)

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|
!
	deallocate(v,h,r,xg,aux,aux1,vecv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: GMRES2D'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      END SUBROUTINE gmres2D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	      END OF GMRES 2D                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



!                      !wwwwwwwwwwwwwwwwwwwwwwww!
!                      !------------------------!
!                      !  Auxiliar subroutines  !
!                      !------------------------!
!                      !wwwwwwwwwwwwwwwwwwwwwwww!



!=====================================================================!
!---------------------------------------------------------------------!
!          AUXILIAR 1: Solution of the upper triangular system        !
!---------------------------------------------------------------------!
!=====================================================================!

      subroutine soln2D(nk,ns,r,hg,v,xg)

!---------------------------------------------------------------------!
!                                                                     !
!     Ce sous-programme resout le systeme triangulaire superieur      !
!     permettant de calculer la solution du probleme de minimisa-     !
!     tion dans l"algorithme gmres. La solution approchee a l'etape   !
!     nk est calculee et stokee dans le vecteur x.                    !
!                                                                     !
!---------------------------------------------------------------------!

!     !**********************************************************!
!     !                    Definitions                           !
!     !**********************************************************!
!     ____________________________________
!     Keys and common parameters 

      implicit none
#     include "common.mpf"
!     ____________________________________
!     Declaration of variables  

      integer:: nk
      integer:: ns
      real*8,dimension(:)   :: r(nk)
      real*8,dimension(:,:) :: hg(nk,nk)
      real*8,dimension(:,:) :: v(ns,nk)
      real*8,dimension(:)   :: xg(ns)
!     ____________________________________
!     Declaration of local variables

      integer:: il
      real*8 :: tem

!     !**********************************************************!
!     !                        Solution                          !
!     !**********************************************************!

      do i=nk,1,-1
         tem = r(i)/hg(i,i)
         r(i)=tem
         do j=i-1,1,-1
            r(j)=r(j)-hg(j,i)*tem
         enddo
      enddo    

      do i=1,nk
         tem = r(i)
         do il=1,ns
            xg(il)=xg(il)+tem*v(il,i)
         enddo
      enddo

      end subroutine

!=====================================================================!
!---------------------------------------------------------------------!
!         AUXILIAR 2: function of the product of two vectors          !
!---------------------------------------------------------------------!
!=====================================================================!

      function prods2D(n,x,y)    

!     !**********************************************************!
!     !                    Definitions                           !
!     !**********************************************************!

      implicit none
!     ____________________________________
!     Declaration of variables  

      integer :: n
      real*8,dimension(:) :: x(n)
      real*8,dimension(:) :: y(n)
!     ____________________________________
!     Declaration of local variables

      real*8::som,prods2D
      integer :: i

!     !**********************************************************!
!     !                        Product: x*y                      !
!     !**********************************************************!

      som=0.0d0
      do i=1,n
         som = som + x(i)*y(i)    
      enddo
      prods2D = som

      return
      end function prods2D

!=====================================================================!
!---------------------------------------------------------------------!
!               AUXILIAR 3:  Construction of vector A*x               !
!---------------------------------------------------------------------!
!=====================================================================!

      subroutine pmat2D(Axg,                               &
                        xg,Am0,Am1,Am2,Am3,Bmv1,Bmv2,Bmv3, &
                        xc,yc,No_cp,nbe,                   &
                        xv,yv,No_vp,nbev)    

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine calculates the resultant vector of the multi-    !
!    plication of the matrix A with the vector xg.                    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <-- Axg    | (N_CELL)   | Matrix-vector multiplication        |  !
!  | --> xg     | (N_CELL)   | Guess function vector               |  !
!  | --> rhs    | (N_CELL)   | right-hand side of the system       |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

!     !**********************************************************!
!     !                    Definitions                           !
!     !**********************************************************!
!     ____________________________________
!     Keys and common parameters 

      implicit none
#     include "common.mpf"
!     ____________________________________
!     Declaration of variables  

      real*8,dimension(:)   :: Axg(N_CELL0)
!     ------------------------------------
      real*8,dimension(:)   :: xg(N_CELL0)
!     ------------------------------------
      real*8,dimension(:)   :: Am0(N_CELL0)
      real*8,dimension(:)   :: Am1(N_CELL0)
      real*8,dimension(:)   :: Am2(N_CELL0)
      real*8,dimension(:)   :: Am3(N_CELL0) 
      real*8,dimension(:)   :: Bmv1(N_CELL0)
      real*8,dimension(:)   :: Bmv2(N_CELL0)
      real*8,dimension(:)   :: Bmv3(N_CELL0)
!     ------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     ____________________________________
!     Declaration of local variables

      real*8,dimension(:) :: phi2D(N_CELL)
      real*8,dimension(:) :: phiv2D(N_VERT)
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
         phi2D(i) = xg(i)
      enddo
!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!     ________________________________________________________
!    |                                                        |
!    |                     Vertex values                      |
!    |________________________________________________________|

!     ------------------------------------------------
!     Vertex interpolation
      call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                           phi2D,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                      phi2D,xc,yc,No_cp,nbe)

!      ________________________________________________________
!     |                                                        |
!     |          Multiplication Ax = A*phi + B*phiv            |
!     |________________________________________________________|

      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
	 Axg(i) =  Am0(i)*phi2D(i)     &
                 + Am1(i)*phi2D(jc1)   &
                 + Am2(i)*phi2D(jc2)   &
                 + Am3(i)*phi2D(jc3)   &
                 + Bmv1(i)*phiv2D(jv1) &
                 + Bmv2(i)*phiv2D(jv2) &
                 + Bmv3(i)*phiv2D(jv3)
      enddo

      return
      end subroutine

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        END OF GMRES SUBROUTINES                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
