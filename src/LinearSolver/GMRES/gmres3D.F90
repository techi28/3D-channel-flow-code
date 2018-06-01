!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!          --   SOLUTION OF SYSTEM BY GMRES (POISSON)  --             !
!                            April 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE gmres3D(phi,phiv,                    &
                         rhs,Gamx,Gamy,Gamz,          &
                         xc,yc,sig,dsig,No_cp,nbe,    &
                         xv,yv,sigv,dsigv,No_vp,nbev, &
                         Hpr,h,etan,                  &
                         Hprv,hv,etav)  

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
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  |    Am1     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  |    AmT     |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  |    AmB     |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  !
!  |____________|____________|_____________________________________|  !
!  |    Bmv1T   |(N_CELL0,NZ)| matrix coeff. vertex 1 top          |  ! 
!  |    Bmv2T   |(N_CELL0,NZ)| matrix coeff. vertex 2 top          |  ! 
!  |    Bmv3T   |(N_CELL0,NZ)| matrix coeff. vertex 3 top          |  ! 
!  |    Bmv1B   |(N_CELL0,NZ)| matrix coeff. vertex 1 bottom       |  ! 
!  |    Bmv2B   |(N_CELL0,NZ)| matrix coeff. vertex 2 bottom       |  ! 
!  |    Bmv3B   |(N_CELL0,NZ)| matrix coeff. vertex 3 bottom       |  !
!  |____________|____________|_____________________________________|  !  
!  |    Newrhs  |(N_CELL0,NZ)| right hand side of the method       |  !
!  |____________|____________|_____________________________________|  !  
!  |    nkr     | Dimension of the Hessenberg matrix               |  !
!  |    Ax      | (N_CELL)  Matrix-vector multiplication           |  !  
!  |____________|__________________________________________________|  !
!                                                                     !
!    Subroutines & functions used:                                    !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - diffusion3D                ( diffusion3D.F90     )        |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - BCcellcenter3D             ( BCcellcenter3D.F90  )        |  !
!  |   - BCvertex3D                 ( BCvertex3D.F90      )        |  !
!  |   - pmat3D                     ( gmres3D.F90         )        |  !
!  |   - soln                       ( gmres3D.F90         )        |  !
!  |   - prods                      ( gmres3D.F90         )        |  !
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

      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      real*8,dimension(:),  allocatable :: xg
      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: funH
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

	allocate(v(ns,nkr),funH(nkr,nkr),r(nkr))
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

!     ________________________________________________________
!    |                                                        |
!    |            Matrix Am & Bm of the diffusion term        |
!    |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0 	
            Am0(i,k)  = 0.0d0 
            Am1(i,k)  = 0.0d0
            Am2(i,k)  = 0.0d0
            Am3(i,k)  = 0.0d0
            AmT(i,k)  = 0.0d0
            AmB(i,k)  = 0.0d0
            Bmv1T(i,k)= 0.0d0 
            Bmv2T(i,k)= 0.0d0 
            Bmv3T(i,k)= 0.0d0 
            Bmv1B(i,k)= 0.0d0 
            Bmv2B(i,k)= 0.0d0 
            Bmv3B(i,k)= 0.0d0 
         enddo
      enddo

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             & 
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!           Residual for the initial guess (no preconditioner)        !
!                                                                     !
!*********************************************************************!

!     -----------------------------------------------                      
!     xg = phi 
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            xg(ii) = phi(i,k)
         enddo
      enddo 
!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg   
      call pmat3D(aux1,xg,ns,                          &
                  Am0,Am1,Am2,Am3,AmT,AmB,             &
                  Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                  xc,yc,sig,dsig,No_cp,nbe,            &
                  xv,yv,sigv,dsigv,No_vp,nbev,         &
                  Hpr,h,etan,                          &
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
         call pmat3D(aux1,aux,ns,                         &
                     Am0,Am1,Am2,Am3,AmT,AmB,             &
                     Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                     xc,yc,sig,dsig,No_cp,nbe,            &
                     xv,yv,sigv,dsigv,No_vp,nbev,         &
                     Hpr,h,etan,                          &
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
            funH(i,j) = tem
            do il=1,ns
                aux(il)=aux(il)-tem*v(il,i)
            enddo
         enddo
!        ------------------------------------------- 
!        H(j+1,j)
         dem= prods(ns,aux,aux)
         dem= dsqrt(dem)
         funH(j+1,j) = dem
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
         hii =funH(i,i)
         hipi=funH(ip,i) 
         gamma= dsqrt(hii*hii+hipi*hipi)
         cosg= hii/gamma
         sing=-hipi/gamma
         do j=i,nk
            hij=funH(i,j)
            hipj=funH(ip,j)
            funH(i,j) = cosg*hij - sing*hipj
            funH(ip,j)= sing*hij + cosg*hipj
         enddo
         raux=r(i)
         r(i) =cosg*raux
         r(ip)=sing*raux
      enddo

!      ________________________________________________________
!     |                                                        |
!     |            Solution of the tridiagonal system          |  
!     |________________________________________________________|

      call soln(nk,ns,r,funH,v,xg)

!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |  
!     |________________________________________________________|

!     -----------------------------------------------                      
!     Matrix-vector product: aux1=A*xg   
      call pmat3D(aux1,xg,ns,                          &
                  Am0,Am1,Am2,Am3,AmT,AmB,             &
                  Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                  xc,yc,sig,dsig,No_cp,nbe,            &
                  xv,yv,sigv,dsigv,No_vp,nbev,         &
                  Hpr,h,etan,                          &
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
         write(*,6)'iter GMRES = ',it,', error= ',res  
         goto 1
      endif

 4    continue
      print*,'iter GMRES = ',it,', error= ',res

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

	deallocate(xg,aux,aux1,vecv)
	deallocate(v,funH,r)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: GMRES 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      END SUBROUTINE gmres3D

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



!=====================================================================!
!---------------------------------------------------------------------!
!          AUXILIAR 1: Solution of the upper triangular system        !
!---------------------------------------------------------------------!
!=====================================================================!

      subroutine soln(nk,ns,r,hg,v,xg)

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

      END

!=====================================================================!
!---------------------------------------------------------------------!
!         AUXILIAR 2: function of the product of two vectors          !
!---------------------------------------------------------------------!
!=====================================================================!

      function prods(n,x,y)    

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

      real*8::som,prods
      integer :: i

!     !**********************************************************!
!     !                        Product: x*y                      !
!     !**********************************************************!

      som=0.0d0
      do i=1,n
         som = som + x(i)*y(i)    
      enddo
      prods = som

      return
      end function prods

!=====================================================================!
!---------------------------------------------------------------------!
!               AUXILIAR 3:  Construction of vector A*x               !
!---------------------------------------------------------------------!
!=====================================================================!

      subroutine pmat3D(Axg,xg,ns,                           &
                        Am0,Am1,Am2,Am3,AmT,AmB,             &
                        Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
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

      implicit none
#     include "common.mpf"
!     ____________________________________
!     Declaration of variables  

      integer :: ns
      real*8,dimension(:)   :: Axg(ns)
      real*8,dimension(:)   :: xg(ns)
!     -------------------------------------
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
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
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
      integer:: ii


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
!      ________________________________________________________
!     |                                                        |
!     |           Multiplication Ax = A*phi + B*phiv           |
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
	    Axg(ii) =   Am0(i,k)*phi(i,k)        &
                      + Am1(i,k)*phi(jc1,k)      &
                      + Am2(i,k)*phi(jc2,k)      &
                      + Am3(i,k)*phi(jc3,k)      &
                      + AmT(i,k)*phi(i,k+1)      &       
                      + AmB(i,k)*phi(i,k-1)      &
!                     ---------------------------
                      + Bmv1T(i,k)*phiv(jv1,k)   &
                      + Bmv2T(i,k)*phiv(jv2,k)   &
                      + Bmv3T(i,k)*phiv(jv3,k)   &
                      + Bmv1B(i,k)*phiv(jv1,k-1) &
                      + Bmv2B(i,k)*phiv(jv2,k-1) &
                      + Bmv3B(i,k)*phiv(jv3,k-1)  
         enddo
      enddo

      return 
      end subroutine

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           END OF GMRES SUBROUTINES                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
