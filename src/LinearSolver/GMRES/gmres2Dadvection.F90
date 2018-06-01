!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               --   SOLUTION OF SYSTEM BY GMRES  --                  !
!                            April 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE gmres2Dadvection(phi2D,rhs,uu,vv,  & 
                                  xc,yc,No_cp,nbe)   

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
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> xg    | (N_CELL)   | Solution & initial guess            |  !
!  | ---> rhs   | (N_CELL)   | Right-hand side of the system       |  !
!  | ---> uu,vv | (N_CELL)   | Velocity profile                    |  ! 
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc | (N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> No_cp | (N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe   | (N_CELL0)  | Tag type cell-center                |  !
!  |____________|____________|_____________________________________|  !
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

      real*8,dimension(:) :: phi2D(N_CELL)
      real*8,dimension(:) :: rhs(N_CELL)
      real*8,dimension(:) :: uu(N_CELL)
      real*8,dimension(:) :: vv(N_CELL)

      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: h
      real*8,dimension(:),  allocatable :: r
      real*8,dimension(:),  allocatable :: xg
      real*8,dimension(:),  allocatable :: aux,aux1,vecv

      data zero/1.d-8/
      integer::il,ind,ip,it,nk,elem
      real*8::zero,rec
      real*8::tem,dem,res,prods
      real*8::gamma,sing,cosg,raux,hij,hii,hipj,hipi

      integer, parameter :: nkr = 20
      integer, parameter :: IterMaxGMRES = 200
      integer :: ns 
                 ns = N_CELL0


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: GMRES adv'
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

      do i=1,N_CELL0	
         xg(i) = phi2D(i)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Residual for the initial guess             |
!     |________________________________________________________|

!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs      
      ind = 1
      call pmat2Dadv(aux,xg,rhs,uu,vv,ind,& 
                     xc,yc,No_cp,nbe)  
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
!      ________________________________________________________
!     |                                                        |
!     |                    GMRES iterations                    |
!     |________________________________________________________|

      it = 0   
1     continue
      it = it + 1 
!     ________________________________________________________
!     Orthogonalisation d'arnoldi  

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
         ind=0 
         call pmat2Dadv(aux1,aux,rhs,uu,vv,ind,& 
                        xc,yc,No_cp,nbe)
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
!     ________________________________________________________
!     Triangularisation et modif. du second membre

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
!     ________________________________________________________
!     Solution of the tridiagonal system    

      call soln(nk,ns,r,h,v,xg)

!     ________________________________________________________
!     Convergence criteria  

!     -----------------------------------------------                      
!     Residual: aux = A*xg-rhs              
      ind = 1
      call pmat2Dadv(aux,xg,rhs,uu,vv,ind,& 
                     xc,yc,No_cp,nbe)  
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)
      res = prods(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
7     format(t10,a24,i5,a9,e10.3)

      if (res.lt.rec) then
         write(*,7) 'Sol GMRES Adv2D: iters =',it,', error =',res 
         goto 4
      elseif (res.gt.1.0d5) then
         write(*,7) 'DIVERGENCE !!!!: iters =',it,', error =',res
            stop
      elseif (it.ge.IterMaxGMRES) then
         write(*,7) 'Non-convergence: iters =',it,', error =',res
         goto 4
      else
         goto 1
      endif

 4    continue

!     ________________________________________________________
!     Final solution

      do i=1,N_CELL0	
         phi2D(i) = xg(i)
      enddo 
!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|
!
	deallocate(v,h,r,xg,aux,aux1,vecv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: GMRES'
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

      SUBROUTINE pmat2Dadv(Axg,xg,rhs,uu,vv,ind, & 
                           xc,yc,No_cp,nbe)    

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine calculates the resultant vector of the multi-    !
!    plication of the matrix A with the vector x:                     !
!                                                                     !
!                      ax =     A*x  if  ind=0                        !
!                      ax = rhs-A*x  if  ind=1                        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <-- Ax     | (N_CELL)   | Matrix-vector multiplication        |  !
!  | --> xg     | (N_CELL)   | Guess function vector               |  !
!  | --> rhs    | (N_CELL)   | right-hand side of the system       |  !
!  | --> ind    | integer    | Tag related to the gmres method     |  ! 
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

      integer :: ind
      real*8,dimension(:)   :: Axg(N_CELL0)
      real*8,dimension(:)   ::  xg(N_CELL0)
      real*8,dimension(:)   :: rhs(N_CELL)
      real*8,dimension(:,:) :: uu(N_CELL)
      real*8,dimension(:,:) :: vv(N_CELL)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 

!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8, dimension(:) :: phi(N_CELL)
      real*8, dimension(:) :: Am0(N_CELL0)
      real*8, dimension(:) :: Am1(N_CELL0)
      real*8, dimension(:) :: Am2(N_CELL0)
      real*8, dimension(:) :: Am3(N_CELL0)
      real*8, dimension(:) :: AmG(N_CELL0)
!     --------------------------------------
      real*8, dimension(:) :: dphidx(N_CELL)
      real*8, dimension(:) :: dphidy(N_CELL)
      real*8 :: sumfx,sumfy,deter
!     --------------------------------------
      real*8, dimension(:) :: C01(N_CELL0)
      real*8, dimension(:) :: C02(N_CELL0)
      real*8, dimension(:) :: C03(N_CELL0)
      real*8 :: u0j,v0j
!     -------------------------------------
      real*8, dimension(:) :: aa(1:3)
      real*8, dimension(:) :: C0j(1:3)
      real*8 :: aa0,aaB,aaT
      real*8 :: C0neg,C0pos
      real*8 :: GF,GF0,GFj,GFT,GFB
      real*8 :: GFpos,GFneg
      real*8 :: dtoVol
!     -------------------------------------
      integer:: jc,jj,jc1,jc2,jc3,elem


!     !**********************************************************!
!     !             Construction of the vector A*x               !
!     !**********************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                   Cell-center values                   |
!     |________________________________________________________|

!     ------------------------------------------------
!     Read the current inside values 
      do i=1,N_CELL0	
         phi(i) = xg(i)
      enddo
!     ------------------------------------------------
!     Boundary conditions of the cell-centers 
      call BCcellcenter2D(phi,xc,yc,No_cp,nbe)
!      ________________________________________________________
!     |                                                        |
!     |         Construction of the Matrix coefficients        |
!     |________________________________________________________|

!     ________________________________________________________
!     Mass flux calculation  

      do i=1,N_CELL0	
         do j=1,3 
	    jc = No_cp(i,j)
!           ---------------------------------
            u0j = 0.5d0*(uu(i)+uu(jc))
            v0j = 0.5d0*(vv(i)+vv(jc))
!           ---------------------------------
            C0j(j) = u0j*dyVV(i,j)-v0j*dxVV(i,j)   
	 enddo
         C01(i) = C0j(1)    
         C02(i) = C0j(2)
         C03(i) = C0j(3)
      enddo
!     ________________________________________________________
!     Gradient of phi  

      do i=1,N_CELL0	
	 sumfx = 0.0d0
	 sumfy = 0.0d0
	 do j=1,3
	    jc = No_cp(i,j)                
	    sumfx = sumfx + dxCC(i,j)*(phi(jc)-phi(i))
	    sumfy = sumfy + dyCC(i,j)*(phi(jc)-phi(i))               
	 enddo
	 deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
         dphidx(i)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
	 dphidy(i)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
!        ------------------------
!        BC gradient
	 if (nbe(i).ne.0) then  
	    do j=1,3
	       nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
	          dphidx(nc) = dphidx(i) 
	          dphidy(nc) = dphidy(i)	
	       endif
            enddo 
         endif  
      enddo
!     ________________________________________________________
!     Advection contributions 

      do i=1,N_CELL0	
         C0j(1) = C01(i)  
         C0j(2) = C02(i)
         C0j(3) = C03(i)
          
         aa0 = 0.0d0
         GF  = 0.0d0
         do j=1,3
!           ---------------------------------
	    jc = No_cp(i,j)
!           ---------------------------------
	    C0pos =  0.5d0*(C0j(j)+abs(C0j(j)))
	    C0neg = -0.5d0*(C0j(j)-abs(C0j(j)))
!           ---------------------------------
	    aa0   =  C0pos + aa0 
	    aa(j) = -C0neg               
!           ---------------------------------
	    GF0 =  dphidx(i)*(xme(i,j)-xc(i)) &
                 + dphidy(i)*(yme(i,j)-yc(i))	 
	    GFj =  dphidx(jc)*(xme(i,j)-xc(jc)) &
                 + dphidy(jc)*(yme(i,j)-yc(jc))
!           ---------------------------------
            GF = GF + C0pos*GF0-C0neg*GFj 
	 enddo
         dtoVol =  dt/areaCell(i)
         Am0(i) =  (1.0d0 + dtoVol*aa0) 
         Am1(i) =  dtoVol*aa(1)   
         Am2(i) =  dtoVol*aa(2)   
         Am3(i) =  dtoVol*aa(3)
	 AmG(i) =  dtoVol*GF  
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Multiplication Ax = A*phi + AmG            |
!     |________________________________________________________|

      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
	 Axg(i) =  Am0(i)*phi(i)   &
                 + Am1(i)*phi(jc1) &
                 + Am2(i)*phi(jc2) &
                 + Am3(i)*phi(jc3) &
                 + AmG(i)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |             Residual case  Ax = rhs - Ax               |
!     |________________________________________________________|

      if (ind.eq.1) then
         do i=1,N_CELL0
            Axg(i) = rhs(i)- Axg(i) 
         enddo
      endif

      return
      end

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      END OF AUXILIAR SUBROUTINES                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
