!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!          GMRES LINEAR SYSTEM SOLVER FOR CELL-CENTER VARIABLES       !
!                              Oct 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE GMRES2Dcc(phi2D,                    &
                           MAm0,MAm1,MAm2,MAm3,Vbm,  & 
                           xc,yc,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the GMRES techni-   !
!    que for different Hesemberg matrix dimension. This system is     !
!    adecuate for the linear systems that only depends on cell-center !
!    values with coefficients "MAm" and righ-hand side "Vbm" .        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi2D | (N_CELL)    | Solution & initial guess           |  !
!  |____________|_____________|____________________________________|  !
!  | ---> MAm0  | (N_CELL0)   | matrix coefficient of element i    |  !
!  | ---> MAm1  | (N_CELL0)   | matrix coeff. horizontal neigborn 1|  ! 
!  | ---> MAm2  | (N_CELL0)   | matrix coeff. horizontal neigborn 2|  ! 
!  | ---> MAm3  | (N_CELL0)   | matrix coeff. horizontal neigborn 3|  ! 
!  | ---> Vbm   | (N_CELL0)   | right hand side of the method      |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - BCcellcenter2D             ( BCcellcenter2D.F90 )         |  !
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

      real*8,dimension(:)   :: phi2D(N_CELL)
!     ------------------------------------
      real*8,dimension(:)   :: MAm0(N_CELL0)
      real*8,dimension(:)   :: MAm1(N_CELL0)
      real*8,dimension(:)   :: MAm2(N_CELL0)
      real*8,dimension(:)   :: MAm3(N_CELL0) 
      real*8,dimension(:)   :: Vbm(N_CELL0)
!     ------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer:: jc1,jc2,jc3
!     ------------------------------------
      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: h
      real*8,dimension(:),  allocatable :: r
      real*8,dimension(:),  allocatable :: xg,aux,aux1,vecv
!     ------------------------------------
      data zero/1.d-8/
      integer:: il,ind,ip,it,nk
      real*8 :: zero,rec
      real*8 :: tem,dem,res,prods2D
      real*8 :: gamma,sing,cosg,raux,hij,hii,hipj,hipi
!     ------------------------------------
      integer :: ns 
      integer :: nkr
                 ns = N_CELL0
                 nkr= DimHess

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: SolGMRES2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

      allocate(v(ns,nkr),h(nkr,nkr),r(nkr))
      allocate(xg(ns),aux(ns),aux1(ns),vecv(ns))

!      ________________________________________________________
!     |                                                        |
!     |             Residual for the initial guess             |
!     |________________________________________________________|

!     ------------------------------------------------
!     Initial guess + Boundary conditions
      do i=1,N_CELL0	
         xg(i) = phi2D(i)
      enddo
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Residual: aux = A*xg-rhs 
      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
	 aux(i) = Vbm(i)- ( MAm0(i)*phi2D(i)     &
                           +MAm1(i)*phi2D(jc1)   &
                           +MAm2(i)*phi2D(jc2)   &
                           +MAm3(i)*phi2D(jc3) )
      enddo 
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)	
      res = prods2D(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
      rec = zero
      if (res.lt.rec) then
         write(*,*) ' '
         write(*,6) '* No GMRES iterations * ',0,', error =',res
         write(*,*) ' '
         goto 4
      endif
      rec = max(rec,res*eps)
!      ________________________________________________________
!     |                                                        |
!     |                      GMRES iterations                  |
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
!        Read the current inside values + BC
         do i=1,N_CELL0	
            phi2D(i) = aux(i)
         enddo
         call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)
!        ------------------------------------------- 
!        Product: aux1 = A*aux
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
	    aux1(i) =  MAm0(i)*phi2D(i)     &
                     + MAm1(i)*phi2D(jc1)   &
                     + MAm2(i)*phi2D(jc2)   &
                     + MAm3(i)*phi2D(jc3)
         enddo
!        ------------------------------------------- 
!        Update: aux = aux1 
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
!     ________________________________________________________
!     Triangularisation et modif. du second membre
5     continue

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
      call soln2D(nk,ns,r,h,v,xg)

!     ________________________________________________________
!     Convergence criteria                  

!     ------------------------------------------------
!     Current guess + Boundary conditions
      do i=1,N_CELL0	
         phi2D(i) = xg(i)
      enddo
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Residual: aux = rhs - A*xg  
      do i=1,N_CELL0
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
	 aux(i) = Vbm(i)- ( MAm0(i)*phi2D(i)     &
                           +MAm1(i)*phi2D(jc1)   &
                           +MAm2(i)*phi2D(jc2)   &
                           +MAm3(i)*phi2D(jc3) )
      enddo  
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)
      res = prods2D(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
6     format(t10,a26,i5,a9,e10.3)  
      if (res.lt.rec) then
         write(*,6) 'Solution GMRES 2D: iters =',it,', error =',res
         goto 4
      elseif (res.gt.1.0d5) then
         write(*,6) 'DIVERGENCE !!!!!!: iters =',it,', error =',res
         stop
      elseif (it.ge.MaxIters) then
         write(*,6) 'Non-convergence!!: iters =',it,', error =',res
         goto 4
      else
         goto 1
      endif

4     continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      deallocate(v,h,r,xg,aux,aux1,vecv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine:SolGMRES2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        END OF GMRES SUBROUTINES                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
