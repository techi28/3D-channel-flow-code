!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!         GMRES LINEAR SYSTEM SOLVER FOR CELL-CENTER VARIABLES        !
!          	          FOR THE VELOCITY                            !
!                        --  Nov 2013  --                             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE GMRES3Dvelocity(phi,                                &
                                 MAm0,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,  & 
                                 xc,yc,sig,dsig,No_cp,nbe,           &
                                 Hpr,h,etan,                         &
                                 tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the GMRES techni-   !
!    que for different Hesemberg matrix dimension. This system is     !
!    adecuate for the linear systems that only depends on cell-center !
!    values with coefficients "MAm" and righ-hand side "Vbm".         !
!    This program is almost identical as the solSOR3D but now we      !
!    can decide between the velocity BC that we want to impose.       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi   |(N_CELL,NZ)  | Solution & initial guess           |  !
!  |____________|_____________|____________________________________|  !
!  | ---> MAm0  |(N_CELL0,NZ) | matrix coefficient of element i    |  !
!  | ---> MAm1  |(N_CELL0,NZ) | matrix coeff. horizontal neigborn 1|  ! 
!  | ---> MAm2  |(N_CELL0,NZ) | matrix coeff. horizontal neigborn 2|  ! 
!  | ---> MAm3  |(N_CELL0,NZ) | matrix coeff. horizontal neigborn 3|  ! 
!  | ---> MAmT  |(N_CELL0,NZ) | matrix coeff. vertical top         |  ! 
!  | ---> MAmB  |(N_CELL0,NZ) | matrix coeff. vertical bottom      |  !
!  | ---> Vbm   |(N_CELL0,NZ) | right hand side of the method      |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> sig   | (NZ)        | sigma value at the cell centers    |  !
!  | ---> dsig  | (NZ)        | = sig(k+1)-sig(k+1)                |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> tagBC | integer     | Tag related to the veclocity BC    |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - BCvelcenter3D               ( BCvelocity.F90 )            |  !
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
!     -------------------------------------
      real*8,dimension(:,:) :: MAm0(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Vbm(N_CELL0,NZ)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     -------------------------------------
      integer :: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer:: jc1,jc2,jc3
!     ------------------------------------
      real*8,dimension(:,:),allocatable :: v
      real*8,dimension(:,:),allocatable :: hh
      real*8,dimension(:),  allocatable :: r
      real*8,dimension(:),  allocatable :: xg,aux,aux1,vecv
!     ------------------------------------
      data zero/1.d-8/
      integer:: ii,il,ind,ip,it,nk
      real*8 :: zero,rec
      real*8 :: tem,dem,res,prods
      real*8 :: gamma,sing,cosg,raux,hij,hii,hipj,hipi
!     ------------------------------------
      integer :: ns 
      integer :: nkr
                 ns = N_CELL0*(NZ-2)
                 nkr= DimHess

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: GMRES3Dvelocity'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

      allocate(v(ns,nkr),hh(nkr,nkr),r(nkr))
      allocate(xg(ns),aux(ns),aux1(ns),vecv(ns))

!      ________________________________________________________
!     |                                                        |
!     |               Residual for the initial guess           |  
!     |________________________________________________________|

!     ------------------------------------------------
!     Initial guess + Boundary conditions
!     -----------------------------------------------                      
!     xg = phi 
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            xg(ii) = phi(i,k)
         enddo
      enddo 
      call BCvelcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan,tagBC)
!     ------------------------------------------------
!     Residual: aux = A*xg-rhs 
      do k=2,NZ-1
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
	    aux(ii) = Vbm(i,k) - ( MAm0(i,k)*phi(i,k)   &
                                  +MAm1(i,k)*phi(jc1,k) &
                                  +MAm2(i,k)*phi(jc2,k) &
                                  +MAm3(i,k)*phi(jc3,k) &
                                  +MAmT(i,k)*phi(i,k+1) &       
                                  +MAmB(i,k)*phi(i,k-1) )       
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
         write(*,6) '* No GMRES iterations * ',0,', error =',res
         write(*,*) ' '
         goto 4
      endif
      rec = max(rec,res*eps)

!      ________________________________________________________
!     |                                                        |
!     |                     GMRES iterations                   |  
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
         do k=2,NZ-1 
            do i=1,N_CELL0
               ii = (k-2)*N_CELL0+i
               phi(i,k) = aux(ii)
            enddo
         enddo 
         call BCvelcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan,tagBC)
!        -------------------------------------------
!        Residual: aux = A*xg-rhs 
         do k=2,NZ-1
            do i=1,N_CELL0
               ii  = (k-2)*N_CELL0+i
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
	       aux1(ii) =  MAm0(i,k)*phi(i,k)   &
                         + MAm1(i,k)*phi(jc1,k) &
                         + MAm2(i,k)*phi(jc2,k) &
                         + MAm3(i,k)*phi(jc3,k) &
                         + MAmT(i,k)*phi(i,k+1) &       
                         + MAmB(i,k)*phi(i,k-1)
            enddo
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
!     ________________________________________________________
!     Solution of the tridiagonal system          
      call soln(nk,ns,r,hh,v,xg)

!     ________________________________________________________
!     Convergence criteria                  

!     ------------------------------------------------
!     Current guess + Boundary conditions
      do k=2,NZ-1 
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            phi(i,k) = xg(ii)
         enddo
      enddo 
      call BCvelcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan,tagBC)
!     ------------------------------------------------
!     Residual: aux = rhs-A*xg 
      do k=2,NZ-1
         do i=1,N_CELL0
            ii = (k-2)*N_CELL0+i
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
	    aux(ii) = Vbm(i,k) - ( MAm0(i,k)*phi(i,k)   &
                                  +MAm1(i,k)*phi(jc1,k) &
                                  +MAm2(i,k)*phi(jc2,k) &
                                  +MAm3(i,k)*phi(jc3,k) &
                                  +MAmT(i,k)*phi(i,k+1) &       
                                  +MAmB(i,k)*phi(i,k-1) )       
         enddo
      enddo
!     ----------------------------------------------                      
!     Error = sqrt(aux*aux)
      res = prods(ns,aux,aux)
      res = dsqrt(res)
!     ----------------------------------------------        
!     Criteria
6     format(t10,a26,i5,a9,e10.3)  
      if (res.lt.rec) then
         write(*,6) 'Solution GMRES 3D: iters =',it,', error =',res
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

      deallocate(v,hh,r,xg,aux,aux1,vecv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine:GMRES3Dvelocity'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           END OF GMRES SUBROUTINES                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
