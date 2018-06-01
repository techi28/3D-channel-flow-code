!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             S.O.R. FOR CELL-CENTER VARIABLES & VELOCITY BC          !
!                             Oct 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SOR3Dccvelocity(phi,                                 &
                                 MAm0,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,   & 
                                 xc,yc,sig,dsig,No_cp,nbe,            &
                                 Hpr,h,etan,                          &
                                 tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the S.O.R. techni-  !
!    que for different relaxion factors: relax. This system is good   !
!    for the linear systems with coefficients "MAm" and righ-hand     !
!    side "Vbm" that only depends on cell-center phi values.          !
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
!  |   - BCvelcenter3D                ( BCvelocity.F90 )           |  !
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
!     --------------------------------------
      real*8  :: errorsys,residu,som
      integer :: it

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: SOR3Dccvelocity'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

      it=0
111   continue
      it=it+1 

!     ________________________________________________________
!     Solution of the system SOR

      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
	    som =  MAm1(i,k)*phi(jc1,k) &
                 + MAm2(i,k)*phi(jc2,k) &
                 + MAm3(i,k)*phi(jc3,k) &
                 + MAmT(i,k)*phi(i,k+1) &       
                 + MAmB(i,k)*phi(i,k-1) 
	    residu = (Vbm(i,k)-som)/MAm0(i,k)-phi(i,k)
	    errorsys = errorsys + abs(residu)
	    phi(i,k) = phi(i,k) + relaxSOR*residu
         enddo
      enddo
!     ------------------------------------------------
!     Boundary conditions
      call BCvelcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,tagBC,Hpr,h,etan)

!     ________________________________________________________
!     Convergence criteria
     
      7 format(t10,a24,i5,a9,e10.3)
      if (errorsys.lt.eps) then
         write(*,7) 'Solution S0R 3D : iters =',it,', error =',errorsys
      elseif (errorsys.gt.1.0d5) then
         write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         stop
      elseif(it.gt.MaxIters) then
         write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
      else
         goto 111
      endif

119   continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine:SOR3Dccvelocity'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      End of S.O.R. Methods 3D                       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
