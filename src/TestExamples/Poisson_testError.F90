!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   TEST OF THE POISSON EQUATION                      !
!                             May 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Poisson_testError(phi,phiv,                      &
                                   SolApprov,SolExactv,SolErrorv, &
                                   xc,yc,sig,dsig,No_cp,nbe,      &
                                   xv,yv,sigv,dsigv,No_vp,nbev,   &     
                                   Hpr,h,etan,                    &
                                   Hprv,hv,etav)              
!---------------------------------------------------------------------!   
!                                                                     !
!     This calculates the error of the poisson equation of the        !
!     variable phi, given the diffusive values, right-hand side       !
!     and the correct boundary conditions.                            !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |    Size   | Description                        |  !  
!  |______________|___________|____________________________________|  !  
!  | <--SolAppro  |(N_CELL,NZ)| Approximation at the cell-centers  |  !
!  | <--SolExact  |(N_CELL,NZ)| Exact solution at the cell-centers |  !
!  | <--SolError  |(N_CELL,NZ)| Errors at the cell-centers         |  !
!  |______________|___________|____________________________________|  !  
!  | <--SolApprov |(N_VERT,NZ)| Approximation at the vertices      |  !
!  | <--SolExactv |(N_VERT,NZ)| Exact solution at the vertices     |  !
!  | <--SolErrorv |(N_VERT,NZ)| Errors at the vertices             |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !  
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbev   |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!        Common variables used                                        !
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:,:) :: SolApprov(N_VERT,NZ-1)
      real*8,dimension(:,:) :: SolExactv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: SolErrorv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     --------------------------------------
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
      real*8,dimension(:,:) :: SolAppro(N_CELL,NZ)
      real*8,dimension(:,:) :: SolExact(N_CELL,NZ)
      real*8,dimension(:,:) :: SolError(N_CELL,NZ)
!     ----------------------------------------      
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8 :: maxErrorA,maxErrorR
      real*8 :: sumErrorA,sumErrorR
      real*8 :: MAXmaxErrorA,MAXmaxErrorR
      real*8 :: SUMsumErrorA,SUMsumErrorR
!     ----------------------------------------
      real*8 :: funSolExam3D,funDiffExam3D
      real*8 :: funExamNSrhsp,funExamNSp
      real*8 :: x,y,z,c

!*********************************************************************!
!                                                                     !
!                                 Error                               !
!                                                                     !
!*********************************************************************!

!     ________________________________________________
!     Constant in the case of the Neumann BC
      c = 0.0d0
#     ifdef Key_NeumannBCp
         z = sig(3)*Hpr(1)-h(1)
         !c = funSolExam3D(xc(1),yc(1),z)-phi(1,3)
         c = funExamNSp(xc(1),yc(1),z)-phi(1,3)
#     endif

!      ________________________________________________________
!     |                                                        |
!     |          Exact solution & Error  (cell-center)         |
!     |________________________________________________________|

      do k=1,NZ 
         do i=1,N_CELL
            x = xc(i)
            y = yc(i)
            z = sig(k)*Hpr(i)-h(i)
            !SolExact(i,k) = funSolExam3D(x,y,z)
            SolExact(i,k) = funExamNSp(x,y,z,time)
            SolAppro(i,k) = phi(i,k)+c
            SolError(i,k) = abs(SolAppro(i,k)-SolExact(i,k))
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |            Exact solution & Error  (Vertex)            |
!     |________________________________________________________|

      do k=1,NZ-1
         do nv=1,N_VERT   
            x = xv(nv)
            y = yv(nv)
            z = sigv(k)*Hprv(nv)-hv(nv)
            !SolExactv(nv,k) = funSolExam3D(x,y,z)
            SolExactv(nv,k) = funExamNSp(x,y,z)
            SolApprov(nv,k) = phiv(nv,k)+c
            SolErrorv(nv,k) = abs(SolApprov(nv,k)-SolExactv(nv,k))
         enddo
      enddo

!*********************************************************************!
!                                                                     !
!                               Norms                                 !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                    Lmax & L2 norms                     |
!     |________________________________________________________|

      maxErrorA = 0.0d0
      maxErrorR = 0.0d0
      sumErrorA = 0.0d0
      sumErrorR = 0.0d0
      do k=2,NZ-1 
         do i=1,N_CELL0 
            ErrorA(i,k) = abs(SolAppro(i,k)-SolExact(i,k))
            ErrorR(i,k) = abs(SolExact(i,k))
            maxErrorA = max(maxErrorA,ErrorA(i,k))
            maxErrorR = max(maxErrorR,ErrorR(i,k))
            sumErrorA = sumErrorA+ErrorA(i,k)**2
            sumErrorR = sumErrorR+ErrorR(i,k)**2
         enddo
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call MAX_parallel(maxErrorA,MAXmaxErrorA)
         call MAX_parallel(maxErrorR,MAXmaxErrorR)
         call SUM_parallel(sumErrorA,SUMsumErrorA)
         call SUM_parallel(sumErrorR,SUMsumErrorR)
         maxErrorA = MAXmaxErrorA
         maxErrorR = MAXmaxErrorR
         sumErrorA = SUMsumErrorA
         sumErrorR = SUMsumErrorR
#     endif	
!     =============== END ================    
!     ====================================
      sumErrorA = dsqrt(sumErrorA/(N_CELL0global*(NZglobal-2)))
      sumErrorR = dsqrt(sumErrorR/(N_CELL0global*(NZglobal-2)))
      maxErrorR = maxErrorA/maxErrorR
      sumErrorR = sumErrorA/sumErrorR
!      ________________________________________________________
!     |                                                        |
!     |                       Display Error                    |
!     |________________________________________________________|

7     format(t26,20a)
8     format(t10,60a)
9     format(t10,a3,t14,e10.3,t26,e10.3,t38,e10.3,t50,e10.3)

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      IF (rang_topo.eq.0) THEN
         write(*,8)'===================================================='
         write(*,8)'____________________________________________________'
         write(*,8)'                                                    '
         write(*,8)'             MPI TEST: Poisson equation             '
         write(*,'(t30,a4,i3)') ' N = ',NN
         write(*,'(t26,a9,f6.3)') ' omega = ',RelaxSOR
#        ifdef KeyMSOR 
         write(*,7) ' Method: MSOR'
#        endif
#        ifdef KeySiOR 
         write(*,7) ' Method: Jacobi'
#        endif
#        ifdef KeyGMRES 
         write(*,7) ' Method: GMRES'
#        endif
#        ifdef KeyPseudoTime 
         write(*,7) ' Method: PseudoTime'
#        endif
         write(*,8)'____________________________________________________'
         write(*,8)'                                                    '
         write(*,8)'          Absolute error          Relative error    '
         write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
         write(*,8)'----------------------------------------------------'
         write(*,9) '   ',maxErrorA,sumErrorA,maxErrorR,sumErrorR
         write(*,8)'----------------------------------------------------'
         write(*,8)'===================================================='
         write(*,*)' '
         write(*,8)'  REMARK: Functions used: funExamNSp,funExamNSrhsp'
         write(*,*)' '
      ENDIF
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
         write(*,8)'===================================================='
         write(*,8)'____________________________________________________'
         write(*,8)'                                                    '
         write(*,8)'               TEST: Poisson equation               '
         write(*,'(t30,a4,i3)') ' N = ',NN
         write(*,'(t26,a9,f6.3)') ' omega = ',RelaxSOR
#        ifdef KeySOR 
         write(*,7) ' Method: SOR'
#        endif
#        ifdef KeySiOR 
         write(*,7) ' Method: Jacobi'
#        endif
#        ifdef KeyGMRES 
         write(*,7) ' Method: GMRES'
#        endif
#        ifdef KeyPseudoTime 
         write(*,7) ' Method: PseudoTime'
#        endif
         write(*,8)'____________________________________________________'
         write(*,8)'                                                    '
         write(*,8)'          Absolute error          Relative error    '
         write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
         write(*,8)'----------------------------------------------------'
         write(*,9) '   ',maxErrorA,sumErrorA,maxErrorR,sumErrorR
         write(*,8)'----------------------------------------------------'
         write(*,8)'===================================================='
         write(*,*)' '
         write(*,8)'  REMARK: Functions used: funExamNSp,funExamNSrhsp'
         write(*,*)' '
#     endif
!     =============== END ================    
!     ==================================== 

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF testPoisson                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
