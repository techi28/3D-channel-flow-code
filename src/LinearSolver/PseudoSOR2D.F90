!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            SOLUTION OF THE 2D POISSON EQUATION BY NEW S.O.R.        !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE PseudoSOR2D(phi2D,phiv2D,        &
                             rhs2D,Gamx2D,Gamy2D, & 
                             xc,yc,No_cp,nbe,     &
                             xv,yv,No_vp,nbev)                  

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the method based    !
!    in the steady-state of a parabolic diffusion equation and using  !
!    the S.O.R. technique for different relaxion factors to solve     !
!    the linear system.                                               !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  | <--> phi2D | (N_CELL)   | Solution & initial guess            |  !
!  | <--> phi2Dv| (N_VERT)   | Solution at the vertex values       |  !
!  |____________|____________|_____________________________________|  !
!  | ---> rhs2D | (N_CELL)   | Right-hand side of the system       |  ! 
!  | ---> Gamx2D| (N_CELL)   | Diffusive coefficients in x         |  !
!  | ---> Gamy2D| (N_CELL)   | Diffusive coefficients in y         |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xc,yc | (N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> No_cp | (N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe   | (N_CELL0)  | Tag type cell-center                |  !
!  |____________|____________|_____________________________________|  !
!  | ---> xv,yv | (N_VERT)   | Coordinates of the vertices         |  !
!  | ---> No_vp | (N_CELL0,3)| Numbering of the 3 cell vertices    |  !
!  | ---> nbev  | (N_VERT)   | Tag type of cell vertex             |  !
!  |____________|____________|_____________________________________|  !
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
      real*8,dimension(:)   :: phiv2D(N_VERT)
!     ----------------------------------------
      real*8,dimension(:)   :: rhs2D(N_CELL)
      real*8,dimension(:)   :: Gamx2D(N_CELL)
      real*8,dimension(:)   :: Gamy2D(N_CELL)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     ----------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:) :: Am0(N_CELL0)
      real*8,dimension(:) :: Am1(N_CELL0)
      real*8,dimension(:) :: Am2(N_CELL0)
      real*8,dimension(:) :: Am3(N_CELL0)
      real*8,dimension(:) :: Bmv1(N_CELL0)
      real*8,dimension(:) :: Bmv2(N_CELL0)
      real*8,dimension(:) :: Bmv3(N_CELL0)
      real*8,dimension(:) :: Newrhs(N_CELL0) 
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     ----------------------------------------
      real*8,dimension(:) :: phin2D(N_CELL0)
      real*8,dimension(:) :: MAm0(N_CELL0)
      real*8,dimension(:) :: MAm1(N_CELL0)
      real*8,dimension(:) :: MAm2(N_CELL0)
      real*8,dimension(:) :: MAm3(N_CELL0)
      real*8,dimension(:) :: Vbm(N_CELL0)  
!     ----------------------------------------
      real*8  :: ErrorConv,som
      integer :: iterN
!     ----------------------------------------
      integer, parameter :: ChooseSolver = 2
!     ---------------------------------------- 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Pseudo SOR 2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |          Boundary condition of the initial guess       |
!     |________________________________________________________|

!     ____________________________________
!     Cell-center BC
      call BCcellcenter2D(phi2D,xc,yc,No_cp,nbe)

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
!            Solution of the system by Pseudo Time Method             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        New Matrix                      |
!     |________________________________________________________|

      do i=1,N_CELL0 	
         MAm0(i) = ttheta*Am0(i)-ttau  
         MAm1(i) = ttheta*Am1(i)
         MAm2(i) = ttheta*Am2(i)
         MAm3(i) = ttheta*Am3(i)
      enddo

!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|

      iterN=0
11    continue
      iterN=iterN+1

!     ________________________________________________________
!     Save current iteration  

      do i=1,N_CELL0
         phin2D(i) = phi2D(i)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                Update vertex values                    |
!     |________________________________________________________|

!     ------------------------------------------------
!     Interpolation 
      call interpolation2D(phiv2D,xv,yv,No_vp,nbev, &
                           phi2D,xc,yc,No_cp,nbe)
!     ------------------------------------------------
!     Vertex BC
      call BCVertex2D(phiv2D,xv,yv,No_vp,nbev, &
                      phi2D,xc,yc,No_cp,nbe)
!      ________________________________________________________
!     |                                                        |
!     |                    New right-hand side                 |
!     |________________________________________________________|

      do i=1,N_CELL0 
!        ------------------------------
!        Cell-center terms
         jc1 = No_cp(i,1)
         jc2 = No_cp(i,2)
         jc3 = No_cp(i,3)
  	 som =   Am0(i)*phi2D(i)   &
               + Am1(i)*phi2D(jc1) &
               + Am2(i)*phi2D(jc2) &
               + Am3(i)*phi2D(jc3)
!        ------------------------------
!        Vertex terms
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
         Newrhs(i) =  rhs2D(i)-( Bmv1(i)*phiv2D(jv1) &
                                +Bmv2(i)*phiv2D(jv2) &
                                +Bmv3(i)*phiv2D(jv3) ) 
!        ------------------------------
!        Final rhs at each iteration
         Vbm(i) = -(1-ttheta)*som -ttau*phi2D(i) + Newrhs(i)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                 Solution of the system                 |
!     |________________________________________________________|

!     ________________________________________________________
!     SOR
      if (ChooseSolver.eq.1) then 
         call SOR2Dcc(phi2D,                      &
                      MAm0,MAm1,MAm2,MAm3,Vbm,    & 
                      xc,yc,No_cp,nbe)

!     ________________________________________________________
!     GMRES
      elseif (ChooseSolver.eq.2) then
         call GMRES2Dcc(phi2D,                    &
                        MAm0,MAm1,MAm2,MAm3,Vbm,  & 
                        xc,yc,No_cp,nbe)
      endif
!      ________________________________________________________
!     |                                                        |
!     |           Error of each convergence iteration          |
!     |________________________________________________________|

      ErrorConv=0.0d0 
      do i=1,N_CELL0 
         ErrorConv = ErrorConv + dabs(phi2D(i)-phin2D(i))**2
      enddo
      ErrorConv = dsqrt(ErrorConv)/float(N_CELL0)

!      ________________________________________________________
!     |                                                        |
!     |          Convergence criteria of the method            |
!     |________________________________________________________|

8     format(t10,a29,i5,a9,e10.3)
        
      if (ErrorConv.lt.epsConv) then
         write(*,*) ' '
         write(*,8) 'Sol Psedo-time Method: iters =',iterN,&
                    ', error =',ErrorConv
         write(*,*) ' '
      elseif (ErrorConv.gt.1.0d5) then
         write(*,*) ' '
         write(*,8) 'DIVERGENCE !!!!: iters =',iterN,&
                    ', error =',ErrorConv
         write(*,*) ' '
         stop
      elseif(iterN.gt.MaxConv) then
         write(*,*) ' '
         write(*,8) 'Non-convergence: iters =',iterN,&
                    ', error =',ErrorConv
         write(*,*) ' '
      else
         goto 11
      endif

19    continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Pseudo SOR 2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	       END OF Pseudo S.O.R                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
