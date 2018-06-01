!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        MODULE:  nsmputilcuda3D                      !
!---------------------------------------------------------------------!
!                     -  subroutine: initialization CUDA              !
!                     -  subroutine: allocate_cuda                    !
!                     -  subroutine: deallocate_cuda                  !
!                     -  subroutine: transfergeometry_cuda            !
!                     -  subroutine: transferInterpo_cuda             !
!                     -  subroutine: transfermatrixA_cuda             !
!                     -  subroutine: transfermatrixB_cuda             !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!    Device variables (processors variables):                         !
!   _______________________________________________________________   !
!  |     Name   |    Size  |          Description                  |  !
!  |____________|__________|_______________________________________|  !
!  | d_dpstar   |(nmax)    | Solution variable: dp                 |  !
!  | d_yg       |(nmax)    | Right hand side of the system         |  !
!  | d_bb       |(nmax)    | Iteration righ-hand side              |  !
!  | d_am       |(2,nmax)  | Matrix coeff. of implicit variables   |  !
!  | d_cm       |(4:6,nmax)| Matrix coeff. of explicit variables   |  !
!  | d_error    |(nx-2)    | Error = SUM(residual**2)              |  !
!  |____________|__________|_______________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!

      MODULE parallelCUDA

      use cudafor

!     =========================================================
      !-----------------.Device--------------------------------
      real*8, dimension(:),   device, allocatable :: d_phiNew
      real*8, dimension(:),   device, allocatable :: d_phi
      real*8, dimension(:),   device, allocatable :: d_phiv
      real*8, dimension(:),   device, allocatable :: d_rhs
!     -----------------
      real*8, dimension(:,:), device, allocatable :: d_Am
      real*8, dimension(:,:), device, allocatable :: d_Bm
!     -----------------
      real*8, dimension(:),   device, allocatable :: d_funfB
      real*8, dimension(:),   device, allocatable :: d_funfBv
!     -----------------
      real*8, dimension(:),   device, allocatable :: d_fundfdnB
      real*8, dimension(:),   device, allocatable :: d_fundfdnBv
!     -----------------
      real*8, dimension(:),   device, allocatable :: d_phivC
      integer,dimension(:),   device, allocatable :: d_Color
      real*8, dimension(:),   device, allocatable :: d_dlCE      
      real*8, dimension(:),   device, allocatable :: d_dlTB
!     -----------------      
      real*8, dimension(:,:), device, allocatable :: d_erreur
!     -----------------
      integer,dimension(:,:), device, allocatable :: d_No_cp
      integer,dimension(:,:), device, allocatable :: d_No_vp
      integer,dimension(:,:), device, allocatable :: d_Surro
      real*8, dimension(:,:), device, allocatable :: d_Weigh
      real*8, dimension(:),   device, allocatable :: d_dzT
      real*8, dimension(:),   device, allocatable :: d_dzB
      integer,dimension(:),   device, allocatable :: d_nbev

!     =========================================================
!     Threads and blocks distribution 
      integer :: THREADSPBLOCK_X,THREADSPBLOCK_Y 
      type(dim3) :: grid,gridv,tBlock

      integer :: shared_mem
      integer :: N_C0,N_C,N_V,N_V0,MaxNumC

      CONTAINS

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.1  |                 Initialization CUDA                    |
!     |______|________________________________________________________|

      SUBROUTINE initialisation_cuda

#     include "cppdefs.h" 
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------   

      integer :: NumHc,NumZc,NumHv,NumZv

      !--------------------------------------------------------------!
      !                                                              !
      !   <--- tBlock  : Dimension of the block                      !
      !   <--- grid    : Block topology (cell-centers points)        !
      !   <--- gridv   : Block topology (vertex points)              !
      !                                                              !
      !--------------------------------------------------------------!

!     ============================
!     SETUP VALUES

      THREADSPBLOCK_X = 1
      THREADSPBLOCK_Y = NZglobal
      
!     ============================
!     CUDA: Domain decomposition
      
      tBlock = dim3(THREADSPBLOCK_X,THREADSPBLOCK_Y,1)
      NumHc = ceiling(real(N_CELL)/tBlock%x)
      NumZc = ceiling(real(NZ)/tBlock%y)
      NumHv = ceiling(real(N_VERT)/tBlock%x)
      NumZv = ceiling(real(NZ-1)/tBlock%y)
      grid  = dim3(NumHc,NumZc,1)
      gridv = dim3(NumHV,NumZv,1)
      
      print*,' '
      print*,'           ---------------------------------------  '
      print*,'                                  ^ |axb  axb  axb| '
      print*,'            BLOCK CONFIGURATION:  | |axb  axb  axb| '
      print*,'                                  D |axb  axb  axb| '
      print*,'                                  | |axb  axb  axb| '
      print*,'                                  v  <---- C ---->  '
      print*,'           ---------------------------------------  '
      print*,'           a = Threads horizontal : ',THREADSPBLOCK_X
      print*,'           b = Threads vertical   : ',THREADSPBLOCK_Y
      print*,'           C = No. horizont blocks: ',grid%x
      print*,'           D = No. vertical blocks: ',grid%y
      print*,'           ---------------------------------------  '
      print*,' '

!     ============================
!     CUDA: Share memory capacity for phi

      shared_mem = THREADSPBLOCK_X*THREADSPBLOCK_Y*8
      if (shared_mem.gt.64*1024) then
         print*,'Shared memory setting for update kernel exceeds 64KB,'
         print*,'please correct it.'
         stop
      endif
      
      END SUBROUTINE initialisation_cuda 
            

!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.2  |               Allocate CUDA variables                  |
!     |______|________________________________________________________|
      
      SUBROUTINE allocate_cuda

#     include "cppdefs.h" 
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------

!     ============================
!     Total number of elements

      N_C0  = N_CELL0*NZ
      N_C   = N_CELL*NZ
      N_V0  = N_VERT*(NZ-1)
      N_V   = N_VERT*NZ

!     ============================
!     CUDA: Allocate memory
      allocate(d_phi(N_C))
      allocate(d_phiNew(N_C))
      allocate(d_rhs(N_C))
      allocate(d_funfB(N_C))
      allocate(d_fundfdnB(N_C))             
      allocate(d_Am(5,N_C0))
      allocate(d_Bm(6,N_C0))
      allocate(d_phiv(N_V0)) 
      allocate(d_funfBv(N_V0))
      allocate(d_fundfdnBv(N_V0))
      allocate(d_phivC(N_V))
      allocate(d_No_cp(N_CELL,3))
      allocate(d_No_vp(N_CELL0,3))
      allocate(d_dzT(NZ-1))
      allocate(d_dzB(NZ-1))
      allocate(d_nbev(N_VERT))
      allocate(d_erreur(grid%x,grid%y))
      allocate(d_dlCE(N_CELL0))
      allocate(d_dlTB(2)) 
#     if defined(KeyMSOR_cuda) || defined(KeyPDMSOR_cuda)
      allocate(d_Color(N_CELL0))
#     endif

      END SUBROUTINE allocate_cuda 
                    
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.3  |              Deallocate CUDA variables                 |
!     |______|________________________________________________________|
      
      SUBROUTINE deallocate_cuda

#     include "cppdefs.h"
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------  
      
!     ============================
!     Deallocate memory (CUDA variables)
      
      deallocate(d_phi(N_C))
      deallocate(d_phiNew(N_C))
      deallocate(d_rhs(N_C))
      deallocate(d_funfB(N_C))
      deallocate(d_fundfdnB(N_C))           
      deallocate(d_Am(5,N_C0))
      deallocate(d_Bm(6,N_C0))
      deallocate(d_phiv(N_V0)) 
      deallocate(d_funfBv(N_V0))
      deallocate(d_fundfdnBv(N_V0))
      deallocate(d_phivC(N_V0))
      deallocate(d_No_cp(N_CELL,3))
      deallocate(d_No_vp(N_CELL0,3))
      deallocate(d_Surro(N_VERT,MaxNumC))
      deallocate(d_Weigh(N_VERT,MaxNumC))
      deallocate(d_dzT(NZ-1))
      deallocate(d_dzB(NZ-1))
      deallocate(d_nbev(N_VERT))
      deallocate(d_erreur(grid%x,grid%y))
      deallocate(d_dlCE(N_CELL0))
      deallocate(d_dlTB(2))
#     if defined(KeyMSOR_cuda) || defined(KeyPDMSOR_cuda)
      deallocate(d_Color(N_CELL0))
#     endif

      END SUBROUTINE deallocate_cuda 
      
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.4  |      Transfer all geometry values (need it once)       |
!     |______|________________________________________________________|

      SUBROUTINE transfergeometry_cuda(No_cp,No_vp,nbev,sig,sigv)

#     include "cppdefs.h"
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------              

      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)
      real*8,dimension(:)    :: sig(NZ)
      real*8,dimension(:)    :: sigv(NZ-1)
      real*8,dimension(:)    :: dsigTB(2)
      real*8, dimension(:), allocatable :: dzT_vec,dzB_vec
      real*8, dimension(:), allocatable :: dlCE_vec        
      real*8 :: dzB,dzT
      integer:: ii                        
      
      !--------------------------------------------------------------!
      !                                                              !
      !   dzT      : top step of the vertical direction (vector)     !
      !   dzTB     : bottom step of the vertical direction (vector)  !
      !   d_No_cp  : three center neighbor indexes                   !
      !   d_No_vp  : three vertex neighbor indexes                   !
      !   d_nbev   : inside or outside tag (vertex)                  !
      !   d_Color  : Color tag of the element (cell-center)          !
      !                                                              !
      !--------------------------------------------------------------!

!     ==================================
!     Distance from edge to ghost center

#     if defined(Key_NeumannBCp) || defined(Key_MixBCp)
         !--------------
         ! Vertical
         dsigTB(1) = sig(NZ)-sig(NZ-1)
         dsigTB(2) = sig(2)-sig(1)
         d_dlTB = dsigTB
         !--------------
         ! Horizontal  
         allocate(dlCE_vec(N_CELL0))
         dlCE_vec = 0.0d0
         do ii=N_CELL0+1,N_CELLexact
            i = No_cp(ii,1)
            j = No_cp(ii,2) 
            dlCE_vec(i) = dlCE(i,j)
         enddo
         d_dlCE = dlCE_vec
         deallocate(dlCE_vec(N_CELL0))
#     endif    

!     ==================================
!     Geometry data
     
      d_No_cp = No_cp
      d_No_vp = No_vp 
      d_nbev  = nbev

!     ==================================
!     Step vertical dzT, dzB  

      allocate(dzT_vec(NZ-1))
      allocate(dzB_vec(NZ-1))          
      do k=1,NZ-1    
         dzT = abs(sigv(k)-sig(k+1))
         dzB = abs(sigv(k)-sig(k))            
         dzT_vec(k) = dzT/(dzT+dzB)
         dzB_vec(k) = dzB/(dzT+dzB)
      enddo
      d_dzT = dzT_vec
      d_dzB = dzB_vec      
      deallocate(dzT_vec(NZ-1))
      deallocate(dzB_vec(NZ-1))        

!     ==================================
!     Color  
     
#     if defined(KeyMSOR_cuda) || defined(KeyPDMSOR_cuda)
         d_Color = ColorCell
#     endif

      END SUBROUTINE transfergeometry_cuda
 
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.5  |      Transfer all geometry values (need it once)       |
!     |______|________________________________________________________|

      SUBROUTINE transferInterpo_cuda

#     include "cppdefs.h"
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------               

      real*8, dimension(:,:), allocatable :: Surro
      real*8, dimension(:,:), allocatable :: Weigh         
      
      !--------------------------------------------------------------!
      !                                                              !
      !   d_Surro  : surrounding index of cell-centers               !
      !   d_Weigh  : weight for each sourounding cell-center         !                  
      !                                                              !
      !--------------------------------------------------------------!
      
!     ==================================
!     Surrounding data for interpolation

!     ------------------------------
!     Maximum number of surrouding neigh
      MaxNumC = 1
      do nv=1,N_VERT
         MaxNumC = max(MaxNumC,Dimsurrounding(nv))
      enddo
!     ------------------------------
!     Allocate device      
      allocate(d_Surro(N_VERT,MaxNumC))
      allocate(d_Weigh(N_VERT,MaxNumC)) 

!     ------------------------------
!     Allocate host
      allocate(Surro(N_VERT,MaxNumC))
      allocate(Weigh(N_VERT,MaxNumC))     
!     ------------------------------ 
!     Assign values     
      do nv=1,N_VERT
         do j=1,MaxNumC
            if (j.le.Dimsurrounding(nv)) then
               Surro(nv,j) = surrounding(nv,j)
               Weigh(nv,j) = weight(nv,j)/dlVsum(nv)
            else
               Surro(nv,j) = 1
               Weigh(nv,j) = 0.0d0            
            endif
         enddo
      enddo 
!     ------------------------------  
!     Transfer to buffer 
      d_Surro = Surro
      d_Weigh = Weigh

!     ------------------------------
!     Deallocate
      deallocate(Surro(N_VERT,MaxNumC))
      deallocate(Weigh(N_VERT,MaxNumC))

      END SUBROUTINE transferInterpo_cuda
 
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.6  |      Transfer matrices coefficients A to buffer        |
!     |______|________________________________________________________|

      SUBROUTINE transfermatrixA_cuda(Am1,Am2,Am3,AmT,AmB)
            
#     include "cppdefs.h"
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------

      real*8,dimension(:,:)  :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:)  :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:)  :: AmB(N_CELL0,NZ)
      real*8, dimension(:,:), allocatable :: A_vec
      integer :: s

      allocate(A_vec(5,N_CELL0*NZ))

      do i=1,N_CELL0
         do k=1,NZ
            s = k+(i-1)*NZ      
            A_vec(1,s) = Am1(i,k)
            A_vec(2,s) = Am2(i,k)
            A_vec(3,s) = Am3(i,k)
            A_vec(4,s) = AmT(i,k)
            A_vec(5,s) = AmB(i,k)
         enddo
      enddo 
      d_Am = A_vec

      deallocate(A_vec(5,N_CELL0*NZ))

      END SUBROUTINE transfermatrixA_cuda
      
!      _______________________________________________________________ 
!     |      |                                                        |
!     | 1.7  |      Transfer matrices coefficients B to buffer        |
!     |______|________________________________________________________|

      SUBROUTINE transfermatrixB_cuda(Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B)
            
#     include "cppdefs.h"
      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------

      real*8,dimension(:,:)  :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv3B(N_CELL0,NZ)
      real*8, dimension(:,:), allocatable :: B_vec
      integer :: s

      allocate(B_vec(6,N_CELL0*NZ)) 

      do i=1,N_CELL0
         do k=1,NZ
            s = k+(i-1)*NZ
            B_vec(1,s) = Bmv1T(i,k)
            B_vec(2,s) = Bmv2T(i,k)
            B_vec(3,s) = Bmv3T(i,k)
            B_vec(4,s) = Bmv1B(i,k)
            B_vec(5,s) = Bmv2B(i,k)
            B_vec(6,s) = Bmv3B(i,k)
         enddo
      enddo 
      d_Bm  = B_vec

      deallocate(B_vec(6,N_CELL0*NZ))

      END SUBROUTINE transfermatrixB_cuda

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!             Method SOR (Successive-Over-Relaxation) by CUDA         !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE MethSOR_CUDA(phi,phiv,rhs,                        &
                              Am1,Am2,Am3,AmT,AmB,                 &
                              Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                              No_cp,No_vp,nbev,sig,sigv,           &
                              funfB,funfBv,                        &
                              fundfdnB,fundfdnBv)

#     include "cppdefs.h"
      USE geometry
      implicit none 

!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|

      real*8,dimension(:,:)  :: phi(N_CELL,NZ)
      real*8,dimension(:,:)  :: phiv(N_VERT,NZ-1)
      real*8,dimension(:,:)  :: rhs(N_CELL,NZ)
      real*8,dimension(:,:)  :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:)  :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:)  :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:)  :: Bmv3B(N_CELL0,NZ) 
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)
      real*8,dimension(:)    :: sig(NZ)
      real*8,dimension(:)    :: sigv(NZ-1)
      real*8,dimension(:,:)  :: funfB(N_CELL,NZ)
      real*8,dimension(:,:)  :: funfBv(N_VERT,NZ-1)
!     ---------------------      
      real*8,dimension(:,:)  :: fundfdnB(N_CELL,NZ)
      real*8,dimension(:,:)  :: fundfdnBv(N_VERT,NZ-1)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
 
      real*8, dimension(:,:), allocatable :: h_erreur
!     ---------------------
      real*8, dimension(:),   allocatable :: phi_vec
      real*8, dimension(:),   allocatable :: rhs_vec
      real*8, dimension(:),   allocatable :: fun_vec
      real*8, dimension(:),   allocatable :: phiv_vec
      real*8, dimension(:),   allocatable :: funv_vec
!     ---------------------
      real*8, dimension(:),   allocatable :: dfdn_vec
      real*8, dimension(:),   allocatable :: dfdnv_vec
!     __________________________________
!     Declaration of local variables
      integer:: iter,s,ss 
      real*8 :: erreur,erreurOld,errorNeum
!     ----------------------------------
      real :: startT,finishT,timeTotal
!      ____________________________________
!     |                                    |
!     |              Parameters            |
!     |____________________________________|

      real*8, parameter :: relaxSORCUDA = 1.2d0
         
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|
 
      allocate(h_erreur(grid%x,grid%y))
!     ---------------------
      allocate(phi_vec(N_CELL*NZ))
      allocate(phiv_vec(N_VERT*(NZ-1)))
      allocate(rhs_vec(N_CELL*NZ))
!     ---------------------  
      allocate(fun_vec(N_CELL*NZ))
      allocate(funv_vec(N_VERT*(NZ-1)))
!     ---------------------
      allocate(dfdn_vec(N_CELL*NZ))       
      allocate(dfdnv_vec(N_VERT*(NZ-1)))

!      ________________________________________________________
!     |                                                        |
!     |       Transfer information from host to device         |
!     |________________________________________________________|

!     =========================================================
!     CUDA: Transfer data from host to device: phi      
      IF ((time-dt).le.1e-06) THEN !<-- Only at initial time
         do i=1,N_CELL
            do k=1,NZ
               s = k+(i-1)*NZ
               phi_vec(s) = phi(i,k)            
            enddo
         enddo                                  
         d_phi    = phi_vec  
         d_phiNew = d_phi
      ENDIF

!     =========================================================
!     CUDA: Transfer data from host to device: rhs
      do i=1,N_CELL
         do k=1,NZ
            s = k+(i-1)*NZ
            rhs_vec(s) = rhs(i,k)
         enddo
      enddo
      d_rhs = rhs_vec

!     =========================================================
!     CUDA: Transfer data from host to device: DIRICHLET functions
#     if defined(Key_DirichletBCp) || defined(Key_MixBCp)      
         do i=1,N_CELL
            do k=1,NZ
               s = k+(i-1)*NZ           
               fun_vec(s) = funfB(i,k)
            enddo
         enddo                             
         do nv=1,N_VERT
            do k=1,NZ-1
               s = k+(nv-1)*(NZ-1)           
               funv_vec(s) = funfBv(nv,k)
            enddo
         enddo 
         d_funfB  = fun_vec
         d_funfBv = funv_vec 
#     endif      

!     =========================================================
!     CUDA: Transfer data from host to device: NEUMANN functions
#     if defined(Key_NeumannBCp) || defined(Key_MixBCp)
         do i=1,N_CELL
            do k=1,NZ
               s = k+(i-1)*NZ           
               dfdn_vec(s) = fundfdnB(i,k)
            enddo
         enddo                             
         do nv=1,N_VERT
            do k=1,NZ-1
               s = k+(nv-1)*(NZ-1)           
               dfdnv_vec(s) = fundfdnBv(nv,k)
            enddo
         enddo
         d_fundfdnB  = dfdn_vec
         d_fundfdnBv = dfdnv_vec
#     endif      

!     =========================================================
!     CUDA: Transfer from host to device: Matrix A
      call transfermatrixA_cuda(Am1,Am2,Am3,AmT,AmB)

!     =========================================================
!     CUDA: Transfer from host to device: Matrix B

#     ifdef KeyFixedFreeSurface
         IF ((time-dt).le.1e-06) THEN   !<-- Only at initial time
         print*,'Transfer of B matrix only at the fisrt step'     
         call transfermatrixB_cuda(Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B)         
         ENDIF
#     else
         call transfermatrixB_cuda(Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B)
#     endif

!      ________________________________________________________
!     |                                                        |
!     |           Boundary conditions & Interpolation          |
!     |________________________________________________________| 

!     ===========================
!     KERNEL: Boundary conditions 
!     ===========================

#     ifdef Key_DirichletBCp
         call kernel_cal_phiDirichletBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)
#     endif
!     -----------------------
#     ifdef Key_NeumannBCp
         call kernel_cal_phiNeumannBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)
#     endif
!     -----------------------
#     ifdef Key_MixBCp
         call kernel_cal_phiMixBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)
#     endif

!     ===========================            
!     KERNEL: Interpolation 
!     =========================== 
                 
      call kernel_cal_phiV<<<gridv,tBlock,shared_mem>>>(MaxNumC,NZ) 

!      ________________________________________________________
!     |                                                        |
!     |                    Begin of iterations                 |
!     |________________________________________________________|

      call cpu_time(startT)
      
      iter=0
500   continue
      iter=iter+1

      erreurOld = erreur
      erreur = 0.0d0

!     ===========================
!     KERNEL: Solution SOR
!     ===========================

!     ________________________
#     ifdef KeyJacobi_cuda
         relaxSOR = 1.0d0 
         call kernel_cal_phiJacobi<<<grid,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR)
         d_phi = d_phiNew
         h_erreur = d_erreur
         erreur   = sum(h_erreur)
#     endif
!     -----------------------
#     ifdef KeyJSOR_cuda
         relaxSOR = 1.1d0 
         call kernel_cal_phiJSOR<<<grid,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR)
         h_erreur = d_erreur
         erreur   = sum(h_erreur)
#     endif
!     -----------------------
#     ifdef KeyMSOR_cuda
         relaxSOR = relaxSORCUDA
         do s=0,N_COLOR-1
            call kernel_cal_phiMSOR<<<grid,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s)
            h_erreur = d_erreur
            erreur = erreur + sum(h_erreur)
         enddo
#     endif        

!     ===========================
!     KERNEL: Boundary conditions
!     ===========================

#     ifdef Key_DirichletBCp
         call kernel_cal_phiDirichletBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)
#     endif
!     -----------------------
#     ifdef Key_NeumannBCp
         call kernel_cal_phiNeumannBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)
#     endif
!     -----------------------
#     ifdef Key_MixBCp
         call kernel_cal_phiMixBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)
#     endif

!     ===========================
!     KERNEL: Interpolation
!     ===========================
                 
      call kernel_cal_phiV<<<gridv,tBlock,shared_mem>>>(MaxNumC,NZ) 

!      ________________________________________________________
!     |                                                        |
!     |            Alternative Error (Neumann stop)            |
!     |________________________________________________________|

      errorNeum = abs(erreur-erreurOld)
!      ________________________________________________________
!     |                                                        |
!     |                 Criteria of convergence                |
!     |________________________________________________________|

!     ________________________________________________________
!     Save multiples simulations varying epsilon  

#     ifdef KeySaveVaryEps 
         write(7100,*),iter,erreur
#     endif
!     ________________________________________________________
!     Criteria
      if ((erreur.lt.eps).or.(errorNeum.lt.1d-5*eps)) then
      !if (erreur.lt.eps) then
         write(*,*) ' '
         write(*,7) 'Pressure SOR-GPU: iters =',iter, &
                    ', error =',erreur,',',errorNeum
         write(*,*) ' '
      elseif (erreur.gt.1.0d10) then
         write(*,*) ' '
         write(*,7) 'DIVERGE SOR-GPU!: iters =',iter, &
                    ', error =',erreur,',',errorNeum
         write(*,*) ' '
      elseif(iter.gt.MaxIters) then
         write(*,*) ' '
         write(*,7) 'Non-conv.SOR-GPU: iters =',iter, &
                    ', error =',erreur,',',errorNeum
         write(*,*) ' '
         goto 2000
      else
         if (0.eq.mod(iter,1000)) then
            write(*,7),'iterCUDA =',iter,', error=',erreur,',',errorNeum
         endif
         goto 500
      endif    
      
!      ________________________________________________________
!     |                                                        |
!     |                    Final of iterations                 |
!     |________________________________________________________|

      call cpu_time(finishT)
      timeTotal = (finishT-startT) 
#     ifdef KeyDisplay 
      print*,' >>>>>> Time of CUDA only iterations=',timeTotal
#     endif

7     format(t9,a26,i5,a9,e10.3,a2,e10.3) 

2000  continue
  
!      ________________________________________________________
!     |                                                        |
!     |           Transfer solution from device to host        |
!     |________________________________________________________|

!     ==========================================
!     CUDA: Transfer data from device to host
      phi_vec = d_phi 
      do i=1,N_CELL
         do k=1,NZ
            s = k+(i-1)*NZ
            phi(i,k) = phi_vec(s)  
         enddo
      enddo
!     ==========================================      
           
!      ________________________________________________________
!     |                                                        |
!     |                      Deallocate                        |
!     |________________________________________________________|

      deallocate(h_erreur(grid%x,grid%y))
!     ---------------------
      deallocate(phi_vec(N_CELL*NZ))
      deallocate(phiv_vec(N_CELL*(NZ-1)))
      deallocate(rhs_vec(N_CELL*NZ))
!     ---------------------      
      deallocate(fun_vec(N_CELL*NZ))
      deallocate(funv_vec(N_CELL*(NZ-1)))
!     ---------------------
      deallocate(dfdn_vec(N_CELL*NZ))
      deallocate(dfdnv_vec(N_CELL*(NZ-1)))
                    
      return
      end       

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                               KERNELS                               !
!---------------------------------------------------------------------!
!                     -  kernel_cal_phiJacobi                         !
!                     -  kernel_cal_phiJSOR                           !
!                     -  kernel_cal_phiMSOR                           !
!                     -  kernel_cal_phiBC                             !
!                     -  kernel_cal_phiV                              !
!                     -  kernel_cal_PDJacobi                          !
!                     -  kernel_cal_PDMCSOR                           !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


!=====================================================================!
!                                                                     !
!               KERNEL 1: kernel_cal_phi  (cell-center)               !
!                                                                     !
!=====================================================================!

!     ================================================================
!     KERNEL 1.1: Jacobi 
!     ================================================================

#     ifdef KeyJacobi_cuda
      attributes(global) subroutine kernel_cal_phiJacobi(N_CELL0,NZ,relaxSOR)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s,blockSize
      real*8  :: residu,dot
      integer :: jc1,jc2,jc3,jv1,jv2,jv3    
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ 
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      

      !d_phiNew(index) = d_phi(index)
      if ((idx.le.(N_CELL0)).and.(idy.lt.NZ).and.(idy.gt.1)) then
         jc1 = idy + (d_No_cp(idx,1)-1)*NZ
         jc2 = idy + (d_No_cp(idx,2)-1)*NZ
         jc3 = idy + (d_No_cp(idx,3)-1)*NZ      
         jv1 = idy + (d_No_vp(idx,1)-1)*(NZ-1)
         jv2 = idy + (d_No_vp(idx,2)-1)*(NZ-1)
         jv3 = idy + (d_No_vp(idx,3)-1)*(NZ-1)
         residu = d_rhs(index)                 & 
                - d_Bm(1,index)*d_phiv(jv1)    &
                - d_Bm(2,index)*d_phiv(jv2)    &
                - d_Bm(3,index)*d_phiv(jv3)    &
                - d_Bm(4,index)*d_phiv(jv1-1)  &
                - d_Bm(5,index)*d_phiv(jv2-1)  &
                - d_Bm(6,index)*d_phiv(jv3-1)  &
                - d_Am(1,index)*d_phi(jc1)     &
                - d_Am(2,index)*d_phi(jc2)     &
                - d_Am(3,index)*d_phi(jc3)     & 
                - d_Am(4,index)*d_phi(index+1) &
                - d_Am(5,index)*d_phi(index-1) & 
                - d_phi(index) 
         d_phiNew(index) = d_phi(index) + relaxSOR*residu
         s_erreur(s_index) = abs(residu)
      endif
!     __________________________________________________
!     Error calculation & assignation

      call syncthreads();
	  blockSize=blockDim%x*blockDim%y
      s = blockSize/2
      do while (s.gt.0)
        if (s_index.le.s) then
            s_erreur(s_index) = s_erreur(s_index) + s_erreur(s_index + s)
        endif
        s = s/2
        call syncthreads();
      enddo
!     ______________
      if (s_index.eq.1) then
        d_erreur(blockIdx%x, blockIdx%y) = s_erreur(1)
      endif      

      end subroutine kernel_cal_phiJacobi
#     endif

!     ================================================================
!     KERNEL 1.2: JSOR      
!     ================================================================
     
#     ifdef KeyJSOR_cuda
      attributes(global) subroutine kernel_cal_phiJSOR(N_CELL0,NZ,relaxSOR)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s
      real*8  :: residu
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ 
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      

      d_phiNew(index) = d_phi(index)
      if ((idx.le.(N_CELL0)).and.(idy.lt.NZ).and.(idy.gt.1)) then
         jc1 = idy + (d_No_cp(idx,1)-1)*NZ
         jc2 = idy + (d_No_cp(idx,2)-1)*NZ
         jc3 = idy + (d_No_cp(idx,3)-1)*NZ      
         jv1 = idy + (d_No_vp(idx,1)-1)*(NZ-1)
         jv2 = idy + (d_No_vp(idx,2)-1)*(NZ-1)
         jv3 = idy + (d_No_vp(idx,3)-1)*(NZ-1)          
         residu = d_rhs(index)                 & 
                - d_Bm(1,index)*d_phiv(jv1)    &
                - d_Bm(2,index)*d_phiv(jv2)    &
                - d_Bm(3,index)*d_phiv(jv3)    &
                - d_Bm(4,index)*d_phiv(jv1-1)  &
                - d_Bm(5,index)*d_phiv(jv2-1)  &
                - d_Bm(6,index)*d_phiv(jv3-1)  &
                - d_Am(1,index)*d_phi(jc1)     &
                - d_Am(2,index)*d_phi(jc2)     &
                - d_Am(3,index)*d_phi(jc3)     & 
                - d_Am(4,index)*d_phi(index+1) &
                - d_Am(5,index)*d_phi(index-1) & 
                - d_phi(index)      
         d_phi(index) = d_phi(index) + relaxSOR*residu
         s_erreur(s_index) = abs(residu)
      endif
!     __________________________________________________
!     Error calculation & assignation

      call syncthreads();
      s = blockDim%x*blockDim%y/2
      do while (s.gt.0)
        if (s_index.le.s) then
            s_erreur(s_index) = s_erreur(s_index) + s_erreur(s_index + s)
        endif
        s = s/2
        call syncthreads();
      enddo
!     ______________
      if (s_index.eq.1) then
        d_erreur(blockIdx%x, blockIdx%y) = s_erreur(1)
      endif      

      end subroutine kernel_cal_phiJSOR
#     endif

!     ================================================================
!     KERNEL 1.3: MultiColorSOR
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR(N_CELL0,NZ,relaxSOR,NumColor)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ      
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      
      if (idx.le.N_CELL0) then
      if (d_Color(idx).eq.NumColor) then
      if ((idy.gt.1).and.(idy.lt.NZ)) then
         jc1 = idy + (d_No_cp(idx,1)-1)*NZ
         jc2 = idy + (d_No_cp(idx,2)-1)*NZ
         jc3 = idy + (d_No_cp(idx,3)-1)*NZ      
         jv1 = idy + (d_No_vp(idx,1)-1)*(NZ-1)
         jv2 = idy + (d_No_vp(idx,2)-1)*(NZ-1)
         jv3 = idy + (d_No_vp(idx,3)-1)*(NZ-1)          
         residu = d_rhs(index)                 & 
                - d_Bm(1,index)*d_phiv(jv1)    &
                - d_Bm(2,index)*d_phiv(jv2)    &
                - d_Bm(3,index)*d_phiv(jv3)    &
                - d_Bm(4,index)*d_phiv(jv1-1)  &
                - d_Bm(5,index)*d_phiv(jv2-1)  &
                - d_Bm(6,index)*d_phiv(jv3-1)  &
                - d_Am(1,index)*d_phi(jc1)     &
                - d_Am(2,index)*d_phi(jc2)     &
                - d_Am(3,index)*d_phi(jc3)     & 
                - d_Am(4,index)*d_phi(index+1) &
                - d_Am(5,index)*d_phi(index-1) & 
                - d_phi(index)      
         d_phi(index) = d_phi(index) + relaxSOR*residu
         s_erreur(s_index) = abs(residu) 
      endif
      endif
      endif
!     __________________________________________________
!     Error calculation & assignation

      call syncthreads();
      s = blockDim%x*blockDim%y/2
      do while (s.gt.0)
        if (s_index.le.s) then
            s_erreur(s_index) = s_erreur(s_index) + s_erreur(s_index + s)
        endif
        s = s/2
        call syncthreads();
      enddo
!     ______________
      if (s_index.eq.1) then
        d_erreur(blockIdx%x, blockIdx%y) = s_erreur(1)
      endif      
      
      end subroutine kernel_cal_phiMSOR
#     endif 

!=====================================================================!
!                                                                     !
!             KERNEL 2: kernel_cal_phiBC  (cell-center)               !
!                                                                     !
!=====================================================================!

!     ================================================================
!     KERNEL 2.1:  DIRICHLET  
!     ================================================================

#     ifdef Key_DirichletBCp

      attributes(global) subroutine kernel_cal_phiDirichletBC(N_CELL0,N_CELL,NZ)
      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: N_CELL
      integer, value :: NZ
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s,jc1
      real*8  :: residu
!     __________________________________________________
!     Index
      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ
!     __________________________________________________
!     d_phiBC calculation (vertical)     
      if ((idx.le.N_CELL0).and.(idy.eq.1)) then      
          d_phi(index) = 2.0d0*d_funfB(index)-d_phi(index+1)
      endif    
      if ((idx.le.N_CELL0).and.(idy.eq.NZ)) then
          d_phi(index) = 2.0d0*d_funfB(index)-d_phi(index-1)         
      endif
!     __________________________________________________
!     d_phiBC calculation (horizontal)                
      call syncthreads();
      if (idx.gt.N_CELL0) then         
	     jc1 = idy + (d_No_cp(idx,1)-1)*NZ 
         d_phi(index) = 2.0d0*d_funfB(index)-d_phi(jc1)                  
      endif
      end subroutine kernel_cal_phiDirichletBC
      
#     endif

!     ================================================================
!     KERNEL 2.2:  NEUMANN  
!     ================================================================

#     ifdef Key_NeumannBCp

      attributes(global) subroutine kernel_cal_phiNeumannBC(N_CELL0,N_CELL,NZ)
      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: N_CELL
      integer, value :: NZ
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s,jc1,nc
      real*8  :: residu        
!     __________________________________________________
!     Index
      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ 
!     __________________________________________________
!     d_phiBC calculation (vertical)
      if ((idx.le.N_CELL0).and.(idy.eq.1)) then
          d_phi(index) = d_dlTB(2)*d_fundfdnB(index) + d_phi(index+1)
      endif    
      if ((idx.le.N_CELL0).and.(idy.eq.NZ)) then
          d_phi(index) = d_dlTB(1)*d_fundfdnB(index) + d_phi(index-1)       
      endif      
!     __________________________________________________
!     d_phiBC calculation (horizontal)                
      call syncthreads();
      if (idx.gt.N_CELL0) then         
	     jc1 = idy + (d_No_cp(idx,1)-1)*NZ
	     nc  = d_No_cp(idx,1)
         d_phi(index) = 2.0d0*d_dlCE(nc)*d_fundfdnB(index) + d_phi(jc1)                  
      endif
      end subroutine kernel_cal_phiNeumannBC
            
#     endif

!     ================================================================
!     KERNEL 2.3:  MIX: NEUMANN + DIRICHLET(TOP)  
!     ================================================================

#     ifdef Key_MixBCp

      attributes(global) subroutine kernel_cal_phiMixBC(N_CELL0,N_CELL,NZ)
      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: N_CELL
      integer, value :: NZ
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s,jc1
      real*8  :: residu        
!     __________________________________________________
!     Index
      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ     
!     __________________________________________________
!     d_phiBC calculation (vertical)
      if ((idx.le.N_CELL0).and.(idy.eq.1)) then
          d_phi(index) = d_dlTB(2)*d_fundfdnB(index) + d_phi(index+1)
      endif    
      if ((idx.le.N_CELL0).and.(idy.eq.NZ)) then
          !d_phi(index) = d_dlTB(1)*d_fundfdnB(index) + d_phi(index-1)
          d_phi(index) = 2.0d0*d_funfB(index)-d_phi(index-1) !<<< Dirichlet  at Top        
      endif      
!     __________________________________________________
!     d_phiBC calculation (horizontal)                
      call syncthreads();
      if (idx.gt.N_CELL0) then         
	     jc1 = idy + (d_No_cp(idx,1)-1)*NZ 
         d_phi(index) = 2.0d0*d_dlCE(idx)*d_fundfdnB(index) + d_phi(jc1)                  
      endif
      end subroutine kernel_cal_phiMixBC
            
#     endif
            
!=====================================================================!
!                                                                     !
!                     KERNEL 3: kernel_cal_phiV (vertex)              !
!                                                                     !
!=====================================================================!

      attributes(global) subroutine kernel_cal_phiV(MaxNumC,NZ)

      implicit none
!     ___________________
!     Variables
      integer, value :: MaxNumC
      integer, value :: NZ
!     ___________________
!     Local variables
      integer :: idx,idy,index,indexv,s_index,s
      real*8  :: som
      integer :: j,jc1,nc

!     -------------------------------------------------
!     Index: [blocks = vertex x (NZ-1)]

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index  =  idy + (idx-1)*(NZ)
      indexv =  idy + (idx-1)*(NZ-1)
!     -------------------------------------------------
!     phivC

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y

      do j=1,MaxNumC
         jc1 = idy + (d_Surro(idx,j)-1)*NZ
         som = som + d_Weigh(idx,j)*d_phi(jc1)
      enddo
      d_phivC(index) = som       
!     -------------------------------------------------
!     Interpolation       

      call syncthreads();
      if ((idy.ge.1).and.(idy.le.NZ-1)) then              
         d_phiv(indexv) = d_dzB(idy)*d_phivC(index+1) &
                         +d_dzT(idy)*d_phivC(index)              
      endif                         
!     -------------------------------------------------
!      Boundary conditions vertex (Dirichlet)

#     ifdef Key_DirichletBCp
         call syncthreads();
         if (idy.eq.1) then      
            d_phiv(indexv) = d_funfBv(indexv)      
         elseif (idy.eq.NZ-1) then     
            d_phiv(indexv) = d_funfBv(indexv)                 
         endif
      
         call syncthreads();
         if (d_nbev(idx).ne.0) then
             d_phiv(indexv) = d_funfBv(indexv)
         endif
#     endif

!     -------------------------------------------------
!      Boundary conditions vertex (Mix: Neumann + Dirichlet(Top))

#     ifdef Key_MixBCp
         call syncthreads();
         if (idy.eq.NZ-1) then     
            d_phiv(indexv) = d_funfBv(indexv)                 
         endif
#     endif       

      end subroutine kernel_cal_phiV
 
!=====================================================================!
!                                                                     !
!                        KERNEL 4: PD-MultiColorSOR                   !
!                                                                     !
!=====================================================================!

!     ================================================================
!     KERNEL 5.1: PD-Jacobi 
!     ================================================================

      attributes(global) subroutine kernel_cal_PDJacobi(N_CELL0,NZ,relaxSOR)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s
      real*8  :: residu
      integer :: jc1,jc2,jc3   
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ 
      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      

      if ((idx.le.(N_CELL0)).and.(idy.lt.NZ).and.(idy.gt.1)) then
         jc1 = idy + (d_No_cp(idx,1)-1)*NZ
         jc2 = idy + (d_No_cp(idx,2)-1)*NZ
         jc3 = idy + (d_No_cp(idx,3)-1)*NZ               
         residu = d_rhs(index)                 & 
                - d_Am(1,index)*d_phi(jc1)     &
                - d_Am(2,index)*d_phi(jc2)     &
                - d_Am(3,index)*d_phi(jc3)     & 
                - d_Am(4,index)*d_phi(index+1) &
                - d_Am(5,index)*d_phi(index-1) & 
                - d_phi(index)      
         d_phiNew(index) = d_phi(index) + relaxSOR*residu
         s_erreur(s_index) = abs(residu)
      endif
!     __________________________________________________
!     Error calculation & assignation

      call syncthreads();
      s = blockDim%x*blockDim%y/2
      do while (s.gt.0)
        if (s_index.le.s) then
            s_erreur(s_index) = s_erreur(s_index) + s_erreur(s_index + s)
        endif
        s = s/2
        call syncthreads();
      enddo
      if (s_index.eq.1) then
        d_erreur(blockIdx%x, blockIdx%y) = s_erreur(1)
      endif      

      end subroutine kernel_cal_PDJacobi


!     ================================================================
!     KERNEL 5.2: PD-MCSOR      
!     ================================================================

      attributes(global) subroutine kernel_cal_PDMCSOR(N_CELL0,NZ,relaxSOR,NumColor)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3  
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index = idy + (idx-1)*NZ      
      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      
      if ((idx.le.(N_CELL0)).and.(idy.lt.NZ).and.(idy.gt.1)) then
      if (d_Color(idx).eq.NumColor) then    
         jc1 = idy + (d_No_cp(idx,1)-1)*NZ
         jc2 = idy + (d_No_cp(idx,2)-1)*NZ
         jc3 = idy + (d_No_cp(idx,3)-1)*NZ              
         residu = d_rhs(index)                 & 
                - d_Am(1,index)*d_phi(jc1)     &
                - d_Am(2,index)*d_phi(jc2)     &
                - d_Am(3,index)*d_phi(jc3)     & 
                - d_Am(4,index)*d_phi(index+1) &
                - d_Am(5,index)*d_phi(index-1) & 
                - d_phi(index)      
         d_phi(index) = d_phi(index) + relaxSOR*residu
         s_erreur(s_index) = abs(residu) 
      endif
      endif
!     __________________________________________________
!     Error calculation & assignation

      call syncthreads();
      s = blockDim%x*blockDim%y/2
      do while (s.gt.0)
        if (s_index.le.s) then
            s_erreur(s_index) = s_erreur(s_index) + s_erreur(s_index + s)
        endif
        s = s/2
        call syncthreads();
      enddo
!     ______________
      if (s_index.eq.1) then
        d_erreur(blockIdx%x, blockIdx%y) = s_erreur(1)
      endif      

      end subroutine kernel_cal_PDMCSOR    
      
      END MODULE parallelCUDA

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          END OF MODULE CUDA                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
