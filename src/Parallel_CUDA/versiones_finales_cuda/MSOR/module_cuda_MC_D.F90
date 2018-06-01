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

      real*8, dimension(:),   device, allocatable :: d_phi
      real*8, dimension(:),   device, allocatable :: d_phiNew
      real*8, dimension(:),   device, allocatable :: d_rhs                          
      real*8, dimension(:,:), device, allocatable :: d_Am
      real*8, dimension(:),   device, allocatable :: d_phiv      
      real*8, dimension(:,:), device, allocatable :: d_Bm                       
      real*8, dimension(:),   device, allocatable :: d_funfB
      real*8, dimension(:),   device, allocatable :: d_funfBv    
      integer,dimension(:,:), device, allocatable :: d_No_cp
      integer,dimension(:,:), device, allocatable :: d_No_vp       
      integer,dimension(:,:), device, allocatable :: d_Surro       
      real*8, dimension(:,:), device, allocatable :: d_Weigh       
      real*8, dimension(:),   device, allocatable :: d_dzT
      real*8, dimension(:),   device, allocatable :: d_dzB 
      integer,dimension(:),   device, allocatable :: d_nbev      
      real*8, dimension(:,:), device, allocatable :: d_erreur
	  real*8, dimension(:,:), device, allocatable :: d_erreurP
	  real*8, dimension(:,:), device, allocatable :: d_erreur_c0P
	  real*8, dimension(:,:), device, allocatable :: d_erreur_c1P
	  real*8, dimension(:,:), device, allocatable :: d_erreur_c2P               
      real*8, dimension(:),   device, allocatable :: d_phivC 
      integer,dimension(:),   device, allocatable :: d_Color 
	  integer,dimension(:), device, allocatable :: d_IndexColor

	  real*8, dimension(:), device, allocatable :: d_erreur_AllColors

!#     ifdef KeyAuxMCSOR
      integer,dimension(:), allocatable :: IniColor
      integer,dimension(:), allocatable :: FinColor
      integer,dimension(:), allocatable :: IndexColor
!#     endif
      

!     =========================================================
!     Threads and blocks distribution 
      integer :: THREADSPBLOCK_X,THREADSPBLOCK_Y 
      type(dim3) :: grid,gridv,tBlock,gridP,gridI
	  type(dim3) :: grid_c0P,grid_c1P,grid_c2P
      integer :: shared_mem
      integer :: N_C0,N_C,N_V,N_V0,MaxNumC,midP,midI,Inie_c,Fine_c 
	  integer :: Elem_c0,Elem_c1,Elem_c2,Tam_erreur_AllColors
	  integer :: Elem_grid_c0P,Elem_grid_c1P,Elem_grid_c2P,Elem_grid_c0P_c1P
                          
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
 
      integer :: NumHc,NumZc,NumHv,NumZv,NumZcP,NumZcI
	  integer :: NumHc_0,NumHc_1,NumHc_2,m,s
	  integer :: pot_i
	  real*8  :: pot_r
       
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


!     Necesitamos THREADSPBLOCK_X * THREADSPBLOCK_Y sea potencia de 2
	  pot_r=log(real(THREADSPBLOCK_X * THREADSPBLOCK_Y))/log(2.0)
	  pot_i=ceiling(pot_r)
	  THREADSPBLOCK_X=ceiling(real(2**pot_i)/THREADSPBLOCK_Y)
      if (THREADSPBLOCK_X*THREADSPBLOCK_Y.gt.1024) then
         print*,'Number of Threads per block exceeds 1024, please correct it.'
         stop
      endif

!     ============================
!     CUDA: Domain decomposition
      
      tBlock = dim3(THREADSPBLOCK_X,THREADSPBLOCK_Y,1)      
      NumHc = ceiling(real(N_CELL)/tBlock%x)
      NumZc = ceiling(real(NZ)/tBlock%y)
      NumHv = ceiling(real(N_VERT)/tBlock%x)
      NumZv = ceiling(real(NZ-1)/tBlock%y)
      grid   = dim3(NumHc,NumZc,1)
      gridv  = dim3(NumHV,NumZv,1)

!     ============================
!     Grid Configuration: Even and Odd

	  !Vertical  
      midP = floor((NZ-1)/2.0d0)  !Mitad Pares
      if (0.eq.mod(midP,2)) then
         midI = midP -1           !Mitad Impares
      else
         midI = midP
      endif

	  NumZcP = ceiling(real(midP)/tBlock%y)
	  NumZcI = ceiling(real(midI)/tBlock%y)
	  gridP   = dim3(NumHc,NumZcP,1)
	  gridI   = dim3(NumHc,NumZcI,1)
!     ============================

!      ________________________________________________________
!     |                                                        |
!     |               Reordering elements by color             |
!     |________________________________________________________|
      
	  allocate(IniColor(0:3)) 
	  allocate(FinColor(0:3))
	  allocate(IndexColor(N_CELL0))

         m = 0
         DO s=0,N_COLOR-1
            IniColor(s) = m+1
            do i=1,N_CELL0
               if (ColorCell(i).eq.s) then
                  m = m+1
                  IndexColor(m) = i
               endif
            enddo
            FinColor(s) = m
         ENDDO
     

	  print*,IniColor(0)
	  print*,IniColor(1)
	  print*,IniColor(2)

!     ============================
!     Grid Configuration: Per cells
      Elem_c0=FinColor(0)-IniColor(0)+1
	  Elem_c1=FinColor(1)-IniColor(1)+1
	  Elem_c2=FinColor(2)-IniColor(2)+1
      NumHc_0 = ceiling(real(Elem_c0)/tBlock%x)
	  NumHc_1 = ceiling(real(Elem_c1)/tBlock%x)
	  NumHc_2 = ceiling(real(Elem_c2)/tBlock%x)
	  grid_c0P=dim3(NumHc_0,NumZcP,1)
	  grid_c1P=dim3(NumHc_1,NumZcP,1)
	  grid_c2P=dim3(NumHc_2,NumZcP,1)

	  Elem_grid_c0P=grid_c0P%x * grid_c0P%y
	  Elem_grid_c1P=grid_c1P%x * grid_c1P%y
	  Elem_grid_c2P=grid_c2P%x * grid_c2P%y
	  Elem_grid_c0P_c1P=Elem_grid_c0P+ Elem_grid_c1P
	  Tam_erreur_AllColors= Elem_grid_c0P_c1P + Elem_grid_c2P

	  print*,Elem_c0
	  print*,Elem_c1
	  print*,Elem_c2

!     ============================

      
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
         print*,'Shared memory setting for update kernel exceeds 64KB, please correct it.'
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
      allocate(d_Am(5,N_C0))
      allocate(d_Bm(6,N_C0))
      allocate(d_phiv(N_V0)) 
      allocate(d_funfBv(N_V0))
      allocate(d_phivC(N_V))      
      allocate(d_No_cp(N_CELL,3))
      allocate(d_No_vp(N_CELL0,3))                                     
      allocate(d_dzT(NZ-1))
      allocate(d_dzB(NZ-1))
      allocate(d_nbev(N_VERT))                                                                    
      allocate(d_erreur(grid%x,grid%y))   
	  allocate(d_erreurP(gridP%x,gridP%y))
	  allocate(d_erreur_c0P(grid_c0P%x,grid_c0P%y))
	  allocate(d_erreur_c1P(grid_c1P%x,grid_c1P%y))
	  allocate(d_erreur_c2P(grid_c2P%x,grid_c2P%y))    
	  allocate(d_erreur_AllColors(Tam_erreur_AllColors))      
#     ifdef KeyMSOR_cuda           
      allocate(d_Color(N_CELL0))
#     endif  
#     ifdef KeyPDMSOR_cuda           
       allocate(d_Color(N_CELL0))
#     endif   

	  allocate(d_IndexColor(N_CELL0))   

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
      deallocate(d_Am(5,N_C0))
      deallocate(d_Bm(6,N_C0))
      deallocate(d_phiv(N_V0)) 
      deallocate(d_funfBv(N_V0))
      deallocate(d_phivC(N_V0))           
      deallocate(d_No_cp(N_CELL,3))
      deallocate(d_No_vp(N_CELL0,3))       
      deallocate(d_Surro(N_VERT,MaxNumC))
      deallocate(d_Weigh(N_VERT,MaxNumC))                               
      deallocate(d_dzT(NZ-1))
      deallocate(d_dzB(NZ-1))
      deallocate(d_nbev(N_VERT))            
      deallocate(d_erreur(grid%x,grid%y))
	  deallocate(d_erreurP(gridP%x,gridP%y))    
	  deallocate(d_erreur_c0P(grid_c0P%x,grid_c0P%y))    
	  deallocate(d_erreur_c1P(grid_c1P%x,grid_c1P%y))    
	  deallocate(d_erreur_c2P(grid_c2P%x,grid_c2P%y))    
	  deallocate(d_erreur_AllColors(Tam_erreur_AllColors))      
	         
#     ifdef KeyMSOR_cuda           
      deallocate(d_Color(N_CELL0))
#     endif  
#     ifdef KeyPDMSOR_cuda           
      deallocate(d_Color(N_CELL0))
#     endif      

      deallocate(d_IndexColor(N_CELL0)) 

	  

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
      real*8, dimension(:), allocatable :: dzT_vec,dzB_vec        
      real*8 :: dzB,dzT                        
      
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
     
#     ifdef KeyMSOR_cuda           
         d_Color = ColorCell
#     endif  
#     ifdef KeyPDMSOR_cuda           
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
!      MaxNumC = 1
!      do nv=1,N_VERT
!         MaxNumC = max(MaxNumC,Dimsurrounding(nv))
!      enddo
      if (MaxNumC.gt.9)then
	      print*,'WARNING: Surrounding is bigger than 9!'
		  stop
      endif
      MaxNumC = 9 
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

      subroutine MethSOR_CUDA(phi,phiv,rhs,                        &
                              Am1,Am2,Am3,AmT,AmB,                 &
                              Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                              No_cp,No_vp,nbev,sig,sigv,           &
                              funfB,funfBv)

#     include "cppdefs.h"
      USE geometry
      implicit none 
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------              

!     __________________________________
!     Declaration of variables
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
      
!     __________________________________ 

      real*8, dimension(:,:), allocatable :: h_erreur
	  real*8, dimension(:,:), allocatable :: h_erreurP
	  real*8, dimension(:,:), allocatable :: h_erreur_c0P
	  real*8, dimension(:,:), allocatable :: h_erreur_c1P
	  real*8, dimension(:,:), allocatable :: h_erreur_c2P
	  real*8, dimension(:), allocatable :: h_erreur_AllColors
      real*8, dimension(:),   allocatable :: phi_vec
      real*8, dimension(:),   allocatable :: rhs_vec
      real*8, dimension(:),   allocatable :: fun_vec
      real*8, dimension(:),   allocatable :: phiv_vec
      real*8, dimension(:),   allocatable :: funv_vec                                                                                       
!     __________________________________
!     Declaration of local variables
      integer:: iter,s,ss 
      real*8 :: erreur
!     ----------------------------------------
      real :: startT,finishT,timeTotal                
       
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________| 
 
    
      allocate(phi_vec(N_CELL*NZ))
      allocate(rhs_vec(N_CELL*NZ))
      allocate(fun_vec(N_CELL*NZ))                  
      allocate(phiv_vec(N_CELL*(NZ-1)))
      allocate(funv_vec(N_CELL*(NZ-1)))
      allocate(h_erreur(grid%x,grid%y)) 
	  allocate(h_erreurP(gridP%x,gridP%y))
	  allocate(h_erreur_c0P(grid_c0P%x,grid_c0P%y))
	  allocate(h_erreur_c1P(grid_c1P%x,grid_c1P%y))     
	  allocate(h_erreur_c2P(grid_c2P%x,grid_c2P%y))

	  allocate(h_erreur_AllColors(Tam_erreur_AllColors))
      
!     =========================================================
!     CUDA: Transfer data from host to device      

      do i=1,N_CELL
         do k=1,NZ
            s = k+(i-1)*NZ
            phi_vec(s) = phi(i,k)            
            rhs_vec(s) = rhs(i,k)
            fun_vec(s) = funfB(i,k)
         enddo
      enddo                             
      do nv=1,N_VERT
         do k=1,NZ-1
            s = k+(nv-1)*(NZ-1)           
            phiv_vec(s) = phiv(nv,k)
            funv_vec(s) = funfBv(nv,k)
         enddo
      enddo

      d_phi     = phi_vec
      d_phiNew  = phi_vec
      d_rhs     = rhs_vec
      d_funfB   = fun_vec
      d_phiv    = funv_vec
      d_funfBv  = funv_vec
      
!     =========================================================
!     CUDA: Transfer matrices & geometry from host to device  

      call transfergeometry_cuda(No_cp,No_vp,nbev,sig,sigv)
      call transferInterpo_cuda      
      call transfermatrixA_cuda(Am1,Am2,Am3,AmT,AmB)
      call transfermatrixB_cuda(Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B)

!     =========================================================
!     CUDA: Transfer IndexColor
      print*,'IndexColor...'
      d_IndexColor=IndexColor
	       
!      ________________________________________________________
!     |                                                        |
!     |                    Begin of iterations                 |
!     |________________________________________________________| 

      call cpu_time(startT)
      
      iter=0
500   continue
      iter=iter+1

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
!         do s=0,N_COLOR-1
!            call kernel_cal_phiMSOR<<<grid,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s)            
!            h_erreur = d_erreur
!            erreur = erreur + sum(h_erreur)
!         enddo
		
!		 do s=0,N_COLOR-1
!		    !Even
!            call kernel_cal_phiMSOR_Even<<<gridP,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s,midP)            
!            h_erreurP = d_erreurP
!            erreur = erreur + sum(h_erreurP)

!			!Odd
!			call kernel_cal_phiMSOR_Odd<<<gridI,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s,midI)            
!            h_erreurP = d_erreurP
!            erreur = erreur + sum(h_erreurP)
!
!         enddo

!         do s=0,N_COLOR-1
!            call kernel_cal_phiMSOR_EvenOdd<<<gridP,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s,midP,midI)            
!            h_erreurP = d_erreurP
!            erreur = erreur + sum(h_erreurP)
!         enddo
		 
		 !do s=0,N_COLOR-1
		    Inie_c=IniColor(0)
			call kernel_cal_phiMSOR_EvenOdd_Color0<<<grid_c0P,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s,midP,midI,Inie_c,Elem_c0)            
            !h_erreur_c0P = d_erreur_c0P
            !erreur = erreur + sum(h_erreur_c0P)

			Inie_c=IniColor(1)
			call kernel_cal_phiMSOR_EvenOdd_Color1<<<grid_c1P,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s,midP,midI,Inie_c,Elem_c1,Elem_grid_c0P)            
            !h_erreur_c1P = d_erreur_c1P
            !erreur = erreur + sum(h_erreur_c1P)

			Inie_c=IniColor(2)
			call kernel_cal_phiMSOR_EvenOdd_Color2<<<grid_c2P,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s,midP,midI,Inie_c,Elem_c2,Elem_grid_c0P_c1P)            
            !h_erreur_c2P = d_erreur_c2P
            !erreur = erreur + sum(h_erreur_c2P)

			h_erreur_AllColors = d_erreur_AllColors
            erreur = sum(h_erreur_AllColors)

         !enddo

#     endif        

!     ===========================
!     KERNEL: Boundary conditions 
!     ===========================

      call kernel_cal_phiBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)

!     ===========================            
!     KERNEL: Interpolation 
!     =========================== 
                 
      call kernel_cal_phiV<<<gridv,tBlock,shared_mem>>>(MaxNumC,NZ) 

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
    
      if (erreur.lt.eps) then
         write(*,*) ' '
         write(*,7) 'Solution MethSOR-CUDA3D: iters =',iter,', error =',erreur
         write(*,*) ' '
      elseif (erreur.gt.1.0d5) then
         write(*,*) ' '
         write(*,7) 'DIVERGENCE !!!!: iters =',iter,', error =',erreur
         write(*,*) ' '
      elseif(iter.gt.MaxIters) then
         write(*,*) ' '
         write(*,7) 'Non-convergence: iters =',iter,', error =',erreur
         write(*,*) ' '
      else
         if (0.eq.mod(iter,1000)) print*, 'iterCUDA =',iter,'Error=',erreur
         goto 500
      endif    
      
!      ________________________________________________________
!     |                                                        |
!     |                    Final of iterations                 |
!     |________________________________________________________| 

      call cpu_time(finishT)
      timeTotal = (finishT-startT)  
      print*,' >>>>>> Time of CUDA only iterations=',timeTotal     

7     format(t10,a32,i8,a9,e10.3) 

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

      deallocate(phi_vec(N_CELL*NZ))
      deallocate(rhs_vec(N_CELL*NZ))
      deallocate(fun_vec(N_CELL*NZ))              
      deallocate(phiv_vec(N_CELL*(NZ-1)))
      deallocate(funv_vec(N_CELL*(NZ-1)))      
      deallocate(h_erreur(grid%x,grid%y))
	  deallocate(h_erreurP(gridP%x,gridP%y))
	  deallocate(h_erreur_c0P(grid_c0P%x,grid_c0P%y))
	  deallocate(h_erreur_c1P(grid_c1P%x,grid_c1P%y))     
	  deallocate(h_erreur_c2P(grid_c2P%x,grid_c2P%y))  
	  deallocate(h_erreur_AllColors(Tam_erreur_AllColors)) 
	  
      deallocate(IniColor(0:3))
	  deallocate(FinColor(0:3))
	  deallocate(IndexColor(N_CELL0)) 
                    
      return
      end       

!=====================================================================!
!                                                                     !
!               KERNEL 1: kernel_cal_phi  (cell-center)               !
!                                                                     !
!=====================================================================!

!     ================================================================
!     KERNEL 1: Jacobi 
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
!     Old
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
      if ((idy.gt.1).and.(idy.lt.NZ).and.(0.eq.mod(idy,2))) then
	     
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
      call syncthreads();
      if ((idy.gt.1).and.(idy.lt.NZ).and.(1.eq.mod(idy,2))) then   
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

!     ================================================================
!     KERNEL 1.4: MultiColorSOR Even
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR_Even(N_CELL0,NZ,relaxSOR,NumColor,midP)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor,midP
!     ___________________
!     Local variables
      integer :: idx,idy,idyl,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idyl   = (blockidx%y-1)*blockdim%y + threadidx%y
	  idy   = idyl*2 !Even
      index = idy + (idx-1)*NZ      
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      

      if (idx.le.N_CELL0) then
      if (d_Color(idx).eq.NumColor) then 
      if (idyl.le.midP) then
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
      
      end subroutine kernel_cal_phiMSOR_Even
#     endif 

!     ================================================================
!     KERNEL 1.5: MultiColorSOR Odd
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR_Odd(N_CELL0,NZ,relaxSOR,NumColor,midI)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor,midI
!     ___________________
!     Local variables
      integer :: idx,idy,idyl,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idyl   = (blockidx%y-1)*blockdim%y + threadidx%y
	  idy   = idyl*2 +1 !Odd
      index = idy + (idx-1)*NZ      
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      

      if (idx.le.N_CELL0) then
      if (d_Color(idx).eq.NumColor) then 
      if (idyl.le.midI) then
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
      
      end subroutine kernel_cal_phiMSOR_Odd
#     endif 

!     ================================================================
!     KERNEL 1.6: MultiColorSOR Even and Odd
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR_EvenOdd(N_CELL0,NZ,relaxSOR,NumColor,midP,midI)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor,midP,midI
!     ___________________
!     Local variables
      integer :: idx,idy,idyl,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idyl   = (blockidx%y-1)*blockdim%y + threadidx%y     
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      

      if (idx.le.N_CELL0) then
      if (d_Color(idx).eq.NumColor) then 
      !if ((idy.gt.1).and.(idy.lt.NZ).and.(0.eq.mod(idy,2))) then
	  if (idyl.le.midP) then   
		 idy=idyl*2;
         index = idy + (idx-1)*NZ 

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
      call syncthreads();
      !if ((idy.gt.1).and.(idy.lt.NZ).and.(1.eq.mod(idy,2))) then   
	  if (idyl.le.midI) then
		 idy=idyl*2+1;
         index = idy + (idx-1)*NZ 
	     
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
         s_erreur(s_index) = s_erreur(s_index) + abs(residu) 
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
        d_erreurP(blockIdx%x, blockIdx%y) = s_erreur(1)
      endif      
      
      end subroutine kernel_cal_phiMSOR_EvenOdd
#     endif 

!     ================================================================
!     KERNEL 1.7.1: MultiColorSOR Even, Odd and Color
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR_EvenOdd_Color0(N_CELL0,NZ,relaxSOR,NumColor,midP,midI,Inie_c,Elem_c0)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor,midP,midI,Inie_c,Elem_c0
!     ___________________
!     Local variables
      integer :: idx,idxl,idy,idyl,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idxl   = (blockidx%x-1)*blockdim%x + threadidx%x
      idyl   = (blockidx%y-1)*blockdim%y + threadidx%y     
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      

      if (idxl.le.Elem_c0) then
      !if (d_Color(idx).eq.NumColor) then 
	  if (idyl.le.midP) then   
		 idy=idyl*2;
		 idx=d_IndexColor(Inie_c+idxl-1)
         index = idy + (idx-1)*NZ 

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
      call syncthreads();
	  if (idyl.le.midI) then
		 idy=idyl*2+1;
		 idx=d_IndexColor(Inie_c+idxl-1)
         index = idy + (idx-1)*NZ 
	     
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
         s_erreur(s_index) = s_erreur(s_index) + abs(residu) 
      endif
      !endif
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
        !d_erreur_c0P(blockIdx%x, blockIdx%y) = s_erreur(1)
		d_erreur_AllColors((blockIdx%x-1) * gridDim%y + blockIdx%y) = s_erreur(1)
      endif      
      
      end subroutine kernel_cal_phiMSOR_EvenOdd_Color0
#     endif 

!     ================================================================
!     KERNEL 1.7.2: MultiColorSOR Even, Odd and Color
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR_EvenOdd_Color1(N_CELL0,NZ,relaxSOR,NumColor,midP,midI,Inie_c,Elem_c1,Elem_grid_c0P)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor,midP,midI,Inie_c,Elem_c1,Elem_grid_c0P
!     ___________________
!     Local variables
      integer :: idx,idxl,idy,idyl,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idxl   = (blockidx%x-1)*blockdim%x + threadidx%x
      idyl   = (blockidx%y-1)*blockdim%y + threadidx%y     
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      

      if (idxl.le.Elem_c1) then
      !if (d_Color(idx).eq.NumColor) then 
	  if (idyl.le.midP) then   
		 idy=idyl*2;
		 idx=d_IndexColor(Inie_c+idxl-1)
         index = idy + (idx-1)*NZ 

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
      call syncthreads();
	  if (idyl.le.midI) then
		 idy=idyl*2+1;
		 idx=d_IndexColor(Inie_c+idxl-1)
         index = idy + (idx-1)*NZ 
	     
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
         s_erreur(s_index) = s_erreur(s_index) + abs(residu) 
      endif
      !endif
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
        !d_erreur_c1P(blockIdx%x, blockIdx%y) = s_erreur(1)
		d_erreur_AllColors((blockIdx%x-1) * gridDim%y + blockIdx%y + Elem_grid_c0P)= s_erreur(1)
      endif      
      
      end subroutine kernel_cal_phiMSOR_EvenOdd_Color1
#     endif 

!     ================================================================
!     KERNEL 1.7.3: MultiColorSOR Even, Odd and Color
!     ================================================================

#     ifdef KeyMSOR_cuda
      attributes(global) subroutine kernel_cal_phiMSOR_EvenOdd_Color2(N_CELL0,NZ,relaxSOR,NumColor,midP,midI,Inie_c,Elem_c2,Elem_grid_c0P_c1P)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: NZ
      real*8,  value :: relaxSOR
      integer, value :: NumColor,midP,midI,Inie_c,Elem_c2,Elem_grid_c0P_c1P
!     ___________________
!     Local variables
      integer :: idx,idxl,idy,idyl,index,s_index,s
      real*8  :: residu,sol
      integer :: jc1,jc2,jc3,jv1,jv2,jv3      
!     ___________________
!     Share variable
      real*8,shared,dimension(*) :: s_erreur

!     __________________________________________________
!     Index

      idxl   = (blockidx%x-1)*blockdim%x + threadidx%x
      idyl   = (blockidx%y-1)*blockdim%y + threadidx%y     
      
!     __________________________________________________
!     Error index

      s_index = (threadidx%x-1)*blockdim%y + threadidx%y
      s_erreur(s_index) = 0.0d0
!     __________________________________________________
!     d_phiNew calculation      
      

      if (idxl.le.Elem_c2) then
      !if (d_Color(idx).eq.NumColor) then 
	  if (idyl.le.midP) then   
		 idy=idyl*2;
		 idx=d_IndexColor(Inie_c+idxl-1)
         index = idy + (idx-1)*NZ 

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
      call syncthreads();
	  if (idyl.le.midI) then
		 idy=idyl*2+1;
		 idx=d_IndexColor(Inie_c+idxl-1)
         index = idy + (idx-1)*NZ 
	     
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
         s_erreur(s_index) = s_erreur(s_index) + abs(residu) 
      endif
      !endif
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
        !d_erreur_c2P(blockIdx%x, blockIdx%y) = s_erreur(1)
		d_erreur_AllColors((blockIdx%x-1) * gridDim%y + blockIdx%y + Elem_grid_c0P_c1P)= s_erreur(1)
      endif      
      
      end subroutine kernel_cal_phiMSOR_EvenOdd_Color2
#     endif 
!=====================================================================!
!                                                                     !
!             KERNEL 3: kernel_cal_phiBC  (cell-center)               !
!                                                                     !
!=====================================================================!

      attributes(global) subroutine kernel_cal_phiBC(N_CELL0,N_CELL,NZ)

      implicit none
!     ___________________
!     Variables
      integer, value :: N_CELL0
      integer, value :: N_CELL
      integer, value :: NZ
!     ___________________
!     Local variables
      integer :: idx,idy,index,s_index,s
      real*8  :: residu
      integer :: jc1          

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
      
      end subroutine kernel_cal_phiBC
            
!=====================================================================!
!                                                                     !
!                     KERNEL 4: kernel_cal_phiV (vertex)              !
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
      integer :: idx,idy,index,indexv
      real*8  :: som
      integer :: j,jc1  

!     -------------------------------------------------
!     Index: [blocks = vertex x (NZ-1)]

      idx   = (blockidx%x-1)*blockdim%x + threadidx%x
      idy   = (blockidx%y-1)*blockdim%y + threadidx%y
      index  =  idy + (idx-1)*(NZ)
      indexv =  idy + (idx-1)*(NZ-1)
!     -------------------------------------------------
!     phivC

!      do j=1,MaxNumC
!         jc1 = idy + (d_Surro(idx,j)-1)*NZ 
!         som = som + d_Weigh(idx,j)*d_phi(jc1) 
!      enddo

      !     Unroll
      som=0.0d0
	  !j=1
	  jc1 = idy + (d_Surro(idx,1)-1)*NZ
	  som = som + d_Weigh(idx,1)*d_phi(jc1)
	  !j=2
	  jc1 = idy + (d_Surro(idx,2)-1)*NZ
	  som = som + d_Weigh(idx,2)*d_phi(jc1)
	  !j=3
	  jc1 = idy + (d_Surro(idx,3)-1)*NZ
	  som = som + d_Weigh(idx,3)*d_phi(jc1)
	  !j=4
	  jc1 = idy + (d_Surro(idx,4)-1)*NZ
	  som = som + d_Weigh(idx,4)*d_phi(jc1)
	  !j=5
	  jc1 = idy + (d_Surro(idx,5)-1)*NZ
	  som = som + d_Weigh(idx,5)*d_phi(jc1)
	  !j=6
	  jc1 = idy + (d_Surro(idx,6)-1)*NZ
	  som = som + d_Weigh(idx,6)*d_phi(jc1)
	  !j=7
	  jc1 = idy + (d_Surro(idx,7)-1)*NZ
	  som = som + d_Weigh(idx,7)*d_phi(jc1)
	  !j=8
	  jc1 = idy + (d_Surro(idx,8)-1)*NZ
	  som = som + d_Weigh(idx,8)*d_phi(jc1)
	  !j=9
	  jc1 = idy + (d_Surro(idx,9)-1)*NZ
	  som = som + d_Weigh(idx,9)*d_phi(jc1)


      d_phivC(index) = som       
!     -------------------------------------------------
!     Interpolation       

      call syncthreads();
      if ((idy.ge.1).and.(idy.le.NZ-1)) then              
         d_phiv(indexv) = d_dzB(idy)*d_phivC(index+1) &
                         +d_dzT(idy)*d_phivC(index)              
      endif                         
!     -------------------------------------------------
!      Boundary conditions vertex

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

      end subroutine kernel_cal_phiV
 
!=====================================================================!
!                                                                     !
!                        KERNEL 5: PD-MultiColorSOR                   !
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