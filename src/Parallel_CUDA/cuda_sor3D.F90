!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!         SOR SOLVERS FOR THE 3D POISSON EQUATION USING CUDA          !
!                               Dic 2017                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE sor3Dcuda(phi,phiv,                    &
                           rhs,Gamx,Gamy,Gamz,          &
                           xc,yc,sig,dsig,No_cp,nbe,    &
                           xv,yv,sigv,dsigv,No_vp,nbev, &
                           Hpr,h,etan,                  &
                           Hprv,hv,etav)                                          

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the Sucessive       !
!    Over-Relaxation technique using a CUDA solver. In this part      !
!    we calculate all those parts that the system solutions already   !
!    needs. For example the calculation of the coefficients AM.       ! 
!    We use different relaxion factors: relaxSOR.                     !
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
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - diffusion3D                ( diffusion3D.F90 )            |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - BCcellcenter3D             ( BCcellcenter3D.F90 )         |  !
!  |   - BCvertex3D                 ( BCvertex3D.F90 )             |  !
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
!     =========  PARALLEL MPI ============
#     ifdef KeyParallel
        USE parallel
#     endif
!     ==================================== 
!     =====================================

!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============
#     ifdef KeyCUDA
	     USE parallelCUDA
#     endif 
!     =====================================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

      USE geometry
      implicit none  
!     -----------------
#     ifdef KeyParallel
!#     include "common.mpf"
#     endif
!     -----------------        
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

      real*8,dimension(:,:) :: phiNew(N_CELL,NZ)
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
      real*8,dimension(:,:) :: Nm(N_CELL0,NZ)
      real*8,dimension(:,:) :: Newrhs(N_CELL0,NZ)
!     --------------------------------------
      integer :: tag,SaveEpsResults
!     ----------------------------------------
      integer :: ii,jj
      real*8  :: x,y,z,zT,zB
      real*8  :: funExamNSp,funExamNSrhsp,NeumanndpdnNS
      real*8  :: nnx,nny,nnz
      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8,dimension(:,:) :: funfBv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: fundfdnB(N_CELL,NZ)
      real*8,dimension(:,:) :: fundfdnBv(N_VERT,NZ-1)
!     ---------------------------------------- 
      integer :: j0,j1,j2,j3,jT,jB,irec
      integer :: iter,jc1,jc2,jc3,jv1,jv2,jv3
      real*8  :: som,residu,errorsys   
!     ----------------------------------------
      real :: start,finish,startT,finishT
      real :: timeTotal,timeOther
      real :: timeBCcc,timeBCvv,timeInter
      real :: timeCUDA
           
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: CUDA-SOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      call cpu_time(startT)
                              
!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                      Initial times                     |
!     |________________________________________________________|

      timeTotal = 0.0d0
      timeOther = 0.0d0
      timeBCcc  = 0.0d0
      timeBCvv  = 0.0d0
      timeInter = 0.0d0
      timeCUDA  = 0.0d0

!      ________________________________________________________
!     |                                                        |
!     |         Special case only testing Poisson              |
!     |________________________________________________________|

#     ifdef KeyTestOnlyPoisson
         do k=1,NZ
            do i=1,N_CELL0 
               x = xc(i)
               y = yc(i) 	
               z = sig(k)*Hpr(i)-h(i)             
               rhs(i,k) = VolPrism(i,k)*Hpr(i)*funExamNSrhsp(x,y,z,time)
            enddo
         enddo
         phi  = 0.0d0
         phiv = 0.0d0
#     endif

!      ________________________________________________________
!     |                                                        |
!     |           Function fB & dfdnB at the boundary          |
!     |________________________________________________________|

      call BCfunction3D(funfB,funfBv,fundfdnB,fundfdnBv, &
                        xc,yc,sig,dsig,No_cp,nbe,        &
                        xv,yv,sigv,dsigv,No_vp,nbev,     &
                        Hpr,h,etan,                      &
                        Hprv,hv,etav) 
                        
!      ________________________________________________________
!     |                                                        |
!     |                   Boundary conditions                  |
!     |________________________________________________________|

      call BCpressure3D(phi,phiv,                        &
                        funfB,funfBv,fundfdnB,fundfdnBv, &
                        xc,yc,sig,dsig,No_cp,nbe,        &
                        xv,yv,sigv,dsigv,No_vp,nbev,     &
                        Hpr,h,etan,                      &
                        Hprv,hv,etav,                    &
                        timeBCcc,timeBCvv,timeInter) 

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

!     ________________________________________________________      
!     Diffusion with Dirichlet BC

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             & 
                       Bmv1T,Bmv2T,Bmv3T,                   &
                       Bmv1B,Bmv2B,Bmv3B,                   &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)

!     ________________________________________________________
!    |                                                        |
!    |                    Update coefficients                 |
!    |________________________________________________________|    

      do k=1,NZ
         do i=1,N_CELL0 	
            Am1(i,k)  = Am1(i,k)/Am0(i,k)
            Am2(i,k)  = Am2(i,k)/Am0(i,k)
            Am3(i,k)  = Am3(i,k)/Am0(i,k)
            AmT(i,k)  = AmT(i,k)/Am0(i,k)
            AmB(i,k)  = AmB(i,k)/Am0(i,k)
            Bmv1T(i,k)= Bmv1T(i,k)/Am0(i,k) 
            Bmv2T(i,k)= Bmv2T(i,k)/Am0(i,k) 
            Bmv3T(i,k)= Bmv3T(i,k)/Am0(i,k) 
            Bmv1B(i,k)= Bmv1B(i,k)/Am0(i,k) 
            Bmv2B(i,k)= Bmv2B(i,k)/Am0(i,k) 
            Bmv3B(i,k)= Bmv3B(i,k)/Am0(i,k) 
            rhs(i,k)  = rhs(i,k)/Am0(i,k)
         enddo
      enddo
  
!*********************************************************************!
!                                                                     !
!             Solution of the system (S.0.R.) 3D by CUDA              !
!                                                                     !
!*********************************************************************!      
 
      call cpu_time(start)
      
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============    
#     ifdef KeyCUDA      
         call MethSOR_CUDA(phi,phiv,rhs,                       &
                           Am1,Am2,Am3,AmT,AmB,                &
                           Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B,&
                           No_cp,No_vp,nbev,sig,sigv,          &
                           funfB,funfBv,                       &
                           fundfdnB,fundfdnBv)                                          
!     ======================================
!     ======================================
!     Only if we want to run sequential 
#      else
         call Jacobi_Sequ(phi,phiv,rhs,                        &
                          Am1,Am2,Am3,AmT,AmB,                 &
                          Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                          funfB,funfBv,                        &
                          xv,yv,sigv,dsigv,No_vp,nbev,         &
                          xc,yc,sig,dsig,No_cp,nbe)                                            
#     endif
!     =====================================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
                    
      call cpu_time(finish)
      timeCUDA = finish-start  
  
!      ________________________________________________________
!     |                                                        |
!     |                   Boundary conditions                  |
!     |________________________________________________________|

      call BCpressure3D(phi,phiv,                        &
                        funfB,funfBv,fundfdnB,fundfdnBv, &
                        xc,yc,sig,dsig,No_cp,nbe,        &
                        xv,yv,sigv,dsigv,No_vp,nbev,     &
                        Hpr,h,etan,                      &
                        Hprv,hv,etav,                    &
                        timeBCcc,timeBCvv,timeInter)        
            
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      call cpu_time(finishT)
      timeTotal = finishT-startT  
      
!      ________________________________________________________
!     |                                                        |
!     |                  Distribution of times                 |
!     |________________________________________________________|

#     ifdef KeyDisplay
      timeOther = timeTotal-timeCUDA
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============ 
#     ifdef KeyCUDA 
         print*,' '
         print*,'        ------------------------------------'
         print*,'        Distribution of the time in SOR-CUDA' 
         print*,'        ------------------------------------'
         print*,'        Time SOR CUDA       :  ',timeCUDA
         print*,'        Time other          :  ',timeOther
         print*,'        ____________________________________'
         print*,'        Total Time          :  ',timeTotal
         print*,'        ------------------------------------'
         print*,' '
         !write(7100,8),time,eps,iter,timeTotal
!     ======================================
!     ======================================         
#     else
         print*,' '
         print*,'        ------------------------------------'
         print*,'        Distribution of the time in SOR-Sequ'
         print*,'        ------------------------------------'
         print*,'        Time Jacobi Seq.    :  ',timeCUDA
         print*,'        Time other          :  ',timeOther
         print*,'        ____________________________________'
         print*,'        Total Time          :  ',timeTotal
         print*,'        ------------------------------------'
         print*,' '
         !write(7100,8),time,eps,iter,timeTotal      
#     endif 
!     =====================================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: CUDA-SOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END
       
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        AUXILIAR 1:  Jacobi sequential               !
!                              Sept 2015                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Jacobi_Sequ(phi,phiv,rhs,                        &
                             Am1,Am2,Am3,AmT,AmB,                 &
                             Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                             funfB,funfBv,                        &
                             xv,yv,sigv,dsigv,No_vp,nbev,         &
                             xc,yc,sig,dsig,No_cp,nbe)

      USE geometry
      implicit none
!     -----------------
#     ifdef KeyParallel
#     include "common.mpf"
#     endif
!     -----------------      
      
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
      real*8,dimension(:,:)  :: funfB(N_CELL,NZ)
      real*8,dimension(:,:)  :: funfBv(N_VERT,NZ-1)
      real*8,dimension(:)    :: xv(N_VERT)
      real*8,dimension(:)    :: yv(N_VERT)
      real*8,dimension(:)    :: sigv(NZ-1)
      real*8,dimension(:)    :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)            
      real*8,dimension(:)    :: xc(N_CELL)
      real*8,dimension(:)    :: yc(N_CELL)
      real*8,dimension(:)    :: sig(NZ)
      real*8,dimension(:)    :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0)       
!     ----------------------------------------
      real*8,dimension(:,:)  :: phiNew(N_CELL,NZ)   
!     ----------------------------------------
      real*8  :: som,residu,errorsys
      integer :: iter,ii,jc1,jc2,jc3,jv1,jv2,jv3

      relaxSOR = 1.0d0
                  
!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|


      iter=0
111   continue
      iter=iter+1 

!      ________________________________________________________
!     |                                                        |
!     |       Solution of the system by SOR at each block      |
!     |________________________________________________________|

      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
!           _________________________________________________
!           Bm(vertex)
            residu = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                               +Bmv2T(i,k)*phiv(jv2,k)   &
                               +Bmv3T(i,k)*phiv(jv3,k)   &
                               +Bmv1B(i,k)*phiv(jv1,k-1) &
                               +Bmv2B(i,k)*phiv(jv2,k-1) &
                               +Bmv3B(i,k)*phiv(jv3,k-1) &
!           _________________________________________________
!           Am(center)
                               +Am1(i,k)*phi(jc1,k)      &
                               +Am2(i,k)*phi(jc2,k)      &
                               +Am3(i,k)*phi(jc3,k)      &
                               +AmT(i,k)*phi(i,k+1)      &       
                               +AmB(i,k)*phi(i,k-1)      &
                               +phi(i,k))
	        errorsys = errorsys + abs(residu)
	        phiNew(i,k) = phi(i,k) + relaxSOR*residu
         enddo
      enddo
      
      !errorsys = errorsys/(N_CELL0*(NZ-2))
    
      do k=2,NZ-1
         do i=1,N_CELL0
             phi(i,k) = phiNew(i,k)
         enddo
      enddo
!     ________________________________________________________
!     Boundary conditions cell-centers
!     ----------                  
!     Vertical
      do i=1,N_CELL0
         phi(i,1) = 2.0d0*funfB(i,1)-phi(i,2)
         phi(i,NZ)= 2.0d0*funfB(i,NZ)-phi(i,NZ-1)
         !if (i.le.10) print*,i,'Hello Serial',phi(i,NZ)
      enddo
!     ----------                 
!     Horizontal       
      do ii=N_CELL0+1,N_CELLexact
	     i = No_cp(ii,1)
         do k=1,NZ     
            phi(ii,k) = 2.0d0*funfB(ii,k)-phi(i,k)                  
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                Update vertex values                    |
!     |________________________________________________________|
 
!     ________________________________________________________
!     Interpolation of the inside vertex points 
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)                         
!     ________________________________________________________
!     Boundary Conditions of the vertex points 
      do nv=1,N_VERT
         phiv(nv,1)    = funfBv(nv,1)
         phiv(nv,NZ-1) = funfBv(nv,NZ-1)
         if (nbev(nv).ne.0) then
             do k=2,NZ-2 
                phiv(nv,k) = funfBv(nv,k)
             enddo
         endif
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |
!     |________________________________________________________|
        
         if (errorsys.lt.eps) then
            write(*,*) ' '
            write(*,7) 'Solution Jacobi 3D : iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
         elseif (errorsys.gt.1.0d5) then
            write(*,*) ' '
            write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
         elseif(iter.gt.MaxIters) then
            write(*,*) ' '
            write(*,7) 'Non-convergence: iters =',iter,&
                    ', error =',errorsys
            write(*,*) ' '
         else
            if (0.eq.mod(iter,100)) then
               print*, 'iterSerial =',iter,'Error=',errorsys
            endif
            goto 111
         endif

119   continue

      7 format(t10,a32,i5,a9,e10.3)

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   End of Jacobi S.O.R. Methods 3D                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
