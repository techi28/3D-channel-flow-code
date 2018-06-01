!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            S.O.R. FOR THE SOLUTION OF THE 3D POISSON EQUATION       !
!                              Dic 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE sor3D(phi,phiv,                    &
                       rhs,Gamx,Gamy,Gamz,          &
                       xc,yc,sig,dsig,No_cp,nbe,    &
                       xv,yv,sigv,dsigv,No_vp,nbev, &
                       Hpr,h,etan,                  &
                       Hprv,hv,etav)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the S.O.R. techni-  !
!    que for different relaxion factors: relax. This system is good   !
!    for the Poisson problem because it consider the cell-center "Am" !
!    and the vertex "Bmv" matrix coefficients.                        !
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
      real*8,dimension(:,:) :: Nm(N_CELL0,NZ)
      real*8,dimension(:,:) :: Newrhs(N_CELL0,NZ)
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      real*8  :: errorsys,residu,som
      real*8  :: errorNeum,errorsysOld
      integer :: iter
!     ----------------------------------------
      integer :: ii,jj
      real*8  :: x,y,z,zT,zB
      real*8  :: fB,dfBdn,Neumanndfdn3D
      real*8  :: funSolExam3D,funDiffExam3D 
      real*8  :: funExamNSp,funExamNSrhsp,NeumanndpdnNS
      real*8  :: nnx,nny,nnz
      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8,dimension(:,:) :: funfBv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: fundfdnB(N_CELL,NZ)
      real*8,dimension(:,:) :: fundfdnBv(N_VERT,NZ-1)
!     ---------------------------------------- 
      real :: start,finish,startT,finishT
      real :: timeTotal,timeOther
      real :: timeBCcc,timeBCvv,timeInter

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: SOR 3D'
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
      
!      ________________________________________________________
!     |                                                        |
!     |         Special case only testing Poisson              |
!     |________________________________________________________|

#     ifdef KeyTestOnlyPoisson
         do k=1,NZ
            do i=1,N_CELL
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
            Nm(i,k)   = 0.0d0
         enddo
      enddo

!     ________________________________________________________      
!     Diffusion with Dirichlet BC
      if (ChooseBoundary.eq.1) then
         call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             & 
                          Bmv1T,Bmv2T,Bmv3T,                   &
                          Bmv1B,Bmv2B,Bmv3B,                   &
                          Gamx,Gamy,Gamz,                      &
                          xc,yc,sig,dsig,No_cp,nbe,            &
                          xv,yv,sigv,dsigv,No_vp,nbev)
!     ________________________________________________________
!     Diffusion with Neumann BC (not used 2017)
      elseif (ChooseBoundary.eq.2) then
         call diffusionNeumann3D(Am0,Am1,Am2,Am3,AmT,AmB,      & 
                                 Bmv1T,Bmv2T,Bmv3T,            &
                                 Bmv1B,Bmv2B,Bmv3B,Nm,         &
                                 Gamx,Gamy,Gamz,               &
                                 xc,yc,sig,dsig,No_cp,nbe,     &
                                 xv,yv,sigv,dsigv,No_vp,nbev)
!        ----------------------------------
         do k=2,NZ-1
            do i=1,N_CELL0
               rhs(i,k) = rhs(i,k) - Nm(i,k)
            enddo
         enddo
      endif

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

!     ________________________________________________________
!    |                                                        |
!    |           rhs with NO vertex iteration update          |
!    |________________________________________________________|

#     ifndef KeyIterationWithVertex
         do k=1,NZ
            do i=1,N_CELL0
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               rhs(i,k) = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                                    +Bmv2T(i,k)*phiv(jv2,k)   &
                                    +Bmv3T(i,k)*phiv(jv3,k)   &
                                    +Bmv1B(i,k)*phiv(jv1,k-1) &
                                    +Bmv2B(i,k)*phiv(jv2,k-1) &
                                    +Bmv3B(i,k)*phiv(jv3,k-1))
            enddo
         enddo
#     endif
            
!*********************************************************************!
!                                                                     !
!                 Solution of the system (S.0.R.) 3D                  !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|

      errorsys = 0.0d0

      iter=0
111   continue
      iter=iter+1 

      errorsysOld = errorsys

!      ________________________________________________________
!     |                                                        |
!     |                Solution of the system SOR              |
!     |________________________________________________________|

      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            residu = rhs(i,k)
!           _________________________________________________
!           Bm(vertex)
#           ifdef KeyIterationWithVertex
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            residu = residu - ( Bmv1T(i,k)*phiv(jv1,k)   &
                               +Bmv2T(i,k)*phiv(jv2,k)   &
                               +Bmv3T(i,k)*phiv(jv3,k)   &
                               +Bmv1B(i,k)*phiv(jv1,k-1) &
                               +Bmv2B(i,k)*phiv(jv2,k-1) &
                               +Bmv3B(i,k)*phiv(jv3,k-1))
#           endif
!           _________________________________________________
!           Am(center)
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
            residu = residu - ( Am1(i,k)*phi(jc1,k)      &
                               +Am2(i,k)*phi(jc2,k)      &
                               +Am3(i,k)*phi(jc3,k)      &
                               +AmT(i,k)*phi(i,k+1)      &       
                               +AmB(i,k)*phi(i,k-1)      &
                               +phi(i,k))                  
            errorsys = errorsys + abs(residu)
            phi(i,k) = phi(i,k) + relaxSOR*residu
         enddo
      enddo

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

!      ________________________________________________________
!     |                                                        |
!     |                        Norm Error                      |
!     |________________________________________________________|

!     ________________________________________________________
!     Original L-2 norm
    
      !errorsys = errorsys/(N_CELL0*(NZ-2))
      
!     ________________________________________________________
!     Alternative Error (Neumann stop)

      errorNeum = abs(errorsys-errorsysOld)

!      ________________________________________________________
!     |                                                        |
!     |     Save multiples simulations varying epsilon         |
!     |________________________________________________________|

#     ifdef KeySaveVaryEps
         call cpu_time(finishT)
         timeTotal = (finishT-startT)  
         write(7100,*),iter,errorsys,timeTotal
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |
!     |________________________________________________________|

      !if ((errorsys.lt.eps).or.(errorNeum.lt.1d-4*eps)) then
      if ((errorsys.lt.eps).or.(errorNeum.lt.1d-0*eps)) then
      !if (errorsys.lt.eps) then
         write(*,*) ' '
         write(*,7) 'Solution S0R 3D : iters =',iter,&
                    ', error =',errorsys,',',errorNeum
         write(*,*) ' '
      elseif (errorsys.gt.1.0d10) then
         write(*,*) ' '
         write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                    ', error =',errorsys,',',errorNeum
         write(*,*) ' '
         stop
      elseif(iter.gt.MaxIters) then
         write(*,*) ' '
         write(*,7) 'Non-convergence: iters =',iter,&
                    ', error =',errorsys,',',errorNeum
         write(*,*) ' '
      else
         if (0.eq.mod(iter,1000)) then
            print*, 'iter=',iter,'Error=',errorsys,errorNeum
         endif
         goto 111
      endif

119   continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      7 format(t10,a24,i5,a9,e10.3,a2,e10.3)
      8 format(e10.3,e10.3,i6,f10.3)
      
!      ________________________________________________________
!     |                                                        |
!     |                  Distribution of times                 |
!     |________________________________________________________|

      call cpu_time(finishT)
      timeTotal = (finishT-startT)  

#     ifdef KeyDisplay 
      timeOther = timeTotal-(timeInter+timeBCcc+timeBCvv)      
      print*,' '
      print*,'        ------------------------------------'
      print*,'           Distribution of the time in SOR'
      print*,'        ------------------------------------'
      print*,'        Time interpolation :  ',timeInter
      print*,'        Time BC cell-center:  ',timeBCcc
      print*,'        Time BC vertex     :  ',timeBCvv
      print*,'        Time other         :  ',timeOther
      print*,'        ____________________________________'
      print*,'        Time total         :  ',timeTotal
      print*,'        ------------------------------------'
      print*,' '
      !write(7100,8),time,eps,iter,timeTotal 
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: SOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   Auxiliar: INTERPOLATION LSM                       !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE interpolation3DLSM(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                    phi,xc,yc,sig,dsig,No_cp,nbe)
       
!---------------------------------------------------------------------!
!                                                                     !
!    This program interpolate the values phi from the center of the   !
!    elements to the vertex of the prism using the LSM technique.     !  
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <-- phiv    |(N_VERT,NZ)| Function phi at the vertices        |  !  
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size    | Description                        |  !  
!  |_____________|____________|____________________________________|  !  
!  | --> phi     |(N_CELL,NZ) |Function phi at the cell center     |  !
!  | --> sigv    |(NZ-1)      |sigma coordinate vetex points       |  !
!  | --> dsigv   |(NZ-1)      |Increment of the sigma vertex points|  !  
!  | --> No_vp   |(N_CELL0,3) |Numbering of cell vertices          |  !
!  | --> nbev    |(N_VERT)    |Type of tag about the kind of vertex|  !
!  |_____________|____________|____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  ! 
!  |--- N_CELL0   | Number of the cell centers inside the domain   |  !
!  |--- N_CELL    | Total number of cell centers                   |  !
!  |--- N_VERT    | Number of the computing vertices               |  !
!  |--- NZ        | Number of points in the vertical direction     |  !  
!  |______________|________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,k       | Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    *   Common variables modified (integers)                         !
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

      real*8, dimension(:,:):: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  

      real*8, dimension(:,:):: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: phivdz0
      real*8,dimension(:),  allocatable :: SumFun      
      real*8 :: dzT,dzB,con1,con2,dxCV,dyCV
      real*8 :: f0,f1,f2,f3,dfundx,dfundy,funpGF
      integer,parameter :: ChooseLSMDist = 1 !(=1 Adding distance weights)
      integer,parameter :: Option = 1

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: interpolation3DLSM'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(phivdz0(N_VERT,NZ),SumFun(N_VERT))

!      ________________________________________________________
!     |                                                        |
!     |                Interpolation inside domain             |  
!     |________________________________________________________|

!     -------------------------------------------------
!     Interpolation vertices (Surrounding option)
      IF (Option.eq.1) THEN
      phivdz0 = 0.0d0
      do k=1,NZ
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then !<--- Only boundary
               do j=1,DimsurroundingLSM(nv)
                  nc = surroundingLSM(nv,j)
                  phivdz0(nv,k)= phivdz0(nv,k) + weightLSM(nv,j)*phi(nc,k)
               enddo
            endif 
         enddo
      enddo
      ELSE     
!     -------------------------------------------------
!     Interpolation vertices (Direct calculation option)
      phivdz0 = 0.0d0
      DO k=1,NZ
          do i=1,N_CELL0
             f0  = phi(i,k)
             f1  = phi(No_cp(i,1),k)
             f2  = phi(No_cp(i,2),k)
             f3  = phi(No_cp(i,3),k)
             dfundx = aGx(i,0)*f0+aGx(i,1)*f1+aGx(i,2)*f2+aGx(i,3)*f3
             dfundy = aGy(i,0)*f0+aGy(i,1)*f1+aGy(i,2)*f2+aGy(i,3)*f3
             !------------  
             do j=1,3
                nv = No_vp(i,j)
                funpGF = f0 + dfundx*(xv(nv)-xc(i))+dfundy*(yv(nv)-yc(i))                  
                if (ChooseLSMDist.eq.1) then
                   phivdz0(nv,k)=phivdz0(nv,k)+funpGF*dlCV(i,j)/dlVsumLSM(nv)
                else
                   phivdz0(nv,k)=phivdz0(nv,k)+funpGF/dlVsumLSM(nv)
                endif
             enddo
          enddo                                      
      ENDDO
      ENDIF 
                  
!     -------------------------------------------------
!     Interpolation vertices at position sigv(k=2:NZ)
      DO k=1,NZ-1        
         dzT = abs(sigv(k)-sig(k+1))
         dzB = abs(sigv(k)-sig(k))
         con1 = dzB/(dzT+dzB)
         con2 = dzT/(dzT+dzB)
         do nv=1,N_VERT
            !if (nbev(nv).ne.0) then !<--- Only boundary
                phiv(nv,k)= con1*phivdz0(nv,k+1) + con2*phivdz0(nv,k)
            !endif
         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(phivdz0,SumFun)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: interpolation3DLSM'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      End of S.O.R. Methods 3D                       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
