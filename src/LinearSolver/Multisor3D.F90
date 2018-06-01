!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            MULTI-COLORING SOR FOR THE 3D POISSON EQUATION           !
!                              Dic 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Multisor3D(phi,phiv,                    &
                            rhs,Gamx,Gamy,Gamz,          &
                            xc,yc,sig,dsig,No_cp,nbe,    &
                            xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,etan,                  &
                            Hprv,hv,etav)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the Multicoloring   !
!    S.O.R. technique. We need the file corresponding to the color    !
!    of each cell as well as the number of colors.                    !
!    We use different relaxion factors: relaxSOR. This system is good !
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

      integer,dimension(:,:):: IniColor(0:3)
      integer,dimension(:,:):: FinColor(0:3)
      integer,dimension(:,:):: IndexColor(N_CELL0)

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
      real*8,dimension(:,:) :: suma(N_CELL0,NZ)
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      real*8  :: errorsys,residu,som,SUMerrorsys
      real*8  :: errorNeum,SUMerrorNeum,errorsysOld
      integer :: s,m,iter
      integer :: SaveEpsResults
!     ----------------------------------------
      integer :: ii,jj,ivert
      integer :: Display
      real*8  :: nnx,nny,nnz,x,y,z,zT,zB
      real*8  :: fB,dfBdn,Neumanndfdn3D
      real*8  :: funSolExam3D,funDiffExam3D
      real*8  :: funExamNSp,funExamNSrhsp,NeumanndpdnNS
      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8,dimension(:,:) :: funfBv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: fundfdnB(N_CELL,NZ)
      real*8,dimension(:,:) :: fundfdnBv(N_VERT,NZ-1)
!     ----------------------------------------
      real :: start,finish,startT,finishT
      real :: timeTotal,timeOther
      real :: timeBCcc,timeBCvv,timeInter
      real :: timeComm
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: elem,Iini
!     ------------------------------------
      real :: timeInterA,timeBCccA,timeBCvvA,timeOtherA
      real :: timeLoopsA,timeCommA,timeTotalA
      real,dimension(Nprocs) :: timeInterV,timeBCccV,timeBCvvV,timeOtherV
      real,dimension(Nprocs) :: timeLoopsV,timeCommV,timeTotalV      
#     endif
!     =============== END ================    
!     ==================================== 
 
!      ____________________________________
!     |                                    |
!     |      Declaration of parameters     |
!     |____________________________________|

      integer, parameter :: BCoption = 2 !=1 Direchlet, 
                                         !=2 Calling BCcellcenter3D

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: MultiSOR 3D'
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
!     |                Initial times & display                 |
!     |________________________________________________________|

      timeTotal = 0.0d0
      timeOther = 0.0d0
      timeBCcc  = 0.0d0
      timeBCvv  = 0.0d0
      timeInter = 0.0d0
      timeComm  = 0.0d0
      
!     ====================================
!     =====    DISPLAY ITERATIONS  =======
      Display = 1
#     ifdef KeyParallel
         if (rang_topo.ne.0) then
            Display = 0
         endif
#     endif	
!     =============== END ================    
!     ====================================

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
!                 Solution of the system (M.S.0.R.) 3D                !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |               Reordering elements by color             |
!     |________________________________________________________|

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
!     |   MSOR: Solution of the system by SOR for each color   |
!     |________________________________________________________|

!     ________________________________________________________
!     Update: Bm(vertex) + Am(old)

      errorsys = 0.0d0
      DO s=0,N_COLOR-1
!        ____________________________________
!        Update color
         do k=2,NZ-1
            do m=IniColor(s),FinColor(s)
               i = IndexColor(m)
               residu = rhs(i,k)
!              ______________________________________________
!              Bm(vertex)
#             ifdef KeyIterationWithVertex
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)
               residu = residu - ( Bmv1T(i,k)*phiv(jv1,k)   &
                                  +Bmv2T(i,k)*phiv(jv2,k)   &
                                  +Bmv3T(i,k)*phiv(jv3,k)   &
                                  +Bmv1B(i,k)*phiv(jv1,k-1) &
                                  +Bmv2B(i,k)*phiv(jv2,k-1) &
                                  +Bmv3B(i,k)*phiv(jv3,k-1))
#              endif
!              ______________________________________________
!              Am(center)
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
         
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call cpu_time(start)
            !call communication3D(phi)
            !call communication3Dtype1(phi)
            call communication3Dtype2(phi) 
            call cpu_time(finish)
            timeComm = timeComm + (finish-start)
#        endif	
!        =============== END ================    
!        ====================================
      ENDDO

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
!     Parallel error

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
#     endif	
!     =============== END ================    
!     ====================================

!     ________________________________________________________
!     Original L-2 norm
    
      !errorsys = errorsys/(N_CELL0global*(NZ-2))

!     ________________________________________________________
!     Alternative Error (Neumann stop)

      errorNeum = abs(errorsys-errorsysOld)

!      ________________________________________________________
!     |                                                        |
!     |     Save multiples simulations varying epsilon         |
!     |________________________________________________________|

#     ifdef KeySaveVaryEps
         SaveEpsResults = 1        
!        ====================================
!        ====================================
#        ifdef KeyParallel
            if (rang_topo.ne.0) then
               SaveEpsResults = 0
            endif
#        endif
!        ====================================   
!        ==================================== 
         if (SaveEpsResults.eq.1) then
            call cpu_time(finishT)
            timeTotal = (finishT-startT)  
            write(7100,*),iter,errorsys,timeTotal
         endif
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |
!     |________________________________________________________|
        
      !if ((errorsys.lt.eps).or.(errorNeum.lt.1d-4*eps)) then
      if ((errorsys.lt.eps).or.(errorNeum.lt.1d-0*eps)) then
         IF (Display.eq.1) THEN
            write(*,*) ' '
            write(*,7) 'Solution M.S.0.R. 3D : iters =',iter,&
                       ', error =',errorsys,',',errorNeum
            write(*,*) ' '
         ENDIF
      elseif (errorsys.gt.1.0d10) then
         IF (Display.eq.1) THEN
            write(*,*) ' '
            write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                       ', error =',errorsys,',',errorNeum
            write(*,*) ' '
         ENDIF
         stop
      elseif(iter.gt.MaxIters) then
         IF (Display.eq.1) THEN
            write(*,*) ' '
            write(*,7) 'Non-convergence: iters =',iter,&
                       ', error =',errorsys,',',errorNeum
            write(*,*) ' '
         ENDIF
      else
         IF (Display.eq.1) THEN
            if (0.eq.mod(iter,1000)) then
               print*, 'iter=',iter,'Error=',errorsys,errorNeum
            endif
         ENDIF
         goto 111
      endif

119   continue

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      7 format(t10,a32,i5,a9,e10.3,a2,e10.3)
      8 format(e10.3,e10.3,i6,f10.3)

!      ________________________________________________________
!     |                                                        |
!     |                  Distribution of times                 |
!     |________________________________________________________|

      call cpu_time(finishT)
      timeTotal = (finishT-startT)   

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
#     ifdef KeyDisplay 
      timeOther = timeTotal-(timeInter+timeBCcc+timeBCvv)
      print*,' '
      print*,'        ------------------------------------'
      print*,'        Distribution of the time in MSOR'
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
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
#     ifdef KeyDisplay 
      timeOther = timeTotal-(timeInter+timeBCcc+timeBCvv+timeComm)
      IF (rang_topo.eq.0) THEN
      print*,'        ------------------------------------'
      print*,'        Distribution of the time in MSOR'
      print*,'        ------------------------------------'
      print*,'        Time interpolation :  ',timeInter
      print*,'        Time BC cell-center:  ',timeBCcc
      print*,'        Time BC vertex     :  ',timeBCvv
      print*,'        Time other         :  ',timeOther
      print*,'        Time communication :  ',timeComm
      print*,'        ____________________________________'
      print*,'        Time total         :  ',timeTotal
      print*,'        ------------------------------------'
      print*,' '
      !write(7100,8),time,eps,iter,timeTotal 
      ENDIF
      call MPI_Barrier(comm3D,code)      
      call MPI_ALLGATHER(timeInter,1,MPI_FLOAT,timeInterV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(timeBCcc, 1,MPI_FLOAT,timeBCccV, 1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(timeBCvv, 1,MPI_FLOAT,timeBCvvV, 1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(timeOther,1,MPI_FLOAT,timeOtherV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(timeComm, 1,MPI_FLOAT,timeCommV, 1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(timeTotal,1,MPI_FLOAT,timeTotalV,1,MPI_FLOAT,comm3D,code)
      timeInterA = sum(timeInterV)/Nprocs 
      timeBCccA  = sum(timeBCccV)/Nprocs
      timeBCvvA  = sum(timeBCvvV)/Nprocs
      timeOtherA = sum(timeOtherV)/Nprocs
      timeCommA  = sum(timeCommV)/Nprocs      
      timeTotalA = sum(timeTotalV)/Nprocs            
      IF (Display.eq.1) THEN 
      write(*,*),'                 Distribution of the time in MSOR'
      write(*,*),'--------------------------------------------------------------------------- '
      write(*,*)'  p |  Interpo. |   BC_Cell | BC_Vertex |    Other  |   Comm.   |   Total   '
      write(*,*)'----|-----------|-----------|-----------|-----------|-----------|-----------'      
      do s=1,Nprocs
         write(*,19),s-1,' |',timeInterV(s),' |',timeBCccV(s),' |',timeBCvvV(s),' |',&
                    timeOtherV(s),' |',timeCommV(s),' |',timeTotalV(s)
      enddo
      19 format(i4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4)      
      write(*,*)'----|-----------|-----------|-----------|-----------|-----------|-----------'
      write(*,18),' Avg |',timeInterA,' |',timeBCccA,' |',timeBCvvA,' |',&
                    timeOtherA,' |',timeCommA,' |',timeTotalA
      18 format(a6,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4)
      write(*,*), ' '
      ENDIF
#     endif     
#     endif
!     =============== END ================    
!     ==================================== 
	         
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: MultiSOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               End of Multicoloring S.O.R. Method 3D                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
