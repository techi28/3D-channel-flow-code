!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           SOLUTION OF THE 3D POISSON EQUATION BY NEW S.O.R.         !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE PDsor3Dcuda(phi,phiv,                    &
                             rhs,Gamx,Gamy,Gamz,          &
                             xc,yc,sig,dsig,No_cp,nbe,    &
                             xv,yv,sigv,dsigv,No_vp,nbev, &
                             Hpr,h,etan,                  &
                             Hprv,hv,etav)                

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
      real*8,dimension(:,:) :: Newrhs(N_CELL0,NZ) 
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      real*8,dimension(:,:) :: phin(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm0(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Vbm(N_CELL0,NZ)  
      real*8,dimension(:,:) :: Newphi(N_CELL,NZ)
!     ----------------------------------------
      real*8  :: ErrorConv,som,residu,errorsys,ErrorPD,SUMErrorPD
      integer :: ii,iterN,it,totalit,elem,m,s
      real*8,dimension(:,:),  allocatable :: M1,M2,M3,MT,MB
      real*8 :: x,y,z,zT,zB,fB,funSolExam3D
      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8 :: SUMErrorConv,SUMerrorsys
!     ----------------------------------------
      real :: start,finish,startT,finishT
      real :: timeTotal,timeOther,timeBCcc,timeBCvv,timeComm,timeInter
!     ----------------------------------------
      integer :: tag,SaveEpsResults,iter,Display,DoThis
!     ----------------------------------------
#     ifdef KeyAuxMCSOR
      integer,dimension(:,:):: IniColor(0:3)
      integer,dimension(:,:):: FinColor(0:3)
      integer,dimension(:,:):: IndexColor(N_CELL0)
#     endif
!     ----------------------------------------
      real*8, dimension(:,:), allocatable :: h_erreur   
      real*8  :: erreur
      integer :: N_CELLNZ                
!     =========================================================  

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Pseudo SOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      timeBCcc  = 0.0d0
      timeBCvv  = 0.0d0
      timeInter = 0.0d0
      timeComm  = 0.0d0

      call cpu_time(startT)

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


!*********************************************************************!
!                                                                     !
!                           CUDA Initialization                       !
!                                                                     !
!*********************************************************************! 

      
!     =========================================================
!     CUDA: Allocate memory  

      N_CELLNZ  = N_CELL*NZ
      
      allocate(d_phi(N_CELLNZ))
      allocate(d_phiNew(N_CELLNZ))                    
      allocate(d_rhs(N_CELLNZ))
      allocate(d_funfB(N_CELLNZ))      
      allocate(d_Am(1:5,N_CELLNZ))
      allocate(d_No_cp(N_CELL,1:3))                   
      allocate(d_erreur(grid%x,grid%y))
      allocate(h_erreur(grid%x,grid%y))         

!     ==========================================================
!     Transfer geometry

#     ifdef KeyMethPDMSOR_cuda
      allocate(d_Color(N_CELL0))
      do i=1,N_CELL0             
         d_Color(i) = ColorCell(i)
      enddo   
#     endif
      d_No_cp = No_cp
      
!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |        Boundary Condition of the initial guess         |
!     |________________________________________________________|

!     ______________________________________________________
!     Function fB at the boundary points
      call cpu_time(start)
!     ___________
!     Vertical
      do i=1,N_CELL0
         x  = xc(i)
         y  = yc(i) 
         zB = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)
         zT = 0.5d0*(sig(NZ-1)+sig(NZ))*Hpr(i)-h(i)
         funfB(i,1) = funSolExam3D(x,y,zB)
         funfB(i,NZ)= funSolExam3D(x,y,zT)
      enddo
!     ___________
!     Horizontal       
      do ii=N_CELL0+1,N_CELLexact
	 i = No_cp(ii,1)
	 j = No_cp(ii,2)
         x = xe(i,j)
         y = ye(i,j)
         do k=1,NZ 
            z = sig(k)*Hpr(i)-h(i)
            funfB(ii,k) = funSolExam3D(x,y,z)          
         enddo
       enddo
!     ______________________________________________________
!     Cell-centers BC   
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
      call cpu_time(finish)
      timeBCcc = finish-start
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

      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             & 
                       Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B, &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)

!*********************************************************************!
!                                                                     !
!                 Solution of the system (S.0.R.) 3D                  !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        New Matrix                      |
!     |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0 	
            MAm0(i,k) = ttheta*Am0(i,k)-ttau  
            MAm1(i,k) = ttheta*Am1(i,k)/MAm0(i,k)
            MAm2(i,k) = ttheta*Am2(i,k)/MAm0(i,k)
            MAm3(i,k) = ttheta*Am3(i,k)/MAm0(i,k)
            MAmT(i,k) = ttheta*AmT(i,k)/MAm0(i,k)
            MAmB(i,k) = ttheta*AmB(i,k)/MAm0(i,k)
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |               Reordering elements by color             |
!     |________________________________________________________|

#     ifdef KeyAuxMCSOR
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
#     endif
!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|

      totalit = 0

      iterN=0
11    continue
      iterN=iterN+1

!     ________________________________________________________
!     Save current iteration at phin 

      do k=1,NZ
         do i=1,N_CELL0	
            phin(i,k) = phi(i,k)
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                Update vertex values                    |
!     |________________________________________________________|
      
!     ------------------------------------------------
!     Interpolation  
      call cpu_time(start)
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)
      call cpu_time(finish)
      timeInter = timeInter + (finish-start)
!     ------------------------------------------------
!     Vertex BC 
      call cpu_time(start)
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav, &
                      phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
      call cpu_time(finish)
      timeBCvv = timeBCvv + (finish-start)
!      ________________________________________________________
!     |                                                        |
!     |                  New right-hand side                   |
!     |________________________________________________________|

      do k=2,NZ-1
         do i=1,N_CELL0 
!           ------------------------------
!           Cell-center terms
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)	
  	        som =   Am0(i,k)*phi(i,k)   &
                  + Am1(i,k)*phi(jc1,k) &
                  + Am2(i,k)*phi(jc2,k) &
                  + Am3(i,k)*phi(jc3,k) &
                  + AmT(i,k)*phi(i,k+1) &       
                  + AmB(i,k)*phi(i,k-1) 
!           ------------------------------
!           Vertex terms
            jv1 = No_vp(i,1)
            jv2 = No_vp(i,2)
            jv3 = No_vp(i,3)
            Newrhs(i,k) = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                                    +Bmv2T(i,k)*phiv(jv2,k)   &
                                    +Bmv3T(i,k)*phiv(jv3,k)   &
                                    +Bmv1B(i,k)*phiv(jv1,k-1) &
                                    +Bmv2B(i,k)*phiv(jv2,k-1) &
                                    +Bmv3B(i,k)*phiv(jv3,k-1) )
!           ------------------------------
!           Final rhs at each iteration
            Vbm(i,k) = (-(1-ttheta)*som -ttau*phi(i,k) + Newrhs(i,k))/MAm0(i,k)
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |Solution of the system by parallel S.0.R. methods (JSOR)|
!     |________________________________________________________|

      
!     =========================================================
!     Transfer data from host to device

      do i=1,N_CELL
         do k=1,NZ
            s = k+(i-1)*NZ
            d_phi(s)   = phi(i,k)
            d_phiNew(s)= phi(i,k)           
            d_rhs(s)   = Vbm(i,k)
            d_funfB(s) = funfB(i,k)             
            if (i.le.N_CELL0) then
               d_Am(1,s) = MAm1(i,k)
               d_Am(2,s) = MAm2(i,k)
               d_Am(3,s) = MAm3(i,k)
               d_Am(4,s) = MAmT(i,k)
               d_Am(5,s) = MAmB(i,k)
            else
               d_Am(1:5,s) = 0.0d0
            endif   
         enddo
      enddo    

      it=0
501   continue
      it=it+1
!     =========================================================
!     CUDA-KERNEL: Solution SOR

#     ifdef KeyMethPDMSOR_cuda
         !erreur = 0.0d0
         !do s=0,N_COLOR-1
         !   call kernel_cal_PDMCSOR<<<grid,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s)
         !   h_erreur = d_erreur
         !   erreur = erreur + sum(h_erreur)
         !enddo
         relaxSOR = 1.0d0         
         call kernel_cal_PDJacobi<<<grid,tBlock,shared_mem>>>(N_CELL0,NZ,relaxSOR,s)
         d_phi = d_phiNew         
         h_erreur = d_erreur
         erreur = erreur + sum(h_erreur)         
#     endif   
  
!     =========================================================
!     CUDA-KERNEL: Boundary conditions

      call kernel_cal_phiBC<<<grid,tBlock,shared_mem>>>(N_CELL0,N_CELL,NZ)

!     =========================================================
!     Criteria of convergence       
      if (erreur.lt.epsConv) then
         if (Display.eq.1) then
            write(*,7)   'Aux. PD-MCSOR: iters =',it,', error =',errorsys
         endif
      elseif (erreur.gt.1.0d5) then
         if (Display.eq.1) then
            write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         endif
         stop
      elseif(it.gt.MaxConv) then
         if (Display.eq.1) then
            write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
         endif
      else
         goto 501
      endif     

2001  continue
                                 
!     =========================================================
!     CUDA: Transfer data from device to host
      do i=1,N_CELL
         do k=1,NZ
            s = k+(i-1)*NZ
            phi(i,k) = d_phi(s)  
         enddo
      enddo   

    
    IF (DoThis.eq.1) THEN
             call AuxJacobi(it,timeBCcc,phi, &
                       funfB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp) 
    ENDIF
    
!      ________________________________________________________
!     |                                                        |
!     |            Error of each convergence iteration         |
!     |________________________________________________________|

      ErrorConv=0.0d0
      do i=1,N_CELL0 
         do k=2,NZ-1
            ErrorConv = ErrorConv +dabs(phi(i,k)-phin(i,k))**2
         enddo
      enddo
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call SUM_parallel(ErrorConv,SUMErrorConv)
         ErrorConv = dsqrt(SUMErrorConv)/float(N_CELL0global*(NZglobal-2))
#     else
         ErrorConv = dsqrt(ErrorConv)/float(N_CELL0*(NZ-2))  
#     endif	
!     =============== END ================    
!     ====================================
!      ________________________________________________________
!     |                                                        |
!     |     Save multiples simulations varying epsilon         |
!     |________________________________________________________|

      totalit = totalit + it

#     ifdef KeySaveVaryEps
         if (Display.eq.1) then
            call cpu_time(finishT)
            timeTotal = (finishT-startT)  
            write(7100,*),iterN,ErrorConv,timeTotal,it
         endif
#     endif
!      ________________________________________________________
!     |                                                        |
!     |          Convergence criteria of the method            |
!     |________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel         
      if (ErrorConv.lt.eps) then
          write(*,*) ' '
          write(*,9) 'Sum of SOR system iterations =',totalit
          write(*,8) 'Sol PDSOR-CUDA Method: iters =',iterN,&
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
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      if (ErrorConv.lt.eps) then
          IF (rang_topo.eq.0) THEN
          write(*,*) ' '
          write(*,9) 'Sum of SOR system iterations =',totalit
          write(*,8) 'Sol Psedo-time Method: iters =',iterN,&
                     ', error =',ErrorConv
          write(*,*) ' '
          ENDIF
      elseif (ErrorConv.gt.1.0d5) then
          IF (rang_topo.eq.0) THEN
          write(*,*) ' '
          write(*,8) 'DIVERGENCE !!!!: iters =',iterN,&
                     ', error =',ErrorConv
          write(*,*) ' '
          ENDIF
          stop
      elseif(iterN.gt.MaxIters) then
          IF (rang_topo.eq.0) THEN
          write(*,*) ' '
          write(*,8) 'Non-convergence: iters =',iterN,&
                     ', error =',ErrorConv
          write(*,*) ' '
          ENDIF
      else
     	  goto 11
      endif
#     endif
!     =============== END ================    
!     ==================================== 

19    continue

7     format(t10,a24,i5,a9,e10.3)
8     format(t10,a32,i5,a9,e10.3)
9     format(t10,a32,i5)

!*********************************************************************!
!                                                                     !
!                             Finalization                            !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |            Deallocate memory (CUDA variables)          |
!     |________________________________________________________| 

      deallocate(d_phi(N_CELLNZ))
      deallocate(d_phiNew(N_CELLNZ))                    
      deallocate(d_rhs(N_CELLNZ))
      deallocate(d_funfB(N_CELLNZ))      
      deallocate(d_Am(1:5,N_CELLNZ))    
      deallocate(d_No_cp(N_CELL,1:3))                          
      deallocate(d_erreur(grid%x,grid%y))      
      deallocate(h_erreur(grid%x,grid%y))
      !---------       
#     ifdef KeyMethMSOR_cuda             
      deallocate(d_Color(N_CELL0))
#     endif 

!      ________________________________________________________
!     |                                                        |
!     |                  Distribution of times                 |
!     |________________________________________________________|

      call cpu_time(finishT)
      timeTotal = (finishT-startT)   

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      timeOther = timeTotal-(timeInter+timeBCcc+timeBCvv)
      print*,' '
      print*,'        -------------------------------------'
      print*,'        Distribution of the time in PDSOR-CUDA'
      print*,'        -------------------------------------'
      print*,'        Time interpolation :  ',timeInter
      print*,'        Time BC cell-center:  ',timeBCcc
      print*,'        Time BC vertex     :  ',timeBCvv
      print*,'        Time other         :  ',timeOther
      print*,'        _____________________________________'
      print*,'        Time total         :  ',timeTotal
      print*,'        -------------------------------------'
      print*,' '
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      timeOther = timeTotal-(timeInter+timeBCcc+timeBCvv+timeComm)
      IF (rang_topo.eq.0) THEN
      print*,'        -------------------------------------'
      print*,'        Distribution of the time in PseudoSOR'
      print*,'        -------------------------------------'
      print*,'        Time interpolation :  ',timeInter
      print*,'        Time BC cell-center:  ',timeBCcc
      print*,'        Time BC vertex     :  ',timeBCvv
      print*,'        Time other         :  ',timeOther
      print*,'        Time communication :  ',timeComm
      print*,'        _____________________________________'
      print*,'        Time total         :  ',timeTotal
      print*,'        -------------------------------------'
      print*,' '
      ENDIF
#     endif
!     =============== END ================    
!     ==================================== 

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: Pseudo SOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	  END OF Pseudo S.O.R. 3D                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



