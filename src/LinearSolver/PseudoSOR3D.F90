!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           SOLUTION OF THE 3D POISSON EQUATION BY NEW S.O.R.         !
!                             April 2013                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE PseudoSOR3D(phi,phiv,                    &
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

!     ====================================
!     =====  START PARALLEL OPTION =======
#     include "cppdefs.h"
#     ifdef KeyParallel
         USE parallel
#     endif
      USE geometry
      implicit none
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
      integer :: tag,SaveEpsResults,iter,Display
!     ----------------------------------------
      integer,dimension(:) :: Index_localBdy(N_CELL0)
      integer :: NNbdy,i0,jj
!     ----------------------------------------
#     ifdef KeyAuxMCSOR
      integer,dimension(:,:):: IniColor(0:3)
      integer,dimension(:,:):: FinColor(0:3)
      integer,dimension(:,:):: IndexColor(N_CELL0)
#     endif
!     ----------------------------------------
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      real :: timeInterA,timeBCccA,timeBCvvA,timeOtherA
      real :: timeLoopsA,timeCommA,timeTotalA
      real,dimension(Nprocs) :: timeInterV,timeBCccV,timeBCvvV,timeOtherV
      real,dimension(Nprocs) :: timeLoopsV,timeCommV,timeTotalV
#     endif
!     =============== END ================    
!     ==================================== 



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

#     ifdef KeyParallel
         i0 = 0
         s = rang_topo +1
         do ii = DimBdyCCdom(s,2),DimBdyCCdom(s,3)
            i = BdyCCdom(ii,2)
            do jj=1,N_CELL0
               if (Index_global(jj).eq.i) then
                 i0 = i0 +1
                 Index_localBdy(i0) = jj
               endif
            enddo
         enddo
         NNbdy = i0
#     endif	

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
 
#     ifdef KeyAuxJacobi
         call AuxJacobi(it,timeBCcc,timeComm,phi, &
                       funfB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp) 
#     endif
#     ifdef KeyAuxJSOR
         call AuxJSOR(it,timeBCcc,timeComm,phi, &
                       funfB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp) 
#     endif
#     ifdef KeyAuxMCSOR
         call AuxMCSOR(it,timeBCcc,timeComm,phi,&
                       funFB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp,&
                       IniColor,FinColor,IndexColor) 
#     endif
#     ifdef KeyAuxPSOR
         call AuxPSOR(it,timeBCcc,timeComm,phi,&
                      funFB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp,&
                      Index_localBdy,NNBdy)
#     endif
  
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
      print*,'          Distribution of the time in PD-SOR'
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
      print*,'          Distribution of the time in PD-SOR'
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
      write(*,*),'                 Distribution of the time in PD-SOR'
      write(*,*),'--------------------------------------------------------------------------- '
      write(*,*),'  p |  Interpo. |   BC_Cell | BC_Vertex |    Other  |   Comm.   |   Total   '
      write(*,*),'----|-----------|-----------|-----------|-----------|-----------|-----------'      
      do s=1,Nprocs
         write(*,201),s-1,' |',timeInterV(s),' |',timeBCccV(s),' |',timeBCvvV(s),' |',&
                    timeOtherV(s),' |',timeCommV(s),' |',timeTotalV(s)
      enddo
      201 format(i4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4)      
      write(*,*),'----|-----------|-----------|-----------|-----------|-----------|-----------'
      write(*,202),' Avg |',timeInterA,' |',timeBCccA,' |',timeBCvvA,' |',&
                    timeOtherA,' |',timeCommA,' |',timeTotalA
      202 format(a6,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4,a2,f10.4)
      write(*,*), ' '
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
!                        AUXILIAR 1:  Jacobi                          !
!                              May 2015                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AuxJacobi(it,timeBCcc,timeComm,psi,&
                           funFB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp)                

#     include "cppdefs.h"
!     ====================================
!     =====  START PARALLEL OPTION =======
#     include "cppdefs.h"
#     ifdef KeyParallel
         USE parallel
#     endif
      USE geometry
      implicit none
!     =============== END ================    
!     ==================================== 

      integer :: totalit
      real    :: timeBCcc,timeComm
      real*8, dimension(:,:) :: psi(N_CELL,NZ)
      real*8, dimension(:,:) :: funfB(N_CELL,NZ)
      real*8, dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: Vbm(N_CELL0,NZ) 
      integer,dimension(:,:) :: No_cp(N_CELL,3)
!     ----------------------------------------
      real*8, dimension(:,:) :: Newpsi(N_CELL,NZ)
      real*8  :: som,residu,errorsys,SUMerrorsys
      real    :: start,finish
      integer :: it,ii,jc1,jc2,jc3,Display

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

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      it=0
111   continue
      it=it+1 
!     ______________________
!     Update internal points
      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
	    residu = Vbm(i,k) -MAm1(i,k)*psi(No_cp(i,1),k) &
                              -MAm2(i,k)*psi(No_cp(i,2),k) &
                              -MAm3(i,k)*psi(No_cp(i,3),k) &
                              -MAmT(i,k)*psi(i,k+1) &       
                              -MAmB(i,k)*psi(i,k-1) &
	                      -psi(i,k)
	    errorsys = errorsys + abs(residu)
	    Newpsi(i,k) = psi(i,k) + residu
         enddo
      enddo
      do k=2,NZ-1
         do i=1,N_CELL0
	    psi(i,k) = Newpsi(i,k)
         enddo
      enddo
!     ______________________
!     Boundary conditions
      call cpu_time(start)
      do i=1,N_CELL0
         psi(i,1) = 2.0d0*funfB(i,1) -psi(i,2)
         psi(i,NZ)= 2.0d0*funfB(i,NZ)-psi(i,NZ-1)
      enddo  
      do ii=N_CELL0+1,N_CELLexact
	 i = No_cp(ii,1)
         do k=1,NZ     
            psi(ii,k) = 2.0d0*funfB(ii,k)-psi(i,k)                  
         enddo
      enddo
      call cpu_time(finish)
      timeBCcc = timeBCcc + (finish-start)

!     ====================================
!     =====  COMMUNICATION & ERROR =======
#     ifdef KeyParallel
         call cpu_time(start)
         !call communication3D(psi)
         !call communication3Dtype1(psi)
         call communication3Dtype2(psi)	 
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
         call cpu_time(finish)
         timeComm = timeComm + (finish-start)
#     endif	
!     =============== END ================    
!     ====================================

!     ______________________
!     Convergence criteria
      if (errorsys.lt.epsConv) then
         if (Display.eq.1) then
            write(*,7) 'Aux. Jacobi: iters =',it,', error =',errorsys
         endif
      elseif (errorsys.gt.1.0d5) then
         if (Display.eq.1) then
            write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         endif
         stop
      elseif(it.gt.MaxConv) then
         if (Display.eq.1) then
            write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
         endif
      else
         if (Display.eq.1) then
            if (0.eq.mod(it,1000)) then
               print*, 'iter=',it,'Error=',errorsys
            endif
         endif
         goto 111
      endif

119   continue

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

7     format(t10,a24,i5,a9,e10.3)
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                         AUXILIAR 2:  JSOR                           !
!                              May 2015                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AuxJSOR(it,timeBCcc,timeComm,psi,&
                         funFB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp)                

!     ====================================
!     =====  START PARALLEL OPTION =======
#     include "cppdefs.h"
#     ifdef KeyParallel
         USE parallel
#     endif
      USE geometry
      implicit none
!     =============== END ================    
!     ====================================  
      integer :: totalit
      real    :: timeBCcc,timeComm
      real*8, dimension(:,:) :: psi(N_CELL,NZ)
      real*8, dimension(:,:) :: funfB(N_CELL,NZ)
      real*8, dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: Vbm(N_CELL0,NZ) 
      integer,dimension(:,:) :: No_cp(N_CELL,3)
!     ----------------------------------------
      real*8  :: som,residu,errorsys,SUMerrorsys
      real    :: start,finish
      integer :: it,ii,jc1,jc2,jc3,Display

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

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      it=0
111   continue
      it=it+1 
!     ______________________
!     Update internal points

      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
	    residu = Vbm(i,k) -MAm1(i,k)*psi(jc1,k) &
                              -MAm2(i,k)*psi(jc2,k) &
                              -MAm3(i,k)*psi(jc3,k) &
                              -MAmT(i,k)*psi(i,k+1) &       
                              -MAmB(i,k)*psi(i,k-1) &
	                      -psi(i,k)
	    errorsys = errorsys + abs(residu)
	    psi(i,k) = psi(i,k) + relaxSOR*residu
         enddo
      enddo
!     ______________________
!     Boundary conditions

      call cpu_time(start)
      do i=1,N_CELL0
         psi(i,1) = 2.0d0*funfB(i,1) -psi(i,2)
         psi(i,NZ)= 2.0d0*funfB(i,NZ)-psi(i,NZ-1)
      enddo  
      do ii=N_CELL0+1,N_CELLexact
	 i = No_cp(ii,1)
         do k=1,NZ     
            psi(ii,k) = 2.0d0*funfB(ii,k)-psi(i,k)                  
         enddo
      enddo
      call cpu_time(finish)
      timeBCcc = timeBCcc + (finish-start)

!     ====================================
!     =====  COMMUNICATION & ERROR =======
#     ifdef KeyParallel
         call cpu_time(start)
         call communication3Dtype2(psi)
         !call communication3Dtype1(psi)
         !call communication3D(psi)
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
         call cpu_time(finish)
         timeComm = timeComm + (finish-start)
#     endif	
!     =============== END ================    
!     ====================================

!     ______________________
!     Convergence criteria
      if (errorsys.lt.epsConv) then
         if (Display.eq.1) then
            write(*,7) 'Aux. JSOR: iters =',it,', error =',errorsys
         endif
      elseif (errorsys.gt.1.0d5) then
         if (Display.eq.1) then
            write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         endif
         stop
      elseif(it.gt.MaxConv) then
         if (Display.eq.1) then
            write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
         endif
      else
         goto 111
      endif 

119   continue

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

7     format(t10,a24,i5,a9,e10.3)
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                         AUXILIAR 3:  MCSOR                          !
!                              May 2015                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AuxMCSOR(it,timeBCcc,timeComm,psi,&
                          funFB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp,&
                          IniColor,FinColor,IndexColor)                

!     ====================================
!     =====  START PARALLEL OPTION =======
#     include "cppdefs.h"
#     ifdef KeyParallel
         USE parallel
#     endif
      USE geometry
      implicit none
!     =============== END ================    
!     ====================================  
      integer :: totalit
      real    :: timeBCcc,timeComm
      real*8, dimension(:,:) :: psi(N_CELL,NZ)
      real*8, dimension(:,:) :: funfB(N_CELL,NZ)
      real*8, dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: Vbm(N_CELL0,NZ) 
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:,:):: IniColor(0:3)
      integer,dimension(:,:):: FinColor(0:3)
      integer,dimension(:,:):: IndexColor(N_CELL0)
!     ----------------------------------------
      real*8  :: som,residu,errorsys,SUMerrorsys
      real    :: start,finish
      integer :: it,ii,jc1,jc2,jc3,Display,m,s

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

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      it=0
111   continue
      it=it+1 
!     ______________________
!     Update internal points

      errorsys = 0.0d0
      DO s=0,N_COLOR-1
!        ____________________________________
!        Update color
         do k=2,NZ-1
            do m=IniColor(s),FinColor(s)
               i = IndexColor(m)
	       residu = Vbm(i,k) -MAm1(i,k)*psi(No_cp(i,1),k) &
                                 -MAm2(i,k)*psi(No_cp(i,2),k) &
                                 -MAm3(i,k)*psi(No_cp(i,3),k) &
                                 -MAmT(i,k)*psi(i,k+1) &       
                                 -MAmB(i,k)*psi(i,k-1) &
	                         -psi(i,k)
	       errorsys = errorsys + abs(residu)
	       psi(i,k) = psi(i,k) + relaxSOR*residu
            enddo
         enddo
!        ____________________________________
!        Parallel communication
!        ====================================
!        =========  COMMUNICATION ===========
#        ifdef KeyParallel
            call cpu_time(start)
            !call communication3Dtype1(psi)
            call communication3Dtype2(psi)
            !call communication3D(psi)
            call cpu_time(finish)
            timeComm = timeComm + (finish-start)
#        endif	
!        =============== END ================    
!        ====================================
      ENDDO
!     ______________________
!     Boundary conditions
      call cpu_time(start)
      do i=1,N_CELL0
         psi(i,1) = 2.0d0*funfB(i,1) -psi(i,2)
         psi(i,NZ)= 2.0d0*funfB(i,NZ)-psi(i,NZ-1)
      enddo  
      do ii=N_CELL0+1,N_CELLexact
	 i = No_cp(ii,1)
         do k=1,NZ     
            psi(ii,k) = 2.0d0*funfB(ii,k)-psi(i,k)                  
         enddo
      enddo
      call cpu_time(finish)
      timeBCcc = timeBCcc + (finish-start)

!     ====================================
!     =====  COMMUNICATION & ERROR =======
#     ifdef KeyParallel
         call cpu_time(start)
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
         call cpu_time(finish)
         timeComm = timeComm + (finish-start)
#     endif	
!     =============== END ================    
!     ====================================

!     ______________________
!     Convergence criteria
      if (errorsys.lt.epsConv) then
         if (Display.eq.1) then
            write(*,7) 'Aux. MCSOR: iters =',it,', error =',errorsys
         endif
      elseif (errorsys.gt.1.0d5) then
         if (Display.eq.1) then
            write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         endif
         stop
      elseif(it.gt.MaxConv) then
         if (Display.eq.1) then
            write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
         endif
      else
         if (Display.eq.1) then
            if (0.eq.mod(it,1000)) then
               print*, 'iter=',it,'Error=',errorsys
            endif
         endif
         goto 111
      endif 

119   continue

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

7     format(t10,a24,i5,a9,e10.3)
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                         AUXILIAR 4:  PSOR                           !
!                              May 2015                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AuxPSOR(it,timeBCcc,timeComm,psi,&
                         funFB,MAm1,MAm2,MAm3,MAmT,MAmB,Vbm,No_cp,&
                         Index_localBdy,NNBdy)                

!     ====================================
!     =====  START PARALLEL OPTION =======
#     include "cppdefs.h"
#     ifdef KeyParallel
         USE parallel
#     endif
      USE geometry
      implicit none
!     =============== END ================    
!     ====================================  
      integer :: totalit
      real    :: timeBCcc,timeComm
      real*8, dimension(:,:) :: psi(N_CELL,NZ)
      real*8, dimension(:,:) :: funfB(N_CELL,NZ)
      real*8, dimension(:,:) :: MAm1(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm2(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAm3(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: MAmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: Vbm(N_CELL0,NZ) 
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: Index_localBdy(N_CELL0)
      integer :: NNbdy
!     ----------------------------------------
      real*8  :: som,residu,errorsys,SUMerrorsys
      real    :: start,finish
      integer :: it,ii,jj,jc1,jc2,jc3,Display,m,s

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

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      it=0
111   continue
      it=it+1 
!     ______________________
!     Update Boundary points
!     ====================================
!     ====================================
#     ifdef KeyParallel
         DO k=2,NZ-1
            do jj=1,NeighborNumber
!              ---------------------------
!              Update
               do ii=1,BlockBdyNumber(jj)
                  i = BlockBdyIndex(ii,jj)
               !do ii = 1,NNbdy
               !   i = Index_localBdy(ii)
                  residu = Vbm(i,k)-( MAm1(i,k)*psi(No_cp(i,1),k) &
                                     +MAm2(i,k)*psi(No_cp(i,2),k) &
                                     +MAm3(i,k)*psi(No_cp(i,3),k) &
                                     +MAmT(i,k)*psi(i,k+1)        &       
                                     +MAmB(i,k)*psi(i,k-1)        &
                                     +psi(i,k))
	              psi(i,k) = psi(i,k) + relaxSOR*residu
               enddo
            enddo
         ENDDO

!        --------------------
!        Communicate
         call cpu_time(start)
         !call communication3D(psi)
         !call communication3Dtype1(psi)
         call communication3Dtype2(psi)
         !call communication3Dcc(psi)
         call cpu_time(finish)
         timeComm = timeComm + (finish-start)

#     endif	
!     =============== END ================    
!     ====================================
!     ______________________
!     Update internal points
      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
	    residu = Vbm(i,k) -MAm1(i,k)*psi(jc1,k) &
                              -MAm2(i,k)*psi(jc2,k) &
                              -MAm3(i,k)*psi(jc3,k) &
                              -MAmT(i,k)*psi(i,k+1) &       
                              -MAmB(i,k)*psi(i,k-1) &
	                      -psi(i,k)
	    errorsys = errorsys + abs(residu)
	    psi(i,k) = psi(i,k) + relaxSOR*residu
         enddo
      enddo
!     ______________________
!     Boundary conditions
      call cpu_time(start)
      do i=1,N_CELL0
         psi(i,1) = 2.0d0*funfB(i,1) -psi(i,2)
         psi(i,NZ)= 2.0d0*funfB(i,NZ)-psi(i,NZ-1)
      enddo  
      do ii=N_CELL0+1,N_CELLexact
	 i = No_cp(ii,1)
         do k=1,NZ     
            psi(ii,k) = 2.0d0*funfB(ii,k)-psi(i,k)                  
         enddo
      enddo
      call cpu_time(finish)
      timeBCcc = timeBCcc + (finish-start)

!     ====================================
!     =====  COMMUNICATION & ERROR =======
#     ifdef KeyParallel
         call cpu_time(start)
         !call communication3Dtype1(psi)
         call communication3Dtype2(psi)
         !call communication3D(psi)
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
         call cpu_time(finish)
         timeComm = timeComm + (finish-start)

#     endif	
!     =============== END ================    
!     ====================================

!     ______________________
!     Convergence criteria
      if (errorsys.lt.epsConv) then
         if (Display.eq.1) then
            write(*,7) 'Aux. PSOR: iters =',it,', error =',errorsys
         endif
      elseif (errorsys.gt.1.0d5) then
         if (Display.eq.1) then
            write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         endif
         stop
      elseif(it.gt.MaxConv) then
         if (Display.eq.1) then
            write(*,7) ' Non-convergence: iters =',it,', error =',errorsys
         endif
      else
         if (Display.eq.1) then
            if (0.eq.mod(it,1000)) then
               print*, 'iter=',it,'Error=',errorsys
            endif
         endif
         goto 111
      endif 

119   continue

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

7     format(t10,a24,i5,a9,e10.3)
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	  END OF Pseudo S.O.R. 3D                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



