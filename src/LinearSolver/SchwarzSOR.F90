!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!    SOLUTION OF THE 3D POISSON EQUATION BY Schwarz Decomposition.    !
!                              May 2015                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SchwarzSOR(phi,phiv,                    &
                            rhs,Gamx,Gamy,Gamz,          &
                            xc,yc,sig,dsig,No_cp,nbe,    &
                            xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,etan,                  &
                            Hprv,hv,etav)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the method based    !
!    on the Schwarz dwecomposition. See Walker2013 as a refenrence    !
!    of the technique.                                                !
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
      real*8 :: ErrorSch,som,residu,errorsys
      real*8 :: SUMErrorSch,SUMerrorsys
      real*8 :: x,y,z,zT,zB,fB,funSolExam3D
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
      integer:: ii,iterN,it,totalit,elem
      integer:: tag,SaveEpsResults,iter,Display
!     ----------------------------------------
      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8,dimension(:,:) :: funfBv(N_VERT,NZ-1)
!     ----------------------------------------
      real :: start,finish,startT,finishT
      real :: timeTotal,timeOther,timeBCcc,timeBCvv,timeComm,timeInter
!     ----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: Schwarz SOR'
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
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |        Boundary Condition of the initial guess         |
!     |________________________________________________________|

!     ______________________________________________________
!     Function fB at the boundary points
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
!     Function fB at the boundary vertex points
!     ___________
!     Vertical
      do nv=1,N_VERT
         x  = xv(nv)
         y  = yv(nv)
         zB = sigv(1)*Hprv(nv)-hv(nv)
         zT = sigv(NZ-1)*Hprv(nv)-hv(nv)
         funfBv(nv,1)    = funSolExam3D(x,y,zB) 
         funfBv(nv,NZ-1) = funSolExam3D(x,y,zT)
      enddo
!     ___________
!     Horizontal
      do k=1,NZ-1 
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
              x = xv(nv)
              y = yv(nv)
              z = sigv(k)*Hprv(nv)-hv(nv)
              funfBv(nv,k) = funSolExam3D(x,y,z) 
            endif
         enddo
      enddo
!     ______________________________________________________
!     Cell-centers BC   
      call cpu_time(start)
      do i=1,N_CELL0
         phi(i,1) = 2.0d0*funfB(i,1) -phi(i,2)
         phi(i,NZ)= 2.0d0*funfB(i,NZ)-phi(i,NZ-1)
      enddo  
      do ii=N_CELL0+1,N_CELLexact
	 i = No_cp(ii,1)
         do k=1,NZ     
            phi(ii,k) = 2.0d0*funfB(ii,k)-phi(i,k)                  
         enddo
      enddo
      call cpu_time(finish)
      timeBCcc = timeBCcc + (finish-start)



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
!                  Solution of the system by Schwarz                  !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|

      totalit = 0

      iterN=0
11    continue
      iterN=iterN+1

!      ________________________________________________________
!     |                                                        |
!     |Solution of the system by parallel S.0.R. methods (JSOR)|
!     |________________________________________________________|

      call AuxSchJacobi(it,timeBCcc,timeBCvv,timeInter,       &
                        phi,phiv,                             &
                        funFB,funfBv,rhs,                     &
                        Am1,Am2,Am3,AmT,AmB,                  &  
                        Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B,  &
                        xc,yc,sig,dsig,No_cp,nbe,             &
                        xv,yv,sigv,dsigv,No_vp,nbev)

!     ====================================
!     =========  COMMUNICATION  ==========
#     ifdef KeyParallel
        call communication3Dtype1(phi)
        !call communication3D(phi)
#     endif
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |            Error of each convergence iteration         |
!     |________________________________________________________|

      ErrorSch=0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            residu = rhs(i,k)-( Bmv1T(i,k)*phiv(No_vp(i,1),k)   &
                               +Bmv2T(i,k)*phiv(No_vp(i,2),k)   &
                               +Bmv3T(i,k)*phiv(No_vp(i,3),k)   &
                               +Bmv1B(i,k)*phiv(No_vp(i,1),k-1) &
                               +Bmv2B(i,k)*phiv(No_vp(i,2),k-1) &
                               +Bmv3B(i,k)*phiv(No_vp(i,3),k-1) &
                               +  Am1(i,k)*phi(No_cp(i,1),k)    &
                               +  Am2(i,k)*phi(No_cp(i,2),k)    &
                               +  Am3(i,k)*phi(No_cp(i,3),k)    &
                               +  AmT(i,k)*phi(i,k+1)           &       
                               +  AmB(i,k)*phi(i,k-1)           &
                               +           phi(i,k))
	    ErrorSch = ErrorSch + abs(residu)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call SUM_parallel(ErrorSch,SUMErrorSch)
         ErrorSch = SUMErrorSch
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
            write(7100,*),iterN,totalit,ErrorSch,timeTotal
         endif
#     endif

!      ________________________________________________________
!     |                                                        |
!     |          Convergence criteria of the method            |
!     |________________________________________________________|

      if (ErrorSch.lt.eps) then
          if (Display.eq.1) then
          write(*,*) ' '
          write(*,9) 'Sum of Schwarz system iter =',totalit
          write(*,8) 'Sol Schwarz Method: iters =',iterN,&
                     ', error =',ErrorSch
          write(*,*) ' '
          endif
      elseif (ErrorSch.gt.1.0d5) then
          if (Display.eq.1) then
          write(*,*) ' '
          write(*,8) 'DIVERGENCE !!!!: iters =',iterN,&
                     ', error =',ErrorSch
          write(*,*) ' '
          endif
          stop
      elseif(iterN.gt.MaxConv) then
          if (Display.eq.1) then
          write(*,*) ' '
          write(*,9) 'Sum of Schwarz system iter =',totalit
          write(*,8) 'Non-convergence: iters =',iterN,&
                     ', error =',ErrorSch
          write(*,*) ' '
          endif
      else
          if (Display.eq.1) then
          write(*,7) 'Schwarz Method: ',iterN,', error =',ErrorSch, &
                     ',  [System: Iters=',it,']'
          endif
     	  goto 11
      endif

19    continue

7     format(t7,a16,i5,a9,e10.3,a18,i5,a1)
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
      print*,'        Distribution of the time in SchwarzSOR'
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
      print*,'        Distribution of the time in SchwarzSOR'
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
         write(*,'(t15,60a)'), '<---- End   subroutine: Schwarz SOR'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                        AUXILIAR 1:  Jacobi                          !
!                              May 2015                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE AuxSchJacobi(it,timeBCcc,timeBCvv,timeInter,       &
                              psi,psiv,                             &
                              funFB,funfBv,rhs,                     &
                              Am1,Am2,Am3,AmT,AmB,                  &  
                              Bmv1T,Bmv2T,Bmv3T,Bmv1B,Bmv2B,Bmv3B,  &
                              xc,yc,sig,dsig,No_cp,nbe,             &
                              xv,yv,sigv,dsigv,No_vp,nbev)

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

      integer :: it
      real    :: timeBCcc,timeBCvv,timeInter
      real*8, dimension(:,:) :: psi(N_CELL,NZ)
      real*8, dimension(:,:) :: psiv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: funfB(N_CELL,NZ)
      real*8, dimension(:,:) :: funfBv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: rhs(N_CELL,NZ)
!     -------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      real*8, dimension(:)   :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!     -------------------------------------
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT) 
!     -------------------------------------
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8, dimension(:,:) :: Bmv3B(N_CELL0,NZ)
!     ----------------------------------------
      real*8, dimension(:,:) :: Newpsi(N_CELL,NZ)
      real*8  :: som,residu,errorsys,SUMerrorsys
      real    :: start,finish
      integer :: ii,jc1,jc2,jc3,Display

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

!     ______________________
!     Interpolation vertex 
      call cpu_time(start)
      call interpolation3D(psiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           psi,xc,yc,sig,dsig,No_cp,nbe)
      call cpu_time(finish)
      timeInter = timeInter + (finish-start)
!     ______________________
!     Boundary conditions v
      call cpu_time(start)
      do nv=1,N_VERT
         psiv(nv,1)    = funfBv(nv,1)
         psiv(nv,NZ-1) = funfBv(nv,NZ-1)
         if (nbev(nv).ne.0) then
            do k=2,NZ-2 
               psiv(nv,k) = funfBv(nv,k)
            enddo
         endif
      enddo
      call cpu_time(finish)
      timeBCvv = timeBCvv + (finish-start)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      it=0
111   continue
      it=it+1 
!     ______________________
!     Update internal points
      errorsys = 0.0d0
      do k=2,NZ-1
         do i=1,N_CELL0
            residu = rhs(i,k)-( Bmv1T(i,k)*psiv(No_vp(i,1),k)   &
                               +Bmv2T(i,k)*psiv(No_vp(i,2),k)   &
                               +Bmv3T(i,k)*psiv(No_vp(i,3),k)   &
                               +Bmv1B(i,k)*psiv(No_vp(i,1),k-1) &
                               +Bmv2B(i,k)*psiv(No_vp(i,2),k-1) &
                               +Bmv3B(i,k)*psiv(No_vp(i,3),k-1) &
                               +  Am1(i,k)*psi(No_cp(i,1),k)    &
                               +  Am2(i,k)*psi(No_cp(i,2),k)    &
                               +  Am3(i,k)*psi(No_cp(i,3),k)    &
                               +  AmT(i,k)*psi(i,k+1)           &       
                               +  AmB(i,k)*psi(i,k-1)           &
                               +           psi(i,k))
	    errorsys = errorsys + abs(residu)
	    psi(i,k) = psi(i,k) + relaxSOR*residu
         enddo
      enddo
      do k=2,NZ-1
         do i=1,N_CELL0
          !psi(i,k) = Newpsi(i,k)
         enddo
      enddo
!     ______________________
!     Boundary conditions c
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
!     ______________________
!     Interpolation vertex 
      call cpu_time(start)
      call interpolation3D(psiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           psi,xc,yc,sig,dsig,No_cp,nbe)
      call cpu_time(finish)
      timeInter = timeInter + (finish-start)
!     ______________________
!     Boundary conditions v
      call cpu_time(start)
      do nv=1,N_VERT
         psiv(nv,1)    = funfBv(nv,1)
         psiv(nv,NZ-1) = funfBv(nv,NZ-1)
         if (nbev(nv).ne.0) then
            do k=2,NZ-2 
               psiv(nv,k) = funfBv(nv,k)
            enddo
         endif
      enddo
      call cpu_time(finish)
      timeBCvv = timeBCvv + (finish-start)

!     ______________________
!     Convergence criteria
      if (errorsys.lt.epsConv) then 
         if (Display.eq.1) then
             !write(*,7) 'Solution Aux: iters =',it,', error =',errorsys
         endif
      elseif (errorsys.gt.1.0d5) then
         if (Display.eq.1) then
         write(*,7) ' DIVERGENCE !!!!: iters =',it,', error =',errorsys
         endif
         stop
      elseif(it.gt.MaxIters) then
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
!---------------------------------------------------------------------!
!                         END OF Schwarz S.O.R.                       !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



