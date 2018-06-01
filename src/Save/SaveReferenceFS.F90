!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  SAVING DATA AT REFERENCE POINTS                    !
!                             March 2016                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

    SUBROUTINE FS_SaveReference(                       &
                      Hprv,etav,ufv,vfv,wfv,pfv,       & 
                      HprvA,etavA,ufvA,vfvA,wfvA,pfvA, &                      
                      h,xc,yc,sig,dsig,No_cp,nbe,      &
                      hv,xv,yv,sigv,dsigv,No_vp,nbev,  & 
                      No_sp) 
        
!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
#     endif
!     =============== END ================    
!     ====================================
      USE geometry
      implicit none
  
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|

      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)     
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)  
!     --------------------------------------      
      real*8, dimension(:)  :: HprvA(N_VERT)
      real*8, dimension(:)  :: etavA(N_VERT)
      real*8,dimension(:,:) :: ufvA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfvA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfvA(N_VERT,NZ-1)     
      real*8,dimension(:,:) :: pfvA(N_VERT,NZ-1)           
!     --------------------------------------
      real*8, dimension(:)  :: h(N_CELL)           
      real*8, dimension(:)  :: hv(N_VERT)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
      integer,dimension(:)  :: nbev(N_VERT)         
!     --------------------------------------
      real*8,dimension(:,:) :: xct(N_CELL,NZ)
      real*8,dimension(:,:) :: yct(N_CELL,NZ)
      real*8,dimension(:,:) :: zct(N_CELL,NZ)
      real*8,dimension(:,:) :: xvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: yvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: zvt(N_VERT,NZ-1)
!     --------------------------------------
      integer,dimension(:)  :: No_sp(N_SPmax)
      integer :: EachStep
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: FS_funu,FS_funv,FS_funw,FS_funp,FS_funeta,FS_funH
      real*8 :: MaxErrorA,sumErrorA,x,y,z
      real*8 :: som,Vol0,Vol,Are,Hadd
      real*8 :: DR,DR0,Min_DR0
      integer:: m,SaveFiles,elem,ii,vert
      integer:: nv1,nv2,Use_range1,Use_range2
      real*8,dimension(:) :: xvR1(N_SPMAX),yvR1(N_SPMAX)
      integer,dimension(:) :: OpenFiles(N_SPMAX)
      integer :: pro,n,s,NumS
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      real*8, dimension(Nprocs) :: DR0V
      integer,dimension(Nprocs) :: nv1V,proV,UseProc
#     endif
!     =============== END ================    
!     ====================================
  
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>> Begin subroutine: FS_SaveReference'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Initial time of the Example                     !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                       Exact solution                   |
!     |________________________________________________________|

      etavA = etav
      HprvA = Hprv
      ufvA  = ufv
      vfvA  = vfv
      wfvA  = wfv
      pfvA  = pfv
         
      do nv=1,N_VERT
!        ---------------------
#        ifdef KeyStandingWave
            etavA(nv) = FS_funeta(xv(nv),yv(nv),time)
            HprvA(nv) = etavA(nv) + hv(nv)
            do k=1,NZ-1
               z = sig(k)*HprvA(nv)-hv(nv)
               ufvA(nv,k) = FS_funu(xv(nv),yv(nv),z,time)
               vfvA(nv,k) = FS_funv(xv(nv),yv(nv),z,time)
               wfvA(nv,k) = FS_funw(xv(nv),yv(nv),z,time) 
               pfvA(nv,k) = FS_funp(xv(nv),yv(nv),z,time)  
            enddo
#        endif
!        ---------------------
#        ifdef KeySolitaryWave
            etavA(nv) = FS_SolitaryWave(xv(nv),time)
            HprvA(nv) = etavA(nv) + hv(nv)
            do k=1,NZ-1
               z = sig(k)*HprvA(nv)-hv(nv)
               ufvA(nv,k) = FS_funu(xv(nv),yv(nv),z,time)
               vfvA(nv,k) = FS_funv(xv(nv),yv(nv),z,time)
               wfvA(nv,k) = FS_funw(xv(nv),yv(nv),z,time) 
               pfvA(nv,k) = FS_funp(xv(nv),yv(nv),z,time)  
            enddo
#        endif

      enddo
!      ________________________________________________________
!     |                                                        |
!     |               Reference points nv calculation          |
!     |________________________________________________________|

      !IF (time.eq.0.0d0) THEN
!        ---------------------
#        ifdef KeyStandingWave
            NumS = 2
            xvR1(1) = 2.25d0
            yvR1(1) = 2.25d0
            xvR1(2) = 0.25d0
            yvR1(2) = 0.25d0
#        endif
!        ---------------------
#        ifdef KeyStaticCylinder
            NumS = 3
            xvR1(1) =  0.00d0
            yvR1(1) =  0.60d0
            xvR1(2) =  0.60d0
            yvR1(2) =  0.00d0
            xvR1(3) = -0.60d0
            yvR1(3) =  0.00d0
#        endif
!        ---------------------
#        ifdef KeyStaticChannel
            NumS = 2
            xvR1(1) =  0.00d0
            yvR1(1) =  0.60d0
            xvR1(2) =  0.60d0
            yvR1(2) =  0.00d0
#        endif
!        ---------------------
#        ifdef KeySolitaryWave
            NumS = 2
            xvR1(1) = 200.0d0
            yvR1(1) = 5.0d0
            xvR1(2) = 400.0d0
            yvR1(2) = 5.0d0
#        endif
!        ---------------------

!      ________________________________________________________
!     |                                                        |
!     |     Find the closest index to the reference points     |
!     |________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel 
         DO s=1,NumS  
            nv1 = 1
            DR0 = sqrt((xvR1(s)-xv(1))**2+(yvR1(s)-yv(1))**2)      
            do nv=2,N_VERT
              DR = sqrt((xvR1(s)-xv(nv))**2+(yvR1(s)-yv(nv))**2)
              if (DR.lt.DR0) then
                 nv1 = nv
                 DR0 = DR
              endif 
            enddo
            No_sp(s) = nv1 
         ENDDO    
!     ====================================
!     ====================================
#     else 
         DO s=1,NumS  
            nv1 = 1
            DR0 = sqrt((xvR1(s)-xv(1))**2+(yvR1(s)-yv(1))**2)      
            do nv=2,N_VERT
               DR = sqrt((xvR1(s)-xv(nv))**2+(yvR1(s)-yv(nv))**2)
               if (DR.lt.DR0) then
                  nv1 = nv
                  DR0 = DR
               endif 
            enddo
            pro = rang_topo
            call MPI_ALLGATHER(pro,1,MPI_INTEGER,proV,1,MPI_INTEGER,comm3D,code)
            call MPI_ALLGATHER(nv1,1,MPI_INTEGER,nv1V,1,MPI_INTEGER,comm3D,code)
            call MPI_ALLGATHER(DR0,1,MPI_DOUBLE_PRECISION,DR0V,1,MPI_DOUBLE_PRECISION,comm3D,code)
            pro = proV(1)
            nv1 = nv1V(1)
            DR0 = DR0V(1)
            do n=2,Nprocs
              if (DR0V(n).lt.DR0) then
                  nv1 = nv1V(n)
                  pro = proV(n)
                  DR0 = DR0V(n)
               endif 
            enddo
            UseProc(s) = pro
            No_sp(s)   = nv1
         ENDDO
#     endif         
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                       Open data files                  |
!     |________________________________________________________|

      OpenFiles(1) = 0
      OpenFiles(2) = 0
      OpenFiles(3) = 0
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         OpenFiles(1) = 1
         OpenFiles(2) = 1
         OpenFiles(3) = 1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         if (rang_topo.eq.UseProc(1))  OpenFiles(1) = 1
         if (rang_topo.eq.UseProc(2))  OpenFiles(2) = 1
         if (rang_topo.eq.UseProc(3))  OpenFiles(3) = 1
#     endif
!     =============== END ================    
!     ====================================

      IF (time.eq.0.0d0.and.OpenFiles(1).eq.1) THEN
         print*,' Saved: output/FS/FS_SolTime1.dat'
         open(9101,file="../output/FS/FS_SolTime1.dat",status='unknown')
         print*,' Saved: output/FS/FS_SolTimeE1.dat'
         open(9103,file="../output/FS/FS_SolTimeE1.dat",status='unknown')
      ENDIF
      IF (time.eq.0.0d0.and.OpenFiles(2).eq.1) THEN
         print*,' Saved: output/FS/FS_SolTime2.dat'
         open(9102,file="../output/FS/FS_SolTime2.dat",status='unknown')
         print*,' Saved: output/FS/FS_SolTimeE2.dat'
         open(9104,file="../output/FS/FS_SolTimeE2.dat",status='unknown')        
      ENDIF
!     --------------------------------------------------------
#     ifdef KeySave1DRef3Points
      IF (time.eq.0.0d0.and.OpenFiles(3).eq.1) THEN
         print*,' Saved: output/FS/FS_SolTime3.dat'
         open(9105,file="../output/FS/FS_SolTime3.dat",status='unknown')
      ENDIF
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                       Saving data                      |
!     |________________________________________________________|
             
      EachStep = floor(1.0d0/(5.0d0*dt)) 
      EachStep = max(EachStep,1)
      EachStep = 10
            
      if (0.eq.mod(nstep,EachStep)) then       
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         IF (rang_topo.eq.UseProc(1)) THEN
         nv = No_sp(1)
         write(9101,*) time,xv(nv),yv(nv),etav(nv),& 
                       ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1),pfv(nv,1:NZ-1)
         write(9103,*) time,xv(nv),yv(nv),etavA(nv),&
                       ufvA(nv,1:NZ-1),vfvA(nv,1:NZ-1),wfvA(nv,1:NZ-1),pfvA(nv,1:NZ-1)
         ENDIF
!        -------------         
         IF (rang_topo.eq.UseProc(2)) THEN
         nv = No_sp(2)
         write(9102,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1),pfv(nv,1:NZ-1)
         write(9104,*) time,xv(nv),yv(nv),etavA(nv),&
                       ufvA(nv,1:NZ-1),vfvA(nv,1:NZ-1),wfvA(nv,1:NZ-1),pfvA(nv,1:NZ-1)
         ENDIF
!        -------------
#        ifdef KeySave1DRef3Points
         IF (rang_topo.eq.UseProc(3)) THEN
         nv = No_sp(3)
         write(9105,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1),pfv(nv,1:NZ-1)
         ENDIF
#        endif
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
         nv = No_sp(1)
         write(9101,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1),pfv(nv,1:NZ-1)
         write(9103,*) time,xv(nv),yv(nv),etavA(nv),&
                       ufvA(nv,1:NZ-1),vfvA(nv,1:NZ-1),wfvA(nv,1:NZ-1),pfvA(nv,1:NZ-1)
!        -------------
         nv = No_sp(2)
         write(9102,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1),pfv(nv,1:NZ-1)
         write(9104,*) time,xv(nv),yv(nv),etavA(nv),&
                       ufvA(nv,1:NZ-1),vfvA(nv,1:NZ-1),wfvA(nv,1:NZ-1),pfvA(nv,1:NZ-1)
!        -------------
#        ifdef KeySave1DRef3Points
         nv = No_sp(3)
         write(9105,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1),vfv(nv,1:NZ-1),wfv(nv,1:NZ-1),pfv(nv,1:NZ-1)
#        endif
#     endif
!     =============== END ================    
!     ====================================
      endif
      !ENDIF
      
!      ________________________________________________________
!     |                                                        |
!     |                   Closing files data                   |
!     |________________________________________________________|
     
      IF (abs(tfin-time).le.1.0d-5.and.OpenFiles(1).eq.1) THEN
            close(9101)
            close(9103)
      ENDIF
      IF (abs(tfin-time).le.1.0d-5.and.OpenFiles(2).eq.1) THEN         
            close(9102)
            close(9104)    
      ENDIF
!     --------------------------------------------------------
#     ifdef KeySave1DRef3Points
      IF (abs(tfin-time).le.1.0d-5.and.OpenFiles(3).eq.1) THEN         
            close(9105)
      ENDIF
#     endif      
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: FS_SaveReference'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END
      

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END OF outsavtec                            !
!                             Jan 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
