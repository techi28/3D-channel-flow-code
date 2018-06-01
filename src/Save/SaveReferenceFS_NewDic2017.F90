!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   OPEN DATA AT REFERENCE POINTS                     !
!                              Dic 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

    SUBROUTINE FS_SaveReferenceOpen

#     include "cppdefs.h"
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
#     endif
!     =============== END ================    
!     ====================================
      implicit none
      integer :: UseThis

      UseThis = 1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if (rang_topo.ne.0)  UseThis = 0
#     endif
!     =============== END ================    
!     ====================================

      IF (UseThis.eq.1) THEN
         print*,'  '
         print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
         print*,'  '
         print*,' Saved: ../output/FS/FS_SolTime1.dat'
         print*,' Saved: ../output/FS/FS_SolTime2.dat'
         print*,' Saved: ../output/FS/FS_SolTimeE1.dat'
         print*,' Saved: ../output/FS/FS_SolTimeE2.dat'   
         open(4101,file="../output/FS/FS_SolTime1.dat",status='unknown')      
         open(4102,file="../output/FS/FS_SolTime2.dat",status='unknown')
         open(4103,file="../output/FS/FS_SolTimeE1.dat",status='unknown')
         open(4104,file="../output/FS/FS_SolTimeE2.dat",status='unknown')
!        -----------------------
#        ifdef KeySave1DRef3Points
         print*,' Saved: ../output/FS/FS_SolTime3.dat'
         open(4105,file="../output/FS/FS_SolTime3.dat",status='unknown')
#        endif
         print*,'  '
         print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
         print*,'  '
      ENDIF

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   CLOSE DATA AT REFERENCE POINTS                    !
!                              Dic 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


    SUBROUTINE FS_SaveReferenceClose

#     include "cppdefs.h"
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
#     endif
!     =============== END ================    
!     ====================================
      implicit none
      integer :: UseThis

      UseThis = 1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if (rang_topo.ne.0)  UseThis = 0
#     endif
!     =============== END ================    
!     ====================================

      IF (UseThis.eq.1) THEN
         close(4101)
         close(4103)
         close(4102)
         close(4104)    
!        -----------------------
#        ifdef KeySave1DRef3Points      
         close(9105)
#        endif 
      ENDIF

      RETURN
      END     
    
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
                      No_sp,nstep) 
        
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
      integer :: nstep
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8, dimension(:) :: xvR1(N_SPMAX),yvR1(N_SPMAX)
      integer,dimension(:) :: OpenFiles(N_SPMAX)
!     --------------------------------------
      real*8 :: FS_funeta,FS_funH     
      real*8 :: FS_funu,FS_funv,FS_funw,FS_funp
      real*8 :: x,y,z
      real*8 :: DR,DR0
      integer :: nv1,pro
      integer :: n,s,NumS
      integer :: UseThis
!     --------------------------------------
      real*8,dimension(:),allocatable :: phiufv
      real*8,dimension(:),allocatable :: phivfv
      real*8,dimension(:),allocatable :: phiwfv     
      real*8,dimension(:),allocatable :: phipfv
      real*8,dimension(:),allocatable :: phiufvA
      real*8,dimension(:),allocatable :: phivfvA
      real*8,dimension(:),allocatable :: phiwfvA   
      real*8,dimension(:),allocatable :: phipfvA
      real*8 :: etavR1,etavAR1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      real*8, dimension(Nprocs) :: DR0V
      integer,dimension(Nprocs) :: nv1V,proV,UseProc
#     endif
!     =============== END ================    
!     ====================================
!      ____________________________________
!     |                                    |
!     |             Parameter              |
!     |____________________________________|

      integer, parameter :: EachStep = 10
    
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

      UseThis = 1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if (rang_topo.ne.0)  UseThis = 0
#     endif
!     =============== END ================    
!     ====================================      
!      ________________________________________________________
!     |                                                        |
!     |                       Exact solution                   |
!     |________________________________________________________|

      do nv=1,N_VERT
         etavA = etav
         HprvA = Hprv
         ufvA  = ufv
         vfvA  = vfv
         wfvA  = wfv
         pfvA  = pfv
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
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            pro = rang_topo
!           -------------------------
!           Results of all processors
            call MPI_ALLGATHER(pro,1,MPI_INTEGER,                    &
                               proV,1,MPI_INTEGER,comm3D,code)
            call MPI_ALLGATHER(nv1,1,MPI_INTEGER,                    &
                               nv1V,1,MPI_INTEGER,comm3D,code)
            call MPI_ALLGATHER(DR0,1,MPI_DOUBLE_PRECISION,           &
                               DR0V,1,MPI_DOUBLE_PRECISION,comm3D,code)
            call MPI_Barrier(comm3D,code)
!           -------------------------
!           Minimum among processors
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
#        endif
!        =============== END ================    
!        ====================================
         No_sp(s) = nv1
      ENDDO
      
!      ________________________________________________________
!     |                                                        |
!     |                       Saving data                      |
!     |________________________________________________________|
             
      IF ((0.eq.mod(nstep,EachStep)).or.(abs(time-tFin).lt.1d-5)) THEN
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel

      allocate(phiufv(NZ-1),  &
               phivfv(NZ-1),  &
               phiwfv(NZ-1),  &   
               phipfv(NZ-1),  &
               phiufvA(NZ-1), &
               phivfvA(NZ-1), &
               phiwfvA(NZ-1), &    
               phipfvA(NZ-1))
               
      DO s=1,NumS
         nv = No_sp(s)
         pro = UseProc(s)
!        -------------------------
!        Re-define variables             
         if (rang_topo.eq.pro) then
            etavR1  = etav(nv)
            phiufv  = ufv(nv,1:NZ-1)
            phivfv  = vfv(nv,1:NZ-1)
            phiwfv  = wfv(nv,1:NZ-1)
            phipfv  = pfv(nv,1:NZ-1)
!           ------
            etavAR1 = etavA(nv)
            phiufvA = ufvA(nv,1:NZ-1)
            phivfvA = vfvA(nv,1:NZ-1)
            phiwfvA = wfvA(nv,1:NZ-1)
            phipfvA = pfvA(nv,1:NZ-1)
         else
            etavR1  = 0.0d0
            phiufv  = 0.0d0
            phivfv  = 0.0d0
            phiwfv  = 0.0d0
            phipfv  = 0.0d0
!           ------
            etavAR1 = 0.0d0
            phiufvA = 0.0d0
            phivfvA = 0.0d0
            phiwfvA = 0.0d0
            phipfvA = 0.0d0
         endif
!        -------------------------
!        Communicate to proc=0
         call MPI_Bcast(etavR1,1,MPI_DOUBLE_PRECISION,pro,comm3D,code)
         call MPI_Bcast(phiufv,1,type_vertical,pro,comm3D,code)        
         call MPI_Bcast(phivfv,1,type_vertical,pro,comm3D,code) 
         call MPI_Bcast(phiwfv,1,type_vertical,pro,comm3D,code) 
         call MPI_Bcast(phipfv,1,type_vertical,pro,comm3D,code)
!        ------
         call MPI_Bcast(etavAR1,1,MPI_DOUBLE_PRECISION,pro,comm3D,code)
         call MPI_Bcast(phiufvA,1,type_vertical,pro,comm3D,code)        
         call MPI_Bcast(phivfvA,1,type_vertical,pro,comm3D,code) 
         call MPI_Bcast(phiwfvA,1,type_vertical,pro,comm3D,code) 
         call MPI_Bcast(phipfvA,1,type_vertical,pro,comm3D,code)
!        -------------------------
!        Save results    
         if (UseThis.eq.1) then
!           ------------
!           Save point 1 
            if (s.eq.1) then
            print*,'    ===> FS_SaveReference used for nstep = ',nstep
            write(4101,*) time,xvR1(s),yvR1(s),etavR1,& 
                          phiufv(1:NZ-1), &
                          phivfv(1:NZ-1), &
                          phiwfv(1:NZ-1), &
                          phipfv(1:NZ-1)
            write(4103,*) time,xvR1(s),yvR1(s),etavAR1,&
                          phiufvA(1:NZ-1), &
                          phivfvA(1:NZ-1), &
                          phiwfvA(1:NZ-1), &
                          phipfvA(1:NZ-1)
!           ------------
!           Save point 2        
            elseif (s.eq.2) then
            write(4102,*) time,xvR1(s),yvR1(s),etavR1,&
                          phiufv(1:NZ-1), &
                          phivfv(1:NZ-1), &
                          phiwfv(1:NZ-1), &
                          phipfv(1:NZ-1)
            write(4104,*) time,xvR1(s),yvR1(s),etavAR1,&
                          phiufvA(1:NZ-1), &
                          phivfvA(1:NZ-1), &
                          phiwfvA(1:NZ-1), &
                          phipfvA(1:NZ-1)
!           ------------
!           Save point 3
#           ifdef KeySave1DRef3Points
            elseif (s.eq.3) then
            write(4105,*) time,xvR1(s),yvR1(s),etavR1,&
                          phiufv(1:NZ-1), &
                          phivfv(1:NZ-1), &
                          phiwfv(1:NZ-1), &
                          phipfv(1:NZ-1)

#           endif
            endif
            call MPI_Barrier(comm3D,code)
         endif


      ENDDO 
      
      deallocate(phiufv,  &
                 phivfv,  &
                 phiwfv,  &   
                 phipfv,  &
                 phiufvA, &
                 phivfvA, &
                 phiwfvA, &    
                 phipfvA)
      
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
         print*,'    ===> FS_SaveReference used for nstep = ',nstep
         nv = No_sp(1)
         write(4101,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1), &
                       vfv(nv,1:NZ-1), &
                       wfv(nv,1:NZ-1), &
                       pfv(nv,1:NZ-1)
         write(4103,*) time,xv(nv),yv(nv),etavA(nv),&
                       ufvA(nv,1:NZ-1), &
                       vfvA(nv,1:NZ-1), &
                       wfvA(nv,1:NZ-1), &
                       pfvA(nv,1:NZ-1)
!        -------------
         nv = No_sp(2)
         write(4102,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1), &
                       vfv(nv,1:NZ-1), &
                       wfv(nv,1:NZ-1), &
                       pfv(nv,1:NZ-1)
         write(4104,*) time,xv(nv),yv(nv),etavA(nv),&
                       ufvA(nv,1:NZ-1), &
                       vfvA(nv,1:NZ-1), &
                       wfvA(nv,1:NZ-1), &
                       pfvA(nv,1:NZ-1)
!        -------------
#        ifdef KeySave1DRef3Points
         nv = No_sp(3)
         write(4105,*) time,xv(nv),yv(nv),etav(nv),&
                       ufv(nv,1:NZ-1), &
                       vfv(nv,1:NZ-1), &
                       wfv(nv,1:NZ-1), &
                       pfv(nv,1:NZ-1)
#        endif
#     endif
!     =============== END ================    
!     ====================================
      ENDIF
 
                     
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
