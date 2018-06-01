!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  WRITE DATA TO DISPLAY AT PARAVIEW                  !
!                             Dic 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_SaveParaview(Hprv,etav,                &
                                 ufv,vfv,wfv,pfv,          &
                                 HprvA,etavA,              &
                                 ufvA,vfvA,wfvA,pfvA,      &
                                 uErrv,vErrv,wErrv,pErrv,  &
                                 xvt,yvt,zvt,No_vp)

!---------------------------------------------------------------------!
!                                                                     !
!    This program writes the main variables to display in the         !
!    program PARAVIEW.                                                ! 
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |     Size    |        Description                 |  !  
!  |____________|_____________|____________________________________|  !
!  | --> ufv    |(N_VERT,NZ-1)| Velocity component u_f             |  !
!  | --> vfv    |(N_VERT,NZ-1)| Velocity component v_f             |  !      
!  | --> wfv    |(N_VERT,NZ-1)| Velocity component w_f             |  !
!  | --> pfv    |(N_VERT,NZ-1)| Pressure of the fluid              |  !
!  |____________|_____________|____________________________________|  !
!  | --> usvA   |(N_VERT,NZ-1)| Velocity component exact           |  !
!  | --> vsvA   |(N_VERT,NZ-1)| Velocity component exact           |  !     
!  | --> wsvA   |(N_VERT,NZ-1)| Velocity component exact           |  !
!  | --> psvA   |(N_VERT,NZ-1)| Pressure of the fluid exact        |  !
!  |____________|_____________|____________________________________|  !
!  | --> xvt    |(N_VERT,NZ-1)| xc at the current time             |  !
!  | --> yvt    |(N_VERT,NZ-1)| yc at the current time             |  !
!  | --> zvt    |(N_VERT,NZ-1)| yc at the current time             |  !
!  |____________|_____________|____________________________________|  !
!  | --> No_vp  |(N_VERT,3 )  | Numbering of cell vertices         |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !  
!  |_____________|_________________________________________________|  ! 
!  | --- N_CELL  | Total number of the cells                       |  !
!  | --- N_CELL0 | Inside number of cells                          |  !
!  | --- N_VERT  | Total number of vertices                        |  !
!  | --- NZ      | Points in the sigma direction                   |  !
!  |   time      | Current simulation time                         |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name         |                 Description                  |  !  
!  |________________|______________________________________________|  !
!  | * i,k          |  Loop counters                               |  !
!  | * SaveCounter  |  Integer counter to save results             |  !    
!  |________________|______________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name          |                 Description                 |  !  
!  |_________________|_____________________________________________|  !
!  | irec            |  Writing file                               |  !
!  | TotalN_VERT     |  Total number de vertex in the 3D domain    |  !
!  | TotalN_CELL0    |  Total number of prism elements             |  !
!  | nv1B,nv2B,nv3B  |  Bottom vertex indices of an element        |  !
!  | nv1T,nv2T,nv3T  |  Top vertex indices of an element           |  !           
!  |_________________|_____________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   ---  Parameters                                                   !
!    *   Common variables modified                                    !
!        Common variables used                                        !
!---------------------------------------------------------------------!

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

      real*8,dimension(:)   :: etav(N_VERT)
      real*8,dimension(:)   :: Hprv(N_VERT)

      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)

      real*8,dimension(:)   :: etavA(N_VERT)
      real*8,dimension(:)   :: HprvA(N_VERT)        

      real*8,dimension(:,:) :: ufvA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfvA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfvA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfvA(N_VERT,NZ-1)

      real*8,dimension(:,:) :: uErrv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vErrv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wErrv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pErrv(N_VERT,NZ-1)  

      real*8,dimension(:,:) :: xvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: yvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: zvt(N_VERT,NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: etav_gl,Hprv_gl
      real*8,dimension(:,:),allocatable :: etavA_gl,HprvA_gl
      real*8,dimension(:,:),allocatable :: xvt_gl,yvt_gl,zvt_gl
      real*8,dimension(:,:),allocatable :: ufv_gl,vfv_gl,wfv_gl
      real*8,dimension(:,:),allocatable :: ufvA_gl,vfvA_gl,wfvA_gl
      real*8,dimension(:,:),allocatable :: pfv_gl
      real*8,dimension(:,:),allocatable :: pfvA_gl
      real*8,dimension(:,:) :: etavt(N_VERT,NZ-1),etavtA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hprvt(N_VERT,NZ-1),HprvtA(N_VERT,NZ-1)
      
!     ----------------------------------------
      integer:: irec,vert,elem
      integer:: nv1B,nv2B,nv3B,nv1T,nv2T,nv3T
      integer:: Dothis,PrintThis
      character*50 filen
      character*50 filenP
      character*50 filenPP
!     ----------------------------------------
      integer :: TotalN_VERT
      integer :: TotalN_ELEM
!     ----------------------------------------

      TotalN_VERT = N_VERTglobal*(NZglobal-1) 
      TotalN_ELEM = N_CELL0global*(NZglobal-2)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: FS_SaveParaview'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                        Formats                         |
!     |________________________________________________________|

 1    format(1(1x,e12.5))
 2    format(3(1x,e12.5))
 3    format(3(1x,i8))
 6    format(7(1x,i8))

 100  format('# vtk DataFile Version 2.0')
 110  format('Unstructured Grid Example')
 120  format('ASCII')
 130  format('DATASET UNSTRUCTURED_GRID')


#     ifndef KeyParallel
!        ________________________________________________________
!       |                                                        |
!       |            Save Mesh with eta (for Matlab)             |
!       |________________________________________________________|
             
         if (SaveCounter.eq.0) then
            open(1000,file='../output/Matlab/PlotMesh/Mesh.txt')
            write(1000,*) N_VERT
            write(1000,*) N_CELL0
            do nv=1,N_VERTglobal
               write(1000,*) xvt(nv,1),yvt(nv,1),etav(nv)       
            enddo
            do i=1,N_CELL0global
               nv1B = No_vp(i,1)
               nv2B = No_vp(i,2)
               nv3B = No_vp(i,3)          
               write(1000,*) nv1B,nv2B,nv3B  
            enddo           
            close(1000)
         endif
         
!        ________________________________________________________
!       |                                                        |
!       |          *********  Write data file  ************      |
!       |________________________________________________________|

!        __________________________________
!       |                                  |
!       |          Open file irec          |
!       |__________________________________| 
 
         irec=60
         filen='../output/Serial/V-    .vtk'
         write(filen(20:23),'(i4.4)') SaveCounter
         open(irec,file=filen)
         write(irec,100)
         write(irec,110)
         write(irec,120)
         write(irec,130)

!        __________________________________
!       |                                  |
!       |           Write points           |
!       |__________________________________| 

         write(irec,*) 'POINTS ',TotalN_VERT,' float'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,2) xvt(nv,k),yvt(nv,k),zvt(nv,k)       
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write cells            |
!       |__________________________________| 

         write(irec,*) 'CELLS ',TotalN_ELEM,(TotalN_ELEM+6*TotalN_ELEM) 
         do k=1,NZ-2
            do i=1,N_CELL0
               nv1B = No_vp(i,1) + N_VERTglobal*(k-1) - 1
               nv2B = No_vp(i,2) + N_VERTglobal*(k-1) - 1
               nv3B = No_vp(i,3) + N_VERTglobal*(k-1) - 1
               nv1T = No_vp(i,1) + N_VERTglobal*k - 1 
               nv2T = No_vp(i,2) + N_VERTglobal*k - 1
               nv3T = No_vp(i,3) + N_VERTglobal*k - 1           
               write(irec,6) 6,nv1B,nv2B,nv3B,nv1T,nv2T,nv3T  
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |         CELL types (prism)       |
!       |__________________________________| 

         write(irec,*) 'CELL_TYPES ',TotalN_ELEM
         do i=1,TotalN_ELEM
            write(irec,3) 13       
         enddo
!        __________________________________
!       |                                  |
!       |           Write scalars          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'POINT_DATA ',TotalN_VERT
         write(irec,*) 'SCALARS eta    float'
         write(irec,*) 'LOOKUP_TABLE default'      
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal        
               write(irec,1) etav(nv) 
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS eta_exact float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal        
               write(irec,1) etavA(nv) 
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS Hpr    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal    
               write(irec,1) Hprv(nv) 
            enddo
         enddo
         !-------- 
         write(irec,*) ' '
         write(irec,*) 'SCALARS Hpr_exact    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal    
               write(irec,1) HprvA(nv) 
            enddo
         enddo  
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS pf    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal        
               write(irec,1) pfv(nv,k) 
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS pf_exact float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal        
               write(irec,1) pfvA(nv,k) 
            enddo
         enddo                  
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'VECTORS vf float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               if (abs(ufv(nv,k)).le.1e-14) ufv(nv,k) = 0.0d0
               if (abs(vfv(nv,k)).le.1e-14) vfv(nv,k) = 0.0d0
               if (abs(wfv(nv,k)).le.1e-14) wfv(nv,k) = 0.0d0
               write(irec,2) ufv(nv,k),vfv(nv,k),wfv(nv,k)        
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'VECTORS vf_exact float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufvA(nv,k),vfvA(nv,k),wfvA(nv,k)        
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Close file irec        |
!       |__________________________________|    

         rewind(irec)
         close(irec)

!        ________________________________________________________
!       |                                                        |
!       |                       Display                          |
!       |________________________________________________________|

!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!#        ifdef KeyDisplay
         write(*,'(t20,a20)')      ' __________________ '
         write(*,'(t20,a20)')      '|Save Results Parav|'
         write(*,'(t20,a20)')      '|__________________|'
         write(*,'(t20,a7,i3)')    '  No.  :',SaveCounter
         write(*,'(t20,a7,e10.3)') '  time :',time
         write(*,'(t20,a20)')      '|__________________|'
         print*,'  '
!#        endif
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#     endif

!*********************************************************************!
!                                                                     !
!                ====================================                 !
!                =====  START PARALLEL OPTION =======                 !
!                                                                     !
!*********************************************************************!

#     ifdef KeyParallel

!        ________________________________________________________
!       |                                                        |
!       |                 Allocate global matrix                 |
!       |________________________________________________________|   

        allocate(xvt_gl(N_VERTglobal,NZglobal-1), &
                 yvt_gl(N_VERTglobal,NZglobal-1), &
                 zvt_gl(N_VERTglobal,NZglobal-1), &
                 etav_gl(N_VERTglobal,NZglobal-1),&
                 Hprv_gl(N_VERTglobal,NZglobal-1),&                
                 ufv_gl(N_VERTglobal,NZglobal-1), &
                 vfv_gl(N_VERTglobal,NZglobal-1), &
                 wfv_gl(N_VERTglobal,NZglobal-1), &
                 pfv_gl(N_VERTglobal,NZglobal-1), & 
                 etavA_gl(N_VERTglobal,NZglobal-1),&
                 HprvA_gl(N_VERTglobal,NZglobal-1),&              
                 ufvA_gl(N_VERTglobal,NZglobal-1), &
                 vfvA_gl(N_VERTglobal,NZglobal-1), &
                 wfvA_gl(N_VERTglobal,NZglobal-1), &
                 pfvA_gl(N_VERTglobal,NZglobal-1))
!        ________________________________________________________
!       |                                                        |
!       |               Reconstruct global matrix                |
!       |________________________________________________________|  

         do k=1,NZ-1
            do nv=1,N_VERT
               etavt(nv,k)  = etav(nv)
               etavtA(nv,k) = etavA(nv)
               Hprvt(nv,k)  = Hprv(nv)
               HprvtA(nv,k) = HprvA(nv)      
            enddo
         enddo
         
         call matgloV(etavt,etav_gl)
         call matgloV(etavtA,etavA_gl)
         call matgloV(Hprvt,Hprv_gl)
         call matgloV(HprvtA,HprvA_gl)         
         call matgloV(xvt,xvt_gl)
         call matgloV(yvt,yvt_gl)
         call matgloV(zvt,zvt_gl)
         call matgloV(ufv,ufv_gl)
         call matgloV(vfv,vfv_gl)
         call matgloV(wfv,wfv_gl)
         call matgloV(pfv,pfv_gl)
         call matgloV(ufvA,ufvA_gl)
         call matgloV(vfvA,vfvA_gl)
         call matgloV(wfvA,wfvA_gl)
         call matgloV(pfvA,pfvA_gl)


         IF (rang_topo.eq.0) THEN
!        ________________________________________________________
!       |                                                        |
!       |            Save Mesh with eta (for Matlab)             |
!       |________________________________________________________|
             
         if (SaveCounter.eq.0) then
            open(1000,file='../output/Matlab/PlotMesh/Mesh.txt')
            write(1000,*) N_VERTglobal
            write(1000,*) N_CELL0global
            do nv=1,N_VERTglobal
               write(1000,*) xvt_gl(nv,1),yvt_gl(nv,1),etav_gl(nv,1)       
            enddo
            do i=1,N_CELL0global
               nv1B = No_vp_global(i,1)
               nv2B = No_vp_global(i,2)
               nv3B = No_vp_global(i,3)          
               write(1000,*) nv1B,nv2B,nv3B  
            enddo           
            close(1000)
         endif  
!        ________________________________________________________
!       |                                                        |
!       |     *********  Write data file  (global) ************  |
!       |________________________________________________________|
       
!        __________________________________
!       |                                  |
!       |          Open file irec          |
!       |__________________________________| 

         irec=60
         filenP='../output/Parallel/V-    .vtk'
         write(filenP(22:25),'(i4.4)') SaveCounter
         open(irec,file=filenP)
         write(irec,100)
         write(irec,110)
         write(irec,120)
         write(irec,130) 
!        __________________________________
!       |                                  |
!       |           Write points           |
!       |__________________________________| 

         write(irec,*) 'POINTS ',TotalN_VERT,' float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) xvt_gl(nv,k),yvt_gl(nv,k),zvt_gl(nv,k)       
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write cells            |
!       |__________________________________| 

         write(irec,*) 'CELLS ',TotalN_ELEM,(TotalN_ELEM+6*TotalN_ELEM) 
         do k=1,NZglobal-2
            do i=1,N_CELL0global
               nv1B = No_vp_global(i,1) + N_VERTglobal*(k-1) - 1
               nv2B = No_vp_global(i,2) + N_VERTglobal*(k-1) - 1
               nv3B = No_vp_global(i,3) + N_VERTglobal*(k-1) - 1
               nv1T = No_vp_global(i,1) + N_VERTglobal*k - 1
               nv2T = No_vp_global(i,2) + N_VERTglobal*k - 1
               nv3T = No_vp_global(i,3) + N_VERTglobal*k - 1           
               write(irec,6) 6,nv1B,nv2B,nv3B,nv1T,nv2T,nv3T  
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |         CELL types (prism)       |
!       |__________________________________| 

         write(irec,*) 'CELL_TYPES ',TotalN_ELEM
         do i=1,TotalN_ELEM
            write(irec,3) 13       
         enddo
!        __________________________________
!       |                                  |
!       |           Write scalars          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'POINT_DATA ',TotalN_VERT                 
         !--------
         write(irec,*) ' '                  
         write(irec,*) 'SCALARS pf     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) pfv_gl(nv,k)      
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS pf_exact     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) pfvA_gl(nv,k)        
            enddo
         enddo
         !--------         
         write(irec,*) 'SCALARS eta     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) etav_gl(nv,k)      
            enddo
         enddo  
         !--------
         write(irec,*) ' '                  
         write(irec,*) 'SCALARS eta_exact     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) etavA_gl(nv,k)      
            enddo
         enddo           
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'VECTORS vf float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufv_gl(nv,k),vfv_gl(nv,k),wfv_gl(nv,k)        
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'VECTORS vf_exact float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufvA_gl(nv,k),vfvA_gl(nv,k),wfvA_gl(nv,k)        
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |          Close file irec         |
!       |__________________________________|  

         rewind(irec)
         close(irec)           

         ENDIF
!        ________________________________________________________
!       |                                                        |
!       | ********** Write data file  (each processor) ********  |
!       |________________________________________________________|

         PrintThis = 0
         IF (PrintThis.eq.1) THEN
         
         TotalN_VERT = N_VERT*(NZ-1) 
         TotalN_ELEM = N_CELL0*(NZ-2)
!        __________________________________
!       |                                  |
!       |          Open file irec          |
!       |__________________________________| 
 
         irec=70
         filenPP='../output/Parallel/VMPI  -    .vtk'
         write(filenPP(24:25),'(i2.2)') rang_topo
         write(filenPP(27:30),'(i4.4)') SaveCounter
         open(irec,file=filenPP)
         write(irec,100)
         write(irec,110)
         write(irec,120)
         write(irec,130) 
!        __________________________________
!       |                                  |
!       |           Write points           |
!       |__________________________________| 

         write(irec,*) 'POINTS ',TotalN_VERT,' float'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,2) xvt(nv,k),yvt(nv,k),zvt(nv,k)       
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write cells            |
!       |__________________________________| 

         write(irec,*) 'CELLS ',TotalN_ELEM,(TotalN_ELEM+6*TotalN_ELEM) 
         do k=1,NZ-2
            do i=1,N_CELL0
               nv1B = No_vp(i,1) + N_VERT*(k-1) - 1
               nv2B = No_vp(i,2) + N_VERT*(k-1) - 1
               nv3B = No_vp(i,3) + N_VERT*(k-1) - 1
               nv1T = No_vp(i,1) + N_VERT*k - 1
               nv2T = No_vp(i,2) + N_VERT*k - 1
               nv3T = No_vp(i,3) + N_VERT*k - 1           
               write(irec,6) 6,nv1B,nv2B,nv3B,nv1T,nv2T,nv3T  
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |         CELL types (prism)       |
!       |__________________________________| 

         write(irec,*) 'CELL_TYPES ',TotalN_ELEM
         do i=1,TotalN_ELEM
            write(irec,3) 13       
         enddo
!        __________________________________
!       |                                  |
!       |           Write scalars          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'POINT_DATA ',TotalN_VERT
         write(irec,*) 'SCALARS pf     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,1) pfv(nv,k)      
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS ps     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,1) pfvA(nv,k)        
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_f float'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,2) ufv(nv,k),vfv(nv,k),wfv(nv,k)        
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_s float'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,2) ufvA(nv,k) ,vfvA(nv,k) ,wfvA(nv,k)        
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |          Close file irec         |
!       |__________________________________|  

         rewind(irec)
         close(irec)  
        
         ENDIF
!        ________________________________________________________
!       |                                                        |
!       |                Deallocate global matrix                |
!       |________________________________________________________|   

        deallocate(xvt_gl,yvt_gl,zvt_gl,&
                   etav_gl,Hprv_gl,&
                   ufv_gl,vfv_gl,wfv_gl,pfv_gl,&
                   etavA_gl,HprvA_gl,&
                   ufvA_gl,vfvA_gl,wfvA_gl,pfvA_gl)

!        ________________________________________________________
!       |                                                        |
!       |                       Display                          |
!       |________________________________________________________|

!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!#        ifdef KeyDisplay
         if(rang_topo.eq.0) then
         write(*,'(t20,a20)')      ' __________________ '
         write(*,'(t20,a20)')      '|Save Results Parav|'
         write(*,'(t20,a20)')      '|__________________|'
         write(*,'(t20,a7,i3)')    '  No.  :',SaveCounter
         write(*,'(t20,a7,e10.3)') '  time :',time
         write(*,'(t20,a20)')      '|__________________|'
         print*,'  '
         endif
!#        endif
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#     endif


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: FS_SaveParaview'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          END OF PARAVIEW                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

