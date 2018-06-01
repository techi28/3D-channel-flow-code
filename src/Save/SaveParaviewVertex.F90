!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  WRITE DATA TO DISPLAY AT PARAVIEW                  !
!                             March 2017                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SaveParaviewVertex(Hprv,etav,                &
                                    ufv,vfv,wfv,pfv,          &
                                    usv,vsv,wsv,psv,          &
                                    uErrv,vErrv,wErrv,pErrv,  &
                                    xvt,yvt,zvt,              &
                                    No_vp)

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
!  | --> ufv    |(N_VERT,NZ-1)| Approx. Velocity component u_f     |  !
!  | --> vfv    |(N_VERT,NZ-1)| Approx. Velocity component v_f     |  !      
!  | --> wfv    |(N_VERT,NZ-1)| Approx. Velocity component w_f     |  !
!  | --> pfv    |(N_VERT,NZ-1)| Approx. Pressure of the fluid      |  !
!  |____________|_____________|____________________________________|  !
!  | --> usv    |(N_VERT,NZ-1)| Exact Velocity component u_s       |  !
!  | --> vsv    |(N_VERT,NZ-1)| Exact Velocity component v_s       |  !     
!  | --> wsv    |(N_VERT,NZ-1)| Exact Velocity component w_s       |  !
!  | --> psv    |(N_VERT,NZ-1)| Exact Pressure of the solid        |  !
!  |____________|_____________|____________________________________|  !
!  | --> uErrv  |(N_VERT,NZ-1)| Error Velocity component u_s       |  !
!  | --> vErrv  |(N_VERT,NZ-1)| Error Velocity component v_s       |  !     
!  | --> wErrv  |(N_VERT,NZ-1)| Error Velocity component w_s       |  !
!  | --> pErrv  |(N_VERT,NZ-1)| Error Pressure of the solid        |  !
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

      real*8,dimension(:)   :: Hprv
      real*8,dimension(:)   :: etav
      
      real*8,dimension(:,:) :: ufv
      real*8,dimension(:,:) :: vfv
      real*8,dimension(:,:) :: wfv
      real*8,dimension(:,:) :: pfv  

      real*8,dimension(:,:) :: usv
      real*8,dimension(:,:) :: vsv
      real*8,dimension(:,:) :: wsv
      real*8,dimension(:,:) :: psv
 
      real*8,dimension(:,:) :: uErrv
      real*8,dimension(:,:) :: vErrv
      real*8,dimension(:,:) :: wErrv
      real*8,dimension(:,:) :: pErrv    

      real*8,dimension(:,:) :: xvt
      real*8,dimension(:,:) :: yvt
      real*8,dimension(:,:) :: zvt
      
      integer,dimension(:,:):: No_vp
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:) :: Magnv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: etavt(N_VERT,NZ-1),etavtA(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Hprvt(N_VERT,NZ-1),HprvtA(N_VERT,NZ-1)
!     ----------------------------------------      
      real*8,dimension(:,:),allocatable :: etav_gl,Hprv_gl
      real*8,dimension(:,:),allocatable :: etavA_gl,HprvA_gl
      real*8,dimension(:,:),allocatable :: xvt_gl,yvt_gl,zvt_gl
      real*8,dimension(:,:),allocatable :: ufv_gl,vfv_gl,wfv_gl,pfv_gl
      real*8,dimension(:,:),allocatable :: usv_gl,vsv_gl,wsv_gl,psv_gl
      real*8,dimension(:,:),allocatable :: uErrv_gl,vErrv_gl,wErrv_gl
      real*8,dimension(:,:),allocatable :: pErrv_gl
      real*8,dimension(:,:),allocatable :: Magnv_gl
      real*8,dimension(:,:),allocatable :: Vortv_gl
      real*8,dimension(:,:),allocatable :: Vortxv_gl,Vortyv_gl,Vortzv_gl       
!     ----------------------------------------
      integer:: irec,vert,elem
      integer:: nv1B,nv2B,nv3B,nv1T,nv2T,nv3T
      integer:: Dothis,PrintThis
      character*50 filen
      character*50 filenP
      character*50 filenPP
!     ----------------------------------------      
      integer :: kmid,nvmid
      real*8  :: cons      
!     ----------------------------------------
      integer :: TotalN_VERT
      integer :: TotalN_ELEM
!     ----------------------------------------
      integer :: Write_eta,Write_etaEx,Write_etaError
      integer :: Write_p,Write_pEx,Write_pError
      integer :: Write_vel,Write_velEx,Write_velError
      integer :: Write_Magni_vel
      integer :: Write_Magni_vort
      integer :: Write_vorticity

      TotalN_VERT = N_VERTglobal*(NZglobal-1) 
      TotalN_ELEM = N_CELL0global*(NZglobal-2)

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: SaveParaviewVertex'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                      Test cases                        |
!     |________________________________________________________|

         Write_eta       = 0
         Write_etaEx     = 0
         Write_etaError  = 0
         Write_p         = 0
         Write_pEx       = 0
         Write_pError    = 0
         Write_vel       = 0
         Write_velEx     = 0
         Write_velError  = 0
         Write_Magni_vel = 0
         Write_Magni_vort= 0
         Write_vorticity = 0
         
!     ------------------------------------
#     if defined(KeyEstuaryGironde)         
         Write_eta       = 1
         Write_p         = 1
         Write_vel       = 1
         Write_Magni_vel = 1
         Write_Magni_vort= 1
         Write_vorticity = 1
!     ------------------------------------         
#     elif defined(KeyStaticCylinder)
         Write_eta       = 1
         Write_p         = 1
         Write_vel       = 1
         Write_Magni_vel = 1
         Write_Magni_vort= 1
         Write_vorticity = 1
!     ------------------------------------         
#     elif defined(KeyStaticChannel)
         Write_p         = 1
         Write_vel       = 1
         Write_Magni_vel = 1
         Write_Magni_vort= 1
         Write_vorticity = 1
!     ------------------------------------
#     elif defined(KeyStandingWave)         
         Write_eta       = 1
         Write_etaEx     = 1
         Write_etaError  = 1
         Write_p         = 1
         Write_pEx       = 1
         Write_pError    = 1
         Write_vel       = 1
         Write_velEx     = 1
         Write_velError  = 1
!     ------------------------------------
#     elif defined(KeyTaylorVortex)
         Write_p         = 1
         Write_pEx       = 1
         Write_pError    = 1
         Write_vel       = 1
         Write_velEx     = 1
         Write_velError  = 1
         Write_Magni_vel = 1
         Write_Magni_vort= 1
         Write_vorticity = 1
!     ------------------------------------
#     elif defined(KeyTestOnlyPoisson)
         Write_p         = 1
         Write_pEx       = 1
         Write_pError    = 1
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                    Other quantities                    |
!     |________________________________________________________|

!     ________________________________________________________
!     Velocity magnitude
      do k=1,NZ-1
         do nv=1,N_VERT
            Magnv(nv,k)=sqrt(ufv(nv,k)**2+vfv(nv,k)**2+wfv(nv,k)**2)
         enddo
      enddo
      
!     ________________________________________________________
!     Difference constant between analytical and numerical  
#     ifdef Key_NeumannBCp
         if (time.gt.0) then
            kmid  = floor(NZ/2.0d0)
            nvmid = floor(N_VERT/2.0d0)
            cons = pfv(nvmid,kmid)-psv(nvmid,kmid)       
            do k=1,NZ-1
               do nv=1,N_VERT
                  !pfv(nv,k)   = pfv(nv,k) - cons
                  !pErrv(nv,k) = abs(pfv(nv,k)-psv(nv,k)-cons)
               enddo
            enddo
         endif      
#     endif

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

	
!*********************************************************************!
!                                                                     !
!                ====================================                 !
!                ==========  SEQUENTIAL =============                 !
!                                                                     !
!*********************************************************************!

#     ifndef KeyParallel

!        ________________________________________________________
!       |                                                        |
!       |            Save Mesh with eta (for Matlab)             |
!       |________________________________________________________|
             
         if (SaveCounter.eq.0) then
            print*,'  '
            print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
            print*,'  '
            print*,' Saved: ../output/Matlab/PlotMesh/Mesh.txt'
            open(1000,file='../output/Matlab/PlotMesh/Mesh.txt')
            write(1000,*) N_VERT
            write(1000,*) N_CELL0
            do nv=1,N_VERT
               write(1000,*) xvt(nv,1),yvt(nv,1),etav(nv)       
            enddo
            do i=1,N_CELL0
               nv1B = No_vp(i,1)
               nv2B = No_vp(i,2)
               nv3B = No_vp(i,3)          
               write(1000,*) nv1B,nv2B,nv3B  
            enddo           
            close(1000)
            print*,'  '
            print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
         endif

!        ________________________________________________________
!       |                                                        |
!       |               Write geometry data files                |
!       |________________________________________________________|

         Dothis = 0
         if (Dothis.eq.1) then
            open(111,file="../output/xv.dat",status='unknown')
            open(112,file="../output/yv.dat",status='unknown')
            open(113,file="../output/zv.dat",status='unknown')
            open(114,file='../output/phiVertex.dat')
            do nv=1,N_VERT
               write(111,6) xvt(nv,1)
               write(112,6) yvt(nv,1)
               write(113,6) zvt(nv,1)
               write(114,6) ufv(nv,3)
            enddo
            close(111)
            close(112)
            close(113)
            close(114)
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
!       |            POINT_DATA            |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'POINT_DATA ',TotalN_VERT
!        __________________________________
!       |                                  |
!       |           Write scalars          |
!       |__________________________________| 

         if (Write_eta.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS eta    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal        
               write(irec,1) etav(nv) 
            enddo
         enddo
         endif
         !--------
         if (Write_p.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS p     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) pfv(nv,k)      
            enddo
         enddo
         endif
         !--------
         if (Write_pEx.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS p_Ex     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) psv(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_pError.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS p_Error     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) pErrv(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_Magni_vel.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS Magni_vel     float'
         write(irec,*) 'LOOKUP_TABLE default'         
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) Magnv(nv,k)
            enddo
         enddo
         endif
         !--------
         if (Write_Magni_Vort.eq.1) then         
         write(irec,*) ' '
         write(irec,*) 'SCALARS Magni_Vort   float'
         write(irec,*) 'LOOKUP_TABLE default'         
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) Vortv(nv,k)
            enddo
         enddo
         endif                 
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 
         
         if (Write_vel.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'VECTORS vel float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufv(nv,k),vfv(nv,k),wfv(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_velEx.eq.1) then        
         write(irec,*) ' '
         write(irec,*) 'VECTORS vel_Ex float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) usv(nv,k),vsv(nv,k),wsv(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_velError.eq.1) then 
         write(irec,*) ' '
         write(irec,*) 'VECTORS vel_Error float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) uErrv(nv,k),vErrv(nv,k),wErrv(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_vorticity.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'VECTORS vorticity float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) Vortxv(nv,k),Vortyv(nv,k),Vortzv(nv,k)        
            enddo
         enddo
         endif          
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

        allocate(etav_gl(N_VERTglobal,NZglobal-1), &
                 Hprv_gl(N_VERTglobal,NZglobal-1), &
                 etavA_gl(N_VERTglobal,NZglobal-1),&
                 HprvA_gl(N_VERTglobal,NZglobal-1),&
                 xvt_gl(N_VERTglobal,NZglobal-1),  &
                 yvt_gl(N_VERTglobal,NZglobal-1),  &
                 zvt_gl(N_VERTglobal,NZglobal-1),  &
                 ufv_gl(N_VERTglobal,NZglobal-1),  &
                 vfv_gl(N_VERTglobal,NZglobal-1),  &
                 wfv_gl(N_VERTglobal,NZglobal-1),  &
                 pfv_gl(N_VERTglobal,NZglobal-1),  &
                 usv_gl(N_VERTglobal,NZglobal-1),  &
                 vsv_gl(N_VERTglobal,NZglobal-1),  &
                 wsv_gl(N_VERTglobal,NZglobal-1),  &
                 psv_gl(N_VERTglobal,NZglobal-1),  &
                 uErrv_gl(N_VERTglobal,NZglobal-1),&
                 vErrv_gl(N_VERTglobal,NZglobal-1),&
                 wErrv_gl(N_VERTglobal,NZglobal-1),&
                 pErrv_gl(N_VERTglobal,NZglobal-1),&
                 Vortv_gl(N_VERTglobal,NZglobal-1),&
                 Vortxv_gl(N_VERTglobal,NZglobal-1),&
                 Vortyv_gl(N_VERTglobal,NZglobal-1),&
                 Vortzv_gl(N_VERTglobal,NZglobal-1),&
                 Magnv_gl(N_VERTglobal,NZglobal-1))
!        ________________________________________________________
!       |                                                        |
!       |               Reconstruct global matrix                |
!       |________________________________________________________|

         call matgloV(xvt,xvt_gl)
         call matgloV(yvt,yvt_gl)
         call matgloV(zvt,zvt_gl)
!        ------------------------
         if (Write_eta.eq.1) then
            do k=1,NZ-1
               do nv=1,N_VERT
                  etavt(nv,k)  = etav(nv)
                  Hprvt(nv,k)  = Hprv(nv)      
               enddo
            enddo
            call matgloV(etavt,etav_gl)
            call matgloV(Hprvt,Hprv_gl)
         endif
         if (Write_p.eq.1) then
            call matgloV(pfv,pfv_gl)
         endif
         if (Write_pEx.eq.1) then
            call matgloV(psv,psv_gl)
         endif
          if (Write_pError.eq.1) then
            call matgloV(pErrv,pErrv_gl)
         endif
         if (Write_Magni_vel.eq.1) then
            call matgloV(Magnv,Magnv_gl)
         endif
         if (Write_Magni_Vort.eq.1) then
            call matgloV(Vortv,Vortv_gl)
         endif
!        ------------------------
         if (Write_vel.eq.1) then
            call matgloV(ufv,ufv_gl)
            call matgloV(vfv,vfv_gl)
            call matgloV(wfv,wfv_gl)
         endif
         if (Write_velEx.eq.1) then                
            call matgloV(usv,usv_gl)
            call matgloV(vsv,vsv_gl)
            call matgloV(wsv,wsv_gl)
         endif
         if (Write_velError.eq.1) then
            call matgloV(uErrv,uErrv_gl)
            call matgloV(vErrv,vErrv_gl)
            call matgloV(wErrv,wErrv_gl)
         endif
         if (Write_vorticity.eq.1) then
            call matgloV(Vortxv,Vortxv_gl)
            call matgloV(Vortyv,Vortyv_gl)
            call matgloV(Vortzv,Vortzv_gl)
         endif    
      
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
!       |            POINT_DATA            |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'POINT_DATA ',TotalN_VERT
!        __________________________________
!       |                                  |
!       |           Write scalars          |
!       |__________________________________| 

         if (Write_eta.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS eta    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal        
               write(irec,1) etav_gl(nv,k) 
            enddo
         enddo
         endif
         !--------
         if (Write_p.eq.1) then        
         write(irec,*) ' '
         write(irec,*) 'SCALARS p     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) pfv_gl(nv,k)      
            enddo
         enddo
         endif
         !--------
         if (Write_pEx.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS p_Ex     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) psv_gl(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_pError.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS p_Error     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) pErrv_gl(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_Magni_vel.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'SCALARS Magni_vel     float'
         write(irec,*) 'LOOKUP_TABLE default'         
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) Magnv_gl(nv,k)
            enddo
         enddo
         endif
         !--------
         if (Write_Magni_Vort.eq.1) then          
         write(irec,*) ' '
         write(irec,*) 'SCALARS Magni_Vort    float'
         write(irec,*) 'LOOKUP_TABLE default'         
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) Vortv_gl(nv,k)
            enddo
         enddo
         endif          
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         if (Write_vel.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'VECTORS vel float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufv_gl(nv,k),vfv_gl(nv,k),wfv_gl(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_velEx.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'VECTORS vel_Ex float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) usv_gl(nv,k),vsv_gl(nv,k),wsv_gl(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_velError.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'VECTORS vel_Error float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) uErrv_gl(nv,k),vErrv_gl(nv,k),wErrv_gl(nv,k)        
            enddo
         enddo
         endif
         !--------
         if (Write_vorticity.eq.1) then
         write(irec,*) ' '
         write(irec,*) 'VECTORS vorticity float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) Vortxv_gl(nv,k),Vortyv_gl(nv,k),Vortzv_gl(nv,k)        
            enddo
         enddo
         endif         
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
!       |            POINT_DATA            |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'POINT_DATA ',TotalN_VERT
!        __________________________________
!       |                                  |
!       |           Write scalars          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'SCALARS eta    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZ-1
            do nv=1,N_VERT        
               write(irec,1) etav(nv) 
            enddo
         enddo
         write(irec,*) ' '         
         write(irec,*) 'SCALARS p     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,1) pfv(nv,k)      
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'VECTORS vel float'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,2) ufv(nv,k),vfv(nv,k),wfv(nv,k)        
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

        deallocate(etav_gl,Hprv_gl,                 &
                   etavA_gl,HprvA_gl,               &
                   xvt_gl,yvt_gl,zvt_gl,            &
                   ufv_gl,vfv_gl,wfv_gl,pfv_gl,     &
                   usv_gl,vsv_gl,wsv_gl,psv_gl,     &
                   uErrv_gl,vErrv_gl,wErrv_gl,      &
                   pErrv_gl,                        &
                   Magnv_gl,                        &
                   Vortv_gl,                        &
                   Vortxv_gl,Vortyv_gl,Vortzv_gl)

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

!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: SaveParaviewVertex'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END OF outsavtec                            !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

