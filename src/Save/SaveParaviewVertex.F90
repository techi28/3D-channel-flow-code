!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  WRITE DATA TO DISPLAY AT PARAVIEW                  !
!                             Dic 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SaveParaviewVertex(alphafv,ufv,vfv,wfv,pfv,         &
                                    alphasv,usv,vsv,wsv,psv,rhosv,   &
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
!  | --> alphafv|(N_VERT,NZ-1)| Fluid control volume               |  !
!  | --> ufv    |(N_VERT,NZ-1)| Velocity component u_f             |  !
!  | --> vfv    |(N_VERT,NZ-1)| Velocity component v_f             |  !      
!  | --> wfv    |(N_VERT,NZ-1)| Velocity component w_f             |  !
!  | --> pfv    |(N_VERT,NZ-1)| Pressure of the fluid              |  !
!  |____________|_____________|____________________________________|  !
!  | --> alphasv|(N_VERT,NZ-1)| Solid control volume               |  !
!  | --> usv    |(N_VERT,NZ-1)| Velocity component u_s             |  !
!  | --> vsv    |(N_VERT,NZ-1)| Velocity component v_s             |  !     
!  | --> wsv    |(N_VERT,NZ-1)| Velocity component w_s             |  !
!  | --> psv    |(N_VERT,NZ-1)| Pressure of the solid              |  !
!  | --> rhosv  |(N_VERT,NZ-1)| Density of the solid               |  !
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

      real*8,dimension(:,:) :: alphafv
      real*8,dimension(:,:) :: ufv
      real*8,dimension(:,:) :: vfv
      real*8,dimension(:,:) :: wfv
      real*8,dimension(:,:) :: pfv  
      real*8,dimension(:,:) :: alphasv
      real*8,dimension(:,:) :: usv
      real*8,dimension(:,:) :: vsv
      real*8,dimension(:,:) :: wsv
      real*8,dimension(:,:) :: psv
      real*8,dimension(:,:) :: rhosv
      real*8,dimension(:,:) :: xvt
      real*8,dimension(:,:) :: yvt
      real*8,dimension(:,:) :: zvt
      integer,dimension(:,:):: No_vp
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: xvt_gl,yvt_gl,zvt_gl
      real*8,dimension(:,:),allocatable :: ufv_gl,vfv_gl,wfv_gl,pfv_gl
      real*8,dimension(:,:),allocatable :: usv_gl,vsv_gl,wsv_gl,psv_gl
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
!       |                     Open file irec                     |
!       |________________________________________________________|
 
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
         write(irec,*) 'SCALARS pf     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               !write(irec,1) pfv(nv,k)      
               write(irec,1) ufv(nv,k) 
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS ps     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               !write(irec,1) psv(nv,k)        
               write(irec,1) usv(nv,k) 
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_f float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufv(nv,k),vfv(nv,k),wfv(nv,k)        
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_s float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) usv(nv,k) ,vsv(nv,k) ,wsv(nv,k)        
            enddo
         enddo
!        ________________________________________________________
!       |                                                        |
!       |                     Close file irec                    |
!       |________________________________________________________|     

         rewind(irec)
         close(irec)
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
                 ufv_gl(N_VERTglobal,NZglobal-1), &
                 vfv_gl(N_VERTglobal,NZglobal-1), &
                 wfv_gl(N_VERTglobal,NZglobal-1), &
                 pfv_gl(N_VERTglobal,NZglobal-1), &
                 usv_gl(N_VERTglobal,NZglobal-1), &
                 vsv_gl(N_VERTglobal,NZglobal-1), &
                 wsv_gl(N_VERTglobal,NZglobal-1), &
                 psv_gl(N_VERTglobal,NZglobal-1))
!        ________________________________________________________
!       |                                                        |
!       |               Reconstruct global matrix                |
!       |________________________________________________________|  

         call matgloV(xvt,xvt_gl)
         call matgloV(yvt,yvt_gl)
         call matgloV(zvt,zvt_gl)
         call matgloV(ufv,ufv_gl)
         call matgloV(vfv,vfv_gl)
         call matgloV(wfv,wfv_gl)
         call matgloV(pfv,pfv_gl)
         call matgloV(usv,usv_gl)
         call matgloV(vsv,vsv_gl)
         call matgloV(wsv,wsv_gl)
!        ________________________________________________________
!       |                                                        |
!       |     *********  Write data file  (global) ************  |
!       |________________________________________________________|

         IF (rang_topo.eq.0) THEN
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
         write(irec,*) 'SCALARS pf     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) ufv_gl(nv,k)      
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS ps     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,1) usv_gl(nv,k)        
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |           Write vectors          |
!       |__________________________________| 

         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_f float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) ufv_gl(nv,k),vfv_gl(nv,k),wfv_gl(nv,k)        
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_s float'
         do k=1,NZglobal-1
            do nv=1,N_VERTglobal
               write(irec,2) usv_gl(nv,k) ,vsv_gl(nv,k) ,wsv_gl(nv,k)        
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

         PrintThis = 1
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
               write(irec,1) ufv(nv,k)      
            enddo
         enddo
         !--------
         write(irec,*) ' '
         write(irec,*) 'SCALARS ps     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do k=1,NZ-1
            do nv=1,N_VERT
               write(irec,1) usv(nv,k)        
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
               write(irec,2) usv(nv,k) ,vsv(nv,k) ,wsv(nv,k)        
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
                   ufv_gl,vfv_gl,wfv_gl,pfv_gl,&
                   usv_gl,vsv_gl,wsv_gl,psv_gl)

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

