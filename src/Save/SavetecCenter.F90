!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  WRITE DATA TO DISPLAY AT TECPLOT                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SavetecCenter(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xvt,yvt,zvt,No_vp)

!---------------------------------------------------------------------!
!                                                                     !
!    This program writes the main variables to display in the         !
!    program tecplot using the cell centers. In the case, the         !
!    flow variables are saved at the cell center. Six vertices        !
!    are used to provide the physical coordinate.                     !
!    Modified by Xin on 28/09/16. Old routine deleted.                !
!---------------------------------------------------------------------!
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |     Size   |        Description                  |  !
!  |____________|____________|_____________________________________|  !
!  | --> alphaf |(N_CELL0,NZ)| Fluid control volume                |  !
!  | --> uf     |(N_CELL0,NZ)| Velocity component u_f              |  !
!  | --> vf     |(N_CELL0,NZ)| Velocity component v_f              |  !
!  | --> wf     |(N_CELL0,NZ)| Velocity component w_f              |  !
!  | --> pf     |(N_CELL0,NZ)| Pressure of the fluid               |  !
!  |____________|____________|_____________________________________|  !
!  | --> alphas |(N_CELL0,NZ)| Solid control volume                |  !
!  | --> us     |(N_CELL0,NZ)| Velocity component u_s              |  !
!  | --> vs     |(N_CELL0,NZ)| Velocity component v_s              |  !
!  | --> ws     |(N_CELL0,NZ)| Velocity component w_s              |  !
!  | --> ps     |(N_CELL0,NZ)| Pressure of the solid               |  !
!  | --> rhos   |(N_CELL0,NZ)| Density of the solid                |  !
!  |____________|____________|_____________________________________|  !
!  | --> xct    |(N_CELL0,NZ)| xc at the current time              |  !
!  | --> yct    |(N_CELL0,NZ)| yc at the current time              |  !
!  | --> zct    |(N_CELL0,NZ)| yc at the current time              |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !
!  |_____________|_________________________________________________|  !
!  | --- N_CELL  | Total number of the cells                       |  !
!  | --- N_CELL0 | Inside number of cells                          |  !
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
!  | DefStructure    |  Tag of structured triangles                |  !
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

!      ________________________________________________________
!     |                                                        |
!     |   Keys, subroutines and common parameters              |
!     |________________________________________________________|

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

!      ________________________________________________________
!     |                                                        |
!     |           Definition of variables                      |
!     |________________________________________________________|

      real*8,dimension(:,:) :: alphafnp
      real*8,dimension(:,:) :: ufnp
      real*8,dimension(:,:) :: vfnp
      real*8,dimension(:,:) :: wfnp
      real*8,dimension(:,:) :: pfnp
      real*8,dimension(:,:) :: alphasnp
      real*8,dimension(:,:) :: usnp
      real*8,dimension(:,:) :: vsnp
      real*8,dimension(:,:) :: wsnp
      real*8,dimension(:,:) :: psnp
      real*8,dimension(:,:) :: rhos
      real*8,dimension(:,:) :: xvt
      real*8,dimension(:,:) :: yvt
      real*8,dimension(:,:) :: zvt
      integer,dimension(:,:):: No_vp
!      ________________________________________________________
!     |                                                        |
!     |           Definition of local variables                |
!     |________________________________________________________|
      real*8,dimension(:,:),allocatable :: xvt_gl,yvt_gl,zvt_gl
      real*8,dimension(:,:),allocatable :: uf_gl,vf_gl,wf_gl,pf_gl
      integer:: irec
      integer:: nv1B,nv2B,nv3B,nv4B,nv1T,nv2T,nv3T,nv4T
      character*50 filen
      integer :: TotalN_VERT
      integer :: TotalN_ELEM
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '----> Begin subroutine: SavetecCenter'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                ====================================                 !
!                ==========  SEQUENTIAL =============                 !
!                                                                     !
!*********************************************************************!

4     format(A80)
5     format(8(1x,i8))
6     format(8(1x,e12.5))

!         ________________________________________________________
!        |                                                        |
!        |                 Vertices and elements                  |
!        |________________________________________________________|
#     ifndef KeyParallel
         TotalN_VERT = N_VERT*(NZ-1)
         TotalN_ELEM = N_CELL0*(NZ-2)
#     else
         TotalN_VERT = N_VERTglobal*(NZglobal-1)
         TotalN_ELEM = N_CELL0global*(NZglobal-2)
#     endif

!         ________________________________________________________
!        |                                                        |
!        |                    Write tecplot file                  |
!        |________________________________________________________|

!         __________________________________
!        |                                  |
!        |          Open file irec          |
!        |__________________________________|

         irec=60
#        ifndef KeyParallel
         filen='../output/Serial/C-     .tec'
         write(filen(20:24),'(i5.5)') SaveCounter
         open(irec,file=filen)
         write(irec,4)'TITLE     = "nsmp 3D prism data"'
         write(irec,4)'VARIABLES = "xc","yc","zc","uf","vf","wf","pf"'
         write(irec,'(a11,i8,a8,i8,a40)')'ZONE N = ',TotalN_VERT,&
                ', E = ', TotalN_ELEM,&
                ', DATAPACKING= BLOCK, ZONETYPE=FEBRICK'
         write(irec,'(a40)') 'VARLOCATION=([4,5,6,7]=CELLCENTERED)'
         write(irec,'(a11,i3,a17,f12.5)') 'StrandID = ',1, &
                     ', SolutionTime = ', time
!         __________________________________
!        |                                  |
!        |         Write vertices           |
!        |__________________________________|

         write(irec,6) ((xvt(i,k),i=1,N_VERT),k=1,NZ-1)
         write(irec,6) ((yvt(i,k),i=1,N_VERT),k=1,NZ-1)
         write(irec,6) ((zvt(i,k),i=1,N_VERT),k=1,NZ-1)
!         __________________________________
!        |                                  |
!        |       Write variable values      |
!        |__________________________________|

         write(irec,6) ((ufnp(i,k),i=1,N_CELL0),k=2,NZ-1)
         write(irec,6) ((vfnp(i,k),i=1,N_CELL0),k=2,NZ-1)
         write(irec,6) ((wfnp(i,k),i=1,N_CELL0),k=2,NZ-1)
         write(irec,6) ((pfnp(i,k),i=1,N_CELL0),k=2,NZ-1)
!         __________________________________
!        |                                  |
!        |  Write interconexions of elements|
!        |__________________________________|
         do k=2,NZ-1
            do i=1,N_CELL0
               nv1B = No_vp(i,1) + N_VERT*(k-2)
               nv2B = No_vp(i,2) + N_VERT*(k-2)
               nv3B = No_vp(i,3) + N_VERT*(k-2)
               nv1T = No_vp(i,1) + N_VERT*(k-1)
               nv2T = No_vp(i,2) + N_VERT*(k-1)
               nv3T = No_vp(i,3) + N_VERT*(k-1)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
            enddo
         enddo
!         __________________________________
!        |                                  |
!        |            Close file            |
!        |__________________________________|

         rewind(irec)
         close(irec)
!         __________________________________
!        |                                  |
!        |            Display file          |
!        |__________________________________|
         write(*,'(t20,a20)')      ' __________________ '
         write(*,'(t20,a20)')      '|   Save Results   |'
         write(*,'(t20,a20)')      '|__________________|'
         write(*,'(t20,a7,i3)')    '  No.  :',SaveCounter
         write(*,'(t20,a7,e10.3)') '  time :',time
         write(*,'(t20,a20)')      '|__________________|'
         print*,'  '
!     =============== END ================
!     ====================================
#        else
         allocate(xvt_gl(N_VERTglobal,NZglobal-1), &
                  yvt_gl(N_VERTglobal,NZglobal-1), &
                  zvt_gl(N_VERTglobal,NZglobal-1), &
                  uf_gl(N_CELL0global,NZglobal), &
                  vf_gl(N_CELL0global,NZglobal), &
                  wf_gl(N_CELL0global,NZglobal), &
                  pf_gl(N_CELL0global,NZglobal))

         call matgloV(xvt,xvt_gl)
         call matgloV(yvt,yvt_gl)
         call matgloV(zvt,zvt_gl)
         call matgloC(ufnp,uf_gl)
         call matgloC(vfnp,vf_gl)
         call matgloC(wfnp,wf_gl)
         call matgloC(pfnp,pf_gl)

     IF(rang_topo .eq. 0) THEN
!         __________________________________
!        |                                  |
!        |            Open file             |
!        |__________________________________|
        filen='../output/Parallel/C-     .tec'
        write(filen(22:26),'(i5.5)') SaveCounter
        open(irec,file=filen)
        write(irec,4)'TITLE     = "nsmp 3D prism data"'
        write(irec,4)'VARIABLES = "xc","yc","zc","uf","vf","wf","pf"'
        write(irec,'(a11,i8,a8,i8,a40)')'ZONE N = ',TotalN_VERT,&
            ', E = ', TotalN_ELEM,&
            ', DATAPACKING= BLOCK, ZONETYPE=FEBRICK'
        write(irec,'(a40)') 'VARLOCATION=([4,5,6,7]=CELLCENTERED)'
        write(irec,'(a11,i3,a17,f12.5)') 'StrandID = ',1, &
            ', SolutionTime = ', time
!         __________________________________
!        |                                  |
!        |         Write vertices           |
!        |__________________________________|

         write(irec,6) ((xvt_gl(i,k),i=1,N_VERTglobal),k=1,NZglobal-1)
         write(irec,6) ((yvt_gl(i,k),i=1,N_VERTglobal),k=1,NZglobal-1)
         write(irec,6) ((zvt_gl(i,k),i=1,N_VERTglobal),k=1,NZglobal-1)
!         __________________________________
!        |                                  |
!        |       Write variable values      |
!        |__________________________________|

         write(irec,6) ((uf_gl(i,k),i=1,N_CELL0global),k=2,NZglobal-1)
         write(irec,6) ((vf_gl(i,k),i=1,N_CELL0global),k=2,NZglobal-1)
         write(irec,6) ((wf_gl(i,k),i=1,N_CELL0global),k=2,NZglobal-1)
         write(irec,6) ((pf_gl(i,k),i=1,N_CELL0global),k=2,NZglobal-1)
!         __________________________________
!        |                                  |
!        |  Write interconexions of elements|
!        |__________________________________|
         do k=2,NZglobal-1
            do i=1,N_CELL0global
               nv1B = No_vp_global(i,1) + N_VERTglobal*(k-2)
               nv2B = No_vp_global(i,2) + N_VERTglobal*(k-2)
               nv3B = No_vp_global(i,3) + N_VERTglobal*(k-2)
               nv1T = No_vp_global(i,1) + N_VERTglobal*(k-1)
               nv2T = No_vp_global(i,2) + N_VERTglobal*(k-1)
               nv3T = No_vp_global(i,3) + N_VERTglobal*(k-1)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
            enddo
         enddo
!         __________________________________
!        |                                  |
!        |            Close file            |
!        |__________________________________|

        rewind(irec)
        close(irec)
!         __________________________________
!        |                                  |
!        |            Display file          |
!        |__________________________________|
         write(*,'(t20,a20)')      ' __________________ '
         write(*,'(t20,a20)')      '|   Save Results   |'
         write(*,'(t20,a20)')      '|__________________|'
         write(*,'(t20,a7,i3)')    '  No.  :',SaveCounter
         write(*,'(t20,a7,e10.3)') '  time :',time
         write(*,'(t20,a20)')      '|__________________|'
         print*,'  '

        ENDIF

         call MPI_Barrier(comm3D,code)

         deallocate(xvt_gl,yvt_gl,zvt_gl,     &
                    uf_gl,vf_gl,wf_gl,pf_gl)
#        endif

!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t8,60a)'), '<---- End   subroutine: SavetecCenter'
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

