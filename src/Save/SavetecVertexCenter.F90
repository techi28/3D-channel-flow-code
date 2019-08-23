!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  WRITE DATA TO DISPLAY AT TECPLOT                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SavetecVC(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                           alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                           xct,yct,zct,No_cp,                 &
                           alphafv,ufv,vfv,wfv,pfv,           &
                           alphasv,usv,vsv,wsv,psv,rhosv,     &
                           xvt,yvt,zvt,No_vp)

!---------------------------------------------------------------------!
!                                                                     !
!    This program writes the main variables to display in the         !
!    program tecplot.                                                 !
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
      real*8,dimension(:,:) :: xct
      real*8,dimension(:,:) :: yct
      real*8,dimension(:,:) :: zct
      integer,dimension(:,:):: No_cp

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
      real*8,dimension(:,:),allocatable :: xct_gl
      real*8,dimension(:,:),allocatable :: yct_gl
      real*8,dimension(:,:),allocatable :: zct_gl
      real*8,dimension(:,:),allocatable :: uf_gl
      real*8,dimension(:,:),allocatable :: vf_gl
      real*8,dimension(:,:),allocatable :: wf_gl
      real*8,dimension(:,:),allocatable :: pf_gl
      real*8,dimension(:,:),allocatable :: xvt_gl
      real*8,dimension(:,:),allocatable :: yvt_gl
      real*8,dimension(:,:),allocatable :: zvt_gl
      real*8,dimension(:,:),allocatable :: ufv_gl
      real*8,dimension(:,:),allocatable :: vfv_gl
      real*8,dimension(:,:),allocatable :: wfv_gl
      real*8,dimension(:,:),allocatable :: pfv_gl
      integer:: irec
      integer:: nv1B,nv2B,nv3B,nv1T,nv2T,nv3T
      character*50 filen
!     ----------------------------------------
      integer :: TotalN_VERT
      integer :: TotalN_ELEM
!     ----------------------------------------
#     ifndef KeyParallel
      TotalN_VERT = N_VERT*(NZ-1)+N_CELL0*(NZ-2)
      TotalN_ELEM = 5*N_CELL0*(NZ-2)
#     else
      TotalN_VERT = N_VERTglobal*(NZglobal-1)+N_CELL0global*(NZglobal-2)
      TotalN_ELEM = 5*N_CELL0global*(NZglobal-2)
#     endif

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: SavetecVC'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                        Formats                         |
!     |________________________________________________________|
4     format(A70)
5     format(8(1x,i8))
6     format(8(1x,e12.5))

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
         filen='../output/Serial/VC-     .tec'
         write(filen(21:25),'(i5.5)') SaveCounter
         open(irec,file=filen)
         write(irec,4)'TITLE     = "nsmp 3D prism data"'
         write(irec,4)'VARIABLES = "x","y","z"'
         write(irec,4)'"uf","vf","wf","pf"'!,"us","vs","ws","ps"'
         write(irec,'(a11,i7,a8,i7,a40)')'ZONE N = ',TotalN_VERT,&
                ', E = ', TotalN_ELEM,&
                ', DATAPACKING= BLOCK, ZONETYPE=FEBRICK'
         write(irec,'(a11,i3,a17,f12.5)') 'StrandID = ',1, &
                     ', SolutionTime = ', time
!        __________________________________
!       |                                  |
!       |         Write vertices           |
!       |__________________________________|

         write(irec,6) (xvt(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (xct(1:N_CELL0,k),k=2,NZ-1)
!        -------------
         write(irec,6) (yvt(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (yct(1:N_CELL0,k),k=2,NZ-1)
!        -------------
         write(irec,6) (zvt(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (zct(1:N_CELL0,k),k=2,NZ-1)
!        __________________________________
!       |                                  |
!       |       Write variable values      |
!       |__________________________________|

         write(irec,6) (ufv(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (ufnp(1:N_CELL0,k),k=2,NZ-1)
!        -------------
         write(irec,6) (vfv(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (vfnp(1:N_CELL0,k),k=2,NZ-1)
!        -------------
         write(irec,6) (wfv(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (wfnp(1:N_CELL0,k),k=2,NZ-1)
!        -------------
         write(irec,6) (pfv(1:N_VERT,k),k=1,NZ-1)
         write(irec,6) (pfnp(1:N_CELL0,k),k=2,NZ-1)
!        __________________________________
!       |                                  |
!       |  Write interconexions of elements|
!       |__________________________________|

         do k=1,NZ-2
            do i=1,N_CELL0
!              -----------------------------------------------------
!              Bottom
               nv1B = No_vp(i,1) + N_VERT*(k-1)
               nv2B = No_vp(i,2) + N_VERT*(k-1)
               nv3B = No_vp(i,3) + N_VERT*(k-1)
               nv1T = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2T = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv3T = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Top
               nv1B = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2B = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv3B = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv1T = No_vp(i,1) + N_VERT*(k)
               nv2T = No_vp(i,2) + N_VERT*(k)
               nv3T = No_vp(i,3) + N_VERT*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Side 1
               nv1B = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2B = No_vp(i,1) + N_VERT*(k-1)
               nv3B = No_vp(i,2) + N_VERT*(k-1)
               nv1T = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2T = No_vp(i,1) + N_VERT*(k)
               nv3T = No_vp(i,2) + N_VERT*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Side 2
               nv1B = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2B = No_vp(i,2) + N_VERT*(k-1)
               nv3B = No_vp(i,3) + N_VERT*(k-1)
               nv1T = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2T = No_vp(i,2) + N_VERT*(k)
               nv3T = No_vp(i,3) + N_VERT*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Side 3
               nv1B = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2B = No_vp(i,3) + N_VERT*(k-1)
               nv3B = No_vp(i,1) + N_VERT*(k-1)
               nv1T = i + N_CELL0*(k-1) + N_VERT*(NZ-1)
               nv2T = No_vp(i,3) + N_VERT*(k)
               nv3T = No_vp(i,1) + N_VERT*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
            enddo
         enddo

!        ________________________________________________________
!       |                                                        |
!       |                     Close file irec                    |
!       |________________________________________________________|

         rewind(irec)
         close(irec)

         write(*,'(t20,a20)')      ' __________________ '
         write(*,'(t20,a20)')      '|   Save Results   |'
         write(*,'(t20,a20)')      '|__________________|'
         write(*,'(t20,a7,i3)')    '  No.  :',SaveCounter
         write(*,'(t20,a7,e10.3)') '  time :',time
         write(*,'(t20,a20)')      '|__________________|'
         print*,'  '

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
!       |                   allocate variable                    |
!       |________________________________________________________|
         allocate(xct_gl(N_CELL0global,NZglobal), &
                  yct_gl(N_CELL0global,NZglobal), &
                  zct_gl(N_CELL0global,NZglobal), &
                  uf_gl(N_CELL0global,NZglobal), &
                  vf_gl(N_CELL0global,NZglobal), &
                  wf_gl(N_CELL0global,NZglobal), &
                  pf_gl(N_CELL0global,NZglobal))

         allocate(xvt_gl(N_VERTglobal,NZglobal-1), &
                  yvt_gl(N_VERTglobal,NZglobal-1), &
                  zvt_gl(N_VERTglobal,NZglobal-1), &
                  ufv_gl(N_VERTglobal,NZglobal-1), &
                  vfv_gl(N_VERTglobal,NZglobal-1), &
                  wfv_gl(N_VERTglobal,NZglobal-1), &
                  pfv_gl(N_VERTglobal,NZglobal-1))
!        ________________________________________________________
!       |                                                        |
!       |                    collect value                       |
!       |________________________________________________________|
            call matgloC(xct,xct_gl)
            call matgloC(yct,yct_gl)
            call matgloC(zct,zct_gl)
            call matgloC(ufnp,uf_gl)
            call matgloC(vfnp,vf_gl)
            call matgloC(wfnp,wf_gl)
            call matgloC(pfnp,pf_gl)

            call matgloV(xvt,xvt_gl)
            call matgloV(yvt,yvt_gl)
            call matgloV(zvt,zvt_gl)
            call matgloV(ufv,ufv_gl)
            call matgloV(vfv,vfv_gl)
            call matgloV(wfv,wfv_gl)
            call matgloV(pfv,pfv_gl)
!        ________________________________________________________
!       |                                                        |
!       |                     Open file irec                     |
!       |________________________________________________________|
         if(rang_topo .eq. 0) then
         irec=60
         filen='../output/Parallel/VC-     .tec'
         write(filen(23:27),'(i5.5)') SaveCounter
         open(irec,file=filen)
         write(irec,4)'TITLE     = "nsmp 3D prism data"'
         write(irec,4)'VARIABLES = "x","y","z"'
         write(irec,4)'"uf","vf","wf","pf"'!,"us","vs","ws","ps"'
         write(irec,'(a11,i7,a8,i7,a40)')'ZONE N = ',TotalN_VERT,&
                ', E = ', TotalN_ELEM,&
                ', DATAPACKING= BLOCK, ZONETYPE=FEBRICK'
         write(irec,'(a11,i3,a17,f12.5)') 'StrandID = ',1, &
                     ', SolutionTime = ', time
!        __________________________________
!       |                                  |
!       |         Write vertices           |
!       |__________________________________|

         write(irec,6) (xvt_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (xct_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        -------------
         write(irec,6) (yvt_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (yct_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        -------------
         write(irec,6) (zvt_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (zct_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        __________________________________
!       |                                  |
!       |       Write variable values      |
!       |__________________________________|

         write(irec,6) (ufv_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (uf_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        -------------
         write(irec,6) (vfv_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (vf_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        -------------
         write(irec,6) (wfv_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (wf_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        -------------
         write(irec,6) (pfv_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (pf_gl(1:N_CELL0global,k),k=2,NZglobal-1)
!        __________________________________
!       |                                  |
!       |  Write interconexions of elements|
!       |__________________________________|

         do k=1,NZglobal-2
            do i=1,N_CELL0global
!              -----------------------------------------------------
!              Bottom
               nv1B = No_vp_global(i,1) + N_VERTglobal*(k-1)
               nv2B = No_vp_global(i,2) + N_VERTglobal*(k-1)
               nv3B = No_vp_global(i,3) + N_VERTglobal*(k-1)
               nv1T = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2T = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv3T = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Top
               nv1B = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2B = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv3B = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv1T = No_vp_global(i,1) + N_VERTglobal*(k)
               nv2T = No_vp_global(i,2) + N_VERTglobal*(k)
               nv3T = No_vp_global(i,3) + N_VERTglobal*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Side 1
               nv1B = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2B = No_vp_global(i,1) + N_VERTglobal*(k-1)
               nv3B = No_vp_global(i,2) + N_VERTglobal*(k-1)
               nv1T = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2T = No_vp_global(i,1) + N_VERTglobal*(k)
               nv3T = No_vp_global(i,2) + N_VERTglobal*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Side 2
               nv1B = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2B = No_vp_global(i,2) + N_VERTglobal*(k-1)
               nv3B = No_vp_global(i,3) + N_VERTglobal*(k-1)
               nv1T = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2T = No_vp_global(i,2) + N_VERTglobal*(k)
               nv3T = No_vp_global(i,3) + N_VERTglobal*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
!              -----------------------------------------------------
!              Side 3
               nv1B = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2B = No_vp_global(i,3) + N_VERTglobal*(k-1)
               nv3B = No_vp_global(i,1) + N_VERTglobal*(k-1)
               nv1T = i + N_CELL0global*(k-1) + N_VERTglobal*(NZ-1)
               nv2T = No_vp_global(i,3) + N_VERTglobal*(k)
               nv3T = No_vp_global(i,1) + N_VERTglobal*(k)
               write(irec,5) nv1B,nv1B,nv1T,nv1T,nv2B,nv3B,nv3T,nv2T
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
!       |                      Display                           |
!       |________________________________________________________|

         write(*,'(t20,a20)')      ' __________________ '
         write(*,'(t20,a20)')      '|   Save Results   |'
         write(*,'(t20,a20)')      '|__________________|'
         write(*,'(t20,a7,i3)')    '  No.  :',SaveCounter
         write(*,'(t20,a7,e10.3)') '  time :',time
         write(*,'(t20,a20)')      '|__________________|'
         print*,'  '

         endif
 !        ________________________________________________________
!       |                                                        |
!       |                 deallocate variable                    |
!       |________________________________________________________|
        call MPI_Barrier(comm3D,code)

        deallocate(xct_gl,yct_gl,zct_gl, &
                   uf_gl,vf_gl,wf_gl,pf_gl)

        deallocate(xvt_gl,yvt_gl,zvt_gl, &
                   ufv_gl,vfv_gl,wfv_gl,pfv_gl)

#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: SavetecVC'
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

