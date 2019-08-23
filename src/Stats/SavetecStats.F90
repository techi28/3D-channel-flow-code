!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  WRITE DATA TO DISPLAY AT TECPLOT                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SavetecStats(suf,svf,swf,spf,  &
                              xvt,yvt,zvt,No_vp)
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

      real*8,dimension(:,:):: suf(N_CELL,NZ)
      real*8,dimension(:,:):: svf(N_CELL,NZ)
      real*8,dimension(:,:):: swf(N_CELL,NZ)
      real*8,dimension(:,:):: spf(N_CELL,NZ)
!     --------------------------------------
      real*8,dimension(:,:) :: xvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: yvt(N_VERT,NZ-1)
      real*8,dimension(:,:) :: zvt(N_VERT,NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
      real*8,dimension(:,:),allocatable :: xvt_gl,yvt_gl,zvt_gl
      real*8,dimension(:,:),allocatable :: uf_gl,vf_gl,wf_gl,pf_gl
!     ----------------------------------------
      integer:: irec
      integer:: nv1B,nv2B,nv3B,nv1T,nv2T,nv3T
      character*50 filen
!     ----------------------------------------
      integer :: TotalN_VERT
      integer :: TotalN_ELEM
!     ----------------------------------------
#     ifndef KeyParallel
      TotalN_VERT = N_VERT*(NZ-1)
      TotalN_ELEM = N_CELL0*(NZ-2)
#     else
      TotalN_VERT = N_VERTglobal*(NZglobal-1)
      TotalN_ELEM = N_CELL0global*(NZglobal-2)
#     endif

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: SavetecStats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                        Formats                         |
!     |________________________________________________________|

4     format(A75)
5     format(8(1x,i8))
6     format(8(1x,e12.5))
7     format(a11,i7,a8,i7,a40)
!*********************************************************************!
!                                                                     !
!                ====================================                 !
!                ==========  SEQUENTIAL =============                 !
!                                                                     !
!*********************************************************************!

#     ifndef KeyParallel
!     ----------------------------------------------------------
!     Allocate vertex values
      allocate(uf_gl(N_VERT,NZ), &
               vf_gl(N_VERT,NZ), &
               wf_gl(N_VERT,NZ), &
               pf_gl(N_VERT,NZ))
!   ------------------------------------------------------------
!    calculate the mean value
         do k=1,NZ
            do i=1,N_CELL
              uf_gl(i,k) = suf(i,k)/scount
              vf_gl(i,k) = svf(i,k)/scount
              wf_gl(i,k) = swf(i,k)/scount
              pf_gl(i,k) = spf(i,k)/scount
            enddo
         enddo
!        ________________________________________________________
!       |                                                        |
!       |                     Open file irec                     |
!       |________________________________________________________|

         irec=60
         filen='../output/Stats/SV-     .tec'
         write(filen(20:24),'(i5.5)') SaveCounter
         open(irec,file=filen)
         write(irec,4)'TITLE     = "nsmp stats data"'
         write(irec,4)'VARIABLES = "xv","yv","zv","us","vs","ws","ps"'
         write(irec,'(a11,i7,a8,i7,a40)')'ZONE N = ',TotalN_VERT,&
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

         write(irec,6) ((uf_gl(i,k),i=1,N_CELL0),k=2,NZ-1)
         write(irec,6) ((vf_gl(i,k),i=1,N_CELL0),k=2,NZ-1)
         write(irec,6) ((wf_gl(i,k),i=1,N_CELL0),k=2,NZ-1)
         write(irec,6) ((pf_gl(i,k),i=1,N_CELL0),k=2,NZ-1)
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
!        ________________________________________________________
!       |                                                        |
!       |                     Close file irec                    |
!       |________________________________________________________|

         rewind(irec)
         close(irec)
         deallocate(uf_gl,vf_gl,wf_gl,pf_gl)
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
                 uf_gl(N_CELL0global,NZglobal), &
                 vf_gl(N_CELL0global,NZglobal), &
                 wf_gl(N_CELL0global,NZglobal), &
                 pf_gl(N_CELL0global,NZglobal))
!        ________________________________________________________
!       |                                                        |
!       |               Reconstruct global matrix                |
!       |________________________________________________________|

         call matgloV(xvt,xvt_gl)
         call matgloV(yvt,yvt_gl)
         call matgloV(zvt,zvt_gl)
         call matgloC(suf,uf_gl)
         call matgloC(svf,vf_gl)
         call matgloC(swf,wf_gl)
         call matgloC(spf,pf_gl)

         IF (rang_topo.eq.0) THEN
!        _________________________________
!       |                                 |
!       |      calculate averaged data    |
!       |_________________________________|
         do k=1,NZglobal
            do i=1,N_CELL0global
              uf_gl(i,k) = uf_gl(i,k)/scount
              vf_gl(i,k) = vf_gl(i,k)/scount
              wf_gl(i,k) = wf_gl(i,k)/scount
              pf_gl(i,k) = pf_gl(i,k)/scount
            enddo
         enddo
!        __________________________________
!       |                                  |
!       |          Open file irec          |
!       |__________________________________|

         irec=60
         filen='../output/Stats/SV-     .tec'
         write(filen(20:24),'(i5.5)') SaveCounter
         open(irec,file=filen)
         write(irec,4)'TITLE     = "nsmp 3D stats data"'
         write(irec,4)'VARIABLES = "x","y","z","us","vs","ws","ps"'
         write(irec,'(a11,i7,a8,i7,a40)')'ZONE N = ',TotalN_VERT,&
                ', E = ', TotalN_ELEM,&
                ', DATAPACKING= BLOCK, ZONETYPE=FEBRICK'
         write(irec,'(a40)') 'VARLOCATION=([4,5,6,7]=CELLCENTERED)'
         write(irec,'(a11,i3,a17,f12.5)') 'StrandID = ',1, &
                     ', SolutionTime = ', time
!        __________________________________
!       |                                  |
!       |         Write vertices           |
!       |__________________________________|

         write(irec,6) (xvt_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (yvt_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
         write(irec,6) (zvt_gl(1:N_VERTglobal,k),k=1,NZglobal-1)
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
!        __________________________________
!       |                                  |
!       |          Close file irec         |
!       |__________________________________|

         rewind(irec)
         close(irec)

         ENDIF

         call MPI_Barrier(comm3D,code)
!        ________________________________________________________
!       |                                                        |
!       |                Deallocate global matrix                |
!       |________________________________________________________|

        deallocate(xvt_gl,yvt_gl,zvt_gl,&
                   uf_gl,vf_gl,wf_gl,pf_gl)

#     endif

!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: SavetecStats'
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

