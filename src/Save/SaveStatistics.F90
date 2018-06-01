!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!=====================================================================!
!---------------------------------------------------------------------!
!             Turbulence Statistics (UNSTRUCTURED MESH)               !
!                     Xin Bai & Miguel Uh Zapata                      !
!                              Jan 2018                               !
!---------------------------------------------------------------------!
!=====================================================================!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!——————————————————————————————————————————————————————————————————---!
!                                                                     !
!    This subroutine calculates the turbulent statistical terms,      !
!    including time averaged velocity, velocity fluctuation,          !
!    Reynolds stress and etc. Currently only valid for rectangle      !
!    structured grid.                                                 !
!                                                                     !
!    SUBROUNTINES:                                                    !
!                                                                     !
!                    [1.1] sm_init                                    !
!                    [1.2] glob_get_mm                                !
!                    [1.3] glob_sample                                !
!                    [1.4] glob_put_mm   [Not used now]               !
!                    [1.5] glob_put_plane_mm                          !
!                    [1.6] update_stats  [Not used now]               !
!                    [1.7] restart_in_stats                           !
!                    [1.8] restart_out_stats                          !
!                    [1.9] SavetecStats                               !
!                                                                     !
!---------------------------------------------------------------------!

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.1]              Initialization Statistics                      !
!                             Jan 2018                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE sm_init

!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         USE geometry
         USE statistics
         implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         USE parallel
         USE geometry
         USE statistics
         implicit none
#     endif
!     =============== END ================
!     ====================================

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: sm_init'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

!     -------------------------------------------------
!     -----    Set the initial values to zero      ----
!     -------------------------------------------------

      sm_uf  = 0.0d0
      sm_vf  = 0.0d0
      sm_wf  = 0.0d0
      sm_pf  = 0.0d0
      sm_u2f = 0.0d0
      sm_v2f = 0.0d0
      sm_w2f = 0.0d0
      sm_uvf = 0.0d0
      sm_vwf = 0.0d0
      sm_wuf = 0.0d0
      sm_p2f = 0.0d0
      sm_ef  = 0.0d0
!     --------------
      sm_n   = 0
      
!     -----------------------------------------------
!     ------Read previously saved momentum stats-----
!     -----------------------------------------------
#     ifdef KeyReStartStatistics
         call glob_get_mm
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: sm_init'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE sm_init

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.2]      Read Global Momentum Stats (for restart)               !
!                             Jan 2018                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE glob_get_mm

!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE statistics
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE statistics
            implicit none
#     endif
!     =============== END ================
!     ====================================
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8, dimension(:),allocatable :: ufaux,vfaux,wfaux
      real*8, dimension(:),allocatable :: u2faux,v2faux,w2faux
      real*8, dimension(:),allocatable :: uvfaux,vwfaux,wufaux
      real*8, dimension(:),allocatable :: pfaux,p2faux
      real*8, dimension(:),allocatable :: efaux
      integer:: m,elem,N_total,kNZ
      integer:: ifile
      character*50 ifiledat
      character*80 title
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: glob_get_mm'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

!     -----------------------------------
!     Total number of elements
#     ifndef KeyParallel
         N_total = N_CELL0*NZ
#     else
         N_total = N_CELL0global*NZglobal
#     endif

!     -----------------------------------
!     Allocate
      allocate(ufaux(N_total),vfaux(N_total),wfaux(N_total),    &
               u2faux(N_total),v2faux(N_total),w2faux(N_total), &
               uvfaux(N_total),vwfaux(N_total),wufaux(N_total), &
               pfaux(N_total), p2faux(N_total),efaux(N_total))

!     -----------------------------------
!     Formats
79    format(18x,i8)
5     format(82a)
35    format(6e12.5)

!      ________________________________________________________
!     |                                                        |
!     |                        Read file                       |
!     |________________________________________________________|

!      __________________________________
!     |                                  |
!     |            Open file             |
!     |__________________________________|

      ifile=605
      ifiledat = '../output/Stats/glob_mm_restart.dat'
      open(ifile,file=ifiledat,status='old')
!      __________________________________
!     |                                  |
!     |              Time                |
!     |__________________________________|

      read(ifile,79) sm_n
!      __________________________________
!     |                                  |
!     |         Fluid variables          |
!     |__________________________________|

      read(ifile,5) title
      read(ifile,35)(ufaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(vfaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(wfaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(u2faux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(v2faux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(w2faux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(uvfaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(vwfaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(wufaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(pfaux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(p2faux(m),m=1,N_total)
      read(ifile,5) title
      read(ifile,35)(efaux(m),m=1,N_total)

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
        do k=1,NZ
           do i=1,N_CELL0
              m = (k-1)*N_CELL0 + i
!             ---------------------------
              sm_uf(i,k) = ufaux(m)
              sm_vf(i,k) = vfaux(m)
              sm_wf(i,k) = wfaux(m)
              sm_pf(i,k) = pfaux(m)
!             ---------------------------
              sm_u2f(i,k) = u2faux(m)
              sm_v2f(i,k) = v2faux(m)
              sm_w2f(i,k) = w2faux(m)
              sm_p2f(i,k) = p2faux(m)
!             ---------------------------
              sm_uvf(i,k) = uvfaux(m)
              sm_vwf(i,k) = vwfaux(m)
              sm_wuf(i,k) = wufaux(m)
              sm_ef(i,k) = efaux(m)
           enddo
        enddo
!     ====================================
!     ============  PARALLEL  ============
#     else
        do k=2,NZ-1
           do i=1,N_CELL0
              elem = index_global(i)
              m = (k-1)*N_CELL0global + elem
!             ---------------------------
              sm_uf(i,k) = ufaux(m)
              sm_vf(i,k) = vfaux(m)
              sm_wf(i,k) = wfaux(m)
              sm_pf(i,k) = pfaux(m)
!             ---------------------------
              sm_u2f(i,k) = u2faux(m)
              sm_v2f(i,k) = v2faux(m)
              sm_w2f(i,k) = w2faux(m)
              sm_p2f(i,k) = p2faux(m)
!             ---------------------------
              sm_uvf(i,k) = uvfaux(m)
              sm_vwf(i,k) = vwfaux(m)
              sm_wuf(i,k) = wufaux(m)
              sm_ef(i,k)  = efaux(m)
           enddo
        enddo
#     endif
!     =============== END ================
!     ====================================

!      __________________________________
!     |                                  |
!     |            Close file            |
!     |__________________________________|

      close(ifile)

!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

!     -----------------------------------
!     Display
#     ifdef KeyParallel
       if(rang_topo .eq. 0) then
           print*, 'Stats restarts successfully read in'
       endif
       call MPI_Barrier(comm3D,code)
#     else
           print*, 'Stats restarts successfully read in'
#     endif

!     -----------------------------------
!     Deallocate
      deallocate(ufaux,vfaux,wfaux,    &
                 u2faux,v2faux,w2faux, &
                 uvfaux,vwfaux,wufaux, &
                 pfaux,p2faux,efaux)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: glob_get_mm'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE glob_get_mm


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.3]               Taking Stats Samples                          !
!                             Jan 2018                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE glob_sample(phiu,phiv,phiw,phip)

!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         USE geometry
         USE statistics
#        ifdef KeyLES
           USE les
#        endif
         implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         USE parallel
         USE geometry
         USE statistics
#        ifdef KeyLES
            USE les
#        endif
         implicit none
#     endif
!     =============== END ================
!     ====================================
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|

      real*8,dimension(:,:) :: phiu(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_CELL,NZ)
      real*8,dimension(:,:) :: phiw(N_CELL,NZ)
      real*8,dimension(:,:) :: phip(N_CELL,NZ)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: tu,tv,tw,tp
      real*8 :: tu2,tv2,tw2,tp2
      real*8 :: tuv,tvw,twu

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: sample_stats '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

      do k=1,NZ
        do i=1,N_CELL0
        tu = phiu(i,k)
        tv = phiv(i,k)
        tw = phiw(i,k)
        tp = phip(i,k)
        tp2 = tp*tp
        tu2 = tu*tu
        tv2 = tv*tv
        tw2 = tw*tw
        tuv = tu*tv
        tvw = tv*tw
        twu = tw*tu
!       -------------------------------
        sm_uf(i,k)  = sm_uf(i,k) + tu
        sm_vf(i,k)  = sm_vf(i,k) + tv
        sm_wf(i,k)  = sm_wf(i,k) + tw
        sm_pf(i,k)  = sm_pf(i,k) + tp
!       -------------------------------
        sm_u2f(i,k) = sm_u2f(i,k)+ tu2
        sm_v2f(i,k) = sm_v2f(i,k)+ tv2
        sm_w2f(i,k) = sm_w2f(i,k)+ tw2
        sm_p2f(i,k) = sm_p2f(i,k)+ tp2
!       -------------------------------
        sm_uvf(i,k) = sm_uvf(i,k)+ tuv
        sm_vwf(i,k) = sm_vwf(i,k)+ tvw
        sm_wuf(i,k) = sm_wuf(i,k)+ twu
!       -------------------------------
#       ifdef KeyLES
        sm_ef(i,k) = sm_ef(i,k)+ eddy(i,k)
#       endif
        enddo
      enddo
!     -----------------------------------
!     Update the sampling counter
      sm_n = sm_n + 1

!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

      if (mod(sm_n,10) .eq. 0) then
#        ifdef KeyParallel
            if(rang_topo .eq. 0) then
            print*, '       =========================================='
            print*, '               STATS SAMPLED',sm_n,'TIMES'
            print*, '       =========================================='
            endif
#        else
            print*, '       =========================================='
            print*, '               STATS SAMPLED',sm_n,'TIMES'
            print*, '       =========================================='
#       endif
     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: sample_stats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE glob_sample

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.4]            Output Global Momentum Stats                     !
!                             Jan 2018                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE glob_put_mm

!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE statistics
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE statistics
            implicit none
#     endif
!     =============== END ================
!     ====================================
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

#     ifdef KeyParallel
      real*8,dimension(:,:),allocatable :: varuf,varvf,varwf
      real*8,dimension(:,:),allocatable :: varu2f,varv2f,varw2f
      real*8,dimension(:,:),allocatable :: varuvf,varvwf,varwuf
      real*8,dimension(:,:),allocatable :: varpf,varp2f,varef
#     endif
      integer:: ifile,elem
      character*50 filedat

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: glob_put_mm'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

79    format('  SampleCounter = ',i8)
35    format(6e12.5)
5     format(t2,80a)

!      ________________________________________________________
!     |         ====================================           |
!     |         ==========  SEQUENTIAL =============           |
!     |________________________________________________________|

#     ifndef KeyParallel
!      __________________________________
!     |                                  |
!     |          Open file file          |
!     |__________________________________|

      ifile=605
      filedat = '../output/Stats/glob_mm_restart.dat'
      open(ifile,file=filedat,status='unknown')
!      __________________________________
!     |                                  |
!     |               Time               |
!     |__________________________________|

      write(ifile,79) sm_n
!      __________________________________
!     |                                  |
!     |         Stats variables          |
!     |__________________________________|

      write(ifile,5) 'uf'
      write(ifile,35)((sm_uf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'vf'
      write(ifile,35)((sm_vf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'wf'
      write(ifile,35)((sm_wf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'u2f'
      write(ifile,35)((sm_u2f(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'v2f'
      write(ifile,35)((sm_v2f(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'w2f'
      write(ifile,35)((sm_w2f(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'uvf'
      write(ifile,35)((sm_uvf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'vwf'
      write(ifile,35)((sm_vwf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'wuf'
      write(ifile,35)((sm_wuf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'pf'
      write(ifile,35)((sm_pf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'p2f'
      write(ifile,35)((sm_p2f(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'ef'
      write(ifile,35)((sm_ef(i,k),i=1,N_CELL0),k=1,NZ)
!      __________________________________
!     |                                  |
!     |            Close file            |
!     |__________________________________|

      close(ifile)
      
#     else

!      ________________________________________________________
!     |         ====================================           |
!     |         ============  PARALLEL =============           |
!     |________________________________________________________|

!      __________________________________
!     |                                  |
!     |            Allocate              |
!     |__________________________________|

      allocate(varuf(N_CELL0global,NZglobal),  &
               varvf(N_CELL0global,NZglobal),  &
               varwf(N_CELL0global,NZglobal),  &
               varu2f(N_CELL0global,NZglobal), &
               varv2f(N_CELL0global,NZglobal), &
               varw2f(N_CELL0global,NZglobal), &
               varuvf(N_CELL0global,NZglobal), &
               varvwf(N_CELL0global,NZglobal), &
               varwuf(N_CELL0global,NZglobal), &
               varpf(N_CELL0global,NZglobal),  &
               varp2f(N_CELL0global,NZglobal), &
               varef(N_CELL0global,NZglobal))

!      __________________________________
!     |                                  |
!     |   Assign local value to global   |
!     |__________________________________|

      varuf  = 0.0d0
      varvf  = 0.0d0
      varwf  = 0.0d0
      varpf  = 0.0d0
      varu2f = 0.0d0
      varv2f = 0.0d0
      varw2f = 0.0d0
      varuvf = 0.0d0
      varvwf = 0.0d0
      varwuf = 0.0d0
      varp2f = 0.0d0
      varef  = 0.0d0

      call matgloC(sm_uf,varuf)
      call matgloC(sm_vf,varvf)
      call matgloC(sm_wf,varwf)
      call matgloC(sm_u2f,varu2f)
      call matgloC(sm_v2f,varv2f)
      call matgloC(sm_w2f,varw2f)
      call matgloC(sm_uvf,varuvf)
      call matgloC(sm_vwf,varvwf)
      call matgloC(sm_wuf,varwuf)
      call matgloC(sm_pf,varpf)
      call matgloC(sm_p2f,varp2f)
      call matgloC(sm_ef,varef)

      IF(rang_topo .EQ. 0) THEN
!      __________________________________
!     |                                  |
!     |          Open file file          |
!     |__________________________________|

      ifile=605
      filedat = '../output/Stats/glob_mm_restart.dat'
      open(ifile,file=filedat,status='unknown')
!      __________________________________
!     |                                  |
!     |               Time               |
!     |__________________________________|

      write(ifile,79) sm_n
!      __________________________________
!     |                                  |
!     |         Stats variables          |
!     |__________________________________|

      write(ifile,5) 'uf'
      write(ifile,35)((varuf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'vf'
      write(ifile,35)((varvf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'wf'
      write(ifile,35)((varwf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'u2f'
      write(ifile,35)((varu2f(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'v2f'
      write(ifile,35)((varv2f(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'w2f'
      write(ifile,35)((varw2f(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'uvf'
      write(ifile,35)((varuvf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'vwf'
      write(ifile,35)((varvwf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'wuf'
      write(ifile,35)((varwuf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'pf'
      write(ifile,35)((varpf(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'p2f'
      write(ifile,35)((varp2f(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'ef'
      write(ifile,35)((varef(i,k),i=1,N_CELL0global),k=1,NZglobal)
!      __________________________________
!     |                                  |
!     |            Close file            |
!     |__________________________________|

      close(ifile)
      
!      __________________________________
!     |                                  |
!     |           Deallocate             |
!     |__________________________________|
      
      deallocate(varuf,varvf,varwf,    &
                 varu2f,varv2f,varw2f, &
                 varuvf,varvwf,varwuf, &
                 varpf,varp2f,varef)
                     
      ENDIF

#     endif
!     =============== END ================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

#     ifdef KeyParallel
!     ----------------------------
!     wait for output complete
      call MPI_Barrier(comm3D,code)
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: glob_put_mm'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE glob_put_mm

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.5]             Output Plane Averaged Stats                     !
!                             Jan 2018                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE glob_put_plane_mm

!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE statistics
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE statistics
            implicit none
#     endif
!     =============== END ================
!     ====================================
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:),allocatable :: Su1
      real*8,dimension(:),allocatable :: Sv1
      real*8,dimension(:),allocatable :: Sw1
      real*8,dimension(:),allocatable :: Su2
      real*8,dimension(:),allocatable :: Sv2
      real*8,dimension(:),allocatable :: Sw2
      real*8,dimension(:),allocatable :: Suv
      real*8,dimension(:),allocatable :: Svw
      real*8,dimension(:),allocatable :: Swu
      real*8 :: u1m,u2m,uvm,v1m,v2m,vwm,w1m,w2m,wum
      real*8 :: urms, vrms, wrms, tkesum
      real*8 :: xrms, yrms, zrms
      real*8 :: Tuvm, Tvwm, Twum, npjx,yplus,yy
      integer :: ns,irec,IDWRITE,kNZ
      integer :: WriteA,WriteB,WriteC 
      character*30 filedat0,filedat1,filedat2
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: output_plane_sm'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      WriteA = 1
      WriteB = 1
      WriteC = 1
!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

!     -----------------------------------
!     Formats
6     format(I3,4(2x,e12.5))

!     -----------------------------------
!     Allocate
      allocate (Su1(NZglobal),Su2(NZglobal),Suv(NZglobal), &
                Sv1(NZglobal),Sv2(NZglobal),Svw(NZglobal), &
                Sw1(NZglobal),Sw2(NZglobal),Swu(NZglobal))

!     -----------------------------------
!     Set the initial values to zero
      do k=1,NZglobal
        Su1(k) = 0.0d0
        Sv1(k) = 0.0d0
        Sw1(k) = 0.0d0
        Su2(k) = 0.0d0
        Sv2(k) = 0.0d0
        Sw2(k) = 0.0d0
        Suv(k) = 0.0d0
        Svw(k) = 0.0d0
        Swu(k) = 0.0d0
      enddo
      
!      ________________________________________________________
!     |                                                        |
!     | Update the averaged plane values (time & horz. space)  |
!     |________________________________________________________|

      ns = sm_n
      
      do k=2,NZ-1
         do i=1,N_CELL0
            Su1(k) = Su1(k) + sm_uf(i,k)/ns
            Sv1(k) = Sv1(k) + sm_vf(i,k)/ns
            Sw1(k) = Sw1(k) + sm_wf(i,k)/ns
            Su2(k) = Su2(k) + sm_u2f(i,k)/ns
            Sv2(k) = Sv2(k) + sm_v2f(i,k)/ns
            Sw2(k) = Sw2(k) + sm_w2f(i,k)/ns
            Suv(k) = Suv(k) + sm_uvf(i,k)/ns
            Svw(k) = Svw(k) + sm_vwf(i,k)/ns
            Swu(k) = Swu(k) + sm_wuf(i,k)/ns
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call MPI_Barrier(comm3D,code)
        call SUM_parallelr(Su1,NZglobal)
        call SUM_parallelr(Sv1,NZglobal)
        call SUM_parallelr(Sw1,NZglobal)
        call SUM_parallelr(Su2,NZglobal)
        call SUM_parallelr(Sv2,NZglobal)
        call SUM_parallelr(Sw2,NZglobal)
        call SUM_parallelr(Suv,NZglobal)
        call SUM_parallelr(Svw,NZglobal)
        call SUM_parallelr(Swu,NZglobal)
#     endif
!     =============== END ================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |          Calculate and output the Statistics           |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         npjx = 1.0d0/N_CELL0global
         if (rang_topo .eq. 0) then
            IDWRITE = 1
         else
            IDWRITE = 0
         endif
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
        npjx = 1.0d0/N_CELL0
        IDWRITE = 1
#     endif
!     =============== END ================
!     ====================================

!     ------------------------------------------------------
      IF (IDWRITE .eq. 1) THEN
!         __________________________________
!        |                                  |
!        |            Open file             |
!        |__________________________________|

         filedat0 ='../output/Stats/sm_     .dat'
         filedat1 ='../output/Stats/srms_     .dat'
         filedat2 ='../output/Stats/stau_     .dat'
         write(filedat0(20:24),'(i5.5)') ns
         write(filedat1(22:26),'(i5.5)') ns
         write(filedat2(22:26),'(i5.5)') ns
         if (WriteA.eq.1) open(604,file=filedat0,status='unknown')
         if (WriteB.eq.1) open(605,file=filedat1,status='unknown')
         if (WriteC.eq.1) open(606,file=filedat2,status='unknown')

!         __________________________________
!        |                                  |
!        |            Write data            |
!        |__________________________________|
         
         tkesum = 0.0d0
         DO k=2,NZglobal-1
!           ---------------------------------
!           Non-dimensional: y+=(Re_tau/delta)*y
            yy = (2.0d0*k-3.0d0)/(2.0d0*NZglobal-4.0d0)*H0
            yplus = Re*yy  !<--- (delta=1)
!           ----------------------
!           [A] Mean Velocity
            u1m = Su1(k)*npjx
            v1m = Sv1(k)*npjx
            w1m = Sw1(k)*npjx
!           ----------------------
!           [B] Root-mean-square velocity fluctuations
            u2m = Su2(k)*npjx
            v2m = Sv2(k)*npjx
            w2m = Sw2(k)*npjx
            xrms = u2m - u1m*u1m
            yrms = v2m - v1m*v1m
            zrms = w2m - w1m*w1m
            urms = dsqrt(xrms)
            vrms = dsqrt(yrms)
            wrms = dsqrt(zrms)
!           ----------------------
!           [C] 
            uvm  = Suv(k)*npjx
            vwm  = Svw(k)*npjx
            wum  = Swu(k)*npjx        
            Tuvm = u1m*v1m - uvm
            Tvwm = v1m*w1m - vwm
            Twum = w1m*u1m - wum
!           ----------------------
!           [D] Turbulent Kinetic Energy
            tkesum = tkesum + 0.5d0*(xrms + yrms + zrms)
!           ----------------------
!           Write files
            if (WriteA.eq.1) write(604,6) k,yplus,u1m,v1m,w1m
            if (WriteB.eq.1) write(605,6) k,yplus,urms,vrms,wrms
            if (WriteC.eq.1) write(606,6) k,yplus,Tuvm,Tvwm,Twum
         ENDDO
!        ------------------------------------------------
!        Close FILE
         !write(604,*) 'TKE =', tkesum/(NZglobal-2)
!         __________________________________
!        |                                  |
!        |            Close file            |
!        |__________________________________|

         if (WriteA.eq.1) close(604)
         if (WriteB.eq.1) close(605)
         if (WriteC.eq.1) close(606)
!         __________________________________
!        |                                  |
!        |             Display              |
!        |__________________________________|
         
         print*,'       ---> Plane Averaged Stats printed!',ns
         print*,' '
      ENDIF
      
!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call MPI_Barrier(comm3D,code)
#     endif
!     =============== END ================
!     ====================================

       deallocate(Su1,Su2,Suv, &
                  Sv1,Sv2,Svw, &
                  Sw1,Sw2,Swu)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: output_plane_sm'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE glob_put_plane_mm

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.6]         Update statictics (cylinder test case)              !
!                             March 2016                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE update_stats(suf,svf,swf,spf,   &
                              ufnp,vfnp,wfnp,pfnp)

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
      real*8,dimension(:,:):: ufnp(N_CELL,NZ)
      real*8,dimension(:,:):: vfnp(N_CELL,NZ)
      real*8,dimension(:,:):: wfnp(N_CELL,NZ)
      real*8,dimension(:,:):: pfnp(N_CELL,NZ)
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: update_stats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

!      ------------------------------------
       scount = scount + 1
!     ----------------------------------
#     ifndef KeyParallel
         print*,'       =================================== '
         print*,'              Field Info Sample Count :',scount
         print*,'       =================================== '
#     else
         if (rang_topo .eq. 0) then
         print*,'       =================================== '
         print*,'              Field Info Sample Count :',scount
         print*,'       =================================== '
         endif
#     endif
!      ------------------------------------
       do k=1,NZ
          do i=1,N_CELL
            suf(i,k) = suf(i,k) + ufnp(i,k)
            svf(i,k) = svf(i,k) + vfnp(i,k)
            swf(i,k) = swf(i,k) + wfnp(i,k)
            spf(i,k) = spf(i,k) + pfnp(i,k)
          enddo
       enddo
!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: update_stats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.7]              READ RESULTS TO RESTART                        !
!                             Mar 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE restart_in_stats(suf,svf,swf,spf)

!---------------------------------------------------------------------!
!                                                                     !
!    This program read all the variables needed to restart a          !
!    new simulation using the final values of the previous si-        !
!    mulation. We need the name of the input restart file de-         !
!    fined in the input subroutine as "filerepin".                    !
!                                                                     !
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

      real*8,dimension(:,:):: suf(N_CELL,NZ)
      real*8,dimension(:,:):: svf(N_CELL,NZ)
      real*8,dimension(:,:):: swf(N_CELL,NZ)
      real*8,dimension(:,:):: spf(N_CELL,NZ)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8, dimension(:),allocatable :: ufaux,vfaux,wfaux,pfaux
      integer :: N_total,m,kNZ,ifile,Dothis
      character*50 filedat
      character*50 title
      logical exist

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: initial_stats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

!     -----------------------------------
!     Check if the restart data exist
      DoThis = 0
      ifile=16
      filedat = '../output/Stats/TimeStats-    .dat'
      write(filedat(27:30),'(i4.4)') SaveCounter
      inquire(file=filedat,exist=exist)
      if (exist) DoThis = 1
      
!      ________________________________________________________
!     |                                                        |
!     |               Reading the re-start data                |
!     |________________________________________________________|

      IF (DoThis.eq.1) THEN

!        -----------------------------------
!        Total number of points
#        ifdef KeyParallel
            N_total = N_CELL0global*NZglobal
#        else
            N_total = N_CELL0*NZ
#        endif

!        -----------------------------------
!        Allocate
         allocate(ufaux(N_total),vfaux(N_total),wfaux(N_total), &
                  pfaux(N_total))

!        -----------------------------------
!        Formats
79       format(16x,i8)
5        format(82a)
35       format(6e12.5)

!        __________________________________
!       |                                  |
!       |      Open file: filerepin        |
!       |__________________________________|

         open(ifile,file=filedat,status='old')
!        __________________________________
!       |                                  |
!       |              Time                |
!       |__________________________________|

         read(ifile,79) scount
!        __________________________________
!       |                                  |
!       |         Fluid variables          |
!       |__________________________________|

         read(ifile,5) title
         read(ifile,35)(ufaux(m),m=1,N_total)
         read(ifile,5) title
         read(ifile,35)(vfaux(m),m=1,N_total)
         read(ifile,5) title
         read(ifile,35)(wfaux(m),m=1,N_total)
         read(ifile,5) title
         read(ifile,35)(pfaux(m),m=1,N_total)
!        __________________________________
!       |                                  |
!       |        Assign final matrices     |
!       |__________________________________|

!        ====================================
!        ===========  SEQUENTIAL  ===========
#        ifndef KeyParallel
            do k=1,NZ
               do i=1,N_CELL0
                  m = (k-1)*N_CELL0 + i
!                 ----------
                  suf(i,k) = ufaux(m)
                  svf(i,k) = vfaux(m)
                  swf(i,k) = wfaux(m)
                  spf(i,k) = pfaux(m)
               enddo
            enddo
#        else
!        ====================================
!        ===========   PARALLEL  ===========
            do k=1,NZ
               do i=1,N_CELL0
                  nc = index_global(i)
#                 ifndef KeyParallelNZ
                     m = (k-1)*N_CELL0global + nc
#                 else
                     kNZ = (NZLayer-1)*(NZ-2) + k
                     m = (kNZ-1)*N_CELL0global + nc
#                 endif
!                 ----------
                  suf(i,k) = ufaux(m)
                  svf(i,k) = vfaux(m)
                  swf(i,k) = wfaux(m)
                  spf(i,k) = pfaux(m)
               enddo
            enddo
#        endif
!        =============== END ================
!        ====================================
!        __________________________________
!       |                                  |
!       |            Close file            |
!       |__________________________________|

         close(ifile)
!        __________________________________
!       |                                  |
!       |            Deallocate            |
!       |__________________________________|

         deallocate(ufaux,vfaux,wfaux,pfaux)

     ELSE
!       -----------------------------------
!       For first time use
#       ifndef KeyParallel
           print*,'   =================================== '
           print*,'          Begin Time Stats Info'
           print*,'   =================================== '
#       else
           if (rang_topo .eq. 0) then
           print*,'   =================================== '
           print*,'          Begin Time Stats Info'
           print*,'   =================================== '
           endif
#       endif
!       -----------------------------------
!       Initialize
        suf = 0.0d0
        svf = 0.0d0
        swf = 0.0d0
        spf = 0.0d0
      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: initial_stats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.8]              WRITE RESULTS TO RESTART                       !
!                             Mar 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE restart_out_stats(suf,svf,swf,spf)

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
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

#      ifdef KeyParallel
        real*8,dimension(:,:),allocatable:: varf1
        real*8,dimension(:,:),allocatable:: varf2
        real*8,dimension(:,:),allocatable:: varf3
        real*8,dimension(:,:),allocatable:: varf4
#      endif
       integer:: ifile
       character*50 filedat

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: restart_out_stats '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

79    format('       scount = ',i8)
35    format(6e12.5)
5     format(t2,80a)

!      ________________________________________________________
!     |         ====================================           |
!     |         ==========  SEQUENTIAL =============           |
!     |________________________________________________________|

#     ifndef KeyParallel
!      __________________________________
!     |                                  |
!     |             Open file            |
!     |__________________________________|

      ifile=16
      filedat = '../output/Stats/TimeStats-    .dat'
      write(filedat(27:30),'(i4.4)') SaveCounter
      open(ifile,file=filedat,status='unknown')
!      __________________________________
!     |                                  |
!     |               Time               |
!     |__________________________________|

      write(ifile,79) scount
!      __________________________________
!     |                                  |
!     |         Fluid variables          |
!     |__________________________________|

      write(ifile,5) 'uf'
      write(ifile,35)((suf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'vf'
      write(ifile,35)((svf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'wf'
      write(ifile,35)((swf(i,k),i=1,N_CELL0),k=1,NZ)
      write(ifile,5) 'pf'
      write(ifile,35)((spf(i,k),i=1,N_CELL0),k=1,NZ)
!      __________________________________
!     |                                  |
!     |            Close file            |
!     |__________________________________|

      close(ifile)

#     endif

!      ________________________________________________________
!     |         ====================================           |
!     |         =====  START PARALLEL OPTION =======           |
!     |________________________________________________________|

#     ifdef KeyParallel
!      __________________________________
!     |                                  |
!     |         Allocate variable        |
!     |__________________________________|

       allocate(varf1(N_CELL0global,NZglobal), &
                varf2(N_CELL0global,NZglobal), &
                varf3(N_CELL0global,NZglobal), &
                varf4(N_CELL0global,NZglobal))

       varf1 = 0.0d0
       varf2 = 0.0d0
       varf3 = 0.0d0
       varf4 = 0.0d0
!      __________________________________
!     |                                  |
!     |           Assign variable        |
!     |__________________________________|

        call matgloC(suf,varf1)
        call matgloC(svf,varf2)
        call matgloC(swf,varf3)
        call matgloC(spf,varf4)

      IF(rang_topo .EQ. 0) THEN
!      __________________________________
!     |                                  |
!     |             Open file            |
!     |__________________________________|

      ifile=16
      filedat = '../output/Matlab/PlotStatistics/TimeStats-    .dat'
      write(filedat(43:46),'(i4.4)') SaveCounter
      open(ifile,file=filedat,status='unknown')
!      __________________________________
!     |                                  |
!     |               Time               |
!     |__________________________________|

      write(ifile,79) scount
!      __________________________________
!     |                                  |
!     |         Fluid variables          |
!     |__________________________________|

      write(ifile,5) 'uf'
      write(ifile,35)((varf1(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'vf'
      write(ifile,35)((varf2(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'wf'
      write(ifile,35)((varf3(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'pf'
      write(ifile,35)((varf4(i,k),i=1,N_CELL0global),k=1,NZglobal)
!      __________________________________
!     |                                  |
!     |            Close file            |
!     |__________________________________|

      close(ifile)
      
      ENDIF
      call MPI_Barrier(comm3D,code)
#     endif

!     =============== END ================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

#     ifdef KeyParallel
      deallocate(varf1,varf2,varf3,varf4)
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: restart_out_stats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   [1.9]           WRITE DATA TO DISPLAY AT TECPLOT                  !
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

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: SavetecStats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |                                                        |
!     |                    Initialization                      |
!     |________________________________________________________|

4     format(A75)
5     format(8(1x,i8))
6     format(8(1x,e12.5))
7     format(a11,i7,a8,i7,a40)

!      ________________________________________________________
!     |         ====================================           |
!     |         ==========  SEQUENTIAL =============           |
!     |________________________________________________________|

#     ifndef KeyParallel
!         _________________________________
!        |                                 |
!        |    Allocate vertex variables    |
!        |_________________________________|

         allocate(uf_gl(N_VERT,NZ), &
                  vf_gl(N_VERT,NZ), &
                  wf_gl(N_VERT,NZ), &
                  pf_gl(N_VERT,NZ))

!         _________________________________
!        |                                 |
!        |    Calculate the mean value     |
!        |_________________________________|

         do k=1,NZ
            do i=1,N_CELL
              uf_gl(i,k) = suf(i,k)/scount
              vf_gl(i,k) = svf(i,k)/scount
              wf_gl(i,k) = swf(i,k)/scount
              pf_gl(i,k) = spf(i,k)/scount
            enddo
         enddo
!         __________________________________
!        |                                  |
!        |          Open file irec          |
!        |__________________________________|

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
!         __________________________________
!        |                                  |
!        |          Close file irec         |
!        |__________________________________|

         rewind(irec)
         close(irec)
         deallocate(uf_gl,vf_gl,wf_gl,pf_gl)
         
#     endif
!      ________________________________________________________
!     |         ====================================           |
!     |         =====  START PARALLEL OPTION =======           |
!     |________________________________________________________|

#     ifdef KeyParallel
!         _________________________________
!        |                                 |
!        |      Allocate global matrix     |
!        |_________________________________|

         allocate(xvt_gl(N_VERTglobal,NZglobal-1), &
                  yvt_gl(N_VERTglobal,NZglobal-1), &
                  zvt_gl(N_VERTglobal,NZglobal-1), &
                  uf_gl(N_CELL0global,NZglobal),   &
                  vf_gl(N_CELL0global,NZglobal),   &
                  wf_gl(N_CELL0global,NZglobal),   &
                  pf_gl(N_CELL0global,NZglobal))
!         _________________________________
!        |                                 |
!        |    Reconstruct global matrix    |
!        |_________________________________|

         call matgloV(xvt,xvt_gl)
         call matgloV(yvt,yvt_gl)
         call matgloV(zvt,zvt_gl)
         call matgloC(suf,uf_gl)
         call matgloC(svf,vf_gl)
         call matgloC(swf,wf_gl)
         call matgloC(spf,pf_gl)

         IF (rang_topo.eq.0) THEN
!         _________________________________
!        |                                 |
!        |      Calculate averaged data    |
!        |_________________________________|

         do k=1,NZglobal
            do i=1,N_CELL0global
              uf_gl(i,k) = uf_gl(i,k)/scount
              vf_gl(i,k) = vf_gl(i,k)/scount
              wf_gl(i,k) = wf_gl(i,k)/scount
              pf_gl(i,k) = pf_gl(i,k)/scount
            enddo
         enddo
!         __________________________________
!        |                                  |
!        |          Open file irec          |
!        |__________________________________|

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
!         __________________________________
!        |                                  |
!        |         Write vertices           |
!        |__________________________________|

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
!         __________________________________
!        |                                  |
!        |          Close file irec         |
!        |__________________________________|

         rewind(irec)
         close(irec)
         ENDIF

         call MPI_Barrier(comm3D,code)          
!         __________________________________
!        |                                  |
!        |      Deallocate global matrix    |
!        |__________________________________|
         
        deallocate(xvt_gl,yvt_gl,zvt_gl,&
                   uf_gl,vf_gl,wf_gl,pf_gl)

#     endif
!      ________________________________________________________
!     |                                                        |
!     |                     Finalization                       |
!     |________________________________________________________|

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----   End subroutine: SavetecStats'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                               END                                   !
!                             Jan 2018                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
