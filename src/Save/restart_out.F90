!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      WRITE RESULTS TO RESTART                       !
!                             Mar 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp, &
                             alphasnp,usnp,vsnp,wsnp,psnp, &
                             zct)

!---------------------------------------------------------------------!
!                                                                     !
!    This program writes all the variables needed to restart a        !
!    new simulation using the final values of the previous si-        !
!    mulation.                                                        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !
!  |______________|___________|____________________________________|  !
!  | --> alphafnp |(N_CELL,NZ)| Fluid control volume at t(n+1)     |  !
!  | --> ufnp     |(N_CELL,NZ)| Velocity component u_f at t(n+1)   |  !
!  | --> vfnp     |(N_CELL,NZ)| Velocity component v_f at t(n+1)   |  !
!  | --> wfnp     |(N_CELL,NZ)| Velocity component w_f at t(n+1)   |  !
!  | --> pfnp     |(N_CELL,NZ)| Pressure of the fluid at t(n+1)    |  !
!  |______________|___________|____________________________________|  !
!  | --> alphasnp |(N_CELL,NZ)| Solid control volume at t(n+1)     |  !
!  | --> usnp     |(N_CELL,NZ)| Velocity component u_s at t(n+1)   |  !
!  | --> vsnp     |(N_CELL,NZ)| Velocity component v_s at t(n+1)   |  !
!  | --> wsnp     |(N_CELL,NZ)| Velocity component w_s at t(n+1)   |  !
!  | --> psnp     |(N_CELL,NZ)| Pressure of the solid at t(n+1)    |  !
!  |______________|___________|____________________________________|  !
!  | --> zct      |(N_CELL,NZ)| zc at t(n+1)                       |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !
!  |_____________|_________________________________________________|  !
!  |--- N_CELL   | Total number of the cells                       |  !
!  |--- NZ       | Points in the sigma direction                   |  !
!  |    time     | Final time calculated                           |  !
!  |    dt       | Time step                                       |  !
!  | SaveCounter | Counting integer of saving files                |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  | * i,k       |  Loop counters                                  |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  |   ifile     |  Writing format                                 |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   ---  Parameters                                                   !
!        Common variables used                                        !
!    *   Common variables modified                                    !
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

      real*8,dimension(:,:):: alphafnp
      real*8,dimension(:,:):: ufnp
      real*8,dimension(:,:):: vfnp
      real*8,dimension(:,:):: wfnp
      real*8,dimension(:,:):: pfnp
      real*8,dimension(:,:):: alphasnp
      real*8,dimension(:,:):: usnp
      real*8,dimension(:,:):: vsnp
      real*8,dimension(:,:):: wsnp
      real*8,dimension(:,:):: psnp
      real*8,dimension(:,:):: zct

!      ________________________________________________________
!     |                                                        |
!     |           Definition of local variables                |
!     |________________________________________________________|
#      ifdef KeyParallel
      real*8,dimension(:,:),allocatable:: alphaf_gl
      real*8,dimension(:,:),allocatable:: uf_gl
      real*8,dimension(:,:),allocatable:: vf_gl
      real*8,dimension(:,:),allocatable:: wf_gl
      real*8,dimension(:,:),allocatable:: pf_gl
#      endif
       integer:: ifile

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: restart_out '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Formats                         |
!     |________________________________________________________|

77    format('         time = ',e12.5)
78    format('           dt = ',e12.5)
79    format('  SaveCounter = ',i8)
80    format('        nstep = ',i8)
35    format(6e12.5)
36    format(i8)
5     format(t2,80a)
!*********************************************************************!
!                                                                     !
!                        Write file (parallel)                        !
!                                                                     !
!*********************************************************************!
!     ====================================
!     ============  PARALLEL =============
#     ifdef KeyParallel
!      ________________________________________________________
!     |                                                        |
!     |                  allocate variable                     |
!     |________________________________________________________|
       allocate(alphaf_gl(N_CELL0global,NZglobal), &
                uf_gl(N_CELL0global,NZglobal), &
                vf_gl(N_CELL0global,NZglobal), &
                wf_gl(N_CELL0global,NZglobal), &
                pf_gl(N_CELL0global,NZglobal))
!     -----------------------------------
             alphaf_gl = 0.0d0
             uf_gl = 0.0d0
             vf_gl = 0.0d0
             wf_gl = 0.0d0
             pf_gl = 0.0d0
!      ________________________________________________________
!     |                                                        |
!     |                    assign variable                     |
!     |________________________________________________________|
        call matgloC(alphafnp,alphaf_gl)
        call matgloC(ufnp,uf_gl)
        call matgloC(vfnp,vf_gl)
        call matgloC(wfnp,wf_gl)
        call matgloC(pfnp,pf_gl)
!      ----------------
      IF(rang_topo .EQ. 0) THEN
!      ________________________________________________________
!     |                                                        |
!     |                  Open file ifile                       |
!     |________________________________________________________|
      ifile=16
      write(filerepout(11:15),'(i4.4)') SaveCounter
      open(ifile,file=filerepout,status='unknown')
!      __________________________________
!     |                                  |
!     |               Time               |
!     |__________________________________|

      write(ifile,77) time
      write(ifile,79) SaveCounter
      write(ifile,80) nstep
!      __________________________________
!     |                                  |
!     |         Fluid variables          |
!     |__________________________________|

      write(ifile,5) 'alphaf'
      write(ifile,35)((alphaf_gl(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'uf'
      write(ifile,35)((uf_gl(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'vf'
      write(ifile,35)((vf_gl(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'wf'
      write(ifile,35)((wf_gl(i,k),i=1,N_CELL0global),k=1,NZglobal)
      write(ifile,5) 'pf'
      write(ifile,35)((pf_gl(i,k),i=1,N_CELL0global),k=1,NZglobal)
!      ________________________________________________________
!     |                                                        |
!     |                      Close file                        |
!     |________________________________________________________|

      close(ifile)
      ENDIF
#     endif

!     =============== END ================
!     ====================================
!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!
#     ifdef KeyParallel
        deallocate(alphaf_gl,uf_gl,vf_gl,wf_gl,pf_gl)
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: restart_out'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END FINALIZATION                            !
!                             Mar 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
