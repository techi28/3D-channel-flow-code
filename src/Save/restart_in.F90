!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      WRITE RESULTS TO RESTART                       !
!                             Mar 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE restart_in(alphafn,ufn,vfn,wfn,pfn, &
                            alphasn,usn,vsn,wsn,psn, &
                            zct)

!---------------------------------------------------------------------!
!                                                                     !
!    This program read all the variables needed to restart a          !
!    new simulation using the final values of the previous si-        !
!    mulation. We need the name of the input restart file de-         !
!    fined in the input subroutine as "filerepin".                    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name     |     Size  |        Description                 |  !
!  |______________|___________|____________________________________|  !
!  | <-- zct      |(N_CELL,NZ)| zc at t(n)                         |  !
!  | <-- alphafn  |(N_CELL,NZ)| Fluid control volume at t(n)       |  !
!  | <-- ufn      |(N_CELL,NZ)| Velocity component u_f at t(n)     |  !
!  | <-- vfn      |(N_CELL,NZ)| Velocity component v_f at t(n)     |  !
!  | <-- wfn      |(N_CELL,NZ)| Velocity component w_f at t(n)     |  !
!  | <-- pfn      |(N_CELL,NZ)| Pressure of the fluid at t(n)      |  !
!  |______________|___________|____________________________________|  !
!  | <-- alphasn  |(N_CELL,NZ)| Solid control volume at t(n)       |  !
!  | <-- usn      |(N_CELL,NZ)| Velocity component u_s at t(n)     |  !
!  | <-- vsn      |(N_CELL,NZ)| Velocity component v_s at t(n)     |  !
!  | <-- wsn      |(N_CELL,NZ)| Velocity component w_s at t(n)     |  !
!  | <-- psn      |(N_CELL,NZ)| Pressure of the solid at t(n)      |  !
!  |______________|___________|____________________________________|  !
!                                                                     !
!    Common parameters & variables used:                              !
!   _______________________________________________________________   !
!  |   Name      |                  Description                    |  !
!  |_____________|_________________________________________________|  !
!  | --- N_CELL  | Total number of the cells                       |  !
!  | --- NZ      | Points in the sigma direction                   |  !
!  |     time    | Final time calculated                           |  !
!  |     dt      | Time step                                       |  !
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
!  |   title     |  Title reading and printing                     |  !
!  |_____________|_________________________________________________|  !
!  |   N_total   |  Length of the vector of total number of element|  !
!  |   m         |  Loop counter of aux vectors                    |  !
!  |   zctaux    | (N_total) zc auxiliar vector                    |  !
!  |   alphafaux | (N_total) Fluid control volume auxiliar vector  |  !
!  |   ufaux     | (N_total) Velocity component u_f auxiliar vector|  !
!  |   vfaux     | (N_total) Velocity component v_f auxiliar vector|  !
!  |   wfaux     | (N_total) Velocity component w_f auxiliar vector|  !
!  |   pfaux     | (N_total) Pressure of the fluid auxiliar vector |  !
!  |   alphasaux | (N_total) Solid control volume auxiliar vector  |  !
!  |   usaux     | (N_total) Velocity component u_s auxiliar vector|  !
!  |   vsaux     | (N_total) Velocity component v_s auxiliar vector|  !
!  |   wsaux     | (N_total) Velocity component w_s auxiliar vector|  !
!  |   psaux     | (N_total) Pressure of the solid auxiliar vector |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   <--  Output variables                                             !
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

      real*8,dimension(:,:):: alphafn
      real*8,dimension(:,:):: ufn
      real*8,dimension(:,:):: vfn
      real*8,dimension(:,:):: wfn
      real*8,dimension(:,:):: pfn
      real*8,dimension(:,:):: alphasn
      real*8,dimension(:,:):: usn
      real*8,dimension(:,:):: vsn
      real*8,dimension(:,:):: wsn
      real*8,dimension(:,:):: psn
      real*8,dimension(:,:):: zct

!      ________________________________________________________
!     |                                                        |
!     |           Definition of local variables                |
!     |________________________________________________________|

      integer :: N_total
      real*8, dimension(:),allocatable :: alphafaux,ufaux,vfaux,wfaux
      real*8, dimension(:),allocatable :: pfaux
      real*8, dimension(:),allocatable :: alphasaux,usaux,vsaux,wsaux
      real*8, dimension(:),allocatable :: psaux
      real*8, dimension(:),allocatable :: zctaux
      real*8, dimension(:),allocatable :: u1aux,u2aux,u3aux,uTaux,uBaux
      real*8, dimension(:),allocatable :: uoaux,voaux,woaux
      integer:: m
      integer:: ifile
      character*80 title

!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: restart_in'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                     Reading the re-start data                       !
!                                                                     !
!*********************************************************************!
#     ifdef KeyParallel
            N_total = N_CELL0global*NZ
#     else
            N_total = N_CELL0*NZ
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                        Allocate                        |
!     |________________________________________________________|


      allocate(alphafaux(N_total),                           &
               ufaux(N_total),vfaux(N_total),wfaux(N_total), &
               pfaux(N_total))

        if(Isolid .eq. 1) then
            allocate(alphasaux(N_total),                      &
                usaux(N_total),vsaux(N_total),wsaux(N_total), &
                psaux(N_total))
        endif

!      ________________________________________________________
!     |                                                        |
!     |                        Formats                         |
!     |________________________________________________________|

77    format(16x,e12.5)
78    format(16x,e12.5)
79    format(16x,i8)
5     format(82a)
35    format(6e12.5)
36    format(i8)


!        __________________________________
!       |                                  |
!       |      Open file: filerepin        |
!       |__________________________________|

         ifile=16
         open(ifile,file=filerepin,status='old')
!        __________________________________
!       |                                  |
!       |              Time                |
!       |__________________________________|

         read(ifile,77) time
         read(ifile,79) SaveCounter
         read(ifile,79) nstep
!        __________________________________
!       |                                  |
!       |         Fluid variables          |
!       |__________________________________|

         read(ifile,5) title
         read(ifile,35)(alphafaux(m),m=1,N_total)
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
!       |         Solid variables          |
!       |__________________________________|
!            if(Isolid .eq. 1) then
!                read(ifile,5) title
!                read(ifile,35)(alphasaux(m),m=1,N_total)
!                read(ifile,5) title
!                read(ifile,35)(usaux(m),m=1,N_total)
!                read(ifile,5) title
!                read(ifile,35)(vsaux(m),m=1,N_total)
!                read(ifile,5) title
!                read(ifile,35)(wsaux(m),m=1,N_total)
!                read(ifile,5) title
!                read(ifile,35)(psaux(m),m=1,N_total)
!            endif
!     ====================================
!     ===========   PARALLEL  ===========
!        __________________________________
!       |                                  |
!       |        Assign final matrices     |
!       |__________________________________|

        do k=1,NZ
           do i=1,N_CELL0
              nc = index_global(i)
              m = (k-1)*N_CELL0global + nc
!             ---------------------------
              alphafn(i,k)= alphafaux(m)
              ufn(i,k) = ufaux(m)
              vfn(i,k) = vfaux(m)
              wfn(i,k) = wfaux(m)
              pfn(i,k) = pfaux(m)
!             ---------------------------
                if(Isolid .eq. 1) then
                    alphasn(i,k)= alphasaux(m)
                    usn(i,k) = usaux(m)
                    vsn(i,k) = vsaux(m)
                    wsn(i,k) = wsaux(m)
                    psn(i,k) = psaux(m)
                endif
           enddo
        enddo

!      __________________________________
!     |                                  |
!     |            Close file            |
!     |__________________________________|

      close(ifile)
!     =============== END ================
!     ====================================
!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                        Deallocate                      |
!     |________________________________________________________|

      deallocate(alphafaux,ufaux,vfaux,wfaux,pfaux)
      if(Isolid .eq. 1) then
          deallocate(alphasaux,usaux,vsaux,wsaux,psaux)
      endif
!     ~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: restart_in'
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
