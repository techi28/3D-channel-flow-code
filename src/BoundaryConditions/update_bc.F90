!  ===================================================================!
!  ===================================================================!
!  ===================================================================!
      SUBROUTINE adjust_tagv(xc,yc,sig,dsig,No_cp,nbe, &
                           xv,yv,sigv,dsigv,No_vp,nbev)

!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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
      integer :: jv10,jv20,jv30,j0,ncount1,ncount2,jc
#     ifdef KeyParallel
      integer :: nvg
#     endif
!     -----------------------------------------
!     For no periodic bc, corner tag need to be fixed
               do nv=1,N_VERT
                  if(nbev(nv) .ne. 0) then
                     if(tagXV(nv)*tagYV(nv) .ne. 0) then
                       if(tagXV(nv) .eq. 2) then
                          tagXV(nv) = 0 ! priority given to y dir
!                          tagYV(nv) = 3 ! corner vertex
                       else
                          tagXV(nv) = 0
                       endif
                  endif
                 endif
               enddo

#      ifdef KeyParallel
               do nv=1,N_VERTglobal
                  if(nbev_global(nv) .ne. 0) then
                     if(tagXV_global(nv)*tagYV_global(nv) .ne. 0) then
                         if(tagXV_global(nv) .eq. 2) then
                           tagXV_global(nv) = 0
                           tagYV_global(nv) = 3
                         else
                           tagXV_global(nv) = 0
                         endif
                  endif
                 endif
               enddo

             if(rang_topo .eq. 1) print*, '       CORNER TAG FIXED!'
#      endif

      RETURN
      END
!  ===================================================================!
!  ===================================================================!
!  ===================================================================!
      SUBROUTINE adjust_bc(xc,yc,sig,dsig,No_cp,nbe, &
                           xv,yv,sigv,dsigv,No_vp,nbev)

!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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
!     -------------------------------------
!     local variables
      integer :: ii,elem,flag,fB,nvg
      integer :: countn
!      ________________________________________________________
!     |                                                        |
!     |             set tag for special ghost points           |
!     |________________________________________________________|
!     In case of a ghost boundary cell, if not periodic, the intersection
!     with its neigbouring inner cell is parrallel to normal vector,
!     hence the cross diffusion contribution is not necessary.
!     Aslo, the gradient of the special ghost point can not be updated using
!     standard boundary condition.
      do i=1,N_CELL
         tagsp(i) = 0
      enddo

      countn = 0

      do i=1,N_CELL0
        if(nbe(i) .ne. 0) then
            do j=1,3
                nc = No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
                  ii = No_cp(nc,3)
                  if(ii .eq. -1) then
                    tagsp(nc) = 1
 !                  ------------
 !                  correction of xe2 vector
                    xe2(i,j) = 0.
                    ye2(i,j) = 0.
#                   ifdef KeyKimChoi
                    dle1(i,j) = 1.0d0
                    dle2(i,j) = 0.
                    dle3(i,j) = 0.0d0
#                   endif
                    countn = countn + 1
                  endif
              endif
            enddo
           endif
       enddo

       print*, '      Special boundary points :', countn
!      ________________________________________________________
!     |                                                        |
!     |                  allocate bc info for NZ bc            |
!     |________________________________________________________|


         bctop   = 1
         bcbot   = 1

!      ________________________________________________________
!     |                                                        |
!     |                  allocate bc info for center           |
!     |________________________________________________________|
       allocate(bcuc(N_CELL),bcvc(N_CELL), &
                bcwc(N_CELL),bcpc(N_CELL))
!      ---------------------------------------------------------
!      assign initial value
       do i=1,N_CELL
          bcuc(i) = -1
          bcvc(i) = -1
          bcwc(i) = -1
          bcpc(i) = -1
       enddo
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component u                  |
!     |________________________________________________________|

!     --------------------------------------------------
!     HORIZONTAL
            DO i=1,N_CELL0
                if (nbe(i).ne.0) then
                    do j=1,3
                        nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        flag = 0
                        ii = No_cp(nc,3)
!                 -------------------------------------
!                  perodic point with pair in same proc
                        if(ii .gt. 0) then
                            fB = 3
                            flag = flag + 1
!                 -------------------------------------
!                  perodic point with pair in diff proc
                        elseif(ii .eq. 0) then
                            fB = -1
                            flag = flag + 1
!                 -------------------------------------
!                  Begin not periodic bc
                        elseif(ii .eq. -1) then
!                 -------------------------------------
!                      assign value for x dir ghost cell
                           if(TagXC(nc) .gt. 0) then
                              if(XBC .eq. 0) then
                                fB = 0
                                flag = flag +1
                             endif
                              if(XBC .eq. 1) then
                                fB = 0
                                flag = flag +1
                              endif
                              if(XBC .eq. 2) then
                                if(TagXC(nc) .eq. 1) then
                                  fB = 2
                                  flag = flag +1
                                elseif(TagXC(nc) .eq. 2) then
                                  fB = 1
                                  flag = flag +1
                                endif
                              endif
                           endif
!                 -------------------------------------
!                      assign value for y dir ghost cell
                           if(TagYC(nc) .gt. 0) then
                              if(YBC .eq. 0) then
                                fB = 0
                                flag = flag +1
                             endif
                              if(YBC .eq. 1) then
                                fB = 2
                                flag = flag +1
                              endif
                           endif
!                 -------------------------------------
!                  End not periodic bc
                        endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR U',flag
                         stop
                       else
                         bcuc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component v                  |
!     |________________________________________________________|
!     --------------------------------------------------
!     HORIZONTAL
    DO i=1,N_CELL0
	    if (nbe(i).ne.0) then
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        flag = 0
                        ii = No_cp(nc,3)
!                 -------------------------------------
!                  perodic point with pair in same proc
                        if(ii .gt. 0) then
                            fB = 3
                            flag = flag + 1
!                 -------------------------------------
!                  perodic point with pair in diff proc
                        elseif(ii .eq. 0) then
                            fB = -1
                            flag = flag + 1
!                 -------------------------------------
!                  Begin not periodic bc
                        elseif(ii .eq. -1) then
!                 -------------------------------------
!                      assign value for x dir ghost cell
                           if(TagXC(nc) .gt. 0) then
                              if(XBC .eq. 0) then
                                fB = 0
                                flag = flag +1
                              elseif(XBC .eq. 1) then
                                fB = 1
                                flag = flag +1
                              elseif(XBC .eq. 2) then
                                if(TagXC(nc) .eq. 1) then
                                  fB = 2
                                  flag = flag +1
                                elseif(TagXC(nc) .eq. 2) then !convective outflow
                                  fB = 1
                                  flag = flag +1
                                endif
                              endif
                           endif
!                 -------------------------------------
!                      assign value for y dir ghost cell
                           if(TagYC(nc) .gt. 0) then
                              if(YBC .eq. 0) then
                                fB = 0
                                flag = flag +1
                              endif
                              if(YBC .eq. 1) then
                                fB = 0
                                flag = flag +1
                              endif
                           endif
!                 -------------------------------------
!                  End not periodic bc
                        endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR V'
                         stop
                       else
                         bcvc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component w                  |
!     |________________________________________________________|
     DO i=1,N_CELL0
	    if (nbe(i).ne.0) then
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        flag = 0
                        ii = No_cp(nc,3)
!                 -------------------------------------
!                  perodic point with pair in same proc
                        if(ii .gt. 0) then
                            fB = 3
                            flag = flag + 1
!                 -------------------------------------
!                  perodic point with pair in diff proc
                        elseif(ii .eq. 0) then
                            fB = -1
                            flag = flag + 1
!                 -------------------------------------
!                  Begin not periodic bc
                        elseif(ii .eq. -1) then
!                 -------------------------------------
!                      assign value for x dir ghost cell
                           if(TagXC(nc) .gt. 0) then
                              if(XBC .eq. 0) then
                                fB = 0
                                flag = flag +1
                              elseif(XBC .eq. 1) then
                                fB = 1
                                flag = flag +1
                              elseif(XBC .eq. 2) then
                                if(TagXC(nc) .eq. 1) then
                                  fB = 2
                                  flag = flag +1
                                elseif(TagXC(nc) .eq. 2) then
                                  fB = 1
                                  flag = flag +1
                                endif
                              endif
                           endif
!                 -------------------------------------
!                      assign value for y dir ghost cell
                           if(TagYC(nc) .gt. 0) then
                              if(YBC .eq. 0) then
                                fB = 0
                                flag = flag +1
                             endif
                              if(YBC .eq. 1) then
                                fB = 1
                                flag = flag +1
                              endif
                           endif

                        endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR W'
                         stop
                       else
                         bcwc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component p                  |
!     |________________________________________________________|
      DO i=1,N_CELL0
	     if (nbe(i).ne.0) then
   	        do j=1,3
	           nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        flag = 0
                        ii = No_cp(nc,3)
!                 -------------------------------------
!                  perodic point with pair in same proc
                        if(ii .gt. 0) then
                            fB = 3
                            flag = flag + 1
!                 -------------------------------------
!                  perodic point with pair in diff proc
                        elseif(ii .eq. 0) then
                            fB = -1
                            flag = flag + 1
!                 -------------------------------------
!                  Begin not periodic bc
                        elseif(ii .eq. -1) then
!                            if(TagXC(nc) .eq. 2) then
!                                fB = 0
!                                flag = flag +1
!                            else
!                    FOR WET CASE, ALL NEWMANN BC ARE USED.
                                fB = 1
                                flag = flag +1
!                            endif
!                 -------------------------------------
!                  End not periodic bc
                        endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR P'
                         stop
                       else
                         bcpc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  allocate bc info for vertex           |
!     |________________________________________________________|
       allocate(bcuv(N_VERT),bcvv(N_VERT), &
                bcwv(N_VERT),bcpv(N_VERT))
!      ---------------------------------------------------------
!      assign initial value
       do nv=1,N_VERT
          bcuv(nv) = -1
          bcvv(nv) = -1
          bcwv(nv) = -1
          bcpv(nv) = -1
       enddo
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component u                  |
!     |________________________________________________________|
        DO nv=1,N_VERT
          if (nbev(nv) .ne. 0) then
            flag = 0
!     --------------------------------------------------
!           assign local info on outside boundary for x dir
            if(TagXV(nv) .gt. 0) then
!     --------------------------------------------------
!              Begin not periodic in x dir
                IF(XPB .eq. 0) THEN
                      if(XBC .eq. 0) then
                        fB = 0
                        flag = flag + 1
                      elseif(XBC .eq. 1) then
                        fB = 0
                        flag = flag + 1
                      elseif(XBC .eq. 2) then
                        if(TagXV(nv) .eq. 1) then
                         fB = 2
                         flag = flag + 1
                        endif
                        if(TagXV(nv) .eq. 2) then
                        fB = 1
                         flag = flag + 1
                        endif
                      endif
!     --------------------------------------------------
!              Begin periodic in x dir
                ELSE
                   fB = 3
                   flag = flag + 1
                ENDIF
             endif
!     --------------------------------------------------
!           assign local info on outside boundary for y dir
            if(TagYV(nv) .gt. 0) then
!     --------------------------------------------------
!               Begin not periodic in y dir
                IF(YPB .eq. 0) THEN
                      if(YBC .eq. 0) then
                        fB = 0
                        flag = flag + 1
                      endif
                      if(YBC .eq. 1) then
                        fB = 1
                        flag = flag + 1
                      endif
!     --------------------------------------------------
!              Begin periodic in x dir
               ELSE
                   fB = 3
                   flag = flag + 1
               ENDIF
             endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR U',flag
                 stop
            else
                 bcuv(nv) = fB
            endif

            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component v                  |
!     |________________________________________________________|
          DO nv=1,N_VERT
            if (nbev(nv) .ne. 0) then
            flag = 0
!     --------------------------------------------------
!           assign local info on outside boundary for x dir
            if(TagXV(nv) .gt. 0) then
!     --------------------------------------------------
!              Begin not periodic in x dir
                IF(XPB .eq. 0) THEN
                      if(XBC .eq. 0) then
                        fB = 0
                        flag = flag + 1
                      elseif(XBC .eq. 1) then
                        fB = 1
                        flag = flag + 1
                      elseif(XBC .eq. 2) then
                        if(TagXV(nv) .eq. 1) then
                         fB = 2
                         flag = flag + 1
                        endif
                        if(TagXV(nv) .eq. 2) then
                        fB = 1
                         flag = flag + 1
                        endif
                      endif
!     --------------------------------------------------
!              Begin periodic in x dir
               ELSE
                   fB = 3
                   flag = flag + 1
               ENDIF
             endif
!     --------------------------------------------------
!           assign local info on outside boundary for y dir
            if(TagYV(nv) .gt. 0) then
!     --------------------------------------------------
!               Begin not periodic in y dir
                IF(YPB .eq. 0) THEN
                      if(YBC .eq. 0) then
                        fB = 0
                        flag = flag + 1
                      endif
                      if(YBC .eq. 1) then
                        fB = 0
                        flag = flag + 1
                      endif
!     --------------------------------------------------
!              Begin periodic in y dir
               ELSE
                   fB = 3
                   flag = flag + 1
               endif
             endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR V'
                 stop
            else
                 bcvv(nv) = fB
            endif

            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component W                  |
!     |________________________________________________________|
          DO nv=1,N_VERT
            if (nbev(nv) .ne. 0) then
            flag = 0
!     --------------------------------------------------
!           assign local info on outside boundary for x dir
            if(TagXV(nv) .gt. 0) then
!     --------------------------------------------------
!              Begin not periodic in x dir
                IF(XPB .eq. 0) THEN
                      if(XBC .eq. 0) then
                        fB = 0
                        flag = flag + 1
                      elseif(XBC .eq. 1) then
                        fB = 1
                        flag = flag + 1
                      elseif(XBC .eq. 2) then
                        if(TagXV(nv) .eq. 1) then
                         fB = 2
                         flag = flag + 1
                        endif
                        if(TagXV(nv) .eq. 2) then
                        fB = 1
                         flag = flag + 1
                        endif
                      endif
!     --------------------------------------------------
!              Begin periodic in x dir
                ELSE
                   fB = 3
                   flag = flag + 1
                ENDIF
             endif
!     --------------------------------------------------
!           assign local info on outside boundary for y dir
            if(TagYV(nv) .gt. 0) then
!     --------------------------------------------------
!               Begin not periodic in y dir
                IF(YPB .eq. 0) THEN
                      if(YBC .eq. 0) then
                        fB = 0
                        flag = flag + 1
                      endif
                      if(YBC .eq. 1) then
                        fB = 1
                        flag = flag + 1
                      endif
!     --------------------------------------------------
!              Begin periodic in y dir
               ELSE
                   fB = 3
                   flag = flag + 1
               ENDIF
             endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR W'
                 stop
            else
                 bcwv(nv) = fB
            endif

            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component P                  |
!     |________________________________________________________|
          DO nv=1,N_VERT
            if (nbev(nv) .ne. 0) then
            flag = 0
!     --------------------------------------------------
!           assign local info on outside boundary for x dir
            if(TagXV(nv) .gt. 0) then
!     --------------------------------------------------
!              Begin not periodic in x dir
                IF(XPB .eq. 0) THEN
                    fB = 1
                    flag = flag +1
                ELSE
                   fB = 3
                   flag = flag + 1
                ENDIF
             endif
!     --------------------------------------------------
!           assign local info on outside boundary for y dir
            if(TagYV(nv) .gt. 0) then
!     --------------------------------------------------
!               Begin not periodic in y dir
                IF(YPB .eq. 0) then
                   fB = 1
                   flag = flag + 1
                ELSE
                   fB = 3
                   flag = flag + 1
                ENDIF
             endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR P'
                 stop
            else
                 bcpv(nv) = fB
            endif
!    -------------------------------------
!    the following is design to ensure the mass conservation
!    with periodic bc, meaning flux_total = 0 exactly by fixing
!    the corner vertex. In old version, if we only have perioidc
!    in one direction, we will encouter problem, though trivial,
!    but it is good to fix this.
#    ifdef KeyTESTpBC
!    Serial mode
#    ifndef KeyParallel
        if(vertex_edge(nv) .gt. 0) then
             bcpv(nv) = 3
        endif
        if(vertex_pair(nv) .gt. 0) then
             bcpv(nv) = 3
        endif
#    else
        nvg = index_globalv(nv)
        if(vpair_global(nvg) .gt. 0) then
             bcpv(nv) = 3
        endif
        if(vedge_global(nvg) .gt. 0) then
             bcpv(nv) = 3
        endif
#    endif
#    endif
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                      FINAL CHECK	               |
!     |________________________________________________________|
#     ifdef KeyTESTpBC
#     ifndef KeyParallel
!     ------------------------
!      check pre-calculated boundary info
      do nv=1,N_VERT
         if(nbev(nv) .ne. 0) then
            if(vertex_edge(nv) .eq. 1) then
               if(bcpv(nv) .ne. 3) then
                print*, 'Serial:Error for edge vertex in bcpv!'
                stop
               endif
           endif
           if(vertex_pair(nv) .gt. 0) then
               if(bcpv(nv) .ne. 3) then
                print*, 'Serial:Error for pair vertex in bcpv!'
                stop
               endif
           endif
        endif
      enddo
#     else
!     ------------------------
!      check pre-calculated boundary info
      do nv=1,N_VERT
        if(nbev(nv) .ne. 0) then
            nvg = index_globalv(nv)
            if(vedge_global(nvg) .eq. 1) then
               if(bcpv(nv) .ne. 3) then
                print*, 'Parallel:Error for edge vertex in bcpv!'
                stop
               endif
           endif
           if(vpair_global(nvg) .gt. 0) then
               if(bcpv(nv) .ne. 3) then
                print*, 'Parallel:Error for pair vertex in bcpv!'
                stop
               endif
           endif
         endif
      enddo
#     endif
#     endif

   RETURN
   END

!  ===================================================================!
!  ===================================================================!
!  ===================================================================!
      SUBROUTINE adjust_bc_BK(xc,yc,sig,dsig,No_cp,nbe, &
                              xv,yv,sigv,dsigv,No_vp,nbev)

!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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
!     -------------------------------------
!     local variables
      integer :: ii,elem,flag,fB,nvg
      integer :: jj,nv1,nv2,v1tag,v2tag
      integer :: countn
!------------------------------------------------------------
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'End subroutine: Update_bc_BK'
#     endif
!     ============================================================
!     IMPORTANT NOTE:
!     THE SUPPORT FOR BLUEKENUE MESH IN CASE OF A NON-STRAIGHT
!     FREE SLIP BOUNDARY CAN GIVE PROBLEM. SHOULD BE VERY CAREFEUL.
!     DEVELOPMENT UNDERWAY. SHOULD BE NO PROBLEM FOR OTHER BC CONDITIONS.
!     ============================================================
!      ________________________________________________________
!     |                                                        |
!     |             set tag for special ghost points           |
!     |________________________________________________________|
!     In case of a ghost boundary cell, if not periodic, the intersection
!     with its neigbouring inner cell is parrallel to normal vector,
!     hence the cross diffusion contribution is not necessary.
      do i=1,N_CELL
         tagsp(i) = 0
      enddo

      countn = 0

      do i=1,N_CELL0
        if(nbe(i) .ne. 0) then
            do j=1,3
                nc = No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
                  ii = No_cp(nc,3)
                  if(ii .eq. -1) then
                    tagsp(nc) = 1
 !                  ------------
 !                  correction of xe2 vector
                    xe2(i,j) = 0.
                    ye2(i,j) = 0.
                    countn = countn + 1
                  endif
              endif
            enddo
           endif
       enddo

       print*, '      Special boundary points :', countn
!      ________________________________________________________
!     |                                                        |
!     |                  allocate bc info for NZ bc            |
!     |________________________________________________________|

         bctop   = 1
         bcbot   = 1

!      ________________________________________________________
!     |                                                        |
!     |                  allocate bc info for center           |
!     |________________________________________________________|
       allocate(bcuc(N_CELL),bcvc(N_CELL), &
                bcwc(N_CELL),bcpc(N_CELL))
!      ---------------------------------------------------------
!      assign initial value
       do i=1,N_CELL
          bcuc(i) = -1
          bcvc(i) = -1
          bcwc(i) = -1
          bcpc(i) = -1
       enddo
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component u                  |
!     |________________________________________________________|

!     --------------------------------------------------
!     HORIZONTAL
            DO i=1,N_CELL0
              if (nbe(i).ne.0) then
                    do j=1,3
                        nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        jj = j+1
                        if(jj .gt. 3) jj = jj-3
                        nv1 = No_vp(i,j)
                        nv2 = No_vp(i,jj)
                        if(nbev(nv1) .gt. nbev(nv2)) then
                            v1tag = nbev(nv2)
                            v2tag = nbev(nv1)
                        else
                            v1tag = nbev(nv1)
                            v2tag = nbev(nv2)
                        endif
                    !  ---------------------------
                           flag = 0
                    !  -----------------------------
                    !   in case of a wall bc vertex

                           if(v2tag .eq. 1) then
                              fB = 0
                              flag = flag + 1
                           elseif(v2tag .eq. 2) then
                              fB = 1
                              flag = flag + 1
                           elseif(v2tag .eq. 3) then
                              fB = 2
                              flag = flag + 1
                           elseif(v2tag .eq. 4) then
                              fB = 1
                              flag = flag + 1
                           endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR U'
                         stop
                       else
                         bcuc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component v                  |
!     |________________________________________________________|
!     --------------------------------------------------
!     HORIZONTAL
    DO i=1,N_CELL0
	    if (nbe(i).ne.0) then
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        jj = j+1
                        if(jj .gt. 3) jj = jj-3
                        nv1 = No_vp(i,j)
                        nv2 = No_vp(i,jj)
                        if(nbev(nv1) .gt. nbev(nv2)) then
                            v1tag = nbev(nv2)
                            v2tag = nbev(nv1)
                        else
                            v1tag = nbev(nv1)
                            v2tag = nbev(nv2)
                        endif
                    !  ---------------------------
                           flag = 0
                    !  -----------------------------
                    !   in case of a wall bc vertex

                           if(v2tag .eq. 1) then
                              fB = 0
                              flag = flag + 1
                           elseif(v2tag .eq. 2) then
                              fB = 1
                              flag = flag + 1
                           elseif(v2tag .eq. 3) then
                              fB = 2
                              flag = flag + 1
                           elseif(v2tag .eq. 4) then
                              fB = 1
                              flag = flag + 1
                           endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR V'
                         stop
                       else
                         bcvc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component w                  |
!     |________________________________________________________|
     DO i=1,N_CELL0
	    if (nbe(i).ne.0) then
   	       do j=1,3
	          nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                        jj = j+1
                        if(jj .gt. 3) jj = jj-3
                        nv1 = No_vp(i,j)
                        nv2 = No_vp(i,jj)
                        if(nbev(nv1) .gt. nbev(nv2)) then
                            v1tag = nbev(nv2)
                            v2tag = nbev(nv1)
                        else
                            v1tag = nbev(nv1)
                            v2tag = nbev(nv2)
                        endif
                    !  ---------------------------
                           flag = 0
                    !  -----------------------------
                    !   in case of a wall bc vertex

                           if(v2tag .eq. 1) then
                              fB = 0
                              flag = flag + 1
                           elseif(v2tag .eq. 2) then
                              fB = 1
                              flag = flag + 1
                           elseif(v2tag .eq. 3) then
                              fB = 2
                              flag = flag + 1
                           elseif(v2tag .eq. 4) then
                              fB = 1
                              flag = flag + 1
                           endif
!                --------------------------------------
!                      check flag and assign new value
                       if(flag .ne. 1) then
                         print*, 'FATAL ERROR FOR BC CENTER FOR W'
                         stop
                       else
                         bcwc(nc) = fB
                       endif

                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component p                  |
!     |________________________________________________________|
      DO i=1,N_CELL0
	     if (nbe(i).ne.0) then
   	        do j=1,3
	           nc=No_cp(i,j)
!                 ====================================
!                 ==========  SEQUENTIAL =============
#                 ifndef KeyParallel
	          if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(nc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
!                 =============== END ================
!                 ====================================
                         bcpc(nc) = 1
                  endif
	          enddo
            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  allocate bc info for vertex           |
!     |________________________________________________________|
       allocate(bcuv(N_VERT),bcvv(N_VERT), &
                bcwv(N_VERT),bcpv(N_VERT))
!      ---------------------------------------------------------
!      assign initial value
       do nv=1,N_VERT
          bcuv(nv) = -1
          bcvv(nv) = -1
          bcwv(nv) = -1
          bcpv(nv) = -1
       enddo
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component u                  |
!     |________________________________________________________|
        DO nv=1,N_VERT
              v1tag = nbev(nv)
          if (v1tag .ne. 0) then
              flag = 0
              if(v1tag .eq. 1) then
                 fB = 0
                 flag = flag + 1
              elseif(v1tag .eq. 2) then
                 fB = 1
                 flag = flag + 1
              elseif(v1tag .eq. 3) then
                 fB = 2
                 flag = flag + 1
              elseif(v1tag .eq. 4) then
                 fB = 1
                 flag = flag + 1
              endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR U',flag
                 stop
            else
                 bcuv(nv) = fB
            endif

            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component v                  |
!     |________________________________________________________|
          DO nv=1,N_VERT
              v1tag = nbev(nv)
          if (v1tag .ne. 0) then
              flag = 0
              if(v1tag .eq. 1) then
                 fB = 0
                 flag = flag + 1
              elseif(v1tag .eq. 2) then
                 fB = 1
                 flag = flag + 1
              elseif(v1tag .eq. 3) then
                 fB = 2
                 flag = flag + 1
              elseif(v1tag .eq. 4) then
                 fB = 1
                 flag = flag + 1
              endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR V'
                 stop
            else
                 bcvv(nv) = fB
            endif

            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component W                  |
!     |________________________________________________________|
          DO nv=1,N_VERT
              v1tag = nbev(nv)
          if (v1tag .ne. 0) then
              flag = 0
              if(v1tag .eq. 1) then
                 fB = 0
                 flag = flag + 1
              elseif(v1tag .eq. 2) then
                 fB = 1
                 flag = flag + 1
              elseif(v1tag .eq. 3) then
                 fB = 2
                 flag = flag + 1
              elseif(v1tag .eq. 4) then
                 fB = 1
                 flag = flag + 1
              endif
!    --------------------------------------
!           check flag and assign new value
            if(flag .ne. 1) then
                 print*, 'FATAL ERROR FOR BC VERTEX FOR W'
                 stop
            else
                 bcwv(nv) = fB
            endif

            endif
          ENDDO
!      ________________________________________________________
!     |                                                        |
!     |                  Velocity component P                  |
!     |________________________________________________________|
          DO nv=1,N_VERT
              v1tag = nbev(nv)
          if (v1tag .ne. 0) then
!    --------------------------------------
!           check flag and assign new value
                 bcpv(nv) = 1
            endif
          ENDDO
!------------------------------------------------------------
#     ifdef KeyDbg
         write(*,'(t15,60a)'), 'End subroutine: Update_bc_BK'
#     endif
   RETURN
   END
