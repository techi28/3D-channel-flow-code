!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               3D BOUNDARY CONDITION FOR VELOCITY FIELD              !
!                               Dec 2015                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      SUBROUTINE calcul_bc(xc,yc,sig,dsig,No_cp,nbe, &
                           xv,yv,sigv,dsigv,No_vp,nbev)
!---------------------------------------------------------------------!
!            Description of Boundary Condition Tags                   !
!---------------------------------------------------------------------!
!     XDIni            = real Initial streamwise value                !
!     XDFin            = real Final streamwise value                  !
!     YDIni            = real Initial spanwise value                  !
!     YDFin            = real Final spanwise value                    !
!     XPB              = 1 Periodic BC in streamwise direction        !
!                      = 0                                            !
!     YPB              = 1 Periodic BC in spanwise direction          !
!                      = 0                                            !
!---------------------------------------------------------------------!
!     TopBC            = 0 No slip wall                               !
!                      = 1 Free slip wall                             !
!                      = 2 Special for Cavity Test case               !
!---------------------------------------------------------------------!
!     BotBC            = 0 No slip wall                               !
!                      = 1 Free slip wall                             !
!---------------------------------------------------------------------!
!     XBC              = 0 No slip wall                               !
!                      = 1 Free slip wall                             !
!                      = 2 Inflow and outflow                         !
!---------------------------------------------------------------------!
!     YBC              = 0 No slip wall                               !
!                      = 1 Free slip wall                             !
!---------------------------------------------------------------------!
!     ChooseBoundary   = 1 Manually set boundary condition            !
!                      = 0 Use input info from BlueKenue (data.txt)   !
!---------------------------------------------------------------------!
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
!     ------------------------------------------
!      LOCAL variable
      integer :: IDISPLAY
      IDISPLAY = 1
#     ifdef KeyParallel
         if(rang_topo .ne. 0) IDISPLAY = 0
#     endif
!     ------------------------------------------
        if(IDISPLAY .eq. 1) then
            if(ChooseBoundary .eq. 0) then
                print*, '      USE BLUEKENUE INPUT BC'
            elseif(ChooseBoundary .eq. 1) then
                print*, '      DOMAIN SIZE: X',XDIni,XDFin
                print*, '      DOMAIN SIZE: Y',YDIni,YDFin
                print*, '      DOMAIN SIZE: Z',SigIni,SigFin
                print*,  '     USE MANNUALLY SET BC'
            else
                print*, '      ERROR IN BC. EXIT'
                stop
            endif
!    -------------------------------------------
!      For Top Cells
         IF (ZPB .eq. 1) then
            print*, '      PERIODIC IN Z DIRECTION'
            print*, '      USE WITH CAUTION'
         ELSE
            if(TopBC .eq. 0) then
                print*, '      TOP BC = NO SLIP WALL'
            elseif(TopBC .eq. 1) then
                print*, '      TOP BC = FREE SLIP WALL'
            elseif(TopBC .eq. 2) then
                print*, '      TOP BC = FIXED VELOCITY'
            else
                print*, '      UNDEFINED TOP BC USED. EXIT'
                stop
            endif
!     ------------------------------------------------
!      For Bottom Cells
            if(BotBC .eq. 0) then
                print*, '      BOT BC = NO SLIP WALL'
            elseif(BotBC .eq. 1) then
                print*, '      BOT BC = FREE SLIP WALL'
            else
                print*, '      UNDEFINED BOT BC USED. EXIT'
                stop
            endif
        ENDIF

        if(ChooseBoundary .eq. 1) then
!    ----------------------------------------------
#    ifndef KeyTESTpBC
        if((XPB .eq. 1) .or. (YPB .eq. 1)) then
            print*, '      ERROR.PERIODIC BC SUPPORT IS TURNED OFF. EXIT'
            print*, '      Check common.mpf and cppdefs.h'
            stop
        endif
#    endif
!     ----------------------------------------------
!       For horizontal cells
       if(XPB*YPB .eq. 1) then
            print*, '      PERIODIC BC IN XY PLANE.'
       else
!      ------------------------------------------------
            IF(XPB .eq. 0) THEN
                if(XBC .eq. 0) then
                    print*, '      XBC = NO SLIP WALL'
                elseif(XBC .eq. 1) then
                    print*, '      XBC = FREE SLIP WALL'
                elseif(XBC .eq. 2) then
                    print*, '      XBC = INFLOW AND OUTFLOW'
                else
                    print*, '      ERROR IN X BC. EXIT'
                    stop
                endif
            ELSE
                print*, '      XBC = PERIODIC'
            ENDIF
!     -------------------------------------------------
            IF(YPB .eq. 0) THEN
                if(YBC .eq. 0) then
                    print*, '      YBC = NO SLIP WALL'
                elseif(YBC .eq. 1) then
                    print*, '      YBC = FREE SLIP WALL'
                else
                    print*, '      ERROR IN Y BC. EXIT'
                    stop
                endif
            ELSE
                print*, '      YBC = PERIODIC'
            ENDIF
       endif ! for non all periodic
       endif ! for manually set bc
     endif ! for IDisplay = 1


      if(ChooseBoundary .eq. 1) then
!     =============================================
!     adjustment of tagxv,tagyv using the bc condition
      call adjust_tagv(xc,yc,sig,dsig,No_cp,nbe, &
                       xv,yv,sigv,dsigv,No_vp,nbev)
!     =============================================
!     adjustment of bc info on ghost cells
      call adjust_bc(xc,yc,sig,dsig,No_cp,nbe, &
                       xv,yv,sigv,dsigv,No_vp,nbev)
!     =============================================
      else
       call adjust_bc_BK(xc,yc,sig,dsig,No_cp,nbe, &
                       xv,yv,sigv,dsigv,No_vp,nbev)
      endif
!     ---------------------------------------------
!     NZ partition
        if(IDISPLAY .eq. 1) then
            if(NZBlock .gt. 1) then
                print*, '      NZ DIVIDED INTO',NZBlock,'PARTS'
            else
                print*, '      NZ NOT PARTITIONED'
            endif

            print*, '      INITIAL BC SET COMPLETED!'
        endif

   RETURN
   END
!  ===================================================================!
!  ===================================================================!
!  ===================================================================!
   SUBROUTINE  BCglobalVC(phiuc,phivc,phiwc,             &
                          phiuv,phivv,phiwv,             &
                          xc,yc,sig,dsig,No_cp,nbe,      &
                          xv,yv,sigv,dsigv,No_vp,nbev,tagdir)
!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
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

      real*8,dimension(:,:) :: phiuc(N_CELL,NZ)
      real*8,dimension(:,:) :: phivc(N_CELL,NZ)
      real*8,dimension(:,:) :: phiwc(N_CELL,NZ)
      real*8,dimension(:,:) :: phiuv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: phivv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: phiwv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
 !     --------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
      integer :: tagdir
!     ---------------------------------------
!     U Boundary Condition
!     ---------------------------------------
          call BCvelcenter3D(phiuc,xc,yc,sig,dsig,No_cp,nbe,TagBCu)
!     ----------------------------------------
!     Interpolation
          call interpolation3D(phiuv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phiuc,xc,yc,sig,dsig,No_cp,nbe,TagBCu,tagdir)
!     ----------------------------------------
!     Vertex Boundary Condition
          call BCvelvertex3D(phiuv,xv,yv,sigv,dsigv,No_vp,nbev,TagBCu,tagdir)

!     ---------------------------------------
!     V Boundary Condition
!     ---------------------------------------
          call BCvelcenter3D(phivc,xc,yc,sig,dsig,No_cp,nbe,TagBCv)
!     ----------------------------------------
!     Interpolation
          call interpolation3D(phivv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phivc,xc,yc,sig,dsig,No_cp,nbe,TagBCv,tagdir)
!     ----------------------------------------
!     Vertex Boundary Condition
          call BCvelvertex3D(phivv,xv,yv,sigv,dsigv,No_vp,nbev,TagBCv,tagdir)

!     ---------------------------------------
!     W Boundary Condition
!     ---------------------------------------
          call BCvelcenter3D(phiwc,xc,yc,sig,dsig,No_cp,nbe,TagBCw)
!     ----------------------------------------
!     Interpolation
          call interpolation3D(phiwv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phiwc,xc,yc,sig,dsig,No_cp,nbe,TagBCw,tagdir)
!     ----------------------------------------
!     Vertex Boundary Condition
          call BCvelvertex3D(phiwv,xv,yv,sigv,dsigv,No_vp,nbev,TagBCw,tagdir)

   RETURN
   END
!  ==================================================================!
!  ==================================================================!
!  ==================================================================!
   SUBROUTINE  BCglobalP(phi,phiv,                      &
                         xc,yc,sig,dsig,No_cp,nbe,      &
                         xv,yv,sigv,dsigv,No_vp,nbev,tagdir)
!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
!     --------------------------------------
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
      integer :: tagdir

!     ---------------------------------------
!     P Boundary Condition
!     ---------------------------------------
          call BCpcenter(phi,xc,yc,sig,dsig,No_cp,nbe)
!     ----------------------------------------
!     Interpolation
          call interpolation3DP(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                                phi,xc,yc,sig,dsig,No_cp,nbe,tagdir)
!     ----------------------------------------
!     Vertex Boundary Condition
          call BCpvertex(phiv,xv,yv,sigv,dsigv,No_vp,nbev)

   RETURN
   END
!  ==================================================================!
!  ==================================================================!
!  ==================================================================!
     SUBROUTINE BCvelcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct cell-center boundary         !
!    condition of the velocity. The tag called "tagBC" is used to     !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | <--- phi    |(N_CELL,NZ)| Function at the cell-center         |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | ---> xc,yc  |(N_CELL)   | Coordinates of the cell centers     |  !
!  | ---> sig    |(NZ)       | sigma value at the cell centers     |  !
!  | ---> dsig   |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | ---> No_cp  |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | ---> nbe    |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer :: tagBC,elem
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z,fB
      real*8 :: inflow,dirchlet
      integer :: ii,jj,kk,flag

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCvelcenter3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      ________________________________________________________
!     |                                                        |
!     |                     FOR ZPB = 0                        |
!     |________________________________________________________|
      IF(ZPB .ne. 1) THEN

!     ________________________________________________________
!     |                                                        |
!     |                       Bottom BC                        |
!     |________________________________________________________|
        if(bcbot .eq. 1) then
            if(BotBC .eq. 0) then
                do i=1,N_CELL0
                    phi(i,1) = -phi(i,2)
                enddo

            elseif(BotBC .eq. 1) then
                if(tagBC .eq. 3) then
                    DO i=1,N_CELL0
                        phi(i,1) = -phi(i,2)
                    enddo
                else
                    DO i=1,N_CELL0
                        phi(i,1) = phi(i,2)
                    enddo
                endif
            endif
        endif
!      ________________________________________________________
!     |                                                        |
!     |                       Top BC                           |
!     |________________________________________________________|
        if(bctop .eq. 1) then
            if(TopBC .eq. 0) then
                do i=1,N_CELL0
                    phi(i,NZ) = -phi(i,NZ-1)
                enddo

            elseif(TopBC .eq. 1) then
                if(tagBC .eq. 3) then
                    DO i=1,N_CELL0
                        phi(i,NZ) = -phi(i,NZ-1)
                    enddo
                else
                    DO i=1,N_CELL0
                        phi(i,NZ) = phi(i,NZ-1)
                    enddo
                endif

            elseif(TopBC .eq. 2) then
                if(tagBC .eq. 1) then
                    DO i=1,N_CELL0
                        phi(i,NZ) = 2.0d0-phi(i,NZ-1)
                    enddo
                else
                    DO i=1,N_CELL0
                        phi(i,NZ) = -phi(i,NZ-1)
                    enddo
                endif
            endif
        endif
!      ________________________________________________________
!     |                                                        |
!     |                  FOR ZPB = 1                           |
!     |________________________________________________________|
      ELSE
            if((bctop + bcbot) .eq. 2) then
                do i=1,N_CELL0
                    phi(i,1)  = phi(i,NZ-1)
                    phi(i,NZ) = phi(i,2)
                enddo
            endif
      ENDIF
!      ________________________________________________________
!     |                                                        |
!     |                   BC for overlapping points            |
!     |________________________________________________________|
!      Communication between procs
!      In this part, the overlapping cells from phi(N_CELL0,2:NZ-1)
#     ifdef KeyParallel
        call MPI_Barrier(comm3D,code)
        call communication3D(phi)
#     endif
!     --------------------------------------------------
!     HORIZONTAL
        DO k=1,NZ
            DO i=1,N_CELL0
                if (nbe(i).ne.0) then
                    do j=1,3
                        nc=No_cp(i,j)
                        if(tagBC .eq. 1) flag = bcuc(nc)
                        if(tagBC .eq. 2) flag = bcvc(nc)
                        if(tagBC .eq. 3) flag = bcwc(nc)
!                 ------------------------------
!                 wall bc + dirchlet bc
                        if(flag .eq. 0) then
!                            x = 0.5*(xc(nc)+xc(i))
!                            y = 0.5*(yc(nc)+yc(i))
!                            z = sig(k)
!                            fB = 2.0d0*dirchlet(x,y,z,time,Re,tagBC)
                            fB = 0.
                            fB = fB-phi(i,k)
                        endif
!                 ------------------------------
!                 newmann bc
                        if(flag .eq. 1) fB = phi(i,k)
!                 ------------------------------
!                 inflow  bc
                        if(flag .eq. 2) then
                            x = 0.5*(xc(nc)+xc(i))
                            y = 0.5*(yc(nc)+yc(i))
                            z = sig(k)
                            fB = 2.0d0*inflow(x,y,z,tagBC)-phi(i,k)
                        endif
!                 ------------------------------
!                 periodic bc
                        if(flag .eq. 3) then
                            ii = No_cp(nc,3)
                            fB = phi(ii,k)
                        endif
!                --------------------------------------
!                  update local value
                        if(flag .ne. -1)  phi(nc,k) = fB
                    enddo
                endif
            ENDDO
        ENDDO
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCvelcenter3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  VELOCITY VERTEX BOUNDARY CONDITION                 !
!                             Nov 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCvelvertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,tagBC,tagdir)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct vertex boundary condition    !
!    of the velocity components. The tag called "tagBC" is used to    !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- phiv   |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | ---> xv,yv  |(N_VERT)   | Coordinates of the vertices         |  !
!  | ---> sigv   |(NZ-1)     | sigma value at the vertices         |  !
!  | ---> dsigv  |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | ---> No_vp  |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | ---> nbev   |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
      integer :: tagBC
      integer :: tagdir
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: x,y,z
      integer :: pair,ncount,j0,nvg,flag
      real*8 :: inflow,dirchlet
      real*8 :: sumdlV
#     ifdef KeyParallel
      real*8,dimension(:,:),allocatable :: phiv_global
      real*8,dimension(:),allocatable   :: phiv2_global
      real*8,dimension(:,:),allocatable :: sumdlV_global
#     else
      integer,dimension(:) :: vpair(1:4)
      real*8,dimension(:),allocatable :: oldphiv
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCvelvertex3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!              Dirichlet Boundary Condition for the velocity          !
!                                                                     !
!*********************************************************************!
!                 ===============
!                 == HORIZONTAL =
!                 ===============
            DO nv=1,N_VERT
                IF(nbev(nv).ne.0) THEN
                    if(tagBC .eq. 1) flag = bcuv(nv)
                    if(tagBC .eq. 2) flag = bcvv(nv)
                    if(tagBC .eq. 3) flag = bcwv(nv)
                    if(flag .eq. 0)  then
                        do k=1,NZ-1
                            x = xv(nv)
                            y = yv(nv)
                            z = sigv(k)
                            phiv(nv,k) = 0.!dirchlet(x,y,z,time,Re,tagBC)
                        enddo
                    endif
                    if(flag .eq. 2) then
                        do k=1,NZ-1
                            x = xv(nv)
                            y = yv(nv)
                            z = sigv(k)
                            phiv(nv,k) = inflow(x,y,z,tagBC)
                        enddo
                    endif
                ENDIF
            ENDDO
!     --------------------------------------------------
!                 ===============
!                 ==  VERTICAL ==
!                 ===============
        IF(tagdir .eq. 3) THEN
        if(ZPB .ne. 1 ) then
!       =======================
!       =========BOT BC =======
!       =======================
            if(bcbot .eq. 1) then
                if(botBC .eq. 0) then
                    do nv=1,N_VERT
                        phiv(nv,1) = 0.
                    enddo

                elseif(botBC .eq. 1) then
                    if(tagBC .eq. 3) then
                        do nv=1,N_VERT
                            phiv(nv,1) = 0.
                        enddo
                    endif
                endif
            endif
!       =======================
!       =========TOP BC =======
!       =======================
            if(bctop .eq. 1) then
                if(topBC .eq. 0) then
                    do nv=1,N_VERT
                        phiv(nv,NZ-1) = 0.
                    enddo

                elseif(topBC .eq. 1) then
                    if(tagBC .eq. 3) then
                        do nv=1,N_VERT
                            phiv(nv,NZ-1) = 0.
                        enddo
                    endif

                elseif(topBC .eq. 2) then
                    if(tagBC .eq. 1) then
                        do nv=1,N_VERT
                            if(nbev(nv).eq.0) then
                                phiv(nv,NZ-1) = 1.0d0
                            endif
                        enddo
                    else
                        do nv=1,N_VERT
                            phiv(nv,NZ-1) = 0.
                        enddo
                    endif
                endif
            endif
        endif
        ENDIF ! END 3D
!     --------------------------------------------------
!     Periodic BC
#     ifdef KeyTESTpBC
        call bc_periodic(phiv,xv,yv,sigv,dsigv,No_vp,nbev,tagBC)
#    endif

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCvelvertex3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!  ===================================================================!
!  ===================================================================!
!  ===================================================================!
!  ===================================================================!
   SUBROUTINE BCpcenter(phi,xc,yc,sig,dsig,No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the boundary condition for the pressure  !
!    in the case of free surface problems at the cell-centers.        !
!                                                                     !
!    In the vertical walls:                                           !
!                              dpdn = 0               (Neumann)       !
!    At the free surface:                                             !
!                                 p = 0               (Dirichlet)     !
!    At the bottom boundary:                                          !
!                            dpdsig = -rho*H*dwdt     (Neumann)       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | <--> phi    |(N_CELL,NZ)| Function at the cell-center         |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> Hpr     |(N_CELL)   | Total depth = h + etan  cell-center |  !
!  | --> h       |(N_CELL)   | still depth             cell-center |  !
!  | --> etan    |(N_CELL)   | free surface            cell-center |  !
!  |_____________|___________|_____________________________________|  !
!  | --> rhof    |(N_CELL,NZ)| Density                             |  !
!  | --> dwdtB   |(N_CELL)   | d(w)/dt at the bottom boundary      |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: fB
      integer :: elem,ii,flag

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCpcenter'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                            Initilaization                           !
!                                                                     !
!*********************************************************************!
!     --------------------------------------------------
!     VERTICAL
            if(ZPB .eq. 0) then
                if(bctop .eq. 1) then
                    do i=1,N_CELL0
                        phi(i,NZ) = phi(i,NZ-1)
                    enddo
                endif

                if(bcbot .eq. 1) then
                    do i=1,N_CELL0
                        phi(i,1) = phi(i,2)
                    enddo
                endif
            else
                if((bctop+bcbot) .eq. 2) then
                    do i=1,N_CELL0
                        phi(i,NZ) = phi(i,2)
                        phi(i,1) = phi(i,NZ-1)
                    enddo
                endif
            endif
!     -------------------------------------------
!      Overlapping cells
#     ifdef KeyParallel
          call MPI_Barrier(comm3D,code)
          call communication3D(phi)
#     endif
!     -------------------------------------------
!      Boundary Cells
        DO k=1,NZ
            DO i=1,N_CELL0
                if (nbe(i).ne.0) then
                    do j=1,3
                        nc=No_cp(i,j)
                        flag = bcpc(nc)
                        if(flag .eq. 0) fB = -phi(i,k)
                        if(flag .eq. 1) fB = phi(i,k)
                        if(flag .eq. 3) then
                            ii = No_cp(nc,3)
                            fB = phi(ii,k)
                        endif
                        if(flag .ne. -1)  phi(nc,k) = fB
                    enddo
                endif
            ENDDO
        ENDDO
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCpcenter'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!   PRESSURE BOUNDARY CONDITIONS FOR FREE SURFACE PROBLEMS (vertex)   !
!                             March 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCpvertex(phiv,xv,yv,sigv,dsigv,No_vp,nbev)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the boundary condition for the pressure  !
!    in the case of free surface problems at the vertex points.       !
!                                                                     !
!    In the vertical walls:                                           !
!                              dpdn = 0               (Neumann)       !
!    At the free surface:                                             !
!                                 p = 0               (Dirichlet)     !
!    At the bottom boundary:                                          !
!                            dpdsig = -rho*H*dwdt     (Neumann)       !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !
!  |_____________|_____________|___________________________________|  !
!  | <-> phiv    |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> Hpr     |(N_CELL)   | Total depth = h + etan  cell-center |  !
!  | --> h       |(N_CELL)   | still depth             cell-center |  !
!  | --> etan    |(N_CELL)   | free surface            cell-center |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!  | --> Hprv    |(N_VERT)   | Total depth = h + etan  vertices    |  !
!  | --> hv      |(N_VERT)   | still depth             vertices    |  !
!  | --> etav    |(N_VERT)   | free surface            vertices    |  !
!  |_____________|___________|_____________________________________|  !
!  | --> rhofv   |(N_VERT,NZ)| Density at the vertices             |  !
!  | --> dwdtB   |(N_VERT)   | d(w)/dt at the bottom boundary      |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
      integer :: pair,ncount,nvg,jv
      real*8 :: sumdlV
      real*8 :: fun_l,fun_g
#     ifdef KeyParallel
      real*8,dimension(:),allocatable :: phiv_global
      real*8,dimension(:),allocatable :: phiv2_global
      real*8,dimension(:),allocatable :: sumdlV_global
#     else
      integer,dimension(:) :: vpair(4)
      real*8,dimension(:),allocatable :: oldphiv
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: BCpvertex'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                              Initialization                         !
!                                                                     !
!*********************************************************************!
!     --------------------------------------------------
!     HORIZONTAL
        DO k=1,NZ-1
          DO nv=1,N_VERT
             IF(nbev(nv).ne.0) THEN
                  if(bcpv(nv) .eq. 0) phiv(nv,k) = 0.
             ENDIF
          ENDDO
        ENDDO
!     -------------------------------------------------
!     Periodic consideration
#     ifdef KeyTESTpBC
       call bc_periodic(phiv,xv,yv,sigv,dsigv,No_vp,nbev,TagBCp)
#    endif

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: BCpvertex'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!*********************************************************************!
!                                                                     !
!                    Periodic BC for the vertex                       !
!                                                                     !
!*********************************************************************!
     SUBROUTINE bc_periodic(phiv,xv,yv,sigv,dsigv,No_vp,nbev,tagBC)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct vertex boundary condition    !
!    of the velocity components. The tag called "tagBC" is used to    !
!    choose between velocity components:                              !
!                          tagBC = 1 is u                             !
!                          tagBC = 2 is v                             !
!                          tagBC = 3 is w                             !
!                          tagBC = 0 is p                             !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !
!  |_____________|_____________|___________________________________|  !
!  | <--- phiv   |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  | ---> xv,yv  |(N_VERT)   | Coordinates of the vertices         |  !
!  | ---> sigv   |(NZ-1)     | sigma value at the vertices         |  !
!  | ---> dsigv  |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | ---> No_vp  |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | ---> nbev   |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8, dimension(:)  :: sigv(NZ-1)
      real*8, dimension(:)  :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
      integer :: tagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
      integer :: pair,ncount,j0,nvg,flag,kk
      real*8 :: sumdlV
#     ifdef KeyParallel
      real*8,dimension(:),allocatable :: phiv2_gl
      real*8,dimension(:),allocatable :: phiv_gl
      real*8,dimension(:),allocatable :: sumdlV_gl
#     else
      integer,dimension(:) :: vpair(1:4)
      real*8,dimension(:),allocatable :: oldphiv
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: bc_periodic'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyTESTpBC
!     ========================================
!     =====          SERIAL MODE    ==========
!     ========================================
#     ifndef KeyParallel
!     ------------------------
!     allocate variable
      allocate(oldphiv(N_VERT))
!     ------------------------
!     start vertical loop
      DO k=1,NZ-1
            ncount = 0
!     -----------------------
!     assign initial value
            do nv=1,N_VERT
                oldphiv(nv) = 0.0d0
                if(vertex_edge(nv) .eq. 1) then
                    ncount = ncount +1
!                    print*, ncount
                    vpair(ncount) = nv
                endif
            enddo
!    --------------------------------------------------------------
          if(XPB*YPB .eq. 1) then
             if(ncount .ne. 4) print*, 'ERROR! EDGE POINT MISSING'
          endif
!    --------------------------------------------------------------
!     exchange periodic info for boundary vertex
            do nv =1,N_VERT
                if(nbev(nv) .ne. 0) then
                    if(vertex_edge(nv) .eq. 1) then
                        sumdlV = 0.0d0
                        do i=1,ncount
                            pair = vpair(i)
                            oldphiv(nv) = oldphiv(nv) + phiv(pair,k)*dlVsum(pair)
                            sumdlV = sumdlV+dlVsum(pair)
                        enddo
                        oldphiv(nv) = oldphiv(nv)/sumdlV
                    elseif(vertex_pair(nv) .gt. 0) then
                        pair = vertex_pair(nv)
                        sumdlV = dlVsum(nv)+dlVsum(pair)
                        oldphiv(nv) = dlVsum(nv)*phiv(nv,k)+dlVsum(pair)*phiv(pair,k)
                        oldphiv(nv) = oldphiv(nv)/sumdlV
                    endif
                endif
            enddo
!   -------------------------------------------------------------------
!    update local periodic vertex info
            do nv=1,N_VERT
                if(nbev(nv) .ne. 0) then
                    if(tagBC .eq. TagBCu) flag = bcuv(nv)
                    if(tagBC .eq. TagBCv) flag = bcvv(nv)
                    if(tagBC .eq. TagBCw) flag = bcwv(nv)
                    if(tagBC .eq. TagBCp) flag = bcpv(nv)
                    if(flag .eq. 3) then
                        phiv(nv,k) = oldphiv(nv)
                    endif
                endif
            enddo
!     ----------------
!     end loop vertical
      ENDDO
!     ---------------------------------------
!     deallocate variable
       deallocate(oldphiv)
!     ========================================
!     =====          PARALLEL MODE  ==========
!     ========================================
#    else
!     ---------------------------------
!     allocate variable
       allocate(phiv_gl(N_VERTglobal),   &
                sumdlV_gl(N_VERTglobal), &
                phiv2_gl(N_VERTglobal))
!     ----------------------------------
!     Begin vertical direction loop
    DO k=1,NZ-1
        do kk=1,NZBlock
!     ----------------------------------
!     Initial global value
            do nv=1,N_VERTglobal
                phiv_gl(nv) = 0.0d0
                phiv2_gl(nv) = 0.0d0
                sumdlV_gl(nv) = 0.0d0
            enddo

            if(NZLayer .eq. kk) then
!    -----------------------------------
!    Assign local value to global
                do nv=1,N_VERT
                    if(nbev(nv) .ne. 0) then
                        nvg = index_globalv(nv)
                        if(vpair_global(nvg) .gt. 0) then
                            phiv_gl(nvg) = phiv(nv,k)*dlVsum(nv)/vcount_global(nvg)
                            sumdlV_gl(nvg) = dlVsum(nv)/vcount_global(nvg)
                        endif
                        if(vedge_global(nvg) .eq. 1) then
                            phiv_gl(nvg) = phiv(nv,k)*dlVsum(nv)/vcount_global(nvg)
                            sumdlV_gl(nvg) = dlVsum(nv)/vcount_global(nvg)
                        endif
                    endif
                enddo
            endif
!    -----------------------------------
!    sum across all procs
                call MPI_Barrier(comm3D,code)
                call SUM_Parallelr(phiv_gl,N_VERTglobal)
                call SUM_Parallelr(sumdlV_gl,N_VERTglobal)
!    -----------------------------------
!       update global info
                do nv=1,N_VERTglobal
                    if(vedge_global(nv) .eq. 1) then
                        sumdlV = 0.0d0
                        do i=1,N_VERTglobal
                            if(vedge_global(i) .eq. 1) then
                                if(i .ne. nv) then
                                    phiv2_gl(nv) = phiv2_gl(nv)+ phiv_gl(i)
                                    sumdlV = sumdlV+sumdlV_gl(i)
                                endif
                            endif
                        enddo
                        phiv2_gl(nv) = phiv2_gl(nv) + phiv_gl(nv)
                        sumdlV = sumdlV + sumdlV_gl(nv)
                        phiv2_gl(nv) = phiv2_gl(nv)/sumdlV
                    elseif(vpair_global(nv) .gt. 0) then
                        pair = vpair_global(nv)
                        sumdlV = sumdlV_gl(nv)+sumdlV_gl(pair)
                        phiv2_gl(nv) = phiv_gl(nv) + phiv_gl(pair)
                        phiv2_gl(nv) = phiv2_gl(nv)/sumdlV
                    endif
                enddo
!    -----------------------------------
!     update local info
            if(NZLayer .eq. kk) then
                do nv=1,N_VERT
                    if(nbev(nv) .ne. 0) then
                        if(tagBC .eq. TagBCu) flag = bcuv(nv)
                        if(tagBC .eq. TagBCv) flag = bcvv(nv)
                        if(tagBC .eq. TagBCw) flag = bcwv(nv)
                        if(tagBC .eq. TagBCp) flag = bcpv(nv)
                        if(flag .eq. 3) then
                            nvg = index_globalv(nv)
                            phiv(nv,k) = phiv2_gl(nvg)
                        endif
                    endif
                enddo
            endif
        enddo
!     ----------------------------------
!     end vertical loop
    ENDDO
!    -----------------------------------------------
!    deallocate variable
        deallocate(phiv_gl,phiv2_gl,sumdlV_gl)
#    endif
!     ========================================
!     =====          PARALLEL MODE  ==========
!     ========================================
#    endif
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: bc_periodic'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
