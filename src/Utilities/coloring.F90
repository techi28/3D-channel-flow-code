!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          COLORING CELLS                             !
!                             Sept 2014                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE coloring(No_cp,No_vp,xv,yv)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the corresponding color of each cell.    !
!    The general rule is that two neighbor cells can not have the     !
!    same color. Aditionally we try to get the minimum number of      !
!    colors.                                                          !
!---------------------------------------------------------------------!
!                                                                     !
!    Common parameters & variables modified:                          !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | * ColorCell | N_CELL0   | color of each cell: 1, 2, 3 or 4    |  !
!  | * N_COLOR   | integer   | Total number of colors used         |  ! 
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <--> No_cp  |(N_CELL,3) | Node No. of surrounding three cell  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!    *   Common variables modified                                    !
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

      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      real*8,dimension(:):: xv(N_VERT)
      real*8,dimension(:):: yv(N_VERT)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer,dimension(:,:), allocatable:: neighborsCC
      integer,dimension(:),   allocatable:: ColorCellaux
      integer :: N0,N1,N2,N3,i1,i2,i3,s,jj
      integer :: nc1,nc2,nc3,c1,c2,c3,c,ielem
      integer :: nv1,nv2,nv3
      integer :: TypeColoring  !=1 (Foward), =2 (Backward)
      TypeColoring = 2

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: coloring'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                        Coloring parallel                            !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(neighborsCC(N_CELL0global,3))
      allocate(ColorCellaux(N_CELL0global))

!      ________________________________________________________
!     |                                                        |
!     |                 Calling the neighbors                  |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         do i=1,N_CELL0global
            do j=1,3
               neighborsCC(i,j) = No_cp_global(i,j)
            enddo
         enddo 
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
         do i=1,N_CELL0global
            do j=1,3
               neighborsCC(i,j) = No_cp(i,j)
            enddo
         enddo
#     endif
!     =============== END ================    
!     ==================================== 

!      ________________________________________________________
!     |                                                        |
!     |             Initial palette: -1 (no color )            |
!     |________________________________________________________|

      do i=1,N_CELL0global
         ColorCellaux(i) = -1
      enddo
      N0 = 0
      N1 = 0
      N2 = 0
      N3 = 0
!      ________________________________________________________
!     |                                                        |
!     |               Considering 4 posible cases              |
!     |________________________________________________________|

         if (TypeColoring.eq.2) then
            ColorCellaux(N_CELL0global) = 0
         else
            ColorCellaux(1) = 0
         endif
         N0 = N0 + 1
         
         DO jj=2,N_CELL0global
         
         !_____________________________________________________
         ! Identify neighbors that has not been colored yet
         if (TypeColoring.eq.2) then
            i = N_CELL0global -(jj-1)
            i1 = neighborsCC(i,1)
            i2 = neighborsCC(i,2)
            i3 = neighborsCC(i,3)
            if (i1.lt.i) i1 = 0
            if (i2.lt.i) i2 = 0
            if (i3.lt.i) i3 = 0
         else
            i = jj
            i1 = neighborsCC(i,1)
            i2 = neighborsCC(i,2)
            i3 = neighborsCC(i,3)
            if (i1.gt.i) i1 = 0
            if (i2.gt.i) i2 = 0
            if (i3.gt.i) i3 = 0
         endif     	
         !_______________________________________________________
         ! CASE 1: All neighbors have no color
         IF ((i1.eq.0).and.(i2.eq.0).and.(i3.eq.0)) THEN
            ColorCellaux(i) = 0
         !_______________________________________________________
         ! CASE 2: One neighbor has been colored
         ELSEIF ((i1.gt.0).and.(i2.eq.0).and.(i3.eq.0)) THEN
            c = ColorCellaux(i1)
            if     (c.eq.0) then 
               ColorCellaux(i)=1
            elseif (c.eq.1) then
               ColorCellaux(i)=0
            elseif (c.eq.2) then
               ColorCellaux(i)=0
            elseif (c.eq.3) then
               ColorCellaux(i)=0
            endif  
         ELSEIF ((i1.eq.0).and.(i2.gt.0).and.(i3.eq.0)) THEN
            c = ColorCellaux(i2)
            if     (c.eq.0) then
               ColorCellaux(i)=1
            elseif (c.eq.1) then
               ColorCellaux(i)=0
            elseif (c.eq.2) then
               ColorCellaux(i)=0
            elseif (c.eq.3) then
               ColorCellaux(i)=0
            endif
         ELSEIF ((i1.eq.0).and.(i2.eq.0).and.(i3.gt.0)) THEN   
            c = ColorCellaux(i3)
            if     (c.eq.0) then
               ColorCellaux(i)=1
            elseif (c.eq.1) then
               ColorCellaux(i)=0 
            elseif (c.eq.2) then
               ColorCellaux(i)=0
            elseif (c.eq.3) then
               ColorCellaux(i)=0
            endif
         !_______________________________________________________
         ! CASE 3: Two neighbors have been colored
         ELSEIF (((i1.eq.0).and.(i2.gt.0).and.(i3.gt.0)).or. &
                ((i1.gt.0).and.(i2.eq.0).and.(i3.gt.0)).or. &
                ((i1.gt.0).and.(i2.gt.0).and.(i3.eq.0))) THEN

            !------------------
            if (i1.eq.0) then
               nc1 = i2
               nc2 = i3
            elseif (i2.eq.0) then
               nc1 = i1
               nc2 = i3
            elseif (i3.eq.0) then
               nc1 = i1
               nc2 = i2
            endif
            !------------------
            c1 = ColorCellaux(nc1);
            c2 = ColorCellaux(nc2);
            if (c1.eq.c2) then
               if (c1.eq.0) then 
                  ColorCellaux(i)=1
               elseif (c1.eq.1) then
                  ColorCellaux(i)=0   
               elseif (c1.eq.2) then
                  ColorCellaux(i)=0  
               elseif (c1.eq.3) then
                  ColorCellaux(i)=0
               endif
            else
               if     ((c1.eq.0).and.(c2.eq.1)) then
                  ColorCellaux(i) = 2
               elseif ((c1.eq.0).and.(c2.eq.2)) then
                  ColorCellaux(i) = 1
               elseif ((c1.eq.0).and.(c2.eq.3)) then
                  ColorCellaux(i) = 1
               elseif ((c1.eq.1).and.(c2.eq.0)) then
                  ColorCellaux(i) = 2
               elseif ((c1.eq.1).and.(c2.eq.2)) then
                  ColorCellaux(i) = 0
               elseif ((c1.eq.1).and.(c2.eq.3)) then
                  ColorCellaux(i) = 0
               elseif ((c1.eq.2).and.(c2.eq.0)) then
                  ColorCellaux(i) = 1
               elseif ((c1.eq.2).and.(c2.eq.1)) then
                  ColorCellaux(i) = 0
               elseif ((c1.eq.2).and.(c2.eq.3)) then
                  ColorCellaux(i) = 0
               elseif ((c1.eq.3).and.(c2.eq.0)) then
                  ColorCellaux(i) = 1
               elseif ((c1.eq.3).and.(c2.eq.1)) then
                  ColorCellaux(i) = 0
               elseif ((c1.eq.3).and.(c2.eq.2)) then
                  ColorCellaux(i) = 0   
               endif 
            endif
         !_______________________________________________________
         ! CASE 4: All neighbors have been colored        
         ELSE
            c1 = ColorCellaux(i1)
            c2 = ColorCellaux(i2)
            c3 = ColorCellaux(i3)

            if (c1.eq.c2) then
               !-----------------
               ! CASE: c1=c2=c3
               if (c1.eq.c3) then
                  if     (c1.eq.0) then
                     ColorCellaux(i) = 1
                  elseif (c1.eq.1) then
                     ColorCellaux(i) = 0
                  elseif (c1.eq.2) then
                     ColorCellaux(i) = 0
                  elseif (c1.eq.3) then
                     ColorCellaux(i) = 0
                  endif
               !-----------------
               ! CASE: (c1=c2)=~c3
               else
                  if     ((c1.eq.0).and.(c3.eq.1)) then
                     ColorCellaux(i) = 2
                  elseif ((c1.eq.0).and.(c3.eq.2)) then
                     ColorCellaux(i) = 1
                  elseif ((c1.eq.0).and.(c3.eq.3)) then
                     ColorCellaux(i) = 1
                  elseif ((c1.eq.1).and.(c3.eq.0)) then
                     ColorCellaux(i) = 2
                  elseif ((c1.eq.1).and.(c3.eq.2)) then
                     ColorCellaux(i) = 0
                  elseif ((c1.eq.1).and.(c3.eq.3)) then
                     ColorCellaux(i) = 0
                  elseif ((c1.eq.2).and.(c3.eq.0)) then
                     ColorCellaux(i) = 1
                  elseif ((c1.eq.2).and.(c3.eq.1)) then
                     ColorCellaux(i) = 0
                  elseif ((c1.eq.2).and.(c3.eq.3)) then
                     ColorCellaux(i) = 0
                  elseif ((c1.eq.3).and.(c3.eq.0)) then
                     ColorCellaux(i) = 1
                  elseif ((c1.eq.3).and.(c3.eq.1)) then
                     ColorCellaux(i) = 0
                  elseif ((c1.eq.3).and.(c3.eq.2)) then
                     ColorCellaux(i) = 0
                  endif
               endif
            else
               !-----------------
               ! CASE: c2=~(c1=c3)
               if (c1.eq.c3) then
                  if     ((c1.eq.0).and.(c2.eq.1)) then
                      ColorCellaux(i) = 2
                  elseif ((c1.eq.0).and.(c2.eq.2)) then
                      ColorCellaux(i) = 1
                  elseif ((c1.eq.0).and.(c2.eq.3)) then
                      ColorCellaux(i) = 1
                  elseif ((c1.eq.1).and.(c2.eq.0)) then
                      ColorCellaux(i) = 2
                  elseif ((c1.eq.1).and.(c2.eq.2)) then
                      ColorCellaux(i) = 0
                  elseif ((c1.eq.1).and.(c2.eq.3)) then
                      ColorCellaux(i) = 0
                  elseif ((c1.eq.2).and.(c2.eq.0)) then
                      ColorCellaux(i) = 1
                  elseif ((c1.eq.2).and.(c2.eq.1)) then
                      ColorCellaux(i) = 0
                  elseif ((c1.eq.2).and.(c2.eq.3)) then
                      ColorCellaux(i) = 0
                  elseif ((c1.eq.3).and.(c2.eq.0)) then
                      ColorCellaux(i) = 1
                  elseif ((c1.eq.3).and.(c2.eq.1)) then
                      ColorCellaux(i) = 0
                  elseif ((c1.eq.3).and.(c2.eq.2)) then
                      ColorCellaux(i) = 0
                  endif
               else
                  !-----------------
                  ! CASE: c1=~(c2=c3)
                  if (c2.eq.c3) then
                      if     ((c1.eq.0).and.(c2.eq.1)) then
                          ColorCellaux(i) = 2
                      elseif ((c1.eq.0).and.(c2.eq.2)) then
                          ColorCellaux(i) = 1
                      elseif ((c1.eq.0).and.(c2.eq.3)) then
                          ColorCellaux(i) = 1
                      elseif ((c1.eq.1).and.(c2.eq.0)) then
                          ColorCellaux(i) = 2
                      elseif ((c1.eq.1).and.(c2.eq.2)) then
                          ColorCellaux(i) = 0
                      elseif ((c1.eq.1).and.(c2.eq.3)) then
                          ColorCellaux(i) = 0
                      elseif ((c1.eq.2).and.(c2.eq.0)) then
                          ColorCellaux(i) = 1
                      elseif ((c1.eq.2).and.(c2.eq.1)) then
                          ColorCellaux(i) = 0
                      elseif ((c1.eq.2).and.(c2.eq.3)) then
                          ColorCellaux(i) = 0
                      elseif ((c1.eq.3).and.(c2.eq.0)) then
                          ColorCellaux(i) = 1
                      elseif ((c1.eq.3).and.(c2.eq.1)) then
                          ColorCellaux(i) = 0
                      elseif ((c1.eq.3).and.(c2.eq.2)) then
                          ColorCellaux(i) = 0
                      endif
                  !-----------------
                  ! CASE: c1=~c2=~c3
                  else
                     ColorCellaux(i) = 6 -(c1+c2+c3);                         
                  endif
               endif        
            endif
         ENDIF
         !_______________________________
         ! Number of each colors
         c = ColorCellaux(i)
         if (c.eq.-1) then
            print*,'Error coloring!!!! No color assigned',i
            print*,i1,ColorCellaux(i1)
            print*,i2,ColorCellaux(i2)
            print*,i3,ColorCellaux(i3)
            stop
         endif
         if (c.eq.0) N0 = N0 + 1
         if (c.eq.1) N1 = N1 + 1
         if (c.eq.2) N2 = N2 + 1
         if (c.eq.3) N3 = N3 + 1
      ENDDO

!      ________________________________________________________
!     |                                                        |
!     |                  Total number of colors                |
!     |________________________________________________________|

      if (N3.gt.0) then
          N_COLOR = 4
      elseif (N2.gt.0) then 
          N_COLOR = 3
      else
          N_COLOR = 2
      endif

!      ________________________________________________________
!     |                                                        |
!     |                 Transfering information                |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         do i=1,N_CELL0
            ielem = Index_global(i)
            if (ielem.le.0) then
               ColorCell(i) = 0
            else
               ColorCell(i) = ColorCellaux(ielem)
            endif
         enddo
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
         do i=1,N_CELL0
            ColorCell(i) = ColorCellaux(i)
         enddo
#     endif
!     =============== END ================    
!     ==================================== 

         do i=1,N_CELL0
            if (ColorCell(i).eq.3) then
               !print*,i
               !print*,xv(No_vp(i,1)),yv(No_vp(i,1))
               !print*,xv(No_vp(i,2)),yv(No_vp(i,2))
               !print*,xv(No_vp(i,3)),yv(No_vp(i,3))
            endif
         enddo
                  
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
   
!      ________________________________________________________
!     |                                                        |
!     |                        Displaying                      |
!     |________________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      if (rang_topo.eq.1) then
      write(*,*) '                                                        '
      write(*,*) '   ==================================================== '
      write(*,*) '                       MPI: COLORING                    '
      write(*,*) '   ==================================================== '
      write(*,*) '                                                        '
      write(*,*) '       -------------------------------------------'
      write(*,*) '       Total number of cells : ',N_CELL0global
      write(*,*) '       -------------------------------------------'
      write(*,*) '       Number of color 1 (blue)  : ',N0
      write(*,*) '       Number of color 2 (red )  : ',N1
      write(*,*) '       Number of color 3 (green) : ',N2
      write(*,*) '       Number of color 4 (yellow): ',N3
      write(*,*) '       -------------------------------------------'
      write(*,*) '       Total number of colors : ',N_COLOR
      write(*,*) '       -------------------------------------------'
      endif
!     ====================================
!     ==========  SEQUENTIAL =============
#     else
      write(*,*) '    ___________________________________________________  '
      write(*,*) '   |                                                   | '
      write(*,*) '   |                      COLORING                     | '
      write(*,*) '   |___________________________________________________| '
      write(*,*) '                                                         '
      write(*,*) '       -------------------------------------------'
      write(*,*) '       Total number of cells : ',N_CELL0global
      write(*,*) '       -------------------------------------------'
      write(*,*) '       Number of color 1 (blue)  : ',N0
      write(*,*) '       Number of color 2 (red )  : ',N1
      write(*,*) '       Number of color 3 (green) : ',N2
      write(*,*) '       Number of color 4 (yellow): ',N3
      write(*,*) '       -------------------------------------------'
      write(*,*) '       Total number of colors : ',N_COLOR
      write(*,*) '       -------------------------------------------'
#     endif
!     =============== END ================    
!     ==================================== 

!      ________________________________________________________
!     |                                                        |
!     |               Saving data (used in Matlab)             |
!     |________________________________________________________|
!                      
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            if (rang_topo.eq.1) then
            print*,'  '
            print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
            print*,'  '
            print*,' Saved: output/Matlab/PlotMultiColors/ColoringData.dat'
            print*,'  '
            print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
            print*,'  '
            open(111,file="../output/Matlab/PlotMultiColors/ColoringData.dat",status='unknown')
            do i=1,N_CELL0global
               nv1 = No_vp_global(i,1)
               nv2 = No_vp_global(i,2)
               nv3 = No_vp_global(i,3)
               write(111,*) ColorCellaux(i),xv_global(nv1),&
                                            xv_global(nv2),&
                                            xv_global(nv3),&
                                            yv_global(nv1),&
                                            yv_global(nv2),&
                                            yv_global(nv3)
            enddo
            close(111)
            endif
!        ====================================
!        ==========  SEQUENTIAL =============
#        else
            print*,'  '
            print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
            print*,'  '
            print*,' Saved: output/Matlab/PlotMulticolors/ColoringData.dat'
            print*,'  '
            print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
            print*,'  '
            open(111,file="../output/Matlab/PlotMultiColors/ColoringData.dat",status='unknown')
            do i=1,N_CELL0
               nv1 = No_vp(i,1)
               nv2 = No_vp(i,2)
               nv3 = No_vp(i,3)
               write(111,2017) ColorCellaux(i),xv(nv1),xv(nv2),xv(nv3),&
                                               yv(nv1),yv(nv2),yv(nv3)
            enddo
            close(111)
#        endif
!        =============== END ================    
!        ====================================     

      2017 format(i3,6(1x,e10.3))

!      ________________________________________________________
!     |                                                        |
!     |                       Deallocate                       |
!     |________________________________________________________|
   
      deallocate(neighborsCC,ColorCellaux)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: coloring'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	      END OF GEOMETRY                         !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
