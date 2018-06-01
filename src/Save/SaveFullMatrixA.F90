!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!  WRITE THE FULL MATRIX A WITH ALL RELATIONS CELL-CENTER & VERTEX    !
!                            May 2015                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE SaveFullMatrixA(No_cp,nbe,No_vp,nbev)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine find all the structure of the matrix A resulting !
!    from the linear system. The relations are given by the cell-     !
!    center conexions and vertex interpolation. The files correspon-  !
!    ding to the 2D and 3D matrix are saved in the files:             !
!             output/Matlab/PlotMatrixA/MatrixA2D.dat'                !
!             output/Matlab/PlotMatrixA/MatrixA2D.dat'                !
!    respectively.                                                    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> No_vp | (N_CELL0,3) | Numbering of the 3 cell vertices   |  !
!  | ---> nbev  | (N_VERT)    | Tag type of cell vertex            |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  |    Am1     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  |    AmT     |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  |    AmB     |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  !
!  |____________|____________|_____________________________________|  !
!  |    Bmv1T   |(N_CELL0,NZ)| matrix coeff. vertex 1 top          |  ! 
!  |    Bmv2T   |(N_CELL0,NZ)| matrix coeff. vertex 2 top          |  ! 
!  |    Bmv3T   |(N_CELL0,NZ)| matrix coeff. vertex 3 top          |  ! 
!  |    Bmv1B   |(N_CELL0,NZ)| matrix coeff. vertex 1 bottom       |  ! 
!  |    Bmv2B   |(N_CELL0,NZ)| matrix coeff. vertex 2 bottom       |  ! 
!  |    Bmv3B   |(N_CELL0,NZ)| matrix coeff. vertex 3 bottom       |  !
!  |____________|____________|_____________________________________|  !  
!  |    Newrhs  |(N_CELL0,NZ)| right hand side of the method       |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - diffusion3D                ( diffusion3D.F90 )            |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - BCcellcenter3D             ( BCcellcenter3D.F90 )         |  !
!  |   - BCvertex3D                 ( BCvertex3D.F90 )             |  !
!  |_______________________________________________________________|  !
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

!     -------------------------------------
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT) 
!     --------------------------------------
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     ----------------------------------------
      integer :: i1,i2,i3,j0,j1,j2,j3,jT,jB,irec
!     ---------------------------------------- 
      integer :: ii,jj,jjT,jjB,ivert
      integer :: Level,maximo,MaxCol2D,MaxCol3D
      integer :: Njv1,Njv2,Njv3,TotalNN

      integer,dimension(:,:),allocatable :: MatrixA2D
      integer,dimension(:,:),allocatable :: MatrixA3D
!     ----------------------------------------

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: FullMatrixA'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

#     ifndef KeyParallel
!     ________________________________________________________
!    |                                                        |
!    |                Dimensions of Matrix A                  |
!    |________________________________________________________|

!     ________________________________________________________
!     Maximum dimension of index for all vertices          

      TotalNN = NZ*N_CELL0
      maximo = 1
      do nv=1,N_VERT
         maximo = max(maximo,Dimsurrounding(nv))
      enddo
      MaxCol2D = 4 + 3*maximo
      MaxCol3D = 6 + 9*maximo 
!     ________________________________________________________
!     allocate         

      allocate(MatrixA2D(N_CELL0,MaxCol2D))
      allocate(MatrixA3D(TotalNN,MaxCol3D))

!     ________________________________________________________
!     Inital index         

      do i=1,N_CELL0 
         do j=1,MaxCol2D
            MatrixA2D(i,j) = 0
         enddo
      enddo
      do k=1,NZ
         do i=1,N_CELL0 
            nc = i + (k-1)*N_CELL0
            do j=1,MaxCol3D
               MatrixA3D(nc,j) = 0
            enddo
         enddo
      enddo
!     ________________________________________________________
!    |                                                        |
!    |                     Fill Matrix A (2D)                 |
!    |________________________________________________________|

      do i=1,N_CELL0 
         !--------------------------------------
         ! Cell-center neighbors
         i1 = No_cp(i,1)
         i2 = No_cp(i,2)
         i3 = No_cp(i,3)	
         if (i1.gt.N_CELL0)  i1 = i
         if (i2.gt.N_CELL0)  i2 = i
         if (i3.gt.N_CELL0)  i3 = i
         MatrixA2D(i,1) = i
         MatrixA2D(i,2) = i1
         MatrixA2D(i,3) = i2
         MatrixA2D(i,4) = i3
         !--------------------------------------
         ! Cell-center of vertex values
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
	 do j=1,Dimsurrounding(jv1)
            jj  = j + 4  
            j0 = surrounding(jv1,j)
            if (j0.gt.N_CELL0)  j0 = i
            MatrixA2D(i,jj) = j0                     
         enddo
	 do j=1,Dimsurrounding(jv2)
            jj  = jj + j
            j0 = surrounding(jv2,j) 
            if (j0.gt.N_CELL0)  j0 = i
            MatrixA2D(i,jj) = j0                    
         enddo
	 do j=1,Dimsurrounding(jv3)
            jj  = jj + j
            j0 = surrounding(jv3,j) 
            if (j0.gt.N_CELL0)  j0 = i
            MatrixA2D(i,jj) = j0               
         enddo
      enddo
!     ________________________________________________________
!    |                                                        |
!    |                   Fill Matrix A (3D)                   |
!    |________________________________________________________|

      do k=1,NZ
         Level = (k-1)*N_CELL0
         do i=1,N_CELL0 
            nc = i  + Level 
            !--------------------------------------
            ! Cell-center neighbors
            i1 = No_cp(i,1)
            i2 = No_cp(i,2)
            i3 = No_cp(i,3)	
            j1 = i1 + Level
            j2 = i2 + Level
            j3 = i3 + Level
            jT = nc + N_CELL0
            jB = nc - N_CELL0
            if (i1.gt.N_CELL0)  j1 = nc
            if (i2.gt.N_CELL0)  j2 = nc
            if (i3.gt.N_CELL0)  j3 = nc
            if (k.eq.1)         jB = nc 
            if (k.eq.NZ)        jT = nc
            MatrixA3D(nc,1) = nc
            MatrixA3D(nc,2) = j1
            MatrixA3D(nc,3) = j2
            MatrixA3D(nc,4) = j3
            MatrixA3D(nc,5) = jT
            MatrixA3D(nc,6) = jB
            !--------------------------------------
            ! Cell-center of vertex values
            jv1 = No_vp(i,1)
            Njv1 = Dimsurrounding(jv1)
	    do j=1,Njv1
               jj  = j   + 6  
               jjT = jj  + Njv1
               jjB = jjT + Njv1
               j0 = surrounding(jv1,j) + Level
               jT = j0 + N_CELL0
               jB = j0 - N_CELL0
               if (j0.gt.N_CELL0)  j0 = nc
               if (k.eq.1)         jB = nc 
               if (k.eq.NZ)        jT = nc
               MatrixA3D(nc,jj)  = j0               
               MatrixA3D(nc,jjT) = jT   
               MatrixA3D(nc,jjB) = jB             
            enddo
            jv2 = No_vp(i,2)
            Njv2 = Dimsurrounding(jv2)
	    do j=1,Njv2
               jj  = j   + 6 + 3*Njv1 
               jjT = jj  + Njv2
               jjB = jjT + Njv2
               j0 = surrounding(jv2,j) + Level
               jT = j0 + N_CELL0
               jB = j0 - N_CELL0
               if (j0.gt.N_CELL0)  j0 = nc
               if (k.eq.1)         jB = nc 
               if (k.eq.NZ)        jT = nc
               MatrixA3D(nc,jj)  = j0               
               MatrixA3D(nc,jjT) = jT   
               MatrixA3D(nc,jjB) = jB            
            enddo
            jv3 = No_vp(i,3)
            Njv3 = Dimsurrounding(jv3)
	    do j=1,Njv3
               jj  = j   + 6 + 3*Njv1 + 3*Njv2
               jjT = jj  + Njv3
               jjB = jjT + Njv3
               j0 = surrounding(jv3,j) + Level
               jT = j0 + N_CELL0
               jB = j0 - N_CELL0
               if (j0.gt.N_CELL0)  j0 = nc
               if (k.eq.1)         jB = nc 
               if (k.eq.NZ)        jT = nc
               MatrixA3D(nc,jj)  = j0               
               MatrixA3D(nc,jjT) = jT   
               MatrixA3D(nc,jjB) = jB               
            enddo
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                   Save in data files                   |
!     |________________________________________________________|
 
      print*,'  '
      print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
      print*,'  '
      print*,' Saved: output/Matlab/PlotMatrixA/MatrixA2D.dat'
      print*,' Saved: output/Matlab/PlotMatrixA/MatrixA3D.dat'
      print*,'  '
      print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'  '
      open(60,file='../output/Matlab/PlotMatrixA/MatrixA2D.dat')
      open(70,file='../output/Matlab/PlotMatrixA/MatrixA3D.dat')
      do i=1,N_CELL0 
            write(60,*) (MatrixA2D(i,j),j=1,MaxCol2D)
      enddo
      do k=1,NZ
         do i=1,N_CELL0 
            nc = i + (k-1)*N_CELL0
            write(70,*) (MatrixA3D(nc,j),j=1,MaxCol3D)
         enddo
      enddo
      close(60)
      close(70)

#     endif


!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      deallocate(MatrixA2D)
      deallocate(MatrixA3D)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: FullMatrixA'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                          End of FullMatrixA                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |             Save matrix A as a data file               |
!     |________________________________________________________|

#     ifdef KeySaveMatrixA
#     ifndef KeyParallel
      print*,'  '
      print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
      print*,'  '
      print*,' Saved: output/Matlab/PlotMatrixA/MatrixA.dat'
      print*,'  '
      print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'  '

3     format(6(1x,i6,1x,e12.5))
      irec = 60
      open(irec,file='../output/Matlab/PlotMatrixA/MatrixA.dat')
      do k=1,NZ
         do i=1,N_CELL0 
            if (No_cp(i,1).gt.N_CELL0) then
               j1 = i
            else
               j1 = No_cp(i,1)
            endif
            if (No_cp(i,2).gt.N_CELL0) then
               j2 = i
            else
               j2 = No_cp(i,2)
            endif
            if (No_cp(i,3).gt.N_CELL0) then
               j3 = i
            else
               j3 = No_cp(i,3)
            endif
            j0 = i  + (k-1)*N_CELL0 	
            j1 = j1 + (k-1)*N_CELL0
            j2 = j2 + (k-1)*N_CELL0
            j3 = j3 + (k-1)*N_CELL0
            jT = i + k*N_CELL0
            jB = i + (k-2)*N_CELL0
            if (k.eq.1) then
               write(irec,3) j0,1.0d0,&
                             j1,0.0d0,&
                             j2,0.0d0,&
                             j3,0.0d0,&
                             jT,0.0d0,&
                             j0,1.0d0   !<----- 
            elseif (k.eq.NZ) then
               write(irec,3) j0,1.0d0,&
                             j1,0.0d0,&
                             j2,0.0d0,&
                             j3,0.0d0,&
                             j0,1.0d0,& !<-----
                             jB,0.0d0
            else
               write(irec,3) j0,Am0(i,k),&
                             j1,Am1(i,k),&
                             j2,Am2(i,k),&
                             j3,Am3(i,k),&
                             jT,AmT(i,k),&
                             jB,AmB(i,k)
            endif
           enddo
        enddo
      close(irec)
#     endif
#     endif
