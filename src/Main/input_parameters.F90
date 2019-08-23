!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  INPUT FILE PARAMETERS & OTHERS                     !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE input_parameters(nbev,hv,h,                     &
                                  No_vp,No_cp,nbe,xv,yv,zbv)

!---------------------------------------------------------------------!
!                                                                     !
!    This program reads the data files:                               !
!                  -   dataParameters.txt                             !
!    corresponding to the input parameters of the problem. This data  !
!    can be changed without need to compile the program again. We     !
!    also perform some operations with this data:                     !
!        1)   Tagging the BC:    -  nbev                              !
!        2)   Depth:             -  hv                                !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !
!  |_____________|___________|_____________________________________|  !
!  |             |           |                                     |  !
!  | --> No_vp   |(N_CELL0,3)| Index number of the 3 cell vertices |  !
!  | --> No_cp   |(N_CELL,3) | Index number of surrounding 3 cells |  !
!  | --> nbe     | N_CELL0   | Type of boundary cell (inside or bc)|  !
!  | --> xv      | N_VERT    | x-coordinate of the vertex          |  !
!  | --> yv      | N_VERT    | y-coordinate of the vertex          |  !
!  | --> zbv     | N_VERT    | The bottom of the river at vertex   |  !
!  |_____________|___________|_____________________________________|  !
!  |             |           |                                     |  !
!  | <-- nbev    | N_VERT    | Type of tag about the kind of vertex|  !
!  | <-- h       | N_CELL    | Depth of the river at cell          |  !
!  | <-- hv      | N_VERT    | Depth of the river at vertex        |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !
!  |______________|________________________________________________|  !
!  |              |                                                |  !
!  |--- N_CELL0   | Number of the cell centers inside the domain   |  !
!  |--- N_CELL    | Total number of the cells                      |  !
!  |--- N_VERT    | Number of the computing vertices               |  !
!  |--- N_CELLmax | Maximum number of cells                        |  !
!  |--- N_VERTmax | maximum Number of vertices                     |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name    |                 Description                       |  !
!  |___________|___________________________________________________|  !
!  |           |                                                   |  !
!  | * IPRESS  | Tag to calculate the pressure                     |  !
!  | * ITURBU  | Tag to calculate the turbulance                   |  !
!  | * ICASE   | Tag to about the release box                      |  !
!  |___________|___________________________________________________|  !
!  | * tFin    | Final simulation time                             |  !
!  | * tIni    | Initial simulation time                           |  !
!  | * dt      | Time step                                         |  !
!  | * dtprint | Step for printing the result output               |  !
!  | * tcomin  | Step for printing                                 |  !
!  |___________|___________________________________________________|  !
!  | ChooseExit|  Choice of the exit format                        |  !
!  |           |  = 1 only tecplot: V-**.tec                       |  !
!  |           |  = 2 tecplot: V-**.tec, C-**.tec, VC-***.tec      |  !
!  |___________|___________________________________________________|  !
!  | * IDEPTH  | Tag to select constant depth                      |  !
!  | * h0      | Constant depth                                    |  !
!  |___________|___________________________________________________|  !
!  | * NmaxCONV| Maximum number of no linear iterations            |  !
!  | * tol     | Tolerance of the linear system                    |  !
!  | * MaxIters| Maximum number of linear system iterations        |  !
!  | * eps     | Tolerance of the linear system                    |  !
!  |___________|___________________________________________________|  !
!  | * NBC     | Number of boundary values                         |  !
!  |___________|___________________________________________________|  !
!  | * ZBV     | Bottom levels at each vertex                      |  !
!  |___________|___________________________________________________|  !
!  | * Iinit0  | Tag to select type of initial condition           |  !
!  | * wl0     | Initial value fot water level                     |  !
!  | * qx0     | Initial value for discharges qx                   |  !
!  | * qy0     | Initial value for discharges qy                   |  !
!  |___________|___________________________________________________|  !
!  | * NC_REF  | Reference cell                                    |  !
!  |___________|___________________________________________________|  !
!  | * nv      | Integer used as the vertex loop counter           |  !
!  | * nc      | Integer used as the cell loop counter             |  !
!  |___________|___________________________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  | IDISPLAY    | Display the data at the windows                 |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!


!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

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
!     ____________________________________
!    |                                    |
!    |     Declaration of variables       |
!    |____________________________________|

      integer,dimension(:)  :: nbev
      real*8 ,dimension(:)  :: h
      real*8 ,dimension(:)  :: hv
!     ------------------------------
      integer,dimension(:,:):: No_vp
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe
      real*8 ,dimension(:)  :: xv
      real*8 ,dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer :: k1,k2,k3,kv1,kv2,kk,nvg
      integer :: count1,count2,count3,count4
      integer:: IDISPLAY = 0
      real*8 :: funh
      character*80 title
      character*80 line
#     ifdef KeyParallel
      integer,dimension(:),allocatable :: ccount_global
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: input_parameters'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Read & write:  input_data.dat                      !
!                                                                     !
!*********************************************************************!

!---------------------------------------------------------------------!
!                                                                     !
!      The file "dataParameters.txt" has the following format:        !
!                                                                     !
!     __________________________________________________________      !
!                             PROBLEM:                                !
!     __________________________________________________________      !
!     Runge-Kutta steps (1 or 2)             | FinRK       : 1        !
!     __________________________________________________________      !
!     Choose the type of save exit file      | ChooseExit : 2         !
!     __________________________________________________________      !
!     Tag to calculate the pressure          | IPRESS  :  1           !
!     Tag to calculate the turbulance        | ITURBU  :  0           !
!     Tag to calculate the input box         | ICASE   :  1           !
!     __________________________________________________________      !
!     Initial time simulation                | tIni    =  3000.0      !
!     Final time simulation                  | tFin    =  3000.0      !
!     Time step                              | dt      =  0.005       !
!     Initial time to save output results    | tInisave=  50.         !
!     Time step to save output results       | dsave   =  50.         !
!     Time step to take stats sample         | dtsample=  50          !
!     Time step to output plane average      | dtplane =  50          !
!     Time limit for the simulation          | tlimit  =  50          !
!     ___________________________________________________________     !
!     Maximum Nolinear iterations (eta)      | NmaxCON =  50          !
!     Tolerance of the nolinear system       | TOL     =  1D-1        !
!     ___________________________________________________________     !
!     LINEAR SOLVER PARAMETERS:                                       !
!     Maximum number of iterations           | MaxIters =  20000      !
!     Tolerance value                        | eps      =  1.0E-05    !
!     SOR relax value                        | relaxSOR =  1.3        !
!     Dim. of Hessemberg matrix (GMRES)      | DimHess  =  20         !
!     Reciprocal of the pesudo-time step     | tau      =  1.0E-04    !
!     Value of the theta-discretization      | theta    =  1.0E-00    !
!     Stop tolerance of the iterations       | epsConv  =  1.0E-06    !
!     Maximum number of iterations           | MaxConv  =  100        !
!     __________________________________________________________      !
!     New system solver parameters:                                   !
!     Reciprocal of the pesudo-time step     | gamma    =  1.0d-2     !
!     Value of the theta-discretization      | theta    =  1.0d0      !
!     Stop tolerance of the iterations       | epsConv  =  1.0d-5     !
!     Maximum number of iterations           | MaxConv  =  100        !
!     __________________________________________________________      !
!     Tag to identify constant depth         | IDEPTH  : 1            !
!     Constant depth                         | h0      = 0.5          !
!     __________________________________________________________      !
!     Initial condition                      | Iinit0 :  1            !
!     Initial value for water level          | WLO    =  9.5          !
!     Initial value for discharges           | QXO    =  0.75         !
!     Initial value for discharges           | QYO    =  0.0          !
!     __________________________________________________________      !
!     REFERENCE CELL                                  :  2881         !
!     __________________________________________________________      !
!                                                                     !
!                                                                     !
!---------------------------------------------------------------------!
#     ifdef KeyParallel
      if(rang_topo .eq. 0) IDISPLAY = 1
#     else
       IDISPLAY = 1
#     endif
      open(13,file='dataParameters.txt',status='OLD')

71    format(A62)
72    format(A25)
73    format(A46,I10)
74    format(A46)
76    format(A46,F15.10)
77    format(A46,E15.5)
78    format(A49,I10)

!      ________________________________________________________
!     |                                                        |
!     |                  Runge-Kutta time steps                |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
!     --------------------------------------------------
!     Read: RK-1 or RK-2
      read(13,78) title,FinRK
      if (IDISPLAY.eq.1) write(*,73) title,FinRK
!      ________________________________________________________
!     |                                                        |
!     |                   Type of saving data                  |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
!     --------------------------------------------------
!     Read: Tag to choose initial simulation
      read(13,78) title,ChooseExit
      if (IDISPLAY.eq.1) write(*,78) title,ChooseExit
!      ________________________________________________________
!     |                                                        |
!     |  Read: Tags to calcule the pressure & solid phase      |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
!     --------------------------------------------------
!     Tag to calculate the pressure
      read(13,73) title,IPRESS
      if (IDISPLAY.eq.1) write(*,73) title,IPRESS
!     --------------------------------------------------
!     Tag to calculate the turbulance
      read(13,73) title,ITURBU
      if (IDISPLAY.eq.1) write(*,73) title,ITURBU
!     --------------------------------------------------
!     Tag about the release box
      read(13,73) title,ICASE
      if (IDISPLAY.eq.1) write(*,73) title,ICASE
!     --------------------------------------------------
!     Tag about the hydrostatic pressure option
      read(13,73) title,IHYDRO
      if (IDISPLAY.eq.1) write(*,73) title,IHYDRO
!      ________________________________________________________
!     |                                                        |
!     |  Read: Time variables                                  |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
!     --------------------------------------------------
!     Read: Initial time
      read(13,76) title,tIni
      if (IDISPLAY.eq.1) write(*,77) title,tIni
!     --------------------------------------------------
!     Read: Final time
      read(13,76) title,tFin
      if (IDISPLAY.eq.1) write(*,77) title,tFin
!     --------------------------------------------------
!     Read: Time Step
      read(13,76) title,dt
      if (IDISPLAY.eq.1) write(*,77) title,dt
!     --------------------------------------------------
!     Read: Time step interval to check CFL
      read(13,73) title,ncfl
      if (IDISPLAY.eq.1) write(*,73) title,ncfl
!     --------------------------------------------------
!     Read: Initial time to start saving results
      read(13,76) title,tInisave
      if (IDISPLAY.eq.1) write(*,77) title,tInisave
!     --------------------------------------------------
!     Read: Time step to save results
      read(13,76) title,dtsave
      if (IDISPLAY.eq.1) write(*,77) title,dtsave
!     --------------------------------------------------
!     Read: Time step to take stats sample
      read(13,73) title,dtsample
      if (IDISPLAY.eq.1) write(*,73) title,dtsample
!     --------------------------------------------------
!     Read: Time step to output averaged stats
      read(13,73) title,dtplane
      if (IDISPLAY.eq.1) write(*,73) title,dtplane
!     --------------------------------------------------
!     Read: Time limit
      read(13,73) title,tlimit
      if (IDISPLAY.eq.1) write(*,73) title,tlimit
!     --------------------------------------------------
!     Read: Reynolds number
      read(13,76) title,Re
      if (IDISPLAY.eq.1) write(*,77) title,Re
80    format(A46,E9.2)
81    format(A47,E9.2)
!      ________________________________________________________
!     |                                                        |
!     |               External body force options              |
!     |________________________________________________________|
      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
      read(13,74) title
      if (IDISPLAY.eq.1) write(*,71) title
!     --------------------------------------------------
!     Read: bdx
      read(13,76) title,bdx
      if (IDISPLAY.eq.1) write(*,77) title,bdx
!     --------------------------------------------------
!     Read: bdy
      read(13,76) title,bdy
      if (IDISPLAY.eq.1) write(*,77) title,bdy
!     --------------------------------------------------
!     Read: bdz
      read(13,76) title,bdz
      if (IDISPLAY.eq.1) write(*,77) title,bdz
!      ________________________________________________________
!     |                                                        |
!     |               Linear system parameters                 |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
      read(13,74) title
      if (IDISPLAY.eq.1) write(*,71) title
!     --------------------------------------------------
!     Read: Maximum number of linear iterations
      read(13,73) title,MaxIters
      if (IDISPLAY.eq.1) write(*,73) title,MaxIters
!     --------------------------------------------------
!     Read: Tolerance of the linear system loop
      read(13,80) title,eps
      if (IDISPLAY.eq.1) write(*,81) title,eps
!     --------------------------------------------------
!     Read: SOR relax value
      read(13,76) title,relaxSOR
      if (IDISPLAY.eq.1) write(*,77) title,relaxSOR
!      ________________________________________________________
!     |                                                        |
!     |  Read: Constant depth                                  |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
!     --------------------------------------------------
!     Read: Simulation duration in second (final time)
      read(13,73) title,IDEPTH
      if (IDISPLAY.eq.1) write(*,73) title,IDEPTH
!     --------------------------------------------------
!     Read: Initial time
      read(13,76) title,h0
      if (IDISPLAY.eq.1) write(*,77) title,h0
!      ________________________________________________________
!     |                                                        |
!     |   Read: Initial condition                              |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
!     --------------------------------------------------
!     Read: Tag initial condition
      read(13,73) title,Iinit0
      if (IDISPLAY.eq.1) write(*,73) title,Iinit0
!     --------------------------------------------------
!     Read: Initial value for water level and discharges
        if (Iinit0.ne.0) then
            read(13,76)  title,wl0
            if (IDISPLAY.eq.1) write(*,77)  title,wl0
            read(13,76)  title,qx0
            if (IDISPLAY.eq.1) write(*,77)  title,qx0
            read(13,76)  title,qy0
            if (IDISPLAY.eq.1) write(*,77)  title,qy0
        else
            wl0 = 0.0d0
            qx0 = 0.0d0
            qy0 = 0.0d0
        endif
!      ________________________________________________________
!     |                                                        |
!     |       Reference cell                                   |
!     |________________________________________________________|

        read(13,71) line
        if (IDISPLAY.eq.1) write(*,71) line
        read(13,73) title,NC_REF
        if (IDISPLAY.eq.1) write(*,73) title,NC_REF
        read(13,71) line
        if (IDISPLAY.eq.1) write(*,71) line
        if (IDISPLAY.eq.1) write(*,*)'      '
!      ________________________________________________________
!     |                                                        |
!     |                         Close                          |
!     |________________________________________________________|

       close(13)


!*********************************************************************!
!                                                                     !
!       Performing and assigning values with the load data            !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |   Tagging the type of vertex (inside or boundary type) |
!     |________________________________________________________|
#     ifdef KeyParallel
!        ----------------------------------------
!        In case of parallel mode
!        Ensure the local vertex bc is the same as global value
         do nv=1,N_VERT
            nvg = index_globalv(nv)
            nbev(nv) = nbev_global(nvg)
         enddo

        allocate(ccount_global(N_CELL0global))

        do i=1,N_CELL0global
             ccount_global(i) = 0
        enddo

        if(NZLayer .eq. 1) then
            do i=1,N_CELL0
                    nc = index_global(i)
                    ccount_global(nc) = ccount_global(nc) +1
            enddo

            do nv=1,N_VERT
                nvg = index_globalv(nv)
                vcount_global(nvg) = vcount_global(nvg) + 1
            enddo
        endif
!     ---------------------------------------------------
!      Obtain the vcount_global for overlapping vertex
        call MPI_Barrier(comm3D,code)
        call SUM_Paralleli(vcount_global,N_VERTglobal)
        call SUM_Paralleli(ccount_global,N_CELL0global)
!     ---------------------------------------------------
!      Check the global count
         IF(IDisplay .eq. 1) THEN
          count1 = 0
          count2 = 0
          count3 = 0
          count4 = 0
          do nv=1,N_VERTglobal
             if(vcount_global(nv) .eq. 1) count1 = count1 +1
             if(vcount_global(nv) .eq. 2) count2 = count2 +1
             if(vcount_global(nv) .eq. 3) count3 = count3 +1
             if(vcount_global(nv) .gt. 3) then
                count4 = count4+1
             endif
          enddo
        print*, '       ==============================='
        print*, '             GLOBAL VERTEX INFO'
        print*, '       ==============================='
        print*, '             TOTAL VERTEX :',N_VERTglobal
        print*, '                 1  PROC  :',count1
        print*, '                 2  PROC  :',count2
        print*, '                 3  PROC  :',count3
        print*, '                >3  PROC  :',count4
        print*, '       ==============================='
!       ------------------------------------------------
!        for center
          count1 = 0
          count2 = 0
          count3 = 0
          do i=1,N_CELL0global
             if(ccount_global(i) .eq. 0) count1 = count1 +1
             if(ccount_global(i) .eq. 1) count2 = count2 +1
             if(ccount_global(i) .gt. 1) count3 = count3 +1
          enddo
        print*, '       ==============================='
        print*, '             GLOBAL CENTER INFO'
        print*, '       ==============================='
        print*, '             TOTAL CENTER :',N_CELL0global
        print*, '                 0  PROC  :',count1
        print*, '                 1  PROC  :',count2
        print*, '                >1  PROC  :',count3
        print*, '       ==============================='

        ENDIF
!      ----------------------------------------------------
!       wait for output to finish
        call MPI_Barrier(comm3D,code)
        deallocate(ccount_global)
#     endif

!      ________________________________________________________
!     |                                                        |
!     |                           Depth                        |
!     |________________________________________________________|

!      __________________________________
!     |                                  |
!     |          At the vertices         |
!     |__________________________________|

!     ------------------------------------
!     Constant depth
      if (iDEPTH.eq.1) then
          do nv=1,N_VERT
             hv(nv) = h0
          enddo
!     ------------------------------------
!     Depth from data
      else
          if (IDISPLAY.eq.1) then
           print*, '      WARNING! DEPTH DATA FROM FUNCTION.'
          endif
	  do nv=1,N_VERT
            hv(nv)= zbv(nv)
            hv(nv) = 0. !<----Function h (depth)
	  enddo
      endif
!      __________________________________
!     |                                  |
!     |        At the cell centers       |
!     |__________________________________|

      do i=1,N_CELL0
         k1=No_vp(i,1)
         k2=No_vp(i,2)
         k3=No_vp(i,3)
         h(i) = (hv(k1)+hv(k2)+hv(k3))/3.0d0
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(h)
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: input_parameters'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF INPUT                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

