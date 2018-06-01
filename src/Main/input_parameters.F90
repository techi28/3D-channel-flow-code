!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  INPUT FILE PARAMETERS & OTHERS                     !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE input_parameters(nbev,hv,h,                     &
                                  No_vp,No_cp,No_wb,No_hb,No_qb, &
                                  nbe,xv,yv,zbv)                                       

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
!  | --> No_wb   | N_WB      | Vertex number of wall boundary      |  !
!  | --> No_qb   | N_QB      | Vertex number of discharge boundary |  !
!  | --> No_hb   | N_HB      | Vertex number of wate surface bound.|  !
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
!  |--- N_WBmax   | Maximum number of N_WB                         |  !
!  |--- N_QBmax   | Maximum number of N_QB                         |  !
!  |--- N_HBmax   | Maximumnumber of  N_HB                         |  !
!  |______________|________________________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name    |                 Description                       |  !  
!  |___________|___________________________________________________|  !   
!  |           |                                                   |  !
!  | * IPRESS  | Tag to calculate the pressure                     |  !
!  | * ISOLID  | Tag to calculate the solid phase                  |  !
!  | * ITURBU  | Tag to calculate the turbulance                   |  !
!  | * ICASE   | Tag to about the release box                      |  !
!  | * IBED    | Tag about the granular bed option                 |  !
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
!  | * N_WB    | Number of the vertices on the wall boundary       |  !
!  | * N_HB    | Number of the vertices on the water level boundary|  !
!  | * N_QB    | Number of the vertices on discharge normal bound. |  !
!  |___________|___________________________________________________|  !
!  | * IRESTART| Tag for the initial simulation                    |  !
!  | filerepout| Name of the output file                           |  ! 
!  | filerepin | Name of the input file                            |  ! 
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
      integer,dimension(:)  :: No_wb 
      integer,dimension(:)  :: No_qb
      integer,dimension(:)  :: No_hb 
      integer,dimension(:)  :: nbe   
      real*8 ,dimension(:)  :: xv
      real*8 ,dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer k1,k2,k3,kv1,kv2,kk,ii
      integer,parameter :: IDISPLAY = 0
      real*8 :: funh
      character*80 title
      character*80 line      

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
!     Tag for initial simulation             | IrestartIN : 0         !
!     Tag for initial simulation             | IrestartOUT: 0         !
!     Name of the restart output file        | restartout             !
!     Name of the restart input file         | restartin              !
!     __________________________________________________________      !
!     Tag to calculate the pressure          | IPRESS  :  1           !
!     Tag to calculate the solid phase       | ISOLID  :  1           !
!     Tag to calculate the turbulance        | ITURBU  :  0           !
!     Tag to calculate the input box         | ICASE   :  1           ! 
!     Tag to calculate the sediment bed      | IBED    :  0           !
!     __________________________________________________________      !
!     Initial time simulation                | tIni    =  3000.0      ! 
!     Final time simulation                  | tFin    =  3000.0      !
!     Time step                              | dt      =  0.005       ! 
!     Initial time to save output results    | tInisave=  50.         !    
!     Time step to save output results       | dsave   =  50.         !
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
!     Time step to save results              | dtsave   =  0.0        !
!     Time step to sample stats              | dtsample =  10         !
!     Time to output averaged results        | dtplane  =  10         !
!     __________________________________________________________      !
!                                                                     !
!                                                                     !
!---------------------------------------------------------------------!

!     ________________________________
#     if defined(KeyEstuaryGironde)
         open(13,file='EstuaryGironde/dataParameters.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStaticCylinder)
         open(13,file='StaticCylinder/dataParameters.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStaticChannel)
         open(13,file='StaticChannel/dataParameters.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStandingWave)
         open(13,file='StandingWave/dataParameters.txt',status='OLD')
!     ________________________________
#     elif defined(KeyTaylorVortex)
         open(13,file='TaylorVortex/dataParameters.txt',status='OLD')
!     ________________________________
#     elif defined(KeyTestOnlyPoisson)
         open(13,file='TestOnlyPoisson/dataParameters.txt',status='OLD')
!     ________________________________
#     else
         open(13,file='dataParameters.txt',status='OLD')
#     endif

71    format(A62)
72    format(A25)
73    format(A46,I10)
74    format(A46)
76    format(A46,F15.10)
77    format(A46,E15.5)
78    format(A49,I10)

!      ________________________________________________________
!     |                                                        |
!     |  Read: Title                                           |
!     |________________________________________________________|

      read(13,71) line
      if (IDISPLAY.eq.1) write(*,71) line
      read(13,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
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
!     |                Re-start time simulation                |
!     |________________________________________________________|

      read(13,71) line	
      if (IDISPLAY.eq.1) write(*,71) line 
!     --------------------------------------------------
!     Read: Tag to choose initial simulation
      read(13,78) title,IrestartIN
      if (IDISPLAY.eq.1) write(*,78) title,IrestartIN
!     --------------------------------------------------
!     Read: Tag to save final simulation
      read(13,78) title,IrestartOUT
      if (IDISPLAY.eq.1) write(*,78) title,IrestartOUT

!     --------------------------------------------------
!     Read:  About the starting file (time) of the simulation
      read(13,'(a36,25a)') title,filerepout
      if (IDISPLAY.eq.1) write(*,'(a36,25a)') title,filerepout
      if (IRESTARTOUT.eq.0) then
          read(13,71) title 
      else
          read(13,'(a36,25a)') title,filerepin
          if (IDISPLAY.eq.1) write(*,'(a36,25a)') title,filerepin
      endif
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
!     Tag to calculate the solid phase   
      read(13,73) title,ISOLID
      if (IDISPLAY.eq.1) write(*,73) title,ISOLID
!     --------------------------------------------------
!     Tag to calculate the turbulance   
      read(13,73) title,ITURBU
      if (IDISPLAY.eq.1) write(*,73) title,ITURBU
!     --------------------------------------------------
!     Tag about the release box   
      read(13,73) title,ICASE
      if (IDISPLAY.eq.1) write(*,73) title,ICASE
!     --------------------------------------------------
!     Tag about the granular bed option   
      read(13,73) title,IBED
      if (IDISPLAY.eq.1) write(*,73) title,IBED
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
!     Read: Initial time to start saving results 
      read(13,76) title,tInisave
      if (IDISPLAY.eq.1) write(*,77) title,tInisave
!     --------------------------------------------------
!     Read: Time step to save results        
      read(13,76) title,dtsave
      if (IDISPLAY.eq.1) write(*,77) title,dtsave

!      ________________________________________________________
!     |                                                        |
!     |          Convergence parameters: Nonlinear loop        |
!     |________________________________________________________|

80    format(A46,E9.2)
81    format(A47,E9.2)

      read(13,71) line	
      if (IDISPLAY.eq.1) write(*,71) line  
!     --------------------------------------------------
!     Read: Maximum number of nonlinear iterations eta 
      read(13,73) title,NmaxCONV
      if (IDISPLAY.eq.1) write(*,73) title,NmaxCONV
!     --------------------------------------------------
!     Read: Tolerance of eta loop                       
      read(13,80) title,tol
      if (IDISPLAY.eq.1) write(*,81) title,tol

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
      !read(13,74) title
      !read(13,*)  eps
      !if (IDISPLAY.eq.1) write(*,77) title,eps
!     --------------------------------------------------
!     Read: SOR relax value
      read(13,76) title,relaxSOR
      if (IDISPLAY.eq.1) write(*,77) title,relaxSOR
!     --------------------------------------------------
!     Read: Heseinberg matrix dimention (GMRES)
      read(13,73) title,DimHess
      if (IDISPLAY.eq.1) write(*,73) title,DimHess
!     --------------------------------------------------
!     Read: Reciprocal of the pesudo-time step                       
      read(13,76) title,ttau
      if (IDISPLAY.eq.1) write(*,77) title,ttau
!     --------------------------------------------------
!     Read: Value of the theta-discretization
      read(13,76) title,ttheta
      if (IDISPLAY.eq.1) write(*,77) title,ttheta
!     --------------------------------------------------
!     Read: Stop tolerance of the iterations
      read(13,80) title,epsConv
      if (IDISPLAY.eq.1) write(*,81) title,epsConv
!     --------------------------------------------------
!     Read: Maximum number of iterations
      read(13,73) title,MaxConv
      if (IDISPLAY.eq.1) write(*,73) title,MaxConv
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

!      ________________________________________________________
!     |                                                        |
!     |   Read: Time Statistics                                |
!     |________________________________________________________|

#     ifdef KeySaveStatistics
!     - - - - - - - - - - - - - - - - - - - - - - - - - 
!     --------------------------------------------------
!     Read: Initial time 
      read(13,76) title,tInistats
      if (IDISPLAY.eq.1) write(*,77) title,tInistats
!     --------------------------------------------------
!     Read: Final time     
      read(13,73) title,dtsample
      if (IDISPLAY.eq.1) write(*,73) title,dtsample
!     --------------------------------------------------
!     Read: Time Step   
      read(13,73) title,dtplane
      if (IDISPLAY.eq.1) write(*,73) title,dtplane
!     --------------------------------------------------
      read(13,71) line	
      if (IDISPLAY.eq.1) write(*,71) line
!     - - - - - - - - - - - - - - - - - - - - - - - - -
#     endif
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

!     --------------------------------------------------
!     Inside vertex                             nbev = 0
      do nv=1,N_VERT
         nbev(nv)=0
      enddo
!     --------------------------------------------------
!     Wall boundary                             nbev = 1
      if (N_WB.ne.0) then
         do i=1,N_WB
            nv=No_wb(i)
            nbev(nv)=1
         enddo
      endif
!     --------------------------------------------------
!     Inflow boundary                           nbev = 2
      if (N_QB.ne.0) then
         do i=1,N_QB
            nv=No_qb(i)
            nbev(nv)=2
         enddo
      endif
!     --------------------------------------------------
!     OutFlow boundary                          nbev = 3
      if (N_HB.ne.0) then
        do i=1,N_HB
           nv=No_hb(i)
           nbev(nv)=3
        enddo
      endif
            
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
	  do nv=1,N_VERT
            hv(nv)= zbv(nv)
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

