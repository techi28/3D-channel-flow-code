!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	     INPUT FILE                               !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE input_data(No_vp,No_cp,No_wb,No_hb,No_qb,No_sp,&
                            nbe,xv,yv,zbv)                                       

!---------------------------------------------------------------------!
!                                                                     !
!    This program reads the data files:                               !
!                  -   input_data.dat                                 !
!    corresponding to the data of the physical domain.                !   
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !   
!  |             |           |                                     |  !
!  | <-- No_vp   |(N_CELL0,3)| Index number of the 3 cell vertices |  !
!  | <-- No_cp   |(N_CELL,3) | Index number of surrounding 3 cells |  !
!  | <-- No_wb   | N_WB      | Vertex number of wall boundary      |  !
!  | <-- No_qb   | N_QB      | Vertex number of discharge boundary |  !
!  | <-- No_hb   | N_HB      | Vertex number of wate surface bound.|  !
!  | <-- No_sp   | N_SP      | Vertex numbering of the sample      |  !
!  | <-- nbe     | N_CELL    | Type of boundary cell (inside or bc)|  !
!  | <-- xv      | N_VERT    | x-coordinate of the vertex          |  !
!  | <-- yv      | N_VERT    | y-coordinate of the vertex          |  !
!  | <-- zbv     | N_VERT    | The bottom of the river at vertex   |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
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

      integer,dimension(:,:):: No_vp
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: No_wb 
      integer,dimension(:)  :: No_qb
      integer,dimension(:)  :: No_hb 
      integer,dimension(:)  :: No_sp
      integer,dimension(:)  :: nbe        
      real*8 ,dimension(:)  :: xv
      real*8 ,dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer,dimension(:),  allocatable :: TagSubdomain
      integer,dimension(:),  allocatable :: NewIndex,OldIndex
      integer,dimension(:),  allocatable :: nbe_old
      integer,dimension(:,:),allocatable :: No_cp_old,No_vp_old
      integer,dimension(:),  allocatable :: IniSub,FinSub
      integer,dimension(:),  allocatable :: NumSubBdy,IniSubBdy,FinSubBdy
      integer,dimension(:),  allocatable :: NumSubInt,IniSubInt,FinSubInt
      integer,dimension(:),  allocatable :: ElementsBdy
      integer,dimension(:),  allocatable :: ElementsInt
      integer :: NumberSubdomains,tag,i0,s,s0,sm,j0
      integer :: Ielem,Jelem,ic1,ic2,ic3,m,ii
      integer :: TotalElementsBdy,TotalElementsInt
!     ---------------
      integer :: k1,k2,k3,kv1,kv2,kk
      integer,parameter :: Idisplay = 1
      real*8 :: funh,maxi
      character*80 title     

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: input'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                     Read & write:  data.txt                         !
!                                                                     !
!*********************************************************************!

!---------------------------------------------------------------------!                                                            
!     The file "data.txt" has the following format:                   !     
!                                                                     !   
!         --------NUMBERING OF CELL VERTICES----                      !
!              1       2     323                                      !
!              2       3     323                                      !
!              3       4     324                                      !
!              ....                                                   !
!         --------COORDINATES OF CELL VERTICES----                    !
!              0.00000      -20.00000                                 !
!              4.68750      -20.00000                                 !
!              9.37500      -20.00000                                 !
!              ....                                                   !
!         --------NUMBERING OF CELL CENTERS----                       !
!              0       2    6401                                      !
!              0    6402       1                                      !
!              0    6403    6402                                      !
!              ....                                                   !
!         --------BOTTOM LEVELS ---------------                       !
!              0.000   0.000   0.000   0.000   0.000   ...            !
!              0.000   0.080   0.310   0.550   0.780   ...            !
!              2.190   2.420   2.660   2.890   3.130   ...            !
!              ....                                                   !
!         --------TYPE OF BOUNDARY CELLS ---------------              !
!              1      1      1      1      1      1    ...            !
!              1      1      1      1      1      1    ...            !
!              2      0      0      0      0      0    ...            !
!              ....                                                   !
!         --------NUMBER OF WALL BOUNDARY CELLS-----------            !
!             638                                                     !
!              2      3      4      5      6      7    ...            !
!              8      9     10     11     12     13    ...            !
!             22     23     24     25     26     27    ...            !
!              ....                                                   !
!        --------NUMBER OF DISCHARGE BOUNDARY CELLS-------            !
!             21                                                      !
!            322    643    964   1285   1606   1927 ...               !
!              1                                                      !
!        --------NUMBER OF WATER BOUNDARY CELLS-------                !
!             21                                                      !
!            642    963   1284   1605   1926   2247 ...               !
!        --------NUMBER OF SAMPLING CELLS-------                      !
!              6                                                      !
!            826    684   2365    698   3711   2113                   !
!                                                                     ! 
!                                                                     !        
!---------------------------------------------------------------------!

71    format(A62)
72    format(A25)
73    format(A46,I10)
74    format(A46)
76    format(A46,F15.10)
77    format(A46,E15.5)

!      ________________________________________________________
!     |                                                        |
!     |                    Read data file                      |
!     |________________________________________________________|

!     ________________________________________________________
!     Open file

!     ________________________________
#     if defined(KeyEstuaryGironde)
         open(22,file='EstuaryGironde/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStaticCylinder)
         open(22,file='StaticCylinder/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStaticChannel)
         open(22,file='StaticChannel/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStandingWave)
         open(22,file='StandingWave/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyTaylorVortex)
         open(22,file='TaylorVortex/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyTestOnlyPoisson)
         open(22,file='TestOnlyPoisson/data.txt',status='OLD')
!     ________________________________
#     else
         open(22,file='data.txt',status='OLD')
#     endif

!     ________________________________________________________
!     Read: Number of cell-center and vertex points 

      read(22,74) title
      read(22,74) title
      read(22,74) title
!     ________________________________________________________
!     Read: Numbering of the cell vertices for each element 
  
      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,88) ((No_vp(i,j),j=1,3),i=1,N_CELL0)      
88    format(3I8)
!     ________________________________________________________
!     Read: Coordinates of cell vertices                    
      
      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,87) (xv(nv),yv(nv),nv=1,N_VERT)
87    format(2F15.5)

!     ________________________________________________________
!     Read: Numbering of the surrounding three cell centers

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,88) ((No_cp(i,j),j=1,3),i=1,N_CELL0)
!     ________________________________________________________
!     Read: Bottom levels                                  
         
      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,86) (zbv(i),i=1,N_VERT)
86    format(10F8.3)
!     ________________________________________________________
!     Read: Type of cell (inside or boundary)               
     
      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,85) (nbe(i),i=1,N_CELL0) 
85    format(20I7)
!     ________________________________________________________
!     Read: Inputting for wall boundary vertices          

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,*)  N_WB

      if (N_WB.gt.N_WBMAX) then
        print*,'-----------------------------------------------'
        print*,'              WARNING!!!!!!!                   '
        print*,'-----------------------------------------------'
        print*,'   No. of wall boundary vertices is greater than'
        print*,'   planned. The program stops.'
        print*,'   N_WB    = ',N_WB
        print*,'   N_WBMAX = ',N_WBMAX
        print*,'-----------------------------------------------'
	    stop
      endif

      if (N_WB.ne.0) then
	 read(22,85) (No_wb(nv),nv=1,N_WB)
      else
         read(22,74) title
         if (IDISPLAY.eq.1) write(*,*)  title
      endif
!     ________________________________________________________
!     Read: Inputting for discharge normal to boundary vert
					 
      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,*)  N_QB

      if (N_QB.gt.N_QBMAX) then
        print*,'-----------------------------------------------'
        print*,'              WARNING!!!!!!!                   '
        print*,'-----------------------------------------------'
        print*,'   No. of Q boundary vertices is greater than  '
        print*,'   planned. The program stops.'
        print*,'   N_QB    = ',N_QB
        print*,'   N_QBMAX = ',N_QBMAX
        print*,'-----------------------------------------------'
        stop
      endif	

      if (N_QB.ne.0) then
	 read(22,85) (No_qb(nv),nv=1,N_QB) 
      else
         read(22,74) title
         if (IDISPLAY.eq.1) write(*,*)  title
      endif
!     ________________________________________________________
!     Read: Inputting for water level boundary vertices   
					
      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,*)  N_HB

      if (N_HB.gt.N_HBMAX) then
        print*,'-----------------------------------------------'
        print*,'              WARNING!!!!!!!                   '
        print*,'-----------------------------------------------'
        print*,'   No. of H boundary vertices is greater than  '
        print*,'   planned. The program stops.'
        print*,'   N_HB    = ',N_HB
        print*,'   N_HBMAX = ',N_HBMAX
        print*,'-----------------------------------------------'
        stop
      endif	

      if (N_HB.ne.0) then
	 read(22,85) (No_hb(nv),nv=1,N_HB)
      else
         read(22,74) title
         if (IDISPLAY.eq.1) write(*,*)  title
      endif 
!     ________________________________________________________
!     Read: Number of sampling cells                        

      read(22,74) title
      if (IDISPLAY.eq.1) write(*,*)  title
      read(22,*) N_SP

      if (N_SP.gt.N_SPMAX) then
        print*,'-----------------------------------------------'
        print*,'              WARNING!!!!!!!                   '
        print*,'-----------------------------------------------'
	print*,'   No. of sample vertices is greater than      '
        print*,'   planned. The program stops.'
        print*,'   N_SP    = ',N_SP
        print*,'   N_SPMAX = ',N_SPMAX
        print*,'-----------------------------------------------'
        stop
      endif	

      if (N_SP.ge.1) then
	 read(22,85) (No_sp(nv),nv=1,N_SP)
      else
         read(22,74) title
         if (IDISPLAY.eq.1) write(*,*)  title
      endif
!     ________________________________________________________
!     Close                          

      close(22)

!      ________________________________________________________
!     |                                                        |
!     |              Displaying all the results                |
!     |________________________________________________________|

      if (Idisplay.eq.1) then 
      write(*,*)'                                              '
      write(*,*)'  Number of inside cell centers      N_CELL0 = ',N_CELL0
      write(*,*)'  Number of vertices                 N_VERT  = ',N_VERT
      write(*,*)'  Number of vertices with walls bdy. N_WB    = ',N_WB
      write(*,*)'  Number of vertices with disch bdy. N_QB    = ',N_QB
      write(*,*)'  Number of vertices with water bdy. N_HB    = ',N_HB
      write(*,*)'  Number of sampling vertices        N_SP    = ',N_SP
      write(*,*)'                                             '
      endif

!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!

      if (Idisplay.ne.1) then
         write(*,'(t12,20a)') ' --- data.txt has been successfully read ---'
      endif
	
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: input'
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

