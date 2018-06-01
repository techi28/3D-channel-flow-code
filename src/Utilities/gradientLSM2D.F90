!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                       GRADIENT USING LSM 2D                         !
!                             Jun 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE grandientLSM2D(dfundx,dfundy, &
                                fun,No_cp,nbe)     
 
!---------------------------------------------------------------------!
!                                                                     !
!    This program approximates the gradient of a function fun using   !
!    the least square technique in 2D (Appendix A).                   !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name      |    Size   | Description                       |  !  
!  |_______________|___________|___________________________________|  !
!  | <-- dfundx    |(N_CELL)   | d(fun)/dx = gradient component x  |  !
!  | <-- dfundy    |(N_CELL)   | d(fun)/dy = gradient component y  |  !  
!  |_______________|___________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> fun     |(N_CELL)   | Function "fun" at the center element|  !
!  | --> nbe     | N_CELL    | Type of boundary cell (inside or bc)|  !
!  | --> No_cp   |(N_CELL,3) | Node No. of surrounding three cell  |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name     |                   Description                    |  !  
!  |____________|__________________________________________________|  ! 
!  |--- N_CELL  |  Total number of the cells                       |  !
!  |--- N_CELL0 |  Number of the cell centers inside the domain    |  !
!  |____________|__________________________________________________|  ! 
!  |  dxCC      |(N_CELL0,3)| = dx from cell center to cell neigb. |  !
!  |  dyCC      |(N_CELL0,3)| = dy from cell center to cell neigb. |  !
!  |  sum_xc2   | N_CELL0   | Sum(xc(i)-xc(neighborn))^2           |  !
!  |  sum_yc2   | N_CELL0   | Sum(yc(i)-yc(neighborn))^2           |  !
!  |  sum_xcyx  | N_CELL0   | Sum(xc(i)-xc(neig))*(yc(i)-yc(neig)) |  !
!  |____________|___________|______________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | deter       |  determinant of the horizontal sums             |  !
!  | sumfx,sumfy |  real variables to add gradient components      |  !
!  |_____________|_________________________________________________|  !   
!  | dfBdn       |  Value of the derivative at the boundary        |  !
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

      real*8 ,dimension(:,:) :: dfundx(N_CELL)
      real*8 ,dimension(:,:) :: dfundy(N_CELL)
      real*8 ,dimension(:,:) :: fun(N_CELL)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0) 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 :: sumfx,sumfy,deter
      real*8 :: dfBdn
      integer:: jc,jc2,elem

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: GradientLSM 2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                         Horizontal gradient                         !
!                                                                     !
!*********************************************************************!

      do i=1,N_CELL0	
!         __________________________________________________
!        |                                                  |
!        |          Inside & boundary cell-centers          |
!        |__________________________________________________|

	  sumfx = 0.0d0
	  sumfy = 0.0d0
	  do j=1,3
	     jc = No_cp(i,j)                
	     sumfx = sumfx + dxCC(i,j)*(fun(jc)-fun(i))
	     sumfy = sumfy + dyCC(i,j)*(fun(jc)-fun(i))     
	  enddo
	  deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
          dfundx(i)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
	  dfundy(i)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
!         __________________________________________________
!        |                                                  |
!        |         Outside cell-centers  (ficticious)       |
!        |__________________________________________________|

!        ----------------------------------
!        Wall 
	 if (nbe(i).eq.1) then  
	    do j=1,3
	       nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                  dfBdn = 0.0d0
	          dfundx(nc) = 2.0d0*dfBdn+dfundx(i) 
		  dfundy(nc) = 2.0d0*dfBdn+dfundy(i)	
	       endif
	    enddo 
	 endif  
!        ----------------------------------
!        Discharge normal 
	 if (nbe(i).eq.2) then  
	    do j=1,3
	       nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                  dfBdn = 0.0d0
	          dfundx(nc) = 2.0d0*dfBdn+dfundx(i)
	          dfundy(nc) = 2.0d0*dfBdn+dfundy(i)	
	      endif
	    enddo 
	 endif  
!        ----------------------------------
!        Water level 
	 if (nbe(i).eq.3) then  
	    do j=1,3
	       nc=No_cp(i,j)
!                  ====================================
!                  ==========  SEQUENTIAL =============
#                  ifndef KeyParallel
	           if (nc.lt.1.OR.nc.gt.N_CELL0) then
!                  ====================================
!                  =====  START PARALLEL OPTION =======
#                  else
                   elem = index_global(nc)
	           if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                  endif	
!                  =============== END ================    
!                  ====================================
                  dfBdn = 0.0d0
	          dfundx(nc) = 2.0d0*dfBdn+dfundx(i)
	          dfundy(nc) = 2.0d0*dfBdn+dfundy(i)	
	       endif
	    enddo 
	 endif  
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication2D(dfundx)
         call communication2D(dfundy)
#     endif	
!     =============== END ================    
!     ====================================

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: GradientLSM 2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	    END OF GRADIENT 2D                        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
