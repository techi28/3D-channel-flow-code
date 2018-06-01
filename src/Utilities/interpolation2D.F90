!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   INTERPOLATION OF MAIN VARIABLES                   !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE interpolation2D(funv,xv,yv,No_vp,nbev,&
                                 fun,xc,yc,No_cp,nbe)       
!---------------------------------------------------------------------!
!                                                                     !
!    This program interpolate the values phi from the center of the   !
!    cell to the vertex of the triangles. It is a 2D interpolation    !  
!    in the horizontal direction. It is important to remark that we   !
!    need the exact values in the boundaries of our domain.           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !
!  | <-- funv    |(N_VERT)   | Function phi at the vertex          |  !  
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> fun     |(N_CELL)   | Function phi at the cell center     |  !
!  | --> No_vp   |(N_CELL0,3)| Numbering of cell vertices          |  !
!  | --> nbev    |(N_VERT)   | Type of tag about the kind of vertex|  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Common parameters used:                                          !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  ! 
!  |--- N_CELL0   | Number of the cell centers inside the domain   |  !
!  |--- N_CELL    | Total number of cell centers                   |  !
!  |--- N_VERT    | Number of the computing vertices               |  !
!  |--- NZ        | Number of points in the vertical direction     |  !  
!  |______________|________________________________________________|  !  
!  |    dlCV      |(CELL,3)  | Distance from the center to vertex  |  !
!  |    dlVsum    | N_VERT   | Total sum of dlCV distances at vert.|  !
!  |    areaCELL  |(N_CELL,3)| Area of each cell                   |  !
!  |    areaVsum  | N_VERT   | Total area related to each vertex.  |  !
!  |______________|__________|_____________________________________|  !
!                                                                     !
!    Common variables modified:                                       !
!   _______________________________________________________________   !
!  |   Name       |                 Description                    |  !  
!  |______________|________________________________________________|  !   
!  | * nv,i,k     | Loop counters: vertices,cells, other           |  !
!  |______________|________________________________________________|  !
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

!      ________________________________________________________
!     |                                                        |
!     |   Keys and common parameters                           |
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
!     |    Declaration of variables                            |
!     |________________________________________________________|

      real*8 ,dimension(:)  :: funv(N_VERT)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  

      real*8 ,dimension(:)  :: fun(N_CELL)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!      ________________________________________________________
!     |                                                        |
!     |    Declaration of local variables                      |
!     |________________________________________________________|

      integer :: nv1,nv2,nv3,nv4
      integer :: nc1,nc2,nc3,nc4,nc5,nc6,nc7,nc8
      real*8  :: dlextra
!     -------------------------
      real*8 ,dimension(:), allocatable :: dfundx,dfundy
      integer,dimension(:), allocatable :: NumFun
      real*8 ,dimension(:), allocatable :: SumFun
      real*8 :: sumfx,sumfy,deter
      real*8 :: funpGF
      integer:: jc,nx,ny
!     -------------------------
      real*8 :: dxCV,dyCV,dlCVout
      real*8 :: f0,f1,f2,f3
      integer:: jj,jv1,jv2,jc1,jc2,jc3
!     -------------------------
      integer,dimension(:) :: dlVsumout(N_VERT) 
      integer :: ChooseLSMDist
      integer :: elem,jvert,kvert

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: interpolation2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                          Initialization                             !
!                                                                     !
!*********************************************************************!

      do nv=1,N_VERT
         funv(nv) = 0.0d0
      enddo

!      ________________________________________________________
!     |                                                        |
!     |  USING NEW TECHNIQUE OF DISTANCE WEIGHTING COLLECTING  |
!     |________________________________________________________|
      
#     ifdef KeyInterpoNew
         do nv=1,N_VERT
            do j=1,Dimsurrounding(nv)
               nc = surrounding(nv,j)
               funv(nv)= funv(nv) + weight(nv,j)*fun(nc) 
            enddo
            funv(nv) = funv(nv)/dlVsum(nv)
         enddo
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                USING DISTANCE WEIGHTING                |
!     |________________________________________________________|
      
#     ifdef KeyInterpoDist
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
!           __________________
!           Inside cells 
            do i=1,N_CELL0
               do j=1,3
                  nv = No_vp(i,j)
                  funv(nv)= funv(nv) + dlCV(i,j)*fun(i) 
               enddo
            enddo
!           __________________
!           Ghost cells
#           ifdef KeyUseInterGhost
            do i=1,N_CELL0
               if (nbe(i).ne.0) then
                  do j=1,3
                     nc = No_cp(i,j)
                     if (nc.lt.1.OR.nc.gt.N_CELL0) then
                        jj=j+1
                        if (jj.gt.3) jj=jj-3
                        jv1 = No_vp(i,j)
                        jv2 = No_vp(i,jj)
                        funv(jv1) = funv(jv1) + dlCV(i,j)*fun(nc)
                        funv(jv2) = funv(jv2) + dlCV(i,jj)*fun(nc)
                     endif 
                  enddo
               endif
            enddo
#           endif
!           __________________
!           Average  
            do nv=1,N_VERT
               funv(nv) = funv(nv)/dlVsum(nv)
            enddo
#        endif
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
!           __________________
!           Inside cells 
            do i=1,N_CELL0
	       do j=1,3
	          nv = No_vp(i,j)
                  funv(nv)= funv(nv) + dlCV(i,j)*fun(i) 
               enddo
            enddo
!           __________________
!           Ghost cells
#           ifdef KeyUseInterGhost
            do i=1,N_CELL0
               if (nbe(i).ne.0) then	
   	          do j=1,3
	             nc = No_cp(i,j)
                     elem = index_global(nc)
	             if (elem.lt.1.OR.elem.gt.N_CELL0global) then
		        jj=j+1
		        if (jj.gt.3) jj=jj-3
		        jv1 = No_vp(i,j)
		        jv2 = No_vp(i,jj)
                        funv(jv1) = funv(jv1) + dlCV(i,j)*fun(nc)
                        funv(jv2) = funv(jv2) + dlCV(i,jj)*fun(nc)
                     endif 
	          enddo
               endif
            enddo
#           endif
!           __________________
!           Overlaping cells 
            do i = N_CELL0+N_CELLghost+1,N_CELL
               elem = index_global(i)
               do j=1,3
                  jvert = No_vp_global(elem,j)
                  do nv=1,N_VERT
                     kvert = index_globalv(nv)  
                     if (jvert.eq.kvert) then
                         funv(nv)= funv(nv) + dlCV(i,j)*fun(i)
                     endif
                  enddo
	       enddo
            enddo     
!           __________________
!           Average  
            do nv=1,N_VERT
               funv(nv) = funv(nv)/dlVsum(nv)
            enddo
#        endif	
!        =============== END ================    
!        ====================================
#     endif
!      ________________________________________________________
!     |                                                        |
!     |         USING AREA WEIGHTING TO INTERPOLATE            |
!     |________________________________________________________|
    
#     ifdef KeyInterpoArea
          do i=1,N_CELL0	
	     do j=1,3
	        nv = No_vp(i,j)
                funv(nv)= funv(nv) + areaCell(i)*fun(i) 
             enddo
          enddo

          do nv=1,N_VERT
             funv(nv) = funv(nv)/areaVsum(nv)
          enddo
#     endif

!      ________________________________________________________
!     |                                                        |
!     |      Using more points in the cell-vertex distance     |
!     |________________________________________________________|
      
#     ifdef KeyInterpoDist2
          do i=1,N_CELL0	
             nv1 = No_vp(i,1)
             nv2 = No_vp(i,2)
             nv3 = No_vp(i,3)
             nc1 = No_cp(i,1)
             nc2 = No_cp(i,2)
             nc3 = No_cp(i,3)
!            ---------------------------------------
!            Vertex 1
             if (nbev(nv1).eq.0) then
                funv(nv1)= funv(nv1) + dlCV(i,1)*fun(i)
                if ((nc2.ge.1).and.(nc2.le.N_CELL0)) then 
                   dlextra    = sqrt((xc(nc2)-xv(nv1))**2 &
                                    +(yc(nc2)-yv(nv1))**2)
                   funv(nv1)  = funv(nv1) + (1.0d0/dlextra)*fun(nc2)
                   dlVsum(nv1)= dlVsum(nv1)+(1.0d0/dlextra)
                endif
             endif
!            ---------------------------------------
!            Vertex 2
             if (nbev(nv2).eq.0) then
                funv(nv2)= funv(nv2) + dlCV(i,2)*fun(i)
                if ((nc3.ge.1).and.(nc3.le.N_CELL0)) then 
                   dlextra    = sqrt((xc(nc3)-xv(nv2))**2 &
                                    +(yc(nc3)-yv(nv2))**2)   
                   funv(nv2)  = funv(nv2) + (1.0d0/dlextra)*fun(nc3)
                   dlVsum(nv2)= dlVsum(nv2)+(1.0d0/dlextra)
                endif
             endif
!            ---------------------------------------
!            Vertex 3
             if (nbev(nv3).eq.0) then
                funv(nv3)= funv(nv3) + dlCV(i,3)*fun(i)
                if ((nc1.ge.1).and.(nc1.le.N_CELL0)) then 
                   dlextra     = sqrt((xc(nc1)-xv(nv3))**2 &
                                     +(yc(nc1)-yv(nv3))**2)
                   funv(nv3)   = funv(nv3)+(1.0d0/dlextra)*fun(nc1)
                   dlVsum(nv3) = dlVsum(nv3)+(1.0d0/dlextra)
                endif
             endif
          enddo

          do nv=1,N_VERT
             if (nbev(nv).eq.0) then
                funv(nv) = funv(nv)/dlVsum(nv)
             endif
          enddo
#     endif

!      ________________________________________________________
!     |                                                        |
!     |            Using the gradient of the LSM               |
!     |________________________________________________________|
      
#     ifdef KeyInterpoLSM

!         --------------------------------------------
!         Allocate
          allocate(dfundx(N_CELL),dfundy(N_CELL),Sumfun(N_VERT))

!         --------------------------------------------
!         Gradient of fun

          !do i=1,N_CELL0
             !sumfx = 0.0d0
             !sumfy = 0.0d0
             !do j=1,3
             !   jc = No_cp(i,j)
             !   sumfx = sumfx + dxCC(i,j)*(fun(jc)-fun(i))
             !   sumfy = sumfy + dyCC(i,j)*(fun(jc)-fun(i))
             !enddo
             !deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
             !dfundx(i)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
             !dfundy(i)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter          
          !enddo
!         --------------------------------------------
!         Adding distance weights 
          ChooseLSMDist = 1
!         --------------------------------------------
!         Contribution to each vertex
          SumFun = 0
          do i=1,N_CELL0
             f0  = fun(i)
             f1  = fun(No_cp(i,1))
             f2  = fun(No_cp(i,2))
             f3  = fun(No_cp(i,3))
             dfundx(i)= aGx(i,0)*f0+aGx(i,1)*f1+aGx(i,2)*f2+aGx(i,3)*f3
             dfundy(i)= aGy(i,0)*f0+aGy(i,1)*f1+aGy(i,2)*f2+aGy(i,3)*f3
             !------------  
             do j=1,3
                nv = No_vp(i,j)
                funpGF =  fun(i) + dfundx(i)*(xv(nv)-xc(i)) &
                                 + dfundy(i)*(yv(nv)-yc(i))                  
                if (ChooseLSMDist.eq.1) then
                   funv(nv)   = funv(nv)   + funpGF*dlCV(i,j)
                   SumFun(nv) = SumFun(nv) + dlCV(i,j)
                else
                   funv(nv)   = funv(nv)   + funpGF
                   SumFun(nv) = SumFun(nv) + 1.0d0
                endif
             enddo
          enddo
!         --------------------------------------------
!         Final average value
          do nv=1,N_VERT
             funv(nv) = funv(nv)/SumFun(nv)
          enddo

!        --------------------------------------------
!        Deallocate
         deallocate(dfundx,dfundy,SumFun)

#     endif

!      ________________________________________________________
!     |                                                        |
!     | Using the gradient of the LSM + Distance + Conditional |
!     |________________________________________________________|
      
#     ifdef KeyInterpoMix

!         --------------------------------------------
!         Allocate
          allocate(dfundx(N_CELL),dfundy(N_CELL),Numfun(N_VERT),Sumfun(N_VERT))

!         --------------------------------------------
!         Gradient of fun

          do i=1,N_CELL0	
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
          enddo

!         --------------------------------------------
!         Number of elements that share each vertex
          do nv=1,N_VERT
             NumFun(nv) = 0
             SumFun(nv) = 0
          enddo
          do i=1,N_CELL0
             do j=1,3
                nv = No_vp(i,j)
                NumFun(nv) =  NumFun(nv) + 1
             enddo
          enddo
!         --------------------------------------------
!         Contribution to each vertex
          do i=1,N_CELL0	
	     do j=1,3
	        nv = No_vp(i,j)
                if (NumFun(nv).le.2) then
                   funpGF =  fun(i) 
                else
                   funpGF =  fun(i) + dfundx(i)*(xv(nv)-xc(i)) &
                                    + dfundy(i)*(yv(nv)-yc(i))
                endif
                funv(nv)   = funv(nv)   + funpGF*dlCV(i,j)
                SumFun(nv) = SumFun(nv) + dlCV(i,j)
             enddo
          enddo

          do nv=1,N_VERT
             funv(nv) = funv(nv)/SumFun(nv)
          enddo

!        --------------------------------------------
!        Deallocate
         deallocate(dfundx,dfundy,NumFun,SumFun)

#     endif


!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: interpolation2D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                                                                     !
!                      	   END OF INTERPOLATION                       !
!                                                                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
