!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         Construct PPE Matrix 3D                     !
!                               Oct 2016                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE    PDE_center(xc,yc,sig,dsig,No_cp,nbe,            &
                               xv,yv,sigv,dsigv,No_vp,nbev)
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
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      real*8, dimension(:)   :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0)
!     --------------------------------------
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
      real*8, dimension(:) :: AA(1:3)
      real*8 :: AA0,AAB,AAT
      real*8 :: co1
      integer:: jc

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: PDE_center'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     do k=2,NZ-1
        do i=1,N_CELL0
!         __________________________________________
!         For XY Plane
           do j =1, 3
               jc = No_cp(i,j)
#              ifdef KeyKimChoi
               co1 = dle1(i,j)/dlCC(i,j)
#              else
               co1 = dsqrt(xe1(i,j)*xe1(i,j)+ye1(i,j)*ye1(i,j))/dlCC(i,j)
#              endif
!              -------------------------
               if(bcpc(jc) .eq. 1) then ! Newmann bc
                   AA(j) = 0.0d0
                   goto 100
               endif
!             --------------
!              principle diffusion
                 AA(j) =  co1*dlVV(i,j)*dsigv(k-1)
100            continue
           enddo
!         __________________________________________
!         Bottom
                 AAB = areaCell(i)/dsig(k-1)
!         __________________________________________
!         Top
                 AAT = areaCell(i)/dsig(k)
!         __________________________________________
!         For no periodic bc
          if(ZPB .eq. 0) then
                 if(k .eq. 2)    then
                    if(bcbot .eq. 1) AAB = 0.0d0
                 endif

                 if(k .eq. NZ-1)  then
                    if(bctop .eq. 1) AAT = 0.0d0
                 endif
          endif


          AA0 = -1.0d0*(AAT+AAB+AA(1)+AA(2)+AA(3))

!          |----------------------------------------------------------|
!          |              Final coefficients contributions            |
!          |----------------------------------------------------------|
	        AmC0(i,k) = AA0
            AmC1(i,k) = AA(1)
            AmC2(i,k) = AA(2)
            AmC3(i,k) = AA(3)
            AmCT(i,k) = AAT
            AmCB(i,k) = AAB
         enddo
      ENDDO

!  ---------------------------------------------------------------
!        Formation of Matrix coefficent for Ax=b
        call get_matrix_coef(AmC0,AmC1,AmC2,AmC3,AmCT,AmCB,   &
                             xc,yc,sig,dsig,No_cp,nbe)
!     ------------------------------------------------------------
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: PDE_center'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF PDE                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

