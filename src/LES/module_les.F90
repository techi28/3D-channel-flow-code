!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               DEFINITION OF THE LES VARIABLES                       !
!                            June 2016                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      MODULE les

#     include "cppdefs.h"
      implicit none
!     ====================================
!     ==========  SEQUENTIAL =============
!#     ifndef KeyParallel
!#        include "common.mpf"
!#     endif
!     =============== END ================
!     ====================================
!      ________________________________________________________
!     |                                                        |
!     |                  Global variables                      |
!     |________________________________________________________|
      real*8, dimension(:,:), allocatable :: walld
      real*8, dimension(:,:), allocatable :: eddy
      integer,dimension(:,:),allocatable :: utau_pair,utau_j
!     ________________________________________________________
!     LES Parameters
      integer :: sg_init = 1  !Initialaztion Flag
      integer,parameter :: sg_model = 1 !Subgrid model identifier
      real*8,parameter :: sg_smag_const = 0.1d0 ! Smagorinsky constant
      real*8,parameter :: sg_smag_aplus = 26.0d0 !VanDriest Damping length
      real*8,parameter :: sg_smag_mason = 2.0d0 !Mason wall matching power
      real*8,parameter :: sg_smag_so = 0.0d0 !Thomas shear constant
      real*8,parameter :: von_Karman = 0.415 !von Karman constant

      CONTAINS

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            SUBROUTINE: DEALLOCATE STATS VARIABLES                   !
!                          JUNE 2016                                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

        SUBROUTINE alloc_les

#	include "common.mpf"

!       --------------------------------------
!       allocate wall normal distance array
        allocate(walld(N_CELL,NZ),eddy(N_CELL,NZ))
        allocate(utau_pair(N_CELL,NZ),utau_j(N_CELL,NZ))
!      ---------------------------------------
!       Set Wall Distance to a very large value
          do k=1,NZ
           do i=1,N_CELL
             walld(i,k) = 1.0E5
             eddy(i,k) = 0.0d0
           enddo
          enddo

	   END SUBROUTINE alloc_les
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            SUBROUTINE: DEALLOCATE STATS VARIABLES                   !
!                          JUNE 2016                                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

        SUBROUTINE dealloc_les

!      ________________________________________________________
!     |                                                        |
!     |                  FLuid variables                       |
!     |________________________________________________________|

	    deallocate(walld,eddy)
        deallocate(utau_pair,utau_j)

	   END SUBROUTINE dealloc_les

       END MODULE les

