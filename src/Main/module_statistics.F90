!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               DEFINITION OF THE STATS VARIABLES                     !
!                            June 2015                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

	MODULE statistics

!---------------------------------------------------------------------!
!                                                                     !
!      SUBROUTINES:     -  alloc_stats_variables                      !
!                       -  dealloc_stats_varaibles                    !
!                                                                     !
!---------------------------------------------------------------------!

	implicit none

!      ________________________________________________________
!     |                                                        |
!     |                  Stats variables                       |
!     |________________________________________________________|

        real*8,dimension(:,:),allocatable :: &
             sm_uf, &      ! Velocity U
             sm_vf, &      ! Velocity V
             sm_wf, &      ! Velocity W
             sm_u2f, &     ! Velocity Products UU
             sm_v2f, &     ! Velocity Products VV
             sm_w2f, &     ! Velocity Products WW
             sm_uvf, &     ! Velocity Products UV
             sm_vwf, &     ! Velocity Products VW
             sm_wuf, &     ! Velocity Products WU
             sm_pf,  &     ! Pressure Term P
             sm_p2f, &     !
             sm_ef         ! Turbulent Viscosity

!      ________________________________________________________
!     |                                                        |
!     |                       Others                           |
!     |________________________________________________________|

!     ----------------------------------------
!     Count of the sample number
      integer :: sm_n

      CONTAINS


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            SUBROUTINE: ALLOCATE STATS VARIABLES                     !
!                            June 2015                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

	SUBROUTINE alloc_stats_variables

#	include "common.mpf"

!      ________________________________________________________
!     |                                                        |
!     |            FLuid Stats variables                       |
!     |________________________________________________________|

	     allocate(sm_uf(N_CELL,NZ),   &
                  sm_vf(N_CELL,NZ),   &
                  sm_wf(N_CELL,NZ),   &
                  sm_u2f(N_CELL,NZ),  &
                  sm_v2f(N_CELL,NZ),  &
                  sm_w2f(N_CELL,NZ),  &
                  sm_uvf(N_CELL,NZ),  &
                  sm_vwf(N_CELL,NZ),  &
                  sm_wuf(N_CELL,NZ),  &
                  sm_pf(N_CELL,NZ),   &
                  sm_p2f(N_CELL,NZ),  &
                  sm_ef(N_CELL,NZ))

!      ---------------------------------------------------------

	END SUBROUTINE alloc_stats_variables


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            SUBROUTINE: DEALLOCATE STATS VARIABLES                   !
!                          JUNE 2015                                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


    SUBROUTINE dealloc_stats_variables

!      ________________________________________________________
!     |                                                        |
!     |                  FLuid variables                       |
!     |________________________________________________________|

	deallocate(sm_uf,sm_vf,sm_wf,     &
               sm_u2f,sm_v2f,sm_w2f,  &
               sm_uvf,sm_vwf,sm_wuf,  &
               sm_pf,sm_p2f,sm_ef)

	END SUBROUTINE dealloc_stats_variables

    END MODULE statistics