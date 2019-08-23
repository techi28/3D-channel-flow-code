!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!         Calculate the pressure coefficient around cylinder          !
!                             March 2016                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      SUBROUTINE update_stats(suf,svf,swf,spf,   &
                              ufnp,vfnp,wfnp,pfnp)

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
!      ________________________________________________________
!     |                                                        |
!     |           Definition of variables                      |
!     |________________________________________________________|
      real*8,dimension(:,:):: suf(N_CELL,NZ)
      real*8,dimension(:,:):: svf(N_CELL,NZ)
      real*8,dimension(:,:):: swf(N_CELL,NZ)
      real*8,dimension(:,:):: spf(N_CELL,NZ)
      real*8,dimension(:,:):: ufnp(N_CELL,NZ)
      real*8,dimension(:,:):: vfnp(N_CELL,NZ)
      real*8,dimension(:,:):: wfnp(N_CELL,NZ)
      real*8,dimension(:,:):: pfnp(N_CELL,NZ)
!     ---------------------------------
!     Local variable
      integer :: Idisplay
!     ----------------------------------
#     ifndef KeyParallel
         Idisplay = 1
#     else
          if(rang_topo .eq. 0) then
            Idisplay = 1
          else
            Idisplay = 0
          endif
#     endif
!    ------------------------------------
       scount = scount + 1

      if(Idisplay .eq. 1) then
         print*,'   =================================== '
         print*,'          Field Info Sample Count :',scount
         print*,'   =================================== '
      endif

       do k=1,NZ
          do i=1,N_CELL
            suf(i,k) = suf(i,k) + ufnp(i,k)
            svf(i,k) = svf(i,k) + vfnp(i,k)
            swf(i,k) = swf(i,k) + wfnp(i,k)
            spf(i,k) = spf(i,k) + pfnp(i,k)
          enddo
       enddo

      RETURN
      END
