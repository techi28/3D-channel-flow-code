
!=====================================================================!
!---------------------------------------------------------------------!
!             Turbulence Stats (UNSTRUCTURED MESH)                    !
!                              Xin Bai                                !
!                             June 2015                               !
!---------------------------------------------------------------------!
!=====================================================================!

!—————————————————————————————————————————————————————————————————————---!
! This subroutine calculates the turbulent statistical terms, including  !
! time averaged velocity, velocity fluctuation, Reynolds stress and etc. !
! Currently only valid for rectangle structured grid.                    !
!========================================================================!
!                        SUBROUNTINE:                                    !
!------------------------------------------------------------------------!
! This subroutine aims to initialize the variables used for LES by using !
! the existing velocity and grid information for original source code.   !
!------------------------------------------------------------------------!
!          NOTE: THE PROGRAMME HAVEN'T BEEN PARRALLELIZED YET.           !
!========================================================================!
!                         1.1 sm_init                                    !
!                         1.2 sm_sample                                  !
!                         1.3 sm_put_mm                                  !
!                         1.4 sm_get_mm                                  !
!                         1.5 sm_put_plane_mm                            !
!------------------------------------------------------------------------!
!      _______________________________________________________________
!     |      |                                                        |
!     | 1.1  |                 Initialization Stats                   |
!     |______|________________________________________________________|


     subroutine sm_init()


#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE stats
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE stats
            implicit none
#     endif
!     =============== END ================
!     ====================================
#     ifdef KeyDbg
        print*, '      ----> Begin subroutine: Initialize_stats '
#     endif
!     -------------------------------------------------
!     -----    Set the initial values to zero      ----
!     -------------------------------------------------
     do k=1,NZ
        do i=1, N_CELL0
            sm_uf(i,k)  = 0.0d0
            sm_vf(i,k)  = 0.0d0
            sm_wf(i,k)  = 0.0d0
            sm_u2f(i,k) = 0.0d0
            sm_v2f(i,k) = 0.0d0
            sm_w2f(i,k) = 0.0d0
            sm_uvf(i,k) = 0.0d0
            sm_vwf(i,k) = 0.0d0
            sm_wuf(i,k) = 0.0d0
            sm_pf(i,k)  = 0.0d0
            sm_p2f(i,k) = 0.0d0
            sm_ef(i,k)  = 0.0d0
        enddo
     enddo
            sm_n = 0
!    ---------------------------------------------------
#     ifdef KeyDbg
          print*, '      <---- End subroutine: Initialize_stats '
#     endif

      return
      end subroutine sm_init


!      _______________________________________________________________
!     |      |                                                        |
!     | 1.2  |                 Taking Stats Samples                   |
!     |______|________________________________________________________|
     subroutine glob_sample (phiu,phiv,phiw,phip)


#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE stats
#     ifdef KeyLES
            USE les
#     endif
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE stats
#     ifdef KeyLES
            USE les
#     endif
            implicit none
#     endif
!     =============== END ================
!     ====================================

!===============================================
!  Latest field variables at Fluids Cell Center
      real*8,dimension(:,:) :: phiu(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_CELL,NZ)
      real*8,dimension(:,:) :: phiw(N_CELL,NZ)
      real*8,dimension(:,:) :: phip(N_CELL,NZ)

!================================================
!==============Local Variables===================
      real*8 :: tu,tv,tw,tp
      real*8 :: tu2,tv2,tw2,tp2
      real*8 :: tuv,tvw,twu
!========================================================
#  ifdef KeyDbg
     print*, '      ----> Begin subroutine: sample_stats '
#  endif
!    -----------------------------------------------
#  ifdef KeyDisplay
     write(*,*) 'Sample Index = ', sm_n+1
#  endif
!     -------------------------------------------------
!      Set the initial values to zero
!     -------------------------------------------------
     do k=1,NZ
        do i=1, N_CELL0
        tu = phiu(i,k)
        tv = phiv(i,k)
        tw = phiw(i,k)
        tp = phip(i,k)
        tp2 = tp * tp
        tu2 = tu * tu
        tv2 = tv * tv
        tw2 = tw * tw
        tuv = tu * tv
        tvw = tv * tw
        twu = tw * tu
        sm_uf(i,k) = sm_uf(i,k)+ tu
        sm_vf(i,k) = sm_vf(i,k)+ tv
        sm_wf(i,k) = sm_wf(i,k)+ tw
        sm_pf(i,k) = sm_pf(i,k)+ tp
        sm_u2f(i,k) = sm_u2f(i,k)+ tu2
        sm_v2f(i,k) = sm_v2f(i,k)+ tv2
        sm_w2f(i,k) = sm_w2f(i,k)+ tw2
        sm_uvf(i,k) = sm_uvf(i,k)+ tuv
        sm_vwf(i,k) = sm_vwf(i,k)+ tvw
        sm_wuf(i,k) = sm_wuf(i,k)+ twu
        sm_p2f(i,k) = sm_p2f(i,k)+ tp2
!       -------------------------------
#       ifdef KeyLES
        sm_ef(i,k) = sm_ef(i,k)+ eddy(i,k)
#       endif
        enddo
     enddo
!    -----------------------------------
!      Update the sampling counter
            sm_n = sm_n+1
!    -----------------------------------
!    output sample times
     if(mod(sm_n,10) .eq. 0) then
#      ifdef KeyParallel
          if(rang_topo .eq. 0) then
            print*, '       =========================================='
            print*, '               STATS SAMPLED',sm_n,'TIMES'
            print*, '       =========================================='
          endif
#      else
           print*, '       =========================================='
           print*, '               STATS SAMPLED',sm_n,'TIMES'
           print*, '       =========================================='
#      endif
     endif

! ------------------------
#    ifdef KeyDbg
         print*, '      <---- End subroutine: sample_stats '
#    endif

      return
      end subroutine glob_sample
!      _______________________________________________________________
!     |      |                                                        |
!     | 1.5  |            Output Plane Averaged Stats                 |
!     |______|________________________________________________________|


     SUBROUTINE glob_put_plane_mm()


#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE stats
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE stats
            implicit none
#     endif
!     =============== END ================
!     ====================================

!===============================================
!--------------LOCAL Variables ------------------
      real*8,dimension(:),allocatable :: Su1
      real*8,dimension(:),allocatable :: Sv1
      real*8,dimension(:),allocatable :: Sw1
      real*8,dimension(:),allocatable :: Su2
      real*8,dimension(:),allocatable :: Sv2
      real*8,dimension(:),allocatable :: Sw2
      real*8,dimension(:),allocatable :: Suv
      real*8,dimension(:),allocatable :: Svw
      real*8,dimension(:),allocatable :: Swu
      real*8 :: u1m,u2m,uvm,v1m,v2m,vwm,w1m,w2m,wum
      real*8 :: urms, vrms, wrms, tkesum
      real*8 :: xrms, yrms, zrms
      real*8 :: Tuvm, Tvwm, Twum, npjx,yplus
      integer :: ns,irec,IDWRITE,kNZ
      character*30 filedat
      character*30 filedat1
      character*30 filedat2
!     ----------------------
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '=======> Begin subroutine: output_plane_sm'
         write(*,*) ''
#     endif

!======================================================================
!      Output data format
6     format(I3,4(2x,e12.5))
!     -------------------------------------------------
!     ALLOCATE THE VARIABLES
      allocate (Su1(NZglobal),Su2(NZglobal),Suv(NZglobal), &
                Sv1(NZglobal),Sv2(NZglobal),Svw(NZglobal), &
                Sw1(NZglobal),Sw2(NZglobal),Swu(NZglobal))
       ns = sm_n
!     -------------------------------------------------
!      Set the initial values to zero
!     -------------------------------------------------
     do k=1,NZglobal
        Su1(k) = 0.0d0
        Sv1(k) = 0.0d0
        Sw1(k) = 0.0d0
        Su2(k) = 0.0d0
        Sv2(k) = 0.0d0
        Sw2(k) = 0.0d0
        Suv(k) = 0.0d0
        Svw(k) = 0.0d0
        Swu(k) = 0.0d0
     enddo
!     -------------------------------------------------
!      Update the averaged plane values
     do k=2,NZ-1
        kNZ = k
       do i=1, N_CELL0
        Su1(kNZ) = Su1(kNZ)+sm_uf(i,k)/ns
        Sv1(kNZ) = Sv1(kNZ)+sm_vf(i,k)/ns
        Sw1(kNZ) = Sw1(kNZ)+sm_wf(i,k)/ns
        Su2(kNZ) = Su2(kNZ)+sm_u2f(i,k)/ns
        Sv2(kNZ) = Sv2(kNZ)+sm_v2f(i,k)/ns
        Sw2(kNZ) = Sw2(kNZ)+sm_w2f(i,k)/ns
        Suv(kNZ) = Suv(kNZ)+sm_uvf(i,k)/ns
        Svw(kNZ) = Svw(kNZ)+sm_vwf(i,k)/ns
        Swu(kNZ) = Swu(kNZ)+sm_wuf(i,k)/ns
       enddo
     enddo
!    -----------------------------------------------------
!    Sum the stats from all processors
#    ifdef KeyParallel
        call MPI_Barrier(comm3D,code)
        call SUM_parallelr(Su1,NZglobal)
        call SUM_parallelr(Sv1,NZglobal)
        call SUM_parallelr(Sw1,NZglobal)
        call SUM_parallelr(Su2,NZglobal)
        call SUM_parallelr(Sv2,NZglobal)
        call SUM_parallelr(Sw2,NZglobal)
        call SUM_parallelr(Suv,NZglobal)
        call SUM_parallelr(Svw,NZglobal)
        call SUM_parallelr(Swu,NZglobal)
#    endif
!    -----------------------------------------------------
!     Calculate and output the Stats
!    -----------------------------------------------------
#    ifdef KeyParallel
       npjx = 1.0d0/N_CELL0global
       if(rang_topo .eq. 0) then
        IDWRITE = 1
        else
        IDWRITE = 0
        endif
#    else
       npjx = 1.0d0/N_CELL0
       IDWRITE = 1
#    endif
!    ------------------------------------------------------
    IF(IDWRITE .eq. 1) THEN
       filedat='../output/Stats/sm-     .dat'
       filedat1='../output/Stats/srms-     .dat'
       filedat2='../output/Stats/stau-     .dat'
       write(filedat(20:24),'(i5.5)') ns
       write(filedat1(22:26),'(i5.5)') ns
       write(filedat2(22:26),'(i5.5)') ns
       open(604,file=filedat,status='unknown')
       open(605,file=filedat1,status='unknown')
       open(606,file=filedat2,status='unknown')
       tkesum = 0.0d0
!   --------------------------------------------------------
       DO k=2,NZglobal-1
!   --------------------------------------
       yplus = Re*(2*k-3)/(2*NZglobal-4)
!   --------------------------------------------------------
!      Mean Velocity
       u1m = Su1(k)*npjx
       v1m = Sv1(k)*npjx
       w1m = Sw1(k)*npjx
       u2m = Su2(k)*npjx
       v2m = Sv2(k)*npjx
       w2m = Sw2(k)*npjx
       uvm = Suv(k)*npjx
       vwm = Svw(k)*npjx
       wum = Swu(k)*npjx
!    ------------------------------------------------------
       Tuvm = u1m*v1m - uvm
       Tvwm = v1m*w1m - vwm
       Twum = w1m*u1m - wum
       xrms = u2m - u1m*u1m
       yrms = v2m - v1m*v1m
       zrms = w2m - w1m*w1m
!      -----------------------------------------
!      Turbulent Kinetic Energy
       tkesum = tkesum + 0.5d0*(xrms + yrms + zrms)
!      -----------------------------------------
!      Velocity Fluctuation
       urms = dsqrt(xrms)
       vrms = dsqrt(yrms)
       wrms = dsqrt(zrms)
       write(604,6) k,yplus,u1m,v1m,w1m
       write(605,6) k,yplus,urms,vrms,wrms
       write(606,6) k,yplus,Tuvm,Tvwm,Twum
       ENDDO
!    ------------------------------------------------
!      Close FILE
       write(604,*) 'TKE =', tkesum/(NZglobal-2)
       close(604)
       close(605)
       close(606)

       print*,'   --- Plane Averaged Stats printed! ---',ns
       ENDIF
!    --------------------------------
#     ifdef KeyParallel
      call MPI_Barrier(comm3D,code)
#     endif
!    ------------------------------------------------
!     DEALLOCATE THE VARIABLES
       deallocate(Su1,Su2,Suv, &
                  Sv1,Sv2,Svw, &
                  Sw1,Sw2,Swu)
!    ------------------------------------------------
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: output_plane_sm'
         write(*,*) ''
#     endif

      RETURN
      END SUBROUTINE glob_put_plane_mm


