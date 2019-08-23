!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                       Calculate the CFL number                      !
!                             March 2016                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
    SUBROUTINE glob_cfl(phiu,phiv,phiw,phip,      &
                        xc,yc,sig,dsig,No_cp,nbe)

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
!    |      Declaration of variables      |
!    |____________________________________|

      real*8, dimension(:,:):: phiu(N_CELL,NZ)
      real*8, dimension(:,:):: phiv(N_CELL,NZ)
      real*8, dimension(:,:):: phiw(N_CELL,NZ)
      real*8, dimension(:,:):: phip(N_CELL,NZ)
!     --------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     -------------------------------------
!     LOCAL VARIABLES
      real*8 :: cfl_x,cfl_y,cfl_z,Vol
      real*8 :: cfl_x_max,cfl_y_max,cfl_z_max
      real*8 :: umax,vmax,wmax
      real*8 :: sumu,sumv,sumw
      real*8 :: aveu,avev,avew
      real*8 :: sumug,sumvg,sumwg
      integer :: IDISPLAY,irec
      character*50 filen
!     ---------------------------------------
      real*8 :: flux,flux_m,flux_l,flux_g,flux_mg
      real*8 :: flux_bc,flux_bcg
      real*8 :: flux_i, flux_o,flux_ig,flux_og
      real*8 :: flux_w,flux_wg
      real*8 :: flux_cfl,cfl,cfl_max
      real*8, dimension(:) :: uface(1:3)
      integer :: jc,elem
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: glob_cfl'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyParallel
      if(rang_topo .eq. 0) then
         IDISPLAY = 1
      else
         IDISPLAY = 0
      endif
#     else
         IDISPLAY = 1
#     endif
!    ------------------------------------------
!     velocity check
      umax = 0.
      vmax = 0.
      wmax = 0.
      sumu = 0.
      sumv = 0.
      sumw = 0.
      sumug = 0.
      sumvg = 0.
      sumwg = 0.
      aveu = 0.
      avev = 0.
      avew = 0.
!    -----------------------------------------
!    ONLY CHECK THE INNER CELLS
      DO k=2,NZ-1
        DO i=1,N_CELL0
!      ---------------------------------------
!         u field
            if(dabs(phiu(i,k)) .gt. dabs(umax)) umax = phiu(i,k)
!      ---------------------------------------
!         v field
            if(dabs(phiv(i,k)) .gt. dabs(vmax)) vmax = phiv(i,k)
!      ----------------------------------------
!        w field
            if(dabs(phiw(i,k)) .gt. dabs(wmax)) wmax = phiw(i,k)
!      ----------------------------------------
!      get average v
            sumu = phiu(i,k)*areaCell(i)*dsig(k) + sumu
            sumv = phiv(i,k)*areaCell(i)*dsig(k) + sumv
            sumw = phiw(i,k)*areaCell(i)*dsig(k) + sumw
        ENDDO
      ENDDO

!    ---------------------------------------------------
!    check the divergence free condition
         cfl = 0.0d0
         flux_l = 0.0d0
         flux_g = 0.0d0
         flux_m = 0.0d0
         flux_mg = 0.0d0
         flux_bc = 0.0d0
         flux_bcg = 0.0d0
         flux_i = 0.0d0
         flux_ig = 0.0d0
         flux_o = 0.0d0
         flux_og = 0.0d0
         flux_w = 0.0d0
         flux_wg = 0.0d0
!   ----------------------------------------------------
!   Loop over all interal cells
         do k=2,NZ-1
          do i=1,N_CELL0
             flux = 0.0d0
             flux_cfl = 0.0d0
             uface(1) = U1FACE(i,k)
             uface(2) = U2FACE(i,k)
             uface(3) = U3FACE(i,k)
             do j=1,3
                flux = flux + uface(j)*dlVV(i,j)*dsig(k-1)
                flux_cfl = flux_cfl + dabs(uface(j)*dlVV(i,j)*dsig(k))
                jc = No_cp(i,j)
!               ---------------------------------------
#                 ifndef KeyParallel
	          if (jc.lt.1.OR.jc.gt.N_CELL0) then
!                 ====================================
!                 =====  START PARALLEL OPTION =======
#                 else
                  elem = index_global(jc)
	          if (elem.lt.1.OR.elem.gt.N_CELL0global) then
#                 endif
                  flux_w = flux_w + uface(j)*dlVV(i,j)*dsig(k)
              endif
             enddo

                flux = flux + UTFACE(i,k)*areaCell(i)
                flux = flux + UBFACE(i,k)*areaCell(i)

                flux_cfl = flux_cfl + dabs(UBFACE(i,k)*areaCell(i)) &
                           + dabs(UTFACE(i,k)*areaCell(i))

                flux_cfl = 0.5d0*dt*flux_cfl/(areaCell(i)*dsig(k))

                if(flux_cfl .gt. cfl) cfl = flux_cfl

                flux_l = flux_l + flux

                if(dabs(flux_m) .lt. dabs(flux)) then
                   flux_m = flux
                endif

                if(k .eq. 2) then
                    if(bcbot .eq. 1) then
                        flux_i = flux_i + UBFACE(i,k)*areaCell(i)
                    endif
                endif

                if(k .eq. NZ-1) then
                    if(bctop .eq. 1) then
                        flux_o = flux_o + UTFACE(i,k)*areaCell(i)
                    endif
                endif

          enddo
         enddo
#     ifndef KeyParallel
          Vol  = TotalVolume*(SigFin-SigIni)
          aveu = sumu/Vol
          avev = sumv/Vol
          avew = sumw/Vol
          cfl_max   = cfl
          cfl_x_max = umax
          cfl_y_max = vmax
          cfl_z_max = wmax
#     else
          Vol  = TotalVolume*(SigFin-SigIni)
          cfl_x = dabs(umax)
          cfl_y = dabs(vmax)
          cfl_z = dabs(wmax)
          call MPI_Barrier(comm3D,code)
          call MAX_parallel(cfl,cfl_max)
          call MAX_parallel(cfl_x,cfl_x_max)
          call MAX_parallel(cfl_y,cfl_y_max)
          call MAX_parallel(cfl_z,cfl_z_max)
          call SUM_parallel(sumu,sumug)
          call SUM_parallel(sumv,sumvg)
          call SUM_parallel(sumw,sumwg)
          aveu = sumug/Vol
          avev = sumvg/Vol
          avew = sumwg/Vol
#     endif

!    -----------------
#    ifndef KeyParallel
           flux_g   = flux_l
           flux_bcg = flux_i+flux_o+flux_w
           flux_ig  = flux_i
           flux_og  = flux_o
           flux_wg  = flux_w
#    else
         call MPI_Barrier(comm3D,code)
         call SUM_parallel(flux_l,flux_g)
         call SUM_parallel(flux_i,flux_ig)
         call SUM_parallel(flux_o,flux_og)
         call SUM_parallel(flux_w,flux_wg)
         call MAX_parallel(flux_m,flux_mg)
         flux_m = flux_mg
         flux_bcg = flux_ig+flux_og+flux_wg
#    endif

        if(IDISPLAY .eq. 1) then
            print*, '        ====================================='
            print*, '        MAX U     =',cfl_x_max
            print*, '        MAX V     =',cfl_y_max
            print*, '        MAX W     =',cfl_z_max
            print*, '        MAX CFL   =',cfl_max
            print*, '        ====================================='
            irec=60
            filen='meanuv.dat'
            open(irec,file=filen,position="append")
            write(irec,'(f12.5,2x,f12.5,2x,f12.5,2x,f12.5)') time,aveu,avev,avew
            close(irec)
        endif

        if(IDISPLAY .eq. 1) then
            print*, '        ====================================='
            print*, '                      MASS FLUX'
            print*, '        -------------------------------------'
            print*, '        Total  flux =', flux_g
            print*, '        Maxim  flux =', flux_m
            print*, '        BC     flux =', flux_bcg
            print*, '        BOT    flux =', flux_ig
            print*, '        TOP    flux =', flux_og
            print*, '        Side   flux =', flux_wg
            print*, '        ====================================='
        endif
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: glob_cfl'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                       Check Periodic BC                             !
!                             March 2016                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
    SUBROUTINE check_pbc(xc,yc,sig,dsig,No_cp,nbe,      &
                         xv,yv,sigv,dsigv,No_vp,nbev)
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
!    |      Declaration of variables      |
!    |____________________________________|
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer,dimension(:)  :: nbev(N_VERT)
!     -------------------------------------
!     LOCAL VARIABLES
      real*8,dimension(:,:) :: phix(N_CELL,NZ)
      real*8,dimension(:,:) :: phiy(N_CELL,NZ)
      integer:: irec,elem
      character*50 filen
!     -------------------------------------
!     Assign initial value equal to cell centre index
      Do k=1,NZ
        do i=1,N_CELL
           phix(i,k) = 0.
           phiy(i,k) = 0.
        enddo
      EndDo


      Do k=2,NZ-1
        do i=1,N_CELL0
           phix(i,k) = xc(i)
           phiy(i,k) = yc(i)
        enddo
      EndDo

#     ifdef KeyParallel
       call BCpcenter(phix,xc,yc,sig,dsig,No_cp,nbe)
       call BCpcenter(phiy,xc,yc,sig,dsig,No_cp,nbe)

       irec=60
       filen='../output/Parallel/tag_  .txt'
       write(filen(24:25),'(i2.2)') rang_topo+1
       open(irec,file=filen)
        k = NZ
        Do nc = N_CELL0+1,N_CELL
            elem = index_global(nc)
            if(elem .gt. 0) then
                write(irec,*),'extra cell',nc,elem,xc_global(elem),yc_global(elem),phix(nc,k),phiy(nc,k)
            else
                write(irec,*),'periodic cell',nc,No_cp(nc,3),xc(nc),yc(nc),phix(nc,k),phiy(nc,k)
            endif
        EndDo

       close(irec)
#     endif


      RETURN
      END
