
!=====================================================================!
!---------------------------------------------------------------------!
!             SUBROUTINE LARGE EDDY SIMULATION (UNSTRUCTURED MESH)    !
!                              Xin Bai                                !
!                             Mar 2015                                !
!---------------------------------------------------------------------!
!=====================================================================!

!_____________________________________________________________________!
!      SUBROUTINES:     1.1)  sgsetup          			              !
!                       1.2)  sgmkwd                     	          !
!                       1.3)  sgmkeddy                      	      !
!                       1.4)  sgsetebc                   	          !
!                       1.5)  sgdamp_mason                    	      !
!—————————————————————————————————————————————————————————————————————!

!========================================================================!
!                        SUBROUNTINE: sgles                              !
!------------------------------------------------------------------------!
! This subroutine aims to initialize the variables used for LES by using !
! the existing velocity and grid information for original source code.   !
!========================================================================!
     subroutine sgles  (phiu,phiv,phiw,                    &
                        xc,yc,sig,dsig,No_cp,nbe,          &
                        xv,yv,sigv,dsigv,No_vp,nbev,       &
                        Hpr,h,eta,etan,                    &
                        Hprv,hv,etav)


#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE les
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE les
            implicit none
#     endif
!     =============== END ================
!     ====================================
!===============================================
      real*8, dimension(:,:):: phiu(N_CELL,NZ)
      real*8, dimension(:,:):: phiv(N_CELL,NZ)
      real*8, dimension(:,:):: phiw(N_CELL,NZ)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!================================================
!==============Local Variables===================
      real*8, dimension(:,:),allocatable :: dudx,dudy,dudz
      real*8, dimension(:,:),allocatable :: dvdx,dvdy,dvdz
      real*8, dimension(:,:),allocatable :: dwdx,dwdy,dwdz
      real*8, dimension(:,:),allocatable :: utau
      real*8 :: Sdudx,Sdudy,Sdudz
      real*8 :: Sdvdx,Sdvdy,Sdvdz
      real*8 :: Sdwdx,Sdwdy,Sdwdz
      real*8 :: sxx, sxy, sxz
      real*8 :: syy, syz, szz
!     -------------------------------------
      real*8 :: nut,DD,Cs,Sij,del,dist,lwplus
      real*8 :: yplusmax,utaumax
      real*8 :: maxeddy,maxeddy_gl
      integer:: IDisplay
!======================================================================
#  ifdef KeyDbg
     write(*,*) '====> Begin: Turbulence model'
     if (sg_model .eq. 1) then
     write(*,*)  'SG_MODEL = Standard Smagrinsky'
     else
     write(*,*) 'Attention: SG_MODEL = unknown!'
     endif
#  endif
!======================================================================
!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
#     ifdef KeyParallel
         if(rang_topo .eq. 0) then
           IDisplay = 1
         else
           IDisplay = 0
         endif
#    else
           IDisplay =  1
#    endif

      if(sg_init .eq. 1) then

!      ----------------------
!       allocate wall distance and eddy viscosity
        call alloc_les
!      -------------------------------------------
!      obtain wall normal distance
        call sgmkwd(xc,yc,sig,dsig,No_cp,nbe,      &
                    xv,yv,sigv,dsigv,No_vp,nbev,   &
                    Hpr,h,eta,etan,                &
                    Hprv,hv,etav)
!     ------------------------------------------
!      print LES parameters
        if(IDisplay .eq. 1) then
                print*, '        ====================================='
                print*, '                    LES Parameters'
                print*, '        -------------------------------------'
                print*, '        sg_smag_const =', sg_smag_const
                print*, '        sg_smag_aplus =', sg_smag_aplus
                print*, '        sg_smag_mason =', sg_smag_mason
                print*, '        sg_smag_so    =', sg_smag_so
                print*, '        von_Karman    =', von_Karman
                print*, '        ====================================='
        endif
!     -------------------------------------------
!      update initialization flag
        sg_init = 0
     endif

!*********************************************************************!
!                                                                     !
!                  Calculate the wall friction stress                 !
!                                                                     !
!*********************************************************************!
      allocate(utau(N_CELL,NZ))
!     ---------------------------------
!     For Open channel flow ONLY!
#     ifdef KeyTESTChannel
         do k=1,NZ
            do i=1,N_CELL
                utau(i,k) = 1.0d0
            enddo
         enddo
#     else
      call sgm_utau(utau,phiu,phiv,phiw,yplusmax,utaumax, &
                    xc,yc,sig,dsig,No_cp,nbe,  &
                    xv,yv,sigv,dsigv,No_vp,nbev)

       if(IDisplay .eq. 1) then
                print*, '        ====================================='
                print*, '                    LES Parameters'
                print*, '        -------------------------------------'
                print*, '        yplus =          ', yplusmax
                print*, '        utau  =          ', utaumax
                print*, '        ====================================='
       endif
!    ----------------------
!     u+ = u/u_tau,make first layer u+ = y+
!#    ifdef KeyTESTChannel
!       do k=2,NZ-1
!         do i=1,N_CELL0
!            phiu(i,2) = phiu(i,2)/dsqrt(2.0d0*phiu(i,2)/(dsig(1)*Re))
!            phiv(i,2) = phiv(i,2)/dsqrt(2.0d0*phiv(i,2)/(dsig(1)*Re))
!         enddo
!        enddo
!#    endif
#     endif
!*********************************************************************!
!                                                                     !
!                      Calculate the gradient                         !
!                                                                     !
!*********************************************************************!
     allocate(dudx(N_CELL,NZ),dudy(N_CELL,NZ),dudz(N_CELL,NZ), &
              dvdx(N_CELL,NZ),dvdy(N_CELL,NZ),dvdz(N_CELL,NZ), &
              dwdx(N_CELL,NZ),dwdy(N_CELL,NZ),dwdz(N_CELL,NZ))

     call fvm_gradient(phiu,phiv,phiw,       &
                       dudx,dudy,dudz,       &
                       dvdx,dvdy,dvdz,       &
                       dwdx,dwdy,dwdz,       &
                  xc,yc,sig,dsig,No_cp,nbe,  &
                  xv,yv,sigv,dsigv,No_vp,nbev)
!*********************************************************************!
!                                                                     !
!                  Calculate the eddy viscosity                       !
!                                                                     !
!*********************************************************************!
      Cs = sg_smag_const
      maxeddy = 0.0d0
      maxeddy_gl = 0.0d0

        DO k=2,NZ-1
         do i=1,N_CELL0
!           ---------------
!           Deardorff Mesh Scale
            DD  = Cs*(areaCell(i)*dsig(k)*Hpr(i))**(1.0d0/3.0d0)
            dist = walld(i,k)
!           ---------------
#           ifndef KeyLESVan
            lwplus = 0.0d0
#           else
            lwplus = dist*utau(i,k)*Re
#           endif
!           --------------------------
!           Near Wall Damping function
            call sgdamp(DD,dist,lwplus,sg_smag_mason,sg_smag_aplus,sg_smag_so,von_Karman)
!           ---------------
!           obtain gradient info
            Sdudx = dudx(i,k)
            Sdvdx = dvdx(i,k)
            Sdwdx = dwdx(i,k)
            Sdudy = dudy(i,k)
            Sdvdy = dvdy(i,k)
            Sdwdy = dwdy(i,k)
            Sdudz = dudz(i,k)
            Sdvdz = dvdz(i,k)
            Sdwdz = dwdz(i,k)
!           ---------------
!           Treat Eddy viscosity as scalar saved in center of the cell
            nut = Hpr(i)*DD
            sxx = (Sdudx+Sdudx)/2.0d0
            sxy = (Sdudy+Sdvdx)/2.0d0
            sxz = (Sdudz+Sdwdx)/2.0d0
            syy = (Sdvdy+Sdvdy)/2.0d0
            syz = (Sdvdz+Sdwdy)/2.0d0
            szz = (Sdwdz+Sdwdz)/2.0d0
            Sij = dsqrt(2.0d0*(sxx*sxx+syy*syy+szz*szz) &
                  + 4.0d0*(sxy*sxy+sxz*sxz+syz*syz))
!          --------------------------------------
!           Assign same value to all coefficient terms
            eddy(i,k) = nut*Sij
            if(eddy(i,k) .gt. maxeddy) maxeddy = eddy(i,k)
         enddo
      ENDDO
!    ---------------------------------------------
!        communicate the eddy visocosity
!        is this correct? or we should set zero on wall?
         call BCpcenter(eddy,xc,yc,sig,dsig,No_cp,nbe)

#        ifdef KeyParallel
           call MAX_parallel(maxeddy,maxeddy_gl)
           if(rang_topo .eq. 0) then
             print*, 'MAX EDDY VISCOSITY =', maxeddy_gl
           endif
#        endif
!    ---------------------------------------------
!        deallocate variables
         deallocate(dudx,dudy,dudz, &
                    dvdx,dvdy,dvdz, &
                    dwdx,dwdy,dwdz)

         deallocate(utau)
!    -----------------------------------------------
#     ifdef KeyDbg
           write(*,*) '<==== Exit: Turbulence Model.'
#     endif
!-------------------------------------------------
      return
      end subroutine sgles

!=============================================================================!
!-----------------------------------------------------------------------------!
!                              END SUBROUTINE sg_setup                        !
!-----------------------------------------------------------------------------!
!=============================================================================!
!  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW!
!  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW!
      SUBROUTINE fvm_gradient(phiu,phiv,phiw, &
                              dudx,dudy,dudz, &
                              dvdx,dvdy,dvdz, &
                              dwdx,dwdy,dwdz, &
                   xc,yc,sig,dsig,No_cp,nbe,  &
                   xv,yv,sigv,dsigv,No_vp,nbev)

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
!    |      Declaration of variables      |
!    |____________________________________|
!     --------------------------------------
      real*8, dimension(:,:) :: phiu(N_CELL,NZ)
      real*8, dimension(:,:) :: phiv(N_CELL,NZ)
      real*8, dimension(:,:) :: phiw(N_CELL,NZ)
!     -------------------------------------
      real*8, dimension(:,:) :: dudx(N_CELL,NZ)
      real*8, dimension(:,:) :: dudy(N_CELL,NZ)
      real*8, dimension(:,:) :: dudz(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdz(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdz(N_CELL,NZ)
!     --------------------------------------
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
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|
      real*8, dimension(:,:),allocatable :: Ui1,Ui2,Ui3,UiT,UiB
      real*8, dimension(:,:),allocatable :: Vi1,Vi2,Vi3,ViT,ViB
      real*8, dimension(:,:),allocatable :: Wi1,Wi2,Wi3,WiT,WiB
      real*8, dimension(:,:),allocatable :: phiuv,phivv,phiwv
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: fvm_gradient'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|
      allocate(Ui1(N_CELL0,NZ),         &
               Ui2(N_CELL0,NZ),         &
               Ui3(N_CELL0,NZ),         &
               UiT(N_CELL0,NZ),         &
               UiB(N_CELL0,NZ),         &
               Vi1(N_CELL0,NZ),         &
               Vi2(N_CELL0,NZ),         &
               Vi3(N_CELL0,NZ),         &
               ViT(N_CELL0,NZ),         &
               ViB(N_CELL0,NZ),         &
               Wi1(N_CELL0,NZ),         &
               Wi2(N_CELL0,NZ),         &
               Wi3(N_CELL0,NZ),         &
               WiT(N_CELL0,NZ),         &
               WiB(N_CELL0,NZ))
!    --------------------------------------------------
      allocate(phiuv(N_VERT,NZ-1),         &
               phivv(N_VERT,NZ-1),         &
               phiwv(N_VERT,NZ-1))
!    --------------------------------------------------
      call BCglobalVC(phiu,phiv,phiw,               &
                      phiuv,phivv,phiwv,             &
                      xc,yc,sig,dsig,No_cp,nbe,      &
                      xv,yv,sigv,dsigv,No_vp,nbev,2)
!*********************************************************************!
!                                                                     !
!                          face value calculation                     !
!                                                                     !
!*********************************************************************!

      do k=1,NZ
         do i=1,N_CELL0
            Ui1(i,k) = 0.0d0
            Ui2(i,k) = 0.0d0
            Ui3(i,k) = 0.0d0
            Vi1(i,k) = 0.0d0
            Vi2(i,k) = 0.0d0
            Vi3(i,k) = 0.0d0
            Wi1(i,k) = 0.0d0
            Wi2(i,k) = 0.0d0
            Wi3(i,k) = 0.0d0
            UiT(i,k) = 0.0d0
            UiB(i,k) = 0.0d0
            ViT(i,k) = 0.0d0
            ViB(i,k) = 0.0d0
            WiT(i,k) = 0.0d0
            WiB(i,k) = 0.0d0
         enddo
      enddo
!    -----------------------------------------------------------
!     GET THE FACE VALUE OF VELOCITY ON THE EDGE FACE
      call calcul_face(Ui1,Ui2,Ui3,UiT,UiB,phiu,phiuv,   &
                         xc,yc,sig,dsig,No_cp,nbe,  &
                         xv,yv,sigv,dsigv,No_vp,nbev,1)
!    -----------------------------------------------------------
      call calcul_face(Vi1,Vi2,Vi3,ViT,ViB,phiv,phivv,   &
                         xc,yc,sig,dsig,No_cp,nbe,  &
                         xv,yv,sigv,dsigv,No_vp,nbev,1)
!    -----------------------------------------------------------
      call calcul_face(Wi1,Wi2,Wi3,WiT,WiB,phiw,phiwv,   &
                         xc,yc,sig,dsig,No_cp,nbe,  &
                         xv,yv,sigv,dsigv,No_vp,nbev,1)
!    -----------------------------------------------------------
!    UPDATE THE GRADIENT USING THE FVM METHOD
      DO k=2,NZ-1
         do i=1,N_CELL0
!           ___________________________________________________
!           Horizontal neighbors
            dudx(i,k) = Ui1(i,k)*xnn(i,1)*dlVV(i,1) &
                      + Ui2(i,k)*xnn(i,2)*dlVV(i,2) &
                      + Ui3(i,k)*xnn(i,3)*dlVV(i,3)
            dvdx(i,k) = Vi1(i,k)*xnn(i,1)*dlVV(i,1) &
                      + Vi2(i,k)*xnn(i,2)*dlVV(i,3) &
                      + Vi3(i,k)*xnn(i,3)*dlVV(i,2)
            dwdx(i,k) = Wi1(i,k)*xnn(i,1)*dlVV(i,1) &
                      + Wi2(i,k)*xnn(i,2)*dlVV(i,2) &
                      + Wi3(i,k)*xnn(i,3)*dlVV(i,3)
            dudy(i,k) = Ui1(i,k)*ynn(i,1)*dlVV(i,1) &
                      + Ui2(i,k)*ynn(i,2)*dlVV(i,2) &
                      + Ui3(i,k)*ynn(i,3)*dlVV(i,3)
            dvdy(i,k) = Vi1(i,k)*ynn(i,1)*dlVV(i,1) &
                      + Vi2(i,k)*ynn(i,2)*dlVV(i,3) &
                      + Vi3(i,k)*ynn(i,3)*dlVV(i,2)
            dwdy(i,k) = Wi1(i,k)*ynn(i,1)*dlVV(i,1) &
                      + Wi2(i,k)*ynn(i,2)*dlVV(i,2) &
                      + Wi3(i,k)*ynn(i,3)*dlVV(i,3)
!           ---------------------------------
            dudx(i,k) = dudx(i,k)/areaCell(i)
            dvdx(i,k) = dvdx(i,k)/areaCell(i)
            dwdx(i,k) = dwdx(i,k)/areaCell(i)
            dudy(i,k) = dudy(i,k)/areaCell(i)
            dvdy(i,k) = dvdy(i,k)/areaCell(i)
            dwdy(i,k) = dwdy(i,k)/areaCell(i)
!          ___________________________________________________
!           Top and Bottom
!           assume uniform spacing in z direction
            dudz(i,k) = (UiT(i,k)-UiB(i,k))/dsig(k)
            dvdz(i,k) = (ViT(i,k)-ViB(i,k))/dsig(k)
            dwdz(i,k) = (WiT(i,k)-WiB(i,k))/dsig(k)
        enddo
      ENDDO
!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|

      deallocate(Ui1,Ui2,Ui3, &
                 Vi1,Vi2,Vi3, &
                 Wi1,Wi2,Wi3)

      deallocate(UiT,UiB, &
                 ViT,ViB, &
                 WiT,WiB)

      deallocate(phiuv,phivv,phiwv)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: fvm_gradient'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW!
!  WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW!
      SUBROUTINE sgm_utau(utau,phiu,phiv,phiw,yplusmax,utaumax, &
                        xc,yc,sig,dsig,No_cp,nbe,  &
                        xv,yv,sigv,dsigv,No_vp,nbev)

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
            USE les
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE les
            implicit none
#     endif
!     =============== END ================
!     ====================================
!     ____________________________________
!    |                                    |
!    |      Declaration of variables      |
!    |____________________________________|
!     --------------------------------------
      real*8, dimension(:,:) :: utau(N_CELL,NZ)
      real*8, dimension(:,:) :: phiu(N_CELL,NZ)
      real*8, dimension(:,:) :: phiv(N_CELL,NZ)
      real*8, dimension(:,:) :: phiw(N_CELL,NZ)
!     --------------------------------------
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
      real*8 :: yplusmax,utaumax
!    ----------------------------------------
!     local variables
      real*8, dimension(:),allocatable :: tauwT,tauwB
      real*8 :: umag,nnx,nny,nnz,dln
      integer :: elem,ii,jj,jv1,jv2
      real*8 :: xv1,yv1,xv2,yv2
      real*8 :: yplus,yplus_g,utau_g
!     ----------------
#     ifdef KeyParallel
      real*8, dimension(:,:),allocatable :: varu,varv,varw,vard
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: sgm_utau'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!*********************************************************************!
!                                                                     !
!                    stress on top and bottom walls                   !
!                                                                     !
!*********************************************************************!
        allocate(tauwT(N_CELL0),tauwB(N_CELL0))

          do i=1,N_CELL0
            tauwT(i) = 0.0d0
            tauwB(i) = 0.0d0
          enddo

          if(ZPB .ne. 1) then
            if(TopBC .eq. 0) then
              do i=1,N_CELL0
                  umag =dsqrt(phiu(i,NZ-1)*phiu(i,NZ-1)+phiv(i,NZ-1)*phiv(i,NZ-1))
                  tauwT(i) = 2.0d0*umag/dsig(NZ-1)
              enddo
            endif

            if(BotBC .eq. 0) then
              do i=1,N_CELL0
                  umag =dsqrt(phiu(i,2)*phiu(i,2)+phiv(i,2)*phiv(i,2))
                  tauwB(i) = 2.0d0*umag/dsig(2)
              enddo
            endif

          endif
!*********************************************************************!
!                                                                     !
!                  global velocity for parallel mode                  !
!                                                                     !
!*********************************************************************!
#        ifdef KeyParallel
            allocate(varu(N_CELL0global,NZglobal),varv(N_CELL0global,NZglobal),&
                varw(N_CELL0global,NZglobal),vard(N_CELL0global,NZglobal))

            call matgloCFull(phiu,varu)
            call matgloCFull(phiv,varv)
            call matgloCFull(phiw,varw)
            call matgloCFull(walld,vard)
#        endif
!*********************************************************************!
!                                                                     !
!                 update the friction stress                          !
!                                                                     !
!*********************************************************************!
           do k=2,NZ-1
              do i=1,N_CELL0
                  ii = utau_pair(i,k) ! minimun wd indicator
                  if(ii .eq. -2) then
                        utau(i,k) = tauwT(i) ! top
                  elseif(ii .eq. -1) then
                        utau(i,k) = tauwB(i) ! bot
                  else
                       j = utau_j(i,k) ! edge indicator
!                 -----------------
#                 ifndef KeyParallel
                       nnx = dxVV(ii,j)
                       nny = dyVV(ii,j)
                       dln = dsqrt(nnx*nnx+nny*nny)
                       nnx = phiu(ii,k)*nnx/dln
                       nny = phiv(ii,k)*nny/dln
                       nnz = phiw(ii,k)
                       umag = dsqrt(nnx*nnx+nny*nny+nnz*nnz)
                       dln = walld(ii,k)
#                 else
                       jj = j+1
                       if(jj .gt. 3) jj = jj-3
                       jv1 = No_vp_global(ii,j)
                       jv2 = No_vp_global(ii,jj)
                       xv1 = xv_global(jv1)
                       yv1 = yv_global(jv1)
                       xv2 = xv_global(jv2)
                       yv2 = yv_global(jv2)
                       nnx = xv2-xv1
                       nny = yv2-yv1
                       dln = dsqrt(nnx*nnx+nny*nny)
                       nnx = varu(ii,k)*nnx/dln
                       nny = varv(ii,k)*nny/dln
                       nnz = varw(ii,k)
                       umag = dsqrt(nnx*nnx+nny*nny+nnz*nnz)
                       dln = vard(ii,k)
#                 endif
                       utau(i,k) = umag/dln
                 endif
              enddo
            enddo
!*********************************************************************!
!                                                                     !
!                 update utau = dsqrt(tau_w/rho)                      !
!                                                                     !
!*********************************************************************!
          utaumax = 0.0d0
          do k=2, NZ-1
             do i=1,N_CELL0
                utau(i,k) = dsqrt(utau(i,k)/Re)
                if(utau(i,k) .gt. utaumax) then
                  utaumax = utau(i,k)
                endif
             enddo
          enddo
!#        ifdef KeyParallel
!        ---------------------------------------------
!         check utau
!          do k=2,NZ-1
!              print*,'z',sig(k),'wd',walld(1,k),'utau',utau(1,k),'procs',rang_topo
!          enddo
!          stop
!#        endif
!*********************************************************************!
!                                                                     !
!                 check first layer yplus value                       !
!                                                                     !
!*********************************************************************!
          yplusmax = 0.0d0
!    --------------------
!    check vertical direction
!    ---------
!    Bot
           k = 2
            do i=1,N_CELL0
              if(utau_pair(i,k) .eq. -1) then
                 yplus = walld(i,k)*utau(i,k)*Re
                if(yplus .gt. yplusmax) then
                  yplusmax =  yplus
                endif
              endif
            enddo
!    --------
!    Top
          k = NZ-1
            do i=1,N_CELL0
              if(utau_pair(i,k) .eq. -2) then
                  yplus = walld(i,k)*utau(i,k)*Re
                  if(yplus .gt. yplusmax) then
                   yplusmax = yplus
                 endif
              endif
            enddo
!    ---------------------------
!     check xy direction
#     ifndef KeyParallel
            do k=2, NZ-1
              do i=1,N_CELL0
              if(utau_pair(i,k) .gt. 0) then
                elem = utau_pair(i,k)
                yplus = walld(elem,k)*utau(elem,k)*Re
                  if(yplus .gt. yplusmax) then
                    yplusmax = yplus
                 endif
              endif
             enddo
            enddo
#    else
             do k=2, NZ-1
               do i=1,N_CELL0
                 if(utau_pair(i,k) .gt. 0) then
                   elem = utau_pair(i,k)
!                  --------------------------------
!                  assume tau(i,k) = tau(elem,k)
                   yplus = vard(elem,k)*utau(i,k)*Re
                     if(yplus .gt. yplusmax) then
                      yplusmax = yplus
                     endif
                 endif
               enddo
             enddo
!            ---------------------------------
!            communicate yplusmax
             yplus_g = 0.0d0
             utau_g = 0.0d0
             call MPI_Barrier(comm3D,code)
             call MAX_Parallel(yplusmax,yplus_g)
             call MAX_Parallel(utaumax,utau_g)
             yplusmax = yplus_g
             utaumax = utau_g
#    endif
!      ________________________________________________________
!     |                                                        |
!     |                         Deallocate                     |
!     |________________________________________________________|
         deallocate(tauwT,tauwB)

#     ifdef KeyParallel
         deallocate(varu,varv,varw,vard)
#     endif

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: sgm_utau'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
!========================================================================!
!                        SUBROUNTINE: sgmkwd                             !
!------------------------------------------------------------------------!
! This subroutine is for set the wall distance for each cell and will be !
! used in the damping function to calculate the turbulent viscosity      !
!------------------------------------------------------------------------!
!          NOTE: THE PROGRAMME HAVEN'T BEEN PARRALLELIZED YET.           !
!========================================================================!

     subroutine sgmkwd (xc,yc,sig,dsig,No_cp,nbe,      &
                        xv,yv,sigv,dsigv,No_vp,nbev,   &
                        Hpr,h,eta,etan,                &
                        Hprv,hv,etav)


#     include "cppdefs.h"
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
            USE geometry
            USE les
            implicit none
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            USE les
            implicit none
#     endif
!     =============== END ================
!     ====================================

!===============================================
!---------------Input Variables----------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!================================================
!==============Local Variables===================
      real*8, dimension(:),allocatable :: dnminx
      real*8, dimension(:),allocatable :: dnminy
      real*8, dimension(:),allocatable :: dnminxy
      integer, dimension(:),allocatable :: u_pair1
      integer, dimension(:),allocatable :: u_pair2
      integer, dimension(:),allocatable :: u_pair3
      integer, dimension(:),allocatable :: uj1,uj2,uj3
      integer, dimension(:,:),allocatable :: TagXC_global
      integer, dimension(:,:),allocatable :: TagYC_global
      real*8  :: dnmint,dnminb,xee,yee
      real*8:: dist,dlx,dly,dls
      real*8:: dnmin1,dnmin2,dnmin
      integer :: ii,jj,jv1,jv2,jv3,jc,elem
!     -------------------------------------------
#  ifdef KeyDbg
     write(*,*) '====> Begin: Set Wall distance!'
#  endif

!   -----------------------------------------------
!   For open channel flow ONLY!
#   ifdef KeyTESTChannel
      do k=1,NZ
         do i=1,N_CELL
            walld(i,k) = dabs(sig(k))
         enddo
      enddo
#   endif
!!------------------------------------------------!
!!       Find wall distance in XY plane           !
!!------------------------------------------------!
!
!      allocate(dnminx(N_CELL),dnminy(N_CELL),dnminxy(N_CELL))
!      allocate(u_pair1(N_CELL),u_pair2(N_CELL),u_pair3(N_CELL))
!      allocate(uj1(N_CELL),uj2(N_CELL),uj3(N_CELL))
!
!       do i=1,N_CELL
!         dnminx(i) = 1.0E5
!         dnminy(i) = 1.0E5
!         u_pair1(i) = 0
!         u_pair2(i) = 0
!         u_pair3(i) = 0
!         uj1(i) = 0
!         uj2(i) = 0
!         uj3(i) = 0
!       enddo
!!     ---------------------
!      if(XPB .ne. 1) then ! no periodic in X direction
!         if(XBC .eq. 0) then ! no slip wall
!             do i=1, N_CELL
!               dnmin1 = dabs(XDIni-xc(i))
!               dnmin2 = dabs(XDFin-xc(i))
!               if(dnmin1 .lt. dnmin2) then
!                 dist = dnmin1
!               else
!                 dist = dnmin2
!               endif
!                dnminx(i) = dist
!             enddo
!         endif
!      endif
!!    ----------------------
!      if(YPB .ne. 1) then ! no periodic in Y direction
!         if(YBC .eq. 0) then ! no slip wall
!             do i=1, N_CELL
!               dnmin1 = dabs(YDIni-yc(i))
!               dnmin2 = dabs(YDFin-yc(i))
!               if(dnmin1 .lt. dnmin2) then
!                 dist = dnmin1
!               else
!                 dist = dnmin2
!               endif
!                dnminy(i) = dist
!             enddo
!         endif
!      endif
!!    -----------------------
!            do i=1,N_CELL
!               if(dnminx(i) .lt. dnminy(i)) then
!                   dnminxy(i) = dnminx(i)
!               else
!                   dnminxy(i) = dnminy(i)
!               endif
!            enddo
!
!!------------------------------------------------!
!!       Find wall distance in vertical           !
!!------------------------------------------------!
!         do k=1, NZ
!           dnmint = 1.0E5
!           dnminb = 1.0E5
!!           ---------------------
!           if(ZPB .ne. 1) then
!              if(TopBC .eq. 0) then
!              dnmint = dabs(sig(k)-sigFin)
!              endif
!!           ---------------------
!              if(BotBC .eq. 0) then
!              dnminb = dabs(sig(k)-sigIni)
!              endif
!           endif
!!           ---------------------
!            if(dnmint .lt. dnminb) then
!              dls = dnmint
!            else
!              dls = dnminb
!            endif
!
!              do i=1, N_CELL
!                if(dnminxy(i) .lt. dls) then
!                   walld(i,k) = dnminxy(i)
!                   utau_pair(i,k) = 1  ! minmum wd in xy dir
!                else
!                   walld(i,k) = dls*Hpr(i)
!                   if(dnmint .lt. dnminb) then
!                      utau_pair(i,k) = -2 ! minmum wd in top z dir
!                   else
!                      utau_pair(i,k) = -1 ! minmum wd in bot z dir
!                   endif
!                endif
!              enddo
!         enddo
!!------------------------------------------------!
!!       assign global bc tag for parallel mode   !
!!------------------------------------------------!
!#  ifdef KeyParallel
!        allocate(TagXC_global(N_CELL0global,2), &
!                 TagYC_global(N_CELL0global,2))
!
!        do i=1, N_CELL0global
!           do j=1,2
!             TagXC_global(i,j) = 0
!             TagYC_global(i,j) = 0
!           enddo
!        enddo
!
!        do i=1,N_CELL0
!          if(nbe(i) .ne. 0) then
!           elem = index_global(i)
!           do j=1,3
!
!             jc = No_cp(i,j)
!             if(TagXC(jc) .gt. 0) then
!               TagXC_global(elem,1) = TagXC(jc)
!               TagXC_global(elem,2) = j
!             endif
!
!             if(TagYC(jc) .gt. 0) then
!               TagYC_global(elem,1) = TagYC(jc)
!               TagYC_global(elem,2) = j
!             endif
!
!            if(TagXC(jc)+TagYC(jc) .eq. 0) then
!                if(index_global(jc) .lt. 0) then ! ghost cell
!                  TagXC_global(elem,2) = j
!                endif
!            endif
!
!           enddo
!          endif
!        enddo
!!  --------------------------------------
!!   commmunicate
!        call MPI_Barrier(comm3D,code)
!        call SUM_Paralleli(TagXC_global(1:N_CELL0global,1),N_CELL0global)
!        call SUM_Paralleli(TagXC_global(1:N_CELL0global,2),N_CELL0global)
!        call SUM_Paralleli(TagYC_global(1:N_CELL0global,1),N_CELL0global)
!        call SUM_Paralleli(TagYC_global(1:N_CELL0global,2),N_CELL0global)
!#  endif
!!------------------------------------------------!
!!                  update u_pair                 !
!!------------------------------------------------!
!! case 1: minmum wd in cylinder
!          if(ICYLINDER .eq. 1) then
!              do i=1, N_CELL0
!                dnmin = 1.0E5
!!               -----------------
!#               ifndef KeyParallel
!                do ii=1, N_CELL0
!                   if(nbe(ii) .ne. 0) then
!                   do j=1,3
!                     jc = No_cp(ii,j)
!                     if(jc .gt. N_CELL0) then
!                        if(TagXC(jc)+TagYC(jc) .eq. 0) then
!                           dlx = xc(i)-xc(ii)
!                           dly = yc(i)-yc(ii)
!                           dls = dsqrt(dlx*dlx+dly*dly)
!                           if(dls .lt. dnmin) then
!                             u_pair1(i) = ii
!                             uj1(i) = j
!                             dnmin = dls
!                           endif
!                        endif
!                     endif
!                    enddo
!                    endif
!                 enddo
!#                else
!!               ---------------------
!                do ii=1, N_CELL0global
!                   if(nbe_global(ii) .ne. 0) then
!                       if(TagXC_global(ii,1)+TagYC_global(ii,1) .eq. 0) then
!                           jv1 = No_vp_global(ii,1)
!                           jv2 = No_vp_global(ii,2)
!                           jv3 = No_vp_global(ii,3)
!                           xee=(xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
!	                   yee=(yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
!                           dlx = xc(i)-xee
!                           dly = yc(i)-yee
!                           dls = dsqrt(dlx*dlx+dly*dly)
!                           if(dls .lt. dnmin) then
!                             u_pair1(i) = ii
!                             uj1(i) = TagXC_global(ii,2)
!                             dnmin = dls
!                          endif
!                       endif
!                   endif
!                 enddo
!#                endif
!               enddo
!          endif
!! case 2: mimum wd in x dir
!          if(XPB .ne. 1) then
!            if(XBC .eq. 0) then
!              do i=1, N_CELL0
!                dnmin = 1.0E5
!!               -----------------
!#               ifndef KeyParallel
!                do ii=1, N_CELL0
!                   if(nbe(ii) .ne. 0) then
!                   do j=1,3
!                     jc = No_cp(ii,j)
!                     if(jc .gt. N_CELL0) then
!                        if(TagXC(jc) .gt. 0) then
!                           dlx = xc(i)-xc(ii)
!                           dly = yc(i)-yc(ii)
!                           dls = dsqrt(dlx*dlx+dly*dly)
!                           if(dls .lt. dnmin) then
!                             u_pair2(i) = ii
!                             uj2(i) = j
!                             dnmin = dls
!                           endif
!                        endif
!                     endif
!                    enddo
!                    endif
!                 enddo
!#                else
!!               ---------------------
!              do ii=1, N_CELL0global
!                   if(nbe_global(ii) .ne. 0) then
!                       if(TagXC_global(ii,1) .gt. 0) then
!                           jv1 = No_vp_global(ii,1)
!                           jv2 = No_vp_global(ii,2)
!                           jv3 = No_vp_global(ii,3)
!                           xee=(xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
!	                   yee=(yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
!                           dlx = xc(i)-xee
!                           dly = yc(i)-yee
!                           dls = dsqrt(dlx*dlx+dly*dly)
!                           if(dls .lt. dnmin) then
!                             u_pair2(i) = ii
!                             uj2(i) = TagXC_global(ii,2)
!                             dnmin = dls
!                          endif
!                       endif
!                   endif
!                 enddo
!#                endif
!               enddo
!             endif
!          endif
!! case 3: mimum wd in y dir
!          if(YPB .ne. 1) then
!            if(YBC .eq. 0) then
!              do i=1, N_CELL0
!                dnmin = 1.0E5
!!               -----------------
!#               ifndef KeyParallel
!                do ii=1, N_CELL0
!                   if(nbe(ii) .ne. 0) then
!                   do j=1,3
!                     jc = No_cp(ii,j)
!                     if(jc .gt. N_CELL0) then
!                        if(TagYC(jc) .gt. 0) then
!                           dlx = xc(i)-xc(ii)
!                           dly = yc(i)-yc(ii)
!                           dls = dsqrt(dlx*dlx+dly*dly)
!                           if(dls .lt. dnmin) then
!                             u_pair3(i) = ii
!                             uj3(i) = j
!                             dnmin = dls
!                           endif
!                        endif
!                     endif
!                    enddo
!                    endif
!                 enddo
!#                else
!!               ---------------------
!              do ii=1, N_CELL0global
!                   if(nbe_global(ii) .ne. 0) then
!                       if(TagYC_global(ii,1) .gt. 0) then
!                           jv1 = No_vp_global(ii,1)
!                           jv2 = No_vp_global(ii,2)
!                           jv3 = No_vp_global(ii,3)
!                           xee=(xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
!	                   yee=(yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
!                           dlx = xc(i)-xee
!                           dly = yc(i)-yee
!                           dls = dsqrt(dlx*dlx+dly*dly)
!                           if(dls .lt. dnmin) then
!                             u_pair3(i) = ii
!                             uj3(i) = TagYC_global(ii,2)
!                             dnmin = dls
!                          endif
!                       endif
!                   endif
!                 enddo
!#                endif
!               enddo
!             endif
!          endif
!!------------------------------------------------!
!!                update utau_pair                !
!!------------------------------------------------!
!          do k=1,NZ
!            do i=1,N_CELL0
!              if(utau_pair(i,k) .gt. 0) then
!                dist = walld(i,k)
!                dlx = dnminx(i)
!                dly = dnminy(i)
!                if(dabs(dlx-dist) .lt. 1.0E-7) then
!                    utau_pair(i,k) = u_pair2(i) ! wd in x dir
!                    utau_j(i,k) = uj2(i)
!                elseif(dabs(dly-dist) .lt. 1.0E-7) then
!                    utau_pair(i,k) = u_pair3(i) ! wd in y dir
!                    utau_j(i,k) = uj3(i)
!                else
!                    utau_pair(i,k) = u_pair1(i)
!                    utau_j(i,k) = uj1(i)
!                endif
!               endif
!             enddo
!           enddo
!!----------------------------------------------------
!!  final check
!         do k=1,NZ
!           do i=1,N_CELL0
!              if(utau_pair(i,k) .gt. 0) then
!                 if(utau_j(i,k) .eq. 0) then
!                  print*, 'FATAL ERROR! UTAU PAIR NOT ASSIGNED'
!                  stop
!                  endif
!              endif
!           enddo
!          enddo
!! ---------------------------------------------------
!!  explicit check for two test cases
!#   ifdef KeyTESTChannel
!        do k=1,NZ
!          do i=1,N_CELL0
!            if(utau_pair(i,k) .ne. -1) then
!                print*, 'ERROR! IN CHANNEL TEST WD SHOULD BE IN VERTICAL'
!                stop
!            endif
!          enddo
!        enddo
!#   endif
!!   ------------------------------------------------
!!   for cylinder case, if top and bottom are free slip walls
!#   ifdef KeyTESTCylinder
!        do k=1,NZ
!           do i=1, N_CELL0
!            if(utau_pair(i,k) .lt. 0) then
!               print*, 'ERROR! SHOULD BE NO VERTICAL WALLS'
!               stop
!            endif
!           enddo
!        enddo
!#   endif
!!   -------------------------------------------------
!!   deallocate variables
!#       ifdef KeyParallel
!          deallocate(TagXC_global,TagYC_global)
!#       endif
!!  -------------------------------------------------
!        deallocate(dnminxy,dnminx,dnminy)
!        deallocate(u_pair1,u_pair2,u_pair3)
!        deallocate(uj1,uj2,uj3)
!----------------------------------------------------
#  ifdef KeyDbg
     write(*,*) '<==== Exit: Set Wall distance!'
#  endif

     return
     end subroutine sgmkwd

!=============================================================================!
!-----------------------------------------------------------------------------!
!                              END SUBROUTINE sg_mkwd                         !
!-----------------------------------------------------------------------------!
!=============================================================================!

!========================================================================!
!                        SUBROUNTINE: sgdamp                             !
!------------------------------------------------------------------------!
! This subroutine is to apply the near wall damping function.            !
!                  a) Mason Damping                                      !
!                  b) Van-Driest Damping                                 !
!------------------------------------------------------------------------!
!========================================================================!
   subroutine sgdamp(del,dist,yplus,                     &
!               ----------------------------------------
                mason,aplus,so,                          &
!               ----------------------------------------
                Karman)

     implicit none
!    ---------------------------------------------------
!    ==================Variable Declaration=============
     real*8 :: del,dist,mason,aplus,so,Karman,yplus
     real*8 :: lwall,n,dn1,dn2
!    ___________________________________________________
!    Mason Damping function
#    ifndef KeyLESVan
     lwall = Karman*dist
     n = mason
     dn1 = (1.0d0/del)**n +(1.0d0/lwall)**n
     dn2 = 1.0d0/dn1
     del = dn2**(1.0d0/n)
!    ____________________________________________________
!    Van Driest Damping function
!    l_sg = Cs*delta*f_damping
#   else
         lwall = Karman*dist ! l_wall*karman
         if(lwall .gt. del) then
           lwall = del
         endif
         del = lwall*(1.0d0-dexp(-yplus/aplus))
#   endif
!    _____________________________________________________
!    (Delta*const)^2
     del = del**2

     return
     endsubroutine sgdamp
!=============================================================================!
!-----------------------------------------------------------------------------!
!                              END SUBROUTINE sg_mkwd                         !
!-----------------------------------------------------------------------------!
!=============================================================================!

