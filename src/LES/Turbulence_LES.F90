!=====================================================================!
!---------------------------------------------------------------------!
!          SUBROUTINE LARGE EDDY SIMULATION (UNSTRUCTURED MESH)       !
!                        Xin Bai   &   Miguel Uh                      !
!                        Mar 2015      Oct 2017                       !
!---------------------------------------------------------------------!
!=====================================================================!
!        NOTE: THE PROGRAMS HAVEN'T BEEN PARRALLELIZED YET.           !
!=====================================================================!

!_____________________________________________________________________!
!                                                                     !
!      MAIN SUBROUTINES:  sgles                                       !
!      AUX. SUBROUTINES:  1.1)  sgmkwd                                !
!                         1.2)  sgdamp                                !
!                                                                     !
!—————————————————————————————————————————————————————————————————————!

!=====================================================================!
!                                                                     !
!                        SUBROUNTINE: sgles                           !
!                                                                     !
!=====================================================================!
!                                                                     !
!   This subroutine aims to initialize the variables used for LES     !
!   by using the existing velocity and grid information for original  !
!   source code.                                                      !
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
!   _______________________________________________________________   !
!  |     Name               |    Size     |   Description          |  !
!  |________________________|_____________|________________________|  !
!  | <--- Gamux,Gamuy,Gamuz |(N_CELL,NZ)  | diffusion coeff. u     |  !
!  | <--- Gamvx,Gamvy,Gamvz |(N_CELL,NZ)  | diffusion coeff. v     |  !
!  | <--- Gamwx,Gamwy,Gamwz |(N_CELL,NZ)  | diffusion coeff. w     |  !
!  |________________________|_____________|________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |    Name    |    Size     |  Description                       |  !
!  |____________|_____________|____________________________________|  !
!  | --> uf     |(N_CELL,NZ)  | u solution at the current RK step  |  !
!  | --> vf     |(N_CELL,NZ)  | v solution at the current RK step  |  !
!  | --> wf     |(N_CELL,NZ)  | w solution at the current RK step  |  !
!  |____________|_____________|____________________________________|  !
!  | -->Hpr_gam |(N_CELL)     | Water depth cell-center (ponderate)|  !
!  | -->Hpr_new |(N_CELL)     | Water depth cell-center (new)      |  !
!  | -->Hpr     |(N_CELL)     | Water depth cell-center (old)      |  !
!  | -->Hprv    |(N_VERT)     | Water depth solution vertex        |  !
!  |____________|_____________|____________________________________|  !
!  | --> h      |(N_CELL)     | Fixed water depth cell-center      |  !
!  | --> hv     |(N_VERT)     | Fixed water depth vertex           |  ! 
!  |____________|_____________|____________________________________|  !
!  | --> xc,yc  |(N_CELL)     | Coordinates of the cell centers    |  !
!  | --> sig    |(NZ)         | Sigma value at the cell centers    |  !
!  | --> dsig   |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | --> No_cp  |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | --> nbe    |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | --> xv,yv  |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | --> sigv   |(NZ-1)       | sigma of the vertex points         |  !
!  | --> dsigv  |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !
!  | --> No_vp  |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | --> nbev   |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
!---------------------------------------------------------------------!

      SUBROUTINE sgles (Gamux,Gamuy,Gamuz,                 &
                        Gamvx,Gamvy,Gamvz,                 &
                        Gamwx,Gamwy,Gamwz,                 &
!                       ------------------------------------
                        uf,vf,wf,                          &
!                       ------------------------------------
                        Hpr_gam,Hpr_new,Hpr,Hprv,          &
!                       ------------------------------------
                        h,hv,                              &
!                       ------------------------------------
                        xc,yc,sig,dsig,No_cp,nbe,          &
                        xv,yv,sigv,dsigv,No_vp,nbev,       &
!                      ------------------------------------
                        No_wb)

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

      real*8, dimension(:,:)  :: Gamux(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamuy(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamuz(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamvx(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamvy(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamvz(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamwx(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamwy(N_CELL,NZ)
      real*8, dimension(:,:)  :: Gamwz(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)    :: Hpr_gam(N_CELL)
      real*8, dimension(:)    :: Hpr_new(N_CELL)
      real*8, dimension(:)    :: Hpr(N_CELL)
      real*8, dimension(:)    :: Hprv(N_VERT)
      real*8, dimension(:,:)  :: uf(N_CELL,NZ)
      real*8, dimension(:,:)  :: vf(N_CELL,NZ)
      real*8, dimension(:,:)  :: wf(N_CELL,NZ)
!     -------------------------------------
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: hv(N_VERT)
!     -------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)

      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     -------------------------------------
      integer,dimension(:)  :: No_wb(N_WB)
      
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8 :: Sdudx,Sdudy,Sdudz
      real*8 :: Sdvdx,Sdvdy,Sdvdz
      real*8 :: Sdwdx,Sdwdy,Sdwdz
      real*8 :: sxx, sxy, sxz
      real*8 :: syy, syz, szz
!     -------------------------------------
      real*8 :: nut,DD,Cs,Sij,del,dist
      integer:: UseInit
!     -------------------------------------
      real*8, dimension(:,:) :: dudx(N_CELL,NZ)
      real*8, dimension(:,:) :: dudy(N_CELL,NZ)
      real*8, dimension(:,:) :: duds(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dvdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dvds(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdx(N_CELL,NZ)
      real*8, dimension(:,:) :: dwdy(N_CELL,NZ)
      real*8, dimension(:,:) :: dwds(N_CELL,NZ)

!     ____________________________________
!    |                                    |
!    |           LES PARAMETERS           |
!    |____________________________________|

!     Initialization Flag: 1= do setup, 0 = done
      integer :: sg_init = 1
!     Subgrid model identifier
      integer :: sg_model = 1
!     Smagorinsky constant
      real*8 :: sg_smag_const = 0.1d0
!     VanDriest Damping length
      real*8 :: sg_smag_aplus = 26.0d0
!     Mason wall matching power
      real*8 :: sg_smag_mason = 2.0d0
!     Thomas shear constant
      real*8 :: sg_smag_so = 0.0d0
!     Von Karman constant
      real*8 :: von_Karman = 0.415

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,*) '====> Begin: Turbulence model'
         if (sg_model .eq. 1) then
         write(*,*) 'SG_MODEL = Standard Smagrinsky'
         else
         write(*,*) 'Attention: SG_MODEL = unknown!'
         endif
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                 Display LES parameters                 |
!     |________________________________________________________|

#     ifdef KeyDisplay
         write(*,*) 'sg_smag_const = ', sg_smag_const
         write(*,*) 'sg_smag_aplus = ', sg_smag_aplus
         write(*,*) 'sg_smag_mason = ', sg_smag_mason
         write(*,*) 'sg_smag_so = ', sg_smag_so
         write(*,*) 'von_Karman = ', von_Karman
#     endif


      call grandientLSM(dudx,dudy,duds,uf,No_cp,nbe,sig) 
      call grandientLSM(dvdx,dvdy,dvds,vf,No_cp,nbe,sig) 
      call grandientLSM(dwdx,dwdy,dwds,wf,No_cp,nbe,sig)

!      ________________________________________________________
!     |                                                        |
!     |                     Wall distance                      |
!     |________________________________________________________|

!     --------------------------------------------------
!     Use general calculation
#     ifdef KeyLESWallGeneral
         UseInit = sg_init
         if (UseInit .eq. 1) then
            do k=1,NZ
               do i=1,N_CELL
                  walld(i,k) = 1.0E5 ! Set Wall Distance to a large value
               enddo
            enddo
            call sgmkwd(UseInit,                        &
                        xc,yc,sig,dsig,No_cp,nbe,       &
                        xv,yv,sigv,dsigv,No_vp,nbev,    &
                        Hpr,h,eta,etan,                 &
                        Hprv,hv,etav,No_wb)
         endif
         sg_init = UseInit
#     endif
!     --------------------------------------------------
!     Use coordinate sigma
#     ifdef KeyLESWallSigma
        do k =1,NZ
           do i=1,N_CELL
              walld(i,k) = dabs(sig(k))
           enddo
        enddo
#     endif

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call communication3Dtype2(walld)
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!                  Calculate the coefficient                          !
!                                                                     !
!*********************************************************************!

        Cs = sg_smag_const
        DO k=1,NZ
         do i=1,N_CELL
            Sdudx = dudx(i,k)+sigmax(i,k)*duds(i,k)
            Sdvdx = dvdx(i,k)+sigmax(i,k)*dvds(i,k)
            Sdwdx = dwdx(i,k)+sigmax(i,k)*dwds(i,k)
            Sdudy = dudy(i,k)+sigmay(i,k)*duds(i,k)
            Sdvdy = dvdy(i,k)+sigmay(i,k)*dvds(i,k)
            Sdwdy = dwdy(i,k)+sigmay(i,k)*dwds(i,k)
            Sdudz = sigmaz(i,k)*duds(i,k)
            Sdvdz = sigmaz(i,k)*dvds(i,k)
            Sdwdz = sigmaz(i,k)*dwds(i,k)
!           ---------------
!           Deardorff Mesh Scale
            DD  = sg_smag_const*(areaCell(i)*dsig(k)*Hpr(i))**(1.0d0/3.0d0)
            dist = walld(i,k)
!           ---------------
!           Near Wall Damping function
            call sgdamp(DD,dist,sg_smag_mason,sg_smag_aplus,sg_smag_so,von_Karman)
!           ---------------
!           Treat Eddy viscosity as scalar saved in center of the cell
            nut = Hpr(i)*DD
!           --------------------------------------------------
!           Same scalar for all coefficients
#           ifdef KeyLESOptEqual
               sxx = (Sdudx+Sdudx)/2.0d0
               sxy = (Sdudy+Sdvdx)/2.0d0
               sxz = (Sdudz+Sdwdx)/2.0d0
               syy = (Sdvdy+Sdvdy)/2.0d0
               syz = (Sdvdz+Sdwdy)/2.0d0
               szz = (Sdwdz+Sdwdz)/2.0d0
               Sij = dsqrt(2.0d0*(sxx*sxx+syy*syy+szz*szz) &
                         + 4.0d0*(sxy*sxy+sxz*sxz+syz*syz))
               !eddy(i,k)  = nut*Sij
               Gamux(i,k) = Gamux(i,k) + nut*Sij
               Gamuy(i,k) = Gamuy(i,k) + nut*Sij
               Gamuz(i,k) = Gamuz(i,k) + nut*Sij
!              ----------
               Gamvx(i,k) = Gamvx(i,k) + nut*Sij
               Gamvy(i,k) = Gamvy(i,k) + nut*Sij
               Gamvz(i,k) = Gamvz(i,k) + nut*Sij
!              ----------
               Gamwx(i,k) = Gamwx(i,k) + nut*Sij
               Gamwy(i,k) = Gamwy(i,k) + nut*Sij
               Gamwz(i,k) = Gamwz(i,k) + nut*Sij
#           endif
!           --------------------------------------------------
!           Different scalar for each coefficient
#           ifdef KeyLESOptDifferent
               Sij = (Sdudx+Sdudx)/2.0d0
               Gamux(i,k) = Gamux(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdudy+Sdvdx)/2.0d0
               Gamuy(i,k) = Gamuy(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdudz+Sdwdx)/2.0d0
               Gamuz(i,k) = Gamuz(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
!              ----------
               Sij = (Sdvdx+Sdudy)/2.0d0
               Gamvx(i,k) = Gamvx(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdvdy+Sdvdy)/2.0d0
               Gamvy(i,k) = Gamvy(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdvdz+Sdwdy)/2.0d0
               Gamvz(i,k) = Gamvz(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
!              ----------
               Sij = (Sdwdx+Sdudz)/2.0d0
               Gamwx(i,k) = Gamwx(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdwdy+Sdvdz)/2.0d0
               Gamwy(i,k) = Gamwy(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
               Sij = (Sdwdz+Sdwdz)/2.0d0
               Gamwz(i,k) = Gamwz(i,k) + nut*dsqrt(2.0d0*Sij*Sij)
#           endif
         enddo
      ENDDO

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      call communication3Dtype2(Gamux)
      call communication3Dtype2(Gamuy)
      call communication3Dtype2(Gamuz)
      call communication3Dtype2(Gamvx)
      call communication3Dtype2(Gamvy)
      call communication3Dtype2(Gamvz)
      call communication3Dtype2(Gamwx)
      call communication3Dtype2(Gamwy)
      call communication3Dtype2(Gamwz)
#     endif
!     =============== END ================
!     ====================================


!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,*) '<==== Exit: Turbulence Model.'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END SUBROUTINE sgles

!=====================================================================!
!---------------------------------------------------------------------!
!                      END SUBROUTINE sg_setup                        !
!---------------------------------------------------------------------!
!=====================================================================!



!=====================================================================!
!                                                                     !
!                   AUXILIAR SUBROUNTINE 1.1) sgmkwd                  !
!                                                                     !
! This is a geometric value and it doesn’t need to be calculated      !
! at each time step. Miguel Uh                                        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
! This subroutine is for set the wall distance for each cell and will !
! be used in the damping function to calculate the turbulent viscosity!
!                                                                     !
!=====================================================================!

      subroutine sgmkwd (flag,                               &
!                        ------------------------------------
                         xc,yc,sig,dsig,No_cp,nbe,           &
                         xv,yv,sigv,dsigv,No_vp,nbev,        &
!                        ------------------------------------
                         Hpr,h,eta,etan,                     &
                         Hprv,hv,etav,No_wb)

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

      integer :: flag
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
      integer,dimension(:) :: No_wb(N_WB)

!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8, dimension(:) :: dnminxy(N_CELL)
      real*8:: dist,dlx,dly,dls
      real*8:: dnmin,dnmin1,dnmin2
      integer:: ii,jj,kk,s
      real*8:: xm1,ym1,xm2,ym2,xm0,ym0,bottom,upper

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,*) '====> Begin: Set Wall distance!'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                           Initialization                            !
!                                                                     !
!*********************************************************************!

!      _______________________________________________
!     |                                               |
!     |  1.1)  Find wall distance in XY plane         |
!     |_______________________________________________|

      dnminxy = 1.0E5

      DO i=1,N_CELL
         dnmin  = 1.0E5
         dnmin1 = 1.0E5
         dnmin2 = 1.0E5
         ii =0
         jj =0
!        ---------------------------------------------
         do j=1,N_WB
            k = No_wb(j)
            dlx = xc(i)-xv(k)
            dly = yc(i)-yv(k)
            dnmin = sqrt(dlx**2+dly**2)
            if(dnmin .lt. dnmin1) then
               dnmin1 = dnmin
               ii = k
            endif
         enddo
!        ---------------------------------------------
         do j=1,N_WB
            if (j.ne.ii) then
               k = No_wb(j)
               dlx = xc(i)-xv(k)
               dly = yc(i)-yv(k)
               dnmin = sqrt(dlx**2+dly**2)
              if(dnmin .lt. dnmin2) then
                 dnmin2 = dnmin
                 jj = k
              endif 
            endif
         enddo
!        ---------------------------------------------
!        Check the validity of minimum wall distance
         if ((ii .eq. 0) .or. (jj .eq. 0)) then
             print*,'Error!!!: Minimum wall distance not found'
         stop
         endif
!        ---------------------------------------------
!        Calculate the distance from cell center to
!        the line connecting two vertices
         xm1 = xv(ii)
         ym1 = yv(ii)
         xm2 = xv(jj)
         ym2 = yv(jj)
         xm0 = xc(i)
         ym0 = yc(i)
         upper = abs((xm2-xm1)*(ym1-ym0)-(xm1-xm0)*(ym2-ym1))
         bottom = sqrt((xm2-xm1)**2+(ym2-ym1)**2)
         dist = upper/bottom
         dnminxy(i) = dist
      ENDDO

!      _______________________________________________
!     |                                               |
!     |  1.2)  Find the global wall distance          |
!     |_______________________________________________|

      DO k=1,NZ
         do i=1,N_CELL
            dist = dnminxy(i)
            dls  = abs(sig(k)-sigIni)*Hpr(i)
            if (dls.lt.dist) then
              dist = dls
            endif
            walld(i,k) = dist
         enddo
      ENDDO

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!    Mark the initialize of the wall distance is done
     flag = 0

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,*) '<==== Exit: Set Wall distance!'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     return
     end subroutine sgmkwd

!=====================================================================!
!                                                                     !
!                       END SUBROUTINE sgmkwd                         !
!                                                                     !
!=====================================================================!



!=====================================================================!
!                   AUXILIAR SUBROUNTINE 1.2) sgdamp                  !
!                                                                     !
!    This subroutine is to apply the near wall damping function.      !
!                  a) Mason Damping                                   !
!                  b) Van-Driest Damping                              !
!                                                                     !
!=====================================================================!

      subroutine sgdamp(del,dist,                           &
!                       ------------------------------------
                        mason,aplus,so,                     &
!                       ------------------------------------
                        Karman)
                        
      implicit none
      real*8 :: del,dist,mason,aplus,so,Karman
      real*8 :: lwall,n,dn1,dn2

!     ___________________________________________________
!     Mason Damping function
#     ifdef KeyLESMason
         lwall = Karman*dist
         n = mason
         dn1 = (1.0d0/del)**n +(1.0d0/lwall)**n
         dn2 = 1.0d0/dn1
         del = dn2**(1.0d0/n)
#     endif

!     ____________________________________________________
!     Van Driest Damping function
#     ifdef KeyLESVan
!        This part haven't been implemented yet. This will require
!        The friction velocity at the near bed.
!        Perhaps for horizontal wall and vertical wall.
#     endif

!     _____________________________________________________
!     (Delta*const)^2
      del = del**2

      return
      end subroutine sgdamp
     
!=====================================================================!
!                                                                     !
!                       END SUBROUTINE sgdamp                         !
!                                                                     !
!=====================================================================!