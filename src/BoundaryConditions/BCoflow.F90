!  ==================================================================!
!  ==================================================================!
!  ==================================================================!
!   Calcualte the outflow variables at the beginning of the simulation
!   later, it will be updated during the loop
   SUBROUTINE  calculate_oflow(uo,vo,wo, &
                               phiu,phiv,phiw,       &    
                               xc,yc,sig,dsig,No_cp,nbe, &
                               xv,yv,sigv,dsigv,No_vp,nbev,flag)
!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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
      real*8,dimension(:,:) :: uo(N_CELL,NZ)
      real*8,dimension(:,:) :: vo(N_CELL,NZ)
      real*8,dimension(:,:) :: wo(N_CELL,NZ)
      real*8,dimension(:,:) :: phiu(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_CELL,NZ)
      real*8,dimension(:,:) :: phiw(N_CELL,NZ)
!     --------------------------------------
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
!     ----------------------------------------
!     local variable
      integer :: flag,elem
      real*8 :: inflow,x,y,z
!     ---------------------------------------
!     Initialization
!    ---------------------------------------
!    allocate indicator for outflow cell
        allocate(tago(N_CELL))
        
        do i=1, N_CELL
          tago(i) = -1.0d0
        enddo
        
        do i=1, N_CELL0
           do j=1,3
            nc = No_cp(i,j)
            if(bcuc(nc) .eq. 4) then
            tago(i) = 1.0d0
            endif
           enddo
       enddo
!     ---------------------------------------
!     assign value
     if(flag .eq. 0) then
      do k=1,NZ
         do i=1,N_CELL
           uo(i,k) = 0.0d0
           vo(i,k) = 0.0d0
           wo(i,k) = 0.0d0
         enddo
      enddo

      DO k=1,NZ
         DO i=1,N_CELL0
	    if (tago(i) .eq. 1) then
               x = xc(i)
               y = yc(i)
               z = sig(k)
               uo(i,k) = inflow(x,y,z,1)
               vo(i,k) = inflow(x,y,z,2)
               wo(i,k) = inflow(x,y,z,3)
            endif
          ENDDO
        ENDDO
     endif

#   ifdef KeyParallel
      call MPI_Barrier(comm3D,code)
      call communication2D(tago)
      call communication3D(uo)
      call communication3D(vo)
      call communication3D(wo)
#   endif

   RETURN
   END
!  ==================================================================!
!  ==================================================================!
!  ==================================================================!
!  update the oflow bc with the latest flow info
!  for mean inflow equals 1.0d0
   SUBROUTINE  update_oflow(uo,vo,wo, &
                            phiu,phiv,phiw,       &    
                            xc,yc,sig,dsig,No_cp,nbe, &
                            xv,yv,sigv,dsigv,No_vp,nbev)
!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |   Keys and common parameters       |
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
      real*8,dimension(:,:) :: uo(N_CELL,NZ)
      real*8,dimension(:,:) :: vo(N_CELL,NZ)
      real*8,dimension(:,:) :: wo(N_CELL,NZ)
      real*8,dimension(:,:) :: phiu(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_CELL,NZ)
      real*8,dimension(:,:) :: phiw(N_CELL,NZ)
!     --------------------------------------
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
!     ----------------------------------------
!     local variable
      real*8  :: fBu,fBv,fBw
!     ---------------------------------------
!     assign value
      DO k=1,NZ
         DO i=1,N_CELL0
            if(tago(i) .gt. 0.0d0) then
              do j=1,3
                nc = No_cp(i,j)
                if(bcuc(nc) .eq. 4) then
                 fBu = dt*(phiu(nc,k)-phiu(i,k))/dlCC(i,j)
                 fBv = dt*(phiv(nc,k)-phiv(i,k))/dlCC(i,j)
                 fBw = dt*(phiw(nc,k)-phiw(i,k))/dlCC(i,j)
                 uo(i,k) = uo(i,k) - fBu
                 vo(i,k) = vo(i,k) - fBv
                 wo(i,k) = wo(i,k) - fBw
                 if(j .eq. 1) then
                    U1FACE(i,k) = uo(i,k)
                 elseif(j .eq. 2) then
                    U2FACE(i,k) = uo(i,k)
                 else
                    U3FACE(i,k) = uo(i,k)
                endif
                goto 100
               endif
             enddo
            endif
100         continue
          ENDDO
        ENDDO

#   ifdef KeyParallel
      call MPI_Barrier(comm3D,code)
      call communication3D(uo)
      call communication3D(vo)
      call communication3D(wo)
#   endif

   RETURN
   END
!  ==================================================================!
!  ==================================================================!
!  ==================================================================!
