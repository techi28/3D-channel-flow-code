!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!      SOLUTION OF THE POISSON EQUATION FOR FREE SURFACE PROBLEMS     !
!                             Nov 2015                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE     bicgstab(phi,phiv,rhs,                  &
                              xc,yc,sig,dsig,No_cp,nbe,      &
                              xv,yv,sigv,dsigv,No_vp,nbev,   &
                              Hpr,h,etan,                    &
                              Hprv,hv,etav)

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

      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ----------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)

!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable :: Am0,Am1,Am2,Am3,AmT,AmB
      real*8 :: Vol
      integer :: flg,IDisplay
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: bicgstab'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IDisplay = 1
#     ifdef KeyParallel
      if(rang_topo .ne. 0) IDisplay = 0
#     endif

      if(IDisplay .eq. 1) print*, '       SOLVE FOR P FIELD'
!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
       flg = 0 ! Flag to normalize the source term and pressure field
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(Am0(N_CELL0,NZ),Am1(N_CELL0,NZ),Am2(N_CELL0,NZ), &
               Am3(N_CELL0,NZ),AmT(N_CELL0,NZ),AmB(N_CELL0,NZ))
!      ________________________________________________________
!     |                                                        |
!     |                     Initial guess                      |
!     |________________________________________________________|
#     ifdef KeyTESTChannel
       flg = 1
#     endif

!     ________________________________________________________
!    |                                                        |
!    |            Matrix Am & Bm of the diffusion term        |
!    |________________________________________________________|
#     ifdef KeyPPECenter
      do k=1,NZ
         do i=1,N_CELL0
            Am0(i,k)  = AmC0(i,k)
            Am1(i,k)  = AmC1(i,k)
            Am2(i,k)  = AmC2(i,k)
            Am3(i,k)  = AmC3(i,k)
            AmT(i,k)  = AmCT(i,k)
            AmB(i,k)  = AmCB(i,k)
         enddo
      enddo
#     endif
!   ---------------------------------------------------------------
      if (flg .eq. 1) call Glob_Norm(rhs,xc,yc,sig,dsig,No_cp,nbe,1)
!*********************************************************************!
!                                                                     !
!                 Solution of the system (BICGSTAB) 3D                !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |           New right-hand side: rhs - Bm(vertex)        |
!     |________________________________________________________|
       call BICGSTAB_SOLVE(Am0,Am1,Am2,Am3,AmT,AmB,         &
                           rhs,phi,                         &
                           xc,yc,sig,dsig,No_cp,nbe,TagBCp)

!     --------------------------------------------------------
       if(flg .eq. 1) call Glob_Norm(phi,xc,yc,sig,dsig,No_cp,nbe,0)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
       call BCglobalP(phi,phiv,                      &
                      xc,yc,sig,dsig,No_cp,nbe,      &
                      xv,yv,sigv,dsigv,No_vp,nbev,3)
!      ________________________________________________________
!     |                                                        |
!     |                         De-allocate                    |
!     |________________________________________________________|

      deallocate(Am0,Am1,Am2,       &
                 Am3,AmT,AmB)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- End   subroutine: bicgstab'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF Poisson for Free surface                  !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

SUBROUTINE get_matrix_coef(Am0,Am1,Am2,Am3,AmT,AmB,         &
                           xc,yc,sig,dsig,No_cp,nbe)

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
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ------------------------------------------
      real*8, dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
!     --------------------------------------
      integer:: jv1,jv2,jv3,elem,ii
      integer :: jc1,jc2,jc3,jc
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin: check PPE matrix'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
!    By modify the coefficients, no need to exchange data on no-periodic boundaries
!    Apply boundary condition, here we set all Newmann bc on all boundaries
!    instead apply bc, here we do check and stop programm if wrong.
!    ---------------------------------------------
!        Top and bottom wall
        IF (ZPB .eq. 0) THEN
            do i=1,N_CELL0
                if(bctop .eq. 1) then
                    if(abs(AmT(i,NZ-1)) .gt. 1.0E-7) then
                        print*, 'Error! Pres Coef for top bc', AmT(i,NZ-1)
                        stop
                    endif
                endif
                if(bcbot .eq. 1) then
                    if(abs(AmB(i,2)) .gt. 1.0E-7) then
                        print*, 'Error! Pres Coef for bot bc', AmB(i,2)
                        stop
                    endif
                endif
            enddo
        ENDIF
!      -------------------------------------------
!       Horizontal boundary
          do i=1,N_CELL0
           if (nbe(i).ne.0) then
                do j=1,3
                nc=No_cp(i,j)
                ii = bcpc(nc)
                IF (ii .eq. 1) THEN
                   if (j .eq. 1) then
                     do k=2,NZ-1
                       if(abs(Am1(i,k)) .gt. 1.0E-7) then
                          print*, 'Error! Pres Coef for xy bc', Am1(i,k)
                          stop
                      endif
                     enddo
                  elseif (j .eq. 2) then
                     do k=2,NZ-1
                       if(abs(Am2(i,k)) .gt. 1.0E-7) then
                          print*, 'Error! Pres Coef for xy bc', Am2(i,k)
                          stop
                      endif
                     enddo
                  elseif (j .eq. 3) then
                     do k=2,NZ-1
                       if(abs(Am3(i,k)) .gt. 1.0E-7) then
                          print*, 'Error! Pres Coef for xy bc', Am3(i,k)
                          stop
                      endif
                     enddo
                 endif
              ENDIF
             enddo
          endif
        enddo

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<----- End: check PPE matrix'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    END OF ALTER MARTIX FOR LINEAR                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BICGSTAB_SOLVE(Am0,Am1,Am2,Am3,AmT,AmB,    &
                                rhs,phi,                    &
                                xc,yc,sig,dsig,No_cp,nbe,TagBC)
!       ____________________________________
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

!     ----------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
!     ----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ------------------------------------------
      real*8, dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8, dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8, dimension(:,:) :: AmB(N_CELL0,NZ)
      integer :: TagBC
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
!    -------------------------------------------
     real*8,dimension(:,:),allocatable :: r0bar,vk,sk,zk
     real*8,dimension(:,:),allocatable :: pk,tk,yk,rk
!     --------------------------------------
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
      integer :: iter,IDISPLAY
      real*8 :: e1norm, b1norm,r1norm
      real*8 :: e2norm, b2norm,r2norm
      real*8 :: rhok,rhokm,r0barvk
      real*8 :: alpha,omegakm
      real*8 :: sum_e1bval,sum_e2bval
      real*8 :: sum_e1resid,sum_e2resid
      real*8 :: resid,maxrhs,enorm
      real*8 :: var1,var2
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
       write(*,'(t5,60a)'), '-----> Begin: BICGSTAB SOLVE'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

     IDISPLAY = 1

#    ifdef KeyParallel
 	  if (rang_topo.ne.0)  IDISPLAY = 0
#    endif

!     --------------------------------------
!      2D Formation
       if(I3D .ne. 1) then
          do k=2,NZ-1
            do i=1,N_CELL0
                Am0(i,k) = Am0(i,k)+AmT(i,k)+AmB(i,k)
                AmT(i,k) = 0.
                AmB(i,k) = 0.
            enddo
          enddo
       endif
!     --------------------------------------
!     Allocate variable for bicgstab solver
       allocate(r0bar(N_CELL,NZ),vk(N_CELL,NZ),sk(N_CELL,NZ),       &
                tk(N_CELL,NZ),yk(N_CELL,NZ),rk(N_CELL,NZ),          &
                pk(N_CELL,NZ),zk(N_CELL,NZ))
!     --------------------------------------
!     Initialize vk, pk
      do k=1, NZ
        do i=1,N_CELL
           r0bar(i,k) = 0.0d0
           rk(i,k) = 0.0d0
           vk(i,k) = 0.0d0
           sk(i,k) = 0.0d0
           tk(i,k) = 0.0d0
           yk(i,k) = 0.0d0
           pk(i,k) = 0.0d0
           zk(i,k) = 0.0d0
        enddo
      enddo
!     -------------------------------------
!      Get bnorm and rnorm
      resid = 0.0d0
      sum_e1bval = 0.0d0
      sum_e1resid = 0.0d0
      sum_e2bval = 0.0d0
      sum_e2resid = 0.0d0
!     -------------------------------------
!     Solving for only Center value
#     ifdef KeyPPECenter
      do k=2, NZ-1
        do i=1, N_CELL0
            jc1 = No_cp(i,1)
            jc2 = No_cp(i,2)
            jc3 = No_cp(i,3)
!           _________________________________________________
!           Bm(vertex)
            resid = rhs(i,k)- ( Am0(i,k)*phi(i,k)     &
                               +Am1(i,k)*phi(jc1,k)   &
                               +Am2(i,k)*phi(jc2,k)   &
                               +Am3(i,k)*phi(jc3,k)   &
                               +AmT(i,k)*phi(i,k+1)   &
                               +AmB(i,k)*phi(i,k-1))

           sum_e2bval = sum_e2bval + rhs(i,k)*rhs(i,k)
           sum_e2resid = sum_e2resid + resid*resid
           sum_e1bval = sum_e1bval + dabs(rhs(i,k))
           sum_e1resid = sum_e1resid + dabs(resid)
           rk(i,k) = resid
           r0bar(i,k) = resid
        enddo
      enddo
#    endif
 !    ----------------------------------------------
 !    For parallazation
#    ifdef KeyParallel
      call MPI_Barrier(comm3D,code)
      call SUM_parallel(sum_e2bval,b2norm)
      call SUM_parallel(sum_e2resid,r2norm)
      call SUM_parallel(sum_e1bval,b1norm)
      call SUM_parallel(sum_e1resid,r1norm)
#    else
      b2norm = sum_e2bval
      r2norm = sum_e2resid
      b1norm = sum_e1bval
      r1norm = sum_e1resid
#    endif
      b2norm = dsqrt(b2norm)
      r2norm = dsqrt(r2norm)
      if (b2norm == 0) b2norm = 1.0d0
      if (b1norm == 0) b1norm = 1.0d0
      e2norm = r2norm/b2norm
      e1norm = r1norm/b1norm
      if((e1norm .lt. eps) .or. (e2norm .lt. eps)) then
       if (IDISPLAY .eq. 1) then
       print*, 'Pressure converged without solve !'
       endif
       goto 999
      endif
!     -----------------------------------------------
!     Set the parameters
      rhokm = 1.0d0
      alpha = 1.0d0
      omegakm = 1.0d0
      iter = 0
!    -----------------------------------------------
666  continue
!    Check convergence
     if (dabs(rhokm) .lt. 1.0e-30) then
      if (IDISPLAY .eq. 1) then
      print*, 'Alogorithm already converged !'
      endif
      goto 999
     endif
     if (dabs(omegakm) .lt. 1.0e-30) then
      if (IDISPLAY .eq. 1) then
      print*, 'Alogorithm already converged !'
      endif
      goto 999
     endif
!    ------------------------------------------------
!     Update rhok
      call bicgstab_innerprod(r0bar,rk,rhok)
!    ------------------------------------------------
!     Update pk
      call update_pk(iter,rhok,rhokm,alpha,omegakm, &
                     pk,rk,vk)
!     -----------------------------------------------
!     DIAG preconditioner
      call diag_precdt(yk,pk,     &
                       Am0,Am1,Am2,Am3,AmT,AmB, &
                       nbe,No_cp)
!     -----------------------------------------------
!     get Apk
      call get_Apk(vk,yk,Am0,Am1,Am2,Am3,AmT,AmB, &
                   xc,yc,sig,dsig,No_cp,nbe,TagBC)
!     ----------------------------------------------
!     Update r0barvk
      call bicgstab_innerprod(r0bar,vk,r0barvk)
!     ----------------------------------------------
!     Check convergence
        if (dabs(r0barvk) .lt. 1.0e-30) then
          if (IDISPLAY .eq. 1) then
          print*, 'Alogorithm already converged !'
          endif
        goto 999
        endif
!    -----------------------------------------------
!     Update alpha
      alpha = rhok/r0barvk
!    -----------------------------------------------
!     Update sk
      call update_sk(alpha,rk,vk,sk)
!    -----------------------------------------------
!     DIAG preconditioner
      call diag_precdt(zk,sk,&
                       Am0,Am1,Am2,Am3,AmT,AmB,&
                       nbe,No_cp)
!    -----------------------------------------------
!     get Apk
      call get_Apk(tk,zk,Am0,Am1,Am2,Am3,AmT,AmB, &
                       xc,yc,sig,dsig,No_cp,nbe,TagBC)
!    -----------------------------------------------
!     Update omegakm
      call bicgstab_innerprod(tk,sk,var1)
      call bicgstab_innerprod(tk,tk,var2)
      omegakm = var1/var2
!    ----------------------------------------------
!      Update rnorm
      call update_xkrk(alpha,omegakm,r1norm,r2norm,maxrhs, &
                       phi,yk,zk,sk,tk,rk)
!    ----------------------------------------------
!    Update enorm
     e2norm = r2norm/b2norm
     e1norm = r1norm/b1norm
     rhokm = rhok
!    ----------------------------------------------
     iter = iter + 1
!    ---------------------------------------------
!    Check convergence
     enorm = 1.0E10
!     if(TagBC .eq. TagBCp) maxrhs = dt*maxrhs
!     if(enorm .gt. maxrhs) enorm = maxrhs
     if(enorm .gt. e1norm) enorm = e1norm
     if(enorm .gt. e2norm) enorm = e2norm
!    ---------------------------------------------
        if (enorm .lt. eps) then
!      	    if (IDISPLAY .eq. 1) then
!      	      print*, 'BiCGstab converged.'
!      	    endif
            goto 888
        elseif (iter .gt. MaxIters) then
            if (IDISPLAY .eq. 1) then
                print*, 'Maximum allowed iteration reached.'
            endif
            goto 888
        else
!           Show convergence critia
            if (0.eq.mod(iter,10)) then
                if (IDISPLAY .eq. 1) then
                    print*, 'iter=',iter,'e1norm=',e1norm,'e2norm=',e2norm!,'maxrhs =', maxrhs
                endif
            endif
            goto 666
        endif

888  continue
      	  if (IDISPLAY .eq. 1) then
           print*, 'iter=',iter,'e1norm=',e1norm,'e2norm=',e2norm
          endif
999  continue
!    ---------------------------------------
!    Deallocate variable
      deallocate(r0bar,vk,sk,tk,yk,pk,rk,zk)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
       write(*,'(t5,60a)'), '<----- End: BICGSTAB SOLVE'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    	END OF BICGSOLVE                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Glob_Norm(field,xc,yc,sig,dsig,No_cp,nbe,tags)
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
      real*8,dimension(:,:) :: field(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8, dimension(:)  :: sig(NZ)
      real*8, dimension(:)  :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer :: tags
      real*8 :: local_sum, glob_sum, avec
      real*8 :: fB,Vol
      integer :: ncount,elem,ii,flag
!     ------------------------------------
      local_sum = 0.0
      glob_sum = 0.0

      do k=2,NZ-1
        do i=1, N_CELL0
          Vol = areaCell(i)*dsig(k)
          if(tags .eq. 1) then
            local_sum = field(i,k) + local_sum
          else
            local_sum = field(i,k)*Vol + local_sum
          endif
         enddo
      enddo

#     ifdef KeyParallel
          call MPI_Barrier(comm3D,code)
          call SUM_parallel(local_sum,glob_sum)
#     else
          glob_sum = local_sum
#     endif

      avec = glob_sum/(TotalVolume*(SigFin-SigIni))

        do k=2,NZ-1
            do i=1, N_CELL0
                if(tags .eq. 1) then
                        field(i,k) = field(i,k) - avec*areaCell(i)*dsig(k)
                else
                        field(i,k) = field(i,k) - avec
                endif
            enddo
        enddo

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    	   END OF Glob_norm                           !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
 SUBROUTINE bicgstab_innerprod(vec1,vec2,prod)
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
      real*8,dimension(:,:) :: vec1(N_CELL,NZ)
      real*8,dimension(:,:) :: vec2(N_CELL,NZ)
      real*8 :: prod
!     ---------------------------------------
      real*8 :: local_sum, glob_sum
      local_sum = 0.0d0
      glob_sum = 0.0d0
!     ---------------------------------------------
      do k=2,NZ-1
        do i=1,N_CELL0
         local_sum = local_sum + vec1(i,k)*vec2(i,k)
        enddo
      enddo
!     --------------------------------------------
!      SUM ACROSS PROCS IF NECCESSARY
#     ifdef KeyParallel
      call MPI_Barrier(comm3D,code)
      call SUM_parallel(local_sum,glob_sum)
#     else
      glob_sum = local_sum
#     endif
!     --------------------------------------------
!     Assign final value
      prod = glob_sum

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  	 END OF bicgstab_innerprod                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
       SUBROUTINE update_pk(iter,rhok,rhokm,alpha,omegakm, &
                                     pk,rk,vk)
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
      real*8,dimension(:,:) :: pk(N_CELL,NZ)
      real*8,dimension(:,:) :: rk(N_CELL,NZ)
      real*8,dimension(:,:) :: vk(N_CELL,NZ)
      integer :: iter
      real*8 :: rhok,rhokm,alpha,omegakm
!     -------------------------------------
      real*8 :: beta
!     --------------------------------------
!     get pk at iter = 0
      if (iter .eq. 0) then
        do k=2,NZ-1
         do i=1, N_CELL0
         pk(i,k) = rk(i,k)
         enddo
        enddo
       else
         beta = (rhok/rhokm)*(alpha/omegakm)
         do k=2,NZ-1
           do i=1,N_CELL0
             pk(i,k) = rk(i,k)+beta*(pk(i,k)-omegakm*vk(i,k))
           enddo
          enddo
        endif
      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!              	 	     END OF update_pk                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      SUBROUTINE diag_precdt(vec1,vec2, &
                             Am0,Am1,Am2,Am3,AmT,AmB, &
                             nbe,No_cp)
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
      real*8,dimension(:,:) :: vec1(N_CELL,NZ)
      real*8,dimension(:,:) :: vec2(N_CELL,NZ)
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ---------------------------------------
!     local variables
      integer :: elem,ii,flag
      real*8 :: fB,resid,omegaf
      integer :: Opt,iter,Kiter
      integer :: jc1,jc2,jc3
      real*8 :: mxm1,mxm2,mxm3,mxmT,mxmB,mdiag
!     -------------------------------------
!     Different preconditioner are used
!     1 : DIAG preconditoner
!     2 : SOR preconditioner
!     3 : Not working seems only working for Cartesian
!     --------------------------------------
      Opt = 1
!     ======================================
!     DIAG preconditioner
!     ======================================
#     ifdef KeyPPECenter
      if(Opt .eq. 1) then
      do k=2,NZ-1
        do i=1, N_CELL0
         vec1(i,k) = vec2(i,k)/Am0(i,k)
        enddo
      enddo
!     =======================================
!     Additive Schwarz pre-conditioner based on use
!     of K (<=10) iterations of SOR
!     with a fixed over-relaxation factor omega.
!     =======================================
      elseif(Opt .eq. 2) then
        omegaf = 1.75
        Kiter = 5
        do k=1,NZ
            do i=1,N_CELL
                vec1(i,k) = 0.0d0
            enddo
        enddo

        do iter=1,Kiter
            do k=2,NZ-1
                do i=1,N_CELL0
                    jc1 = No_cp(i,1)
                    jc2 = No_cp(i,2)
                    jc3 = No_cp(i,3)
                    resid =   Am0(i,k)*vec1(i,k)     &
                            + Am1(i,k)*vec1(jc1,k)   &
                            + Am2(i,k)*vec1(jc2,k)   &
                            + Am3(i,k)*vec1(jc3,k)   &
                            + AmT(i,k)*vec1(i,k+1)   &
                            + AmB(i,k)*vec1(i,k-1)
                    resid = -resid + vec2(i,k)
                    vec1(i,k) = vec1(i,k) + omegaf*resid/Am0(i,k)
                enddo
            enddo
        enddo
!     =======================================
!     Additive Schwarz pre-conditioner based on
!    "Incomplete Poisson" preconditioner
!     =======================================
      elseif(Opt .eq. 3) then
       do k=1,NZ
         do i=1,N_CELL
            vec1(i,k) = 0.0d0
         enddo
       enddo

      do k=2,NZ-1
          do i=1,N_CELL0
              jc1 = No_cp(i,1)
              jc2 = No_cp(i,2)
              jc3 = No_cp(i,3)
              if(jc1 .le. N_CELL0) then
                 mxm1 = 1.0d0/Am0(jc1,k)
              else
                 mxm1 = 0.0d0
              endif
              if(jc2 .le. N_CELL0) then
                 mxm2 = 1.0d0/Am0(jc2,k)
              else
                 mxm2 = 0.0d0
              endif
              if(jc3 .le. N_CELL0) then
                 mxm3 = 1.0d0/Am0(jc3,k)
              else
                 mxm3 = 0.0d0
              endif
              if(k .eq. 2) then
                mxmB = 0.0d0
              else
                mxmB = 1.0d0/Am0(i,k-1)
              endif
              if(k .eq. NZ-1) then
                mxmT = 0.0d0
              else
                mxmT = 1.0d0/Am0(i,k+1)
              endif
              mdiag = 1.0d0 + 0.5d0*(mxm1*mxm1  &
                         + mxm2*mxm2 +mxm3*mxm3 &
                         + mxmT*mxmT +mxmB*mxmB)
              vec1(i,k) =  mdiag*vec2(i,k)     &
                           +mxm1*vec2(jc1,k)   &
                           +mxm2*vec2(jc2,k)   &
                           +mxm3*vec2(jc3,k)   &
                           +mxmT*vec2(i,k+1)   &
                           +mxmB*vec2(i,k-1)
          enddo
      enddo

      endif
#    endif

     RETURN
     END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  	   END OF bicg_precdt                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
      SUBROUTINE get_Apk(vec1,vec2,Am0,Am1,Am2,Am3,AmT,AmB, &
                       xc,yc,sig,dsig,No_cp,nbe,TagBC)
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
      real*8,dimension(:,:) :: vec1(N_CELL,NZ)
      real*8,dimension(:,:) :: vec2(N_CELL,NZ)
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
!     -----------------------------------------
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
      integer :: TagBC
!    -----------------------------------------
      real*8 :: rsum
      integer :: jc1,jc2,jc3
!    -----------------------------------------------
!    set boundary condition for vec1 as it will be
!    used to multiply with the pressure coefficient
     if(TagBC .eq. TagBCp) then
      call BCpcenter(vec2,xc,yc,sig,dsig,No_cp,nbe)
     else
      call BCvelcenter3D(vec2,xc,yc,sig,dsig,No_cp,nbe,TagBC)
     endif

#    ifdef KeyPPECenter
        do k=2,NZ-1
            do i=1, N_CELL0
                jc1 = No_cp(i,1)
                jc2 = No_cp(i,2)
                jc3 = No_cp(i,3)
                rsum =   Am0(i,k)*vec2(i,k)     &
                        +Am1(i,k)*vec2(jc1,k)   &
                        +Am2(i,k)*vec2(jc2,k)   &
                        +Am3(i,k)*vec2(jc3,k)   &
                        +AmT(i,k)*vec2(i,k+1)   &
                        +AmB(i,k)*vec2(i,k-1)

                vec1(i,k) = rsum
            enddo
        enddo
#    endif

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  	 	END OF get_Apk                        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
 	SUBROUTINE update_sk(alpha,rk,vk,sk)
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
      real*8,dimension(:,:) :: sk(N_CELL,NZ)
      real*8,dimension(:,:) :: rk(N_CELL,NZ)
      real*8,dimension(:,:) :: vk(N_CELL,NZ)
      real*8 :: alpha

         do k=2,NZ-1
           do i=1,N_CELL0
             sk(i,k) = rk(i,k)-alpha*vk(i,k)
           enddo
          enddo

      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  	   END OF update_sk                           !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
	SUBROUTINE update_xkrk(alpha,omegakm,r1norm,r2norm,maxrhs,phi,yk,zk,sk,tk,rk)
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
      real*8,dimension(:,:) :: yk(N_CELL,NZ)
      real*8,dimension(:,:) :: zk(N_CELL,NZ)
      real*8,dimension(:,:) :: sk(N_CELL,NZ)
      real*8,dimension(:,:) :: tk(N_CELL,NZ)
      real*8,dimension(:,:) :: rk(N_CELL,NZ)
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8 :: alpha,omegakm
      real*8 :: r1norm,r2norm,maxrhs,rhsmax
      real*8 :: sum_r1norm,sum_r2norm
!     ---------------------------------------
      real*8 :: resid,xnew
      r1norm = 0.0d0
      r2norm = 0.0d0
      rhsmax = 0.0d0
      maxrhs = 0.0d0
!     ---------------------------------------
         do k=2,NZ-1
           do i=1,N_CELL0
             xnew = phi(i,k)+alpha*yk(i,k)+omegakm*zk(i,k)
             resid = sk(i,k) - omegakm*tk(i,k)
             r2norm = r2norm + resid*resid
             r1norm = r1norm + dabs(resid)
             phi(i,k) = xnew
             rk(i,k) = resid
             if(dabs(resid) .gt. rhsmax) rhsmax = dabs(resid)
           enddo
          enddo
!    ----------------------------------------
#     ifdef KeyParallel
      call MPI_Barrier(comm3D,code)
      call SUM_parallel(r2norm,sum_r2norm)
      call SUM_parallel(r1norm,sum_r1norm)
      call MAX_parallel(rhsmax,maxrhs)
#     else
      sum_r2norm = r2norm
      sum_r1norm = r1norm
      maxrhs = rhsmax
#     endif

      r2norm = dsqrt(sum_r2norm)
      r1norm = sum_r1norm
      RETURN
      END
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  	 END OF update_xkrk                           !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
