!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   SOLUTION OF THE POISSON EQUATION                  !
!                             March 2017                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Poisson(phi,phiv,                      &
                         rhs,Gamx,Gamy,Gamz,            &
!                        -------------------------------
                         Hpr,etan,                      &
                         Hprv,etav,                     &
!                        -------------------------------
                         h,hv,                          &
!                        -------------------------------
                         xc,yc,sig,dsig,No_cp,nbe,      &
                         xv,yv,sigv,dsigv,No_vp,nbev)
 
!---------------------------------------------------------------------!
!                                                                     !
!     This subroutine solves the poisson equation of the variable     !
!     phi, given the diffusive values, right-hand side and the        !
!     correct boundary conditions.                                    !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |   Name     |   Size     | Description                         |  !
!  |____________|____________|_____________________________________|  !
!  | <--phi     |(N_CELL,NZ) | Cell-center solution at t(n+1)      |  !
!  | <--phiv    |(N_VERT,NZ) | Cell-vertex solution at t(n+1)      |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !
!  |____________|____________|_____________________________________|  !
!  | --> Gam    |(N_CELL,NZ) | Diffusive coefficient of the press. |  !
!  | --> rhs    |(N_CELL,NZ) | right-hand side                     |  !
!  |____________|____________|_____________________________________|  !
!  | --> phiOld |(N_CELL,NZ) | Cell-center solution at t(n)        |  !
!  | --> xc,yc  |(N_CELL)    | Coordinates of the cell centers     |  !
!  | --> sig    |(NZ)        | Sigma value at the cell centers     |  !
!  | --> dsig   |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_cp  |(N_CELL,3)  | Numbering of surrounding three cells|  !
!  | --> nbe    |(N_CELL)    | Tag: Type of cell (inside or bc)    |  !
!  |____________|____________|_____________________________________|  !
!  | --> phiOldv|(N_VERT,NZ) | Cell-vertex solution at t(n)        |  !
!  | --> xv,yv  |(N_VERT)    | Coordinates of the cell vertices    |  !
!  | --> sigv   |(NZ-1)      | sigma of the vertex points          |  !
!  | --> dsigv  |(NZ-1)      | Increment = sig(k+1)-sig(k)         |  !
!  | --> No_vp  |(N_CELL0,3) | Numbering of the cell vertices      |  !
!  | --> nbev   |(N_VERT)    | Tag: Type of vertex (inside or bc)  |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !
!  |____________|____________|_____________________________________|  !
!  |    phiT    |(N_CELL)| Dirichlet boundary condition top        |  !
!  |    phiB    |(N_CELL)| Dirichlet boundary condition bottom     |  !
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |   - sor2D                      ( sor2D.F90     )              |  !
!  |   - sor3D                      ( sor3D.F90     )              |  !
!  |   - gmres2D                    ( gmres2D.F90   )              |  !
!  |   - gmres3D                    ( gmres3D.F90   )              |  !
!  |   - sorNew2D                   ( sorNew2D.F90  )              |  !
!  |   - sorNew3D                   ( sorNew3D.F90  )              |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!        Common variables used                                        !
!    *   Common variables modified                                    !
!---------------------------------------------------------------------!

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
!     ----------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
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
!     ----------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
!     ----------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)       
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

!     --------------------------------------      
      real*8,dimension(:,:) :: phiOld(N_CELL,NZ)
      real*8,dimension(:,:) :: rhsOld(N_CELL,NZ)
      real*8,dimension(:,:) :: phiOldv(N_VERT,NZ-1) 
      real*8 :: Vol
            
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- Begin subroutine: Poisson'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                 Solution of the Poisson Equation 3D                 !
!                                                                     !
!*********************************************************************!
      
         phiOld  = phi
         phiOldv = phiv
!        ________________________________________________________
!       |                                                        |
!       |                     Approximation                      |
!       |________________________________________________________|

!        ________________________________________________________
!        Jacobi
#        ifdef KeyJacobi 
            phi  = phiOld
            phiv = phiOldv
            call  Jacobi(phi,phiv,                    &
                         rhs,Gamx,Gamy,Gamz,          &
                         xc,yc,sig,dsig,No_cp,nbe,    &
                         xv,yv,sigv,dsigv,No_vp,nbev, &
                         Hpr,h,etan,                  &
                         Hprv,hv,etav)
!           ====================================
!           ====================================
#           ifdef KeyParallel
               call MPI_Barrier(comm3D,code)
#           endif
!           ====================================
!           ====================================
#        endif                  
!        ________________________________________________________
!        S.O.R. Method 
#        ifdef KeySOR 
            phi  = phiOld
            phiv = phiOldv
            call sor3D(phi,phiv,                    &
                       rhs,Gamx,Gamy,Gamz,          &
                       xc,yc,sig,dsig,No_cp,nbe,    &
                       xv,yv,sigv,dsigv,No_vp,nbev, &
                       Hpr,h,etan,                  &
                       Hprv,hv,etav)
#        endif
!        ________________________________________________________
!        Simultaneous Over-Relaxation Method  (SiOR)
#        ifdef KeySiOR 
            phi  = phiOld
            phiv = phiOldv
            call sior3D(phi,phiv,                    &
                        rhs,Gamx,Gamy,Gamz,          &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        Hpr,h,etan,                  &
                        Hprv,hv,etav)
!           ====================================
!           ====================================
#           ifdef KeyParallel
               call MPI_Barrier(comm3D,code)
#           endif
!           ====================================
!           ====================================
#        endif
!        ________________________________________________________
!        Multicoloring Over-Relaxation Method  (MSOR)
#        ifdef KeyMultiSOR 
            phi  = phiOld
            phiv = phiOldv
            call Multisor3D(phi,phiv,                    &
                            rhs,Gamx,Gamy,Gamz,          &
                            xc,yc,sig,dsig,No_cp,nbe,    &
                            xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,etan,                  &
                            Hprv,hv,etav)
!           ====================================
!           ====================================
#           ifdef KeyParallel
               call MPI_Barrier(comm3D,code)
#           endif
!           ====================================
!           ====================================
#        endif
!        ________________________________________________________
!        Simultaneous Over-Relaxation Method  (JSOR)
#        ifdef KeyJSOR  
            phi  = phiOld
            phiv = phiOldv
            call Jsor3D(phi,phiv,                    &
                        rhs,Gamx,Gamy,Gamz,          &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        Hpr,h,etan,                  &
                        Hprv,hv,etav)
!           ====================================
!           ====================================
#           ifdef KeyParallel
               call MPI_Barrier(comm3D,code)
#           endif
!           ====================================
!           ====================================
#        endif
!        ________________________________________________________
!        Partition Over-Relaxation Method  (PSOR)
#        ifdef KeyPSOR 
            phi  = phiOld
            phiv = phiOldv
            call psor3D(phi,phiv,                    &
                        rhs,Gamx,Gamy,Gamz,          &
                        xc,yc,sig,dsig,No_cp,nbe,    &
                        xv,yv,sigv,dsigv,No_vp,nbev, &
                        Hpr,h,etan,                  &
                        Hprv,hv,etav)
!           ====================================
!           ====================================
#           ifdef KeyParallel
               call MPI_Barrier(comm3D,code)
#           endif
!           ====================================
!           ====================================
#        endif
!        ________________________________________________________
!        GMRES Method 
#        ifdef KeyGMRES
            phi  = phiOld
            phiv = phiOldv
            call gmres3D(phi,phiv,                    &
                         rhs,Gamx,Gamy,Gamz,          &
                         xc,yc,sig,dsig,No_cp,nbe,    &
                         xv,yv,sigv,dsigv,No_vp,nbev, &
                         Hpr,h,etan,                  &
                         Hprv,hv,etav)
#        endif  
!        _________________________________________________________
!        PD-SOR: Pseudo Time Method 
#        ifdef KeyPseudoTime
            phi  = phiOld
            phiv = phiOldv
            call PseudoSOR3D(phi,phiv,                    &
                             rhs,Gamx,Gamy,Gamz,          &
                             xc,yc,sig,dsig,No_cp,nbe,    &
                             xv,yv,sigv,dsigv,No_vp,nbev, &
                             Hpr,h,etan,                  &
                             Hprv,hv,etav)
#        endif
!        ________________________________________________________
!        Schwarz SOR Method 
#        ifdef KeySchwarzSOR
            phi  = phiOld
            phiv = phiOldv
            call SchwarzSOR(phi,phiv,                    &
                            rhs,Gamx,Gamy,Gamz,          &
                            xc,yc,sig,dsig,No_cp,nbe,    &
                            xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,etan,                  &
                            Hprv,hv,etav)
#        endif
!        ________________________________________________________
!        CUDA METHODS: SOR

!        *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!        =========  PARALLEL CUDA ============
#        ifdef KeyCUDA
            phi  = phiOld
            phiv = phiOldv
!           ___________________________
!           MSOR & Jacobi (2017)             
#           if defined(KeyMSOR_cuda) || defined(KeyJacobi_cuda)
               call sor3Dcuda(phi,phiv,                    &
                              rhs,Gamx,Gamy,Gamz,          &
                              xc,yc,sig,dsig,No_cp,nbe,    &
                              xv,yv,sigv,dsigv,No_vp,nbev, &
                              Hpr,h,etan,                  &
                              Hprv,hv,etav)
#           endif
!           ___________________________
!           JSOR
#           ifdef KeyJSOR_cuda
               call sor3Dcuda(phi,phiv,                    &
                              rhs,Gamx,Gamy,Gamz,          &
                              xc,yc,sig,dsig,No_cp,nbe,    &
                              xv,yv,sigv,dsigv,No_vp,nbev, &
                              Hpr,h,etan,                  &
                              Hprv,hv,etav)
#           endif
!           ___________________________
!           PD-SOR (Pseudo Time Method)
#           ifdef KeyPDMSOR_cuda
               call PDsor3Dcuda(phi,phiv,                    &
                                rhs,Gamx,Gamy,Gamz,          &
                                xc,yc,sig,dsig,No_cp,nbe,    &
                                xv,yv,sigv,dsigv,No_vp,nbev, &
                                Hpr,h,etan,                  &
                                Hprv,hv,etav)
#           endif
#        endif
!        =====================================
!        *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*   
                              
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t5,60a)'), '<---- End   subroutine: Poisson'
      write(*,*) ' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                             END OF Poisson                          !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
