!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!            MULTI-COLORING SOR FOR THE 3D POISSON EQUATION           !
!                              Sept 2014                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE Multisor3D(phi,phiv,                    &
                            rhs,Gamx,Gamy,Gamz,          &
                            xc,yc,sig,dsig,No_cp,nbe,    &
                            xv,yv,sigv,dsigv,No_vp,nbev, &
                            Hpr,h,etan,                  &
                            Hprv,hv,etav)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine solves a linear system using the Multicoloring   !
!    S.O.R. technique. We need the file corresponding to the color    !
!    of each cell as well as the number of colors.                    !
!    We use different relaxion factors: relaxSOR. This system is good !
!    for the Poisson problem because it consider the cell-center "Am" !
!    and the vertex "Bmv" matrix coefficients.                        !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | <--> phi   |(N_CELL,NZ)  | Solution & initial guess           |  !
!  | <--> phiv  |(N_VERT,NZ-1)| Solution at the vertex values      |  !
!  |____________|_____________|____________________________________|  !
!  | ---> rhs   |(N_CELL,NZ)  | Right-hand side of the system      |  !
!  | ---> Gamx  |(N_CELL,NZ)  | Diffusive coefficient in x         |  ! 
!  | ---> Gamy  |(N_CELL,NZ)  | Diffusive coefficient in y         |  ! 
!  | ---> Gamz  |(N_CELL,NZ)  | Diffusive coefficient in z         |  ! 
!  |____________|_____________|____________________________________|  !
!  | ---> xc,yc | (N_CELL)    | Coordinates of the cell centers    |  !
!  | ---> sig   | (NZ)        | sigma value at the cell centers    |  !
!  | ---> dsig  | (NZ)        | = sig(k+1)-sig(k+1)                |  !
!  | ---> No_cp | (N_CELL,3)  | Numbering of surrounding 3 cell-cen|  !
!  | ---> nbe   | (N_CELL0)   | Tag type cell-center               |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv | (N_VERT)    | Coordinates of the vertices        |  !
!  | ---> sigv  | (NZ-1)      | sigma value at the vertices        |  !
!  | ---> dsigv | (NZ-1)      | = sigv(k+1)-sigv(k)                |  !
!  | ---> No_vp | (N_CELL0,3) | Numbering of the 3 cell vertices   |  !
!  | ---> nbev  | (N_VERT)    | Tag type of cell vertex            |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |     Name   |    Size    |  Description                        |  !  
!  |____________|____________|_____________________________________|  !
!  |    Am0     |(N_CELL0,NZ)| matrix coefficient of element i     |  !
!  |    Am1     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 1 |  ! 
!  |    Am2     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 2 |  ! 
!  |    Am3     |(N_CELL0,NZ)| matrix coeff. horizontal neigborn 3 |  ! 
!  |    AmT     |(N_CELL0,NZ)| matrix coeff. vertical top          |  ! 
!  |    AmB     |(N_CELL0,NZ)| matrix coeff. vertical bottom       |  !
!  |____________|____________|_____________________________________|  !
!  |    Bmv1T   |(N_CELL0,NZ)| matrix coeff. vertex 1 top          |  ! 
!  |    Bmv2T   |(N_CELL0,NZ)| matrix coeff. vertex 2 top          |  ! 
!  |    Bmv3T   |(N_CELL0,NZ)| matrix coeff. vertex 3 top          |  ! 
!  |    Bmv1B   |(N_CELL0,NZ)| matrix coeff. vertex 1 bottom       |  ! 
!  |    Bmv2B   |(N_CELL0,NZ)| matrix coeff. vertex 2 bottom       |  ! 
!  |    Bmv3B   |(N_CELL0,NZ)| matrix coeff. vertex 3 bottom       |  !
!  |____________|____________|_____________________________________|  !  
!  |    Newrhs  |(N_CELL0,NZ)| right hand side of the method       |  !  
!  |____________|____________|_____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - diffusion3D                ( diffusion3D.F90 )            |  !
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - BCcellcenter3D             ( BCcellcenter3D.F90 )         |  !
!  |   - BCvertex3D                 ( BCvertex3D.F90 )             |  !
!  |_______________________________________________________________|  !
!                                                                     !
!   --->  Input variables                                             !
!   <---  Output variables                                            !
!                                                                     !
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
!     -------------------------------------
      real*8,dimension(:,:) :: rhs(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8,dimension(:,:) :: Gamz(N_CELL,NZ)
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
      real*8, dimension(:)  :: etan(N_CELL)
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:) :: phiNew(N_CELL,NZ)
      real*8,dimension(:,:) :: Am0(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am1(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am2(N_CELL0,NZ)
      real*8,dimension(:,:) :: Am3(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmT(N_CELL0,NZ)
      real*8,dimension(:,:) :: AmB(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3T(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv1B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv2B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Bmv3B(N_CELL0,NZ)
      real*8,dimension(:,:) :: Nm(N_CELL0,NZ)
      real*8,dimension(:,:) :: Newrhs(N_CELL0,NZ)
      real*8,dimension(:,:) :: suma(N_CELL0,NZ)
      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     --------------------------------------
      real*8  :: errorsys,residu,som,SUMerrorsys
      real*8  :: errorNeum,SUMerrorNeum
      integer :: iter
!     ----------------------------------------
      integer :: jj,ivert
      real*8  :: nnx,nny,nnz,x,y,z
      real*8  :: dfBdn,Neumanndfdn3D
!     ---------------------------------------- 
      integer :: s
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: elem,Iini
#     endif
!     =============== END ================    
!     ====================================  

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: PSOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |                     Initial guess                      |
!     |________________________________________________________|

!     ______________________________________________________
!     Cell-centers BC   
!     ----------------------------------------------
!     Boundary Condition of the cell-center points  
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3D(phi)
#     endif	
!     =============== END ================    
!     ====================================
!     ______________________________________________________
!     Vertex interpolation & BC
!     ------------------------------------------------
!     Interpolation of the inside vertex points 
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)
!     ------------------------------------------------
!     Boundary Conditions of the vertex points 
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                      phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
!     ________________________________________________________
!    |                                                        |
!    |            Matrix Am & Bm of the diffusion term        |
!    |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0 	
            Am0(i,k)  = 0.0d0 
            Am1(i,k)  = 0.0d0
            Am2(i,k)  = 0.0d0
            Am3(i,k)  = 0.0d0
            AmT(i,k)  = 0.0d0
            AmB(i,k)  = 0.0d0
            Bmv1T(i,k)= 0.0d0 
            Bmv2T(i,k)= 0.0d0 
            Bmv3T(i,k)= 0.0d0 
            Bmv1B(i,k)= 0.0d0 
            Bmv2B(i,k)= 0.0d0 
            Bmv3B(i,k)= 0.0d0 
         enddo
      enddo

!     ________________________________________________________      
!     Diffusion with Dirichlet BC
      call diffusion3D(Am0,Am1,Am2,Am3,AmT,AmB,             & 
                       Bmv1T,Bmv2T,Bmv3T,                   &
                       Bmv1B,Bmv2B,Bmv3B,                   &
                       Gamx,Gamy,Gamz,                      &
                       xc,yc,sig,dsig,No_cp,nbe,            &
                       xv,yv,sigv,dsigv,No_vp,nbev)

      
!*********************************************************************!
!                                                                     !
!                 Solution of the system (M.S.0.R.) 3D                !
!                                                                     !
!*********************************************************************!

!     ________________________________________________________
!    |                                                        |
!    |                    Update coefficients                 |
!    |________________________________________________________|

      do k=1,NZ
         do i=1,N_CELL0 	
            Am1(i,k)  = Am1(i,k)/Am0(i,k)
            Am2(i,k)  = Am2(i,k)/Am0(i,k)
            Am3(i,k)  = Am3(i,k)/Am0(i,k)
            AmT(i,k)  = AmT(i,k)/Am0(i,k)
            AmB(i,k)  = AmB(i,k)/Am0(i,k)
            Bmv1T(i,k)= Bmv1T(i,k)/Am0(i,k) 
            Bmv2T(i,k)= Bmv2T(i,k)/Am0(i,k) 
            Bmv3T(i,k)= Bmv3T(i,k)/Am0(i,k) 
            Bmv1B(i,k)= Bmv1B(i,k)/Am0(i,k) 
            Bmv2B(i,k)= Bmv2B(i,k)/Am0(i,k) 
            Bmv3B(i,k)= Bmv3B(i,k)/Am0(i,k) 
            rhs(i,k)  = rhs(i,k)/Am0(i,k)
         enddo
      enddo

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call communication3D(rhs)
#     endif	
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |             Initial calculations of the loop           |
!     |________________________________________________________|


      iter=0
111   continue
      iter=iter+1 
!      ________________________________________________________
!     |                                                        |
!     |   MSOR: Solution of the system by SOR for each color   |
!     |________________________________________________________|


!     ____________________________________________________
!     Old terms: Bm(vertex) + Am(old) 
      do i=1,N_CELL0
         jv1 = No_vp(i,1)
         jv2 = No_vp(i,2)
         jv3 = No_vp(i,3)
         do k=2,NZ-1
            suma(i,k) = rhs(i,k)-( Bmv1T(i,k)*phiv(jv1,k)   &
                                  +Bmv2T(i,k)*phiv(jv2,k)   &
                                  +Bmv3T(i,k)*phiv(jv3,k)   &
                                  +Bmv1B(i,k)*phiv(jv1,k-1) &
                                  +Bmv2B(i,k)*phiv(jv2,k-1) &
                                  +Bmv3B(i,k)*phiv(jv3,k-1))
         enddo
      enddo

!     ____________________________________________________
!     Am(center) 
      errorsys = 0.0d0
      DO s=0,N_COLOR-1
         do i=1,N_CELL0
            if (ColorCell(i).eq.s) then
               jc1 = No_cp(i,1)
               jc2 = No_cp(i,2)
               jc3 = No_cp(i,3)
               do k=2,NZ-1
                  residu = suma(i,k)-(Am1(i,k)*phi(jc1,k)      &
                                     +Am2(i,k)*phi(jc2,k)      &
                                     +Am3(i,k)*phi(jc3,k)      &
                                     +AmT(i,k)*phi(i,k+1)      &       
                                     +AmB(i,k)*phi(i,k-1)      &
                                     +phi(i,k))
	          errorsys = errorsys + abs(residu)
	          phi(i,k) = phi(i,k) + relaxSOR*residu
               enddo
            endif
         enddo

!        _____________________________________________________
!        Parallel communication

!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call communication3Dtype1(phi)
#        endif	
!        =============== END ================    
!        ====================================

      ENDDO

!     _____________________________________________________
!     Boundary conditions cell-centers
      call BCcellcenter3D(phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)

!     ________________________________________________________
!     Parallel error

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call SUM_parallel(errorsys,SUMerrorsys)
         errorsys = SUMerrorsys
#     endif	
!     =============== END ================    
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |                Update vertex values                    |
!     |________________________________________________________|
 
!     ________________________________________________________
!     Interpolation of the inside vertex points 
      call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                           phi,xc,yc,sig,dsig,No_cp,nbe)

!     ________________________________________________________
!     Boundary Conditions of the vertex points 
      call BCVertex3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,Hprv,hv,etav,&
                      phi,xc,yc,sig,dsig,No_cp,nbe,Hpr,h,etan)
!      ________________________________________________________
!     |                                                        |
!     |                  Convergence criteria                  |
!     |________________________________________________________|
        
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         if (errorsys.lt.eps) then
            write(*,*) ' '
            write(*,7) 'Solution M.S.0.R. 3D : iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
         elseif (errorsys.gt.1.0d5) then
            write(*,*) ' '
            write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
            !stop
         elseif(iter.gt.MaxIters) then
            write(*,*) ' '
            write(*,7) 'Non-convergence: iters =',iter,&
                    ', error =',errorsys
            write(*,*) ' '
         else
            if (0.eq.mod(iter,1000)) then
               print*, 'iter=',iter,'Error=',errorsys
            endif
            goto 111
         endif
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         if (errorsys.lt.eps) then
            IF (rang_topo.eq.0) THEN
            write(*,*) ' '
            write(*,*) '        PROCESSOR:',rang_topo
            write(*,7) 'Solution M.S.0.R. 3D : iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
            ENDIF
         elseif (errorsys.gt.1.0d5) then
            IF (rang_topo.eq.0) THEN
            write(*,*) ' '
            write(*,7) 'DIVERGENCE !!!!: iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
            ENDIF
            stop
         elseif(iter.gt.MaxIters) then
            IF (rang_topo.eq.0) THEN
            write(*,*) ' '
            write(*,7) 'Non-convergence: iters =',iter,&
                       ', error =',errorsys
            write(*,*) ' '
            ENDIF
         else
            IF (rang_topo.eq.0) THEN
            if (0.eq.mod(iter,1000)) then
               print*, 'iter=',iter,'Error=',errorsys
            endif
            ENDIF
            goto 111
         endif
#     endif
!     =============== END ================    
!     ==================================== 

119   continue

      7 format(t10,a32,i5,a9,e10.3)
      6 format(t10,a14,i4,a9,e10.3)

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: MultiSOR 3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               End of Multicoloring S.O.R. Method 3D                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
