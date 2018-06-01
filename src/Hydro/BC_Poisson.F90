!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           BOUNDARY CONDITIONS OF THE PRESSURE (POISSON EQN)         !
!                               Dic 2017                              !
!           SUBROUTINES:                                              !
!                           -  BCfunction3D                           !
!                           -  BCpressure3D                           !
!                           -  NeumannVertexBC (not good so far)      !
!                                                                     !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCfunction3D(funfB,funfBv,fundfdnB,fundfdnBv, &
                              xc,yc,sig,dsig,No_cp,nbe,        &
                              xv,yv,sigv,dsigv,No_vp,nbev,     &
                              Hpr,h,etan,                      &
                              Hprv,hv,etav)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine calculates the known functions at the boundary   !
!    for the pressure that it is applied to the Poisson equation.     !
!    It depends on Dirichlet, Neumann and Mix BC. They are selected   !
!    directly from the cppdefs.h file using the following keys:       !
!                                                                     !
!                   #       define Key_DirichletBCp                   !
!                   #       define Key_NeumannBCp                     !
!                   #       define Key_MixBCp                         !
!                                                                     !
!     REMARK: This way to define the function is only for efficiency. !
!             It belongs to the BCpressure3D; however, we only need   !
!             them once. Thus it is better to calculate separeted.    !  
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name     |    Size     |  Description                     |  !  
!  |______________|_____________|__________________________________|  !
!  | --> funfB    |(N_CELL,NZ)  | Solution & initial guess         |  !
!  | --> funfBv   |(N_VERT,NZ-1)| Solution at the vertex values    |  !
!  | --> fundfdnB |(N_CELL,NZ)  | Solution & initial guess         |  !
!  | --> fundfdnBv|(N_VERT,NZ-1)| Solution at the vertex values    |  !
!  |______________|_____________|__________________________________|  !
!  | ---> rhs     |(N_CELL,NZ)  | Right-hand side of the system    |  !
!  | ---> Gamx    |(N_CELL,NZ)  | Diffusive coefficient in x       |  ! 
!  | ---> Gamy    |(N_CELL,NZ)  | Diffusive coefficient in y       |  ! 
!  | ---> Gamz    |(N_CELL,NZ)  | Diffusive coefficient in z       |  ! 
!  |______________|_____________|__________________________________|  !
!  | ---> xc,yc   | (N_CELL)    | Coordinates of the cell centers  |  !
!  | ---> sig     | (NZ)        | sigma value at the cell centers  |  !
!  | ---> dsig    | (NZ)        | = sig(k+1)-sig(k+1)              |  !
!  | ---> No_cp   | (N_CELL,3)  | Numbering of surrounding 3 cell-c|  !
!  | ---> nbe     | (N_CELL0)   | Tag type cell-center             |  !
!  |______________|_____________|__________________________________|  !
!  | ---> xv,yv   | (N_VERT)    | Coordinates of the vertices      |  !
!  | ---> sigv    | (NZ-1)      | sigma value at the vertices      |  !
!  | ---> dsigv   | (NZ-1)      | = sigv(k+1)-sigv(k)              |  !
!  | ---> No_vp   | (N_CELL0,3) | Numbering of the 3 cell vertices |  !
!  | ---> nbev    | (N_VERT)    | Tag type of cell vertex          |  !
!  |______________|_____________|__________________________________|  !
!                                                                     !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - NeumannVertexBC            ( BC_Poisson.F90 )             |  !
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

      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8,dimension(:,:) :: funfBv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: fundfdnB(N_CELL,NZ)
      real*8,dimension(:,:) :: fundfdnBv(N_VERT,NZ-1)
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

      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     ----------------------------------------
      integer :: ii,jj
      real*8  :: x,y,z,zT,zB
      real*8  :: fB,dfBdn,Neumanndfdn3D
      real*8  :: funSolExam3D,funDiffExam3D 
      real*8  :: funExamNSp,funExamNSrhsp,NeumanndpdnNS
      real*8  :: nnx,nny,nnz

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: BCfunction3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

      funfB     = 0.0d0         
      funfBv    = 0.0d0
      fundfdnB  = 0.0d0         
      fundfdnBv = 0.0d0
      
!      ________________________________________________________
!     |                                                        |
!     |         Function fB at the boundary (Dirichlet)        |
!     |________________________________________________________|

#     if defined(Key_DirichletBCp) || defined(Key_MixBCp)      
!        ______________________________________________________
!        Cell-center points
!        ___________
!        Vertical
         do i=1,N_CELL0
            x  = xc(i)
            y  = yc(i) 
            zB = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)
            zT = 0.5d0*(sig(NZ-1)+sig(NZ))*Hpr(i)-h(i)
            funfB(i,1) = funExamNSp(x,y,zB,time)
            funfB(i,NZ)= funExamNSp(x,y,zT,time)               
         enddo
!        ___________
!        Horizontal       
         do ii=N_CELL0+1,N_CELLexact
            i = No_cp(ii,1)
            j = No_cp(ii,2)
            x = xe(i,j)
            y = ye(i,j)
            do k=1,NZ 
               z = sig(k)*Hpr(i)-h(i)
               funfB(ii,k) = funExamNSp(x,y,z,time)            
            enddo
         enddo
!        ______________________________________________________
!        Vertex points
!        ___________
!        Vertical
         do nv=1,N_VERT
            x  = xv(nv)
            y  = yv(nv)
            zB = sigv(1)*Hprv(nv)-hv(nv)
            zT = sigv(NZ-1)*Hprv(nv)-hv(nv)
            funfBv(nv,1)    = funExamNSp(x,y,zB,time) 
            funfBv(nv,NZ-1) = funExamNSp(x,y,zT,time)                  
         enddo
!        ___________
!        Horizontal
         do k=1,NZ-1 
            do nv=1,N_VERT
               !if (nbev(nv).ne.0) then
               x = xv(nv)
               y = yv(nv)
               z = sigv(k)*Hprv(nv)-hv(nv)
               funfBv(nv,k) = funExamNSp(x,y,z,time)
               !endif
            enddo
         enddo
!        ______________________________________________________
!        SPECIAL CASE: Standing Wave         
#        ifndef KeyFixedFreeSurface
            funfB  = 0.0d0  !<-- Only used at the top
            funfBv = 0.0d0  !<-- Only used at the top
#        endif 
         
#     endif

!      ________________________________________________________
!     |                                                        |
!     |         Function fB at the boundary (Neumann)          |
!     |________________________________________________________|

#     if defined(Key_NeumannBCp) || defined(Key_MixBCp)
!        ______________________________________________________
!        Cell-center points
!        ___________
!        Vertical
         do i=1,N_CELL0
            x  = xc(i)
            y  = yc(i) 
            zB = 0.5d0*(sig(1)+sig(2))*Hpr(i)-h(i)
            zT = 0.5d0*(sig(NZ-1)+sig(NZ))*Hpr(i)-h(i)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0*Hpr(i)
            fundfdnB(i,1) = NeumanndpdnNS(x,y,zB,time,nnx,nny,nnz)
            nnz =  1.0d0*Hpr(i)
            fundfdnB(i,NZ)= NeumanndpdnNS(x,y,zT,time,nnx,nny,nnz)               
         enddo
!        ___________
!        Horizontal       
         do ii=N_CELL0+1,N_CELLexact
            i = No_cp(ii,1)
            j = No_cp(ii,2)
            x = xe(i,j)
            y = ye(i,j)
            nnx = normxc(i,j)
            nny = normyc(i,j)
            nnz = 0.0d0                  
            do k=1,NZ 
               z = sig(k)*Hpr(i)-h(i)
               fundfdnB(ii,k) = NeumanndpdnNS(x,y,z,time,nnx,nny,nnz)            
            enddo
         enddo
!        ______________________________________________________
!        Vertex points (not used in the calculation april 2017)
!        ___________
!        Vertical
         do nv=1,N_VERT
            x  = xv(nv)
            y  = yv(nv)
            zB = sigv(1)*Hprv(nv)-hv(nv)
            zT = sigv(NZ-1)*Hprv(nv)-hv(nv)
            nnx =  0.0d0
            nny =  0.0d0
            nnz = -1.0d0         
            fundfdnBv(nv,1) = NeumanndpdnNS(x,y,zB,time,nnx,nny,nnz)
            nnz =  1.0d0         
            fundfdnBv(nv,NZ-1) = NeumanndpdnNS(x,y,zT,time,nnx,nny,nnz)                  
         enddo
!        ___________
!        Horizontal
         do k=1,NZ-1 
            do nv=1,N_VERT
               if (nbev(nv).ne.0) then         
                  x = xv(nv)
                  y = yv(nv)
                  z = sigv(k)*Hprv(nv)-hv(nv)
                  nnx = normxv(nv)
                  nny = normyv(nv)
                  nnz = 0.0d0               
                  fundfdnBv(nv,k) = NeumanndpdnNS(x,y,z,time,nnx,nny,nnz)
               endif
            enddo
         enddo
!        ______________________________________________________
!        SPECIAL CASE: Standing Wave         
#        ifndef KeyFixedFreeSurface
            fundfdnB  = 0.0d0
            fundfdnBv = 0.0d0
            do i=1,N_CELL0
               fundfdnB(i,1) = -rho_f*Hpr(i)*dw1dt(i)              
            enddo 
#        endif 
         
#     endif  


!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: BCfunction3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           BOUNDARY CONDITIONS OF THE PRESSURE (POISSON EQN)         !
!                               Nov 2017                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE BCpressure3D(phi,phiv,                        &
                              funfB,funfBv,fundfdnB,fundfdnBv, &
                              xc,yc,sig,dsig,No_cp,nbe,        &
                              xv,yv,sigv,dsigv,No_vp,nbev,     &
                              Hpr,h,etan,                      &
                              Hprv,hv,etav,                    &
                              timeBCcc,timeBCvv,timeInter)                

!---------------------------------------------------------------------!
!                                                                     !
!    This subroutine applies the boundary conditions for the pressure.!
!    It is applied to the Poisson equation. We have three options     !
!    to choose: Dirichlet, Neumann and Mix. They can be selected      !
!    directly from the cppdefs.h file using the following keys:       !
!                                                                     !
!                   #       define Key_DirichletBCp                   !
!                   #       define Key_NeumannBCp                     !
!                   #       define Key_MixBC                          !
!                                                                     !
!     Notice that an important requirent for an efficient way to      !
!     implement the boundaries are the predefinitions of the functions!
!     used at the boundary: funfB,funfBv,fundfdnB,fundfdnBv.          !  
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |     Name     |    Size     |  Description                     |  !  
!  |______________|_____________|__________________________________|  !
!  | <--> phi     |(N_CELL,NZ)  | Solution & initial guess         |  !
!  | <--> phiv    |(N_VERT,NZ-1)| Solution at the vertex values    |  !
!  |______________|_____________|__________________________________|  !
!  | --> funfB    |(N_CELL,NZ)  | Solution & initial guess         |  !
!  | --> funfBv   |(N_VERT,NZ-1)| Solution at the vertex values    |  !
!  | --> fundfdnB |(N_CELL,NZ)  | Solution & initial guess         |  !
!  | --> fundfdnBv|(N_VERT,NZ-1)| Solution at the vertex values    |  !
!  |______________|_____________|__________________________________|  !
!  | ---> rhs     |(N_CELL,NZ)  | Right-hand side of the system    |  !
!  | ---> Gamx    |(N_CELL,NZ)  | Diffusive coefficient in x       |  ! 
!  | ---> Gamy    |(N_CELL,NZ)  | Diffusive coefficient in y       |  ! 
!  | ---> Gamz    |(N_CELL,NZ)  | Diffusive coefficient in z       |  ! 
!  |______________|_____________|__________________________________|  !
!  | ---> xc,yc   | (N_CELL)    | Coordinates of the cell centers  |  !
!  | ---> sig     | (NZ)        | sigma value at the cell centers  |  !
!  | ---> dsig    | (NZ)        | = sig(k+1)-sig(k+1)              |  !
!  | ---> No_cp   | (N_CELL,3)  | Numbering of surrounding 3 cell-c|  !
!  | ---> nbe     | (N_CELL0)   | Tag type cell-center             |  !
!  |______________|_____________|__________________________________|  !
!  | ---> xv,yv   | (N_VERT)    | Coordinates of the vertices      |  !
!  | ---> sigv    | (NZ-1)      | sigma value at the vertices      |  !
!  | ---> dsigv   | (NZ-1)      | = sigv(k+1)-sigv(k)              |  !
!  | ---> No_vp   | (N_CELL0,3) | Numbering of the 3 cell vertices |  !
!  | ---> nbev    | (N_VERT)    | Tag type of cell vertex          |  !
!  |______________|_____________|__________________________________|  !
!                                                                     !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - interpolation3D            ( interpolation3D.F90 )        |  !
!  |   - NeumannVertexBC            ( BC_Poisson.F90 )             |  !
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
      real*8,dimension(:,:) :: funfB(N_CELL,NZ)
      real*8,dimension(:,:) :: funfBv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: fundfdnB(N_CELL,NZ)
      real*8,dimension(:,:) :: fundfdnBv(N_VERT,NZ-1)
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
!     --------------------------------------
      real :: timeBCcc,timeBCvv,timeInter      
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      integer:: jc1,jc2,jc3
      integer:: jv1,jv2,jv3
!     ----------------------------------------
      integer :: ii,jj,nvp,i0,elem
      real*8  :: x,y,z,zT,zB
      real*8  :: fB,dfBdn,Neumanndfdn3D
      real*8  :: funSolExam3D,funDiffExam3D 
      real*8  :: funExamNSp,funExamNSrhsp,NeumanndpdnNS
      real*8  :: nnx,nny,nnz
!     ---------------------------------------- 
      real :: start,finish,startT,finishT

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '----> Begin subroutine: BCpressure3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |             Boundary conditions (Dirichlet)            |
!     |________________________________________________________|
     
#     ifdef Key_DirichletBCp
!         ____________________________________
!        |                                    |
!        |           Cell-centers BC          |
!        |____________________________________|
         call cpu_time(start)     
!        ----------                  
!        Vertical
         do i=1,N_CELL0
            phi(i,1) = 2.0d0*funfB(i,1)-phi(i,2)
            phi(i,NZ)= 2.0d0*funfB(i,NZ)-phi(i,NZ-1)
         enddo
!        ----------                 
!        Horizontal       
         do ii=N_CELL0+1,N_CELLexact
            i = No_cp(ii,1)
            do k=1,NZ     
               phi(ii,k) = 2.0d0*funfB(ii,k)-phi(i,k)
            enddo
         enddo
         call cpu_time(finish)
         timeBCcc = timeBCcc + (finish-start) 
!         ____________________________________
!        |                                    |
!        |               Vertex BC            |
!        |____________________________________| 
#        ifdef KeyIterationWithVertex                           
!        --------------------------------------
!        Vertex interpolation  
         call cpu_time(start)
         call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              phi,xc,yc,sig,dsig,No_cp,nbe)
         call cpu_time(finish)
         timeInter = timeInter + (finish-start) 
            
!        --------------------------------------
!        Boundary conditions
         call cpu_time(start)
         do nv=1,N_VERT
            phiv(nv,1)    = funfBv(nv,1)
            phiv(nv,NZ-1) = funfBv(nv,NZ-1)
            if (nbev(nv).ne.0) then
               do k=2,NZ-2 
                  phiv(nv,k) = funfBv(nv,k)
               enddo
            endif
         enddo
         call cpu_time(finish)
         timeBCvv = timeBCvv + (finish-start)
#        endif
  
#     endif /* Key_DirichletBCp */

!      ________________________________________________________
!     |                                                        |
!     |              Boundary conditions (Neumann)             |
!     |________________________________________________________|
      
#     ifdef Key_NeumannBCp 
!         ____________________________________
!        |                                    |
!        |           Cell-centers BC          |
!        |____________________________________|

         call cpu_time(start)     
!        ----------                  
!        Vertical
         do i=1,N_CELL0
            phi(i,1) = dsig(1)*fundfdnB(i,1) + phi(i,2)
            phi(i,NZ)= dsig(NZ-1)*fundfdnB(i,NZ) + phi(i,NZ-1)
         enddo
!        ----------                 
!        Horizontal       
         do ii=N_CELL0+1,N_CELLexact
            i = No_cp(ii,1)
            j = No_cp(ii,2)
            do k=1,NZ     
               phi(ii,k) = 2.0d0*dlCE(i,j)*fundfdnB(ii,k) +phi(i,k)
            enddo
         enddo
         call cpu_time(finish)
         timeBCcc = timeBCcc + (finish-start)
!         ____________________________________
!        |                                    |
!        |               Vertex BC            |
!        |____________________________________| 

#        ifdef KeyIterationWithVertex
!        --------------------------------------
!        Vertex interpolation 
         call cpu_time(start)
         call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              phi,xc,yc,sig,dsig,No_cp,nbe)                             
         call cpu_time(finish)
         timeInter = timeInter + (finish-start)
         
!        --------------------------------------
!        Boundary conditions
!        --------                                                          
#        ifdef KeyUseNeumannBCvertex
            call NeumannVertexBC(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                 phi,xc,yc,sig,dsig,No_cp,nbe,     &
                                 Hprv,hv,etav,                     &
                                 Hpr,h,etan) 
#        endif
!        -------- 
#        endif
                            
#     endif /* Key_NeumannBCp */

!      ________________________________________________________
!     |                                                        |
!     |  Boundary conditions ( Mix: Neumann + Dirichlet(top) ) |
!     |________________________________________________________|
      
#     ifdef Key_MixBCp
!         ____________________________________
!        |                                    |
!        |           Cell-centers BC          |
!        |____________________________________|

         call cpu_time(start) 
!        ----------                  
!        Vertical
         do i=1,N_CELL0
            phi(i,1) = dsig(1)*fundfdnB(i,1) + phi(i,2)
            phi(i,NZ)= 2.0d0*funfB(i,NZ)-phi(i,NZ-1) !<<< Dirichlet at the top           
         enddo
!        ----------                 
!        Horizontal       
         do ii=N_CELL0+1,N_CELLexact
            i = No_cp(ii,1)
            j = No_cp(ii,2)
            do k=1,NZ     
               phi(ii,k) = 2.0d0*dlCE(i,j)*fundfdnB(ii,k) + phi(i,k)
            enddo
         enddo
!        ----------          
         call cpu_time(finish)
         timeBCcc = timeBCcc + (finish-start)
!         ____________________________________
!        |                                    |
!        |               Vertex BC            |
!        |____________________________________|

#        ifdef KeyIterationWithVertex
!        --------------------------------------
!        Vertex interpolation 
         call cpu_time(start)
         call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              phi,xc,yc,sig,dsig,No_cp,nbe)
         call cpu_time(finish)
         timeInter = timeInter + (finish-start)
         
!        --------------------------------------
!        Boundary conditions
         call cpu_time(start)
         do nv=1,N_VERT
            phiv(nv,NZ-1) = funfBv(nv,NZ-1) !<<<< Dirichlet at the top
         enddo
         call cpu_time(finish)
         timeBCvv = timeBCvv + (finish-start)
#        endif
            
#     endif /* Key_MixBCp */

!      ________________________________________________________
!     |                                                        |
!     |             Boundary conditions (Periodic)             |
!     |________________________________________________________|
     
#     ifdef Key_PeriodicBCp
!         ____________________________________
!        |                                    |
!        |           Cell-centers BC          |
!        |____________________________________|

         call cpu_time(start)     
!        ----------                  
!        Vertical
         do i=1,N_CELL0
            phi(i,1) = phi(i,2)
            phi(i,NZ)= phi(i,NZ-1)
         enddo
!        ----------                 
!        Horizontal
       
         DO ii=N_CELL0+1,N_CELLexact
            i0 = No_cp(ii,3)
!           ====================================
!           ==========  SEQUENTIAL =============
#           ifndef KeyParallel
            if (i0.ge.1.AND.i0.le.N_CELL0) then
               do k=1,NZ     
                  phi(ii,k) = phi(i0,k)
               enddo
            endif
!           ====================================
!           ===========  PARALLEL ==============
#           else
               elem = index_global(i0)
               if (elem.ge.1.AND.elem.le.N_CELL0global) then
                  do k=1,NZ
                     phi(ii,k) = phi(i0,k)
                  enddo
               endif
#           endif
!           ====================================
!           ====================================
         ENDDO

#        ifdef KeyNoHacer
            DO i=1,N_CELL0
               do j=1,3
                  nc=No_cp(i,j)
!                 ====================================
#                 ifndef KeyParallel
!                 ==========  SEQUENTIAL =============
                  if (nc.lt.1.OR.nc.gt.N_CELL0) then
                     i0 = No_cp(nc,3)
                     if (i0.ge.1.AND.i0.le.N_CELL0) then
                        do k=1,NZ
                           phi(nc,k)  = phi(i0,k)
                        enddo
                     endif
                  endif
!                 ====================================
#                 else
!                 =====  START PARALLEL OPTION =======
                  elem = index_global(nc)
                  if (elem.lt.1.OR.elem.gt.N_CELL0global) then
                     i0 = No_cp(nc,3)
                     elem = index_global(i0)
                     if (elem.ge.1.AND.elem.le.N_CELL0global) then
                        do k=1,NZ
                           phi(nc,k)  = phi(i0,k)
                        enddo
                     endif
                  endif
!                 =============== END ================
#                 endif
!                 ====================================  
               enddo                                
            ENDDO
#        endif

!        ---------- 
         call cpu_time(finish)
         timeBCcc = timeBCcc + (finish-start) 
         
!         ____________________________________
!        |                                    |
!        |               Vertex BC            |
!        |____________________________________|

#        ifdef KeyIterationWithVertex
!        --------------------------------------
!        Vertex interpolation 
         call cpu_time(start)
         call interpolation3D(phiv,xv,yv,sigv,dsigv,No_vp,nbev,&
                              phi,xc,yc,sig,dsig,No_cp,nbe)
         call cpu_time(finish)
         timeInter = timeInter + (finish-start)
         
!        --------------------------------------
!        Boundary conditions
         call cpu_time(start)
!        ====================================
!        ====== SEQUENTIAL & PARALLEL* ======
         DO nv=1,N_VERT
            if (TagPeriodicBC(nv,1).eq.1) then
               nvp = TagPeriodicBC(nv,2)
               do k=1,NZ-1 
                  phiv(nv,k) = phiv(nvp,k)
               enddo
            endif   
         ENDDO
!        =============== END ================
!        ====================================
         call cpu_time(finish)
         timeBCvv = timeBCvv + (finish-start)
#        endif


#     endif /* Key_PeriodicBCp */

!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<---- End   subroutine: BCpressure3D'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  NEUMANN VERTEX BOUNDARY CONDITION                  !
!           *** This is not implemented successfully so far ***       !
!                               Oct 2017                              !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE NeumannVertexBC(phiv,xv,yv,sigv,dsigv,No_vp,nbev, &
                                 phi,xc,yc,sig,dsig,No_cp,nbe,     &
                                 Hprv,hv,etav,                     &
                                 Hpr,h,eta)

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the correct vertex boundary condition    !
!    using Newumann type.                                             !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |     Size    | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <-- phiv    |(N_VERT,NZ-1)| Function at the vertices          |  !
!  |_____________|_________________________________________________|  !
!  | * i,j,nv,k  |   Loop counters                                 |  !    
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !  
!  | --> xc,yc   |(N_CELL)   | Coordinates of the cell centers     |  !
!  | --> sig     |(NZ)       | sigma value at the cell centers     |  !
!  | --> dsig    |(NZ)       | = sig(k+1)-sig(k+1)                 |  !
!  | --> No_cp   |(N_CELL,3) | Numbering of surrounding 3 cell-cent|  !
!  | --> nbe     |(N_CELL0)  | Tag type cell-center                |  !
!  |_____________|___________|_____________________________________|  !
!  | --> xv,yv   |(N_VERT)   | Coordinates of the vertices         |  !
!  | --> sigv    |(NZ-1)     | sigma value at the vertices         |  !
!  | --> dsigv   |(NZ-1)     | = sigv(k+1)-sigv(k)                 |  !
!  | --> No_vp   |(N_VERT,3) | Numbering of the 3 cell vertices    |  !
!  | --> nbev    |(N_VERT)   | Tag type of cell vertex             |  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!    Local variables:                                                 !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !   
!  | xee,yee     |  Perpendicular intersection of each edge        |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   <->  Input and output variables                                   !
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

      real*8,dimension(:,:) :: phiv(N_VERT,NZ-1)
      real*8, dimension(:)  :: xv(N_VERT)
      real*8, dimension(:)  :: yv(N_VERT)
      real*8,dimension(:,:) :: sigv(NZ-1)
      real*8,dimension(:,:) :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)  
!     ---------------------------------
      real*8,dimension(:,:) :: phi(N_CELL,NZ)
      real*8, dimension(:)  :: xc(N_CELL)
      real*8, dimension(:)  :: yc(N_CELL)
      real*8,dimension(:,:) :: sig(NZ)
      real*8,dimension(:,:) :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     --------------------------------------
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: hv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)      
!     --------------------------------------
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: h(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
 
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8 ,dimension(:,:):: dphidx(N_CELL,NZ)
      real*8 ,dimension(:,:):: dphidy(N_CELL,NZ)
      real*8 :: sumfx,sumfy,deter
!     ---------------------------------
      real*8 :: dnn,nnx,nny,nnz
      real*8 :: xnn1,ynn1,xnn2,ynn2
      real*8 :: x1,y1,x2,y2,z1,z2,z3
!     -----------------------------------
      real*8 :: f0,f1,f2,h1,h2,deno
      real*8 :: x,y,z,fB,dfBdn,dfBdnv
      real*8 :: NeumanndpdnNS
      real*8 :: funSolExam2D,funSolExam3D
!     ---------------------------------
      integer :: jv1,jv2,jv3,jj,jc
      integer :: ic1,ic2,tag
      integer, parameter :: OptionNeuBC=1

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: VertexBC'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                  Vertex Neumann BC approximation                    !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |           Approximation of the gradient 2D             |
!     |________________________________________________________|

      dphidx = 0.0d0
      dphidy = 0.0d0
      do k=1,NZ
         do i=1,N_CELL0	
            sumfx = 0.0d0
            sumfy = 0.0d0
            do j=1,3
               jc = No_cp(i,j)                
               sumfx = sumfx + dxCC(i,j)*(phi(jc,k)-phi(i,k))
               sumfy = sumfy + dyCC(i,j)*(phi(jc,k)-phi(i,k))
            enddo
            deter = sum_xc2(i)*sum_yc2(i)-sum_xcyc(i)*sum_xcyc(i)
            dphidx(i,k)=(sum_yc2(i)*sumfx-sum_xcyc(i)*sumfy)/deter
            dphidy(i,k)=(sum_xc2(i)*sumfy-sum_xcyc(i)*sumfx)/deter
         enddo
      enddo

!      ________________________________________________________
!     |                                                        |
!     |                        Neumann BC                      |  
!     |________________________________________________________|

      IF (OptionNeuBC.eq.1) THEN
         tag=0
!        ---------------------------------------------------------  
!        Neumann BC approximation  
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
            tag = tag + 1
!              ---------------------------
!              Exact Neumann BC
               x = xv(nv)
               y = yv(nv)
               nnx = normxv(nv)
               nny = normyv(nv)
               do k=1,NZ-1
                  z = sigv(k)*Hprv(nv)-hv(nv)
                  nnz =  0.0d0
                  dfBdn = -NeumanndpdnNS(x,y,z,time,nnx,nny,nnz)
                  ! It is negative beacuse we approx the derivative in
                  ! the opposite direction of the outnormal                   
!                 ---------------------------
!                 Function approximations
                  ic1 = icxn1(nv)
                  ic2 = icxn2(nv)
                  f1 = phi(ic1,k) + dphidx(ic1,k)*(xn1(nv)-xc(ic1)) &
                                  + dphidy(ic1,k)*(yn1(nv)-yc(ic1))   
                  f2 = phi(ic2,k) + dphidx(ic2,k)*(xn2(nv)-xc(ic2)) &
                                  + dphidy(ic2,k)*(yn2(nv)-yc(ic2))
!                 ---------------------------
!                 Approx at the vertex point
                  h1 = 1.0d0*dn(nv)
                  h2 = 2.0d0*dn(nv)
                  deno = h1*(h2**2)-h2*(h1**2)
                  f0 =(-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)/(h2**2-h1**2)
                  phiv(nv,k) = f0 
               enddo
            endif           
         enddo
         !if (tag.ge.1) print*,'Hello',time,tag,ic1,ic2,f0
      ELSE
         
!      ________________________________________________________
!     |                                                        |
!     |                        Neumann BC                      |  
!     |________________________________________________________|

      tag = 0
      DO i=1,N_CELL0
!        --------------------------------------------------
!        Length dnn in the normal direction
         dnn = 0.55*min(dlVV(i,1),dlVV(i,2),dlVV(i,3))
         do j=1,3
            nv=No_vp(i,j)
            if (nbev(nv).ne.0) then 
!              _____________________________________________
!              Unit normal n in the out direction
               nnx = normxv(nv)
               nny = normyv(nv)
               !print*,i,'nv=',nv,'nnx,nny=',nnx,nny
!              _____________________________________________
!              Points in the normal line (inside direction)
               xnn1 = xv(nv)+dnn*(-nnx)
               ynn1 = yv(nv)+dnn*(-nny)
               xnn2 = xv(nv)+2.0d0*dnn*(-nnx)
               ynn2 = yv(nv)+2.0d0*dnn*(-nny)
               !print*,'            nv=',nv,'(xnn2,ynn2)=',xnn2,ynn2
!              _____________________________________________
!              Is the point inside the triangle?
               jv1 = No_vp(i,1)
               jv2 = No_vp(i,2)
               jv3 = No_vp(i,3)

               x1 = dxVV(i,1)
               y1 = dyVV(i,1)
               x2 = xnn2-xv(jv1)
               y2 = ynn2-yv(jv1)
               z1 = x1*y2-y1*x2

               x1 = dxVV(i,2)
               y1 = dyVV(i,2)
               x2 = xnn2-xv(jv2)
               y2 = ynn2-yv(jv2)
               z2 = x1*y2-y1*x2

               x1 = dxVV(i,3)
               y1 = dyVV(i,3)
               x2 = xnn2-xv(jv3)
               y2 = ynn2-yv(jv3)
               z3 = x1*y2-y1*x2
               !print*,'            nv=',nv,'z1,z2,z3',z1,z2,z3
!              _____________________________________________
!              If YES, approximate the vertex point
               if ((z1.ge.0).and.(z2.ge.0).and.(z3.ge.0)) then
                  !print*,'            nv=',nv,'YES'
                  do k=1,NZ-1
!                    ---------------------------
!                    Exact Neumann BC
                     x = xv(nv)
                     y = yv(nv)
                     z = sigv(k)*Hprv(nv)-hv(nv)
                     nnz =  0.0d0
                     dfBdn = -NeumanndpdnNS(x,y,z,time,nnx,nny,nnz)
!                    ---------------------------
!                    Function approximations
                     f1 = phi(i,k) + dphidx(i,k)*(xnn1-xc(i)) &
                                   + dphidy(i,k)*(ynn1-yc(i))   
                     f2 = phi(i,k) + dphidx(i,k)*(xnn2-xc(i)) &
                                   + dphidy(i,k)*(ynn2-yc(i))
!                    ---------------------------
!                    Approx at the vertex point
                     h1 = 1.0d0*dnn
                     h2 = 2.0d0*dnn
                     deno = (h1*h2*h2-h2*h1*h1)
                     f0 = (-(h1**2)*f2+(h2**2)*f1-dfBdn*deno)&
                           /(h2*h2-h1*h1)
                     phiv(nv,k) = f0
                  enddo                   
                  !print*,'            nv=',nv,'fn1,fn2=',f1,f2
                  !print*,'------------------------'
               else
                  tag = 1
                  !print*,'            nv=',nv,'NO',xv(nv),yv(nv),nnx,nny
                  !print*,'------------------------'
               endif           
            endif
         enddo
      ENDDO
      !if (tag.eq.1) print*,'Hello',time,tag,z1,z2,z3
      
      ENDIF
       
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      <----   End subroutine: VertexBC'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                      	   END OF INTERPOLATION                       !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
