!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                    ERROR OF FREE-SURFACE EXAMPLES                   !
!                              March 2014                             !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE FS_TimeError(Hpr,eta,Hprv,etav,               &
                              ufnp,vfnp,wfnp,pfnp,             &
                              ufv,vfv,wfv,pfv,                 &
                              h,xc,yc,sig,dsig,No_cp,nbe,      &
                              hv,xv,yv,sigv,dsigv,No_vp,nbev,  &
                              No_sp)                  

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the exact solution and errors of free    !
!    surface time problems. An important remark is that we display    !
!    more variables than we actually have.                            !
!                                                                     !
!---------------------------------------------------------------------!
!    NO Output  variables:                                            !
!   _______________________________________________________________   !
!  |     Name    |    Size     | Description                       |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- ufnpE  |(N_CELL,NZ)  | u Exact solution cell-center      |  !  
!  | <--- Erruf  |(N_CELL,NZ)  | u Error cell-center               |  !    
!  | <--- ufv    |(N_VERT,NZ-1)| u Approximate solution vertex     |  !
!  | <--- ufvE   |(N_VERT,NZ-1)| u Exact solution vertex           |  ! 
!  | <--- Errufv |(N_VERT,NZ-1)| u Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- vfnpE  |(N_CELL,NZ)  | v Exact solution cell-center      |  !  
!  | <--- Errvf  |(N_CELL,NZ)  | v Error cell-center               |  !    
!  | <--- vfv    |(N_VERT,NZ-1)| v Approximate solution vertex     |  !
!  | <--- vfvE   |(N_VERT,NZ-1)| v Exact solution vertex           |  ! 
!  | <--- Errvfv |(N_VERT,NZ-1)| v Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- wfnpE  |(N_CELL,NZ)  | w Exact solution cell-center      |  !  
!  | <--- Errwf  |(N_CELL,NZ)  | w Error cell-center               |  !    
!  | <--- wfv    |(N_VERT,NZ-1)| w Approximate solution vertex     |  !
!  | <--- wfvE   |(N_VERT,NZ-1)| w Exact solution vertex           |  ! 
!  | <--- Errwfv |(N_VERT,NZ-1)| w Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!  | <--- pfnpE  |(N_CELL,NZ)  | p Exact solution cell-center      |  !  
!  | <--- Errpf  |(N_CELL,NZ)  | p Error cell-center               |  !    
!  | <--- pfv    |(N_VERT,NZ-1)| p Approximate solution vertex     |  !
!  | <--- pfvE   |(N_VERT,NZ-1)| p Exact solution vertex           |  ! 
!  | <--- Errpfv |(N_VERT,NZ-1)| p Error vertex                    |  !  
!  |_____________|_____________|___________________________________|  !
!                                                                     !
!    Input variables:                                                 !
!   _______________________________________________________________   !
!  |     Name   |    Size     |  Description                       |  !  
!  |____________|_____________|____________________________________|  !
!  | ---> ufnp  |(N_CELL,NZ)  | u Old solution of the equation     |  !
!  | ---> vfnp  |(N_CELL,NZ)  | v Old solution of the equation     |  !  
!  | ---> wfnp  |(N_CELL,NZ)  | w Old solution of the equation     |  !  
!  | ---> pfnp  |(N_CELL,NZ)  | p Old solution of the equation     |  !    
!  | ---> xc,yc |(N_CELL)     | Coordinates of the cell centers    |  !
!  | ---> sig   |(NZ)         | Sigma value at the cell centers    |  !
!  | ---> dsig  |(NZ)         | Increment = sig(k+1)-sig(k)        |  !
!  | ---> No_cp |(N_CELL,3)   | Numbering of surrounding 3 cells   |  !
!  | ---> nbe   |(N_CELL)     | Tag: Type of cell (inside or bc)   |  !
!  |____________|_____________|____________________________________|  !
!  | ---> xv,yv |(N_VERT)     | Coordinates of the cell vertices   |  !
!  | ---> sigv  |(NZ-1)       | sigma of the vertex points         |  !
!  | ---> dsigv |(NZ-1)       | Increment = sigv(k+1)-sigv(k)      |  !  
!  | ---> No_vp |(N_CELL0,3)  | Numbering of the cell vertices     |  !
!  | ---> nbev  |(N_VERT)     | Tag: Type of vertex (inside or bc) |  !
!  |____________|_____________|____________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !  
!  |   - interpolation3D             ( interpolation3D.F90 )       |  !
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
      real*8, dimension(:)  :: Hpr(N_CELL)
      real*8, dimension(:)  :: eta(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
!     --------------------------------------    
      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)  
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)   
!     -------------------------------------
      real*8, dimension(:)  :: h(N_CELL)
      real*8,dimension(:)   :: xc(N_CELL)
      real*8,dimension(:)   :: yc(N_CELL)
      real*8,dimension(:)   :: sig(NZ)
      real*8,dimension(:)   :: dsig(NZ)
      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0) 
!     -------------------------------------
      real*8, dimension(:)  :: hv(N_VERT)
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     --------------------------------------
      integer,dimension(:)  :: No_sp(N_SPmax)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

!     --------------------------------------
      real*8, dimension(:)  :: etaE(N_CELL)
      real*8, dimension(:)  :: HprE(N_CELL)
      real*8, dimension(:)  :: etavE(N_VERT)
      real*8, dimension(:)  :: HprvE(N_VERT)       
      real*8,dimension(:,:) :: ufnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Erruf(N_CELL,NZ) 
      real*8,dimension(:,:) :: ufvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errvf(N_CELL,NZ) 
      real*8,dimension(:,:) :: vfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errvfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errwf(N_CELL,NZ) 
      real*8,dimension(:,:) :: wfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errwfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errpf(N_CELL,NZ) 
      real*8,dimension(:,:) :: pfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errpfv(N_VERT,NZ-1)
!     --------------------------------------
      real*8,dimension(:) :: Err2D(N_CELL)
      real*8,dimension(:) :: Err2Dv(N_VERT)
!     -------------------------------------
      real*8,dimension(:,:) :: phiA(N_CELL,NZ)
      real*8,dimension(:,:) :: phiE(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8,dimension(:) :: SaveErrorMax(4)
      real*8,dimension(:) :: SaveErrorSum(4)
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
      real*8 :: MAXmaxErrorA,MAXmaxErrorR
      real*8 :: SUMsumErrorA,SUMsumErrorR
      real*8 :: Peakvalue,PeakvalueEx,PeakError
!     --------------------------------------
      real*8 :: x,y,z,som,som1,som2,som3,som4,c1
      real*8 :: TimeExample2D,TimeExample3D     
      real*8 :: FS_funu,FS_funv,FS_funw,FS_funp,FS_funeta
      integer:: DisplayThis,s,ss
!     --------------------------------------
      real*8,dimension(:,:) :: aux(N_VERT,NZ)

!*********************************************************************!
!                                                                     !
!                             Initialization                          !
!                                                                     !
!*********************************************************************!

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: TestTimeErrorFS'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      IF ((time.ge.tInisave).and.&
          (1d-08.ge.dtsave-(time-LastTimeSave))) THEN

!      ________________________________________________________
!     |                                                        |
!     |            Exact solution & Absolute Errors            |
!     |________________________________________________________|

!     _________________________________________________________
!     Cell-center

      do i=1,N_CELL
         x = xc(i)
         y = yc(i)
         etaE(i)  = FS_funeta(x,y,time)
         HprE(i)  = etaE(i) + h(i)
         do k=1,NZ 
            z = sig(k)*HprE(i)-h(i)
!           -----------------------------------
            ufnpE(i,k) = FS_funu(x,y,z,time)
            vfnpE(i,k) = FS_funv(x,y,z,time)
            wfnpE(i,k) = FS_funw(x,y,z,time) 
            pfnpE(i,k) = FS_funp(x,y,z,time)
            !pfnpE(i,k) = etaE(i) 
!           -----------------------------------
            Erruf(i,k) = abs(ufnp(i,k)-ufnpE(i,k))
            Errvf(i,k) = abs(vfnp(i,k)-vfnpE(i,k))
            Errwf(i,k) = abs(wfnp(i,k)-wfnpE(i,k))
            Errpf(i,k) = abs(pfnp(i,k)-pfnpE(i,k))
            !Errpf(i,k) = abs(eta(i)-etaE(i))
         enddo
      enddo
!     -----------------------------------
!     Error BC
      do k=1,NZ  
         do i=1,N_CELL0
            if (nbe(i).ne.0) then
               Erruf(i,k) = 0.0d0
               Errvf(i,k) = 0.0d0
               Errwf(i,k) = 0.0d0 
               Errpf(i,k) = 0.0d0 
            endif
            if (k.eq.1)  Erruf(i,k) = 0.0d0
            if (k.eq.NZ) Erruf(i,k) = 0.0d0
         enddo
      enddo  

!     _________________________________________________________
!     Vertices
!     -----------------------------------
!     Exact vertex  
      do nv=1,N_VERT
         x = xv(nv)
         y = yv(nv)
         etavE(nv) = FS_funeta(x,y,time) 
         HprvE(nv) = etavE(nv) + hv(nv)
         do k=1,NZ-1
            z = sig(k)*HprvE(nv)-hv(nv)
            ufvE(nv,k) = FS_funu(x,y,z,time)
            vfvE(nv,k) = FS_funv(x,y,z,time)
            wfvE(nv,k) = FS_funw(x,y,z,time) 
            pfvE(nv,k) = FS_funp(x,y,z,time) 
         enddo
      enddo
!     -----------------------------------
!     Error (interpolation)

      DO k=1,NZ-1   
         do nv=1,N_VERT                 
             som1 = 0.0d0
             som2 = 0.0d0
             som3 = 0.0d0
             som4 = 0.0d0                                               
	         do j=1,Dimsurrounding(nv)
                nc = surrounding(nv,j)
                c1 = weight(nv,j)
                som1 = som1 + c1*(Erruf(nc,k)+Erruf(nc,k+1))
                som2 = som2 + c1*(Errvf(nc,k)+Errvf(nc,k+1))  
                som3 = som3 + c1*(Errwf(nc,k)+Errwf(nc,k+1))  
                som4 = som4 + c1*(Errpf(nc,k)+Errpf(nc,k+1))                                                                                                          
             enddo                       
             Errufv(nv,k) = 0.5d0*som1/dlVsum(nv)
             Errvfv(nv,k) = 0.5d0*som2/dlVsum(nv)
             Errwfv(nv,k) = 0.5d0*som3/dlVsum(nv)
             Errpfv(nv,k) = 0.5d0*som4/dlVsum(nv)
         enddo
      ENDDO
      
!     -----------------------------------
!     Error BC
      do k=1,NZ-1  
         do nv=1,N_VERT
            if (nbev(nv).ne.0) then
               Errufv(nv,k) = 0.0d0
               Errvfv(nv,k) = 0.0d0
               Errwfv(nv,k) = 0.0d0 
               Errpfv(nv,k) = 0.0d0 
            endif
            if (k.eq.1)    Errufv(nv,k) = 0.0d0
            if (k.eq.NZ-1) Errufv(nv,k) = 0.0d0
         enddo
      enddo
           
!      ________________________________________________________
!     |                                                        |
!     |                      Norm Errors                       |
!     |________________________________________________________|

      DisplayThis = 1
#     ifdef KeyParallel
       if (rang_topo.ne.0) then
          DisplayThis = 0
       endif
#     endif

      if (DisplayThis.eq.1) then
      write(*,8)'===================================================='
      write(*,8)'                                                    '
      write(*,8)'           TEST: Free Surface problem               '
      write(*,6)' Time = ',time,' (dt=',dt,')'
      write(*,7)' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      endif

      DO s=1,5
!        ________________________________________________________
!        Norm errors 
         maxErrorA = 0.0d0
         maxErrorR = 0.0d0
         sumErrorA = 0.0d0
         sumErrorR = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
               if (s.eq.1) then
                  phiA(i,k) = ufnp(i,k)
                  phiE(i,k) = ufnpE(i,k)
               elseif (s.eq.2) then
                  phiA(i,k) = vfnp(i,k)
                  phiE(i,k) = vfnpE(i,k)
               elseif (s.eq.3) then
                  phiA(i,k) = wfnp(i,k)
                  phiE(i,k) = wfnpE(i,k)
               elseif (s.eq.4) then
                  phiA(i,k) = pfnp(i,k)
                  phiE(i,k) = pfnpE(i,k)
               elseif (s.eq.5) then
                  phiA(i,k) = eta(i)
                  phiE(i,k) = etaE(i)
               endif
               ErrorA(i,k) = abs(phiA(i,k)-phiE(i,k))
               ErrorR(i,k) = abs(phiE(i,k))
               maxErrorA = max(maxErrorA,ErrorA(i,k))
               maxErrorR = max(maxErrorR,ErrorR(i,k))
               sumErrorA = sumErrorA+ErrorA(i,k)**2
               sumErrorR = sumErrorR+ErrorR(i,k)**2
            enddo
         enddo
!        ====================================
!        =====  START PARALLEL OPTION =======
#        ifdef KeyParallel
            call MAX_parallel(maxErrorA,MAXmaxErrorA)
            call MAX_parallel(maxErrorR,MAXmaxErrorR)
            call SUM_parallel(sumErrorA,SUMsumErrorA)
            call SUM_parallel(sumErrorR,SUMsumErrorR)
            maxErrorA = MAXmaxErrorA
            maxErrorR = MAXmaxErrorR
            sumErrorA = SUMsumErrorA
            sumErrorR = SUMsumErrorR
#        endif	
!        =============== END ================    
!        ====================================
         sumErrorA = dsqrt(sumErrorA/(N_CELL0global*(NZglobal-2)))
         sumErrorR = dsqrt(sumErrorR/(N_CELL0global*(NZglobal-2)))
        
         if (maxErrorR.le.1e-8) then
            maxErrorR = 0.0d0
         else
            maxErrorR = maxErrorA/maxErrorR
         endif
         if (sumErrorR.le.1e-8) then
            sumErrorR = 0.0d0
         else
            sumErrorR = sumErrorA/sumErrorR
         endif
!        ________________________________________________________
!        Save norm errors 
         SaveErrorMax(s) = maxErrorA
         SaveErrorSum(s) = sumErrorA 
!        ________________________________________________________
!        Display solution 
         if (DisplayThis.eq.1) then
            if (s.eq.1) write(*,9) ' u ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.2) write(*,9) ' v ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.3) write(*,9) ' w ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.4) write(*,9) ' p ',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
            if (s.eq.5) write(*,9) 'eta',maxErrorA,sumErrorA,&
                                         maxErrorR,sumErrorR
         endif
      ENDDO
      if (DisplayThis.eq.1) then
      write(*,8)'===================================================='
      write(*,*)' '
      endif                        
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      5  format(t12,a4,i3,a18,e10.3)
      6  format(t14,a7,f8.4,a5,f8.4,a2)
      7  format(t25,a4,i3)
      8  format(t5,60a)
      9  format(t5,a3,t9,e10.3,t21,e10.3,t34,e10.3,t46,e10.3)

      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |                Save Errors (FS_ErrTime.dat)            |
!     |________________________________________________________|
           
      write(9100,*) time,SaveErrorMax(1),SaveErrorSum(1),&
                         SaveErrorMax(2),SaveErrorSum(2),&
                         SaveErrorMax(3),SaveErrorSum(3),&
                         SaveErrorMax(4),SaveErrorSum(4)
                         
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: TestTimeErrorFS'
         write(*,*) ''
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   END OF TEST NAVIER-STOKES EQUATION                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
