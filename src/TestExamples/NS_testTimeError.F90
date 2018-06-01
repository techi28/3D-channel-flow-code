!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        ERROR TIME EQUATION                      !
!                              Nov 2013                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                  ERROR OF THE NAVIER-STOKES EQUATION                !
!                            March 2017                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE NS_testTimeError(ufnp,ufnpE,Erruf,ufv,ufvE,Errufv,   &
                                 vfnp,vfnpE,Errvf,vfv,vfvE,Errvfv,    &
                                 wfnp,wfnpE,Errwf,wfv,wfvE,Errwfv,    &
                                 pfnp,pfnpE,Errpf,pfv,pfvE,Errpfv,    &
!                                -------------------------------------- 
                                 Hprn,etan,                           &
                                 Hprv,etav,                           &
                                 h,hv,                                &  
!                                --------------------------------------
                                 xc,yc,sig,dsig,No_cp,nbe,            &
                                 xv,yv,sigv,dsigv,No_vp,nbev)               

!---------------------------------------------------------------------!
!                                                                     !
!    This program calculates the exact solution and errors of the     !
!    time problems. An important remark is that we display more       !
!    varaibles than we actually have.                                 !
!                                                                     !
!---------------------------------------------------------------------!
!    Output  variables:                                               !
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

      real*8,dimension(:,:) :: ufnp(N_CELL,NZ)             
      real*8,dimension(:,:) :: ufnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Erruf(N_CELL,NZ) 
      real*8,dimension(:,:) :: ufv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: ufvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errufv(N_VERT,NZ-1)

      real*8,dimension(:,:) :: vfnp(N_CELL,NZ)     
      real*8,dimension(:,:) :: vfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errvf(N_CELL,NZ) 
      real*8,dimension(:,:) :: vfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: vfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errvfv(N_VERT,NZ-1)

      real*8,dimension(:,:) :: wfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errwf(N_CELL,NZ) 
      real*8,dimension(:,:) :: wfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: wfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errwfv(N_VERT,NZ-1)

      real*8,dimension(:,:) :: pfnp(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnpE(N_CELL,NZ)
      real*8,dimension(:,:) :: Errpf(N_CELL,NZ) 
      real*8,dimension(:,:) :: pfv(N_VERT,NZ-1)
      real*8,dimension(:,:) :: pfvE(N_VERT,NZ-1)
      real*8,dimension(:,:) :: Errpfv(N_VERT,NZ-1)
!     --------------------------------------
      real*8, dimension(:)  :: Hprn(N_CELL)
      real*8, dimension(:)  :: etan(N_CELL)
      real*8, dimension(:)  :: Hprv(N_VERT)
      real*8, dimension(:)  :: etav(N_VERT)
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
!     -------------------------------------
      real*8,dimension(:)   :: xv(N_VERT)
      real*8,dimension(:)   :: yv(N_VERT)
      real*8,dimension(:)   :: sigv(NZ-1)
      real*8,dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:):: No_vp(N_CELL0,3)
      integer,dimension(:)  :: nbev(N_VERT)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real*8,dimension(:,:) :: ufnpA(N_CELL,NZ)
      real*8,dimension(:,:) :: vfnpA(N_CELL,NZ)
      real*8,dimension(:,:) :: wfnpA(N_CELL,NZ)
      real*8,dimension(:,:) :: pfnpA(N_CELL,NZ)      
!     -------------------------------------      
      real*8,dimension(:) :: Err2D(N_CELL)
      real*8,dimension(:) :: Err2Dv(N_VERT)
!     -------------------------------------
      real*8,dimension(:,:) :: phiA(N_CELL,NZ)
      real*8,dimension(:,:) :: phiE(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorA(N_CELL,NZ)
      real*8,dimension(:,:) :: ErrorR(N_CELL,NZ)
      real*8,dimension(:) :: SaveErrorMaxA(4)
      real*8,dimension(:) :: SaveErrorMaxR(4)
      real*8,dimension(:) :: SaveErrorSumA(4)      
      real*8,dimension(:) :: SaveErrorSumR(4)      
      real*8 :: MaxErrorA,MaxErrorR
      real*8 :: sumErrorA,sumErrorR
      real*8 :: MAXmaxErrorA,MAXmaxErrorR
      real*8 :: SUMsumErrorA,SUMsumErrorR
      real*8 :: Peakvalue,PeakvalueEx,PeakError
!     --------------------------------------
      real*8 :: x,y,z
      real*8 :: TimeExample2D,TimeExample3D     
      real*8 :: funExamNSu,funExamNSv,funExamNSw
      real*8 :: funExamNSp,funExamNSrhsp,cons,cons1,cons2
      integer:: kmid1,imid1,kmid2,imid2
      integer:: DisplayThis,s
!     --------------------------------------

!*********************************************************************!
!                                                                     !
!                             Initialization                          !
!                                                                     !
!*********************************************************************!


!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '>>>>> Begin subroutine: NS_TestTimeError'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      
!*********************************************************************!
!                                                                     !
!                        Navier-Stokes solution                       !
!                                                                     !
!*********************************************************************!
!      ________________________________________________________
!     |                                                        |
!     |            Exact solution & Absolute Errors            |
!     |________________________________________________________|

      ufnpA = ufnp
      vfnpA = vfnp
      wfnpA = wfnp            
      pfnpA = pfnp
      
!     _________________________________________________________
!     Difference constant between analytical and numerical
#     ifdef Key_NeumannBCp
         kmid1 = floor(NZ/2.0d0)
         imid1 = floor(N_CELL0/2.0d0)
         cons1 = pfnp(imid1,kmid1)-pfnpE(imid1,kmid1)
         !kmid2 = floor(NZ/2.0d0)
         !imid2 = floor(N_CELL0/4.0d0)
         !cons2 = pfnp(imid2,kmid2)-pfnpE(imid2,kmid2)
         cons  = cons1!max(cons1,cons2)
         !print*,cons1,cons2         
         do k=1,NZ 
            do i=1,N_CELL
               pfnpA(i,k) = pfnp(i,k) - cons
            enddo
         enddo      
#     endif

!     _________________________________________________________
!     Cell-center
      do k=1,NZ 
         do i=1,N_CELL
            x = xc(i)
            y = yc(i)
            z = sig(k)*Hprn(i)-h(i)
!           -----------------------------------
!           Exact values 
            ufnpE(i,k) = funExamNSu(x,y,z,time)
            vfnpE(i,k) = funExamNSv(x,y,z,time)
            wfnpE(i,k) = funExamNSw(x,y,z,time) 
            pfnpE(i,k) = funExamNSp(x,y,z,time) 
!           -----------------------------------
!           Error 
            Erruf(i,k) = abs(ufnpA(i,k)-ufnpE(i,k))
            Errvf(i,k) = abs(vfnpA(i,k)-vfnpE(i,k))
            Errwf(i,k) = abs(wfnpA(i,k)-wfnpE(i,k))
            Errpf(i,k) = abs(pfnpA(i,k)-pfnpE(i,k))
         enddo
      enddo 

!     _________________________________________________________
!     Vertices
  
      do k=1,NZ-1  
         do nv=1,N_VERT
            x = xv(nv)
            y = yv(nv)
            z = sigv(k)*Hprv(nv)-hv(nv)
            ufvE(nv,k) = funExamNSu(x,y,z,time)
            vfvE(nv,k) = funExamNSv(x,y,z,time)
            wfvE(nv,k) = funExamNSw(x,y,z,time) 
            pfvE(nv,k) = funExamNSp(x,y,z,time)
            !---------
!           Error 
            Errufv(nv,k) = abs(ufv(nv,k)-ufvE(nv,k))
            Errvfv(nv,k) = abs(vfv(nv,k)-vfvE(nv,k))
            Errwfv(nv,k) = abs(wfv(nv,k)-wfvE(nv,k))
            Errpfv(nv,k) = abs(pfv(nv,k)-pfvE(nv,k))    
         enddo
      enddo     
!      ________________________________________________________
!     |                                                        |
!     |                      Norm Errors                       |
!     |________________________________________________________|

      DisplayThis = 1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
       if (rang_topo.ne.0) then
          DisplayThis = 0
       endif
#     endif
!     =============== END ================    
!     ====================================

      IF (1d-08.ge.dtsave-(time-LastTimeSave)) THEN
      if (DisplayThis.eq.1) then
      write(*,8)'===================================================='
      write(*,8)'                                                    '
      write(*,8)'           TEST: Navier-Stokes problem              '
      write(*,6)' Time = ',time,', (dt= ',dt,')'
      write(*,7)' N = ',NN
      write(*,8)'____________________________________________________'
      write(*,8)'                                                    '
      write(*,8)'          Absolute error          Relative error    '
      write(*,8)'     Max norm    L-2 norm     Max norm    L-2 norm  '
      write(*,8)'----------------------------------------------------'
      endif
      ENDIF

      DO s=1,4
!        ________________________________________________________
!        Norm errors 
         maxErrorA = 0.0d0
         maxErrorR = 0.0d0
         sumErrorA = 0.0d0
         sumErrorR = 0.0d0
         do k=2,NZ-1 
            do i=1,N_CELL0 
               if (s.eq.1) then
                  phiA(i,k) = ufnpA(i,k)
                  phiE(i,k) = ufnpE(i,k)
               elseif (s.eq.2) then
                  phiA(i,k) = vfnpA(i,k)
                  phiE(i,k) = vfnpE(i,k)
               elseif (s.eq.3) then
                  phiA(i,k) = wfnpA(i,k)
                  phiE(i,k) = wfnpE(i,k)
               elseif (s.eq.4) then
                  phiA(i,k) = pfnpA(i,k)
                  phiE(i,k) = pfnpE(i,k)
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
         sumErrorR = sumErrorA/sumErrorR
         maxErrorR = maxErrorA/maxErrorR
!        ________________________________________________________
!        Save norm errors 
         SaveErrorMaxA(s) = maxErrorA
         SaveErrorMaxR(s) = maxErrorR         
         SaveErrorSumA(s) = sumErrorA
         SaveErrorSumR(s) = sumErrorR 
      ENDDO              
!      ________________________________________________________
!     |                                                        |
!     |                  Display Norm Errors                   |
!     |________________________________________________________|

      IF (1d-08.ge.dtsave-(time-LastTimeSave)) THEN
      if (DisplayThis.eq.1) then                          
            write(*,9) ' u ',SaveErrorMaxA(1),SaveErrorSumA(1),&         
                             SaveErrorMaxR(1),SaveErrorSumR(1)
            write(*,9) ' v ',SaveErrorMaxA(2),SaveErrorSumA(2),&         
                             SaveErrorMaxR(2),SaveErrorSumR(2)
            write(*,9) ' w ',SaveErrorMaxA(3),SaveErrorSumA(3),&         
                             SaveErrorMaxR(3),SaveErrorSumR(3)
            write(*,9) ' p ',SaveErrorMaxA(4),SaveErrorSumA(4),&         
                             SaveErrorMaxR(4),SaveErrorSumR(4)
      write(*,8)'===================================================='
      write(*,*)' '
      endif
      ENDIF

!      ________________________________________________________
!     |                                                        |
!     |             Save Norm Errors (Matlab view)             |
!     |________________________________________________________|

      write(8100,10) time,SaveErrorSumA(1),SaveErrorSumA(2),&
                          SaveErrorSumA(3),SaveErrorSumA(4),&
                          SaveErrorSumR(1),SaveErrorSumR(2),&
                          SaveErrorSumR(3),SaveErrorSumR(4)
                          
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!

      6  format(t16,a7,f8.4,a6,f7.5,a1)
      7  format(t25,a4,i3)
      8  format(t5,60a)
      9  format(t5,a3,t9,e10.3,t21,e10.3,t34,e10.3,t46,e10.3)
      10 format(9(1x,e12.5))
      11 format(f7.5,4(1x,e12.5))    

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t15,60a)'), '<<<<< End   subroutine: NS_TestTimeError'
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
