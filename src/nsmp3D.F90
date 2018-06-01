!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	    MAIN PROGRAM                              !
!       Solution of the Navier-Stokes Multi-phase equations           !
!                      Miguel Angel Uh Zapata                         !
!                  Last modification: Nov 2017                        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!
      PROGRAM nsmp3D
!
!---------------------------------------------------------------------!
!                                                                     !
!     This program find the solution of the Navier-Stokes equations   !
!     corresponding to the free surface two-phase flow problem in 3D. !
!     In this program all the defined and common variables are modi-  !
!     fied.                                                           !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !
!  |_____________|_________________________________________________|  !
!  | nstep       | Current number of time step                     |  !
!  | tcpu        | cpu time of the simulation                      |  !
!  | tt(1:2)     | user and system simulation time                 |  !
!  | nt          | Approx number of time steps to display time     |  !
!  | time0       | Initial time of the current simulation          |  !
!  | usepas      | Elapsed time for one point & one time step      |  !
!  | idhr        | Number of hours of the simulation               |  !
!  | idmin       | Number of minutes of the simulation             |  !
!  | idsec       | Number of seconds of the simulation             |  !
!  |_____________|_________________________________________________|  !
!                                                                     !
!    Subroutines used:                                                !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |  STRUCTURAL:                                                  |  !
!  |                      ° common.mpf                             |  !
!  |                      ° cppdefs.h                              |  !
!  |                      ° definition_variables                   |  !
!  |                      ° alloc_variables                        |  !
!  |                      ° dealloc_variables                      |  !
!  |                      ° interfaces.F90                         |  !
!  |                                                               |  !
!  |  INITIALIZATION:                                              |  !
!  |                      ° input_data.F90                         |  !
!  |                      ° input_parameters.F90                   |  !
!  |                      ° input_initial.F90                      |  !
!  |                      ° geometry.F90                           |  !
!  |                                                               |  !
!  |  TIME LOOP UPDATES:                                           |  !
!  |                      ° hydro.F90                              |  !
!  |                                                               |  !
!  |  RE-START & SAVING:                                           |  !
!  |                      ° restart_out.F90                        |  !
!  |                      ° restart_in.F90                         |  !
!  |                                                               |  !
!  |_______________________________________________________________|  !
!                                                                     !
!---------------------------------------------------------------------!
!*********************************************************************!
!                                                                     !
!                            Definitions                              !
!                                                                     !
!*********************************************************************!

!     ____________________________________
!    |                                    |
!    |     Keys and common parameters     |
!    |____________________________________|

#     include "cppdefs.h"

!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============
#     ifdef KeyCUDA
          USE parallelCUDA
#     endif
!     =============== END ================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
          USE parallel
#     endif
!     =============== END ================
!     ====================================

      USE variables
      USE geometry
      USE interfaces
      USE statistics
      implicit none

!     ____________________________________
!    |                                    |
!    |      Declaration of variables      |
!    |____________________________________|

      real,dimension(2) :: tt
      real   :: tcpu,MAXtcpu,start_nsmp,finish_nsmp
      real*8 :: nt,usepas
      real*8 :: time0
      integer:: idmin,idhr,idsec
      integer:: m,SaveFiles
      integer:: nRK
      character*80 title

!*********************************************************************!
!                                                                     !
!                         I) Initialization                           !
!                                                                     !
!*********************************************************************!

!      _____________________________________________________________
!     |      |                                                      |
!     | I.1  |          Start time simulation & display             |
!     |______|______________________________________________________|

      call cpu_time(start_nsmp)

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      print*,'                                                            '
      print*,'                                                            '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'------------------------------------------------------------'
      print*,'            __________________________________              '
      print*,'           |                                  |             '
      print*,'           |    ***  PROGRAM: NSMP3D   ***    |             '
      print*,'           |__________________________________|             '
      print*,'                                                            '
      print*,'------------------------------------------------------------'
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                            '
      print*,'                                                            '
      print*,'    ___________________________________________________     '
      print*,'   |                                                   |    '
      print*,'   |                  INITIALIZATION                   |    '
      print*,'   |___________________________________________________|    '
      print*,'                                                            '
#     endif
!     =============== END ================
!     ====================================

!     ____________________________________
!     WARNING
#     ifndef KeyFixedFreeSurface
#     ifdef KeyMaximumEta
         print*,'    WARNING-REMARK: We are using a numerical trick '
         print*,'                --> Setting a maximum water level'
#     endif
#     endif

!      _____________________________________________________________
!     |      |                                                      |
!     | I.2  |     Read: Number of cell-center and vertex points    |
!     |______|______________________________________________________|

!     ________________________________
#     if defined(KeyEstuaryGironde)
         open(25,file='EstuaryGironde/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStaticCylinder)
         open(25,file='StaticCylinder/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStaticChannel)
         open(25,file='StaticChannel/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyStandingWave)
         open(25,file='StandingWave/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyTaylorVortex)
         open(25,file='TaylorVortex/data.txt',status='OLD')
!     ________________________________
#     elif defined(KeyTestOnlyPoisson)
         open(25,file='TestOnlyPoisson/data.txt',status='OLD')
!     ________________________________
#     else
         open(22,file='data.txt',status='OLD')
#     endif

      read(25,*) title
      read(25,*) N_VERTglobal
      read(25,*) N_CELL0global

      close(25)

!     _____________________________________________________________
!    |      |                                                      |
!    | I.3  |               Defition of parameters                 |
!    |______|______________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         NZ      = NZglobal
         N_VERT  = N_VERTglobal
         N_CELL0 = N_CELL0global
         N_CELL  = N_CELL0 + N_CELLghostMAX
         print*,'                                                       '
         print*,'       DOMAIN:                                         '
         print*,'          Number of elements  =',N_CELL0Global
         print*,'          Number of vertices  =',N_VERTGlobal
         print*,'          Number of NZ points =',NZglobal-1
         print*,'                                                       '
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!                         P) Parallelization                          !
!                                                                     !
!*********************************************************************!

!     _____________________________________________________________
!    |      |                                                      |
!    | P.1  |            Parallel structure MPI                    |
!    |______|______________________________________________________|

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
!       __________________________________
!       Distribution
        call initialisation_mpi
        call parallel_input_global
        call parallel_distribution
!       __________________________________
!       Topology
        call parallel_topology
        call parallel_neighbors
!       __________________________________
!       Defition of parameters
        call parallel_parameters
!       __________________________________
!       Parallel utilities
        call parallel_index
        call parallel_shareindex
        call parallel_type
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!           1)  Allocate variables, input & output data files         !
!                                                                     !
!*********************************************************************!

!      _____________________________________________________________
!     |      |                                                      |
!     | 1.1  |                Allocate variables                    |
!     |______|______________________________________________________|

      call alloc_variables
      call alloc_geometry
!     -----------------------
#     ifdef KeySaveStatistics
      call alloc_stats_variables
#     endif

!      _____________________________________________________________
!     |      |                                                      |
!     | 1.2  |          Input domain data (read data.txt)           |
!     |______|______________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
        call input_data(No_vp,No_cp,No_wb,No_hb,No_qb,No_sp,&
                        nbe,xv,yv,zbv)
!        __________________________
!        Index re-ordering
#        ifdef KeyOrderingIndex
            call ReOrderingIndex(No_vp,No_cp,nbe)
#        endif
!        __________________________
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        call parallel_input_local(No_vp,No_cp,No_wb,No_qb,No_hb,No_sp,&
                                  nbe,xv,yv,zbv)
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.3  |       Input parameters (read dataParameters.txt)       |
!     |______|________________________________________________________|

      call input_parameters(nbev,hv,h,                     &
                            No_vp,No_cp,No_wb,No_hb,No_qb, &
                            nbe,xv,yv,zbv)

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.4  |      Open global output files: time, error, iters      |
!     |______|________________________________________________________|

      SaveFiles = 0
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         SaveFiles = 1
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         if (rang_topo.eq.0) then
            SaveFiles = 1
         endif
#     endif
!     =============== END ================
!     ====================================

      IF (SaveFiles.eq.1) THEN
         print*,'  '
         print*,' wwwwwwwwwwwwwwwwww  OUTPUT FILE  wwwwwwwwwwwwwwwwwwww'
         print*,'  '
         print*,' Saved: ../output/Matlab/PlotConvergence/ItersTime.dat'
         open(7100,file='../output/Matlab/PlotConvergence/ItersTime.dat')
         !---------------------------
#        ifdef KeyTaylorVortex
         print*,' Saved: ../output/NS/ErrorTime.dat'
         open(8100,file="../output/NS/ErrorTime.dat",status='unknown')
#        endif
         !---------------------------
#        ifdef KeyStandingWave
         print*,' Saved: ../output/FS/FS_ErrTime.dat'
         open(9100,file="../output/FS/FS_ErrTime.dat",status='unknown')
#        endif
         print*,'  '
         print*,' wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
         print*,'  '
      ENDIF

!     ------------------------------------
!     For the new version (Dic 2017)
!#     ifndef KeyFixedFreeSurface
!#        ifdef KeySave1DReference
!            call FS_SaveReferenceOpen
!#        endif
!#     endif

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.4  |            Coloring (used for MultiSOR)                |
!     |______|________________________________________________________|

#     ifdef KeyMultiSOR
         call coloring(No_cp,No_vp,xv,yv)
#     endif
#     ifdef KeyAuxMCSOR
         call coloring(No_cp,No_vp,xv,yv)
#     endif

!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============
#     ifdef KeyCUDA
#        ifdef KeyMSOR_cuda
            call coloring(No_cp,No_vp,xv,yv)
#        endif
#        ifdef KeyPDMSOR_cuda
            call coloring(No_cp,No_vp,xv,yv)
#        endif
#     endif
!     =====================================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

!     _____________________________________________________________
!    |      |                                                      |
!    | P.2  |            Parallel structure CUDA                   |
!    |______|______________________________________________________|

!     REMARK: It needs coloring before to apply the fastest version

!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============
#     ifdef KeyCUDA
         print*,'    ___________________________________________________  '
         print*,'   |                                                   | '
         print*,'   |   * * *  ======> PARALLEL CUDA  <=======  * * *   | '
         print*,'   |___________________________________________________| '
         print*,'                                                         '
         call initialisation_cuda
         call allocate_cuda
#     endif
!     =============== END ================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

!*********************************************************************!
!                                                                     !
!                  2) Geometry & Initial conditions                   !
!                                                                     !
!*********************************************************************!

!      _______________________________________________________________
!     |      |                                                        |
!     | 2.1  |            Geometry variables of the mesh              |
!     |______|________________________________________________________|

      call calcul_geometry(xc,yc,sig,dsig,No_cp,nbe,h,  &
                           xv,yv,sigv,dsigv,No_vp,nbev, &
                           ic1tec,ic2tec,ic3tec)

!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============
#     ifdef KeyCUDA
      call transfergeometry_cuda(No_cp,No_vp,nbev,sig,sigv)
      call transferInterpo_cuda
#     endif
!     =====================================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

!     _______________________________________________________________
!     Save full matrix structure
#     ifdef KeySaveMatrixA
         call SaveFullMatrixA(No_cp,nbe,No_vp,nbev)
#     endif
!      _______________________________________________________________
!     |      |                                                        |
!     | 2.2  |              Input initial conditions                  |
!     |______|________________________________________________________|

!     ___________________________________________
!     Initial time, time step & save counter
      time  = 0.0d0
      nstep = 0
      SaveCounter = 0
      scount = 0

!     ___________________________________________
!     Initial values of the main variables
      call input_initial(alphafn,ufn,vfn,wfn,pfn,viscof,rhof,   &
                         alphasn,usn,vsn,wsn,psn,viscos,rhos,   &
                         alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv, &
                         alphasv,usv,vsv,wsv,psv,viscosv,rhosv, &
                         etan,etav,                             &
                         Hprn,Hprv,                             &
                         h,hv,                                  &
                         xc,yc,sig,dsig,No_cp,nbe,              &
                         xv,yv,sigv,dsigv,No_vp,nbev,           &
                         Heaviside,mask)

!     ___________________________________________
!     Vorticity for initial values
      call Vorticity(ufn,vfn,wfn,               &
                     xc,yc,sig,dsig,No_cp,nbe,  &
                     xv,yv,sigv,dsigv,No_vp,nbev)

!     ___________________________________________
!     Initial time averaged stats
#     ifdef KeySaveStatistics
         call sm_init
         !call restart_in_stats(sm_uf,sm_vf,sm_wf,sm_pf)
#     endif

!      _______________________________________________________________
!     |      |                                                        |
!     | 2.3  |            Exact solutions & errors                    |
!     |______|________________________________________________________|

!     ________________________________________________________
!     Navier-Stokes error: Taylor vortex
#     ifdef KeyTaylorVortex
         call NS_testTimeError(ufn,usn,uErr,ufv,usv,uErrv,     &
                               vfn,vsn,vErr,vfv,vsv,vErrv,     &
                               wfn,wsn,wErr,wfv,wsv,wErrv,     &
                               pfn,psn,pErr,pfv,psv,pErrv,     &
!                              --------------------------------
                               Hprn,etan,                      &
                               Hprv,etav,                      &
!                              --------------------------------
                               h,hv,                           &
!                              --------------------------------
                               xc,yc,sig,dsig,No_cp,nbe,       &
                               xv,yv,sigv,dsigv,No_vp,nbev)
#     endif

!     ________________________________________________________
!     Free Surface error: Standing wave problem
#     ifdef KeyStandingWave
         call FS_TimeError(Hprn,etan,Hprv,etav,             &
!                          --------------------------------
                           ufnp,vfnp,wfnp,pfnp,             &
                           ufv,vfv,wfv,pfv,                 &
!                          --------------------------------
                           h,xc,yc,sig,dsig,No_cp,nbe,      &
                           hv,xv,yv,sigv,dsigv,No_vp,nbev,  &
!                          --------------------------------
                           No_sp)
#     endif
!      _______________________________________
!     |                                       |
!     |       Coordinates (xt,yt,zt)          |
!     |_______________________________________|

!     -------------
!     Cell-Center
      do i=1,N_CELL
         do k=1,NZ
            xct(i,k) = xc(i)
            yct(i,k) = yc(i)
            zct(i,k) = sig(k)*Hpr(i)-h(i)
         enddo
      enddo
!     -------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            xvt(nv,k) = xv(nv)
            yvt(nv,k) = yv(nv)
            zvt(nv,k) = sigv(k)*Hprv(nv)-hv(nv)
         enddo
      enddo

!      _______________________________________________________________
!     |      |                                                        |
!     | 2.4  |               Re-start solution (n)                    |
!     |______|________________________________________________________|

!     ________________________________________________________
!     Re-start time

      IF (IrestartIN.eq.1) THEN
!        ---------------------------------------------------
!        Read main variables
         call restart_in(alphafn,ufn,vfn,wfn,pfn,  &
                         alphasn,usn,vsn,wsn,psn,  &
                         zct)

!        ---------------------------------------------------
!        Update the free surface
      print*,'free surface restart not yet'
!         do i=1,N_CELL
!            etan(i) = zct(i,NZ)
!            Hprn(i) = etan(i) + h(i)
!         enddo
!      ________________________________________________________
!     |                                                        |
!     |                    Vertex values                       |
!     |________________________________________________________|
!      not yet BC
      call interpolation3D(ufv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           ufn,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(vfv,xv,yv,sigv,dsigv,No_vp,nbev, &
                          vfn,xc,yc,sig,dsig,No_cp,nbe)
      call interpolation3D(wfv,xv,yv,sigv,dsigv,No_vp,nbev, &
                           wfn,xc,yc,sig,dsig,No_cp,nbe)
!        ---------------------------------------------------
!        Display re-start initial time
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
         print*,'    ==================================================  '
         print*,'                  Re-start simulation                   '
         print*,'    ==================================================  '
         print*,'                                                        '
         write(*,'(t17,a13,f12.5)') '       time = ', time
         write(*,'(t17,a13,f12.5)') '         dt = ', dt
         write(*,'(t17,a13,i8)')    'SaveCounter = ', SaveCounter
         print*,'                                                        '
         print*,'    ==================================================  '
         print*,'                                                        '
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
         IF (rang_topo.eq.0) THEN
         print*,'    ==================================================  '
         print*,'                  Re-start simulation                   '
         print*,'            WARNING!! Parallel not done yet             '
         print*,'    ==================================================  '
         print*,'                                                        '
         write(*,'(t17,a13,f12.5)') '       time = ', time
         write(*,'(t17,a13,f12.5)') '         dt = ', dt
         write(*,'(t17,a13,i8)')    'SaveCounter = ', SaveCounter
         print*,'                                                        '
         print*,'    ==================================================  '
         print*,'                                                        '
         ENDIF
#        endif
!        =============== END ================
!        ====================================
      ENDIF

!      _______________________________________________________________
!     |      |                                                        |
!     | 2.5  |             Save initial condition values              |
!     |______|________________________________________________________|

!      _______________________________________
!     |                                       |
!     |        Save at reference point        |
!     |_______________________________________|

#     ifndef KeyFixedFreeSurface
#     ifdef KeySave1DReference
         call FS_SaveReference(Hprv,etav,ufv,vfv,wfv,pfv,       &
                               HprvA,etavA,ufvA,vfvA,wfvA,pfvA, &
                               h,xc,yc,sig,dsig,No_cp,nbe,      &
                               hv,xv,yv,sigv,dsigv,No_vp,nbev,  &
                               No_sp,nstep)
#     endif
#     endif
!      _______________________________________
!     |                                       |
!     |         Save Paraview/Tecplot         |
!     |_______________________________________|

      uErrv = 0
      vErrv = 0
      wErrv = 0
      pErrv = 0

!     -------------------------------------------------------
!     Last time when the results were saved
      LastTimeSave = time

!     -------------------------------------------------------
!     Save Tecplot file cell centers
      if ((ChooseExit.eq.1).or.(ChooseExit.eq.2)) then
         call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,       &
                            alphasv,usv,vsv,wsv,psv,rhosv, &
                            xvt,yvt,zvt,No_vp)
      endif
!     -------------------------------------------------------
!     Save Paraview file at vertex points
      if (ChooseExit.eq.3) then
!        ---------------------
#        ifdef KeyStandingWave
            call FS_SaveParaview(Hprv,etav,               &
                                 ufv,vfv,wfv,pfv,         &
                                 HprvA,etavA,             &
                                 ufvA,vfvA,wfvA,pfvA,     &
                                 uErrv,vErrv,wErrv,pErrv, &
                                 xvt,yvt,zvt,             &
                                 No_vp)
!        ---------------------
#        else
            call SaveParaviewVertex(Hprv,etav,                &
                                    ufv,vfv,wfv,pfv,          &
                                    usv,vsv,wsv,psv,          &
                                    uErrv,vErrv,wErrv,pErrv,  &
                                    xvt,yvt,zvt,              &
                                    No_vp)
#        endif
!        ---------------------
      endif

!*********************************************************************!
!                                                                     !
!              3) Main structure: time step simulations               !
!                                                                     !
!*********************************************************************!

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         print*,'                                                        '
         print*,'    ___________________________________________________ '
         print*,'   |                                                   |'
         print*,'   |                TIME SIMULATIONS                   |'
         print*,'   |___________________________________________________|'
         print*,'                                                        '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.0) THEN
         print*,'                                                        '
         print*,'   ==================================================== '
         print*,'                    MPI: TIME SIMULATIONS               '
         print*,'   ==================================================== '
         print*,'                                                        '
         ENDIF
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.1  |          Time loop: Initial values: (n+1)=(n)          |
!     |______|________________________________________________________|

      Hprnp = Hprn
      etanp = etan
!     -------
      alphafnp = alphafn
      alphasnp = alphasn
!     -------
      ufnp = ufn
      vfnp = vfn
      wfnp = wfn
      pfnp = pfn
!     -------
      usnp = usn
      vsnp = vsn
      wsnp = wsn
      psnp = psn

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.2  |         Time loop:  Beggining of the new time step     |
!     |______|________________________________________________________|

10    continue

      time  = time + dt
      nstep = nstep + 1

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifndef KeyDisplay
!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         write(*,'(a8,i8)')  'Step : ',nstep
         write(*,'(a8,f10.5)')'Time = ',time
         print*,'  '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.0) THEN
         write(*,'(a8,i8)')  'Step : ',nstep
         write(*,'(a8,f10.5)')'Time = ',time
         print*,'  '
         ENDIF
#     endif
!     =============== END ================
!     ====================================
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!      ________________________________________________________
!     |      |                                                 |
!     | 3.2.1|  Runge-Kutta: Initial known values: (known)=(n) |
!     |______|_________________________________________________|

      eta = etan
      Hpr = Hprn
!     -------
      alphaf = alphafn
      alphas = alphasn
!     -------
      uf = ufn
      vf = vfn
      wf = wfn
      pf = pfn
!     -------
      us = usn
      vs = vsn
      ws = wsn
      ps = psn

!      ________________________________________________________
!     |      |                                                 |
!     | 3.2.2|  Runge-Kutta: The two step loop                 |
!     |______|_________________________________________________|

      DO nRK=1,FinRK

!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         IF (FinRK==2) THEN
!        ====================================
!        ==========  SEQUENTIAL =============
#        ifndef KeyParallel
            print*,' '
            write(*,'(t10,a8,i3)')'RK step:',nRK
!        ====================================
!        =====  START PARALLEL OPTION =======
#        else
            IF (rang_topo.eq.0) THEN
            print*,' '
            write(*,'(t10,a8,i3)')'RK step:',nRK
            ENDIF
#        endif
!        =============== END ================
!        ====================================
         ENDIF
!        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!         _______________________________________
!        |                                       |
!        | ************************************* |
!        |            Update variables           |
!        |    Navier-Stokes equations with FS    |
!        | ************************************* |
!        |_______________________________________|

         call Hydro(ufnp,vfnp,wfnp,pfnp,             &
                    uf,vf,wf,pf,                     &
                    ufv,vfv,wfv,pfv,                 &
!                   ----------------------------------
                    Hprnp,etanp,                     &
                    Hpr,eta,                         &
                    Hprv,etav,                       &
!                   ----------------------------------
                    h,hv,                            &
!                   ----------------------------------
                    xc,yc,sig,dsig,No_cp,nbe,        &
                    xv,yv,sigv,dsigv,No_vp,nbev,     &
!                   ----------------------------------
                    No_wb,No_qb,No_hb,No_sp,         &
!                   ----------------------------------
                    nRK)

!         _______________________________________
!        |                                       |
!        |   Update RK2 values:(known)=(n+1)     |
!        |_______________________________________|

         eta = etanp
         Hpr = Hprnp
!        -------
         alphaf = alphafnp
         alphas = alphasnp
!        -------
         uf = ufnp
         vf = vfnp
         wf = wfnp
         pf = pfnp
!        -------
         us = usnp
         vs = vsnp
         ws = wsnp
         ps = psnp

      ENDDO

!      ________________________________________________________
!     |      |                                                 |
!     | 3.2.3|  Runge-Kutta: Final solution of the RK2 formula |
!     |______|_________________________________________________|

      IF (FinRK==2) THEN
         etanp = 0.5d0*(etan+etanp)
         Hprnp = 0.5d0*(Hprn+Hprnp)
!        -------
         alphafnp = 0.5d0*(alphafn+alphafnp)
         alphasnp = 0.5d0*(alphasn+alphasnp)
!        -------
         ufnp = 0.5d0*(ufn+ufnp)
         vfnp = 0.5d0*(vfn+vfnp)
         wfnp = 0.5d0*(wfn+wfnp)
         pfnp = 0.5d0*(pfn+pfnp)
!        -------
         usnp = 0.5d0*(usn+usnp)
         vsnp = 0.5d0*(vsn+vsnp)
         wsnp = 0.5d0*(wsn+wsnp)
         psnp = 0.5d0*(psn+psnp)
      ENDIF

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.3  |     Time loop: Update the variables & input box        |
!     |______|________________________________________________________|

      Hprn=Hprnp
      etan=etanp
!     -------
      alphafn = alphafnp
      alphasn = alphasnp
!     -------
      ufn = ufnp
      vfn = vfnp
      wfn = wfnp
      pfn = pfnp
!     -------
      usn = usnp
      vsn = vsnp
      wsn = wsnp
      psn = psnp

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.4  |     Time loop:  Errors during simulation               |
!     |______|________________________________________________________|

!     ________________________________________________________
!     Navier-Stokes error: Taylor vortex
#     ifdef KeyTaylorVortex
         call NS_testTimeError(ufnp,usnp,uErr,ufv,usv,uErrv,     &
                               vfnp,vsnp,vErr,vfv,vsv,vErrv,     &
                               wfnp,wsnp,wErr,wfv,wsv,wErrv,     &
                               pfnp,psnp,pErr,pfv,psv,pErrv,     &
!                              -----------------------------------
                               Hpr,etan,                         &
                               Hprv,etav,                        &
!                              -----------------------------------
                               h,hv,                             &
!                              -----------------------------------
                               xc,yc,sig,dsig,No_cp,nbe,         &
                               xv,yv,sigv,dsigv,No_vp,nbev)
#     endif

!     ________________________________________________________
!     Free surface error: standing wave
#     ifdef KeyStandingWave
         call FS_TimeError(Hpr,eta,Hprv,etav,               &
                           ufnp,vfnp,wfnp,pfnp,             &
                           ufv,vfv,wfv,pfv,                 &
                           h,xc,yc,sig,dsig,No_cp,nbe,      &
                           hv,xv,yv,sigv,dsigv,No_vp,nbev,  &
                           No_sp)
#     endif

!     ________________________________________________________
!     Test only Poisson equation (exact=psv)
#     ifdef KeyTestOnlyPoisson
          call Poisson_testError(pfnp,pfv,                      &
                                 pfv,psv,pErrv,                 &
                                 xc,yc,sig,dsig,No_cp,nbe,      &
                                 xv,yv,sigv,dsigv,No_vp,nbev,   &
                                 Hpr,h,eta,                     &
                                 Hprv,hv,etav)
          goto 9999
#     endif

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.5  |    Time loop:  Saving results during simulation        |
!     |______|________________________________________________________|

!      _______________________________________
!     |                                       |
!     |       Coordinates (xt,yt,zt)          |
!     |_______________________________________|

!     -------------
!     Cell-Center
      do i=1,N_CELL
         do k=1,NZ
            xct(i,k) = xc(i)
            yct(i,k) = yc(i)
            zct(i,k) = sig(k)*Hpr(i)-h(i)
         enddo
      enddo
!     -------------
!     Vertex
      do nv=1,N_VERT
         do k=1,NZ-1
            xvt(nv,k) = xv(nv)
            yvt(nv,k) = yv(nv)
            zvt(nv,k) = sigv(k)*Hprv(nv)-hv(nv)
         enddo
      enddo

!      _______________________________________
!     |                                       |
!     |            Save Statistics            |
!     |_______________________________________|

#     ifdef KeySaveStatistics
!     ________________________________________
!     Obtain time average data for cylinder case
      !if (mod(nstep,dtsample).eq.0) then
      !   call update_stats(sm_uf,sm_vf,sm_wf,sm_pf, &
      !                     ufnp,vfnp,wfnp,pfnp)
      !endif
!     ________________________________________
!     Stats for Open Channel Case
      IF (time.gt.tInistats) THEN
         if (mod(nstep,dtsample).eq.0) then
            call glob_sample(ufnp,vfnp,wfnp,pfnp)
         endif
         if (mod(nstep,dtplane).eq.0) then
            call glob_put_plane_mm
         endif
      ENDIF
#     endif

!      _______________________________________
!     |                                       |
!     |        Save at reference point        |
!     |_______________________________________|

#     ifndef KeyFixedFreeSurface
#     ifdef KeySave1DReference
         call FS_SaveReference(Hprv,etav,ufv,vfv,wfv,pfv,       &
                               HprvA,etavA,ufvA,vfvA,wfvA,pfvA, &
                               h,xc,yc,sig,dsig,No_cp,nbe,      &
                               hv,xv,yv,sigv,dsigv,No_vp,nbev,  &
                               No_sp,nstep)
#     endif
#     endif

!      _______________________________________
!     |                                       |
!     |         Save Paraview/Tecplot         |
!     |_______________________________________|

      IF ((time.ge.tInisave).and.&
         (1d-08.ge.dtsave-(time-LastTimeSave))) THEN
!        ---------------------------------------------------
!        Last time saved
         SaveCounter  = SaveCounter+1
         LastTimeSave = SaveCounter*dtsave
!        ---------------------------------------------------
!        Save Tecplot file vertices
         if ((ChooseExit.eq.1).or.(ChooseExit.eq.2)) then
            call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,       &
                               alphasv,usv,vsv,wsv,psv,rhosv, &
                               xvt,yvt,zvt,No_vp)
         endif
!        ---------------------------------------------------
!        Save Paraview file vertices
         if (ChooseExit.eq.3) then
!           ---------------------
#           ifdef KeyStandingWave
            call FS_SaveParaview(Hprv,etav,                &
                                 ufv,vfv,wfv,pfv,          &
                                 HprvA,etavA,              &
                                 ufvA,vfvA,wfvA,pfvA,      &
                                 uErrv,vErrv,wErrv,pErrv,  &
                                 xvt,yvt,zvt,              &
                                 No_vp)
!           ---------------------
#           else
            call SaveParaviewVertex(Hprv,etav,                &
                                    ufv,vfv,wfv,pfv,          &
                                    usv,vsv,wsv,psv,          &
                                    uErrv,vErrv,wErrv,pErrv,  &
                                    xvt,yvt,zvt,              &
                                    No_vp)
#           endif
!           ---------------------
         endif
!        ---------------------------------------------------
!        Save Tecplot file for stats
#        ifdef KeySaveStatistics
            !call SavetecStats(sm_uf,sm_vf,sm_wf,sm_pf,xvt,yvt,zvt,No_vp)
            !call restart_out_stats(sm_uf,sm_vf,sm_wf,sm_pf)
#        endif
!        ---------------------------------------------------
!        Save Re-start data
         if (IrestartOUT.eq.1) then
            call restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp,    &
                             alphasnp,usnp,vsnp,wsnp,psnp,    &
                             zct)
         endif
      ENDIF

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.6  |    Time loop:  Criteria to finish time simulations     |
!     |______|________________________________________________________|

      if ((tfin-time).le.1d-08) goto 9999
      goto 10

9999  continue


!*********************************************************************!
!                                                                     !
!                          F) Finalization                            !
!                                                                     !
!*********************************************************************!

!      _______________________________________________________________
!     |      |                                                        |
!     | F.1  |                Saving final results                    |
!     |______|________________________________________________________|

      if ((tfin-LastTimeSave).ge.1d-8) then
!        ---------------------------------------------------
!        Last time saved
         SaveCounter  = SaveCounter+1
!        ---------------------------------------------------
!        Save Tecplot file vertices
         if ((ChooseExit.eq.1).or.(ChooseExit.eq.2)) then
            print*,'       --- Tecplot final time solution printed ---'
            call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,           &
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
         endif
!        ---------------------------------------------------
!        Save Paraview file vertices
         if (ChooseExit.eq.3) then
            print*,'       --- Paraview final time solution printed ---'
!           ---------------------
#           ifdef KeyStandingWave
            call FS_SaveParaview(Hprv,etav,               &
                                 ufv,vfv,wfv,pfv,         &
                                 HprvA,etavA,             &
                                 ufvA,vfvA,wfvA,pfvA,     &
                                 uErrv,vErrv,wErrv,pErrv, &
                                 xvt,yvt,zvt,             &
                                 No_vp)
!           ---------------------
#           else
            call SaveParaviewVertex(Hprv,etav,                &
                                    ufv,vfv,wfv,pfv,          &
                                    usv,vsv,wsv,psv,          &
                                    uErrv,vErrv,wErrv,pErrv,  &
                                    xvt,yvt,zvt,              &
                                    No_vp)
#           endif
!           ---------------------
         endif
!        ---------------------------------------------------
!        Save Re-start data
         if (IrestartOUT.eq.1) then
             call restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp,       &
                              alphasnp,usnp,vsnp,wsnp,psnp,       &
                              zct)
         endif
      endif

!      _______________________________________________________________
!     |      |                                                        |
!     | F.2  |               Close time error file                    |
!     |______|________________________________________________________|

      IF (SaveFiles.eq.1) THEN
!        ---------------------------
         close(7100)
!        ---------------------------
#        ifdef KeyTaylorVortex
            close(8100)
#        endif
!        ---------------------------
#        ifndef KeyStandingWave
            close(9100)
#        endif
      ENDIF
!     ------------------------------
!     For the new version (Dic 2017)
!#     ifndef KeyFixedFreeSurface
!#        ifdef KeySave1DReference
!            call FS_SaveReferenceClose
!#        endif
!#     endif

!      _______________________________________________________________
!     |      |                                                        |
!     | F.3  |               Final simulation time                    |
!     |______|________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         if (FinRK==1) then
            print*,' '
            print*,'          ---------------------------------------'
            print*,'          ------------ Using Euler --------------'
            print*,'          ---------------------------------------'
            print*,' '
         elseif (FinRK==2) then
            print*,' '
            print*,'          ---------------------------------------'
            print*,'          ------------ Using RK-2 ---------------'
            print*,'          ---------------------------------------'
            print*,' '
         endif
         print*,'                                                        '
         print*,'    ___________________________________________________ '
         print*,'   |                                                   |'
         print*,'   |                   FINALIZATION                    |'
         print*,'   |___________________________________________________|'
         print*,'                                                        '
         print*,'     -------------------------------------------------- '
         print*,'                      Simulation Time                   '
         print*,'     -------------------------------------------------- '
         print*,'                                                        '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         IF (rang_topo.eq.1) THEN
         if (FinRK==1) then
            print*,' '
            print*,'          ======================================='
            print*,'          ============ Using Euler =============='
            print*,'          ======================================='
            print*,' '
         elseif (FinRK==2) then
            print*,' '
            print*,'          ======================================='
            print*,'          ============ Using RK-2 ==============='
            print*,'          ======================================='
            print*,' '
         endif
         print*,'                                                        '
         print*,'   ==================================================== '
         print*,'                    MPI: FINALIZATION                   '
         print*,'   ==================================================== '
         print*,'                                                        '
         print*,'   ____________________________________________________ '
         print*,'                                                        '
         print*,'                      Simulation Time                   '
         print*,'   ____________________________________________________ '
         print*,'                                                        '
         ENDIF
#     endif
!     =============== END ================
!     ====================================

      call cpu_time(finish_nsmp)
!     -------------------------------------------------------
!     Calculating the simulation time
      tcpu  = finish_nsmp-start_nsmp !etime(tt)
      idhr  = tcpu/3600
      idmin = tcpu/60-idhr*60
      idsec = tcpu-(idhr*3600+idmin*60)
!     -------------------------------------------------------
!     Print CPU Time per point
      nt=(tfin-time0)/dt+1
      usepas=tcpu/nt/(N_CELL*NZ)*1000

!     -------------------------------------------------------
!     Display elapsed time (CPU)

916   format(t6,'PROC:',i3,'  Elapsed time : ',f10.3,' sec CPU (',i2,':',i2,':',i2,')')
917   format(t6,'  Elapsed time : ',f10.3,' sec CPU (',i2,':',i2,':',i2,')')
918   format(t6,'          user : ',f10.3,' sec')
919   format(t6,'        system : ',f10.3,' sec')
920   format(t6,'Time per point : ',f10.3,' msec CPU /pas/point ')

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         write(*,917) tcpu,idhr,idmin,idsec
         write(*,918) tt(1)
         write(*,919) tt(2)
         print*,'  '
         write(*,920) usepas
         print*,'                                                        '
         print*,'     -------------------------------------------------- '
         print*,'                                                        '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
         write(*,916) rang_topo,tcpu,idhr,idmin,idsec
         call MPI_ALLREDUCE(tcpu,MAXtcpu,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                          comm3D,code)
         IF (rang_topo.eq.0) THEN
            print*,' '
            print*,'    Maximum time = ',MAXtcpu
            print*,' '
            print*,'   ==================================================== '
            print*,'                                                        '
         ENDIF
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________________________
!     |      |                                                        |
!     | F.4  |             Free memory: deallocating                  |
!     |______|________________________________________________________|

      call dealloc_geometry
      call dealloc_variables
!     ----------------------
#     ifdef KeySaveStatistics
      call dealloc_stats_variables
#     endif

!      _______________________________________________________________
!     |      |                                                        |
!     | F.5  |                 Closing programs                       |
!     |______|________________________________________________________|

8888  continue

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call finalisation_mpi
#     endif
!     =============== END ================
!     ====================================

!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
!     =========  PARALLEL CUDA ============
#     ifdef KeyCUDA
         call deallocate_cuda
#     endif
!     =============== END ================
!     *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
      print*,'                                                            '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'------------------------------------------------------------'
      print*,'                      ***  END  ***                         '
      print*,'------------------------------------------------------------'
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                            '
      print*,'                                                            '
!     ====================================
!     =====  START PARALLEL OPTION =======
#     else
      IF (rang_topo.eq.0) THEN
      print*,'                                                            '
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'============================================================'
      print*,'                  ***  END MPI PROGRAM ***                  '
      print*,'============================================================'
      print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
      print*,'                                                            '
      print*,'                                                            '
      ENDIF
#     endif
!     =============== END ================

      END PROGRAM

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                            END of nsmp3D                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
