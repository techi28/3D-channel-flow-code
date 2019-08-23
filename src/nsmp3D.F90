!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	    MAIN PROGRAM                              !
!       Solution of the Navier-Stokes Multi-phase equations           !
!                      Miguel Angel Uh Zapata                         !
!                  Last modification: Aug 2015                        !
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
!*********************************************************************!
!                                                                     !
!                            Definitions                              !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |            Keys, subroutines and parameters            |
!     |________________________________________________________|

#     include "cppdefs.h"
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
#     endif
            USE variables
            USE geometry
            USE interfaces
            USE stats
            implicit none
!     =============== END ================
!     ====================================

!      ________________________________________________________
!     |                                                        |
!     |           Definition of local variables                |
!     |________________________________________________________|

      real,dimension(2) :: tt
      real ::tcpu,MAXtcpu,tstart,tfinish,tnow,tpst
      real*8 :: nt,usepas
      real*8 :: time0,t1,t2
      real*8 :: x,y,z,UT,VT
      real*8 :: tInistats
      integer:: idmin,idhr,idsec
      integer:: m,SaveFiles,elem,ii
!     --------------------------
      integer:: nRK,nflag,IDisplay
!     --------------------------
      character*80 title
!     --------------------------
      nflag = 0
      IDisplay = 1
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

!*********************************************************************!
!                                                                     !
!                         0) Parallelization                          !
!                                                                     !
!*********************************************************************!
      call cpu_time(tstart)
!      _____________________________________________________________
!     |      |                                                      |
!     | 0.0  |     Read: Number of cell-center and vertex points    |
!     |______|______________________________________________________|

      open(25,file='data.txt',status='old')
      read(25,*) title
      read(25,*) N_VERTglobal
      read(25,*) N_CELL0global
      close(25)
!     _____________________________________________________________
!    |      |                                                      |
!    | 0.1  |               Defition of parameters                 |
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
         print*,'          Number of NZ points =',NZglobal
         print*,'                                                       '
#     endif
!     =============== END ================
!     ====================================

!     _____________________________________________________________
!    |      |                                                      |
!    | 0.2  |                Parallel structure                    |
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

         if(rang_topo .ne.0) IDisplay = 0
#     endif
!     =============== END ================
!     ====================================

!*********************************************************************!
!                                                                     !
!          1)           Variables & input data files                  !
!                                                                     !
!*********************************************************************!

!      _____________________________________________________________
!     |      |                                                      |
!     | 1.1  |                  Allocate variables                  |
!     |______|______________________________________________________|

      call alloc_variables
      call alloc_geometry
#     ifdef KeyTESTChannel
      tInistats = 60.0d0
      call alloc_stats_variables
#     endif
!      _____________________________________________________________
!     |      |                                                      |
!     | 1.2  |                    Input data                        |
!     |______|______________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
        call input(No_vp,No_cp,nbe,nbev,xv,yv,zbv)
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
        if(rang_topo .eq. 1) print *,"lecture de Parallel_In"
        call parallel_input_local(No_vp,No_cp,nbe,nbev,xv,yv,zbv)
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.3  |                 Input parameters                       |
!     |______|________________________________________________________|

      call input_parameters(nbev,hv,h,                     &
                            No_vp,No_cp,nbe,xv,yv,zbv)


!*********************************************************************!
!                                                                     !
!                           2) Initialization                         !
!                                                                     !
!*********************************************************************!

!      _______________________________________________________________
!     |      |                                                        |
!     | 2.1  |            Geometry variables of the mesh              |
!     |______|________________________________________________________|

      call calcul_geometry(xc,yc,sig,dsig,No_cp,nbe,h, &
                           xv,yv,sigv,dsigv,No_vp,nbev,&
                           ic1tec,ic2tec,ic3tec)
!      _______________________________________________________________
!     |      |                                                        |
!     | 2.2  |             Boundary Condition                         |
!     |______|________________________________________________________|

      call calcul_bc(xc,yc,sig,dsig,No_cp,nbe, &
                           xv,yv,sigv,dsigv,No_vp,nbev)
!      _______________________________________________________________
!     |       |                                                        |
!     | 2.2a  |             Pressure Martix(Ax=b)                      |
!     |______ |________________________________________________________|
#     ifdef KeyPPECenter
       call PDE_center(xc,yc,sig,dsig,No_cp,nbe,    &
                       xv,yv,sigv,dsigv,No_vp,nbev)
#     endif

!#     ifdef KeyTESTpBC
!      call check_pbc(xc,yc,sig,dsig,No_cp,nbe, &
!                     xv,yv,sigv,dsigv,No_vp,nbev)
!#     endif

!      _______________________________________________________________
!     |      |                                                        |
!     | 2.3  |                 Initial conditions                     |
!     |______|________________________________________________________|

!     ___________________________________________
!     Initial time and step

      time  = 0.0d0
      nstep = 0
      scount = 0
!     ___________________________________________
!     Initial values of the main variables

      call initial(alphafn,ufn,vfn,wfn,pfn,viscof,rhof,   &
                   alphasn,usn,vsn,wsn,psn,viscos,rhos,   &
                   alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv, &
                   alphasv,usv,vsv,wsv,psv,viscosv,rhosv, &
                   etan,etav,Hpr,Hprv,                    &
                   xct,yct,zct,                           &
                   xvt,yvt,zvt,                           &
                   xc,yc,sig,dsig,No_cp,nbe,              &
                   xv,yv,sigv,dsigv,No_vp,nbev,           &
                   h,hv)
!    -------------------
!      Initial time averaged stats
        if(IStats .eq. 1) then
            call initial_stats(suf,svf,swf,spf)
        endif
!      _______________________________________________________________
!     |      |                                                        |
!     | 2.4  |               Re-start solution (n)                    |
!     |______|________________________________________________________|
!     ________________________________________________________
!     Save counter
      SaveCounter = 0
!    ------------------------------------------------------------------
        call user_initial(ufn,vfn,wfn,pfn,            &
                         ufv,vfv,wfv,pfv,             &
                         xc,yc,sig,dsig,No_cp,nbe,    &
                         xv,yv,sigv,dsigv,No_vp,nbev, &
                         Hpr,h,etan,                  &
                         Hprv,hv,etav,                &
                         xct,yct,zct,                 &
                         xvt,yvt,zvt)
!    ==================================================================
#       ifdef KeyTESTChannel
         call  sm_init()
#       endif
!      _______________________________________________________________
!     |      |                                                        |
!     | 2.5  |             Save initial condition values              |
!     |______|________________________________________________________|

!     -------------------------------------------------------
!     Last time when the results were saved
      LastTimeSave = time

      if (ChooseExit.eq.1) then
!        ---------------------------------------------------
!        Save Tecplot file vertex
         call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,              &
                            alphasv,usv,vsv,wsv,psv,rhosv,        &
                            xvt,yvt,zvt,No_vp)
      elseif(ChooseExit.eq.2) then
!        ---------------------------------------------------
!        Save Tecplot file cell centers
         call SavetecCenter(alphafn,ufn,vfn,wfn,pfn,              &
                            alphasn,usn,vsn,wsn,psn,rhos,         &
                            xvt,yvt,zvt,No_vp)

      elseif (ChooseExit.eq.3) then
!        ---------------------------------------------------
!        Save Tecplot file cell centers and vertex
            call SavetecVC(alphafn,ufn,vfn,wfn,pfn,      &
                           alphasn,usn,vsn,wsn,psn,rhos, &
                           xct,yct,zct,No_cp,                 &
                           alphafv,ufv,vfv,wfv,pfv,           &
                           alphasv,usv,vsv,wsv,psv,rhosv,     &
                           xvt,yvt,zvt,No_vp)
      elseif (ChooseExit.eq.4) then
!        ---------------------------------------------------
!        Save paraview file vertex
         call SaveParaviewVertex(alphafv,ufv,vfv,wfv,pfv,         &
                                 alphasv,usv,vsv,wsv,psv,rhosv,   &
                                 xvt,yvt,zvt,No_vp)
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
        if(IDisplay .eq. 1) then
            print*,'                                                        '
            print*,'   ==================================================== '
            print*,'                    MPI: TIME SIMULATIONS               '
            print*,'   ==================================================== '
            print*,'                                                        '
        endif
#     endif
!     =============== END ================
!     ====================================

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.1  |          Time loop:  Initial values: (n+1)=(n)         |
!     |______|________________________________________________________|

      call LoopTimeInitial(etanp,                              &
                           alphafnp,ufnp,vfnp,wfnp,pfnp,       &
                           alphasnp,usnp,vsnp,wsnp,psnp,       &
                           xct,yct,zct,                        &
                           xvt,yvt,zvt,                        &
                           etan,                               &
                           alphafn,ufn,vfn,wfn,pfn,            &
                           alphasn,usn,vsn,wsn,psn,            &
                           xc,yc,sig,Hpr,h,                    &
                           xv,yv,sigv,Hprv,hv,                 &
                           No_cp,nbe)
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.2  |         Time loop:  Beggining of the new time step     |
!     |______|________________________________________________________|

      flag_ab = 0

10    continue
!     ----------------------------------
!     get time usage per step
#        ifdef KeyParallel
            call MPI_Barrier(comm3D,code)
            t1 = MPI_Wtime()
#        endif

      time  = time + dt
      nstep = nstep + 1
!     --------------------------
!     Display current time step
        if(IDisplay .eq. 1) then
            print*,'  '
            write(*,'(a8,i8)')  'Step : ',nstep
            write(*,'(a8,f8.4)')'Time = ',time
            print*,'  '
        endif
!      ________________________________________________________
!     |      |                                                 |
!     | 3.2.1|  Runge-Kutta: Initial known values: (known)=(n) |
!     |______|_________________________________________________|

      do i=1,N_CELL
!        ----------------------------
!        Free surface
         eta(i)=etan(i)
         do k=1,NZ
!           -------------------------
!           Volume Fraction
            alphaf(i,k) = alphafn(i,k)
            alphas(i,k) = alphasn(i,k)
!           -------------------------
!           Velocity
            uf(i,k) = ufn(i,k)
            vf(i,k) = vfn(i,k)
            wf(i,k) = wfn(i,k)
            us(i,k) = usn(i,k)
            vs(i,k) = vsn(i,k)
            ws(i,k) = wsn(i,k)
!           -------------------------
!           Pressure
            pf(i,k) = pfn(i,k)
            ps(i,k) = psn(i,k)
        enddo
      enddo

!      ________________________________________________________
!     |      |                                                 |
!     | 3.2.2|  Runge-Kutta: The two step loop                 |
!     |______|_________________________________________________|

        DO nRK=1,FinRK

            IF (IDisplay .eq. 1) THEN
                if(flag_ab .eq. 0) then
                    write(*,'(t15,a8,i3)')'RK step:',nRK
                else
                    write(*,*)'         ====Adams-Bashforth 2nd order==='
                endif
            ENDIF
!         _______________________________________
!        |                                       |
!        |            Update variables           |
!        |_______________________________________|
         call hydro(ufnp,ufn,uf,ufv,                          &
                    vfnp,vfn,vf,vfv,                          &
                    wfnp,wfn,wf,wfv,                          &
                    pfnp,pfn,pf,pfv,                          &
!                   -------------------------------------------
                    rhof,viscof,                              &
                    rhofv,viscofv,                            &
!                   -------------------------------------------
                    xc,yc,sig,dsig,No_cp,nbe,                 &
                    xv,yv,sigv,dsigv,No_vp,nbev,              &
!                   --------------------------------------------
                    Hpr,h,etanp,etan,                         &
                    Hprv,hv,etav,                             &
!                   --------------------------------------------
                    flag_ab)
!         _______________________________________
!        |                                       |
!        |   Update RK2 values:(known)=(n+1)     |
!        |_______________________________________|
       if(FinRK .eq. 2) then
         do i=1,N_CELL
!           ----------------------------
!           Free surface
            eta(i) = etanp(i)
            do k=1,NZ
!              -------------------------
!              Volume Fraction
               alphaf(i,k) = alphafnp(i,k)
               alphas(i,k) = alphasnp(i,k)
!              -------------------------
!              Velocity
               uf(i,k) = ufnp(i,k)
               vf(i,k) = vfnp(i,k)
               wf(i,k) = wfnp(i,k)
               us(i,k) = usnp(i,k)
               vs(i,k) = vsnp(i,k)
               ws(i,k) = wsnp(i,k)
!              -------------------------
!              Pressure
               pf(i,k) = pfnp(i,k)
               ps(i,k) = psnp(i,k)
           enddo
         enddo
       endif

      ENDDO
!     _________________________________________________________
!     In the case of only one-step case
        IF (FinRK .eq. 1) THEN
            if (IDisplay.eq.1) then
#           ifndef KeyImplicit
                if(flag_ab .eq. 0) then
                    print*,' '
                    print*,'          ======================================='
                    print*,'          ============ Using RK-1 ==============='
                    print*,'          ======================================='
                    print*,' '
                endif
#           else
                    print*,' '
                    print*,'          ======================================='
                    print*,'          ============ Using CN-2 ==============='
                    print*,'          ======================================='
                    print*,' '
#           endif
            endif
!      ________________________________________________________
!     |      |                                                 |
!     | 3.2.3|  Runge-Kutta: Final solution of the RK2 formula |
!     |______|_________________________________________________|
      ELSE
          do i=1,N_CELL
             do k=1,NZ
!              -------------------------
!              Volume Fraction
                alphafnp(i,k) = 0.5d0*(alphafn(i,k)+alphafnp(i,k))
                alphasnp(i,k) = 0.5d0*(alphasn(i,k)+alphasnp(i,k))
!              -------------------------
!              Velocity
                ufnp(i,k) = 0.5d0*(ufn(i,k)+ufnp(i,k))
                vfnp(i,k) = 0.5d0*(vfn(i,k)+vfnp(i,k))
                wfnp(i,k) = 0.5d0*(wfn(i,k)+wfnp(i,k))
                usnp(i,k) = 0.5d0*(usn(i,k)+usnp(i,k))
                vsnp(i,k) = 0.5d0*(vsn(i,k)+vsnp(i,k))
                wsnp(i,k) = 0.5d0*(wsn(i,k)+wsnp(i,k))
             enddo
          enddo
      ENDIF
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.3  |     Time loop: Update the variables & input box        |
!     |______|________________________________________________________|

      call LoopTimeUpdate(etan,                                &
                          alphafn,ufn,vfn,wfn,pfn,             &
                          alphasn,usn,vsn,wsn,psn,             &
                          xct,yct,zct,                         &
                          xvt,yvt,zvt,                         &
                          etanp,                               &
                          alphafnp,ufnp,vfnp,wfnp,pfnp,        &
                          alphasnp,usnp,vsnp,wsnp,psnp,        &
                          xc,yc,sig,Hpr,h,                     &
                          xv,yv,sigv,Hprv,hv,                  &
                          No_cp,nbe)
!      _______________________________________________________________
!     |       |                                                        |
!     | 3.4b  |      Obtain time average data for cylinder case        |
!     |______ |________________________________________________________|
        if(IStats .eq. 1) then
            if(mod(nstep,dtsample) .eq. 0) then
                call update_stats(suf,svf,swf,spf,   &
                    ufnp,vfnp,wfnp,pfnp)
            endif
        endif
!      _______________________________________________________________
!     |       |                                                        |
!     | 3.4c  |                Stats for Open Channel Case             |
!     |______ |________________________________________________________|
#     ifdef KeyTESTChannel
       IF(time .gt. tInistats) then
         if(mod(nstep,dtsample) .eq. 0) then
            call glob_sample(ufnp,vfnp,wfnp,pfnp)
         endif

         if(mod(nstep,dtplane) .eq. 0) then
            call glob_put_plane_mm()
         endif
       ENDIF
#     endif
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.5  |    Time loop:  Saving results during simulation        |
!     |______|________________________________________________________|

      if(mod(nstep,ncfl) .eq. 0) then


      if(IDisplay .eq. 1) then
         print*,' '
         write(*,'(t15,a10,i6)')'CHECK CFL at:',nstep
      endif
!     --------------------------------------------------
!     check the velocity field
         call glob_cfl(ufnp,vfnp,wfnp,pfnp,           &
                       xc,yc,sig,dsig,No_cp,nbe)
      endif

      if ((time.ge.tInisave).and.&
         (1d-08.ge.dtsave-(time-LastTimeSave))) then
!        ---------------------------------------------------
!        Last time saved
         SaveCounter = SaveCounter+1
         LastTimeSave = time
!        ---------------------------------------------------
!        Save Tecplot file vertices
         if (ChooseExit.eq.1) then
            call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,           &
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save Tecplot file cell-centers
         elseif (ChooseExit.eq.2) then
            call SavetecCenter(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save Tecplot file cell-centers & vertex
         elseif (ChooseExit.eq.3) then
            call SavetecVC(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                           alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                           xct,yct,zct,No_cp,                 &
                           alphafv,ufv,vfv,wfv,pfv,           &
                           alphasv,usv,vsv,wsv,psv,rhosv,     &
                           xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save paraview file vertex
         elseif (ChooseExit.eq.4) then
            call SaveParaviewVertex(alphafv,ufv,vfv,wfv,pfv,      &
                                    alphasv,usv,vsv,wsv,psv,rhosv,&
                                    xvt,yvt,zvt,No_vp)
         endif
!        ---------------------------------------------------
!        Save Tecplot file for stats
         if(IStats .eq. 1) then
            call SavetecStats(suf,svf,swf,spf, &
                              xvt,yvt,zvt,No_vp)
         endif

      endif

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.6  |    Time loop:  Criteria to finish time simulations     |
!     |______|________________________________________________________|
#       ifdef KeyParallel
         t2 = MPI_Wtime()
         t2 = t2 -t1
         call MPI_Barrier(comm3D,code)
         call MAX_parallel(t2, tnsmp)
         if(IDisplay .eq. 1) then
            print*, 'Time per step    =',tnsmp
            print*,'  '
         endif
#       endif

      if ((tfin-time).le.1d-08) goto 9999
      call cpu_time(tnow)
      tpst = tnow - tstart
      if(tpst .gt. tlimit) goto 9999
      goto 10

9999  continue


!*********************************************************************!
!                                                                     !
!                          4) Finalization                            !
!                                                                     !
!*********************************************************************!

!      _______________________________________________________________
!     |      |                                                        |
!     | 4.1  |                Saving final results                    |
!     |______|________________________________________________________|

      if ((tfin-LastTimeSave).ge.1d-8) then
!        ---------------------------------------------------
!        Save Tecplot file vertices
        if ((ChooseExit.eq.1)) then
            if(IDisplay .eq. 1) then
                print*,'   --- Tecplot final time solution printed for V ---'
            endif
            call SavetecVertex(alphafv,ufv,vfv,wfv,pfv,           &
                               alphasv,usv,vsv,wsv,psv,rhosv,     &
                               xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save Tecplot file cell centers
         elseif (ChooseExit.eq.2) then
            if(IDisplay .eq. 1) then
                print*,'   --- Tecplot final time solution printed for C ---'
            endif
            call SavetecCenter(alphafnp,ufnp,vfnp,wfnp,pfnp,      &
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                               xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save Tecplot file cell centers & vertex
        elseif (ChooseExit.eq.3) then
            if(IDisplay .eq. 1) then
                print*,'   --- Tecplot final time solution printed for VC ---'
            endif
            call SavetecVC(alphafnp,ufnp,vfnp,wfnp,pfnp,&
                alphasnp,usnp,vsnp,wsnp,psnp,rhos, &
                xct,yct,zct,No_cp,                 &
                alphafv,ufv,vfv,wfv,pfv,           &
                alphasv,usv,vsv,wsv,psv,rhosv,     &
                xvt,yvt,zvt,No_vp)
!        ---------------------------------------------------
!        Save paraview file cell centers & vertex
        elseif (ChooseExit.eq.4) then
            if(IDisplay .eq. 1) then
                print*,'   --- Paraview final time solution printed for V ---'
            endif
            call SaveParaviewVertex(alphafv,ufv,vfv,wfv,pfv,      &
                                    alphasv,usv,vsv,wsv,psv,rhosv,&
                                    xvt,yvt,zvt,No_vp)
        endif
!        ---------------------------------------------------
!        Save Tecplot file for stats file
            if(IStats .eq. 1) then
                call SavetecStats(suf,svf,swf,spf, &
                    xvt,yvt,zvt,No_vp)
                if(IDisplay .eq. 1) then
                    print*,'   --- Stats final time solution printed ---'
                    print*,'   --- Final Stats restart file data saved! ---'
                endif
            endif

         endif
!      _______________________________________________________________
!     |      |                                                        |
!     | 4.3  |                  Final simulation time                 |
!     |______|________________________________________________________|

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
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

!     -------------------------------------------------------
!     Calculating the simulation time
      call cpu_time(tfinish)
      tcpu = tfinish - tstart
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
!     | 4.4  |             Free memory: deallocating                  |
!     |______|________________________________________________________|

      call dealloc_variables
      call dealloc_geometry
#     ifdef KeyTESTChannel
      call dealloc_stats_variables
#     endif
!      _______________________________________________________________
!     |      |                                                        |
!     | 4.5  |                 Closing programns                      |
!     |______|________________________________________________________|

8888  continue

!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         call finalisation_mpi
#     endif
!     =============== END ================
!     ====================================

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
!                            END OF nsmp3D                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
