!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               ---  TEST MPI COMMUNICATION 3D ---                    !
!                      Miguel Angel Uh Zapata                         !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

       SUBROUTINE TestMPInsmp3D(No_cp,nbe)

!---------------------------------------------------------------------!
!                                                                     !
!     Programe to test definitions and comunications of the MPI       !
!     formulation. I define a function which sum its three neighbors  !
!     values and itself. I calculate the sum of all independent       !
!     elements and use an MPI command to collect the values in all    !
!     processors. The final result is compare with the global solu-   !
!     tion after a determinated number of iterations.                 !
!                                                                     !
!---------------------------------------------------------------------!  

!*********************************************************************!
!                                                                     !
!                        Definition of variables                      !
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

      integer,dimension(:,:):: No_cp(N_CELL,3)
      integer,dimension(:)  :: nbe(N_CELL0)
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      real,dimension(2) :: tt
      real :: tcpu,MAXtcpu
      integer:: idmin,idhr,idsec
!     --------------------------------
      integer :: jj,kk
      integer :: iter,itermax
      integer :: nc1,nc2,nc3
      real*8  :: f1,f2,f3
      real*8  :: suma,sumaglob,SumaNoMPI,SumaMPI
!     --------------------------------
      real*8 ,dimension(:,:), allocatable :: fun1_global,fun2_global
      real*8 ,dimension(:,:), allocatable :: fun1,fun2,fun2old
      real*8, dimension(:),   allocatable :: xc_global,yc_global 
      real*8, dimension(:),   allocatable :: xc_local,yc_local
      integer :: s,ielem,jv1,jv2,jv3,Display
!     ----------------------------------------
      real :: start,finish,startC,finishC
      real :: TimeGlobal,TimeLocal,TimeComm,TimeInner   
!     ----------------------------------------
      real :: SpL,SpI,SpLA,SpIA
      real :: TimeGlobalA,TimeLocalA,TimeCommA,TimeInnerA
      integer :: EpL,EpI,EpLA,EpIA,PercCommA
      real*8, dimension(Nprocs) :: SumaNoMPIV,SumaMPIV
      real,   dimension(Nprocs) :: TimeGlobalV,TimeLocalV,TimeCommV,TimeInnerV 
      real,   dimension(Nprocs) :: SpLV,SpIV
      integer,dimension(Nprocs) :: EpLV,EpIV  
      integer,dimension(Nprocs) :: PercCommV
!     ----------------------------------------      
      
      itermax = 1000
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (rang_topo.eq.0) THEN
         write(*,*) ' '     
         write(*,'(t6,60a)'), '>>>>> Begin subroutine: TestMPI3D'
      ENDIF
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     ====================================
!     =====    DISPLAY ITERATIONS  =======
      Display = 1
#     ifdef KeyParallel
         if (rang_topo.ne.0) then
            Display = 0
         endif
#     endif	
!     =============== END ================    
!     ====================================
      TimeComm = 0.0d0
	
!*********************************************************************!
!                                                                     !
!                      Functions: global and local                    !
!                                                                     !
!*********************************************************************!
!      __________________________________
!     |                                  |
!     |             Global               |
!     |__________________________________| 	
      
!     ___________________________________
!     allocate
      allocate(fun1_global(N_CELL0Global,NZglobal))
      allocate(fun2_global(N_CELL0Global,NZglobal))
      allocate(xc_global(N_CELL0global))
      allocate(yc_global(N_CELL0global))
!     ___________________________________
!     global cell-center
      do i=1,N_CELL0global
	 jv1=No_vp_global(i,1)
	 jv2=No_vp_global(i,2)
	 jv3=No_vp_global(i,3)              
	 xc_global(i)=(xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
	 yc_global(i)=(yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
      enddo
!     __________________________________
!     Initial function          
      do k=1,NZglobal
         do i=1,N_CELL0Global
            fun1_global(i,k) = k*dsin(xc_global(i)+yc_global(i))**2 
         enddo
      enddo      
!      __________________________________
!     |                                  |
!     |              Local               |
!     |__________________________________| 

!     ___________________________________
!     allocate
      allocate(fun1(N_CELL,NZ),fun2(N_CELL,NZ),fun2old(N_CELL,NZ))
      allocate(xc_local(N_CELL))
      allocate(yc_local(N_CELL))
!     ___________________________________
!     local cell-center
      do i=1,N_CELL
         ielem = Index_global(i)
         xc_local(i) = xc_global(ielem)
         yc_local(i) = yc_global(ielem)
      enddo
!     __________________________________
!     Initial function       
      do k=1,NZ
         do i=1,N_CELL
            fun1(i,k) = k*dsin(xc_local(i)+yc_local(i))**2  
         enddo
      enddo
      
!*********************************************************************!
!                                                                     !
!                    Serial (done for each processor)                 !
!                                                                     !
!*********************************************************************!


      IF (Display.eq.1) THEN
      print*,'                                                        '
      print*,'    ___________________________________________________ '
      print*,'   |                                                   |'
      print*,'   |TESTING... EXAMPLE MPI 3D UNSTRUCTURED GRID: GLOBAL|'
      print*,'   |___________________________________________________|'
      print*,'                                                        '
      ENDIF

      call cpu_time(start)     
      
      iter = 0
11    continue
      iter = iter + 1
!     ___________________________________
!     Calculations      
      sumaNoMPI = 0.0d0
      do k=1,NZglobal      
         do i=1,N_CELL0global
            if (nbe_global(i).eq.0) then
                nc1 = No_cp_global(i,1)
                nc2 = No_cp_global(i,2)
                nc3 = No_cp_global(i,3)
                f1  = fun1_global(nc1,k)
                f2  = fun1_global(nc2,k)
                f3  = fun1_global(nc3,k)
                fun2_global(i,k) = dsin(fun1_global(i,k) + f1 + f2 + f3)
            else
                fun2_global(i,k) = fun1_global(i,k)
            endif
            sumaNoMPI = sumaNoMPI + fun2_global(i,k)
         enddo
      enddo
!     ___________________________________
!     Update       
      do k=1,NZglobal
         do i=1,N_CELL0global
            fun1_global(i,k) = fun2_global(i,k)
         enddo
      enddo      
!     ___________________________________
!     Criteria
      if (iter.lt.itermax) then
         goto 11
      else
         goto 12
      endif
12   continue      
      
      call cpu_time(finish)
      timeGlobal = finish-start      

!*********************************************************************!
!                                                                     !
!                              Parallel                               !
!                                                                     !
!*********************************************************************!


      IF (Display.eq.1) THEN
      print*,'    ___________________________________________________ '
      print*,'   |                                                   |'
      print*,'   |TESTING... EXAMPLE MPI 3D UNSTRUCTURED GRID: LOCAL |'
      print*,'   |___________________________________________________|'
      print*,'                                                        '
      ENDIF

      call cpu_time(start)			

      iter = 0
 111  continue
      iter = iter + 1
!     ___________________________________
!     Calculations
      suma = 0.0d0
      do k=1,NZ      
         do i=1,N_CELL0
            if (nbe(i).eq.0) then
                nc1 = No_cp(i,1)
                nc2 = No_cp(i,2)
                nc3 = No_cp(i,3)
                f1  = fun1(nc1,k)
                f2  = fun1(nc2,k)
                f3  = fun1(nc3,k)
                fun2(i,k) = dsin(fun1(i,k) + f1 + f2 + f3) 
            else
                fun2(i,k) = fun1(i,k)
            endif
            suma = suma + fun2(i,k)
         enddo
      enddo 
!     __________________________________
!     Communication
      call cpu_time(startC)
      !call communication3D(fun2)
      !call communication3Dtype1(fun2)      
      call communication3Dtype2(fun2)
      call SUM_parallel(suma,sumaMPI)  
      call cpu_time(finishC)
      TimeComm = TimeComm + (finishC-startC)
!     ___________________________________
!     Update               
      do k=1,NZ      
         do i=1,N_CELL
            fun1(i,k) = fun2(i,k)  
         enddo
      enddo
!     ___________________________________
!     Criteria
      if (iter.lt.itermax) then
         goto 111
      else
         goto 112
      endif

112   continue

      call cpu_time(finish)
      TimeLocal = finish-start
      TimeInner = TimeLocal - TimeComm
      
!*********************************************************************!
!                                                                     !
!                           Finalization                              !
!                                                                     !
!*********************************************************************!
!      __________________________________
!     |                                  |
!     |            Deallocate            |
!     |__________________________________|

      deallocate(fun1_global,fun2_global)
      deallocate(xc_global,yc_global)
      deallocate(fun1,fun2,fun2old)
      deallocate(xc_local,yc_local)
      
!      __________________________________
!     |                                  |
!     |         Display results          |
!     |__________________________________|

      call MPI_ALLGATHER(sumaNoMPI,1,MPI_DOUBLE_PRECISION,sumaNoMPIV,1,MPI_DOUBLE_PRECISION,comm3D,code)
      call MPI_ALLGATHER(sumaMPI,1,MPI_DOUBLE_PRECISION,sumaMPIV,1,MPI_DOUBLE_PRECISION,comm3D,code)     
      call MPI_ALLGATHER(TimeGLobal,1,MPI_FLOAT,TimeGLobalV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(TimeLocal,1,MPI_FLOAT,TimeLocalV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(TimeInner,1,MPI_FLOAT,TimeInnerV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(TimeComm,1,MPI_FLOAT,TimeCommV,1,MPI_FLOAT,comm3D,code)
      TimeGlobalA = sum(TimeGLobalV)/Nprocs 
      TimeLocalA  = sum(TimeLocalV)/Nprocs
      TimeInnerA  = sum(TimeInnerV)/Nprocs
      TimeCommA   = sum(TimeCommV)/Nprocs 
      SpL = TimeGLobalA/TimeLocal
      SpI = TimeGLobalA/TimeInner
      EpL = floor(SpL/Nprocs*100)
      EpI = floor(SpI/Nprocs*100)     
      call MPI_ALLGATHER(SpL,1,MPI_FLOAT,SpLV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(SpI,1,MPI_FLOAT,SpIV,1,MPI_FLOAT,comm3D,code)
      call MPI_ALLGATHER(EpL,1,MPI_INTEGER,EpLV,1,MPI_INTEGER,comm3D,code)
      call MPI_ALLGATHER(EpI,1,MPI_INTEGER,EpIV,1,MPI_INTEGER,comm3D,code)  
      SpLA =  sum(SpLV)/Nprocs
      SpIA =  sum(SpIV)/Nprocs
      EpLA = floor(SpLA/Nprocs*100)
      EpIA = floor(SpIA/Nprocs*100)      
      
      IF (Display.eq.1) THEN
      write(*,*),' '
      write(*,*),'====================================================================================== '
      write(*,*),'                            Results in all procerssors                                '
      write(*,*),'                               Iterations=',iterMAX
      write(*,*),'====================================================================================== '
      write(*,*),' '
      write(*,*),'                          CALCULATIONS & TIME COMPARISON                              '
      write(*,*),'--------------------------------------------------------------------------------------'         
      write(*,*),'  p  |    sumaNoMPI   |    sumaMPI     | Time Glob  | Time Local :   Inner     Comm'
      write(*,*),'-----|----------------|----------------|------------|---------------------------------'
      19 format(I4,a3,f14.4,a3,f14.4,a3,f10.5,a3,f10.5,a3,f10.5,f10.5)
      do s=1,Nprocs
         write(*,19),s-1,' |',sumaNoMPIV(s),' |',sumaMPIV(s),' |',&
               TimeGLobalV(s),' |',TimeLocalV(s),' :',TimeInnerV(s),TimeCommV(s)      
      enddo
      write(*,*),'---------------------------------------------------------------------------------------'
      18 format(a7,f14.4,a3,f14.4,a3,f10.5,a3,f10.5,a3,f10.5,f10.5)
      write(*,18),' Avg |',sumaNoMPIV(1),' |',sumaMPIV(1),' |',&
               TimeGLobalA,' |',TimeLocalA,' :',TimeInnerA,TimeCommA
      write(*,*),' '   
      write(*,*),'                               SPEEDUP & EFFICIENCY                                    '
      write(*,*),'---------------------------------------------------------------------------------------'      
      write(*,*),'  p  | Time Glob  | Time Local     Sp     Ep  |  Time Inner   Sp    Ep   | % of commun.'
      write(*,*),'-----|------------|---------------------------|--------------------------|-------------'
      17 format(I4,a3,f10.5,a3,f10.5,f7.2,a2,i5,a5,f10.5,f6.2,a2,i4,a4,i4,a2)      
      do s=1,Nprocs 
         PercCommV(s) = floor(TimeCommV(s)/TimeLocalV(s)*100)
         write(*,17),s-1,' |',TimeGLobalA,' |',&
                    TimeLocalV(s),SpLV(s),'x',EpLV(s),'\% | ',&
                    TimeInnerV(s),SpIV(s),'x',EpIV(s),'\% | ',PercCommV(s),'\%'
      enddo
      PercCommA   = ceiling(1.0d0*sum(PercCommV)/Nprocs)
      write(*,*),'---------------------------------------------------------------------------------------'
      16 format(a7,f10.5,a3,f10.5,f7.2,a2,i5,a5,f10.5,f6.2,a2,i4,a4,i4,a2) 
      write(*,16),' Avg |',TimeGlobalA,' |',TimeLocalA,SpLA,'x',EpLA,'\% | ',&
                         TimeInnerA,SpIA,'x',EpIA,'\% | ',PercCommA,'\%' 
      write(*,*),'======================================================================================= '
      write(*,*),' '
      ENDIF

!      __________________________________
!     |                                  |
!     |    Displaying simulation time    |
!     |__________________________________|

917   format(t6,'  Elapsed time : ',f10.3,' sec CPU (',i2,':',i2,':',i2,')')
918   format(t6,'          user : ',f10.3,' sec')
919   format(t6,'        system : ',f10.3,' sec') 

!     -------------------------------------------------------
!     Calculating the simulation time     
      tcpu  = etime(tt)
      call MPI_ALLREDUCE(tcpu,MAXtcpu,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                          comm3D,code)
      tcpu = MAXtcpu                    
      idhr  = tcpu/3600
      idmin = tcpu/60-idhr*60
      idsec = tcpu-(idhr*3600+idmin*60)

!     -------------------------------------------------------
!     Displaying 
      IF (Display.eq.1) THEN
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'            Maximum Simulation Time in Parallel         '
      print*,'     ================================================== '
      print*,'                                                        '
      write(*,917), tcpu,idhr,idmin,idsec
      write(*,918), tt(1)
      write(*,919), tt(2)   
      print*,'                                                        '
      print*,'     ================================================== '
      print*,'                                                        '
      ENDIF
      
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (Display.eq.1) THEN
         write(*,'(t6,60a)'), '<<<<< End   subroutine: TestMPI3D'
         print*,' '
      ENDIF
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	         
      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                           END of TestMPI                            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
