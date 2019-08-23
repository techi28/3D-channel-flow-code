!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             MODULE PARALLEL 3D PROBLEM (UNSTRUCTURED MESH)          !
!                      Miguel Angel Uh Zapata                         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!      SUBROUTINES:     1.1)  initialisation_mpi                      !
!                       1.2)  parallel_neighbors                      !
!                       1.3)  parallel_topology                       !
!                       1.4)  parallel_index                          !
!                       1.5)  parallel_shareindex                     !
!                       1.6)  finalisation_mpi                        !
!                                                                     !
!---------------------------------------------------------------------!

      MODULE parallel

      implicit none
#     include "mpif.h"
#     include "common.mpf"
!     ------------------------------------
      integer :: code,rang_topo
      integer :: comm3D,comm2
      integer :: Nprocs,Nindex,Nedges
      integer :: NindexNZ,NedgesNZ
      integer :: sNZ,NZLayer
      integer :: TotalSubdomains
      integer :: NeighborNumber
      integer :: NextraMax,NshareMax
      integer :: MPIOption
!     ------------------------------------
      integer,dimension(:), allocatable :: edges
      integer,dimension(:), allocatable :: indexval
      integer,dimension(:), allocatable :: Neighbors
!     -------------------------------------
      integer,dimension(:), allocatable :: edgesNZ
      integer,dimension(:), allocatable :: indexvalNZ
      integer,dimension(:), allocatable :: NeighborsNZ
!     ------------------------------
      integer,dimension(:), allocatable :: type_send
      integer,dimension(:), allocatable :: type_recv
      integer,dimension(:), allocatable :: type_sendVEC
      integer,dimension(:), allocatable :: type_recvVEC
!      integer,dimension(:), allocatable :: type_BdyCCVEC
      integer :: type_blocV,type_blocC,type_bloc2D
!     ------------------------------------
      integer,dimension(:,:), allocatable :: CCdom
      integer,dimension(:,:), allocatable :: DimCCdom
      integer,dimension(:,:), allocatable :: ExtraCCdom
      integer,dimension(:,:), allocatable :: DimExtraCCdom
      integer,dimension(:,:), allocatable :: VVdom
      integer,dimension(:,:), allocatable :: DimVVdom
      integer,dimension(:,:), allocatable :: BdyCCdom
      integer,dimension(:,:), allocatable :: DimBdyCCdom
!     ------------------------------
      integer,dimension(:,:), allocatable :: IniExtraFrom
      integer,dimension(:,:), allocatable :: FinExtraFrom
      integer,dimension(:,:), allocatable :: DimExtraFrom
!     ------------------------------
      integer,dimension(:),   allocatable :: Index_global
      integer,dimension(:),   allocatable :: Index_globalv
      integer,dimension(:),   allocatable :: IniExtraIndex_local
      integer,dimension(:,:), allocatable :: BlockBdyIndex
      integer,dimension(:),   allocatable :: BlockBdyNumber
!      integer,dimension(:,:), allocatable :: BlockBdyCCIndex
!      integer,dimension(:),   allocatable :: BlockBdyCCNumber
      integer,dimension(:,:), allocatable :: InterCCdom
      integer,dimension(:,:), allocatable :: DimInterCCdom
      integer,dimension(:,:), allocatable :: BoundCCdom
      integer,dimension(:,:), allocatable :: DimBoundCCdom
      real*8, dimension(:),   allocatable :: SharePhi
      real*8, dimension(:,:), allocatable :: SharePhiNew
      integer,dimension(:), allocatable :: index_neigbour
!     ------------------------------
      real*8,dimension(:), allocatable :: xc_global
      real*8,dimension(:), allocatable :: yc_global
      integer,dimension(:),   allocatable :: TagParallel
      integer,dimension(:,:), allocatable :: No_vp_global
      integer,dimension(:,:), allocatable :: No_cp_global
      integer,dimension(:),   allocatable :: nbe_global
!     -----------------------------------------------
      real*8, dimension(:),   allocatable :: xv_global
      real*8, dimension(:),   allocatable :: yv_global
      real*8, dimension(:),   allocatable :: zbv_global
      integer,dimension(:), allocatable :: nbev_global
      integer,dimension(:), allocatable :: vcount_global
      integer,dimension(:), allocatable :: TagXV_global
      integer,dimension(:), allocatable :: TagYV_global
      integer :: Nglobalghost
!     -------------------------------
#     ifdef KeyTESTpBC
      integer ::  NtotalMax
      integer,dimension(:), allocatable :: vpair_global
      integer,dimension(:), allocatable :: vedge_global
      integer,dimension(:), allocatable :: ghost_local
      integer,dimension(:), allocatable :: aux_new
      integer,dimension(:,:), allocatable :: ghost_index
      integer,dimension(:,:), allocatable :: ghostCCdom_aux
      integer,dimension(:,:), allocatable :: ghostCCdom
      integer,dimension(:,:), allocatable :: DimghostCCdom
      integer,dimension(:,:), allocatable :: DimghostFrom
      integer,dimension(:,:), allocatable :: InighostFrom
      integer,dimension(:,:), allocatable :: FinghostFrom
      integer,dimension(:,:), allocatable :: BlockBdyghostIndex
      integer,dimension(:),   allocatable :: BlockBdyghostNumber
      integer,dimension(:),   allocatable :: InighostIndex_local
      integer,dimension(:), allocatable :: type_send_ghost
      integer,dimension(:), allocatable :: type_recv_ghost
#     endif

      CONTAINS
!      _______________________________________________________________
!     |      |                                                        |
!     | 1.1  |                 Initialization MPI                     |
!     |______|________________________________________________________|

      SUBROUTINE initialisation_mpi

      !--------------------------------------------------------------!
      !                                                              !
      !   ---> code      : name of the MPI program                   !
      !   <--- Nprocs    : total number of processors                !
      !   <--- rang_topo : rank (name) of each processor             !
      !                                                              !
      !--------------------------------------------------------------!

!     __________________________
!     Initialisation de MPI
      call MPI_INIT(code)
!     __________________________
!     Number of processors
      call MPI_COMM_SIZE(MPI_COMM_WORLD,Nprocs,code)
!     __________________________
!     Rang of each processor
      call MPI_COMM_RANK(MPI_COMM_WORLD,rang_topo,code)

!      ________________________________________________________
!     |                                                        |
!     |          Elements in the vertical direction: NZ        |
!     |________________________________________________________|
       NZ = NZglobal
!     __________________________
!     Display global parameters
        IF (rang_topo.eq.0) THEN
            print*,'                                                         '
            print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
            print*,'========================================================='
            print*,'           __________________________________            '
            print*,'          |                                  |           '
            print*,'          |  ***  MPI PROGRAM: NSMP3D   ***  |           '
            print*,'          |__________________________________|           '
            print*,'                                                         '
            print*,'========================================================='
            print*,'wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww'
            print*,'                                                         '
            print*,'   ====================================================  '
            print*,'                   MPI: INITIALIZATION                   '
            print*,'   ====================================================  '
            print*,'                                                         '
            print*,'       NUMBER OF MPI PROCESSORS: ', Nprocs
            print*,'                                                         '
            print*,'   ____________________________________________________  '
            print*,'                                                         '
            print*,'       GLOBAL DOMAIN:                                    '
            print*,'          Global number of elements  =',N_CELL0Global
            print*,'          Global number of vertices  =',N_VERTGlobal
            print*,'          Global number of NZ points =',NZglobal
            print*,'          Local  number of NZ points =',NZ
            print*,'                                                         '
        ENDIF

      END SUBROUTINE initialisation_mpi
!      _______________________________________________________________
!     |       |                                                        |
!     | 1.1a  |                 parallel_ghost                         |
!     |______ |________________________________________________________|

       SUBROUTINE parallel_ghost

        integer :: ielem,ii,nc,cc,kk,jj,jv1,jv2,jv3,i0
        integer :: jv10,jv20,jv30,ncount,vcount,jv1p,jv2p
        real*8 :: xmin,ymin,xmax,ymax,dxl,dyl
        real*8 :: xee,yee,dxody,dxody2
        real*8 :: xcc,ycc
!    ---------------------------------------
!     Global ghost cell count
!     To use in save tec center(parallel mode)
                ii = 0
                DO i= 1,N_CELL0global
                    ielem = i
                    if (nbe_global(ielem).gt.0) then
                        do j=1,3
                            nc = No_cp_global(ielem,j)
                            if ((nc.lt.1).or.(nc.gt.N_CELL0global)) then
                                ii=ii+1
                            endif
                        enddo
                    endif
                ENDDO

            Nglobalghost = ii
!     -------------------------------------
!     Obtain global vertex bc in case of no periodic bc
#     ifndef KeyTESTpBC
            do i= 1,N_CELL0global
                if (nbe_global(i) .gt. 0) then
                    do j=1,3
                        nc = No_cp_global(i,j)
                        if ((nc.lt.1).or.(nc.gt.N_CELL0global)) then
                            jj=j+1
                            if (jj.gt.3) jj=jj-3
                            jv1 = No_vp_global(i,j)
                            jv2 = No_vp_global(i,jj)
                            nbev_global(jv1) = 1    ! set global boundary index
                            nbev_global(jv2) = 1    ! set global boundary index
                        endif
                    enddo
                endif
            enddo
#     else
!     -------------------------------------
!     To use in Newmann BC
             DO nv=1,N_VERTglobal
              vpair_global(nv) = 0
              vedge_global(nv) = 0
             ENDDO
!       --------------------------------------------------------------
             do i=1,N_CELL0global
                jv1=No_vp_global(i,1)
                jv2=No_vp_global(i,2)
                jv3=No_vp_global(i,3)
                xc_global(i)=(xv_global(jv1)+xv_global(jv2)+xv_global(jv3))/3.0d0
                yc_global(i)=(yv_global(jv1)+yv_global(jv2)+yv_global(jv3))/3.0d0
              enddo
!     ----------------------------------------------
            allocate(ghostCCdom(1:Nglobalghost,5))

            ii = 0
            vcount = 0
            do i= 1,N_CELL0global
                ielem = i
                if (nbe_global(i).gt.0) then
                    do j=1,3
                        nc = No_cp_global(ielem,j)
                        if ((nc.lt.1).or.(nc.gt.N_CELL0global)) then
                        ii=ii+1
                        ghostCCdom(ii,1) = i !------> ii is from index i
                        ghostCCdom(ii,2) = j  !------> ii is the j neibouger of nc
                        ghostCCdom(ii,3) = i
!                 -------------------------------------------
!                 Intersection on the triangle boundary edge
                        jj=j+1
                        if (jj.gt.3) jj=jj-3
                        jv1 = No_vp_global(i,j)
                        jv2 = No_vp_global(i,jj)
                        nbev_global(jv1) = 1    ! set global boundary index
                        nbev_global(jv2) = 1    ! set global boundary index
                        if (dabs(yv_global(jv2)-yv_global(jv1)).lt.1.0E-7) then
                            yee = yv_global(jv2)
                            xee = xc_global(i)
                        else
                            dxody  = xv_global(jv2)-xv_global(jv1)
                            dxody  = dxody/(yv_global(jv2)-yv_global(jv1))
                            dxody2 = dxody*dxody
                            yee= (yc_global(i)-(xv_global(jv2)-xc_global(i))*dxody &
                                +yv_global(jv2)*dxody2)/(1.0d0+dxody2)
                            xee= xv_global(jv2)+(yee-yv_global(jv2))*dxody
                        endif
                        xcc=2.0d0*xee-xc_global(i)
                        ycc=2.0d0*yee-yc_global(i)
!          -------------------------------------------
!                      Periodic BC Implementation
                       xmin = XDIni
                       xmax = XDFin
                       ymin = YDIni
                       ymax = YDFin
                       dxl = xmax - xmin
                       dyl = ymax - ymin
                       jv1p = 0
                       jv2p = 0
                       cc = 0
!           --------------------------------------------
                        IF(XPB .EQ. 1) THEN
                            if (xcc.gt. xmax) then
                                dxl = -1.0d0*dxl
                                dyl = 0.
                            elseif(xcc .lt. xmin) then
                                dxl = dxl
                                dyl = 0.
                            endif
                        ENDIF
!          ----------------------------------------------
                        IF(YPB .EQ. 1) THEN
                            if (ycc .gt. ymax) then
                                dxl = 0.
                                dyl = -1.0d0*dyl
                            elseif(ycc .lt. ymin) then
                                dxl = 0.
                                dyl = dyl
                            endif
                        ENDIF
!          ================================================
                        IF(dxl*dyl .gt. 1.0E-7) THEN
                           vcount = vcount + 1
                        ELSE
!           -----------------------------------------------
!            Pair for jv1
                            ncount = 0
                            xee = xv_global(jv1) + dxl
                            yee = yv_global(jv1) + dyl
                            do nv =1,N_VERTglobal
                                if(dabs(xee-xv_global(nv)) .lt. 1.0E-7) then
                                    if(dabs(yee-yv_global(nv)) .lt. 1.0E-7) then
                                        jv1p = nv
                                        ncount = ncount +1
                                    endif
                                endif
                            enddo
                            if(ncount .ne. 1) print*, 'ERROR!VERTEX PAIR FOR JV1'

!           -----------------------------------------------
!            Pair for jv2
                            ncount = 0
                            xee = xv_global(jv2) + dxl
                            yee = yv_global(jv2) + dyl
                            do nv =1,N_VERTglobal
                                if(dabs(xee-xv_global(nv)) .lt. 1.0E-7) then
                                    if(dabs(yee-yv_global(nv)) .lt. 1.0E-7) then
                                        jv2p = nv
                                        ncount = ncount+1
                                    endif
                                endif
                            enddo
                            if(ncount .ne. 1) print*, 'ERROR!VERTEX PAIR FOR JV2'
!           -----------------------------------------------
                            ncount = 0
                            do i0 =1, N_CELL0global
                                if(nbe_global(i0) .ne. 0) then

                                    jv10 = No_vp_global(i0,1)
                                    jv20 = No_vp_global(i0,2)
                                    jv30 = No_vp_global(i0,3)

                                    if(jv10*jv20 .eq. jv1p*jv2p) then
                                        if(jv10 .eq. jv1p) then
                                            cc = jv30
                                            ghostCCdom(ii,3) = i0
                                            ncount = ncount+1
                                        elseif(jv10 .eq. jv2p) then
                                            cc = jv30
                                            ghostCCdom(ii,3) = i0
                                            ncount = ncount+1
                                        endif
                                    endif

                                    if(jv30*jv20 .eq. jv1p*jv2p) then
                                        if(jv30 .eq. jv1p) then
                                            cc = jv10
                                            ghostCCdom(ii,3) = i0
                                            ncount = ncount+1
                                        elseif(jv30 .eq. jv2p) then
                                            cc = jv10
                                            ghostCCdom(ii,3) = i0
                                            ncount = ncount+1
                                        endif
                                    endif

                                    if(jv10*jv30 .eq. jv1p*jv2p) then
                                        if(jv10 .eq. jv1p) then
                                            cc = jv20
                                            ghostCCdom(ii,3) = i0
                                            ncount = ncount+1
                                        elseif(jv10 .eq. jv2p) then
                                            cc = jv20
                                            ghostCCdom(ii,3) = i0
                                            ncount = ncount+1
                                        endif
                                    endif
                                endif
                            enddo
!                   ------------------------------------------------------
                       if (ncount .ne. 1) print*, 'ERROR FIND PERIDOIC CELL'
!                   ------------------------------------------------------
                   ENDIF
               endif
            enddo
         endif
       enddo
#       endif
!*********************************************************************!
!                                                                     !
!                       Assign TagXV,TagYV info                       !
!                                                                     !
!*********************************************************************!
            cc = 0
            DO nv=1,N_VERTglobal
                IF (nbev_global(nv).ne. 0) THEN
                    ncount = 0
                    xee = xv_global(nv)
                    yee = yv_global(nv)
                    if(dabs(xee - XDIni) .lt. 1.0d-7) then
                        tagXV_global(nv) = 1
                        ncount = ncount +1
                    endif

                    if(dabs(xee - XDFin) .lt. 1.0d-7) then
                        tagXV_global(nv) = 2
                        ncount = ncount +1
                    endif

                    if(dabs(yee - YDIni) .lt. 1.0d-7) then
                        tagYV_global(nv) = 1
                        ncount = ncount +1
                    endif

                    if(dabs(yee - YDFin) .lt. 1.0d-7) then
                        tagYV_global(nv) = 2
                        ncount = ncount +1
                    endif
!         _______________________
!         periodic pair for the vertex
#         ifdef KeyTESTpBC
!           ____________________________________
!            Boundary Point
                    if(ncount .eq. 1) then
                        do i=1,N_VERTglobal
                            if (nbev_global(i) .ne. 0) then
!                    _________________________________________________
!                       periodic pair in y direction
                        if (YPB .eq. 1) then
                            if(dabs(xee - xv_global(i)) .lt. 1.0d-7) then
                                if(dabs(yee - yv_global(i)) .gt. (YDFin-YDIni-1.0d-7)) then
                                    vpair_global(nv) = i
                                endif
                            endif
                        endif
!                   _____________________________________________________
!                       periodic pair in x direction
                        if (XPB .eq. 1) then
                            if(dabs(yee - yv_global(i)) .lt. 1.0d-7) then
                                if(dabs(xee - xv_global(i)) .gt. (XDFin-XDIni-1.0d-7)) then
                                    vpair_global(nv) = i
                                endif
                            endif
                        endif
                    endif
                enddo
            endif ! end ncount = 1
!           _______________________________________
!           Corner Point
            if(ncount .eq. 2) then
                if(XPB*YPB .eq. 1) then
                    vedge_global(nv) = 1
                else
                    vedge_global(nv) = 0
                    do i=1,N_VERTglobal
                        if (nbev_global(i) .ne. 0) then
!                    _________________________________________________
!                       periodic pair in y direction
                     if (YPB .eq. 1) then
                        if(dabs(xee - xv_global(i)) .lt. 1.0d-7) then
                            if(abs(yee - yv_global(i)) .gt. (YDFin-YDIni-1.0d-7)) then
                             vpair_global(nv) = i
                            endif
                        endif
                     endif
!                   _____________________________________________________
!                       periodic pair in x direction
                        if (XPB .eq. 1) then
                            if(dabs(yee - yv_global(i)) .lt. 1.0d-7) then
                                if(dabs(xee - xv_global(i)) .gt. (XDFin-XDIni-1.0d-7)) then
                                    vpair_global(nv) = i
                                endif
                            endif
                        endif
                    endif
                    enddo
                endif ! end XPB*YPB = 1
            endif ! end ncount = 2
!           ___________________
#           endif
            if (ncount .eq. 0) cc = cc + 1
            if (ncount .gt. 2) print*, 'ERROR! Global Vertex Coordinate!',ncount
           ENDIF ! end nbev .ne. 0
       ENDDO

!      ----------------------------------------------
         if(rang_topo .eq. 0) then
            print*, '        ======================================='
#        ifdef KeyTESTpBC
            print*, '                     WALL CELL :', vcount
#        endif
            print*, '        ======================================='
          endif
!*********************************************************************!
!                                                                     !
!                       End TagXV,TagYV info                          !
!                                                                     !
!*********************************************************************!

       call MPI_Barrier(comm3D,code)

       END SUBROUTINE parallel_ghost

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.2  |               Processors topology                      |
!     |______|________________________________________________________|

      SUBROUTINE parallel_topology

      !---------------------------------------------------------------!
      !                                                               !
      !  Division of the domain in subdomains by a graph topology.    !
      !                                                               !
      ! ---> MPI_COMM_WORLD : the communicator group we are using.    !
      ! ---> nbprocs: the number of processors.                       !
      ! ---> index  : array of integers describing node degrees       !
      ! ---> edges  : array of integers describing graph edges        !
      ! --->   0    : if we don't want to order processes in the group!
      ! <--- comm3D : the communicator which represents the graph     !
      !                                                               !
      !---------------------------------------------------------------!
      call MPI_GRAPH_CREATE(MPI_COMM_WORLD,Nprocs,indexval,edges,0,comm3D,code)
!     __________________________
!     Display topology
     !IF (rang_topo.eq.1) THEN
     !print*,'   ____________________________________________________  '
     !print*,'                                                         '
     !print*,'       MPI GRAPH TOPOLOGY:                               '
     !print*,'          Number of indexes =',Nindex
     !print*,'          Number of edges   =',Nedges
     !print*,'          Index values      =',indexval(0:Nindex-1)
     !print*,'          Edge values       =',edges(0:Nedges-1)
     !print*,'                                                        '
     !ENDIF

      END SUBROUTINE parallel_topology

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.3  |                    Neighbors                           |
!     |______|________________________________________________________|

       SUBROUTINE parallel_neighbors

       integer, dimension(:),allocatable :: tag
       integer :: ii,ielem,kelem

      !---------------------------------------------------------------!
      !                                                               !
      !  Calculation of the neighbors:                                !
      !  <--- neighbourNumber : is number of neighbors of the "s" proc!
      !  <--- neighbour       : is the array neighbors of "s" proc.   !
      !                                                               !
      !---------------------------------------------------------------!

!     _________________________________________________________________
!     Neighbors calculation
!     --------------------------------
!     Number of neighborns
      call MPI_Graph_neighbors_count(comm3D,rang_topo,NeighborNumber,code)
!     --------------------------------
!     Array of neighborns
      allocate(Neighbors(NeighborNumber))
      call MPI_Graph_neighbors(comm3D,rang_topo,NeighborNumber,Neighbors,code)
!     ----------------------------------
!     Select MPI communication option
!      MPIOption = 1 ! first send then recv
       MPIOption = 2 ! Non-block send and recv
       allocate(index_neigbour(NeighborNumber))

        if(MPIOption .eq. 2) then

            allocate(tag(NeighborNumber))

            do k=1,NeighborNumber
                tag(k) = Neighbors(k)
            enddo

            do k=1,NeighborNumber-1
                do ii=k+1,NeighborNumber
                    ielem= tag(k)
                    kelem= tag(ii)
                    if(ielem .gt. kelem) then
                        tag(k) = kelem
                        tag(ii) = ielem
                    endif
                enddo
            enddo

            do k=1,NeighborNumber
               do ii= 1, NeighborNumber
                 if(tag(k) .eq. Neighbors(ii)) then
                    index_neigbour(k) = ii
                    goto 100
                 endif
                enddo
100            continue
            enddo

            deallocate(tag)
        endif

      END SUBROUTINE parallel_neighbors
!      _______________________________________________________________
!     |      |                                                        |
!     | 1.4  |           Assignation of global index                  |
!     |______|________________________________________________________|

      SUBROUTINE parallel_index

!---------------------------------------------------------------------!
!                                                                     !
!     Global index values corresponding to local index values for     !
!     each subdomain.                                                 !
!                                                                     !
!     IMPORTANT: It is really important to remark the order of the    !
!                overlaping cell-centers are located. I will include  !
!                the ghost points and later the extra points.         !
!                                                                     !
!---------------------------------------------------------------------!

      integer,dimension(Nprocs) :: NNeigh_vec,N_CELL0_vec,N_VERT_vec
      integer :: Iini,Ifin,i0,s
      integer :: k,ii,nv

      s = sNZ

      allocate(Index_globalv(N_VERT))
      allocate(Index_global(N_CELL))
!     ________________________________________________________
!     Global index of vertex points in each subdomain

      nv = 0
      Iini = DimVVdom(s,2)
      Ifin = DimVVdom(s,3)
      do i=Iini,Ifin
         nv = nv + 1
         Index_globalv(nv) = VVdom(i,1)
      enddo
!      ________________________________________________________
!     Global index of cell-center points in each subdomain
!     -----------------
!     Interior
      i0 = 0
      Iini = DimCCdom(s,2)
      Ifin = DimCCdom(s,3)
      do i=Iini,Ifin
         i0 = i0 + 1
         Index_global(i0) = CCdom(i,1)
      enddo
!     -----------------
!     Ghost (No global index needed, temporal = -10)
      do i=1,N_CELLghost
         i0 = i0 + 1
         Index_global(i0) = -10
      enddo
!     -----------------
!     extra overlaping
      Iini = DimExtraCCdom(s,2)
      Ifin = DimExtraCCdom(s,3)
      do i=Iini,Ifin
         i0 = i0 + 1
         Index_global(i0) = ExtraCCdom(i,3)
      enddo
!     _______________________________________________________
!     Display subdomain information
!
!     IF (rang_topo.eq.0) THEN
!        print*,'   ____________________________________________________  '
!        print*,'                                                         '
!        print*,'       NEIGHBORS:                                        '
!        print*,'                                                         '
!        print*,'    p  | N_VERT | N_CELL0| N_CELLghost | N_CELL | No. Neigh & Neighbors'
!        print*,'   ----|--------|--------|-------------|--------|-----------------------'
!     ENDIF
      call MPI_Barrier(comm3D,code)
!     write(*,9),rang_topo,' |',N_VERT,' |',N_CELL0,' |',N_CELLghost,' |',N_CELL, ' |',&
!                NeighborNumber,': ',Neighbors(1:NeighborNumber)
!     9 format(I6,a3,I6,a3,I6,a3,I6,a3,I6,a3,I2,a3,10I4)


      END SUBROUTINE parallel_index

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.5  |              Index to communication                    |
!     |______|________________________________________________________|

      SUBROUTINE parallel_shareindex

      !---------------------------------------------------------------!
      !                                                               !
      !       Assign the vertex and cell-centers to each domain       !
      !                                                               !
      !   <--- IniExtraIndex_local: Initial local index to received   !
      !   <--- BlockBdyIndex      : Index of bdy. elem. for neighbor  !
      !   <--- BlockBdyNumber     : Number of bdy. elem for neighborn !
      !   <--- BlockIntIndex      : Index interior elem. in the block !
      !   <--- BlockIntNumber     : Number interior elem. in the block!
      !   <--- SharePhi           : Just allocate this vector         !
      !                                                               !
      !---------------------------------------------------------------!

      integer :: Iini,Ifin,i0,k0,ielem,jelem,kelem,ivert,s,sm,suma,m
      integer :: k,ii,jj,tag
      integer :: nv1,nv2,nv3
      character*80 filen


!     ________________________________________________________
!     Initial index of extra cells by neighborn

      allocate(IniExtraIndex_local(NeighborNumber))

      s  = rang_topo + 1
!     _________________________________
!     Share the extra cell info
      suma = N_CELL0 + N_CELLghost + 1

      do i=1,NeighborNumber
          sm = Neighbors(i) + 1
          IniExtraIndex_local(i) = suma
          suma = suma + DimExtraFrom(s,sm)
      enddo

!     ________________________________________________________
!     Share index to send
      NextraMax = 0
      do s=1,TotalSubdomains
         NextraMax = max(NextraMax,DimExtraCCdom(s,1))
      enddo

      allocate(BlockBdyIndex(NextraMax,NeighborNumber))
      allocate(BlockBdyNumber(NeighborNumber))

      s = rang_topo + 1

      do k=1,NeighborNumber
         sm = Neighbors(k) + 1
         Iini = IniExtraFrom(sm,s)
         Ifin = FinExtraFrom(sm,s)
         i0 = 0
         do j=Iini,Ifin
            i0 = i0 + 1
            jelem = ExtraCCdom(j,3)
            do ii=1,N_CELL0
               ielem = Index_global(ii)
               if (ielem.eq.jelem) then
                   BlockBdyIndex(i0,k) = ii
               endif
            enddo
         enddo
         BlockBdyNumber(k) = i0
      enddo

!     ---------------------------------------------------------------------
#     ifdef KeyDbgPBC
            filen='../output/PBC/extra_send-  .txt'
            write(filen(26:27),'(i2.2)') rang_topo
            open(210,file=filen)
            write(210,*),'Proc =',s
            do k=1,NeighborNumber
                sm = Neighbors(k) + 1
                i0 = BlockBdyNumber(k)
                write(210,*), 'Neighbor proc',sm,'value send count=',i0,DimExtraFrom(sm,s)
                do ii=1,i0
                    ielem = BlockBdyIndex(ii,k)
                    jelem = Index_global(ielem)
                    write(210,*), ii, BlockBdyIndex(ii,k),jelem,TagParallel(jelem)+1
                enddo
                write(210,*), 'Neighbor proc',sm, 'value recv count =',DimExtraFrom(s,sm),&
                    'start from:',IniExtraIndex_local(k)
                do i=1,DimExtraFrom(s,sm)
                    ielem = i + IniExtraIndex_local(k) -1
                    jelem = Index_global(ielem)
                    write(210,*), i,ielem,jelem,TagParallel(jelem)+1
                enddo
            enddo
            close(210)
#     endif
!     ________________________
!     Share the ghost cell
#     ifdef KeyTESTpBC
!     ________________________________________________________
!     Initial index of ghost cells by neighborn
      allocate(InighostIndex_local(NeighborNumber))
!     _________________________________
!     Share the extra cell info
      suma = N_CELL0
      do i=1,NeighborNumber
          sm = Neighbors(i) + 1
          InighostIndex_local(i) = suma + InighostFrom(s,sm)
      enddo

      NextraMax = 0
      do s=1,TotalSubdomains
         NextraMax = max(NextraMax,DimghostCCdom(s,1))
      enddo

      allocate(BlockBdyghostIndex(NextraMax,NeighborNumber))
      allocate(BlockBdyghostNumber(NeighborNumber))

      s = rang_topo + 1
      do k=1,NeighborNumber
         sm = Neighbors(k) + 1
         i0 = 0
!       ---------------------------------
!       For ghost cell
         Iini = DimghostCCdom(sm,2)
         Ifin = DimghostCCdom(sm,3)
         do j=Iini,Ifin
          if(ghostCCdom(j,5) .eq. s) then
            i0 = i0 + 1
            kelem = ghostCCdom(j,3)
            jelem = CCdom(kelem,1)
            if(CCdom(kelem,2) .ne.s) then
                print*, 'FATAL ERROR! GHOST NOT MATCH'
            endif
            do ii=1,N_CELL0
               ielem = Index_global(ii)
               if (ielem.eq.jelem) then
                   BlockBdyghostIndex(i0,k) = ii
               endif
            enddo
           endif
         enddo
!       ---------------------------------
         BlockBdyghostNumber(k) = i0
      enddo
!     ---------------------------------------------------------------------
!#     ifdef KeyDbgPBC
!            filen='../output/PBC/ghost_send-  .txt'
!            write(filen(26:27),'(i2.2)') rang_topo
!            open(160,file=filen)
!            write(160,*),'Proc =',s
!            do k=1,NeighborNumber
!                sm = Neighbors(k) + 1
!                i0 = BlockBdyghostNumber(k)
!                write(160,*), 'Neighbor proc',sm,'value send count=',i0,DimghostFrom(sm,s)
!                do ii=1,i0
!                    ielem = BlockBdyghostIndex(ii,k)
!                    jelem = Index_global(ielem)
!                    write(160,*), ii, ielem,jelem,TagParallel(jelem)+1
!                enddo
!                write(160,*), 'Neighbor proc',sm, 'value recv count =',DimghostFrom(s,sm),&
!                    'start from:',InighostIndex_local(k)
!                do ii=1,DimghostFrom(s,sm)
!                    ielem = ii -1 + InighostIndex_local(k) - N_CELL0 + DimghostCCdom(s,2) -1
!                    jelem = ghostCCdom(ielem,3)
!                    kelem = CCdom(jelem,1)
!                    write(160,*),ii,ielem,kelem,TagParallel(kelem)+1
!                enddo
!            enddo
!            close(160)
!#     endif
#     endif

!     -----------------------------------------
!     Virtual boundary index of the neardest cell-center values
!      NextraMax = 0
!      do s=1,TotalSubdomains
!         NextraMax = max(NextraMax,DimExtraCCdom(s,1))
!      enddo
!      allocate(BlockBdyCCIndex(NextraMax,NeighborNumber))
!      allocate(BlockBdyCCNumber(NeighborNumber))
!
!      s  = rang_topo + 1
!      do k=1,NeighborNumber
!         sm = Neighbors(k) + 1
!         Iini = DimBdyCCdom(s,2)
!         Ifin = DimBdyCCdom(s,3)
!         i0 = 0
!         do j=Iini,Ifin
!            if (BdyCCdom(s,3).eq.sm) then
!               i0 = i0 + 1
!               BlockBdyCCIndex(i0,k)= BdyCCdom(s,1)
!            endif
!         enddo
!         BlockBdyCCNumber(k) = i0
!      enddo
!     ________________________________________________________
!     Free memory ExtraCCdom global
      deallocate(ExtraCCdom)
      deallocate(DimExtraCCdom)

      END SUBROUTINE parallel_shareindex
!      _______________________________________________________________
!     |      |                                                        |
!     | 1.5  |              Index to communication                    |
!     |______|________________________________________________________|

      SUBROUTINE parallel_shareindex_NZ

      !---------------------------------------------------------------!
      !                                                               !
      !       Assign the vertex and cell-centers to each domain       !
      !                                                               !
      !   <--- IniExtraIndex_local: Initial local index to received   !
      !   <--- BlockBdyIndex      : Index of bdy. elem. for neighbor  !
      !   <--- BlockBdyNumber     : Number of bdy. elem for neighborn !
      !   <--- BlockIntIndex      : Index interior elem. in the block !
      !   <--- BlockIntNumber     : Number interior elem. in the block!
      !   <--- SharePhi           : Just allocate this vector         !
      !                                                               !
      !---------------------------------------------------------------!

      integer :: Iini,Ifin,i0,k0,ielem,jelem,kelem,ivert,sm,suma,m
      integer :: k,ii,jj,tag,smNZ,smLayer
      integer :: nv1,nv2,nv3
      character*80 filen


!     ________________________________________________________
!     Initial index of extra cells by neighborn

      allocate(IniExtraIndex_local(NeighborNumber))
!     _________________________________
!     Share the extra cell info
        suma = N_CELL0 + N_CELLghost + 1

        do i=1,NeighborNumber
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer) then
                IniExtraIndex_local(i) = suma
                suma = suma + DimExtraFrom(sNZ,smNZ)
            elseif(smLayer .gt. NZLayer) then
                if((sm-rang_topo-1) .eq. TotalSubdomains) then
                    IniExtraIndex_local(i) = NZ-1 ! send phi(N_CELL,NZ-1)
                else
                    IniExtraIndex_local(i) = 2  ! send phi(N_CELL,2)
                endif
            elseif(smLayer .lt. NZLayer) then
                if((rang_topo-sm+1) .eq. TotalSubdomains) then
                    IniExtraIndex_local(i) = 2 ! send phi(N_CELL,2)
                else
                    IniExtraIndex_local(i) = NZ-1 ! send phi(N_CELL,NZ-1)
                endif
            endif
        enddo
!     ________________________________________________________
!     Share index to send
      NextraMax = N_CELL

      allocate(BlockBdyIndex(NextraMax,NeighborNumber))
      allocate(BlockBdyNumber(NeighborNumber))

        do k=1,NeighborNumber
            sm = Neighbors(k) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer) then
                Iini = IniExtraFrom(smNZ,sNZ)
                Ifin = FinExtraFrom(smNZ,sNZ)
                i0 = 0
                do j=Iini,Ifin
                    i0 = i0 + 1
                    jelem = ExtraCCdom(j,3)
                    do ii=1,N_CELL0
                        ielem = Index_global(ii)
                        if (ielem.eq.jelem) then
                            BlockBdyIndex(i0,k) = ii
                        endif
                    enddo
                enddo
                BlockBdyNumber(k) = i0
            else
                BlockBdyNumber(k) = N_CELL
                do j=1,N_CELL
                    BlockBdyIndex(j,k) = j
                enddo
            endif
        enddo

#     ifdef KeyTESTpBC
!     ________________________
!     Share the ghost cell
!     ________________________________________________________
!     Initial index of ghost cells by neighborn
      allocate(InighostIndex_local(NeighborNumber))
!     _________________________________
!     Share the extra cell info
      suma = N_CELL0
        do i=1,NeighborNumber
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer) then
                InighostIndex_local(i) = suma + InighostFrom(sNZ,smNZ)
            else
                InighostIndex_local(i) = 0
            endif
        enddo

      NextraMax = 0
      do i=1,TotalSubdomains
         NextraMax = max(NextraMax,DimghostCCdom(i,1))
      enddo

      allocate(BlockBdyghostIndex(NextraMax,NeighborNumber))
      allocate(BlockBdyghostNumber(NeighborNumber))


        do k=1,NeighborNumber
            sm = Neighbors(k) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer) then
                i0 = 0
                Iini = DimghostCCdom(smNZ,2)
                Ifin = DimghostCCdom(smNZ,3)
                do j=Iini,Ifin
                    if(ghostCCdom(j,5) .eq. sNZ) then
                        i0 = i0 + 1
                        kelem = ghostCCdom(j,3)
                        jelem = CCdom(kelem,1)
                        if(CCdom(kelem,2) .ne. sNZ) then
                            print*, 'FATAL ERROR! GHOST NOT MATCH'
                        endif
                        do ii=1,N_CELL0
                            ielem = Index_global(ii)
                            if (ielem.eq.jelem) then
                                BlockBdyghostIndex(i0,k) = ii
                            endif
                        enddo
                    endif
                enddo
                BlockBdyghostNumber(k) = i0
            else
                BlockBdyghostNumber(k) = 0
            endif
        enddo
#     endif

!     =========================================================
!     THE FOLLOWING HAS NO USE IN THE CURRENT VERSION.
!     THE ORIGINAL FORMAT IS MODIFED ACCORDING TO NZ PARTITION.
!     BUT NOT USED.
!     =========================================================
!     -----------------------------------------
!     Virtual boundary index of the neardest cell-center values
!      NextraMax = 0
!      do s=1,TotalSubdomains
!         NextraMax = max(NextraMax,DimExtraCCdom(s,1))
!      enddo
!      allocate(BlockBdyCCIndex(NextraMax,NeighborNumber))
!      allocate(BlockBdyCCNumber(NeighborNumber))
!
!      s  = rang_topo + 1
!        do k=1,NeighborNumber
!            sm = Neighbors(k) + 1
!            smNZ = mod(sm-1,TotalSubdomains) + 1
!            smLayer = (sm-1)/TotalSubdomains + 1
!            if(smLayer .eq. NZLayer) then
!                Iini = DimBdyCCdom(sNZ,2)
!                Ifin = DimBdyCCdom(sNZ,3)
!                i0 = 0
!                do j=Iini,Ifin
!                    if (BdyCCdom(sNZ,3).eq.smNZ) then
!                        i0 = i0 + 1
!                        BlockBdyCCIndex(i0,k)= BdyCCdom(sNZ,1)
!                    endif
!                enddo
!                BlockBdyCCNumber(k) = i0
!            else
!                BlockBdyCCNumber(k) = 0
!        enddo
!     ________________________________________________________
!     Free memory ExtraCCdom global

      deallocate(ExtraCCdom)
      deallocate(DimExtraCCdom)

      END SUBROUTINE parallel_shareindex_NZ

!      _______________________________________________________________
!     |      |                                                        |
!     | 1.6  |                 Finalization MPI                       |
!     |______|________________________________________________________|

      SUBROUTINE finalisation_mpi

      !---------------------------------------------------------------!
      !                                                               !
      !   Finalization of definitions and the MPI version             !
      !                                                               !
      !---------------------------------------------------------------!

!     ________________________________________________________
!     Free memory
!     __________________
!     Index and edges
      deallocate(TagParallel)
      deallocate(indexval,edges)
      deallocate(Neighbors)
!     __________________
!     Global variables
!     ---------------------
      deallocate(VVdom)
      deallocate(DimVVdom)
      deallocate(BdyCCdom)
      deallocate(DimBdyCCdom)
      deallocate(CCdom)
      deallocate(DimCCdom)
!     ---------------------
      deallocate(IniExtraFrom)
      deallocate(FinExtraFrom)
      deallocate(DimExtraFrom)
!     ---------------------
      deallocate(No_vp_global,  &
                 No_cp_global,  &
                 nbe_global,    &
                 xv_global,     &
                 yv_global,     &
                 zbv_global)
!    -----------------------
      deallocate(nbev_global,  &
                 vcount_global,&
                 TagXV_global, &
                 TagYV_global)
!      -----------------------
!      deallocate local variable
      deallocate(xc_global)
      deallocate(yc_global)
!     __________________
!     Local variables
      deallocate(Index_global)
      deallocate(Index_globalv)
      deallocate(IniExtraIndex_local)
      deallocate(BlockBdyIndex)
      deallocate(BlockBdyNumber)
!      deallocate(BlockBdyCCIndex)
!      deallocate(BlockBdyCCNumber)
      deallocate(InterCCdom)
      deallocate(DimInterCCdom)
      deallocate(BoundCCdom)
      deallocate(DimBoundCCdom)
      deallocate(SharePhi)
      deallocate(SharePhiNew)
      deallocate(index_neigbour)
!     __________________
!     Communication variables
      deallocate(type_recv)
      deallocate(type_send)
      deallocate(type_recvVEC)
      deallocate(type_sendVEC)
!     ____________________
!     For periodic bc
#     ifdef KeyTESTpBC
      deallocate(vpair_global, &
                 vedge_global)
      deallocate(DimghostFrom)
      deallocate(InighostFrom)
      deallocate(FinghostFrom)
      deallocate(ghostCCdom)
      deallocate(DimghostCCdom)
      deallocate(ghost_local)
      deallocate(InighostIndex_local)
      deallocate(BlockBdyghostIndex)
      deallocate(BlockBdyghostNumber)
      deallocate(type_recv_ghost)
      deallocate(type_send_ghost)
#     endif
!     ________________________________________________________
!     Free block vertex type
      call MPI_TYPE_FREE(type_blocV,code)
      call MPI_TYPE_FREE(type_blocC,code)
      call MPI_TYPE_FREE(type_bloc2D,code)
!     ________________________________________________________
!     Free code
      call MPI_FINALIZE(code)


      END SUBROUTINE finalisation_mpi

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   MODULE PARALLEL (COMMUNICATION)                   !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!      SUBROUTINES:    3.1)  parallel_type                            !
!                      3.2)  communication2D                          !
!                      3.3)  communication3Dtype1                     !
!                      3.4)  communication3D                          !
!---------------------------------------------------------------------!

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.1  |    Definition of the type of communication vectors     |
!     |______|________________________________________________________|

      SUBROUTINE parallel_type

      !---------------------------------------------------------------!
      !                                                               !
      ! Specifying the vector type:                                   !
      !                                                               !
      !   <--- type_send(i)   : vector of extra elements to send to   !
      !                         each neighborn i (continuous type).   !
      !   <--- type_recv(i)   : vector of extra elements to receive to!
      !                         each neighborn i (continuous type).   !
      !   <--- type_sendVEC(i): vector of extra elements to send to   !
      !                         each neighborn i (vector type).       !
      !   <--- type_recvVEC(i): vector of extra elements to receive to!
      !                         each neighborn i (vector type).       !
      !                                                               !
      !---------------------------------------------------------------!

      integer :: number_send,number_recv
      integer :: i,s,sm,number_BdyCC

      allocate(type_send(NeighborNumber))
      allocate(type_recv(NeighborNumber))
      allocate(type_sendVEC(NeighborNumber))
      allocate(type_recvVEC(NeighborNumber))
!      allocate(type_BdyCCVEC(NeighborNumber))
!     ________________________________________
!     Continuous type: size of sent array
      do i=1,NeighborNumber
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_send = DimExtraFrom(sm,s)
!        --------------------
!        Continuous type
         call MPI_TYPE_CONTIGUOUS(number_send,&
                                  MPI_DOUBLE_PRECISION,type_send(i),code)
         call MPI_TYPE_COMMIT(type_send(i),code)
!        --------------------
!        Vector type
         call MPI_TYPE_VECTOR(NZ-2,number_send,N_CELL,&
                              MPI_DOUBLE_PRECISION,type_sendVEC(i),code)
         call MPI_TYPE_COMMIT(type_sendVEC(i),code)
!        --------------------
      enddo
!     ________________________________________
!     Continuous type: size of received array
      do i=1,NeighborNumber
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_recv = DimExtraFrom(s,sm)
!        --------------------
!        Continuous type
         call MPI_TYPE_CONTIGUOUS(number_recv,&
              MPI_DOUBLE_PRECISION,type_recv(i),code)
         call MPI_TYPE_COMMIT(type_recv(i),code)
!        --------------------
!        Vector type
         call MPI_TYPE_VECTOR(NZ-2,number_recv,N_CELL,&
              MPI_DOUBLE_PRECISION,type_recvVEC(i),code)
         call MPI_TYPE_COMMIT(type_recvVEC(i),code)
      enddo
!    For Periodic BC
#    ifdef KeyTESTpBC
      allocate(type_send_ghost(NeighborNumber))
      allocate(type_recv_ghost(NeighborNumber))
!     ________________________________________
!     Continuous type: size of sent array
      do i=1,NeighborNumber
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_send = DimghostFrom(sm,s)
!        --------------------
!        Continuous type
         call MPI_TYPE_CONTIGUOUS(number_send,&
                                  MPI_DOUBLE_PRECISION,type_send_ghost(i),code)
         call MPI_TYPE_COMMIT(type_send_ghost(i),code)
      enddo
!     ________________________________________
!     Continuous type: size of received array
      do i=1,NeighborNumber
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_recv = DimghostFrom(s,sm)
!        --------------------
!        Continuous type
         call MPI_TYPE_CONTIGUOUS(number_recv,&
              MPI_DOUBLE_PRECISION,type_recv_ghost(i),code)
         call MPI_TYPE_COMMIT(type_recv_ghost(i),code)
      enddo

#    endif
!     ________________________________________
!     Continuous type: share bdy: Only Cell-center(Vector type)

!      do i=1,NeighborNumber
!         s  = rang_topo + 1
!         sm = Neighbors(i) + 1
!         number_BdyCC = 0
!         do j=DimBdyCCdom(s,2),DimBdyCCdom(s,3)
!            if (BdyCCdom(s,2).eq.sm) then
!               number_BdyCC = number_BdyCC + 1
!            endif
!         enddo
!         call MPI_TYPE_VECTOR(NZ-2,number_BdyCC,N_CELL,&
!              MPI_DOUBLE_PRECISION,type_BdyCCVEC(i),code)
!         call MPI_TYPE_COMMIT(type_BdyCCVEC(i),code)
!      enddo

!     ________________________________________
!     Vertex global block type
      call MPI_TYPE_VECTOR(NZ-1,N_VERTglobal,N_VERTglobal,&
           MPI_DOUBLE_PRECISION,type_blocV,code)
      call MPI_TYPE_COMMIT(type_blocV,code)

!     ________________________________________
!     Cell-center global continuous type
      call MPI_TYPE_VECTOR(NZ,N_CELL0global,N_CELL0global,&
                          MPI_DOUBLE_PRECISION,type_blocC,code)
      call MPI_TYPE_COMMIT(type_blocC,code)
!     ________________________________________
!     Cell-center global continuous type
      call MPI_TYPE_CONTIGUOUS(N_CELL0global,&
                            MPI_DOUBLE_PRECISION,type_bloc2D,code)
      call MPI_TYPE_COMMIT(type_bloc2D,code)
!     ________________________________________
!     Allocate auxiliar variable to communicate

      NshareMax = 0
      do i=1,NeighborNumber
         NshareMax = max(NshareMax,BlockBdyNumber(i))
      enddo
#     ifdef KeyTESTpBC
      do i=1,NeighborNumber
         NshareMax = max(NshareMax,BlockBdyghostNumber(i))
      enddo
#     endif
      allocate(SharePhi(NshareMax))
      allocate(SharePhiNew(N_CELL-2,NZ))

!     =============================================
!     ORGINIALLY WE HAVE TYPE2 COMMUNICATION, WHICH
!     DIVIDE THE NZ INTO TWO PARTS. WITH NZ PARTITION
!     THIS NO LONGER NEEDED!
!     =============================================

      END SUBROUTINE parallel_type
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.1  |    Definition of the type of communication vectors     |
!     |______|________________________________________________________|

      SUBROUTINE parallel_type_NZ

      !---------------------------------------------------------------!
      !                                                               !
      ! Specifying the vector type:                                   !
      !                                                               !
      !   <--- type_send(i)   : vector of extra elements to send to   !
      !                         each neighborn i (continuous type).   !
      !   <--- type_recv(i)   : vector of extra elements to receive to!
      !                         each neighborn i (continuous type).   !
      !   <--- type_sendVEC(i): vector of extra elements to send to   !
      !                         each neighborn i (vector type).       !
      !   <--- type_recvVEC(i): vector of extra elements to receive to!
      !                         each neighborn i (vector type).       !
      !                                                               !
      !---------------------------------------------------------------!

      integer :: number_send,number_recv
      integer :: i,s,sm,number_BdyCC
      integer :: smNZ,smLayer

      allocate(type_send(NeighborNumber))
      allocate(type_recv(NeighborNumber))
      allocate(type_sendVEC(NeighborNumber))
      allocate(type_recvVEC(NeighborNumber))
!      allocate(type_BdyCCVEC(NeighborNumber))
!     ________________________________________
!     Continuous type: size of sent array
        do i=1,NeighborNumber
            s  = rang_topo + 1
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer)  then
                number_send = DimExtraFrom(smNZ,sNZ)
            else
                number_send = N_CELL
            endif
!        --------------------
!        Continuous type
            call MPI_TYPE_CONTIGUOUS(number_send,&
                MPI_DOUBLE_PRECISION,type_send(i),code)
            call MPI_TYPE_COMMIT(type_send(i),code)
!        --------------------
!        Vector type
            call MPI_TYPE_VECTOR(NZ-2,number_send,N_CELL,&
                MPI_DOUBLE_PRECISION,type_sendVEC(i),code)
            call MPI_TYPE_COMMIT(type_sendVEC(i),code)

        enddo
!     ________________________________________
!     Continuous type: size of received array
        do i=1,NeighborNumber
            s  = rang_topo + 1
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer)  then
                number_recv = DimExtraFrom(sNZ,smNZ)
            else
                number_recv = N_CELL
            endif
!        --------------------
!        Continuous type
            call MPI_TYPE_CONTIGUOUS(number_recv,&
                MPI_DOUBLE_PRECISION,type_recv(i),code)
            call MPI_TYPE_COMMIT(type_recv(i),code)
!        --------------------
!        Vector type
            call MPI_TYPE_VECTOR(NZ-2,number_recv,N_CELL,&
                MPI_DOUBLE_PRECISION,type_recvVEC(i),code)
            call MPI_TYPE_COMMIT(type_recvVEC(i),code)
        enddo
!    ---------------
!    For Periodic BC
#    ifdef KeyTESTpBC
      allocate(type_send_ghost(NeighborNumber))
      allocate(type_recv_ghost(NeighborNumber))
!     ________________________________________
!     Continuous type: size of sent array
        do i=1,NeighborNumber
            s  = rang_topo + 1
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer)  then
                number_send = DimghostFrom(smNZ,sNZ)
            else
                number_send = 0
            endif
!        --------------------
!        Continuous type
            call MPI_TYPE_CONTIGUOUS(number_send,&
                MPI_DOUBLE_PRECISION,type_send_ghost(i),code)
            call MPI_TYPE_COMMIT(type_send_ghost(i),code)
        enddo
!     ________________________________________
!     Continuous type: size of received array
        do i=1,NeighborNumber
            s  = rang_topo + 1
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer)  then
                number_recv = DimghostFrom(sNZ,smNZ)
            else
                number_recv = 0
            endif
!        --------------------
!        Continuous type
            call MPI_TYPE_CONTIGUOUS(number_recv,&
                MPI_DOUBLE_PRECISION,type_recv_ghost(i),code)
            call MPI_TYPE_COMMIT(type_recv_ghost(i),code)
        enddo

#    endif
!     ________________________________________
!     Continuous type: share bdy: Only Cell-center(Vector type)
!        do i=1,NeighborNumber
!            s  = rang_topo + 1
!            sm = Neighbors(i) + 1
!            smNZ = mod(sm-1,TotalSubdomains) + 1
!            smLayer = (sm-1)/TotalSubdomains + 1
!            number_BdyCC = 0
!            if(smLayer .eq. NZLayer)  then
!                do j=DimBdyCCdom(sNZ,2),DimBdyCCdom(sNZ,3)
!                    if (BdyCCdom(sNZ,2).eq.smNZ) then
!                        number_BdyCC = number_BdyCC + 1
!                    endif
!                enddo
!            endif
!            call MPI_TYPE_VECTOR(NZ-2,number_BdyCC,N_CELL,&
!                 MPI_DOUBLE_PRECISION,type_BdyCCVEC(i),code)
!            call MPI_TYPE_COMMIT(type_BdyCCVEC(i),code)
!      enddo

!     ________________________________________
!     Vertex global block type
      call MPI_TYPE_VECTOR(NZ-1,N_VERTglobal,N_VERTglobal,     &
           MPI_DOUBLE_PRECISION,type_blocV,code)
      call MPI_TYPE_COMMIT(type_blocV,code)
!     ________________________________________
!     Cell-center global continuous type
      call MPI_TYPE_VECTOR(NZ,N_CELL0global,N_CELL0global,     &
                           MPI_DOUBLE_PRECISION,type_blocC,code)
      call MPI_TYPE_COMMIT(type_blocC,code)
!     ________________________________________
!     Cell-center global continuous type
      call MPI_TYPE_CONTIGUOUS(N_CELL0global,&
                            MPI_DOUBLE_PRECISION,type_bloc2D,code)
      call MPI_TYPE_COMMIT(type_bloc2D,code)
!     ________________________________________
!     Allocate auxiliar variable to communicate
      NshareMax = 0
      do i=1,NeighborNumber
         NshareMax = max(NshareMax,BlockBdyNumber(i))
      enddo

#     ifdef KeyTESTpBC
      do i=1,NeighborNumber
         NshareMax = max(NshareMax,BlockBdyghostNumber(i))
      enddo
#     endif

      allocate(SharePhi(NshareMax))
      allocate(SharePhiNew(N_CELL,NZ-2))
!     =============================================
!     ORGINIALLY WE HAVE TYPE2 COMMUNICATION, WHICH
!     DIVIDE THE NZ INTO TWO PARTS. WITH NZ PARTITION
!     THIS NO LONGER NEEDED!
!     =============================================

      END SUBROUTINE parallel_type_NZ

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.2  |        Communication continuos type 2D                 |
!     |______|________________________________________________________|

       SUBROUTINE communication2D(phi)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable ar with one entry phi(N_CELL).                    !
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,elem,Iini
      integer :: smNZ,smLayer,sm
!     ______________________________
!     Send information
        do i=1,NeighborNumber
!        -----------
                do j=1,BlockBdyNumber(i)
                    elem = BlockBdyIndex(j,i)
                    SharePhi(j) = phi(elem)
                enddo
!        -----------
                call MPI_SEND(SharePhi(1),1,type_send(i),Neighbors(i),&
                    etiquette,comm3D,code)
        enddo

      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
        do i=1,NeighborNumber
                Iini = IniExtraIndex_local(i)
                call MPI_RECV(phi(Iini),1,type_recv(i),Neighbors(i),&
                    etiquette,comm3D,statut,code)
        enddo

      END SUBROUTINE communication2D

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.3  |         Communication 3D: vector type  1               |
!     |______|________________________________________________________|

       SUBROUTINE communication3Dtype1(phi)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                !
      !    Difference with communication3D:                           !
      !    Use vector type instead of continous type.                 !
      !    No communication of ghost points.                          !
      !---------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini
      integer :: sm,smNZ,smLayer

!     ______________________________
!     Send information
      do i=1,NeighborNumber
!        -----------
         do j=1,BlockBdyNumber(i)
            elem = BlockBdyIndex(j,i)
            do k=2,NZ-1
               SharePhiNew(j,k-1) = phi(elem,k)
            enddo
         enddo
!        -----------
         call MPI_SEND(SharePhiNew(1,1),1,type_sendVEC(i),Neighbors(i),&
                       etiquette,comm3D,code)

      enddo

      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
      do i=1,NeighborNumber
         Iini = IniExtraIndex_local(i)
         call MPI_RECV(phi(Iini,2),1,type_recvVEC(i),Neighbors(i),&
                       etiquette,comm3D,statut,code)
      enddo

      END SUBROUTINE communication3Dtype1
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.6  |        Communication 3D: NZ continuos (type 0)         |
!     |______|________________________________________________________|

       SUBROUTINE communication3D(phi)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                !
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini,iproc
      integer :: s,sm,smNZ,smLayer

      DO k=1,NZ
!      ========================
!      ====Sequence SENDRECV===
!      ========================
       if(MPIOption .eq. 1) then
!     ______________________________
!     Send information
         do i=1,NeighborNumber
#           ifdef KeyParallelNZ
            sm = Neighbors(i) + 1
            smNZ = mod(sm-1,TotalSubdomains) + 1
            smLayer = (sm-1)/TotalSubdomains + 1
            if(smLayer .eq. NZLayer) then
#           endif
!           -----------
            do j=1,BlockBdyNumber(i)
               elem = BlockBdyIndex(j,i)
               SharePhi(j) = phi(elem,k)
            enddo
!           -----------
            call MPI_SEND(SharePhi(1),1,type_send(i),Neighbors(i),&
                          etiquette,comm3D,code)
#           ifdef KeyParallelNZ
            endif
#           endif
         enddo

      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
         do i=1,NeighborNumber
            Iini = IniExtraIndex_local(i)
            call MPI_RECV(phi(Iini,k),1,type_recv(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo

#    ifdef KeyTESTpBC
      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Send information
        do i=1,NeighborNumber
            Iini = InighostIndex_local(i)

                do j=1,BlockBdyghostNumber(i)
                    elem = BlockBdyghostIndex(j,i)
                    SharePhi(j) = phi(elem,k)
                enddo

        call MPI_SEND(SharePhi(1),1,type_send_ghost(i),Neighbors(i),&
                      etiquette,comm3D,code)
        enddo

      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
         do i=1,NeighborNumber
            Iini = InighostIndex_local(i)
            call MPI_RECV(phi(Iini,k),1,type_recv_ghost(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo
#     endif
!      ========================
!      ====Parallel SENDRECV===
!      ========================
      else
!     ______________________________
!     Send information
         do iproc=1,NeighborNumber
            i = index_neigbour(iproc)
            Iini = IniExtraIndex_local(i)

!           -----------
            do j=1,BlockBdyNumber(i)
               elem = BlockBdyIndex(j,i)
               SharePhi(j) = phi(elem,k)
            enddo

            call MPI_SENDRECV(SharePhi(1), 1, type_send(i), Neighbors(i), 2*k, &
                              phi(Iini,k), 1, type_recv(i), Neighbors(i), 2*k, &
                              comm3D,statut,code)
#           ifdef KeyParallelNZ
            endif
#           endif
         enddo

#    ifdef KeyTESTpBC
         call MPI_Barrier(comm3D,code)
!     ______________________________
!     Send information
         do iproc=1,NeighborNumber
            i = index_neigbour(iproc)
            Iini = InighostIndex_local(i)

                do j=1,BlockBdyghostNumber(i)
                    elem = BlockBdyghostIndex(j,i)
                    SharePhi(j) = phi(elem,k)
                enddo

                call MPI_SENDRECV(SharePhi(1),1,type_send_ghost(i),Neighbors(i),2*k,&
                                  phi(Iini,k),1,type_recv_ghost(i),Neighbors(i),2*k,&
                                  comm3D,statut,code)
        enddo
#     endif
       endif
      ENDDO

      END SUBROUTINE communication3D
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.7  |        Communication NZ                                |
!     |______|________________________________________________________|

       SUBROUTINE communicationNZ(phi)

      !------------------------------------------------------------------------!
      !                                                                        !
      !    The communication is done at the vertical overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                         !
      !------------------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,ksend,krecv,elem,Iini,iproc
      integer :: s,sm,smNZ,smLayer
!     -----------------------------
!     For NZ partition communication
        if(MPIOption .eq. 1) then
!           ---------------------
!           FIRST SEND
            do i=1,NeighborNumber
                sm = Neighbors(i) + 1
                smNZ = mod(sm-1,TotalSubdomains) + 1
                smLayer = (sm-1)/TotalSubdomains + 1
                if(smLayer .ne. NZLayer) then
                    ksend =IniExtraIndex_local(i)
                    do j=1,BlockBdyNumber(i)
                        elem = BlockBdyIndex(j,i)
                        SharePhi(j) = phi(elem,ksend)
                    enddo

                    call MPI_SEND(SharePhi(1), 1, type_send(i), &
                        Neighbors(i), etiquette,                &
                        comm3D,statut,code)
                endif
            enddo

            call MPI_Barrier(comm3D,code)
!           ------------------------
!           THEN RECV
            do i=1,NeighborNumber
                sm = Neighbors(i) + 1
                smNZ = mod(sm-1,TotalSubdomains) + 1
                smLayer = (sm-1)/TotalSubdomains + 1
                if(smLayer .ne. NZLayer) then
                    if(IniExtraIndex_local(i) .eq. 2) then
                        krecv = 1
                    elseif(IniExtraIndex_local(i) .eq. NZ-1) then
                        krecv = NZ
                    endif

                    call MPI_RECV(phi(1,krecv), 1, type_recv(i), &
                        Neighbors(i), etiquette,                 &
                        comm3D,statut,code)
                endif
            enddo
!     ========================
!     ===NON-BLOCK SENDRECV===
!     ========================
        else
            do i=1,NeighborNumber
                sm = Neighbors(i) + 1
                smNZ = mod(sm-1,TotalSubdomains) + 1
                smLayer = (sm-1)/TotalSubdomains + 1
                if(smLayer .ne. NZLayer) then

                    if(IniExtraIndex_local(i) .eq. 2) then
                        ksend = 2
                        krecv = 1
                    elseif(IniExtraIndex_local(i) .eq. NZ-1) then
                        ksend = NZ-1
                        krecv = NZ
                    endif

                    do j=1,BlockBdyNumber(i)
                        elem = BlockBdyIndex(j,i)
                        SharePhi(j) = phi(elem,ksend)
                    enddo

                    call MPI_SENDRECV(SharePhi(1), 1, type_send(i), &
                        Neighbors(i), etiquette,                    &
                        phi(1,krecv), 1, type_recv(i),              &
                        Neighbors(i), etiquette,                    &
                        comm3D,statut,code)


                endif
            enddo
        endif

      END SUBROUTINE communicationNZ
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.5  |     Communication 3D only nearest cc(vector type)      |
!     |______|________________________________________________________|

!       SUBROUTINE communication3Dcc(phi)
!
!      !---------------------------------------------------------------!
!      !                                                               !
!      !    The communication is done at the overlaping elements of a  !
!      !    variable phi with two entry phi(N_CELL,NZ).                !
!      !                                                               !
!      !---------------------------------------------------------------!
!
!      real*8,dimension(:,:):: phi
!      integer,parameter  :: etiquette=100
!      integer,dimension(MPI_STATUS_SIZE) :: statut
!      integer :: i,j,k,elem,Iini
!
!!     ______________________________
!!     Send information
!      do i=1,NeighborNumber
!!        -----------
!         do j=1,BlockBdyCCNumber(i)
!            elem = BlockBdyCCIndex(j,i)
!            do k=2,NZ-1
!               SharePhiNew(j,k-1) = phi(elem,k)
!            enddo
!         enddo
!!        -----------
!         call MPI_SEND(SharePhiNew(1,1),1,type_BdyCCVEC(i),Neighbors(i),&
!                       etiquette,comm3D,code)
!      enddo
!      call MPI_Barrier(comm3D,code)
!!     ______________________________
!!     Receive information
!      do i=1,NeighborNumber
!         Iini = IniExtraIndex_local(i)
!         call MPI_RECV(phi(Iini,2),1,type_BdyCCVEC(i),Neighbors(i),&
!                       etiquette,comm3D,statut,code)
!      enddo
!
!
!      END SUBROUTINE communication3Dcc

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!           MODULE PARALLEL (ESPECIAL FUNCTIONS: min,max,sum)         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!      SUBROUTINES:     3.1)  MIN_parallel                            !
!                       3.2)  MAX_parallel                            !
!                       3.3)  SUM_parallel                            !
!                                                                     !
!---------------------------------------------------------------------!

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.1  |              MIN value all processors                  |
!     |______|________________________________________________________|

       SUBROUTINE MIN_parallel(value,valuemin)

      !---------------------------------------------------------------!
      !                                                               !
      !    Return the minimum value among all subdomains:             !
      !    <--- valuemin = min{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------!

       real*8 :: value,valuemin

       call MPI_ALLREDUCE(value,valuemin,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
                          comm3D,code)

       END SUBROUTINE MIN_parallel

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.2  |               MAX value all processors                 |
!     |______|________________________________________________________|

       SUBROUTINE MAX_parallel(value,valuemax)

      !---------------------------------------------------------------!
      !                                                               !
      !    Return the maximum value among all subdomains:             !
      !    <--- valuemax = max{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------!

       real*8 :: value,valuemax

       call MPI_ALLREDUCE(value,valuemax,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
                          comm3D,code)

       END SUBROUTINE MAX_parallel

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.2  |               MAX value all processors                 |
!     |______|________________________________________________________|

       SUBROUTINE MAX_paralleli(value)

      !---------------------------------------------------------------!
      !                                                               !
      !    Return the maximum value among all subdomains:             !
      !    <--- valuemax = max{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------!

       integer :: value,valuemax

       call MPI_ALLREDUCE(value,valuemax,1,MPI_INTEGER,MPI_MAX,&
                          comm3D,code)

       value = valuemax

       END SUBROUTINE MAX_paralleli

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.3  |                 SUM value all processors               |
!     |______|________________________________________________________|

       SUBROUTINE SUM_parallel(value,valuesum)

      !---------------------------------------------------------------!
      !                                                               !
      !    Return the addiction of a value among all subdomains:      !
      !    <--- valuesum = sum{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------!

       real*8 :: value,valuesum

       call MPI_ALLREDUCE(value,valuesum,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
                          comm3D,code)

       END SUBROUTINE SUM_parallel
!      _______________________________________________________________
!     |      |                                                        |
!     | 3.4  |            SUM value all processors (size n)           |
!     |______|________________________________________________________|

       SUBROUTINE SUM_paralleli(value,n)

      !---------------------------------------------------------------!
      !                                                               !
      !    Return the addiction of a value among all subdomains:      !
      !    <--- valuesum = sum{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------!
       integer :: n
       integer,dimension(:) :: value(n)
       integer,dimension(:) :: valuesum(n)

       call MPI_ALLREDUCE(value,valuesum,n,MPI_INTEGER,MPI_SUM,&
                          comm3D,code)

       do i=1,n
          value(i) = valuesum(i)
       enddo

       END SUBROUTINE SUM_paralleli
      !---------------------------------------------------------------!
!     |      |                                                        |
!     | 3.5  |           SUM value all processors (size(n1,n2)        |
!     |______|________________________________________________________|

       SUBROUTINE SUM_parallelr(value,n)

      !---------------------------------------------------------------!
      !                                                               !
      !    Return the addiction of a value among all subdomains:      !
      !    <--- valuesum = sum{value}_proc                            !
      !                                                               !
      !---------------------------------------------------------------!
       integer :: n
       real*8,dimension(:) :: value(n)
       real*8,dimension(:) :: valuesum(n)
       call MPI_ALLREDUCE(value,valuesum,n,MPI_DOUBLE_PRECISION,MPI_SUM,&
                            comm3D,code)

         do i=1,n
            value(i) = valuesum(i)
         enddo

       END SUBROUTINE SUM_parallelr

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                   MODULE PARALLEL (GLOBAL MATRIX)                   !
!                              Jul 2014                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!---------------------------------------------------------------------!
!                                                                     !
!      SUBROUTINES:    4.1)  matgloV                                  !
!                      4.2)  matgloC                                  !
!                      4.3)  matgloVFull                              !
!                      4.4)  matgloCFull                              !
!                                                                     !
!---------------------------------------------------------------------!

!      _______________________________________________________________
!     |      |                                                        |
!     | 4.1  |           Building the global matrix (vertices)        |
!     |______|________________________________________________________|

      SUBROUTINE matgloV(phiv,phiv_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phiv_global(N_VERTglobal)                   !
      !    using all the subdomain matrices phiv(N_VERT).             !
      !    Global data saved only on proc# 0                          !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phiv,phiv_global
      integer,parameter :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,vert,ii,Iini,Ifin
      integer :: sm,smNZ,smLayer,kvert
      real*8,dimension(:,:),allocatable :: arglo


      allocate(arglo(N_VERTglobal,NZ-1))
!      ___________________________________________
!     |                                           |
!     |                ONE PROCESSOR              |
!     |___________________________________________|

      IF (Nprocs.eq.1) THEN
         do k=1,NZ-1
            do nv=1,N_VERT
               vert = Index_globalv(nv)
               phiv_global(vert,k) = phiv(nv,k)
            enddo
         enddo
!      ___________________________________________
!     |                                           |
!     |          MORE THAN ONE PROCESSOR          |
!     |___________________________________________|

      ELSE
!       ___________________________________
!       SENDIND INFORMATION
        if (rang_topo.eq.0) then
            do k=1,NZ-1
                do nv=1,N_VERT
                    vert = Index_globalv(nv)
                    phiv_global(vert,k) = phiv(nv,k)
                enddo
            enddo
        else
            do k=1,NZ-1
                do nv=1,N_VERT
                    arglo(nv,k) = phiv(nv,k)
                enddo
            enddo

            call MPI_SEND(arglo(1,1),1,type_blocV,0,&
                etiquette,comm3D,code)
        endif
!       ___________________________________
!       RECEIVING INFORMATION
        if (rang_topo.eq.0) then
	    do s=1,Nprocs-1

	       call MPI_RECV(arglo(1,1),1,type_blocV,s,&
     		             etiquette,comm3D,statut,code)

               Iini = DimVVdom(s+1,2)
               Ifin = DimVVdom(s+1,3)

               nv = 0
               do i=Iini,Ifin
                  nv = nv + 1
                  vert = VVdom(i,1)
                  do k=1,NZglobal-1
                     phiv_global(vert,k) = arglo(nv,k)
                  enddo
               enddo
	    enddo
        endif
      ENDIF

      deallocate(arglo)

      END SUBROUTINE matgloV

!      _______________________________________________________________
!     |      |                                                        |
!     | 4.2  |   Building the global matrix (cell-centers)            |
!     |______|________________________________________________________|

      SUBROUTINE matgloC(phi,phi_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phi_global(N_CELL0global)                   !
      !    using all the subdomain matrices phi(N_CELL0).             !
      !    Global data saved only on proc# 0                          !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phi,phi_global
      integer,parameter :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,elem,ii,Iini,Ifin
      integer :: sm,smNZ,smLayer,kc
      real*8,dimension(:,:),allocatable :: argloC


      allocate(argloC(N_CELL0global,NZ))
!      ___________________________________________
!     |                                           |
!     |                ONE PROCESSOR              |
!     |___________________________________________|

      IF (Nprocs.eq.1) THEN
         do k=1,NZ
            do i=1,N_CELL0
               elem = Index_global(i)
               phi_global(elem,k) = phi(i,k)
            enddo
         enddo
!      ___________________________________________
!     |                                           |
!     |          MORE THAN ONE PROCESSOR          |
!     |___________________________________________|

      ELSE
!           ___________________________________
!           SENDIND INFORMATION
            if (rang_topo.eq.0) then
                do k=1,NZ
                    do i=1,N_CELL0
                        elem = Index_global(i)
                        phi_global(elem,k) = phi(i,k)
                    enddo
                enddo
            else
                do k=1,NZ
                    do i=1,N_CELL0
                        argloC(i,k) = phi(i,k)
                    enddo
                enddo

                    call MPI_SEND(argloC(1,1),1,type_blocC,0,&
                        etiquette,comm3D,code)
 !                   print*,'send from proc#',rang_topo,'to proc# 0'
            endif
!           ___________________________________
!           RECEIVING INFORMATION
        if (rang_topo.eq.0) then
!       =====================
!       NO NZ Partition
#       ifndef KeyParallelNZ
            do s=1,Nprocs-1

                call MPI_RECV(argloC(1,1),1,type_blocC,s,&
                              etiquette,comm3D,statut,code)

                Iini = DimCCdom(s+1,2)
                Ifin = DimCCdom(s+1,3)
                nc = 0
                do i=Iini,Ifin
                    nc = nc + 1
                    elem = CCdom(i,1)
                    do k=1,NZ
                        phi_global(elem,k) = argloC(nc,k)
                    enddo
                enddo
            enddo
!       =====================
!        NZ Partition
#        else
            do s=1,Nprocs-1
                sm = s + 1
                smNZ = mod(sm-1,TotalSubdomains) + 1
                smLayer = (sm-1)/TotalSubdomains + 1

                call MPI_RECV(argloC(1,1),1,type_blocC,s,&
                    etiquette,comm3D,statut,code)

!                print*, 'recv from proc#',s

                if(smNZ .ne. sNZ) then
                    Iini = DimCCdom(smNZ,2)
                    Ifin = DimCCdom(smNZ,3)
                    nc = 0
                    do i=Iini,Ifin
                        nc = nc + 1
                        elem = CCdom(i,1)
                        do k=1,NZ
                            kc = (NZ-2)*(smLayer-1) + k
                            phi_global(elem,kc) = argloC(nc,k)
                        enddo
                    enddo
                else
                    do k=1,NZ
                        kc = (NZ-2)*(smLayer-1) + k
                        do i=1,N_CELL0
                            elem = Index_global(i)
                            phi_global(elem,kc) = argloC(i,k)
                        enddo
                    enddo
                endif
            enddo
#           endif
!   =========END NZ PARTITION=======
        endif
      ENDIF

      deallocate(argloC)

      END SUBROUTINE matgloC
!      _______________________________________________________________
!     |      |                                                        |
!     | 4.1  |           Building the global matrix (vertices)        |
!     |______|________________________________________________________|

      SUBROUTINE matgloVFull(phiv,phiv_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phiv_global(N_VERTglobal)                   !
      !    using all the subdomain matrices phiv(N_VERT).             !
      !    Global data saved only on all procs                        !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phiv,phiv_global

      if(Nprocs .gt. 1) then
           do k= 1,NZglobal-1
            do i=1,N_VERTglobal
               phiv_global(i,k) = 0.
            enddo
           enddo
!         ------------------------------
!          collect all info on proc#0
          call matgloV(phiv,phiv_global)
          call MPI_Barrier(comm3D,code)
          do k=1,NZglobal-1
            call SUM_parallelr(phiv_global(1:N_VERTglobal,k),N_VERTglobal)
          enddo
      endif

      END SUBROUTINE matgloVFull

!      _______________________________________________________________
!     |      |                                                        |
!     | 4.2  |   Building the global matrix (cell-centers)            |
!     |______|________________________________________________________|

      SUBROUTINE matgloCFull(phi,phi_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phi_global(N_CELL0global)                   !
      !    using all the subdomain matrices phi(N_CELL0).             !
      !    Global data saved only on all procs                        !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phi,phi_global

      if(Nprocs .gt. 1) then
           do k= 1,NZglobal
            do i=1,N_CELL0global
               phi_global(i,k) = 0.
            enddo
           enddo
!         ------------------------------
!          collect all info on proc#0
          call matgloC(phi,phi_global)
          call MPI_Barrier(comm3D,code)
!         ------------------------------
!          copy info to all procs
          do k=1,NZglobal
            call SUM_parallelr(phi_global(1:N_CELL0global,k),N_CELL0global)
          enddo
      endif

      END SUBROUTINE matgloCFull
!      _______________________________________________________________
!     |      |                                                        |
!     | 4.5  |   Building the global matrix (cell-centers)            |
!     |______|________________________________________________________|

      SUBROUTINE matgloC_2D(phi,phi_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phi_global(N_CELL0global)                   !
      !    using all the subdomain matrices phi(N_CELL0).             !
      !                                                               !
      !---------------------------------------------------------------!

      real*8, dimension(:) :: phi,phi_global
      integer,parameter :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,elem,ii,Iini,Ifin
      real*8,dimension(:),allocatable :: argloC


      allocate(argloC(N_CELL0global))
!      ___________________________________________
!     |                                           |
!     |                ONE PROCESSOR              |
!     |___________________________________________|

      IF (Nprocs.eq.1) THEN
            do i=1,N_CELL0
               elem = Index_global(i)
               phi_global(elem) = phi(i)
            enddo
!      ___________________________________________
!     |                                           |
!     |          MORE THAN ONE PROCESSOR          |
!     |___________________________________________|

      ELSE
!           ___________________________________
!           SENDIND INFORMATION
            if (rang_topo.eq.0) then
               do i=1,N_CELL0
                  elem = Index_global(i)
                  phi_global(elem) = phi(i)
               enddo
            else
               do i=1,N_CELL0
                  argloC(i) = phi(i)
               enddo
#           ifndef KeyParallelNZ
               call MPI_SEND(argloC(1),1,type_bloc2D,0,&
                             etiquette,comm3D,code)
#           else
               if(NZLayer .eq. 1) then
                 call MPI_SEND(argloC(1),1,type_bloc2D,0,&
                             etiquette,comm3D,code)
               endif
#           endif
            endif
!           ___________________________________
!           RECEIVING INFORMATION
            if (rang_topo.eq.0) then
               do s=1,TotalSubdomains-1
                  call MPI_RECV(argloC(1),1,type_bloc2D,s,&
                                etiquette,comm3D,statut,code)
                  Iini = DimCCdom(s+1,2)
                  Ifin = DimCCdom(s+1,3)
                  nc = 0
                  do i=Iini,Ifin
                     nc = nc + 1
                     elem = CCdom(i,1)
                     phi_global(elem) = argloC(nc)
                  enddo
               enddo
            endif
      ENDIF

      deallocate(argloC)

      END SUBROUTINE matgloC_2D

      END MODULE parallel

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                         END MODULE PARALLEL                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
