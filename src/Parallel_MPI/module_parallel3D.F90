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
      integer :: TotalSubdomains
      integer :: NeighborNumber
      integer :: NextraMax,NshareMax
      integer :: N_VERTmax
!     ------------------------------------
      integer,dimension(:), allocatable :: edges
      integer,dimension(:), allocatable :: indexval
      integer,dimension(:), allocatable :: Neighbors
!     ------------------------------
      integer,dimension(:), allocatable :: type_send
      integer,dimension(:), allocatable :: type_recv
      integer,dimension(:), allocatable :: type_sendVEC
      integer,dimension(:), allocatable :: type_recvVEC
      integer,dimension(:), allocatable :: type_BdyCCVEC
      integer :: type_bloc1
      integer :: type_blocC
      integer :: type_bloc1D
      integer :: type_continuousC
      integer :: type_VERTsend
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
      integer,dimension(:),   allocatable :: Local_index_Inter
      integer,dimension(:),   allocatable :: Local_index_Bound
      integer,dimension(:),   allocatable :: Index_global
      integer,dimension(:),   allocatable :: Index_globalv
      integer,dimension(:),   allocatable :: IniExtraIndex_local
      integer,dimension(:,:), allocatable :: BlockBdyIndex
      integer,dimension(:),   allocatable :: BlockBdyNumber
      integer,dimension(:,:), allocatable :: BlockBdyCCIndex
      integer,dimension(:),   allocatable :: BlockBdyCCNumber
      integer,dimension(:),   allocatable :: BlockIntIndex
      integer,dimension(:,:), allocatable :: InterCCdom
      integer,dimension(:,:), allocatable :: DimInterCCdom
      integer,dimension(:,:), allocatable :: BoundCCdom
      integer,dimension(:,:), allocatable :: DimBoundCCdom
      integer :: BlockIntNumber,NN_Inter,NN_Bound
      real*8, dimension(:),   allocatable :: SharePhi
      real*8, dimension(:,:), allocatable :: SharePhiNew
!     ------------------------------
      integer,dimension(:),   allocatable :: wb_aux,qb_aux,hb_aux,sp_aux
!     ------------------------------
      integer,dimension(:),   allocatable :: TagParallel
!     ------------------------------
      integer,dimension(:,:), allocatable :: No_vp_global
      integer,dimension(:,:), allocatable :: No_cp_global
      integer,dimension(:),   allocatable :: nbe_global
      integer,dimension(:),   allocatable :: No_wb_global
      integer,dimension(:),   allocatable :: No_qb_global
      integer,dimension(:),   allocatable :: No_hb_global
      integer,dimension(:),   allocatable :: No_sp_global
      real*8, dimension(:),   allocatable :: xv_global
      real*8, dimension(:),   allocatable :: yv_global
      real*8, dimension(:),   allocatable :: zbv_global
!     ------------------------------
      integer :: NZdiv,div,NZres
      real*8, dimension(:,:), allocatable :: SharePhiDiv1,SharePhiDiv2
      integer,dimension(:),   allocatable :: type_sendDiv1,type_sendDiv2
      integer,dimension(:),   allocatable :: type_recvDiv1,type_recvDiv2
      integer :: type_vertical

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
      print*,'          Global number of NZ points =',NZglobal-1
      print*,'                                                         '
      ENDIF

      END SUBROUTINE initialisation_mpi

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

      s = rang_topo + 1

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

     IF (rang_topo.eq.0) THEN
        print*,'   ____________________________________________________  '
        print*,'                                                         '
        print*,'       NEIGHBORS:                                        '
        print*,'                                                         '
        print*,'    p  | N_VERT | N_CELL0| N_CELL | No. Neigh & Neighbors'
        print*,'   ----|--------|--------|--------|----------------------'
     ENDIF
     call MPI_Barrier(comm3D,code)
     write(*,9),rang_topo,' |',N_VERT,' |',N_CELL0,' |',N_CELL, ' |',&
                NeighborNumber,': ',Neighbors(1:NeighborNumber)
     9 format(I6,a3,I6,a3,I6,a3,I6,a3,I2,a3,10I4)


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

      integer :: Iini,Ifin,i0,k0,ielem,jelem,ivert,s,sm,suma
      integer :: k,ii,jj,tag
      integer :: nv1,nv2,nv3

      s = rang_topo + 1

!     ________________________________________________________
!     Initial index of extra cells by neighborn

      allocate(IniExtraIndex_local(NeighborNumber))

      s  = rang_topo + 1
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

!     -----------------------------------------
!     Virtual boundary cell-center values (index)

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

!     -----------------------------------------
!     Virtual boundary index of the neardest cell-center values

      allocate(BlockBdyCCIndex(NextraMax,NeighborNumber))
      allocate(BlockBdyCCNumber(NeighborNumber))

      s  = rang_topo + 1
      do k=1,NeighborNumber
         sm = Neighbors(k) + 1
         Iini = DimBdyCCdom(s,2)
         Ifin = DimBdyCCdom(s,3)
         i0 = 0
         do j=Iini,Ifin
            if (BdyCCdom(s,3).eq.sm) then
               i0 = i0 + 1
               BlockBdyCCIndex(i0,k)= BdyCCdom(s,1)
            endif
         enddo
         BlockBdyCCNumber(k) = i0
      enddo

!     -----------------------------------------
!     Internal cell-center values (Index)

      allocate(BlockIntIndex(N_CELL0))

      i0 = 0
      do i=1,N_CELL0
         tag = 0
         do k=1,NeighborNumber
            do jj=1,BlockBdyNumber(k)
               ii = BlockBdyIndex(jj,k)
               if (i.eq.ii) then
                  tag = 1
               endif
            enddo
         enddo
         if (tag.eq.0) then
             i0 = i0 + 1
             BlockIntIndex(i0) = i
         endif
      enddo
      BlockIntNumber = i0

!     ________________________________________________________
!     Internal and Boundary elements

      NN_Inter = DimInterCCdom(s,1)
      NN_Bound = DimBoundCCdom(s,1)
      allocate(Local_index_Inter(NN_Inter))
      allocate(Local_index_Bound(NN_Bound))

      Iini = DimBoundCCdom(s,2)
      Ifin = DimBoundCCdom(s,3)
      i0 = 0
      k0 = 0
      do i=1,N_CELL0
         ielem = Index_global(i)
         tag = 0
         do jj = Iini,Ifin
            jelem = BoundCCdom(jj,1)
            if (ielem.eq.jelem) tag = 1
         enddo
         if (tag.eq.0) then
            i0 = i0 + 1
            Local_index_Inter(i0) = i
         else
            k0 = k0 + 1
            Local_index_Bound(k0) = i
         endif
      enddo
!     ________________________________________________________
!     Free memory ExtraCCdom global

      deallocate(ExtraCCdom)
      deallocate(DimExtraCCdom)

      END SUBROUTINE parallel_shareindex

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
                 No_wb_global,  &
                 No_qb_global,  &
                 No_hb_global,  &
                 No_sp_global,  &
                 xv_global,     &
                 yv_global,     &
                 zbv_global)
!     __________________
!     Local variables
      deallocate(Index_global)
      deallocate(Index_globalv)
!     ---------------------
      deallocate(Local_index_Inter)
      deallocate(Local_index_Bound)
!     ---------------------
      deallocate(IniExtraIndex_local)
      deallocate(BlockBdyIndex)
      deallocate(BlockBdyNumber)
      deallocate(BlockBdyCCIndex)
      deallocate(BlockBdyCCNumber)
      deallocate(InterCCdom)
      deallocate(DimInterCCdom)
      deallocate(BoundCCdom)
      deallocate(DimBoundCCdom)
      deallocate(BlockIntIndex)
      deallocate(SharePhi)
      deallocate(SharePhiNew)
!     __________________
!     Communication variables
      deallocate(type_recv)
      deallocate(type_send)
      deallocate(type_recvVEC)
      deallocate(type_sendVEC)
      deallocate(type_BdyCCVEC)
!     __________________
      deallocate(SharePhiDiv1,SharePhiDiv2)
      deallocate(type_recvDiv1,type_recvDiv2)
      deallocate(type_sendDiv1,type_sendDiv2)
!     ________________________________________________________
!     Free block vertex type
      call MPI_TYPE_FREE(type_bloc1,code)
      call MPI_TYPE_FREE(type_blocC,code)
      call MPI_TYPE_FREE(type_bloc1D,code)
      call MPI_TYPE_FREE(type_VERTsend,code)
      call MPI_TYPE_FREE(type_continuousC,code)

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
!                      3.1)  communication2D                          !
!                      3.2)  communication3Dtype1                     !
!                      3.4)  communication3Dtype2                     !
!                      3.5)  communication3D                          !
!                      3.6)  communication3Dcc                        !
!                                                                     !
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
      allocate(type_BdyCCVEC(NeighborNumber))
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
         call MPI_TYPE_VECTOR(NZ,number_send,N_CELL,&
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
         call MPI_TYPE_VECTOR(NZ,number_recv,N_CELL,&
              MPI_DOUBLE_PRECISION,type_recvVEC(i),code)
         call MPI_TYPE_COMMIT(type_recvVEC(i),code)
      enddo

!     ________________________________________
!     Continuous type: share bdy: Only Cell-center(Vector type)

      do i=1,NeighborNumber
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_BdyCC = 0
         do j=DimBdyCCdom(s,2),DimBdyCCdom(s,3)
            if (BdyCCdom(s,2).eq.sm) then
               number_BdyCC = number_BdyCC + 1
            endif
         enddo
         call MPI_TYPE_VECTOR(NZ,number_BdyCC,N_CELL,&
              MPI_DOUBLE_PRECISION,type_BdyCCVEC(i),code)
         call MPI_TYPE_COMMIT(type_BdyCCVEC(i),code)
      enddo

!     ________________________________________
!     Vertex global block type  (Matglob)
      call MPI_TYPE_VECTOR(NZglobal-1,N_VERTglobal,N_VERTglobal,&
           MPI_DOUBLE_PRECISION,type_bloc1,code)
      call MPI_TYPE_COMMIT(type_bloc1,code)

!     ________________________________________
!     Vertex global block type  (Matglob 1D)
      call MPI_TYPE_VECTOR(1,N_VERTglobal,N_VERTglobal,&
           MPI_DOUBLE_PRECISION,type_bloc1D,code)
      call MPI_TYPE_COMMIT(type_bloc1D,code)
!     ________________________________________
!     Vertex global block type maximum N_VERT (Matglob)
      call MPI_ALLREDUCE(N_VERT,N_VERTmax,1, &
                         MPI_DOUBLE_PRECISION,MPI_MAX,comm3D,code)
!    -------------------
      call MPI_TYPE_CONTIGUOUS((NZglobal-1)*N_VERTmax,&
                                MPI_DOUBLE_PRECISION,type_VERTsend,code)
      call MPI_TYPE_COMMIT(type_VERTsend,code)

!     ________________________________________
!     Cell-center global continuous type
      call MPI_TYPE_CONTIGUOUS(N_CELL0global,&
                               MPI_DOUBLE_PRECISION,type_continuousC,code)
      call MPI_TYPE_COMMIT(type_continuousC,code)
!     ________________________________________
!     Cell-center global continuous type
      call MPI_TYPE_VECTOR(NZglobal,N_CELL0global,N_CELL0global,&
                          MPI_DOUBLE_PRECISION,type_blocC,code)
      call MPI_TYPE_COMMIT(type_blocC,code)
!     ________________________________________
!     Allocate auxiliar variable to communicate

      NshareMax = 0
      do i=1,NeighborNumber
         NshareMax = max(NshareMax,BlockBdyNumber(i))
      enddo

      allocate(SharePhi(NshareMax))
      allocate(SharePhiNew(N_CELL,NZ))

!     =================================================================
!     For communication type 2
!     __________
!     Block size
      NZdiv = 16
      div = floor(1.0d0*NZglobal/NZdiv)
      NZres = NZglobal-(div-1)*NZdiv
!     __________
!     Block type
      allocate(type_sendDiv1(NeighborNumber))
      allocate(type_recvDiv1(NeighborNumber))
      allocate(type_sendDiv2(NeighborNumber))
      allocate(type_recvDiv2(NeighborNumber))
      do i=1,NeighborNumber
         s  = rang_topo + 1
         sm = Neighbors(i) + 1
         number_send = DimExtraFrom(sm,s)
         number_recv = DimExtraFrom(s,sm)
         call MPI_TYPE_VECTOR(NZdiv,number_send,N_CELL,MPI_DOUBLE_PRECISION,type_sendDiv1(i),code)
         call MPI_TYPE_VECTOR(NZdiv,number_recv,N_CELL,MPI_DOUBLE_PRECISION,type_recvDiv1(i),code)
         call MPI_TYPE_VECTOR(NZres,number_send,N_CELL,MPI_DOUBLE_PRECISION,type_sendDiv2(i),code)
         call MPI_TYPE_VECTOR(NZres,number_recv,N_CELL,MPI_DOUBLE_PRECISION,type_recvDiv2(i),code)
         call MPI_TYPE_COMMIT(type_sendDiv1(i),code)
         call MPI_TYPE_COMMIT(type_recvDiv1(i),code)
         call MPI_TYPE_COMMIT(type_sendDiv2(i),code)
         call MPI_TYPE_COMMIT(type_recvDiv2(i),code)
      enddo
!     __________
!     Share
      allocate(SharePhiDiv1(N_CELL,NZdiv))
      allocate(SharePhiDiv2(N_CELL,NZres))


!     =================================================================
!     For vertical communication type (FS_SaveReference)
!     vector of NZ-2 elements separated by 1

      call MPI_TYPE_VECTOR(NZ-1,1,1,MPI_DOUBLE_PRECISION,&
                           type_vertical,code)
      call MPI_TYPE_COMMIT(type_vertical,code)

      END SUBROUTINE parallel_type

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
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini

!     ______________________________
!     Send information
      do i=1,NeighborNumber
!        -----------
         do j=1,BlockBdyNumber(i)
            elem = BlockBdyIndex(j,i)
            do k=1,NZ
               SharePhiNew(j,k) = phi(elem,k)
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
         call MPI_RECV(phi(Iini,1),1,type_recvVEC(i),Neighbors(i),&
                       etiquette,comm3D,statut,code)
      enddo


      END SUBROUTINE communication3Dtype1

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.4  | Communication 3D: vector type 2 (for high-resolutions) |
!     |______|________________________________________________________|

       SUBROUTINE communication3DtypeD(phi)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZdiv).             !
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini
      integer :: s,ks,kIni,kFin

!     ______________________________
!     Send information
      DO s=1,div-1
         kIni = 1 + (s-1)*(NZdiv)
         kFin = s*(NZdiv)
         do i=1,NeighborNumber
!           -----------
            do j=1,BlockBdyNumber(i)
               ks = 0
               do k = kIni,kFin
                  ks = ks + 1
                  SharePhiDiv1(j,ks) = phi(BlockBdyIndex(j,i),k)
               enddo
            enddo
!           -----------
            call MPI_SEND(SharePhiDiv1(1,1),1,type_sendDiv1(i),Neighbors(i),&
                          etiquette,comm3D,code)
         enddo

         call MPI_Barrier(comm3D,code)
!        ______________________________
!        Receive information
         do i=1,NeighborNumber
            Iini = IniExtraIndex_local(i)
            call MPI_RECV(phi(Iini,kIni),1,type_recvDiv1(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo
      ENDDO
      IF (div.gt.1) THEN
         kIni = 1 + (div-1)*(NZdiv)
         kFin = NZglobal
         do i=1,NeighborNumber
!           -----------
            do j=1,BlockBdyNumber(i)
               ks = 0
               do k=kIni,kFin
                  ks = ks + 1
                  SharePhiDiv2(j,ks) = phi(BlockBdyIndex(j,i),k)
               enddo
            enddo
!           -----------
            call MPI_SEND(SharePhiDiv2(1,1),1,type_sendDiv2(i),Neighbors(i),&
                          etiquette,comm3D,code)
         enddo
         call MPI_Barrier(comm3D,code)
!        ______________________________
!        Receive information
         do i=1,NeighborNumber
            Iini = IniExtraIndex_local(i)
            call MPI_RECV(phi(Iini,kIni),1,type_recvDiv2(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo
      ENDIF

      END SUBROUTINE communication3DtypeD

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.5  |        Communication 3D: NZ continuos (type 0)         |
!     |______|________________________________________________________|

       SUBROUTINE communication3Dtype2(phi)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                !
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini

!     ______________________________
!     Send information
      DO k=1,NZ
         do i=1,NeighborNumber
!           -----------
            do j=1,BlockBdyNumber(i)
               elem = BlockBdyIndex(j,i)
               SharePhi(j) = phi(elem,k)
            enddo
!           -----------
            call MPI_SEND(SharePhi(1),1,type_send(i),Neighbors(i),&
                          etiquette,comm3D,code)
         enddo
      ENDDO

      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
      DO k=1,NZ
         do i=1,NeighborNumber
            Iini = IniExtraIndex_local(i)
            call MPI_RECV(phi(Iini,k),1,type_recv(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo
      ENDDO

      END SUBROUTINE communication3Dtype2

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
      integer :: i,j,k,elem,Iini

!     ______________________________
!     Send information
      DO k=1,NZ
         do i=1,NeighborNumber
!           -----------
            do j=1,BlockBdyNumber(i)
               elem = BlockBdyIndex(j,i)
               SharePhi(j) = phi(elem,k)
            enddo
!           -----------
            call MPI_SEND(SharePhi(1),1,type_send(i),Neighbors(i),&
                          etiquette,comm3D,code)
         enddo
      ENDDO

      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
      DO k=1,NZ
         do i=1,NeighborNumber
            Iini = IniExtraIndex_local(i)
            call MPI_RECV(phi(Iini,k),1,type_recv(i),Neighbors(i),&
                          etiquette,comm3D,statut,code)
         enddo
      ENDDO

      END SUBROUTINE communication3D

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.5  |     Communication 3D only nearest cc(vector type)      |
!     |______|________________________________________________________|

       SUBROUTINE communication3Dcc(phi)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable phi with two entry phi(N_CELL,NZ).                !
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:,:):: phi
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: i,j,k,elem,Iini

!     ______________________________
!     Send information
      do i=1,NeighborNumber
!        -----------
         do j=1,BlockBdyCCNumber(i)
            elem = BlockBdyCCIndex(j,i)
            do k=1,NZ
               SharePhiNew(j,k) = phi(elem,k)
            enddo
         enddo
!        -----------
         call MPI_SEND(SharePhiNew(1,1),1,type_BdyCCVEC(i),Neighbors(i),&
                       etiquette,comm3D,code)
      enddo
      call MPI_Barrier(comm3D,code)
!     ______________________________
!     Receive information
      do i=1,NeighborNumber
         Iini = IniExtraIndex_local(i)
         call MPI_RECV(phi(Iini,1),1,type_BdyCCVEC(i),Neighbors(i),&
                       etiquette,comm3D,statut,code)
      enddo


      END SUBROUTINE communication3Dcc

!      _______________________________________________________________
!     |      |                                                        |
!     | 3.6  |           Communication 2D vertical                    |
!     |______|________________________________________________________|

       SUBROUTINE communication2Dvertical(phi,pro)

      !---------------------------------------------------------------!
      !                                                               !
      !    The communication is done at the overlaping elements of a  !
      !    variable ar with one entry phi(N_CELL).                    !
      !                                                               !
      !---------------------------------------------------------------!

      real*8,dimension(:):: phi
      integer :: pro
      integer,parameter  :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut

!     ______________________________
!     Send & Receive information
      if (pro.ne.0) then
         call MPI_SENDRECV(phi(1),1,type_vertical,0,etiquette,   &
                           phi(1),1,type_vertical,pro,etiquette, &
                           comm3D,statut,code)
      endif

      END SUBROUTINE communication2Dvertical


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
!     | 3.4  |     SUM value all processors in the Z-direction        |
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

      call MPI_ALLREDUCE(value,valuesum,n, &
                         MPI_DOUBLE_PRECISION,MPI_SUM,comm3D,code)

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
!                                                                     !
!---------------------------------------------------------------------!
!      _______________________________________________________________
!     |      |                                                        |
!     | 4.1.2|  Building the global matrix (vertices) - Version 2017  |
!     |______|________________________________________________________|

      SUBROUTINE matgloV(phiv,phiv_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phiv_global(N_VERTglobal)                   !
      !    using all the subdomain matrices phiv(N_VERT).             !
      !                                                               !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phiv,phiv_global
      integer,parameter :: etiquette=200
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,vert,ii,Iini,Ifin,jj
      real*8,dimension(:),allocatable :: arLoc
      real*8,dimension(:),allocatable :: arGlo
      integer :: NumLoc,NumGlo

      NumLoc = N_VERTmax
      NumGlo = N_VERTmax*Nprocs
      allocate(arLoc(NumLoc))
      allocate(arGlo(NumGlo))
      arLoc = 0.0d0
      arGlo = 0.0d0

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
        DO k=1,NZ-1
!          ___________________________________
!          SAVE INFORMATION FOR EACH PROCESSOR
           jj = 0
           do nv=1,N_VERT
              jj = jj + 1
              arLoc(jj) = phiv(nv,k)
           enddo
!          ___________________________________
!          SENDIND & RECEIVING INFORMATION
           call MPI_GATHER(arLoc,NumLoc,MPI_DOUBLE_PRECISION, &
                           arGlo,NumLoc,MPI_DOUBLE_PRECISION, &
                           0,comm3D,code)
!          ___________________________________
!          GLOBAL ASSIGNATION
           if (rang_topo.eq.0) then
              do s=0,Nprocs-1
                  Iini = DimVVdom(s+1,2)
                  Ifin = DimVVdom(s+1,3)
                  nv = s*N_VERTmax
                  do i=Iini,Ifin
                     nv = nv + 1
                     vert = VVdom(i,1)
                     phiv_global(vert,k) = arGlo(nv)
                  enddo
              enddo
           endif
        ENDDO
      ENDIF
!     ___________________________________
      deallocate(arLoc)
      deallocate(arGlo)

      END SUBROUTINE matgloV

!      _______________________________________________________________
!     |      |                                                        |
!     |4.1.2 |     Building the global matrix (original-version)      |
!     |______|________________________________________________________|

      SUBROUTINE matgloV_2017(phiv,phiv_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phiv_global(N_VERTglobal)                   !
      !    using all the subdomain matrices phiv(N_VERT).             !
      !                                                               !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phiv,phiv_global
      integer,parameter :: etiquette=100
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,vert,ii,Iini,Ifin
      real*8,dimension(:,:),allocatable :: arglo


      allocate(arglo(N_VERTglobal,NZglobal-1))
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
	   call MPI_SEND(arglo(1,1),1,type_bloc1,0,&
                         etiquette,comm3D,code)
        endif
!       ___________________________________
!       RECEIVING INFORMATION
        if (rang_topo.eq.0) then
	    do s=1,Nprocs-1
	       call MPI_RECV(arglo(1,1),1,type_bloc1,s,&
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

      END SUBROUTINE matgloV_2017

!      _______________________________________________________________
!     |      |                                                        |
!     | 4.1.3|  Building the global matrix (vertices) - Version 2018  |
!     |______|________________________________________________________|

      SUBROUTINE matgloV_2018(phiv,phiv_global)

      !---------------------------------------------------------------!
      !                                                               !
      !    Reconstruction of the global matrix at cell-centeres:      !
      !                   phiv_global(N_VERTglobal)                   !
      !    using all the subdomain matrices phiv(N_VERT).             !
      !                                                               !
      !---------------------------------------------------------------!

      real*8, dimension(:,:) :: phiv,phiv_global
      !integer,parameter :: etiquette=200
      integer,dimension(20) :: etiquette
      integer,dimension(MPI_STATUS_SIZE) :: statut
      integer :: s,i,k,vert,ii,Iini,Ifin
!      real*8,dimension(:,:),allocatable :: arglo
      real*8,dimension(:),allocatable :: arglo
      integer :: NumLoc,NumGlo

      NumLoc = N_VERTmax
      NumGlo = N_VERTmax*Nprocs

      do i=1,20
          etiquette(i) = 200 + i
      enddo


      allocate(arglo(NumGlo))
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
        DO k=1,NZ-1
           if (rang_topo.eq.0) then
              do nv=1,N_VERT
                 vert = Index_globalv(nv)
                 phiv_global(vert,k) = phiv(nv,k)
              enddo
           endif

           do s=1,Nprocs-1
!          ___________________________________
!          SENDIND INFORMATION PROC: s
           if (rang_topo.eq.s) then
              do nv=1,N_VERT
                 arglo(nv) = phiv(nv,k)
              enddo
              call MPI_SEND(arglo(1),1,type_bloc1D,0,&
                            etiquette(rang_topo),comm3D,code)
!          ___________________________________
!          RECEIVING INFORMATION  PROC: 0
           elseif (rang_topo.eq.0) then
               call MPI_RECV(arglo(1),1,type_bloc1D,s,&
                            etiquette(s),comm3D,statut,code)
               Iini = DimVVdom(s+1,2)
               Ifin = DimVVdom(s+1,3)
               nv = 0
               do i=Iini,Ifin
                  nv = nv + 1
                  vert = VVdom(i,1)
                  phiv_global(vert,k) = arglo(nv)
               enddo
           endif
           enddo
        ENDDO
      ENDIF

      deallocate(arglo)

      END SUBROUTINE matgloV_2018

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


      allocate(argloC(N_CELL0global,NZglobal))
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
            do s=1,Nprocs-1

                call MPI_RECV(argloC(1,1),1,type_blocC,s,&
                              etiquette,comm3D,statut,code)

                Iini = DimCCdom(s+1,2)
                Ifin = DimCCdom(s+1,3)
                nc = 0
                do i=Iini,Ifin
                    nc = nc + 1
                    elem = CCdom(i,1)
                    do k=1,NZglobal
                        phi_global(elem,k) = argloC(nc,k)
                    enddo
                enddo
            enddo
        endif
      ENDIF

      deallocate(argloC)

      END SUBROUTINE matgloC

!      _______________________________________________________________
!     |      |                                                        |
!     | 4.3  |   Building the global matrix (cell-centers) NO YET !!! |
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
               call MPI_SEND(argloC(1),1,type_continuousC,0,&
                             etiquette,comm3D,code)
            endif
!           ___________________________________
!           RECEIVING INFORMATION
            if (rang_topo.eq.0) then
               do s=1,Nprocs-1
                  call MPI_RECV(argloC(1),1,type_continuousC,s,&
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
