!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        PARALLEL PARAMETERS                          !
!                      Miguel Angel Uh Zapata                         !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE parallel_parameters

!---------------------------------------------------------------------!
!                                                                     !
!     Calculate the local parameters corresponding to the ones ini-   !
!     tially was read in the input data of the sequential version.    !
!     It also calculate the number of ghost cell-centers and finally  !
!     the total number of cell-centers. At the end auxiliar vectors   !
!     are used to storage the local number of boundary elements.      !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Parameters:                                                      !
!   _______________________________________________________________   !
!  |      Name     |                Description                    |  !
!  |_______________|_______________________________________________|  !
!  |               |                                               |  !
!  |  NZ           | Number of vertical variables (vertex points)  |  !
!  |  N_VERT       | Number of vertices in each subdomain          |  !
!  |  N_CELL0      | Number of cell-centers in each subdomain      |  !
!  |  N_CELLghost  | Number of ghost cell-centers                  |  !
!  |  N_CELLexact  | Number of interior + ghost cell-centers       |  !
!  |  N_CELLextra  | Number of overlaping cell-centers (extra)     |  !
!  |  N_CELL       | Total number of cell-centers                  |  !
!  |_______________|_______________________________________________|  !
!                                                                     !
!    We need variables:                                               !
!   _______________________________________________________________   !
!  |                                                               |  !
!  |   - CCdom             ( parallel_distribution.F90 )           |  !
!  |   - DimVVdom          ( parallel_distribution.F90 )           |  !
!  |   - VVdom             ( parallel_distribution.F90 )           |  !
!  |   - DimCCdom          ( parallel_distribution.F90 )           |  !
!  |   - DimExtraCCdom     ( parallel_distribution.F90 )           |  !
!  |   - No_cp_global      ( parallel_input.F90        )           |  !
!  |   - nbe_global        ( parallel_input.F90        )           |  !
!  |   - nbev_global       ( parallel_input.F90        )           |  !
!  |_______________________________________________________________|  !
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
!    |   Declaration of local variables   |
!    |____________________________________|

      integer :: ii,Iini,Ifin,ielem,s
      integer :: ivert,jvert
      integer :: nv1,nv2,nv3,nv4

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: parallel_parameters'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Initialization                           !
!                                                                     !
!*********************************************************************!
        s = rang_topo + 1
!      ________________________________________________________
!     |                                                        |
!     |              N_CELL0, N_VERT, N_CELLextra              |
!     |________________________________________________________|

      N_VERT  = DimVVdom(s,1)
      N_CELL0 = DimCCdom(s,1)
      N_CELLextra = DimExtraCCdom(s,1)

!      ________________________________________________________
!     |                                                        |
!     |    Exact number of ghost points:  N_CELLghost          |
!     |________________________________________________________|

      ii = 0
      Iini = DimCCdom(s,2)
      Ifin = DimCCdom(s,3)
        do i=Iini,Ifin
            ielem  = CCdom(i,1)
            if (nbe_global(ielem).gt.0) then
                do j=1,3
                    nc = No_cp_global(ielem,j)
                    if ((nc.lt.1).or.(nc.gt.N_CELL0global)) then
                        ii=ii+1
                    endif
                enddo
            endif
        enddo

       N_CELLghost = ii
!      ________________________________________________________
!     |                                                        |
!     |              Total number of cells: N_CELL             |
!     |________________________________________________________|

      N_CELLexact = N_CELL0 + N_CELLghost
      N_CELL      = N_CELL0 + N_CELLghost + N_CELLextra

!      ________________________________________________________
!     |                                                        |
!     |       Total number of boundary vertex: N_BC            |
!     |________________________________________________________|
      ii = 0
      Iini = DimVVdom(s,2)
      Ifin = DimVVdom(s,3)
        do i=Iini,Ifin
            ielem  = VVdom(i,1)
            if (nbev_global(ielem).gt.0)  ii=ii+1
        enddo

      N_BC = ii
!*********************************************************************!
!                                                                     !
!                             Finalization                            !
!                                                                     !
!*********************************************************************!

#     ifdef KeyDisplayParallel
        IF (rang_topo.eq.1) THEN
            print*,'_________________________________________________________'
            print*,'                                                         '
            print*,'       ===========================================       '
            print*,'               LOCAL SUBDOMAIN PARAMETERS                '
            print*,'       ===========================================       '
            print*,'                                                         '
        ENDIF
        print*,'SUBDOMAIN:',s-1
        print*,'   N_VERT      =',N_VERT
        print*,'   N_CELL0     =',N_CELL0
        print*,'   N_CELLextra =',N_CELLextra
        print*,'   N_CELLghost =',N_CELLghost
        print*,'   N_CELL      =',N_CELL
        print*,'   N_BC        =',N_BC
        print*,'   --------------------------------------'
        print*,'   '
#     endif
      call MPI_Barrier(comm3D,code)
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> End   subroutine: parallel_parameters'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                     END OF parallel_parameters                      !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

