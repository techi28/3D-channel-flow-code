!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        INPUT PARALLEL FILE                          !
!                      Miguel Angel Uh Zapata                         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE parallel_input_local(No_vp,No_cp,No_wb,No_qb,No_hb,No_sp,&
                                      nbe,xv,yv,zbv)

!---------------------------------------------------------------------!
!                                                                     !
!     Transfer information from global variables to local             !
!     variables in each subdomain. These variables correspond to      !
!     the ones read in data.txt.                                      !
!---------------------------------------------------------------------!
!                                                                     !
!    Variables:                                                       !
!   _______________________________________________________________   !
!  |        Name      |             Description                    |  !  
!  |__________________|____________________________________________|  !   
!  |                  |                                            |  !
!  | <--- No_vp       |  Coming from: No_vp_global                 |  !
!  | <--- No_cp       |  Coming from: No_cp_global                 |  !
!  | <--- No_wp       |  Coming from: No_wp_global                 |  !
!  | <--- No_qp       |  Coming from: No_qp_global                 |  !
!  | <--- No_hp       |  Coming from: No_hp_global                 |  !
!  | <--- No_sp       |  Coming from: No_sp_global                 |  !
!  | <--- nbe         |  Coming from: nbe_global                   |  !
!  | <--- xv          |  Coming from: xv_global                    |  !
!  | <--- yv          |  Coming from: yv_global                    |  !
!  | <--- zbv         |  Coming from: zbv_global                   |  !
!  |__________________|____________________________________________|  ! 
!                                                                     !
!   REMARK: This are not the final variables, the good ones will be   !
!           without the word "_local", but they need to be define     !
!           after the definition of local variables.                  !
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

      integer,dimension(:,:):: No_vp  
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: No_wb 
      integer,dimension(:)  :: No_qb
      integer,dimension(:)  :: No_hb 
      integer,dimension(:)  :: No_sp 
      integer,dimension(:)  :: nbe 
      real*8, dimension(:)  :: xv
      real*8, dimension(:)  :: yv  
      real*8 ,dimension(:)  :: zbv
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer :: ielem,jelem,kelem
      integer :: ivert,jvert,kvert
      integer :: jv1,jv2,jv3
      real*8  :: tag

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: parallel_input_local'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                        Tranfer data variables                       !
!                                                                     !
!*********************************************************************!

!      ________________________________________________________
!     |                                                        |
!     |                    No_cp(N_CELL,3)                     |
!     |________________________________________________________|

!     ______________________________________
!     Independent cell-centers (inside + BC)

      do i=1,N_CELL0
         ielem = Index_global(i)
         do j=1,3
             Jelem = No_cp_global(ielem,j)
             if (Jelem.eq.0) then
                No_cp(i,j) = 0
             else
                tag=0
                do k=1,N_CELL
                   kelem = Index_global(k)
                   if (Jelem.eq.kelem) then
                       No_cp(i,j) = k
                       tag=1
                   endif
                enddo
                if (tag.eq.0) then
                    print*,'Processor:',rang_topo,'Error!!!: Missing cell-centers'
                    print*,'i=',i,',j',j,',No_cp_global=',Jelem
                    stop
                endif
             endif
         enddo
      enddo
!     ______________________________________
!     For the overlaping elements: No_cp = 1

      do i=N_CELL0+1,N_CELL
         No_cp(i,1) = 0
         No_cp(i,2) = 0
         No_cp(i,3) = 0
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                      No_vp(N_CELL0,3)                  |
!     |________________________________________________________|

      do i=1,N_CELL0
         ielem = Index_global(i)
         do j=1,3
             Jvert = No_vp_global(ielem,j)
             tag = 0
             do nv=1,N_VERT
                kvert = Index_globalv(nv)
                if (Jvert.eq.kvert) then
                   No_vp(i,j) = nv
                   tag = 1
                endif
             enddo
             if (tag.eq.0) then
                print*,'Processor:',rang_topo,'Error!!!: Missing vertices'
                print*,'i=',i,',j',j,',No_vp_global=',Jvert
                stop
             endif
         enddo
      enddo
!      ________________________________________________________
!     |                                                        |
!     |                        nbe(N_CELL0)                    |
!     |________________________________________________________|

      do i=1,N_CELL0
         ielem  = Index_global(i)
         nbe(i) = nbe_global(ielem)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |          xv(N_VERT), yv(N_VERT), zbv(N_VERT)           |
!     |________________________________________________________|

      do nv=1,N_VERT
         ivert = index_globalv(nv)
         xv(nv)  = xv_global(ivert)
         yv(nv)  = yv_global(ivert)
         zbv(nv) = zbv_global(ivert)
      enddo

!      ________________________________________________________
!     |                                                        |
!     |          No_wb(N_WB), No_qb(N_QB), No_hb(N_HB)         |
!     |________________________________________________________|

      do nv=1,N_WB
         No_wb(nv) = wb_aux(nv)
      enddo
      do nv=1,N_QB
         No_qb(nv) = qb_aux(nv)
      enddo
      do nv=1,N_HB
         No_hb(nv) = hb_aux(nv)
      enddo
!      ________________________________________________________
!     |                                                        |
!     |          Sample vertex points: No_sp(N_SP)             |
!     |________________________________________________________|

      if (N_SP.ge.1) then
         do nv=1,N_SP
            No_sp(nv) = sp_aux(nv)
         enddo
      endif


!*********************************************************************!
!                                                                     !
!                             Finalization                            !
!                                                                     !
!*********************************************************************!

      deallocate(wb_aux,qb_aux,hb_aux,sp_aux)

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> End   subroutine: parallel_input_local'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF INPUT                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

