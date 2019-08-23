!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        INPUT PARALLEL FILE                          !
!                      Miguel Angel Uh Zapata                         !
!                             Jun 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE parallel_input_local(No_vp,No_cp,nbe,nbev,xv,yv,zbv)

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
      integer,dimension(:)  :: nbe
      integer,dimension(:)  :: nbev
      real*8, dimension(:)  :: xv
      real*8, dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer :: ielem,jelem,kelem,ii,jj
      integer :: ivert,jvert,kvert
      integer :: jv1,jv2,jv3,s0,j0
      real*8  :: tag,xvert,yvert,x1vert,y1vert
      real*8  :: deltax,deltay
      integer :: tag_local,tag_global ! corner tag
      integer :: countpv,countpv2,countpv4,countcp
#     ifdef KeyDbgPBC
      character*80 filen
#     endif

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
#    ifdef KeyDbgPBC
         filen='../output/PBC/index_global-  .txt'
         write(filen(28:29),'(i2.2)') rang_topo
         open(180,file=filen)
         write(180,*),'Proc =',rang_topo,'N_CELL0=',N_CELL0,'N_CELL=',N_CELL
         do i =1,N_CELL
            ielem = Index_global(i)
             if(ielem .gt. 0) then
             write(180,*),i,ielem,rang_topo,TagParallel(ielem)
             else
             write(180,*),i,ielem,rang_topo
             endif
         enddo
        close(180)
#     endif
!      ________________________________________________________
!     |                                                        |
!     |                    No_cp(N_CELL,3)                     |
!     |________________________________________________________|

!     ______________________________________
!     Independent cell-centers (inside + BC)

      do i=1,N_CELL0
         ielem = Index_global(i)
         countcp = 0
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
                       tag= tag +1
                   endif
                enddo
                if (tag.eq.0) then
                    print*,'Processor:',rang_topo,'Error!!!: Missing cell-centers'
                    print*,'i=',i,',j',j,',No_cp_global=',Jelem
                    stop
                elseif(tag .ge. 2) then
                    countcp = countcp +1
                endif
             endif
         enddo
         if(countcp .gt. 0) print*,'subdomain',rang_topo, 'cp assigned twice:',countcp
      enddo
!     ______________________________________
!     For the ghost elements: No_cp = 0
      do i=N_CELL0+1,N_CELL
         No_cp(i,1) = -1
         No_cp(i,2) = -1
         No_cp(i,3) = -1
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

