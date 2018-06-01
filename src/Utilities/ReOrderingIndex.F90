!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	     INPUT FILE                               !
!                             Jul 2014                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!


      SUBROUTINE ReOrderingIndex(No_vp_new,No_cp_new,nbe_new)                                       

!---------------------------------------------------------------------!
!                                                                     !
!    This program changes the order of the mesh indexing given        !
!    by BlueKenue program according to the technique of block         !
!    domain decomposition made by chaco program.                      !   
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
!    Output variables:                                                !
!   _______________________________________________________________   !
!  |     Name    |    Size   | Description                         |  !  
!  |_____________|___________|_____________________________________|  !   
!  |             |           |                                     |  !
!  | <-- No_vp   |(N_CELL0,3)| Index number of the 3 cell vertices |  !
!  | <-- No_cp   |(N_CELL,3) | Index number of surrounding 3 cells |  !
!  | <-- nbe     | N_CELL0   | Type of boundary cell (inside or bc)|  !
!  |_____________|___________|_____________________________________|  !
!                                                                     !
!                                                                     !
!    Local parameters & variables:                                    !
!   _______________________________________________________________   !
!  |   Name      |                 Description                     |  !  
!  |_____________|_________________________________________________|  !           
!  |_____________|_________________________________________________|  !
!                                                                     !
!   -->  Input variables                                              !
!   <--  Output variables                                             !
!   ---  Parameters                                                   !
!    -   Common variables used                                        !
!    *   Common variables modified                                    !
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
            !--------------            
            integer,dimension(:,:):: No_vp_new(N_CELL0,3)
            integer,dimension(:,:):: No_cp_new(N_CELL,3)
            integer,dimension(:)  :: nbe_new(N_CELL0)
            !--------------
            integer,dimension(:,:) :: No_vp_old(N_CELL0,3)
            integer,dimension(:,:) :: No_cp_old(N_CELL,3)
            integer,dimension(:)   :: nbe_old(N_CELL0)
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
            USE parallel
            USE geometry
            implicit none
            !--------------            
            integer,dimension(:,:):: No_vp_new(N_CELL0global,3)
            integer,dimension(:,:):: No_cp_new(N_CELL0global,3)
            integer,dimension(:)  :: nbe_new(N_CELL0global)
            !--------------
            integer,dimension(:,:) :: No_vp_old(N_CELL0global,3)
            integer,dimension(:,:) :: No_cp_old(N_CELL0global,3)
            integer,dimension(:)   :: nbe_old(N_CELL0global)
#     endif
!     =============== END ================    
!     ====================================  
!     ____________________________________
!    |                                    |
!    |   Declaration of local variables   |
!    |____________________________________|

      integer,dimension(:),  allocatable :: TagSubdomain
      integer,dimension(:),  allocatable :: NewIndex,OldIndex      
      integer,dimension(:),  allocatable :: IniSub,FinSub
      integer,dimension(:),  allocatable :: NumSubBdy,IniSubBdy,FinSubBdy
      integer,dimension(:),  allocatable :: NumSubInt,IniSubInt,FinSubInt
      integer,dimension(:),  allocatable :: ElementsBdy
      integer,dimension(:),  allocatable :: ElementsInt
      integer :: NCELLS
      integer :: NumberSubdomains,tag,i0,s,s0,sm,j0
      integer :: Ielem,Jelem,ic1,ic2,ic3,m,ii
      integer :: TotalElementsBdy,TotalElementsInt
!     ---------------
      integer :: k1,k2,k3,kv1,kv2,kk
      real*8 :: funh
      character*80 title

!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*, '      ----> Begin subroutine: ReOrderingIndex'
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!*********************************************************************!
!                                                                     !
!                            Inicialization                           !
!                                                                     !
!*********************************************************************!

!     ====================================
!     ==========  SEQUENTIAL =============
#     ifndef KeyParallel
         NCELLS = N_CELL0
         print*,' '
         print*,'     -----------------------------------------------'
         print*,'               *** Index re-ordering ***            '
         print*,'     -----------------------------------------------'
         print*,' '
#     endif
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         NCELLS = N_CELL0global
         if (rang_topo.eq.0) then
         print*,' '
         print*,'     -----------------------------------------------'
         print*,'             *** MPI: Index re-ordering ***         '
         print*,'     -----------------------------------------------'
         print*,' '
         endif
#     endif
!     =============== END ================    
!     ====================================
	
!     ________________________________________________________
!     Allocate

      allocate(TagSubdomain(NCELLS))
      allocate(NewIndex(NCELLS)) 
      allocate(OldIndex(NCELLS)) 
      allocate(IniSub(NumberSubdomains)) 
      allocate(FinSub(NumberSubdomains)) 
      allocate(ElementsBdy(NCELLS),&
               ElementsInt(NCELLS))
      allocate(NumSubBdy(NumberSubdomains),&
               IniSubBdy(NumberSubdomains),&
               FinSubBdy(NumberSubdomains)) 
      allocate(NumSubInt(NumberSubdomains),&
               IniSubInt(NumberSubdomains),&
               FinSubInt(NumberSubdomains)) 

89    format(I8)

!      ________________________________________________________
!     |                                                        |
!     |          Index re-ordering according to Chaco          |
!     |________________________________________________________|
	
!     ________________________________________________________
!     Tags corresponding for each subdomain
      open(12,file='dataChaco.txt',status='old')
      do i=1,NCELLS
         read(12,89) TagSubdomain(i)
      enddo
      close(12)
!     ________________________________________________________
!     Number of subdomains
      NumberSubdomains = 1
      do i=2,NCELLS
         tag = 0
         do j=1,i-1
             if (TagSubdomain(i).eq.TagSubdomain(j)) tag=1
         enddo
         if (tag.eq.0) NumberSubdomains = NumberSubdomains + 1
      enddo
!     ________________________________________________________
!     New cell-center indexes 
      i0 = 0
      do s=1,NumberSubdomains
         IniSub(s) = i0+1 !<----Initial index (subdomain)
         do i=1,NCELLS
            if (TagSubdomain(i).eq.(s-1)) then
                i0 = i0 + 1
                NewIndex(i)  = i0 
                OldIndex(i0) = i 
            endif
         enddo
         FinSub(s) = i0   !<---- Final index (subdomain)
      enddo
!     ________________________________________________________
!     New index asignation to each element

      No_vp_old = No_vp_new 
      No_cp_old = No_cp_new
      nbe_old   = nbe_new
      do i0=1,NCELLS
         i = OldIndex(i0)
         do j=1,3
            No_vp_new(i0,j) = No_vp_old(i,j)
            No_cp_new(i0,j) = NewIndex(No_cp_old(i,j))
         enddo
         nbe_new(i0) = nbe_old(i)
      enddo  
!      ________________________________________________________
!     |                                                        |
!     |     Index re-ordering according to boundary elemets    |
!     |________________________________________________________|
	
!     _________________________________________________________
!     Boundary cell-centers                
      j0 = 0
      i0 = 0
      do s=1,NumberSubdomains-1           
         IniSubBdy(s) = j0 + 1     !<----Initial index boundary (subdomain)
         IniSubInt(s) = i0 + 1     !<----Initial index interior (subdomain)
!        -------------------------------------- 
         do i=IniSub(s),FinSub(s)   
            ic1 = No_cp_new(i,1)
            ic2 = No_cp_new(i,2)
            ic3 = No_cp_new(i,3)
            tag = 0
            do sm=s+1,NumberSubdomains  !<--- Only the subdomains with bigger indices       
               do j=IniSub(sm),FinSub(sm)
                  if ((j.eq.ic1).or.(j.eq.ic2).or.(j.eq.ic3)) then
                     tag = 1
                  endif
               enddo
            enddo
            if (tag.eq.1) then 
               j0 = j0 + 1
               ElementsBdy(j0) = i    ! index of Subdomain
            else
               i0 = i0 + 1
               ElementsInt(i0) = i     
            endif
         enddo
!        --------------------------------------
         FinSubBdy(s) = j0       !<---- Final index boundary (subdomain)
         FinSubInt(s) = i0       !<---- Final index interior (subdomain)  
      enddo

      s=NumberSubdomains
!     -----------------
      IniSubBdy(s) = 0
      FinSubBdy(s) = 0
!     -----------------
      IniSubInt(s) = i0 + 1
      do i=IniSub(s),FinSub(s) 
         i0 = i0 + 1
         ElementsInt(i0) = i
      enddo
      FinSubInt(s) = i0
!     -----------------
      TotalElementsBdy = j0
      TotalElementsInt = i0

!     ________________________________________________________
!     Re-index considering boundary elements

!     --------------------------
!     New index: boundary first
      do i0=1,TotalElementsBdy
         i = ElementsBdy(i0)
         OldIndex(i0) = i
         NewIndex(i)  = i0 
      enddo 
!     --------------------------
!     New index: internal second
      do ii=1,TotalElementsInt
         i0 = ii + TotalElementsBdy
         i  = ElementsInt(ii)
         OldIndex(i0) = i
         NewIndex(i)  = i0
      enddo 
!     --------------------------
!     Re-index
      No_vp_old = No_vp_new 
      No_cp_old = No_cp_new
      nbe_old   = nbe_new
      do i0=1,NCELLS
         i = OldIndex(i0)
         do j=1,3
            No_vp_new(i0,j) = No_vp_old(i,j)
            No_cp_new(i0,j) = NewIndex(No_cp_old(i,j))
         enddo
         nbe_new(i0) = nbe_old(i)
      enddo  

!     ________________________________________________________
!     Deallocate

      deallocate(TagSubdomain,&
                 NewIndex,OldIndex)
      deallocate(IniSub,FinSub,&
                 ElementsBdy,ElementsInt,&
                 NumSubBdy,IniSubBdy,FinSubBdy,&
                 NumSubInt,IniSubInt,FinSubInt)

!*********************************************************************!
!                                                                     !
!                            Finalization                             !
!                                                                     !
!*********************************************************************!
	
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
           print*,'      <----   End subroutine: ReOrderingIndex'
           print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	RETURN
	END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	   END OF INPUT                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

