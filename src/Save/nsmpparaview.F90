!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!             --   SAVE SOLUTION FOR PARAVIEW DISPLAY  --             !
!                      Miguel Angel Uh Zapata                         !
!                            Feb 2014                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      subroutine outsavParaview(alphaf,uf,wf,pf,                &
                                energftn,dissipftn,effstress,   &
                                alphas,us,ws,ps,                &
                                energstn,dissipstn,             &
                                x,dx,hpr,h,sig,dsig,nrec)

!     *****************************************************************
!     !                                                               !
!     !                         Definitions                           !
!     !                                                               !
!     *****************************************************************
!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

#     include "cppdefs.h"
      use interf
!     ====================================
!     =====  START PARALLEL OPTION =======
#     ifdef KeyParallel
         use parallel
         implicit none
#     else
         implicit none
#        include "common.mpf"
#     endif
!     ====================================
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|

      real*8,dimension(:,:)::alphaf,uf,wf,pf,energftn,dissipftn
      real*8,dimension(:,:)::alphas,us,ws,ps,energstn,dissipstn,effstress
      real*8,dimension(:)  ::x,dx,hpr,h
      real*8,dimension(:)  ::sig,dsig
      integer :: nrec
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|

      real*8,dimension(:,:),allocatable:: xt,zt, &
        ust,wst,pst,energft,dissipft,efft,       &
        uft,wft,pft,energst,dissipst
      real*8,dimension(nxglobal,nzglobal)::xt_gl,zt_gl,               &
        alphaf_gl,uft_gl,wft_gl,pft_gl,energft_gl,dissipft_gl,efft_gl,&
        alphas_gl,ust_gl,wst_gl,pst_gl,energst_gl,dissipst_gl

      character*11 filen
      integer:: irec,i,j

!     *****************************************************************
!     !                                                               !
!     !                        Initialisation                         !
!     !                                                               !
!     *****************************************************************
!      ________________________________________________________
!     |                                                        |
!     |                         Allocate                       |
!     |________________________________________________________|

      allocate(ust(ndimx,ndimz),wst(ndimx,ndimz),pst(ndimx,ndimz),&
             xt(ndimx,ndimz),zt(ndimx,ndimz),energft(ndimx,ndimz),&
             dissipft(ndimx,ndimz),efft(ndimx,ndimz),uft(ndimx,ndimz),&
             wft(ndimx,ndimz),pft(ndimx,ndimz),energst(ndimx,ndimz),&
             dissipst(ndimx,ndimz))

!      ________________________________________________________
!     |                                                        |
!     |                         Formats                        |
!     |________________________________________________________|

 2    format(1(1x,e12.5))
 3    format(3(1x,e12.5))

 100  format('# vtk DataFile Version 2.0')
 110  format('Sample rectilinear grid')
 120  format('ASCII')
 130  format('DATASET RECTILINEAR_GRID')
 140  format('DIMENSIONS ',i5,1x,i5,' 1')
 150  format('X_COORDINATES ',i5,' float')
 160  format('Y_COORDINATES ',i5,' float')
 170  format('Z_COORDINATES 1 float')
 180  format('0.0')
 190  format('POINT_DATA ',i7)

!      ________________________________________________________
!     |                                                        |
!     |                  Physiscal coordinates                 |
!     |________________________________________________________|

      do i=1,nx
         do j=1,nz
            xt(i,j) = x(i)
            zt(i,j) = sig(j)*hpr(i)-h(i)
         enddo
      enddo

!     fluid
      call uvpakepsnodetec(uft,wft,pft, &
                           efft,energft,dissipft,&
                           uf,wf,pf,&
                           effstress,energftn,dissipftn,&
                           xt,zt,dx,dsig)
!     solide
      call uvpakepsnodetec(ust,wst,pst, &
                           efft,energst,dissipst,&
                           us,ws,ps,&
                           effstress,energstn,dissipstn,&
                           xt,zt,dx,dsig)

!      ________________________________________________________
!     |                                                        |
!     |                       File name                        |
!     |________________________________________________________|

      irec=60
      nrec = nrec+1
      filen='P-    .vtk'
      write(filen(3:6),'(i4.4)') nrec

!      ________________________________________________________
!     |                                                        |
!     |                       Parallel                         |
!     |________________________________________________________|

#     ifdef KeyParallel
         call matglo(uft,uft_gl)
         call matglo(wft,wft_gl)
         call matglo(pft,pft_gl)
         call matglo(alphaf,alphaf_gl)
         call matglo(energft,energft_gl)
         call matglo(dissipft,dissipft_gl)
         call matglo(ust,ust_gl)
         call matglo(wst,wst_gl)
         call matglo(pst,pst_gl)
         call matglo(efft,efft_gl)
         call matglo(alphas,alphas_gl)
         call matglo(energst,energst_gl)
         call matglo(dissipst,dissipst_gl)
         call matglo(xt,xt_gl)
         call matglo(zt,zt_gl)

         if(rang_topo.eq.0) then 
            open(irec,FILE=filen)
            write(irec,100)
            write(irec,110)
            write(irec,120)      
            write(irec,130)
            write(irec,140) nxglobal,nzglobal
!           ___________________________________________________
!           Coordinates
            write(irec,150) nxglobal
            do i=1,nxglobal
               write(irec,2) xt_gl(i,1)     
            enddo
            write(irec,160) nzglobal
            do j=1,nzglobal
               write(irec,2) zt_gl(1,j)     
            enddo
            write(irec,170) 
            write(irec,180) 
!           ___________________________________________________
!           Scalars
            write(irec,190) nxglobal*nzglobal
            write(irec,*) 'SCALARS pf     float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) pft_gl(i,j)    
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS alphaf float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) alphaf_gl(i,j)   
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS kf     float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) energft_gl(i,j)  
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS ef     float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) dissipft_gl(i,j) 
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS ps     float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) pst_gl(i,j)
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS alphas float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) alphas_gl(i,j) 
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS ks     float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) energst_gl(i,j)   
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS es     float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) dissipst_gl(i,j)   
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'SCALARS eff    float'
            write(irec,*) 'LOOKUP_TABLE default'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,2) efft_gl(i,j) 
               enddo
            enddo
            write(irec,*) ' '
!           ___________________________________________________
!           Vectors
            write(irec,*) 'VECTORS velocity_f float'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,3) uft_gl(i,j),wft_gl(i,j),0.0d0       
               enddo
            enddo
            write(irec,*) ' '
            write(irec,*) 'VECTORS velocity_s float'
            do j=1,nzglobal
               do i=1,nxglobal
                  write(irec,3) ust_gl(i,j),wst_gl(i,j),0.0d0       
               enddo
            enddo
            rewind(irec)
            close(irec)
         endif
!      ________________________________________________________
!     |                                                        |
!     |                         Serial                         |
!     |________________________________________________________|

#     else
         open(irec,FILE=filen)
         write(irec,100)
         write(irec,110)
         write(irec,120)
         write(irec,130)
         write(irec,140) nx,nz
!        ___________________________________________________
!        Coordinates
         write(irec,150) nx
         do i=1,nx
            write(irec,2) xt(i,1)
         enddo
         write(irec,160) nz
         do j=1,nz
            write(irec,2) zt(1,j)
         enddo
         write(irec,170)
         write(irec,180)
!        ___________________________________________________
!        Scalars
         write(irec,190) nx*nz
         write(irec,*) 'SCALARS pf     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) pft(i,j)    
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS alphaf float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) alphaf(i,j)   
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS kf     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) energft(i,j)  
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS ef     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) dissipft(i,j) 
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS ps     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) pst(i,j)
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS alphas float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) alphas(i,j) 
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS ks     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) energst(i,j)   
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS es     float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) dissipst(i,j)   
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'SCALARS eff    float'
         write(irec,*) 'LOOKUP_TABLE default'
         do j=1,nz
            do i=1,nx
               write(irec,2) efft(i,j) 
            enddo
         enddo
         write(irec,*) ' '
!        ___________________________________________________
!        Vectors
         write(irec,*) 'VECTORS velocity_f float'
         do j=1,nz
            do i=1,nx
               write(irec,3) uft(i,j),wft(i,j),0.0d0       
            enddo
         enddo
         write(irec,*) ' '
         write(irec,*) 'VECTORS velocity_s float'
         do j=1,nz
            do i=1,nx
               write(irec,3) ust(i,j),wst(i,j),0.0d0       
            enddo
         enddo
         rewind(irec)
         close(irec)
#     endif

!     *****************************************************************
!     !                                                               !
!     !                           Finalization                        !
!     !                                                               !
!     *****************************************************************
!      ________________________________________________________
!     |                                                        |
!     |                       Deallocate                       |
!     |________________________________________________________|

      deallocate(ust,wst,pst,xt,zt,energft,dissipft,&
                 efft,uft,wft,pft,energst,dissipst)

      return
      end

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                                 End                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
