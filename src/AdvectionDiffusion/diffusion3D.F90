!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                       DIFFUSION 3D Stress                           !
!                            Apr 2016                                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE diffusion3DNEW(rhsu,rhsv,rhsw,                   &
                                phiu,phiv,phiw,                   &
                                phiuv,phivv,phiwv,                &
                                Gamx,Gamy,Gamz,                   &
                                xc,yc,sig,dsig,No_cp,nbe,         &
                                xv,yv,sigv,dsigv,No_vp,nbev,Option)

!*********************************************************************!
!                                                                     !
!                           Definitions                               !
!                                                                     !
!*********************************************************************!
!      ____________________________________
!     |                                    |
!     |     Keys and common parameters     |
!     |____________________________________|

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
!      ____________________________________
!     |                                    |
!     |      Declaration of variables      |
!     |____________________________________|

      real*8, dimension(:,:) :: phiu(N_CELL,NZ)
      real*8, dimension(:,:) :: phiv(N_CELL,NZ)
      real*8, dimension(:,:) :: phiw(N_CELL,NZ)
      real*8, dimension(:,:) :: phiuv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: phivv(N_VERT,NZ-1)
      real*8, dimension(:,:) :: phiwv(N_VERT,NZ-1)
!     --------------------------------------
      real*8, dimension(:,:) :: rhsu(N_CELL,NZ)
      real*8, dimension(:,:) :: rhsv(N_CELL,NZ)
      real*8, dimension(:,:) :: rhsw(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:,:) :: Gamx(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamy(N_CELL,NZ)
      real*8, dimension(:,:) :: Gamz(N_CELL,NZ)
!     --------------------------------------
      real*8, dimension(:)   :: xc(N_CELL)
      real*8, dimension(:)   :: yc(N_CELL)
      real*8, dimension(:)   :: sig(NZ)
      real*8, dimension(:)   :: dsig(NZ)
      integer,dimension(:,:) :: No_cp(N_CELL,3)
      integer,dimension(:)   :: nbe(N_CELL0)
!     --------------------------------------
      real*8, dimension(:)   :: xv(N_VERT)
      real*8, dimension(:)   :: yv(N_VERT)
      real*8, dimension(:)   :: sigv(NZ-1)
      real*8, dimension(:)   :: dsigv(NZ-1)
      integer,dimension(:,:) :: No_vp(N_CELL0,3)
      integer,dimension(:)   :: nbev(N_VERT)
      integer :: Option
!      ____________________________________
!     |                                    |
!     |   Declaration of local variables   |
!     |____________________________________|
      real*8, dimension(:,:),allocatable :: dudx,dudy,duds
      real*8, dimension(:,:),allocatable :: dvdx,dvdy,dvds
      real*8, dimension(:,:),allocatable :: dwdx,dwdy,dwds
      real*8, dimension(:) :: dudn(1:3),dvdn(1:3),dwdn(1:3)
      real*8 :: dudT,dudB,dvdT,dvdB,dwdT,dwdB
      real*8 :: funA,funB,co1,co2,co3
      real*8 :: diff1,diff2,diff3,coef
      real*8 :: sourceu,sourcev,sourcew,Vol
      integer :: jc,jv1,jv2,jj
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     ifdef KeyDbg
         write(*,'(t22,60a)'), '----> Begin subroutine: diffusion3D'
         print*,' '
#     endif
!     ------------------------------------------------------------

        allocate(dudx(N_CELL,NZ),dudy(N_CELL,NZ), duds(N_CELL,NZ), &
                 dvdx(N_CELL,NZ),dvdy(N_CELL,NZ), dvds(N_CELL,NZ), &
                 dwdx(N_CELL,NZ),dwdy(N_CELL,NZ), dwds(N_CELL,NZ))

        call grandientLSM(dudx,dudy,duds,phiu,xc,yc,sig,dsig,No_cp,nbe)
        call grandientLSM(dvdx,dvdy,dvds,phiv,xc,yc,sig,dsig,No_cp,nbe)
        call grandientLSM(dwdx,dwdy,dwds,phiw,xc,yc,sig,dsig,No_cp,nbe)


#     ifdef KeyKimChoi
!      ________________________________________________________
!     |                                                        |
!     |                     Main routine                       |
!     |________________________________________________________|

     do k=2, NZ-1
        do i=1, N_CELL0
          Vol = areaCell(i)*dsigv(k-1)
!         ----------------------
!         HORIZONTAL FACE
          do j=1,3
            jc = No_cp(i,j)
            jj = j+1
            if(jj .gt. 3)  jj= jj-3
            jv1 = No_vp(i,j)
            jv2 = No_vp(i,jj)
            co1 = dle1(i,j)
            co2 = dle2(i,j)
            co3 = dle3(i,j)
!           -----------------
!           dudn
            diff1 = (phiu(jc,k) -phiu(i,k))/dlCC(i,j)
            funA = phiuv(jv1,k)
            funB = phiuv(jv2,k)
            diff2 = (funB-funA)/dlVV(i,j)
            diff3 = (xmo(i,j)-xme(i,j))*(dudx(jc,k)-dudx(i,k)) &
                   +(ymo(i,j)-yme(i,j))*(dudy(jc,k)-dudy(i,k))
            dudn(j) = diff1*co1 + diff2*co2 + diff3*co3
!           ----------------
!           dvdn
            diff1 =  (phiv(jc,k) -phiv(i,k))/dlCC(i,j)
            funA = phivv(jv1,k)
            funB = phivv(jv2,k)
            diff2 = (funB-funA)/dlVV(i,j)
            diff3 = (xmo(i,j)-xme(i,j))*(dvdx(jc,k)-dvdx(i,k)) &
                   +(ymo(i,j)-yme(i,j))*(dvdy(jc,k)-dvdy(i,k))
            dvdn(j) = diff1*co1 + diff2*co2 + diff3*co3
!           ----------------
!           dwdn
            diff1 =  (phiw(jc,k) -phiw(i,k))/dlCC(i,j)
            funA = phiwv(jv1,k)
            funB = phiwv(jv2,k)
            diff2 = (funB-funA)/dlVV(i,j)
            diff3 = (xmo(i,j)-xme(i,j))*(dwdx(jc,k)-dwdx(i,k)) &
                   +(ymo(i,j)-yme(i,j))*(dwdy(jc,k)-dwdy(i,k))
            dwdn(j) = diff1*co1 + diff2*co2 + diff3*co3
          enddo
!         ----------------------
!         HORIZONTAL FACE
           IF(I3D .EQ. 1) THEN
           dudT = (phiu(i,k+1)-phiu(i,k))/dsig(k)
           dudB = (phiu(i,k-1)-phiu(i,k))/dsig(k-1)
           dVdT = (phiv(i,k+1)-phiv(i,k))/dsig(k)
           dVdB = (phiv(i,k-1)-phiv(i,k))/dsig(k-1)
           dwdT = (phiw(i,k+1)-phiw(i,k))/dsig(k)
           dwdB = (phiw(i,k-1)-phiw(i,k))/dsig(k-1)
           ELSE
            dudT = 0.
            dudB = 0.
            dVdT = 0.
            dVdB = 0.
            dwdT = 0.
            dwdB = 0.
           ENDIF
!          -------------------
!          calculate extra source term by mass flux
           sourceu =  dudn(1)*dlVV(i,1)*dsigv(k-1)+ &
                      dudn(2)*dlVV(i,2)*dsigv(k-1)+ &
                      dudn(3)*dlVV(i,3)*dsigv(k-1)+ &
                      (dudT+dudB)*areaCell(i)

           sourcev =  dvdn(1)*dlVV(i,1)*dsigv(k-1)+ &
                      dvdn(2)*dlVV(i,2)*dsigv(k-1)+ &
                      dvdn(3)*dlVV(i,3)*dsigv(k-1)+ &
                      (dvdT+dvdB)*areaCell(i)

           sourcew =  dwdn(1)*dlVV(i,1)*dsigv(k-1)+ &
                      dwdn(2)*dlVV(i,2)*dsigv(k-1)+ &
                      dwdn(3)*dlVV(i,3)*dsigv(k-1)+ &
                      (dwdT+dwdB)*areaCell(i)
!          -----------------------------------------------
!          update the right hand side term
           rhsu(i,k) = rhsu(i,k) + Gamx(i,k)*sourceu/Vol
           rhsv(i,k) = rhsv(i,k) + Gamy(i,k)*sourcev/Vol
           rhsw(i,k) = rhsw(i,k) + Gamz(i,k)*sourcew/Vol
         enddo
        enddo
#     else
!*********************************************************************!
!                                                                     !
!                  Initialization of the subroutine                   !
!                                                                     !
!*********************************************************************!
!     Option = 1 ! for explicit scheme
!     Option = 2 ! for semi-implicit scheme
!      ________________________________________________________
!     |                                                        |
!     |                     Main routine                       |
!     |________________________________________________________|

     do k=2, NZ-1
        do i=1, N_CELL0
            do j=1,3
                jc = No_cp(i,j)
                co1 = dsqrt(xe1(i,j)*xe1(i,j)+ye1(i,j)*ye1(i,j))/dlCC(i,j)
                co2 = xe2(i,j)
                co3 = ye2(i,j)

                if(Option .eq. 2) then
                    co1 = 0.5d0*co1
                    co2 = co2
                    co3 = co3
                endif
!              -----------------
!              dudn
                diff1 = phiu(jc,k) -phiu(i,k)
                funA = dudx(i,k)
                funB = dudx(jc,k)
                diff2 = 0.5d0*(funB+funA)
                funA = dudy(i,k)
                funB = dudy(jc,k)
                diff3 = 0.5d0*(funB+funA)
                dudn(j) = diff1*co1 + diff2*co2 + diff3*co3
!               ----------------
!               dvdn
                diff1 = phiv(jc,k) -phiv(i,k)
                funA = dvdx(i,k)
                funB = dvdx(jc,k)
                diff2 = 0.5d0*(funB+funA)
                funA = dvdy(i,k)
                funB = dvdy(jc,k)
                diff3 = 0.5d0*(funB+funA)
                dvdn(j) = diff1*co1 + diff2*co2 + diff3*co3
!               ----------------
!                dwdn
                diff1 = phiw(jc,k) -phiw(i,k)
                funA = dwdx(i,k)
                funB = dwdx(jc,k)
                diff2 = 0.5d0*(funB+funA)
                funA = dwdy(i,k)
                funB = dwdy(jc,k)
                diff3 = 0.5d0*(funB+funA)
                dwdn(j) = diff1*co1 + diff2*co2 + diff3*co3
            enddo
!         ----------------------
!         VERTICAL FACE
           if(Option .eq. 2) then
             coef = 0.5d0
           else
             coef = 1.0d0
           endif

           if(I3D .ne. 1) coef = 0.

           dudT = coef*(phiu(i,k+1)-phiu(i,k))/dsig(k)
           dudB = coef*(phiu(i,k-1)-phiu(i,k))/dsig(k-1)
           dVdT = coef*(phiv(i,k+1)-phiv(i,k))/dsig(k)
           dVdB = coef*(phiv(i,k-1)-phiv(i,k))/dsig(k-1)
           dwdT = coef*(phiw(i,k+1)-phiw(i,k))/dsig(k)
           dwdB = coef*(phiw(i,k-1)-phiw(i,k))/dsig(k-1)
!          -------------------
!          calculate extra source term by mass flux
           sourceu =  dudn(1)*dlVV(i,1)*dsigv(k-1)+ &
                      dudn(2)*dlVV(i,2)*dsigv(k-1)+ &
                      dudn(3)*dlVV(i,3)*dsigv(k-1)+ &
                      (dudT+dudB)*areaCell(i)

           sourcev =  dvdn(1)*dlVV(i,1)*dsigv(k-1)+ &
                      dvdn(2)*dlVV(i,2)*dsigv(k-1)+ &
                      dvdn(3)*dlVV(i,3)*dsigv(k-1)+ &
                      (dvdT+dvdB)*areaCell(i)

           sourcew =  dwdn(1)*dlVV(i,1)*dsigv(k-1)+ &
                      dwdn(2)*dlVV(i,2)*dsigv(k-1)+ &
                      dwdn(3)*dlVV(i,3)*dsigv(k-1)+ &
                      (dwdT+dwdB)*areaCell(i)
!          -----------------------------------------------
!          update the right hand side term
           Vol = areaCell(i)*dsigv(k-1)
           rhsu(i,k) = rhsu(i,k) + Gamx(i,k)*sourceu/Vol
           rhsv(i,k) = rhsv(i,k) + Gamy(i,k)*sourcev/Vol
           rhsw(i,k) = rhsw(i,k) + Gamz(i,k)*sourcew/Vol
         enddo
        enddo
#     endif

!*********************************************************************!
!                                                                     !
!                   Finalization of the subroutine                    !
!                                                                     !
!*********************************************************************!
         deallocate(dudx,dudy,duds,dvdx,dvdy,dvds,dwdx,dwdy,dwds)

#     ifdef KeyDbg
         write(*,'(t22,60a)'), '<---- End   subroutine: diffusion3DNew'
         print*,' '
#     endif
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                      	END OF DIFFUSION3DNEW                         !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

