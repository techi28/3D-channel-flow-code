!    -----------------------------------------------------------
     call fvm_gradient(phiu,phiv,phiw, &
                  dudx,dudy,dudz, &
                  dvdx,dvdy,dvdz, &
                  dwdx,dwdy,dwdz, &
                  xc,yc,sig,dsig,No_cp,nbe,  &
                  xv,yv,sigv,dsigv,No_vp,nbev)
!   ------------------------------------------------------------

     do k=2,NZ-1
        do i=1,N_CELL0
           coef = Gamx(i,k)
           Vol = areaCell(i)*disgv(k-1)
           AB(1) = 0.0d0
           AB(2) = 0.0d0
           AB(3) = 0.0d0
!         __________________________________________
!         For XY Plane
           do j =1, 3
               jc = No_cp(i,j)
!             -------------
!              aa(j) refers to stress tensor Tij
               AA(1) = 2.0d0*(phiu(jc,k)-phiu(i,k))/(xc(jc)-xc(i))
               AA(2) = (phiu(jc,k)-phiu(i,k))/(yc(jc)-yc(i)) + & 
                       (phiv(jc,k)-phiv(i,k))/(xc(jc)-xc(i))
               AB(j) = AA(1)*xnn(i,j)+AA(2)*ynn(i,j)
               AB(j) = AB(j)*dlVV(i,j)*dsigv(k-1)
               dwdxB = 0.5d0*(phiw(i,k)+phiw(i,k-1))*
           enddo
!         __________________________________________
!         Bottom
                        phiw
                 AAB = (phiu(i,k)-phiu(i,k-1))/dsigv(k-1) + &
                        dwdxB
                 AAB = -1.0d0*AAB*areaCell(i)
!         __________________________________________
!         Top
                 AAT = (phiu(i,k+1)-phiu(i,k)/dsigv(k) + &
                       0.5d0*(dwdx(i,k) + dwdx(i,k+1))
                 AAT = 1.0d0*AAT*areaCell(i)
!         ___________________________________________
!         total coef for centre

                 AA0 = AB(1)+AB(2)+AB(3)+AAB+AAT
                 AA0 = coef*AA0/Vol
                 rhsu(i,k) = rhsu(i,k)-AA0
          enddo
        enddo


     do k=2,NZ-1
        do i=1,N_CELL0
           coef = Gamx(i,k)
           Vol = areaCell(i)*disgv(k-1)
           AB(1) = 0.0d0
           AB(2) = 0.0d0
           AB(3) = 0.0d0
!         __________________________________________
!         For XY Plane
           do j =1, 3
               jc = No_cp(i,j)
!             -------------
!              aa(j) refers to stress tensor Tij
               AA(1) = (phiv(jc,k)-phiv(i,k))/(xc(jc)-xc(i)) + &
                       (phiu(jc,k)-phiu(i,k))/(yc(jc)-yc(i))
               AA(2) = 2.0d0*(phiv(jc,k)-phiv(i,k))/(yc(jc)-yc(i))
               AB(j) = AA(1)*xnn(i,j)+AA(2)*ynn(i,j)
               AB(j) = AB(j)*dlVV(i,j)*dsigv(k-1)
           enddo
!         __________________________________________
!         Bottom
                 AAB = (phiv(i,k)-phiv(i,k-1))/dsigv(k-1) + &
                        0.5d0*(dwdy(i,k) + dwdy(i,k-1))
                 AAB = -1.0d0*AAB*areaCell(i)
!         __________________________________________
!         Top
                 AAT = (phiv(i,k+1)-phiv(i,k)/dsigv(k) + &
                       0.5d0*(dwdy(i,k) + dwdy(i,k+1))
                 AAT = 1.0d0*AAT*areaCell(i)
!         ___________________________________________
!         total coef for centre

                 AA0 = AB(1)+AB(2)+AB(3)+AAB+AAT
                 AA0 = coef*AA0/Vol
                 rhsv(i,k) = rhsu(i,k)-AA0
          enddo
        enddo


     do k=2,NZ-1
        do i=1,N_CELL0
           coef = Gamx(i,k)
           Vol = areaCell(i)*disgv(k-1)
           AB(1) = 0.0d0
           AB(2) = 0.0d0
           AB(3) = 0.0d0
!         __________________________________________
!         For XY Plane
           do j =1, 3
               jc = No_cp(i,j)
!             -------------
!              aa(j) refers to stress tensor Tij
               AA(1) = (phiw(jc,k)-phiw(i,k))/(xc(jc)-xc(i)) + &
                       dudz
               AA(2) = (phiw(jc,k)-phiw(i,k))/(yc(jc)-yc(i)) + &
                       dvdz
               AB(j) = AA(1)*xnn(i,j)+AA(2)*ynn(i,j)
               AB(j) = AB(j)*dlVV(i,j)*dsigv(k-1)
           enddo
!         __________________________________________
!         Bottom
                 AAB = 2.0d0*(phiw(i,k)-phiw(i,k-1))/dsigv(k-1)
                 AAB = -1.0d0*AAB*areaCell(i)
!         __________________________________________
!         Top
                 AAT = 2.0d0*(phiw(i,k+1)-phiw(i,k)/dsigv(k) 
                 AAT = 1.0d0*AAT*areaCell(i)
!         ___________________________________________
!         total coef for centre

                 AA0 = AB(1)+AB(2)+AB(3)+AAB+AAT
                 AA0 = coef*AA0/Vol
                 rhsw(i,k) = rhsw(i,k)-AA0
          enddo
        enddo
