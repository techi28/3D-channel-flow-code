!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        MODULE: interfaces                           !
!                 Declaration of all the subroutines                  !
!                             Jan 2013                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      MODULE INTERFACES
      interface

!---------------------------------------------------------------------!
!     This MODULE contains the following subroutines:                 !
!       1)  input.F90                                                 !
!---------------------------------------------------------------------!

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                             Input.F90 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE input_data(No_vp,No_cp,No_wb,No_hb,No_qb,No_sp,&
                            nbe,xv,yv,zbv)                  

      integer,dimension(:,:):: No_vp
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: No_wb 
      integer,dimension(:)  :: No_qb
      integer,dimension(:)  :: No_hb
      integer,dimension(:)  :: No_sp
      integer,dimension(:)  :: nbe 
      real*8 ,dimension(:)  :: xv
      real*8 ,dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv


     END SUBROUTINE input_data

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                      input_parameters.F90 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE input_parameters(nbev,hv,h,                     &
                                  No_vp,No_cp,No_wb,No_hb,No_qb, &
                                  nbe,xv,yv,zbv)                  

      integer,dimension(:)  :: nbev    
      real*8 ,dimension(:)  :: h
      real*8 ,dimension(:)  :: hv
      integer,dimension(:,:):: No_vp
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: No_wb 
      integer,dimension(:)  :: No_qb
      integer,dimension(:)  :: No_hb
      integer,dimension(:)  :: nbe   
      real*8 ,dimension(:)  :: xv
      real*8 ,dimension(:)  :: yv
      real*8 ,dimension(:)  :: zbv

     END SUBROUTINE input_parameters


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                      parallel_transfer.F90 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE parallel_input_local(No_vp,No_cp,No_wb,No_qb,No_hb,No_sp, &
                                      nbe,xv,yv,zbv)

      integer,dimension(:,:):: No_cp
      integer,dimension(:,:):: No_vp  
      integer,dimension(:)  :: No_wb 
      integer,dimension(:)  :: No_qb
      integer,dimension(:)  :: No_hb 
      integer,dimension(:)  :: No_sp
      integer,dimension(:)  :: nbe  
      real*8, dimension(:)  :: xv
      real*8, dimension(:)  :: yv  
      real*8 ,dimension(:)  :: zbv

     END SUBROUTINE parallel_input_local

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                             geometry.F90 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE calcul_geometry(xc,yc,sig,dsig,No_cp,nbe,h,  &
                                 xv,yv,sigv,dsigv,No_vp,nbev, &
                                 ic1tec,ic2tec,ic3tec)

      real*8, dimension(:)  :: xc
      real*8, dimension(:)  :: yc  
      real*8, dimension(:)  :: sig
      real*8, dimension(:)  :: dsig 
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe   
      real*8, dimension(:)  :: h      
      real*8, dimension(:)  :: xv
      real*8, dimension(:)  :: yv 
      real*8, dimension(:)  :: sigv
      real*8, dimension(:)  :: dsigv
      integer,dimension(:,:):: No_vp
      integer,dimension(:)  :: nbev  
      integer,dimension(:)  :: ic1tec
      integer,dimension(:)  :: ic2tec
      integer,dimension(:)  :: ic3tec
 
      END SUBROUTINE calcul_geometry

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                             initial.F90                             !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

      SUBROUTINE input_initial(alphafn,ufn,vfn,wfn,pfn,viscof,rhof,  &
                         alphasn,usn,vsn,wsn,psn,viscos,rhos,        &
                         alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv,      &
                         alphasv,usv,vsv,wsv,psv,viscosv,rhosv,      &
                         etan,etav,                                  &
                         Hpr,Hprv,                                   &
                         h,hv,                                       &
                         xc,yc,sig,dsig,No_cp,nbe,                   &
                         xv,yv,sigv,dsigv,No_vp,nbev,                &
                         Heaviside,mask)

      real*8,dimension(:,:) :: alphafn
      real*8,dimension(:,:) :: ufn,vfn,wfn,pfn
      real*8,dimension(:,:) :: rhof,viscof      
      real*8,dimension(:,:) :: alphasn
      real*8,dimension(:,:) :: usn,vsn,wsn,psn
      real*8,dimension(:,:) :: rhos,viscos
      real*8,dimension(:,:) :: alphafv,ufv,vfv,wfv,pfv,viscofv,rhofv
      real*8,dimension(:,:) :: alphasv,usv,vsv,wsv,psv,viscosv,rhosv
      real*8,dimension(:)   :: etan,etav   
      real*8,dimension(:)   :: Hpr,Hprv
      real*8,dimension(:)   :: h,hv   
      real*8,dimension(:)   :: xc,yc,sig,dsig
      real*8,dimension(:)   :: xv,yv,sigv,dsigv
      integer,dimension(:,:):: No_cp,No_vp
      integer,dimension(:)  :: nbe,nbev
      real*8,dimension(:)   :: Heaviside
      real*8,dimension(:,:) :: mask

      END SUBROUTINE input_initial


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                           Re-start & saved data                     !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |    restart_out.F90                                     |
!     |________________________________________________________|

      SUBROUTINE restart_out(alphafnp,ufnp,vfnp,wfnp,pfnp, & 
                             alphasnp,usnp,vsnp,wsnp,psnp, &
                             zct)

      real*8,dimension(:,:):: alphafnp
      real*8,dimension(:,:):: ufnp
      real*8,dimension(:,:):: vfnp
      real*8,dimension(:,:):: wfnp
      real*8,dimension(:,:):: pfnp   
      real*8,dimension(:,:):: alphasnp
      real*8,dimension(:,:):: usnp
      real*8,dimension(:,:):: vsnp
      real*8,dimension(:,:):: wsnp
      real*8,dimension(:,:):: psnp
      real*8,dimension(:,:):: zct

      END SUBROUTINE restart_out

!      ________________________________________________________
!     |                                                        |
!     |    restart_in.F90                                      |
!     |________________________________________________________|

      SUBROUTINE restart_in(alphafn,ufn,vfn,wfn,pfn, & 
                            alphasn,usn,vsn,wsn,psn, &
                            zct)

      real*8,dimension(:,:):: alphafn
      real*8,dimension(:,:):: ufn
      real*8,dimension(:,:):: vfn
      real*8,dimension(:,:):: wfn
      real*8,dimension(:,:):: pfn   
      real*8,dimension(:,:):: alphasn
      real*8,dimension(:,:):: usn
      real*8,dimension(:,:):: vsn
      real*8,dimension(:,:):: wsn
      real*8,dimension(:,:):: psn
      real*8,dimension(:,:):: zct

      END SUBROUTINE restart_in

!      ________________________________________________________
!     |                                                        |
!     |    SavetecVertex.F90                                   |
!     |________________________________________________________|

      SUBROUTINE SavetecVertex(alphafv,ufv,vfv,wfv,pfv,       &
                               alphasv,usv,vsv,wsv,psv,rhosv, &
                               xvt,yvt,zvt,No_vp)

      real*8,dimension(:,:) :: alphafv
      real*8,dimension(:,:) :: ufv
      real*8,dimension(:,:) :: vfv
      real*8,dimension(:,:) :: wfv
      real*8,dimension(:,:) :: pfv  
      real*8,dimension(:,:) :: alphasv
      real*8,dimension(:,:) :: usv
      real*8,dimension(:,:) :: vsv
      real*8,dimension(:,:) :: wsv
      real*8,dimension(:,:) :: psv
      real*8,dimension(:,:) :: rhosv
      real*8,dimension(:,:) :: xvt
      real*8,dimension(:,:) :: yvt
      real*8,dimension(:,:) :: zvt
      integer,dimension(:,:):: No_vp

      END SUBROUTINE SavetecVertex

!      ________________________________________________________
!     |                                                        |
!     |    SaveParaviewVertex.F90                              |
!     |________________________________________________________|

      SUBROUTINE SaveParaviewVertex(Hprv,etav,                 &
                                    ufv,vfv,wfv,pfv,           &
                                    usv,vsv,wsv,psv,           &
                                    uErrv,vErrv,wErrv,pErrv,   &
                                    xvt,yvt,zvt,               &
                                    No_vp)

      real*8,dimension(:)   :: Hprv
      real*8,dimension(:)   :: etav
      
      real*8,dimension(:,:) :: ufv
      real*8,dimension(:,:) :: vfv
      real*8,dimension(:,:) :: wfv
      real*8,dimension(:,:) :: pfv  

      real*8,dimension(:,:) :: usv
      real*8,dimension(:,:) :: vsv
      real*8,dimension(:,:) :: wsv
      real*8,dimension(:,:) :: psv
      
      real*8,dimension(:,:) :: uErrv
      real*8,dimension(:,:) :: vErrv      
      real*8,dimension(:,:) :: wErrv
      real*8,dimension(:,:) :: pErrv      
      
      real*8,dimension(:,:) :: xvt
      real*8,dimension(:,:) :: yvt
      real*8,dimension(:,:) :: zvt
      
      integer,dimension(:,:):: No_vp

      END SUBROUTINE SaveParaviewVertex

!      ________________________________________________________
!     |                                                        |
!     |    SavetecCenter.F90                                   |
!     |________________________________________________________|

      SUBROUTINE SavetecCenter(alphafnp,ufnp,vfnp,wfnp,pfnp,     &
                               alphasnp,usnp,vsnp,wsnp,psnp,rhos,&
                               xct,yct,zct,No_cp,                &
                               ic1tec,ic2tec,ic3tec)

      real*8,dimension(:,:) :: alphafnp
      real*8,dimension(:,:) :: ufnp
      real*8,dimension(:,:) :: vfnp
      real*8,dimension(:,:) :: wfnp
      real*8,dimension(:,:) :: pfnp  
      real*8,dimension(:,:) :: alphasnp
      real*8,dimension(:,:) :: usnp
      real*8,dimension(:,:) :: vsnp
      real*8,dimension(:,:) :: wsnp
      real*8,dimension(:,:) :: psnp
      real*8,dimension(:,:) :: rhos
      real*8,dimension(:,:) :: xct
      real*8,dimension(:,:) :: yct
      real*8,dimension(:,:) :: zct
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: ic1tec
      integer,dimension(:)  :: ic2tec
      integer,dimension(:)  :: ic3tec

      END SUBROUTINE SavetecCenter

!      ________________________________________________________
!     |                                                        |
!     |    SavetecVetexCenter.F90                              |
!     |________________________________________________________|

      SUBROUTINE SavetecVC(alphafnp,ufnp,vfnp,wfnp,pfnp,&
                                     alphasnp,usnp,vsnp,wsnp,psnp,rhos,&
                                     xct,yct,zct,No_cp,                &
                                     alphafv,ufv,vfv,wfv,pfv,          &
                                     alphasv,usv,vsv,wsv,psv,rhosv,    &
                                     xvt,yvt,zvt,No_vp)

      real*8,dimension(:,:) :: alphafnp
      real*8,dimension(:,:) :: ufnp
      real*8,dimension(:,:) :: vfnp
      real*8,dimension(:,:) :: wfnp
      real*8,dimension(:,:) :: pfnp  
      real*8,dimension(:,:) :: alphasnp
      real*8,dimension(:,:) :: usnp
      real*8,dimension(:,:) :: vsnp
      real*8,dimension(:,:) :: wsnp
      real*8,dimension(:,:) :: psnp
      real*8,dimension(:,:) :: rhos
      real*8,dimension(:,:) :: xct
      real*8,dimension(:,:) :: yct
      real*8,dimension(:,:) :: zct
      integer,dimension(:,:):: No_cp

      real*8,dimension(:,:) :: alphafv
      real*8,dimension(:,:) :: ufv
      real*8,dimension(:,:) :: vfv
      real*8,dimension(:,:) :: wfv
      real*8,dimension(:,:) :: pfv  
      real*8,dimension(:,:) :: alphasv
      real*8,dimension(:,:) :: usv
      real*8,dimension(:,:) :: vsv
      real*8,dimension(:,:) :: wsv
      real*8,dimension(:,:) :: psv
      real*8,dimension(:,:) :: rhosv
      real*8,dimension(:,:) :: xvt
      real*8,dimension(:,:) :: yvt
      real*8,dimension(:,:) :: zvt
      integer,dimension(:,:):: No_vp

      END SUBROUTINE SavetecVC

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                             Time loop                               !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |    LoopTimeInitial.F90                                 |
!     |________________________________________________________|


      SUBROUTINE LoopTimeInitial(Hprnp,etanp,                   &
                                 alphafnp,ufnp,vfnp,wfnp,pfnp,  &
                                 alphasnp,usnp,vsnp,wsnp,psnp,  &                                 
                                 xct,yct,zct,                   & 
                                 xvt,yvt,zvt,                   &    
                                 Hprn,etan,                     &  
                                 alphafn,ufn,vfn,wfn,pfn,       &
                                 alphasn,usn,vsn,wsn,psn,       &
                                 xc,yc,sig,Hpr,h,               &
                                 xv,yv,sigv,Hprv,hv,            &
                                 No_cp,nbe)  

      real*8,dimension(:)  :: Hprn,etan  
      real*8,dimension(:,:):: alphafn,ufn,vfn,wfn,pfn
      real*8,dimension(:,:):: alphasn,usn,vsn,wsn,psn
      real*8,dimension(:,:):: xct,yct,zct  
      real*8,dimension(:,:):: xvt,yvt,zvt   
      real*8,dimension(:)  :: Hprnp,etanp    
      real*8,dimension(:,:):: alphafnp,ufnp,vfnp,wfnp,pfnp
      real*8,dimension(:,:):: alphasnp,usnp,vsnp,wsnp,psnp
      real*8,dimension(:)  :: xc,yc,sig,Hpr,h   
      real*8,dimension(:)  :: xv,yv,sigv,Hprv,hv  
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe    

      END SUBROUTINE LoopTimeInitial

!      ________________________________________________________
!     |                                                        |
!     |     LoopTimeFinal.F90                                  |
!     |________________________________________________________|

      SUBROUTINE LoopTimeUpdate(Hprn,etan,                    &
                                alphafn,ufn,vfn,wfn,pfn,      &
                                alphasn,usn,vsn,wsn,psn,      & 
                                xct,yct,zct,                  &  
                                xvt,yvt,zvt,                  &                                 
                                mask,                         & 
                                Hprnp,etanp,                  &
                                alphafnp,ufnp,vfnp,wfnp,pfnp, &
                                alphasnp,usnp,vsnp,wsnp,psnp, &
                                xc,yc,sig,Hpr,h,              &
                                xv,yv,sigv,Hprv,hv,           &
                                No_cp,nbe)    

      real*8,dimension(:)  :: Hprn,etan  
      real*8,dimension(:,:):: alphafn,ufn,vfn,wfn,pfn
      real*8,dimension(:,:):: alphasn,usn,vsn,wsn,psn
      real*8,dimension(:,:):: xct,yct,zct  
      real*8,dimension(:,:):: xvt,yvt,zvt   
      real*8,dimension(:,:):: mask    
      real*8,dimension(:)  :: Hprnp,etanp    
      real*8,dimension(:,:):: alphafnp,ufnp,vfnp,wfnp,pfnp
      real*8,dimension(:,:):: alphasnp,usnp,vsnp,wsnp,psnp
      real*8,dimension(:)  :: xc,yc,sig,Hpr,h   
      real*8,dimension(:)  :: xv,yv,sigv,Hprv,hv   
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe    

      END SUBROUTINE LoopTimeUpdate


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                             interpolation.F90                       !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |    interpolation.F90                                   |
!     |________________________________________________________|

      SUBROUTINE interpolation(phiv,phi,&
                               sig,dsig,sigv,dsigv,No_vp,nbev) 

      real*8 ,dimension(:,:) :: phiv
      real*8 ,dimension(:,:) :: phi
      real*8, dimension(:,:) :: sigv
      real*8, dimension(:,:) :: dsigv
      real*8, dimension(:,:) :: sig
      real*8, dimension(:,:) :: dsig
      integer,dimension(:,:) :: No_vp
      integer,dimension(:)   :: nbev

      END SUBROUTINE interpolation

!      ________________________________________________________
!     |                                                        |
!     |    interpolationETA.F90                                |
!     |________________________________________________________|

      SUBROUTINE interpolationEta(etav,eta,&
                                  No_vp,dlCV,dlVsum,areaCELL,areaVsum)    

      real*8 ,dimension(:)  :: etav
      real*8 ,dimension(:)  :: eta
      integer,dimension(:,:):: No_vp

      END SUBROUTINE interpolationEta

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                         Nonlinear convergence                       !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                   GrandientLSM.F90                     |
!     |________________________________________________________|


      SUBROUTINE grandientLSM(dfundx,dfundy,dfundsig,        &
                              fun,No_cp,nbe,                 &
                              sig)     

      real*8 ,dimension(:,:) :: dfundx
      real*8 ,dimension(:,:) :: dfundy
      real*8 ,dimension(:,:) :: dfundsig
      real*8 ,dimension(:,:) :: fun
      integer,dimension(:,:) :: No_cp
      integer,dimension(:)   :: nbe
      real*8, dimension(:)   :: sig


      END SUBROUTINE GrandientLSM


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                               Test case                             !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                        test.F90                        |
!     |________________________________________________________|

      SUBROUTINE test(fun,funv,                         &
                      xc,yc,sig,dsig,                   &
                      xv,yv,sigv,dsigv,                 &
                      No_cp,nbe,                        &
                      No_vp,nbev)  
  
      real*8,dimension(:,:) :: fun
      real*8,dimension(:,:) :: funv
      real*8,dimension(:)   :: xc,yc,sig,dsig
      real*8,dimension(:)   :: xv,yv,sigv,dsigv
      integer,dimension(:,:):: No_cp
      integer,dimension(:,:):: No_vp
      integer,dimension(:)  :: nbe
      integer,dimension(:)  :: nbev

      END SUBROUTINE test

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                               Advection                             !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                    resolvesystem.F90                   |
!     |________________________________________________________|


      SUBROUTINE solvesystem(xg,Am0,Am1,Am2,Am3,AmT,AmB,bm,&
                               No_cp,nbe)

      real*8, dimension(:,:) ::  xg
      real*8, dimension(:,:) ::  Am0
      real*8, dimension(:,:) ::  Am1
      real*8, dimension(:,:) ::  Am2
      real*8, dimension(:,:) ::  Am3
      real*8, dimension(:,:) ::  AmT
      real*8, dimension(:,:) ::  AmB
      real*8, dimension(:,:) ::  bm
      integer,dimension(:,:) ::  No_cp
      integer,dimension(:)   ::  nbe

      END SUBROUTINE solvesystem

!      ________________________________________________________
!     |                                                        |
!     |                     massflux.F90                       |
!     |________________________________________________________|

      SUBROUTINE massflux(C01,C02,C03,C0T,C0B,     &
                          uu,vv,ww,                &
                          No_cp,nbe,dsigv)

      real*8 ,dimension(:,:) :: C01
      real*8 ,dimension(:,:) :: C02
      real*8 ,dimension(:,:) :: C03
      real*8 ,dimension(:,:) :: C0T
      real*8 ,dimension(:,:) :: C0B
      real*8 ,dimension(:,:) :: uu
      real*8 ,dimension(:,:) :: vv
      real*8 ,dimension(:,:) :: ww
      integer,dimension(:,:) :: No_cp
      integer,dimension(:)   :: nbe 
      real*8, dimension(:)   :: dsigv

      END SUBROUTINE massflux

!      ________________________________________________________
!     |                                                        |
!     |                     advection.F90                      |
!     |________________________________________________________|

      SUBROUTINE advection(Am0,Am1,Am2,Am3,AmT,AmB,bm,         &
                           uu,vv,ww,                           &
                           phi,xc,yc,sig,No_cp,nbe,            &
                           sigv,dsigv)     

      real*8 ,dimension(:,:) :: Am0
      real*8 ,dimension(:,:) :: Am1             
      real*8 ,dimension(:,:) :: Am2            
      real*8 ,dimension(:,:) :: Am3            
      real*8 ,dimension(:,:) :: AmT
      real*8 ,dimension(:,:) :: AmB
      real*8 ,dimension(:,:) :: bm
      real*8 ,dimension(:,:) :: uu
      real*8 ,dimension(:,:) :: vv
      real*8 ,dimension(:,:) :: ww
      real*8 ,dimension(:,:) :: phi
      real*8, dimension(:)   :: xc
      real*8, dimension(:)   :: yc
      real*8, dimension(:)   :: sig
      integer,dimension(:,:) :: No_cp
      integer,dimension(:)   :: nbe 
      real*8, dimension(:)   :: sigv
      real*8, dimension(:)   :: dsigv

      END SUBROUTINE advection

!      ________________________________________________________
!     |                                                        |
!     |                  boundaryCondition.F90                 |
!     |________________________________________________________|

      SUBROUTINE boundaryCondition(phi,                        &
                                   phiT,phiB,                  &                                  
                                   xc,yc,sig,                  &
                                   xv,yv,sigv,                 &
                                   No_cp,nbe,                  &        
                                   No_vp,nbev)      

      real*8,dimension(:,:) :: phi
      real*8,dimension(:)   :: phiT
      real*8,dimension(:)   :: phiB
      real*8,dimension(:)   :: xc,yc,sig
      real*8,dimension(:)   :: xv,yv,sigv
      integer,dimension(:,:):: No_cp
      integer,dimension(:)  :: nbe
      integer,dimension(:,:):: No_vp
      integer,dimension(:)  :: nbev

      END SUBROUTINE boundaryCondition


!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!                               Utilyties                             !
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                              .F90                      |
!     |________________________________________________________|


      end interface

      END MODULE INTERFACES

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                                  END                                !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!



