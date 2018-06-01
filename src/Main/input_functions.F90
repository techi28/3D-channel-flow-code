!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 1 !        AUXILIAR: Function eta (free surface) & h (depth)        !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |            Free surface: eta (TestOnlyPoisson)         |
!     |________________________________________________________|

      function funEta(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"    
      real*8 :: x,y,funEta
      real*8  :: pi,delta 
      integer :: i
      pi = 4.0d0*datan(1.d0)
      
!     ----------------------------------------
!     KeyTestOnlyPoisson
#     if  defined(KeyTestOnlyPoisson) 
         !funEta = delta*dsin(x/8)*dsin(y/8) + 2.0d0*pi 
!        ------------------------
         !funEta = 0.1d0*dsin(2*pi*x)*dsin(4*pi*y)
!        ------------------------
!        Paper ParallelPoisson3D 
         !funEta = 0.25d0
!        ------------------------
!        Original Poisson2014
         !funEta = 0.1d0*dsin(2*pi*x)*dsin(4*pi*y)+0.25d0
!        ------------------------
!        ParCFD2015 Montreal, Canada
         !funEta = 0.05d0*dsin(2*pi*x)*dsin(4*pi*y)+0.25d0 
!        ------------------------
         !funEta = 2*pi + 0.1d0*(2*pi)*exp(-0.5d0*(x-pi)**2-0.5d0*(y-pi)**2)
!        ------------------------
!        Test Poisson only (Dic 2017)
         funEta = 0.0d0
#     endif      
!     ----------------------------------------      
!     KeyTaylorVortex: Function Paper IJNMF_2017  
#     if  defined(KeyTaylorVortex)   
         delta  = 0.25d0*(2*pi)
         funEta = 2*pi !+ delta*cos(x)*cos(y)
#     endif 

      return
      end function funEta
!      ________________________________________________________
!     |                                                        |
!     |                Depth: h (TestOnlyPoisson)              |
!     |________________________________________________________|

      function funh(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"          
      real*8 :: x,y,funh,pi,delta
      integer :: i
      pi = 4.0d0*datan(1.d0)
      
!     ----------------------------------------
!     KeyTestOnlyPoisson
#     if  defined(KeyTestOnlyPoisson)    
         !funh = -0.2d0*(sin(3.14*(x+y))-4.0d0)
!        ------------------------
!        Original Poisson2014
         !funh = -0.2d0*(sin(3.14*(x+y))-4.0d0)
!        ------------------------
!        ParCFD2015 Montreal, Canada
         !funh = -0.1d0*(sin(3.14*(x+y))+(x+y)-3.0d0)
!        ------------------------
         !funh  = 0.0d0 + 0.1d0*(2*pi)*exp(-0.5d0*(x-pi)**2-0.5d0*(y-pi)**2)
!        ------------------------
!        Test Poisson only (Dic 2017)
         funh = 2.0d0
#     endif
!     ----------------------------------------      
!     KeyTaylorVortex: Function Paper IJNMF_2017  
#     if  defined(KeyTaylorVortex)
         delta = 0.25d0*(2*pi)
         funh  = 0.0d0 !+ delta*cos(x)*cos(y)
#     endif

      return
      end function funh
      
!      ________________________________________________________
!     |                                                        |
!     |                Funtion Free Surface                    |
!     |________________________________________________________|

      function funFreeSurface(x,y,t)    

      implicit none
      real*8 :: x,y,t,funFreeSurface
      real*8 :: pi,L,D,A,k,TT,s
!     ----------------------------------------

      pi = 4.0d0*datan(1.d0)
      L  = 20.0d0  ! Basin length
      D  = 10.0d0  ! Basin depth
      A  = 0.1d0   ! Amplitude of the standing wave
      k  = pi/D
      TT = 3.59d0  ! Wave period approximation
      s  = 2*pi/TT ! s^2 = g*k*tanh(k*D)

      funFreeSurface  = A*dcos(k*pi)*dcos(s*t)

      return
      end function funFreeSurface     
 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 2 !              AUXILIAR: Solution N-S equations                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |                Exact velocity & pressure               |
!     |________________________________________________________|

!     ________________________________________________________
!     Velocity component u

      function funExamNSu(x,y,z,t)    
!     ------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSu
      real*8 :: c 
!     ------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      funExamNSu =-0.5d0*c*(dsqrt(3.0d0)*dcos(x)*dsin(y)*dsin(z)&
                                       + dsin(x)*dcos(y)*dcos(z))
      return
      end function funExamNSu
      
!     ________________________________________________________
!     Velocity component v

      function funExamNSv(x,y,z,t)    
!     ------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"      
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSv
      real*8 :: c 
!     ------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      funExamNSv = 0.5d0*c*(dsqrt(3.0d0)*dsin(x)*dcos(y)*dsin(z)&
                                       - dcos(x)*dsin(y)*dcos(z))
      return
      end function funExamNSv
      
!     ________________________________________________________
!     Velocity component w

      function funExamNSw(x,y,z,t)    
!     ------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"            
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSw
      real*8 :: c 
!     ------------------------
!     Function
      c = dexp(-3.0d0*t/Re)
      funExamNSw = c*dcos(x)*dcos(y)*dsin(z)

      return
      end function funExamNSw
      
!     ________________________________________________________
!     Pressure p

      function funExamNSp(x,y,z,t)    
!     ------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"       
#     include "common.mpf"
      real*8 :: x,y,z,t,funExamNSp
      real*8 :: uu,vv,ww,c,d 
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
      real*8 :: funtnp,funtn 
      
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)      
!     ------------------------
!     Function at time (n+1) = tt
      c = dexp(-3.0d0*t/Re)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funtnp = -0.5d0*(uu**2+vv**2+ww**2)
!     ------------------------
!     Function at time (n) = t - dt
      c = dexp(-3.0d0*(t-dt)/Re)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funtn = -0.5d0*(uu**2+vv**2+ww**2)
!     ------------------------
#     ifdef Key_ProjectionSecond
         funExamNSp = funtnp-funtn
#     else
         funExamNSp = funtnp
#     endif
!     ------------------------
#     ifdef KeyTestOnlyPoisson
         funExamNSp = sin(0.5d0*pi*x)*sin(0.5d0*pi*y)*sin(0.5d0*pi*z)
#     endif      

      return
      end function funExamNSp
      
!      ________________________________________________________
!     |                                                        |
!     |        Function Dp = pf^(n+1)-pf^(n) (pressure)        |
!     |________________________________________________________|

      function funExamNSDp(x,y,z,tt)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"     
#     include "common.mpf"
      real*8 :: x,y,z,tt,funExamNSDp
      real*8 :: uu,vv,ww,c,d,funtnp,funtn 
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
!     ----------------------------------------
!     Function at time (n+1) = tt
      c = dexp(-3.0d0*tt/Re)
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funtnp = -0.5d0*(uu**2+vv**2+ww**2)
!     ----------------------------------------
!     Function at time (n) = tt - dt
      c = dexp(-3.0d0*(tt-dt)/Re)
      uu =-0.5d0*c*(d*cosx*siny*sinz+sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz-cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      funtn = -0.5d0*(uu**2+vv**2+ww**2)
!     ----------------------------------------
!     Final function
      funExamNSDp = funtnp-funtn

      return
      end function funExamNSDp
      
!      ________________________________________________________
!     |                                                        |
!     |        Neumman boundary condition of NS equation "p"   |
!     |________________________________________________________|

      function NeumanndpdnNS(x,y,z,t,nnx,nny,nnz)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"      
#     include "common.mpf"
      real*8 :: NeumanndpdnNS
      real*8 :: x,y,z,t,nnx,nny,nnz
      real*8 :: uu,vv,ww,pp
      real*8 :: c,d 
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
      real*8 :: dudx,dudy,dudz
      real*8 :: dvdx,dvdy,dvdz
      real*8 :: dwdx,dwdy,dwdz
      real*8 :: dpdx,dpdy,dpdz      
      
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)            
!     ------------------------
!     Functions (u,v,w) and p
      c = dexp(-3.0d0*t/Re)
      uu =-0.5d0*c*(d*cosx*siny*sinz + sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz - cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
      pp = -0.5d0*(uu**2+vv**2+ww**2)
!     ------------------------
!     Derivatives of (u,v,w)
      dudx =-0.5d0*c*(-d*sinx*siny*sinz + cosx*cosy*cosz)
      dudy =-0.5d0*c*( d*cosx*cosy*sinz - sinx*siny*cosz)
      dudz =-0.5d0*c*( d*cosx*siny*cosz - sinx*cosy*sinz)
      dvdx = 0.5d0*c*( d*cosx*cosy*sinz + sinx*siny*cosz)
      dvdy = 0.5d0*c*(-d*sinx*siny*sinz - cosx*cosy*cosz)
      dvdz = 0.5d0*c*( d*sinx*cosy*cosz + cosx*siny*sinz)
      dwdx = c*(-sinx*cosy*sinz)
      dwdy = c*(-cosx*siny*sinz)
      dwdz = c*cosx*cosy*cosz
!     ------------------------
!     Derivatives of p
      dpdx = -uu*dudx-vv*dvdx-ww*dwdx
      dpdy = -uu*dudy-vv*dvdy-ww*dwdy
      dpdz = -uu*dudz-vv*dvdz-ww*dwdz
!     ------------------------
#     ifdef KeyTestOnlyPoisson
      dpdx = 0.5d0*cos(0.5d0*x)*sin(0.5d0*y)*sin(0.5d0*z)
      dpdy = 0.5d0*sin(0.5d0*x)*cos(0.5d0*y)*sin(0.5d0*z)
      dpdz = 0.5d0*sin(0.5d0*x)*sin(0.5d0*y)*cos(0.5d0*z)
#     endif                     
!     ------------------------
#     ifdef KeyStaticCylinder 
      dpdx = 0.0d0
      dpdy = 0.0d0
      dpdz = 0.0d0                
#     endif
!     ------------------------
#     ifdef KeyStaticChannel 
      dpdx = 0.0d0
      dpdy = 0.0d0
      dpdz = 0.0d0                
#     endif 
!     ------------------------
!     Function dp/dn 
      NeumanndpdnNS = nnx*dpdx + nny*dpdy + nnz*dpdz

      return
      end function NeumanndpdnNS
            
!      ________________________________________________________
!     |                                                        |
!     |     Exact right-hand side of the Poisson problem       |
!     |________________________________________________________|

      function funExamNSrhsp(x,y,z,tt)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"      
#     include "common.mpf"
      real*8 :: x,y,z,tt,funExamNSrhsp
      real*8 :: c,d 
      real*8 :: uu,vv,ww      
      real*8 :: sinx,siny,sinz
      real*8 :: cosx,cosy,cosz
      real*8 :: dudx,dudy,dudz
      real*8 :: dvdx,dvdy,dvdz
      real*8 :: dwdx,dwdy,dwdz
!     ----------------------------------------
!     Function
      c = dexp(-3.0d0*tt/Re)
      d = dsqrt(3.0d0)
      sinx = dsin(x)
      siny = dsin(y)
      sinz = dsin(z)
      cosx = dcos(x)
      cosy = dcos(y)
      cosz = dcos(z)
!     ------------------------
      uu =-0.5d0*c*(d*cosx*siny*sinz + sinx*cosy*cosz)
      vv = 0.5d0*c*(d*sinx*cosy*sinz - cosx*siny*cosz)
      ww = c*cosx*cosy*sinz
!     ------------------------
      dudx =-0.5d0*c*(-d*sinx*siny*sinz + cosx*cosy*cosz)
      dudy =-0.5d0*c*( d*cosx*cosy*sinz - sinx*siny*cosz)
      dudz =-0.5d0*c*( d*cosx*siny*cosz - sinx*cosy*sinz)
      dvdx = 0.5d0*c*( d*cosx*cosy*sinz + sinx*siny*cosz)
      dvdy = 0.5d0*c*(-d*sinx*siny*sinz - cosx*cosy*cosz)
      dvdz = 0.5d0*c*( d*sinx*cosy*cosz + cosx*siny*sinz)
      dwdx = -c*sinx*cosy*sinz
      dwdy = -c*cosx*siny*sinz
      dwdz =  c*cosx*cosy*cosz
!     ------------------------
      funExamNSrhsp = 3.0d0*(uu**2+vv**2+ww**2)   &
                       -( dudx**2+dudy**2+dudz**2 &
                         +dvdx**2+dvdy**2+dvdz**2 &
                         +dwdx**2+dwdy**2+dwdz**2 ) 

!     ------------------------
#     ifdef KeyTestOnlyPoisson
      funExamNSrhsp = -3.0d0*(pi**2)*0.25d0* &
                       sin(0.5d0*pi*x)*sin(0.5d0*pi*y)*dsin(0.5d0*pi*z)
#     endif 

      return
      end function funExamNSrhsp 
 
!      ________________________________________________________
!     |                                                        |
!     |         Inflow velocity of the cylinder problem        |
!     |________________________________________________________|
 
      function funInflow(x,y,z,t)    
!     ------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
#     include "common.mpf"
      real*8 :: x,y,z,t,funInflow
      real*8 :: c
      integer:: ChooseInflowProfile
      
!     ------------------------
!     Choose profile
#     if   defined(KeyInflowConstant)
         ChooseInflowProfile = 1
#     elif defined(KeyInflowPoiseuille)
         ChooseInflowProfile = 2
#     else
         ChooseInflowProfile = 1
#     endif

!     ------------------------
!     Constant profile
      IF (ChooseInflowProfile.eq.1) THEN
         funInflow = Uinflow
!     ------------------------
!     Poiseuille profile
      ELSEIF (ChooseInflowProfile.eq.2) THEN   
         funInflow = Uinflow*(1.0d0-(z/H0)**2)
!     ------------------------
!     Other profile
      ELSEIF (ChooseInflowProfile.eq.3) THEN 
         if (z.ge.0.0d0) then
            funInflow = 0.0d0
         else
            funInflow = Uinflow*(abs(z)**(1.0d0/3.0d0))
         endif
      ENDIF   
      
      return
      end function funInflow
 
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!        AUXILIAR: Functions of the Standing Waves problem            !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

!      ________________________________________________________
!     |                                                        |
!     |           Standing Wave:  Exact free surface           |
!     |________________________________________________________|

      function FS_funeta(x,y,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: x,y,t,FS_funeta,nn
      real*8 :: L,W,h,A,k,omega
      real*8 :: kx,ky,TT
      real*8, parameter :: pi = 3.14159265359d0
      real*8, parameter :: g  = 9.80665d0 
!     ----------------------------------------
      L  = 10.0d0 ! m
      W  = 10.0d0 ! m
      h  = 10.0d0 ! m
      A  = 0.1d0  ! m     
!     ----------------------------------------
!     Function 2D
#     ifdef KeyStanWaveExample2D
         nn = 2.0d0
         k  = 2*pi/(nn*L)   
         omega = dsqrt(g*k*tanh(k*h)) 
         FS_funeta = A*dcos(k*x)*dcos(omega*t)
#     endif
!     ----------------------------------------
!     Function 3D
#     ifdef KeyStanWaveExample3D
         kx = pi/L   
         ky = pi/W
         TT  = 3.01d0
         omega = 2*pi/TT
         FS_funeta = A*dcos(kx*x)*dcos(ky*y)*dcos(omega*t)
#     endif      

      return
      end function FS_funeta

!      ________________________________________________________
!     |                                                        |
!     |           Standing Wave:  Exact velocity               |
!     |________________________________________________________|

!     ________________________________________________________
!     Velocity component u

      function FS_funu(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: x,y,z,t,FS_funu
      real*8 :: L,W,h,A,k,omega,c,nn
      real*8 :: kx,ky,TT
      real*8, parameter :: pi = 3.14159265359d0
      real*8, parameter :: g  = 9.80665d0
!     ----------------------------------------      
      L  = 10.0d0     
      W  = 10.0d0    
      h  = 10.0d0      
      A  = 0.1d0
!     ----------------------------------------
!     Function 2D
#     ifdef KeyStanWaveExample2D
         nn = 2.0d0       
         k  = 2.0d0*pi/(nn*L)   
         omega = dsqrt(g*k*tanh(k*h)) 
         c = A*omega*dcosh(k*(h+z))/dsinh(k*h)
         FS_funu = c*dsin(k*x)*dsin(omega*t)
#     endif
!     ----------------------------------------
!     Function 3D
#     ifdef KeyStanWaveExample3D
         kx = pi/L   
         ky = pi/W
         k  = 0.44d0 ! = sqrt(kx**2 +ky**2)        
         TT = 3.01d0   
         omega = 2*pi/TT
         c = (A*g*kx/omega)*dcosh(k*(h+z))/dcosh(k*h)
         FS_funu = c*dsin(kx*x)*dcos(ky*y)*dsin(omega*t)
#     endif 
      
      return
      end function FS_funu

!     ________________________________________________________
!     Velocity component v

      function FS_funv(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: x,y,z,t,FS_funv
      real*8 :: L,W,h,A,k,omega,c,nn
      real*8 :: kx,ky,TT
      real*8, parameter :: pi = 3.14159265359d0
      real*8, parameter :: g  = 9.80665d0
!     ----------------------------------------
      L  = 10.0d0     
      W  = 10.0d0    
      h  = 10.0d0      
      A  = 0.1d0 
!     ----------------------------------------
!     Function 2D
#     ifdef KeyStanWaveExample2D
         FS_funv = 0.0d0
#     endif
!     ----------------------------------------
!     Function 3D
#     ifdef KeyStanWaveExample3D
         kx = pi/L   
         ky = pi/W
         k  = 0.44d0 ! = sqrt(kx**2 +ky**2)      
         TT = 3.01d0 
         omega = 2*pi/TT 
         c = (A*g*ky/omega)*dcosh(k*(h+z))/dcosh(k*h)
         FS_funv = c*dcos(kx*x)*dsin(ky*y)*dsin(omega*t)        
#     endif

      return
      end function FS_funv

!     ________________________________________________________
!     Velocity component w

      function FS_funw(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: x,y,z,t,FS_funw
      real*8 :: L,W,h,A,k,omega,c,nn
      real*8 :: kx,ky,TT
      real*8, parameter :: pi = 3.14159265359d0
      real*8, parameter :: g  = 9.80665d0 
!     ----------------------------------------      
      L  = 10.0d0     
      W  = 10.0d0    
      h  = 10.0d0      
      A  = 0.1d0
!     ----------------------------------------
!     Function 2D
#     ifdef KeyStanWaveExample2D
         nn = 2.0d0       
         k  = 2*pi/(nn*L)
         omega = dsqrt(g*k*tanh(k*h))
         c = - A*omega*dsinh(k*(h+z))/dsinh(k*h)
         FS_funw = c*dcos(k*x)*dsin(omega*t)
#     endif
!     ----------------------------------------
!     Function 3D
#     ifdef KeyStanWaveExample3D
         kx = pi/L   
         ky = pi/W
         k  = 0.44d0      
         TT = 3.01d0 
         omega = 2*pi/TT
         c = - (A*g*k/omega)*dsinh(k*(h+z))/dcosh(k*h)
         FS_funw = c*dcos(kx*x)*dcos(ky*y)*dsin(omega*t)
#     endif

      return
      end function FS_funw

!      ________________________________________________________
!     |                                                        |
!     |         Standing Wave:  Exact dynamic pressure         |
!     |________________________________________________________|

      function FS_funp(x,y,z,t)
#     include "cppdefs.h"
!     ----------------------------------------
!     Declaration of variables         
      implicit none
      real*8 :: x,y,z,t,FS_funp
      real*8 :: L,W,h,A,k,omega,c,eta,rho,nn
      real*8 :: kx,ky,TT
      real*8, parameter :: pi = 3.14159265359d0
      real*8, parameter :: g  = 9.80665d0
!     ----------------------------------------
      L  = 10.0d0
      W  = 10.0d0
      h  = 10.0d0
      A  = 0.1d0
!     ----------------------------------------
!     Function 2D
#     ifdef KeyStanWaveExample2D
         rho= 1.0d0
         nn = 2.0d0
         k  = 2*pi/(nn*L)
         omega = dsqrt(g*k*tanh(k*h))
         c   = A*rho*g*dcosh(k*(h+z))/dcosh(k*h)
         eta = A*dcos(k*x)*dcos(omega*t)
         FS_funp = -rho*g*eta + c*dcos(k*x)*dcos(omega*t)
#     endif
!     ----------------------------------------
!     Function 3D
#     ifdef KeyStanWaveExample3D
         rho= 1.0d0
         kx = pi/L   
         ky = pi/W
         k  = 0.44d0      
         TT = 3.01d0 
         omega = 2*pi/TT
         c = - (A*g*k/omega)*dsinh(k*(h+z))/dcosh(k*h)
         eta = A*dcos(kx*x)*dcos(ky*y)*dcos(omega*t)        
         FS_funp = -rho*g*eta + c*dcos(kx*x)*dcos(ky*y)*dcos(omega*t)
#     endif

      return
      end function FS_funp

!      ________________________________________________________
!     |                                                        |
!     |           Standing Wave:  Initial Water Depth          |
!     |________________________________________________________|

      function FS_funH(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: x,y,FS_funH
      real*8 :: xx,yy,fun,rr,theta
      real*8, parameter :: pi = 3.14159265359d0
      real*8, parameter :: g  = 9.80665d0
!     ----------------------------------------
      FS_funH = 10.d0
!     ----------------------------------------
      return
      end function FS_funH
      
!      ________________________________________________________
!     |                                                        |
!     |                      Incident wave                     |
!     |________________________________________________________|

      function FS_waveBC(t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: t,FS_waveBC
      real*8 :: WL0,WL1,tt0,tt1
      integer :: II
      real*8, parameter :: height = 4.64
      real*8, parameter :: period = 1.0
!     ----------------------------------------
!     Linear function with period and heigh d
      WL0 = height
      WL1 = 0.0d0
      II = int(t/period) + 1
      tt0 = (II-1)*period
      tt1 = II*period
      FS_waveBC = (WL1-WL0)/(tt1-tt0)*(t-tt0) + wl0
!     ----------------------------------------
      return
      end function FS_waveBC
      
!      ________________________________________________________
!     |                                                        |
!     |                      Solitary wave                     |
!     |________________________________________________________|

      function FS_SolitaryWave(x,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "cppdefs.h"
      real*8 :: x,t,FS_SolitaryWave
      real*8 :: H,d,x0,c,eta
      real*8, parameter :: g = 9.80665
!     ----------------------------------------
!     Linear function with period and heigh d
      H  = 2.0d0  ! wave amplitud
      d  = 10.0d0 ! water depth
      x0 = 80.0d0 
      c  = sqrt(g*(d+H));
      eta = H/cosh(sqrt((0.75d0)*H/(d**3))*(x-x0-c*t))**2
      FS_SolitaryWave = eta
!     ----------------------------------------
      return
      end function FS_SolitaryWave
      
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 3 !           AUXILIAR: Functions of channel flow                   !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!

#     ifdef KeyStaticChannel
!     ________________________________________________________
!     Velocity component u
      function funu_ChannelFlow(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funu_ChannelFlow
      real*8 :: c,umax,uran,drand48
      integer :: iseed
      umax = 18.0d0
      uran = 10.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z.lt.0.0d0).or.(z.gt.1.0d0)) then
       c = 0.0d0
      else
!     ----------
#     ifdef KeyBCtopNoSlip
         c = umax*z*(1.0d0-z) + uran*2.0d0*(drand48-0.5d0)
#     else
         c = umax*z*(2.0d0-z) + uran*2.0d0*(drand48-0.5d0)
#     endif
!     ----------
      endif
      funu_ChannelFlow = c
      return
      end function funu_ChannelFlow
!     ________________________________________________________
!     Velocity component v

      function funv_ChannelFlow(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funv_ChannelFlow
      real*8 :: c, uran,drand48
      integer :: iseed
      uran = 10.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z.le.0.0d0).or.(z.ge.1.0d0)) then
       c = 0.0d0
      else
       c = uran*2.0d0*(drand48-0.5d0)
      endif
      funv_ChannelFlow = c
      return
      end function funv_ChannelFlow
!     ________________________________________________________
!     Velocity component w

      function funw_ChannelFlow(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funw_ChannelFlow
      real*8 :: c, uran,drand48
      integer :: iseed
      uran = 10.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z.le.0.0d0).or.(z.ge. 1.0d0)) then
       c = 0.0d0
      else
       c = uran*2.0d0*(drand48-0.5d0)
      endif
      funw_ChannelFlow = c
      return
      end function funw_ChannelFlow
!     ________________________________________________________
!     Pressure

      function funp_ChannelFlow(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funp_ChannelFlow
      real*8 :: c, uran,drand48
      integer :: iseed
      uran = 8.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z.le.0.0d0).or.(z.ge.1.0d0)) then
       c = 0.0d0
      else
       c = uran*2.0d0*(drand48-0.5d0)
      endif
      funp_ChannelFlow = 0.0d0 !<--- NOT RANDOM
      return
      end function funp_ChannelFlow

#     endif

!      ________________________________________________________
!     |                                                        |
!     |                       TestPBC3D                        |
!     |________________________________________________________|

      function TestPBC3D(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,TestPBC3D
      real*8 :: x0,y0,z0,uu,vv,ww,rx,ry,rz,ra,c1,c2,L0,Gam
!     __________________________________________
!     Example: Pure advection (need to fix) 3D
         c2 = 0.5d0*(1.0d0/dsqrt(2.0d0))
         x0 = -c2
         y0 = -c2
         z0 = -c2
         rx =  x-x0
         ry =  y-y0
         rz =  z-z0
         ra = dsqrt(rx**2+ry**2+rz**2)
         if (ra.le.0.25d0) then
            TestPBC3D = dcos(2.0d0*pi*ra)**2
         else
            TestPBC3D = 0.0d0
         endif
      end function TestPBC3D

      
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 3 !           AUXILIAR: Functions of 2D examples                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                         Exam2D                         |
!     |________________________________________________________|

      function funSolExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,funSolExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funSolExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funSolExam2D = (x**2-1)*(y**2-1)
      elseif (FunctionExample.eq.2) then
          funSolExam2D = dsin(pi*x)*dsin(pi*y)
      endif 

      return
      end function funSolExam2D
!      ________________________________________________________
!     |                                                        |
!     |                     Gradient Exam2D                    |
!     |________________________________________________________|
!     ________________________________________________________
!     dfdx 

      function dfdxExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,dfdxExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdxExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdxExam2D = (2.0d0*x)*(y**2-1)
      elseif (FunctionExample.eq.2) then
          dfdxExam2D = pi*dcos(pi*x)*dsin(pi*y) 
      endif 

      return
      end function dfdxExam2D
!     ________________________________________________________
!     dfdy

      function dfdyExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,dfdyExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.1) then
          dfdyExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdyExam2D = (x**2-1)*(2.0d0*y)
      elseif (FunctionExample.eq.2) then
          dfdyExam2D = pi*dsin(pi*x)*dcos(pi*y)  
      endif 

      return
      end function dfdyExam2D
!      ________________________________________________________
!     |                                                        |
!     |         Neumman boundary condition of Exam2D           |
!     |________________________________________________________|

      function Neumanndfdn2D(x,y,nnx,nny)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,nnx,nny,Neumanndfdn2D
      real*8 :: dfdx,dfdy
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)
          dfdy = (x**2-1)*(2.0d0*y)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y) 
          dfdy = pi*dsin(pi*x)*dcos(pi*y)  
      endif 
      Neumanndfdn2D = nnx*dfdx+nny*dfdy

      return
      end function Neumanndfdn2D
!      ________________________________________________________
!     |                                                        |
!     |               Advection term of Exam2D                 |
!     |________________________________________________________|

      function funAdvExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,funAdvExam2D
      real*8 :: dfdx,dfdy
!     ----------------------------------------
!     Function
      if  (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)
          dfdy = (x**2-1)*(2.0d0*y)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y) 
          dfdy = pi*dsin(pi*x)*dcos(pi*y)  
      endif 
      funAdvExam2D = dfdx+dfdy

      return
      end function funAdvExam2D
!      ________________________________________________________
!     |                                                        |
!     |               Diffusion term of Exam2D                 |
!     |________________________________________________________|

      function funDiffExam2D(x,y)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,funDiffExam2D
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funDiffExam2D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funDiffExam2D = 2.0d0*((y*y-1)+(x*x-1))
      elseif (FunctionExample.eq.2) then
          funDiffExam2D = -2.0d0*pi*pi*dsin(pi*x)*dsin(pi*y) 
      endif 

      return
      end function funDiffExam2D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 4 !              AUXILIAR: Functions of 3D examples                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                         Exam3D                         |
!     |________________________________________________________|

      function funSolExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,funSolExam3D,c     
      integer, parameter :: n = 1
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funSolExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funSolExam3D = -(x**2-1)*(y**2-1)*(z**2-1)
      elseif (FunctionExample.eq.2) then
          !funSolExam3D = dsin(pi*x)*dsin(pi*y)*dsin(pi*z)
          !funSolExam3D = dsin(pi*(x+1)/2)*dsin(pi*(y+1)/2)*dsin(pi*(z+1)/2)
          funSolExam3D = dsin(n*pi*x)*dsin(n*pi*y)*dsin(n*pi*z)/(n*n)
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          funSolExam3D = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function funSolExam3D

!      ________________________________________________________
!     |                                                        |
!     |                Diffusion term of Exam3D                |
!     |           Poisson right-hand side examples             |
!     |________________________________________________________|

      function funDiffExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,funDiffExam3D
      integer, parameter :: n = 1
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          funDiffExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          funDiffExam3D = -2.0d0*( (y**2-1)*(z**2-1) &
                                  +(x**2-1)*(z**2-1) &
                                  +(x**2-1)*(y**2-1) ) 
      elseif (FunctionExample.eq.2) then
          !funDiffExam3D = -3.0d0*pi*pi*dsin(pi*x)*dsin(pi*y)*dsin(pi*z)
          !funDiffExam3D = -3.0d0*pi*pi/4*dsin(pi*(x+1)/2)*dsin(pi*(y+1)/2)*dsin(pi*(z+1)/2)
          funDiffExam3D = -3.0d0*pi*pi*dsin(n*pi*x)*dsin(n*pi*y)*dsin(n*pi*z)
      elseif (FunctionExample.eq.3) then
          funDiffExam3D = 0.0d0
      endif 

      return
      end function funDiffExam3D

!      ________________________________________________________
!     |                                                        |
!     |                   Gradient of Exam3D                   |
!     |________________________________________________________|
!     ________________________________________________________
!     dfdx                 
      function dfdxExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,dfdxExam3D,c
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdxExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdxExam3D = (2.0d0*x)*(y**2-1)*(z**2-1)
      elseif (FunctionExample.eq.2) then
          dfdxExam3D = pi*dcos(pi*x)*dsin(pi*y)*dsin(pi*z) 
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdxExam3D = c*sqrt(2.0d0)*pi*(2.0d0*cosh(pi*sqrt(2.0d0)*x)-cosh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function dfdxExam3D
!     ________________________________________________________
!     dfdy

      function dfdyExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,dfdyExam3D,c
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdyExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdyExam3D = (x**2-1)*(2.0d0*y)*(z**2-1)
      elseif (FunctionExample.eq.2) then
          dfdyExam3D = pi*dsin(pi*x)*dcos(pi*y)*dsin(pi*z)  
      elseif (FunctionExample.eq.3) then
          c = pi*dcos(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdyExam3D = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function dfdyExam3D
!     ________________________________________________________
!     dfdz

      function dfdzExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,dfdzExam3D,c
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdzExam3D = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdzExam3D = (x**2-1)*(y**2-1)*(2.0d0*z)
      elseif (FunctionExample.eq.2) then
          dfdzExam3D = pi*dsin(pi*x)*dsin(pi*y)*dcos(pi*z)  
      elseif (FunctionExample.eq.3) then
          c = pi*dsin(pi*y)*pi*dcos(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdzExam3D = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 

      return
      end function dfdzExam3D
!      ________________________________________________________
!     |                                                        |
!     |          Neumman boundary condition of Exam3D          |
!     |________________________________________________________|

      function Neumanndfdn3D(x,y,z,nnx,nny,nnz)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,nnx,nny,nnz,Neumanndfdn3D,c
      real*8 :: dfdx,dfdy,dfdz
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
          dfdz = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)*(z**2-1)
          dfdy = (x**2-1)*(2.0d0*y)*(z**2-1)
          dfdz = (x**2-1)*(y**2-1)*(2.0d0*z)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y)*dsin(pi*z)
          dfdy = pi*dsin(pi*x)*dcos(pi*y)*dsin(pi*z)
          dfdz = pi*dsin(pi*x)*dsin(pi*y)*dcos(pi*z)
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdx = c*sqrt(2.0d0)*pi*(2.0d0*cosh(pi*sqrt(2.0d0)*x)-cosh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dcos(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdy = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dsin(pi*y)*pi*dcos(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdz = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      elseif (FunctionExample.eq.5) then
          dfdx = 0.0d0
          dfdy = 0.0d0
          dfdz = 0.0d0
      endif 
      Neumanndfdn3D = nnx*dfdx+nny*dfdy+nnz*dfdz

      return
      end function Neumanndfdn3D
!      ________________________________________________________
!     |                                                        |
!     |                Advection term of Exam3D                |
!     |________________________________________________________|

      function funAdvExam3D(x,y,z)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,funAdvExam3D,c
      real*8 :: dfdx,dfdy,dfdz
!     ----------------------------------------
!     Function
      if (FunctionExample.eq.0) then
          dfdx = 0.0d0
          dfdy = 0.0d0
          dfdz = 0.0d0
      elseif (FunctionExample.eq.1) then
          dfdx = (2.0d0*x)*(y**2-1)*(z**2-1)
          dfdy = (x**2-1)*(2.0d0*y)*(z**2-1)
          dfdz = (x**2-1)*(y**2-1)*(2.0d0*z)
      elseif (FunctionExample.eq.2) then
          dfdx = pi*dcos(pi*x)*dsin(pi*y)*dsin(pi*z)
          dfdy = pi*dsin(pi*x)*dcos(pi*y)*dsin(pi*z)
          dfdz = pi*dsin(pi*x)*dsin(pi*y)*dcos(pi*z)
      elseif (FunctionExample.eq.3) then
          c = dsin(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdx = c*sqrt(2.0d0)*pi*(2.0d0*cosh(pi*sqrt(2.0d0)*x)-cosh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dcos(pi*y)*dsin(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdy = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
          c = pi*dsin(pi*y)*pi*dcos(pi*z)/sinh(pi*sqrt(2.0d0))
          dfdz = c*(2.0d0*sinh(pi*sqrt(2.0d0)*x)+sinh(pi*sqrt(2.0d0)*(1-x)))
      endif 
      funAdvExam3D = dfdx+dfdy+dfdz

      return
      end function funAdvExam3D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
! 5 !               AUXILIAR: Time problems solutions                 !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                           2D                           |
!     |________________________________________________________|

      function TimeExample2D(x,y,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,t,TimeExample2D
      real*8 :: x0,y0,uu,vv,rx,ry,ra,c1,c2,L0,Gam,cx,cy
!     __________________________________________
!     Example: RK2
      if (RUNtestRK2approx.eq.1) then
          TimeExample2D = dsin(t)
      endif
!     __________________________________________
!     Example: Pure advection 2D
      if (RUNtestAdvEqn.eq.1) then
         x0 = -0.5d0*dcos(2.0d0*pi*t)
         y0 = -0.5d0*dsin(2.0d0*pi*t)
         rx = x-x0
         ry = y-y0
         ra = dsqrt(rx**2+ry**2) 
         if (ra.le.0.25d0) then 
            TimeExample2D = dcos(2.0d0*pi*ra)**2
         else
            TimeExample2D = 0.0d0
         endif
      endif
!     __________________________________________
!     Example: Pure diffusion 2D
      if (RUNtestDiffEqn.eq.1) then
         Gam =  5.0d-2 
         L0  =  7.5d-2
         c2  = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1  =  1.0d0/(4.0d0*pi*(Gam*t+L0**2))         
         rx = x-0.0d0
         ry = y-0.0d0
         TimeExample2D = c1*dexp(c2*(rx**2+ry**2)) 
      endif
!     __________________________________________
!     Example: advection-diffusion 2D
      if (RUNtestAdvDiffEqn.eq.1) then
         Gam =  3.0d-2 
         L0  =  7.5d-2
         uu  =  5.0d0
         vv  =  5.0d0
         x0  = -0.5d0
         y0  = -0.5d0
         c2  = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1  =  1.0d0/(4.0d0*pi*(Gam*t+L0**2))         
         rx = x-x0-uu*t
         ry = y-y0-vv*t
         TimeExample2D = c1*dexp(c2*(rx**2+ry**2))            
      endif

      return
      end function TimeExample2D

!      ________________________________________________________
!     |                                                        |
!     |                           3D                           |
!     |________________________________________________________|

      function TimeExample3D(x,y,z,t)    
!     ----------------------------------------
!     Declaration of variables         
      implicit none
#     include "common.mpf"
      real*8 :: x,y,z,t,TimeExample3D
      real*8 :: x0,y0,z0,uu,vv,ww,rx,ry,rz,ra,c1,c2,L0,Gam
!     __________________________________________
!     Example: RK2
      if (RUNtestRK2approx.eq.1) then         
          TimeExample3D = dsin(t)
      endif
!     __________________________________________
!     Example: Pure advection (need to fix) 3D
      if (RUNtestAdvEqn.eq.1) then
         c2 = 0.5d0*(1.0d0/dsqrt(2.0d0))
         x0 = -c2*dcos(2.0d0*pi*t)
         y0 =  c2*dcos(2.0d0*pi*t)
         z0 = -c2*dcos(2.0d0*pi*t)
         rx =  x-x0
         ry =  y-y0
         rz =  z-z0      
         ra = dsqrt(rx**2+ry**2+rz**2)
         if (ra.le.0.25d0) then 
            TimeExample3D = dcos(2.0d0*pi*ra)**2
         else
            TimeExample3D = 0.0d0
         endif
      endif
!     __________________________________________
!     Example: Pure diffusion 3D
      if (RUNtestDiffEqn.eq.1) then
         Gam=  5.0d-2 
         L0 =  7.5d-2
         c2 = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1 =  1.0d0/dsqrt((4.0d0*pi*(Gam*t+L0**2))**3)
         rx =  x-0.0d0
         ry =  y-0.0d0
         rz =  z-0.0d0                     
         TimeExample3D = c1*dexp(c2*(rx**2+ry**2+rz**2))  
      endif
!     __________________________________________
!     Example: advection-diffusion 3D
      if (RUNtestAdvDiffEqn.eq.4) then
         Gam=  3.0d-2 
         L0 =  7.5d-2
         uu =  5.0d0
         vv =  5.0d0
         ww =  5.0d0
         x0 = -0.5d0
         y0 = -0.5d0
         z0 = -0.5d0
         c2 = -1.0d0/(4.0d0*(Gam*t+L0**2))
         c1 =  1.0d0/dsqrt((4.0d0*pi*(Gam*t+L0**2))**3)
         rx =  x-x0-uu*t
         ry =  y-y0-vv*t
         rz =  z-z0-ww*t                     
         TimeExample3D = c1*dexp(c2*(rx**2+ry**2+rz**2))  
      endif

      return
      end function TimeExample3D

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        End of auxiliar functions                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
