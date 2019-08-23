!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!               Functions for User Defined Initial Field              !
!                              Dec 2015                               !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!      ________________________________________________________
!     |                                                        |
!     |                     Exact velocity                     |
!     |________________________________________________________|
!     fOR CHANNEL FLOW CASE
!     ________________________________________________________
!     Velocity component u

      function funu(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funu
      real*8 :: c,umax,uran,drand48
      integer :: iseed
      umax = 18.0d0
      uran = 8.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z .lt. 0.0d0) .or. (z .gt. 1.0d0)) then
       c = 0.0d0
      else
       c = umax*z*(2.0d0-z)+uran*2.0d0*(drand48-0.5d0)
      endif
      funu = c
      return
      end function funu
!     ________________________________________________________
!     Velocity component v

      function funv(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funv
      real*8 :: c, uran,drand48
      integer :: iseed
      uran = 10.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z .lt. 0.0d0) .or. (z .gt. 1.0d0)) then
       c = 0.0d0
      else
       c = uran*2.0d0*(drand48-0.5d0)
      endif
      funv = c
      return
      end function funv
!     ________________________________________________________
!     Velocity component w

      function funw(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funw
      real*8 :: c, uran,drand48
      integer :: iseed
      uran = 10.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z .lt. 0.0d0) .or. (z .gt. 1.0d0)) then
       c = 0.0d0
      else
       c = uran*2.0d0*(drand48-0.5d0)
      endif
      funw = c
      return
      end function funw
!     ________________________________________________________
!     Pressure

      function funp(x,y,z)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,funp
      real*8 :: c, uran,drand48
      integer :: iseed
      uran = 2.0d0
!     ----------------------------------------
!     Function
      call random_seed(iseed)
      call random_number(drand48)
      if ((z .lt. 0.0d0) .or. (z .gt. 1.0d0)) then
       c = 0.0d0
      else
       c = uran*2.0d0*(drand48-0.5d0)
      endif
      funp = 0.
      return
      end function funp
!      ________________________________________________________
!     |                                                        |
!     |                     Exact velocity                     |
!     |________________________________________________________|
      function inflow(x,y,z,flag)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,inflow
      integer :: flag
      if(flag .eq. 2) then
         inflow = 1.0d0      ! for u inflow
      else
         inflow = 0.0d0      ! the other two direction
      endif
      return
      end function inflow
!      ________________________________________________________
!     |                                                        |
!     |                     Exact velocity                     |
!     |________________________________________________________|
      function dirchlet(x,y,z,t,Re,flag)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t,Re,dirchlet
      real*8 :: funu2,funv2,funw2
      integer :: flag
!     -------------------
          if(flag .eq. 1) then
             dirchlet = funu2(x,y,z,t,Re)      ! for u
          elseif(flag .eq. 2) then
             dirchlet = funv2(x,y,z,t,Re)      ! for v
          elseif(flag .eq. 3) then
             dirchlet = funw2(x,y,z,t,Re)      ! for w
          else
            print*, '    ERROR! VELOCITY DIRECTION NOT SPECIFICED'
            stop
          endif
      return
      end function dirchlet
!      ________________________________________________________
!     |                                                        |
!     |                     Exact velocity                     |
!     |________________________________________________________|
!     fOR Green Taylor Vortex Case
!     ________________________________________________________
!     Velocity component u

      function funu2(x,y,z,t,Re)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t,funu2,Re
      real*8 :: c,drand48
      integer :: iseed
      real*8 :: pi = 3.14159265359d0
       call random_seed(iseed)
       call random_number(drand48)
       c = -1.0d0*dcos(pi*x)*dsin(pi*y)*dexp(-2.0d0*pi*pi*t/Re)
      funu2 =  c
      return
      end function funu2
!     ________________________________________________________
!     Velocity component v

      function funv2(x,y,z,t,Re)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t,funv2,Re
      real*8 :: c,drand48
      integer :: iseed
      real*8 :: pi = 3.14159265359d0
       call random_seed(iseed)
       call random_number(drand48)
      c = dsin(pi*x)*dcos(pi*y)*dexp(-2.0d0*pi*pi*t/Re)
      funv2 = c
      return
      end function funv2
!     ________________________________________________________
!     Velocity component w

      function funw2(x,y,z,t,Re)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t,funw2,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       c = 0.
      funw2 = 0.0d0
      return
      end function funw2
!     ________________________________________________________
!     Pressure

      function funp2(x,y,z,t,Re)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t,funp2,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       c = -0.25d0*(dcos(2.0*pi*x)+dcos(2.0*pi*y))*dexp(-4.0d0*pi*pi*t/Re)
      funp2 = c
      return
      end function funp2
!     ________________________________________________________
!     Source term

      function funs(x,y,z,t,Re)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t,funs,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       c = pi*pi*(dcos(2.0*pi*x)+dcos(2.0*pi*y))*dexp(-4.0d0*pi*pi*t/Re)
      funs = c
      return
      end function funs
!     ________________________________________________________
!     UFSTAR

      function funustar(x,y,z,t1,Re,t2)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t1,t2,t3,funustar,Re
      real*8 :: un,vn,uconv,udiff
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
      integer :: DoThis
       DoThis = 3
       t3 = t1-t2
       vn = dsin(pi*x)*dcos(pi*y)*dexp(-2.0d0*pi*pi*t3/Re)
       un = -1.0d0*dcos(pi*x)*dsin(pi*y)*dexp(-2.0d0*pi*pi*t3/Re)
       uconv =  pi*un*(dsin(pi*x)*dsin(pi*y)*dexp(-2.0d0*pi*pi*t3/Re)) &
                     -pi*vn*(dcos(pi*x)*dcos(pi*y)*dexp(-2.0d0*pi*pi*t3/Re))
       udiff = -2.0d0*pi*pi*un/Re
       c = un
       if (DoThis .eq. 1) c = c + udiff*t2
       if (DoThis .eq. 2) c = c - uconv*t2
       if (DoThis .eq. 3) c = c + udiff*t2 -uconv*t2
       funustar = c
      return
      end function funustar

!     ________________________________________________________
!     VFSTAR

      function funvstar(x,y,z,t1,Re,t2)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t1,t2,t3,funvstar,Re
      real*8 :: un,vn,vconv,vdiff
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
      integer :: DoThis
       DoThis = 3
       t3 = t1-t2
       vn = dsin(pi*x)*dcos(pi*y)*dexp(-2.0d0*pi*pi*t3/Re)
       un = -1.0d0*dcos(pi*x)*dsin(pi*y)*dexp(-2.0d0*pi*pi*t3/Re)
       vconv =  pi*un*(dcos(pi*x)*dcos(pi*y)*dexp(-2.0d0*pi*pi*t3/Re)) &
                     -pi*vn*(dsin(pi*x)*dsin(pi*y)*dexp(-2.0d0*pi*pi*t3/Re))
       vdiff = -2.0d0*pi*pi*vn/Re
       c = vn
       if (DoThis .eq. 1) c = c + vdiff*t2
       if (DoThis .eq. 2) c = c - vconv*t2
       if (DoThis .eq. 3) c = c + vdiff*t2 -vconv*t2
       funvstar = c
      return
      end function funvstar
!     ________________________________________________________
!     du/dx

      function fundudx(x,y,z,t1,Re,t2)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t1,t2,t3,fundudx,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       t3 = t1-t2
       c =  dcos(x)*dcos(y)*dexp(-2.0d0*t3/Re)
       fundudx = c
      return
      end function fundudx
!     ________________________________________________________
!     du/dy

      function fundudy(x,y,z,t1,Re,t2)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t1,t2,t3,fundudy,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       t3 = t1-t2
       c =  -dsin(x)*dsin(y)*dexp(-2.0d0*t3/Re)
       fundudy = c
      return
      end function fundudy
!     ________________________________________________________
!     dv/dx

      function fundvdx(x,y,z,t1,Re,t2)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t1,t2,t3,fundvdx,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       t3 = t1-t2
       c =  dsin(x)*dsin(y)*dexp(-2.0d0*t3/Re)
       fundvdx = c
      return
      end function fundvdx
!     ________________________________________________________
!     dv/dx

      function fundvdy(x,y,z,t1,Re,t2)
!     ----------------------------------------
!     Declaration of variables
      implicit none
      real*8 :: x,y,z,t1,t2,t3,fundvdy,Re
      real*8 :: c
      real*8 :: pi = 3.14159265359d0
       t3 = t1-t2
       c =  -1.0d0*dcos(x)*dcos(y)*dexp(-2.0d0*t3/Re)
       fundvdy = c
      return
      end function fundvdy

!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
!---------------------------------------------------------------------!
!                        End of auxiliar functions                    !
!---------------------------------------------------------------------!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww!
