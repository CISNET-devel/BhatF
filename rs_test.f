      implicit real*8(a-h,o-z)

      parameter(nmax=30)
      dimension a(nmax,nmax),w(nmax),z(nmax,nmax),fv1(nmax),fv2(nmax)

      n=3

      a(1,1)=1.; a(1,2)=0.; a(1,3)=-1.
      a(2,1)=0.; a(2,2)=4.; a(2,3)=-0.
      a(3,1)=-1.; a(3,2)=0.; a(3,3)=2.

      call rs(nmax,n,a,w,matz,z,fv1,fv2,ierr)

      print*,w(1:3)
      stop
      end
