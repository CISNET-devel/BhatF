c     --- back transformations 

      subroutine btrafo(npar,u,p)
      implicit real*8(a-h,o-z)
      parameter(nmax=100)
      dimension u(nmax),p(nmax)
      COMMON /BOUNDS/SU(nmax),SE(NMAX),S1(nmax),S2(nmax)
      
c      this assumes logit transformations of parameters
c      that are bounded from below by a=0 and from above by b=constant

      do i=1,npar
         
         rho=dexp(u(i))
         p(i)=(rho*s2(i)+s1(i))/(1.d0+rho)

      enddo

      return
      end


         
