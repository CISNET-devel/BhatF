      double precision function randgs (iseed, xmean, sd)
      implicit real*8(a-h,o-z)

c generate a normally distributed random number, i.e., generate random
c numbers with a gaussian distribution.  these random numbers are not
c exceptionally good -- especially in the tails of the distribution,
c but this implementation is simple and suitable for most applications.
c see r. w. hamming, numerical methods for scientists and engineers,
c mcgraw-hill, 1962, pages 34 and 389.

c             input arguments --
c xmean  the mean of the gaussian distribution.
c sd     the standard deviation of the gaussian function
c          exp (-1/2 * (x-xmean)**2 / sd**2)

      randgs = -6.d0
      do 10 i=1,12
        randgs = randgs + ran2(iseed)
 10   continue

      randgs = xmean + sd*randgs

      return
      end


