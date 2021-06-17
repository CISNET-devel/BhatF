      implicit real*8(a-h,o-z)
      external alngam
      
      print*,'input:'
      read*,x
      print*,alngam(x)
      stop
      end
