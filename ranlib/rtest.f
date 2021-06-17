      implicit real*8(a-h,o-z)
      integer iseed1,iseed2

      PRINT*,'iseed1, iseed2:'
      READ(5,*) iseed1,iseed2
      print*,'setall'
      CALL SETALL(ISEED1,ISEED2)

      sa=10.
      sr=sa*5.

      do i=1,1000
         print*,gengam(real(sa),real(sr))
      enddo
      stop
      end
