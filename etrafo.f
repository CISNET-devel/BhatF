c --- eigen transformations 

      subroutine etrafo(npar,eu,z)
      implicit real*8(a-h,o-z)
      parameter(nmax=100)
      dimension eu(nmax),z(nmax,nmax)
      COMMON /BOUNDS/SU(NMAX),SE(NMAX),S1(NMAX),S2(NMAX)
      
C ***   INVERT RELATIONSHIP S^{-1}

	DO I=1,NPAR
	 SU(I)=0.D0
	 DO J=1,NPAR
	  SU(I)=SU(I)+EU(J)*Z(J,I)
         ENDDO

C ***   CHECK BOUNDS

c	IF(SU(I).LE.S1(I)) PRINT*,'problem <:',i,su(i),s1(i)
c	IF(SU(I).GE.S2(I)) PRINT*,'problem >:',i,su(i),s2(i)

	ENDDO
	
      return
      end


