      dimension ivn(10)
      dimension x(10), a(10),b(10)

      do i=1,10
         a(i)=i
         b(i)=i**2
      enddo
      
      ivn(1)=10
      ivn(2)=5
      ivn(3)=1

      x(ivn(1:3))=b(1:3)/a(1:3)
      print*,x
      stop
      end
