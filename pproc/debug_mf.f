      Subroutine Debug(dQ, drYi)
  
      include 'header'

      Real*8 dQ(My,Mz,Mx,5), drYi(My,Mz,Mx,1)

      integer I, J, K

* Checking rhs value for species equation at all boundaries
      if(xrank .EQ. Px-1)then
        print *, maxval(drYi(:,:,Mx-1,1)),minval(drYi(:,:,Mx-1,1)),'x'
      endif

      if(yrank .EQ. Py-1)then
        print *, maxval(drYi(My-1,:,:,1)),minval(drYi(My-1,:,:,1)),'y'
      endif

      if(zrank .EQ. Pz-1)then
        print *, maxval(drYi(:,Mz-1,:,1)),minval(drYi(:,Mz-1,:,1)),'z'
      endif

      return

      End Subroutine Debug      
