* Adding MMS source terms (radial characteristics)
* x - boundary
      Do I = 1,Mx,Mx-1
       if(  ((xrank.EQ.Px-1).AND.(I.EQ.Mx)) .OR.
     .      ((xrank.EQ.   0).AND.(I.EQ.1))  )then
        Do J = 1,My
          Do K = 1,Mz
           r = rad(J,K,I)

           a2 = sin(r) *
     .          (gamma-1D0)*cos(r)
           a5 = cos(r) *
     .          (sqrt(gamma)+sin(r)) *
     .          (1D0+sqrt(gamma)*(2D0+sin(r)))

           d1 = (1D0/gamma)*(a2 + 5.0D-1*a5)
           d2 = 5.0D-1 * a5
           d3 = d2 / (dsqrt(gamma)*(2D0 + sin(r)) )

           Qrhs(J,K,I,4) = Qrhs(J,K,I,4) + d1
           Qrhs(J,K,I,5) = Qrhs(J,K,I,5)
     .                   + 5.0D-1*sin(r)*sin(r)*d1
     .                   + d2/(gamma-1D0)
     .                   + (2D0+sin(r))*sin(r)*d3
           Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + (xcor(I)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)
           Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + (ycor(J)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)
           Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + (zcor(K)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)

          End Do
        End Do
       endif
      End Do

* y - boundary
      Do J = 1,My,My-1
       if(  ((yrank.EQ.Py-1).AND.(J.EQ.My)) .OR.
     .      ((yrank.EQ.   0).AND.(J.EQ.1))  )then

        if(xrank .EQ. 0)then                            ! Ensuring boundaries are not evaluated twice
          I_start = 2; I_end = Mx
        elseif(xrank .EQ. Px-1)then
          I_start = 1; I_end = Mx-1
        elseif(Px .EQ. 1)then
          I_start = 2; I_end = Mx-1
        else
          I_start = 1; I_end = Mx
        endif

        Do I = I_start,I_end
         Do K = 1,Mz
           r = rad(J,K,I)

           a2 = sin(r) *
     .          (gamma-1D0)*cos(r)
           a5 = cos(r) *
     .          (sqrt(gamma)+sin(r)) *
     .          (1D0+sqrt(gamma)*(2D0+sin(r)))

           d1 = (1D0/gamma)*(a2 + 5.0D-1*a5)
           d2 = 5.0D-1 * a5
           d3 = d2 / (dsqrt(gamma)*(2D0 + sin(r)) )

           Qrhs(J,K,I,4) = Qrhs(J,K,I,4) + d1
           Qrhs(J,K,I,5) = Qrhs(J,K,I,5)
     .                   + 5.0D-1*sin(r)*sin(r)*d1
     .                   + d2/(gamma-1D0)
     .                   + (2D0+sin(r))*sin(r)*d3
           Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + (xcor(I)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)
           Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + (ycor(J)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)
           Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + (zcor(K)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)

         End Do
        End Do
       endif
      End Do


* z - boundary
      Do K = 1,Mz,Mz-1
       if(  ((zrank.EQ.Pz-1).AND.(K.EQ.Mz)) .OR.
     .      ((zrank.EQ.   0).AND.(K.EQ.1))  )then

        if(xrank .EQ. 0)then                            ! Ensuring boundaries are not evaluated twice
          I_start = 2; I_end = Mx
        elseif(xrank .EQ. Px-1)then
          I_start = 1; I_end = Mx-1
        elseif(Px .EQ. 1)then
          I_start = 2; I_end = Mx-1
        else
          I_start = 1; I_end = Mx
        endif

        if(yrank .EQ. 0)then
          J_start = 2; J_end = My
        elseif(yrank .EQ. Py-1)then
          J_start = 1; J_end = My-1
        elseif(Py .EQ. 1)then
          J_start = 2; J_end = My-1
        else
          J_start = 1; J_end = My
        endif

        Do I = I_start,I_end
          Do J = J_start,J_end
           r = rad(J,K,I)

           a2 = sin(r) *
     .          (gamma-1D0)*cos(r)
           a5 = cos(r) *
     .          (sqrt(gamma)+sin(r)) *
     .          (1D0+sqrt(gamma)*(2D0+sin(r)))

           d1 = (1D0/gamma)*(a2 + 5.0D-1*a5)
           d2 = 5.0D-1 * a5
           d3 = d2 / (dsqrt(gamma)*(2D0 + sin(r)) )

           Qrhs(J,K,I,4) = Qrhs(J,K,I,4) + d1
           Qrhs(J,K,I,5) = Qrhs(J,K,I,5)
     .                   + 5.0D-1*sin(r)*sin(r)*d1
     .                   + d2/(gamma-1D0)
     .                   + (2D0+sin(r))*sin(r)*d3
           Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + (xcor(I)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)
           Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + (ycor(J)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)
           Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + (zcor(K)/r) *
     .                     (sin(r)*d1 + (2D0+sin(r))*d3)

          End Do
        End Do
       endif
      End Do

      if(xrank.EQ.0)then
        print *, maxval(Qrhs(:,:,1,1)), minval(Qrhs(:,:,1,1)),1
        print *, maxval(Qrhs(:,:,1,2)), minval(Qrhs(:,:,1,2)),2
        print *, maxval(Qrhs(:,:,1,3)), minval(Qrhs(:,:,1,3)),3
        print *, maxval(Qrhs(:,:,1,4)), minval(Qrhs(:,:,1,4)),4
        print *, maxval(Qrhs(:,:,1,5)), minval(Qrhs(:,:,1,5)),5
      elseif(xrank.EQ.Px-1)then
        print *, maxval(Qrhs(:,:,Mx,1)), minval(Qrhs(:,:,Mx,1)),1
        print *, maxval(Qrhs(:,:,Mx,2)), minval(Qrhs(:,:,Mx,2)),2
        print *, maxval(Qrhs(:,:,Mx,3)), minval(Qrhs(:,:,Mx,3)),3
        print *, maxval(Qrhs(:,:,Mx,4)), minval(Qrhs(:,:,Mx,4)),4
        print *, maxval(Qrhs(:,:,Mx,5)), minval(Qrhs(:,:,Mx,5)),5
      endif
      if(yrank.EQ.0) then
        print *, maxval(Qrhs(1,:,:,1)), minval(Qrhs(1,:,:,1)),1
        print *, maxval(Qrhs(1,:,:,2)), minval(Qrhs(1,:,:,2)),2
        print *, maxval(Qrhs(1,:,:,3)), minval(Qrhs(1,:,:,3)),3
        print *, maxval(Qrhs(1,:,:,4)), minval(Qrhs(1,:,:,4)),4
        print *, maxval(Qrhs(1,:,:,5)), minval(Qrhs(1,:,:,5)),5
      elseif(yrank.EQ.Py-1)then
        print *, maxval(Qrhs(My,:,:,1)), minval(Qrhs(My,:,:,1)),1
        print *, maxval(Qrhs(My,:,:,2)), minval(Qrhs(My,:,:,2)),2
        print *, maxval(Qrhs(My,:,:,3)), minval(Qrhs(My,:,:,3)),3
        print *, maxval(Qrhs(My,:,:,4)), minval(Qrhs(My,:,:,4)),4
        print *, maxval(Qrhs(My,:,:,5)), minval(Qrhs(My,:,:,5)),5
      endif
      if(zrank.EQ.0)then
        print *, maxval(Qrhs(:,1,:,1)), minval(Qrhs(:,1,:,1)),1
        print *, maxval(Qrhs(:,1,:,2)), minval(Qrhs(:,1,:,2)),2
        print *, maxval(Qrhs(:,1,:,3)), minval(Qrhs(:,1,:,3)),3
        print *, maxval(Qrhs(:,1,:,4)), minval(Qrhs(:,1,:,4)),4
        print *, maxval(Qrhs(:,1,:,5)), minval(Qrhs(:,1,:,5)),5
      elseif(zrank.EQ.Pz-1)then
        print *, maxval(Qrhs(:,Mz,:,1)), minval(Qrhs(:,Mz,:,1)),1
        print *, maxval(Qrhs(:,Mz,:,2)), minval(Qrhs(:,Mz,:,2)),2
        print *, maxval(Qrhs(:,Mz,:,3)), minval(Qrhs(:,Mz,:,3)),3
        print *, maxval(Qrhs(:,Mz,:,4)), minval(Qrhs(:,Mz,:,4)),4
        print *, maxval(Qrhs(:,Mz,:,5)), minval(Qrhs(:,Mz,:,5)),5
      endif

* Adding MMS source terms (x,y,z- characteristics)
* x - boundary
      Do I = 1,Mx,Mx-1
       if(  ((xrank.EQ.Px-1).AND.(I.EQ.Mx)) .OR.
     .      ((xrank.EQ.   0).AND.(I.EQ.1))  )then
        Do J = 1,My
          Do K = 1,Mz

           if(xrank .EQ. Px-1)then
            Lsource(1) = 0D0
            Lsource(2) = 0D0
            Lsource(3) = 0D0
            Lsource(4) = 0D0
            Lsource(5) = cos(xcor(I))*(sin(xcor(I))+dsqrt(gamma))
     .                 * ( sin(ycor(J))*sin(zcor(K))
     .             + (2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))
     .               *dsqrt(gamma) )
           else
            Lsource(1) = cos(xcor(I))*(sin(xcor(I))-dsqrt(gamma))
     .                 * ( sin(ycor(J))*sin(zcor(K))
     .             - (2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))
     .               *dsqrt(gamma) )
            Lsource(2) = 0D0
            Lsource(3) = 0D0
            Lsource(4) = 0D0
            Lsource(5) = 0D0
           endif

           d(J,K,I,1) = d(J,K,I,1) - (1D0/gamma) *
     .          ( Lsource(2) + 5.0D-1*(Lsource(1)+Lsource(5)) )
           d(J,K,I,2) = d(J,K,I,2) - 5.0D-1*(Lsource(5)+Lsource(1))
           d(J,K,I,3) = d(J,K,I,3) -
     .                 (5.0D-1/
     .     ((2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))*dsqrt(gamma)))
     .                 *(Lsource(5)-Lsource(1))
           d(J,K,I,4) = d(J,K,I,4) - Lsource(3)
           d(J,K,I,5) = d(J,K,I,5) - Lsource(4)
          End Do
        End Do
       endif
      End Do

*y-boundary
      Do J = 1,My,My-1
       if(  ((yrank.EQ.Py-1).AND.(J.EQ.My)) .OR.
     .      ((yrank.EQ.   0).AND.(J.EQ.1))  )then
        Do I = 1,Mx
         Do K = 1,Mz

            if(yrank .EQ. Py-1)then
             Msource(1) = 0D0
             Msource(2) = 0D0
             Msource(3) = 0D0
             Msource(4) = 0D0
             Msource(5) = cos(ycor(J))*(sin(ycor(J))+dsqrt(gamma))
     .                 * ( sin(xcor(I))*sin(zcor(K))
     .             + (2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))
     .               *dsqrt(gamma) )
            else
             Msource(1) = cos(ycor(J))*(sin(ycor(J))-dsqrt(gamma))
     .                 * ( sin(xcor(I))*sin(zcor(K))
     .             - (2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))
     .               *dsqrt(gamma) )
             Msource(2) = 0D0
             Msource(3) = 0D0
             Msource(4) = 0D0
             Msource(5) = 0D0
            endif

           e(J,K,I,1) = e(J,K,I,1) - (1D0/gamma) *
     .          (Msource(2) + 5.0D-1*(Msource(1)+Msource(5)) )
           e(J,K,I,2) = e(J,K,I,2) - 5.0D-1*(Msource(5)+Msource(1))
           e(J,K,I,3) = e(J,K,I,3) -
     .                 (5.0D-1/
     .     ((2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))*dsqrt(gamma)))
     .                 *(Msource(5)-Msource(1))
           e(J,K,I,4) = e(J,K,I,4) - Msource(3)
           e(J,K,I,5) = e(J,K,I,5) - Msource(4)

         End Do
        End Do

c        print*, maxval(e(J,:,:,1)),minval(e(J,:,:,1)),'e1'
c        print*, maxval(e(J,:,:,2)),minval(e(J,:,:,2)),'e2'
c        print*, maxval(e(J,:,:,3)),minval(e(J,:,:,3)),'e3'
c        print*, maxval(e(J,:,:,4)),minval(e(J,:,:,4)),'e4'
c        print*, maxval(e(J,:,:,5)),minval(e(J,:,:,5)),'e5'

       endif
      End Do

* z - boundary
      Do K = 1,Mz,Mz-1
       if(  ((zrank.EQ.Pz-1).AND.(K.EQ.Mz)) .OR.
     .      ((zrank.EQ.   0).AND.(K.EQ.1))  )then
        Do I = 1,Mx
          Do J = 1,My

            if(zrank .EQ. Pz-1)then
              Nsource(1) = 0D0
              Nsource(2) = 0D0
              Nsource(3) = 0D0
              Nsource(4) = 0D0
              Nsource(5) = cos(zcor(K))*(sin(zcor(K))+dsqrt(gamma))
     .                   * ( sin(ycor(J))*sin(xcor(I))
     .               + (2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))
     .                 *dsqrt(gamma) )
           else
              Nsource(1) = cos(zcor(K))*(sin(zcor(K))-dsqrt(gamma))
     .                   * ( sin(ycor(J))*sin(xcor(I))
     .               - (2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))
     .                 *dsqrt(gamma) )
              Nsource(2) = 0D0
              Nsource(3) = 0D0
              Nsource(4) = 0D0
              Nsource(5) = 0D0
           endif

           f(J,K,I,1) = f(J,K,I,1) - (1D0/gamma) *
     .          ( Nsource(2) + 5.0D-1*(Nsource(1)+Nsource(5)) )
           f(J,K,I,2) = f(J,K,I,2) - 5.0D-1*(Nsource(5)+Nsource(1))
           f(J,K,I,3) = f(J,K,I,3) -
     .                 (5.0D-1/
     .     ((2D0+sin(xcor(I))*sin(ycor(J))*sin(zcor(K)))*dsqrt(gamma)))
     .                 *(Nsource(5)-Nsource(1))
           f(J,K,I,4) = f(J,K,I,4) - Nsource(3)
           f(J,K,I,5) = f(J,K,I,5) - Nsource(4)

          End Do
         End Do

        print*, maxval(f(:,K,:,1)),minval(f(:,K,:,1)),'f1'
        print*, maxval(f(:,K,:,2)),minval(f(:,K,:,2)),'f2'
        print*, maxval(f(:,K,:,3)),minval(f(:,K,:,3)),'f3'
        print*, maxval(f(:,K,:,4)),minval(f(:,K,:,4)),'f4'
        print*, maxval(f(:,K,:,5)),minval(f(:,K,:,5)),'f5'

        endif
       End Do

