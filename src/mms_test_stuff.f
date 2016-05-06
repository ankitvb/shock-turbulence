cWarning: Adding source terms for MMS , TO BE REMOVED
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + sin(2*ycor(J))
     .                     + cos(zcor(K))*sin(ycor(J))
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + sin(2*zcor(K))
     .                     + cos(ycor(J))*sin(zcor(K))
            Qrhs(J,K,I,4) = Qrhs(J,K,I,4) + cos(zcor(K))
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + fac2*P(J,K,I)*cos(zcor(K))
     .                      + (3.0/4.0)*sin(zcor(K))*sin(2*zcor(K))

          End Do
        End Do
      End Do

* MMS source terms (4/18/08)
* rho, P const, u,v,w varying

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
           fac2 = gamma / (gamma - 1D0)
           P = 1.0D2
           Qrhs(J,K,I,4) = Qrhs(J,K,I,4) + cos(xcor(I)+Pi/6D0)
     .   + cos(ycor(J)+Pi/3D0) + cos(zcor(K)+Pi/2D0)
           Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + sin(2*xcor(I)+Pi/3D0)
     .   + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)
     .   + sin(xcor(I)+Pi/6D0)*cos(zcor(K)+Pi/2D0)
           Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + sin(2*ycor(J)+2D0*Pi/3D0)
     .   + sin(ycor(J)+Pi/3D0)*cos(xcor(I)+Pi/6D0)
     .   + sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
           Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + sin(2*zcor(K)+Pi)
     .   + sin(zcor(K)+Pi/2D0)*cos(xcor(I)+Pi/6D0)
     .   + sin(zcor(K)+Pi/2D0)*cos(ycor(J)+Pi/3D0)
           Qrhs(J,K,I,5) = Qrhs(J,K,I,5)
     .   + fac2*P(J,K,I)*(cos(xcor(I)+Pi/6D0)
     .   + cos(ycor(J)+Pi/3D0) + cos(zcor(K)+Pi/2D0) )
     .   + (3D0/2D0)
     .   *(sin(xcor(I)+Pi/6D0)*sin(xcor(I)+Pi/6D0)*cos(xcor(I)+Pi/6D0)
     .   + sin(ycor(J)+Pi/3D0)*sin(ycor(J)+Pi/3D0)*cos(ycor(J)+Pi/3D0)
     .   + sin(zcor(K)+Pi/2D0)*sin(zcor(K)+Pi/2D0)*cos(zcor(K)+Pi/2D0))
     .   + 5.0D-1 *
     .   + (sin(ycor(J)+Pi/3D0)*sin(ycor(J)+Pi/3D0)*cos(xcor(I)+Pi/6D0)
     .   +  sin(zcor(K)+Pi/2D0)*sin(zcor(K)+Pi/2D0)*cos(xcor(I)+Pi/6D0)
     .   +  sin(xcor(I)+Pi/6D0)*sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)
     .   +  sin(zcor(K)+Pi/2D0)*sin(zcor(K)+Pi/2D0)*cos(ycor(J)+Pi/3D0)
     .   +  sin(xcor(I)+Pi/6D0)*sin(xcor(I)+Pi/6D0)*cos(zcor(K)+Pi/2D0)
     .   +  sin(ycor(J)+Pi/3D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0))
          End Do
        End Do
      End Do

* MMS source terms(4/20/08)
* u, v, w constant,P, rho varying

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
          fac2 = gamma / (gamma - 1D0)
          Qrhs(J,K,I,4) = Qrhs(J,K,I,4)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,1) = Qrhs(J,K,I,1)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,2) = Qrhs(J,K,I,2)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,3) = Qrhs(J,K,I,3)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + (3D0/2D0 + fac2)
     .   *(
     .      cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
     .    )
          End Do
        End Do
      End Do

* MMS source terms (4/25/08)
* u,v,w,P constant, rho varying

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
          Qrhs(J,K,I,4) = Qrhs(J,K,I,4)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,1) = Qrhs(J,K,I,1)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,2) = Qrhs(J,K,I,2)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,3) = Qrhs(J,K,I,3)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + (3D0/2D0)
     .   *(
     .      cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
     .    )
          End Do
        End Do
      End Do

* MMS cancelling terms(4/25/08)
* rho, u, v, w constant, P varying

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
          Qrhs(J,K,I,4) = Qrhs(J,K,I,4)
          Qrhs(J,K,I,1) = Qrhs(J,K,I,1)
     .    + cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,2) = Qrhs(J,K,I,2)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,3) = Qrhs(J,K,I,3)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
          Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + (fac2)
     .   *(
     .      cos(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*cos(ycor(J)+Pi/3D0)*sin(zcor(K)+Pi/2D0)
     .    + sin(xcor(I)+Pi/6D0)*sin(ycor(J)+Pi/3D0)*cos(zcor(K)+Pi/2D0)
     .    )
          End Do
        End Do
      End Do

* Adding MMS source terms for viscous part only
* 5/28/08

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz

           Qrhs(J,K,I,4) = 0D0
           Qrhs(J,K,I,1) = Qrhs(J,K,I,1)
     .     +       4D0*(sin(xcor(I))*cos(ycor(J))*cos(zcor(K)))
           Qrhs(J,K,I,2) = Qrhs(J,K,I,2)
     .     +       4D0*(cos(xcor(I))*sin(ycor(J))*cos(zcor(K)))
           Qrhs(J,K,I,3) = Qrhs(J,K,I,3)
     .     +       4D0*(cos(xcor(I))*cos(ycor(J))*sin(zcor(K)))
           Qrhs(J,K,I,5) = Qrhs(J,K,I,5)
     .     + (4D0/3D0)*(sin(xcor(I))*cos(ycor(J))*cos(zcor(K)))**2D0
     .     +           (cos(xcor(I))*sin(ycor(J))*cos(zcor(K)))**2D0
     .     +           (cos(xcor(I))*cos(ycor(J))*sin(zcor(K)))**2D0
     .     +           (sin(xcor(I))*cos(ycor(J))*cos(zcor(K)))**2D0
     .     + (4D0/3D0)*(cos(xcor(I))*sin(ycor(J))*cos(zcor(K)))**2D0
     .     +           (cos(xcor(I))*cos(ycor(J))*sin(zcor(K)))**2D0
     .     +           (sin(xcor(I))*cos(ycor(J))*cos(zcor(K)))**2D0
     .     +           (cos(xcor(I))*sin(ycor(J))*cos(zcor(K)))**2D0
     .     + (4D0/3D0)*(cos(xcor(I))*cos(ycor(J))*sin(zcor(K)))**2D0
     .     + (2D0/3D0)*(sin(xcor(I))*cos(ycor(J))*cos(zcor(K)))**2D0
     .     + (2D0/3D0)*(cos(xcor(I))*sin(ycor(J))*cos(zcor(K)))**2D0
     .     + (2D0/3D0)*(cos(xcor(I))*cos(ycor(J))*sin(zcor(K)))**2D0
     .     -       4D0*(sin(xcor(I))*sin(ycor(J))*cos(zcor(K)))**2D0
     .     -       4D0*(sin(xcor(I))*cos(ycor(J))*sin(zcor(K)))**2D0
     .     -       4D0*(cos(xcor(I))*sin(ycor(J))*sin(zcor(K)))**2D0

          End Do
        End Do
      End Do

