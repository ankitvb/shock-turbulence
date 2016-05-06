************************************************************************
* Subroutine: bc_rhs                                                   *
*                                                                      *
* Function: Set the boundary conditions at the end of each             *
*           RK4 substep for shock-acoustic wave interaction            *
*                                                                      *
*           Left boundary : Supersonic inflow, all variables           *
*                           specified                                  *
*           Right boundary: One incoming acoustic wave, subsonic bdy.  *
*                                                                      *
************************************************************************

      Subroutine bc_RHS(Q, Qrhs, time)

      Include 'header'

      Real*8 Q(Ny,Nz,Mx,5), Qrhs(Ny,Nz,Mx,5)
      Real*8 time
      Real*8 U(Nx,Nz,My,3), dF(Nx,Nz,My,5)
      Real*8 Av, Ae, psi, kx, ky, U1, M1
      Real*8 tmp2

      Real*8 fac1, fac2

      Integer I, J, K, M, JJ

* Left boundary - supersonic flow, all variables specified
      Av = 2.5D-2
      Ae = 2.5D-2

c      psi = atan(1D0) ! 45 deg
      psi = 1.30899693899575D0 ! 75 deg

      ky = 2D0
      kx = ky / tan(psi)

      M1 = 1.5D0
      U1 = M1

* Inlet boundary : Supersonic inflow with disturbances
      if(xrank .EQ. 1)then
        I = 1
        Do J = 1,My
          Do K = 1,Mz

            JJ = yrank * My + J

            Q(J,K,1,4) = 1D0  + 1D0 * Ae
     .        * cos( ky*ycor(JJ) - U1*kx*time )
            Q(J,K,1,1) = Q(J,K,1,4)
     .        * ( U1 + U1 * Av * sin(psi)
     .           * cos( ky*ycor(JJ) - U1*kx*time ) )
            Q(J,K,1,2) = -Q(J,K,1,4)
     .        * ( U1 * Av * cos(psi)
     .            * cos(ky*ycor(JJ) - U1*kx*time ) )
            Q(J,K,1,3) = 0D0
            tmp2 = 0.714286D0
            Q(J,K,1,5) = tmp2/ (gamma - 1D0)
     .                   + (5D-1 / Q(J,K,1,4)) *
     .                   (   Q(J,K,1,1) * Q(J,K,1,1)
     .                     + Q(J,K,1,2) * Q(J,K,1,2) )

          End Do
        End Do
      endif

      Return

      End Subroutine bc_RHS

