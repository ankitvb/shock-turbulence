************************************************************************
* Subroutine: Set_BC                                                   *
*                                                                      *
* Function: Set non reflecting boundary conditions for spherical       *
*           shock-turbulence interaction. Subsonic boundary,           *
*           outgoing waves only, incoming wave amplitudes set to zero  *
*           Sponge damping applied to relax the flow to steady ref     *                                                          
*           state                                                      *
*           Assumes origin of coordinate system is at center of domain *
************************************************************************

      Subroutine bc_RHS(Q, Qrhs, time)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8 P(My,Mz,Mx), U(My,Mz,Mx,3)
      Real*8 Ur(My,Mz,Mx), rad(My,Mz,Mx), dF(My,Mz,Mx,9)
      Real*8 Vs(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 L(My,Mz,Mx,5)
      Real*8 x(My,Mz,Mx), y(My,Mz,Mx), z(My,Mz,Mx)
      Real*8 time, r

      Real*8 fac1, fac2, fac3, fac, Ms
      Real*8 tmp2, d1, d2, d3, a2, a5
      Real*8 dPdr, dUrdr, drhodr

      Integer I, J, K, M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ
      Integer I_start, I_end, J_start, J_end
     
*  Compute specific volume (avoid division in computations), pressure,
*  temperature and transverse velocity
*
*   Vs  = 1 / rho
*
*    P  = (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*           gamma1     P        1
*    T  = ---------   ---  ------------
*         gamma1 - 1  rho  (Y1 + eps*Y2)
*

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
             P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                          ( Q(J,K,I,1) * Q(J,K,I,1)
     .                          + Q(J,K,I,2) * Q(J,K,I,2)
     .                          + Q(J,K,I,3) * Q(J,K,I,3) ) )
             T(J,K,I)   = fac2 * Vs(J,K,I) * P(J,K,I)

             U(J,K,I,1) = Q(J,K,I,1) * Vs(J,K,I)
             U(J,K,I,2) = Q(J,K,I,2) * Vs(J,K,I)
             U(J,K,I,3) = Q(J,K,I,3) * Vs(J,K,I)

             II = xrank * Mx + I
             JJ = yrank * My + J
             KK = zrank * Mz + K          

             IX = REAL(Nx/2 - II) + 0.5D0
             IY = REAL(Ny/2 - JJ) + 0.5D0
             IZ = REAL(Nz/2 - KK) + 0.5D0

             rad(J,K,I) = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

             Ur(J,K,I) = (xcor(I)/rad(J,K,I)) * U(J,K,I,1)
     .                 + (ycor(J)/rad(J,K,I)) * U(J,K,I,2)
     .                 + (zcor(K)/rad(J,K,I)) * U(J,K,I,3) 
          End Do
        End Do
      End Do

* Precomputing derivatives for later use
      Call DFDX(P,  dF(1,1,1,5), 0, 1.0D0)
      Call DFDX(Ur, dF(1,1,1,6), 0, 1.0D0)
      Call DFDX(Q(1,1,1,4), dF(1,1,1,7), 0, 1.0D0)
      Call DFDY(P,  dF(1,1,1,1), 0, 1.0D0)
      Call DFDY(Ur, dF(1,1,1,2), 0, 1.0D0) 
      Call DFDY(Q(1,1,1,4), dF(1,1,1,8), 0, 1.0D0)
      Call DFDZ(P,  dF(1,1,1,3), 0, 1.0D0)
      Call DFDZ(Ur, dF(1,1,1,4), 0, 1.0D0)
      Call DFDZ(Q(1,1,1,4), dF(1,1,1,9), 0, 1.0D0)

* Calculating the wave amplitude variation for incoming/outgoing waves
* x - boundary
      Do I = 1,Mx,Mx-1
       if(  ((xrank.EQ.Px-1).AND.(I.EQ.Mx)) .OR. 
     .      ((xrank.EQ.   0).AND.(I.EQ.1))  )then
        Do J = 1,My
          Do K = 1,Mz
           dPdr   = (xcor(I)/rad(J,K,I)) * dF(J,K,I,5)
     .            + (ycor(J)/rad(J,K,I)) * dF(J,K,I,1)
     .            + (zcor(K)/rad(J,K,I)) * dF(J,K,I,3)
           dUrdr  = (xcor(I)/rad(J,K,I)) * dF(J,K,I,6)
     .            + (ycor(J)/rad(J,K,I)) * dF(J,K,I,2)
     .            + (zcor(K)/rad(J,K,I)) * dF(J,K,I,4)
           drhodr = (xcor(I)/rad(J,K,I)) * dF(J,K,I,7)
     .            + (ycor(J)/rad(J,K,I)) * dF(J,K,I,8)
     .            + (zcor(K)/rad(J,K,I)) * dF(J,K,I,9)

           tmp2 = dsqrt(gamma*Vs(J,K,I)*P(J,K,I))                ! Sound speed
 
           L(J,K,I,1) = 0D0                                      ! Incoming characteristic (acoustic)
           L(J,K,I,2) = Ur(J,K,I) *                              ! Outgoing characteristic (entropy)
     .                  (tmp2*tmp2*drhodr - dPdr)
           L(J,K,I,5) = (Ur(J,K,I) + tmp2) *                     ! Outgoing characteristic (acoustic)
     .                  (dPdr + Q(J,K,I,4) * tmp2 * dUrdr)

           d1 = (1D0/(tmp2*tmp2)) *
     .          ( L(J,K,I,2) + 5.0D-1*(L(J,K,I,1)+L(J,K,I,5)) )
           d2 = 5.0D-1*(L(J,K,I,5)+L(J,K,I,1))
           d3 = (5.0D-1/tmp2)*(L(J,K,I,5)-L(J,K,I,1))/Q(J,K,I,4)

           Qrhs(J,K,I,4) = -d1
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d1
     .                   - d2/(gamma-1D0) 
     .                   - Q(J,K,I,4)*Ur(J,K,I)*d3
           Qrhs(J,K,I,1) = -(xcor(I)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3) 
           Qrhs(J,K,I,2) = -(ycor(J)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
           Qrhs(J,K,I,3) = -(zcor(K)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
          
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
           dPdr  = (xcor(I)/rad(J,K,I)) * dF(J,K,I,5)
     .           + (ycor(J)/rad(J,K,I)) * dF(J,K,I,1)
     .           + (zcor(K)/rad(J,K,I)) * dF(J,K,I,3)
           dUrdr = (xcor(I)/rad(J,K,I)) * dF(J,K,I,6)
     .           + (ycor(J)/rad(J,K,I)) * dF(J,K,I,2)
     .           + (zcor(K)/rad(J,K,I)) * dF(J,K,I,4)
           drhodr = (xcor(I)/rad(J,K,I)) * dF(J,K,I,7)
     .            + (ycor(J)/rad(J,K,I)) * dF(J,K,I,8)
     .            + (zcor(K)/rad(J,K,I)) * dF(J,K,I,9)

           tmp2 = sqrt(gamma*Vs(J,K,I)*P(J,K,I))                 ! Sound speed

           L(J,K,I,1) = 0D0                                      ! Incoming characteristic (acoustic)
           L(J,K,I,2) = Ur(J,K,I) *                              ! Outgoing characteristic (entropy)
     .                  (tmp2*tmp2*drhodr - dPdr)
           L(J,K,I,5) = (Ur(J,K,I) + tmp2) *                     ! Outgoing characteristic (acoustic)
     .                  (dPdr + Q(J,K,I,4) * tmp2 * dUrdr)

           d1 = (1D0/(tmp2*tmp2)) *
     .          ( L(J,K,I,2) + 5.0D-1*(L(J,K,I,1)+L(J,K,I,5)) )
           d2 = 5.0D-1*(L(J,K,I,5)+L(J,K,I,1))
           d3 = (5.0D-1/tmp2)*(L(J,K,I,5)-L(J,K,I,1))/Q(J,K,I,4)

           Qrhs(J,K,I,4) = -d1
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d1
     .                   - d2/(gamma-1D0)
     .                   - Q(J,K,I,4)*Ur(J,K,I)*d3
           Qrhs(J,K,I,1) = -(xcor(I)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
           Qrhs(J,K,I,2) = -(ycor(J)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
           Qrhs(J,K,I,3) = -(zcor(K)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
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
           dPdr  = (xcor(I)/rad(J,K,I)) * dF(J,K,I,5)
     .           + (ycor(J)/rad(J,K,I)) * dF(J,K,I,1)
     .           + (zcor(K)/rad(J,K,I)) * dF(J,K,I,3)
           dUrdr = (xcor(I)/rad(J,K,I)) * dF(J,K,I,6)
     .           + (ycor(J)/rad(J,K,I)) * dF(J,K,I,2)
     .           + (zcor(K)/rad(J,K,I)) * dF(J,K,I,4)
           drhodr = (xcor(I)/rad(J,K,I)) * dF(J,K,I,7)
     .            + (ycor(J)/rad(J,K,I)) * dF(J,K,I,8)
     .            + (zcor(K)/rad(J,K,I)) * dF(J,K,I,9)

           tmp2 = sqrt(gamma*Vs(J,K,I)*P(J,K,I))                 ! Sound speed

           L(J,K,I,1) = 0D0                                      ! Incoming characteristic (acoustic)
           L(J,K,I,2) = Ur(J,K,I) *                              ! Outgoing characteristic (entropy)
     .                  (tmp2*tmp2*drhodr - dPdr)
           L(J,K,I,5) = (Ur(J,K,I) + tmp2) *                     ! Outgoing characteristic (acoustic)
     .                  (dPdr + Q(J,K,I,4) * tmp2 * dUrdr)

           d1 = (1D0/(tmp2*tmp2)) *
     .          ( L(J,K,I,2) + 5.0D-1*(L(J,K,I,1)+L(J,K,I,5)) )
           d2 = 5.0D-1*(L(J,K,I,5)+L(J,K,I,1))
           d3 = (5.0D-1/tmp2)*(L(J,K,I,5)-L(J,K,I,1))/Q(J,K,I,4)

           Qrhs(J,K,I,4) = -d1
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d1
     .                   - d2/(gamma-1D0)
     .                   - Q(J,K,I,4)*Ur(J,K,I)*d3
           Qrhs(J,K,I,1) = -(xcor(I)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
           Qrhs(J,K,I,2) = -(ycor(J)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
           Qrhs(J,K,I,3) = -(zcor(K)/rad(J,K,I)) *
     .                     (Ur(J,K,I)*d1 + Q(J,K,I,4)*d3)
          End Do
        End Do
       endif 
      End Do

      End Subroutine bc_RHS
