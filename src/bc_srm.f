************************************************************************
* Subroutine: Set_BC                                                   *
*                                                                      *
* Function: Set boundary conditions for spherical Richtmeyer-Meshkov   *
*                                                                      *
*           Left boundary : Symmetric boundary                         *
*                                                                      *
*                        d f                                           *
*                  v_n,  --- = 0   for all quantities f                *
*                        d n       & velocity normal to face           *
*                                                                      *
*           Right boundary: Subsonic boundary, outgoing waves only     *
*                           incoming wave amplitudes set to zero       *
*                                                                      *
************************************************************************

      Subroutine bc_RHS(Q, Qrhs, Y, time)

      Include 'header'

 
      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8 P(My,Mz,Mx), gamma_eff(My,Mz,Mx), U(My,Mz,Mx,3)
      Real*8 Ur(My,Mz,Mx), rad(My,Mz,Mx), dF(My,Mz,Mx,6)
      Real*8 Y(My,Mz,Mx,1), Vs(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 L(My,Mz,Mx,5)
      Real*8 time

      Common /THERMO/ gamma_eff  

      Real*8 fac1, fac2, fac3, fac, Ms
      Real*8 tmp2, d1, d2, d3
      Real*8 dPdr, dUrdr
      Real*8 base_pressure
      Real*8 P_postshock, P_preshock, P_ratio

      Integer I,J,K,M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

      Ms = 1.24D0 ! Initial shock mach number  
      base_pressure = 0.71429D0

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

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz

             fac1 = gamma_eff(J,K,I) - 1.0D0
             fac2 = gamma_eff(J,K,I) / fac1
             fac3 = gamma1 / (gamma1 - 1D0)
             fac3 = fac3 *
     .       1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

            Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
             P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
             T(J,K,I)   = fac3 * Vs(J,K,I) * P(J,K,I)
             U(J,K,I,1) = Q(J,K,I,1) * Vs(J,K,I)
             U(J,K,I,2) = Q(J,K,I,2) * Vs(J,K,I)
             U(J,K,I,3) = Q(J,K,I,3) * Vs(J,K,I)
             Ur(J,K,I)  = dsqrt(U(J,K,I,1) * U(J,K,I,1) 
     .                       + U(J,K,I,2) * U(J,K,I,2)
     .                       + U(J,K,I,3) * U(J,K,I,3) ) 

              II = rank * Mx + I
              IX = REAL(II) - 1D0
              IY = REAL(J)  - 1D0
              IZ = REAL(K)  - 1D0

              rad(J,K,I) = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

          End Do
        End Do
      End Do
  
* Precomputing derivatives for later use
      Call DFDX_sym(P,  dF(1,1,1,5), 0, 1.0D0)
      Call DFDX_sym(Ur, dF(1,1,1,6), 0, 1.0D0)
      Call DFDY_sym(P,  dF(1,1,1,1), 0, 1.0D0)
      Call DFDY_sym(Ur, dF(1,1,1,2), 0, 1.0D0) 
      Call DFDZ_sym(P,  dF(1,1,1,3), 0, 1.0D0)
      Call DFDZ_sym(Ur, dF(1,1,1,4), 0, 1.0D0)

c      Call Compute_Max(P, P_postshock) 
c      P_preshock  = base_pressure
c      P_ratio     = 0.99D0 * P_postshock / P_preshock

* Shock Mach number
c           Ms = sqrt((gamma1-1D0)/(2D0*gamma1) +
c     .               (gamma1+1D0)/(2D0*gamma1)*P_ratio)

c      if(rank .EQ. (Px-1))print *, Ms

* Make velocities normal to inner faces zero
* 
* Plane x = 0 ; u_x = 0
* Plane y = 0 ; u_y = 0
* Plane z = 0 ; u_z = 0
*

      if(xrank .EQ. 0)then
       Do J = 1,My
         Do K = 1,Mz
          Q(J,K,1,1) = 0D0
        End Do
       End Do
      endif

      if(yrank .EQ. 0)then
       Do I = 1,Mx
         Do K = 1,Mz
           Q(1,K,I,2) = 0D0
         End Do
       End Do
      endif

      if(zrank .EQ. 0)then
       Do I = 1,Mx
         Do J = 1,My
           Q(J,1,I,3) = 0D0
        End Do
       End Do
      endif

* Calculating the wave amplitude variation for incoming acoustic waves
* x - boundary
      if(xrank .EQ. Px-1)then
        I = Mx
        Do J = 1,My
          Do K = 1,Mz
           JJ = yrank * My + J
           KK = zrank * Mz + K
           if((JJ. EQ. 1) .AND. (KK .EQ. 1))then
             dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
             dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
           elseif((KK .EQ. 1) .AND. (JJ .GT. 1))then
             dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .             + (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
             dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .             + (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
           elseif((JJ .EQ. 1) .AND. (KK .GT. 1))then
             dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .             + (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
             dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .             + (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
           else
             dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .             + (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
     .             + (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
             dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .             + (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
     .             + (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
           endif

           tmp2 = dsqrt(gamma_eff(J,K,I)*Vs(J,K,I)*P(J,K,I))         ! Sound speed
  
           L(J,K,I,1) = 0D0                                      ! Incoming characteristic
           L(J,K,I,2) = 0D0                                      ! Incoming characteristic
           L(J,K,I,5) = (Ur(J,K,I) + tmp2) *                     ! Outgoing characteristic 
     .                    (dPdr + Q(J,K,I,4) * tmp2 * dUrdr)


           d1 = (1D0/(tmp2*tmp2)) *
     .          ( L(J,K,I,2) + 5.0D-1*(L(J,K,I,1)+L(J,K,I,5)) )
           d2 = 5.0D-1*(L(J,K,I,5)+L(J,K,I,1))
           d3 = 5.0D-1*(L(J,K,I,5)-L(J,K,I,1))

           Qrhs(J,K,I,4) = -d1
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)
     .                   - d2/(gamma_eff(J,K,I)-1D0) 
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

* y - boundary
      if(yrank .EQ. Py-1)then
        J = My
        Do I = 1,Mx
          Do K = 1,Mz
             II = xrank * Mx + I
             KK = zrank * Mz + K
             if((II. EQ. 1) .AND. (KK .EQ. 1))then
               dPdr  = (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
               dUrdr = (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
             elseif((KK .EQ. 1) .AND. (II .GT. 1))then
               dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
               dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
             elseif((II .EQ. 1) .AND. (KK .GT. 1))then
               dPdr  = (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
               dUrdr = (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
             else
               dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
               dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
             endif

            tmp2 = dsqrt(gamma_eff(J,K,I)*Vs(J,K,I)*P(J,K,I))         ! Sound speed
 
           L(J,K,I,1) = 0D0                                      ! Incoming characteristic
           L(J,K,I,2) = 0D0                                      ! Incoming characteristic
           L(J,K,I,5) = (Ur(J,K,I) + tmp2) *                     ! Outgoing characteristic
     .                    (dPdr + Q(J,K,I,4) * tmp2 * dUrdr)

           d1 = (1D0/(tmp2*tmp2)) *
     .          ( L(J,K,I,2) + 5.0D-1*(L(J,K,I,1)+L(J,K,I,5)) )
           d2 = 5.0D-1*(L(J,K,I,5)+L(J,K,I,1))
           d3 = 5.0D-1*(L(J,K,I,5)-L(J,K,I,1))

           Qrhs(J,K,I,4) = -d1
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)
     .                   - d2/(gamma_eff(J,K,I)-1D0)
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

* z - boundary
      if(zrank .EQ. Pz-1)then
        K = Mz
        Do I = 1,Mx
          Do J = 1,My
             II = xrank * Mx + I 
             JJ = yrank * My + J
             if((JJ. EQ. 1) .AND. (II .EQ. 1))then
               dPdr  = (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
               dUrdr = (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
             elseif((II .EQ. 1) .AND. (JJ .GT. 1))then
               dPdr  = (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
               dUrdr = (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
             elseif((JJ .EQ. 1) .AND. (II .GT. 1))then
               dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
               dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
             else
               dPdr  = (rad(J,K,I)/xcor(I)) * dF(J,K,I,5)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,1)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,3)
               dUrdr = (rad(J,K,I)/xcor(I)) * dF(J,K,I,6)
     .               + (rad(J,K,I)/ycor(J)) * dF(J,K,I,2)
     .               + (rad(J,K,I)/zcor(K)) * dF(J,K,I,4)
             endif

           tmp2 = dsqrt(gamma_eff(J,K,I)*Vs(J,K,I)*P(J,K,I))         ! Sound speed

           L(J,K,I,1) = 0D0                                      ! Incoming characteristic
           L(J,K,I,2) = 0D0                                      ! Incoming characteristic
           L(J,K,I,5) = (Ur(J,K,I) + tmp2) *                     ! Outgoing characteristic
     .                    (dPdr + Q(J,K,I,4) * tmp2 * dUrdr)

           d1 = (1D0/(tmp2*tmp2)) *
     .          ( L(J,K,I,2) + 5.0D-1*(L(J,K,I,1)+L(J,K,I,5)) )
           d2 = 5.0D-1*(L(J,K,I,5)+L(J,K,I,1))
           d3 = 5.0D-1*(L(J,K,I,5)-L(J,K,I,1))

           Qrhs(J,K,I,4) = -d1
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)
     .                   - d2/(gamma_eff(J,K,I)-1D0)
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

      End Subroutine bc_RHS
