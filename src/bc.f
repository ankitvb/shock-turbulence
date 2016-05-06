************************************************************************
* Subroutine: Set_BC                                                   *
*                                                                      *
* Function: Set non reflecting boundary conditions for spherical       *
*           shock-turbulence interaction. Subsonic boundary,           *
*           outgoing waves only, incoming wave amplitudes set to zero  *
*           Sponge damping applied to relax the flow to steady ref     *                                                          
*           state                                                      *
************************************************************************

      Subroutine bc_RHS(Q, Qrhs, time)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8 P(My,Mz,Mx), U(My,Mz,Mx,3), Ur(My,Mz,Mx)
      Real*8 dF(My,Mz,Mx,3)
      Real*8 Vs(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 L(5), M(5), N(5)
      Real*8 time
 
      Real*8 Qconvx(My,Mz,Mx,5), Qconvy(My,Mz,Mx,5), Qconvz(My,Mz,Mx,5)
      Real*8 Qvisc(My,Mz,Mx,5)

      Common /RHSV/ Qconvx, Qconvy, Qconvz, Qvisc

      Real*8 Qh(My,Mz,Mx,5)

      Common /RHSVH/ Qh

      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)

      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz

      Real*8 fac1, fac2, fac3, fac, Ms
      Real*8 tmp2
      Real*8 d(My,Mz,Mx,5), e(My,Mz,Mx,5), f(My,Mz,Mx,5)

      Integer I, J, K
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

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
*           gamma1     P  
*    T  = ---------   ---
*         gamma1 - 1  rho 
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

             Ur(J,K,I) = dsqrt(U(J,K,I,1)*U(J,K,I,1)
     .                       + U(J,K,I,2)*U(J,K,I,2)
     .                       + U(J,K,I,3)*U(J,K,I,3))

              II = xrank * Mx + I
              JJ = yrank * My + J
              KK = zrank * Mz + K          

              IX = REAL(Nx/2 - II) + 0.5D0
              IY = REAL(Ny/2 - JJ) + 0.5D0
              IZ = REAL(Nz/2 - KK) + 0.5D0

          End Do
        End Do
      End Do

* Precomputing derivatives for later use
      Call DFDX(P, dF(1,1,1,1), 0, 1.0D0)
      Call DFDY(P, dF(1,1,1,2), 0, 1.0D0)
      Call DFDZ(P, dF(1,1,1,3), 0, 1.0D0)

* Faces

* x - boundary
      Do I = 1,Mx,Mx-1
       if(  ((xrank.EQ.Px-1).AND.(I.EQ.Mx)) .OR. 
     .      ((xrank.EQ.   0).AND.(I.EQ.1))  )then

        Do J = 1,My
          Do K = 1,Mz
           tmp2 = dsqrt(gamma*Vs(J,K,I)*P(J,K,I))                ! Sound speed
  
           if(xrank .EQ. Px-1)then
             L(1) = 0D0                                          ! Incoming characteristic
             L(2) = 0D0                                          ! Incoming characteristic
             L(5) = (U(J,K,I,1) + tmp2) *                        ! Outgoing characteristic 
     .              (dF(J,K,I,1) + tmp2*Q(J,K,I,4)*dux_dx(J,K,I))
             L(3) = U(J,K,I,1) * duy_dx(J,K,I)
             L(4) = U(J,K,I,1) * duz_dx(J,K,I)
           else !(xrank .EQ. 0)
             L(1) = (U(J,K,I,1) - tmp2) *                        ! Outgoing characteristic
     .              (dF(J,K,I,1) - tmp2*Q(J,K,I,4)*dux_dx(J,K,I))
             L(2) = 0D0                                          ! Incoming characteristic
             L(5) = 0D0                                          ! Incoming characteristic
             L(3) = U(J,K,I,1) * duy_dx(J,K,I)
             L(4) = U(J,K,I,1) * duz_dx(J,K,I)
           endif


           d(J,K,I,1) = (1D0/(tmp2*tmp2)) *
     .          ( L(2) + 5.0D-1*(L(1)+L(5)) )
           d(J,K,I,2) = 5.0D-1*(L(5)+L(1))
           d(J,K,I,3) = (5.0D-1/(Q(J,K,I,4)*tmp2))*(L(5)-L(1))
           d(J,K,I,4) = L(3)
           d(J,K,I,5) = L(4)

           Qrhs(J,K,I,4) = -d(J,K,I,1) + Qconvy(J,K,I,4) 
     .                   + Qconvz(J,K,I,4) 
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
     .                   - d(J,K,I,2)/(gamma-1D0) 
     .                   - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
     .                   - Q(J,K,I,3)*d(J,K,I,5)
     .                   + Qconvy(J,K,I,5) + Qconvz(J,K,I,5)
     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*d(J,K,I,1) 
     .                     - Q(J,K,I,4)*d(J,K,I,3) 
     .                   + Qconvy(J,K,I,1) + Qconvz(J,K,I,1) 
     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*d(J,K,I,1) 
     .                     - Q(J,K,I,4)*d(J,K,I,4)
     .                   + Qconvy(J,K,I,2) + Qconvz(J,K,I,2)
     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*d(J,K,I,1) 
     .                     - Q(J,K,I,4)*d(J,K,I,5)
     .                   + Qconvy(J,K,I,3) + Qconvz(J,K,I,3)
     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
         End Do
        End Do 
       endif
      End Do 

* y - boundary
      Do J = 1,My,My-1
       if(  ((yrank.EQ.Py-1).AND.(J.EQ.My)) .OR.
     .      ((yrank.EQ.   0).AND.(J.EQ.1))  )then

        Do I = 1,Mx
         Do K = 1,Mz
           tmp2 = dsqrt(gamma*Vs(J,K,I)*P(J,K,I))                ! Sound speed

           if(yrank .EQ. Py-1)then
             M(1) = 0D0                                          ! Incoming characteristic
             M(2) = 0D0                                          ! Incoming characteristic
             M(5) = (U(J,K,I,2) + tmp2) *                        ! Outgoing characteristic
     .              (dF(J,K,I,2) + tmp2*Q(J,K,I,4)*duy_dy(J,K,I))
             M(3) = U(J,K,I,2) * dux_dy(J,K,I)
             M(4) = U(J,K,I,2) * duz_dy(J,K,I) 
           else !(yrank .EQ. 0)
             M(1) = (U(J,K,I,2) - tmp2) *                        ! Outgoing characteristic
     .              (dF(J,K,I,2) - tmp2*Q(J,K,I,4)*duy_dy(J,K,I))
             M(2) = 0D0                                          ! Incoming characteristic
             M(5) = 0D0                                          ! Incoming characteristic
             M(3) = U(J,K,I,2) * dux_dy(J,K,I)
             M(4) = U(J,K,I,2) * duz_dy(J,K,I)
           endif

           e(J,K,I,1) = (1D0/(tmp2*tmp2)) *
     .          ( M(2) + 5.0D-1*(M(1)+M(5)) )
           e(J,K,I,2) = 5.0D-1*(M(5)+M(1))
           e(J,K,I,3) = (5.0D-1/(Q(J,K,I,4)*tmp2))*(M(5)-M(1))
           e(J,K,I,4) = M(3)
           e(J,K,I,5) = M(4)

           Qrhs(J,K,I,4) = -e(J,K,I,1) + Qconvx(J,K,I,4) 
     .                   + Qconvz(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
     .                   - e(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
     .                   - Q(J,K,I,3)*e(J,K,I,5)
     .                   + Qconvx(J,K,I,5) + Qconvz(J,K,I,5)
     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*e(J,K,I,1) 
     .                     - Q(J,K,I,4)*e(J,K,I,3)
     .                   + Qconvx(J,K,I,2) + Qconvz(J,K,I,2)
     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*e(J,K,I,1) 
     .                     - Q(J,K,I,4)*e(J,K,I,4)
     .                   + Qconvx(J,K,I,1) + Qconvz(J,K,I,1)
     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,5)
     .                   + Qconvx(J,K,I,3) + Qconvz(J,K,I,3)
     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
          End Do
        End Do
       endif
      End Do

* z - boundary
      Do K = 1,Mz,Mz-1
       if(  ((zrank.EQ.Pz-1).AND.(K.EQ.Mz)) .OR.
     .      ((zrank.EQ.   0).AND.(K.EQ.1))  )then

        Do I = 1,Mx
          Do J = 1,My
           tmp2 = dsqrt(gamma*Vs(J,K,I)*P(J,K,I))               ! Sound speed

           if(zrank .EQ. Pz-1)then
             N(1) = 0D0                                         ! Incoming characteristic
             N(2) = 0D0                                         ! Incoming characteristic
             N(5) = (U(J,K,I,3) + tmp2) *                       ! Outgoing characteristic
     .              (dF(J,K,I,3) + tmp2*Q(J,K,I,4)*duz_dz(J,K,I))
             N(3) = U(J,K,I,3) * dux_dz(J,K,I)
             N(4) = U(J,K,I,3) * duy_dz(J,K,I)
           else
             N(1) = (U(J,K,I,3) - tmp2) *                       ! Outgoing characteristic
     .              (dF(J,K,I,3) - tmp2*Q(J,K,I,4)*duz_dz(J,K,I))
             N(2) = 0D0                                         ! Incoming characteristic
             N(5) = 0D0                                         ! Incoming characteristic
             N(3) = U(J,K,I,3) * dux_dz(J,K,I)
             N(4) = U(J,K,I,3) * duy_dz(J,K,I)
           endif

           f(J,K,I,1) = (1D0/(tmp2*tmp2)) *
     .          ( N(2) + 5.0D-1*(N(1)+N(5)) )
           f(J,K,I,2) = 5.0D-1*(N(5)+N(1))
           f(J,K,I,3) = (5.0D-1/(Q(J,K,I,4)*tmp2))*(N(5)-N(1))
           f(J,K,I,4) = N(3)
           f(J,K,I,5) = N(4)

           Qrhs(J,K,I,4) = -f(J,K,I,1) 
     .                   + Qconvx(J,K,I,4) + Qconvy(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
     .                   - f(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
     .                   - Q(J,K,I,2)*f(J,K,I,5)
     .                   + Qconvx(J,K,I,5) + Qconvy(J,K,I,5)
     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*f(J,K,I,1) 
     .                     - Q(J,K,I,4)*f(J,K,I,3)
     .                   + Qconvx(J,K,I,3) + Qconvy(J,K,I,3)
     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*f(J,K,I,1) 
     .                     - Q(J,K,I,4)*f(J,K,I,4)
     .                   + Qconvx(J,K,I,1) + Qconvy(J,K,I,1)
     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*f(J,K,I,1) 
     .                     - Q(J,K,I,4)*f(J,K,I,5)
     .                   + Qconvx(J,K,I,2) + Qconvy(J,K,I,2)
     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
          End Do
        End Do
       endif 
      End Do

* Edges

* x-y edges
      if(((xrank.EQ.0).AND.(yrank.EQ.0)).OR.
     .   ((xrank.EQ.0).AND.(yrank.EQ.Py-1)).OR.
     .   ((xrank.EQ.Px-1).AND.(yrank.EQ.0)).OR.
     .   ((xrank.EQ.Px-1).AND.(yrank.EQ.Py-1)))then
        
        if((xrank.EQ.0).AND.(yrank.EQ.0))then
          I = 1; J = 1
        endif
        if((xrank.EQ.0).AND.(yrank.EQ.Py-1))then
          I = 1; J = My
        endif
        if((xrank.EQ.Px-1).AND.(yrank.EQ.0))then
          I = Mx; J = 1
        endif
        if((xrank.EQ.Px-1).AND.(yrank.EQ.Py-1))then
          I = Mx; J = My
        endif

        Do K = 1,Mz
           Qrhs(J,K,I,4) = -d(J,K,I,1) - e(J,K,I,1)
     .                   + Qconvz(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
     .                   - d(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
     .                   - Q(J,K,I,3)*d(J,K,I,5)
     .                     -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
     .                   - e(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
     .                   - Q(J,K,I,3)*e(J,K,I,5)
     .                   + Qconvz(J,K,I,5)
     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*d(J,K,I,1)
     .                     - Q(J,K,I,4)*d(J,K,I,3)
     .                     - U(J,K,I,1)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,4)
     .                   + Qconvz(J,K,I,1)
     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*d(J,K,I,1)
     .                     - Q(J,K,I,4)*d(J,K,I,4)
     .                     - U(J,K,I,2)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,3)
     .                   + Qconvz(J,K,I,2)
     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*d(J,K,I,1)
     .                     - Q(J,K,I,4)*d(J,K,I,5)
     .                     - U(J,K,I,3)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,5)
     .                   + Qconvz(J,K,I,3)
     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
        End Do

      endif

* y-z edges
      if(((yrank.EQ.0).AND.(zrank.EQ.0)).OR.
     .   ((yrank.EQ.0).AND.(zrank.EQ.Pz-1)).OR.
     .   ((yrank.EQ.Py-1).AND.(zrank.EQ.0)).OR.
     .   ((yrank.EQ.Py-1).AND.(zrank.EQ.Pz-1)))then

        if((yrank.EQ.0).AND.(zrank.EQ.0))then
          J = 1; K = 1
        endif
        if((yrank.EQ.0).AND.(zrank.EQ.Pz-1))then
          J = 1; K = Mz
        endif
        if((yrank.EQ.Py-1).AND.(zrank.EQ.0))then
          J = My; K = 1
        endif
        if((yrank.EQ.Py-1).AND.(zrank.EQ.Pz-1))then
          J = My; K = Mz
        endif

        Do I = 1,Mx
           Qrhs(J,K,I,4) = -e(J,K,I,1) - f(J,K,I,1)
     .                   + Qconvx(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
     .                   - e(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
     .                   - Q(J,K,I,3)*e(J,K,I,5)
     .                   -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
     .                   - f(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
     .                   - Q(J,K,I,2)*f(J,K,I,5)
     .                   + Qconvx(J,K,I,5) 
     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,3)
     .                     - U(J,K,I,2)*f(J,K,I,1)
     .                     - Q(J,K,I,4)*f(J,K,I,5)
     .                   + Qconvx(J,K,I,2)
     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,4)
     .                     - U(J,K,I,1)*f(J,K,I,1)
     .                     - Q(J,K,I,4)*f(J,K,I,4)
     .                   + Qconvx(J,K,I,1)
     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,5)
     .                     - U(J,K,I,3)*f(J,K,I,1)
     .                     - Q(J,K,I,4)*f(J,K,I,3)
     .                   + Qconvx(J,K,I,3)
     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
        End Do

      endif

* z-x edges
      if(((zrank.EQ.0).AND.(xrank.EQ.0)).OR.
     .   ((zrank.EQ.0).AND.(xrank.EQ.Px-1)).OR.
     .   ((zrank.EQ.Pz-1).AND.(xrank.EQ.0)).OR.
     .   ((zrank.EQ.Pz-1).AND.(xrank.EQ.Px-1)))then

        if((zrank.EQ.0).AND.(xrank.EQ.0))then
          K = 1; I = 1
        endif
        if((zrank.EQ.0).AND.(xrank.EQ.Px-1))then
          K = 1; I = Mx
        endif
        if((zrank.EQ.Pz-1).AND.(xrank.EQ.0))then
          K = Mz; I = 1
        endif
        if((zrank.EQ.Pz-1).AND.(xrank.EQ.Px-1))then
          K = Mz; I = Mx
        endif

        Do J = 1,My
           Qrhs(J,K,I,4) = -f(J,K,I,1) - d(J,K,I,4)
     .                   + Qconvy(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
     .                   - f(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
     .                   - Q(J,K,I,2)*f(J,K,I,5)
     .                   -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
     .                   - d(J,K,I,2)/(gamma-1D0)
     .                   - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
     .                   - Q(J,K,I,3)*d(J,K,I,5)
     .                   + Qconvy(J,K,I,5)
     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*f(J,K,I,1)
     .                     - Q(J,K,I,4)*f(J,K,I,3)
     .                     - U(J,K,I,3)*d(J,K,I,1)
     .                     - Q(J,K,I,4)*d(J,K,I,5)
     .                   + Qconvy(J,K,I,3)
     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*f(J,K,I,1)
     .                     - Q(J,K,I,4)*f(J,K,I,4)
     .                     - U(J,K,I,1)*d(J,K,I,1)
     .                     - Q(J,K,I,4)*d(J,K,I,3)
     .                   + Qconvy(J,K,I,1)
     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*f(J,K,I,1)
     .                     - Q(J,K,I,4)*f(J,K,I,5)
     .                     - U(J,K,I,2)*d(J,K,I,1)
     .                     - Q(J,K,I,4)*d(J,K,I,4)
     .                   + Qconvy(J,K,I,2)
     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
        End Do

      endif

* Corners
      if(((xrank.EQ.0).OR.(xrank.EQ.Px-1)).AND.
     .   ((yrank.EQ.0).OR.(yrank.EQ.Py-1)).AND.
     .   ((zrank.EQ.0).OR.(zrank.EQ.Pz-1)))then
    
        if(xrank.EQ.0) then
          I = 1
        else 
          I = Mx
        endif
        if(yrank.EQ.0) then
          J = 1
        else
          J = My
        endif
        if(zrank.EQ.0) then
          K = 1
        else
          K = Mz
        endif

        Qrhs(J,K,I,4) = -d(J,K,I,1) - e(J,K,I,1) - f(J,K,I,1)
        Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
     .                - d(J,K,I,2)/(gamma-1D0)
     .                - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
     .                - Q(J,K,I,3)*d(J,K,I,5)
     .                  -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
     .                - e(J,K,I,2)/(gamma-1D0)
     .                - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
     .                - Q(J,K,I,3)*e(J,K,I,5)
     .                  -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
     .                - f(J,K,I,2)/(gamma-1D0)
     .                - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
     .                - Q(J,K,I,2)*f(J,K,I,5)
     .                + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
        Qrhs(J,K,I,1) = - U(J,K,I,1)*d(J,K,I,1)
     .                  - Q(J,K,I,4)*d(J,K,I,3)
     .                  - U(J,K,I,1)*e(J,K,I,1)
     .                  - Q(J,K,I,4)*e(J,K,I,4)
     .                  - U(J,K,I,1)*f(J,K,I,1)
     .                  - Q(J,K,I,4)*f(J,K,I,4)
     .                + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
        Qrhs(J,K,I,2) = - U(J,K,I,2)*d(J,K,I,1)
     .                  - Q(J,K,I,4)*d(J,K,I,4)
     .                  - U(J,K,I,2)*e(J,K,I,1)
     .                  - Q(J,K,I,4)*e(J,K,I,3)
     .                  - U(J,K,I,2)*f(J,K,I,1)
     .                  - Q(J,K,I,4)*f(J,K,I,5)
     .                + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
        Qrhs(J,K,I,3) = - U(J,K,I,3)*d(J,K,I,1)
     .                  - Q(J,K,I,4)*d(J,K,I,5)
     .                  - U(J,K,I,3)*e(J,K,I,1)
     .                  - Q(J,K,I,4)*e(J,K,I,5)
     .                  - U(J,K,I,3)*f(J,K,I,1)
     .                  - Q(J,K,I,4)*f(J,K,I,3)
     .                + Qvisc(J,K,I,3)  + Qh(J,K,I,3)

      endif

      End Subroutine bc_RHS
