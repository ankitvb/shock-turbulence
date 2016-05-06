************************************************************************
* Subroutine: Set_BC                                                   *
*                                                                      *
* Function: Set non reflecting boundary conditions for spherical       *
*           shock-turbulence interaction. Subsonic boundary,           *
*           outgoing waves only, incoming wave amplitudes set to zero  *
*           Sponge damping applied to relax the flow to steady ref     *                                                          
*           state                                                      *
************************************************************************

      Subroutine bc_RHS(Q, Qrhs, Y, Yrhs, time)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8 Y(My,Mz,Mx,1), Yrhs(My,Mz,Mx,1)
      Real*8 Yconvx(My,Mz,Mx,1), Yconvy(My,Mz,Mx,1)
      Real*8 Yconvz(My,Mz,Mx,1), Ydiff(My,Mz,Mx,1)

      Common /SPECIES_RHSV/ Yconvx, Yconvy, Yconvz, Ydiff

      Real*8 gamma_eff(My,Mz,Mx)
      Real*8 P(My,Mz,Mx), U(My,Mz,Mx,3), Ur(My,Mz,Mx)
      Real*8 dF(My,Mz,Mx,3)
      Real*8 Vs(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 L(6), M(6), N(6)
      Real*8 time
 
      Common /THERMO/ gamma_eff

      Real*8  dYdx(My,Mz,Mx), dYdy(My,Mz,Mx), dYdz(My,Mz,Mx)

      Common /SPECIES_DERS/ dYdx, dYdy, dYdz

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
      Real*8 d(My,Mz,Mx,6), e(My,Mz,Mx,6), f(My,Mz,Mx,6)

      Integer I, J, K

* Make velocities normal to inner faces zero
*
* Plane x = 0 ; u_x = 0
* Plane y = 0 ; u_y = 0
* Plane z = 0 ; u_z = 0
*

      Do J = 1,My
       Do K = 1,Mz
        if(xrank .EQ. 0)then
          Q(J,K,1,1) = 0D0
        endif
       End Do
      End Do

      Do I = 1,Mx
       Do J = 1,My
        if(zrank .EQ. 0)then
         Q(J,1,I,3) = 0D0
        endif
       End Do
      End Do

      Do I = 1,Mx
       Do K = 1,Mz
        if(yrank .EQ. 0)then
         Q(1,K,I,2) = 0D0
        endif
       End Do
      End Do

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

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac1 = gamma_eff(J,K,I) - 1.0D0
            fac2 = gamma_eff(J,K,I) / fac1
            fac3 = gamma1 / (gamma1 - 1D0)
            fac3 = fac3 *
     .      1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

            Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
             P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) ) )
             T(J,K,I) = fac3 * Vs(J,K,I)  * P(J,K,I)

             U(J,K,I,1) = Q(J,K,I,1) * Vs(J,K,I)
             U(J,K,I,2) = Q(J,K,I,2) * Vs(J,K,I)
             U(J,K,I,3) = Q(J,K,I,3) * Vs(J,K,I)

             Ur(J,K,I) = dsqrt(U(J,K,I,1)*U(J,K,I,1)
     .                       + U(J,K,I,2)*U(J,K,I,2)
     .                       + U(J,K,I,3)*U(J,K,I,3))
          End Do
        End Do
      End Do

* Precomputing derivatives for later use
      Call DFDX_sym(P, dF(1,1,1,1), 0, 1.0D0)
      Call DFDY_sym(P, dF(1,1,1,2), 0, 1.0D0)
      Call DFDZ_sym(P, dF(1,1,1,3), 0, 1.0D0)

* Faces

* x - boundary
       if(xrank.EQ.Px-1) then 
        I = Mx
        Do J = 1,My
          Do K = 1,Mz
           tmp2 = dsqrt(gamma_eff(J,K,I)*Vs(J,K,I)*P(J,K,I))   ! Sound speed
  
           L(1) = 0D0                                          ! Incoming characteristic
           L(2) = 0D0                                          ! Incoming characteristic
           L(5) = (U(J,K,I,1) + tmp2) *                        ! Outgoing characteristic 
     .            (dF(J,K,I,1) + tmp2*Q(J,K,I,4)*dux_dx(J,K,I))
           L(3) = 0D0 !U(J,K,I,1) * duy_dx(J,K,I)                   ! Incoming characteristic
           L(4) = 0D0 !U(J,K,I,1) * duz_dx(J,K,I)                   ! Incoming characteristic
           L(6) = 0D0 !U(J,K,I,1) * dYdx(J,K,I)                     ! Incoming characteristic

           d(J,K,I,1) = (1D0/(tmp2*tmp2)) *
     .          ( L(2) + 5.0D-1*(L(1)+L(5)) )
           d(J,K,I,2) = 5.0D-1*(L(5)+L(1))
           d(J,K,I,3) = (5.0D-1/(Q(J,K,I,4)*tmp2))*(L(5)-L(1))
           d(J,K,I,4) = L(3)
           d(J,K,I,5) = L(4)
           d(J,K,I,6) = L(6)

           Qrhs(J,K,I,4) = -d(J,K,I,1) !+ Qconvy(J,K,I,4) 
c     .                   + Qconvz(J,K,I,4) 
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
     .                   - d(J,K,I,2)/(gamma_eff(J,K,I)-1D0) 
     .                   - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
     .                   - Q(J,K,I,3)*d(J,K,I,5)
c     .                   + Qconvy(J,K,I,5) + Qconvz(J,K,I,5)
c     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*d(J,K,I,1) 
     .                     - Q(J,K,I,4)*d(J,K,I,3) 
c     .                   + Qconvy(J,K,I,1) + Qconvz(J,K,I,1) 
c     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*d(J,K,I,1) 
     .                     - Q(J,K,I,4)*d(J,K,I,4)
c     .                   + Qconvy(J,K,I,2) + Qconvz(J,K,I,2)
c     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*d(J,K,I,1) 
     .                     - Q(J,K,I,4)*d(J,K,I,5)
c     .                   + Qconvy(J,K,I,3) + Qconvz(J,K,I,3)
c     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
           Yrhs(J,K,I,1) = - Q(J,K,I,4) * d(J,K,I,6) 
c     .                   + Yconvy(J,K,I,1) + Yconvz(J,K,I,1)
c     .                   + Ydiff(J,K,I,1)
         End Do
        End Do 
       endif

* y - boundary
       if(yrank.EQ.Py-1) then
        J = My 
        Do I = 1,Mx
         Do K = 1,Mz
           tmp2 = dsqrt(gamma_eff(J,K,I)*Vs(J,K,I)*P(J,K,I))   ! Sound speed

           M(1) = 0D0                                          ! Incoming characteristic
           M(2) = 0D0                                          ! Incoming characteristic
           M(5) = (U(J,K,I,2) + tmp2) *                        ! Outgoing characteristic
     .            (dF(J,K,I,2) + tmp2*Q(J,K,I,4)*duy_dy(J,K,I))
           M(3) = 0D0 !U(J,K,I,2) * dux_dy(J,K,I)                   ! Incoming characteristic
           M(4) = 0D0 !U(J,K,I,2) * duz_dy(J,K,I)                   ! Incoming characteristic
           M(6) = 0D0 !U(J,K,I,2) * dYdy(J,K,I)                     ! Incoming characteristic

           e(J,K,I,1) = (1D0/(tmp2*tmp2)) *
     .          ( M(2) + 5.0D-1*(M(1)+M(5)) )
           e(J,K,I,2) = 5.0D-1*(M(5)+M(1))
           e(J,K,I,3) = (5.0D-1/(Q(J,K,I,4)*tmp2))*(M(5)-M(1))
           e(J,K,I,4) = M(3)
           e(J,K,I,5) = M(4)
           e(J,K,I,6) = M(6)

           Qrhs(J,K,I,4) = -e(J,K,I,1) !+ Qconvx(J,K,I,4) 
c     .                   + Qconvz(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
     .                   - e(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
     .                   - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
     .                   - Q(J,K,I,3)*e(J,K,I,5)
c     .                   + Qconvx(J,K,I,5) + Qconvz(J,K,I,5)
c     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*e(J,K,I,1) 
     .                     - Q(J,K,I,4)*e(J,K,I,3)
c     .                   + Qconvx(J,K,I,2) + Qconvz(J,K,I,2)
c     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*e(J,K,I,1) 
     .                     - Q(J,K,I,4)*e(J,K,I,4)
c     .                   + Qconvx(J,K,I,1) + Qconvz(J,K,I,1)
c     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*e(J,K,I,1)
     .                     - Q(J,K,I,4)*e(J,K,I,5)
c     .                   + Qconvx(J,K,I,3) + Qconvz(J,K,I,3)
c     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
           Yrhs(J,K,I,1) = - Q(J,K,I,4) * e(J,K,I,6)
c     .                   + Yconvx(J,K,I,1) + Yconvz(J,K,I,1)
c     .                   + Ydiff(J,K,I,1)
          End Do
        End Do
       endif

* z - boundary
       if(zrank.EQ.Pz-1) then
        K = Mz 
        Do I = 1,Mx
          Do J = 1,My
           tmp2 = dsqrt(gamma_eff(J,K,I)*Vs(J,K,I)*P(J,K,I))   ! Sound speed
  
           N(1) = 0D0                                          ! Incoming characteristic
           N(2) = 0D0                                          ! Incoming characteristic
           N(5) = (U(J,K,I,3) + tmp2) *                        ! Outgoing characteristic
     .            (dF(J,K,I,3) + tmp2*Q(J,K,I,4)*duz_dz(J,K,I))
           N(3) = 0D0 !U(J,K,I,3) * dux_dz(J,K,I)                   ! Incoming characteristic
           N(4) = 0D0 !U(J,K,I,3) * duy_dz(J,K,I)                   ! Incoming characteristic
           N(6) = 0D0 !U(J,K,I,3) * dYdz(J,K,I)                     ! Incoming characteristic

           f(J,K,I,1) = (1D0/(tmp2*tmp2)) *
     .          ( N(2) + 5.0D-1*(N(1)+N(5)) )
           f(J,K,I,2) = 5.0D-1*(N(5)+N(1))
           f(J,K,I,3) = (5.0D-1/(Q(J,K,I,4)*tmp2))*(N(5)-N(1))
           f(J,K,I,4) = N(3)
           f(J,K,I,5) = N(4)
           f(J,K,I,6) = N(6)

           Qrhs(J,K,I,4) = -f(J,K,I,1) 
c     .                   + Qconvx(J,K,I,4) + Qconvy(J,K,I,4)
           Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
     .                   - f(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
     .                   - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
     .                   - Q(J,K,I,2)*f(J,K,I,5)
c     .                   + Qconvx(J,K,I,5) + Qconvy(J,K,I,5)
c     .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
           Qrhs(J,K,I,3) = - U(J,K,I,3)*f(J,K,I,1) 
     .                     - Q(J,K,I,4)*f(J,K,I,3)
c     .                   + Qconvx(J,K,I,3) + Qconvy(J,K,I,3)
c     .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
           Qrhs(J,K,I,1) = - U(J,K,I,1)*f(J,K,I,1) 
     .                     - Q(J,K,I,4)*f(J,K,I,4)
c     .                   + Qconvx(J,K,I,1) + Qconvy(J,K,I,1)
c     .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
           Qrhs(J,K,I,2) = - U(J,K,I,2)*f(J,K,I,1) 
     .                     - Q(J,K,I,4)*f(J,K,I,5)
c     .                   + Qconvx(J,K,I,2) + Qconvy(J,K,I,2)
c     .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
           Yrhs(J,K,I,1) = - Q(J,K,I,4) * f(J,K,I,6)
c     .                   + Yconvx(J,K,I,1) + Yconvy(J,K,I,1)
c     .                   + Ydiff(J,K,I,1)
          End Do
        End Do
       endif 

* Edges

* x-y edges
!       if((xrank.EQ.Px-1).AND.(yrank.EQ.Py-1))then
!         I = Mx; J = My
!          Do K = 1,Mz
!            Qrhs(J,K,I,4) = -d(J,K,I,1) - e(J,K,I,1)
!      .                   + Qconvz(J,K,I,4)
!            Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
!      .                   - d(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                   - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
!      .                   - Q(J,K,I,3)*d(J,K,I,5)
!      .                     -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
!      .                   - e(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                   - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
!      .                   - Q(J,K,I,3)*e(J,K,I,5)
!      .                   + Qconvz(J,K,I,5)
!      .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
!            Qrhs(J,K,I,1) = - U(J,K,I,1)*d(J,K,I,1)
!      .                     - Q(J,K,I,4)*d(J,K,I,3)
!      .                     - U(J,K,I,1)*e(J,K,I,1)
!      .                     - Q(J,K,I,4)*e(J,K,I,4)
!      .                   + Qconvz(J,K,I,1)
!      .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
!            Qrhs(J,K,I,2) = - U(J,K,I,2)*d(J,K,I,1)
!      .                     - Q(J,K,I,4)*d(J,K,I,4)
!      .                     - U(J,K,I,2)*e(J,K,I,1)
!      .                     - Q(J,K,I,4)*e(J,K,I,3)
!      .                   + Qconvz(J,K,I,2)
!      .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
!            Qrhs(J,K,I,3) = - U(J,K,I,3)*d(J,K,I,1)
!      .                     - Q(J,K,I,4)*d(J,K,I,5)
!      .                     - U(J,K,I,3)*e(J,K,I,1)
!      .                     - Q(J,K,I,4)*e(J,K,I,5)
!      .                   + Qconvz(J,K,I,3)
!      .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
!            Yrhs(J,K,I,1) = - Q(J,K,I,4) * d(J,K,I,6)
!      .                     - Q(J,K,I,4) * e(J,K,I,6)
!      .                   + Yconvz(J,K,I,1) + Ydiff(J,K,I,1)
!         End Do
!       endif
! 
! * y-z edges
!       if((yrank.EQ.Py-1).AND.(zrank.EQ.Pz-1))then
!         J = My; K = Mz
!         Do I = 1,Mx
!            Qrhs(J,K,I,4) = -e(J,K,I,1) - f(J,K,I,1)
!      .                   + Qconvx(J,K,I,4)
!            Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
!      .                   - e(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                   - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
!      .                   - Q(J,K,I,3)*e(J,K,I,5)
!      .                   -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
!      .                   - f(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                   - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
!      .                   - Q(J,K,I,2)*f(J,K,I,5)
!      .                   + Qconvx(J,K,I,5) 
!      .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
!            Qrhs(J,K,I,2) = - U(J,K,I,2)*e(J,K,I,1)
!      .                     - Q(J,K,I,4)*e(J,K,I,3)
!      .                     - U(J,K,I,2)*f(J,K,I,1)
!      .                     - Q(J,K,I,4)*f(J,K,I,5)
!      .                   + Qconvx(J,K,I,2)
!      .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
!            Qrhs(J,K,I,1) = - U(J,K,I,1)*e(J,K,I,1)
!      .                     - Q(J,K,I,4)*e(J,K,I,4)
!      .                     - U(J,K,I,1)*f(J,K,I,1)
!      .                     - Q(J,K,I,4)*f(J,K,I,4)
!      .                   + Qconvx(J,K,I,1)
!      .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
!            Qrhs(J,K,I,3) = - U(J,K,I,3)*e(J,K,I,1)
!      .                     - Q(J,K,I,4)*e(J,K,I,5)
!      .                     - U(J,K,I,3)*f(J,K,I,1)
!      .                     - Q(J,K,I,4)*f(J,K,I,3)
!      .                   + Qconvx(J,K,I,3)
!      .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
!            Yrhs(J,K,I,1) = - Q(J,K,I,4) * e(J,K,I,6)
!      .                     - Q(J,K,I,4) * f(J,K,I,6)
!      .                     + Yconvx(J,K,I,1) + Ydiff(J,K,I,1)
!         End Do
!       endif
! 
! * z-x edges
!       if((zrank.EQ.Pz-1).AND.(xrank.EQ.Px-1))then
!         K = Mz; I = Mx
!         Do J = 1,My
!            Qrhs(J,K,I,4) = -f(J,K,I,1) - d(J,K,I,4)
!      .                   + Qconvy(J,K,I,4)
!            Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
!      .                   - f(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                   - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
!      .                   - Q(J,K,I,2)*f(J,K,I,5)
!      .                   -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
!      .                   - d(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                   - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
!      .                   - Q(J,K,I,3)*d(J,K,I,5)
!      .                   + Qconvy(J,K,I,5)
!      .                   + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
!            Qrhs(J,K,I,3) = - U(J,K,I,3)*f(J,K,I,1)
!      .                     - Q(J,K,I,4)*f(J,K,I,3)
!      .                     - U(J,K,I,3)*d(J,K,I,1)
!      .                     - Q(J,K,I,4)*d(J,K,I,5)
!      .                   + Qconvy(J,K,I,3)
!      .                   + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
!            Qrhs(J,K,I,1) = - U(J,K,I,1)*f(J,K,I,1)
!      .                     - Q(J,K,I,4)*f(J,K,I,4)
!      .                     - U(J,K,I,1)*d(J,K,I,1)
!      .                     - Q(J,K,I,4)*d(J,K,I,3)
!      .                   + Qconvy(J,K,I,1)
!      .                   + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
!            Qrhs(J,K,I,2) = - U(J,K,I,2)*f(J,K,I,1)
!      .                     - Q(J,K,I,4)*f(J,K,I,5)
!      .                     - U(J,K,I,2)*d(J,K,I,1)
!      .                     - Q(J,K,I,4)*d(J,K,I,4)
!      .                   + Qconvy(J,K,I,2)
!      .                   + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
!            Yrhs(J,K,I,1) = - Q(J,K,I,4) * f(J,K,I,6)
!      .                     - Q(J,K,I,4) * d(J,K,I,6)
!      .                     + Yconvy(J,K,I,1) + Ydiff(J,K,I,1)
!         End Do
!       endif
! 
! * Corners
!       if( (xrank.EQ.Px-1).AND.
!      .    (yrank.EQ.Py-1).AND.
!      .    (zrank.EQ.Pz-1) )then
! 
!         I = Mx
!         J = My
!         K = Mz
! 
!         Qrhs(J,K,I,4) = -d(J,K,I,1) - e(J,K,I,1) - f(J,K,I,1)
!         Qrhs(J,K,I,5) = -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*d(J,K,I,1)
!      .                - d(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                - Q(J,K,I,1)*d(J,K,I,3) - Q(J,K,I,2)*d(J,K,I,4)
!      .                - Q(J,K,I,3)*d(J,K,I,5)
!      .                  -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*e(J,K,I,1)
!      .                - e(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                - Q(J,K,I,2)*e(J,K,I,3) - Q(J,K,I,1)*e(J,K,I,4)
!      .                - Q(J,K,I,3)*e(J,K,I,5)
!      .                  -5.0D-1*Ur(J,K,I)*Ur(J,K,I)*f(J,K,I,1)
!      .                - f(J,K,I,2)/(gamma_eff(J,K,I)-1D0)
!      .                - Q(J,K,I,3)*f(J,K,I,3) - Q(J,K,I,1)*f(J,K,I,4)
!      .                - Q(J,K,I,2)*f(J,K,I,5)
!      .                + Qvisc(J,K,I,5)  + Qh(J,K,I,5)
!         Qrhs(J,K,I,1) = - U(J,K,I,1)*d(J,K,I,1)
!      .                  - Q(J,K,I,4)*d(J,K,I,3)
!      .                  - U(J,K,I,1)*e(J,K,I,1)
!      .                  - Q(J,K,I,4)*e(J,K,I,4)
!      .                  - U(J,K,I,1)*f(J,K,I,1)
!      .                  - Q(J,K,I,4)*f(J,K,I,4)
!      .                + Qvisc(J,K,I,1)  + Qh(J,K,I,1)
!         Qrhs(J,K,I,2) = - U(J,K,I,2)*d(J,K,I,1)
!      .                  - Q(J,K,I,4)*d(J,K,I,4)
!      .                  - U(J,K,I,2)*e(J,K,I,1)
!      .                  - Q(J,K,I,4)*e(J,K,I,3)
!      .                  - U(J,K,I,2)*f(J,K,I,1)
!      .                  - Q(J,K,I,4)*f(J,K,I,5)
!      .                + Qvisc(J,K,I,2)  + Qh(J,K,I,2)
!         Qrhs(J,K,I,3) = - U(J,K,I,3)*d(J,K,I,1)
!      .                  - Q(J,K,I,4)*d(J,K,I,5)
!      .                  - U(J,K,I,3)*e(J,K,I,1)
!      .                  - Q(J,K,I,4)*e(J,K,I,5)
!      .                  - U(J,K,I,3)*f(J,K,I,1)
!      .                  - Q(J,K,I,4)*f(J,K,I,3)
!      .                + Qvisc(J,K,I,3)  + Qh(J,K,I,3)
!         Yrhs(J,K,I,1) = - Q(J,K,I,4) * d(J,K,I,6)
!      .                  - Q(J,K,I,4) * e(J,K,I,6)
!      .                  - Q(J,K,I,4) * f(J,K,I,6) 
!      .                  + Ydiff(J,K,I,1)
!       endif

      Return

      End Subroutine bc_RHS
