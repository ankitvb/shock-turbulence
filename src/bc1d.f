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
      Real*8 L(5)
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
      Real*8 d(My,Mz,Mx,5)

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

          End Do
        End Do
      End Do

* Precomputing derivatives for later use
      Call DFDX(P, dF(1,1,1,1), 0, 1.0D0)
      Call DFDX(Q(1,1,1,4), dF(1,1,1,2), 0, 1.0D0)

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
             L(2) = U(J,K,I,1) *                                 ! Outgoing characteristic
     .              (tmp2*tmp2*dF(J,K,I,2)-dF(J,K,I,1))
             L(5) = (U(J,K,I,1) + tmp2) *                        ! Outgoing characteristic 
     .              (dF(J,K,I,1) + tmp2*Q(J,K,I,4)*dux_dx(J,K,I))
             L(3) = U(J,K,I,1) * duy_dx(J,K,I)
             L(4) = U(J,K,I,1) * duz_dx(J,K,I)
           else !(xrank .EQ. 0)
             L(1) = (U(J,K,I,1) - tmp2) *                        ! Outgoing characteristic
     .              (dF(J,K,I,1) - tmp2*Q(J,K,I,4)*dux_dx(J,K,I))
             L(2) = 0D0                                          ! Incoming characteristic
             L(5) = 0D0                                          ! Incoming characteristic
             L(3) = 0D0 !U(J,K,I,1) * duy_dx(J,K,I)
             L(4) = 0D0 !U(J,K,I,1) * duz_dx(J,K,I)
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

      End Subroutine bc_RHS
