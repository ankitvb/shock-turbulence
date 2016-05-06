************************************************************
* Routines pertaining to hyperviscosity                    *
************************************************************
    
************************************************************
* Subroutine Compute_Contraction                           *
*                                                          *
* Computes the contraction of a 2nd order tensor           *
* according to the formula Sij*Sji                         *
************************************************************

      Subroutine Compute_Contraction(F, Fc)
 
      Include 'header'

      Real*8 F(My,Mz,Mx,9), Fc(My,Mz,Mx)
       
      Integer I,J,K

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Fc(J,K,I) = dsqrt(5.0D-1*(F(J,K,I,1)*F(J,K,I,1) 
     .                + F(J,K,I,5)*F(J,K,I,5) + F(J,K,I,9)*F(J,K,I,9)
     .                + 2*(F(J,K,I,2)*F(J,K,I,4) + F(J,K,I,3)*F(J,K,I,7)
     .                +    F(J,K,I,6)*F(J,K,I,8))))
          End Do
        End Do
      End Do

      Return

      End 

************************************************************
* Subroutine Compute_Laplacian                             *
*                                                          *
* Computes the laplacian of a scalar field                 *
************************************************************

      Subroutine Compute_Laplacian(S, d2S)

      Include 'header'

      Real*8 S(My,Mz,Mx), dS(My,Mz,Mx), d2S(My,Mz,Mx)
      Real*8 temp(My,Mz,Mx,2)      
  
      Integer I,J,K

      Call DFDX_sym(S, dS, 0, 1D0)
      Call DFDX_asym(dS, temp(1,1,1,1), 0, 1D0)
      Call DFDY_sym(S, dS, 0, 1D0)
      Call DFDY_asym(dS, temp(1,1,1,2), 0, 1D0)
      Call DFDZ_sym(S, dS, 0, 1D0)
      Call DFDZ_asym(dS, d2S, 0, 1D0)

c      Call D2FDX2_sym(S, dS, temp(1,1,1,1), 0)
c      Call D2FDY2_sym(S, dS, temp(1,1,1,2), 0)
c      Call D2FDZ2_sym(S, dS, d2S, 0)

      d2S = d2S + temp(:,:,:,1) + temp(:,:,:,2)

      Return
    
      End

***********************************************************
* Subroutine Compute_hyper                                *
*                                                         *
* Computes the nonlinear hyperviscosity based on          *
* Cook and Cabot model (JCP, 2004)                        *
***********************************************************
      Subroutine Compute_hyper(Q, Y)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), U(My,Mz,Mx,3)
      Real*8 Y(My,Mz,Mx,1), gamma_eff(My,Mz,Mx)
      Common /THERMO/ gamma_eff

      Real*8 P(My,Mz,Mx),   Vs(My,Mz,Mx)
      Real*8 T(My,Mz,Mx),   e(My,Mz,Mx)
      Real*8 Sij(My,Mz,Mx,9), S(My,Mz,Mx), d2S(My,Mz,Mx)
      Real*8 d4S(My,Mz,Mx), d4Sg(My,Mz,Mx)
      Real*8 d4e(My,Mz,Mx), d4eg(My,Mz,Mx)
      Real*8 d2Y(My,Mz,Mx), d4Y(My,Mz,Mx)
      Real*8 eta4(My,Mz,Mx,4), eta4g(My,Mz,Mx,4)
      Real*8 beta_hyper(My,Mz,Mx)
      Real*8 visc_hyper(My,Mz,Mx), k_hyper(My,Mz,Mx)
      Real*8 D_hyper(My,Mz,Mx)
      Real*8 C_mu4, C_beta4, C_k4, C_D, C_Y

* Locals
      Real*8  fac1, fac2, fac3, cs, tmp1, tmp2
      Integer I, J, K

      Common /HYPER/ beta_hyper, visc_hyper, k_hyper, D_hyper

      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)
      Real*8 dTdx(My,Mz,Mx),   dTdy(My,Mz,Mx),   dTdz(My,Mz,Mx)

      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz,
     .              dTdx, dTdy, dTdz

* Coefficients based on Cook and Cabot (JCP, 2004)
      C_mu4   = 2.0D-3
      C_beta4 = 1.0D0
      C_k4    = 1.0D-2
      C_D     = 1.0D-2
      C_Y     = 100D0

* Setting up the strain matrix
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            U(J,K,I,1) = Q(J,K,I,1) / Q(J,K,I,4)
            U(J,K,I,2) = Q(J,K,I,2) / Q(J,K,I,4)
            U(J,K,I,3) = Q(J,K,I,3) / Q(J,K,I,4)
          End Do
        End Do
      End Do

       Sij(:,:,:,1) = dux_dx
       Sij(:,:,:,4) = duy_dx
       Sij(:,:,:,7) = duz_dx

       Sij(:,:,:,2) = dux_dy
       Sij(:,:,:,5) = duy_dy
       Sij(:,:,:,8) = duz_dy

       Sij(:,:,:,3) = dux_dz
       Sij(:,:,:,6) = duy_dz
       Sij(:,:,:,9) = duz_dz

* Making it symmetric
      Sij(:,:,:,2) = 5.0D-1*(Sij(:,:,:,2) + Sij(:,:,:,4))
      Sij(:,:,:,3) = 5.0D-1*(Sij(:,:,:,3) + Sij(:,:,:,7))
      Sij(:,:,:,6) = 5.0D-1*(Sij(:,:,:,6) + Sij(:,:,:,8))
      Sij(:,:,:,4) = Sij(:,:,:,2)
      Sij(:,:,:,7) = Sij(:,:,:,3)
      Sij(:,:,:,8) = Sij(:,:,:,6)

* Computing P, T and e required for hyperconductivity
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             fac1 = gamma_eff(J,K,I) - 1.0D0
             fac2 = gamma_eff(J,K,I) / fac1
             fac3 = gamma1 / (gamma1 - 1D0)
             fac3 = fac3 *
     .       1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

             Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
              P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                         ( Q(J,K,I,1) * Q(J,K,I,1)
     .                         + Q(J,K,I,2) * Q(J,K,I,2)
     .                         + Q(J,K,I,3) * Q(J,K,I,3) ) )
              T(J,K,I) = fac3 * Vs(J,K,I) * P(J,K,I)
              e(J,K,I) = Vs(J,K,I)* ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
          End Do
        End Do
      End Do

* Computing the magnitude of the strain tensor
      Call Compute_Contraction(Sij, S)

* Applying biharmonic operator to the strain magnitude and energy 
      Call Compute_Laplacian(S,   d2S)
      Call Compute_Laplacian(d2S, d4S)

      Call Compute_Laplacian(e, d2S)
      Call Compute_Laplacian(d2S, d4e)

      Call Compute_Laplacian(Y(1,1,1,1), d2Y)
      Call Compute_Laplacian(d2Y, d4Y)

      d4S = abs(d4S)
      d4e = abs(d4e)
      d4Y = abs(d4Y)

* Filtering the field with Gaussian filter (Cook & Cabot, JCP 2004)
c      Call Gaussian_filter(d4S, d4Sg, 1)
c      Call Gaussian_filter(d4e, d4eg, 1)

* Calculating the hyperviscosities
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            cs          = dsqrt(gamma_eff(J,K,I) * P(J,K,I) * Vs(J,K,I))
            fac3          = (Q(J,K,I,4)/T(J,K,I)) * cs
            eta4(J,K,I,1) = Q(J,K,I,4) * d4S(J,K,I)
            eta4(J,K,I,2) =       fac3 * d4e(J,K,I)
            eta4(J,K,I,3) =         cs * d4Y(J,K,I)

            if(Y(J,K,I,1) .LT. 1D0) then
               tmp1 = 0D0
            else
               tmp1 = 1D0
            endif
            if(Y(J,K,I,1) .LT. 0D0) then
               tmp2 = 0D0
            else
               tmp2 = 1D0
            endif
            eta4(J,K,I,4) = cs *
     .                    ( ( Y(J,K,I,1) - 1D0 ) * tmp1
     .                    - Y(J,K,I,1) * (1D0 - tmp2) )
          End Do
        End Do
      End Do

* Filtering the field with Gaussian filter (Cook & Cabot, JCP 2004)
      Call Gaussian_filter(eta4(1,1,1,1), eta4g(1,1,1,1), 4)

* Calculating the hyperviscosities
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            beta_hyper(J,K,I) = C_beta4 * eta4g(J,K,I,1) * (DX**6D0)
            visc_hyper(J,K,I) = C_mu4   * eta4g(J,K,I,1) * (DX**6D0)
               k_hyper(J,K,I) = C_k4    * eta4g(J,K,I,2) * (DX**5D0)
               D_hyper(J,K,I) = C_D     * eta4g(J,K,I,3) * (DX**5D0)
     .                        + C_Y     * eta4g(J,K,I,4) *  DX  
          End Do
        End Do
      End Do

      Return

      End



