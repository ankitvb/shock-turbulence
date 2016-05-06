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

      Call D2FDX2(S, dS, temp(1,1,1,1), 0)
      Call D2FDY2(S, dS, temp(1,1,1,2), 0)
      Call D2FDZ2(S, dS, d2S, 0)

      d2S = d2S + temp(:,:,:,1) + temp(:,:,:,2)

      Return
    
      End

***********************************************************
* Subroutine Compute_coeff1                               *
*                                                         *
* Computes the dynamic hyperviscosity coefficients        *
* based on secondary dilatation based sensor              *
*                                                         *
*                                                         *
*                      dil^2                              *
*  C_beta =    ---------------------                      *
*              dil^2 + omega^2 + eps                      *
*                                                         *
***********************************************************
      Subroutine Compute_coeff1(Q, C_mu, C_beta, C_k)

      Include 'header'

      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4), temp(My,Mz,Mx)
      Real*8  U(My,Mz,Mx,3)

      Real*8 C_mu(My,Mz,Mx), C_beta(My,Mz,Mx)
      Real*8 C_k(My,Mz,Mx)
      Real*8 omega2(My,Mz,Mx)
      Real*8 eps, tmp, fac, fac2

      Integer I, J, K

*  Convert momenta to velocities
      Do I = 1,Mx
        Do J = 1,My
           Do K = 1,Mz
            tmp = 1.0D0 / Q(J,K,I,4)

            U(J,K,I,1) = Q(J,K,I,1) * tmp
            U(J,K,I,2) = Q(J,K,I,2) * tmp
            U(J,K,I,3) = Q(J,K,I,3) * tmp
           End Do
        End Do
      End Do

*  x-vorticity
*
*        d u_z   d u_y
*  w_x = ----- - -----
*         d y     d z

      Call DFDY(U(1,1,1,3), Grad(1,1,1,1), 0,  1.0D0)
      Call DFDZ(U(1,1,1,2), Grad(1,1,1,1), 1, -1.0D0)

*  y-vorticity
*
*        d u_x   d u_z
*  w_y = ----- - -----
*         d z     d x

      Call DFDX(U(1,1,1,3), Grad(1,1,1,2), 0, -1.0D0)
      Call DFDZ(U(1,1,1,1), Grad(1,1,1,2), 1, 1.0D0)

*  z-vorticity
*
*        d u_y   d u_x
*  w_z = ----- - -----
*         d x     d y

      Call DFDX(U(1,1,1,2), Grad(1,1,1,3), 0, 1.0D0)
      Call DFDY(U(1,1,1,1), Grad(1,1,1,3), 1, -1.0D0)

*  Dilatation
*
*        d u_x   d u_y   d u_z
*  Dil = ----- + ----- + -----
*         d x     d y     d z

      Grad(:,:,:,4) = 0D0

      Call DFDX(U(1,1,1,1), Grad(1,1,1,4), 0, 1d0)

      Call DFDY(U(1,1,1,2), temp, 0, 1d0)
      Grad(:,:,:,4) = Grad(:,:,:,4) + temp

      Call DFDZ(U(1,1,1,3), temp, 0, 1d0)
      Grad(:,:,:,4) = Grad(:,:,:,4) + temp

      eps  = 1.0D-7
      fac2 = 10D0 

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
                    fac   = 5.0D-1 * 
     .              (1.0D0 - tanh(2.5D0 + fac2 * Grad(J,K,I,4)))

            omega2(J,K,I) = Grad(J,K,I,1)*Grad(J,K,I,1)
     .                    + Grad(J,K,I,2)*Grad(J,K,I,2)
     .                    + Grad(J,K,I,3)*Grad(J,K,I,3)
 
            C_beta(J,K,I) = fac * Grad(J,K,I,4)*Grad(J,K,I,4) 
     .                    / ( Grad(J,K,I,4)*Grad(J,K,I,4)
     .                    + omega2(J,K,I) + eps ) + 1D-20
          End Do
        End Do
      End Do           
 
      C_mu(:,:,:) = 2.0D-3
      C_k(:,:,:)  = 1.0D-2

      End Subroutine Compute_coeff1

***********************************************************
* Subroutine Compute_hyper                                *
*                                                         *
* Computes the nonlinear hyperviscosity based on          *
* Cook and Cabot model (JCP, 2004)                        *
***********************************************************
      Subroutine Compute_hyper(Q)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), U(My,Mz,Mx,3)
      Real*8 P(My,Mz,Mx),   Vs(My,Mz,Mx)
      Real*8 T(My,Mz,Mx),   e(My,Mz,Mx)
 
      Real*8 Sij(My,Mz,Mx,9), S(My,Mz,Mx), d2S(My,Mz,Mx)
      Real*8 d4S(My,Mz,Mx), d4Sg(My,Mz,Mx)
      Real*8 d4e(My,Mz,Mx), d4eg(My,Mz,Mx)
      Real*8 eta4(My,Mz,Mx,3), beta_hyper(My,Mz,Mx)
      Real*8 eta4g(My,Mz,Mx,3)
      Real*8 visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)
      Real*8 C_mu4, C_beta4, C_k4
      Real*8 C_mu_dyn(My,Mz,Mx), C_beta_dyn(My,Mz,Mx)
      Real*8 C_k_dyn(My,Mz,Mx)

      Common /HYPER_COEFF/ C_beta_dyn

      Real*8  fac1, fac2, fac3

      Common /HYPER/ beta_hyper, visc_hyper, k_hyper,
     .               Sij, eta4, C_mu4, C_beta4

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

      Call Compute_coeff1(Q, C_mu_dyn, C_beta_dyn, C_k_dyn)

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
      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
              P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                         ( Q(J,K,I,1) * Q(J,K,I,1)
     .                         + Q(J,K,I,2) * Q(J,K,I,2)
     .                         + Q(J,K,I,3) * Q(J,K,I,3) ) )
              T(J,K,I) = fac2 * Vs(J,K,I) * P(J,K,I)
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

      d4S = abs(d4S)
      d4e = abs(d4e)

* Filtering the field with Gaussian filter (Cook & Cabot, JCP 2004)
c      Call Gaussian_filter(d4S, d4Sg, 1)
c      Call Gaussian_filter(d4e, d4eg, 1)

* Calculating the hyperviscosities
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac3          = (Q(J,K,I,4)/T(J,K,I)) 
     .                    * dsqrt(gamma * P(J,K,I) * Vs(J,K,I))
            eta4(J,K,I,1) = Q(J,K,I,4) * d4S(J,K,I)
            eta4(J,K,I,3) = Q(J,K,I,4) * d4S(J,K,I)
            eta4(J,K,I,2) =       fac3 * d4e(J,K,I)

          End Do
        End Do
      End Do

* Filtering the field with Gaussian filter (Cook & Cabot, JCP 2004)
      Call Gaussian_filter(eta4, eta4g, 3)

* Calculating the hyperviscosities
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            beta_hyper(J,K,I) = C_beta_dyn(J,K,I) *
     .                          eta4g(J,K,I,3) * (DX**6D0)
            visc_hyper(J,K,I) = C_mu_dyn(J,K,I) *
     .                          eta4g(J,K,I,1) * (DX**6D0)
               k_hyper(J,K,I) = C_k_dyn(J,K,I) * 
     .                          eta4g(J,K,I,2) * (DX**5D0)
          End Do
        End Do
      End Do

      Return

      End



