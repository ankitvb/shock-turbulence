****************************************************************
* Calculating flow statistics                                  *
*                                                              *
****************************************************************
      Subroutine Compute_Stats

      Include 'header'
   
      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)
c      Real*8  S1(40), S2(40)

      Common  /RKWORK/ Q, Grad

      Real*8  P(My,Mz,Mx),   Vs(My,Mz,Mx)
      Real*8  T(My,Mz,Mx)
      Real*8  U(My,Mz,Mx,3)
      Real*8  U_temp(My,Mz,Mx)

* Locals
      Integer I,J,K
      Real*8  temp1(My,Mz,Mx),temp2(My,Mz,Mx),temp3(My,Mz,Mx)
      Real*8  tmp, var1, var2, var3
      Real*8  fac1, fac2
      Real*8  lambda_x, lambda_y, lambda_z
      Real*8  skew_x, skew_y, skew_z
      Real*8  u_rms, mu_avg, rho_avg, P_avg

*  Compute specific volume  pressure,
*  temperature
*
*   Vs  = 1 / rho
*
*    P  = (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*           gamma     P
*    T  = ---------  ---
*         gamma - 1  rho
*
      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
              T(J,K,I)   = fac2 * Vs(J,K,I) * P(J,K,I)
          End Do
        End Do
      End Do

*     Compute RMS Values of Density, Pressure and Temperature     
*                 
*     f_rms =  dsqrt( < (f - f_mean)^2 > ) 
*
      temp1 = Q(:,:,:,4)
      Call Compute_RMS(temp1, rho_rms)
      Call Compute_RMS(P, P_rms)
      Call Compute_RMS(T, T_rms) 

*    Compute skewness and flatness
*
*                      <f^3>
*    Skewness(f) =  -----------       
*                    <f^2>^3/2
*                      <f^4>
*    Flatness(f) =  -----------
*                     <f^2>^2 

*    Dilatation skewness 
 
      temp1 = (Grad(:,:,:,4)**3d0)    
      temp2 = (Grad(:,:,:,4)**2d0)
      
      Call Compute_Mean(temp1, var1)
      Call Compute_Mean(temp2, var3)
      Call Compute_RMS(Grad(1,1,1,4), var2)
      
      dilatation_skewness = var1/dsqrt(var3**3d0)
      dilatation_square   = var2
 
*    Dilatation flatness

      temp1 = (Grad(:,:,:,4)**4d0)

      Call Compute_Mean(temp1, var1)

      dilatation_flatness = var1/(var3**2d0)

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

*    Velocity Derivative skewness and Taylor microscale
*
*             d u_x
*    temp1 = -------   
*              d x         
*
      Call DFDX_asym(U(1,1,1,1), temp1, 0, 1d0)
      U_temp = U(:,:,:,1)**2d0

      Call Compute_Mean(U_temp, var1)
      temp3 = temp1**2d0        
      Call Compute_Mean(temp3, var2)
      lambda_x = dsqrt(var1/var2)  

      temp2 = temp1**3d0
      Call Compute_Mean(temp2, var1)

      skew_x = var1/(var2**1.5d0)

*
*             d u_y
*    temp1 = -------
*              d y
*
      Call DFDY_asym(U(1,1,1,2), temp1, 0, 1d0)
      U_temp = U(:,:,:,2)**2d0

      Call Compute_Mean(U_temp, var1)
      temp3 = temp1**2d0
      Call Compute_Mean(temp3, var2)
      lambda_y = dsqrt(var1/var2)

      temp2 = temp1**3d0
      Call Compute_Mean(temp2, var1)

      skew_y = var1/(var2**1.5d0)

*
*             d u_z
*    temp1 = -------
*              d z
*
      Call DFDZ_asym(U(1,1,1,3), temp1, 0, 1d0)
      U_temp = U(:,:,:,3)**2d0

      Call Compute_Mean(U_temp, var1)
      temp3 = temp1**2d0
      Call Compute_Mean(temp3, var2)
      lambda_z = dsqrt(var1/var2)

      temp2 = temp1**3d0
      Call Compute_Mean(temp2, var1)

      skew_z = var1/(var2**1.5d0)

      velocity_derivative_skewness = (skew_x + skew_y + skew_z)/3d0
      taylor_microscale = (lambda_x + lambda_y + lambda_z)/3d0

*     
*     Microscale Reynolds number
*                  rho  *  u_rms * lambda
*     Re_lambda =  ------------------------
*                          mu
*

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             temp1(J,K,I) = ((fac1 * T(J,K,I)) ** vis_exp ) * Re_inv
          End Do
        End Do
      End Do

      Call Compute_Mean(Q(1,1,1,4), rho_avg) 
      Call Compute_Mean(temp1, mu_avg)
      u_rms     = dsqrt((2D0/3D0)*KE)
      Re_lambda = (rho_avg * u_rms * taylor_microscale)/mu_avg

      Return
 
      End 
