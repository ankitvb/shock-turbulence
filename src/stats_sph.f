****************************************************************
* Calculating flow statistics                                  *
*                                                              *
****************************************************************
      Subroutine Compute_ShockTurb_Stats

      Include 'header'
   
      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)

      Common  /RKWORK/ Q, Grad

      Real*8  P(My,Mz,Mx),   Vs(My,Mz,Mx)
      Real*8  T(My,Mz,Mx)
      Real*8  U(My,Mz,Mx,3)
      Real*8  U_temp(My,Mz,Mx)

* Locals
      Integer I,J,K,M
      Real*8  temp1(My,Mz,Mx),temp2(My,Mz,Mx),temp3(My,Mz,Mx)
      Real*8  tmp, var1(N_avg), var2(N_avg), var3(N_avg)
      Real*8  fac1, fac2
      Real*8  skew_x(N_avg), skew_y(N_avg), skew_z(N_avg)
      Real*8  u_rms(N_avg), mu_avg(N_avg), r(N_avg)
      Real*8  rho_shell_avg(N_avg), P_shell_avg(N_avg)

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

      Call Compute_Shell_Avg(P, P_shell_avg)
      Call Compute_Shell_Avg(Q(1,1,1,4), rho_shell_avg)

*     Compute RMS Values of Density, Pressure and Temperature     
*                 
*     f_rms =  dsqrt( < (f - f_mean)^2 > ) 
*
      temp1 = Q(:,:,:,4)
      Call Compute_shell_RMS(temp1, rho_shell_rms)
      Call Compute_shell_RMS(P, P_shell_rms)
      Call Compute_shell_RMS(T, T_shell_rms) 

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
      
      Call Compute_shell_avg(temp1, var1)
      Call Compute_shell_avg(temp2, var3)
      Call Compute_shell_RMS(Grad(1,1,1,4), var2)
      
      dilatation_shell_skewness = var1/dsqrt(var3**3d0)
      dilatation_shell_square   = var2
 
*    Dilatation flatness

      temp1 = (Grad(:,:,:,4)**4d0)

      Call Compute_shell_avg(temp1, var1)

      dilatation_shell_flatness = var1/(var3**2d0)

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
      Call DFDX(U(1,1,1,1), temp1, 0, 1d0)
      U_temp = U(:,:,:,1)**2d0

      Call Compute_shell_avg(U_temp, var1)
      temp3 = temp1**2d0        
      Call Compute_shell_avg(temp3, var2)

      temp2 = temp1**3d0
      Call Compute_shell_avg(temp2, var1)

      Do M = 1,N_avg      
        skew_x(M) = var1(M)/(var2(M)**1.5d0)
      End Do

*
*             d u_y
*    temp1 = -------
*              d y
*
      Call DFDY(U(1,1,1,2), temp1, 0, 1d0)
      U_temp = U(:,:,:,2)**2d0

      Call Compute_shell_avg(U_temp, var1)
      temp3 = temp1**2d0
      Call Compute_shell_avg(temp3, var2)

      temp2 = temp1**3d0
      Call Compute_shell_avg(temp2, var1)

      Do M = 1,N_avg
        skew_y(M) = var1(M)/(var2(M)**1.5d0)
      End Do

*
*             d u_z
*    temp1 = -------
*              d z
*
      Call DFDZ(U(1,1,1,3), temp1, 0, 1d0)
      U_temp = U(:,:,:,3)**2d0

      Call Compute_shell_avg(U_temp, var1)
      temp3 = temp1**2d0
      Call Compute_shell_avg(temp3, var2)

      temp2 = temp1**3d0
      Call Compute_shell_avg(temp2, var1)

      Do M = 1,N_avg
        skew_z(M) = var1(M)/(var2(M)**1.5d0)
      End Do

      velocity_derivative_shell_skewness = 
     .                    (skew_x + skew_y + skew_z)/3d0

* Setting bin levels
* There are N/4 +4 levels, spaced Delta * sqrt(3) apart in r-space

      r(1) = DX * sqrt(3D0) / 2D0

      Do M = 2,N_avg
        r(M) = r(M-1) + DX * sqrt(3D0)
      End Do


      if(rank.EQ.0)then
        Do M = 1, N_avg
          Write(stats2_unit,6) r(M),rho_shell_rms(M), P_shell_rms(M),
     .                         rho_shell_avg(M), P_shell_avg(M)
          Write(stats1_unit,7) velocity_derivative_shell_skewness(M),
     .                         dilatation_shell_skewness(M), 
     .                         dilatation_shell_flatness(M),
     .                         dilatation_shell_square(M)
        End Do
      endif

 6    Format(1x, 5(ES13.5))
 7    Format(1x, 4(ES13.5))

      Return
 
      End 
