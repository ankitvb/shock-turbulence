************************************************************************
*  Subroutine: Sponge_Forcing                                          *
*                                                                      *
*  Function:   It adds on the sponge forcing terms to the right hand   *
*              time derivative vector.                                 *
************************************************************************
      Subroutine Sponge_Forcing(Q, Qrhs, Y, Yrhs)

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8  Y(My,Mz,Mx,1), Yrhs(My,Mz,Mx,1)
      Real*8  Qref(My,Mz,Mx,5), Q_aux(My,Mz,Mx,5), Yref(My,Mz,Mx,1)
      Real*8  P(My,Mz,Mx), Vs(My,Mz,Mx)
      Real*8  sigma_x(My,Mz,Mx), sigma_y(My,Mz,Mx), sigma_z(My,Mz,Mx)
      Real*8  tmp, Ms, base_pressure, c0
      Real*8  rho_ref, P_ref, vel_ref, Y_ref      

      Real*8  xc, yc, zc

      Integer I, J, K, L, M
      Integer gI(Mx), gJ(My), gK(Mz)
      Real*8  A_sponge, s_exp
      Real*8  frac_sponge_x, frac_sponge_y, frac_sponge_z
      Integer I_sponge_start, J_sponge_start, K_sponge_start
      Integer I_sponge_length, J_sponge_length, K_sponge_length

      Integer status(MPI_STATUS_SIZE), req(2)

      base_pressure = 0.71429D0
      Ms            = 1.8D0

      A_sponge      = 4.0D0                        ! Sponge coefficient
      s_exp         = 3.0D0                        ! Power of sponge polynomial
      frac_sponge_x = 0.1D0                        ! Fraction of domain in x treated with sponge BC
      frac_sponge_y = 0.1D0                        ! Fraction of domain in y treated with sponge BC
      frac_sponge_z = 0.1D0                        ! Fraction of domain in z treated with sponge BC

      vel_ref = 0D0
      rho_ref = (gamma1+1D0)*Ms*Ms/((gamma1-1D0)*Ms*Ms+2D0)  
      P_ref   = ((2D0*gamma1/(gamma1+1D0))*(Ms*Ms-1D0)+1D0)*
     .            base_pressure      
      Y_ref   = 1.0D0 

      I_sponge_start  = Nx - CEILING(frac_sponge_x*Nx)
      J_sponge_start  = Ny - CEILING(frac_sponge_y*Ny)
      K_sponge_start  = Nz - CEILING(frac_sponge_z*Nz)

      I_sponge_length = Nx - I_sponge_start
      J_sponge_length = Ny - J_sponge_start
      K_sponge_length = Nz - K_sponge_start

* Setting reference state
      fac1 = gamma1 - 1D0      

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
           Vs(J,K,I) = 1D0 / Q(J,K,I,4)
            P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                      ( Q(J,K,I,1) * Q(J,K,I,1)
     .                      + Q(J,K,I,2) * Q(J,K,I,2)
     .                      + Q(J,K,I,3) * Q(J,K,I,3) ) )
          End Do
        End Do
      End Do

      Qref(:,:,:,1:3) =  vel_ref 
      Qref(:,:,:,4)   =  rho_ref
      Qref(:,:,:,5)   =  P_ref / (gamma1-1D0)
      Yref(:,:,:,1)   =  Y_ref

!      print *, vel_ref, rho_ref, P_ref, Y_ref

* Map local indices to global ones and setting sponge strength
      Do I = 1,Mx
       Do J = 1,My
        Do K = 1,Mz
         gI(I) = xrank * Mx + I
         gJ(J) = yrank * My + J
         gK(K) = zrank * Mz + K

         if(gI(I).GE.I_sponge_start)then
           sigma_x(J,K,I) = A_sponge*
     .     ( REAL(gI(I)-I_sponge_start)/
     .       REAL(I_sponge_length) )**s_exp
         else
           sigma_x(J,K,I) = 0D0
         endif

         if(gJ(J).GE.J_sponge_start)then
           sigma_y(J,K,I) = A_sponge*
     .     ( REAL(gJ(J)-J_sponge_start)/
     .       REAL(J_sponge_length) )**s_exp
         else
           sigma_y(J,K,I) = 0D0
         endif

         if(gK(K).GE.K_sponge_start)then
           sigma_z(J,K,I) = A_sponge*
     .     ( REAL(gK(K)-K_sponge_start)/
     .       REAL(K_sponge_length) )**s_exp
         else
           sigma_z(J,K,I) = 0D0
         endif

        End Do
       End Do
      End Do

       Do I = 1,Mx
        Do J = 1,My
         Do K = 1,Mz
          Do L = 1,5
           Qrhs(J,K,I,L) = Qrhs(J,K,I,L)
     .           - ( sigma_x(J,K,I) + sigma_y(J,K,I) + sigma_z(J,K,I) )
     .          *  ( Q(J,K,I,L) - Qref(J,K,I,L) )
          End Do
          Yrhs(J,K,I,1) = Yrhs(J,K,I,1)
     .      - ( sigma_x(J,K,I) + sigma_y(J,K,I) + sigma_z(J,K,I) ) 
     .      *  ( Q(J,K,I,4) * Y(J,K,I,1) - Q(J,K,I,4) * Yref(J,K,I,1) )
         End Do
        End Do
       End Do

      End Subroutine Sponge_Forcing

