************************************************************************
*  Subroutine: Sponge_Forcing                                          *
*                                                                      *
*  Function:   It adds on the sponge forcing terms to the right hand   *
*              time derivative vector.                                 *
************************************************************************
      Subroutine Sponge_Forcing(Q, Qrhs)

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8  Qref(My,Mz,Mx,5), Q_aux(My,Mz,Mx,5)
      Real*8  P(My,Mz,Mx), Vs(My,Mz,Mx)
      Real*8  sigma_x(My,Mz,Mx)
      Real*8  tmp, P0

      Real*8  xc, yc, zc

      Integer I, J, K, L, M
      Integer gI(Mx)
      Real*8  A_sponge, s_exp
      Real*8  frac_sponge_x, frac_sponge_y, frac_sponge_z
      Integer I_sponge_start_left, I_sponge_start_right
      Integer I_sponge_length

      Integer status(MPI_STATUS_SIZE), req(2)

      A_sponge      = 5.0D0                        ! Sponge coefficient
      s_exp         = 2.0D0                        ! Power of sponge polynomial
      frac_sponge_x = 0.4D0                        ! Fraction of domain in x treated with sponge BC

      I_sponge_start_left  = CEILING((frac_sponge_x*Nx)/2D0)
      I_sponge_start_right = Nx - I_sponge_start_left

      I_sponge_length = I_sponge_start_left

* Setting reference state
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

      P0 = 0.71463D0

      Qref(:,:,:,1)   = 5.75D-1
      Qref(:,:,:,2:3) = 0D0
      Qref(:,:,:,4)   = 1.000482D0
      Qref(:,:,:,5)   = P0 / (gamma-1D0) + 5.0D-1 *
     .                  (1D0/1.000482)*3.30625D-1 

* Map local indices to global ones and setting sponge strength
      Do I = 1,Mx
       Do J = 1,My
        Do K = 1,Mz
         gI(I) = xrank * Mx + I

         if( gI(I).LE.I_sponge_start_left)then
           sigma_x(J,K,I) = A_sponge*
     .     ( REAL(I_sponge_start_left-gI(I)+1)/
     .       REAL(I_sponge_length) )**s_exp
         elseif(gI(I).GE.I_sponge_start_right)then
           sigma_x(J,K,I) = A_sponge*
     .     ( REAL(gI(I)-I_sponge_start_right)/
     .       REAL(I_sponge_length) )**s_exp
         else
           sigma_x(J,K,I) = 0D0
         endif

        End Do
       End Do
      End Do

       Do I = 1,Mx
        Do J = 1,My
         Do K = 1,Mz
          Do L = 1,5
           Qrhs(J,K,I,L) = Qrhs(J,K,I,L)
     .           - ( sigma_x(J,K,I) )
     .           *  ( Q(J,K,I,L) - Qref(J,K,I,L) )
          End Do
         End Do
        End Do
       End Do

      End Subroutine Sponge_Forcing

