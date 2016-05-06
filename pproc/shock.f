      Subroutine Add_Shock_Expanding(Q)

      include 'header'

*  Globals
      Real*8  Q(My,Mz,Mx,5), P(My,Mz,Mx)
      Real*8  U(My,Mz,Mx,3), Vs(My,Mz,Mx)

*  Locals
      Integer      I, J, K, L, M, N
      Real*8       fac1

      Integer      II, JJ, KK
      Real*8       IX, IY, IZ

      Real*8 fac, base_pressure
      Real*8 rad, eps

      eps = 0.1D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            II = xrank * Mx + I
            JJ = yrank * My + J
            KK = zrank * Mz + K

            IX = REAL(II - Nx/2) - 0.5D0
            IY = REAL(JJ - Ny/2) - 0.5D0
            IZ = REAL(KK - Nz/2) - 0.5D0

            rad = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

! Compute primitive variables from turbulence field
            Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
            U(J,K,I,1) = Q(J,K,I,1) / Q(J,K,I,4)
            U(J,K,I,2) = Q(J,K,I,2) / Q(J,K,I,4)
            U(J,K,I,3) = Q(J,K,I,3) / Q(J,K,I,4)

            P(J,K,I) = (gamma - 1D0)*(Q(J,K,I,5) - 5.0D-1*Vs(J,K,I) *
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) ) )

! Add shock to primitive variables
            if(     (abs(IX).EQ.0.5D0)
     .         .AND.(abs(IY).EQ.0.5D0)
     .         .AND.(abs(IZ).EQ.0.5D0)  )then
              P(J,K,I)  = P(J,K,I) + (gamma-1D0) * Q(J,K,I,4) * 25D0 *
     .                 * 0.1528415451D0 * exp(-rad*rad/(eps*eps)) ! Factor 4 for 512 grid, 1 for 256
     .                 / (eps * eps * eps)  
            elseif( (abs(IX).LT.2D0)
     .         .AND.(abs(IY).LT.2D0)
     .         .AND.(abs(IZ).LT.2D0)  )then
              P(J,K,I)  = P(J,K,I) + (gamma-1D0) * Q(J,K,I,4) * 25D0 *
     .                 * 0.1528415451D0 * exp(-rad*rad/(eps*eps)) ! Factor 2 for 512 grid, 1 for 256
     .                 / (eps * eps * eps)
            elseif( (abs(IX).LT.3D0)
     .         .AND.(abs(IY).LT.3D0)
     .         .AND.(abs(IZ).LT.3D0)  )then
              P(J,K,I)  = P(J,K,I) + (gamma-1D0) * Q(J,K,I,4)* 25D0 *
     .                 * 0.1528415451D0 * exp(-rad*rad/(eps*eps)) ! Factor 1.5 for 512 grid, 1 for 256
     .                 / (eps * eps * eps)
            else 
              P(J,K,I)  = P(J,K,I) + (gamma-1D0) * Q(J,K,I,4)* 25D0 *
     .                 * 0.1528415451D0 * exp(-rad*rad/(eps*eps)) 
     .                 / (eps * eps * eps)
            endif

! And recompute conservative variables
            Q(J,K,I,1)   = Q(J,K,I,4) * U(J,K,I,1)
            Q(J,K,I,2)   = Q(J,K,I,4) * U(J,K,I,2)
            Q(J,K,I,3)   = Q(J,K,I,4) * U(J,K,I,3)

            Q(J,K,I,5) = ( P(J,K,I)/(gamma - 1D0) ) +
     .        (5D-1/Q(J,K,I,4)) * ( Q(J,K,I,1)*Q(J,K,I,1)
     .                          +   Q(J,K,I,2)*Q(J,K,I,2)
     .                          +   Q(J,K,I,3)*Q(J,K,I,3) )

          End Do
        End Do
      End Do

      Return

      End Subroutine Add_Shock_Expanding

***********************************************************************************

      Subroutine Add_Shock_Converging(Q)

      include 'header'

*  Globals
      Real*8  Q(My,Mz,Mx,5), P(My,Mz,Mx)
      Real*8  U(My,Mz,Mx,3), Vs(My,Mz,Mx)
      Real*8  rad(My,Mz,Mx)

*  Locals
      Integer      I, J, K, L, M, N
      Real*8       fac1

      Integer      II, JJ, KK
      Real*8       IX, IY, IZ

      Real*8 fac, base_pressure, Ms, Ms1, r_shock

      Ms  = 4.0D0
      Ms1 = 4.0D0!/sqrt(3D0)
      r_shock = 0.75D0 * (x_leng / 2D0)
 
      Do I = 1,Mx
       Do J = 1,My
         Do K = 1,Mz

          II = Mx * xrank + I
          JJ = My * yrank + J
          KK = Mz * zrank + K

          IX = REAL(II - Nx/2) - 0.5D0
          IY = REAL(JJ - Ny/2) - 0.5D0
          IZ = REAL(KK - Nz/2) - 0.5D0

          rad(J,K,I) = sqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

! Compute primitive variables from turbulence field
          Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
          U(J,K,I,1) = Q(J,K,I,1) / Q(J,K,I,4)
          U(J,K,I,2) = Q(J,K,I,2) / Q(J,K,I,4)
          U(J,K,I,3) = Q(J,K,I,3) / Q(J,K,I,4)

          P(J,K,I) = (gamma - 1D0)*(Q(J,K,I,5) - 5.0D-1*Vs(J,K,I) *
     .                     ( Q(J,K,I,1) * Q(J,K,I,1)
     .                     + Q(J,K,I,2) * Q(J,K,I,2)
     .                     + Q(J,K,I,3) * Q(J,K,I,3) ) )

! Add shock to primitive variables
          Q(J,K,I,4) = Q(J,K,I,4) + Q(J,K,I,4) * 5.0D-1* 
     .       (1D0 + tanh(50D0*(rad(J,K,I)-r_shock)/3D0) ) *       ! factor 25 for 512 grid, 50 for 256
     .       ( (gamma+1D0)*Ms*Ms/((gamma-1D0)*Ms*Ms+2D0) - 1D0)

          P(J,K,I) = P(J,K,I) +  P(J,K,I) * 5.0D-1 *
     .       (1D0 + tanh(50D0*(rad(J,K,I)-r_shock)/3D0)) *        ! factor 25 for 512 grid, 50 for 256
     .       (2D0*gamma/(gamma+1D0)) * (Ms1*Ms1 - 1D0)

! And recompute conservative variables
          Q(J,K,I,1)   = Q(J,K,I,4) * U(J,K,I,1)
          Q(J,K,I,2)   = Q(J,K,I,4) * U(J,K,I,2)
          Q(J,K,I,3)   = Q(J,K,I,4) * U(J,K,I,3)

          Q(J,K,I,5) = ( P(J,K,I)/(gamma - 1D0) ) +
     .      (5D-1/Q(J,K,I,4)) * ( Q(J,K,I,1)*Q(J,K,I,1)
     .                        +   Q(J,K,I,2)*Q(J,K,I,2)
     .                        +   Q(J,K,I,3)*Q(J,K,I,3) )
         End Do
       End Do
      End Do

      Return
 
      End Subroutine Add_Shock_Converging
************************************************************************
      Subroutine Generate_Shock_Surface

      Include 'header'
      Include 'mpif.h'

      Real*8 beta_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper 

      Real*8 rad(My,Mz,Mx)
      Real*8 beta_threshold, beta_max
      Real*8 xc_local, yc_local, zc_local, n_local,
     .       rad_local
      Real*8 xc, yc, zc, n, rad_avg, rad_skew, rad_skew_rms
      
      Integer I, J, K
      Integer ierr

      Call Compute_Max(beta_hyper, beta_max)

      beta_threshold = 0.1D0 * beta_max 
     
      xc_local = 0D0; yc_local = 0D0; zc_local = 0D0
      n_local = 0D0
      xc = 0D0; yc = 0D0; zc = 0D0; n = 0D0 

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            if(beta_hyper(J,K,I) .GT. beta_threshold)then
               xc_local = xc_local + xcor(I)
               yc_local = yc_local + ycor(J)
               zc_local = zc_local + zcor(K)
                n_local = n_local  + 1D0
            endif
          End Do
        End Do
      End Do
   
      Call MPI_ALLREDUCE(xc_local, xc, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(yc_local, yc, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(zc_local, zc, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(n_local, n, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)

      xc = xc/n; yc = yc/n; zc = zc/n

      rad_local = 0D0
      rad_avg   = 0D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            If(beta_hyper(J,K,I).GT.beta_threshold)then
              rad(J,K,I) = sqrt((xcor(I)-xc)**2d0 + (ycor(J)-yc)**2d0
     .              + (zcor(K)-zc)**2d0)
              rad_local = rad_local + rad(J,K,I)
            Endif
          End Do
        End Do
      End Do 
   
      Call MPI_ALLREDUCE(rad_local, rad_avg, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)
      
      rad_avg = rad_avg/n

      rad_skew     = 0d0
      rad_skew_rms = 0d0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            if(beta_hyper(J,K,I).GT.beta_threshold)then
              rad_skew = rad_skew + (rad(J,K,I) - rad_avg)**2d0
            endif
          End Do
        End Do
      End Do

      Call MPI_ALLREDUCE(rad_skew, rad_skew_rms, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)

      rad_skew_rms = sqrt(rad_skew_rms)/(n*rad_avg)

      if(rank .EQ. 0)then
        OPEN(UNIT=7,FILE="shock_stats",FORM="FORMATTED",
     .       STATUS="UNKNOWN", POSITION="APPEND")
        WRITE(7,*) rad_avg, rad_skew_rms
        CLOSE(7) 
      endif

      Return

      End Subroutine Generate_Shock_Surface

