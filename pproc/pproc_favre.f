**********************************************************
* Post-processing stuff                                  *
**********************************************************
      Subroutine Compute_Favre_Stats
 
      Include 'header'
      Include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)

      Common /RKWORK/ Q, Grad

      Real*8 Y(My,Mz,Mx,1)

      Common /SPECIES/ Y     

      Real*8 beta_hyper(My,Mz,Mx), visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx), D_hyper(My,Mz,Mx)
      
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper, D_hyper

      Real*8 gamma_eff(My,Mz,Mx), gradrho(My,Mz,Mx,3)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 S(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3), Ur(My,Mz,Mx), rad(My,Mz,Mx) 
      Real*8 omr(My,Mz,Mx)
      Real*8 omrp2(My,Mz,Mx), omtp2(My,Mz,Mx), vort(My,Mz,Mx)
      Real*8 urp2(My,Mz,Mx), utp2(My,Mz,Mx), tke(My,Mz,Mx)      
      Real*8 Yp2(My,Mz,Mx)
 
      Real*8 buf(My,Mz,Mx,7)
      Real*8 fac1, fac2, tmp
 
      Common /THERMO/ gamma_eff

      Real*8 r(N_avg+1)
      Real*8 rho_shell_avg(N_avg), P_shell_avg(N_avg)
      Real*8 Ur_shell_avg(N_avg), Y_shell_avg(N_avg)
      
      Real*8 tke_shell_avg(N_avg), vort_shell_avg(N_avg)
      Real*8 omr_shell_avg(N_avg), Yp2_shell_avg(N_avg) 
      Real*8 omrp2_shell_avg(N_avg), omtp2_shell_avg(N_avg) 
      Real*8 urp2_shell_avg(N_avg), utp2_shell_avg(N_avg) 

      Real*8 beta_shell_avg(N_avg), visc_shell_avg(N_avg)
      Real*8 k_shell_avg(N_avg), D_shell_avg(N_avg)

      Real*8 tke_total, hmix

      Integer I, J, K, M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

      Integer status_arr(MPI_STATUS_SIZE,6), ierr, req(6), rreq(6)
      Integer status(MPI_STATUS_SIZE)

*  Compute specific volume (avoid division in computations), pressure,
*  temperature and transverse velocity
*
*   Vs  = 1 / rho
*
*    P  = rho * (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*           gamma1     P        1
*    T  = ---------   ---  ------------ * rho
*         gamma1 - 1  rho  (Y1 + eps*Y2)
*
*   Ux5 = u_y
*
*  Flux = (rho u_x) u_y
* 
* For Favre Averaged stats 

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             fac1 = gamma_eff(J,K,I) - 1.0D0
             fac2 = gamma_eff(J,K,I) / fac1
             fac3 = gamma1 / (gamma1 - 1D0)
             fac3 = fac3 *
     .       1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

             Y(J,K,I,1) = Q(J,K,I,4) * Y(J,K,I,1)
             Vs(J,K,I)  = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)  = fac1 *  Q(J,K,I,4) *
     .                           ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
              T(J,K,I)  = fac3 * P(J,K,I)
              S(J,K,I)  = (1/gamma_eff(J,K,I))
     .                   *log(P(J,K,I))-log(Q(J,K,I,4))
          End Do
        End Do
      End Do

*  Convert momenta to velocities
      Do I = 1,Mx
       Do J = 1,My
         Do K = 1,Mz
          tmp = 1.0D0 / Q(J,K,I,4)
          U(J,K,I,1) = Q(J,K,I,1) * tmp
          U(J,K,I,2) = Q(J,K,I,2) * tmp
          U(J,K,I,3) = Q(J,K,I,3) * tmp

          II = xrank * Mx + I
          JJ = yrank * My + J
          KK = zrank * Mz + K

          IX = DBLE(II-1) 
          IY = DBLE(JJ-1)
          IZ = DBLE(KK-1)

          rad(J,K,I) = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

          if((II.EQ.1).AND.(JJ.EQ.1).AND.(KK.EQ.1))then
            Ur(J,K,I) = 0D0
          else
            Ur(J,K,I) = Q(J,K,I,4) * 
     .                  (xcor(I)/rad(J,K,I)) * U(J,K,I,1)
     .                + (ycor(J)/rad(J,K,I)) * U(J,K,I,2)
     .                + (zcor(K)/rad(J,K,I)) * U(J,K,I,3)
          endif

          if((II.EQ.1).AND.(JJ.EQ.1).AND.(KK.EQ.1))then
            omr(J,K,I) = 0D0
          else
            omr(J,K,I) = Q(J,K,I,4) *
     .                   Grad(J,K,I,1) * (IX*DX/rad(J,K,I))
     .                 + Grad(J,K,I,2) * (IY*DY/rad(J,K,I))
     .                 + Grad(J,K,I,3) * (IZ*DZ/rad(J,K,I))
          endif

        End Do
       End Do
      End Do

! Compute shell averages
      r(1) = 0D0

      Do M = 2,N_avg+1
        if(M .LE. 4)then
          r(M) = r(M-1) + DX * 2d0 * sqrt(3D0)
        else
          r(M) = r(M-1) + DX * sqrt(3D0)
        endif
      End Do

      Call Compute_Shell_Avg(Q(1,1,1,4), rho_shell_avg)
      Call Compute_Shell_Avg(P, P_shell_avg)
      Call Compute_Shell_Avg(Ur, Ur_shell_avg)
      Call Compute_Shell_Avg(Y, Y_shell_avg)
      Call Compute_Shell_Avg(omr, omr_shell_avg) 

* Computing favre averages of mean quantities:
* \tilde(P), \tilde(Ur), \tilde(Y), \tilde(omr)
      Do M = 1,N_avg
        P_shell_avg(M)   = P_shell_avg(M) / rho_shell_avg(M)
        Ur_shell_avg(M)  = Ur_shell_avg(M) / rho_shell_avg(M)
        Y_shell_avg(M)   = Y_shell_avg(M) / rho_shell_avg(M)
        omr_shell_avg(M) = omr_shell_avg(M) / rho_shell_avg(M)
      End Do

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz

            Ur(J,K,I)  = Ur(J,K,I)  / Q(J,K,I,4)
            omr(J,K,I) = omr(J,K,I) / Q(J,K,I,4) 
            Y(J,K,I,1) = Y(J,K,I,1) / Q(J,K,I,4)

            Do M = 1,N_avg
              if((rad(J,K,I).GE.r(M)).AND.(rad(J,K,I).LT.r(M+1)))then
               II = xrank * Mx + I
               JJ = yrank * My + J
               KK = zrank * Mz + K

               IX = DBLE(II-1)
               IY = DBLE(JJ-1)
               IZ = DBLE(KK-1)
 
               urp2(J,K,I)     = (Ur(J,K,I) - Ur_shell_avg(M))**2d0
               if((II.GT.1).AND.(JJ.GT.1).AND.(KK.GT.1))then
                 tmp1            = U(J,K,I,1) -
     .                             Ur(J,K,I) * (IX*DX/rad(J,K,I))
                 tmp2            = U(J,K,I,2) -
     .                             Ur(J,K,I) * (IY*DY/rad(J,K,I))
                 tmp3            = U(J,K,I,3) -
     .                             Ur(J,K,I) * (IZ*DZ/rad(J,K,I))
               else
                 tmp1 = 0D0; tmp2 = 0D0; tmp3 = 0D0
               endif

               utp2(J,K,I)     = tmp1**2d0 + tmp2**2d0 + tmp3**2d0

               omrp2(J,K,I)    = (omr(J,K,I) - omr_shell_avg(M))**2d0

               if((II.GT.1).AND.(JJ.GT.1).AND.(KK.GT.1))then
                 tmp1            = Grad(J,K,I,1) -
     .                             omr(J,K,I) * (IX*DX/rad(J,K,I))
                 tmp2            = Grad(J,K,I,2) -
     .                             omr(J,K,I) * (IY*DY/rad(J,K,I))
                 tmp3            = Grad(J,K,I,3) -
     .                             omr(J,K,I) * (IZ*DZ/rad(J,K,I))
               else
                 tmp1 = 0D0; tmp2 = 0D0; tmp3 = 0D0
               endif

               omtp2(J,K,I)    = tmp1**2d0 + tmp2**2d0 + tmp3**2d0

               tke(J,K,I)      = utp2(J,K,I) + urp2(J,K,I)
               vort(J,K,I)     = omtp2(J,K,I) + omrp2(J,K,I)

               Yp2(J,K,I)      = (Y(J,K,I,1) - Y_shell_avg(M))**2D0
              endif
            End Do
          End Do
        End Do
      End Do

      Call Compute_Shell_Avg(omrp2, omrp2_shell_avg) 
      Call Compute_Shell_Avg(omtp2, omtp2_shell_avg) 

      Call Compute_Shell_Avg(urp2, urp2_shell_avg) 
      Call Compute_Shell_Avg(utp2, utp2_shell_avg) 

!      Call Compute_Shell_Avg(tke, tke_shell_avg)
!      Call Compute_Shell_Avg(vort, vort_shell_avg)
      Call Compute_Shell_Avg(Yp2, Yp2_shell_avg)

      Call Compute_Shell_Avg(beta_hyper, beta_shell_avg)
      Call Compute_Shell_Avg(visc_hyper, visc_shell_avg)
      Call Compute_Shell_Avg(k_hyper, k_shell_avg)
      Call Compute_Shell_Avg(D_hyper, D_shell_avg)

      if(rank.EQ.0)then
          OPEN(UNIT=5,FILE="avg_stats_favre",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          OPEN(UNIT=6,FILE="srm_stats_favre",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
        Do M = 1, N_avg
          WRITE(5,6) rho_shell_avg(M), P_shell_avg(M),
     .               Ur_shell_avg(M), Y_shell_avg(M) 
          WRITE(6,7) urp2_shell_avg(M),  utp2_shell_avg(M),
     .               omrp2_shell_avg(M), omtp2_shell_avg(M), 
     .               Yp2_shell_avg(M)
        End Do
        CLOSE(5); CLOSE(6)
      endif

 6    Format(1x, 4(ES13.5))
 7    Format(1x, 5(ES13.5))

      End Subroutine Compute_Favre_Stats
