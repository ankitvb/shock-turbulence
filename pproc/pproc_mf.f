**********************************************************
* Post-processing stuff                                  *
**********************************************************
      Subroutine Compute_SRM_Stats
 
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

      Real*8 Y_pdf(N_pdf,2), rho_pdf(N_pdf,2)
      Real*8 tke_pdf(N_pdf,2), vort_pdf(N_pdf,2)
 
      Real*8 tke_total, hmix
      Real*8 tke_total_local, tke_total_global
      Real*8 vort_total_local, vort_total_global
      Real*8 count_local, count_global

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
*    P  = (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*           gamma1     P        1
*    T  = ---------   ---  ------------
*         gamma1 - 1  rho  (Y1 + eps*Y2)
*
*   Ux5 = u_y
*
*  Flux = (rho u_x) u_y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             fac1 = gamma_eff(J,K,I) - 1.0D0
             fac2 = gamma_eff(J,K,I) / fac1
             fac3 = gamma1 / (gamma1 - 1D0)
             fac3 = fac3 *
     .       1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

             Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
              T(J,K,I)   = fac3 * Vs(J,K,I) * P(J,K,I)
              S(J,K,I) = (1/gamma_eff(J,K,I))
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
            Ur(J,K,I) = (xcor(I)/rad(J,K,I)) * U(J,K,I,1)
     .                + (ycor(J)/rad(J,K,I)) * U(J,K,I,2)
     .                + (zcor(K)/rad(J,K,I)) * U(J,K,I,3)
          endif

          if((II.EQ.1).AND.(JJ.EQ.1).AND.(KK.EQ.1))then
            omr(J,K,I) = 0D0
          else
            omr(J,K,I) = Grad(J,K,I,1) * (IX*DX/rad(J,K,I))
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

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
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

      tke_total_local  = 0D0; vort_total_local  = 0D0
      tke_total_global = 0D0; vort_total_global = 0D0
      count_local = 0D0; count_global = 0D0

      Do I = 1,Mx
       Do J = 1,My
        Do K = 1,Mz
         if(Y(J,K,I,1)*(1D0-Y(J,K,I,1)).GT.0.01D0)then
           tke_total_local  = tke_total_local + tke(J,K,I)/Q(J,K,I,4)
           vort_total_local = vort_total_local + vort(J,K,I)/Q(J,K,I,4)
           count_local      = count_local + 1D0
         endif
        End Do
       End Do
      End Do

      Call MPI_ALLREDUCE(tke_total_local, tke_total_global, 1,
     .                   MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(vort_total_local, vort_total_global, 1,
     .                   MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(count_local, count_global, 1,
     .                   MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      tke_total_global  = tke_total_global / count_global
      vort_total_global = vort_total_global / count_global

      if(rank.EQ.0)then
        OPEN(UNIT=5,FILE="eps_stats",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
        WRITE(5,6) tke_total_global, vort_total_global
        CLOSE(5)
      endif


      Call Compute_Shell_Avg(omrp2, omrp2_shell_avg) 
      Call Compute_Shell_Avg(omtp2, omtp2_shell_avg) 

      Call Compute_Shell_Avg(urp2, urp2_shell_avg) 
      Call Compute_Shell_Avg(utp2, utp2_shell_avg) 

      Call Compute_Shell_Avg(tke, tke_shell_avg)
      Call Compute_Shell_Avg(vort, vort_shell_avg)
      Call Compute_Shell_Avg(Yp2, Yp2_shell_avg)

      Call Compute_Shell_Avg(beta_hyper, beta_shell_avg)
      Call Compute_Shell_Avg(visc_hyper, visc_shell_avg)
      Call Compute_Shell_Avg(k_hyper, k_shell_avg)
      Call Compute_Shell_Avg(D_hyper, D_shell_avg)

      if(rank.EQ.0)then
          OPEN(UNIT=5,FILE="avg_stats",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          OPEN(UNIT=6,FILE="srm_stats",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          OPEN(UNIT=8,FILE="hyper_stats",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
        Do M = 1, N_avg
          WRITE(5,6) rho_shell_avg(M), P_shell_avg(M),
     .               Ur_shell_avg(M), Y_shell_avg(M) 
          WRITE(6,7) urp2_shell_avg(M),  utp2_shell_avg(M),
     .               omrp2_shell_avg(M), omtp2_shell_avg(M), 
     .               Yp2_shell_avg(M)
          WRITE(8,6) beta_shell_avg(M), visc_shell_avg(M),
     .               k_shell_avg(M), D_shell_avg(M)
        End Do
        CLOSE(5); CLOSE(6); CLOSE(8)
      endif

      buf(:,:,:,1) = Q(:,:,:,4)
      buf(:,:,:,2) = Ur
      buf(:,:,:,3) = P
      buf(:,:,:,4) = Y(:,:,:,1)

      Call Print_Diagonal(buf,4)
      Call Print_Line_X(buf,4)
      Call Print_Line_Y(buf,4)
      Call Print_Line_Z(buf,4)

      Call Compute_PDF(Y(1,1,1,1), Y_pdf)
      Call Compute_PDF(tke, tke_pdf)
      Call Compute_PDF(vort, vort_pdf)

!      if(rank.EQ.0)then
!        Do M=1,N_pdf
!          print *, Y_pdf(M,1), tke_pdf(M,1), vort_pdf(M,1)
!        End Do
!      endif 

      if(rank.EQ.0)then
          OPEN(UNIT=8,FILE="pdf_stats",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
       Do M = 1, N_pdf
          WRITE(8,8) Y_pdf(M,1), Y_pdf(M,2),
     .               tke_pdf(M,1), tke_pdf(M,2), 
     .               vort_pdf(M,1), vort_pdf(M,2)
       End Do
        CLOSE(8)
      endif

 6    Format(1x, 4(ES13.5))
 7    Format(1x, 5(ES13.5))
 8    Format(1x, 6(ES13.5))

      End Subroutine Compute_SRM_Stats

*******************************************************************
* Subroutine : Print_Diagonal                                     *
*                                                                 *
* Function : Writes the diagonal of an array                      *
*            Assumes Mx = My = Mz, Px = Py = Pz                   *
*******************************************************************
      Subroutine Print_Diagonal(F, N)

      Include 'header'
      Include 'mpif.h'

      Integer I, J, K, L, M, N

      Real*8 F(My,Mz,Mx,N), F_aux(Mx,N)
      Real*8 F_diag(Nx,N), temp(Mx,N)

      Integer status(MPI_STATUS_SIZE), ierr, req(2)

      M = 1

* Put all the diagonal elements of all diagonal processes into temporary array
      if((xrank.EQ.yrank) .AND. (xrank.EQ.zrank))then

       Do I = 1,Mx
        Do J = 1,My
         Do K = 1,Mz
           if((I .EQ. J) .AND. (I .EQ. K))then
            Do L = 1,N
             temp(M,L) = F(J,K,I,L)
            End Do
            M = M + 1
           endif
         End Do
        End Do
       End Do

* Send all diagonal temp arrays to root process which receives them in order
       if(rank .NE. 0)then
        Call MPI_ISEND(temp, Mx*N, MPI_REAL8, 0, rank,
     .                MPI_COMM_WORLD, req(1), ierr)
       else

        F_diag(1:Mx,1:N) = temp

        Do M = 1,Px-1
         L = (Px*Px + Px + 1)*M
         Call MPI_RECV(F_aux, Mx*N, MPI_REAL8, L, L,
     .                MPI_COMM_WORLD, status, ierr)

         F_diag(M*Mx+1:(M+1)*Mx,1:N) = F_aux
        End Do
       endif

      endif

* Write to output
      if(rank .EQ. 0)then
       Do I = 1,Nx
        OPEN(UNIT=stats3_unit,FILE="d_trace",FORM="FORMATTED",
     .       STATUS="UNKNOWN", POSITION="APPEND")
        WRITE(stats3_unit,8) F_diag(I,1), F_diag(I,2), F_diag(I,3),
     .  F_diag(I,4)
        CLOSE(stats3_unit)
       End Do
      endif

 8    Format(1x,4(ES13.5)) 

      End Subroutine Print_Diagonal

*******************************************************************
* Subroutine : Print_Line_X                                       *
*                                                                 *
* Function : Writes a line of an array along x                    *
*            Assumes Py & Pz >= 2                                 *
*******************************************************************
      Subroutine Print_Line_X(F,N)
      Include 'header'
      Include 'mpif.h'

      Integer I, J, K, L, M, N

      Real*8 F(My,Mz,Mx,N), F_aux(My,Mz,Mx,N)

      Integer status(MPI_STATUS_SIZE), ierr, req(2)

      if((yrank.EQ.0).AND.(zrank.EQ.0))then
        Call MPI_ISEND(F(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, 0, xrank,
     .                 MPI_X_COMM, req(1), ierr)

        if(xrank.EQ.0)then
          OPEN(UNIT=16,FILE="x_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Px-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_X_COMM, status, ierr)
            Do I = 1, Mx
             write(16,7) F_aux(1,1,I,1), F_aux(1,1,I,2),
     .                   F_aux(1,1,I,3), F_aux(1,1,I,4)
            End Do
          End Do
          CLOSE(16)
        endif
      endif

      Call MPI_BARRIER(MPI_X_COMM, ierr)

 7    Format(1x, 4(ES13.5))

      End Subroutine Print_Line_X

*******************************************************************
* Subroutine : Print_Line_Y                                       *
*                                                                 *
* Function : Writes a line of an array along y                    *
*            Assumes Px & Pz >= 2                                 *
*******************************************************************
      Subroutine Print_Line_Y(F,N)
      
      Include 'header'
      Include 'mpif.h'
      
      Integer I, J, K, L, M, N
      
      Real*8 F(My,Mz,Mx,N), F_aux(My,Mz,Mx,N)

      Integer status(MPI_STATUS_SIZE), ierr, req(2)
      
      if((xrank.EQ.0).AND.(zrank.EQ.0))then
        Call MPI_ISEND(F(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, 0, yrank,
     .                 MPI_Y_COMM, req(1), ierr)
     
        if(yrank .EQ. 0)then
          OPEN(UNIT=17,FILE="y_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Py-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_Y_COMM, status, ierr)
            Do J = 1, My 
             write(17,7) F_aux(J,1,1,1), F_aux(J,1,1,2),
     .                        F_aux(J,1,1,3), F_aux(J,1,1,4)
            End Do            
          End Do
          CLOSE(17)
        endif
      endif
      
      Call MPI_BARRIER(MPI_Y_COMM, ierr)
      
 7    Format(1x, 4(ES13.5))
 
      End Subroutine Print_Line_Y
      
*******************************************************************
* Subroutine : Print_Line_Z                                       *
*                                                                 *
* Function : Writes a line of an array along z                    *
*            Assumes Px & Py >= 2                                 *
*******************************************************************
      Subroutine Print_Line_Z(F,N)

      Include 'header'
      Include 'mpif.h'

      Integer I, J, K, L, M, N

      Real*8 F(My,Mz,Mx,N), F_aux(My,Mz,Mx,N)

      Integer status(MPI_STATUS_SIZE), ierr, req(2)

      if((xrank.EQ.0).AND.(yrank.EQ.0))then
        Call MPI_ISEND(F(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, 0, zrank,
     .                 MPI_Z_COMM, req(1), ierr)

        if(zrank .EQ. 0)then
          OPEN(UNIT=18,FILE="z_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Pz-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_Z_COMM, status, ierr)
            Do K = 1, Mz
             write(18,7) F_aux(1,K,1,1), F_aux(1,K,1,2),
     .                            F_aux(1,K,1,3), F_aux(1,K,1,4)
            End Do
          End Do
          CLOSE(18)
        endif
      endif

      Call MPI_BARRIER(MPI_Z_COMM, ierr)

 7    Format(1x, 4(ES13.5))

      End Subroutine Print_Line_Z


