**********************************************************
* Post-processing stuff                                  *
**********************************************************
      Subroutine Compute_SRM_Stats(step)
 
      Include 'header'
      Include 'mpif.h'

      Integer step      

      Real*8 Q(My,Mz,Mx,5)

      Common /RKWORK/ Q

      Real*8 Y(My,Mz,Mx,1)

      Common /SPECIES/ Y     

      Real*8 gamma_eff(My,Mz,Mx), gradrho(My,Mz,Mx,3)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 S(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3), Ur(My,Mz,Mx), rad(My,Mz,Mx) 
      Real*8 buf(My,Mz,Mx,7)
      Real*8 fac1, fac2, tmp
 
      Common /THERMO/ gamma_eff

      Real*8 r(N_avg+1)
      Real*8 Ur_shell_avg(N_avg), Y_shell_avg(N_avg)
      Real*8 urp2(My,Mz,Mx), utp2(My,Mz,Mx), tke(My,Mz,Mx)
      Real*8 tke_total, hmix

      Integer I, J, K, M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

      Real*8 beta_hyper(My,Mz,Mx),visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx), D_hyper(My,Mz,Mx)
      
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper, D_hyper

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
        End Do
       End Do
      End Do

! Compute shell averages
      r(1) = 0D0

      Do M = 2,N_avg+1
        r(M) = r(M-1) + DX * sqrt(3D0)
      End Do

      Call Compute_Shell_Avg(Ur, Ur_shell_avg)
      Call Compute_Shell_Avg(Y, Y_shell_avg)

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

               tke(J,K,I)      = utp2(J,K,I) + urp2(J,K,I)
              endif
            End Do
          End Do
        End Do
      End Do

! Compute tke and mixing zone width
      call Compute_Sum(tke, tke_total)
      tke_total = tke_total / DBLE(Nx*Ny*Nz)
      hmix = 0D0
      Do M = 1,N_avg
        if((Y_shell_avg(M).GT.0D0).AND.(Y_shell_avg(M).LT.1D0))then
           hmix = hmix + 8D0 * Y_shell_avg(M) * (1D0 - Y_shell_avg(M)) 
     .                 * DX * sqrt(3D0)
        endif
      End Do

      if(rank .EQ. 0)then
        OPEN(UNIT=7,FILE="mix_stats",FORM="FORMATTED",
     .       STATUS="UNKNOWN", POSITION="APPEND")
        WRITE(7,7) hmix, tke_total
        CLOSE(7)
      endif

c      call Compute_Gradient(Q(1,1,1,4), gradrho(1,1,1,1))
c      buf(:,:,:,1) = gradrho(:,:,:,1) !dsqrt(gradrho(:,:,:,1)**2D0 + 
c     .  gradrho(:,:,:,2)**2D0 + gradrho(:,:,:,3)**2D0)       

c       Call Print_Line_X(buf,1)

 6    Format(1x, 7(ES13.5))
 7    Format(1x, 2(ES13.5))

      End Subroutine Compute_SRM_Stats

***********************************************************************

      Subroutine Fix_Pressure(Q,Y)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), Y(My,Mz,Mx,1)
      Real*8 gamma_eff(My,Mz,Mx) 
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 fac1, fac2, fac3, base_pressure

      Common /THERMO/ gamma_eff

      Integer step
      Integer I,J,K

      base_pressure = 0.71429D0 

* Pressure fix
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac1 = gamma_eff(J,K,I) - 1.0D0
            fac2 = gamma_eff(J,K,I) / fac1
            fac3 = gamma1 / (gamma1 - 1D0)
            fac3 = fac3 *
     .      1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

            Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
             P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) ) )

             if(P(J,K,I) .LT. 0D0) P(J,K,I) = base_pressure

             Q(J,K,I,5) = (1D0/fac1)*P(J,K,I) + 5.0D-1*Vs(J,K,I)*
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) ) 

          End Do
        End Do
      End Do

      End Subroutine Fix_Pressure

********************************************************************************
      Subroutine Fix_Mass_Fraction(Y)

      Include 'header'

      Real*8 Y(My,Mz,Mx,1)

      Integer I, J, K, II, JJ, KK

*  Mass Fraction fix
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             II = xrank*Mx + I
             JJ = yrank*My + J
             KK = zrank*Mz + K
             if (Y(J,K,I,1).LT.0D0) Y(J,K,I,1) = 0D0 
c             if (Y(J,K,I,1).GT.1D0) Y(J,K,I,1) = 1D0 
          End Do
        End Do
      End Do

      End Subroutine Fix_Mass_Fraction

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
         Call MPI_IRECV(F_aux, Mx*N, MPI_REAL8, L, L,
     .                MPI_COMM_WORLD, req(2), ierr)

*        Make sure all communication is complete
         Call MPI_WAIT(req(2), status, ierr)

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
     .  F_diag(I,4),F_diag(I,5), F_diag(I,6), F_diag(I,7)
        CLOSE(stats3_unit)
       End Do
      endif

 8    Format(1x,7(ES13.5)) 

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
          OPEN(UNIT=16,FILE="x_trace_conv",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Px-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_X_COMM, status, ierr)
            Do I = 1, Mx
             write(16,7) F_aux(My,Mz,I,1)!, F_aux(My,Mz,I,2),
c     .                        F_aux(My,Mz,I,3), F_aux(My,Mz,I,4)
            End Do
          End Do
          CLOSE(16)
        endif
      endif

      Call MPI_BARRIER(MPI_X_COMM, ierr)

 7    Format(1x, 1(ES13.5))

      End Subroutine Print_Line_X

