**********************************************************
* Debug stuff                                            *
**********************************************************
      Subroutine Compute_Shockturb_Stats   
 
      Include 'header'

      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)

      Common  /RKWORK/ Q, Grad

      Real*8 beta_hyper(My,Mz,Mx), visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)
      
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3), buf(My,Mz,Mx,6)
      Real*8 Ur(My,Mz,Mx), omr(My,Mz,Mx), omt(My,Mz,Mx)
      Real*8 omrp2(My,Mz,Mx), omtp2(My,Mz,Mx) 
      Real*8 fac1, fac2, tmp, tmp1, tmp2, tmp3
      Real*8 rad(My,Mz,Mx)
      Real*8 temp(My,Mz,Mx,2), gradrho(My,Mz,Mx,3)

      Real*8 urp2(My,Mz,Mx), utp2(My,Mz,Mx), tke(My,Mz,Mx)
      Real*8 urp(My,Mz,Mx), utp(My,Mz,Mx)
      Real*8 durp_dr2(My,Mz,Mx), dutp_dr2(My,Mz,Mx)
      Real*8 Pp2(My,Mz,Mx), rhop2(My,Mz,Mx), pvort(My,Mz,Mx)
 
      Real*8 omr_shell_avg(N_avg), omt_shell_avg(N_avg)
      Real*8 omrp2_shell_avg(N_avg), omtp2_shell_avg(N_avg)
      Real*8 Dil_shell_avg(N_avg)
      Real*8 Ur_shell_avg(N_avg), Ut_shell_avg(N_avg)
      Real*8 rho_shell_avg(N_avg)
      Real*8 urp2_shell_avg(N_avg), utp2_shell_avg(N_avg)
      Real*8 P_shell_avg(N_avg), r(N_avg+1), temp_shell(N_avg,2)

      Real*8 lambda_ur_r(N_avg), lambda_ut_r(N_avg)
      Real*8 tke_shell_avg(N_avg), Pp2_shell_avg(N_avg)
      Real*8 rhop2_shell_avg(N_avg)

      Real*8 beta_shell_avg(N_avg), visc_shell_avg(N_avg)
      Real*8 k_shell_avg(N_avg)


      Integer I, J, K, M
      Integer req(4)
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

* Setting bin levels
* There are N/4 + 4 levels, spaced Delta * sqrt(3) apart in r-space

      r(1) = DX * sqrt(3D0) / 2D0 - 1d-4

      Do M = 2,N_avg+1
        if(M .LE. 4)then
          r(M) = r(M-1) + DX * 2d0 * sqrt(3D0)
        else
          r(M) = r(M-1) + DX * sqrt(3D0)
        endif
      End Do

*  Compute specific volume (avoid division in computations), pressure,
*  temperature and streamwise velocity
*
*  Vs = 1 / rho
*
*  P  = (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*         gamma     P
*  T  = ---------  ---
*       gamma - 1  rho
*

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
             P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) ) )
             T(J,K,I) = fac2 * Vs(J,K,I)  * P(J,K,I)
             k_hyper(J,K,I) = ((fac1 * T(J,K,I)) ** vis_exp ) * Re_inv
          End Do
        End Do
      End Do

*  Convert momenta to velocities
      Do I = 1,Mx
       Do J = 1,My
         Do K = 1,Mz
          II = xrank * Mx + I
          JJ = yrank * My + J
          KK = zrank * Mz + K

          IX = DBLE(II - Nx/2) - 0.5D0
          IY = DBLE(JJ - Ny/2) - 0.5D0
          IZ = DBLE(KK - Nz/2) - 0.5D0

          rad(J,K,I) = 
     .    dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

          tmp = 1.0D0 / Q(J,K,I,4)
          U(J,K,I,1) = Q(J,K,I,1) * tmp
          U(J,K,I,2) = Q(J,K,I,2) * tmp
          U(J,K,I,3) = Q(J,K,I,3) * tmp

          Ur(J,K,I)  = U(J,K,I,1) * (IX*DX/rad(J,K,I))
     .               + U(J,K,I,2) * (IY*DY/rad(J,K,I))
     .               + U(J,K,I,3) * (IZ*DZ/rad(J,K,I))
          
          omr(J,K,I) = Grad(J,K,I,1) * (IX*DX/rad(J,K,I))
     .               + Grad(J,K,I,2) * (IY*DY/rad(J,K,I))
     .               + Grad(J,K,I,3) * (IZ*DZ/rad(J,K,I))
          pvort(J,K,I) = sqrt(Grad(J,K,I,1)*Grad(J,K,I,1)
     .                 +      Grad(J,K,I,2)*Grad(J,K,I,2)
     .                 +      Grad(J,K,I,3)*Grad(J,K,I,3))/
     .                        Q(J,K,I,4)
        End Do
       End Do
      End Do

      Call Compute_Shell_Avg(Q(1,1,1,4), rho_shell_avg)
      Call Compute_Shell_Avg(P, P_shell_avg)
      Call Compute_Shell_Avg(Ur, Ur_shell_avg)
      Call Compute_Shell_Avg(Grad(1,1,1,4),Dil_shell_avg)
      Call Compute_Shell_Avg(omr, omr_shell_avg)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Do M = 1,N_avg
              if((rad(J,K,I).GE.r(M)).AND.(rad(J,K,I).LT.r(M+1)))then
               II = xrank * Mx + I
               JJ = yrank * My + J
               KK = zrank * Mz + K

               IX = DBLE(II - Nx/2) - 0.5D0
               IY = DBLE(JJ - Ny/2) - 0.5D0
               IZ = DBLE(KK - Nz/2) - 0.5D0

               urp2(J,K,I)     = (Ur(J,K,I) - Ur_shell_avg(M))**2d0
               urp(J,K,I)      = Ur(J,K,I) - Ur_shell_avg(M)
               tmp1            = U(J,K,I,1) - 
     .                           Ur(J,K,I) * (IX*DX/rad(J,K,I))
               tmp2            = U(J,K,I,2) - 
     .                           Ur(J,K,I) * (IY*DY/rad(J,K,I))
               tmp3            = U(J,K,I,3) - 
     .                           Ur(J,K,I) * (IZ*DZ/rad(J,K,I))
               utp2(J,K,I)     = tmp1**2d0 + tmp2**2d0 + tmp3**2d0
               utp(J,K,I)      = dsqrt(utp2(J,K,I))

               omrp2(J,K,I)    = (omr(J,K,I) - omr_shell_avg(M))**2d0
               tmp1            = Grad(J,K,I,1) -
     .                           omr(J,K,I) * (IX*DX/rad(J,K,I))
               tmp2            = Grad(J,K,I,2) -
     .                           omr(J,K,I) * (IY*DY/rad(J,K,I))
               tmp3            = Grad(J,K,I,3) -
     .                           omr(J,K,I) * (IZ*DZ/rad(J,K,I))
               omtp2(J,K,I)    = tmp1**2d0 + tmp2**2d0 + tmp3**2d0

               tke(J,K,I)      = utp2(J,K,I) + urp2(J,K,I)
              endif
            End Do
          End Do
        End Do
      End Do          

      Call Compute_Shell_RMS(Ur, urp2_shell_avg)
      urp2_shell_avg = urp2_shell_avg**2d0
      Call Compute_Shell_RMS(omr, omrp2_shell_avg)
      omrp2_shell_avg = omrp2_shell_avg**2d0
      Call Compute_Shell_RMS(P, Pp2_shell_avg)
      Pp2_shell_avg = Pp2_shell_avg**2d0
      Call Compute_Shell_RMS(Q(1,1,1,4), rhop2_shell_avg)
      rhop2_shell_avg = rhop2_shell_avg**2d0

      Call Compute_Shell_Avg(utp2, utp2_shell_avg)
      Call Compute_Shell_Avg(omtp2, omtp2_shell_avg)

      Call DFDX(urp, buf(1,1,1,1), 0, 1.0D0)
      Call DFDY(urp, buf(1,1,1,2), 0, 1.0D0)
      Call DFDZ(urp, buf(1,1,1,3), 0, 1.0D0)

      Call DFDX(utp, buf(1,1,1,4), 0, 1.0D0)
      Call DFDY(utp, buf(1,1,1,5), 0, 1.0D0)
      Call DFDZ(utp, buf(1,1,1,6), 0, 1.0D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            II = xrank * Mx + I
            JJ = yrank * My + J
            KK = zrank * Mz + K

            IX = DBLE(II - Nx/2) - 0.5D0
            IY = DBLE(JJ - Ny/2) - 0.5D0
            IZ = DBLE(KK - Nz/2) - 0.5D0

            durp_dr2(J,K,I) = ((IX*DX/rad(J,K,I))*buf(J,K,I,1)
     .                      + (IY*DY/rad(J,K,I))*buf(J,K,I,2)
     .                      + (IZ*DZ/rad(J,K,I))*buf(J,K,I,3))**2D0
            dutp_dr2(J,K,I) = ((IX*DX/rad(J,K,I))*buf(J,K,I,4)
     .                      + (IY*DY/rad(J,K,I))*buf(J,K,I,5)
     .                      + (IZ*DZ/rad(J,K,I))*buf(J,K,I,6))**2D0
          End Do
        End Do
      End Do

      Call Compute_Shell_Avg(durp_dr2, lambda_ur_r)
      Call Compute_Shell_Avg(dutp_dr2, lambda_ut_r)

      Do M = 1,N_avg
        lambda_ur_r(M) = dsqrt(urp2_shell_avg(M)/lambda_ur_r(M))
        lambda_ut_r(M) = dsqrt(utp2_shell_avg(M)/lambda_ut_r(M))
      End Do

      Call Compute_Shell_Avg(beta_hyper, beta_shell_avg)
      Call Compute_Shell_Avg(visc_hyper, visc_shell_avg) 
      Call Compute_Shell_Avg(k_hyper, k_shell_avg)

c      if(rank.EQ.0)then
c          OPEN(UNIT=5,FILE="avg_stats",FORM="FORMATTED",
c     .         STATUS="UNKNOWN", POSITION="APPEND")
c          OPEN(UNIT=6,FILE="vel_stats",FORM="FORMATTED",
c     .         STATUS="UNKNOWN", POSITION="APPEND")
c          OPEN(UNIT=7,FILE="turb_stats",FORM="FORMATTED",
c     .         STATUS="UNKNOWN", POSITION="APPEND")
c          OPEN(UNIT=8,FILE="hyper_stats",FORM="FORMATTED",
c     .         STATUS="UNKNOWN", POSITION="APPEND")
c        Do M = 1, N_avg
c          WRITE(5,*) rho_shell_avg(M), P_shell_avg(M),
c     .               Ur_shell_avg(M) 
c          WRITE(6,*) urp2_shell_avg(M),  utp2_shell_avg(M),
c     .               omrp2_shell_avg(M), omtp2_shell_avg(M),
c     .               Pp2_shell_avg(M),   rhop2_shell_avg(M)
c          WRITE(7,*) lambda_ur_r(M), lambda_ut_r(M) 
c          WRITE(8,*) beta_shell_avg(M), visc_shell_avg(M),
c     .               k_shell_avg(M)
c        End Do
c        CLOSE(5); CLOSE(6); CLOSE(7); CLOSE(8)
c      endif

      Call Compute_Gradient(Q(1,1,1,4), gradrho)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            buf(J,K,I,1) = dsqrt(gradrho(J,K,I,1)**2D0
     .                         + gradrho(J,K,I,2)**2D0
     .                         + gradrho(J,K,I,3)**2D0)
          End Do
        End Do
      End Do
     
c      buf(:,:,:,1) = Q(:,:,:,4)
      buf(:,:,:,2) = P
      buf(:,:,:,3) = Ur
      buf(:,:,:,4) = Grad(:,:,:,4) 

      Call Print_Diagonal(buf(1,1,1,1),1)
      Call Print_Line_X(buf(1,1,1,1),1)
      Call Print_Line_Y(buf(1,1,1,1),1)
      Call Print_Line_Z(buf(1,1,1,1),1)

      Return
     
 5    Format(1x, 3(ES13.5))
 6    Format(1x, 6(ES13.5))
 7    Format(1x, 2(ES13.5))
 8    Format(1x, 3(ES13.5))

      End Subroutine Compute_Shockturb_Stats 

*******************************************************************
* Subroutine : Print_Diagonal                                     *
*                                                                 *
* Function : Writes the diagonal of an array                      *
*            Assumes Mx = My = Mz, Px = Py = Pz                   *
*******************************************************************
      Subroutine Print_Diagonal(F,N)

      Include 'header'
      Include 'mpif.h'

      Integer I, J, K, L, M, N

      Real*8 F(My,Mz,Mx,N), F_aux(Mx,N)
      Real*8 F_diag(Nx,N), temp(Mx,N)

      Integer status(MPI_STATUS_SIZE), ierr, req(2)

      L = 1
       
* Put all the diagonal elements of all diagonal processes into temporary array
      if((xrank.EQ.yrank) .AND. (xrank.EQ.zrank))then

       Do I = 1,Mx
        Do J = 1,My
         Do K = 1,Mz
           if((I .EQ. J) .AND. (I .EQ. K))then
            Do M = 1,N
              temp(L,M) = F(J,K,I,M)
            End Do
            L = L + 1
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
       OPEN(UNIT=15,FILE="diag_trace",FORM="FORMATTED",
     .      STATUS="UNKNOWN", POSITION="APPEND")
       Do I = 1,Nx
        write(15,7) F_diag(I,1)!, F_diag(I,2), F_diag(I,3),
c     .                       F_diag(I,4) 
       End Do
       CLOSE(15)
      endif

 7    Format(1x, 1(ES13.5))

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

      if((yrank.EQ.(Py/2-1)).AND.(zrank.EQ.(Pz/2-1)))then
        Call MPI_ISEND(F(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, 0, xrank,
     .                 MPI_X_COMM, req(1), ierr)

        if(xrank .EQ. 0)then
          OPEN(UNIT=16,FILE="x_trace",FORM="FORMATTED",
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

      if((xrank.EQ.(Px/2-1)).AND.(zrank.EQ.(Pz/2-1)))then
        Call MPI_ISEND(F(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, 0, yrank,
     .                 MPI_Y_COMM, req(1), ierr)

        if(yrank .EQ. 0)then
          OPEN(UNIT=17,FILE="y_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Py-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_Y_COMM, status, ierr)
            Do J = 1, My
             write(17,7) F_aux(J,Mz,Mx,1)!, F_aux(J,Mz,Mx,2),
c     .                        F_aux(J,Mz,Mx,3), F_aux(J,Mz,Mx,4)
            End Do
          End Do
          CLOSE(17)
        endif
      endif

      Call MPI_BARRIER(MPI_Y_COMM, ierr)

 7    Format(1x, 1(ES13.5))

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

      if((xrank.EQ.(Px/2-1)).AND.(yrank.EQ.(Py/2-1)))then
        Call MPI_ISEND(F(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, 0, zrank,
     .                 MPI_Z_COMM, req(1), ierr)

        if(zrank .EQ. 0)then
          OPEN(UNIT=18,FILE="z_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Pz-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_Z_COMM, status, ierr)
            Do K = 1, Mz
             write(18,7) F_aux(My,K,Mx,1)!, F_aux(My,K,Mx,2),
c     .                            F_aux(My,K,Mx,3), F_aux(My,K,Mx,4)
            End Do
          End Do
          CLOSE(18)
        endif
      endif

      Call MPI_BARRIER(MPI_Z_COMM, ierr)

 7    Format(1x, 1(ES13.5))

      End Subroutine Print_Line_Z






