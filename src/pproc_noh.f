**********************************************************
* Debug stuff                                            *
**********************************************************
      Subroutine Debug(Q, step)   
 
      Include 'header'
      Include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Qf(My,Mz,Mx,5)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 S(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3) 
      Real*8 Grad(My,Mz,Mx), buf(My,Mz,Mx,6)
      Real*8 fac1, fac2, tmp
      Real*8 time

      Integer step
      Integer I,J,K

      Real*8 beta_hyper(My,Mz,Mx),visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)
      
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

      Real*8 Q_aux(My,Mz,Mx,5), P_aux(My,Mz,Mx)
      Integer status(MPI_STATUS_SIZE), ierr, req(6)

      Integer II, JJ, KK
      Real*8  IX, IY, IZ

      Real*8  rad(My,Mz,Mx)
      Integer N, N_sum
      Real*8  rho_avg, rho_sum, rho2_avg, rho2_sum

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
             P(J,K,I) = fac1 * ( Q(J,K,I,5) - 5.0D-1 *Vs(J,K,I) *
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)))
             T(J,K,I) = fac2 * Vs(J,K,I)  * P(J,K,I)
             S(J,K,I) = (1/gamma)*log(P(J,K,I)) - log(Q(J,K,I,4))
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
        End Do
       End Do
      End Do

      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, xrank,
     .               MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(P, My*Mz*Mx, MPI_REAL8, 0, xrank+14,
     .               MPI_X_COMM, req(2), ierr)

      if(xrank .EQ. 0)then
        Do M = 0,(Px-1)
          Call MPI_RECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M,
     .                 MPI_X_COMM, status, ierr)
          Call MPI_RECV(P_aux, My*Mz*Mx, MPI_REAL8, M, M+14,
     .                 MPI_X_COMM, status, ierr)

          Do I = 1, Mx
           if((yrank .EQ. 0).AND.(zrank .EQ. 0))then
             write(stats1_unit,6) Q_aux(1,1,I,4),
     .                         Q_aux(1,1,I,1)/Q_aux(1,1,I,4),
     .                         P_aux(1,1,I)
           endif
          End Do
        End Do
      endif

      Call MPI_BARRIER(MPI_X_COMM, ierr)

      buf(:,:,:,1) = Q(:,:,:,4)
      buf(:,:,:,2) = Q(:,:,:,1)/Q(:,:,:,4)
      buf(:,:,:,3) = P

      Call Print_Diagonal(buf,3)

* Computing RMS and average density in the shell 0.15 < r < 0.17
      rho_avg  = 0D0
      rho2_avg = 0D0
      rho_sum  = 0D0
      rho2_sum = 0D0
      N_sum    = 0
      N        = 0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            II = xrank * Mx + I
            JJ = yrank * My + J
            KK = zrank * Mz + K

            IX = REAL(II-1)
            IY = REAL(JJ-1)
            IZ = REAL(KK-1)

            rad(J,K,I) = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

            if((rad(J,K,I).GE.0.15D0).AND.(rad(J,K,I).LE.0.17D0))then
              N_sum    = N_sum + 1
              rho_sum  = rho_sum  + Q(J,K,I,4)
              rho2_sum = rho2_sum + Q(J,K,I,4)*Q(J,K,I,4)
            endif
          End Do
        End Do
      End Do

      Call MPI_ALLREDUCE(rho_sum, rho_avg, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(rho2_sum, rho2_avg, 1, MPI_REAL8, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(N_sum, N, 1, MPI_INTEGER, MPI_SUM,
     .                   MPI_COMM_WORLD, ierr)

      rho_avg  = rho_avg / REAL(N)
      rho2_avg = rho2_avg / REAL(N)

      If(rank .EQ. 0)then
        Print *, step, rho_avg, rho2_avg
      End If

 6    Format(1x, 3(ES16.8))
 7    Format(1x, 3(ES16.8))

      End Subroutine Debug

***********************************************************************

      Subroutine Fix_Pressure(Q)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5) 
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 fac1, fac2, base_pressure

      Integer step
      Integer I,J,K

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      base_pressure = fac1 * 1.0D-5

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
       
* Pressure fix

             if(P(J,K,I) .LT. 0D0) P(J,K,I) = base_pressure

             Q(J,K,I,5) = (1D0/fac1)*P(J,K,I) + 5.0D-1*Vs(J,K,I)*
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) ) 

          End Do
        End Do
      End Do

      End Subroutine Fix_Pressure

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
       Do I = 1,Nx
        write(stats2_unit,7) F_diag(I,1), F_diag(I,2), F_diag(I,3)
       End Do
      endif

 7    Format(1x, 3(ES16.8))

      End Subroutine Print_Diagonal

