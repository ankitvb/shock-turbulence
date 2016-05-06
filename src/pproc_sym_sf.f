**********************************************************
* Post-processing stuff                                  *
**********************************************************
      Subroutine Debug(Q, step)   
 
      Include 'header'
      Include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 S(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3), Ur(My,Mz,Mx), rad(My,Mz,Mx) 
      Real*8 buf(My,Mz,Mx,7)
      Real*8 fac1, fac2, tmp
      Real*8 time
 
      Integer step
      Integer I, J, K, M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

      Real*8 beta_hyper(My,Mz,Mx),visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)
      
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

      Real*8 Q_aux(My,Mz,Mx,5), P_aux(My,Mz,Mx)
      Real*8 beta_aux(My,Mz,Mx), k_aux(My,Mz,Mx)
      Integer status_arr(MPI_STATUS_SIZE,6), ierr, req(6), rreq(6)
      Integer status(MPI_STATUS_SIZE)
      Logical flag

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

          IX = REAL(II-1) 
          IY = REAL(JJ-1)
          IZ = REAL(KK-1)

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

      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, xrank,
     .               MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(P, My*Mz*Mx, MPI_REAL8, 0, xrank+14,
     .               MPI_X_COMM, req(2), ierr)
      Call MPI_ISEND(beta_hyper, My*Mz*Mx, MPI_REAL8, 0, xrank+42,
     .               MPI_X_COMM, req(3), ierr)
      Call MPI_ISEND(k_hyper, My*Mz*Mx, MPI_REAL8, 0, xrank+56,
     .               MPI_X_COMM, req(4), ierr)

      if(xrank .EQ. 0)then
        Do M = 0,(Px-1)
          Call MPI_IRECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M,
     .                 MPI_X_COMM, rreq(1), ierr)
          Call MPI_IRECV(P_aux, My*Mz*Mx, MPI_REAL8, M, M+14,
     .                 MPI_X_COMM, rreq(2), ierr)
          Call MPI_IRECV(beta_aux, My*Mz*Mx, MPI_REAL8, M, M+42,
     .                 MPI_X_COMM, rreq(3), ierr)
          Call MPI_IRECV(k_aux, My*Mz*Mx, MPI_REAL8, M, M+56,
     .                 MPI_X_COMM, rreq(4), ierr)

          Call MPI_WAITALL(4, rreq, status_arr, ierr)


          Do I = 1, Mx
           if((yrank .EQ. 0).AND.(zrank .EQ. 0))then 
             write(stats2_unit,6) Q_aux(1,1,I,4),
     .                         Q_aux(1,1,I,1)/Q_aux(1,1,I,4),
     .                         P_aux(1,1,I),
     .                         P_aux(1,1,I)/Q_aux(1,1,I,4)**gamma
           endif
          End Do
        End Do
      endif
 
      Call MPI_BARRIER(MPI_X_COMM, ierr)

      buf(:,:,:,1) = Q(:,:,:,4)
      buf(:,:,:,2) = Ur
      buf(:,:,:,3) = P
      buf(:,:,:,4) = P/Q(:,:,:,4)**gamma

      Call Print_Diagonal(buf, 5)

 6    Format(1x, 4(ES13.5))
 7    Format(1x, 3(ES16.8))

      End Subroutine Debug 

***********************************************************************

      Subroutine Fix_Pressure(Q)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 fac1, fac2, fac3, base_pressure

      Integer step
      Integer I,J,K

      base_pressure = 0.71429D0

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
        write(stats3_unit,8) F_diag(I,1), F_diag(I,2), F_diag(I,3),
     .                       F_diag(I,4)
       End Do
      endif

 8    Format(1x,4(ES13.5)) 

      End Subroutine Print_Diagonal





