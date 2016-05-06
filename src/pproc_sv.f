************************************************************************
*  Subroutine: Write_Yavg                                              *
*                                                                      *
*  Function:   It writes out the y-averaged flow field data            *
*              which is needed for post-processing.                    *
*                                                                      *
*              The data are recorded in a total of Px sections in      *
*              descending order of processor rank.                     *
************************************************************************
      Subroutine Write_Yavg(time, step)

      Include 'header'
      Include 'mpif.h'

      Integer step 
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)

      Common  /RKWORK/ Q, Grad

      Real*8  QW(Mz,Mx,2), QW_aux(Mz,Mx,2)

      Common  QW, QW_aux

      Real*8  Q_aux(My,Mz,Mx,5), Grad_aux(My,Mz,Mx,4)
      Real*8  kin_e(My,Mz,Mx),     vortz2(My,Mz,Mx)
      Real*8  kin_e_aux(My,Mz,Mx), vortz2_aux(My,Mz,Mx)
      Real*8  U(My,Mz,Mx,3), U_aux(My,Mz,Mx,2)

*  Locals
      Integer      I, J, K, L, M, II, JJ
      Integer      status(MPI_STATUS_SIZE), ierr, req(3)

* Calculating fluctuation kinetic energy and mean square z-vorticity
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            II = xrank * Mx + I
            if(II .LT. 75)then
              U(J,K,I,1) = (Q(J,K,I,1)/Q(J,K,I,4)) - 1.5D0
            else if(II .GE. 75)then
              U(J,K,I,1) = (Q(J,K,I,1)/Q(J,K,I,4)) - (1.5D0
     .                                     * 0.537037D0)
            endif

            U(J,K,I,2) = Q(J,K,I,2)/Q(J,K,I,4)
            U(J,K,I,3) = Q(J,K,I,3)/Q(J,K,I,4)

            kin_e(J,K,I)    =  U(J,K,I,1) * U(J,K,I,1)
     .                      +  U(J,K,I,2) * U(J,K,I,2)
            vortz2(J,K,I) = 5.0D-1 *
     .                       (Grad(J,K,I,3) * Grad(J,K,I,3))
          End Do
        End Do
      End Do

*  Send data to rank 0 processor

      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0,
     .                             xrank, MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(Grad(1,1,1,1), My*Mz*Mx*4, MPI_REAL8, 0,
     .                             xrank, MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(vortz2(1,1,1), My*Mz*Mx, MPI_REAL8, 0,
     .                             xrank, MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(kin_e(1,1,1),  My*Mz*Mx, MPI_REAL8, 0,
     .                             xrank, MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(U(1,1,1,1),    My*Mz*Mx*2, MPI_REAL8, 0,
     .                             xrank, MPI_X_COMM, req(1), ierr)

*  Rank 0 processor writes the data to file
      If (rank .EQ. 0) Then
        Do M = 0, Px-1

          Call MPI_RECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M,
     .                  M, MPI_X_COMM, status, ierr)
          Call MPI_RECV(Grad_aux(1,1,1,1), My*Mz*Mx*4, MPI_REAL8, M,
     .                  M, MPI_X_COMM, status, ierr)
          Call MPI_RECV(vortz2_aux(1,1,1), My*Mz*Mx, MPI_REAL8, M,
     .                  M, MPI_X_COMM, status, ierr)
          Call MPI_RECV(kin_e_aux(1,1,1), My*Mz*Mx, MPI_REAL8, M,
     .                  M, MPI_X_COMM, status, ierr)
          Call MPI_RECV(U_aux(1,1,1,1), My*Mz*Mx*2, MPI_REAL8, M,
     .                  M, MPI_X_COMM, status, ierr)

c          Do I = 1, Mx
c            II = M * Mx + I
c            Do J = 1,Ny
c             Write(stats1_unit,52) time, gxcor(II),
c     .       vortz2_aux(J,1,I) , kin_e_aux(J,1,I),
c     .       U_aux(J,1,I,1),U_aux(J,1,I,2)
c             Write(stats2_unit,53) time,gxcor(II),
c     .       Q_aux(J,1,I,4), Grad_aux(J,1,I,3)
c            End Do
c          End Do
        End Do

 52     Format(x,6(ES13.6))
 53     Format(x,6(ES13.6))

      End If

      Call MPI_BARRIER(MPI_X_COMM, ierr)

      Return

      End

