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

        if(xrank .EQ. 0)then
          OPEN(UNIT=16,FILE="x_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do M = 0,(Px-1)
            Call MPI_RECV(F_aux(1,1,1,1), My*Mz*Mx*N, MPI_REAL8, M, M,
     .                   MPI_X_COMM, status, ierr)
            Do I = 1, Mx
             write(16,7) F_aux(1,1,I,1), F_aux(1,1,I,2),
     .                   F_aux(1,1,I,3), F_aux(1,1,I,4),
     .                   F_aux(1,1,I,5)
            End Do
          End Do
          CLOSE(16)
        endif
      endif

      Call MPI_BARRIER(MPI_X_COMM, ierr)

 7    Format(1x, 5(ES13.5))

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
     .                   F_aux(J,1,1,3), F_aux(J,1,1,4),
     .                   F_aux(J,1,1,5)
            End Do
          End Do
          CLOSE(17)
        endif
      endif

      Call MPI_BARRIER(MPI_Y_COMM, ierr)

 7    Format(1x, 5(ES13.5))

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
     .                   F_aux(1,K,1,3), F_aux(1,K,1,4),
     .                   F_aux(1,K,1,5)
            End Do
          End Do
          CLOSE(18)
        endif
      endif

      Call MPI_BARRIER(MPI_Z_COMM, ierr)

 7    Format(1x, 5(ES13.5))

      End Subroutine Print_Line_Z

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
        write(15,7) F_diag(I,1), F_diag(I,2), F_diag(I,3),
     .              F_diag(I,4), F_diag(I,5)
       End Do
       CLOSE(15)
      endif

 7    Format(1x, 5(ES13.5))

      End Subroutine Print_Diagonal

