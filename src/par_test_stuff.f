* Testing All-to-All

      Do I = 1,8
        temp_send(I) = 8D0*REAL(xrank) + REAL(I)
      End Do

      temp_recv = 0D0

      Call MPI_AlltoAll(temp_send, 1, MPI_REAL8, temp_recv, 1,
     .                  MPI_REAL8, MPI_X_COMM, ierr)

      Do I = 1,8
       print *, xrank, temp_send(I), temp_recv(I)
      End Do

* Testing derivative routines

      Call DFDX(Q(1,1,1,1), Q(1,1,1,2), 0, 1.0D0)

      Do I = 1,Mx
        Do J = 1,Ny
          Do K = 1,Nz
            Q(J,K,I,2) = Q(J,K,I,2) - cos(xcor(I))
          End Do
        End Do
      End Do

      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, xrank,
     .               MPI_X_COMM, req(1), ierr)

      if(xrank .EQ. 0)then
        Do M = 0,(Px-1)
          Call MPI_RECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M,
     .                 MPI_X_COMM, status, ierr)
          Do I = 1, Mx
           print *, Q_aux(1,1,I,1), Q_aux(1,1,I,2)
          End Do
        End Do
      endif

      Call MPI_BARRIER(MPI_X_COMM, ierr)

* Testing filters
      Call MPI_ISEND(Qf(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, yrank,
     .               MPI_Y_COMM, req(1), ierr)
      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, yrank+4,
     .               MPI_Y_COMM, req(1), ierr)

      if(yrank .EQ. 0)then
        Do M = 0,(Py-1)
          Call MPI_RECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M,
     .                 MPI_Y_COMM, status, ierr)
          Call MPI_RECV(QW(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M+4,
     .                 MPI_Y_COMM, status, ierr)

          Do J = 1, My
           if((xrank .EQ. 0) .AND. (zrank .EQ. 0))then
            print *, QW(J,1,1,4), Q_aux(J,1,1,4)
           end if
          End Do
        End Do
      endif

      Call MPI_BARRIER(MPI_Y_COMM, ierr)


      Call MPI_ISEND(Qf(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, zrank,
     .               MPI_Z_COMM, req(1), ierr)
      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, zrank+4,
     .               MPI_Z_COMM, req(1), ierr)


      if(zrank .EQ. 0)then
        Do M = 0,(Pz-1)
          Call MPI_RECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M,
     .                 MPI_Z_COMM, status, ierr)
          Call MPI_RECV(QW(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M+4,
     .                 MPI_Z_COMM, status, ierr)

          Do K = 1, Mz
           if((xrank .EQ. 0) .AND. (yrank .EQ. 0))then
            print *, QW(1,K,1,4), Q_aux(1,K,1,4)
           end if
          End Do
        End Do
      endif

      Call MPI_BARRIER(MPI_Z_COMM, ierr)



      Call MPI_ISEND(Qf(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, xrank,
     .               MPI_X_COMM, req(1), ierr)
      Call MPI_ISEND(Q(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, 0, xrank+4,
     .               MPI_X_COMM, req(1), ierr)


      if(xrank .EQ. 0)then
        Do M = 0,(Px-1)
          Call MPI_RECV(Q_aux(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M,
     .                 MPI_X_COMM, status, ierr)
          Call MPI_RECV(QW(1,1,1,1), My*Mz*Mx*5, MPI_REAL8, M, M+4,
     .                 MPI_X_COMM, status, ierr)

          Do I = 1, Mx
           if((yrank .EQ. 0) .AND. (zrank .EQ. 0))then
            print *, QW(1,1,I,4), Q_aux(1,1,I,4)
           end if
          End Do
        End Do
      endif

      Call MPI_BARRIER(MPI_X_COMM, ierr)





 
