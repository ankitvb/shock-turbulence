************************************************************************
*  Subroutine: Write_Grid                                              *
*                                                                      *
*  Function:   It writes the grid coordinates to files.                *
************************************************************************
      Subroutine Write_Grid

      Include 'header'
      Include 'mpif.h'

      Integer step
      Real*8  time

      Real*8 xout(Mx,My,Mz), yout(Mx,My,Mz), zout(Mx,My,Mz)

*  Locals
      Character*10 index
      Character*25 filename
      Integer      I, J, K, L, M, N
      Real*8       fac1
      Integer      status(MPI_STATUS_SIZE), ierr
      Integer      req(2)
      Integer      gsizes(3), start_indices(3)
      Integer      lsizes(3), ndims
      Integer      count1, count2, count3
      Integer      fh, filetype, local_array_size
      Integer(kind=MPI_OFFSET_KIND)      disp
      Integer      intSize, realSize

      Integer      Ngrid
      Integer      nints, nreals
      Integer      sizeofRecord, sizeofBuffer
      Real*8       buf(4)

      Ngrid = 1            ! Single grid

      gsizes(1) = Nx
      gsizes(2) = Ny
      gsizes(3) = Nz

      count1 = 4
      count2 = 1
      count3 = 1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            xout(I,J,K) = xcor(I)
            yout(I,J,K) = ycor(J)
            zout(I,J,K) = zcor(K)
          End Do
        End Do
      End Do

! Generating file name
      filename = TRIM(ADJUSTL(vis_path)) // '/' //
     .           TRIM(ADJUSTL(prefix)) // '.grid.dat'

      Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,
     .                   MPI_MODE_CREATE + MPI_MODE_WRONLY,
     .                   MPI_INFO_NULL, fh, ierr)
      if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file'

      Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
      Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

      disp = 0

! Writing header info
! Number of grids
      nints = 1
      sizeofRecord = nints*intSize
      sizeofBuffer = sizeofRecord

      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, Ngrid, nints, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset
      disp = disp + sizeofBuffer

! Grid sizes
      nints = 3
      sizeofRecord = nints*intSize
      sizeofBuffer = sizeofRecord !+ 2*intSize

      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, gsizes, nints, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset
      disp = disp + sizeofBuffer

! Writing grids
      ndims = 3

      lsizes(1) = Mx
      lsizes(2) = My
      lsizes(3) = Mz

      start_indices(1) = xrank * Mx
      start_indices(2) = yrank * My
      start_indices(3) = zrank * Mz

      local_array_size = Mx * My * Mz 

      Call MPI_TYPE_CREATE_SUBARRAY(ndims,gsizes,lsizes,start_indices,
     .                        MPI_ORDER_FORTRAN, MPI_REAL8,
     .                        filetype, ierr)

      Call MPI_TYPE_COMMIT(filetype, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, filetype,
     .                      "native", MPI_INFO_NULL, ierr)
      If(ierr.NE.MPI_SUCCESS)Stop 'Error setting file view'

      Call MPI_FILE_WRITE_ALL(fh, xout, local_array_size,
     .                        MPI_REAL8, status, ierr)
      Call MPI_FILE_WRITE_ALL(fh, yout, local_array_size,
     .                        MPI_REAL8, status, ierr)
      Call MPI_FILE_WRITE_ALL(fh, zout, local_array_size,
     .                        MPI_REAL8, status, ierr)

      if(ierr.NE.MPI_SUCCESS)
     . stop 'Error during collective write of solution file.'

      Call MPI_FILE_CLOSE(fh, ierr)

      End Subroutine Write_Grid

************************************************************************
*  Subroutine: Write_Visual                                            *
*                                                                      *
*  Function:   It writes the values of conservative variables to file  *
*              q1 in the following sequence:                           *
*              [ rho u_x,  rho u_y,  rho u_z,  rho,  rho e_t]          *
*              and writes the values of vorticity, dilatation and      *
*              passive scalar to file q2 in the following sequence:    *
*              [ w_x,  w_y,  w_z,  dil,  P ]                           *
************************************************************************
      Subroutine Write_Visual(time, step)

      Include 'header'
      Include 'mpif.h'

      Integer step
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)
      Real*8  Qout(Mx,My,Mz,5)

      Common  /RKWORK/ Q, Grad

*  Locals
      Character*10 index
      Character*25 filename
      Integer      I, J, K, L, M, N
      Real*8       fac1
      Integer      status(MPI_STATUS_SIZE), ierr
      Integer      req(2)
      Integer      gsizes(4), start_indices(4)
      Integer      lsizes(4), ndims
      Integer      count1, count2, count3
      Integer      fh, filetype, local_array_size
      Integer(kind=MPI_OFFSET_KIND)      disp
      Integer      intSize, realSize

      Integer      Ngrid
      Integer      nints, nreals
      Integer      sizeofRecord, sizeofBuffer
      Real*8       buf(4)

      Ngrid = 1            ! Single grid

      gsizes(1) = Nx
      gsizes(2) = Ny
      gsizes(3) = Nz
      gsizes(4) = 5

      count1 = 4
      count2 = 1
      count3 = 1

! Rearranging data in I,J,K-orientation as required by plot3D format
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qout(I,J,K,1) = Q(J,K,I,4)
            Qout(I,J,K,2) = Q(J,K,I,1)
            Qout(I,J,K,3) = Q(J,K,I,2)
            Qout(I,J,K,4) = Q(J,K,I,3)
            Qout(I,J,K,5) = Q(J,K,I,5)
          End Do
        End Do
      End Do

! Generating file name
      WRITE(index, '(I10)') cue_vis
      filename = TRIM(ADJUSTL(vis_path)) // '/' //
     .           TRIM(ADJUSTL(prefix))//
     .         '.'// TRIM(ADJUSTL(index)) //'.dat'

      Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,
     .                   MPI_MODE_CREATE + MPI_MODE_WRONLY,
     .                   MPI_INFO_NULL, fh, ierr)
      if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file'

      Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
      Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

      disp = 0

! Writing header info
! Number of grids
      nints = 1
      sizeofRecord = nints*intSize
      sizeofBuffer = sizeofRecord 

      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, Ngrid, nints, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset
      disp = disp + sizeofBuffer

! Grid sizes
      nints = 3
      sizeofRecord = nints*intSize 
      sizeofBuffer = sizeofRecord
      
      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, gsizes, nints, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset 
      disp = disp + sizeofBuffer

! Write Ma, alfa, Re & time 
      buf(1) = 0D0         ! Ma
      buf(2) = 0D0         ! alpha
      buf(3) = Re          ! Re
      buf(4) = time        ! time

      nreals = 4
      sizeofRecord = nreals*realSize
      sizeofBuffer = sizeofRecord
      
      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, buf, nreals, MPI_REAL8, 
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset
      disp = disp + sizeofBuffer 

! Writing individual grid sizes

      ndims = 4

      lsizes(1) = Mx
      lsizes(2) = My
      lsizes(3) = Mz
      lsizes(4) = 5

      start_indices(1) = xrank * Mx
      start_indices(2) = yrank * My
      start_indices(3) = zrank * Mz
      start_indices(4) = 0

      local_array_size = Mx * My * Mz * 5

      Call MPI_TYPE_CREATE_SUBARRAY(ndims,gsizes,lsizes,start_indices,
     .                        MPI_ORDER_FORTRAN, MPI_REAL8,
     .                        filetype, ierr)

      Call MPI_TYPE_COMMIT(filetype, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, filetype,
     .                      "native", MPI_INFO_NULL, ierr)
      If(ierr.NE.MPI_SUCCESS)Stop 'Error setting file view'

      Call MPI_FILE_WRITE_ALL(fh, Qout, local_array_size,
     .                        MPI_REAL8, status, ierr)
      if(ierr.NE.MPI_SUCCESS)
     . stop 'Error during collective write of solution file.'

      Call MPI_FILE_CLOSE(fh, ierr)

      Return

      End

************************************************************************
*  Subroutine: Write_Func                                              *
*                                                                      *
*  Function:   Writes the function file in Plot3D format.              *
*              In particular the Grad array [omx, omy, omz, dil]       *
************************************************************************
      Subroutine Write_Func

      Include 'header'
      Include 'mpif.h'

      Integer step
      Real*8  time

      Real*8 Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)
      Real*8 Y(My,Mz,Mx,1)
      Real*8 Grad_out(Mx,My,Mz,5)

      Common  /RKWORK/ Q, Grad
      Common /SPECIES/ Y

*  Locals
      Character*10 index
      Character*25 filename
      Integer      I, J, K, L, M, N
      Real*8       fac1
      Integer      status(MPI_STATUS_SIZE), ierr
      Integer      req(2)
      Integer      gsizes(4), start_indices(4)
      Integer      lsizes(4), ndims
      Integer      count1, count2, count3
      Integer      fh, filetype, local_array_size
      Integer(kind=MPI_OFFSET_KIND)      disp
      Integer      intSize, realSize

      Integer      Ngrid
      Integer      nints, nreals
      Integer      sizeofRecord, sizeofBuffer

      Ngrid = 1            ! Single grid

      gsizes(1) = Nx
      gsizes(2) = Ny
      gsizes(3) = Nz
      gsizes(4) = 5

      count1 = 4
      count2 = 1
      count3 = 1

! Rearranging data to I,J,K-orientation as required by Plot3D format
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Grad_out(I,J,K,1:4) = Grad(J,K,I,1:4)
            Grad_out(I,J,K,5)   = Y(J,K,I,1)
          End Do
        End Do
      End Do

! Generating file name
      WRITE(index, '(I10)') cue_vis
      filename = TRIM(ADJUSTL(vis_path)) // '/' //
     .           TRIM(ADJUSTL(prefix))// 'f' //
     .         '.'// TRIM(ADJUSTL(index)) //'.dat'

      Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,
     .                   MPI_MODE_CREATE + MPI_MODE_WRONLY,
     .                   MPI_INFO_NULL, fh, ierr)
      if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file'

      Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
      Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

      disp = 0

! Writing header info
! Number of grids
      nints = 1
      sizeofRecord = nints*intSize
      sizeofBuffer = sizeofRecord 

      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, Ngrid, nints, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset
      disp = disp + sizeofBuffer

! Grid sizes
      nints = 4
      sizeofRecord = nints*intSize
      sizeofBuffer = sizeofRecord 

      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, gsizes, nints, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
      endif

! Update file offset
      disp = disp + sizeofBuffer

! Writing grids
      ndims = 4

      lsizes(1) = Mx
      lsizes(2) = My
      lsizes(3) = Mz
      lsizes(4) = 5

      start_indices(1) = xrank * Mx
      start_indices(2) = yrank * My
      start_indices(3) = zrank * Mz
      start_indices(4) = 0

      local_array_size = Mx * My * Mz * 5

      Call MPI_TYPE_CREATE_SUBARRAY(ndims,gsizes,lsizes,start_indices,
     .                        MPI_ORDER_FORTRAN, MPI_REAL8,
     .                        filetype, ierr)

      Call MPI_TYPE_COMMIT(filetype, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, filetype,
     .                      "native", MPI_INFO_NULL, ierr)
      If(ierr.NE.MPI_SUCCESS)Stop 'Error setting file view'

      Call MPI_FILE_WRITE_ALL(fh, Grad_out, local_array_size,
     .                        MPI_REAL8, status, ierr)
      if(ierr.NE.MPI_SUCCESS)
     . stop 'Error during collective write of solution file.'

      Call MPI_FILE_CLOSE(fh, ierr)

      End Subroutine Write_Func

************************************************************************
*  Subroutine: Write_Restart                                           *
*                                                                      *
*  Function:   All processors sends data to Rank 0 processor.  Rank 0  *
*              processor receives them and write them to restart file  *
*              at the end of a specific number of time steps or when   *
*              numerical instability may occur.  The restart file is   *
*              used to continue a run at a later time.                 *
*                                                                      *
*              The data are recorded in a total of Px sections in      *
*              descending order of processor rank.                     *
************************************************************************
      Subroutine Write_Restart(time, step)

      Include 'header'
      Include 'mpif.h'

      Integer step
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), Y(My,Mz,Mx,1)
      Real*8  Q_write(My,Mz,Mx,6)

      Common  /RKWORK/ Q
      Common  /SPECIES/ Y

*  Locals
      Character*10 index
      Character*25 filename
      Integer      I, J, K, L, M, N
      Integer      Nx_read, Ny_read, Nz_read, Px_read
      Integer      status(MPI_STATUS_SIZE), ierr
      Integer      gsizes(4), start_indices(4)
      Integer      lsizes(4), ndims
      Integer      count1, count2, count3
      Integer      fh, filetype, local_array_size
      Integer(kind=MPI_OFFSET_KIND) disp
      Integer intSize, realSize

      Q_write(:,:,:,1:5) = Q(:,:,:,1:5)
      Q_write(:,:,:,6)   = Y(:,:,:,1)

      gsizes(1) = Ny
      gsizes(2) = Nz
      gsizes(3) = Nx
      gsizes(4) = 6

      count1 = 4
      count2 = 1
      count3 = 1

! Generating file name
      WRITE(index, '(I10)') cue_restart
      filename = TRIM(ADJUSTL(restart_path)) // '/' // 
     .           TRIM(ADJUSTL(prefix)) //'.restart.'// 
     .           TRIM(ADJUSTL(index))

      Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,
     .                   MPI_MODE_CREATE + MPI_MODE_WRONLY,
     .                   MPI_INFO_NULL, fh, ierr)
      if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file'

      Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
      Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

! Only processor 0 writes the header
      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, gsizes, count1, MPI_INTEGER,
     .                     MPI_STATUS_IGNORE, ierr)
        Call MPI_FILE_WRITE(fh, time, count2, MPI_REAL8,
     .                     MPI_STATUS_IGNORE, ierr)
        Call MPI_FILE_WRITE(fh, step, count3, MPI_INTEGER,
     .                     MPI_STATUS_IGNORE, ierr)
      endif

      ndims = 4

      lsizes(1) = My
      lsizes(2) = Mz
      lsizes(3) = Mx
      lsizes(4) = 6 

      start_indices(1) = yrank * My
      start_indices(2) = zrank * Mz
      start_indices(3) = xrank * Mx
      start_indices(4) = 0

      local_array_size = Mx * My * Mz * 6

      disp = 5*intSize + realSize

      Call MPI_TYPE_CREATE_SUBARRAY(ndims,gsizes,lsizes,start_indices,
     .                        MPI_ORDER_FORTRAN, MPI_REAL8,
     .                        filetype, ierr)

      Call MPI_TYPE_COMMIT(filetype, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, filetype,
     .                      "native", MPI_INFO_NULL, ierr)
      If(ierr.NE.MPI_SUCCESS)Stop 'Error setting file view'

      Call MPI_FILE_WRITE_ALL(fh, Q_write, local_array_size,
     .                        MPI_REAL8, status, ierr)
      if(ierr.NE.MPI_SUCCESS)
     . stop 'Error during collective write of solution file.'

      Call MPI_FILE_CLOSE(fh, ierr)
      if(ierr.NE.MPI_SUCCESS)Stop 'Error closing file'

      Return

      End

************************************************************************
*  Subroutine: Read_Restart                                            *
*                                                                      *
*  Function:   Rank 0 processor reads data from restart file and sends *
*              them to all other processors.                           *
************************************************************************
      Subroutine Read_Restart(time, step)

      Include 'header'
      Include 'mpif.h'

      Integer step
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), Y(My,Mz,Mx,1)
      Real*8  Q_read(My,Mz,Mx,6)

      Common  /RKWORK/ Q
      Common  /SPECIES/ Y

*  Locals
      Character*10 index
      Character*25 filename
      Integer      I, J, K, L, M, N
      Integer      Nx_read, Ny_read, Nz_read, Px_read
      Integer      status(MPI_STATUS_SIZE), ierr
      Integer      gsizes(4), start_indices(4)
      Integer      lsizes(4), ndims
      Integer      count1, count2, count3
      Integer      fh, filetype, local_array_size
      Integer(kind=MPI_OFFSET_KIND) disp
      Integer intSize, realSize

      gsizes(1) = Ny
      gsizes(2) = Nz
      gsizes(3) = Nx
      gsizes(4) = 6 

      count1 = 4
      count2 = 1
      count3 = 1

! Generating file name
      WRITE(index, '(I10)') cue_restart
      filename = TRIM(ADJUSTL(restart_path)) // '/' //
     .           TRIM(ADJUSTL(prefix)) //'.restart.'//
     .           TRIM(ADJUSTL(index))

      Call MPI_FILE_OPEN(MPI_COMM_WORLD, filename,
     .                   MPI_MODE_RDONLY,
     .                   MPI_INFO_NULL, fh, ierr)
     
!      print *, fh, filename, ierr, MPI_SUCCESS  
!      if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file'
 
      Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
      Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

      Call MPI_FILE_READ(fh, gsizes, count1, MPI_INTEGER,
     .                   MPI_STATUS_IGNORE, ierr)
      Call MPI_FILE_READ(fh, time, count2, MPI_REAL8,
     .                   MPI_STATUS_IGNORE, ierr)
      Call MPI_FILE_READ(fh, step, count3, MPI_INTEGER,
     .                   MPI_STATUS_IGNORE, ierr)

      ndims = 4

      lsizes(1) = My
      lsizes(2) = Mz
      lsizes(3) = Mx
      lsizes(4) = 6

      start_indices(1) = yrank * My
      start_indices(2) = zrank * Mz
      start_indices(3) = xrank * Mx
      start_indices(4) = 0

      local_array_size = Mx * My * Mz * 6

      disp = 5*intSize + realSize

      Call MPI_TYPE_CREATE_SUBARRAY(ndims,gsizes,lsizes,start_indices,
     .                        MPI_ORDER_FORTRAN, MPI_REAL8,
     .                        filetype, ierr)

      Call MPI_TYPE_COMMIT(filetype, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL8, filetype,
     .                      "native", MPI_INFO_NULL, ierr)
      If(ierr.NE.MPI_SUCCESS)Stop 'Error setting file view'

      Call MPI_FILE_READ_ALL(fh, Q_read, local_array_size,
     .                        MPI_REAL8, status, ierr)
      if(ierr.NE.MPI_SUCCESS)
     . stop 'Error during collective read of solution file.'

      Call MPI_FILE_CLOSE(fh, ierr)
      if(ierr.NE.MPI_SUCCESS)Stop 'Error closing file'

      Q(:,:,:,1:5) = Q_read(:,:,:,1:5)
      Y(:,:,:,1)   = Q_read(:,:,:,6)

      Return

      End


************************************************************************
*  Subroutine: Output                                                  *
*                                                                      *
*  Function:   It records the instantaneous flow variables to restart  *
*              and/or data files.                                      *
************************************************************************
      Subroutine Output(time, step)

      Include 'header'

      Integer step
      Real*8  time

      Real*8 Q(My,Mz,Mx,5)

      Common /RKWORK/ Q

      Real*8 Y(My,Mz,Mx,1)
 
      Common /SPECIES/ Y

*  Compute Vorticity, Enstrophy and dilatation
      Call Compute_Grad
      Call Compute_KE(Q);   Call Compute_Mom(Q)
      Call Compute_Mass(Q,Y); Call Compute_Energy(Q)
      Call Compute_Entropy(Q,Y)

* Computing hyperviscosity      
      Call Compute_gamma(Y)
      Call Compute_hyper(Q,Y)

*  Compute flow statistics
!       Call Compute_Vorticity_Budget(Q,Y)
       Call Compute_SRM_Stats
!       Call Compute_Favre_Stats
!       Call Generate_Shock_Surface

*  Write conservative variables, vorticity and dilatation
!         Call Write_Visual_Silo(time, step)
!         cue_vis = cue_vis + 1

      Return

      End

