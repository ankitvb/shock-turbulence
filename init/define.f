************************************************************************
*  Subroutine: Define_Param                                            *
*                                                                      *
*  Function:   Rank 0 processor reads in all the input and then        *
*              broadcasts the information to all other processors.     *
************************************************************************
      Subroutine Define_Param

      Include 'header'
      Include 'mpif.h'

*  Locals
      Character*80 comments
      Integer      ierr

*  Rank 0 processor reads input parameters
      If (rank .EQ. 0) Then
         open(unit=in_unit, file='srm.input.1')

*  Input/Output path and filename
        Read(in_unit,1) comments
        Read(in_unit,5) vis_path
        Read(in_unit,5) restart_path
        Read(in_unit,5) stats_path
        Read(in_unit,5) probe_path
        Read(in_unit,5) source_path
        Read(in_unit,5) prefix
        Read(in_unit,*)

*  Grid parameters
        Read(in_unit,1) comments
        Read(in_unit,*) Grid_X
        Read(in_unit,*) Grid_Y
        Read(in_unit,*) x_start
        Read(in_unit,*) y_start
        Read(in_unit,*) z_start
        Read(in_unit,*) period_x
        Read(in_unit,*) period_y
        Read(in_unit,*) period_z
        Read(in_unit,*)

        x_leng = period_x !* (Nx-1)/Nx
        y_leng = period_y !* (Ny-1)/Ny
        z_leng = period_z !* (Nz-1)/Nz

*  Flow parameters
        Read(in_unit,1) comments
        Read(in_unit,*) Re
        Read(in_unit,*) Pr
        Read(in_unit,*) Sc
        Read(in_unit,*) gamma1
        Read(in_unit,*) gamma2
        Read(in_unit,*) epsilon
        Read(in_unit,*) vis_exp
        Read(in_unit,*) bulk_ratio
        Read(in_unit,*) T_ratio
        Read(in_unit,*) hyper_param
        Read(in_unit,*)

*  Time advancement parameters
        Read(in_unit,1) comments
        Read(in_unit,*) march_mode
        Read(in_unit,*) Tstep
        Read(in_unit,*) cfl_cal
        Read(in_unit,*) Nsteps
        Read(in_unit,*) run_time
        Read(in_unit,*) comp_stats
        Read(in_unit,*)

*  Input/Output parameters
        Read(in_unit,1) comments
        Read(in_unit,*) out_grid
        Read(in_unit,*) out_init
        Read(in_unit,*) out_incr_x
        Read(in_unit,*) out_incr_y
        Read(in_unit,*) out_incr_z
        Read(in_unit,*) out_restart
        Read(in_unit,*) out_vis
        Read(in_unit,*) out_probe
        Read(in_unit,*) out_stats
        Read(in_unit,*) out_source

      End If

*  Rank 0 processor broadcasts input parameters to all other processors

*  Input/Output path and filename
      Call MPI_BCAST(vis_path,     20, MPI_CHARACTER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(restart_path, 20, MPI_CHARACTER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(stats_path,   20, MPI_CHARACTER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(probe_path,   20, MPI_CHARACTER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(source_path,  20, MPI_CHARACTER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(prefix,       20, MPI_CHARACTER, 0, 
     .               MPI_COMM_WORLD, ierr)

*  Grid parameters
      Call MPI_BCAST(Grid_X,            1, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Grid_Y,            1, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(x_start,            1, MPI_REAL8,   0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(y_start,            1, MPI_REAL8,   0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(z_start,            1, MPI_REAL8,   0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(x_leng,            1, MPI_REAL8,   0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(y_leng,           1, MPI_REAL8,   0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(z_leng,            1, MPI_REAL8,   0, 
     .               MPI_COMM_WORLD, ierr)

*  Flow parameters
      Call MPI_BCAST(Re,         1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(Pr,         1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(Sc,         1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(gamma1,     1, MPI_REAL8,   0, MPI_COMM_WORLD,
     .               ierr)
      Call MPI_BCAST(gamma2,     1, MPI_REAL8,   0, MPI_COMM_WORLD,
     .               ierr)
      Call MPI_BCAST(epsilon,    1, MPI_REAL8,   0, MPI_COMM_WORLD,
     .               ierr)
      Call MPI_BCAST(vis_exp,    1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(bulk_ratio, 1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(T_ratio,    1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(hyper_param,1, MPI_INTEGER, 0, MPI_COMM_WORLD,
     .               ierr)

*  Time advancement parameters
      Call MPI_BCAST(march_mode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(Tstep,      1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(cfl_cal,    1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(Nsteps,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(run_time,   1, MPI_REAL8,   0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(comp_stats, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)

*  Input/Output parameters
      Call MPI_BCAST(out_grid,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_init,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_incr_x,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_incr_y,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_incr_z,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_restart, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_vis,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_probe,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_stats,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)
      Call MPI_BCAST(out_source,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     .               ierr)

*  Define universal constants
      Pi  = 4.0D0 * atan(1.0D0)
      eye = (0.0D0, 1.0D0)

*  Define inverse of Reynolds number (allow inviscid calculation)
      If (Re .GT. 0.0D0) Then
        Re_inv = 1.0D0 / Re
      Else
        Re_inv = 0.0D0
      End If

      If (Pr .GT. 0.0D0) Then
        Pr_inv = 1.0D0 / Pr
      Else
        Pr_inv = 0.0D0
      End If

*  Define grid size of visualization file
      If (Mod(Nx,out_incr_x) .EQ. 0) Then
        Nx_write = Nx / out_incr_x

      Else
        Nx_write = Nx / out_incr_x + 1

      End If

      If (Mod(Ny,out_incr_y) .EQ. 0) Then
        Ny_write = Ny / out_incr_y

      Else
        Ny_write = Ny / out_incr_y + 1

      End If

      If (Mod(Nz,out_incr_z) .EQ. 0) Then
        Nz_write = Nz / out_incr_z

      Else
        Nz_write = Nz / out_incr_z + 1

      End If

*  Set output cue vector
      cue_vis = 1
      cue_restart = 1

 1    Format(A80)
 5    Format(A20)

      Return

      End


************************************************************************
*  Subroutine: Check_Param                                             *
*                                                                      *
*  Function:   It performs some basic checking on input parameters.    *
************************************************************************
      Subroutine Check_Param

      Include 'header'

*  Number of grid points in streamwise and transverse directions must 
*  be integral multiples of total number of processors
      If (Mod(Nx,Px) .NE. 0) Stop 'Check_Param: Problem with Nx!'
      If (Mod(Ny,Py) .NE. 0) Stop 'Check_Param: Problem with Ny!'
      If (Mod(Nz,Pz) .NE. 0) Stop 'Check_Param: Problem with Nz!'

*  Check CFL for accuracy
c      If ( (march_mode .EQ. 1  .OR.  march_mode .EQ. 3) .AND. 
c     .         cfl_cal .GT. resv_limit )
c     .  Stop 'Check_Param: cfl_cal too large!'

      Return

      End


************************************************************************
*  Subroutine: Write_Info                                              *
*                                                                      *
*  Function:   It writes out the input parameters to an output file.   *
************************************************************************
      Subroutine Write_Info

      Include 'header'

*  Locals
      Integer M, N

      Write(out_unit,5)  'Visualization through ', vis_path,     prefix
      Write(out_unit,5)  'Restart Input through ', restart_path, prefix
      Write(out_unit,5)  'Stats Output  through ', stats_path,   prefix
      Write(out_unit,5)  'Probe Output  through ', probe_path,   prefix
      Write(out_unit,5)  'Source Output through ', source_path,  prefix
      Write(out_unit,*)
      Write(out_unit,10) 'Grid Size:'
      Write(out_unit,15) 'Px  = ', Px
      Write(out_unit,15) 'Py  = ', Py
      Write(out_unit,15) 'Pz  = ', Pz
      Write(out_unit,15) 'Nx  = ', Nx
      Write(out_unit,15) 'Ny  = ', Ny
      Write(out_unit,15) 'Nz  = ', Nz
      Write(out_unit,*)
      Write(out_unit,10) 'Grid Parameters:'
      Write(out_unit,15) 'Grid_X            = ', Grid_X
      Write(out_unit,15) 'Grid_Y            = ', Grid_Y
      Write(out_unit,20) 'x_start           = ', x_start
      Write(out_unit,20) 'y_start           = ', y_start
      Write(out_unit,20) 'z_start           = ', z_start
      Write(out_unit,20) 'x_leng            = ', x_leng
      Write(out_unit,20) 'y_leng            = ', y_leng
      Write(out_unit,20) 'z_leng            = ', z_leng
      Write(out_unit,*)
      Write(out_unit,10) 'Flow Parameters:'
      Write(out_unit,20) 'Re         = ', Re
      Write(out_unit,20) 'Pr         = ', Pr
      Write(out_unit,20) 'Sc         = ', Sc
      Write(out_unit,20) 'gamma1     = ', gamma1
      Write(out_unit,20) 'gamma2     = ', gamma2
      Write(out_unit,20) 'epsilon    = ', epsilon
      Write(out_unit,20) 'vis_exp    = ', vis_exp
      Write(out_unit,20) 'bulk_ratio = ', bulk_ratio
      Write(out_unit,20) 'T_ratio    = ', T_ratio
      Write(out_unit,15) 'hyper_param= ', hyper_param
      Write(out_unit,*)
      Write(out_unit,10) 'Time Advancement Parameters:'
      Write(out_unit,15) 'march_mode = ', march_mode
      Write(out_unit,20) 'Tstep      = ', Tstep
      Write(out_unit,20) 'cfl_cal    = ', cfl_cal
      Write(out_unit,15) 'Nsteps     = ', Nsteps
      Write(out_unit,20) 'run_time   = ', run_time
      Write(out_unit,15) 'comp_stats = ', comp_stats
      Write(out_unit,*)
      Write(out_unit,10) 'Input/Output Parameters:'
      Write(out_unit,15) 'out_grid    = ', out_grid
      Write(out_unit,15) 'out_init    = ', out_init
      Write(out_unit,15) 'out_incr_x  = ', out_incr_x
      Write(out_unit,15) 'out_incr_y  = ', out_incr_y
      Write(out_unit,15) 'out_incr_z  = ', out_incr_z
      Write(out_unit,15) 'out_restart = ', out_restart
      Write(out_unit,15) 'out_vis     = ', out_vis
      Write(out_unit,15) 'out_probe   = ', out_probe
      Write(out_unit,15) 'out_stats   = ', out_stats
      Write(out_unit,15) 'out_source  = ', out_source
      Write(out_unit,*)
      Write(out_unit,*)

 5    Format(1X, 3A)
 10   Format(1X, A)
 15   Format(1X, A, I7)
 20   Format(1X, A, E15.8)
 25   Format(1X, A, I2, A, I2, A, 2E15.8)
 30   Format(1X, A, I2, A, I2, A, E15.8)

      Return

      End


************************************************************************
*  Function: Leng_Char                                                 *
*                                                                      *
*  It returns the real length of a character string which is defined   *
*  as the number of characters before the first space.                 *
************************************************************************
      Function Leng_Char(string)

      Character*(*) string
      Integer Leng_Char, l

      l = 1

 1    If (string(l:l) .EQ. ' ')  Go To 5

      l = l + 1
      Go To 1

 5    Continue

      Leng_Char = l - 1

      Return

      End
