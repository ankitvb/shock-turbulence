************************************************************************
*  Program: Shock_Cell                                                 *
*  Author:  Calvin Lui  (based on Jonathan Freund's code)              *
*  Date:    August 14, 2001                                            *
*                                                                      *
*  This code solves the three-dimensional, compressible Navier-Stokes  *
*  equations in its conservation form in Cartesian coordinate system   *
*  on parallel machines using Message Passing Interface (MPI).         *
*                                                                      *
*  The spatial discretization scheme is as follows:                    *
*  streamwise direction - optimized sixth order "spectral-like" scheme *
*  transverse direction - optimized sixth order "spectral-like" scheme *
*  spanwise direction   - Fourier spectral scheme                      *
*                                                                      *
*  The time advancement method is two-step (5-stage and 6-stage)       *
*  fourth order Low-dissipation-and-dispersion Runge-Kutta scheme.     *
*                                                                      *
*  The components of the conservative variable vector are arranged in  *
*  the following sequence:                                             *
*    [ rho u_x,  rho u_y,  rho u_z,  rho,  rho e_t ]                   *
************************************************************************
      Program Main

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), dQ(My,Mz,Mx,5), QW(My,Mz,Mx,5)
      Common /RKWORK/  Q, dQ, QW

      Real*8  Y(My,Mz,Mx,1), drYi(My,Mz,Mx,1), rYW(My,Mz,Mx,1)
      Common /SPECIES/ Y, drYi, rYW

      Real*8  gamma_eff(My,Mz,Mx)
      Common /THERMO/  gamma_eff

      Real*8  Qf(My,Mz,Mx,5), Yf(My,Mz,Mx,1), Q_aux(My,Mz,Mx,5)

*  Locals
      Integer ierr, req(3)
      Integer status(MPI_STATUS_SIZE)
      Integer step, step_count
      Real*8  cfl_stab, cfl_resv, time

*  Locals for timing
      Real*8  start_time, end_time, mpi_node_time
      Real*8  Tic_time, Ttotal_time, Trk_time, Tavg_time
 
      Integer L, I, J, K, M, II

       Call Initialize_MPI(Px,Py,Pz)

*  Start timing
      mpi_start_time = MPI_WTIME()
      Call MPI_BCAST(start_time, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

      mpi_node_time = start_time

*  Define problem parameters
      Call Define_Param
      If (rank .EQ. Px-1) Call Check_Param

      If (rank .EQ. 0) Call Write_Info

*  Generate grid
      Call Setup_Grid

*  Initialize derivatives
      Call Setup_Derivatives_Sym
      Call Setup_Derivatives_Asym

*  Initialize filters
      Call Init_Filters_Sym
      Call Init_Filters_Asym

*  Read in restart file
      cue_restart = 82 
      Call Read_Restart(time, step)
      cue_restart = cue_restart + 1

!      Call Gaussian_Filter(Q, Qf, 5)
!      Call Gaussian_Filter(Y, Yf, 1)
!      Q = Qf
!      Y = Yf

*  Initialize starting time of simulation for statistics purpose
      sim_start_time = time

*  Open output file and write out problem parameters
      OPEN(UNIT=out_unit,FILE="output",FORM="FORMATTED",
     .      STATUS="UNKNOWN", POSITION="APPEND")
      If (rank .EQ. 0) Write(out_unit,5) 'Restart file: time = ', time, 
     .                                                 'step = ', step
      CLOSE(out_unit)

*  Write out initial fields
      cue_vis = 1 
      If (out_init .EQ. 1) Call Output(time, 0)
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

* Computing hyperviscosity terms
      Call Compute_gamma(Y)
      Call Compute_hyper(Q, Y)

*  Compute time step subject to a constant CFL criterion
      If (march_mode .EQ. 1)then
         Call Compute_CFL(Q, cfl_stab, cfl_resv, step)
         Call Compute_Tstep(Q, Y, cfl_resv, step) 
      End If

      Call Compute_KE(Q);   Call Compute_Mom(Q)
      Call Compute_Mass(Q,Y); Call Compute_Energy(Q)
      Call Compute_Entropy(Q,Y)

*  Write column headers
      If (rank .EQ. 0) Then
       OPEN(UNIT=out_unit,FILE="output",FORM="FORMATTED",
     .      STATUS="UNKNOWN", POSITION="APPEND")
       WRITE(out_unit,10) 'n ', 't   ', 'dt    ', 'CFL_resv  ',
     .    'KE     ', 'MaxDil  ','Enstrophy ',
     .       'MomX  ', 'MomY  ','MomZ  ','Mass  ','Energy ','Entropy'
       WRITE(out_unit,11) step, time, Tstep, cfl_resv*Tstep,KE, MaxDil,
     .  Enstrophy, MomX, MomY, MomZ, Mass_t, e_t, Entropy
       CLOSE(out_unit)
      End If

*  Determine initial startup time
      end_time = MPI_WTIME()
      Tic_time = end_time - start_time

*  Time advancement loop

************************************************************************
*  Constant CFL mode: computes with a specified number of time steps   *
*                     (Nsteps)                                         *
************************************************************************
      If (march_mode .EQ. 1) Then
        Do step_count = 1, Nsteps

*  Compute time step subject to a constant CFL criterion
          Call Compute_CFL(Q, cfl_stab, cfl_resv, step)
          Call Compute_Tstep(Q, Y, cfl_resv, step)

*  Advance the flow
          Call LDDRK4(time, step)

*  Write out flow field data if necessary
          Call Output(time, step)

*  Write out time stepping information
          If (rank .EQ. 0) Then
            OPEN(UNIT=out_unit,FILE="output",FORM="FORMATTED",
     .      STATUS="UNKNOWN", POSITION="APPEND")
            WRITE(out_unit,11) step, time, Tstep, cfl_resv*Tstep, KE,
     .       MaxDil,Enstrophy, MomX, MomY, MomZ, Mass_t, e_t, Entropy
            CLOSE(out_unit)
          End If
 
        End Do

************************************************************************
*  Constant time step mode: computes with a specified number of time   *
*                           steps (Nsteps)                             *
************************************************************************
      Else If (march_mode .EQ. 2) Then
        Do step_count = 1, Nsteps

*  Compute CFL
          Call Compute_CFL(Q, cfl_stab, cfl_resv, step)
          cfl_stab = cfl_stab * Tstep
          cfl_resv = cfl_resv * Tstep

*  Check for numerical stability
c          If (cfl_stab .GT. stab_limit) Then
c            Call Write_Restart(time, step, rank)
c            Stop 'Main: High CFL, numerical instability is probable!'
c          End If

*  Advance the flow
          Call LDDRK4(time, step)

*  Write out flow field data if necessary
          Call Output(time, step)

*  Write out time stepping information
          If (rank .EQ. 0) Then
            Write(out_unit,11) step, time, Tstep, cfl_resv*Tstep, KE,
     .       MaxDil,Enstrophy, MomX, MomY, MomZ, Mass_t, e_t, Entropy
            Write(stats2_unit,7) time, KE*(2d0/3d0), rho_rms, P_rms,
     .                  T_rms
            Write(stats1_unit,7) taylor_microscale, Re_lambda,
     .                          velocity_derivative_skewness,
     .                         dilatation_skewness, dilatation_flatness,
     .                         dilatation_square

          End If

        End Do

************************************************************************
*  Constant CFL mode: computes with a specified duration of real time  *
*                     (run_time)                                       *
************************************************************************
      Else If (march_mode .EQ. 3) Then
        step_count = 0

*  Check if run time is up
        Do while (mpi_node_time - start_time .LT. run_time)
          step_count = step_count + 1

*  Compute time step subject to a constant CFL criterion
          Call Compute_CFL(Q, cfl_stab, cfl_resv, step)
          Tstep = cfl_cal / cfl_resv

*  Write out time stepping information
          If (rank .EQ. 0) 
     .      Write(out_unit,10)   'n = ', step,
     .                           't = ', time, 
     .                         'd t = ', Tstep,
     .                    'CFL_stab = ', cfl_stab * Tstep, 
     .                    'CFL_resv = ', cfl_cal

*  Advance the flow
          Call LDDRK4(time, step)

*  Write out flow field data if necessary
          Call Output(time, step)

*  Check the clock
          mpi_node_time = MPI_WTIME()
          Call MPI_BCAST(mpi_node_time, 1, MPI_REAL8, 0, 
     .                   MPI_COMM_WORLD, ierr)

        End Do

************************************************************************
*  Constant time step mode: computes with a specified duration of real *
*                           time (run_time)                            *
************************************************************************
      Else If (march_mode .EQ. 4) Then
        step_count = 0

*  Check if run time is up
        Do while (mpi_node_time - start_time .LT. run_time)
          step_count = step_count + 1

*  Compute CFL
          Call Compute_CFL(Q, cfl_stab, cfl_resv, step)
          cfl_stab = cfl_stab * Tstep
          cfl_resv = cfl_resv * Tstep

*  Check for numerical stability
          If (cfl_stab .GT. stab_limit) Then
            Call Write_Restart(time, step)
            Stop 'Main: High CFL, numerical instability is probable!'
          End If

*  Check for numerical accuracy
          If (cfl_resv .GT. resv_limit) Then
            Call Write_Restart(time, step)
            Stop 'Main: High CFL, numerical inaccuracy is probable!'
          End If

*  Write out time stepping information
          If (rank .EQ. 0) 
     .      Write(out_unit,10)   'n = ', step,
     .                           't = ', time, 
     .                         'd t = ', Tstep,
     .                    'CFL_stab = ', cfl_stab, 
     .                    'CFL_resv = ', cfl_resv

*  Advance the flow
          Call LDDRK4(time, step)

*  Write out flow field data if necessary
          Call Output(time, step)

*  Check the clock
          mpi_node_time = MPI_WTIME()
          Call MPI_BCAST(mpi_node_time, 1, MPI_REAL8, 0, 
     .                   MPI_COMM_WORLD, ierr)

        End Do

      Else
        Stop 'Main: march_mode is not well-defined!'

      End If

* Write Restart file at exit
!      Call Write_Restart(time, step)
!      Call Write_Visual_Silo(time, step)

*  Determine final timing
      end_time = MPI_WTIME()

      Ttotal_time = end_time - start_time
      Trk_time    = Ttotal_time - Tic_time

      If (march_mode .EQ. 1  .OR.  march_mode .EQ. 2) Then
        Tavg_time = Trk_time / Dble(Nsteps)

      Else
        Tavg_time = Trk_time / Dble(step_count)

      End If

*  Write out timing information
      If (rank .EQ. 0) Then
        OPEN(UNIT=out_unit,FILE="output",FORM="FORMATTED",
     .      STATUS="UNKNOWN", POSITION="APPEND")

        Write(out_unit,15) 'Main: IC Time (total) = ', Tic_time, 
     .                     ' wall sec'
        Write(out_unit,15) 'Main: RK Time (total) = ', Trk_time, 
     .                     ' wall sec'
        Write(out_unit,15) 'Main: Avg time / step = ', Tavg_time, 
     .                     ' wall sec'
        Write(out_unit,15) 'Main: Total time      = ', Ttotal_time, 
     .                     ' wall sec'

        CLOSE(out_unit)
      End If

*  Close Message Passage Interface
      Call MPI_FINALIZE(ierr)

 5    Format(1X, A, E15.8, 2X, A, I8)
 6    Format(1x, a7,4(a10))
 7    Format(1x, 6(ES13.5))
 10   Format(1x,a5,12(a13))
 11   Format(1x,i5,12(ES13.5))
 15   Format(1X, A, ES15.8, A)

      Stop

      End
