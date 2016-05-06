**********************************************************************
* Code to generate spherical shock with density interface            *
* Output restart file                                                *
**********************************************************************
      Program genfield      

      include 'header'
      Include 'mpif.h'

      Integer step
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), P(My,Mz,Mx)
      Real*8  Q_write(My,Mz,Mx,6)
      Real*8  Y(My,Mz,Mx,1)
      Real*8  gamma(My,Mz,Mx)
      Real*8  rad(My,Mz,Mx), r_interface(My,Mz,Mx)
      Real*8  S1(40), S2(40), S3(40), S4(40)

      Common  Q_aux, Y_aux, S1, S2, S3, S4

* MPI vars
      Integer      status(MPI_STATUS_SIZE), ierr
      Integer      req(2)
      Integer      gsizes(4), start_indices(4)
      Integer      lsizes(4), ndims
      Integer      count1, count2, count3
      Integer      fh, filetype, local_array_size
      Integer(kind=MPI_OFFSET_KIND) disp
      Integer intSize, realSize

      Real*8 start_time, end_time, Tic_time

*  Locals
      Character*15 concat
      Integer      I, J, K, L, M, N
      Integer      II, JJ, KK
      Real*8       IX, IY, IZ
      
      Real*8 fac, fac1, fac2, base_pressure, Ms, num, den
      Real*8 r_shock, eps

* Perturbed interface     
      Integer LMAX
  
      Parameter(LMAX=16)
      Real*8 Y_lm(LMAX,LMAX), Y_l0(LMAX)
      Real*8 A_lm(LMAX,LMAX), phase_lm(LMAX,LMAX), phase_l0(LMAX)
      Real*8 A_l(LMAX)
      Real*8 cost, theta, phi, R0, A, sum

*  Initialize Message Passing Interface
      Call Initialize_MPI(Px,Py,Pz)

*  Define problem parameters
      Call Define_Param
      If (rank .EQ. Px-1) Call Check_Param

* Setting up grid
      Call Setup_Grid
 
      step = 0
      time = 0d0

********************************************************************************************
* Setting up random phase shift and amplitude for each spherical harmonic mode
      phase_lm(:,:) = 0D0
      phase_l0(:)   = 0D0
 
      if(rank .EQ. 0)then
        Do L = 1,LMAX
          phase_l0(L) = 2 * Pi * RAND()
          A_l(L) = RAND()
          Do M = 1,L
            phase_lm(L,M) = 2 * Pi * RAND()
            A_lm(L,M) = RAND()
          End Do
        End Do
      endif

c      if(rank .EQ. 0)then
c        Do l = 1,LMAX
c          print *, phase_l0(l),l,'0'
c          Do m = 1,l
c            print *, phase_lm(l,m),l,m
c          End Do
c        End Do
c      endif

* Sending it to all processors
      Call MPI_BCAST(phase_lm, LMAX*LMAX, MPI_REAL8, 0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(phase_l0, LMAX, MPI_REAL8, 0,
     .               MPI_COMM_WORLD, ierr)

********************************************************************************************
* Setting up perturbed material interface
      R0 = 1.0D0                          ! Initial radius of material interface (from 0.25)
      A  = 2.0D-2                        ! Initial amplitude of perturbation
      r_interface(:,:,:) = R0             ! Default interface location

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            II = Mx * xrank + I
            JJ = My * yrank + J
            KK = Mz * zrank + K
            IX = REAL(II) - 1.0D0
            IY = REAL(JJ) - 1.0D0
            IZ = REAL(KK) - 1.0D0

            rad(J,K,I) = sqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

            if((II.EQ.1).AND.(JJ.EQ.1).AND.(KK.EQ.1))then
              cost = 1D0
            else
              cost = IZ*DZ/rad(J,K,I)
            endif

            if(cost .GT. 1d0)then                   ! Making sure argument is within bounds
              cost = 1d0
            elseif(cost .LT. -1d0)then
              cost = -1d0
            endif

            theta = dacos(cost) 
            phi   = datan(IY*DY/(IX*DX+1D-8))

            Call Spherical_Harmonics(Y_lm, Y_l0, phase_lm, phase_l0,
     .                   theta, phi, LMAX)

            sum = 0D0
            Do L = 64,64
!              sum = sum + A * Y_l0(L)
!              Do M = 1,L
                sum = sum +  A *
     .                dcos(DBLE(L)*phi + phase_l0(L)) *
     .                dcos(DBLE(L)*theta + phase_l0(L)) !* Y_lm(L,M)
!              End Do
            End Do

            r_interface(J,K,I) = r_interface(J,K,I) + sum
          End Do
        End Do
      End Do

********************************************************************************

      Ms = 2.0D0 ! Shock mach number (1.5 previously)
      base_pressure = 0.71429D0
      eps = 0.1D0

*  Generating Initial field (converging shock, smooth interface)
      r_shock            = 1.1D0             ! Initial radius of shock (from 0.5)

      Do I = 1,Mx
       Do J = 1,My
         Do K = 1,Mz

          II = Mx * xrank + I
          JJ = My * yrank + J
          KK = Mz * zrank + K
          IX = REAL(II) - 1.0D0
          IY = REAL(JJ) - 1.0D0
          IZ = REAL(KK) - 1.0D0

          rad(J,K,I) = sqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)
          if((II.EQ.1).AND.(JJ.EQ.1).AND.(KK.EQ.1))then
            fac2 = 0D0
          else
            fac2 = 1D0 / rad(J,K,I) 
          endif

          Y(J,K,I,1) = 5.0D-1 *
     .       (1D0 + tanh(200D0*(rad(J,K,I)-r_interface(J,K,I))/2D0) )
c           Y(J,K,I,1) = 1D0

          fac  = 1D0/(Y(J,K,I,1) + epsilon*(1D0-Y(J,K,I,1)))
          fac1 = (2D0/(gamma1+1D0))*(Ms*Ms-1D0)/Ms

          Q(J,K,I,4) = fac + 5.0D-1*
     .       (1D0 + tanh(200D0*(rad(J,K,I)-r_shock)/2D0) ) *
     .       (2.44D0 - fac)                                   !( (gamma1+1D0)*Ms*Ms/((gamma1-1D0)*Ms*Ms+2D0) - fac)
  
          Q(J,K,I,1) = 1D-50 !Q(J,K,I,4) * 5.0D-1*
c     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock)/1D0)) *
c     .       (2D0/(gamma1+1D0)) * ((Ms*Ms - 1D0)/Ms )
c     .      * (-IX*DX * fac2) 

          Q(J,K,I,2) = 1D-50 !Q(J,K,I,4) * 5.0D-1*
c     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock)/1D0)) *
c     .       (2D0/(gamma1+1D0)) * ((Ms*Ms - 1D0)/Ms)
c     .      * (-IY*DY * fac2)

          Q(J,K,I,3) = 1D-50 !Q(J,K,I,4) * 5.0D-1*
c     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock)/1D0)) *
c     .       (2D0/(gamma1+1D0)) * ((Ms*Ms - 1D0)/Ms)
c     .      * (-IZ*DZ * fac2)

          P(J,K,I) = base_pressure + base_pressure * 5.0D-1 *
     .       (1D0 + tanh(200D0*(rad(J,K,I)-r_shock)/2D0)) *
     .       (3.52D0  - 1)                                     !(2D0*gamma1/(gamma1+1D0)) * (Ms*Ms - 1D0)

         End Do
       End Do
      End Do

******************************************************************************
*  Compute effective gamma from mass fractions
*
*                Y1 gm1          Y2 gm2
*               --------- + eps ---------
*                gm1 - 1         gm2 - 1
*  gamma_eff = --------------------------
*                  Y1              Y2
*               --------- + eps ---------
*                gm1 - 1         gm2 - 1

      Do I = 1,Mx
       Do J = 1,My
        Do K = 1,Mz
          fac1 = gamma1 - 1D0
          fac2 = gamma2 - 1D0

          num = ( Y(J,K,I,1) * gamma1 / fac1 )
     .        + ( epsilon * (1D0 - Y(J,K,I,1)) * gamma2 / fac2 )
          den = ( Y(J,K,I,1) / fac1 )
     .        + ( epsilon * (1D0 - Y(J,K,I,1)) / fac2 )

          gamma(J,K,I) = num / den
        End Do
       End Do
      End Do

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             fac1 = 1.0D0/(gamma(J,K,I) - 1.0D0)
             Q(J,K,I,5) = fac1 * P(J,K,I)
             Q(J,K,I,5) = Q(J,K,I,5) + (5.0D-1/Q(J,K,I,4))
     .                                       *(Q(J,K,I,1)*Q(J,K,I,1)
     .                               +         Q(J,K,I,2)*Q(J,K,I,2)
     .                               +         Q(J,K,I,3)*Q(J,K,I,3))
          End Do
        End Do
      End Do

      Q_write(:,:,:,1:5) = Q(:,:,:,1:5)
      Q_write(:,:,:,6)   = Y(:,:,:,1)

      gsizes(1) = Ny
      gsizes(2) = Nz
      gsizes(3) = Nx
      gsizes(4) = 6

      count1 = 4
      count2 = 1
      count3 = 1

      Call MPI_FILE_OPEN(MPI_COMM_WORLD, 'srm.restart.0',
     .                   MPI_MODE_CREATE + MPI_MODE_WRONLY,
     .                   MPI_INFO_NULL, fh, ierr)
      if(ierr.NE.MPI_SUCCESS) Stop 'Error opening file'

      Call MPI_TYPE_SIZE(MPI_INTEGER, intSize, ierr)
      Call MPI_TYPE_SIZE(MPI_REAL8, realSize, ierr)

      if(rank .EQ. 0)then
        Call MPI_FILE_WRITE(fh, gsizes, count1, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
        Call MPI_FILE_WRITE(fh, time, count2, MPI_REAL8,
     .                      MPI_STATUS_IGNORE, ierr)
        Call MPI_FILE_WRITE(fh, step, count3, MPI_INTEGER,
     .                      MPI_STATUS_IGNORE, ierr)
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

* For scaling
      start_time = MPI_WTIME()

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

      Call MPI_TYPE_FREE(filetype, ierr)

      Call MPI_FILE_CLOSE(fh, ierr)
      if(ierr.NE.MPI_SUCCESS)Stop 'Error closing file'

      Call MPI_FINALIZE(ierr)

* For scaling
      end_time = MPI_WTIME()
      Tic_time = end_time - start_time

      if(rank .EQ. 0)then
       print *, Tic_time
      endif

      Stop

      End

