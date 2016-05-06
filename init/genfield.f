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
      Real*8  Q(My,Mz,Mx,5), P(My,Mz,Mx), S(My,Mz,Mx)
      Real*8  Q_write(My,Mz,Mx,6)
      Real*8  Y(My,Mz,Mx,1)
      Real*8  gamma(My,Mz,Mx), buf(My,Mz,Mx,5)
      Real*8  rad(My,Mz,Mx), r_interface(My,Mz,Mx)

      Common  Q_aux, Y_aux

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
      Real*8 r_shock, eps, base_entropy

* Chisnell shock parameters
      Real*8 alpha, qq, B, D, F, Vs, eta, q1, alpha1
      Real*8 tol, xi, v, ff, dfdv, uus, rsr, PPs, cs, c0

* Perturbed interface     
      Integer LMAX
  
      Parameter(LMAX=32)
      Real*8 Y_lm(LMAX,LMAX), Y_l0(LMAX)
      Real*8 A_lm(LMAX,LMAX), phase_lm(LMAX,LMAX), phase_l0(LMAX)
      Real*8 A_l(LMAX)
      Real*8 cost, sint, theta, phi, psi, R0, A, sum
      Real*8 myrand

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
 
c      if(rank .EQ. 0)then
c        Do L = 1,LMAX
c          phase_l0(L) = 2D0 * Pi * RAND()
c          A_l(L) = RAND()
c          Do M = 1,L
c            phase_lm(L,M) = 2D0 * Pi * RAND()
c            A_lm(L,M) = RAND()
c          End Do
c        End Do
c      endif

c      if(rank .EQ. 0)then
c        Do l = 1,LMAX
c          print *, phase_l0(l),l,'0'
c          Do m = 1,l
c            print *, phase_lm(l,m),l,m
c          End Do
c        End Do
c
c      endif

* Sending it to all processors
c      Call MPI_BCAST(phase_lm, LMAX*LMAX, MPI_REAL8, 0,
c     .               MPI_COMM_WORLD, ierr)
c      Call MPI_BCAST(phase_l0, LMAX, MPI_REAL8, 0,
c     .               MPI_COMM_WORLD, ierr)

********************************************************************************************
* Setting up perturbed material interface
      R0 = 0.9D0                          ! Initial radius of material interface (from 0.9)
      A  = 2.75D-2                        ! Initial amplitude of perturbation (from 2.75D-2)
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
              sint = 0D0
            else
              cost = IZ*DZ/rad(J,K,I)
              sint = dsqrt(1D0 - cost*cost)
            endif

            if(cost .GT. 1d0)then                   ! Making sure argument is within bounds
              cost = 1D0
              sint = 0D0
            elseif(cost .LT. -1d0)then
              cost = -1D0
              sint = 0D0
            endif

            theta = dacos(cost) 
c            theta = datan(IX*DX/(IZ*DZ+1D-8))
            phi   = datan(IY*DY/(IX*DX+1D-8))
            psi   = datan(IX*DX/(IY*DY+1D-8))

            myrand = (DBLE(MOD(FLOOR(theta*1D2),10))/10D0)*2D0*Pi

!            Call Spherical_Harmonics(Y_lm, Y_l0, phase_lm, phase_l0,
!     .                   theta, phi, LMAX)

            sum = 0D0
            Do L = 40,40!20,20 !32,32
!              sum = sum + A * Y_l0(L)
!              Do M = 16,16
                M = IDNINT(sint*DBLE(L))
                if(MOD(M,8).NE.0)then
                  if(MOD(M+1,8).EQ.0)then
                    M = M+1
                  elseif(MOD(M+2,8).EQ.0)then
                    M = M+2
                  elseif(MOD(M+3,8).EQ.0)then
                   M = M+3
                  elseif(MOD(M+4,8).EQ.0)then
                    M = M+4
                  elseif(MOD(M+5,8).EQ.0)then
                    M = M+5
                  elseif(MOD(M+6,8).EQ.0)then
                    M = M+6
                  elseif(MOD(M+7,8).EQ.0)then
                    M = M+7
                  endif
                endif
                
                sum = sum +  A *
     .               dcos(DBLE(L)*theta)*
     .               dcos(DBLE(M)*phi)
!     .              +        A * 0.1 *
!     .                dcos(4D0*DBLE(L)*theta) *
!     .                dcos(4D0*DBLE(M)*phi)     
!              End Do
            End Do

            r_interface(J,K,I) = r_interface(J,K,I) + sum
          End Do
        End Do
      End Do

********************************************************************************

      Ms = 1.8D0 ! Shock mach number 
      base_pressure = 0.71429D0
      base_entropy  = base_pressure
      eps = 0.1D0

*  Generating Initial field (converging shock, smooth interface)
      r_shock            = 0.93D0             ! Initial radius of shock (from 0.93)

      Do I = 1,Mx
       Do J = 1,My
         Do K = 1,Mz

          II = Mx * xrank + I
          JJ = My * yrank + J
          KK = Mz * zrank + K
          IX = DBLE(II) - 1.0D0
          IY = DBLE(JJ) - 1.0D0
          IZ = DBLE(KK) - 1.0D0

          rad(J,K,I) = sqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)
          if((II.EQ.1).AND.(JJ.EQ.1).AND.(KK.EQ.1))then
            fac2 = 0D0
          else
            fac2 = 1D0 / rad(J,K,I) 
          endif

          Y(J,K,I,1) = 5.0D-1 *
     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_interface(J,K,I))) )
c           Y(J,K,I,1) = 1D0

          fac  = 1D0/(Y(J,K,I,1) + epsilon*(1D0-Y(J,K,I,1)))
          fac1 = (2D0/(gamma1+1D0))*(Ms*Ms-1D0)/Ms
          c0 = dsqrt(gamma1*base_pressure/1D0)

          Q(J,K,I,4) = fac + 5.0D-1*
     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock)) ) *
     .       ( (gamma1+1D0)*Ms*Ms/((gamma1-1D0)*Ms*Ms+2D0) - fac) !( (gamma1+1D0)/(gamma1-1D0) - fac)

          P(J,K,I) = base_pressure + base_pressure * 5.0D-1 *
     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock))) *
     .       ((2D0*gamma1/(gamma1+1D0))*(Ms*Ms-1D0) + 1D0 - 1D0) !(c0*c0*(2D0/(gamma1+1D0)*Ms*Ms)  - 1D0) 

c          Q(J,K,I,4) = (P(J,K,I)/S(J,K,I))**(1D0/gamma1)

          cs = dsqrt(gamma1*P(J,K,I)/Q(J,K,I,4))

          Q(J,K,I,1) = -c0 * 5.0D-1 *
     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock))) *
     .       (2D0/(gamma1+1D0))*(Ms*Ms-1)/Ms * IX*DX/rad(J,K,I) ! (2D0/(gamma1+1D0))*Ms
          Q(J,K,I,2) = -c0 * 5.0D-1 *
     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock))) *
     .       (2D0/(gamma1+1D0))*(Ms*Ms-1)/Ms * IY*DY/rad(J,K,I) ! (2D0/(gamma1+1D0))*Ms
          Q(J,K,I,3) = -c0 * 5.0D-1 *
     .       (1D0 + tanh(100D0*(rad(J,K,I)-r_shock))) *
     .       (2D0/(gamma1+1D0))*(Ms*Ms-1)/Ms * IZ*DZ/rad(J,K,I) ! (2D0/(gamma1+1D0))*Ms

          if((xrank.EQ.0).AND.(yrank.EQ.0).AND.(zrank.EQ.0))then
             Q(1,1,1,1:3) = 0D0
          endif
         End Do
       End Do
      End Do

******************************************************************************
* Computing Chisnell's post-shock solution
      alpha  = 0.71716D0
      qq     = -9.86D0
      B      = 12.643D0
      D      = 12.65D0
      F      = -3.91D0
      Vs     = 0.5976D0
      eta    = 0.2312D0
      q1     = qq/Vs
      alpha1 = alpha/Vs
      tol    = 1D-3

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

          if((rad(J,K,I).GT.(1.01D0*r_shock)))then   ! 1.04*r_shock for tanh profile
            xi = rad(J,K,I)/r_shock       ! Begin iteration
            v  = 0.01D0
            ff = 10D0
            Do while(dabs(ff) > tol)
              ff   = (v**alpha)*((v+q1)/(1D0+q1))**(-F) - 1D0/xi
              dfdv = alpha*(v**(alpha-1D0))*((v+q1)/(1D0+q1))**(-F) 
     .             - F*(v**alpha)/(1D0+q1)*((v+q1)/(1D0+q1))**(-F-1D0)
             v    = v - ff/dfdv
            End Do
            uus = xi * v
            rsr = ((alpha1-v)/(alpha1-1D0))**(-eta)*((v+q1)/(1D0+q1))
     .            **(-D)
            PPs = v**(2D0*(1D0-alpha))*((v+q1)/(1D0+q1))**(B+D+2D0*F) ! End iteration

         ! Setting post-shock values
c            Q(J,K,I,4)   = Q(J,K,I,4)/rsr
c            P(J,K,I)     = P(J,K,I)*PPs
c            Q(J,K,I,1:3) = Q(J,K,I,1:3)*uus
          endif
         End Do
       End Do
      End Do

      buf(:,:,:,1) = P 
      buf(:,:,:,2) = Q(:,:,:,4)
      buf(:,:,:,3) = Q(:,:,:,1)
      buf(:,:,:,4) = Q(:,:,:,2)
      buf(:,:,:,5) = Q(:,:,:,3)

      Call Print_Line_X(buf,5)
      Call Print_Line_Y(buf,5)
      Call Print_Line_Z(buf,5)
      Call Print_Diagonal(buf,5)

******************************************************************************
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

         End Do
       End Do
      End Do

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

