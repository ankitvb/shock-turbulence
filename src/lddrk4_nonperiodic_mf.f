************************************************************************
*  Subroutine: Compute_CFL                                             *
*                                                                      *
*  Function:   It computes the CFL number at every grid point and      *
*              selects its maximum value over the whole domain.        *
************************************************************************
      Subroutine Compute_CFL(Q, cfl_stab_all, cfl_resv_all, step)

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), cfl_stab_all, cfl_resv_all
      Real*8  gamma_eff(My,Mz,Mx)

      Common /THERMO/ gamma_eff

*  Globals
      Real*8  W1(My,Mz,Mx), W2(My,Mz,Mx)

      Common  W1, W2 

*  Locals
      Integer I, J, K, ierr, step
      Real*8  fac, tmp1, tmp2, ux, uy, uz
      Real*8  cfl_stab_proc, cfl_resv_proc

      Integer  IDAMAX

      External IDAMAX

*  Compute CFL number at every grid point
      Do I = 1, Mx
        Do K = 1, Mz
          Do J = 1, My
            fac = gamma_eff(J,K,I) * (gamma_eff(J,K,I) - 1.0D0)
            tmp1 = 1.0D0 / Q(J,K,I,4)
            tmp2 = dsqrt( fac * tmp1 * ( Q(J,K,I,5) - 5.0D-1 * tmp1 * 
     .                    ( Q(J,K,I,1) * Q(J,K,I,1) 
     .                    + Q(J,K,I,2) * Q(J,K,I,2) 
     .                    + Q(J,K,I,3) * Q(J,K,I,3) ) ) )
            ux   = dabs(Q(J,K,I,1) * tmp1)
            uy   = dabs(Q(J,K,I,2) * tmp1)
            uz   = dabs(Q(J,K,I,3) * tmp1)

            W1(J,K,I) = ( ( tmp2 + ux ) * kx_max / DX
     .                  + ( tmp2 + uy ) * ky_max / DY
     .                  + ( tmp2 + uz ) * kz_max / DZ )
            W2(J,K,I) = ( ( tmp2 + ux ) * kx_hug / DX
     .                  + ( tmp2 + uy ) * ky_hug / DY
     .                  + ( tmp2 + uz ) * kz_hug / DZ )

          End Do
        End Do
      End Do

*  Find out maximum CFL
      cfl_stab_proc  = maxval(W1)
      cfl_resv_proc  = maxval(W2)

      Call MPI_ALLREDUCE(cfl_stab_proc, cfl_stab_all, 1, MPI_REAL8, 
     .                   MPI_MAX, MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(cfl_resv_proc, cfl_resv_all, 1, MPI_REAL8, 
     .                   MPI_MAX, MPI_COMM_WORLD, ierr)
      
   
      End

************************************************************************
* Subroutine: Compute_Tstep                                            *
*                                                                      *
* Function: Compute time step ensuring stability in presence of        *
*           additional terms from subgrid models                       *
************************************************************************

      Subroutine Compute_Tstep(Q, Y, cfl_resv, step)

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), cfl_resv
      Real*8  Y(My,Mz,Mx,1)
      Real*8  gamma_eff(My,Mz,Mx)

      Common /THERMO/ gamma_eff

*  Globals
      Real*8  W3(My,Mz,Mx), W4(My,Mz,Mx)
      Real*8  W5(My,Mz,Mx), W6(My,Mz,Mx)
      Real*8  T(My,Mz,Mx)

      Common  W3, W4, W5, W6

*  Locals
      Integer I, J, K, ierr, step
      Real*8  t_hvisc, t_hbeta, t_hvisc_all, t_hbeta_all
      Real*8  t_hk, t_hk_all, t_hD, t_hD_all
      Real*8  del, fac, tmp1, tmp2

      Real*8 beta_hyper(My,Mz,Mx), visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx), D_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper, visc_hyper, k_hyper, D_hyper

* Default acoustic timestep if no model
      Tstep = cfl_cal / cfl_resv
    
c      del = (DX*DY*DZ)**(1D0/3D0)
      del = DX

      if((hyper_param .EQ. 1) .AND. (step .GT. 0)) then

* Calculating viscous stable time step due to hyperviscosity

        Do I = 1, Mx
          Do K = 1, Mz
            Do J = 1, My
               fac  = gamma_eff(J,K,I) * (gamma_eff(J,K,I) - 1.0D0)
               fac1 = Y(J,K,I,1) + (epsilon*(gamma1-1D0)/(gamma2-1D0))
     .                * (1D0 - Y(J,K,I,1))
               tmp1 = 1.0D0 / Q(J,K,I,4)
               tmp2 = dsqrt( fac * tmp1 * ( Q(J,K,I,5) - 5.0D-1 * tmp1 *
     .                   ( Q(J,K,I,1) * Q(J,K,I,1)
     .                   + Q(J,K,I,2) * Q(J,K,I,2)
     .                   + Q(J,K,I,3) * Q(J,K,I,3) ) ) )
               T(J,K,I)  = gamma1 * tmp1 *
     .                           ( Q(J,K,I,5) - 5.0D-1 * tmp1 *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )

               T(J,K,I)  = T(J,K,I) / fac1

               W3(J,K,I) = Q(J,K,I,4)*(del**2D0)/visc_hyper(J,K,I)
               W4(J,K,I) = Q(J,K,I,4)*(del**2D0)/beta_hyper(J,K,I)
               W5(J,K,I) = Q(J,K,I,4)*((del*tmp2)**2D0) /
     .                     (k_hyper(J,K,I) * T(J,K,I))
               W6(J,K,I) = (del**2D0)/D_hyper(J,K,I)
            End Do
          End Do
       End Do
     
        t_hvisc = minval(W3)
        t_hbeta = minval(W4)
        t_hk    = minval(W5)
        t_hD    = minval(W6)

        Call MPI_ALLREDUCE(t_hvisc, t_hvisc_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)
        Call MPI_ALLREDUCE(t_hbeta, t_hbeta_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)
        Call MPI_ALLREDUCE(t_hk, t_hk_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)
        Call MPI_ALLREDUCE(t_hD, t_hD_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)

        Tstep = min(cfl_cal/cfl_resv, 0.2D0*t_hvisc_all,
     .              0.2D0*t_hbeta_all, 0.2D0*t_hk_all, 0.2D0*t_hD_all)
      endif

      End

************************************************************************
* Subroutine: Compute_KE                                               *
*                                                                      *
* Function: Computes the average integrated Kinetic energy for the     *
*           fluid field                                                *
***********************************************************************
      subroutine Compute_KE(Q)      

      include 'header'

      real*8, dimension(My,Mz,Mx,5) :: Q
      real*8, dimension(My,Mz,Mx)   :: Q_temp
      real*8, dimension(My,Mz,Mx,3) :: U
      real*8  tmp
*  Locals
      integer I, J, K, L
      Q_temp = 0d0
   
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

c     KE per unit mass     
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz 
            Do L = 1, 3
              Q_temp(J,K,I) = Q_temp(J,K,I) + Q(J,K,I,4)
     .                      * U(J,K,I,L) * U(J,K,I,L) 
            End Do
          End Do
        End Do
      End Do

      call Compute_Sum(Q_temp, KE)

      KE = 0.5d0 * KE /(Nx*Ny*Nz) !* x_leng*y_leng*z_leng

      end
************************************************************************
* Subroutine: Compute_Mom                                              *
*                                                                      *
* Function: Computes the integrated x,y,and z momenta for the          *
*           fluid field                                                *
***********************************************************************
      subroutine Compute_Mom(Q)

      include 'header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(My,Mz,Mx)      :: P

      P =  Q(:,:,:,1)
      call Compute_Sum(P, MomX)
      MomX = MomX / (Nx*Ny*Nz) !* x_leng*y_leng*z_leng

      P =  Q(:,:,:,2)
      call Compute_Sum(P, MomY)
      MomY = MomY / (Nx*Ny*Nz) !* x_leng*y_leng*z_leng

      P =  Q(:,:,:,3)
      call Compute_Sum(P, MomZ)
      MomZ = MomZ / (Nx*Ny*Nz) !* x_leng*y_leng*z_leng  

      end

************************************************************************
* Subroutine: Compute_Mass                                             *
*                                                                      *
* Function: Computes the integrated mass for the                       *
*           fluid field                                                *
***********************************************************************
      subroutine Compute_Mass(Q, Y)

      include 'header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(My,Mz,Mx,1)    :: Y
      real*8, dimension(My,Mz,Mx)      :: P

      P = Q(:,:,:,4)
      call Compute_Sum(P, Mass_t)
      P = Q(:,:,:,4) * Y(:,:,:,1)
      call Compute_Sum(P, Mass_Yt)

      Mass_t = Mass_t / (Nx*Ny*Nz) !* x_leng*y_leng*z_leng
      Mass_Yt = Mass_Yt / (Nx*Ny*Nz)

      end

************************************************************************
* Subroutine: Compute_Energy                                           *
*                                                                      *
* Function: Computes the integrated total energy for the               *
*           fluid field                                                *
************************************************************************
      subroutine Compute_Energy(Q)

      include 'header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(My,Mz,Mx)      :: P

      P =  Q(:,:,:,5)
      call Compute_Sum(P, e_t)
      e_t = e_t /(Nx*Ny*Nz) !* x_leng*y_leng*z_leng

      end

************************************************************************
* Subroutine: Compute_Entropy                                          *
*                                                                      *
* Function: Computes the integrated total entropy for the              *
*           fluid field                                                *
*                                                                      * 
* s = ln(p/rho^gamma)                                                  *
************************************************************************
      Subroutine Compute_Entropy(Q,Y)

      include 'header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(My,Mz,Mx,1)    :: Y
      real*8, dimension(My,Mz,Mx)      :: gamma_eff
      real*8, dimension(My,Mz,Mx)      :: P
      real*8, dimension(My,Mz,Mx)      :: Vs
      real*8, dimension(My,Mz,Mx,2)    :: S
      
      Real*8 fac1

      Common /THERMO/ gamma_eff

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             fac1 = gamma_eff(J,K,I) - 1.0D0
             Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
              S(J,K,I,1) = Q(J,K,I,4)*log(P(J,K,I)*
     .                    (Vs(J,K,I)**gamma_eff(J,K,I)))
              S(J,K,I,2) = Q(J,K,I,4)*S(J,K,I,1)*S(J,K,I,1) 
          End Do
        End Do
      End Do

      Call Compute_Sum(S(1,1,1,1), Entropy)
      Call Compute_Sum(S(1,1,1,2), Entropy2)

      Entropy  = Entropy/(Nx*Ny*Nz)
      Entropy2 = Entropy2/(Nx*Ny*Nz)
 
      Return

      End

************************************************************************
* Subroutine: Compute_gamma                                            *
*                                                                      *
* Function: Compute effective gamma from mass fractions                *
*                                                                      *
*                Y1 gm1          Y2 gm2                                *
*               --------- + eps ---------                              *
*                gm1 - 1         gm2 - 1                               *
*  gamma_eff = --------------------------                              *
*                  Y1              Y2                                  *
*               --------- + eps ---------                              *
*                gm1 - 1         gm2 - 1                               *
************************************************************************
      Subroutine Compute_gamma(Y)

      Include 'header'

* Globals
      Real*8 Y(My,Mz,Mx,1), gamma_eff(My,Mz,Mx)

      Common /THERMO/ gamma_eff

* Locals
      Real*8 fac1, fac2, num, den

      Integer I,J,K

      Do I = 1,Mx
       Do J = 1,My
        Do K = 1,Mz
          fac1 = gamma1 - 1D0
          fac2 = gamma2 - 1D0

          num = ( Y(J,K,I,1) * gamma1 / fac1 )
     .        + ( epsilon * (1D0 - Y(J,K,I,1)) * gamma2 / fac2 )
          den = ( Y(J,K,I,1) / fac1 )
     .        + ( epsilon * (1D0 - Y(J,K,I,1)) / fac2 )

          gamma_eff(J,K,I) = num / den
        End Do
       End Do
      End Do

      End Subroutine Compute_gamma

************************************************************************
*  Subroutine: Compute_rY                                              *
*                                                                      *
*  Function:   It computes rho*Y from Y                                *
************************************************************************
      Subroutine Compute_rY(Q,Y,rY)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), Y(My,Mz,Mx,1), rY(My,Mz,Mx,1)

      Integer I,K,J

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            rY(J,K,I,1) = Q(J,K,I,4) * Y(J,K,I,1)
          End Do
        End Do
      End Do

      End Subroutine Compute_rY


************************************************************************
*  Subroutine: Compute_Y                                               *
*                                                                      *
*  Function:   It computes Y from rho*Y                                *
************************************************************************
      Subroutine Compute_Y(Q,Y,rY)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5), Y(My,Mz,Mx,1), rY(My,Mz,Mx,1)

      Integer I,K,J

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Y(J,K,I,1) = rY(J,K,I,1) / Q(J,K,I,4)
          End Do
        End Do
      End Do

      End Subroutine Compute_Y

************************************************************************
*  Subroutine: LDDRK4                                                  *
*                                                                      *
*  Function:   It advances the flow in time using two-step (5-stage,   *
*              6-stage) low-dissipation-and-dispersion fourth order    *
*              Runge-Kutta scheme.                                     *
************************************************************************
      Subroutine LDDRK4(time,step)

      Include 'header'

      Integer step
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), dQ(My,Mz,Mx,5), QW(My,Mz,Mx,5)
      Real*8  Qf(My,Mz,Mx,5)
      Real*8  Y(My,Mz,Mx,1), rY(My,Mz,Mx,1), rYf(My,Mz,Mx,1)
      Real*8  drYi(My,Mz,Mx,1), rYW(My,Mz,Mx,1)
      Real*8  Yf(My,Mz,Mx,1)
 
      Common  /RKWORK/  Q, dQ, QW
      Common  /SPECIES/ Y, drYi, rYW

*  Locals
      Integer I, J, K, L, M

*  Each Runge-Kutta substep involves the following operations:
*  - update time variable
*  - compute dQ/dt according to the Navier-Stokes equation
*  - time advance solutions
*  - filter fields if necessary
*  - set boundary conditions

*  Step 1: 5-stage RK4

      QW  = 0.0D0
      rYW = 0.0D0

      Do K = 1, N_stage1_rk4

        Call Compute_gamma(Y)
        Call RHS(Q, dQ, Y)
        Call Hyper_RHS(Q, dQ, Y)
        Call Species_RHS(Q, dQ, Y, drYi)
        Call bc_RHS(Q, dQ, Y, drYi, time+c1_rk4(K)*Tstep)
        Call Sponge_Forcing(Q, dQ, Y, drYi)
        Call Compute_rY(Q,Y,rY)

        QW = alpha1_rk4(K) * QW + Tstep * dQ
        Q  = Q + beta1_rk4(K) * QW

        rYW = alpha1_rk4(K) * rYW + Tstep * drYi
        rY  = rY + beta1_rk4(K) * rYW

c        if(step .GT. 10) then
           Call Compact_filter(Q, Qf, 5)
           Call Compact_filter_X_sym(rY(1,1,1,1),rY(1,1,1,1))
           Call Compact_filter_Y_sym(rY(1,1,1,1),rY(1,1,1,1))
           Call Compact_filter_Z_sym(rY(1,1,1,1),rYf(1,1,1,1))
           Q  = Qf   ! Setting field to filtered value
           rY = rYf
c        endif

        Call Compute_Y(Q,Y,rY)
        Call Fix_Pressure(Q, Y)
c        Call Fix_Density(Q)

      End Do

      time = time + Tstep
 
*  Step 2: 6-stage RK4

      QW  = 0.0D0
      rYW = 0.0D0
 
      Do K = 1, N_stage1_rk4

       Call Compute_gamma(Y)
       Call RHS(Q, dQ, Y)
       Call Hyper_RHS(Q, dQ, Y)
       Call Species_RHS(Q, dQ, Y, drYi)
       Call bc_RHS(Q, dQ, Y, drYi, time+c2_rk4(K)*Tstep)
       Call Sponge_Forcing(Q, dQ, Y, drYi)
       Call Compute_rY(Q,Y,rY)

        QW = alpha2_rk4(K) * QW + Tstep * dQ
         Q = Q + beta2_rk4(K)* QW

        rYW = alpha2_rk4(K) * rYW + Tstep * drYi
         rY = rY + beta2_rk4(K)* rYW

* Filter flow variables
c        if(step .GT. 10) then
           Call Compact_filter(Q, Qf, 5)
           Call Compact_filter_X_sym(rY(1,1,1,1),rY(1,1,1,1))
           Call Compact_filter_Y_sym(rY(1,1,1,1),rY(1,1,1,1))
           Call Compact_filter_Z_sym(rY(1,1,1,1),rYf(1,1,1,1))
           Q = Qf   ! Setting field to filtered value
           rY = rYf
c        endif

        Call Compute_Y(Q,Y,rY)
        Call Fix_Pressure(Q, Y)
c        Call Fix_Density(Q)

      End Do

      time = time + Tstep
      step = step + 1

      call Compute_Mom(Q)
      call Compute_KE(Q)
      call Compute_Mass(Q,Y)
      call Compute_Energy(Q)
      call Compute_Entropy(Q,Y)

      Return

      End

