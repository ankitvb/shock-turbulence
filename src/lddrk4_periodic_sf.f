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

*  Globals
      Real*8  W1(My,Mz,Mx), W2(My,Mz,Mx)
c      Real*8  S1(40), S2(40) 

      Common  W1, W2 

*  Locals
      Integer I, J, K, ierr, step
      Real*8  fac, tmp1, tmp2, ux, uy, uz
      Real*8  cfl_stab_proc, cfl_resv_proc

      Integer  IDAMAX

      External IDAMAX

      fac = gamma * (gamma - 1.0D0)

*  Compute CFL number at every grid point
      Do I = 1, Mx
        Do K = 1, Mz
          Do J = 1, My
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

      Subroutine Compute_Tstep(Q, cfl_resv, step)

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), cfl_resv

*  Globals
      Real*8  W3(My,Mz,Mx), W4(My,Mz,Mx)
      Real*8  W5(My,Mz,Mx), T(My,Mz,Mx)
c      Real*8  S3(40), S4(40)

      Common  W3, W4, W5

*  Locals
      Integer I, J, K, ierr, step
      Real*8  t_hvisc, t_hbeta, t_hvisc_all, t_hbeta_all
      Real*8  t_hk, t_hk_all
      Real*8  del, fac, tmp1, tmp2

      Real*8 beta_hyper(My,Mz,Mx), visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

* Default acoustic timestep if no model
      Tstep = cfl_cal / cfl_resv
    
c      del = (DX*DY*DZ)**(1D0/3D0)
      del = DX
      fac = gamma * (gamma - 1.0D0)

      if((hyper_param .EQ. 1).AND.(step.GT.1)) then

* Calculating viscous stable time step due to hyperviscosity

        Do I = 1, Mx
          Do K = 1, Mz
            Do J = 1, My
               tmp1 = 1.0D0 / Q(J,K,I,4)
               tmp2 = dsqrt( fac * tmp1 * ( Q(J,K,I,5) - 5.0D-1 * tmp1 *
     .                   ( Q(J,K,I,1) * Q(J,K,I,1)
     .                   + Q(J,K,I,2) * Q(J,K,I,2)
     .                   + Q(J,K,I,3) * Q(J,K,I,3) ) ) )
               T(J,K,I)  = gamma * tmp1 * ( Q(J,K,I,5) - 5.0D-1 * tmp1 *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )

               W3(J,K,I) = Q(J,K,I,4)*(del**2D0)/visc_hyper(J,K,I)
               W4(J,K,I) = Q(J,K,I,4)*(del**2D0)/beta_hyper(J,K,I)
               W5(J,K,I) = Q(J,K,I,4)*((del*tmp2)**2D0) /
     .                     ((k_hyper(J,K,I)) * T(J,K,I))
            End Do
          End Do
       End Do
     
        t_hvisc = minval(W3)
        t_hbeta = minval(W4)
        t_hk    = minval(W5)

        Call MPI_ALLREDUCE(t_hvisc, t_hvisc_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)
        Call MPI_ALLREDUCE(t_hbeta, t_hbeta_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)
        Call MPI_ALLREDUCE(t_hk, t_hk_all, 1, MPI_REAL8,
     .                       MPI_MIN, MPI_COMM_WORLD, ierr)

        Tstep = min(cfl_cal/cfl_resv, 0.2*t_hvisc_all,
     .              0.2*t_hbeta_all,  0.2*t_hk_all )
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
      subroutine Compute_Mass(Q)

      include 'header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(My,Mz,Mx)      :: P

      P =  Q(:,:,:,4)
      call Compute_Sum(P, Mass_t)
      Mass_t = Mass_t / (Nx*Ny*Nz) !* x_leng*y_leng*z_leng

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
      Subroutine Compute_Entropy(Q)

      include 'header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(My,Mz,Mx)      :: P
      real*8, dimension(My,Mz,Mx)      :: Vs
      real*8, dimension(My,Mz,Mx,2)    :: S
      
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
              S(J,K,I,1) = Q(J,K,I,4)*log(P(J,K,I)*(Vs(J,K,I)**gamma))
              S(J,K,I,2) = Q(J,K,I,4)*S(J,K,I,1)*S(J,K,I,1) 
          End Do
        End Do
      End Do

      Call Compute_Sum(S(1,1,1,1), Entropy)
c      Call Compute_Sum(S(1,1,1,2), Entropy2)

      Entropy  = Entropy/(Nx*Ny*Nz)
c      Entropy2 = Entropy2/(Nx*Ny*Nz)
 
      Return

      End

************************************************************************
*  Subroutine: LDDRK4                                                  *
*                                                                      *
*  Function:   It advances the flow in time using two-step (5-stage,   *
*              6-stage) low-dissipation-and-dispersion fourth order    *
*              Runge-Kutta scheme.                                     *
************************************************************************
      Subroutine LDDRK4(time,step)

      Include 'header'
      Include 'part_header'

      Integer step
      Real*8  time

*  Globals
      Real*8  Q(My,Mz,Mx,5), dQ(My,Mz,Mx,5), QW(My,Mz,Mx,5)
      Real*8  Qf(My,Mz,Mx,5)
 
c      Real*8  S1(40), S2(40), S3(40)

      Common  /RKWORK/ Q, dQ, QW

*  Locals
      Integer I, J, K, L, M

*  Each Runge-Kutta substep involves the following operations:
*  - update time variable
*  - compute dQ/dt according to the Navier-Stokes equation
*  - time advance solutions
*  - filter fields if necessary
*  - set boundary conditions

*  Step 1: 5-stage RK4

      QW = 0.0D0
      if(Npart>0) partW = 0D0

      Do K = 1, N_stage1_rk4

        Call RHS(Q, dQ)
        Call Hyper_RHS(Q, dQ)

*   CHECK DQ AFTER FIRST TIME STEP
     
        QW = alpha1_rk4(K) * QW + Tstep * dQ
        Q  = Q + beta1_rk4(K) * QW

* Filter flow variables
        Call Compact_filter(Q, Qf, 5)
        Q = Qf   ! Setting field to filtered value

        Call Fix_Pressure(Q)

* Time stepping the particles
        if(part_param .EQ. 1)then
          Call Part_rhs(Q)
 
          if(Npart>0)then
            partW = alpha2_rk4(K) * partW + Tstep * dpart
            part(:,1:3) = part(:,1:3) + beta2_rk4(K) * partW
          endif

          Call Reorder_Part
        endif
 
      End Do

      time = time + Tstep
 
*  Step 2: 6-stage RK4

      QW = 0.0D0
      if(Npart>0) partW = 0D0
 
      Do K = 1, N_stage2_rk4
  
         Call RHS(Q, dQ)
         Call Hyper_RHS(Q, dQ)

         QW = alpha2_rk4(K) * QW + Tstep * dQ
          Q = Q + beta2_rk4(K)* QW

* Filter flow variables
         Call Compact_filter(Q, Qf, 5)
         Q = Qf   ! Setting field to filtered value
 
        Call Fix_Pressure(Q)

* Time stepping the particles
        if(part_param .EQ. 1)then
          Call Part_rhs(Q)

          if(Npart>0)then
            partW = alpha2_rk4(K) * partW + Tstep * dpart
            part(:,1:3) = part(:,1:3) + beta2_rk4(K) * partW
          endif
  
          Call Reorder_Part
        endif

      End Do

c      if(part_param .EQ. 1)then
c        if(Npart>0)then
c          Do I = 1,Npart
c            if(part(I,NpDof).EQ.1D0)then
c              write(spec_unit,6) time, part(I,4), part(I,5), 
c     .                 part(I,6), part(I,7)
c            endif
c          End Do
c        endif
c      endif


      time = time + Tstep
      step = step + 1

      call Compute_Mom(Q)
      call Compute_KE(Q)
      call Compute_Mass(Q)
      call Compute_Energy(Q)
      call Compute_Entropy(Q)

      Return

 6    Format(1x, 5(ES13.5))

      End

***********************************************************************

      Subroutine Fix_Pressure(Q)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 fac1, fac2, base_pressure

      Integer step
      Integer I,J,K

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      base_pressure = 0.71463D0

* Pressure fix
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )

             if(P(J,K,I) .LT. 0D0) P(J,K,I) = base_pressure

             Q(J,K,I,5) = (1D0/fac1)*P(J,K,I) + 5.0D-1*Vs(J,K,I)*
     .                       ( Q(J,K,I,1) * Q(J,K,I,1)
     .                       + Q(J,K,I,2) * Q(J,K,I,2)
     .                       + Q(J,K,I,3) * Q(J,K,I,3) )

          End Do
        End Do
      End Do

      End Subroutine Fix_Pressure

***********************************************************************

      Subroutine Fix_Density(Q)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5)
      Real*8 base_density

      Integer I,J,K

      base_density = 1D-6

*  Density fix
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             if(Q(J,K,I,4) .LT. 0D0) then
              Q(J,K,I,4) = base_density
             endif
          End Do
        End Do
      End Do

      End Subroutine Fix_Density


