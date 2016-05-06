************************************************************
* Routines pertaining to hyperviscosity                    *
************************************************************
* Subroutine Compute_hyper                                *
*                                                         *
* Computes the nonlinear hyperviscosity based on          *
* Cook and Cabot model (JCP, 2004)                        *
***********************************************************
      Subroutine Compute_hyper(Q)

      Include 'header'

      Real*8 Q(My,Mz,Mx,5)
      Real*8 Sij(My,Mz,Mx,9), S(My,Mz,Mx), d2S(My,Mz,Mx)
      Real*8 d4S(My,Mz,Mx), d4Sg(My,Mz,Mx)
      Real*8 d2Y(My,Mz,Mx), d4Y(My,Mz,Mx), d4Yg(My,Mz,Mx)
      Real*8 eta4(My,Mz,Mx,3), beta_hyper(My,Mz,Mx)
      Real*8 eta4g(My,Mz,Mx,3)
      Real*8 visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx), D_hyper(My,Mz,Mx)
      Real*8 C_mu4, C_beta4, C_k4, C_D4, C_Y4

      Real*8  fac1, fac2, fac3
      Real*8  tmp1, tmp2, cs
   
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper,
     .               Sij, eta4, C_mu4, C_beta4

* Coefficients based on Cook and Cabot (JCP, 2004)
      C_mu4   = 2.0D-3
      C_beta4 = 1.0D0
      C_k4    = 1.0D-2
      C_D4    = 1.0D-2
      C_Y4    = 100D0

* Applying biharmonic operator to the strain magnitude and energy 
      Call D2FDX2(Q(1,1,1,1), S, d2S, 0)
      Call D2FDX2(d2S, S, d4S, 0)
      Call D2FDX2(Y(1,1,1,1), S, d2Y, 0)
      Call D2FDX2(d2Y, S, d4Y, 0)

      d4S = abs(d4S)
      d4Y = abs(d4Y)

* Calculating the hyperviscosities
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            cs = dsqrt(gamma_eff(J,K,I) * P(J,K,I) * Vs(J,K,I))
            eta4(J,K,I,1) = d4S(J,K,I)
            eta4(J,K,I,2) = d4Y(J,K,I) 

            if(Y(J,K,I,1) .LT. 1D0) then
               tmp1 = 0D0
            else
               tmp1 = 1D0
            endif
            if(Y(J,K,I,1) .LT. 0D0) then
               tmp2 = 0D0
            else
               tmp2 = 1D0
            endif
            eta4(J,K,I,3) = cs *
     .                    ( ( Y(J,K,I,1) - 1D0 ) * tmp1
     .                    - Y(J,K,I,1) * (1D0 - tmp2) )
          End Do
        End Do
      End Do

* Filtering the field with Gaussian filter (Cook & Cabot, JCP 2004)
      Call Gaussian_filter(eta4, eta4g, 3)

* Calculating the hyperviscosities
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            beta_hyper(J,K,I) = C_beta4 * eta4g(J,K,I,1) * (DX**5D0)
            visc_hyper(J,K,I) = C_mu4   * eta4g(J,K,I,1) * (DX**5D0)
            D_hyper(J,K,I)    = C_D4    * eta4g(J,K,I,2) * (DX**5D0)
     .                        + C_Y4    * eta4g(J,K,I,3) *  DX 
          End Do
        End Do
      End Do

      Return

      End

************************************************************************
*  Subroutine: Hyper_RHS                                               *
*                                                                      *
*  Function:   It computes the right hand side vector (time derivative *
*              of conservative variables) associated with              *
*              hyperviscosity/hyperconductivity/hyperdiffusivity terms *
*              in conservative form.                                   *
*                                                                      *
*  Author: Ankit Bhagatwala                                            *
*                                                                      *
************************************************************************
      Subroutine Hyper_RHS(Q, Qrhs)

      Include 'header'
      Include 'mpif.h'

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5), Qh(My,Mz,Mx,5)

      Common /RHSVH/ Qh

*  Globals
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 durdr(My,Mz,Mx)
      Real*8 dila(My,Mz,Mx)

      Real*8 Flux(My,Mz,Mx), FF(My,Mz,Mx)
      Real*8 dFx(My,Mz,Mx), dFy(My,Mz,Mx), dFz(My,Mz,Mx)
      Real*8 dEx(My,Mz,Mx), dEy(My,Mz,Mx), dEz(My,Mz,Mx)

      Real*8 S(My,Mz,Mx,3), rad(Mx)

* Hyperviscosity vars
      Real*8 visc_hyper(My,Mz,Mx), beta_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

*  Locals
      Integer I, II, J, K, L
      Real*8 IX
      Real*8  fac1, fac2
     
      Integer ierr

      if(hyper_param .EQ. 0)then
         Return
      endif

      Qh = 0D0

* Calcuating hyper coefficient terms
      Call Compute_hyper(Q)

c      visc_hyper(:,:,:) = 1d0
c      beta_hyper(:,:,:) = 1d0
c      k_hyper(:,:,:)    = 0d0

      Do I = 1,Mx
        II     = xrank * Mx + I
        IX     = REAL(II - Nx/2) - 0.5D0
        rad(I) = dabs(IX*DX)
      End Do

*  d u_r          1   d            d u_r
*  ----- = ... + --- ---[ r^2 2 nu ----- ]
*   d t          r^2 d r            d r

       Call DFDX(Q(1,1,1,1), durdr, 0, 1D0)

       Do I = 1,Mx
         Do J = 1,My
           Do K = 1,Mz
             FF(J,K,I) = rad(I) * rad(I) * 2D0 *
     .                   visc_hyper(J,K,I) * durdr(J,K,I)
           End Do
         End Do
       End Do

       Call DFDX(FF, Flux, 0, 1D0)

       Do I = 1,Mx
         Do J = 1,My
           Do K = 1,Mz
             Qh(J,K,I,1) = (1D0/(rad(I)*rad(I))) * Flux(J,K,I)
           End Do
         End Do
       End Do

*  d u_r          1   d                        d u_r   2 u_r
*  ----- = ... + --- ---[ r^2 (beta - 2/3nu) ( ----- + ----- )]
*   d t          r^2 d r                        d r      r

       Do I = 1,Mx
         Do J = 1,My
           Do K = 1,Mz
             FF(J,K,I) = rad(I) * rad(I) *
     .               (beta_hyper(J,K,I) - (2D0/3D0)*visc_hyper(J,K,I))
     .                  * (durdr(J,K,I) + 2D0*Q(J,K,I,1)/rad(I))
           End Do
         End Do
       End Do

       Call DFDX(FF, Flux, 0, 1D0)

       Do I = 1,Mx
         Do J = 1,My
           Do K = 1,Mz
             Qh(J,K,I,1) = Qh(J,K,I,1)
     .                   + (1D0/(rad(I)*rad(I)))*Flux(J,K,I)
           End Do
         End Do
       End Do

* Adding Hyper terms to Navier-Stokes RHS
      Qrhs = Qrhs + Qh

 7    Format(1x, 4(ES13.5))

      Return

      End

*************************************************************************
* Hyper_RHS_test                                                        *
*************************************************************************
      Subroutine Hyper_RHS_test(Q, Qh)

      include 'header'
      include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Qh(My,Mz,Mx,5)
      Real*8 rad(Mx)
      Integer I, J, K

      Do I = 1,Mx
        II     = xrank * Mx + I
        IX     = REAL(II - Nx/2) - 0.5D0
        rad(I) = dabs(IX*DX)
      End Do

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qh(J,K,I,1) = Qh(J,K,I,1)  
     .    - (1D0/3D0)*dsin(rad(I))*(2D0/(rad(I)*rad(I)) - 7D0)
     .    - (16D0/3D0)*(1D0/rad(I))*dcos(rad(I))
          End Do
        End Do
      End Do

      Do M = 0,Px-1
        if(xrank.EQ.M)then
          OPEN(UNIT=16,FILE="x_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do I = 1,Mx
            WRITE(16,7) rad(I), Qh(1,1,I,1)
         End Do
          CLOSE(16)
        endif
        Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      End Do

      Stop

 7    Format(1x, 4(ES13.5))

      Return

      End

