**********************************************************
* Debug stuff                                            *
**********************************************************
      Subroutine Debug(Q, step)   
 
      Include 'header'
      Include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Qf(My,Mz,Mx,5)
      Real*8 Q_aux(My,Mz,Mx,5)
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 S(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3), buf(My,Mz,Mx,4)
      Real*8 Ur(My,Mz,Mx) 
      Real*8 Grad(My,Mz,Mx,3)
      Real*8 fac1, fac2, tmp, tmp1, tmp2, tmp3

      Integer step
      Integer I, J, K, II, M
      Integer ierr, req(4), status(MPI_STATUS_SIZE)

      Real*8 beta_hyper(My,Mz,Mx),visc_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)
      
      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

*  Compute specific volume (avoid division in computations), pressure,
*  temperature and streamwise velocity
*
*  Vs = 1 / rho
*
*  P  = (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*         gamma     P
*  T  = ---------  ---
*       gamma - 1  rho
*

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
             P(J,K,I) = Q(J,K,I,5) * Q(J,K,I,4)**gamma
             S(J,K,I) = P(J,K,I)/Q(J,K,I,4)**gamma
          End Do
        End Do
      End Do

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

      buf(:,:,:,1) = Q(:,:,:,4)
      buf(:,:,:,2) = U(:,:,:,1)
      buf(:,:,:,3) = P
      buf(:,:,:,4) = S

      Call Print_Line_X(buf,4)

      End Subroutine Debug

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

      base_pressure = 1.013D0

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
