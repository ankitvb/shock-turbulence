*******************************************************************
* Subroutine: Species_RHS                                         *                                                 
*                                                                 *
* Function: Computes the right hand side of the species           *
*           evolution equation                                    *
*                                                                 *
* Author: Ankit Bhagatwala                                        *
*******************************************************************        
      Subroutine Species_RHS(Q, Qrhs, Y, Yrhs)

      Include 'header'

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8  Y(My,Mz,Mx,1), Yrhs(My,Mz,Mx,1)
      Real*8  Yconvx(My,Mz,Mx,1), Yconvy(My,Mz,Mx,1)
      Real*8  Yconvz(My,Mz,Mx,1), Ydiff(My,Mz,Mx,1)

      Common /SPECIES_RHSV/ Yconvx, Yconvy, Yconvz, Ydiff

      Real*8  Flux(My,Mz,Mx), FF(My,Mz,Mx)
      
      Real*8  dYdx(My,Mz,Mx), dYdy(My,Mz,Mx), dYdz(My,Mz,Mx)

      Common /SPECIES_DERS/ dYdx, dYdy, dYdz

      Real*8 gamma_eff(My,Mz,Mx)
      Common /THERMO/ gamma_eff

* Hyperviscosity vars
      Real*8 visc_hyper(My,Mz,Mx), beta_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx), D_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper, visc_hyper, k_hyper, D_hyper

      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 fac1, fac2, fac3
  
*  Locals
      Integer I, J, K, L

      if(hyper_param .EQ. 0)then
         Return
      endif

      Yrhs = 0D0; Yconvx = 0D0; Yconvy = 0D0; Yconvz = 0D0
      Ydiff = 0D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac1 = gamma_eff(J,K,I) - 1.0D0
            fac2 = gamma_eff(J,K,I) / fac1
            fac3 = gamma1 / (gamma1 - 1D0)
            fac3 = fac3 *
     .      1D0 / ( Y(J,K,I,1) + epsilon * (1D0 - Y(J,K,I,1)) )

            Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
             P(J,K,I)   = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                        ( Q(J,K,I,1) * Q(J,K,I,1)
     .                        + Q(J,K,I,2) * Q(J,K,I,2)
     .                        + Q(J,K,I,3) * Q(J,K,I,3) ) )
             T(J,K,I)   = fac3 * Vs(J,K,I)  * P(J,K,I)
          End Do
        End Do
      End Do

  
* Precomputing mass fraction derivatives to be used later     
      Call DFDX_sym(Y(1,1,1,1), dYdx, 0, 1.0D0)
      Call DFDY_sym(Y(1,1,1,1), dYdy, 0, 1.0D0)
      Call DFDZ_sym(Y(1,1,1,1), dYdz, 0, 1.0D0)

* 
* Flux = (rho u_x) Y_i
*
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,1) * Y(J,K,I,1) 
          End Do
        End Do
      End Do

* d (rho Y_i)          d 
* ----------- = ... - ---[ (rho u_x) Y_i]
*     d t             d x

      Call DFDX_asym(Flux, Yconvx(1,1,1,1), 1, -1.0D0)

*
* Flux = (rho u_y) Y_i
*
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,2) * Y(J,K,I,1)
          End Do
        End Do
      End Do

* d (rho Y_i)          d
* ----------- = ... - ---[ (rho u_y) Y_i]
*     d t             d y

      Call DFDY_asym(Flux, Yconvy(1,1,1,1), 1, -1.0D0)

*
* Flux = (rho u_z) Y_i
*
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,3) * Y(J,K,I,1)
          End Do
        End Do
      End Do

* d (rho Y_i)          d
* ----------- = ... - ---[ (rho u_z) Y_i]
*     d t             d z

      Call DFDZ_asym(Flux, Yconvz(1,1,1,1), 1, -1.0D0)

*                    d Y_i             d Y_1        d Y_2
* Flux =   rho ( D_i ----- - Y_i [ D_1 -----  + D_2 ----- ] )
*                     d x               d x          d x

*  d (rho Y_i)          d              d Y_i             d Y_1        d Y_2
* -----------  = ... + --- ( rho ( D_i ----- - Y_i [ D_1 -----  + D_2 ----- ] ) )
*      d t             d x              d x               d x          d x

* d (rho E_t)          d
* ----------- = ... - ---( h_1 J_1 + h_2 J_2 ) 
*      d t            d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,4) * D_hyper(J,K,I) * dYdx(J,K,I)
     .                  + Q(J,K,I,4) * dYdx(J,K,I) * Re1_inv * Sc_inv
            FF(J,K,I) = 
     .      (1D0 - epsilon*gamma2*(gamma1-1D0)/(gamma1*(gamma2-1D0)))
     .              * T(J,K,I) * Flux(J,K,I) 
          End Do
        End Do
      End Do 
 
      Call DFDX_asym(Flux, Ydiff(1,1,1,1), 1, 1.0D0)
      Call DFDX_asym(FF, Qrhs(1,1,1,5), 1, 1.0D0)

*                    d Y_i             d Y_1        d Y_2
* Flux =   rho ( D_i ----- - Y_i [ D_1 -----  + D_2 ----- ] )
*                     d y               d y          d y

*  d (rho Y_i)          d              d Y_i             d Y_1        d Y_2
* -----------  = ... + --- ( rho ( D_i ----- - Y_i [ D_1 -----  + D_2 ----- ] ) )
*      d t             d y              d y               d y          d y

* d (rho E_t)          d
* ----------- = ... - ---( h_1 J_1 + h_2 J_2 )
*      d t            d y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,4) * D_hyper(J,K,I) * dYdy(J,K,I)
     .                  + Q(J,K,I,4) * dYdy(J,K,I) * Re1_inv * Sc_inv
            FF(J,K,I) = 
     .      (1D0 - epsilon*gamma2*(gamma1-1D0)/(gamma1*(gamma2-1D0)))
     .              * T(J,K,I) * Flux(J,K,I)
          End Do
        End Do
      End Do
     
      Call DFDY_asym(Flux, Ydiff(1,1,1,1), 1, 1.0D0)
      Call DFDY_asym(FF, Qrhs(1,1,1,5), 1, 1.0D0)

*                    d Y_i             d Y_1        d Y_2
* Flux =   rho ( D_i ----- - Y_i [ D_1 -----  + D_2 ----- ] )
*                     d z               d z          d z

*  d (rho Y_i)          d              d Y_i             d Y_1        d Y_2
* -----------  = ... + --- ( rho ( D_i ----- - Y_i [ D_1 -----  + D_2 ----- ] ) )
*      d t             d z              d z               d z          d z

* d (rho E_t)          d
* ----------- = ... - ---( h_1 J_1 + h_2 J_2 )
*      d t            d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,4) * D_hyper(J,K,I) * dYdz(J,K,I)
     .                  + Q(J,K,I,4) * dYdz(J,K,I) * Re1_inv * Sc_inv
            FF(J,K,I) = 
     .      (1D0 - epsilon*gamma2*(gamma1-1D0)/(gamma1*(gamma2-1D0)))
     .              * T(J,K,I) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDZ_asym(Flux, Ydiff(1,1,1,1), 1, 1.0D0)
      Call DFDZ_asym(FF, Qrhs(1,1,1,5), 1, 1.0D0)

      Yrhs = Yconvx + Yconvy + Yconvz + Ydiff

      End Subroutine Species_RHS
