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

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5), Qh(My,Mz,Mx,5)

*  Globals
      Real*8 Vs(My,Mz,Mx), P(My,Mz,Mx), T(My,Mz,Mx)
      Real*8 U(My,Mz,Mx,3)
      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)
      Real*8 dTdx(My,Mz,Mx), dTdy(My,Mz,Mx), dTdz(My,Mz,Mx)
      Real*8 dila(My,Mz,Mx)

      Real*8 Flux(My,Mz,Mx), FF(My,Mz,Mx) 
      Real*8 dFx(My,Mz,Mx), dFy(My,Mz,Mx), dFz(My,Mz,Mx)
      Real*8 dEx(My,Mz,Mx), dEy(My,Mz,Mx), dEz(My,Mz,Mx)
      
      Real*8 S(My,Mz,Mx,3)

      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz,
     .              dTdx, dTdy, dTdz

* Hyperviscosity vars
      Real*8 visc_hyper(My,Mz,Mx), beta_hyper(My,Mz,Mx)
      Real*8 k_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper, visc_hyper, k_hyper

*  Locals
      Integer I, J, K, L
      Real*8  fac1, fac2

      if(hyper_param .EQ. 0)then
         Return
      endif
 
      Qh = 0D0

* Calcuating hyper coefficient terms
      Call Compute_hyper(Q)

c      visc_hyper(:,:,:) = 1d0
c      beta_hyper(:,:,:) = 0d0
c      k_hyper(:,:,:)    = 0d0

*  U1 = u_x
*  U2 = u_y
*  U3 = u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            U(J,K,I,1) = Q(J,K,I,1) / Q(J,K,I,4)
            U(J,K,I,2) = Q(J,K,I,2) / Q(J,K,I,4)
            U(J,K,I,3) = Q(J,K,I,3) / Q(J,K,I,4)
          End Do
        End Do
      End Do

* Computing dilatation and shear strains once and for all

*         d u_x    d u_y   d u_z 
* dila =  ----- +  ----- + -----
*          d x      d y     d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             dila(J,K,I) = dux_dx(J,K,I) + duy_dy(J,K,I) + duz_dz(J,K,I)
          End Do
        End Do
      End Do

*       d u_x   d u_y 
* S1 =  ----- + ----- 
*        d y     d x   

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            S(J,K,I,1) = dux_dy(J,K,I) + duy_dx(J,K,I)
          End Do
        End Do
      End Do

*       d u_x   d u_z
* S2 =  ----- + -----
*        d z     d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            S(J,K,I,2) = dux_dz(J,K,I) + duz_dx(J,K,I)
          End Do
        End Do
      End Do

*       d u_y   d u_z
* S3 =  ----- + -----
*        d z     d y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            S(J,K,I,3) = duy_dz(J,K,I) + duz_dy(J,K,I)
          End Do
        End Do
      End Do

*  Compute specific volume (avoid division in computations), pressure,
*  temperature 
*
*   Vs  = 1 / rho
*
*    P  = (gamma - 1) *
*
*                     1
*     { (rho e_t) - -----  [ (rho u_x)^2 + (rho u_y)^2 + (rho u_z)^2 ] }
*                   2 rho
*
*           gamma     P
*    T  = ---------  ---
*         gamma - 1  rho
*

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
              T(J,K,I)   = fac2 * Vs(J,K,I) * P(J,K,I)
          End Do
        End Do
      End Do

* Computing temperature gradients once and for all
      Call DFDX_sym(T, dTdx, 0, 1.0D0)
      Call DFDY_sym(T, dTdy, 0, 1.0D0)
      Call DFDZ_sym(T, dTdz, 0, 1.0D0)

*                 d u_x              2         d u_x    d u_y   d u_z
* Flux = mu_h ( 2 ----- ) + (mu_bh - - mu_h) ( ----- +  ----- + ----- )
*                  d x               3          d x      d y     d z

*                     d u_x              2        d u_x    d u_y   d u_z
* FF = u_x ( mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) )
*                      d x               3         d x      d y     d z

* d (rho u_x)          d            d u_x              2        d u_x    d u_y   d u_z
* ----------- = ... + --- (mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) )
*     d t             d x            d x               3         d x      d y     d z

* d (rho e_t)          d                   d u_x              2        d u_x    d u_y   d u_z
* ----------- = ... + --- ( u_x ( mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) ) )
*     d t             d x                   d x               3         d x      d y     d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = 2.0D0 * visc_hyper(J,K,I) * dux_dx(J,K,I)
     .   + (beta_hyper(J,K,I) - (2D0/3D0)*visc_hyper(J,K,I))*dila(J,K,I)
            FF(J,K,I) = U(J,K,I,1) * Flux(J,K,I)
          End Do
        End Do
      End Do
   
      Call DFDX_sym(Flux, dFx, 0, 1.0D0)   
      Call DFDX_asym(FF, dEx, 0, 1.0D0)

*               d u_x   d u_y 
* Flux = mu_h ( ----- + ----- )
*                d y     d x

*                   d u_x   d u_y
* FF = u_x ( mu_h ( ----- + ----- ) )
*                    d y     d x

* d (rho u_x)          d          d u_x   d u_y
* ----------- = ... + ---( mu_h ( ----- + ----- ) )
*     d t             d y          d y     d x

* d (rho e_t)          d                 d u_x   d u_y
* ----------- = ... + --- ( u_x ( mu_h ( ----- + ----- ) ) )
*     d t             d y                 d y     d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = visc_hyper(J,K,I) * S(J,K,I,1)
              FF(J,K,I) = U(J,K,I,1) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDY_asym(Flux, dFy, 0, 1.0D0)
      Call DFDY_asym(FF,   dEy, 0, 1.0D0)

*               d u_x   d u_z
* Flux = mu_h ( ----- + ----- )
*                d z     d x

*                   d u_x   d u_z
* FF = u_x ( mu_h ( ----- + ----- ) )
*                    d z     d x

* d (rho u_x)           d           d u_x   d u_z
* ----------- = ... +  --- ( mu_h ( ----- + ----- ) )
*    d t               d z           d z     d x

* d (rho e_t)          d                 d u_x   d u_z
* ----------- = ... + --- ( u_x ( mu_h ( ----- + ----- ) ) )
*     d t             d z                 d z     d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = visc_hyper(J,K,I) * S(J,K,I,2)
              FF(J,K,I) = U(J,K,I,1) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDZ_asym(Flux, dFz, 0, 1.0D0)
      Call DFDZ_asym(FF,   dEz, 0, 1.0D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Qh(J,K,I,1) = Qh(J,K,I,1) + dFx(J,K,I)
     .                   + dFy(J,K,I)  + dFz(J,K,I)
             Qh(J,K,I,5) = Qh(J,K,I,5) + dEx(J,K,I)
     .                   + dEy(J,K,I)  + dEz(J,K,I)
          End Do
        End Do
      End Do

*                 d u_y              2        d u_x    d u_y   d u_z
* Flux = mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- )
*                  d y               3         d x      d y     d z

*                     d u_y              2        d u_x    d u_y   d u_z
* FF = u_y ( mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) )
*                      d y               3         d x      d y     d z

* d (rho u_y)          d            d u_y              2        d u_x    d u_y   d u_z
* ----------- = ... + --- (mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ))
*     d t             d y            d y               3         d x      d y     d z

* d (rho e_t)          d                   d u_y              2        d u_x    d u_y   d u_z
* ----------- = ... + --- ( u_y ( mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) ) )
*     d t             d y                   d y               3         d x      d y     d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = 2.0D0 * visc_hyper(J,K,I) * duy_dy(J,K,I)
     .   + (beta_hyper(J,K,I) - (2D0/3D0)*visc_hyper(J,K,I))*dila(J,K,I)
            FF(J,K,I) = U(J,K,I,2) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDY_sym(Flux, dFy, 0, 1.0D0)
      Call DFDY_asym(FF,   dEy, 0, 1.0D0)

*               d u_x   d u_y
* Flux = mu_h ( ----- + ----- )
*                d y     d x

*                   d u_x   d u_y
* FF = u_y ( mu_h ( ----- + ----- ) )
*                    d y     d x

* d (rho u_y)          d           d u_x   d u_y
* ----------- = ... + --- ( mu_h ( ----- + ----- ) )
*     d t             d x           d y     d x

* d (rho e_t)          d                d u_x   d u_y
* ----------- = ... + ---( u_y ( mu_h ( ----- + ----- ) ) )
*    d t              d x                d y     d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = visc_hyper(J,K,I) * S(J,K,I,1)
               FF(J,K,I) = U(J,K,I,2) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDX_asym(Flux, dFx, 0, 1.0D0)
      Call DFDX_asym(FF, dEx, 0, 1.0D0)
      
*               d u_y   d u_z
* Flux = mu_h ( ----- + ----- )
*                d z     d y

*                   d u_y   d u_z
* FF = u_y ( mu_h ( ----- + ----- ) )
*                    d z     d y

* d (rho u_y)          d           d u_y   d u_z
* ----------- = ... + --- ( mu_h ( ----- + ----- ) )
*     d t             d z           d z     d y

* d (rho e_t)          d                d u_y   d u_z
* ----------- = ... + ---( u_y ( mu_h ( ----- + ----- ) ) )
*    d t              d z                d z     d y


      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = visc_hyper(J,K,I) * S(J,K,I,3)
               FF(J,K,I) = U(J,K,I,2) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDZ_asym(Flux, dFz, 0, 1.0D0)
      Call DFDZ_asym(FF,   dEz, 0, 1.0D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Qh(J,K,I,2) = Qh(J,K,I,2) + dFx(J,K,I)
     .                   + dFy(J,K,I)  + dFz(J,K,I)
             Qh(J,K,I,5) = Qh(J,K,I,5) + dEx(J,K,I)
     .                   + dEy(J,K,I)  + dEz(J,K,I)
          End Do
        End Do
      End Do

*                 d u_z              2        d u_x    d u_y   d u_z
* Flux = mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- )
*                  d z               3         d x      d y     d z

*                     d u_z              2        d u_x    d u_y   d u_z
* FF = u_z ( mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) )
*                      d z               3         d x      d y     d z

* d (rho u_z)          d            d u_z              2        d u_x    d u_y   d u_z
* ----------- = ... + --- (mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ))
*     d t             d z            d z               3         d x      d y     d z

* d (rho e_t)          d                   d u_z              2        d u_x    d u_y   d u_z
* ----------- = ... + --- ( u_z ( mu_h ( 2 ----- ) + (mu_bh - - mu_h) (----- +  ----- + ----- ) ) )
*     d t             d z                   d z               3         d x      d y     d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = 2.0D0 * visc_hyper(J,K,I) * duz_dz(J,K,I)
     .   + (beta_hyper(J,K,I) - (2D0/3D0)*visc_hyper(J,K,I))*dila(J,K,I)
            FF(J,K,I) = U(J,K,I,3) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDZ_sym(Flux, dFz, 0, 1.0D0)
      Call DFDZ_asym(FF,   dEz, 0, 1.0D0)

*               d u_z   d u_x
* Flux = mu_h ( ----- + ----- )
*                d x     d z

*                   d u_z   d u_x
* FF = u_z ( mu_h ( ----- + ----- ) )
*                    d x     d z

* d (rho u_z)          d           d u_z   d u_x
* ----------- = ... + --- ( mu_h ( ----- + ----- ) )
*     d t             d x           d x     d z

* d (rho e_t)          d                d u_z   d u_x
* ----------- = ... + ---( u_z ( mu_h ( ----- + ----- ) ) )
*    d t              d x                d x     d z


      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = visc_hyper(J,K,I) * S(J,K,I,2)
               FF(J,K,I) = U(J,K,I,3) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDX_asym(Flux, dFx, 0, 1.0D0)
      Call DFDX_asym(FF, dEx, 0, 1.0D0)

*               d u_y   d u_z
* Flux = mu_h ( ----- + ----- )
*                d z     d y

*                   d u_y   d u_z
* FF = u_z ( mu_h ( ----- + ----- ) )
*                    d z     d y

* d (rho u_z)          d           d u_y   d u_z
* ----------- = ... + --- ( mu_h ( ----- + ----- ) )
*     d t             d y           d z     d y

* d (rho e_t)          d                d u_y   d u_z
* ----------- = ... + ---( u_z ( mu_h ( ----- + ----- ) ) )
*    d t              d y                d z     d y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
              Flux(J,K,I) = visc_hyper(J,K,I) * S(J,K,I,3)
                FF(J,K,I) = U(J,K,I,3) * Flux(J,K,I)
          End Do
        End Do
      End Do

      Call DFDY_asym(Flux, dFy, 0, 1.0D0)
      Call DFDY_asym(FF,   dEy, 0, 1.0D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Qh(J,K,I,3) = Qh(J,K,I,3) + dFx(J,K,I)
     .                   + dFy(J,K,I)  + dFz(J,K,I)
             Qh(J,K,I,5) = Qh(J,K,I,5) + dEx(J,K,I)
     .                   + dEy(J,K,I)  + dEz(J,K,I)
          End Do
        End Do
      End Do

* Heat conduction terms
 
* d (rho e_t)          d     d T      d     d T      d     d T                
* ----------- = ... + ---( k --- ) + ---( k --- ) + ---( k --- )                      
*    d t              d x    d x     d y    d y     d z    d z             

      Do I = 1,Mx 
        Do J = 1,My
          Do K = 1,Mz
              Flux(J,K,I) = k_hyper(J,K,I) * dTdx(J,K,I)
          End Do
        End Do
      End Do

      Call DFDX_sym(Flux, dFx, 0, 1.0D0)
      
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
              Flux(J,K,I) = k_hyper(J,K,I) * dTdy(J,K,I)
          End Do
        End Do
      End Do
 
      Call DFDY_sym(Flux, dFy, 0, 1.0D0)
      
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
              Flux(J,K,I) = k_hyper(J,K,I) * dTdz(J,K,I)
          End Do
        End Do
      End Do

      Call DFDZ_sym(Flux, dFz, 0, 1.0D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Qh(J,K,I,5) = Qh(J,K,I,5) + dFx(J,K,I) + dFy(J,K,I)
     .                   + dFz(J,K,I)
          End Do      
        End Do
      End Do

* Adding Hyper terms to Navier-Stokes RHS
      Qrhs = Qrhs + Qh

      Return


      End

