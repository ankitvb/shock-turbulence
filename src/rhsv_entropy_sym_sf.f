************************************************************************
*  Subroutine: RHS                                                     *
*                                                                      *
*  Function:   It computes the right hand side vector (time derivative *
*              of conservative variables) by making calls to           *
*              Subroutine RHS_x and RHS_yz which computes all the      *
*              x-derivatives and yz-derivatives respectively.          *
************************************************************************
      Subroutine RHS(Q, Qrhs)

      Include 'header'

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5), Ux(My,Mz,Mx,5)

      Ux = Q

*  Code up terms involving x-derivatives
      Call RHS_x(Ux, Qrhs) 

*  Code up terms involving y- and z-derivatives
      Call RHS_yz(Q, Ux, Qrhs)

      Return

      End


************************************************************************
*  Subroutine: RHS_x                                                   *
*                                                                      *
*  Function:   It computes all the x-derivatives in the right hand     *
*              side vector.  To reduce the data stride in the          *
*              computation, data come into this subroutine in the form *
*              Q(My,Mz,Mx).                                            *
************************************************************************
      Subroutine RHS_x(Q, Qrhs)

      Include 'header'

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)

*  Globals
      Real*8   P(My,Mz,Mx),   Vs(My,Mz,Mx),  T(My,Mz,Mx)
      Real*8  FF(My,Mz,Mx), Flux(My,Mz,Mx), dF(My,Mz,Mx)

      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)
      Real*8 dTdx(My,Mz,Mx),   dTdy(My,Mz,Mx),   dTdz(My,Mz,Mx)

      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz,
     .              dTdx, dTdy, dTdz

*  Assign blank common statement to save memory requirement
c      Real*8  S1(40), S2(40), S3(40), S4(40), S5(40), S6(40)

      Common  P, Vs, T, FF, Flux, dF

*  Locals
      Integer I, J, K, L
      Real*8  fac1, fac2

      Qrhs = 0D0

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
*  dF = u_x

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Vs(J,K,I) = 1.0D0 / Q(J,K,I,4)
             P(J,K,I) = Q(J,K,I,5) * (Q(J,K,I,4)**gamma)
             T(J,K,I) = fac2 * Vs(J,K,I)  * P(J,K,I)
            dF(J,K,I) = Q(J,K,I,1) * Vs(J,K,I)
          End Do
        End Do
      End Do

*  Convection terms

*  FF = u_x

      Do I = 1, Mx
        Do J = 1, My
          Do K = 1, Mz
            FF(J,K,I) = dF(J,K,I)
          End Do
        End Do
      End Do

*  Flux = rho *  u_x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,4) * FF(J,K,I)
          End Do
        End Do
      End Do
      
*  d rho          d
*  ----- = ... - --- [ rho (u_x) ]
*   d t          d x

      Call DFDX_asym(Flux, Qrhs(1,1,1,4), 1, -1.0D0)

*  Flux = (rho s) u_x
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,5) * FF(J,K,I)
          End Do
        End Do
      End Do
 
*  d (rho s)            d
*  ----------- = ... - --- [ (rho s) u_x  ]
*      d t             d x

      Call DFDX_asym(Flux, Qrhs(1,1,1,5), 1, -1.0D0)

*           
*  FF =  u_x
*         
      Do I = 1, Mx
        Do J = 1, My
          Do K = 1, Mz
            FF(J,K,I) = dF(J,K,I)
          End Do
        End Do
      End Do

*                  
*  Flux = (rho u_x) u_x
*                

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,1) * FF(J,K,I)
          End Do
        End Do
      End Do

*  d (rho u_x)         1   d                   
*  ----------- = ... - -  --- [ (rho u_x) u_x ]
*      d t             2  d x                

      Call DFDX_sym(Flux, Qrhs(1,1,1,1), 1, -5.0D-1)

*  d (rho u_x)         d P
*  ----------- = ... - ---
*      d t             d x

      Call DFDX_sym(P, Qrhs(1,1,1,1), 1, -1.0D0)

                     
*  Flux = (rho u_y) u_x
*                   

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,2) * FF(J,K,I)
          End Do
        End Do
      End Do

*  d (rho u_y)         1  d                     
*  ----------- = ... - - --- [ (rho u_y) u_x ]
*      d t             2 d x                   

      Call DFDX_asym(Flux, Qrhs(1,1,1,2), 1, -5.0D-1)

*  Flux = (rho u_z) u_x
*                        

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,3) * FF(J,K,I)
          End Do
        End Do
      End Do

*  d (rho u_z)         1  d                    
*  ----------- = ... - - --- [ (rho u_z) u_x) ]
*      d t             2 d x                  

      Call DFDX_asym(Flux, Qrhs(1,1,1,3), 1, -5.0D-1)

*       d (rho u_x)
*  FF = -----------
*           d x

      Call DFDX_asym(Q(1,1,1,1), FF, 0, 1.0D0)

*  Convert momenta to velocities, compute viscosity and its derivative

*  Q1 = u_x
*  Q2 = u_y
*  Q3 = u_z
*
*       mu   1
*  Q5 = -- = -- ( (gamma - 1) T ) ** n
*       Re   Re
*
*         1  d  mu 
*  Flux = -- -----
*         Re  d T

      fac2 = vis_exp - 1.0D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Q(J,K,I,1)  = dF(J,K,I)
            Q(J,K,I,2)  = Q(J,K,I,2) * Vs(J,K,I)
            Q(J,K,I,3)  = Q(J,K,I,3) * Vs(J,K,I)
            Q(J,K,I,5)  = ((fac1 * T(J,K,I)) ** vis_exp ) * Re_inv
            Flux(J,K,I) = ((fac1 * T(J,K,I)) ** fac2    ) * Re_inv
     .                  * vis_exp * fac1
          End Do
        End Do
      End Do

*  Continuation of skew-symmetric formulation of x-convection terms:
*
*  d (rho u_x)         1      d
*  ----------- = ... - - u_x --- (rho u_x)
*      d t             2     d x
*
*  d (rho u_y)         1      d
*  ----------- = ... - - u_y --- (rho u_x)
*      d t             2     d x
*
*  d (rho u_z)         1      d
*  ----------- = ... - - u_z --- (rho u_x)
*      d t             2     d x
    
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = -5.0D-1 * FF(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2 * Q(J,K,I,1)
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2 * Q(J,K,I,2)
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2 * Q(J,K,I,3)
          End Do
        End Do
      End Do

*  Conduction terms

*       d T
*  FF = ---
*       d x
*
*       d^2 T
*  dF = -----
*       d x^2

      Call D2FDX2_sym(T, FF, dF, 1)

      dTdx = FF
      
*  d (rho e_t)         1    mu  d^2 T     1  d  mu     d T  
*  ----------- = ... + -- [ --  ----- + ( -- ----- ) ( --- )^2 ]
*      d t             Pr   Re  d x^2     Re  d T      d x
 
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz

         fac1 = Pr_inv 

         Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + fac1 * ( Q(J,K,I,5) * dF(J,K,I)
     .                  + Flux(J,K,I) * FF(J,K,I) * FF(J,K,I) )
          End Do
        End Do
      End Do

*         d u_x
*  Flux = -----
*          d x
*
*         d^2 u_x
*  dF   = -------
*          d x^2

      Call D2FDX2_asym(Q(1,1,1,1), Flux, dF, 1)

      dux_dx = Flux

*  d (rho u_x)              mu   mu_b   4   d^2 u_x
*  ----------- = ... +      -- ( ---- + - ) -------
*      d t                  Re    mu    3    d x^2
*
*  d (rho e_t)              mu   mu_b   4   d^2 u_x
*  ----------- = ... + u_x  -- ( ---- + - ) -------
*      d t                  Re    mu    3    d x^2

      fac1 = 4.0D0 / 3.0D0 + bulk_ratio

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = fac1 * Q(J,K,I,5) * dF(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,1) * fac2

          End Do
        End Do
      End Do

*       d u_y
*  T  = -----
*        d x
*
*        d^2 u_y
*  dF =  -------
*         d x^2

      Call D2FDX2_sym(Q(1,1,1,2), T, dF, 1)

      duy_dx = T

*  d (rho u_y)              mu  d^2 u_y
*  ----------- = ... +      --  -------
*      d t                  Re   d x^2
*
*  d (rho e_t)              mu  d^2 u_y
*  ----------- = ... + u_y  --  -------
*      d t                  Re   d x^2
*
*  Store up x-derivatives for further calls by RHS_yz
*
*       d u_y
*  Q2 = -----
*        d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = Q(J,K,I,5) * dF(J,K,I)
    
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,2) * fac2

            Q(J,K,I,2) = T(J,K,I)
          End Do
        End Do
      End Do

*       d u_z
*  Vs = -----
*        d x
*
*       d^2 u_z
*  dF = -------
*        d x^2

      Call D2FDX2_sym(Q(1,1,1,3), Vs, dF, 1)

      duz_dx = Vs

*  d (rho u_z)              mu  d^2 u_z
*  ----------- = ... +      --  -------
*      d t                  Re   d x^2
*
*  d (rho e_t)              mu  d^2 u_z
*  ----------- = ... + u_z  --  -------
*      d t                  Re   d x^2
*
*  Store up x-derivatives for further calls by RHS_yz
*
*       d u_z
*  Q3 = -----
*        d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = Q(J,K,I,5) * dF(J,K,I)

            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,3) * fac2

            Q(J,K,I,3) = Vs(J,K,I)
          End Do
        End Do
      End Do

*  Continuation of skew-symmetric formulation of x-convection terms:
*
*  d (rho u_x)         1           d u_x
*  ----------- = ... - - (rho u_x) -----
*      d t             2            d x
*
*  d (rho u_y)         1           d u_y
*  ----------- = ... - - (rho u_x) -----
*      d t             2            d x
*
*  d (rho u_z)         1           d u_z
*  ----------- = ... - - (rho u_x) -----
*      d t             2            d x
*
*  Store up x-derivatives for further calls by RHS_yz
*
*       d u_x
*  Q1 = -----
*        d x
*
*        d T
*  Q4 =  ---
*        d x

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = - 5.0D-1 * Q(J,K,I,4) * Q(J,K,I,1)

            Q(J,K,I,1) = Flux(J,K,I)
            Q(J,K,I,4) = FF(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2 * Q(J,K,I,1)
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2 * Q(J,K,I,2)
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2 * Q(J,K,I,3)
          End Do
        End Do
      End Do

      Return

      End


************************************************************************
*  Subroutine: RHS_yz                                                  *
*                                                                      *
*  Function:   It computes all the y- and z-derivatives in the right   *
*              hand side vector.  To reduce the data stride in the     *
*              computation of y-derivatives, data come into this       *
*              subroutine in the form Q(My,Mz,Mx).                     *
************************************************************************
      Subroutine RHS_yz(Q, Ux, Qrhs)

      Include 'header'
     
      Real*8  Q(My,Mz,Mx,5), Ux(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      
*  Globals
      Real*8    P(My,Mz,Mx),   Vs(My,Mz,Mx)
      Real*8    T(My,Mz,Mx),  vis(My,Mz,Mx)
      Real*8   dF(My,Mz,Mx), Flux(My,Mz,Mx)
      Real*8  Tdy(My,Mz,Mx),  Tdz(My,Mz,Mx)
      Real*8  Udz(My,Mz,Mx),  dFz(My,Mz,Mx)
      Real*8  Udx(My,Mz,Mx),  Udy(My,Mz,Mx) 

      Common  P, Vs, T, vis, dF, Flux, Tdz,
     .        Udz, dFz, Udx, Udy

      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)
      Real*8 dTdx(My,Mz,Mx),   dTdy(My,Mz,Mx),   dTdz(My,Mz,Mx)

      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz,
     .              dTdx, dTdy, dTdz
 
*  Locals
      Integer I, J, K, L
      Real*8  fac1, fac2, fac3, dil, fac_X, fac_Y, fac_Z
      Real*8  tmp1, tmp2, tmp3

*  Compute specific volume (avoid division in computations), pressure,  
*  temperature and transverse velocity
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
*   Ux5 = u_y
*
*  Flux = (rho u_x) u_y

      fac1 = gamma - 1.0D0
      fac2 = gamma / fac1

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Vs(J,K,I)   = 1.0D0 / Q(J,K,I,4)
              P(J,K,I)   = Q(J,K,I,5) * (Q(J,K,I,4)**gamma)
              T(J,K,I)   = fac2 * Vs(J,K,I) * P(J,K,I)
             Ux(J,K,I,5) = Q(J,K,I,2) * Vs(J,K,I)
           Flux(J,K,I)   = Q(J,K,I,1) * Ux(J,K,I,5)

          End Do
        End Do
      End Do

* y-convection terms

*       d (rho u_y)
*  dF = -----------
*           d y

      Call DFDY_asym(Q(1,1,1,2), dF, 0, 1.0D0)

*  d rho         d (rho u_y)
*  ----- = ... - -----------
*   d t              d y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,4) = Qrhs(J,K,I,4) - dF(J,K,I)
          End Do
        End Do
      End Do

*  d (rho u_x)         1 d (rho u_x) u_y
*  ----------- = ... - - ---------------
*      d t             2       d y

      Call DFDY_asym(Flux, Qrhs(1,1,1,1), 1, -5.0D-1)

*  Flux = (rho u_y) u_y 

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,2) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho u_y)         1 d (rho u_y) u_y
*  ----------- = ... - - ---------------
*      d t             2       d y

      Call DFDY_sym(Flux, Qrhs(1,1,1,2), 1, -5.0D-1)

*  d (rho u_y)         d P
*  ----------- = ... - ---
*      d t             d y

      Call DFDY_sym(P, Qrhs(1,1,1,2), 1, -1.0D0)

*  Flux = (rho u_z) u_y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,3) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho u_z)         1 d (rho u_z) u_y
*  ----------- = ... - - ---------------
*      d t             2       d y

      Call DFDY_asym(Flux, Qrhs(1,1,1,3), 1, -5.0D-1)

*  Flux = (rho s) u_y

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
           Flux(J,K,I) =  Q(J,K,I,5) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho s)             d (rho s) u_y
*  ----------- = ... - -------------------
*      d t                  d y

      Call DFDY_asym(Flux, Qrhs(1,1,1,5), 1, -1.0D0)

*         d (rho u_z)
*  dFz =  -----------
*            d z

      Call DFDZ_asym(Q(1,1,1,3), dFz, 0, 1.0D0)

*  z-convection terms

*  d rho         d (rho u_z)
*  ----- = ... - -----------
*   d t              d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,4) = Qrhs(J,K,I,4) - dFz(J,K,I)
          End Do
        End Do
      End Do

*  Ux5 = u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Ux(J,K,I,5) = Q(J,K,I,3) * Vs(J,K,I)
          End Do
        End Do
      End Do

*  Flux = (rho u_x) u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,1) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho u_x)         1 d (rho u_x) u_z
*  ----------- = ... - -  ---------------
*      d t             2        d z

      Call DFDZ_asym(Flux, Qrhs(1,1,1,1), 1, -5.0D-1)

*  Flux = (rho u_y) u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,2) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho u_y)         1  d (rho u_y) u_z
*  ----------- = ... - -  ---------------
*      d t             2        d z

      Call DFDZ_asym(Flux, Qrhs(1,1,1,2), 1, -5.0D-1)

*  Flux = (rho u_z) u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             Flux(J,K,I) = Q(J,K,I,3) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho u_z)         1  d (rho u_z) u_z
*  ----------- = ... - -  ---------------
*      d t             2        d z

      Call DFDZ_sym(Flux, Qrhs(1,1,1,3), 1, -5.0D-1)

*  d (rho u_z)         d P
*  ----------- = ... - ---
*      d t             d z

      Call DFDZ_sym(P, Qrhs(1,1,1,3), 1, -1.0D0)

*  Flux = (rho s) u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Flux(J,K,I) = Q(J,K,I,5) * Ux(J,K,I,5)
          End Do
        End Do
      End Do

*  d (rho s)             d (rho s) u_z
*  ----------- = ... - -------------------
*      d t                  d z

      Call DFDZ_asym(Flux, Qrhs(1,1,1,5), 1, -1.0D0)

*  Convert momenta to velocities

*  Q1 = u_x
*  Q2 = u_y
*  Q3 = u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Q(J,K,I,1) = Q(J,K,I,1) * Vs(J,K,I)
            Q(J,K,I,2) = Q(J,K,I,2) * Vs(J,K,I)
            Q(J,K,I,3) = Q(J,K,I,3) * Vs(J,K,I)
          End Do
        End Do
      End Do

*  Continuation of skew-symmetric formulation of y and z-convection terms:
*
*  d (rho u_x)         1      d
*  ----------- = ... - - u_x --- (rho u_y)
*      d t             2     d y
*
*  d (rho u_y)         1      d
*  ----------- = ... - - u_y --- (rho u_y)
*      d t             2     d y
* 
*  d (rho u_z)         1      d
*  ----------- = ... - - u_z --- (rho u_y)
*      d t             2     d y
*
*  d (rho u_x)         1      d
*  ----------- = ... - - u_x --- (rho u_z)
*      d t             2     d z
*
*  d (rho u_y)         1      d
*  ----------- = ... - - u_y --- (rho u_z)
*      d t             2     d z
*
*  d (rho u_z)         1      d
*  ----------- = ... - - u_z --- (rho u_z)
*      d t             2     d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = -5.0D-1*dF(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2 * Q(J,K,I,1)
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2 * Q(J,K,I,2)
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2 * Q(J,K,I,3)
          End Do
        End Do
      End Do

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = -5.0D-1*dFz(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2 * Q(J,K,I,1)
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2 * Q(J,K,I,2)
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2 * Q(J,K,I,3)
          End Do
        End Do
      End Do
      

*  Heat conduction terms
*
*        d T
*  Tdy = ---
*        d y
*
*        d^2 T
*  Ux5 = -----
*        d y^2

      Call D2FDY2_sym(T, Tdy, Ux(1,1,1,5), 1)

      dTdy = Tdy

*         d T
*  Tdz  = ---
*         d z
*
*         d^2 T
*  Flux = -----
*         d z^2

      Call D2FDZ2_sym(T, Tdz, Flux, 1)

      dTdz = Tdz

*  Compute viscosity and its derivative

*        mu   1
*  vis = -- = -- ( (gamma - 1) T ) ** n
*        Re   Re
*
*        1  d  mu 
*  Udz = -- -----
*        Re  d T

      fac2 = vis_exp - 1.0D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            vis(J,K,I) = ((fac1 * T(J,K,I)) ** vis_exp ) * Re_inv
            Udz(J,K,I) = ((fac1 * T(J,K,I)) ** fac2    ) * Re_inv
     .                 * vis_exp * fac1
          End Do
        End Do
      End Do

*  d (rho e_t)         mu  1  d^2 T   mu  1  d^2 T
*  ----------- = ... + --  -- ----- + --  -- ----- 
*      d t             Re  Pr d y^2   Re  Pr d z^2
*
*                        1  d  mu   1    d T
*                    + ( -- ----- ) -- ( --- )^2
*                        Re  d T    Pr   d y
*
*                        1  d  mu   1    d T  
*                    + ( -- ----- ) -- ( --- )^2
*                        Re  d T    Pr   d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = Pr_inv 

            Qrhs(J,K,I,5) =Qrhs(J,K,I,5) + fac2 * ( vis(J,K,I) *
     .                     ( Ux(J,K,I,5) + Flux(J,K,I) )
     .                    + Udz(J,K,I) * ( Tdy(J,K,I) * Tdy(J,K,I)
     .                                   + Tdz(J,K,I) * Tdz(J,K,I) ) )
          End Do
        End Do
      End Do

*        1  d  mu  d T
*  Ux4 = -- -----  ---
*        Re  d T   d x
*
*        1  d  mu  d T
*  Tdy = -- -----  ---
*        Re  d T   d y
*
*        1  d  mu  d T
*  Tdz = -- -----  ---
*        Re  d T   d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
              Ux(J,K,I,4) = Udz(J,K,I) * Ux(J,K,I,4)
              Tdy(J,K,I)  = Udz(J,K,I) * Tdy(J,K,I)
              Tdz(J,K,I)  = Udz(J,K,I) * Tdz(J,K,I)
          End Do
        End Do
      End Do

*        d rho
*  Udz = -----
*         d y

      Call DFDY_sym(Q(1,1,1,4), Udz, 0, 1.0D0)

*        d u_x
*   P  = -----
*         d y
*
*        d^2 u_x
*  Ux5 = -------
*         d y^2

      Call D2FDY2_sym(Q(1,1,1,1), P, Ux(1,1,1,5), 1)

      dux_dy = P

*         d u_y
*  Flux = -----
*          d y
*
*         d^2 u_y
*   dF  = -------
*          d y^2

      Call D2FDY2_asym(Q(1,1,1,2), Flux, dF, 1)

      duy_dy = Flux

*        d u_z
*  Udz = -----
*         d y
*
*        d^2 u_z
*   T  = -------
*         d y^2

      Call D2FDY2_sym(Q(1,1,1,3), Udz, T, 1)

      duz_dy = Udz

*
*         d u_x
*  Udx =  -----
*          d z

      Call DFDZ_sym(Q(1,1,1,1), Udx, 0, 1D0) 

      dux_dz = Udx

*
*         d u_y
*  Udy =  -----
*          d z
   
      Call DFDZ_sym(Q(1,1,1,2), Udy, 0, 1D0)

      duy_dz = Udy

*
*         d u_z
*  dFz =  -----
*          d z

      Call DFDZ_asym(Q(1,1,1,3), dFz, 0, 1D0)      

      duz_dz = dFz

*  Continuation of skew-symmetric formulation of y and z-convection terms:
*
*  d (rho u_x)         1           d u_x
*  ----------- = ... - - (rho u_y) -----
*      d t             2            d y
*
*  d (rho u_y)         1           d u_y
*  ----------- = ... - - (rho u_y) -----
*      d t             2            d y
*
*  d (rho u_z)         1           d u_z
*  ----------- = ... - - (rho u_y) -----
*      d t             2            d y
*
*  d (rho u_x)         1           d u_x
*  ----------- = ... - - (rho u_z) -----
*      d t             2            d z
*
*  d (rho u_y)         1           d u_y
*  ----------- = ... - - (rho u_z) -----
*      d t             2            d z
*
*  d (rho u_z)         1           d u_z
*  ----------- = ... - - (rho u_z) -----
*      d t             2            d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = - 5.0D-1 * Q(J,K,I,4) * Q(J,K,I,2)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2 *    P(J,K,I)
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2 * Flux(J,K,I)
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2 *  Udz(J,K,I)
          End Do
        End Do
      End Do

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac2 = - 5.0D-1 * Q(J,K,I,4) * Q(J,K,I,3)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac2 * Udx(J,K,I)
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac2 * Udy(J,K,I)
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac2 * dFz(J,K,I)
          End Do
        End Do
      End Do


*  Viscous terms

*  d (rho u_x)           mu  d^2 u_x
*  ----------- = ... +   --  -------
*      d t               Re   d y^2
*
*  d (rho u_y)           mu   mu_b   4   d^2 u_y
*  ----------- = ... +   -- ( ---- + - ) -------
*      d t               Re    mu    3    d y^2
*
*  d (rho u_z)           mu  d^2 u_z
*  ----------- = ... +   --  -------
*      d t               Re   d y^2
*
*  d (rho e_t)               mu  d^2 u_x
*  ----------- = ... + u_x [ --  ------- ]
*      d t                   Re   d y^2
*
*                            mu   mu_b   4   d^2 u_y
*                    + u_y [ -- ( ---- + - ) ------- ]
*                            Re    mu    3    d y^2
*
*                            mu  d^2 u_z
*                    + u_z [ --  ------- ]
*                            Re   d y^2


      fac1 = bulk_ratio + 4.0D0 / 3.0D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz

            fac_X =        vis(J,K,I) * Ux(J,K,I,5)
            fac_Y = fac1 * vis(J,K,I) * dF(J,K,I)
            fac_Z =        vis(J,K,I) *  T(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac_X
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac_Y
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac_Z
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,1) * fac_X
     .                                    + Q(J,K,I,2) * fac_Y
     .                                    + Q(J,K,I,3) * fac_Z
          End Do
        End Do
      End Do

*       d u_z
*  dF = -----
*        d z
*
*       d^2 u_z
*  T  = -------
*        d z^2

      Call D2FDZ2_asym(Q(1,1,1,3), dF, T, 0)

      dF = duz_dz

*         d    d u_z
*  Ux5 = --- ( ----- )
*        d y    d z

      Call DFDY_sym(dF, Ux(1,1,1,5), 0, 1.0D0)

*        d u_x   d u_y   d u_z
*  dil = ----- + ----- + -----
*         d x     d y     d z
*
*  d (rho u_x)               1  d  mu d T     d u_x     mu_b   2 
*  ----------- = ... +       -- ----- --- [ 2 ----- + ( ---- - - ) dil ]
*      d t                   Re  d T  d x      d x       mu    3
*
*  d (rho u_y)               1  d  mu d T     d u_y     mu_b   2
*  ----------- = ... +       -- ----- --- [ 2 ----- + ( ---- - - ) dil ]
*      d t                   Re  d T  d y      d y       mu    3
*
*                            mu   mu_b   1    d    d u_z
*                    +       -- ( ---- + - ) --- ( ----- )
*                            Re    mu    3   d y    d z
*
*  d (rho u_z)               1  d  mu d T     d u_z     mu_b   2
*  ----------- = ... +       -- ----- --- [ 2 ----- + ( ---- - - ) dil ]
*      d t                   Re  d T  d z      d z       mu    3
*
*                            mu   mu_b   4   d^2 u_z
*                    +       -- ( ---- + - ) -------
*                            Re    mu    3    d z^2
*
*  d (rho e_t)               1  d  mu d T     d u_x     mu_b   2 
*  ----------- = ... + u_x   -- ----- --- [ 2 ----- + ( ---- - - ) dil ]
*      d t                   Re  d T  d x      d x       mu    3 
*
*                            1  d  mu d T     d u_y     mu_b   2
*                    + u_y { -- ----- --- [ 2 ----- + ( ---- - - ) dil ]
*                            Re  d T  d y      d y       mu    3
*
*                            mu   mu_b   1    d    d u_z
*                          + -- ( ---- + - ) --- ( ----- ) }
*                            Re    mu    3   d y    d z
*
*                            1  d  mu d T     d u_z     mu_b   2
*                    + u_z { -- ----- --- [ 2 ----- + ( ---- - - ) dil ]
*                            Re  d T  d z      d z       mu    3
*
*                              mu   mu_b   4   d^2 u_z
*                            + -- ( ---- + - ) ------- }
*                              Re    mu    3    d z^2
*
*                      mu       d u_x           d u_y
*                    + -- [ 2 ( ----- )^2 + 2 ( ----- )^2 
*                      Re        d x             d y
*
*                               d u_z         mu_b   2
*                         + 2 ( ----- )^2 + ( ---- - - ) dil^2 ]
*                                d z           mu    3

      fac1 = bulk_ratio - 2.0D0 / 3.0D0
      fac2 = bulk_ratio + 4.0D0 / 3.0D0
      fac3 = bulk_ratio + 1.0D0 / 3.0D0

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            tmp1  = Ux(J,K,I,1)
            dil   = tmp1 + Flux(J,K,I) + dF(J,K,I)

            fac_X = Ux(J,K,I,4) * (2.0D0 * tmp1        + fac1 * dil)
            fac_Y = Tdy(J,K,I)  * (2.0D0 * Flux(J,K,I) + fac1 * dil)
     .            + fac3 * vis(J,K,I) * Ux(J,K,I,5)
            fac_Z = Tdz(J,K,I)  * (2.0D0 *   dF(J,K,I) + fac1 * dil)
     .            + fac2 * vis(J,K,I) *  T(J,K,I)

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac_X
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac_Y
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac_Z

            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,1) * fac_X
     .                                    + Q(J,K,I,2) * fac_Y
     .                                    + Q(J,K,I,3) * fac_Z
     .                        + vis(J,K,I) *
     .                          ( 2.0D0 * ( tmp1 * tmp1
     .                                  + Flux(J,K,I) * Flux(J,K,I)
     .                                  +   dF(J,K,I) *   dF(J,K,I) )
     .                                  + fac1 * dil * dil )
          End Do
        End Do
      End Do

*       d u_x
*  T  = -----
*        d z
*
*       d^2 u_x
*  dF = -------
*        d z^2

      Call D2FDZ2_sym(Q(1,1,1,1), T, dF, 0)

      T = dux_dz

*         d u_y
*  Flux = -----
*          d z
*
*         d^2 u_y
*  Ux5  = -------
*          d z^2

      Call D2FDZ2_sym(Q(1,1,1,2), Flux, Ux(1,1,1,5), 0)

      Flux = duy_dz

*  d (rho u_x)               1  d  mu d T   d u_y   d u_x
*  ----------- = ... +       -- ----- --- ( ----- + ----- )
*      d t                   Re  d T  d y    d x     d y
*
*                            1  d  mu d T   d u_x   d u_z
*                    +       -- ----- --- ( ----- + ----- )
*                            Re  d T  d z    d z     d x
*
*                            mu d^2 u_x
*                    +       -- -------
*                            Re  d z^2
*
*  d (rho u_y)               1  d  mu d T   d u_y   d u_x
*  ----------- = ... +       -- ----- --- ( ----- + ----- )
*      d t                   Re  d T  d x    d x     d y
*
*                            1  d  mu d T   d u_z   d u_y
*                    +       -- ----- --- ( ----- + ----- )
*                            Re  d T  d z    d y     d z
*
*                            mu d^2 u_y
*                    +       -- -------
*                            Re  d z^2
*
*  d (rho u_z)               1  d  mu d T   d u_x   d u_z 
*  ----------- = ... +       -- ----- --- ( ----- + ----- )
*      d t                   Re  d T  d x    d z     d x
*
*                            1  d  mu d T   d u_z   d u_y 
*                    +       -- ----- --- ( ----- + ----- )
*                            Re  d T  d y    d y     d z
*
*  d (rho e_t)               1  d  mu d T   d u_y   d u_x
*  ----------- = ... + u_x [ -- ----- --- ( ----- + ----- )
*      d t                   Re  d T  d y    d x     d y
*
*                            1  d  mu d T   d u_x   d u_z
*                          + -- ----- --- ( ----- + ----- )
*                            Re  d T  d z    d z     d x
*
*                            mu d^2 u_x
*                          + -- ------- ]
*                            Re  d z^2
*
*                            1  d  mu d T   d u_y   d u_x
*                    + u_y [ -- ----- --- ( ----- + ----- )
*                            Re  d T  d x    d x     d y
*
*                            1  d  mu d T   d u_z   d u_y 
*                          + -- ----- --- ( ----- + ----- )
*                            Re  d T  d z    d y     d z
*
*                            mu d^2 u_y
*                          + -- ------- ]
*                            Re  d z^2
*
*                            1  d  mu d T   d u_x   d u_z 
*                    + u_z [ -- ----- --- ( ----- + ----- )
*                            Re  d T  d x    d z     d x
*
*                            1  d  mu d T   d u_z   d u_y 
*                          + -- ----- --- ( ----- + ----- ) ]
*                            Re  d T  d y    d y     d z
*
*                      mu     d u_y   d u_x         d u_x   d u_z
*                    + -- [ ( ----- + ----- )^2 + ( ----- + ----- )^2 
*                      Re      d x     d y           d z     d x
*
*                             d u_z   d u_y
*                         + ( ----- + ----- )^2 ]
*                              d y     d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            tmp1  =  Ux(J,K,I,2) +    P(J,K,I)
            tmp2  =  Ux(J,K,I,3) +    T(J,K,I)
            tmp3  = Udz(J,K,I)   + Flux(J,K,I)

            fac_X = Tdy(J,K,I)  * tmp1 + Tdz(J,K,I) * tmp2 + vis(J,K,I)
     .                * dF(J,K,I)
            fac_Y = Ux(J,K,I,4) * tmp1 + Tdz(J,K,I) * tmp3 + vis(J,K,I)
     .                * Ux(J,K,I,5)
            fac_Z = Ux(J,K,I,4) * tmp2 + Tdy(J,K,I) * tmp3

            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac_X
            Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac_Y
            Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac_Z
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,1) * fac_X
     .                                    + Q(J,K,I,2) * fac_Y
     .                                    + Q(J,K,I,3) * fac_Z
     .                        + vis(J,K,I) * ( tmp1 * tmp1 + tmp2 * tmp2
     .                                       + tmp3 * tmp3 )
          End Do
        End Do
      End Do

*         d    d u_x
*  Ux4 = --- ( ----- )
*        d y    d x

      Call DFDY_sym(Ux(1,1,1,1), Ux(1,1,1,4), 0, 1.0D0)

*         d    d u_x
*  Ux5 = --- ( ----- )
*        d z    d x

      Call DFDZ_sym(Ux(1,1,1,1), Ux(1,1,1,5), 0, 1.0D0)

*         d    d u_y
*  Tdy = --- ( ----- )
*        d y    d x

      Call DFDY_asym(Ux(1,1,1,2), Tdy, 0, 1.0D0)

*         d    d u_z
*  Tdz = --- ( ----- )
*        d z    d x

      Call DFDZ_asym(Ux(1,1,1,3), Tdz, 0, 1.0D0)

*         d    d u_y
*   T  = --- ( ----- )
*        d y    d z

      Call DFDY_asym(Flux, T, 0, 1.0D0)

*  d (rho u_x)           mu   mu_b   1      d    d u_y
*  ----------- = ... +   -- ( ---- + - ) [ --- ( ----- )
*      d t               Re    mu    3     d y    d x
*
*                                           d    d u_z 
*                                        + --- ( ----- ) ]
*                                          d z    d x
*
*  d (rho u_y)           mu   mu_b   1    d    d u_x
*  ----------- = ... +   -- ( ---- + - ) --- ( ----- )
*      d t               Re    mu    3   d y    d x
*
*  d (rho u_z)           mu   mu_b   1      d    d u_x
*  ----------- = ... +   -- ( ---- + - ) [ --- ( ----- )
*      d t               Re    mu    3     d z    d x
*
*                                           d    d u_y 
*                                        + --- ( ----- ) ]
*                                          d y    d z
*
*  d (rho e_t)                mu   mu_b   1      d    d u_y
*  ----------- = ... + u_x { -- ( ---- + - ) [ --- ( ----- )
*      d t                    Re    mu    3     d y    d x
*
*                                               d    d u_z 
*                                            + --- ( ----- ) ] }
*                                              d z    d x
*
*                            mu   mu_b   1    d    d u_x
*                    + u_y { -- ( ---- + - ) --- ( ----- ) }
*                            Re    mu    3   d y    d x
*
*                            mu   mu_b   1      d    d u_x
*                    + u_z { -- ( ---- + - ) [ --- ( ----- )
*                            Re    mu    3     d z    d x
*
*                                               d    d u_y 
*                                            + --- ( ----- ) ] }
*                                              d y    d z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
             fac_X = fac3 * vis(J,K,I) * ( Tdy(J,K,I)  + Tdz(J,K,I) )
             fac_Y = fac3 * vis(J,K,I) *   Ux(J,K,I,4)
             fac_Z = fac3 * vis(J,K,I) * ( Ux(J,K,I,5) +   T(J,K,I) )

             Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + fac_X
             Qrhs(J,K,I,2) = Qrhs(J,K,I,2) + fac_Y
             Qrhs(J,K,I,3) = Qrhs(J,K,I,3) + fac_Z
             Qrhs(J,K,I,5) = Qrhs(J,K,I,5) + Q(J,K,I,1) * fac_X
     .                                     + Q(J,K,I,2) * fac_Y
     .                                     + Q(J,K,I,3) * fac_Z

          End Do
        End Do
      End Do

*  Convert velocities to momenta

*  Q1 = rho u_x
*  Q2 = rho u_y
*  Q3 = rho u_z

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Q(J,K,I,1) = Q(J,K,I,1) * Q(J,K,I,4)
            Q(J,K,I,2) = Q(J,K,I,2) * Q(J,K,I,4)
            Q(J,K,I,3) = Q(J,K,I,3) * Q(J,K,I,4)
          End Do
        End Do
      End Do

      Return

      End
