************************************************************************
*  Subroutine: Compute_Grad                                            *
*                                                                      *
*  Function:   It computes the dilatation and vorticity fields from a  *
*              solution field of conservative variables.  The spectral *
*              z-derivative in this routine has not been de-aliased,   *
*              hence, may contain aliasing error.  But this aliasing   *
*              error has no effect on the time evolution of the        *
*              solution field.  Output variable Grad contains the      *
*              vorticity and dilatation field variables in the         *
*              following sequence: [ w_x,  w_y,  w_z,  Dil ]           *
************************************************************************
      Subroutine Compute_Grad

      Include 'header'

*  Globals
      Real*8  Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4), temp(My,Mz,Mx)

      Common  /RKWORK/ Q, Grad

      Real*8  U(My,Mz,Mx,3) 
      Real*8  W(My,Mz,Mx)
      Real*8  Omega2(My,Mz,Mx)

      Common  U, Ut_x, dF_x

      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)
      Real*8 dTdx(My,Mz,Mx), dTdy(My,Mz,Mx), dTdz(My,Mz,Mx)
       
      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz,
     .              dTdx, dTdy, dTdz

*  Locals
      Integer I, J, K
      Real*8  tmp

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

*  x-vorticity
*
*        d u_z   d u_y
*  w_x = ----- - -----
*         d y     d z

      Call DFDY(U(1,1,1,3), duz_dy, 0, 1.0D0)
      Call DFDZ(U(1,1,1,2), duy_dz, 0, 1.0D0)

      Grad(:,:,:,1) = duz_dy - duy_dz

*  y-vorticity
*
*        d u_x   d u_z
*  w_y = ----- - -----
*         d z     d x

      Call DFDX(U(1,1,1,3), duz_dx, 0, 1.0D0)
      Call DFDZ(U(1,1,1,1), dux_dz, 0, 1.0D0)

      Grad(:,:,:,2) = dux_dz - duz_dx

*  z-vorticity
*
*        d u_y   d u_x
*  w_z = ----- - -----
*         d x     d y

      Call DFDX(U(1,1,1,2), duy_dx, 0, 1.0D0)
      Call DFDY(U(1,1,1,1), dux_dy, 0, 1.0D0)

      Grad(:,:,:,3) = duy_dx - dux_dy

*  Dilatation
*
*        d u_x   d u_y   d u_z
*  Dil = ----- + ----- + -----
*         d x     d y     d z
 
      Grad(:,:,:,4) = 0D0

      Call DFDX(U(1,1,1,1), dux_dx, 0, 1.0D0)
      Call DFDY(U(1,1,1,2), duy_dy, 0, 1.0D0)
      Call DFDZ(U(1,1,1,3), duz_dz, 0, 1.0D0)

      Grad(:,:,:,4) = dux_dx + duy_dy + duz_dz

      W = dabs(Grad(:,:,:,4))

      call Compute_Max(W, MaxDil)

*  Enstrophy
*
*  Enstrophy = volume_integral(0.5*(w_x^2+w_y^2+w_z^2)*dV)
*
      Omega2 = (Grad(:,:,:,1)**2d0 + Grad(:,:,:,2)**2d0
     .          + Grad(:,:,:,3)**2d0)     

      call Compute_Sum(Omega2, Enstrophy)

      Enstrophy = 0.5D0* Enstrophy/(Nx*Ny*Nz)!*period_x*period_y*period_z

      Return

      End

************************************************************************
* Subroutine: Compute_Gradient
*
* Computes the gradient of a scalar field:
*   dF  dF  dF
* { --, --, -- }
*   dx  dy  dz
*
* In: F: Scalar field
*
* Out: Fgrad: gradient field of F (itself a vector of length 3)
************************************************************************
      Subroutine Compute_Gradient(F, Fgrad)

      Include 'header'

      Real*8 F(My,Mz,Mx), Fgrad(My,Mz,Mx,3)
      Integer I, J, K, L

      Call DFDX(F, Fgrad(1,1,1,1), 0, 1D0)
      Call DFDY(F, Fgrad(1,1,1,2), 0, 1D0)
      Call DFDZ(F, Fgrad(1,1,1,3), 0, 1D0)

      End Subroutine Compute_Gradient
************************************************************************

