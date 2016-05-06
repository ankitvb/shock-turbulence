************************************************************************
* Subroutine: Set_BC                                                   *
*                                                                      *
* Function: Set the boundary conditions at the end of each             *
*           RK4 substep for the 3D Noh                                 *
*                                                                      *
*           Left boundary : Symmetric boundary                         *
*                                                                      *
*                        d f                                           *
*                  v_n,  --- = 0   for all quantities f                *
*                        d n       & velocity normal to face           *
*                                                                      *
*           Right boundary: Supersonic boundary, all quantities        *
*                           specified                                  *
*                                                                      *
************************************************************************

      Subroutine bc_RHS(Q, Qrhs, time)

      Include 'header'

 
      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8 P(My,Mz,Mx)
      Real*8 time

      Real*8 fac1, fac, rad

      Integer I,J,K,M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ

      base_pressure = (gamma - 1D0) * 1.0D-5 


* Make velocities normal to inner faces zero
* 
* Plane x = 0 ; u_x = 0
* Plane y = 0 ; u_y = 0
* Plane z = 0 ; u_z = 0
*

      Do J = 1,My
       Do K = 1,Mz
        if(xrank .EQ. 0)then
          Q(J,K,1,1) = 0D0
        endif
       End Do
      End Do

      Do I = 1,Mx
       Do J = 1,My
        if(zrank .EQ. 0)then
         Q(J,1,I,3) = 0D0
        endif
       End Do
      End Do

      Do I = 1,Mx
       Do K = 1,Mz
        if(yrank .EQ. 0)then
         Q(1,K,I,2) = 0D0
        endif
       End Do
      End Do
   
* Outer faces, supersonic inflow
* x- boundary
      Do J = 1,My
        Do K = 1,Mz
            if(xrank. EQ. Px-1)then

              JJ = My * yrank + J
              KK = Mz * zrank + K

              IX = REAL(Nx) - 1D0
              IY = REAL(JJ) - 1D0
              IZ = REAL(KK) - 1D0 

              fac = 1D0/dsqrt(IX*IX + IY*IY + IZ*IZ)

              rad = DX / fac

              Q(J,K,Mx,4) = (1D0 + time/rad) ** 2D0

              Q(J,K,Mx,1) = - IX * fac * Q(J,K,Mx,4)
              Q(J,K,Mx,2) = - IY * fac * Q(J,K,Mx,4)
              Q(J,K,Mx,3) = - IZ * fac * Q(J,K,Mx,4)

              P(J,K,Mx)   = base_pressure

              Q(J,K,Mx,5) = ( P(J,K,Mx)/(gamma - 1D0) ) +
     .             ( 5D-1 / Q(J,K,Mx,4) )
     .                   * ( Q(J,K,Mx,1)*Q(J,K,Mx,1)
     .                   +   Q(J,K,Mx,2)*Q(J,K,Mx,2)
     .                   +   Q(J,K,Mx,3)*Q(J,K,Mx,3) )
            endif
 
        End Do
      End Do

* y- boundary
      Do I = 1,Mx
        Do K = 1,Mz
           if(yrank .EQ. Py-1)then

            II = Mx * xrank + I
            KK = Mz * zrank + K

            IX = REAL(II) - 1D0
            IY = REAL(Ny) - 1D0
            IZ = REAL(KK) - 1D0
          
            fac = 1D0/dsqrt(IX*IX + IY*IY + IZ*IZ)
 
            rad = DX / fac

            Q(My,K,I,4) = (1D0 + time/rad) ** 2D0
 
            Q(My,K,I,1) = - IX * fac * Q(My,K,I,4) 
            Q(My,K,I,2) = - IY * fac * Q(My,K,I,4) 
            Q(My,K,I,3) = - IZ * fac * Q(My,K,I,4)

            P(My,K,I)   = base_pressure

            Q(My,K,I,5) = ( P(My,K,I)/(gamma - 1D0) ) +
     .           ( 5D-1 / Q(My,K,I,4) )
     .               * ( Q(My,K,I,1)*Q(My,K,I,1)
     .               +   Q(My,K,I,2)*Q(My,K,I,2)
     .               +   Q(My,K,I,3)*Q(My,K,I,3) )
          endif
        End Do
      End Do
     
* z- boundary
      Do I = 1,Mx
        Do J = 1,My
          if(zrank .EQ. Pz-1)then     

            II = Mx * xrank + I
            JJ = My * yrank + J

            IX = REAL(II) - 1D0
            IY = REAL(JJ) - 1D0
            IZ = REAL(Nz) - 1D0

            fac = 1D0/dsqrt(IX*IX + IY*IY + IZ*IZ)

            rad = DX / fac

            Q(J,Mz,I,4) = (1D0 + time/rad) ** 2D0

            Q(J,Mz,I,1) = - IX * fac * Q(J,Mz,I,4)
            Q(J,Mz,I,2) = - IY * fac * Q(J,Mz,I,4)
            Q(J,Mz,I,3) = - IZ * fac * Q(J,Mz,I,4)

            P(J,Mz,I)   = base_pressure

            Q(J,Mz,I,5) = ( P(J,Mz,I)/(gamma - 1D0) ) +
     .          ( 5D-1 / Q(J,Mz,I,4) )
     .               * ( Q(J,Mz,I,1)*Q(J,Mz,I,1)
     .               +   Q(J,Mz,I,2)*Q(J,Mz,I,2)
     .               +   Q(J,Mz,I,3)*Q(J,Mz,I,3) )
          endif
        End Do
      End Do


      End Subroutine bc_RHS
