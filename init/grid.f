************************************************************************
*  Subroutine: Setup_Grid                                              *
*                                                                      *
*  Function:   It generates a computational grid whose nodes are       *
*              uniformly spaced.  It is responsible for computing and  *
*              saving the metric terms in global storage which will be *
*              used by the derivative routines.                        *
*              It also outputs the computational grid into formatted   *
*              files.                                                  *
*                                                                      *
*  Output (metric terms):                                              *
*  s1 = 1 / (dx_i/ds)                                                  *
*  s2 = s1 * s1                                                        *
************************************************************************
       Subroutine Setup_Grid

       Include 'header'
       Include 'mpif.h'

      Integer Leng_Char

*  Locals
      Integer I, J, K, ierr

*  Grid generation
*
*  x-grid:
      Call Uniform_X(Nx, x_leng, gxcor)

*  Compute metric terms for x-grid
      sx(1) = 1d0 / x_leng
      sx(2) = sx(1) * sx(1)

*  Assign xcor values to corresponding processors
      Do I = 1, Mx
        xcor(I) = gxcor(xrank*Mx+I)
      End Do

*  Compute grid spacings for x-grid
      DX = gxcor(2) - gxcor(1)

*  Output x-grid information
      If (rank .EQ. 0) Then
        Open(unit = tmp_unit,
     .       file = prefix(1:Leng_Char(prefix))//'.grid.x')

        Do I = 1, Nx
          Write(tmp_unit,10) I, gxcor(I), sx(1), sx(2)
        End Do

        Close(tmp_unit)
      End If

*  y_grid:
      Call Uniform_Y(Ny, y_leng, gycor)

*  Compute metric terms for y-grid
      sy(1) = 1d0 / y_leng
      sy(2) = sy(1) * sy(1)

*  Assign ycor values to corresponding processors
      Do J = 1, My
        ycor(J) = gycor(yrank*My+J)
      End Do

*  Compute grid spacings for y-grid
      DY = ycor(2) - ycor(1)
      
*  Output y-grid information
      If (rank .EQ. 0) Then
        Open(unit = tmp_unit,
     .       file = prefix(1:Leng_Char(prefix))//'.grid.y')

        Do J = 1, Ny
          Write(tmp_unit,10) J, gycor(J), sy(1), sy(2)
        End Do

        Close(tmp_unit)
      End If

*  z-grid:
      Call Uniform_Z(Nz, z_leng, gzcor)

*  Compute metric terms for z-grid
      sz(1) = 1d0 / z_leng
      sz(2) = sz(1) * sz(1)

*  Assign zcor values to corresponding processors
      Do K = 1, Mz
        zcor(K) = gzcor(zrank*Mz+K)
      End Do

*  Compute grid spacing for z-grid
      DZ = zcor(2) - zcor(1)

*  Output z-grid information
      If (rank .EQ. 0) Then
        Open(unit = tmp_unit,
     .       file = prefix(1:Leng_Char(prefix))//'.grid.z')

        Do K = 1, Nz
          Write(tmp_unit,10) K, gzcor(K), sz(1), sz(2)
        End Do

        Close(tmp_unit)
      End If

 10   Format(I5, 3(1X, E20.12))

      Return

      End


************************************************************************
      Subroutine Uniform_X(Nx, x_leng, gxcor)

      implicit none
      Integer Nx
      Real*8  x_leng, gxcor(Nx)

*  Locals
      Integer I

*  Generate x-mesh
      Do I = 1, Nx
        gxcor(I) = x_leng * Dble(I - 1) / Dble(Nx - 1)
      End Do

      Return

      End


************************************************************************
      Subroutine Uniform_Y(Ny, y_leng, gycor)

      implicit none
      Integer Ny
      Real*8  y_leng, gycor(Ny)

*  Locals
      Integer J

*  Generate y-mesh
      Do J = 1, Ny
        gycor(J)   = y_leng * Dble(J - 1) / Dble(Ny - 1)
      End Do

      Return

      End


************************************************************************
      Subroutine Uniform_Z(Nz, z_leng, gzcor)

      implicit none
      Integer Nz
      Real*8  z_leng, gzcor(Nz)

*  Locals
      Integer K

*  Generate z-mesh
      Do K = 1, Nz
        gzcor(K) = z_leng * Dble(K - 1) / Dble(Nz - 1)
      End Do

      Return

      End
