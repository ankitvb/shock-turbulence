************************************************************************
*  Subroutine: Init_Filters                                            *
*                                                                      *
*  Function:   It is responsible for defining the left hand matrix and *
*              performing LDU decomposition on them.                   *
*                                                                      *
*  Written by Daniel J. Bodony (bodony@stanford.edu)                   *
************************************************************************
      Subroutine Init_Filters_asym
      Include 'header'
      Include 'LES_header'

*  Locals
      Integer I, J, K, N
      Real*8  Work_x(Nx,Nx), Work_y(Ny,Ny), Work_z(Nz,Nz)

*  x- filter is non-periodic eighth order Pade scheme 
*  with reduced order of accuracy near the boundaries
*  y- and z- filters are periodic 
*
*  The nodes are sixth order accurate
*
*  beta * F(I-2) + alpha * F(I-1) + F(I) + 
*
*  beta * F(I+2) + alpha * F(I+1) 
*
*          [ f(I+3) + f(I-3) ]        [ f(I+2) + f(I-2) ]
*    = a3 --------------------- + a2 ---------------------
*                   2                          2 
*
*            [ f(I+1) + f(I-1) ]
*      + a1 --------------------- + a0 f(I)
*                     2 
*
*  where
*
*  alpha =   0.663204892041650 E-00
*   beta =   0.164867523042040 E-00
*    a0  =   0.25 * (2 + 3 * alpha)
*    a1  =   (1/16) * (9 + 16 * alpha + 10 * beta)
*    a2  =   0.25 * (alpha + 4 * beta)
*    a3  =   (1/16) * (6 * beta - 1)


*  x-filter
      Do N = 1,LotX
       Do I = 1, Mx
        ldu_filt_x_asym_A(N,I) =  beta_filt4
        ldu_filt_x_asym_B(N,I) = alpha_filt4
        ldu_filt_x_asym_C(N,I) =       1.0D0
        ldu_filt_x_asym_D(N,I) = alpha_filt4
        ldu_filt_x_asym_E(N,I) =  beta_filt4
       End Do
      End Do

*  y-filter
      Do N = 1,LotY
       Do J = 1, My
        ldu_filt_y_asym_A(N,J) =  beta_filt4
        ldu_filt_y_asym_B(N,J) = alpha_filt4
        ldu_filt_y_asym_C(N,J) =       1.0D0
        ldu_filt_y_asym_D(N,J) = alpha_filt4
        ldu_filt_y_asym_E(N,J) =  beta_filt4
       End Do
      End Do

*  z-filter
      Do N = 1,LotZ
       Do K = 1, Mz
        ldu_filt_z_asym_A(N,K) =  beta_filt4
        ldu_filt_z_asym_B(N,K) = alpha_filt4
        ldu_filt_z_asym_C(N,K) =       1.0D0
        ldu_filt_z_asym_D(N,K) = alpha_filt4
        ldu_filt_z_asym_E(N,K) =  beta_filt4
       End Do 
      End Do

*  The first and last nodes are filtered explicitly
*                               
*    f'(1)   =  f(1) 
*                               
*    f'(N)   =  f(N) 
*                   
      if(xrank .EQ. 0)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,1) = 0.0D0
        ldu_filt_x_asym_B(N,1) = 0.0D0
        ldu_filt_x_asym_C(N,1) = 1.0D0
        ldu_filt_x_asym_D(N,1) = 0.0D0
        ldu_filt_x_asym_E(N,1) = 0.0D0
       End Do
      endif


      if(xrank .EQ. Px-1)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,Mx) = 0.0D0
        ldu_filt_x_asym_B(N,Mx) = 0.0D0
        ldu_filt_x_asym_C(N,Mx) = 1.0D0
        ldu_filt_x_asym_D(N,Mx) = 0.0D0
        ldu_filt_x_asym_E(N,Mx) = 0.0D0
   
       End Do
      endif


      if(yrank .EQ. 0)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,1) = 0.0D0
        ldu_filt_y_asym_B(N,1) = 0.0D0
        ldu_filt_y_asym_C(N,1) = 1.0D0
        ldu_filt_y_asym_D(N,1) = 0.0D0
        ldu_filt_y_asym_E(N,1) = 0.0D0
       End Do
      endif

      if(yrank .EQ. Py-1)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,My) = 0.0D0
        ldu_filt_y_asym_B(N,My) = 0.0D0
        ldu_filt_y_asym_C(N,My) = 1.0D0
        ldu_filt_y_asym_D(N,My) = 0.0D0
        ldu_filt_y_asym_E(N,My) = 0.0D0
       End Do
      endif

      if(zrank .EQ. 0)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,1) = 0.0D0
        ldu_filt_z_asym_B(N,1) = 0.0D0
        ldu_filt_z_asym_C(N,1) = 1.0D0
        ldu_filt_z_asym_D(N,1) = 0.0D0
        ldu_filt_z_asym_E(N,1) = 0.0D0
       End Do
      endif

      if(zrank .EQ. Pz-1)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,Mz) = 0.0D0
        ldu_filt_z_asym_B(N,Mz) = 0.0D0
        ldu_filt_z_asym_C(N,Mz) = 1.0D0
        ldu_filt_z_asym_D(N,Mz) = 0.0D0
        ldu_filt_z_asym_E(N,Mz) = 0.0D0
       End Do
      endif

*  The second and second to last nodes use use one-sided stencils
*  and are fourth order accurate
*
*  bndy_alpha_filt F(1) + F(2) + bndy_alpha_filt F(3)
*
*      b                 b
*    = - f(1) + a f(2) + - f(3)
*      2                 2
*
*
*  bndy_alpha_filt F(N-2) + F(N-1) + bndy_alpha_filt F(N)
*
*      b                     b
*    = - f(N-2) + a f(N-1) + - f(N)
*      2                     2
*
* where a = 0.5D0 + bndy_alpha_filt
*       b = a

      if(xrank .EQ. 0)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,2)     = 0.0D0
        ldu_filt_x_asym_B(N,2)     = alpha_filt4
        ldu_filt_x_asym_C(N,2)     = 1.0D0 - beta_filt4
        ldu_filt_x_asym_D(N,2)     = alpha_filt4
        ldu_filt_x_asym_E(N,2)     = beta_filt4
       End Do
      endif

      if(xrank .EQ. Px-1)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,Mx-1)     = 0.0D0
        ldu_filt_x_asym_B(N,Mx-1)     = bndy_alpha_filt
        ldu_filt_x_asym_C(N,Mx-1)     = 1.0D0
        ldu_filt_x_asym_D(N,Mx-1)     = bndy_alpha_filt
        ldu_filt_x_asym_E(N,Mx-1)     = 0.0D0
       End Do
      endif

      if(yrank .EQ. 0)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,2)     = 0.0D0
        ldu_filt_y_asym_B(N,2)     = alpha_filt4
        ldu_filt_y_asym_C(N,2)     = 1.0D0 - beta_filt4
        ldu_filt_y_asym_D(N,2)     = alpha_filt4
        ldu_filt_y_asym_E(N,2)     = beta_filt4
       End Do
      endif

      if(yrank .EQ. Py-1)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,My-1)     = 0.0D0
        ldu_filt_y_asym_B(N,My-1)     = bndy_alpha_filt
        ldu_filt_y_asym_C(N,My-1)     = 1.0D0
        ldu_filt_y_asym_D(N,My-1)     = bndy_alpha_filt
        ldu_filt_y_asym_E(N,My-1)     = 0.0D0
       End Do
      endif

      if(zrank .EQ. 0)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,2)     = 0.0D0
        ldu_filt_z_asym_B(N,2)     = alpha_filt4
        ldu_filt_z_asym_C(N,2)     = 1.0D0 - beta_filt4
        ldu_filt_z_asym_D(N,2)     = alpha_filt4
        ldu_filt_z_asym_E(N,2)     = beta_filt4
       End Do
      endif

      if(zrank .EQ. Pz-1)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,Mz-1)     = 0.0D0
        ldu_filt_z_asym_B(N,Mz-1)     = bndy_alpha_filt
        ldu_filt_z_asym_C(N,Mz-1)     = 1.0D0
        ldu_filt_z_asym_D(N,Mz-1)     = bndy_alpha_filt
        ldu_filt_z_asym_E(N,Mz-1)     = 0.0D0
       End Do
      endif

*  The third and third to last nodes are fourth order accurate
*
*  bndy_alpha_filt * F(2) + F(3) + bndy_alpha_filt * F(4)
*
*      = (c / 2) * f(1) + (b / 2) * f(2) + a * f(3)
*        (c / 2) * f(5) + (b / 2) * f(4)
*
*  bndy_alpha_filt * F(N-3) + F(N-2) + bndy_alpha_filt * F(N-1)
*
*      = (c / 2) * f(N-4) + (b / 2) * f(N-3) + a * f(N-2)
*        (c / 2) * f(N-1) + (b / 2) * f(N)
*
* where
*
*    a = (1 / 8) * ( 5 + 6 * bndy_alpha_filt)
*    b = (1 / 2) * ( 1 + 2 * bndy_alpha_filt)
*    c = (1 / 8) * (-1 + 2 * bndy_alpha_filt)

      if(xrank .EQ. 0)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,3)     = beta_filt4
        ldu_filt_x_asym_B(N,3)     = alpha_filt4
        ldu_filt_x_asym_C(N,3)     = 1.0D0
        ldu_filt_x_asym_D(N,3)     = alpha_filt4
        ldu_filt_x_asym_E(N,3)     = beta_filt4
       End Do
      endif

      if(xrank .EQ. Px-1)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,Mx-2)     = 0.0D0
        ldu_filt_x_asym_B(N,Mx-2)     = bndy_alpha_filt
        ldu_filt_x_asym_C(N,Mx-2)     = 1.0D0
        ldu_filt_x_asym_D(N,Mx-2)     = bndy_alpha_filt
        ldu_filt_x_asym_E(N,Mx-2)     = 0.0D0
       End Do
      endif

      if(yrank .EQ. 0)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,3)     = beta_filt4
        ldu_filt_y_asym_B(N,3)     = alpha_filt4
        ldu_filt_y_asym_C(N,3)     = 1.0D0
        ldu_filt_y_asym_D(N,3)     = alpha_filt4
        ldu_filt_y_asym_E(N,3)     = beta_filt4
       End Do
      endif

      if(yrank .EQ. Py-1)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,My-2)     = 0.0D0
        ldu_filt_y_asym_B(N,My-2)     = bndy_alpha_filt
        ldu_filt_y_asym_C(N,My-2)     = 1.0D0
        ldu_filt_y_asym_D(N,My-2)     = bndy_alpha_filt
        ldu_filt_y_asym_E(N,My-2)     = 0.0D0
       End Do
      endif

      if(zrank .EQ. 0)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,3)     = beta_filt4
        ldu_filt_z_asym_B(N,3)     = alpha_filt4
        ldu_filt_z_asym_C(N,3)     = 1.0D0
        ldu_filt_z_asym_D(N,3)     = alpha_filt4
        ldu_filt_z_asym_E(N,3)     = beta_filt4
       End Do
      endif

      if(zrank .EQ. Pz-1)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,Mz-2)     = 0.0D0
        ldu_filt_z_asym_B(N,Mz-2)     = bndy_alpha_filt
        ldu_filt_z_asym_C(N,Mz-2)     = 1.0D0
        ldu_filt_z_asym_D(N,Mz-2)     = bndy_alpha_filt
        ldu_filt_z_asym_E(N,Mz-2)     = 0.0D0
       End Do
      endif

*  The fourth and fourth to last nodes are sixth order accurate
*
*  bndy_alpha_filt * F(3) + F(4) + bndy_alpha_filt * F(5)
*
*      = (d / 2) * f(1) + (c / 2) * f(2) + (b / 2) * f(3) + a * f(4)
*        (d / 2) * f(7) + (c / 2) * f(6) + (b / 2) * f(5)
*
*  bndy_alpha_filt * F(N-3) + F(N-2) + bndy_alpha_filt * F(N-1)
*
*      = (d / 2) * f(N-6) + (c / 2) * f(N-5) + (b / 2) * f(N-4) + a * f(N-3)
*        (d / 2) * f(N)   + (c / 2) * f(N-1) + (b / 2) * f(N-2)
*
* where
*
*    a = (1 / 16) * (11 + 10 * alpha_filt4 - 10 * beta_filt4)
*    b = (1 / 32) * (15 + 34 * alpha_filt4 + 30 * beta_filt4)
*    c = (1 / 16) * (-3 + 6  * alpha_filt4 + 26 * beta_filt4)
*    d = (1 / 32) * (1  - 2  * alpha_filt4 + 2  * beta_filt4)

      if(xrank .EQ. 0)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,4)  =  beta_filt4 
        ldu_filt_x_asym_B(N,4)  = alpha_filt4
        ldu_filt_x_asym_C(N,4)  =       1.0D0
        ldu_filt_x_asym_D(N,4)  = alpha_filt4
        ldu_filt_x_asym_E(N,4)  =  beta_filt4
       End Do
      endif

      if(xrank .EQ. Px-1)then
       Do N = 1,LotX
        ldu_filt_x_asym_A(N,Mx-3)  =  beta_filt4
        ldu_filt_x_asym_B(N,Mx-3)  = alpha_filt4
        ldu_filt_x_asym_C(N,Mx-3)  =       1.0D0
        ldu_filt_x_asym_D(N,Mx-3)  = alpha_filt4
        ldu_filt_x_asym_E(N,Mx-3)  =  beta_filt4
       End Do
      endif

      if(yrank .EQ. 0)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,4)  =  beta_filt4
        ldu_filt_y_asym_B(N,4)  = alpha_filt4
        ldu_filt_y_asym_C(N,4)  =       1.0D0
        ldu_filt_y_asym_D(N,4)  = alpha_filt4
        ldu_filt_y_asym_E(N,4)  =  beta_filt4
       End Do
      endif

      if(yrank .EQ. Py-1)then
       Do N = 1,LotY
        ldu_filt_y_asym_A(N,My-3)  =  beta_filt4
        ldu_filt_y_asym_B(N,My-3)  = alpha_filt4
        ldu_filt_y_asym_C(N,My-3)  =       1.0D0
        ldu_filt_y_asym_D(N,My-3)  = alpha_filt4
        ldu_filt_y_asym_E(N,My-3)  =  beta_filt4
       End Do
      endif

      if(zrank .EQ. 0)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,4)  =  beta_filt4
        ldu_filt_z_asym_B(N,4)  = alpha_filt4
        ldu_filt_z_asym_C(N,4)  =       1.0D0
        ldu_filt_z_asym_D(N,4)  = alpha_filt4
        ldu_filt_z_asym_E(N,4)  =  beta_filt4
       End Do
      endif

      if(zrank .EQ. Pz-1)then
       Do N = 1,LotZ
        ldu_filt_z_asym_A(N,Mz-3)  =  beta_filt4
        ldu_filt_z_asym_B(N,Mz-3)  = alpha_filt4
        ldu_filt_z_asym_C(N,Mz-3)  =       1.0D0
        ldu_filt_z_asym_D(N,Mz-3)  = alpha_filt4
        ldu_filt_z_asym_E(N,Mz-3)  =  beta_filt4
       End Do
      endif

      Return

      End

*********************************************************************
      Subroutine Compact_Filter_X_asym(F, F0)

      Include 'header'
      Include 'LES_header'
      Include 'mpif.h'

      Real*8  F(My,Mz,Mx), F0(My,Mz,Mx)
      Real*8  FP(My,Mz,4), FN(My,Mz,4)
      Real*8  W(LotX, Mx), W2(My,Mz,Mx)
      Real*8  temp(My,Mz,4,2)
      Integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)

*  Locals
      Integer I, J, K, N

* Putting first 4 and last 4 data planes into temp arrays for sending
      temp(:,:,:,1) = F(:,:,1:4)
      temp(:,:,:,2) = F(:,:,Mx-3:Mx)

* Taking care of boundary nodes
      if(xrank .EQ. 0)then
        source_right_x = MPI_PROC_NULL
        dest_left_x    = MPI_PROC_NULL
      elseif(xrank .EQ. Px-1)then
        source_left_x  = MPI_PROC_NULL
        dest_right_x   = MPI_PROC_NULL
      endif

* Transferring data from neighbouring processors
      Call MPI_IRECV(FP(1,1,1),   My*Mz*4, MPI_REAL8,source_right_x,
     .               source_right_x,                                 ! Recv last 4 planes
     .               MPI_XYZ_COMM, req(1), ierr)

      Call MPI_IRECV(FN(1,1,1),   My*Mz*4, MPI_REAL8,source_left_x,
     .               source_left_x,                                  ! Recv first 4 planes
     .               MPI_XYZ_COMM, req(2), ierr)

      Call MPI_ISEND(temp(1,1,1,2), My*Mz*4, MPI_REAL8,dest_right_x,
     .               rank,                                           ! Send last 4 planes
     .               MPI_XYZ_COMM, req(3), ierr)

      Call MPI_ISEND(temp(1,1,1,1), My*Mz*4, MPI_REAL8,dest_left_x,
     .               rank,                                           ! Send first 4 planes
     .               MPI_XYZ_COMM, req(4), ierr)


* Compute interior in the meantime
      Do J = 1, My
        Do K = 1, Mz
          N =  (J-1) * Mz + K
          Do I = 5, Mx-4
            W(N,I) =     a4_filt4 * ( F(J,K,I+4) + F(J,K,I-4) )
     .                 + a3_filt4 * ( F(J,K,I+3) + F(J,K,I-3) )
     .                 + a2_filt4 * ( F(J,K,I+2) + F(J,K,I-2) )
     .                 + a1_filt4 * ( F(J,K,I+1) + F(J,K,I-1) )
     .                 + a0_filt4 *   F(J,K,I)
          End Do
        End Do
      End Do


*  Wait to recieve planes
      Call MPI_WAITALL(4, req, status_arr, ierr)

      Do J = 1, My
        Do K = 1, Mz

          N =  (J-1) * Mz + K

          if(xrank .EQ. 0)then
           W(N,1) = 0D0
 
           W(N,2) = a4_filt4 * ( F(J,K,6) - F(J,K,4) )
     .            + a3_filt4 * ( F(J,K,5) - F(J,K,3) )
     .            + a2_filt4 * ( F(J,K,4) - F(J,K,2) )
     .            + a1_filt4 * ( F(J,K,1) + F(J,K,3) )
     .            + a0_filt4 *   F(J,K,2)

           W(N,3) = a4_filt4 * ( F(J,K,7) - F(J,K,3) )
     .            + a3_filt4 * ( F(J,K,6) - F(J,K,2) )
     .            + a2_filt4 * ( F(J,K,1) + F(J,K,5) )
     .            + a1_filt4 * ( F(J,K,2) + F(J,K,4) )
     .            + a0_filt4 *   F(J,K,3)

           W(N,4) = a4_filt4 * ( F(J,K,8) - F(J,K,2) )
     .            + a3_filt4 * ( F(J,K,1) + F(J,K,7) )
     .            + a2_filt4 * ( F(J,K,2) + F(J,K,6) )
     .            + a1_filt4 * ( F(J,K,3) + F(J,K,5) )
     .            + a0_filt4 *   F(J,K,4)
          else
           W(N,1) =  a4_filt4 * ( F(J,K,5) + FP(J,K,1) )
     .             + a3_filt4 * ( F(J,K,4) + FP(J,K,2) )
     .             + a2_filt4 * ( F(J,K,3) + FP(J,K,3) )
     .             + a1_filt4 * ( F(J,K,2) + FP(J,K,4) )
     .             + a0_filt4 *   F(J,K,1)


           W(N,2) =  a4_filt4 * ( F(J,K,6) + FP(J,K,2) )
     .             + a3_filt4 * ( F(J,K,5) + FP(J,K,3) )
     .             + a2_filt4 * ( F(J,K,4) + FP(J,K,4) )
     .             + a1_filt4 * ( F(J,K,3) + F(J,K,1) )
     .             + a0_filt4 *   F(J,K,2)

           W(N,3) =  a4_filt4 * ( F(J,K,7) + FP(J,K,3) )
     .             + a3_filt4 * ( F(J,K,6) + FP(J,K,4) )
     .             + a2_filt4 * ( F(J,K,5) + F(J,K,1)  )
     .             + a1_filt4 * ( F(J,K,4) + F(J,K,2)  )
     .             + a0_filt4 *   F(J,K,3)

           W(N,4) =  a4_filt4 * ( F(J,K,8) + FP(J,K,4) )
     .             + a3_filt4 * ( F(J,K,7) + F(J,K,1)  )
     .             + a2_filt4 * ( F(J,K,6) + F(J,K,2)  )
     .             + a1_filt4 * ( F(J,K,5) + F(J,K,3)  )
     .             + a0_filt4 *   F(J,K,4)
          endif
        End Do
      End Do

      Do J = 1, My
        Do K = 1, Mz

          N =  (J-1) * Mz + K

          if(xrank .EQ. Px-1)then
           W(N,Mx)   = bndy_a1_closure1_filt * F(J,K,Mx)
     .               + bndy_a2_closure1_filt * F(J,K,Mx-1)
     .               + bndy_a3_closure1_filt * F(J,K,Mx-2)
     .               + bndy_a4_closure1_filt * F(J,K,Mx-3)
     .               + bndy_a5_closure1_filt * F(J,K,Mx-4)

           W(N,Mx-1) = bndy_a1_closure2_filt * F(J,K,Mx)
     .               + bndy_a2_closure2_filt * F(J,K,Mx-1)
     .               + bndy_a3_closure2_filt * F(J,K,Mx-2)
     .               + bndy_a4_closure2_filt * F(J,K,Mx-3)
     .               + bndy_a5_closure2_filt * F(J,K,Mx-4)

           W(N,Mx-2) = bndy_a2_closure3_filt * F(J,K,Mx-4)
     .               + bndy_a1_closure3_filt * F(J,K,Mx-3)
     .               + bndy_a0_closure3_filt * F(J,K,Mx-2)
     .               + bndy_a1_closure3_filt * F(J,K,Mx-1)
     .               + bndy_a2_closure3_filt * F(J,K,Mx)

           W(N,Mx-3) = bndy_a3_closure4_filt * F(J,K,Mx-6)
     .               + bndy_a2_closure4_filt * F(J,K,Mx-5)
     .               + bndy_a1_closure4_filt * F(J,K,Mx-4)
     .               + bndy_a0_closure4_filt * F(J,K,Mx-3)
     .               + bndy_a1_closure4_filt * F(J,K,Mx-2)
     .               + bndy_a2_closure4_filt * F(J,K,Mx-1)
     .               + bndy_a3_closure4_filt * F(J,K,Mx)
          else
           W(N,Mx-3) =   a4_filt4 * ( FN(J,K,1)   + F(J,K,Mx-7) )
     .                 + a3_filt4 * ( F(J,K,Mx)   + F(J,K,Mx-6) )
     .                 + a2_filt4 * ( F(J,K,Mx-1) + F(J,K,Mx-5) )
     .                 + a1_filt4 * ( F(J,K,Mx-2) + F(J,K,Mx-4) )
     .                 + a0_filt4 *   F(J,K,Mx-3)

           W(N,Mx-2) =   a4_filt4 * ( FN(J,K,2) + F(J,K,Mx-6)   )
     .                 + a3_filt4 * ( FN(J,K,1) + F(J,K,Mx-5)   )
     .                 + a2_filt4 * ( F(J,K,Mx) + F(J,K,Mx-4)   )
     .                 + a1_filt4 * ( F(J,K,Mx-1) + F(J,K,Mx-3) )
     .                 + a0_filt4 *   F(J,K,Mx-2)

           W(N,Mx-1) =   a4_filt4 * ( FN(J,K,3) + F(J,K,Mx-5) )
     .                 + a3_filt4 * ( FN(J,K,2) + F(J,K,Mx-4) )
     .                 + a2_filt4 * ( FN(J,K,1) + F(J,K,Mx-3) )
     .                 + a1_filt4 * ( F(J,K,Mx) + F(J,K,Mx-2) )
     .                 + a0_filt4 *   F(J,K,Mx-1)

           W(N,Mx)   =   a4_filt4 * ( FN(J,K,4) + F(J,K,Mx-4) )
     .                 + a3_filt4 * ( FN(J,K,3) + F(J,K,Mx-3) )
     .                 + a2_filt4 * ( FN(J,K,2) + F(J,K,Mx-2) )
     .                 + a1_filt4 * ( FN(J,K,1) + F(J,K,Mx-1) )
     .                 + a0_filt4 *   F(J,K,Mx)
          endif
        End Do
      End Do

*  Solve the matrix system to obtain the derivative dF
      Call pentadiagonal(ldu_filt_x_asym_A, ldu_filt_x_asym_B,
     .                   ldu_filt_x_asym_C, ldu_filt_x_asym_D,
     .                   ldu_filt_x_asym_E, W,
     .                   Mx, LotX, 0, MPI_X_COMM)

* Converting 2-D array to 3-D for copying
      Do N = 0,LotX-1
       Do I = 1,Mx
         J = N / Mz + 1
         K = MOD(N,Mz) + 1
         W2(J,K,I) = W(N+1,I)
         F0(J,K,I) = W2(J,K,I)
       End Do
      End Do

      Return

      End

*******************************************************************
      Subroutine Compact_Filter_Y_asym(F, F0)

      Include 'header'
      Include 'LES_header'
      Include 'mpif.h'

      Real*8  F(My,Mz,Mx), F0(My,Mz,Mx)
      Real*8  FP(4,Mz,Mx), FN(4,Mz,Mx)
      Real*8  W(LotY, My), W2(My,Mz,Mx)
      Real*8  temp(4,Mz,Mx,2)
      Integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)

*  Locals
      Integer I, J, K, N

* Putting first 4 and last 4 data planes into temp arrays for sending
      temp(:,:,:,1) = F(1:4,:,:)
      temp(:,:,:,2) = F(My-3:My,:,:)

* Taking care of boundary nodes
      if(yrank .EQ. 0)then
        source_right_y = MPI_PROC_NULL
        dest_left_y    = MPI_PROC_NULL
      elseif(yrank .EQ. Py-1)then
        source_left_y  = MPI_PROC_NULL
        dest_right_y   = MPI_PROC_NULL
      endif

* Transferring data from neighbouring processors
      Call MPI_IRECV(FP(1,1,1),   Mz*Mx*4, MPI_REAL8,source_right_y,
     .               source_right_y,                                 ! Recv last 4 planes
     .               MPI_XYZ_COMM, req(1), ierr)

      Call MPI_IRECV(FN(1,1,1),   Mz*Mx*4, MPI_REAL8,source_left_y,
     .               source_left_y,                                 ! Recv last 4 planes
     .               MPI_XYZ_COMM, req(2), ierr)

      Call MPI_ISEND(temp(1,1,1,2), Mz*Mx*4, MPI_REAL8,dest_right_y,
     .               rank,                                           ! Send last 4 planes
     .               MPI_XYZ_COMM, req(3), ierr)

      Call MPI_ISEND(temp(1,1,1,1), Mz*Mx*4, MPI_REAL8,dest_left_y,
     .               rank,                                           ! Send first 4 planes
     .               MPI_XYZ_COMM, req(4), ierr)


* Compute interior in the meantime
      Do I = 1, Mx
        Do K = 1, Mz
          N =  (I-1) * Mz + K
          Do J = 5, My-4
            W(N,J) =     a4_filt4 * ( F(J+4,K,I) + F(J-4,K,I) )
     .                 + a3_filt4 * ( F(J+3,K,I) + F(J-3,K,I) )
     .                 + a2_filt4 * ( F(J+2,K,I) + F(J-2,K,I) )
     .                 + a1_filt4 * ( F(J+1,K,I) + F(J-1,K,I) )
     .                 + a0_filt4 * F(J,K,I)
          End Do
        End Do
      End Do

*  Wait to recieve planes
      Call MPI_WAITALL(4, req, status_arr, ierr)

      Do I = 1, Mx
        Do K = 1, Mz

          N =  (I-1) * Mz + K

          if(yrank .EQ. 0)then
           W(N,1) = 0D0

           W(N,2) = a4_filt4 * ( F(6,K,I) - F(4,K,I) )
     .            + a3_filt4 * ( F(5,K,I) - F(3,K,I) )
     .            + a2_filt4 * ( F(4,K,I) - F(2,K,I) )
     .            + a1_filt4 * ( F(1,K,I) + F(3,K,I) )
     .            + a0_filt4 *   F(2,K,I)

           W(N,3) = a4_filt4 * ( F(7,K,I) - F(3,K,I) )
     .            + a3_filt4 * ( F(6,K,I) - F(2,K,I) )
     .            + a2_filt4 * ( F(1,K,I) + F(5,K,I) )
     .            + a1_filt4 * ( F(2,K,I) + F(4,K,I) )
     .            + a0_filt4 *   F(3,K,I)

           W(N,4) = a4_filt4 * ( F(8,K,I) - F(2,K,I) )
     .            + a3_filt4 * ( F(1,K,I) + F(7,K,I) )
     .            + a2_filt4 * ( F(2,K,I) + F(6,K,I) )
     .            + a1_filt4 * ( F(3,K,I) + F(5,K,I) )
     .            + a0_filt4 *   F(4,K,I)
          else
           W(N,1) =  a4_filt4 * ( F(5,K,I) + FP(1,K,I) )
     .             + a3_filt4 * ( F(4,K,I) + FP(2,K,I) )
     .             + a2_filt4 * ( F(3,K,I) + FP(3,K,I) )
     .             + a1_filt4 * ( F(2,K,I) + FP(4,K,I) )
     .             + a0_filt4 *   F(1,K,I)


           W(N,2) =  a4_filt4 * ( F(6,K,I) + FP(2,K,I) )
     .             + a3_filt4 * ( F(5,K,I) + FP(3,K,I) )
     .             + a2_filt4 * ( F(4,K,I) + FP(4,K,I) )
     .             + a1_filt4 * ( F(3,K,I) + F(1,K,I) )
     .             + a0_filt4 *   F(2,K,I)

           W(N,3) =  a4_filt4 * ( F(7,K,I) + FP(3,K,I) )
     .             + a3_filt4 * ( F(6,K,I) + FP(4,K,I) )
     .             + a2_filt4 * ( F(5,K,I) + F(1,K,I) )
     .             + a1_filt4 * ( F(4,K,I) + F(2,K,I) )
     .             + a0_filt4   * F(3,K,I)

           W(N,4) =  a4_filt4 * ( F(8,K,I) + FP(4,K,I) )
     .             + a3_filt4 * ( F(7,K,I) + F(1,K,I) )
     .             + a2_filt4 * ( F(6,K,I) + F(2,K,I) )
     .             + a1_filt4 * ( F(5,K,I) + F(3,K,I) )
     .             + a0_filt4   * F(4,K,I)
          endif
        End Do
      End Do
       
      Do I = 1, Mx
        Do K = 1, Mz

          N =  (I-1) * Mz + K

         if(yrank .EQ. Py-1)then
          W(N,My)   = bndy_a1_closure1_filt * F(My,K,I)
     .              + bndy_a2_closure1_filt * F(My-1,K,I)
     .              + bndy_a3_closure1_filt * F(My-2,K,I)
     .              + bndy_a4_closure1_filt * F(My-3,K,I)
     .              + bndy_a5_closure1_filt * F(My-4,K,I)

          W(N,My-1) = bndy_a1_closure2_filt * F(My,K,I)
     .              + bndy_a2_closure2_filt * F(My-1,K,I)
     .              + bndy_a3_closure2_filt * F(My-2,K,I)
     .              + bndy_a4_closure2_filt * F(My-3,K,I)
     .              + bndy_a5_closure2_filt * F(My-4,K,I)

          W(N,My-2) = bndy_a2_closure3_filt * F(My-4,K,I)
     .              + bndy_a1_closure3_filt * F(My-3,K,I)
     .              + bndy_a0_closure3_filt * F(My-2,K,I)
     .              + bndy_a1_closure3_filt * F(My-1,K,I)
     .              + bndy_a2_closure3_filt * F(My  ,K,I)

          W(N,My-3) = bndy_a3_closure4_filt * F(My-6,K,I)
     .              + bndy_a2_closure4_filt * F(My-5,K,I)
     .              + bndy_a1_closure4_filt * F(My-4,K,I)
     .              + bndy_a0_closure4_filt * F(My-3,K,I)
     .              + bndy_a1_closure4_filt * F(My-2,K,I)
     .              + bndy_a2_closure4_filt * F(My-1,K,I)
     .              + bndy_a3_closure4_filt * F(My  ,K,I)
         else
          W(N,My-3) = a4_filt4 * ( FN(1,K,I) + F(My-7,K,I) )
     .              + a3_filt4 * ( F(My,K,I) + F(My-6,K,I) )
     .              + a2_filt4 * ( F(My-1,K,I) + F(My-5,K,I) )
     .              + a1_filt4 * ( F(My-2,K,I) + F(My-4,K,I) )
     .              + a0_filt4 *   F(My-3,K,I)

          W(N,My-2) = a4_filt4 * ( FN(2,K,I) + F(My-6,K,I) )
     .              + a3_filt4 * ( FN(1,K,I) + F(My-5,K,I) )
     .              + a2_filt4 * ( F(My,K,I) + F(My-4,K,I) )
     .              + a1_filt4 * ( F(My-1,K,I) + F(My-3,K,I) )
     .              + a0_filt4 *   F(My-2,K,I)

          W(N,My-1) = a4_filt4 * ( FN(3, K,I) + F(My-5,K,I) )
     .              + a3_filt4 * ( FN(2, K,I) + F(My-4,K,I) )
     .              + a2_filt4 * ( FN(1, K,I) + F(My-3,K,I) )
     .              + a1_filt4 * ( F(My,K,I) + F(My-2,K,I) )
     .              + a0_filt4 * F(My-1,K,I)

          W(N,  My) = a4_filt4 * ( FN(4,K,I) + F(My-4,K,I) )
     .              + a3_filt4 * ( FN(3,K,I) + F(My-3,K,I) )
     .              + a2_filt4 * ( FN(2,K,I) + F(My-2,K,I) )
     .              + a1_filt4 * ( FN(1,K,I) + F(My-1,K,I) )
     .              + a0_filt4 *  F(My,K,I)
         endif
        End Do
      End Do

*  Solve the matrix system to obtain the derivative dF
      Call pentadiagonal(ldu_filt_y_asym_A, ldu_filt_y_asym_B,
     .                   ldu_filt_y_asym_C, ldu_filt_y_asym_D,
     .                   ldu_filt_y_asym_E, W,
     .                   My, LotY, 0, MPI_Y_COMM)

* Converting 2-D array to 3-D for copying
      Do N = 0,LotY-1
       Do J = 1,My
         I = N / Mz + 1
         K = MOD(N,Mz) + 1
         W2(J,K,I) = W(N+1,J)
         F0(J,K,I) = W2(J,K,I)
       End Do
      End Do

      Return

      End

******************************************************************
      Subroutine Compact_Filter_Z_asym(F, F0)

      Include 'header'
      Include 'LES_header'
      Include 'mpif.h'

      Real*8  F(My,Mz,Mx), F0(My,Mz,Mx)
      Real*8  FP(My,4,Mx), FN(My,4,Mx)
      Real*8  W(LotZ, Mz), W2(My,Mz,Mx)
      Real*8  temp(My,4,Mx,2)
      Integer status_arr(MPI_STATUS_SIZE,4), ierr, req(4)

*  Locals
      Integer I, J, K, N

* Putting first 4 and last 4 data planes into temp arrays for sending
      temp(:,:,:,1) = F(:,1:4,:)
      temp(:,:,:,2) = F(:,Mz-3:Mz,:)

* Taking care of boundary nodes
      if(zrank .EQ. 0)then
        source_right_z = MPI_PROC_NULL
        dest_left_z    = MPI_PROC_NULL
      elseif(zrank .EQ. Pz-1)then
        source_left_z  = MPI_PROC_NULL
        dest_right_z   = MPI_PROC_NULL
      endif

* Transferring data from neighbouring processors
      Call MPI_IRECV(FP(1,1,1),    Mx*My*4, MPI_REAL8,source_right_z,
     .               source_right_z,                                 ! Recv last 4 planes
     .               MPI_XYZ_COMM, req(1), ierr)

      Call MPI_IRECV(FN(1,1,1),    Mx*My*4, MPI_REAL8,source_left_z,
     .               source_left_z,                                  ! Recv first 4 planes
     .               MPI_XYZ_COMM, req(2), ierr)

      Call MPI_ISEND(temp(1,1,1,2), Mx*My*4, MPI_REAL8,dest_right_z,
     .               rank,                                           ! Send last 4 planes
     .               MPI_XYZ_COMM, req(3), ierr)

      Call MPI_ISEND(temp(1,1,1,1), Mx*My*4, MPI_REAL8,dest_left_z,
     .               rank,                                           ! Send first 4 planes
     .               MPI_XYZ_COMM, req(4), ierr)

* Compute interior in the meantime
      Do I = 1, Mx
        Do J = 1, My
          N =  (I-1) * My + J
          Do K = 5, Mz-4
            W(N,K) =     a4_filt4 * ( F(J,K+4,I) + F(J,K-4,I) )
     .                 + a3_filt4 * ( F(J,K+3,I) + F(J,K-3,I) )
     .                 + a2_filt4 * ( F(J,K+2,I) + F(J,K-2,I) )
     .                 + a1_filt4 * ( F(J,K+1,I) + F(J,K-1,I) )
     .                 + a0_filt4 * F(J,K,I)
          End Do
        End Do
      End Do

*  Wait to recieve planes
      Call MPI_WAITALL(4, req, status_arr, ierr)

      Do I = 1, Mx
        Do J = 1, My

          N =  (I-1) * My + J

          if(zrank .EQ. 0)then
           W(N,1) = 0D0

           W(N,2) = a4_filt4 * ( F(J,6,I) - F(J,4,I) )
     .              + a3_filt4 * ( F(J,5,I) - F(J,3,I) )
     .              + a2_filt4 * ( F(J,4,I) - F(J,2,I) )
     .              + a1_filt4 * ( F(J,1,I) + F(J,3,I) )
     .              + a0_filt4 *   F(J,2,I)

           W(N,3) = a4_filt4 * ( F(J,7,I) - F(J,3,I) )
     .              + a3_filt4 * ( F(J,6,I) - F(J,2,I) )
     .              + a2_filt4 * ( F(J,1,I) + F(J,5,I) )
     .              + a1_filt4 * ( F(J,2,I) + F(J,4,I) )
     .              + a0_filt4 *   F(J,3,I)

           W(N,4) = a4_filt4 * ( F(J,8,I) - F(J,2,I) )
     .              + a3_filt4 * ( F(J,1,I) + F(J,7,I) )
     .              + a2_filt4 * ( F(J,2,I) + F(J,6,I) )
     .              + a1_filt4 * ( F(J,3,I) + F(J,5,I) )
     .             + a0_filt4 *   F(J,4,I)
          else
           W(N,1) = a4_filt4 * ( F(J,5,I) + FP(J,1,I) )
     .            + a3_filt4 * ( F(J,4,I) + FP(J,2,I) )
     .            + a2_filt4 * ( F(J,3,I) + FP(J,3,I) )
     .            + a1_filt4 * ( F(J,2,I) + FP(J,4,I) )
     .            + a0_filt4 *   F(J,1,I)

           W(N,2) = a4_filt4 * ( F(J,6,I) + FP(J,2,I) )
     .            + a3_filt4 * ( F(J,5,I) + FP(J,3,I) )
     .            + a2_filt4 * ( F(J,4,I) + FP(J,4,I) )
     .            + a1_filt4 * ( F(J,3,I) + F(J,1,I) )
     .            + a0_filt4 *   F(J,2,I)

           W(N,3) = a4_filt4 * ( F(J,7,I) + FP(J,3,I) )
     .            + a3_filt4 * ( F(J,6,I) + FP(J,4,I) )
     .            + a2_filt4 * ( F(J,5,I) + F(J,1,I) )
     .            + a1_filt4 * ( F(J,4,I) + F(J,2,I) )
     .            + a0_filt4   * F(J,3,I)

           W(N,4) = a4_filt4 * ( F(J,8,I) + FP(J,4,I) )
     .            + a3_filt4 * ( F(J,7,I) + F(J,1,I) )
     .            + a2_filt4 * ( F(J,6,I) + F(J,2,I) )
     .            + a1_filt4 * ( F(J,5,I) + F(J,3,I) )
     .            + a0_filt4   * F(J,4,I)
          endif
        End Do
      End Do

      Do I = 1, Mx
        Do J = 1, My

          N =  (I-1) * My + J

          if(zrank .EQ. Pz-1)then
           W(N,Mz)   = bndy_a1_closure1_filt * F(J,Mz,I)
     .               + bndy_a2_closure1_filt * F(J,Mz-1,I)
     .               + bndy_a3_closure1_filt * F(J,Mz-2,I)
     .               + bndy_a4_closure1_filt * F(J,Mz-3,I)
     .               + bndy_a5_closure1_filt * F(J,Mz-4,I)

           W(N,Mz-1) = bndy_a1_closure2_filt * F(J,Mz,I)
     .               + bndy_a2_closure2_filt * F(J,Mz-1,I)
     .               + bndy_a3_closure2_filt * F(J,Mz-2,I)
     .               + bndy_a4_closure2_filt * F(J,Mz-3,I)
     .               + bndy_a5_closure2_filt * F(J,Mz-4,I)

           W(N,Mz-2) = bndy_a2_closure3_filt * F(J,Mz-4,I)
     .               + bndy_a1_closure3_filt * F(J,Mz-3,I)
     .               + bndy_a0_closure3_filt * F(J,Mz-2,I)
     .               + bndy_a1_closure3_filt * F(J,Mz-1,I)
     .               + bndy_a2_closure3_filt * F(J,Mz  ,I)

           W(N,Mz-3) = bndy_a3_closure4_filt * F(J,Mz-6,I)
     .               + bndy_a2_closure4_filt * F(J,Mz-5,I)
     .               + bndy_a1_closure4_filt * F(J,Mz-4,I)
     .               + bndy_a0_closure4_filt * F(J,Mz-3,I)
     .               + bndy_a1_closure4_filt * F(J,Mz-2,I)
     .               + bndy_a2_closure4_filt * F(J,Mz-1,I)
     .               + bndy_a3_closure4_filt * F(J,Mz  ,I)
          else
           W(N,Mz-3) = a4_filt4 * ( FN(J,1,I) + F(J,Mz-7,I) )
     .               + a3_filt4 * ( F(J,Mz,I) + F(J,Mz-6,I) )
     .               + a2_filt4 * ( F(J,Mz-1,I) + F(J,Mz-5,I) )
     .               + a1_filt4 * ( F(J,Mz-2,I) + F(J,Mz-4,I) )
     .               + a0_filt4 *   F(J,Mz-3,I)

           W(N,Mz-2) = a4_filt4 * ( FN(J,2,I) + F(J,Mz-6,I) )
     .               + a3_filt4 * ( FN(J,1,I) + F(J,Mz-5,I) )
     .               + a2_filt4 * ( F(J,Mz,I) + F(J,Mz-4,I) )
     .               + a1_filt4 * ( F(J,Mz-1,I) + F(J,Mz-3,I) )
     .               + a0_filt4 *   F(J,Mz-2,I)

           W(N,Mz-1) = a4_filt4 * ( FN(J,3,I) + F(J,Mz-5,I) )
     .               + a3_filt4 * ( FN(J,2,I) + F(J,Mz-4,I) )
     .               + a2_filt4 * ( FN(J,1,I) + F(J,Mz-3,I) )
     .               + a1_filt4 * ( F(J,Mz,I) + F(J,Mz-2,I) )
     .               + a0_filt4 *   F(J,Mz-1,I)

           W(N,Mz) = a4_filt4 * ( FN(J,4,I) + F(J,Mz-4,I) )
     .             + a3_filt4 * ( FN(J,3,I) + F(J,Mz-3,I) )
     .             + a2_filt4 * ( FN(J,2,I) + F(J,Mz-2,I) )
     .             + a1_filt4 * ( FN(J,1,I) + F(J,Mz-1,I) )
     .             + a0_filt4 *   F(J,Mz,I)
          endif
        End Do
      End Do

*  Solve the matrix system to obtain the derivative dF
      Call pentadiagonal(ldu_filt_z_asym_A, ldu_filt_z_asym_B,
     .                   ldu_filt_z_asym_C, ldu_filt_z_asym_D,
     .                   ldu_filt_z_asym_E, W,
     .                   Mz, LotZ, 0, MPI_Z_COMM)

* Converting 2-D array to 3-D for copying
      Do N = 0,LotZ-1
       Do K = 1,Mz
         I = N / My + 1
         J = MOD(N,My) + 1
         W2(J,K,I) = W(N+1,K)
         F0(J,K,I) = W2(J,K,I)
       End Do
      End Do

      Return

      End
