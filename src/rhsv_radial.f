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

      Real*8  Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)

*  Code up terms involving x-derivatives
      Call RHS_r(Q, Qrhs) 

      Return

      End

************************************************************************
*  Subroutine: RHS_r                                                   *
*                                                                      *
*  Function:   It computes all the x-derivatives in the right hand     *
*              side vector.  To reduce the data stride in the          *
*              computation, data come into this subroutine in the form *
*              Q(My,Mz,Mx).                                            *
************************************************************************
      Subroutine RHS_r(Q, Qrhs)

      Include 'header'
      Include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
      Real*8 Qconvr(My,Mz,Mx,5)
      Real*8 Qvisc(My,Mz,Mx,5) 

      Common /RHSV/ Qconvr, Qvisc

*  Globals
      Real*8   P(My,Mz,Mx),   Vs(My,Mz,Mx),  T(My,Mz,Mx)
      Real*8  FF(My,Mz,Mx), Flux(My,Mz,Mx), dF(My,Mz,Mx)
      Real*8  rad(Mx), rad_aux(Mx)

*  Locals
      Integer I, II, J, K, L, M
      Real*8 IX
      Real*8  fac1, fac2
      Integer ierr, req(2), status(MPI_STATUS_SIZE)

      Qrhs = 0D0; Qconvr = 0D0
      Qvisc = 0D0

      Do I = 1,Mx
        II     = xrank * Mx + I
        IX     = REAL(II - Nx/2) - 0.5D0
        rad(I) = dabs(IX*DX)
      End Do

*  Q(:,:,:,1) = u_r
*  Q(:,:,:,2) = 0
*  Q(:,:,:,3) = 0
*  Q(:,:,:,4) = rho
*  Q(:,:,:,5) = P/rho^gamma

*  d rho         d rho
*  ----- = - u_r -----
*   d t           d r

      Call DFDX(Q(1,1,1,4), Qrhs(1,1,1,4), 1, -1D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,4) = Qrhs(J,K,I,4) * Q(J,K,I,1)
          End Do
        End Do
      End Do

*  d rho         rho  d
*  ----- = ... - --- --- [ u_r r^2]
*   d t          r^2 d r

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            FF(J,K,I) = Q(J,K,I,1) * rad(I) * rad(I)
          End Do
        End Do
      End Do

      Call DFDX(FF, Flux, 0, 1D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,4) = Qrhs(J,K,I,4) - Flux(J,K,I)  
     .                    * (Q(J,K,I,4)/(rad(I)*rad(I)))
          End Do
        End Do
      End Do

*  d u_r         d u_r
*  ----- = - u_r ----- 
*   d t           d r

      Call DFDX(Q(1,1,1,1), Qrhs(1,1,1,1), 1, -1D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz        
            Qrhs(J,K,I,1) = Q(J,K,I,1) * Qrhs(J,K,I,1)
          End Do
        End Do
      End Do

*  d u_r          1  d P
*  ----- = ... - --- ---
*   d t          rho d r

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            P(J,K,I) = Q(J,K,I,5) * Q(J,K,I,4)**gamma 
          End Do
        End Do
      End Do

      Call DFDX(P, FF, 0, 1D0)
       
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) 
     .                    - (FF(J,K,I)/Q(J,K,I,4))  
          End Do
        End Do
      End Do

* d (P/rho^gamma)          d (P/rho^gamma)
* --------------- = - u_r  ---------------
*      d t                      d r

      Call DFDX(Q(1,1,1,5), Qrhs(1,1,1,5), 1, -1D0)

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) * Q(J,K,I,1)
          End Do
        End Do
      End Do

      Return

      End

*************************************************************************
* Radial Test                                                           *
*************************************************************************
      Subroutine RHS_test(Q, Qrhs)

      include 'header'
      include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Qrhs(My,Mz,Mx,5)
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
            Qrhs(J,K,I,1) = Qrhs(J,K,I,1) + dcos(rad(I))*dsin(rad(I))
     .                    - dsin(rad(I))
            Qrhs(J,K,I,4) = Qrhs(J,K,I,4) + dcos(rad(I))
     .                    + 2D0*dsin(rad(I))/rad(I)
            Qrhs(J,K,I,5) = Qrhs(J,K,I,5) - dsin(rad(I))*dsin(rad(I))
          End Do
        End Do
      End Do


      Do M = 0,Px-1
        if(xrank.EQ.M)then
          OPEN(UNIT=16,FILE="x_trace",FORM="FORMATTED",
     .         STATUS="UNKNOWN", POSITION="APPEND")
          Do I = 1,Mx
            WRITE(16,7) rad(I), Qrhs(1,1,I,4),
     .                  Qrhs(1,1,I,1), Qrhs(1,1,I,5)
         End Do
          CLOSE(16)
        endif
        Call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      End Do

      Stop

 7    Format(1x, 4(ES13.5))

      Return

      End
