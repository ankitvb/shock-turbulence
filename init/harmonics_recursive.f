********************************************************
* Subroutine : Spherical_Harmonics                     *
*                                                      *
* Function : Returns normalized spherical harmonics    *
*            function of specified degree and order    *
*            at a polar and azimuthal location on the  *
*            sphere                                    *
********************************************************      
      Subroutine Spherical_Harmonics(Y_lm, cost, phi, LMAX)

      Include 'header'

      Integer LMAX
      Real*8 Y_lm(LMAX,LMAX,2), P_lm(LMAX,LMAX,2)
      Integer l, m, n
      Real*8 cost, phi

      Integer lmm_fac, lpm_fac
      Real*8 fac

      P_lm = 0D0

      Call Legendre_Function(P_lm, cost, LMAX)

      Do l = 1,LMAX
        Do m = 1,l
          Call Factorial(l-m,   lmm_fac)
          Call Factorial(l+m-2, lpm_fac)
          fac = dsqrt( (2D0*DBLE(l)+1D0)*DBLE(lmm_fac) /
     .               (4D0*Pi)*DBLE(lpm_fac) )
          Y_lm(l,m,1) = fac * P_lm(l,m,1) * dcos(DBLE(m-1)*phi)
          fac = dsqrt( (2D0*DBLE(l)+1D0)*DBLE(lpm_fac) /
     .               (4D0*Pi)*DBLE(lmm_fac) )
          Y_lm(l,m,2) = fac * P_lm(l,m,2) * dcos(DBLE(m-1)*phi)
        End Do
      End Do


      Return

      End Subroutine Spherical_Harmonics

********************************************************
* Subroutine : Legendre_Function                       *
*                                                      *
* Function : Returns Legendre function of specified    *
*            degree and order for a given x            *
*            (cos(theta))                              *
*                            m                         *
*            P(l+1,m+1,1) = P                          *
*                            l                         *
*                            -m                        *
*            P(l+1,m+1,2) = P                          *
*                            l                         *
********************************************************
      Recursive Subroutine Legendre_Function(P, x, LMAX)

      Real*8  x
      Integer LMAX
      Real*8 P(LMAX,LMAX,2)
      Integer l, m, lmm_fac, lpm_fac
      Real*8 fac

      P(1,1,1) = 1D0             ! P_0,0

      P(2,1,1) = x               ! P_1,0
      P(2,2,1) = -dsqrt(1-x*x)   ! P_1,1
      P(2,2,2) = -0.5D0*P(2,2,1) ! P_1,-1
     
      Do l = 3,LMAX
        P(l,1,1) = ( (2D0*DBLE(l-2)+1D0)*x*P(l-1,1,1)-
     .             DBLE(l-2)*P(l-2,1,1) )/DBLE(l-1)
        Do m = 2,l
          if((x-1).LE.1d-3)then
            P(l,m,1) = 0D0 
          else
            P(l,m,1) = ( DBLE(l-m+1)*x*P(l,m-1,1) - 
     .                 DBLE(l+m-3)*P(l-1,m-1,1) )/dsqrt(1D0-x*x)
          endif

          if(MOD(m-1,2).EQ.0) then 
            fac = 1D0 
          else
            fac = -1D0
          endif
          Call Factorial(l-m,   lmm_fac)
          Call Factorial(l+m-2, lpm_fac)
          P(l,m,2) = fac * DBLE(lmm_fac)/DBLE(lpm_fac)*P(l,m,1)
        End Do
      End Do  
 
      Return

      End Subroutine Legendre_Function

********************************************************
* Subroutine : Factorial                               *
*                                                      *
* Function : Returns factorial of an integer n         *
********************************************************
      Recursive Subroutine Factorial(n, n_fac)

      Integer n, n_fac

      if(n .EQ. 0)then
        n_fac = 1
      else
        Call Factorial(n-1, n_fac)
        n_fac = n * n_fac
      endif

      Return

      End Subroutine Factorial
