********************************************************
* Subroutine : Spherical_Harmonics                     *
*                                                      *
* Function : Returns normalized spherical harmonics    *
*            function of specified degree and order    *
*            at a polar and azimuthal location on the  *
*            sphere                                    *
********************************************************      
      Subroutine Spherical_Harmonics(Y_lm, Y_l0, phase_lm, phase_l0,
     .                theta, phi, LMAX)

      Include 'header'

      Integer LMAX
      Real*8 Y_lm(LMAX,LMAX), Y_l0(LMAX)
      Real*8 phase_lm(LMAX,LMAX), phase_l0(LMAX)
      Integer l, m
      Real*8 cost, theta, phi, P_lm
      Real*8 m_dbl

      Do l = 1,LMAX
        cost = dcos(theta) !+phase_l0(l))
        Call getSphPlm(l, 0, cost, P_lm)
        Y_l0(l) = P_lm
        Do m = 1,l
          Call getSphPlm(l, m, cost, P_lm)                  ! Calling GSL routine for the associated Legendre function
          m_dbl = DBLE(m)
          Y_lm(l,m) = P_lm * dcos(m_dbl*phi) !+phase_lm(l,m))  ! Computing spherical harmonics from ALF
        End Do
      End Do

      Return

      End Subroutine Spherical_Harmonics
