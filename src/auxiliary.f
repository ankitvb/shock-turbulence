************************************************************************
* This file contains auxiliary routine(s) which do not inherently belong
* to a specific part of the program.
************************************************************************

************************************************************************
* Subroutine: Compute_Mean
*
* Computes the mean value of the elements in F, which can be distributed
* over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_mean: mean value
************************************************************************
      Subroutine Compute_Mean(F, F_mean)

      Include 'header'
      Include 'mpif.h'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_mean

*  Locals
      Integer I, J, K, ierr
      Real*8  F_sum

*  Sum the F-field over the entire domain
      F_sum = zero

*  ..Sum over the local domain
      Do I = 1,Mx
         Do J = 1,My
            Do K = 1,Mz
                F_sum = F_sum + F(J,K,I)
            End Do
         End Do
      End Do

*  ..Sum over all processors
      Call MPI_ALLREDUCE(F_sum, F_mean, 1, MPI_REAL8, MPI_SUM, 
     .                   MPI_COMM_WORLD, ierr)

*  Compute the average of F-field
      F_mean = F_mean / Dble(Ny*Nz*Nx)

      Return

      End

************************************************************************
* Subroutine: Compute_RMS
*
* Computes the rms value of the elements in F, which can be distributed
* over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_rms: rms value
************************************************************************
      Subroutine Compute_RMS(F, F_RMS)

      Include 'header'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_rms

*  Locals
      Integer I, J, K, ierr
      Real*8  Ftemp(My,Mz,Mx)
      Real*8  F_mean

      Call Compute_Mean(F, F_mean)
      
*  ..Modify the local domain
      Do I = 1,Mx
         Do J = 1,My
            Do K = 1,Mz
c               Get fluctuating quantities for the subdomain                
                Ftemp(J,K,I) = F(J,K,I) - F_mean
c               Square them
                Ftemp(J,K,I) = Ftemp(J,K,I) * Ftemp(J,K,I)       
            End Do
         End Do
      End Do

c     Take mean of the squared fluctuating quantities
      Call Compute_Mean(Ftemp, F_rms)

c     Take square root to get the RMS quantity
      F_rms = dsqrt(F_rms)

      Return
    
      End

************************************************************************
* Subroutine: Compute_shell_RMS
*
* Computes the rms value of the elements in F, which can be distributed
* over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_rms: rms value
************************************************************************
      Subroutine Compute_shell_RMS(F, F_RMS)

      Include 'header'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_rms(N_avg)

*  Locals
      Integer I, J, K, M, ierr
      Real*8  Ftemp(My,Mz,Mx)
      Real*8  F_mean(N_avg)

      Integer II, JJ, KK
      Real*8  IX, IY, IZ
      Real*8  rad, r(N_avg+1)


      Call Compute_Shell_Avg(F, F_mean)

* Setting bin levels
* There are N/4 + 4 levels, spaced Delta * sqrt(3) apart in r-space

      r(1) = DX * sqrt(3D0) / 2D0 - 1d-4

      Do M = 2,N_avg+1
        if(M .LE. 4)then
          r(M) = r(M-1) + DX * 2d0 * sqrt(3D0)
        else
          r(M) = r(M-1) + DX * sqrt(3D0)
        endif
      End Do

*  ..Modify the local domain
      Do I = 1,Mx
         Do J = 1,My
            Do K = 1,Mz
                II = xrank * Mx + I
                JJ = yrank * My + J
                KK = zrank * Mz + K

                IX = DBLE(II - Nx/2) - 0.5D0
                IY = DBLE(JJ - Ny/2) - 0.5D0
                IZ = DBLE(KK - Nz/2) - 0.5D0

                rad = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)
                Do M = 1,N_avg
                  if((rad.GE.r(M)).AND.(rad.LT.r(M+1)))then
                    Ftemp(J,K,I) = F(J,K,I) - F_mean(M) ! Get fluctuating quantities for the subdomain
                    Ftemp(J,K,I) = Ftemp(J,K,I) * Ftemp(J,K,I)
                  endif
                End Do
            End Do
         End Do
      End Do

c     Take mean of the squared fluctuating quantities
      Call Compute_Shell_Avg(Ftemp, F_rms)

c     Take square root to get the RMS quantity
      F_rms = dsqrt(F_rms)

      Return

      End


************************************************************************
* Subroutine: Compute_Sum
*
* Computes the sum of the elements in F, which can be distributed
* over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_sum: sum
************************************************************************
      Subroutine Compute_Sum(F, F_sum)

      Include 'header'
      Include 'mpif.h'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_sum, F_temp

*  Locals
      Integer I, J, K, ierr

*  Sum the F-field over the entire domain
      F_temp = zero

*  ..Sum over the local domain
      Do I = 1,Mx
         Do J = 1,My
             Do K = 1,Mz
                F_temp = F_temp + F(J,K,I)
             End Do
         End Do 
      End Do

*  ..Sum over all processors
      Call MPI_ALLREDUCE(F_temp, F_sum, 1, MPI_REAL8, MPI_SUM, 
     .                   MPI_COMM_WORLD, ierr)

      Return

      End

************************************************************************
* Subroutine: Compute_Zavg
*
* Computes the 2D, z-averaged form of a 3D F, which can be distributed
* over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_zavg: z-averaged F
************************************************************************
      Subroutine Compute_Zavg(F, F_zavg)

      Include 'header'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_zavg(My,Mx)
*  Locals
      Integer I, J, K

      do I = 1, Mx
         do J = 1, My
            F_zavg(J,I) = sum(F(J,:,I))
         enddo
      enddo
      F_zavg = F_zavg / Mz

      Return
      End

************************************************************************
* Subroutine: Compute_shell_avg
*
* Computes the shell-averaged form of a 3D spherical F, which can be distributed
* over several nodes. The data is binned appropriately
*
* In: F: array with dimensions my, mz, mx 
*
* !!!!! WARNING : ROUTINE ASSUMES Nx = Ny = Nz !!!!!!!!!!!!     
*
* Out: F_zavg: z-averaged F in radial bins
************************************************************************
      Subroutine Compute_shell_avg(F, F_shell_avg)

      Include 'header'
      Include 'mpif.h'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_shell_avg(N_avg)
      Real*8  F_local_sum(N_avg)

*  Locals
      Integer I, J, K, M
      Integer II, JJ, KK
      Real*8  IX, IY, IZ
      Real*8  rad, r(N_avg+1)
      Integer count_local(N_avg), count_sum(N_avg)
      Integer ierr

* Setting bin levels
* There are N/4 + 4 levels, spaced Delta * sqrt(3) apart in r-space

      r(1) = DX * sqrt(3D0) / 2D0 - 1d-4

      Do M = 2,N_avg+1
        if(M .LE. 4)then
          r(M) = r(M-1) + DX * 2d0 * sqrt(3D0)
        else
          r(M) = r(M-1) + DX * sqrt(3D0)
        endif 
      End Do

* Initializing shell avg for each processor
      F_local_sum(:) = 0D0
      F_shell_avg(:) = 0D0
      count_local(:) = 0
      count_sum(:)   = 0

* Computing the average for each processor
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            II = xrank * Mx + I
            JJ = yrank * My + J
            KK = zrank * Mz + K

            IX = DBLE(II - Nx/2) - 0.5D0
            IY = DBLE(JJ - Ny/2) - 0.5D0
            IZ = DBLE(KK - Nz/2) - 0.5D0

            rad = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

            Do M = 1,N_avg      
             if((rad.GE.r(M)).AND.(rad.LT.r(M+1)))then
              F_local_sum(M) = F_local_sum(M) + F(J,K,I)
              count_local(M) = count_local(M) + 1
             endif
            End Do
          End Do
        End Do
      End Do 

* Summing up over all processors for each bin
      Call MPI_ALLREDUCE(F_local_sum, F_shell_avg, N_avg, MPI_REAL8,
     .                   MPI_SUM, MPI_COMM_WORLD, ierr)
      Call MPI_ALLREDUCE(count_local, count_sum, N_avg, MPI_INTEGER,
     .                   MPI_SUM, MPI_COMM_WORLD, ierr)

* Computing shell average
      Do M = 1,N_avg
        F_shell_avg(M) = F_shell_avg(M) / DBLE(count_sum(M))
      End Do     

      Return

      End

************************************************************************
* Subroutine: Compute_Max
*
* Computes the maximum value of the elements in F, which can be
* distributed over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_max: max value
************************************************************************
      Subroutine Compute_Max(F, F_max)

      Include 'header'
      Include 'mpif.h'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_max

*  Locals
      Integer I, ierr
      Real*8  F_localmax

*  Compute the maximum value on this processor
      F_localmax = maxval(F)

*  Find the maximum value over all processors
      Call MPI_ALLREDUCE(F_localmax, F_max, 1, MPI_REAL8, MPI_MAX, 
     .                   MPI_COMM_WORLD, ierr)

      End


************************************************************************
* Subroutine: Compute_Min
*
* Computes the minimum value of the elements in F, which can be
* distributed over several nodes.
*
* In: F: array with dimensions ny, nz, mx
*     mx, ny, nz: array dimensions
*
* Out: F_min: min value
************************************************************************
      Subroutine Compute_Min(F, F_min)

      Include 'header'
      Include 'mpif.h'

*  Procedure parameters
      Real*8  F(My,Mz,Mx), F_min

*  Locals
      Integer I, ierr
      Real*8  F_localmin

*  Compute the minimum value on this processor
      F_localmin = minval(F)

*  Find the minimum value over all processors
      Call MPI_ALLREDUCE(F_localmin, F_min, 1, MPI_REAL8, MPI_MIN, 
     .                   MPI_COMM_WORLD, ierr)

      End


************************************************************************
*  Subroutine: INT2STRING
*
*  Converts a non-negative integer to a string.
*
*  In:
*    number: integer >= 0
*
*  Out :
*    snumber: 25-character-long string containing the number in the
*             front and the rest padded with blanks.
************************************************************************
      subroutine int2string(number, snumber)
      implicit none

*  Procedure parameters
      character snumber*25
      integer number

*  Locals
      character srev*25
      integer inumber, digit, ndig, slen, i, imax

*  Check input value
      inumber = number
      if (inumber < 0) STOP '*** invalid integer'

*  Initialize strings with blanks
      do i = 1, 25
        snumber(i:i) = ' '
      end do

*  Add each digit consecutively to the string.
      i = 1
      digit = mod(inumber, 10)
      do while (inumber > 9)
        srev(i:i) = char(ichar('0')+digit)
        inumber = (inumber - digit) / 10
        i = i+1
        digit = mod(inumber, 10)
      end do
      srev(i:i) = char(ichar('0')+digit)
      imax = i
      do i = 1, imax
        snumber(i:i) = srev(imax+1-i:imax+1-i)
      end do

      end


************************************************************************
* Subroutine REAL2STRING
*
* Takes a real value r and writes its equivalent into the string s in
* the form xxxx.yyy, where ldigits is the number of digits before the
* decimal point and rdigits is the number of digits after it.
* r must be non-negative and smaller than 10^ldigits, otherwise a
* blank string is returned.
*
* In: r (real value)
* Out: s (string of length slen)
*
* NB: The rounding does not always work hypercorrect.
************************************************************************
      subroutine real2string(r, s)
        implicit none

*  Local Parameters
        integer slen, ldigits, rdigits
        parameter (slen = 15, ldigits = 4, rdigits = 3)

*  Procedure parameters
        real*8 r
        character*(slen) s

*  Locals
        integer i, pos, curdig
        real*8 x, curval

*  Initialize string with blanks and copy input value to local variable
        do i = 1, slen
          s(i:i) = ' '
        end do
        x = r

*  Round
        if ((10**rdigits)*x - 10d0*floor((10**(rdigits-1))*x) >= 5d-1)
     .  then
          x = x + 5d-1 * 10d0**(-rdigits)
        end if

*  Check input value
        if (x < 0 .or. x >= 10**(ldigits)) return

*  Copy numeric value from real to string
        pos = 1
        curval = 10**(ldigits-1)
        do i = (ldigits-1), -rdigits, -1
          curdig = aint(x / curval)
          s(pos:pos) = char(ichar('0') + curdig)
          x = x - curdig * curval
          curval = curval / 10d0
          if (i == 0) then
            pos = pos + 1
            s(pos:pos) = '.'
          end if
          pos = pos + 1
        end do

      end

************************************************************************
* Subroutine: CROSS_PRODUCT
*
* Computes the cross-product c = a x b
*
* In: a,b: vectors of length 3 and type real*8
*
* Out: c: cross-product of a and b (itself a vector of length 3)
************************************************************************
      Subroutine Cross_Product(A,B,C)

      Include 'header'

      Real*8 A(My,Mz,Mx,3), B(My,Mz,Mx,3), C(My,Mz,Mx,3)
      Integer I, J, K, L

      Do I = 1,Mx
        Do J = 1,Ny
          Do K = 1,Nz
            Do L = 1,3
              C(J,K,I,1) = A(J,K,I,2)*B(J,K,I,3) - A(J,K,I,3)*B(J,K,I,2)
              C(J,K,I,2) = A(J,K,I,3)*B(J,K,I,1) - A(J,K,I,1)*B(J,K,I,3)
              C(J,K,I,3) = A(J,K,I,1)*B(J,K,I,2) - A(J,K,I,2)*B(J,K,I,1)
            End Do
          End Do
        End Do
      End Do

      End Subroutine Cross_Product
*************************************************************************
