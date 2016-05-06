c----------------------------------------------------------------------
      subroutine Init_Part

      include 'header'
      include 'part_header'
      include 'mpif.h'
      
      real*8  Q(My,Mz,Mx,5)
      integer I, J, K, M, N, pcount, count
      integer II, JJ, KK
      real*8  IX, IY, IZ
      real*8,  dimension(:,:), allocatable  :: W
      real*8,  dimension(:), allocatable    :: T, p
      integer, dimension(:), allocatable    :: WStag
      integer, dimension(3)                 :: coords
      real*8  rad, tag
      integer ndims

      Common /RKWORK/ Q

      Integer  status(MPI_STATUS_SIZE), req(2), ierr

* Rank 0 processor allocates particles to other processors
* Computing radius at every point in the domain (origin at center) 
      ndims = 3
      if(rank == 0)then
        Np(:) = 0
        Do M = 0,Px*Py*Pz-1
         Do I = 1,Mx
          Do J = 1,My
           Do K = 1,Mz
            Call MPI_CART_COORDS(MPI_XYZ_COMM, M, ndims, coords, ierr)
            II =  coords(3) * Mx + I
            JJ =  coords(2) * My + J
            KK =  coords(1) * Mz + K

            IX = REAL(II - Nx/2) - 0.5D0 
            IY = REAL(JJ - Ny/2) - 0.5D0 
            IZ = REAL(KK - Nz/2) - 0.5D0 

            rad = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

            if((rad .GE. (0.15D0*x_leng)) .AND.      ! Converging 0.1 for expanding shock
     .         (rad .LE. (0.1525D0*x_leng)))then
                  Np(M) = Np(M) + 1
            endif
           End Do
          End Do
         End Do
        End Do
      endif

* Debug stuff
      if(rank .EQ. 0)then
        Do M = 0,Px*Py*Pz-1
          print *, Np(M), M
        End Do
      endif

* Broadcast number of particles on each processor
* Initialize Np, Nl, Nr
      Call MPI_BCAST(Np, Px*Py*Pz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      Nxlo = 0; Nxli = 0; Nxro = 0; Nxri = 0
      Nylo = 0; Nyli = 0; Nyro = 0; Nyri = 0
      Nzlo = 0; Nzli = 0; Nzro = 0; Nzri = 0

      Npart = Np(rank)
      NpMax = maxval(Np)

* Every processor allocates particle arrays
      if(Npart>0)then
        allocate(part(Npart,NpDof),STAT=ierr)
        if(ierr .ne. 0) Stop '*** Particle array allocation failed!'
        allocate(dpart(Npart,3),STAT=ierr)
        if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
        allocate(partW(Npart,3),STAT=ierr)
        if(ierr .ne. 0) Stop '*** Particle array allocation failed!'
      endif

* Processor 0 initializes particles with their tags and sends them to
* their respective processors
      if(rank .EQ. 0)then
        tag = 0D0
        Do M = 0,Px*Py*Pz-1
          if(Np(M) > 0)then
            allocate(W(Np(M),NpDof),STAT=ierr)
            if(ierr .ne. 0) STOP '*** Particle array allocation failed!'

            Do N = 1,Np(M)
              tag = tag + 1D0
              W(N,NpDof) = tag
            End Do
            Call MPI_SSEND(W, Np(M)*NpDof, MPI_REAL8, M, M,
     .                     MPI_COMM_WORLD, ierr)

            if(ALLOCATED(W)) deallocate(W)            
          endif
        End Do
      endif

* All processors recieve particle arrays from root processor
      if(Npart > 0)then
        Call MPI_RECV(part, Np(rank)*NpDof, MPI_REAL8, 0, rank,
     .                MPI_COMM_WORLD, status, ierr)
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
      if(Npart > 0)then
* Processors update particles with local field information
        count = 0
        Do I = 1,Mx
          Do J = 1,My
            Do K = 1,Mz
              II =  xrank * Mx + I
              JJ =  yrank * My + J
              KK =  zrank * Mz + K

              IX = REAL(II - Nx/2) - 0.5D0
              IY = REAL(JJ - Ny/2) - 0.5D0
              IZ = REAL(KK - Nz/2) - 0.5D0

              rad = dsqrt((IX*DX)**2D0+(IY*DY)**2D0+(IZ*DZ)**2D0)

              if((rad .GE. (0.15D0*x_leng)) .AND.
     .           (rad .LE. (0.1525D0*x_leng)))then
                  count = count + 1
                  part(count,1) = xcor(I)
                  part(count,2) = ycor(J)
                  part(count,3) = zcor(K)
                  part(count,4) = Q(J,K,I,4)
              endif
            End Do
          End Do
        End Do    
      
      endif 

      Return
      
      End subroutine Init_Part
         
c-----------------------------------------------------------------
      subroutine interpolate(Q)

      include 'header'
      include 'part_header'

      real*8, dimension(My,Mz,Mx,5)    :: Q
      real*8, dimension(:,:), pointer  :: Qp
      real*8, dimension(My,Mz,NpDof-1) :: Qxbuf
      real*8, dimension(Mz,Mx,NpDof-1) :: Qybuf
      real*8, dimension(My,Mx,NpDof-1) :: Qzbuf      
      real*8, dimension(Mz,NpDof-1)    :: Qxybuf
      real*8, dimension(Mx,NpDof-1)    :: Qyzbuf
      real*8, dimension(My,NpDof-1)    :: Qxzbuf
      real*8, dimension(NpDof-1)       :: Qxyzbuf

      COMMON /MIDPART/ Qp /MIDPART_BUF/ Qxbuf, Qybuf, Qzbuf,
     .                                  Qxybuf, Qyzbuf, Qxzbuf,
     .                                  Qxyzbuf

      real*8, dimension(My,Mz,Mx,2) :: pGrad

      Common /PART_GRAD/ pGrad
      
      real*8, dimension(My+1,Mz+1,Mx+1,NpDof-1):: W     !Work array for Q
      real*8, dimension(:), pointer            :: x,y,z
      integer,dimension(:), allocatable        :: ix,iy,iz
      real*8, dimension(:), allocatable        :: xs,ys,zs !Scaled for trilin. int.
           
      integer I, J, K, L, n, err

      if(Npart > 0)then
        allocate(ix(Npart),iy(Npart),iz(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'
        allocate(xs(Npart),ys(Npart),zs(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'

* Converting first three DOFs of W to velocities from momenta, fifth DOF to P
        Do I = 1,Mx
          Do J = 1,My
            Do K = 1,Mz
              W(J,K,I,1) = Q(J,K,I,1)/Q(J,K,I,4)
              W(J,K,I,2) = Q(J,K,I,2)/Q(J,K,I,4)
              W(J,K,I,3) = Q(J,K,I,3)/Q(J,K,I,4)
              W(J,K,I,4) = Q(J,K,I,4)
              W(J,K,I,5) = (gamma-1D0) *
     .                     ( Q(J,K,I,5) - (5.0D-1 / Q(J,K,I,4)) *
     .                     ( Q(J,K,I,1) * Q(J,K,I,1)
     .                     + Q(J,K,I,2) * Q(J,K,I,2)
     .                     + Q(J,K,I,3) * Q(J,K,I,3) ) )
              W(J,K,I,6) = pGrad(J,K,I,1)
              W(J,K,I,7) = pGrad(J,K,I,2)
            End Do
          End Do
        End Do

        do J = 1, My
          do K = 1, Mz
            W(J,K,Mx+1,1) = Qxbuf(J,K,1)/Qxbuf(J,K,4)
            W(J,K,Mx+1,2) = Qxbuf(J,K,2)/Qxbuf(J,K,4)
            W(J,K,Mx+1,3) = Qxbuf(J,K,3)/Qxbuf(J,K,4)
            W(J,K,Mx+1,4) = Qxbuf(J,K,4)
            W(J,K,Mx+1,5) = (gamma-1D0) *
     .                     ( Qxbuf(J,K,5) - (5.0D-1/Qxbuf(J,K,4)) *
     .                     ( Qxbuf(J,K,1) * Qxbuf(J,K,1)
     .                     + Qxbuf(J,K,2) * Qxbuf(J,K,2)
     .                     + Qxbuf(J,K,3) * Qxbuf(J,K,3) ) )
            W(J,K,Mx+1,6) = Qxbuf(J,K,6)
            W(J,K,Mx+1,7) = Qxbuf(J,K,7)
          enddo
        enddo
        do K = 1, Mz
          do I = 1, Mx
            W(My+1,K,I,1) = Qybuf(K,I,1)/Qybuf(K,I,4)
            W(My+1,K,I,2) = Qybuf(K,I,2)/Qybuf(K,I,4)
            W(My+1,K,I,3) = Qybuf(K,I,3)/Qybuf(K,I,4)
            W(My+1,K,I,4) = Qybuf(K,I,4)
            W(My+1,K,I,5) = (gamma-1D0) *
     .                     ( Qybuf(K,I,5) - (5.0D-1/Qybuf(K,I,4)) *
     .                     ( Qybuf(K,I,1) * Qybuf(K,I,1)
     .                     + Qybuf(K,I,2) * Qybuf(K,I,2)
     .                     + Qybuf(K,I,3) * Qybuf(K,I,3) ) )
            W(My+1,K,I,6) = Qybuf(K,I,6)
            W(My+1,K,I,7) = Qybuf(K,I,7)
          enddo
        enddo
        do J = 1, My
          do I = 1, Mx
            W(J,Mz+1,I,1) = Qzbuf(J,I,1)/Qzbuf(J,I,4)
            W(J,Mz+1,I,2) = Qzbuf(J,I,2)/Qzbuf(J,I,4)
            W(J,Mz+1,I,3) = Qzbuf(J,I,3)/Qzbuf(J,I,4)
            W(J,Mz+1,I,4) = Qzbuf(J,I,4)
            W(J,Mz+1,I,5) = (gamma-1D0) *
     .                     ( Qzbuf(J,I,5) - (5.0D-1/Qzbuf(J,I,4)) *
     .                     ( Qzbuf(J,I,1) * Qzbuf(J,I,1)
     .                     + Qzbuf(J,I,2) * Qzbuf(J,I,2)
     .                     + Qzbuf(J,I,3) * Qzbuf(J,I,3) ) )
            W(J,Mz+1,I,6) = Qzbuf(J,I,6)
            W(J,Mz+1,I,7) = Qzbuf(J,I,7)
          enddo
        enddo
 
        do K = 1,Mz
          W(My+1,K,Mx+1,1) = Qxybuf(K,1)/Qxybuf(K,4)
          W(My+1,K,Mx+1,2) = Qxybuf(K,2)/Qxybuf(K,4)
          W(My+1,K,Mx+1,3) = Qxybuf(K,3)/Qxybuf(K,4)
          W(My+1,K,Mx+1,4) = Qxybuf(K,4)
          W(My+1,K,Mx+1,5) = (gamma-1D0)*
     .                     ( Qxybuf(K,5) - (5.0D-1/Qxybuf(K,4)) *
     .                     ( Qxybuf(K,1) * Qxybuf(K,1)
     .                     + Qxybuf(K,2) * Qxybuf(K,2)
     .                     + Qxybuf(K,3) * Qxybuf(K,3) ) )
          W(My+1,K,Mx+1,6) = Qxybuf(K,6)
          W(My+1,K,Mx+1,7) = Qxybuf(K,7)
        enddo

        do I = 1,Mx
          W(My+1,Mz+1,I,1) = Qyzbuf(I,1)/Qyzbuf(I,4)
          W(My+1,Mz+1,I,2) = Qyzbuf(I,2)/Qyzbuf(I,4)
          W(My+1,Mz+1,I,3) = Qyzbuf(I,3)/Qyzbuf(I,4)
          W(My+1,Mz+1,I,4) = Qyzbuf(I,4)
          W(My+1,Mz+1,I,5) = (gamma-1D0)*
     .                     ( Qyzbuf(I,5) - (5.0D-1/Qyzbuf(I,4)) *
     .                     ( Qyzbuf(I,1) * Qyzbuf(I,1)
     .                     + Qyzbuf(I,2) * Qyzbuf(I,2)
     .                     + Qyzbuf(I,3) * Qyzbuf(I,3) ) )
          W(My+1,Mz+1,I,6) = Qyzbuf(I,6)
          W(My+1,Mz+1,I,7) = Qyzbuf(I,7)
        enddo

        do J = 1,My
          W(J,Mz+1,Mx+1,1) = Qxzbuf(J,1)/Qxzbuf(J,4)
          W(J,Mz+1,Mx+1,2) = Qxzbuf(J,2)/Qxzbuf(J,4)
          W(J,Mz+1,Mx+1,3) = Qxzbuf(J,3)/Qxzbuf(J,4)
          W(J,Mz+1,Mx+1,4) = Qxzbuf(J,4)
          W(J,Mz+1,Mx+1,5) = (gamma-1D0)*
     .                     ( Qxzbuf(J,5) - (5.0D-1/Qxzbuf(J,4)) *
     .                     ( Qxzbuf(J,1) * Qxzbuf(J,1)
     .                     + Qxzbuf(J,2) * Qxzbuf(J,2)
     .                     + Qxzbuf(J,3) * Qxzbuf(J,3) ) )
          W(J,Mz+1,Mx+1,6) = Qxzbuf(J,6)
          W(J,Mz+1,Mx+1,7) = Qxzbuf(J,7)
        enddo

        W(My+1,Mz+1,Mx+1,1) = Qxyzbuf(1)/Qxyzbuf(4)
        W(My+1,Mz+1,Mx+1,2) = Qxyzbuf(2)/Qxyzbuf(4)
        W(My+1,Mz+1,Mx+1,3) = Qxyzbuf(3)/Qxyzbuf(4)
        W(My+1,Mz+1,Mx+1,4) = Qxyzbuf(4)
        W(My+1,Mz+1,Mx+1,5) = (gamma-1D0)*
     .                     ( Qxyzbuf(5) - (5.0D-1/Qxyzbuf(4)) *
     .                     ( Qxyzbuf(1) * Qxyzbuf(1)
     .                     + Qxyzbuf(2) * Qxyzbuf(2)
     .                     + Qxyzbuf(3) * Qxyzbuf(3) ) )
        W(My+1,Mz+1,Mx+1,6) = Qxyzbuf(6)
        W(My+1,Mz+1,Mx+1,7) = Qxyzbuf(7)

        x => part(:,1); y => part(:,2); z => part(:,3)

        ix = FLOOR((x-xcor(1))/DX)+1 
        iy = FLOOR((y-ycor(1))/DY)+1
        iz = FLOOR((z-zcor(1))/DZ)+1

        xs = (x - xcor(ix))/DX
        ys = (y - ycor(iy))/DY
        zs = (z - zcor(iz))/DZ

        do n = 1,Npart
          Qp(n,:) = W(iy(n),iz(n),ix(n),:)*(1D0-ys(n))*(1D0-zs(n))
     .          *(1D0-xs(n))
     .          + W(iy(n)+1,iz(n),ix(n),:)*ys(n)*(1D0-zs(n))*(1D0-xs(n))
     .          + W(iy(n),iz(n)+1,ix(n),:)*(1D0-ys(n))*zs(n)*(1D0-xs(n))
     .          + W(iy(n),iz(n),ix(n)+1,:)*(1D0-ys(n))*(1D0-zs(n))*xs(n)
     .          + W(iy(n)+1,iz(n),ix(n)+1,:)*ys(n)*(1D0-zs(n))*xs(n)
     .          + W(iy(n),iz(n)+1,ix(n)+1,:)*(1D0-ys(n))*zs(n)*xs(n)
     .          + W(iy(n)+1,iz(n)+1,ix(n),:)*ys(n)*zs(n)*(1D0-xs(n))
     .          + W(iy(n)+1,iz(n)+1,ix(n)+1,:)*ys(n)*zs(n)*xs(n)
        enddo
      
        DEALLOCATE(ix,iy,iz,xs,ys,zs,STAT=err)
          if(err .ne. 0) STOP '*** Particle array deallocation failed!'
      
      endif

      end subroutine interpolate

c-----------------------------------------------------------------
      subroutine Prepare_Buffer(Q)

      include 'header'
      include 'part_header'
      include 'mpif.h'

      real*8, dimension(My,Mz,Mx,5)       :: Q
      real*8, dimension(:,:), pointer     :: Qp
      real*8, dimension(My,Mz,Mx,NpDof-1) :: Q2
      real*8, dimension(My,Mz,NpDof-1)    :: Qxbuf
      real*8, dimension(Mz,Mx,NpDof-1)    :: Qybuf
      real*8, dimension(My,Mx,NpDof-1)    :: Qzbuf
      real*8, dimension(Mz,NpDof-1)       :: Qxybuf
      real*8, dimension(Mx,NpDof-1)       :: Qyzbuf
      real*8, dimension(My,NpDof-1)       :: Qxzbuf
      real*8, dimension(NpDof-1)          :: Qxyzbuf

      COMMON /MIDPART/ Qp /MIDPART_BUF/ Qxbuf, Qybuf, Qzbuf,
     .                                  Qxybuf, Qyzbuf, Qxzbuf,
     .                                  Qxyzbuf
      
      real*8, dimension(My,Mz,NpDof-1)  :: Wx  
      real*8, dimension(Mz,Mx,NpDof-1)  :: Wy
      real*8, dimension(My,Mx,NpDof-1)  :: Wz
      real*8, dimension(Mz,NpDof-1)     :: Wxy
      real*8, dimension(Mx,NpDof-1)     :: Wyz
      real*8, dimension(My,NpDof-1)     :: Wxz
      real*8, dimension(NpDof-1)        :: Wxyz

      real*8, dimension(My,Mz,Mx,2) :: pGrad

      Common /PART_GRAD/ pGrad

      Integer x_source, y_source, z_source
      Integer x_dest, y_dest, z_dest               
      Integer xy_source, yz_source, xz_source, xyz_source
      Integer xy_dest, yz_dest, xz_dest, xyz_dest
      Integer status_arr(MPI_STATUS_SIZE,3), status(MPI_STATUS_SIZE)
      Integer coords(3), ierr, req(4)
      integer I, J, K, L

      Q2(:,:,:,1:5) = Q(:,:,:,1:5)
      Q2(:,:,:,6)   = pGrad(:,:,:,1)
      Q2(:,:,:,7)   = pGrad(:,:,:,2)

      do L = 1,NpDof-1 
        do K = 1, Mz
          do J = 1, My
             Wx(J,K,L) = Q2(J,K,1,L)
          enddo
          Wxy(K,L) = Q2(1,K,1,L)
        enddo
        do I = 1,Mx
           do K = 1,Mz
              Wy(K,I,L) = Q2(1,K,I,L)
           enddo
           Wyz(I,L) = Q2(1,1,I,L)
        enddo
        do J = 1,My
           do I = 1,Mx
              Wz(J,I,L) = Q2(J,1,I,L)
           enddo
           Wxz(J,L) = Q2(J,1,1,L)
        enddo
        Wxyz(L) = Q2(1,1,1,L)   
      enddo

      

* Face, Edge and Corner data
      x_dest = xrank-1
      if(xrank.EQ.0) x_dest = Px-1
      y_dest = yrank-1
      if(yrank.EQ.0) y_dest = Py-1
      z_dest = zrank-1
      if(zrank.EQ.0) z_dest = Pz-1

      x_source = xrank+1
      if(xrank.EQ.Px-1) x_source = 0
      y_source = yrank+1
      if(yrank.EQ.Py-1) y_source = 0
      z_source = zrank+1
      if(zrank.EQ.Pz-1) z_source = 0 

      coords(1) = zrank; coords(2) = yrank-1; coords(3) = xrank-1
      if(yrank.EQ.0) coords(2) = Py-1
      if(xrank.EQ.0) coords(3) = Px-1
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, xy_dest, ierr)    
      coords(1) = zrank-1; coords(2) = yrank-1; coords(3) = xrank     
      if(zrank.EQ.0) coords(1) = Pz-1
      if(yrank.EQ.0) coords(2) = Py-1
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, yz_dest, ierr)
      coords(1) = zrank-1; coords(2) = yrank; coords(3) = xrank-1
      if(zrank.EQ.0) coords(1) = Pz-1
      if(xrank.EQ.0) coords(3) = Px-1
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, xz_dest, ierr)
      coords(1) = zrank-1; coords(2) = yrank-1; coords(3) = xrank-1
      if(xrank.EQ.0) coords(3) = Px-1
      if(yrank.EQ.0) coords(2) = Py-1
      if(zrank.EQ.0) coords(1) = Pz-1
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, xyz_dest, ierr)      

      coords(1) = zrank; coords(2) = yrank+1; coords(3) = xrank+1
      if(yrank.EQ.Py-1) coords(2) = 0
      if(xrank.EQ.Px-1) coords(3) = 0
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, xy_source, ierr)
      coords(1) = zrank+1; coords(2) = yrank+1; coords(3) = xrank
      if(zrank.EQ.Pz-1) coords(1) = 0
      if(yrank.EQ.Py-1) coords(2) = 0
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, yz_source, ierr)
      coords(1) = zrank+1; coords(2) = yrank; coords(3) = xrank+1
      if(zrank.EQ.Pz-1) coords(1) = 0
      if(xrank.EQ.Px-1) coords(3) = 0
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, xz_source, ierr)
      coords(1) = zrank+1; coords(2) = yrank+1; coords(3) = xrank+1
      if(xrank.EQ.Px-1) coords(3) = 0
      if(yrank.EQ.Py-1) coords(2) = 0
      if(zrank.EQ.Pz-1) coords(1) = 0
      Call MPI_CART_RANK(MPI_XYZ_COMM, coords, xyz_source, ierr)

* All processors send face data
        Call MPI_ISEND(Wx, My*Mz*(NpDof-1), MPI_REAL8, x_dest,
     .                 xrank, MPI_X_COMM, req(1), ierr)
        Call MPI_ISEND(Wy, Mx*Mz*(NpDof-1), MPI_REAL8, y_dest,
     .                 yrank, MPI_Y_COMM, req(2), ierr)
        Call MPI_ISEND(Wz, My*Mx*(NpDof-1), MPI_REAL8, z_dest,
     .                 zrank, MPI_Z_COMM, req(3), ierr)

c All processors receive it in buffer array
        Call MPI_RECV(Qxbuf, My*Mz*(NpDof-1), MPI_REAL8, x_source,
     .                x_source, MPI_X_COMM, status, ierr)
        Call MPI_RECV(Qybuf, Mx*Mz*(NpDof-1), MPI_REAL8, y_source,
     .                y_source, MPI_Y_COMM, status, ierr)
        Call MPI_RECV(Qzbuf, My*Mx*(NpDof-1), MPI_REAL8, z_source,
     .                z_source, MPI_Z_COMM, status, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_ISEND(Wxy, Mz*(NpDof-1), MPI_REAL8, xy_dest,
     .               rank, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_ISEND(Wyz, Mx*(NpDof-1), MPI_REAL8, yz_dest,
     .               rank, MPI_COMM_WORLD, req(2), ierr)
      Call MPI_ISEND(Wxz, My*(NpDof-1), MPI_REAL8, xz_dest,
     .               rank, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_ISEND(Wxyz, NpDof-1, MPI_REAL8, xyz_dest,
     .               rank, MPI_COMM_WORLD, req(4), ierr)

* All processors recieve edge and corner data
      Call MPI_RECV(Qxybuf, Mz*(NpDof-1), MPI_REAL8, xy_source,
     .              xy_source, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(Qyzbuf, Mx*(NpDof-1), MPI_REAL8, yz_source,
     .              yz_source, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(Qxzbuf, My*(NpDof-1), MPI_REAL8, xz_source,
     .              xz_source, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(Qxyzbuf, NpDof-1, MPI_REAL8, xyz_source,
     .              xyz_source, MPI_COMM_WORLD, status, ierr)
 
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      end subroutine Prepare_Buffer

c-----------------------------------------------------------------
      Subroutine Prepare_Grad

      include 'header'
    
      real*8, dimension(My,Mz,Mx,2) :: pGrad
      real*8, dimension(My,Mz,Mx,4) :: Grad
 
      Common /PART_GRAD/ pGrad
       
      Real*8 dux_dx(My,Mz,Mx), duy_dx(My,Mz,Mx), duz_dx(My,Mz,Mx)
      Real*8 dux_dy(My,Mz,Mx), duy_dy(My,Mz,Mx), duz_dy(My,Mz,Mx)
      Real*8 dux_dz(My,Mz,Mx), duy_dz(My,Mz,Mx), duz_dz(My,Mz,Mx)

      Common /GRAD/ dux_dx, duy_dx, duz_dx,
     .              dux_dy, duy_dy, duz_dy,
     .              dux_dz, duy_dz, duz_dz

      Integer I, J, K

      Grad(:,:,:,1) = duz_dy - duy_dz

      Grad(:,:,:,2) = dux_dz - duz_dx

      Grad(:,:,:,3) = duy_dx - dux_dy

      Grad(:,:,:,4) = dux_dx + duy_dy + duz_dz

      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            pGrad(J,K,I,1) = dsqrt(Grad(J,K,I,1)**2d0 +
     .                       Grad(J,K,I,2)**2d0+ Grad(J,K,I,3)**2d0)
             pGrad(J,K,I,2) = Grad(J,K,I,4)
          End Do
        End Do
      End Do
       
      End Subroutine Prepare_Grad      

c-----------------------------------------------------------------
      subroutine Part_rhs(Q)

      include 'header'
      include 'part_header'

      real*8, dimension(My,Mz,Mx,5)   :: Q
      real*8, dimension(:,:), pointer :: Qp
      real*8, dimension(My,Mz,5)      :: Qxbuf
      real*8, dimension(Mz,Mx,5)      :: Qybuf
      real*8, dimension(My,Mx,5)      :: Qzbuf

      COMMON /MIDPART/ Qp /MIDPART_BUF/ Qxbuf, Qybuf, Qzbuf
 
      real*8, dimension(My,Mz,Mx,2) :: pGrad

      Common /PART_GRAD/ pGrad             
            
      integer i, j, k, L
      integer err

c Prepare the buffer vector for each processor
      call Prepare_Grad
      call Prepare_Buffer(Q)

* Allocating Qp, DOFs are u, v, w, rho, P
      if(Npart>0)then
        allocate(Qp(Npart,NpDof-1),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'

        call interpolate(Q)
 
        dpart = 0D0
        dpart(:,1:3) = Qp(:,1:3)

        part(:,4:NpDof-1) = Qp(:,4:NpDof-1)

        deallocate(Qp, STAT=err)
        if(err .ne. 0) STOP '*** Particle array deallocation failed!'
      endif

      return

      end subroutine Part_rhs

*************************************************************************************
* Subroutine: Reorder_Part                                                          *
*                                                                                   *
* Function: To move particles from one processor to another after each              *
*           particle time step                                                      *
*************************************************************************************       
      subroutine Reorder_Part

      include 'header'
      include 'part_header'
      include 'mpif.h'

      
      integer, dimension(:), allocatable   :: ix, iy, iz
      real*8,  dimension(:), allocatable   :: x, y, z
      real*8,  dimension(:,:), allocatable :: Wx, WxLI, WxLO, WxRI, WxRO
      real*8,  dimension(:,:), allocatable :: Wy, WyLI, WyLO, WyRI, WyRO
      real*8,  dimension(:,:), allocatable :: Wz, WzLI, WzLO, WzRI, WzRO
      integer, dimension(:), allocatable   :: index_Wx, index_WxLO, 
     .                                        index_WxRO
      integer, dimension(:), allocatable   :: index_Wy, index_WyLO,
     .                                        index_WyRO
      integer, dimension(:), allocatable   :: index_Wz, index_WzLO,
     .                                        index_WzRO
      integer, dimension(3)                :: coords
      integer  n, L, I, J, K, M 
      integer  IP, JP, KP
      integer  irank, irank_left, irank_right, irank_lbdy, irank_rbdy
      integer  count_Wx, count_WxLO, count_WxRO
      integer  count_Wy, count_WyLO, count_WyRO
      integer  count_Wz, count_WzLO, count_WzRO
      integer  Npart_old
      integer  status(MPI_STATUS_SIZE)
      integer  status_arr(MPI_STATUS_SIZE,2), ierr, req(2), err
      integer  PxMax 

      PxMax = 256

!---------------------------------------------------------------------------------------------!
!                                          X                                                  !
!---------------------------------------------------------------------------------------------!                
      
      if(Npart>0)then
        allocate(ix(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'
        allocate(x(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'    

         x           = part(:,1)
         ix          = FLOOR((x-xcor(1))/DX)+1
         Nxlo(rank)  = COUNT(MASK=(ix<1))
         Nxro(rank)  = COUNT(MASK=(ix>Mx))
         n           = COUNT(MASK=(ix>=1 .and. ix<=Mx)) 
      else
         n = 0
         Nxlo(rank) = 0
         Nxro(rank) = 0
      endif

      if(rank>0)then
        Call MPI_ISEND(Nxlo(rank),1, MPI_INTEGER, 0, PxMax+rank,
     .                 MPI_COMM_WORLD, req(1), ierr)
        Call MPI_ISEND(Nxro(rank),1, MPI_INTEGER, 0, 2*PxMax+rank,
     .                 MPI_COMM_WORLD, req(2), ierr)
      endif
      

      if(rank == 0)then
        Do M = Px*Py*Pz-1, 1, -1
          Call MPI_RECV(Nxlo(M), 1, MPI_INTEGER, M, PxMax+M,
     .                  MPI_COMM_WORLD, status, ierr)
          Call MPI_RECV(Nxro(M), 1, MPI_INTEGER, M, 2*PxMax+M,
     .                  MPI_COMM_WORLD, status, ierr)
        Enddo
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c-------------------------------------------------------------
c Prepare and BCAST Nli, Nri, Nro, Nlo, updated Np
      if(rank==0)then
        Do IP = 0,Px-1
          Do JP = 0,Py-1
            Do KP = 0,Pz-1
              coords(1) = KP; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank,ierr)
              coords(1) = KP; coords(2) = JP; coords(3) = IP-1
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_left,ierr)
              coords(1) = KP; coords(2) = JP; coords(3) = IP+1
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_right,ierr)
              coords(1) = KP; coords(2) = JP; coords(3) = 0 
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_lbdy,ierr)
              coords(1) = KP; coords(2) = JP; coords(3) = Px-1
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_rbdy,ierr)

              if(IP>0)then
                 Nxli(irank) = Nxro(irank_left)
              else
                 Nxli(irank) = Nxro(irank_rbdy)
              endif

              if(IP<Px-1)then
                 Nxri(irank) = Nxlo(irank_right)
              else
                 Nxri(irank) = Nxlo(irank_lbdy)
              endif
         
              Np(irank) = Np(irank) 
     .                  - Nxro(irank) - Nxlo(irank) 
     .                  + Nxri(irank) + Nxli(irank)
            Enddo
          Enddo
        Enddo 
      endif

      Call MPI_BCAST(Np, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)
     
      Call MPI_BCAST(Nxlo, Px*Py*Pz, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nxro, Px*Py*Pz, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nxli, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nxri, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)

         
      NpMax = maxval(Np)
      Npart_old = Npart
      Npart = Np(rank)

c      print *, xrank, Np(rank), Nxli(rank), Nxlo(rank),
c     .                Nxri(rank), Nxro(rank)

c Allocate work arrays
      if(n>0)then
      allocate(Wx(n,NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      
      if(Nxro(rank)>0)then
      allocate(WxRO(Nxro(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxlo(rank)>0)then
      allocate(WxLO(Nxlo(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxri(rank)>0)then
      allocate(WxRI(Nxri(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxli(rank)>0)then      
      allocate(WxLI(Nxli(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif

      
c Allocate index arrays
      if(n>0)then
      allocate(index_Wx(n),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxro(rank)>0)then
      allocate(index_WxRO(Nxro(rank)),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxlo(rank)>0)then
      allocate(index_WxLO(Nxlo(rank)),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif

c Pick particles to retain and to send out
      count_Wx = 0; 
      count_WxRO = 0; count_WxLO = 0
      
      if(Npart_old > 0)then
       do i = 1, Npart_old
         if(ix(i)>=1 .and. ix(i)<=Mx)then
            count_Wx = count_Wx + 1
            index_Wx(count_Wx) = i
         elseif(ix(i)<1)then
            count_WxLO = count_WxLO + 1
            index_WxLO(count_WxLO) = i
         elseif(ix(i)>Mx)then
            count_WxRO = count_WxRO + 1
            index_WxRO(count_WxRO) = i
         endif
       enddo
      endif

      if(count_Wx>0) Wx     = part(index_Wx,:)
      if(count_WxRO>0) WxRO = part(index_WxRO,:)
      if(count_WxLO>0) WxLO = part(index_WxLO,:)

      
c Send the left-travelling particles in x- direction
      if(Nxlo(rank) > 0)then
        if (xrank > 0 ) then
         Call MPI_ISEND(WxLO, Npdof*Nxlo(rank), MPI_REAL8, xrank-1,
     .                  xrank, MPI_X_COMM, req(1), ierr)
        else
         WxLO(:,1) = WxLO(:,1) + x_leng
         Call MPI_ISEND(WxLO, Npdof*Nxlo(rank), MPI_REAL8, Px-1,
     .                  xrank, MPI_X_COMM, req(1), ierr)
        endif
      endif

      if(Nxri(rank) > 0)then
        if (xrank < Px-1 ) then
           Call MPI_RECV(WxRI, Npdof*Nxri(rank), MPI_REAL8, xrank+1,
     .                   xrank+1, MPI_X_COMM, status, ierr)
        else
           Call MPI_RECV(WxRI, Npdof*Nxri(rank), MPI_REAL8, 0,
     .                   0, MPI_X_COMM, status, ierr)
        endif
      endif

c Send the right-travelling particles in x- direction
      if(Nxro(rank) > 0)then
       if (xrank < Px-1 ) then
        Call MPI_ISEND(WxRO, Npdof*Nxro(rank), MPI_REAL8, xrank+1,
     .                 xrank+PxMax, MPI_X_COMM, req(2), ierr)
       else
        WxRO(:,1) = WxRO(:,1) - x_leng
        Call MPI_ISEND(WxRO, Npdof*Nxro(rank), MPI_REAL8, 0,
     .                 xrank+PxMax, MPI_X_COMM, req(2), ierr)
       endif
      endif

      if(Nxli(rank) > 0)then
       if (xrank > 0 ) then
           Call MPI_RECV(WxLI, Npdof*Nxli(rank), MPI_REAL8, xrank-1,
     .                   xrank-1+PxMax, MPI_X_COMM, status, ierr)
       else
           Call MPI_RECV(WxLI, Npdof*Nxli(rank), MPI_REAL8, Px-1,
     .                   Px-1+PxMax, MPI_X_COMM, status, ierr)
       endif
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c Allocate new particle arrays
      if(Npart_old > 0) deallocate(part,dpart)

      if(Npart>0)then
        allocate(part(Npart,NpDof),dpart(Npart,3),STAT=ierr)        
        if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
        L = 0
        if(n>0)then
           part(L+1:n,:) = Wx
           L = n
        endif
        if(Nxri(rank)>0)then
           part(L+1:L+Nxri(rank),:) = WxRI ! Right travelling particles
           L = L+Nxri(rank)
        endif
        if(Nxli(rank)>0)then
           part(L+1:L+Nxli(rank),:) = WxLI ! Left travelling particles
           L = L+Nxli(rank)
        endif
      endif

c Dellocate work arrays
      if(allocated(Wx)) deallocate(Wx,STAT=ierr)
      if(allocated(WxRO)) deallocate(WxRO,STAT=ierr)
      if(allocated(WxLO)) deallocate(WxLO,STAT=ierr)
      if(allocated(WxRI)) deallocate(WxRI,STAT=ierr)
      if(allocated(WxLI)) deallocate(WxLI,STAT=ierr)

c-------------------------------------------------------------------------------------
c Done with reordering part. Start reordering partW
c-------------------------------------------------------------------------------------

c Allocate work arrays
      if(n>0)then
      allocate(Wx(n,3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxro(rank)>0)then
      allocate(WxRO(Nxro(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxlo(rank)>0)then
      allocate(WxLO(Nxlo(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxri(rank)>0)then
      allocate(WxRI(Nxri(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nxli(rank)>0)then      
      allocate(WxLI(Nxli(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif

c Reorder Work array for RK4, partW
      if(count_Wx>0) Wx     = partW(index_Wx,:)
      if(count_WxRO>0) WxRO = partW(index_WxRO,:)
      if(count_WxLO>0) WxLO = partW(index_WxLO,:)
      
c Send the left-travelling particles
      if(Nxlo(rank) > 0)then
        if (xrank > 0 ) then
         Call MPI_ISEND(WxLO, 3*Nxlo(rank), MPI_REAL8, xrank-1,
     .                 xrank, MPI_X_COMM, req(1), ierr)
        else
         Call MPI_ISEND(WxLO, 3*Nxlo(rank), MPI_REAL8, Px-1,
     .                 xrank, MPI_X_COMM, req(1), ierr)
        endif
      endif
      if(Nxri(rank) > 0)then
        if (xrank < Px-1 ) then
           Call MPI_RECV(WxRI, 3*Nxri(rank), MPI_REAL8, xrank+1,
     .                   xrank+1, MPI_X_COMM, status, ierr)
        else
           Call MPI_RECV(WxRI, 3*Nxri(rank), MPI_REAL8, 0,
     .                   0, MPI_X_COMM, status, ierr)
        endif
      endif
      
c Send the right-travelling particles
      if(Nxro(rank) > 0)then
       if (xrank < Px-1 ) then
        Call MPI_ISEND(WxRO, 3*Nxro(rank), MPI_REAL8, xrank+1,
     .                 xrank+PxMax, MPI_X_COMM, req(2), ierr)
       else
        Call MPI_ISEND(WxRO, 3*Nxro(rank), MPI_REAL8, 0,
     .                 xrank+PxMax, MPI_X_COMM, req(2), ierr)
       endif
      endif
      
      if(Nxli(rank) > 0)then
       if (xrank > 0 ) then
           Call MPI_RECV(WxLI, 3*Nxli(rank), MPI_REAL8, xrank-1,
     .                   xrank-1+PxMax, MPI_X_COMM, status, ierr)
       else
           Call MPI_RECV(WxLI, 3*Nxli(rank), MPI_REAL8, Px-1,
     .                   Px-1+PxMax, MPI_X_COMM, status, ierr)
       endif
      endif
      
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c Allocate new particle arrays
      if(Npart_old > 0) deallocate(partW)
      
      if(Npart>0)then

      allocate(partW(Npart,3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      
        L = 0
        if(n>0)then
           partW(L+1:n,:) = Wx
           L = n
        endif
        if(Nxri(rank)>0)then
           partW(L+1:L+Nxri(rank),:) = WxRI(1:Nxri(rank),:) ! Right travelling particles
           L = L+Nxri(rank)
        endif
        if(Nxli(rank)>0)then
           partW(L+1:L+Nxli(rank),:) = WxLI(1:Nxli(rank),:) ! Left travelling particles
           L = L+Nxli(rank)
        endif

      endif

c Dellocate work arrays
      if(allocated(Wx)) deallocate(Wx,STAT=ierr)
      if(allocated(WxRO)) deallocate(WxRO,STAT=ierr)
      if(allocated(WxLO)) deallocate(WxLO,STAT=ierr)
      if(allocated(WxRI)) deallocate(WxRI,STAT=ierr)
      if(allocated(WxLI)) deallocate(WxLI,STAT=ierr)
                                                                  
c Deallocate index arrays
      if(allocated(index_Wx)) deallocate(index_Wx)
      if(allocated(index_WxRO)) deallocate(index_WxRO)
      if(allocated(index_WxLO)) deallocate(index_WxLO)
      if(allocated(x)) deallocate(x)
      if(allocated(ix)) deallocate(ix)

!-------------------------------------------------------------------------------------------!
!                                          Y                                                !
!-------------------------------------------------------------------------------------------!                
      if(Npart>0)then
        allocate(iy(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'
        allocate(y(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'     

        y           = part(:,2)
        iy          = FLOOR((y-ycor(1))/DY)+1
        Nylo(rank)  = COUNT(MASK=(iy<1))
        Nyro(rank)  = COUNT(MASK=(iy>My))
        n           = COUNT(MASK=(iy>=1 .and. iy<=My)) 
      else
         n = 0
         Nylo(rank) = 0
         Nyro(rank) = 0
      endif

      if(rank>0)then
        Call MPI_ISEND(Nylo(rank),1, MPI_INTEGER, 0, PxMax+rank,
     .                 MPI_COMM_WORLD, req(1), ierr)
        Call MPI_ISEND(Nyro(rank),1, MPI_INTEGER, 0, 2*PxMax+rank,
     .                 MPI_COMM_WORLD, req(2), ierr)
      endif
      
      if(rank == 0)then
        Do M = Px*Py*Pz-1, 1, -1
          Call MPI_RECV(Nylo(M), 1, MPI_INTEGER, M, PxMax+M,
     .                  MPI_COMM_WORLD, status, ierr)
          Call MPI_RECV(Nyro(M), 1, MPI_INTEGER, M, 2*PxMax+M,
     .                  MPI_COMM_WORLD, status, ierr)
        Enddo
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c Prepare and BCAST Nli, Nri, Nro, Nlo, updated Np      
      if(rank==0)then
        Do IP = 0,Px-1
          Do JP = 0,Py-1
            Do KP = 0,Pz-1
              coords(1) = KP; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank,ierr)
              coords(1) = KP; coords(2) = JP-1; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_left,ierr)
              coords(1) = KP; coords(2) = JP+1; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_right,ierr)
              coords(1) = KP; coords(2) = 0; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_lbdy,ierr)
              coords(1) = KP; coords(2) = Py-1; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_rbdy,ierr)

              if(JP>0)then
                 Nyli(irank) = Nyro(irank_left)
              else
                 Nyli(irank) = Nyro(irank_rbdy)
              endif

              if(JP<Py-1)then
                 Nyri(irank) = Nylo(irank_right)
              else
                 Nyri(irank) = Nylo(irank_lbdy)
              endif

              Np(irank) = Np(irank)
     .                  - Nyro(irank) - Nylo(irank)
     .                  + Nyri(irank) + Nyli(irank)
            Enddo
          Enddo
        Enddo
      endif
      
      Call MPI_BCAST(Np, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)      

      Call MPI_BCAST(Nylo, Px*Py*Pz, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nyro, Px*Py*Pz, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nyli, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nyri, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)     

      NpMax = maxval(Np)
      Npart_old = Npart
      Npart = Np(rank)

c      print *, yrank, Np(rank), Nyli(rank), Nylo(rank),
c     .                Nyri(rank), Nyro(rank)
      
c Allocate work arrays
      if(n>0)then
      allocate(Wy(n,NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      
      if(Nyro(rank)>0)then
      allocate(WyRO(Nyro(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nylo(rank)>0)then
      allocate(WyLO(Nylo(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nyri(rank)>0)then
      allocate(WyRI(Nyri(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nyli(rank)>0)then      
      allocate(WyLI(Nyli(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif            

c Allocate index arrays
      if(n>0)then
      allocate(index_Wy(n),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      
      if(Nyro(rank)>0)then
      allocate(index_WyRO(Nyro(rank)),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nylo(rank)>0)then
      allocate(index_WyLO(Nylo(rank)),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif      

c Pick particles to retain and to send out
      count_Wy   = 0; 
      count_WyRO = 0; count_WyLO = 0
      
      if(Npart_old > 0)then
       do i = 1, Npart_old
         if(iy(i)>=1 .and. iy(i)<=My)then
            count_Wy = count_Wy + 1
            index_Wy(count_Wy) = i      
         elseif(iy(i)<1)then
            count_WyLO = count_WyLO + 1
            index_WyLO(count_WyLO) = i
         elseif(iy(i)>My)then
            count_WyRO = count_WyRO + 1
            index_WyRO(count_WyRO) = i
         endif
       enddo
      endif    
      
      if(count_Wy>0) Wy     = part(index_Wy,:)   
      if(count_WyRO>0) WyRO = part(index_WyRO,:)
      if(count_WyLO>0) WyLO = part(index_WyLO,:)               

c Send the left-travelling particles in y- direction
      if(Nylo(rank) > 0)then
        if (yrank > 0 ) then
         Call MPI_ISEND(WyLO, Npdof*Nylo(rank), MPI_REAL8, yrank-1,
     .                  yrank, MPI_Y_COMM, req(1), ierr)
        else
         WyLO(:,2) = WyLO(:,2) + y_leng*Ny/(Ny-1) 
         Call MPI_ISEND(WyLO, Npdof*Nylo(rank), MPI_REAL8, Py-1,
     .                  yrank, MPI_Y_COMM, req(1), ierr)
        endif
      endif

      if(Nyri(rank) > 0)then
        if (yrank < Py-1 ) then
           Call MPI_RECV(WyRI, Npdof*Nyri(rank), MPI_REAL8, yrank+1,
     .                   yrank+1, MPI_Y_COMM, status, ierr)
        else
           Call MPI_RECV(WyRI, Npdof*Nyri(rank), MPI_REAL8, 0,
     .                   0, MPI_Y_COMM, status, ierr)
        endif
      endif

c Send the right-travelling particles in y- direction
      if(Nyro(rank) > 0)then
       if (yrank < Py-1 ) then
        Call MPI_ISEND(WyRO, Npdof*Nyro(rank), MPI_REAL8, yrank+1,
     .                 yrank+PxMax, MPI_Y_COMM, req(2), ierr)
       else
        WyRO(:,2) = WyRO(:,1) - y_leng*Ny/(Ny-1)
        Call MPI_ISEND(WyRO, Npdof*Nyro(rank), MPI_REAL8, 0,
     .                 yrank+PxMax, MPI_Y_COMM, req(2), ierr)
       endif
      endif

      if(Nyli(rank) > 0)then
       if (yrank > 0 ) then
           Call MPI_RECV(WyLI, Npdof*Nyli(rank), MPI_REAL8, yrank-1,
     .                   yrank-1+PxMax, MPI_Y_COMM, status, ierr)
       else
           Call MPI_RECV(WyLI, Npdof*Nyli(rank), MPI_REAL8, Py-1,
     .                   Py-1+PxMax, MPI_Y_COMM, status, ierr)
       endif
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c Allocate new particle arrays
      if(Npart_old > 0) deallocate(part,dpart)

      if(Npart>0)then
        allocate(part(Npart,NpDof),dpart(Npart,3),STAT=ierr)        
        if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
        L = 0
        if(n>0)then
           part(L+1:n,:) = Wy
           L = n
        endif
        if(Nyri(rank)>0)then
           part(L+1:L+Nyri(rank),:) = WyRI ! Right travelling particles
           L = L+Nyri(rank)
        endif
        if(Nyli(rank)>0)then
           part(L+1:L+Nyli(rank),:) = WyLI ! Left travelling particles
           L = L+Nyli(rank)
        endif
      endif

c Dellocate work arrays
      if(allocated(Wy)) deallocate(Wy,STAT=ierr)
      if(allocated(WyRO)) deallocate(WyRO,STAT=ierr)
      if(allocated(WyLO)) deallocate(WyLO,STAT=ierr)
      if(allocated(WyRI)) deallocate(WyRI,STAT=ierr)
      if(allocated(WyLI)) deallocate(WyLI,STAT=ierr)

c-------------------------------------------------------------------------------------
c Done with reordering part. Start reordering partW
c-------------------------------------------------------------------------------------
c Allocate work arrays
      if(n>0)then
      allocate(Wy(n,3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nyro(rank)>0)then
      allocate(WyRO(Nyro(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nylo(rank)>0)then
      allocate(WyLO(Nylo(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nyri(rank)>0)then
      allocate(WyRI(Nyri(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nyli(rank)>0)then      
      allocate(WyLI(Nyli(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif

c Reorder Work array for RK4, partW
      if(count_Wy>0) Wy     = partW(index_Wy,:)
      if(count_WyRO>0) WyRO = partW(index_WyRO,:)
      if(count_WyLO>0) WyLO = partW(index_WyLO,:)
      
c Send the left-travelling particles
      if(Nylo(rank) > 0)then
        if (yrank > 0 ) then
         Call MPI_ISEND(WyLO, 3*Nylo(rank), MPI_REAL8, yrank-1,
     .                 yrank, MPI_Y_COMM, req(1), ierr)
        else
         Call MPI_ISEND(WyLO, 3*Nylo(rank), MPI_REAL8, Py-1,
     .                 yrank, MPI_Y_COMM, req(1), ierr)
        endif
      endif
      if(Nyri(rank) > 0)then
        if (yrank < Py-1 ) then
           Call MPI_RECV(WyRI, 3*Nyri(rank), MPI_REAL8, yrank+1,
     .                   yrank+1, MPI_Y_COMM, status, ierr)
        else
           Call MPI_RECV(WyRI, 3*Nyri(rank), MPI_REAL8, 0,
     .                   0, MPI_Y_COMM, status, ierr)
        endif
      endif
      
c Send the right-travelling particles
      if(Nyro(rank) > 0)then
       if (yrank < Py-1 ) then
        Call MPI_ISEND(WyRO, 3*Nyro(rank), MPI_REAL8, yrank+1,
     .                 yrank+PxMax, MPI_Y_COMM, req(2), ierr)
       else
        Call MPI_ISEND(WyRO, 3*Nyro(rank), MPI_REAL8, 0,
     .                 yrank+PxMax, MPI_Y_COMM, req(2), ierr)
       endif
      endif
      
      if(Nyli(rank) > 0)then
       if (yrank > 0 ) then
           Call MPI_RECV(WyLI, 3*Nyli(rank), MPI_REAL8, yrank-1,
     .                   yrank-1+PxMax, MPI_Y_COMM, status, ierr)
       else
           Call MPI_RECV(WyLI, 3*Nyli(rank), MPI_REAL8, Py-1,
     .                   Py-1+PxMax, MPI_Y_COMM, status, ierr)
       endif
      endif
      
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c Allocate new particle arrays
      if(Npart_old > 0) deallocate(partW)
      
      if(Npart>0)then

      allocate(partW(Npart,3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      
        L = 0
        if(n>0)then
           partW(L+1:n,:) = Wy
           L = n
        endif
        if(Nyri(rank)>0)then
           partW(L+1:L+Nyri(rank),:) = WyRI(1:Nyri(rank),:) ! Right travelling particles
           L = L+Nxri(xrank)
        endif
        if(Nyli(rank)>0)then
           partW(L+1:L+Nyli(rank),:) = WyLI(1:Nyli(rank),:) ! Left travelling particles
           L = L+Nyli(rank)
        endif

      endif

c Dellocate work arrays
      if(allocated(Wy)) deallocate(Wy,STAT=ierr)
      if(allocated(WyRO)) deallocate(WyRO,STAT=ierr)
      if(allocated(WyLO)) deallocate(WyLO,STAT=ierr)
      if(allocated(WyRI)) deallocate(WyRI,STAT=ierr)
      if(allocated(WyLI)) deallocate(WyLI,STAT=ierr)
                                                               
c Deallocate index arrays
      if(allocated(index_Wy)) deallocate(index_Wy)
      if(allocated(index_WyRO)) deallocate(index_WyRO)
      if(allocated(index_WyLO)) deallocate(index_WyLO)
      if(allocated(y)) deallocate(y)
      if(allocated(iy)) deallocate(iy)

!--------------------------------------------------------------------------------------------!
!                                          Z                                                 !
!--------------------------------------------------------------------------------------------!                 

      if(Npart>0)then
        allocate(iz(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'
        allocate(z(Npart),STAT=err)
        if(err .ne. 0) STOP '*** Particle array allocation failed!'    

         z           = part(:,3)
         iz          = FLOOR((z-zcor(1))/DZ)+1
         Nzlo(rank)  = COUNT(MASK=(iz<1))
         Nzro(rank)  = COUNT(MASK=(iz>Mz))
         n           = COUNT(MASK=(iz>=1 .and. iz<=Mz))
      else
         n = 0
         Nzlo(rank) = 0
         Nzro(rank) = 0
      endif 

      if(rank>0)then
        Call MPI_ISEND(Nzlo(rank),1, MPI_INTEGER, 0, PxMax+rank,
     .                 MPI_COMM_WORLD, req(1), ierr)
        Call MPI_ISEND(Nzro(rank),1, MPI_INTEGER, 0, 2*PxMax+rank,
     .                 MPI_COMM_WORLD, req(2), ierr)
      endif            

      if(rank == 0)then
        Do M = Px*Py*Pz-1, 1, -1
          Call MPI_RECV(Nzlo(M), 1, MPI_INTEGER, M, PxMax+M,
     .                  MPI_COMM_WORLD, status, ierr)
          Call MPI_RECV(Nzro(M), 1, MPI_INTEGER, M, 2*PxMax+M,
     .                  MPI_COMM_WORLD, status, ierr)          
        Enddo
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c Prepare and BCAST Nli, Nri, Nro, Nlo, updated Np                  
      if(rank==0)then
        Do IP = 0,Px-1
          Do JP = 0,Py-1
            Do KP = 0,Pz-1
              coords(1) = KP; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank,ierr)
              coords(1) = KP-1; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_left,ierr)
              coords(1) = KP+1; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_right,ierr)
              coords(1) = 0; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_lbdy,ierr)
              coords(1) = Pz-1; coords(2) = JP; coords(3) = IP
              Call MPI_CART_RANK(MPI_XYZ_COMM, coords, irank_rbdy,ierr)

              if(KP>0)then
                 Nzli(irank) = Nzro(irank_left)
              else
                 Nzli(irank) = Nzro(irank_rbdy)
              endif

              if(KP<Pz-1)then
                 Nzri(irank) = Nzlo(irank_right)
              else
                 Nzri(irank) = Nzlo(irank_lbdy)
              endif

              Np(irank) = Np(irank)
     .                  - Nzro(irank) - Nzlo(irank)
     .                  + Nzri(irank) + Nzli(irank)
            Enddo
          Enddo
        Enddo
      endif

      Call MPI_BCAST(Np, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)      


      Call MPI_BCAST(Nzlo, Px*Py*Pz, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nzro, Px*Py*Pz, MPI_INTEGER, 0, 
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nzli, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)
      Call MPI_BCAST(Nzri, Px*Py*Pz, MPI_INTEGER, 0,
     .               MPI_COMM_WORLD, ierr)           
     
      NpMax = maxval(Np)
      Npart_old = Npart
      Npart = Np(rank)

c Allocate work arrays
      if(n>0)then
      allocate(Wz(n,NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif                       

      if(Nzro(rank)>0)then
      allocate(WzRO(Nzro(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzlo(rank)>0)then
      allocate(WzLO(Nzlo(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzri(rank)>0)then
      allocate(WzRI(Nzri(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzli(rank)>0)then      
      allocate(WzLI(Nzli(rank),NpDof),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif            

c Allocate index arrays
      if(n>0)then
      allocate(index_Wz(n),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
                  
      if(Nzro(rank)>0)then
      allocate(index_WzRO(Nzro(rank)),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzlo(rank)>0)then
      allocate(index_WzLO(Nzlo(rank)),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif            
            
c Pick particles to retain and to send out
      count_Wz = 0;       
      count_WzRO = 0; count_WzLO = 0            
      
      if(Npart_old > 0)then
       do i = 1, Npart_old
         if(iz(i)>=1 .and. iz(i)<=Mz)then
            count_Wz = count_Wz + 1
            index_Wz(count_Wz) = i       
         elseif(iz(i)<1)then
            count_WzLO = count_WzLO + 1
            index_WzLO(count_WzLO) = i
         elseif(iz(i)>Mz)then
            count_WzRO = count_WzRO + 1
            index_WzRO(count_WzRO) = i            
         endif
       enddo
      endif    
                        
      if(count_Wz>0) Wz     = part(index_Wz,:)      
      if(count_WzRO>0) WzRO = part(index_WzRO,:)
      if(count_WzLO>0) WzLO = part(index_WzLO,:)      

c Send the left-travelling particles in z- direction
      if(Nzlo(rank) > 0)then
        if (zrank > 0 ) then
         Call MPI_ISEND(WzLO, Npdof*Nzlo(rank), MPI_REAL8, zrank-1,
     .                  zrank, MPI_Z_COMM, req(1), ierr)
        else
         WzLO(:,1) = WzLO(:,1) + z_leng 
         Call MPI_ISEND(WzLO, Npdof*Nzlo(rank), MPI_REAL8, Pz-1,
     .                  zrank, MPI_Z_COMM, req(1), ierr)
        endif
      endif

      if(Nzri(rank) > 0)then
        if (zrank < Pz-1 ) then
           Call MPI_RECV(WzRI, Npdof*Nzri(rank), MPI_REAL8, zrank+1,
     .                   zrank+1, MPI_Z_COMM, status, ierr)
        else
           Call MPI_RECV(WzRI, Npdof*Nzri(rank), MPI_REAL8, 0,
     .                   0, MPI_Z_COMM, status, ierr)
        endif
      endif

c Send the right-travelling particles in z- direction
      if(Nzro(rank) > 0)then
       if (zrank < Pz-1 ) then
        Call MPI_ISEND(WzRO, Npdof*Nzro(rank), MPI_REAL8, zrank+1,
     .                 zrank+PxMax, MPI_Z_COMM, req(2), ierr)
       else
        WzRO(:,1) = WzRO(:,1) - z_leng
        Call MPI_ISEND(WzRO, Npdof*Nzro(rank), MPI_REAL8, 0,
     .                 zrank+PxMax, MPI_Z_COMM, req(2), ierr)
       endif
      endif

      if(Nzli(rank) > 0)then
       if (zrank > 0 ) then
           Call MPI_RECV(WzLI, Npdof*Nzli(rank), MPI_REAL8, zrank-1,
     .                   zrank-1+PxMax, MPI_Z_COMM, status, ierr)
       else
           Call MPI_RECV(WzLI, Npdof*Nzli(rank), MPI_REAL8, Pz-1,
     .                   Pz-1+PxMax, MPI_Z_COMM, status, ierr)
       endif
      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      
c Allocate new particle arrays
      if(Npart_old > 0) deallocate(part,dpart)

      if(Npart>0)then
        allocate(part(Npart,NpDof),dpart(Npart,3),STAT=ierr)        
        if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
        L = 0
        if(n>0)then
           part(L+1:n,:) = Wz
           L = n
        endif
        if(Nzri(rank)>0)then
           part(L+1:L+Nzri(rank),:) = WzRI ! Right travelling particles
           L = L+Nzri(rank)
        endif
        if(Nzli(rank)>0)then
           part(L+1:L+Nzli(rank),:) = WzLI ! Left travelling particles
           L = L+Nzli(rank)
        endif
      endif
      
c Dellocate work arrays
      if(allocated(Wz)) deallocate(Wz,STAT=ierr)
      if(allocated(WzRO)) deallocate(WzRO,STAT=ierr)
      if(allocated(WzLO)) deallocate(WzLO,STAT=ierr)
      if(allocated(WzRI)) deallocate(WzRI,STAT=ierr)
      if(allocated(WzLI)) deallocate(WzLI,STAT=ierr)

c-------------------------------------------------------------------------------------
c Done with reordering part. Start reordering partW
c-------------------------------------------------------------------------------------

c Allocate work arrays
      if(n>0)then
      allocate(Wz(n,3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzro(rank)>0)then
      allocate(WzRO(Nzro(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzlo(rank)>0)then
      allocate(WzLO(Nzlo(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzri(rank)>0)then
      allocate(WzRI(Nzri(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif
      if(Nzli(rank)>0)then      
      allocate(WzLI(Nzli(rank),3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      endif

c Reorder Work array for RK4, partW
      if(count_Wz>0) Wz     = partW(index_Wz,:)
      if(count_WzRO>0) WzRO = partW(index_WzRO,:)
      if(count_WzLO>0) WzLO = partW(index_WzLO,:)
      
c Send the left-travelling particles
      if(Nzlo(rank) > 0)then
        if (zrank > 0 ) then
         Call MPI_ISEND(WzLO, 3*Nzlo(rank), MPI_REAL8, zrank-1,
     .                 zrank, MPI_Z_COMM, req(1), ierr)
        else
         Call MPI_ISEND(WzLO, 3*Nzlo(rank), MPI_REAL8, Pz-1,
     .                 zrank, MPI_Z_COMM, req(1), ierr)
        endif
      endif
      if(Nzri(rank) > 0)then
        if (zrank < Pz-1 ) then
           Call MPI_RECV(WzRI, 3*Nzri(rank), MPI_REAL8, zrank+1,
     .                   zrank+1, MPI_Z_COMM, status, ierr)
        else
           Call MPI_RECV(WzRI, 3*Nzri(rank), MPI_REAL8, 0,
     .                   0, MPI_Z_COMM, status, ierr)
        endif
      endif

c Send the right-travelling particles
      if(Nzro(rank) > 0)then
       if (zrank < Pz-1 ) then
        Call MPI_ISEND(WzRO, 3*Nzro(rank), MPI_REAL8, zrank+1,
     .                 zrank+PxMax, MPI_Z_COMM, req(2), ierr)
       else
        Call MPI_ISEND(WzRO, 3*Nzro(rank), MPI_REAL8, 0,
     .                 zrank+PxMax, MPI_Z_COMM, req(2), ierr)
       endif
      endif
      
      if(Nzli(rank) > 0)then
       if (zrank > 0 ) then
           Call MPI_RECV(WzLI, 3*Nzli(rank), MPI_REAL8, zrank-1,
     .                   zrank-1+PxMax, MPI_Z_COMM, status, ierr)
       else
           Call MPI_RECV(WzLI, 3*Nzli(rank), MPI_REAL8, Pz-1,
     .                   Pz-1+PxMax, MPI_Z_COMM, status, ierr)
       endif
      endif
      
      Call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

c Allocate new particle arrays
      if(Npart_old > 0) deallocate(partW)
      
      if(Npart>0)then

      allocate(partW(Npart,3),STAT=ierr)
      if(ierr .ne. 0) STOP '*** Particle array allocation failed!'
      
        L = 0
        if(n>0)then
           partW(L+1:n,:) = Wz
           L = n
        endif
        if(Nzri(rank)>0)then
           partW(L+1:L+Nzri(rank),:) = WzRI(1:Nzri(rank),:) ! Right travelling particles
           L = L+Nzri(rank)
        endif
        if(Nzli(rank)>0)then
           partW(L+1:L+Nzli(rank),:) = WzLI(1:Nzli(rank),:) ! Left travelling particles
           L = L+Nzli(rank)
        endif

      endif

c Dellocate work arrays
      if(allocated(Wz)) deallocate(Wz,STAT=ierr)
      if(allocated(WzRO)) deallocate(WzRO,STAT=ierr)
      if(allocated(WzLO)) deallocate(WzLO,STAT=ierr)
      if(allocated(WzRI)) deallocate(WzRI,STAT=ierr)
      if(allocated(WzLI)) deallocate(WzLI,STAT=ierr)
                                                         
c Deallocate index arrays
      if(allocated(index_Wz)) deallocate(index_Wz)
      if(allocated(index_WzRO)) deallocate(index_WzRO)
      if(allocated(index_WzLO)) deallocate(index_WzLO)
      if(allocated(z)) deallocate(z)
      if(allocated(iz)) deallocate(iz)

      return

      end subroutine Reorder_Part
