subroutine Initialize_MPI(Px,Py,Pz)
!  include 'header'

  include 'mpif.h'
  integer :: ierr
  integer :: rank, xrank, yrank, zrank, mpi_size
  integer :: source_left_x, dest_left_x, source_right_x, dest_right_x
  integer :: source_left_y, dest_left_y, source_right_y, dest_right_y
  integer :: source_left_z, dest_left_z, source_right_z, dest_right_z
  integer :: MPI_XYZ_COMM, MPI_X_COMM, MPI_Y_COMM, MPI_Z_COMM

  Common /MPI/ rank, xrank, yrank, zrank, mpi_size,&
               source_left_x, dest_left_x, source_right_x,&
               dest_right_x, source_left_y, dest_left_y,&
               source_right_y, dest_right_y, source_left_z,&
               dest_left_z, source_right_z, dest_right_z,&
               MPI_XYZ_COMM, MPI_X_COMM, MPI_Y_COMM,&
               MPI_Z_COMM

  integer :: Px, Py, Pz
  !intarr - number of processors in each directions
  !coords - returned coordinates of current processor
  integer, dimension(3) :: intarr, coords
  logical :: reorder
  logical, dimension(3) :: tmplog
  integer :: color, key

!  write(*,*)"MPI Setup"
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
!Check to make sure px,py,pz add up to correct number
  if (mpi_size.ne.(Px*Py*Pz)) stop ' Number or processors numprocs is not equal to Px*Py*Pz'

!Set the dimensions of the cartesian arrangement, made to be same order as breakup order
  intarr(1)=Pz
  intarr(2)=Py
  intarr(3)=Px

!Set periodicity of points
  tmplog(1) = .true.
  tmplog(2) = .true.
  tmplog(3) = .false.
  reorder   = .false.

  call MPI_CART_CREATE(MPI_COMM_WORLD, 3, intarr, tmplog, reorder, MPI_XYZ_COMM,ierr)
  call MPI_CART_GET(MPI_XYZ_COMM, 3, intarr, tmplog, coords, ierr)

  ! x-direction
  call MPI_CART_SHIFT(MPI_XYZ_COMM, 2, -1, source_left_x, dest_left_x, ierr)
  call MPI_CART_SHIFT(MPI_XYZ_COMM, 2, 1, source_right_x, dest_right_x, ierr)

  ! y-direction
  call MPI_CART_SHIFT(MPI_XYZ_COMM, 1, -1, source_left_y, dest_left_y, ierr)
  call MPI_CART_SHIFT(MPI_XYZ_COMM, 1, 1, source_right_y, dest_right_y, ierr)

  ! z-direction
  call MPI_CART_SHIFT(MPI_XYZ_COMM, 0, -1, source_left_z, dest_left_z, ierr)
  call MPI_CART_SHIFT(MPI_XYZ_COMM, 0, 1, source_right_z, dest_right_z, ierr)


  color = coords(1)*Px*Py+coords(2)*Px
  key = coords(3)
  call MPI_COMM_SPLIT(MPI_XYZ_COMM, color, key, MPI_X_COMM, ierr)
  call MPI_COMM_RANK(MPI_X_COMM, xrank, ierr)

  color = coords(1)*Px*Py+coords(3)
  key = coords(2)
  call MPI_COMM_SPLIT(MPI_XYZ_COMM, color, key, MPI_Y_COMM, ierr)
  call MPI_COMM_RANK(MPI_Y_COMM, yrank, ierr)
  
  color = coords(2)*Px+coords(3)
  key = coords(1)
  call MPI_COMM_SPLIT(MPI_XYZ_COMM, color, key, MPI_Z_COMM, ierr)
  call MPI_COMM_RANK(MPI_Z_COMM, zrank, ierr)

  !TODO: Potential add shift operations here

!  if((xrank.EQ.yrank).AND.(xrank.EQ.zrank))then
!   write(*,"(i4 i4 i4 i4 i4 i4 i4)"),rank,zrank,yrank,xrank
!  endif

end subroutine Initialize_MPI
