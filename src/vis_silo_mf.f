***********************************************************************
* Subroutine: Write_Visual_Silo                                       *
*                                                                     *
* Function: To write visualization file in the Visit native Silo      *
*           format                                                    *
***********************************************************************
      Subroutine Write_Visual_Silo(time, step)

      include 'mpif.h'
      include 'silo.inc'
      include 'header'
      
      Parameter (num_groups = 8)  ! Number of output files (= # of I/O channels)
      Parameter (nfields = 11)       ! Number of output variables 

      Real*8       Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4)
      Real*8       W(My+1,Mz+1,Mx+1,nfields)
      Real*8       xx(Mx+1), yy(My+1), zz(Mz+1)
      Real*8       time
      Integer      step
      Integer      dbfile, ierr, err
      Integer      optlistid
      Character*128 rank_str, r_str, g_id_str
      Character*128 group_rank_str, index_str(5)
      Character*128 dirname, groupfilename, rootfilename, varname
      Character*128 meshfilenames(Px*Py*Pz), varfilenames(Px*Py*Pz)
      Character*128 varfnames(Px*Py*Pz)
      Character*20  vis_path_C
      Integer      lgroupfilename, lrootfilename
      Integer      lmeshfilenames(Px*Py*Pz), lvarfilenames(Px*Py*Pz)
      Integer      lvarfnames(Px*Py*Pz)
      Integer      meshtypes(Px*Py*Pz), vartypes(Px*Py*Pz)
      Integer      nmesh, nvars, nvarfiles, N

      Common /RKWORK/ Q, Grad

      Real*8 Y(My,Mz,Mx,1)

      Common /SPECIES/ Y

      Real*8 gamma_eff(My,Mz,Mx)

      Common /THERMO/ gamma_eff

      Real*8 beta_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper

      Integer group_size
      Integer group_rank, group_id, g_id

      Integer ldims(3), ndims
 
      group_size = mpi_size/num_groups
      group_rank = MOD(rank, group_size)
      group_id   = FLOOR(REAL(rank/group_size))

      ldims = (/My+1,Mz+1,Mx+1/)

      ndims     = 3
      nmesh     = Px*Py*Pz
      nvars     = nfields 
      nvarfiles = Px*Py*Pz 

c      ierr = dbshowerrors_(DB_ALL)

      Call MakePaddedArrays(W, xx, yy, zz, nfields)

**************************************************************************************
      WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_vis/100))
      WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_vis,100)/10))
      WRITE(index_str(3), '(I128)') MOD(cue_vis,10)

      dirname = TRIM(ADJUSTL(index_str(1))) // 
     .          TRIM(ADJUSTL(index_str(2))) //
     .          TRIM(ADJUSTL(index_str(3)))

      vis_path_C = TRIM(ADJUSTL(vis_path)) // CHAR(0)
      Call Make_Dir(cue_vis, vis_path_C)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

* Writing rank as a five character string
      WRITE(index_str(1), '(I128)') FLOOR(REAL(rank/10000))
      WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(rank,10000)/1000))
      WRITE(index_str(3), '(I128)') FLOOR(REAL(MOD(rank,1000)/100))
      WRITE(index_str(4), '(I128)') FLOOR(REAL(MOD(rank,100)/10))
      WRITE(index_str(5), '(I128)') MOD(rank,10)

      rank_str = TRIM(ADJUSTL(index_str(1))) //  
     .           TRIM(ADJUSTL(index_str(2))) // 
     .           TRIM(ADJUSTL(index_str(3))) //
     .           TRIM(ADJUSTL(index_str(4))) //
     .           TRIM(ADJUSTL(index_str(5)))

*----------------------------------- Local data (mesh + vars) file ---------------------------------------*

* Creating mesh filenames for each group
      WRITE(index_str(1), '(I128)') FLOOR(REAL(group_id/100))
      WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(group_id,100)/10))
      WRITE(index_str(3), '(I128)') MOD(group_id,10)

      groupfilename = TRIM(ADJUSTL(vis_path))     // '/' //
     .                TRIM(ADJUSTL(dirname))      // '/' //
     .                TRIM(ADJUSTL(prefix))       // '.' //
     .                TRIM(ADJUSTL(index_str(1))) //
     .                TRIM(ADJUSTL(index_str(2))) //
     .                TRIM(ADJUSTL(index_str(3)))   

      lgroupfilename = LEN_TRIM(groupfilename)

* Writing mesh and variable data to file           
      Do N = 0, group_size-1
        if(group_rank .EQ. N)then
          if(N .EQ. 0)then                                           ! group root creates group mesh file
           ierr = dbcreate(groupfilename, lgroupfilename, DB_CLOBBER,
     .                     DB_LOCAL, "LMesh", 5, DB_PDB, dbfile)
            if(dbfile .EQ. -1) then
              print *,'Could not create silo file !'
            endif
          else
            ierr = dbopen(groupfilename, lgroupfilename, DB_PDB, 
     .                    DB_APPEND, dbfile)                         ! every other process opens the group mesh file
            if(dbfile .EQ. -1)then
              print *,'Could not open silo file !'
            endif
          endif
         
          err  = dbmkdir(dbfile, rank_str,5, ierr)                   ! each process creates and changes to its
          ierr = dbsetdir(dbfile, rank_str, 5)                       ! own subdirectory within the file

          err  = dbputqm(dbfile, "quadmesh", 8, "yc", 2,             ! each process writes to the group mesh file 
     .           "zc", 2, "xc", 2, yy, zz, xx, ldims,
     .           ndims, DB_DOUBLE, DB_COLLINEAR, DB_F77NULL, ierr) 

          err = dbputqv1(dbfile, "rho", 3, "quadmesh", 8,          ! All processes write to group var file
     .          W(1,1,1,4),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "U", 1, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,1),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "V", 1, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,2),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "W", 1, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,3),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "P", 1, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,5),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "omx", 3, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,6),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "omy", 3, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,7),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "omz", 3, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,8),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "dil", 3, "quadmesh", 8,            ! All processes write to group var file
     .          W(1,1,1,9),  ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "gamma", 5, "quadmesh", 8,
     .          W(1,1,1,10), ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
          err = dbputqv1(dbfile, "Y", 1, "quadmesh", 8,
     .          W(1,1,1,11), ldims, ndims, DB_F77NULL, 0,
     .          DB_DOUBLE, DB_NODECENT, DB_F77NULL, ierr)
 
          ierr = dbclose(dbfile)
        endif

        Call MPI_BARRIER(MPI_COMM_WORLD, ierr)                       ! barrier to synchronize file access
      End Do     

*---------------------------------- Master file -----------------------------------------------*

* Global root creates the master silo file
* Generating master file name
      WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_vis/100))
      WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_vis,100)/10))
      WRITE(index_str(3), '(I128)') MOD(cue_vis,10)
      
      rootfilename = TRIM(ADJUSTL(vis_path)) // '/' //
     .               TRIM(ADJUSTL(prefix))   // '.' //
     .               TRIM(ADJUSTL(index_str(1))) //
     .               TRIM(ADJUSTL(index_str(2))) //
     .               TRIM(ADJUSTL(index_str(3))) // 
     .               '.silo'  

      lrootfilename = LEN_TRIM(rootfilename)

      if(rank .EQ. 0)then
        ierr = dbcreate(rootfilename, lrootfilename, DB_CLOBBER,
     .                  DB_LOCAL,"Master File",11, DB_PDB, dbfile)
        if(dbfile .EQ. -1) then
          print *,'Could not create silo file !'
        endif

* Generating mesh and var filenames
        Do N = 0, mpi_size-1

          g_id = FLOOR(REAL(N/group_size))
          WRITE(index_str(1), '(I128)') FLOOR(REAL(g_id/100))
          WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(g_id,100)/10))
          WRITE(index_str(3), '(I128)') MOD(g_id,10)
          g_id_str = TRIM(ADJUSTL(index_str(1))) //
     .               TRIM(ADJUSTL(index_str(2))) //
     .               TRIM(ADJUSTL(index_str(3))) 

          WRITE(index_str(1), '(I128)') FLOOR(REAL(N/10000))
          WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(N,10000)/1000))
          WRITE(index_str(3), '(I128)') FLOOR(REAL(MOD(N,1000)/100))
          WRITE(index_str(4), '(I128)') FLOOR(REAL(MOD(N,100)/10))
          WRITE(index_str(5), '(I128)') MOD(N,10)

          r_str    = TRIM(ADJUSTL(index_str(1))) //
     .               TRIM(ADJUSTL(index_str(2))) //
     .               TRIM(ADJUSTL(index_str(3))) //
     .               TRIM(ADJUSTL(index_str(4))) //
     .               TRIM(ADJUSTL(index_str(5)))
          
          meshfilenames(N+1)  = TRIM(ADJUSTL(dirname))  // '/' //
     .                          TRIM(ADJUSTL(prefix))   // '.' //
     .                          TRIM(ADJUSTL(g_id_str)) // ':'//
     .                          TRIM(ADJUSTL(r_str))    //'/quadmesh'

          varfilenames(N+1)   = TRIM(ADJUSTL(dirname))  // '/' // 
     .                          TRIM(ADJUSTL(prefix))   // '.' // 
     .                          TRIM(ADJUSTL(g_id_str)) // ':' //
     .                          TRIM(ADJUSTL(r_str))    // '/'

          lmeshfilenames(N+1) = LEN_TRIM(meshfilenames(N+1))
          lvarfilenames(N+1)  = LEN_TRIM(varfilenames(N+1))

               meshtypes(N+1) = DB_QUAD_RECT
               vartypes(N+1)  = DB_QUADVAR
        End Do

* Saving step and simulation time as options
        ierr = dbmkoptlist(2, optlistid)
        ierr = dbaddiopt(optlistid, DBOPT_CYCLE, step)
        ierr = dbadddopt(optlistid, DBOPT_DTIME, time)

* Set max string length
        oldlen = dbget2dstrlen()
        err    = dbset2dstrlen(128)

* Write the multimesh object
        err   = dbputmmesh(dbfile, "quadmesh", 8, nmesh, meshfilenames,
     .                     lmeshfilenames, meshtypes, optlistid, ierr)

* Write the multivar objects
        Do M = 1, nvars
          Call GetVarName(M, varname)
          varname  = TRIM(ADJUSTL(varname))
          lvarname = LEN_TRIM(varname)

          Do N = 0,mpi_size-1
            varfnames(N+1) = TRIM(ADJUSTL(varfilenames(N+1))) //
     .                       TRIM(ADJUSTL(varname))
            lvarfnames(N+1) = lvarfilenames(N+1) + lvarname
          End Do

          err = dbputmvar(dbfile, varname, lvarname, nvarfiles,
     .                      varfnames, lvarfnames, vartypes,
     .                      DB_F77NULL, ierr)
        End Do

* Cleaning up
        err = dbset2dstrlen(oldlen)
        err = dbfreeoptlist(optlistid)

        ierr = dbclose(dbfile)
      endif

      End Subroutine Write_Visual_Silo

****************************************************************************
* Subroutine: GetVarName                                                   *
*                                                                          *
* Function: Return variable name corresponding to its number               *
*                                                                          *
****************************************************************************
      Subroutine GetVarName(M, varname)

      Integer M
      Character*128 varname
  
      Select Case (M)      
        Case(1)
          varname = "rho"
        Case(2)
          varname = "U"
        Case(3)
          varname = "V"
        Case(4)
          varname = "W"
        Case(5)
          varname = "P"
        Case(6)
          varname = "omx"
        Case(7)
          varname = "omy"
        Case(8)
          varname = "omz"
        Case(9)
          varname = "dil"
        Case(10)
          varname = "gamma"
        Case(11)
          varname = "Y"

      End Select
        
      Return

      End Subroutine GetVarName

******************************************************************************
* Subroutine: Write_Part                                                     *
*                                                                            *
* Function: Write particle data to file                                      *
*                                                                            *
******************************************************************************
      Subroutine Write_Part(time, step)

      Include 'mpif.h'
      Include 'silo.inc'
      Include 'header'
      Include 'part_header'

      Integer step, count, M, I
      Real*8 time

      Character*128 index_str(5)
      Character*128 partfilename
      Integer lpartfilename
      Integer dbfile, optlistid

      Real*8, dimension(:,:), pointer     :: partall
      Real*8, dimension(:,:), pointer     :: recvbuf
      Integer Npartall, ndims

      Integer status(MPI_STATUS_SIZE), ierr, err, req(3)

      ndims = 3

      WRITE(index_str(1), '(I128)') FLOOR(REAL(cue_vis/100))
      WRITE(index_str(2), '(I128)') FLOOR(REAL(MOD(cue_vis,100)/10))
      WRITE(index_str(3), '(I128)') MOD(cue_vis,10)

      partfilename = TRIM(ADJUSTL(vis_path)) // '/' //
     .               'part.' //
     .               TRIM(ADJUSTL(index_str(1))) //
     .               TRIM(ADJUSTL(index_str(2))) //
     .               TRIM(ADJUSTL(index_str(3))) //
     .               '.silo'

      lpartfilename = LEN_TRIM(partfilename)

* All processes send particle data
      if(Npart>0)then
        Call MPI_ISEND(part, Np(rank)*NpDof, MPI_REAL8, 0, rank,
     .                 MPI_COMM_WORLD, req(1), ierr)
      endif

* Process # 0 recieves particle data in order
      if(rank .EQ. 0)then
        Npartall = sum(Np)
        allocate(partall(Npartall,NpDof),STAT=ierr)
        if(ierr .ne. 0) Stop '*** Particle array allocation failed!'

        count = 1

        Do M = 0, Px*Py*Pz-1
          if(Np(M)>0)then
            allocate(recvbuf(Np(M),NpDof),STAT=ierr)
            if(ierr .ne. 0) Stop '*** Particle array allocation failed!'

            Call MPI_RECV(recvbuf, Np(M)*NpDof, MPI_REAL8, M,
     .                    M, MPI_COMM_WORLD, status, ierr)
            partall(count:count+Np(M)-1,:) = recvbuf(1:Np(M),:)
            count = count + Np(M)

            deallocate(recvbuf, STAT=ierr)
            if(ierr .ne. 0) Stop '*** recvbuf deallocation failed!'
          endif
        End Do

      endif

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      if(rank .EQ. 0)then
        ierr = dbcreate(partfilename, lpartfilename, DB_CLOBBER,
     .                  DB_LOCAL,"Particle Data",13, DB_PDB, dbfile)
        if(dbfile .EQ. -1) then
          print *,'Could not create silo file !'
        endif

        ierr = dbmkoptlist(2, optlistid)
        ierr = dbaddiopt(optlistid, DBOPT_CYCLE, step)
        ierr = dbadddopt(optlistid, DBOPT_DTIME, time)

        ierr = dbputpm(dbfile, "pointmesh", 9, ndims, partall(1,1),
     .                 partall(1,2), partall(1,3), Npartall, DB_DOUBLE,
     .                 optlistid, err)

        ierr = dbputpv1(dbfile, "rho", 3, "pointmesh", 9,
     .               partall(1,4), Npartall, DB_DOUBLE, DB_F77NULL, err)
        ierr = dbputpv1(dbfile, "P", 1, "pointmesh", 9,
     .               partall(1,5), Npartall, DB_DOUBLE, DB_F77NULL, err)
        ierr = dbputpv1(dbfile, "om", 2, "pointmesh", 9,
     .               partall(1,6), Npartall, DB_DOUBLE, DB_F77NULL, err)
        ierr = dbputpv1(dbfile, "dil", 3, "pointmesh", 9,
     .               partall(1,7), Npartall, DB_DOUBLE, DB_F77NULL, err)

        err = dbfreeoptlist(optlistid)

        ierr = dbclose(dbfile)

        deallocate(partall,STAT=ierr)
        if(ierr .ne. 0) Stop '*** Particle array deallocation failed!'

      endif

      End Subroutine Write_Part

********************************************************************************
* Subroutine: MakePaddedArrays                                                 *
*                                                                              *
* Function: To generate local arrays padded by one extra node in all           *
*           directions. This node contains redundant data from neighbouring    *
*           processors. These ghost zones are neeeded by VisIt to plot         *
*           seamlessly without domain decomposition artifacts                  *
*           For boundary nodes the data is simply copied from the penultimate  *
*           grid points                                                        *
********************************************************************************
      Subroutine MakePaddedArrays(W, xx, yy, zz, nfields)

      include 'header'
      include 'mpif.h'

      Real*8 Q(My,Mz,Mx,5), Grad(My,Mz,Mx,4), Vs(My,Mz,Mx)
      Real*8 W(My+1,Mz+1,Mx+1,nfields), Q2(My,Mz,Mx,nfields)
      Real*8 xx(Mx+1), yy(My+1), zz(Mz+1)
      Real*8 fac1

      Common /RKWORK/ Q, Grad

      Real*8 Y(My,Mz,Mx,1)

      Common /SPECIES/ Y

      Real*8 beta_hyper(My,Mz,Mx)

      Common /HYPER/ beta_hyper

      Real*8 gamma_eff(My,Mz,Mx)

      Common /THERMO/ gamma_eff


      Real*8 Wx(My,Mz,nfields), Wy(Mz,Mx,nfields), Wz(My,Mx,nfields)
      Real*8 Wxy(Mz,nfields), Wyz(Mx,nfields), Wxz(My,nfields)
      Real*8 Wxyz(nfields)
   
      Real*8 Qxbuf(My,Mz,nfields), Qybuf(Mz,Mx,nfields), 
     .       Qzbuf(My,Mx,nfields)
      Real*8 Qxybuf(Mz,nfields), Qyzbuf(Mx,nfields), Qxzbuf(My,nfields)
      Real*8 Qxyzbuf(nfields)

      Integer nfields
      Integer I, J, K, L

      Integer x_source, y_source, z_source
      Integer x_dest, y_dest, z_dest
      Integer xy_source, yz_source, xz_source, xyz_source
      Integer xy_dest, yz_dest, xz_dest, xyz_dest
      Integer status_arr(MPI_STATUS_SIZE,3), status(MPI_STATUS_SIZE)
      Integer coords(3), ierr, req(4)

* Padding grid arrays 
      xx(1:Mx)              = xcor(1:Mx)
      yy(1:My)              = ycor(1:My)
      zz(1:Mz)              = zcor(1:Mz)
      xx(Mx+1)              = xx(Mx) + DX
      yy(My+1)              = yy(My) + DY
      zz(Mz+1)              = zz(Mz) + DZ

* Converting to primitive variables
      Do I = 1,Mx
        Do J = 1,My
          Do K = 1,Mz
            fac1 = gamma_eff(J,K,I) - 1.0D0
            Vs(J,K,I)   = 1D0/Q(J,K,I,4)
            Q2(J,K,I,1) = Q(J,K,I,1)/Q(J,K,I,4)
            Q2(J,K,I,2) = Q(J,K,I,2)/Q(J,K,I,4)
            Q2(J,K,I,3) = Q(J,K,I,3)/Q(J,K,I,4)
            Q2(J,K,I,4) = Q(J,K,I,4)
            Q2(J,K,I,5) = fac1 * ( Q(J,K,I,5) - 5.0D-1 * Vs(J,K,I) *
     .                           ( Q(J,K,I,1) * Q(J,K,I,1)
     .                           + Q(J,K,I,2) * Q(J,K,I,2)
     .                           + Q(J,K,I,3) * Q(J,K,I,3) ) )
            Q2(J,K,I,6:9) = Grad(J,K,I,1:4)
            Q2(J,K,I,10)  = gamma_eff(J,K,I) !beta_hyper(J,K,I)
            Q2(J,K,I,11)  = Y(J,K,I,1)
          End Do
        End Do
      End Do 


* Padding data arrays
      do L = 1,nfields
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
        Call MPI_ISEND(Wx, My*Mz*nfields, MPI_REAL8, x_dest,
     .                 xrank, MPI_X_COMM, req(1), ierr)
        Call MPI_ISEND(Wy, Mx*Mz*nfields, MPI_REAL8, y_dest,
     .                 yrank, MPI_Y_COMM, req(2), ierr)
        Call MPI_ISEND(Wz, My*Mx*nfields, MPI_REAL8, z_dest,
     .                 zrank, MPI_Z_COMM, req(3), ierr)

c All processors receive it in buffer array
        Call MPI_RECV(Qxbuf, My*Mz*nfields, MPI_REAL8, x_source,
     .                x_source, MPI_X_COMM, status, ierr)
        Call MPI_RECV(Qybuf, Mx*Mz*nfields, MPI_REAL8, y_source,
     .                y_source, MPI_Y_COMM, status, ierr)
        Call MPI_RECV(Qzbuf, My*Mx*nfields, MPI_REAL8, z_source,
     .                z_source, MPI_Z_COMM, status, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      Call MPI_ISEND(Wxy, Mz*nfields, MPI_REAL8, xy_dest,
     .               rank, MPI_COMM_WORLD, req(1), ierr)
      Call MPI_ISEND(Wyz, Mx*nfields, MPI_REAL8, yz_dest,
     .               rank, MPI_COMM_WORLD, req(2), ierr)
      Call MPI_ISEND(Wxz, My*nfields, MPI_REAL8, xz_dest,
     .               rank, MPI_COMM_WORLD, req(3), ierr)
      Call MPI_ISEND(Wxyz, nfields, MPI_REAL8, xyz_dest,
     .               rank, MPI_COMM_WORLD, req(4), ierr)

* All processors recieve edge and corner data
      Call MPI_RECV(Qxybuf, Mz*nfields, MPI_REAL8, xy_source,
     .              xy_source, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(Qyzbuf, Mx*nfields, MPI_REAL8, yz_source,
     .              yz_source, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(Qxzbuf, My*nfields, MPI_REAL8, xz_source,
     .              xz_source, MPI_COMM_WORLD, status, ierr)
      Call MPI_RECV(Qxyzbuf, nfields, MPI_REAL8, xyz_source,
     .              xyz_source, MPI_COMM_WORLD, status, ierr)

      Call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      W(1:My,1:Mz,1:Mx,1:nfields) = Q2
 
      do J = 1, My
        do K = 1, Mz
          if(xrank < Px-1)then
            W(J,K,Mx+1,1:nfields) = Qxbuf(J,K,1:nfields)
          else
            W(J,K,Mx+1,1:nfields) = W(J,K,Mx,1:nfields)
          endif
        enddo
      enddo
      do K = 1, Mz
        do I = 1, Mx
          if(yrank < Py-1)then
            W(My+1,K,I,1:nfields) = Qybuf(K,I,1:nfields)
          else
            W(My+1,K,I,1:nfields) = W(My,K,I,1:nfields)
          endif
        enddo
      enddo
      do J = 1, My
        do I = 1, Mx
          if(zrank < Pz-1)then
            W(J,Mz+1,I,1:nfields) = Qzbuf(J,I,1:nfields)
          else
            W(J,Mz+1,I,1:nfields) = W(J,Mz,I,1:nfields)
          endif
        enddo
      enddo

      do K = 1,Mz
        if((xrank < Px-1).AND.(yrank < Py-1))then
          W(My+1,K,Mx+1,1:nfields) = Qxybuf(K,1:nfields)
        else
          W(My+1,K,Mx+1,1:nfields) = W(My,K,Mx,1:nfields) 
        endif
      enddo

      do I = 1,Mx
        if((yrank < Py-1).AND.(zrank < Pz-1))then
          W(My+1,Mz+1,I,1:nfields) = Qyzbuf(I,1:nfields)
        else
          W(My+1,Mz+1,I,1:nfields) = W(My,Mz,I,1:nfields)
        endif
      enddo

      do J = 1,My
        if((zrank < Pz-1).AND.(xrank < Px-1))then
          W(J,Mz+1,Mx+1,1:nfields) = Qxzbuf(J,1:nfields)
        else
          W(J,Mz+1,Mx+1,1:nfields) = W(J,Mz,Mx,1:nfields)
        endif
      enddo

      if((xrank < Px-1).AND.(yrank < Py-1).AND.(zrank < Pz-1))then
        W(My+1,Mz+1,Mx+1,1:nfields) = Qxyzbuf(1:nfields)
      else
        W(My+1,Mz+1,Mx+1,1:nfields) = W(My,Mz,Mx,1:nfields)
      endif
      
      Return

      End Subroutine MakePaddedArrays 
