************************************************************************
*  Subroutine: Open_File.f                                             *
*                                                                      *
*  Function:   It is a set of routines to open Fortran files using an  *
*              optional version number.                                *
************************************************************************
      Subroutine Open_File(u, path, name, f, s, verno, rewind)

      Character     f, s
      Character*(*) path, name
      Integer       u, verno, rewind

*  Locals
      Character*40 filename, ff, ss
      Integer      I

*  External routines
      Character*40 MakeFileName
      Integer      nextVersion, maxVersion

*  Get the filename
      If (verno .EQ. 0) Then
        I = nextVersion(path, name)

      Else If (verno .EQ. -1) Then
        I = maxVersion(path, name)

      Else
        I = verno

      End If

      filename = MakeFileName(path, name, I)

*  Determine file type
      If (f .EQ. 'F') ff = 'formatted'
      If (f .EQ. 'U') ff = 'unformatted'

      If (s .EQ. 'U') ss = 'unknown'
      If (s .EQ. 'O') ss = 'old'
      If (s .EQ. 'N') ss = 'new'

*  Open file
      If (rewind .EQ. 0) Then
        Open(unit=u, file = filename, status = ss, form = ff, 
     .       position = 'append')

      Else
        Open(unit=u, file = filename, status = ss, form = ff)

      End If

      Return

      End
