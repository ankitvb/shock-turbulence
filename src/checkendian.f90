PROGRAM checkEndianness
IMPLICIT NONE

IF(isLittleEndian())THEN
 WRITE(*,FMT="(A)")"is Little Endian"
ELSE
 WRITE(*,FMT="(A)")"is Big Endian"
ENDIF

CONTAINS

FUNCTION isLittleEndian() RESULT(whatzit)
IMPLICIT NONE

LOGICAL :: whatzit

INTEGER(KIND = 1)  :: byteArray(1:2)
INTEGER(KIND = 2) :: shortInteger

shortInteger = 4660
byteArray    = TRANSFER(shortInteger, byteArray)

IF( byteArray(1) == 52 )THEN
 whatzit = .TRUE.
ELSE
 whatzit = .FALSE.
ENDIF
END FUNCTION isLittleEndian

END PROGRAM checkEndianness
