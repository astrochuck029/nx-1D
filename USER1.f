      SUBROUTINE USER1(GG,T,C,Q,QD,R,TIME,DT,IT,KODE,NOCON,MAXNO,
     + ICONV,DUM1,DUM2,DTP,TF)

      IMPLICIT NONE

C======================================================================
C  USER1 ARGUMENTS
C======================================================================

      DOUBLE PRECISION T(*),TIME,DT

      REAL GG(*),C(*),Q(*),QD(*),R(*)
      REAL DUM1,DUM2,DTP,TF

      INTEGER IT,KODE,NOCON,MAXNO,ICONV(*)

C======================================================================
C  TMG COMMON BLOCKS
C======================================================================

      REAL TDMAX,PRTFLG,PARAMS(80000)
      REAL GRAV,GV(3),TABS,RGAS
      REAL PSTD,TSTD,SIGMA

      INTEGER IRUN,IR(1)
      INTEGER MAXN1,MAXN2

      COMMON/TDMAX/TDMAX
      COMMON/PRTFLG/PRTFLG
      COMMON/IRUN/IRUN,IR
      COMMON/MAXNOQ/MAXN1,MAXN2
      COMMON/PARAMS/PARAMS
      COMMON/GRAV/GRAV,GV,TABS,RGAS,PSTD,TSTD,SIGMA

      SAVE

C======================================================================
C  GROUP VARIABLES
C======================================================================

      INTEGER GID
      INTEGER GLENGTH

      CHARACTER*7  SNAME
      CHARACTER*80 LNAME

      LOGICAL FIRST

      DATA FIRST/.TRUE./

C======================================================================
C  INITIALIZATION
C======================================================================

      IF (FIRST) THEN

         PRINT *,' '
         PRINT *,'=========================================='
         PRINT *,' CFRP USER1 STAGE 3A'
         PRINT *,' LONGNAME ONLY'
         PRINT *,'=========================================='

C        External NX group name
         LNAME = 'Plate_3D'

         PRINT *,' '
         PRINT *,'CALLING LONGNAME'
         PRINT *,'EXTERNAL GROUP = ',LNAME

C        Map external group name -> internal TMG name
         CALL LONGNAME(SNAME,LNAME,GID,GLENGTH,2)

         PRINT *,' '
         PRINT *,'LONGNAME RETURNED'
         PRINT *,'INTERNAL NAME   = ',SNAME
         PRINT *,'GROUP ID        = ',GID
         PRINT *,'NAME LENGTH     = ',GLENGTH
         PRINT *,' '

         PRINT *,'=========================================='
         PRINT *,' LONGNAME TEST COMPLETE'
         PRINT *,'=========================================='
         PRINT *,' '

         FIRST=.FALSE.

      ENDIF

      RETURN
      END