MODULE CDBGlobals
!-----------------------------------------------------------------------
! This module conatins variables required to read or write a CDB
!
! By:   David Rudeen
!       GRAM, Inc.
!
! Modification History
! 12/19/02 - DKR: Original version
!-----------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER :: nev=25, nhv=51, nep = 71
CHARACTER*8   HisName(nhv), NAMXYZ(3), EleName(nev), blockName
CHARACTER*80  HEAD
CHARACTER*255 outCDBFilename, inCDBFilename
INTEGER       IDBout, NOUT, NumHisVar, LDUM, IERRDB, NDIM, NUMNP, &
              NELBLK, NELX, NELY, NELZ, MAXLNK, NEB, NUMPRP, NUMATR, &
	          NUMEL, NumEleVar, outCDBFileID, inCDBFileID, &
	          NINFO, numCDBsteps, numParam, maxNUMPRP
INTEGER       IDBin, NQAREC, NELBLKin, NUMATRin, NUMPRPin, NUMNPS, LNPSNL, &
              NUMESS, LGELEL, LGNODL, DRSblock, reposBlock, wellBLock1, wellBlock2
LOGICAL,     DIMENSION(:,:), ALLOCATABLE :: ISEVOK, ISPROK
LOGICAL,     DIMENSION(:),   ALLOCATABLE :: ISATOK
INTEGER,     DIMENSION(:,:), ALLOCATABLE :: Link
INTEGER,     DIMENSION(:),   ALLOCATABLE :: EleMap, IDELB, NUMLNK, NUMELB
REAL(8),     DIMENSION(:),   ALLOCATABLE :: Xn, Yn, ELEVAT, THICK, DEL_Y, DEL_X
REAL(8),     DIMENSION(:,:), ALLOCATABLE :: EleVal, PROP
CHARACTER*80,DIMENSION(:),   ALLOCATABLE :: INFO
CHARACTER*8 ,DIMENSION(:),   ALLOCATABLE :: NAMELB, PRPname, ParamName, NAMATR
CHARACTER*8 ,DIMENSION(:,:), ALLOCATABLE :: QAREC

REAL        RMEM(1)  !apg must be REAL
REAL(8)     HisVal(nhv), ZN(1)
CHARACTER*1 CMEM(1)
LOGICAL     WHOLE, CDBinput

END



