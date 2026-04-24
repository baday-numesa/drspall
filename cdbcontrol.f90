SUBROUTINE inputCDBdata
!-----------------------------------------------------------------------
! Read input CDB data
!
! By:   David Rudeen
!
! Modification History
!       4/23/03 DKR: Original Version
!-----------------------------------------------------------------------
USE CDBglobals
USE GLOBALS

IMPLICIT None

INTEGER i, ierr, istat, idum, j

IDBIN=10
CALL DBIOPEN(IDBIN,inCDBFilename,idum,idum,IERR)
IF(IERR.NE.0) CALL QAABORT('PROBLEM WITH INPUT CDB')


CALL DBISIZES(IDBin, NQAREC, NINFO, NDIM, NUMNP, NELBLKin, NUMEL, &
              NELX, NELY, NELZ, MAXLNK, NUMATRin, NUMPRPin, NUMNPS, &
              LNPSNL,NUMESS, LGELEL, LGNODL, IERR)

NELBLK = NELBLKin+neb
ALLOCATE (INFO   (NINFO+1),    STAT=istat)
ALLOCATE (QAREC  (4,NQAREC+1), STAT=istat)
ALLOCATE (NUMELB (NELBLK),     STAT=istat)
ALLOCATE (NAMELB (NELBLK),     STAT=istat)
ALLOCATE (IDELB  (NELBLK),     STAT=istat)
ALLOCATE (NUMLNK (NELBLK),     STAT=istat)
! for now NUMPRP max number of properties allowed
! changed to actual number later
NUMPRP = Max(2*NUMPRPin, NUMPRPin+numParam)
maxNUMPRP = NUMPRP
ALLOCATE (PRPname(NUMPRP),       STAT=istat)
ALLOCATE (ISPROK (NELBLK,NUMPRP),STAT=istat)
ALLOCATE (PROP   (NELBLK,NUMPRP),STAT=istat)

CALL DBITITLE(IDBin, HEAD, IERR)

CALL DBIINFO (IDBin, NINFO, INFO, IERR)

CALL DBIQAREC(IDBin, NQAREC, QAREC, IERR)

CALL DBIELBLK(IDBin, NELBLKin ,NAMELB,IDELB,IERR)

CALL DBINELB(IDBin, NELBLKin, NUMELB,IERR)

!zero elements in all blocks
!    - new blocks created later for new grid
DO i = 1, NELBLKin
 NUMELB(i) = 0.0
ENDDO

DO i = 1, NUMPRPin
  CALL DBIPROP(IDBin, PRPname(i),i,ISPROK(1,i),PROP(1,i),idum,idum,IERR)
ENDDO

RETURN
END






SUBROUTINE createCDBHeader
!-----------------------------------------------------------------------
! Generates Output CDB Header records
!
! By:   David Rudeen
!
! Modification History
!       4/23/03 DKR: Original Version
!-----------------------------------------------------------------------
USE Globals
USE CDBglobals

IMPLICIT None

INTEGER i, ierr, istat, idum

! use un-connected output CDB
IDBOUT=11
CALL DBOOPEN(IDBOUT,-1,99,outCDBFilename,IERR)
IF(IERR.NE.0) CALL QAABORT('PROBLEM WITH OUTPUT CDB')


CALL DBOTITLE(IDBout,HEAD,ierrDB)



! if qarec not allocated, none are
if(.NOT.ALLOCATED(QAREC)) then
  ALLOCATE (QAREC(4,NQAREC+1),STAT=istat)
  ALLOCATE (INFO(NINFO),STAT=istat)
endif


!add QA records
CALL QAMAKREC(NQAREC,QAREC)

DO i = 1, NQAREC
  CALL DBOQAREC(IDBout,'ADD',QAREC(1,i),IERR)
  IF(IERR.NE.0) CALL QAABORT('PROBLEM WITH DBOQAREC')
ENDDO


CALL createCDBgrid


!add information records
IF(.NOT.CDBinput)THEN
  INFO(1) = 'RANDOM SEED ...'
  INFO(3) = 'SI - TIME(sec)'
ENDIF
WRITE(INFO(2),'(3I10)')NELX, NELY, NELZ

CALL DBOINFO(IDBout,NINFO,INFO,ierrDB)

! close out header section of CDB
CALL DBOHEAD(IDBOUT,idum,idum,IERR)
IF(IERR.NE.0) CALL QAABORT('PROBLEM WITH DBOHEAD')

RETURN
END



SUBROUTINE createCDBVariables
!-----------------------------------------------------------------------
! Generates CDB history and element variables names and sends them to
! the CDB file
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------
USE Globals
USE CDBGlobals
integer idum

NumHisVar = nhv

HisName(1)  = 'PUMPRS'
HisName(2)  = 'BOTPRS'
HisName(3)  = 'CAVPRS'
HisName(4)  = 'DRILLRAD'
HisName(5)  = 'CAVRAD'
HisName(6)  = 'TENSRAD'
HisName(7)  = 'CUTRAD'
HisName(8)  = 'WBSUPVEL'
HisName(9)  = 'FLUIDVEL'
HisName(10) = 'MUDEJVEL'
HisName(11) = 'WASWELL'
HisName(12) = 'WASEJCT'
HisName(13) = 'CUTMASMX'
HisName(14) = 'GASINJ'
HisName(15) = 'WELLGAS'
HisName(16) = 'GASEJCT'
HisName(17) = 'GASPOSN'
HisName(18) = 'WASPOSN'
HisName(19) = 'CPUTIME'
HisName(20) = 'RUNSTEP'
HisName(21) = 'VOLSTORE'
HisName(22) = 'GASTORE'
HisName(23) = 'WASTORE'
HisName(24) = 'WASINJ'
HisName(25) = 'GASCAV'
HisName(26) = 'SWELLGAS'
HisName(27) = 'SREPOGAS'
HisName(28) = 'GASTOTAL'
HisName(29) = 'GASFROMW'
HisName(30) = 'CUTMASS'
HisName(31) = 'SPLMASS'
HisName(32) = 'TOTMASS'
HisName(33) = 'CUTVOLEQ'
HisName(34) = 'SPLVOLEQ'
HisName(35) = 'TOTVOLEQ'
HisName(36) = 'CUTRUVOL'
HisName(37) = 'CUTRUMAS'
HisName(38) = 'PUMPRATE'
HisName(39) = 'SHEARRAD'
HisName(40) = 'NOZLVEL'
HisName(41) = 'WBUPVEL'
HisName(42) = 'FLUIDTIM'
HisName(43) = 'SWELLWAS'
HisName(44) = 'WASFROMR'
HisName(45) = 'WASTOTAL'
HisName(46) = 'PITGAIN'
HisName(47) = 'MUDEJCT'
HisName(48) = 'SPLVOL2'
HisName(49) = 'SPLMAS2'
HisName(50) = 'BEDDEPTH'
HisName(51) = 'FORCHRAT'


CALL DBOVRNAM(IDBout,'HIS','ADD',NumHisVar,HisName, &
     idum,idum,idum,IERRDB)

EleName(1)  = 'POREPRS'
EleName(2)  = 'RADEFSTR'
EleName(3)  = 'TANEFSTR'
EleName(4)  = 'POREVEL'
EleName(5)  = 'RADELSTR'
EleName(6)  = 'TANELSTR'
EleName(7)  = 'RADSPSTR'
EleName(8)  = 'TANSPSTR'
EleName(9)  = 'FLUDSTRT'
EleName(10) = 'FLUDSTOP'
EleName(11) = 'FAILSTRT'
EleName(12) = 'SUPRVEL'
!EleName(13) = 'FORCHRAT'


EleName(13) = 'WELLPRS'
EleName(14) = 'WELLVEL'
EleName(15) = 'WELLGSMS'
EleName(16) = 'WELLWSMS'
EleName(17) = 'WELLRHO'
EleName(18) = 'WELLWSVF'
EleName(19) = 'WELLGSVF'
EleName(20) = 'WELLSAVF'
EleName(21) = 'WELLWSMF'
EleName(22) = 'WELLGSMF'
EleName(23) = 'WELLMDMF'
EleName(24) = 'WELLVOL'
EleName(25) = 'COORD'



CALL DBOVRNAM(IDBout,'ELE','ADD',NumEleVar,EleName, &
     idum,idum,idum,IERRDB)


ALLOCATE(ISEVOK(NELBLK,numEleVar), STAT=istat)


!setup and write element variable existance table
do i = 1,numEleVar
  do j = 1, NELBLK
    ISEVOK(j,i) = .FALSE.
  enddo
enddo

!repository variables
do i = 1,12
  ISEVOK(reposBlock,i) = .TRUE.
enddo
ISEVOK(reposBlock,25) = .TRUE.

!wellbore Variables
do i = 13,25
 if(firstWellZone == 1)then
   ISEVOK(WellBlock1,i) = .TRUE.
 endif
 ISEVOK(WellBlock2,i) = .TRUE.
enddo


do i = 1, numEleVar
  CALL DBOEVROK(IDBout,eleName(i), i, ISEVOK(1,i), idum,idum, IERRDB)
enddo

RETURN
END


SUBROUTINE createCDBgrid
!-----------------------------------------------------------------------
! Adds repository grid to output CDB assuming 2-D four node elements.
!
! By:   David Rudeen
!
! Modification History
!       12/14/02 DKR: Original Version
!        1/13/03 DKR: added wellbore grid and data
!-----------------------------------------------------------------------

USE CDBGlobals
USE GLOBALS

Implicit None
Integer nodesPerRow, npc, nrn, nwn, nn, nn2, nn3, i, j, block, &
        istat, idmx, ierr, nrz
Real(8) temp
Logical ALLOCATED


NDIM  = 2
NUMNP = 2*(numReposZones+1) + 2*(numWellZones+1)
nodesPerRow = numReposZones+1
NUMEL  = numReposZones + numWellZones
NUMATR = 4
numEleVar = nev
NAMXYZ(1) = 'X'
NAMXYZ(2) = 'Y'
NAMXYZ(3) = 'Z'
NELBLK    = NELBLKin + neb

ALLOCATE (xn    (NUMNP),          STAT=istat)
ALLOCATE (yn    (NUMNP),          STAT=istat)
ALLOCATE (EleVal(NUMEL,numEleVar),STAT=istat)
ALLOCATE (Link  (4,NUMEL),        STAT=istat)
ALLOCATE (EleMap(NUMEL),          STAT=istat)
ALLOCATE (ELEVAT(NUMEL),          STAT=istat)
ALLOCATE (THICK (NUMEL),          STAT=istat)
ALLOCATE (DEL_X (NUMEL),          STAT=istat)
ALLOCATE (DEL_Y (NUMEL),          STAT=istat)
ALLOCATE (ISATOK(NELBLK),        STAT=istat)

! repository
nrn=2*(numReposZones+1)
xn(1) = initialCavityRadius
xn(1+nodesPerRow) = xn(1)
yn(1) = -0.5
yn(1+nodesPerRow) = 0.5

DO i = 2,numReposZones+1
  xn(i) = reposRadiusH(i)
  xn(i+nodesPerRow) = xn(i)
  yn(i) = -0.5
  yn(i+nodesPerRow) = 0.5

END DO

DEL_X(1) = reposDR(1)
DEL_Y(1) = 1.0
THICK(1) = 2.0*Pi*reposRadius(1)
ELEVAT(1)= repositoryTop - 0.5*repositoryThickness
do i = 2, numReposZones
  DEL_X(i) = reposDR(i)
  DEL_Y(i) = 1.0
  THICK(i) = 2.0*Pi*reposRadius(i)
  ELEVAT(i)= ELEVAT(1)
enddo



!wellbore one long string

npc = numWellZones +1
xn(nrn+1) = 0.0
xn(nrn+1+npc) = 0.0
yn(nrn+1) = 1.0
yn(nrn+1+npc) = 2.0

DO i = 2,numWellZones+1
  xn(nrn+i)     = wellPos(i)-0.5*wellzoneSize(i)
  xn(nrn+i+npc) = xn(nrn+i)
  yn(nrn+i)     = 1.0
  yn(nrn+i+npc) = 2.0
END DO

CALL DBOXYZ(IDBout,NDIM,NUMNP,NAMXYZ,XN,YN,ZN,IERRDB)

nrz = numReposZones
DEL_X(nrz+1) = wellZoneSize(1)
DEL_Y(nrz+1) = 1.0
THICK(nrz+1) = wellArea(1)

DO i = 2, numWellZones
  DEL_X(nrz+i) = wellZoneSize(i)
  DEL_Y(nrz+i) = 1.0
  THICK(nrz+i) = wellArea(i)
END DO

ELEVAT(nrz+1)= SurfaceElevation - 0.5*wellZoneSize(1)

Do i = 2, wellBottomIndex-1
  ELEVAT(nrz+i)= ELEVAT(nrz+i-1) - 0.5*(wellZoneSize(i)+wellZoneSize(i-1))
END DO

i = wellBottomIndex
ELEVAT(nrz+i)   = ELEVAT(nrz+i-1)
ELEVAT(nrz+i+1) = ELEVAT(nrz+i)

DO i = wellBottomIndex+2, numWellZones
  ELEVAT(nrz+i) = ELEVAT(nrz+i-1) + 0.5*(wellZoneSize(i)+wellZoneSize(i-1))
ENDDO


NELX      = NUMEL
NELY      = 1
NELZ      = 0
MAXLNK    = 4

if(.NOT.ALLOCATED(IDELB))then
  ALLOCATE(IDELB  (NELBLK), STAT=istat)
  ALLOCATE(NAMELB (NELBLK), STAT=istat)
  ALLOCATE(PRPname(NUMPRP), STAT=istat)
  ALLOCATE(NUMLNK (NELBLK), STAT=istat)
  ALLOCATE(NUMELB (NELBLK), STAT=istat)
  NUMPRP = NUMPRPin + numParam
  maxNUMPRP = NUMPRP
  ALLOCATE(ISPROK(NELBLK,NUMPRP), STAT=istat)
  ALLOCATE(PROP  (NELBLK,NUMPRP), STAT=istat)

endif



do i = 1, NELBLK
  NUMLNK(i) = 4
enddo

if(nelblkin >0)then
  idmx = 0
  do i = 1, NELBLKin
    idmx = max(idmx,IDELB(i))
 enddo
 do  i = NELBLKin+1, NELBLK
    idmx = idmx + 1
   IDELB(i) = idmx
 enddo
else
  do i = 1, NELBLK
   IDELB(i) = i
  enddo
endif


!NAMELB(1) = 'GLOBAL'
!NAMELB(2) = 'DR_SPALL'
!NAMELB(3) = 'REFCON'
!NAMELB(4) = 'BLOWOUT'
!NAMELB(5) = 'BRINESAL'
DRSblock   = NELBLKin+1
reposBlock = NELBLKin+2
wellBlock1 = NELBLKin+3
wellBlock2 = NELBLKin+4
NAMELB(DRSblock)   = 'DATAUSED'
NAMELB(reposBlock) = 'REPOS'
NAMELB(wellBlock1) = 'DOWN_WB'
NAMELB(wellBlock2) = 'UP_WB'

!NUMELB(1) = 0
!NUMELB(2) = 0
!NUMELB(3) = 0
!NUMELB(4) = 0
!NUMELB(5) = 0
NUMELB(DRSblock)   = 0
NUMELB(reposBlock) = numReposZones
NUMELB(wellBlock1) = 2*numWellZones1
NUMELB(wellBlock2) = numWellZones - NUMELB(wellBLock1)

CALL DBOELBLK(IDBout, NELBLK,NAMELB,IDELB,ierrDB)


DO i = 1, numReposZones
  EleMap(i) = i
  LINK(1,i) = i
  LINK(2,i) = i+1
  LINK(3,i) = i+1+nodesPerRow
  LINK(4,i) = i+nodesPerRow
END DO


!along stretched out wellbore
Do i = 1, numWellzones
  j  = i + numReposZones
  eleMap(j) = j
  LINK(1,j) = nrn+i
  LINK(2,j) = nrn+i+1
  LINK(3,j) = nrn+i+1+npc
  LINK(4,j) = nrn+i+npc
END DO


CALL DBOLINK(IDBout,NELBLK,NELX,NELY,NELZ,NUMELB,MAXLNK,NUMLNK, &
             LINK,IERRDB)

CALL DBOMAP(IDBout,NUMEL, eleMap, IERRDB)

DO i = 1, NELBLK
 ISATOK(i) = .FALSE.
END DO
ISATOK(reposBlock) = .TRUE.
ISATOK(wellBlock1) = .TRUE.
ISATOK(wellBLock2) = .TRUE.

CALL DBOATTR(IDBout, 'ADD','ELEVAT',ISATOK, ELEVAT, IERR)
CALL DBOATTR(IDBout, 'ADD','THICK', ISATOK, THICK, IERR)
CALL DBOATTR(IDBout, 'ADD','DEL_X', ISATOK, DEL_X, IERR)
CALL DBOATTR(IDBout, 'ADD','DEL_Y', ISATOK, DEL_Y, IERR)


! Merge properties used with those from input CDB. Orignal values maintained in
! original block. Values used this run are stored in block NELBLKin+1

!Initialize existence table for new blocks and properties
!NUMPRP - 2nd dimension. Will be redefine to actual number of properties in DEFPROP
do i = NELBLKin+1, NELBLK
  do j = 1, NUMPRP
    isprok(i,j) = .FALSE.
  enddo
enddo
do i = 1, NELBLK
  do j = NUMPRPin+1, NUMPRP
    isprok(i,j) = .FALSE.
  enddo
enddo


NUMPRP = NUMPRPin


! define and add properties to CDB - STRING MUST BE 8 CHARACTERS - BLANK FILL
CALL DEFPROP('SURFELEV', surfaceElevation,      DRSblock)
CALL DEFPROP('REPOSTOP', repositoryTop,         DRSblock)
CALL DEFPROP('REPOSTCK', repositoryThickness,   DRSblock)
CALL DEFPROP('DRZTCK  ', dRZThickness,          DRSblock)
CALL DEFPROP('DRZPERM ', dRZPerm,               DRSblock)
CALL DEFPROP('REPOTRAD', repositoryOuterRadius, DRSblock)
CALL DEFPROP('REPIPRES', repositoryInitialPressure,DRSblock)
CALL DEFPROP('FFPORPRS', farFieldPorePressure,  DRSblock)
CALL DEFPROP('FFSTRESS', farfieldStress,        DRSblock)
CALL DEFPROP('REPIPOR ', repositoryInitialPorosity,DRSblock)
CALL DEFPROP('REPIPERM', repositoryInitialPerm, DRSblock)
CALL DEFPROP('FRCHBETA', forchBeta,             DRSblock)
CALL DEFPROP('BIOTBETA', biotBeta,              DRSblock)
CALL DEFPROP('POISSRAT', poissonsRatio,         DRSblock)
CALL DEFPROP('COHESION', cohesion,              DRSblock)
CALL DEFPROP('FRICTANG', frictionAngle,         DRSblock)
CALL DEFPROP('TENSLSTR', tensileStrength,       DRSblock)
CALL DEFPROP('CHARLEN ',  Lt,                   DRSblock)
CALL DEFPROP('PARTDIAM', particleDiameter,      DRSblock)
CALL DEFPROP('GASBSDEN', gasBaseDensity,        DRSblock)
CALL DEFPROP('GASVISCO', gasViscosity,          DRSblock)
CALL DEFPROP('INITMDEN', initialMudDensity,     DRSblock)
CALL DEFPROP('MUDVISCO', mudViscosity,          DRSblock)
CALL DEFPROP('PIPEROUG', wallRoughness(1),      DRSblock)
CALL DEFPROP('ANNUROUG', wallRoughness(2),      DRSblock)
CALL DEFPROP('MUDSOLMX', mudSolidsmax,          DRSblock)
CALL DEFPROP('MUDSOLVE', mudSolidsViscosityExponent,DRSblock)
CALL DEFPROP('BITDIAM ', bitDiameter,           DRSblock)
CALL DEFPROP('PIPEDIAM', pipeDiameter,          DRSblock)
CALL DEFPROP('COLRDIAM', collarDiameter,        DRSblock)
CALL DEFPROP('PIPEID  ', pipeInsideDiameter,    DRSblock)
CALL DEFPROP('COLRLNGT', collarLength,          DRSblock)
CALL DEFPROP('EXITPLEN', exitPipeLength,        DRSblock)
CALL DEFPROP('EXITPDIA', exitPipeDiameter,      DRSblock)
CALL DEFPROP('DRILRATE', drillingRate,          DRSblock)
CALL DEFPROP('BITABOV ', initialBitAboveRepository,DRSblock)
CALL DEFPROP('MUDPRATE', mudPumpRate,           DRSblock)
CALL DEFPROP('MAXPPRES', maxPumpPressure,       DRSblock)
CALL DEFPROP('DDZTHICK', dDZThickness,          DRSblock)
CALL DEFPROP('DDZPERM ', dDZPerm,               DRSblock)
CALL DEFPROP('STPDVOLR', stopDrillingExitVolRate,DRSblock)
CALL DEFPROP('STPPVOLR', stopPumpingExitVolRate,DRSblock)
CALL DEFPROP('STPDTIME', stopDrillingTime,      DRSblock)
CALL DEFPROP('REPODR  ', initialReposZoneSize,  DRSblock)
CALL DEFPROP('WELLDZ  ', initialWellZoneSize,   DRSblock)
CALL DEFPROP('REPODDR ', GrowthRate,            DRSblock)
CALL DEFPROP('WELLDDZ ', wellGrowthRate,        DRSblock)

CALL DEFPROP('GEOMEXP ', dble(geomExponent),    DRSblock)  !apg real(,8)
if(allowFluidization == 'Y')then
  temp = 1.
elseif(allowFluidization == 'A')then
  temp = 2.
else
  temp = 0.
endif
CALL DEFPROP('ALLOWFLD', temp,                  DRSblock)
CALL DEFPROP('WELLSTAB', wellStabilityFactor,   DRSblock)
CALL DEFPROP('REPOSTAB', reposStabilityFactor,  DRSblock)
CALL DEFPROP('MASSDIFF', massDiffusionFactor,   DRSblock)
CALL DEFPROP('MOMDIFF ' , momentumDiffusionFactor,DRSblock)
temp = dble(ValidationTestCase)+0.1*dble(validationSubCase)  !apg real(,4)
CALL DEFPROP('VALIDTC ' ,temp,DRSblock)

CALL DEFPROP('PI      ', Pi,                  DRSblock)
CALL DEFPROP('REFPRES ', AtmosphericPressure, DRSblock)
CALL DEFPROP('GRAVACC ', gravity,             DRSblock)
CALL DEFPROP('RGAS    ', GasConstant,         DRSblock)
CALL DEFPROP('TREPO   ', ReposTemp,           DRSblock)
CALL DEFPROP('H2OCOMP ', WaterCompressibility,DRSblock)
CALL DEFPROP('WASTDENS', WasteDensity,        DRSblock)
CALL DEFPROP('SALTDENS', SaltDensity,         DRSblock)
CALL DEFPROP('SHAPFAC ', ShapeFactor,         DRSblock)
CALL DEFPROP('TENSVEL ', TensileVelocity,     DRSblock)
CALL DEFPROP('BITNZNO ', BitNozzleNumber,     DRSblock)
CALL DEFPROP('BITNZDIA', BitNozzleDiameter,   DRSblock)
CALL DEFPROP('CHOKEFF ', ChokeEfficiency,     DRSblock)
CALL DEFPROP('CAVRAD0 ', initialCavityRadius, DRSblock)
CALL DEFPROP('MINCHVEL', minCharVel         , DRSblock)
temp = minNumLt
CALL DEFPROP('MINNUMLT', temp               , DRSblock)



do i = 1, NUMPRP
 CALL DBOPROP(IDBout,'ADD',PRPname(i),ISPROK(1,i),PROP(1,i),IERR)
enddo

RETURN
END




SUBROUTINE DEFPROP(name, value,block)
!-----------------------------------------------------------------------
! Defines Property Arrays for CDB output. DRSPALL variables must be stored
! in REAL(4) variables for storage on CDB.
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------
Use Globals
Use CDBGlobals
Implicit none
Real(8) value
Integer block, i
Character*8 name


do i = 1, NUMPRPin
  IF( name == PRPname(i))Then
    ! name exist - set exixtence table for block blk
    ISPROK(block,i) = .TRUE.
    PROP  (block,i) = value
    go to 10
  ENDIF
enddo
! add drspall parameter to property list
NUMPRP = NUMPRP + 1

if(NUMPRP > maxNUMPRP)then
 write(diagnosticFileId,*)' ERROR: exceeded property array size'
 write(*,*)' ERROR: exceeded property array size'
 call QAABORT ('exceeded property array size')  !apg was STOP
endif

PRPname(NUMPRP) = name
ISPROK (block,NUMPRP) = .TRUE.
PROP   (block,NUMPRP) = value


10 continue


RETURN
END


SUBROUTINE WriteHisToCDB
!-----------------------------------------------------------------------
! Writes history only record to CDB. DRSPALL variables must be stored
! in REAL(4) variables for storage on CDB.
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------
USE GLOBALS
USE CDBGlobals
integer idum

HisVal(1)  = curPumpPressure
HisVal(2)  = curWellBottomPressure
HisVal(3)  = curCavityPressure
HisVal(4)  = curDrilledRadius
HisVal(5)  = curCavityRadius
HisVal(6)  = curTensileRadius
HisVal(7)  = cuttingsRadiusMax
HisVal(8)  = curWasteBoundaryPoreVelocity
HisVal(9)  = curFluidizationVelocity
HisVal(10) = curMudEjectionVelocity
HisVal(11) = curWasteInWell
HisVal(12) = curWasteEjected
HisVal(13) = cuttingsMassMax
HisVal(14) = curGasInjected
HisVal(15) = curGasInWell
HisVal(16) = curGasEjected
HisVal(17) = curGasPosInWell
HisVal(18) = curWastePosInWell
CALL CPU_TIME(cpuEnd)
cpuTime = cpuEnd/60. - cpuBegin
HisVal(19) = cpuTime
HisVal(20) = RunIndex
HisVal(21) = volStore
HisVal(22) = gasStore
HisVal(23) = wasteStore
HisVal(24) = wasteInjected
hisVal(25) = cavityGasMass
HisVal(26) = sumWellGasMass
hisVal(27) = sumReposGasMass
HisVal(28) = cavityGasMass+sumReposGasMass+gasStore+curGasEjected+sumWellGasMass
HisVal(29) = totalGasFromWaste
cutMass = cutVol*wasteDensity
totMass = TotVol*wasteDensity
splMass = (totVol-cutvol)*wasteDensity
HisVal(30) = cutMass
HisVal(31) = splMass
HisVal(32) = totMass
! assumes uncompacted waste porosity=0.85
cutvoleq = cutVol/(1.0-uncompactedWastePorosity)
splvoleq = (totVol-cutVol)/(1.0-uncompactedWastePorosity)
totvoleq = TotVol/(1.0-uncompactedWastePorosity)
HisVal(33) = cutvoleq
HisVal(34) = splvoleq
HisVal(35) = totvoleq
HisVal(36) = cutTrueVol/(1.0-uncompactedWastePorosity)
HisVal(37) = cutTrueVol*wasteDensity
HisVal(38) = mudPumpRate
HisVal(39) = reposRadiusH(maxShearFailedIndex+1)
HisVal(40) = wellDeltaVInt(wellBottomIndex)
HisVal(41) = wellV        (wellBottomIndex+1)
HisVal(42) = fluidizationTime(firstIntactZone-1)
HisVal(43) = sumWellWasteMass
HisVal(44) = totalWasteFromRepos
HisVal(45) = cavityWasteMass+wasteStore+curWasteEjected+sumWellWasteMass
HisVal(46) = pitGain
HisVal(47) = mudMassEjected
splvol2  = splVol/(1.0-uncompactedWastePorosity)
splmass2 = splvol*wasteDensity
HisVal(48) = splvol2
HisVal(49) = splmass2
beddepth = curTensileRadius-curCavityRadius
HisVal(50) = beddepth
HisVal(51) = maxForchRatio

DO i = 1, NumHisVar
  CALL DBOVAR(IDBout,'HIS',HisName(i),i,HisVal(i),idum,idum,IERRDB)
END DO


RETURN
END



SUBROUTINE WriteEleToCDB
!-----------------------------------------------------------------------
! Writes whole step (element and History variables) to CDB.
! DRSPALL variables must be transfered in REAL(4) variables for storage
! on CDB.
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------

USE GLOBALS
USE CDBGlobals
integer idum

!IF(numReposZones .NE. NUMEL) call QAABORT ('NUMEL problem')  !apg was STOP

! add spatial dependent variables to CDB for current time
if(firstintactzone > 1)then
  do i = 1,firstintactZone-1
    DO j= 3, 8
      EleVal(i,j) = 0.
    ENDDO
      EleVal(i,1)  = reposPres(0)
      if(tensileStrength > 1.01d7)then  !apg e7
        EleVal(i,2)  = 0.0
      else
        EleVal(i,2)  = -TensileStrength
      endif
      EleVal(i,5)  = reposPres(0)
      EleVal(i,9)  = FluidStartTime(i)
      EleVal(i,10) = FluidStopTime (i)
      EleVal(i,11) = FailStartTime (i)
      EleVal(i,12) = 0.0
      EleVal(i,25) = reposRadius(i)
  enddo
endif

do i =firstintactZone,numReposZones
    EleVal(i,1) = reposPres(i)
    EleVal(i,2) = radEffStress(i)
    EleVal(i,3) = tanEffStress(i)
    EleVal(i,4) = poreVelocity(i)
    EleVal(i,5) = radElasticStress(i)
    EleVal(i,6) = tanElasticStress(i)
    EleVal(i,7) = radSeepageStress(i)
    EleVal(i,8) = tanSeepageStress(i)
    EleVal(i,9) = FluidStartTime(i)
    EleVal(i,10) = FluidStopTime(i)
    EleVal(i,11) = FailStartTime(i)
    EleVal(i,12) = superficialVelocity(i)
    EleVal(i,25) = reposRadius(i)
END DO

do i = 1, numWellZones
   j = numreposZones+i
   eleVal(j,13) = wellPres     (i)
   eleVal(j,14) = wellV        (i)
   eleVal(j,15) = wellGasMass  (i)
   eleVal(j,16) = wellWasteMass(i)
   eleVal(j,17) = wellRho(i)
   eleVal(j,18) = wellWasteVol (i)/wellVol(i)
   eleVal(j,19) = wellGasVol   (i)/wellVol(i)
   eleVal(j,20) = wellSaltVol  (i)/wellVol(i)
   eleVal(j,21) = wellWasteMass(i)/(wellRho(i)*wellVol(i))
   eleVal(j,22) = wellGasMass  (i)/(wellRho(i)*wellVol(i))
   eleVal(j,23) = wellMudMass  (i)/(wellRho(i)*wellVol(i))
   eleVal(j,24) = wellVol      (i)
   eleVal(j,25) = wellPos      (i)

enddo

DO i = 1, NumEleVar
  CALL DBOVAR(IDBout,'ELE',EleName(i),i,EleVal(1,i),idum,idum,IERRDB)
END DO


RETURN
END

