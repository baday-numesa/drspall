PROGRAM DRSPALL
!------------------------------------------------------------------------------
! Original Preface:
!
! Program 'DR_Spall' is written as part of the work related to
! Analysis Plan AP-048 entitled:
!
! ANALYSIS PLAN FOR DEVELOPMENT OF IMPROVED SPALL FAILURE PHYSICAL
! MODEL AND COMPUTATIONAL METHODOLOGY
!
! The work is performed by John F. Schatz Research & Consulting, Inc.
! under consulting contract AX-9033 to Sandia National Laboratories.

! The purpose of the Analysis Plan is to create and test an improved
! computational tool for calculating the volume of waste subject to
! material failure during an inadvertent drilling intrusion of WIPP.
! To do so, the following initial objectives are set:
! (a) convert the existing 'GasOut' code to Fortran to allow operation
! on either Intel-based personal computers (PCs) or on DEC Alpha
! hardware (under Windows NT) and be more compatible with existing
! Sandia National Laboratories (SNL) performance assessment (PA)
! software systems,
! (b) investigate the possible replacement of the existing GasOut
! finite difference scheme with a second-order accurate fully explicit
! method, for the reduction of numerical errors in areas of high
! gradients,
! (c) add an improved wellbore model, and
! (d) improve the waste tensile failure model to allow smoother and
! more realistic stress relaxation.  Finally, allowance is made for
! modification and testing of the program as required.
!
! The program is written using Digital Visual Fortran Professional
! Edition and Microsoft Developer Studio.  Within Visual Fortran, the
! 'QuickWin' option is used.  QuickWin is a Digital-enabled subset of
! the full Windows API that supports pixel-based graphics, real-
! coordinate graphics, text windows, character fonts, user-defined
! menus, mouse events, and editing (select/copy/paste) of text,
! graphics, or both.  Using QuickWin avoids some of the complexities of
! the full Windows API that are not necessary for this project.
!
! Addendum:
!
! This particular implementation for the Alpha OpenVMS platform by
! David Rudeen (GRAM, Inc) is the WIPP Qualified version of the software.
! An interface to the WIPP CAMDAT Database (CDB) has been implemented
! for I/O and dynamic memory has been added. All Windows based graphics
! and I/O have been removed.
!
!apg Version 1.21  05/xx/13  *** Migrate from Alpha OpenVMS to Solaris ***
!apg Modified by Amy Gilkey and Dave Rudeen (see !apg or !dkr comments)
!apg    All real variables are now double precision.
!apg    All error STOPs were changed to calls to QAABORT, 
!apg     which returns an error status to the system.
!apg    Minor code cleanup, such as allowing Unix file names
!apg     or using "trim" to print strings.
!apg    Minor cosmetic changes may be visible in the output.
!apg    The validation file names are now prefixed by the output camdat
!apg     file name rather than fixed names; the rest of name is unchanged.
!apg    The header for the validation files is shorter.  It contains the QAPAGE
!apg     line and the file name (new) only.
!    Version 1.22  08/2014  Modified by Todd Zeitler (!trz) and Dwayne Kicker
!       Correct an error in the implementation of the finite difference equations
!       in accordance with SPR 13-001. Change three Forchterm equations
!       (first cell, interior cell, last cell). Correct the derivation of the
!       constant zone size equations. Remove the variable zone size component
!       of DRSPALL and run the code exclusively with a constant zone size.
!       See Change Control Form for detailed changes.
!
!-----------------------------------------------------------------------------------------

Use Globals
Use CDBGLobals
Implicit NONE

Integer(2)  ierr
Integer(2) j !trz

CALL CPU_TIME(cpuBegin)
cpubegin = cpuBegin/60.0

!Platform is software  version (VMS or PC)
!machine is hardware  used for execution (VMS or PC)
platform = 'VMS'
machine  = 'VMS'


FIX = 1
ENG = 2


!Initialize Step control location variables
stepControl = 1
stepControlName(1)='Trial'
stepControlName(2)='Repos'
stepControlName(3)='Tensil'
stepControlName(4)='Wellbor'
stepControlName(5)='MaxDelt'


! For each file to be saved/created, initialize a file identifier
parameterFileID    = 100
diagnosticFileID   = 108
outCDBFileID       = 109
inCDBFileID        = 110

! test Case #1
chanValidationFileID = 111

! test case #4,
couplingValidationFileID        = 120
stressValidationFileID          = 121
fluidizationValidationFileID    = 122
expulsionValidationFileID       = 123
fluidizationTimeValidationFileID = 124 !trz

! test case #5,
wellboreValidationFileID		= 125

summaryFileID = diagnosticFileID


! Input and initialization control
CALL ALPINIT


!run control - time loop
Call RunLoop


WRITE(diagnosticFileID, '(//,'' Final Summary'',/)')
Call WriteSummaryFile

CALL CPU_TIME(cpuEnd)
cpuTime = cpuEnd/60. - cpuBegin
write(6,*)' cpu_time=',cpuTime

write(diagnosticFileID,*)' cpu_time=',cpuTime

Call CloseRunFiles

CALL DBOCLOSE(IDBout,IERR)
IF(IERR.NE.0) CALL QAABORT('PROBLEM WITH DBOCLOSE')

STOP 'DRSPALL Normal Completion'
END




      SUBROUTINE ALPINIT
!-----------------------------------------------------------------------
!
!    ADAPTED FROM     *** PANEL v4.0 Feb 1998 ***
!   --
!   --External software used:
!   --   CAMDAT_LIB package (read/write CAMDAT database)
!   --   CAMCON_LIB package (QA routines, etc.)
!   --   SUPES package (dynamic memory, FORTRAN extensions)
!   --
!   --Runs on DEC Alpha Open VMS AXP Version 7.3
!
!-----------------------------------------------------------------------

USE GLOBALS
USE CDBGlobals

CHARACTER*255 FILES(11)
CHARACTER*132 DEFNAM(8)
CHARACTER*8   status
INTEGER i  !apg NEW

!apg CALL IQAERRUNI(diagnosticFileID)  !apg set after OPEN

CALL QASETUP( 'DRSPALL',' ',' ','John Schatz',' ')

CALL QABANNER (0,' ',' ',' ')

! Setup to read filenames from command line
 CALL FILDFNAM('Input Control file',     'in', 'req','-user',  ' ')
 CALL FILDFNAM('Diagnostic file',        'out','req','-debug', ' ')
 CALL FILDFNAM('Input CAMDAT file',      'in', 'opt','-input', ' ')
 CALL FILDFNAM('Output CAMDAT file',     'out','req','-output',' ')

 CALL FILRDNAMS(FILES,IERR)
 CALL FILWRNAMS(0,FILES)

 IF(IERR.NE.0) THEN
   CALL FILWRNAMS(0,FILES)
!   CALL FILWRNAMS(diagnosticFileID,FILES)
   CALL QAABORT(' Incorrect file assignments')
 END IF

 fullparameterFilename = FILES(1)
 diagnosticFilename    = FILES(2)
 inCDBFilename         = FILES(3)
 outCDBFilename        = FILES(4)
 i = INDEX (outCDBFilename, '.', .TRUE.)  !apg
 validationFilePrefix = outCDBFilename(:i-1)  !apg

! IF no CDB files specified than set switch used to bypass
 CDBinput   = .TRUE.
 CDBoutput  = .TRUE.
 IF(inCDBFilename  .EQ. '')CDBinput  = .FALSE.

 OPEN(parameterFileID,  FILE=fullparameterFilename,FORM='FORMATTED', &
      STATUS='OLD', ACTION = 'READ') !dkr

 if(machine == 'PC')then
   status = 'REPLACE'
 else
   status = 'replace'
 endif


 OPEN(diagnosticFileID, FILE=diagnosticFilename,   FORM='FORMATTED', &
      STATUS=status, RECL=2000) !dkr

 CALL IQAERRUNI(diagnosticFileID)  !apg V1.22 set after OPEN
 CALL QAPAGE   (diagnosticFileID,' ')
 CALL QABANNER (diagnosticFileID,'DRSPALL ',' ',' ')
 CALL QADOEDIS (diagnosticFileID,'*')
 CALL FILWRNAMS(diagnosticFileID,FILES)



!number of extra CDB blocks for DRSPALL stuff
 neb  = 4
!number of parameters added to CDB as properties
 numParam = nep
 NOUT = diagnosticFileID

! Initialize CAMSUPES dynamic memory
   CALL MDINIT (RMEM)
   CALL MCINIT (CMEM)
   CALL DBSETUP(RMEM,CMEM,IERR)



 IF(CDBinput)THEN
   CALL inputCDBdata

 ELSE
   NINFO    = 3
   NQAREC   = 0
   NELBLKin = 0
   NUMPRPin = 0
   NUMPRP   = numParam
   NUMATRin = 0
   HEAD = 'WIPP DRSPALL OUTPUT DATA'
 ENDIF

 CALL LoadDefaultParameters
 CALL readControlFile
 CALL checkParameterBounds




!open validation test case files
If (validationTestCase > 0)then
  Call OpenTCFiles
endif

NOUT = diagnosticFileID


! Sets up repository and Wellbore models
! moved from RUNCALC
   Call Reset(RUNNING)



! open  output CDBs

   CALL createCDBHeader

   CALL createCDBVariables

   DEALLOCATE (INFO,    STAT=istat)
   DEALLOCATE (QAREC,   STAT=istat)
   DEALLOCATE (NUMELB,  STAT=istat)
   DEALLOCATE (IDELB,   STAT=istat)
   DEALLOCATE (NUMLNK,  STAT=istat)
   DEALLOCATE (PRPname, STAT=istat)
   DEALLOCATE (ISPROK,  STAT=istat)
   DEALLOCATE (PROP,    STAT=istat)
   DEALLOCATE (ISEVOK,  STAT=istat)
   DEALLOCATE (xn,      STAT=istat)
   DEALLOCATE (yn,      STAT=istat)
   DEALLOCATE (Link ,   STAT=istat)
   DEALLOCATE (EleMap,  STAT=istat)
   DEALLOCATE (ELEVAT,  STAT=istat)
   DEALLOCATE (THICK,   STAT=istat)
   DEALLOCATE (DEL_X,   STAT=istat)
   DEALLOCATE (DEL_Y,   STAT=istat)
   DEALLOCATE (ISATOK,  STAT=istat)
   DEALLOCATE (NAMELB,  STAT=istat)


 Call writeParameters(diagnosticFileID)

RETURN
END






SUBROUTINE readControlFile
!-----------------------------------------------------------------------
! Reads input control parameters
!
! Format is very similar to PC Version so data can be easily moved between
! PC and VMS versions. If CDB flag is set, then only a minimal set of
! control parameters are read that would be constant across vectors
! in a PA analysis.
!
! By:   David K.Rudeen
!       GRAM, Inc
!       dkrudee@sandia.gov
!
! Modification History
!       12/17/02 DKR: Original
!-----------------------------------------------------------------------

 USE Globals
 USE CDBGlobals
 Implicit None
 Character string*132, substring*20
 Character*20 readcfa
 Character*80 runString
 Integer      errorFlag
 Real(8)      READCF, temp

!----------------------------------
! Echo input to diagnostic file
!----------------------------------
      WRITE (diagnosticFileID,1000)
 1000 FORMAT(/,' Input Echo',/ &
              ,' ----------')

    5 Read  (parameterFileID, '(a)',END=10) string
        Write (diagnosticFileID,'(1x,a)') TRIM(string)
        GO TO 5
   10 Rewind(parameterFileID)



!-----------------------------------------------------------------
! Read parameters from input control file.  Logic
! will read data copied from input parameter file for PC Version
!-----------------------------------------------------------------

      WRITE (diagnosticFileID,1001)  !apg V1.22
 1001 FORMAT(/,' Input Echo (with numeric values)',/ &  !apg V1.22
              ,' --------------------------------')  !apg V1.22

!************
! Repository
!************

 CALL FINDKW(parameterFileID,'REPOSITORY',errorFlag)
 IF(errorFlag > 0)THEN
  write(diagnosticFileID,*)'ERROR: Unable to find Keyword', 'REPOS'
  call QAABORT ('ERROR: Unable to find Keyword')  !apg was STOP
 ENDIF

 surfaceElevation          = READCF(parameterFileID,surfaceElevation,'surfaceElevation')
 repositoryTop             = READCF(parameterFileID,repositoryTop,'repositoryTop')
 repositoryThickness0      = READCF(parameterFileID,repositoryThickness,'repositoryThickness')
 dRZThickness              = READCF(parameterFileID,dRZThickness,'dRZThickness')
 dRZPerm                   = READCF(parameterFileID,dRZPerm,'dRZPerm')
 repositoryOuterRadius     = READCF(parameterFileID,repositoryOuterRadius,'repositoryOuterRadius')
 repositoryInitialPressure = READCF(parameterFileID,repositoryInitialPressure,'repositoryInitialPressure')
 farFieldPorePressure      = repositoryInitialPressure
 farFieldStress            = READCF(parameterFileID,farFieldStress,'farFieldStress')

!************
! Waste
!************

 CALL FINDKW(parameterFileID,'WASTE',errorFlag)
 IF(errorFlag > 0)THEN
  write(diagnosticFIleID,*)'ERROR: Unable to find Keyword', 'WASTE'
  call QAABORT ('ERROR: Unable to find Keyword')  !apg was STOP
 ENDIF

 repositoryInitialPorosity = READCF(parameterFileID,repositoryInitialPorosity,'repositoryInitialPorosity')
 repositoryInitialPerm     = READCF(parameterFileID,repositoryInitialPerm,'repositoryInitialPerm')
 forchBeta                 = READCF(parameterFileID,forchBeta,'forchBeta')
 biotBeta                  = READCF(parameterFileID,biotBeta,'biotBeta')
 poissonsRatio             = READCF(parameterFileID,poissonsRatio,'poissonsRatio')
 cohesion                  = READCF(parameterFileID,cohesion,'cohesion')
 frictionAngle             = READCF(parameterFileID,frictionAngle,'frictionAngle')
 tensileStrength           = READCF(parameterFileID,tensileStrength,'tensileStrength')
 Lt                        = READCF(parameterFileID,Lt,'Lt')
 particleDiameter          = READCF(parameterFileID,particleDiameter,'particleDiameter')
 gasViscosity              = READCF(parameterFileID,gasViscosity,'gasViscosity')
 gasBaseDensity = atmosphericPressure/(gasConstant*repostemp)

!************
! Mud
!************

 CALL FINDKW(parameterFileID,'MUD',errorFlag)
 IF(errorFlag > 0)THEN
  write(diagnosticFIleID,*)'ERROR: Unable to find Keyword', 'MUD'
  call QAABORT ('ERROR: Unable to find Keyword')  !apg was STOP
 ENDIF

 initialMudDensity           = READCF(parameterFileID,initialMudDensity,'initialMudDensity')
 mudViscosity                = READCF(parameterFileID,mudViscosity,'mudViscosity')
 wallRoughness(1)            = READCF(parameterFileID,wallRoughness(1),'wallRoughness(1)')
 wallRoughness(2)            = READCF(parameterFileID,wallRoughness(2),'wallRoughness(2)')
 mudSolidsMax                = READCF(parameterFileID,mudSolidsMax,'mudSolidsMax')
 mudSolidsViscosityExponent  = READCF(parameterFileID,mudSolidsViscosityExponent,'mudSolidsViscosityExponent')
! chokeEfficiency             = READCF(parameterFileID,chokeEfficiency,'chokeEfficiency')

!*******************
! Wellbore/Drilling
!*******************

 CALL FINDKW(parameterFileID,'WELL',errorFlag)
 IF(errorFlag > 0)THEN
  write(diagnosticFIleID,*)'ERROR: Unable to find Keyword', 'WELLbore'
  call QAABORT ('ERROR: Unable to find Keyword')  !apg was STOP
 ENDIF

 bitDiameter                 = READCF(parameterFileID,bitDiameter,'bitDiameter')
 originalbitdiameter         = bitDiameter
 pipeDiameter                = READCF(parameterFileID,pipeDiameter,'pipeDiameter')
 collarDiameter              = READCF(parameterFileID,collarDiameter,'collarDiameter')
 pipeInsideDiameter          = READCF(parameterFileID,pipeInsideDiameter,'pipeInsideDiameter')
 collarLength                = READCF(parameterFileID,collarLength,'collarLength')
 exitPipeLength              = READCF(parameterFileID,exitPipeLength,'exitPipeLength')
 exitPipeDiameter            = READCF(parameterFileID,exitPipeDiameter,'exitPipeDiameter')
 drillingRate                = READCF(parameterFileID,drillingRate,'drillingRate')
 initialBitAboveRepository   = READCF(parameterFileID,initialBitAboveRepository,'initialBitAboveRepository')
 mudPumpRate                 = READCF(parameterFileID,mudPumpRate,'mudPumpRate')
 maxPumpPressure             = READCF(parameterFileID,maxPumpPressure,'maxPumpPressure')
 dDZThickness                = READCF(parameterFileID,dDZThickness,'dDZThickness')
 dDZPerm                     = READCF(parameterFileID,dDZPerm,'dDZPerm')
 stopDrillingExitVolRate     = READCF(parameterFileID,stopDrillingExitVolRate,'stopDrillingExitVolRate')
 stopPumpingExitVolRate      = READCF(parameterFileID,stopPumpingExitVolRate,'stopPumpingExitVolRate')
 stopDrillingTime            = READCF(parameterFileID,stopDrillingTime,'stopDrillingTime')

!**************
!Computational
!**************

 CALL FINDKW(parameterFileID,'COMPUTATIONAL',errorFlag)
 IF(errorFlag > 0)THEN
  write(diagnosticFIleID,*)'ERROR: Unable to find Keyword', 'COMPutational'
  call QAABORT ('ERROR: Unable to find Keyword')  !apg was STOP
 ENDIF

 geometry                = READCFA(parameterFileID,'geometry')
 allowFluidization       = READCFA(parameterFileID,'allowFluidization')
 maxTime                 = READCF (parameterFileID,maxTime,'maxTime')
 initialReposZoneSize    = READCF (parameterFileID,initialReposZoneSize,'initialReposZoneSize')
    string = READCFA(parameterFileID,'radius,growthRate')
    read(string,*)reposRadius1, growthRate
 if (growthRate < 1.0 .or. growthRate > 1.0) then !trz
   call QAABORT ('Growth rate must be 1.0') !trz
 endif !trz
 initialWellZoneSize     = READCF (parameterFileID,initialWellZoneSize,'initialWellZoneSize')
 wellGrowthRate          = READCF (parameterFileID,wellGRowthRate,'wellGRowthRate')
 if (wellGrowthRate < 1.0 .or. wellGrowthRate > 1.0) then !trz
   call QAABORT ('Well growth rate must be 1.0') !trz
 endif !trz
 firstWellZone           = READCF (parameterFileID,dble(firstWellZone),'firstWellZone')  !apg real(,8)
 wellStabilityFactor     = READCF (parameterFileID,wellStabilityFactor,'wellStabilityFactor')
 reposStabilityFactor    = READCF (parameterFileID,reposStabilityFactor,'reposStabilityFactor')
 massDiffusionFactor     = READCF (parameterFileID,massDiffusionFactor,'massDiffusionFactor')
 momentumDiffusionFactor = READCF (parameterFileID,momentumDiffusionFactor,'momentumDiffusionFactor')
 validationTestCase  = 0


!************
! Optional Validation
!************

 CALL FINDKW(parameterFileID,'VALIDATION',errorFlag)
 IF(errorFlag > 0)THEN
   runString = ' STANDARD RUN'
 ELSE

! parse into test case and subcase numbers
   temp = READCF(parameterFileID,dfloat(validationTestCase),'ValidationTestCase')
   validationTestCase = temp
   validationSubcase  = 10.0*(temp -dfloat(validationTestCase)+0.01)
   if(validationSubcase == 0)then
     Write(runString,*) ' VALIDATION TEST CASE: ', validationTestCase
   else
     Write(runString,*) ' VALIDATION TEST CASE: ', validationTestCase, &
                        ' SUBCASE: ', validationSubcase
   endif


 ENDIF
!apg V1.22 write(diagnosticFIleID,'(/,a,/)')runString
!apg V1.22 write(*,'(/,a,/)')runString

!*******************************
! optional initial cavity radius
!*******************************

    CALL FINDKW(parameterFileID,'INITIAL CAVITY RADIUS ',errorFlag)
    if(errorFlag <= 0)then
      backspace parameterFileID
      inputCavityRadius = READCF (parameterFileID,inputCavityRadius,'inputCavityRadius')
    endif

!*******************************
! optional min Characteristic Velocity
!*******************************

    CALL FINDKW(parameterFileID,'MINIMUM CHAR ',errorFlag)
    if(errorFlag <= 0)then
      backspace parameterFileID
      minCharVel = READCF (parameterFileID,minCharVel,'minCharVel')
    endif


!*******************************
! optional min Number of zones per Lt
!*******************************

    CALL FINDKW(parameterFileID,'MINIMUM NUMB ',errorFlag)
    if(errorFlag <= 0)then
      backspace parameterFileID
      minNumLt = READCF (parameterFileID,minCharVel,'minCharVel')
    endif




!************
!parameters
!************

 CALL FINDKW(parameterFileID,'PARAMETER',errorFlag)
 IF(errorFlag > 0)THEN
  write(diagnosticFIleID,*)' Default parameter values will be used'
 ELSE
   Pi                   =READCF(parameterFileID,Pi,'Pi')
   AtmosphericPressure  =READCF(parameterFileID,AtmosphericPressure,'AtmosphericPressure')
   gravity              =READCF(parameterFileID,gravity,'gravity')
   GasConstant          =READCF(parameterFileID,GasConstant,'GasConstant')
   ReposTemp            =READCF(parameterFileID,ReposTemp,'ReposTemp')
   gasBaseDensity= atmosphericPressure/(gasConstant*repostemp)
   WaterCompressibility =READCF(parameterFileID,WaterCompressibility,'WaterCompressibility')
   WasteDensity         =READCF(parameterFileID,WasteDensity,'WasteDensity')
   SaltDensity          =READCF(parameterFileID,SaltDensity,'SaltDensity')
   ShapeFactor          =READCF(parameterFileID,ShapeFactor,'ShapeFactor')
   TensileVelocity      =READCF(parameterFileID,TensileVelocity,'TensileVelocity')
   BitNozzleNumber      =READCF(parameterFileID,BitNozzleNumber,'BitNozzleNumber')
   BitNozzleDiameter    =READCF(parameterFileID,BitNozzleDiameter,'BitNozzleDiameter')
   ChokeEfficiency      =READCF(parameterFileID,ChokeEfficiency,'ChokeEfficiency')

 ENDIF

 write(diagnosticFIleID,'(/,a,/)')runString  !apg V1.22 move to end
 write(*,'(/,a,/)')runString  !apg V1.22 move to end

RETURN
END






REAL(8) FUNCTION DBVALUE(ELBNAM,PROPNAM,DFVALUE)
!-----------------------------------------------------------------------
! Function for retrieving specified property from specified element
! on CDB. If not there, use default.
!
! By:   Jim Garner
!
! Modification History
!       12/17/02 DKR: Code converted to F90 syntax
!        4/23/03 DKR: modified for properties previously read into
!                     PROP array
!-----------------------------------------------------------------------
USE CDBGlobals
IMPLICIT NONE
INTEGER IP, ID, I
CHARACTER*8 PNAME(90),ELBNAM,PROPNAM
REAL(8) DFVALUE

! find property name
  DO IP=1,NUMPRPin
   IF(PROPNAM.EQ.PRPNAME(IP)) GO TO 120
  END DO

! use default
  GO TO 170


120 CONTINUE
! find element block name
      DO I=1,NELBLKin
        IF(ELBNAM.EQ.NAMELB(I)) THEN
          IF(.NOT.ISPROK(I,IP)) GO TO 170
          DBVALUE=PROP(I,IP)
          WRITE(6,   10000) ELBNAM,PROPNAM,DBVALUE
!apg V1.22          WRITE(NOUT,10000) ELBNAM,PROPNAM,DBVALUE
10000     FORMAT(' FROM CDB, ',A8,1X,A8,'       :',1PE13.4)
          GO TO 999
        END IF
      END DO

! Use default

  170 DBVALUE=DFVALUE
      WRITE(6,   10010) PROPNAM,DFVALUE
      WRITE(NOUT,10010) PROPNAM,DFVALUE

10010 FORMAT(' PROPNAM ',A,' NOT FOUND!! SET TO DEFAULT:', &
        1PE12.3)

  999 RETURN
      END






SUBROUTINE FINDKW(iu,kw,ierr)
!-----------------------------------------------------------------------
!     Finds record beginning with string kw on unit iu.
!     iu must be formatted file
!
!     Created: 12/8/00
!
!     By:      David K. Rudeen
!              UNM/NMERI
!
! Routines Called
!     located in CAMCON_LIB
!     STRCMPRS  - removes extra blanks from string and returns length
!     STRUPCASE - converts string to upper case
!
!-----------------------------------------------------------------------
USE Globals  !apg V1.22
IMPLICIT NONE

INTEGER iu, len, kwlen, LNBLNK, ierr
CHARACTER kw*(*), string*80

REWIND(iu)
ierr = 0

    CALL STRCMPRS (kw,kwlen)
    CALL STRUPCASE(kw)

100 READ(iu,'(a)',END=900)string
    CALL STRCMPRS (string,len)
    CALL STRUPCASE(string)
    IF(len >= kwlen .and. string(1:kwlen) .EQ. kw(1:kwlen)) THEN
      GO TO 990
    ELSE
      GO TO 100
    ENDIF

900 ierr = 1

990 CONTINUE  !apg V1.22
    if (ierr == 0) WRITE (diagnosticFileID,'(a)')  trim(kw)  !apg V1.22
    RETURN
    END



REAL(8) FUNCTION READCF(iu,default,varname)
!-----------------------------------------------------------------------
! Reads Real data after delimiter
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------
USE GLOBALS
IMPLICIT NONE
CHARACTER*80 string
CHARACTER*16 string1, string2
CHARACTER*8  property, material
CHARACTER*(*) varname
INTEGER iu, len, i
REAL(8) DBVALUE
REAL(8) default
LOGICAL fromCDB


READ(iu,'(a)')string
CALL STRUPCASE(string)
CALL GetDataStrings(string,string1, string2)


! check for number or CDB (material property) pair
 Call STRCMPRS(string1,len)
 IF( len > 0)THEN

   if(string1(1:7) == 'DEFAULT')THEN
    READCF = default
    WRITE (diagnosticFileID,111)  trim(varname), READCF, 'DEFAULT'  !apg V1.22
111 format (2x,a, t32,' :', 1pe14.6, :, 3x,'!!', a)  !apg V1.22
112 format (2x,a, t32,' : ', a, :, 3x,'!!', a)  !apg V1.22

   elseif(string1(1:1) >= 'A' .and. string1(1:1) <= 'Z') then
     ! CDB materail property
       material = string1

       CALL STRCMPRS(string2,len)
       IF(len >0)THEN
         property = string2
         READCF = DBVALUE(material, property, default)
         WRITE (diagnosticFileID,111)  trim(varname), READCF, trim(material)//' '//trim(property)  !apg V1.22

       else
         write(diagnosticFileID,*)' ERROR: Could not find property ',varname
         call QAABORT (' ERROR: Could not find property ')  !apg was STOP
       endif

   else

     READ (string1,*) READCF
     WRITE (diagnosticFileID,112)  trim(varname), trim(string1)  !apg V1.22

   endif

 ELSE
   write(diagnosticFileID,*)' ERROR: Could not find data on input record for ',varname
   call QAABORT (' ERROR: Could not find data on input record')  !apg was STOP

 ENDIF

RETURN
END



SUBROUTINE getDataStrings(string,string1, string2)
!-----------------------------------------------------------------------
! parses  string into two data strings separated by delimiter
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------
Implicit None
character*(*) string, string1, string2
Character*80  stringt
integer i1,ierr,i, lens, j
CHARACTER*1  delim(4)
DATA delim/':', '=', ',' ,' '/

 string1 = ''
 string2 = ''

 call STRCMPRS(string,lens)

! find delimiter between description and data
 i1= 0
 i = 1
 DO WHILE (i1 == 0 .AND. i <= 4)
  i1 = INDEX(string,delim(i))
  i = i + 1
 ENDDO

 if(i1 > 0) then

   stringt = string(i1+1:lens)

! use blank as deliminter
   do i = 1,lens
     do j = 1, 4
       if(stringt(i:i) == delim(j)) then
         stringt(i:i) = ' '
       endif
     enddo
   enddo

   CALL STRCMPRS(stringt,lens)

   i1 = INDEX(stringt,' ')

   if(i1 >0)then
     string1 = stringt(1:i1-1)

     if(i1+1 <= lens) then
       string2 = stringt(i1+1:lens)
       lens=len(string2)
       if(lens > 0)then

         i1 = INDEX(string2,' ')

         if(i1 >0)then
          string2(i1:lens) = ''
         endif
       endif
     endif
   endif
 endif

RETURN
END




CHARACTER*(*) FUNCTION READCFA(iu,varname)
!-----------------------------------------------------------------------
! Reads ASCII text string data and extracts data after delimiter
!
! By:   David Rudeen
!
! Modification History
!       12/19/02 DKR: Original Version
!-----------------------------------------------------------------------


USE GLOBALS
IMPLICIT NONE
CHARACTER*80 string
CHARACTER*1  delim(4)
CHARACTER*(*) varname
INTEGER iu,len, i, ifmt, i1

DATA delim/':','=','-',';'/


 READ(iu,'(a)')string
 CALL STRCMPRS (string,len)
 CALL STRUPCASE(string)

 i1= 0
 i = 1
 DO WHILE (i1 == 0 .AND. i <= 4)
  i1 = INDEX(string,delim(i))
  i = i + 1
 ENDDO


 IF( i1 > 0)THEN
   i1 = i1+1
   CALL STRCMPRS (string(i1:),len)

   READCFA = string(i1:)
   WRITE (diagnosticFileID,112)  trim(varname), trim(READCFA)  !apg V1.22
111 format (2x,a, t32,' :', 1pe14.6, :, 3x,'!!', a)  !apg V1.22
112 format (2x,a, t32,' : ', a, :, 3x,'!!', a)  !apg V1.22

 ELSE
   write(diagnosticFileID,*) &
       ' ERROR: Could not find delimiter for variable=',Trim(varname)
   call QAABORT ('ERROR: Could not find delimiter for variable')  !apg was STOP

 ENDIF

RETURN
END

!-----------------------------------------------------------------------------------------

Subroutine OpenTCFiles
!-----------------------------------------------------------------------
! Opens special validation test case output files
!
! By:   David Rudeen
!
! Modification History
!       6/03 DKR: Original Version
!-----------------------------------------------------------------------

Use Globals
Implicit None

Character status*8, blank*1

if(machine == 'PC')THEN
  status = 'REPLACE'
else
  status = 'replace'
endif


if(validationTestCase == 1)then

!apg  if(validationSubcase >0) then
!apg    WRITE(TC1ChanFileName ,"('DRS_TC1',i1,'_chan.dat')")validationSubcase
!apg  else
!apg    TC1ChanFileName ='DRS_TC1'//'_chan.dat'
!apg  endif
  TC1ChanFileName = trim(validationFilePrefix)//'_chan.dat'  !apg
  Open(chanValidationFileID, FILE= TC1chanFilename, RECL=2048, FORM='FORMATTED', STATUS=status)


 CALL QAPAGE  (chanValidationFileID,' ')
 write(chanValidationFileID,'(1x,a)') trim(TC1ChanFileName)  !apg new
!apg CALL QABANNER(chanValidationFileID,' ',' ',' ')
!apg CALL QADOEDIS(chanValidationFileID,'*')

endif



if(validationTestCase == 4)then
!apg  if(validationSubcase >0)then
!apg   WRITE(TC4CoupleFileName ,"('DRS_TC4',i1,'_coupling.dat')")    validationSubcase
!apg   WRITE(TC4StressFileName ,"('DRS_TC4',i1,'_stress.dat')")      validationSubcase
!apg   WRITE(TC4FluidFileName  ,"('DRS_TC4',i1,'_fluidization.dat')")validationSubcase
!apg   WRITE(TC4EjectFileName  ,"('DRS_TC4',i1,'_expulsion.dat')")   validationSubcase
!apg  else
!apg   TC4CoupleFileName ='DRS_TC4'//'_coupling.dat'
!apg   TC4StressFileName ='DRS_TC4'//'_stress.dat'
!apg   TC4FluidFileName  ='DRS_TC4'//'_fluidization.dat'
!apg   TC4EjectFileName  ='DRS_TC4'//'_expulsion.dat'
!apg  endif
  TC4CoupleFileName = trim(validationFilePrefix)//'_coupling.dat'  !apg
  TC4StressFileName = trim(validationFilePrefix)//'_stress.dat'  !apg
  TC4FluidFileName  = trim(validationFilePrefix)//'_fluidization.dat'  !apg
  TC4EjectFileName  = trim(validationFilePrefix)//'_expulsion.dat'  !apg
  TC4FluidTimeFileName  = trim(validationFilePrefix)//'_fluidization_time.dat' !trz

 Open(couplingValidationFileID,     FILE= TC4couplefilename, RECL=2048, FORM='FORMATTED', STATUS=status)
 Open(stressValidationFileID,       FILE= TC4stressfilename, RECL=2048, FORM='FORMATTED', STATUS=status)
 Open(fluidizationValidationFileID, FILE= TC4Fluidfilename,  RECL=2048, FORM='FORMATTED', STATUS=status)
 Open(expulsionValidationFileID,    FILE= TC4Ejectfilename,  RECL=2048, FORM='FORMATTED', STATUS=status)
 Open(fluidizationTimeValidationFileID, FILE= TC4FluidTimefilename,  RECL=2048, FORM='FORMATTED', &
                                        STATUS=status) !trz

  blank = " "


  ! coupling file
  CALL QAPAGE  (couplingValidationFileID,blank)
  write(couplingValidationFileID,'(1x,a)') trim(TC4couplefilename)  !apg new
!apg  CALL QABANNER(couplingValidationFileID,blank,blank,blank)
!apg  CALL QADOEDIS(couplingValidationFileID,'*')


  ! stress file
  CALL QAPAGE  (stressValidationFileID,blank)
  write(stressValidationFileID,'(1x,a)') trim(TC4stressfilename)  !apg new
!apg  CALL QABANNER(stressValidationFileID,blank,blank,blank)
!apg  CALL QADOEDIS(stressValidationFileID,'*')


  ! fluidization file
  CALL QAPAGE  (fluidizationvalidationFileID,blank)
  write(fluidizationValidationFileID,'(1x,a)') trim(TC4Fluidfilename)  !apg new
!apg  CALL QABANNER(fluidizationvalidationFileID,blank,blank,blank)
!apg  CALL QADOEDIS(fluidizationvalidationFileID,'*')


  ! expulsion file
  CALL QAPAGE  (expulsionValidationFileID,blank)
  write(expulsionValidationFileID,'(1x,a)') trim(TC4Ejectfilename)  !apg new
!apg  CALL QABANNER(expulsionValidationFileID,blank,blank,blank)
!apg  CALL QADOEDIS(expulsionValidationFileID,'*')

  ! fluidization time file
  CALL QAPAGE  (fluidizationTimeValidationFileID,blank)
  write(fluidizationTimeValidationFileID,'(1x,a)') trim(TC4FluidTimefilename)  !trz

endif





if(validationTestCase == 5)then
!apg  if(validationSubcase > 0)then
!apg    write(TC5WellFileName,"('DRS_TC5',i1,'_wellbore.dat')")validationSubcase
!apg  else
!apg    TC5WellFileName='DRS_TC5'//'_wellbore.dat'
!apg  endif
  TC5WellFileName = trim(validationFilePrefix)//'_wellbore.dat'  !apg
  Open(wellboreValidationFileID,    FILE= TC5Wellfilename, RECL=2048, FORM='FORMATTED', status = status)

  CALL QAPAGE  (wellboreValidationFileID,' ')
  write(wellboreValidationFileID,'(1x,a)') trim(TC5wellfilename)  !apg new
!apg  CALL QABANNER(wellboreValidationFileID,' ',' ',' ')
!apg  CALL QADOEDIS(wellboreValidationFileID,'*')

endif


RETURN
END



!---------------------------------------------------
! These routine satify externals used in PC version
! Retained in VMS  for compatibility
!----------------------------------------------------
Subroutine PlotControl
return
end
Subroutine RunShutDown
return
end
Subroutine SetupPlots(i)
return
end
