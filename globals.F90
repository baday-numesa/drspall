Module Globals
!-------------------------------------------------------------------------
!
! File Globals.f90 -- contains global parameters, variables, and constants
!
!--------------------------------------------------------------------------

!Use DFWin
!Use DFLib
!Use GraphicsRecords
Implicit NONE


! Input Parameters
! ----------------

! Repository
Real(8) surfaceElevation, repositoryTop, repositoryThickness, dRZThickness, &
  dRZPerm, uncompactedWastePorosity, repositoryHeight, repositoryThickness0, &
  repositoryOuterRadius, repositoryInitialPressure, farFieldPorePressure, &
  farfieldStress, reposRadius1, growthRate

! Waste
Real(8) repositoryInitialPorosity, repositoryInitialPerm, forchBeta
Real(8) biotBeta, poissonsRatio, cohesion, frictionAngle, tensileStrength, &
        particleDiameter, fluidizationVelocity, Lt

! Gas
Real(8) gasBaseDensity, gasViscosity

! Mud
Real(8) initialMudDensity, mudViscosity, wallRoughness(2), mudSolidsMax,  &
        mudSolidsViscosityExponent

! Wellbore/Drilling
Real(8) bitDiameter, pipeDiameter, collarDiameter, pipeInsideDiameter, collarLength, &
  drillingRate, initialBitAboveRepository, mudPumpRate, dDZThickness, dDZPerm, &
  stopDrillingExitVolRate, stopPumpingExitVolRate, stopDrillingTime, timeToFloor, &
  maxPumpPressure, initialPumpRate, pitGain

! Computational
Character(1) geometry, allowFluidization
Real(8) maxTime, initialReposZoneSize, initialWellZoneSize, &
  wellStabilityFactor, reposStabilityFactor

! APG Plot axes minimums / maximums (left here for compatibility with PC version)
Real(8) minPlotRadius, maxPlotRadius, minPlotWellPos, maxPlotWellPos, &
  minPlotPressure, maxPlotPressure, minPlotStress, maxPlotStress, &
  minPlotVelocity, maxPlotVelocity, &
  minPlotWellPressure, maxPlotWellPressure, &
  minPlotWellVelocity, maxPlotWellVelocity, &
  minPlotWellRho, maxPlotWellRho, &
  minPlotWellFraction, maxPlotWellFraction
Logical newPlotRadius, newPlotWellPos, &
  newPlotPressure, newPlotStress, newPlotVelocity, &
  newPlotWellPressure, newPlotWellVelocity, newPlotWellRho, newPlotWellFraction

! Internal Parameters and Variables
! ---------------------------------

! Repository array
!dkr Logical(1)
Logical , DIMENSION(:), ALLOCATABLE :: &
  tensileFailureStarted, tensileFailureCompleted, shearFailed, &
  fluidizationStarted, fluidizationCompleted
REAL(8), DIMENSION(:), ALLOCATABLE :: &
  reposRadius, reposRadiusH, reposPres, porosity, permeability, forchratio, &
  poreVelocity, superficialVelocity, reposFactor, pL, &
  psi, dCoeff, radEffStress, tanEffStress, shearStress, &
  fractionTensileFailed, tensileFailureTime, meanEffStress, shearStrength, &
  fractionFluidized, fluidizationTime, FluidStartTime, FluidStopTime, FailStartTime, &
  radElasticStress, tanElasticStress, radSeepageStress, tanSeepageStress, &
  reposDR, reposDRH, reposVol, reposGasMass, drillingFailure
Real(8) , DIMENSION(:), ALLOCATABLE :: aa, bb, cc, rr, gam


! Repository non-array
!dkr Integer(2)
Integer geomExponent
!dkr Integer(4)
! 2026 BC-revision: stress problem is posed on the intact-waste domain with
! the radial-stress inner BC applied at the failed/intact interface.
!   firstFailedZone : inner edge of the failed/fluidizing annulus. Advances
!                     on fluidization completion.
!   firstIntactZone : inner edge of the intact waste; stress anchor.
!                     Advances in Lt-batch-sized jumps on batch completion.
!   batchEndZone    : outer cell of the active Lt tensile-failure batch; 0
!                     when no batch is active (surfaceFailureAllowed=.TRUE.).
! When no failed annulus exists, firstIntactZone == firstFailedZone.
Integer numReposZones, maxTensileFailedIndex, maxShearFailedIndex, &
           firstFailedZone, firstIntactZone, &
           interface1, interface2, lastFailedZone, firstInternalFailedZone, &
           minNumLt
Real(8) reposZoneSize, wasteBoundaryPoreVelocity, sumReposGasMass, &
        totalGasFromWaste, totalWasteFromRepos, maxForchRatio

!dkr Integer(4), Logical(1)
Integer batchEndZone
Logical surfaceFailureAllowed

! Cavity
Real(8) cavityPres, initialCavityRadius, initialCavityArea, initialCavityVol, initialCavityGasVol, &
  cavityGasMass, cavityWasteVol, cavityWasteMass, deltaGasIntoWell, deltaWasteIntoWell, &
  deltaSaltIntoWell, deltaGasFromWaste, exitPoreVelocity, exitPoreArea, &
  inputCavityRadius, cutMass, splMass, totMass, splvoleq, cutvoleq, totvoleq, &
  cutVol, splVol, splvol2, TotVol, cutTrueVol, splmass2, beddepth, minCharVel

! Wellbore array
REAL(8), DIMENSION(:), ALLOCATABLE :: &
  wellPos, wellPres, wellRhoV, wellVol, WellArea, wellZoneSize, &
  wellV, wellRho, wellAreaInt, wellDeltaVInt, wellFactor, &
  wellMudMass, wellGasMass, wellWasteMass, wellSaltMass, &
  wellMudVol, wellGasVol, wellWasteVol, wellSaltVol, hydraulicDia
REAL(8), DIMENSION(:), ALLOCATABLE :: &
  wellMassOrig, wellPresOrig, wellMudMassOrig, &
  wellGasMassOrig, wellWasteMassOrig, wellSaltMassOrig
REAL(8), DIMENSION(:), ALLOCATABLE :: &
  wellGRhoStar, wellGRhoVStar, wellGRhoT, wellGRhoVT, &
  wellVolInt, momDiff, massDiff

! Wellbore non-array
!dkr Logical(1), Integer(4)
Logical repositoryPenetrated, firstPenetration, drilling, pumping, exitPipe, &
           internalfailure

Integer numWellZones, wellBottomIndex, numExitZones, firstWellZone, &
   numWellZones1, numWellZones2, numWellZones3, numWellZones4

Real(8) wellLength, bitAboveRepository, dRZDiameter, dRZArea, wellGrowthRate,  &
  bitArea, bitFlowArea, pipeArea, collarAnnulusArea, pipeAnnulusArea, &
  exitPipeLength, exitPipeDiameter, exitPipeArea, wellDepth, &
  mudMassEjected, gasMassEjected, wasteMassEjected, saltMassEjected, &
  totalMudInWell, totalGasInjected, totalGasInWell, totalWasteInWell, totalSaltInWell, &
  cuttingsMassMax, cuttingsRadiusMax, equivDrilledCavityRadius, &
  volStore, wasteStore, gasStore, wasteInjected,originalBitDiameter, &
  sumWellGasMass,   sumWellWasteMass, Z1, Z2, Z3, Z4

! General array
REAL(8), DIMENSION(:), ALLOCATABLE :: invPorosity

! General non-array
Real(8) invGasViscosity, invWellZoneSize, invWasteDensity, invSaltDensity, invWaterCompressibility, &
  invInitialMudDensity, invGasBaseDensity, invInitialCavityGasVol, invPi, invAtmosphericPressure, &
  invReposZoneSize, timeOfPenetration

! Computational Parameters and Constants
!---------------------------------------

! Run management
Real(8) runTime, deltaTime
Double Precision cpuTime, cpuBegin,cpuEnd, &  !apg REAL
        cpuT0, cpuT1, cpuT2, cpuT3, cpuWBF, cpuWSF, cpuWSS

!dkr Integer(4), Logical(1)
Integer runIndex
Logical stopRun, isRunning, pauseRun, CDBoutput, textOutput
Character*3 platform, machine
Character*132 root

! Window unit constants
!dkr Integer(4)
Integer ParameterUnit, TimeUnit, RepositorySpatialPressureUnit, RepositorySpatialStressUnit, &
  RepositorySpatialVelocityUnit, WellboreSpatialPressureUnit, WellboreSpatialMassUnit, &
  ResultsSummaryUnit

! File logic
!Type(T_OpenFilename) fName
!Logical(1) fileLoaded, fileSaved, fileNamed
Character(255) fullParameterFilename, pressureFileName, radEffStressFileName, &
  tanEffStressFileName, poreVelocityFileName, timeFileName, summaryFileName, outParamsFileName,&
   diagnosticFileName, validationFilePrefix
!apg Added validationFilePrefix for files below
Character(255) TC1chanFileName, TC4couplefilename, TC4stressfilename, &
   TC4Fluidfilename, TC4Ejectfilename, TC4FluidTimefilename, TC5Wellfilename

!dkr Integer(4)
Integer parameterFileID, pressureFileID, radEffStressFileID, tanEffStressFileID, &
  poreVelocityFileID, timeFileID, summaryFileID, outParamsFileID, diagnosticFileID, &
  chanValidationFileID, &
  !test case #4
  couplingValidationFileID, stressValidationFileID, fluidizationvalidationFileID, &
  expulsionValidationFileID, fluidizationTimeValidationFileID, &
  !test case #5
  wellboreValidationFileID

! Character-based user input parameter values
Character(30) gridItems(2,44)
Character(10) gridUnits(2,44)
Character(14) gridValues(2,44)

! Parameter input grids
!dkr Integer(2)
Integer numGridRows, numGridCols, rowHeight, margin
Integer gridWidth, colWidth(3)

! DKR Graphics records moved to seperate module for portatibility with VMS


! Plot and save times
!dkr Real(4), Integer(4)
Real(8) spatialPlotTime(1000), spatialSaveTime(0:1000), &  !apg Real
   timeSaveTime(0:2000), chanSaveTime (10)
Integer spatialPlotIndex, spatialSaveIndex, timeSaveIndex, &
   chanSaveIndex

! APG for time plots
Integer ixtplt
!dkr Real(4) times 3
Real(8) finaltplt  !apg Real
Integer, Parameter :: maxtplt = 1000
Real(8) tpltTime(0:maxtplt)  !apg Real
Real(8) tpltPumpPressure(0:maxtplt), tpltWellBottomPressure(0:maxtplt), tpltCavityPressure(0:maxtplt), &  !apg Real
  tpltDrilledRadius(0:maxtplt), tpltCavityRadius(0:maxtplt), tpltTensileRadius(0:maxtplt), &
  tpltWasteBoundaryPoreVelocity(0:maxtplt), tpltFluidizationVelocity(0:maxtplt), tpltMudEjectionVelocity(0:maxtplt), &
  tpltWasteInWell(0:maxtplt), tpltWasteEjected(0:maxtplt), &
  tpltGasInjected(0:maxtplt), tpltGasInWell(0:maxtplt), tpltGasEjected(0:maxtplt), &
  tpltGasPosInWell(0:maxtplt), tpltWastePosInWell(0:maxtplt)

! for plotting and saving
! time plots
!dkr Real(4)
Real(8) curTime, &  !apg Real
  curPumpPressure, curWellBottomPressure, curCavityPressure, &
  curDrilledRadius, curCavityRadius, curTensileRadius, &
  curWasteBoundaryPoreVelocity, curFluidizationVelocity, curMudEjectionVelocity, &
  curWasteInWell, curWasteEjected, &
  curGasInjected, curGasInWell, curGasEjected, &
  curGasPosInWell, curWastePosInWell

!dkr Real(4) times 2
Real(8) prevTime  !apg Real
Real(8) prevPumpPressure, prevWellBottomPressure, prevCavityPressure, &  !apg Real
  prevDrilledRadius, prevCavityRadius, prevTensileRadius, &
  prevWasteBoundaryPoreVelocity, prevFluidizationVelocity, prevMudEjectionVelocity, &
  prevWasteInWell, prevWasteEjected, &
  prevGasInjected, prevgasInWell, prevGasEjected, &
  prevGasPosInWell, prevWastePosInWell

! summary display
!dkr Real(4)
Real(8) displayedTime, &  !apg Real
  displayedPumpPressure, displayedBitAbove, &
  displayedWellBottomPressure, displayedCavityPressure, displayedWasteBoundaryPoreVelocity, &
  displayedFluidizationVelocity, displayedEquivDrilledCavityRadius, displayedCavityRadius, &
  displayedTensileFailedRadius, displayedShearFailedRadius, displayedGasInjected, &
  displayedMudInWell, displayedGasInWell, displayedWasteInWell,displayedSaltInWell, &
  displayedMudVelocity, displayedMudVolRate, displayedMudEjected, displayedGasEjected, &
  displayedWasteEjected, displayedSaltEjected, &
  displayedMudExitFraction, displayedGasExitFraction, displayedWasteExitFraction, displayedSaltExitFraction, &
  displayedGasPosInWell, displayedWastePosInWell, displayedradEffStress, &
  displayedFirstFailedZone, displayedFirstIntactZone

Integer cellControl(5), oldCellControl,stepControl, oldStepControl
Integer DisplayedCell
Character*8 stepControlName(5),displayedControl


! Validation Test Case specific parameters
!dkr Integer(2)
Integer validationTestCase, validationSubcase
Real(8)    characteristicTime

!test case #4
Real(8) stressSaveTime, stressSaveDelta, fluidizationSaveTime, fluidizationSaveDelta,&
		initialGasInRepository

!DKR - parameters switched to variables. set in reset
REAL(8) gravity, MassDiffusionFactor, MomentumDiffusionFactor, chokeEfficiency, &
        WaterCompressibility, GasConstant, ReposTemp, WasteDensity, saltDensity, &
        BitNozzleNumber,BitNozzleDiameter, Pi, AtmosphericPressure, shapeFactor, &
        TensileVelocity
!dkr Integer(2)
Integer FIX, ENG

! Constants
!dkr Interger(2) times 3
Integer, Parameter :: FULL_RESET = 1
Integer, Parameter :: RUNNING    = 2
Integer, Parameter :: REOPENING  = 3


end

!-----------------------------------------------------------------------------------------
