!-----------------------------------------------------------------------------------------
!
! File Parameters.f90, contains routines for setting up run parameters
!
!-----------------------------------------------------------------------------------------

Subroutine LoadDefaultParameters

Use Globals
Implicit NONE



! Repository
surfaceElevation          = 1037.3
repositoryTop             = 385.31
repositoryThickness0      = 0.0
dRZThickness              = 0.85
dRZPerm                   = 1.0d-15  !apg e-
repositoryOuterRadius     = 19.2
repositoryInitialPressure = 14.5d6  !apg e6
farFieldPorePressure      = repositoryInitialPressure
farfieldStress            = 14.8d6  !apg e6



! Waste
repositoryInitialPorosity = 0.575
repositoryInitialPerm     = 1.7d-13  !apg e-
forchBeta                 = 1.15d-6  !apg e-
biotBeta                  = 1.0
poissonsRatio             = 0.35
cohesion                  = 18.85*6895.0
frictionAngle             = 45.0
tensileStrength           = 1.0*6895.0
Lt                        = 0.02
particleDiameter          = 0.001
gasViscosity              = 10.0d-6  !apg e-



! Mud
initialMudDensity         = 1210.0
mudViscosity              = 0.01
wallRoughness(1)          = 0.5d-5
wallRoughness(2)          = 0.5d-5
mudSolidsmax              = 0.6
mudSolidsViscosityExponent= -1.8


! Wellbore/Drilling
bitDiameter               = 0.3112
pipeDiameter              = 0.1143
collarDiameter            = 0.2032
pipeInsideDiameter        = 0.1005
collarLength              = 182.9
exitPipeLength            = 15.24
exitPipeDiameter          = 0.0508
drillingRate              = 0.005
initialBitAboveRepository = 0.15
mudPumpRate               = 0.021443  ! where 0.00631 = 100 gal/min
maxPumpPressure           = 27.5d6  !apg e6    !4000 psi
dDZThickness              = 0.15
dDZPerm                   = 1.0d-14  !apg e-
stopDrillingExitVolRate   = 1000.0
stopPumpingExitVolRate    = 1000.0
stopDrillingTime          = 1000.0


! Computational
geometry                = 'S'
allowFluidization       = 'Y'
maxTime                 = 900.0
initialReposZoneSize    = 0.002
reposRadius1            = 0.5
growthRate              = 1.01
initialWellZoneSize     = 2.0
firstWellZone           = 2
wellStabilityFactor     = 0.1
reposStabilityFactor    = 5.0
massDiffusionFactor     = 0.02
momentumDiffusionFactor = 0.002
ValidationTestCase      = 0
inputCavityRadius       = 0.


!dkr - constants moved from global parameters
Pi                   = 3.141592653589793D0
AtmosphericPressure  = 101300.0d0
gravity              = 9.8067d0
GasConstant          = 4124.0d0
ReposTemp            = 300.0d0
gasBaseDensity       = atmosphericPressure/(gasConstant*repostemp)  !0.082
WasteDensity         = 2650.0d0
SaltDensity          = 2201.0d0
TensileVelocity      = 1000.0d0
BitNozzleNumber      = 3.0d0
BitNozzleDiameter    = 0.4375d0*0.0254d0 !=1.1113e-2
ChokeEfficiency      = 0.9d0
WaterCompressibility = 3.1d-10
ShapeFactor          = 1.0d0
TensileVelocity      = 1000.0d0
minCharVel           = 1.0d-6  !apg e-
minNumLt             = 5

uncompactedWastePorosity = 0.85
RepositoryHeight         = 3.96

FIX = 1
ENG = 2



! APG Plot axes minimum / maximum
! left here for compatibility with PC Version
minPlotRadius        = 0.0
maxPlotRadius        = 10.0
minPlotPressure      = 12.0d6  !apg e6
maxPlotPressure      = 16.0d6  !apg e6
minPlotStress        = -2.0d6  !apg e6
maxPlotStress        = 8.0d6  !apg e6
minPlotVelocity      = 0.0
maxPlotVelocity      = 2.0
minPlotWellPos       = 0.0
maxPlotWellPos       = 1500.0  !modify for validation test 2 and 5
minPlotWellPressure  = 0.0
maxPlotWellPressure  = 50.0d6  !apg e6  !modify for validation test 2
minPlotWellVelocity  = 0.0
maxPlotWellVelocity  = 10.0
minPlotWellRho       = 0.0
maxPlotWellRho       = 2000.0
minPlotWellFraction  = 0.0
maxPlotWellFraction  = 1.0
newPlotRadius        = .TRUE.
newPlotPressure      = .TRUE.
newPlotStress        = .TRUE.
newPlotVelocity      = .TRUE.
newPlotWellPos       = .TRUE.
newPlotWellRho       = .TRUE.
newPlotWellPressure  = .TRUE.
newPlotWellVelocity  = .TRUE.
newPlotWellFraction  = .TRUE.


return
end

!*****************************************************************************
Subroutine checkParameterBounds

Use Globals
Implicit NONE

Real(8) temp
Integer ie, boundcheck

!dkr changed for QE0110
IF(validationTestCase > 0)RETURN


! Repository
ie = ie+ boundcheck('surfaceElevation',      surfaceElevation,    900.0D0, 1100.0D0)
ie = ie+ boundcheck('repositoryTop',         repositoryTop,       300.0D0, 500.0D0)
! 0. is a flag for calculating repository thickness, so it is acceptable
if(repositoryThickness <= 0.0D0) then
else
  temp = repositoryThickness
  ie = ie+ boundcheck('repositoryThickness',   temp,   1.0D0, 2.0D0)
endif
ie = ie+ boundcheck('dRZThickness',          dRZThickness,          &
                     0.1D0, 1.0D0)
ie = ie+ boundcheck('dRZPerm',               dRZPerm,               &
                     1.0D-21, 1.0D-12)
ie = ie+ boundcheck('repositoryOuterRadius', repositoryOuterRadius, 10.0D0, 50.0D0)
!dkr changed for QE0110
ie = ie+ boundcheck('repositoryInitialPressure', repositoryInitialPressure, &
                     6.0D6, 18.0D6)
ie = ie+ boundcheck('farfieldStress',             farfieldStress,   &
                     14.0D6, 18.0D6)



! Waste
ie = ie+ boundcheck('repositoryInitialPorosity', repositoryInitialPorosity, &
                     0.20D0, 0.85D0)
!dkr changed for QE0110
ie = ie+ boundcheck('repositoryInitialPerm',     repositoryInitialPerm, &
                     1.0D-17, 1.0D-11)
!dkr changed for QE0110
! RM2019 - changed to be a proper inequality comparison for REAL
if (ABS(forchBeta) > 0.0D0) then
  temp = ABS(forchBeta)
  ie = ie+ boundcheck('forchBeta',     temp,            1D-12, 1.0D-4)
endif
ie = ie+ boundcheck('biotBeta',        biotBeta,         0.5D0, 1.0D0)
ie = ie+ boundcheck('poissonsRatio',   poissonsRatio,    0.2D0, 0.45D0)
ie = ie+ boundcheck('cohesion',        cohesion,         1.0D4, 1.0D7)
ie = ie+ boundcheck('frictionAngle',   frictionAngle,    30.0D0, 60.0D0)
ie = ie+ boundcheck('tensileStrength', tensileStrength,  1.0D4, 6.91D6)
ie = ie+ boundcheck('Lt', Lt,  -0.0001D0, 0.1D2)
!dkr changed for QE0110
ie = ie+ boundcheck('particleDiameter',particleDiameter, 1.0D-3, 1.0D0)
ie = ie+ boundcheck('gasViscosity' ,   gasViscosity,     1.0D-6, 1.0D-4)

! Mud
ie = ie+ boundcheck('initialMudDensity',initialMudDensity, 1000.0D0, 1500.0D0)
ie = ie+ boundcheck('mudViscosity',     mudViscosity,      1.0D-4, 1.0D-1)
ie = ie+ boundcheck('wallRoughness1',    wallRoughness(1),     1.0D-14, 1.0d-2)
ie = ie+ boundcheck('wallRoughness2',    wallRoughness(2),     1.0D-14, 1.0D-2)
ie = ie+ boundcheck('mudSolidsmax',     mudSolidsmax,      0.5D0, 0.8D0)
ie = ie+ boundcheck('mudSolidsViscosityExponent',mudSolidsViscosityExponent, &
                     -2.0D0, -1.2D0)

! Wellbore/Drilling
ie = ie+ boundcheck('bitDiameter',       bitDiameter,       0.25D0, 1.32D0)
ie = ie+ boundcheck('pipeDiameter',      pipeDiameter,      0.1D0, 0.15D0)
ie = ie+ boundcheck('collarDiameter',    collarDiameter,    0.15D0, 0.25D0)
ie = ie+ boundcheck('pipeInsideDiameter',pipeInsideDiameter,0.08D0, 0.13D0)
ie = ie+ boundcheck('collarLength',      collarLength,       0.0D0, 250.0D0)
ie = ie+ boundcheck('exitPipeLength',    exitPipeLength,     0.0D0, 25.0D0)
ie = ie+ boundcheck('exitPipeDiameter',  exitPipeDiameter,   0.0254D0, 0.3048D0)

ie = ie+ boundcheck('drillingRate',      drillingRate,      0.00D0, 0.01D0)
ie = ie+ boundcheck('initialBitAboveRepository',initialBitAboveRepository, &
                     0.01D0, 1.0D0)
ie = ie+ boundcheck('mudPumpRate',       mudPumpRate,       0.001D0, 0.1D0)
ie = ie+ boundcheck('maxPumpPressure',   maxPumpPressure,   1.0D5, 100.0D6)
ie = ie+ boundcheck('dDZThickness',      dDZThickness,      0.1D0, 0.3D0)
ie = ie+ boundcheck('dDZPerm',           dDZPerm,           1.0D-21, 1.0D-12)
ie = ie+ boundcheck('stopDrillingExitVolRate',stopDrillingExitVolRate, &
                     0.05D0, 2000.D0)
ie = ie+ boundcheck('stopPumpingExitVolRate', stopPumpingExitVolRate, &
                     1.0D0, 2000.D0)



if(ie > 0)then
  write(*,*)' ERROR: Parameter Bounds Exceeded'
  write(diagnosticFIleID,*)' ERROR: Parameter Bounds Exceeded'
  call QAABORT ('Parameter Bounds Exceeded')  !apg was STOP
endif

return
end
!*****************************************************************************

Integer Function boundCheck(varname, var, minval, maxval)

USE globals
Implicit None
Character*(*) varname
Real(8) var
Double Precision minval, maxval  !apg REAL

!DKR changed for QA0110
!if( var < 0.9999*minval .OR. var > 1.0001*maxval)then
if( var < minval-ABS(0.001*minval) .OR. var > maxval+ABS(0.001*maxval))then
  boundcheck = 1
  write(diagnosticFileID,*) &
     ' ERROR:', Trim(varname),'=',var,' Outside parameter bounds, (',minval,',',maxval,')'
else
  boundcheck = 0
endif

return
end



