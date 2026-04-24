!-----------------------------------------------------------------------------------------
!
! File SetupCalc.f90, contains routines for setup calculations
!
!-----------------------------------------------------------------------------------------

Subroutine SetupWellbore

Use Globals
Implicit NONE

Integer    i, j, numZones, NZ  !apg Integer(4)
Real(8)    wellMudVolOld, rhoIntLocal, r8, Z

!function for number of zones to cover Z with growth
NZ(Z) = log(1.0-(Z/initialWellZoneSize)*(1.0-wellGrowthRate))/log(wellGrowthRate)

! initialize

mudMassEjected     = 0.0
gasMassEjected     = 0.0
wasteMassEjected   = 0.0
firstPenetration   = .FALSE.
totalMudInWell     = 0.0
totalGasInWell     = 0.0
totalGasInjected   = 0.0
totalWasteInWell   = 0.0
totalSaltInWell    = 0.0
totalGasFromWaste  = 0.0
totalWasteFromRepos= 0.0

!JFS3 automatically account for validation test case 2 (coal) by
!setting penetration, drilling and pumping rates in the input file

if (initialBitAboveRepository > 0.0) then
  repositoryPenetrated = .FALSE.
else
  repositoryPenetrated = .TRUE.
endif

if(drillingRate > 0.0) then
  drilling = .TRUE.
  timeToFloor = (repositoryThickness+initialBitAboveRepository)/drillingrate
  timeOfPenetration = initialBitAboveRepository/drillingRate
else
  drilling = .FALSE.
  timeToFloor = 1000.
  timeofPenetration = 0.0
endif


if (mudPumpRate > 0.0) then
  pumping  = .TRUE.
else
  pumping  = .FALSE.
endif

!write(diagnosticFileId,*)' Drilling time to floor (s): ', timeToFloor

! Initialize wellbore geometry

!****************************************************************************************
!	must set bit diameter to effective bit diameter in test case #1
!	see VVP for explanation
!	also retain original bit diameter for calculating dimensionless plotting parameter zeta
!	D.Lord, 01.22.2003

IF (validationTestcase .EQ. 1)THEN
!	originalBitDiameter = bitDiameter
	bitDiameter=2*sqrt(originalBitDiameter*repositoryThickness)
ENDIF
!**************************************************************************************

dRZDiameter = bitDiameter
dRZArea     = Pi*(0.5*dRZDiameter)**2
pipeArea    = Pi*0.25*pipeInsideDiameter**2
bitArea     = Pi*0.25*bitDiameter**2

Call CavityMaxes

bitFlowArea       = BitNozzleNumber*Pi*0.25*BitNozzleDiameter**2
collarAnnulusArea = bitArea-Pi*0.25*collarDiameter**2
pipeAnnulusArea   = bitArea-Pi*0.25*pipeDiameter**2

bitAboveRepository = initialBitAboveRepository
wellDepth          = surfaceElevation-repositoryTop

Z1= 0.5*(surfaceElevation-repositoryTop)
Z2= 0.5*collarLength
Z3= 0.5*(surfaceElevation-repositoryTop-collarlength)
Z4= 0.5*exitPipeLength

if(wellGrowthRate < 1.00001) then

  numWellZones1 = Z1/initialWellZoneSize
  numWellZones2 = Z2/initialWellZoneSize
  numWellZones3 = Z3/initialWellZoneSize
  numWellZones4 = Z4/initialWellZoneSize

else
  numWellZones1 = NZ(Z1)
  numWellZones2 = NZ(Z2)
  numWellZones3 = NZ(Z3)
  numWellZones4 = NZ(Z4)
endif


! add 1 psuedo zone for cavity at well bottom
numwellZones    = 2*(numwellzones1+numwellzones2+numwellzones3+numwellzones4)+1
wellBottomIndex = 2*numwellzones1+ 1




! special initialization if running validation test case 2 (coal)
!IF( validationTestcase .EQ. 2)THEN
!JFS3 use 12" pipe and 1/2 area, but not needed, is set in parameters
!  pipeAnnulusArea   = 0.5*7.26E-2
!  collarAnnulusArea = 0.35*7.26E-2
!  make zone 1 the cavity
!  wellBottomIndex = 1
!  numWellZones = numWellZones/2
!ENDIF



Call ZoneWellbore



wellLength= wellPos(numwellZones)+ 0.5*wellZoneSize(numWellZones)

!firstWellZone = 1 => entire wellbore domain
!firstWellZone > 1 => sets up full problem, pumps mud into bottom of well
!                     calculates and outputs from wellBottomIndex up
if(firstWellZone > 1) then
  firstWellZone = wellBottomIndex
endif



! Initialize densities and pressures
do i = 1, numWellZones
  wellGasMass  (i) = 0.0
  wellWasteMass(i) = 0.0
  wellSaltMass (i) = 0.0
  wellGasVol   (i) = 0.0
  wellWasteVol (i) = 0.0
  wellSaltVol  (i) = 0.0

! approximate starting pressures for viscous flow
!dkr - account for psuedo cell at well bottom
  if (i < wellBottomIndex) then
    !this is a very rough estimate
    wellPres(i) = 2.0*initialMudDensity*Gravity*wellDepth

  elseif (i == wellBottomIndex) then
    wellPres(i) = AtmosphericPressure + initialMudDensity*Gravity*welldepth
  elseif( i <= numWellZones- numExitZones) then
    wellPres(i) = AtmosphericPressure  &
                + initialMudDensity*Gravity*(wellLength-WellPos(i)-exitPipeLength)
  else
    ! horizontal exit pipe
    wellPres(i) = AtmosphericPressure
  end if

!dkr
!  wellRho    (i) = initialMudDensity*Exp(WaterCompressibility*(wellPres(i)-AtmosphericPressure))
  wellRho    (i) = initialMudDensity*(1.0d0+WaterCompressibility*(wellPres(i)-AtmosphericPressure))
  wellMudVol (i) = wellVol(i)
  wellMudMass(i) = wellRho(i)*wellVol(i)

end do

! initialize for test case 5
if(validationTestCase == 5) then
  wellRho (1) = initialMudDensity
  wellPres(1) = 6.0d6*mudPumpRate/0.02018 + atmosphericPressure
  Do i = 2, wellbottomIndex-1

    ! two step iteration
    wellRho (i) = wellRho (i-1)
    wellpres(i) = wellpres(i-1) + wellrho(i)*gravity*wellZoneSize(i)
    wellRho (i) = initialMudDensity*(1.0d0+WaterCompressibility*(wellPres(i)-AtmosphericPressure))
    wellPres(i) = wellpres(i-1) + wellrho(i)*gravity*wellZoneSize(i)
  enddo

    i = numWellZones+1
    wellpres(i) = atmosphericPressure
    wellRho (i) = initialMudDensity

  do i = numWellZones, wellBottomIndex+1, -1
    wellRho (i) = wellRho (i+1)
    wellpres(i) = wellpres(i+1) + wellrho(i)*gravity*wellZoneSize(i)
    wellRho (i) = initialMudDensity*(1.0d0+WaterCompressibility*(wellPres(i)-AtmosphericPressure))
    wellPres(i) = wellpres(i+1) + wellrho(i)*gravity*wellZoneSize(i)
  enddo

    i = wellBottomIndex
    wellRho (i) = wellRho (i+1)
    wellpres(i) = wellpres(i+1)

endif



! Initialize velocities assuming no viscous effects and steady mass flow rate
wellV   (1) = mudPumpRate/(Pi*0.25*pipeInsideDiameter**2)
wellRhoV(1) = wellRho(1)*wellV(1)
wellV   (0) = wellV   (1)
wellRhoV(0) = wellRhoV(1)
wellRho (0) = wellRho (1)
do i = 2, numWellZones
!dkr changed from vol ratio to area
  wellV   (i) = wellV(1)*wellArea(1)/wellArea(i)
  wellRhoV(i) = wellRho(i)*wellV(i)
end do

! cavity at well bottom - mixing celll, no flow
wellV   (wellBottomIndex) = 0.0d0
wellRhoV(wellBottomIndex) = 0.0d0



! interface velocity
wellDeltaVInt(1) = wellV(1)
do i = 2, numWellZones
!  wellDeltaVInt(i) = 0.5*(wellV(i-1)+wellV(i))
  wellDeltaVInt(i) = wellV(i-1)+wellfactor(i)*(wellV(i)-wellV(i-1))
end do



!choke calculation
i = wellBottomIndex
!  rhoIntLocal = wellRho (i-1)+wellfactor(i)*(wellRho (i)-wellRho(i-1))
  rhoIntLocal = wellRho (i-1)
  r8          = 2.0d0*(wellPres(i-1)-wellPres(i)) &
                   /rhoIntLocal
  if (r8 > 0.0) then
    wellDeltaVInt(i) =  DSqrt(r8) *ChokeEfficiency
  else
    wellDeltaVInt(i) = -DSqrt(-r8)*ChokeEfficiency
  end if




! Initialize salt in well, reinitialize other quantities
!dkr account for psuedo cell at wbi
  do i = wellBottomIndex+1, numWellZones
! RM2019 - changed to be a proper inequality comparison for REAL
    if (ABS(wellV(i)) > 0.0) then
      wellSaltMass(i) = drillingRate*bitArea*SaltDensity*wellZoneSize(i)/wellV(i)
    else
      wellSaltMass(i) = 0.0
    endif
    totalSaltInWell = totalSaltInWell+wellSaltMass(i)
    wellSaltVol (i) = wellSaltMass(i)/SaltDensity
    wellMudVolOld   = wellMudVol(i)
    wellMudVol  (i) = wellMudVolOld-wellSaltVol(i)
    wellMudMass (i) = wellMudMass(i)*wellMudVol(i)/wellMudVolOld
    wellRho     (i) = (wellMudMass(i)+wellSaltMass(i))/wellVol(i)
    wellRhoV    (i) = wellRho(i)*wellV(i)
  end do
  i=wellBottomIndex
    wellSaltMass(i) = wellsaltMass(i+1)*wellVol(i)/wellVol(i+1)
    totalSaltInWell = totalSaltInWell+wellSaltMass(i)
    wellSaltVol (i) = wellSaltMass(i)/SaltDensity
    wellMudVolOld   = wellMudVol(i)
    wellMudVol  (i) = wellMudVolOld-wellSaltVol(i)
    wellMudMass (i) = wellMudMass(i)*wellMudVol(i)/wellMudVolOld
    wellRho     (i) = (wellMudMass(i)+wellSaltMass(i))/wellVol(i)
    wellRhoV    (i) = 0.0



! reassign data if running validation test case 2 (air but no salt or mud)
IF(validationTestcase .EQ. 2) THEN
  !JFS3
  i = wellBottomIndex
  gasBaseDensity = 1.0
  wasteDensity   = 1400.0
  wellPres    (i) = repositoryInitialPressure
  wellRho     (i) = gasBaseDensity*wellPres(i)/AtmosphericPressure
  wellGasVol  (i) = wellVol(i)
  wellGasMass (i) = wellRho(i)*wellVol(i)
  totalGasInWell  = wellGasMass(i)
  wellV       (i) = 0.0
  wellSaltMass(i) = 0.0
  totalSaltInWell = 0.0
  wellSaltVol (i) = 0.0
  wellMudVol  (i) = 0.0
  wellMudMass (i) = 0.0
  wellRhoV    (i) = 0.0
!  wellRho     (0) = wellRho(1)

!JFS3 estimate initial gradient for compressible gas
  do i = wellBottomIndex+1, numWellZones
    wellPres    (i) = wellPres(i-1) - wellRho(i-1)*gravity*wellZoneSize(i)
    wellRho     (i) = gasBaseDensity*wellPres(i)/AtmosphericPressure
    wellGasVol  (i) = wellVol(i)
    wellGasMass (i) = wellRho(i)*wellVol(i)
    totalGasInWell = totalGasInWell+wellGasMass(i)
    wellV       (i) = 0.0
    wellSaltMass(i) = 0.0
    wellSaltVol (i) = 0.0
    wellMudVol  (i) = 0.0
    wellMudMass (i) = 0.0
    wellRhoV    (i) = 0.0
  enddo

  wellDeltaVInt(wellBottomIndex)  = 0.0
  wellPres(numWellZones+1) = AtmosphericPressure
  wellRho (numWellZones+1) = gasBaseDensity

elseif(validationTestCase == 5)then
  do i = 1, numWellZones
    wellRho     (i) = initialMudDensity*(1.0d0+WaterCompressibility*(wellPres(i)-AtmosphericPressure))
    wellGasVol  (i) = 0.0
    wellGasMass (i) = 0.0
    wellMudVol  (i) = wellVol(i)
    wellMudMass (i) = wellRho(i)*wellVol(i)
    wellSaltMass(i) = 0.
    totalSaltInWell = 0.
    wellSaltVol (i) = 0.
  end do
  wellPres(numWellZones+1) = AtmosphericPressure
  wellRho (numWellZones+1) = initialMudDensity

endif



! Initialize mud in well
do i = 1, numWellZones
  totalMudInWell = totalMudInWell+wellMudMass(i)
end do



! arbitrary boundary assignments to avoid zero divides
wellMassOrig     (firstWellZone-1) = 1.0
wellGasMassOrig  (firstWellZone-1) = 1.0
wellWasteMassOrig(firstWellZone-1) = 1.0
wellSaltMassOrig (firstWellZone-1) = 1.0
wellDeltaVInt    (firstWellZone)   = 0.0
wellMassOrig     (numWellZones+1)  = 1.0
wellGasMassOrig  (numWellZones+1)  = 1.0
wellWasteMassOrig(numWellZones+1)  = 1.0
wellSaltMassOrig (numWellZones+1)  = 1.0
wellDeltaVInt    (numWellZones+1)  = 0.0
wellRho          (numWellZones+1)  = 1.0
wellRhoV         (numWellZones+1)  = 1.0
wellRho          (0) = wellRho (1)
wellRhoV         (0) = 1.0



! repository cavity
cavityPres          = repositoryInitialPressure
initialCavityGasVol = initialCavityVol*repositoryInitialPorosity
cavityGasMass       = initialCavityGasVol*gasBaseDensity*cavityPres/AtmosphericPressure
cavityWasteVol      = initialCavityVol-initialCavityGasVol
cavityWasteMass     = WasteDensity*cavityWasteVol

deltaGasIntoWell   = 0.0
deltaWasteIntoWell = 0.0
deltaSaltIntoWell  = 0.0
deltaGasFromWaste  = 0.0
volStore           = 0.0
wasteStore         = 0.0
gasStore           = 0.0



! reassign values if running validation test case 2 (coal)
IF(validationTestCase .EQ. 2) THEN
  initialCavityGasVol = initialCavityVol
  cavityGasMass       = initialCavityGasVol*gasBaseDensity*cavityPres/AtmosphericPressure
  cavityWasteVol      = initialCavityVol-initialCavityGasVol
  cavityWasteMass     = WasteDensity*cavityWasteVol

ELSEIF(validationTestCase == 1)THEN
  cavityPres = 0.0

ENDIF

return
end

!-----------------------------------------------------------------------------------------

Subroutine CavityMaxes

Use Globals
Implicit NONE

Real(8) drilledVolMax, cuttingsAreaMax

!pre penetration
initialCavityArea = Pi*0.25*bitDiameter**2
drilledVolMax     = initialCavityArea*repositoryThickness
cuttingsMassMax   = drilledVolMax*(1.0-repositoryInitialPorosity)*WasteDensity
cuttingsAreaMax   = initialCavityArea+Pi*bitDiameter*repositoryThickness

!Original logic - Cavity area = bottom of wellbore circular area
IF(inputCavityRadius <= 0.0)then
  !if penetrated, use wellbore surface area including bottom circle
  if(initialBitAboveRepository < 0.0) Then

    if(initialBitAboveRepository > -0.99*repositoryThickness )then
       initialCavityArea = initialCavityArea &
                        -Pi*bitDiameter*initialBitAboveRepository

   !if near bottom don't include bottom circle
    else
       initialCavityArea = Pi*bitDiameter*repositoryThickness
    endif
  endif

else
!Initial cavity Radius is specified by user
!implemented for TestCase 2 - coalbed methane
!**************************************************************
!IF (validationTestCase .EQ. 2) THEN
!JFS3 use 1.2 ft hole, made cylindrical, for the field case
!  initialCavityRadius = 0.366 !cycle 1
!  initialCavityRadius = 0.425 !cycle 2
!*************************************************************
  IF (geometry == 'S') THEN
    initialCavityArea = 2.0d0*Pi*inputCavityRadius**2
  ELSE
    if(initialBitAboveRepository < 0.0)then
      initialCavityArea = 2.0d0*Pi*inputCavityRadius*(-initialBitAboveRepository)
    else
      initialCavityArea = 2.0d0*Pi*inputCavityRadius*repositoryThickness
    endif
  ENDIF

endif

if (geometry == 'S') then

  initialCavityRadius = DSqrt(initialCavityArea/(2.0d0*Pi))
  initialCavityVol    = (2.0d0/3.0d0)*Pi*initialCavityRadius**3
  cuttingsRadiusMax   = DSqrt(cuttingsAreaMax/(2.0d0*Pi))

else

  initialCavityRadius = initialCavityArea/(2.0d0*Pi*repositoryThickness)
  initialCavityVol  = Pi*repositoryThickness*initialCavityRadius**2
  cuttingsRadiusMax = cuttingsAreaMax/(2.0d0*Pi*repositoryThickness)
end if

equivDrilledCavityRadius = initialCavityRadius

return
end
!-----------------------------------------------------------------------------------------

Subroutine ZoneWellbore

Use Globals
Implicit NONE

Integer i, j, numZones, collarIndex, exitPipeIndex  !apg Integer(4)
Real(8) diffRatio, dzh, Dia

! add exit pipe (50 ft 2-10 in diameter)
if(exitPipeLength > 0.00001)then
  numExitZones = 2*numWellZones4
  exitPipeArea = 0.25*Pi*exitPipeDiameter**2

else
  exitPipeLength = 0.0
  numExitZones   = 0
  exitPipeArea   = pipeAnnulusArea
  exitPipeDiameter = 2.0d0*sqrt(exitPipeArea/Pi)
endif

CALL allocateWellArrays

!-------------------
! Down inside pipe
!-------------------
  if(wellGrowthRate < 1.00001)then
    wellZoneSize(1) = Z1/numwellZones1
  else
    wellZoneSize(1) = Z1*(1.0-wellGrowthRate)/(1.0-wellGrowthRate**numwellzones1)
  endif
  wellZoneSize(0) = wellZoneSize(1)
  wellVol     (0) = pipeArea*wellZoneSize(1)
  wellPos     (0) = -0.5*wellZoneSize(1)
  wellAreaInt (1) = pipeArea
  wellVolInt  (1) = pipeArea*wellZoneSize(1)

!grow zone size first half
do i = 2, numWellZones1
  wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
enddo

!reverse growth rate
wellgrowthrate = 1.0d0/wellGrowthRate
i =  numWellZones1 + 1
wellZoneSize(i) = wellzoneSize(i-1)

!reduce zones size second half
do j = 2, numWellZones1
  i = j + numwellZones1
  wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
enddo

do i = 1, 2*numwellzones1
  dzh = 0.5*(wellZoneSize(i-1)+wellZoneSize(i))
  wellPos     (i) = wellPos(i-1)+dzh
  wellArea    (i) = pipeArea
  wellVol     (i) = pipeArea*wellZoneSize(i)
  hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*pipeInsideDiameter)
enddo

!reverse growth rate
wellgrowthrate = 1.0d0/wellGrowthRate

!-------------------
!wellBottomIndex
!-------------------
i = wellBottomIndex
!wellZoneSize(i) = 0.5*bitDiameter
wellZoneSize(i) = wellzoneSize(i-1)
dzh = 0.5*(wellzoneSize(i) + wellZoneSize(i-1))
wellPos     (i) = wellPos(i-1)+dzh
!wellArea    (i) = 0.5*(pipeArea+collarAnnulusArea)
wellArea    (i) = collarAnnulusArea
wellVol     (i) = wellArea(i)*wellZoneSize(i)
!hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*bitDiameter)
hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*(bitDiameter+collarDiameter))

numzones = wellbottomIndex
i        = numzones + 1

!-------------------
!collar annulus
!-------------------

if(numwellZones2 > 0)then
  if(wellGrowthRate < 1.00001)then
    wellZoneSize(i) = Z2/numwellZones2
  else
    wellZoneSize(i) = Z2*(1.0-wellGrowthRate)/(1.0-wellGrowthRate**numwellzones2)
  endif

!grow zone size first half
  do j = 2, numWellZones2
    i = numzones+j
    wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
  enddo

!reverse growth rate
  wellgrowthrate = 1.0d0/wellGrowthRate
  numzones = numzones + numWellZones2
  i = numzones+1
  wellZoneSize(i) = wellzoneSize(i-1)

!reduce zones size second half
  do j = 2, numWellZones2
    i = numZones + j
    wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
  enddo

  numzones = numzones - numWellZones2
  do j = 1, 2*numwellzones2
    i = numzones+j
    dzh = 0.5*(wellZoneSize(i)+wellZoneSize(i-1))
    wellPos (i) = wellPos(i-1)+dzh
    wellArea(i) = collarAnnulusArea
    wellVol (i) = wellArea(i)*wellZoneSize(i)
    hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*(bitDiameter+collarDiameter))
  enddo

!reverse growth rate
  wellgrowthrate = 1.0d0/wellGrowthRate
  numzones = numzones + 2*numwellZones2
  i = numzones +1
endif
  collarIndex = i


!-------------------
! above collar
!-------------------
if(wellGrowthRate < 1.00001)then
  wellZoneSize(i) = Z3/numwellZones3
else
  wellZoneSize(i) = Z3*(1.0-wellGrowthRate)/(1.0-wellGrowthRate**numwellzones3)
endif

!grow zone size first half
do j = 2, numWellZones3
  i = numzones+j
  wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
enddo

!reverse growth rate
wellgrowthrate = 1.0d0/wellGrowthRate
numzones = numzones + numWellZones3
i = numzones+1
wellZoneSize(i) = wellzoneSize(i-1)


!reduce zones size second half
do j = 2, numWellZones3
  i = numZones + j
  wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
enddo

numzones = numzones - numWellZones3
do j = 1, 2*numwellzones3
  i = numzones+j
  dzh = 0.5*(wellZoneSize(i)+wellZoneSize(i-1))
  wellPos (i) = wellPos(i-1)+dzh
  wellArea(i) = pipeAnnulusArea
  wellVol (i) = wellArea(i)*wellZoneSize(i)
  hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*(bitDiameter+pipeDiameter))
enddo

wellgrowthrate = 1.0d0/wellGrowthRate
numzones = numzones + 2*numwellZones3

!-------------------
!exit pipe (if there is one)
!-------------------
if(exitPipeLength > 0.0)then
  i = numzones +1
  exitPipeIndex = i
  if(wellGrowthRate < 1.00001)then
    wellZoneSize(i) = Z4/numwellZones4
  else
    wellZoneSize(i) = Z4*(1.0-wellGrowthRate)/(1.0-wellGrowthRate**numwellzones4)
  endif

  !grow zone size first half
  do j = 2, numWellZones4
    i = numzones+j
    ! Position
    wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
  enddo

  !reverse growth rate
  wellgrowthrate = 1.0d0/wellGrowthRate

  numzones = numzones + numWellZones4
  i = numzones+1
  wellZoneSize(i) = wellzoneSize(i-1)

  !reduce zone size second half
  do j = 2, numWellZones4
    i = numZones + j
    wellZoneSize(i) = wellzoneSize(i-1)*wellGrowthRate
  enddo

  numzones = numzones - numWellZones4
  do j = 1, 2*numwellzones4
    i = numzones+j
    dzh = 0.5*(wellZoneSize(i)+wellZoneSize(i-1))
    wellPos (i) = wellPos(i-1)+dzh
    wellArea(i) = exitPipeArea
    wellVol (i) = wellArea(i)*wellZoneSize(i)
    hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*exitPipeDiameter)
  enddo
endif

wellPos     (numWellZones+1) = wellPos     (numWellZones)+ wellZoneSize(numWellZones)
!dkr changed for QA0110
wellZoneSize(numWellZones+1) = wellZoneSize(numWellZones)



! smooth transition at collar
if(numwellzones2 > 0)Then
  wellArea(collarIndex-1) = (2d0/3d0)*collarAnnulusArea+(1d0/3d0)*pipeAnnulusArea  !dkr v1.22, was .3333/.6667
  wellArea(collarIndex)   = (1d0/3d0)*collarAnnulusArea+(2d0/3d0)*pipeAnnulusArea  !dkr v1.22, was .6667/.3333
  wellVol (collarIndex-1) = wellArea(collarIndex-1)*wellZoneSize(collarIndex-1)
  wellVol (collarIndex)   = wellArea(collarIndex)  *wellZoneSize(collarIndex)
  Dia = 2.0d0*SQRT((0.5d0*bitDiameter)**2 -wellArea(collarIndex-1)/Pi)
  hydraulicDia(collarIndex-1) = 4.0d0*wellArea(collarIndex-1)/(Pi*(Dia+bitDiameter))
  Dia = 2.0d0*SQRT((0.5d0*bitDiameter)**2 -wellArea(collarIndex)/Pi)
  hydraulicDia(collarIndex)   = 4.0d0*wellArea(collarIndex)  /(Pi*(Dia+bitDiameter))
endif

! modified for Version QA0110 for restart with cylindrical geometry
IF(validationTestCase .EQ. 2)THEN !coal
  wellVol (wellBottomIndex) = initialCavityVol
  wellarea(wellBottomIndex) = initialCavityArea

!JFS3 don't need this, taken care of by pipe annulus and collar annulus assignments
!  wellArea(1) = wellVol(1)/wellZoneSize(1) !??????????????????
!  do i = 2, numWellZones
!    wellArea(i) = pipeArea
!    wellVol (i) = wellArea(i)*wellZoneSize(i)
!    hydraulicDia(i) = 4.0d0*wellArea(i)/(Pi*1.27e-2)
! enddo

ENDIF

wellzoneSize(0) = wellzoneSize(1)
wellArea    (0) = wellArea(1)
wellVol     (0) = wellVol (1)
do i = 1, numWellZones
  wellFactor (i) = wellzoneSize(i-1) / (wellZoneSize(i-1)  +wellZoneSize(i))
  wellAreaInt(i) = wellArea(i-1)+wellfactor(i)*(wellArea(i)-wellArea(i-1))
  wellVolInt (i) = wellVol (i-1)+wellfactor(i)*(wellVol (i)-wellVol (i-1))
end do
wellRho(0) = wellRho(1)

! wellVolInt only used in artificial mass and momentum diffusion
! setting to zero at interfaces limits it to interior cells only
wellVolInt(wellBottomIndex)   = 0.
wellVolInt(wellBottomIndex+1) = 0.
wellVolInt(1)                 = 0.
wellVolInt(numWellZones+1)    = 0.



! override for choking at bit
!IF(validationtestCase .EQ. 2)THEN !no bit choking for coal
!  wellVolInt (wellbottomIndex+1) = wellVol (wellbottomIndex+1)
!ELSE
  wellAreaInt(wellBottomIndex)  = bitFlowArea
!ENDIF



! exit
wellAreaInt(wellBottomIndex+1) = wellArea(wellbottomIndex+1)
wellAreaInt(numWellZones+1) = wellAreaInt(numWellZones)




! numerical diffusion
do i = 1, numWellZones
! RM2019 - changed to be a proper inequality comparison for REAL
  if(ABS(wellVolInt(i)) > 0.0)then
    diffRatio  = wellAreaInt(i)*wellZoneSize(i)/wellVolInt(i)
  else
    diffRatio  = wellArea(i)*wellZoneSize(i)/wellVol(i)
  endif
  momDiff (i) = diffRatio*MomentumDiffusionFactor
  massDiff(i) = diffRatio*MassDiffusionFactor
end do

momDiff (wellBottomIndex)  = 0.0
massDiff(wellBottomIndex)  = 0.0

return
end

!-----------------------------------------------------------------------------------------

Subroutine SetupRepository

Use Globals
Implicit NONE

Integer i, nrz  !apg Integer(4)
Real(8) reposRadius1a, rad, dr, reposRadius2, minNumCells

reposRadius2= 1.0
minNumCells = minNumLt

! Geometry
if (geometry == 'S') then
  geomExponent = 3
else
  geomExponent = 2
end if

numReposZones = repositoryOuterRadius/initialReposZoneSize + 1 !trz
reposZoneSize = (repositoryOuterRadius-initialCavityRadius)/(numReposZones - 1) !trz

!estimate number zones with growth rate for dynamic memory allocation
IF( reposRadius1 < (repositoryOuterRadius-initialCavityRadius) .and. GrowthRate > 1.00001) then
    call QAABORT ('Growth Rate should be 1.0 ')
!dkr 7/23 use brute force
  nrz = (reposRadius1)/reposZoneSize+1
  dr  = reposZoneSize
  rad = reposRadius1
  do while (rad < repositoryOuterRadius)
   nrz = nrz + 1
   rad = rad + dr
   if (rad > reposRadius2 .or. (rad <reposRadius2 .and. dr < Lt/minNumCells) )Then
     dr = dr * growthrate
   endif
  enddo
  numReposZones = nrz+10
ELSE
  growthRate = 1.0d0
ENDIF
 
CALL allocateReposArrays


fluidizationTime(0) = 0.0

!JFS3
fluidizationWaitZone = 0
surfaceFailureAllowed = .True.

reposRadius (0) = initialCavityRadius - 0.5*ReposZoneSize
reposRadius (1) = initialCavityRadius + 0.5*ReposZoneSize
reposRadiusH(0) = initialCavityRadius - ReposZoneSize
reposRadiusH(1) = initialCavityRadius
reposDR (1)     = ReposZoneSize
reposDRH(1)     = ReposZoneSize
reposRadius1a   = reposRadius1+initialCavityRadius


i = 1
do while (reposRadiusH(i) < repositoryOuterRadius-0.00001)
!dkr
  i = i + 1
  if(i > numReposZones)then
    write(6,*)
    write(diagnosticFileId,*)' Error estimating number of repository zones.'
    write(diagnosticFileId,*)' Check input parameters for consistancy.'
    call QAABORT (' numReposZones exceeded')  !apg was STOP
  endif

!trz  if(reposRadius(i-1) > reposRadius1a)then
!trz    reposDR (i) = reposDR(i-1)*GrowthRate
!trz    if(reposRadius(i-1) < reposRadius2) Then
!trz      reposDR(i)= MIN(LT/minNumCells, reposDR (i))
!trz    endif
!trz
!trz  else
   reposDR (i) = reposDR(i-1)
!trz  endif
  reposDRH(i) = 0.5*(reposDR(i) + reposDR(i-1))

  reposRadius (i) = reposRadius(i-1) + reposDRH(i)
  reposRadiusH(i) = reposRadius(i) - 0.5*reposDR(i)

end do
numReposZones = i-1

do i = 1 , numReposZones

    reposFactor(i) = reposDR(i-1) /(reposDR(i-1)+reposDR(i))

    if (geometry == 'S') then
      reposVol(i) = (2d0/3d0)*Pi*(reposRadiusH(i+1)**3 &
                                  -reposRadiusH(i)**3)  !dkr v1.22, was .666667
    else
      reposVol(i) = Pi*repositoryThickness*(reposRadiusH(i+1)**2 &
                                           -reposRadiusH(i)**2)
    end if
enddo

reposFactor(1) = 0.5
reposFactor(numReposZones) = 0.5

! Initialize other array variables (except stresses, which are done last)

do i = 1, numReposZones
  reposPres   (i) = repositoryInitialPressure
  porosity    (i) = repositoryInitialPorosity
  invPorosity (i) = 1.0d0/porosity(i)
  permeability(i) = repositoryInitialPerm
  tensileFailureStarted  (i) = .FALSE.
  tensileFailureCompleted(i) = .FALSE.
  shearFailed            (i) = .FALSE.
  fluidizationStarted    (i) = .FALSE.
  fluidizationCompleted  (i) = .FALSE.
  poreVelocity           (i) = 0.0

  !JFS3
  superficialVelocity    (i) = 0.0

  fractionTensileFailed  (i) = 0.0
  fractionFluidized      (i) = 0.0
  tensileFailureTime     (i) = 2.0*Pi*reposRadius(i)/TensileVelocity
  fluidStartTime         (i) = 0.0
  fluidStopTime          (i) = 0.0
  psi                    (i) = 0.0
  drillingFailure        (i) = 1.0
end do

! Initialize failure-related variables.
maxTensileFailedIndex = 0
maxShearFailedIndex   = 0
! 2026 BC-revision: firstIntactZone (renamed from prior usage) is now the
! cavity wall / first non-fluidized cell. firstIntactZone (NEW) is the inner
! edge of the intact-waste region used by CalculateWasteStresses; equals
! firstFailedZone until tensile failure has completed in at least one cell.
firstFailedZone       = 1
firstIntactZone       = 1
fluidizationVelocity  = 1.0 !(arbitrary)

! initalize boundary pressures
reposPres(0)               = repositoryInitialPressure
reposPres(numReposZones+1) = farFieldPorePressure

! Initialize pseudopressure
do i = 0, numReposZones+1
  psi(i) = reposPres(i)**2/gasViscosity
end do

! Initialize stresses (done last because variable values above are required).
radEffStress(0) = 0.0
tanEffStress(0) = 0.0


! Initial stress state in waste
Call CalculateWasteStresses


! Boundary variables for plotting
poreVelocity(firstFailedZone) = 0.0
wasteBoundaryPoreVelocity     = 0.0


! new formula works in either geometry
initialGasInRepository = 0.0
do i = 1, numReposZones
	initialGasInRepository = initialGasInRepository +&
		(repositoryInitialPressure*reposvol(i)*repositoryInitialPorosity) &
	   *(gasBaseDensity/AtmosphericPressure)
end do

return
end

!-----------------------------------------------------------------------------------------

Subroutine Reset(status)

Use Globals
Implicit NONE

Integer status  !apg Integer(2)

! General
invGasViscosity         = 1.0d0/gasViscosity
invWasteDensity         = 1.0d0/wasteDensity
invSaltDensity          = 1.0d0/saltDensity
invWaterCompressibility = 1.0d0/WaterCompressibility
invInitialMudDensity    = 1.0d0/initialMudDensity
invGasBaseDensity       = 1.0d0/gasBaseDensity
invPi                   = 1.0d0/Pi
invAtmosphericPressure  = 1.0d0/AtmosphericPressure

!if not specified, calulate repository thickness from porosity and initial height
if(repositoryThickness0 <= 0.0) then
  repositoryThickness = RepositoryHeight &
                          *(1.0d0-uncompactedWastePorosity) &
                          /(1.0d0-repositoryInitialPorosity)
else
  repositoryThickness = repositoryThickness0
endif


Call SetupWellbore


invWellZoneSize        = 1.0d0/initialwellZoneSize
invInitialCavityGasVol = 1.0d0/initialCavityGasVol


Call SetupRepository

invReposZoneSize       = 1.0d0/reposZoneSize

! Initialize run management
runTime   = 0.0
runIndex  = 0
deltaTime = 0.0001


Call CalculateTimestep


! left for compatibility with PC version
Call SetupPlots(status)

! Reinitialize run save files
  Call SetSpatialSaveTimes
  Call SetTimeSaveTimes
  spatialSaveIndex = 1
  timeSaveIndex    = 1
  chanSaveIndex    = 1


if(validationTestCase > 0)then
 Call writeTCFileHeaders
endif


if (status .NE. RUNNING) Call CloseRunFiles

return
end

!-----------------------------------------------------------------------------------------

Subroutine SetSpatialSaveTimes

Use Globals
Implicit NONE

Integer i  !apg Integer(2)

IF(validationTestCase .EQ. 1)THEN

!***********************************************************
! Calculate save times for Chan short time solutions
! -D.Lord, 12/16/2002
!***********************************************************

  characteristicTime = (repositoryInitialPorosity * gasViscosity * initialCavityRadius**2)/&
                    (repositoryInitialPerm * farFieldPorePressure)

  chanSaveTime(1) =  0.01 * characteristicTime
  chanSaveTime(2) =  0.10 * characteristicTime
  chanSaveTime(3) =  1.00 * characteristicTime
  chanSaveTime(4) =  10.0 * characteristicTime
  chanSaveTime(5) = 100.0 * characteristicTime


  spatialSaveTime(0) =  0.0
  spatialSaveTime(1) =  0.01 * characteristicTime
  spatialSaveTime(2) =  0.10 * characteristicTime
  spatialSaveTime(3) =  1.00 * characteristicTime
  spatialSaveTime(4) =  10.0 * characteristicTime
  spatialSaveTime(5) = 100.0 * characteristicTime
  do i = 6, 200
    spatialSaveTime(i) = spatialSaveTime(i-1)+50.0*characteristicTime
  end do

ELSEIF(validationTestCase .EQ. 2) THEN
  spatialSaveTime(0) =  0.0
  spatialSaveTime(1) =  0.0001
  spatialSaveTime(2) =  0.0002
  spatialSaveTime(3) =  0.0005
  spatialSaveTime(4) =  0.001
  spatialSaveTime(5) =  0.002
  spatialSaveTime(6) =  0.005
  spatialSaveTime(7) =  0.01
  spatialSaveTime(8) =  0.02
  spatialSaveTime(9) =  0.05
  spatialSaveTime(10) =  0.10
  spatialSaveTime(11) =  0.20
  spatialSaveTime(12) =  0.50
  spatialSaveTime(13) =  1.0
  spatialSaveTime(14) =  2.0
  spatialSaveTime(15) =  5.0
  spatialSaveTime(16) = 10.0
  do i = 17, 200
    spatialSaveTime(i) = spatialSaveTime(i-1)+10.0
  end do

ELSE
  spatialSaveTime(0) =  0.0
  spatialSaveTime(1) =  0.001
  spatialSaveTime(2) =  0.01
  spatialSaveTime(3) =  0.10
  spatialSaveTime(4) =  0.50
  spatialSaveTime(5) =  1.0
  spatialSaveTime(6) =  5.0
  spatialSaveTime(7) = 10.0
  do i = 8, 200
    spatialSaveTime(i) = spatialSaveTime(i-1)+10.0
  end do

ENDIF

!********************************
! test case #4
if (validationTestCase == 4) then
	stressSaveDelta = 10.0
	fluidizationSaveDelta = 1.0 !trz 10.0
endif	
!********************************

return
end

!-----------------------------------------------------------------------------------------
Subroutine SetTimeSaveTimes

Use Globals
Implicit NONE

Integer i  !apg Integer(2)

timeSaveTime(0)  =  0.0
timeSaveTime(1)  =  0.0001
timeSaveTime(2)  =  0.001
timeSaveTime(3)  =  0.01
timeSaveTime(4)  =  0.1
timeSaveTime(5)  =  1.0
timeSaveTime(6)  =  2.0
timeSaveTime(7)  =  4.0
timeSaveTime(8)  =  6.0
timeSaveTime(9)  =  8.0
timeSaveTime(10) =  10.0
do i = 11, 2000
  timeSaveTime(i) = timeSaveTime(i-1)+5.0
end do

return
end

!-----------------------------------------------------------------------------------------

SUBROUTINE allocateReposArrays

USE GLOBALS

INTEGER ialloc_status

if(allocated(tensileFailureStarted)) then
  DEALLOCATE (tensileFailureStarted  ,STAT=ialloc_status)
  DEALLOCATE (tensileFailureCompleted,STAT=ialloc_status)
  DEALLOCATE (shearFailed            ,STAT=ialloc_status)
  DEALLOCATE (fluidizationStarted    ,STAT=ialloc_status)
  DEALLOCATE (fluidizationCompleted  ,STAT=ialloc_status)

  DEALLOCATE (reposRadius ,STAT=ialloc_status)
  DEALLOCATE (reposRadiusH,STAT=ialloc_status)
  DEALLOCATE (reposDR     ,STAT=ialloc_status)
  DEALLOCATE (reposDRH    ,STAT=ialloc_status)
  DEALLOCATE (reposVol    ,STAT=ialloc_status)
  DEALLOCATE (reposGasMass,STAT=ialloc_status)

  DEALLOCATE (reposPres   ,STAT=ialloc_status)
  DEALLOCATE (porosity    ,STAT=ialloc_status)
  DEALLOCATE (permeability,STAT=ialloc_status)
  DEALLOCATE (forchRatio,  STAT=ialloc_status)
  DEALLOCATE (poreVelocity,STAT=ialloc_status)
  DEALLOCATE (psi         ,STAT=ialloc_status)
  DEALLOCATE (dCoeff      ,STAT=ialloc_status)
  DEALLOCATE (radEffStress,STAT=ialloc_status)
  DEALLOCATE (tanEffStress,STAT=ialloc_status)
  DEALLOCATE (shearStress ,STAT=ialloc_status)
  DEALLOCATE (fractionTensileFailed,STAT=ialloc_status)
  DEALLOCATE (tensileFailureTime   ,STAT=ialloc_status)
  DEALLOCATE (fractionFluidized    ,STAT=ialloc_status)
  DEALLOCATE (fluidizationTime     ,STAT=ialloc_status)
  DEALLOCATE (fluidStartTime       ,STAT=ialloc_status)
  DEALLOCATE (fluidStopTime        ,STAT=ialloc_status)
  DEALLOCATE (failStartTime        ,STAT=ialloc_status)
  DEALLOCATE (drillingFailure      ,STAT=ialloc_status)

  DEALLOCATE (invPorosity          ,STAT=ialloc_status)

  DEALLOCATE (radElasticStress,STAT=ialloc_status)
  DEALLOCATE (tanElasticStress,STAT=ialloc_status)
  DEALLOCATE (radSeepageStress,STAT=ialloc_status)
  DEALLOCATE (tanSeepageStress,STAT=ialloc_status)

  DEALLOCATE (aa ,STAT=ialloc_status)
  DEALLOCATE (bb ,STAT=ialloc_status)
  DEALLOCATE (cc ,STAT=ialloc_status)
  DEALLOCATE (rr ,STAT=ialloc_status)
  DEALLOCATE (gam,STAT=ialloc_status)
endif

nrz = numReposZones+10
ALLOCATE (tensileFailureStarted  (0:nrz),STAT=ialloc_status)
ALLOCATE (tensileFailureCompleted(0:nrz),STAT=ialloc_status)
ALLOCATE (shearFailed            (0:nrz),STAT=ialloc_status)
ALLOCATE (fluidizationStarted    (0:nrz),STAT=ialloc_status)
ALLOCATE (fluidizationCompleted  (0:nrz),STAT=ialloc_status)
ALLOCATE (superficialVelocity    (0:nrz),STAT=ialloc_status)

ALLOCATE (reposRadius (0:nrz),STAT=ialloc_status)
ALLOCATE (reposRadiusH(0:nrz),STAT=ialloc_status)
ALLOCATE (reposDR     (0:nrz),STAT=ialloc_status)
ALLOCATE (reposDRH    (0:nrz),STAT=ialloc_status)
ALLOCATE (reposFactor (0:nrz),STAT=ialloc_status)
ALLOCATE (reposVol    (0:nrz),STAT=ialloc_status)
ALLOCATE (reposGasMass(0:nrz),STAT=ialloc_status)
ALLOCATE (pL          (0:nrz),STAT=ialloc_status)
ALLOCATE (reposPres   (0:nrz),STAT=ialloc_status)
ALLOCATE (porosity    (0:nrz),STAT=ialloc_status)
ALLOCATE (permeability(0:nrz),STAT=ialloc_status)
ALLOCATE (forchRatio  (0:nrz),STAT=ialloc_status)
ALLOCATE (poreVelocity(0:nrz),STAT=ialloc_status)

ALLOCATE (psi         (0:nrz),STAT=ialloc_status)
ALLOCATE (dCoeff      (0:nrz),STAT=ialloc_status)
ALLOCATE (radEffStress(0:nrz),STAT=ialloc_status)
ALLOCATE (tanEffStress(0:nrz),STAT=ialloc_status)
ALLOCATE (shearStress (0:nrz),STAT=ialloc_status)

ALLOCATE (fractionTensileFailed(0:nrz),STAT=ialloc_status)
ALLOCATE (tensileFailureTime   (0:nrz),STAT=ialloc_status)
ALLOCATE (fractionFluidized    (0:nrz),STAT=ialloc_status)
ALLOCATE (fluidizationTime     (0:nrz),STAT=ialloc_status)
ALLOCATE (fluidStartTime       (0:nrz),STAT=ialloc_status)
ALLOCATE (fluidStopTime        (0:nrz),STAT=ialloc_status)
ALLOCATE (failStartTime        (0:nrz),STAT=ialloc_status)
ALLOCATE (drillingFailure      (0:nrz),STAT=ialloc_status)

ALLOCATE (invPorosity          (0:nrz),STAT=ialloc_status)

ALLOCATE (shearStrength   (0:nrz),STAT=ialloc_status)
ALLOCATE (meanEffStress   (0:nrz),STAT=ialloc_status)
ALLOCATE (radElasticStress(0:nrz),STAT=ialloc_status)
ALLOCATE (tanElasticStress(0:nrz),STAT=ialloc_status)
ALLOCATE (radSeepageStress(0:nrz),STAT=ialloc_status)
ALLOCATE (tanSeepageStress(0:nrz),STAT=ialloc_status)

ALLOCATE (aa (0:nrz),STAT=ialloc_status)
ALLOCATE (bb (0:nrz),STAT=ialloc_status)
ALLOCATE (cc (0:nrz),STAT=ialloc_status)
ALLOCATE (rr (0:nrz),STAT=ialloc_status)
ALLOCATE (gam(0:nrz),STAT=ialloc_status)

RETURN
END

!-----------------------------------------------------------------------------------------

SUBROUTINE allocateWellArrays

USE GLOBALS

INTEGER ialloc_status

if(allocated(wellPos)) THEN
  DEALLOCATE (wellPos           ,STAT=ialloc_status)
  DEALLOCATE (wellPres          ,STAT=ialloc_status)
  DEALLOCATE (wellRhoV          ,STAT=ialloc_status)
  DEALLOCATE (wellVol           ,STAT=ialloc_status)
  DEALLOCATE (WellArea          ,STAT=ialloc_status)
  DEALLOCATE (wellV             ,STAT=ialloc_status)
  DEALLOCATE (wellRho           ,STAT=ialloc_status)
  DEALLOCATE (wellAreaInt       ,STAT=ialloc_status)
  DEALLOCATE (wellMudMass       ,STAT=ialloc_status)
  DEALLOCATE (wellGasMass       ,STAT=ialloc_status)
  DEALLOCATE (wellWasteMass     ,STAT=ialloc_status)
  DEALLOCATE (wellDeltaVInt     ,STAT=ialloc_status)
  DEALLOCATE (wellSaltMass      ,STAT=ialloc_status)
  DEALLOCATE (wellMudVol        ,STAT=ialloc_status)
  DEALLOCATE (wellGasVol        ,STAT=ialloc_status)
  DEALLOCATE (wellWasteVol      ,STAT=ialloc_status)
  DEALLOCATE (wellSaltVol       ,STAT=ialloc_status)
  DEALLOCATE (wellMassOrig      ,STAT=ialloc_status)
  DEALLOCATE (wellPresOrig      ,STAT=ialloc_status)
  DEALLOCATE (wellMudMassOrig   ,STAT=ialloc_status)
  DEALLOCATE (wellGasMassOrig   ,STAT=ialloc_status)
  DEALLOCATE (wellWasteMassOrig ,STAT=ialloc_status)
  DEALLOCATE (wellSaltMassOrig  ,STAT=ialloc_status)
  DEALLOCATE (wellGRhoStar      ,STAT=ialloc_status)
  DEALLOCATE (wellGRhoVStar     ,STAT=ialloc_status)
  DEALLOCATE (wellGRhoT         ,STAT=ialloc_status)
  DEALLOCATE (wellGRhoVT        ,STAT=ialloc_status)
  DEALLOCATE (wellVolInt        ,STAT=ialloc_status)
  DEALLOCATE (momDiff           ,STAT=ialloc_status)
  DEALLOCATE (massDiff          ,STAT=ialloc_status)
  DEALLOCATE (hydraulicDia      ,STAT=ialloc_status)
endif

nwz = numWellZones+10
ALLOCATE (wellZoneSize      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellPos           (0:nwz),STAT=ialloc_status)
ALLOCATE (wellPres          (0:nwz),STAT=ialloc_status)
ALLOCATE (wellRhoV          (0:nwz),STAT=ialloc_status)
ALLOCATE (wellVol           (0:nwz),STAT=ialloc_status)
ALLOCATE (WellArea          (0:nwz),STAT=ialloc_status)
ALLOCATE (wellV             (0:nwz),STAT=ialloc_status)
ALLOCATE (wellRho           (0:nwz),STAT=ialloc_status)
ALLOCATE (wellAreaInt       (0:nwz),STAT=ialloc_status)
ALLOCATE (wellMudMass       (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGasMass       (0:nwz),STAT=ialloc_status)
ALLOCATE (wellWasteMass     (0:nwz),STAT=ialloc_status)
ALLOCATE (wellDeltaVInt     (0:nwz),STAT=ialloc_status)
ALLOCATE (wellSaltMass      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellMudVol        (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGasVol        (0:nwz),STAT=ialloc_status)
ALLOCATE (wellWasteVol      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellSaltVol       (0:nwz),STAT=ialloc_status)
ALLOCATE (wellMassOrig      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellPresOrig      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellMudMassOrig   (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGasMassOrig   (0:nwz),STAT=ialloc_status)
ALLOCATE (wellWasteMassOrig (0:nwz),STAT=ialloc_status)
ALLOCATE (wellSaltMassOrig  (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGRhoStar      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGRhoVStar     (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGRhoT         (0:nwz),STAT=ialloc_status)
ALLOCATE (wellGRhoVT        (0:nwz),STAT=ialloc_status)
ALLOCATE (wellVolInt        (0:nwz),STAT=ialloc_status)
ALLOCATE (momDiff           (0:nwz),STAT=ialloc_status)
ALLOCATE (massDiff          (0:nwz),STAT=ialloc_status)
ALLOCATE (hydraulicDia      (0:nwz),STAT=ialloc_status)
ALLOCATE (wellFactor        (0:nwz),STAT=ialloc_status)

RETURN
END
!-----------------------------------------------------------------------------------------



