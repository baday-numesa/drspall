!-----------------------------------------------------------------------------------------
!
! File WellboreCalc.f90, contains routines for the wellbore flow calculation
!
!-----------------------------------------------------------------------------------------

Subroutine CalculateWellboreFlow

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) massTot, invMassRatio, r8, rhoIntLocal

! velocity and pressure
Call WellVelocityAndPressure



! Save certain original values for later use
do i = firstWellZone, numWellZones
  wellMassOrig(i) = wellRho(i)*wellVol(i)
  wellPresOrig(i) = wellPres(i)
end do
wellPresOrig(numWellZones+1) = atmosphericPressure

! correct masses for any cumulative roundoff error
do i = firstWellZone, numWellZones
  massTot = wellMudMass(i)+wellGasMass(i)+wellWasteMass(i)+wellSaltMass(i)

!DKR
  IF(massTot .lt. 0.0)then
   write(6,*)
   write(diagnosticFIleId,*)' i=',i,' massTot=',massTot
   call QAABORT ('massTot <= 0.0')  !apg was STOP
  ENDIF

  wellMudMassOrig  (i) = wellMudMass  (i)
  wellGasMassOrig  (i) = wellGasMass  (i)
  wellWasteMassOrig(i) = wellWasteMass(i)
  wellSaltMassOrig (i) = wellSaltMass (i)
end do
wellMudMassOrig  (0) = 0.0



! propagate drillbit
if ((drilling) &
  .AND. (runTime < stopDrillingTime) &
  .AND. ((wellV(numWellZones)*wellAreaInt(numWellZones) < stopDrillingExitVolRate) &
  .AND. (-bitAboveRepository < repositoryThickness) &
  .AND. (runTime < TimeToFloor) &
  .OR. (runTime < 5.0))) then
  bitAboveRepository = bitAboveRepository-drillingRate*deltaTime
else
  drilling =.FALSE.
end if

if ((bitAboveRepository <= 0.0) .AND. (.NOT.repositoryPenetrated)) then
  repositoryPenetrated = .TRUE.
  firstPenetration     = .TRUE.
end if

! prepare for convection, no expansion yet
do i = firstWellZone+1, numWellZones
  ! interface velocity account for variable zone size
  !wellDeltaVInt(i) = 0.5d0*(wellV(i-1)+wellV(i))
  wellDeltaVInt(i) = wellV(i-1) + wellfactor(i)*(wellV(i)-wellV(i-1))
end do


         !-----------------
         !cavity interface
         !-----------------

i = wellBottomIndex
wellDeltaVInt(i+1) = wellV(i+1)


! use Bernoulli choke flow for bit
!  rhoIntLocal = 0.5d0*(wellRho (i-1)+wellRho (i))
  rhoIntLocal = wellRho (i-1)+wellFactor(i)*(wellRho(i)-wellRho(i-1))
  rhoIntLocal = wellRho(i-1)
  r8          = 2.0d0*(wellPres(i-1)-wellPres(i)) &
                   /rhoIntLocal
  if (r8 > 0.0) then
    wellDeltaVInt(i) = DSqrt(r8)*ChokeEfficiency
  else
    wellDeltaVInt(i) = -DSqrt(-r8)*ChokeEfficiency
  end if


IF(validationTestCase == 2 .or. firstWellZone == i)THEN
  wellDeltaVInt(i) = 0.0

ENDIF

! limit to sonic
if (wellDeltaVInt(i) > 1000.0) then
    wellDeltaVInt(i) = 1000.0
else if (wellDeltaVInt(i) < -1000.0) then
         wellDeltaVInt(i) = -1000.0
end if



! mass convection
Call WellMassConvection

! momentum convection
Call WellMomentumConvection

! mass sources
Call WellMassSources

! momentum sources
Call WellMomentumSources

! mass diffusion
Call WellMassDiffusion

! momentum diffusion
Call WellMomentumDiffusion


!gas mass balance term
sumWellGasMass   = 0.0
sumWellWasteMass = 0.0
do i = firstWellZone, numWellZones
  sumWellGasMass   = sumWellGasMass   + wellGasMass  (i)
  sumWellWasteMass = sumWellWasteMass + wellWasteMass(i)
enddo

return
end

!-----------------------------------------------------------------------------------------

Subroutine WellVelocityAndPressure

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) a, b, c, wellAvailVol, mudRho, temp, volGas0, gasPress, gasRho, gasVF

do i = firstWellZone, numWellZones

  ! velocity
  wellV(i) = wellRhoV(i)/wellRho(i)

  ! pressure
  wellWasteVol(i) = wellWasteMass(i)*invWasteDensity
  wellSaltVol (i) = wellSaltMass (i)*invSaltDensity
  wellAvailVol    = wellVol(i)-wellWasteVol(i)-wellSaltVol(i)

! RM2019 - changed to be a proper inequality comparison for REAL
  if (ABS(wellGasMass(i)) > 0.0) then
  ! gas present

    temp    = AtmosphericPressure*WaterCompressibility
    VolGas0 = wellGasMass(i)*invGasBaseDensity
    a = 1.0d0-temp
    b = wellMudMass(i)*invInitialMudDensity &
           +temp*volGas0 -a*wellAvailVol
    c = -temp*volGas0*wellAvailVol

!DKR
    temp = b**2-4.0d0*a*c
    IF(temp < 0) THEN
      call QAABORT ('SQRT(-x) wellGasVol')  !apg was STOP
    ENDIF

    wellGasVol(i) = 0.5d0*(-b+DSqrt(temp))/a

!DKR
    if(wellgasVol(i) < 0.0)Then
      write(diagnosticFileID,'(a,1pe11.4)')' WARN: wellgasvol=', wellgasVol(i)
      write(diagnosticFileID,'(a,i5)')' i= ', i
      write(diagnosticFileID,'(a,1pe11.4)')' wellmudmass= ', wellmudmass(i)
      write(diagnosticFileID,'(a,1pe11.4)')' wellgasmass= ', wellgasmass(i)
      write(diagnosticFileID,'(a,1pe11.4)') ' a= ', a
      write(diagnosticFileID,'(a,1pe11.4)') ' b= ', b
      write(diagnosticFileID,'(a,1pe11.4)') ' c= ', c
      write(diagnosticFileID,'(a,1pe11.4)') ' sqrt=', DSqrt(temp)

       wellGasVol (i) = 0.0d0
       wellGasMass(i) = 0.0d0
    endif

    wellMudVol(i) = wellAvailVol-wellGasVol(i)

  else
    ! no gas present
    wellMudVol(i) = wellAvailVol

  end if


  gasVF = wellgasVol(i)/wellvol(i)

! if gas vol. frac. is > .95 use gas EOS else use mud
  IF(validationTestCase .EQ. 2 .or. gasVF > 0.95) THEN

!    wellPres(i) = (wellRho(i)/gasBaseDensity)*AtmosphericPressure
    if(wellgasVol(i) > 0.0)then
      gasRho      = wellGasMass(i)/wellGasVol(i)
      wellPres(i) = (gasRho/gasBaseDensity)*AtmosphericPressure
    else
      wellPres(i) = 1.0d5
    endif


  ELSE
    mudRho = wellMudMass(i)/wellMudVol(i)
!   wellPres(i) = AtmosphericPressure &
!                +DLog(mudRho*invInitialMudDensity)*invWaterCompressibility
    wellPres(i) = AtmosphericPressure &
                 +(mudRho*invInitialMudDensity-1.0d0)*invWaterCompressibility

!dkr added gas pressure check
    if(wellGasVol(i) > 0.0)then
      gasRho   = wellGasMass(i)/wellGasVol(i)
      gasPress = AtmosphericPressure*invGasBaseDensity* gasRho
    endif
  ENDIF

! DKR: This causes problem zone consistancy - particularly at small DZ
!  if (wellPres(i) < 0.5*AtmosphericPressure) then
!    wellPres(i) = AtmosphericPressure
!  end if

end do

return
end

!-----------------------------------------------------------------------------------------

Subroutine WellMassConvection

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) valRightEdge, valLeftEdge, deltaMudMass, deltaGasMass, &
  deltaWasteMass, deltaSaltMass, r8, rhoL, rhoR


! general convection
do i = firstWellZone, numWellZones

!DKR - upwind
  valLeftEdge  = deltaTime*wellAreaInt(i)  *wellDeltaVInt(i)
  valRightEdge = deltaTime*wellAreaInt(i+1)*wellDeltaVInt(i+1)


  deltaMudMass   = 0.0
  deltaGasMass   = 0.0
  deltaWasteMass = 0.0
  deltaSaltMass  = 0.0

  if (valLeftEdge > 0.0) then
    r8 = 1.0d0/wellMassOrig(i-1)
    valLeftedge    = valLeftedge*wellRho(i-1)
    deltaMudMass   = valLeftEdge*wellMudMassOrig  (i-1)*r8
    deltaGasMass   = valLeftEdge*wellGasMassOrig  (i-1)*r8
    deltaWasteMass = valLeftEdge*wellWasteMassOrig(i-1)*r8
    deltaSaltMass  = valLeftEdge*wellSaltMassOrig (i-1)*r8
  else
    r8 = 1.0d0/wellMassOrig(i)
    valLeftedge    = valLeftedge*wellRho(i)
    deltaMudMass   = valLeftEdge*wellMudMassOrig  (i)*r8
    deltaGasMass   = valLeftEdge*wellGasMassOrig  (i)*r8
    deltaWasteMass = valLeftEdge*wellWasteMassOrig(i)*r8
    deltaSaltMass  = valLeftEdge*wellSaltMassOrig (i)*r8
  end if

  if (valRightEdge < 0.0) then
    r8 = 1.0d0/wellMassOrig(i+1)
    valRightedge   = valRightedge*wellRho(i+1)
    deltaMudMass   = deltaMudMass  -valRightEdge*wellMudMassOrig  (i+1)*r8
    deltaGasMass   = deltaGasMass  -valRightEdge*wellGasMassOrig  (i+1)*r8
    deltaWasteMass = deltaWasteMass-valRightEdge*wellWasteMassOrig(i+1)*r8
    deltaSaltMass  = deltaSaltMass -valRightEdge*wellSaltMassOrig (i+1)*r8
  else
    r8 = 1.0d0/wellMassOrig(i)
    valRightedge   = valRightedge*wellRho(i)
    deltaMudMass   = deltaMudMass  -valRightEdge*wellMudMassOrig  (i)*r8
    deltaGasMass   = deltaGasMass  -valRightEdge*wellGasMassOrig  (i)*r8
    deltaWasteMass = deltaWasteMass-valRightEdge*wellWasteMassOrig(i)*r8
    deltaSaltMass  = deltaSaltMass -valRightEdge*wellSaltMassOrig (i)*r8
  end if

  wellGRhoStar (i) = wellMassOrig     (i)-valRightEdge +valLeftEdge
  wellMudMass  (i) = wellMudMassOrig  (i)+deltaMudMass
  wellGasMass  (i) = wellGasMassOrig  (i)+deltaGasMass
  wellWasteMass(i) = wellWasteMassOrig(i)+deltaWasteMass
  wellSaltMass (i) = wellSaltMassOrig (i)+deltaSaltMass

  if(validationTestCase .EQ. 2)wellMudMass (i) = 0.0

  if (wellGasMass  (i) < 0.0) wellGasMass  (i) = 0.0
  if (wellWasteMass(i) < 0.0) wellWasteMass(i) = 0.0
  if (wellSaltMass (i) < 0.0) wellSaltMass (i) = 0.0

end do

return
end

!-----------------------------------------------------------------------------------------

Subroutine WellMomentumConvection

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) rhoVInt(0:2000),valRightEdge, valLeftEdge


! general convection
do i = firstWellZone, numWellZones
!dkr - upwind

  if(wellDeltaVInt(i) > 0.0 .and. i > firstWellZone) Then
    valLeftEdge = deltaTime*wellAreaInt(i)  *wellDeltaVInt(i)*wellRhoV(i-1)
  else
    valLeftEdge = deltaTime*wellAreaInt(i)  *wellDeltaVInt(i)*wellRhoV(i)
  endif

  if(wellDeltaVInt(i+1) < 0.0 .and. i < numwellZones) Then
    valRightEdge = deltaTime*wellAreaInt(i+1)*wellDeltaVInt(i+1)*wellRhoV(i+1)
  else
    valRightEdge = deltaTime*wellAreaInt(i+1)*wellDeltaVInt(i+1)*wellRhoV(i)
  endif

  wellGRhoVStar(i) = wellVol(i)*wellRhoV(i) -valRightEdge +valLeftEdge
end do

wellGRhoVStar(wellbottomIndex) = 0.0d0


return
end

!-----------------------------------------------------------------------------------------

Subroutine WellMassSources

Use Globals
Implicit NONE

Integer i, firstUncutZone  !apg Integer(4)
Real(8) pumpedMudMass, ejectedTotalMass, ejectedMassRatio, deltaMudEjected, deltaGasEjected, &
  deltaWasteEjected, deltaSaltEjected, lenDam, lenDRZ, lenTot, effPerm,  &
  deltaGasFromFluidization, deltaWasteFromFluidization, equivDrilledCavityArea, a, b, c, &
  deltaGasFromDrilling, deltaWasteFromDrilling, deltaSaltFromDrilling, curDensity, &
  deltaVolFromFluidization, deltaVolFromDrilling, deltaVol, deltaGas, deltaWaste, &
  cavityCharTime, storeLostRatio, temp, oldPres, charVel, mudRho, curGasDen, &
  forchNumber, DvolStore, DGasStore, DWasteStore

Data DvolStore, DGasStore, DWasteStore/0.,0.,0./
Data oldPres / -1.0/
Data firstUncutZone /1/



! initialize
do i = firstWellZone, numWellZones
  wellGRhoT(i) = wellGRhoStar(i)
end do

deltaVolFromFluidization   = 0.0
deltaGasFromFluidization   = 0.0
deltaWasteFromFluidization = 0.0
deltaVolFromDrilling       = 0.0
deltaGasFromDrilling       = 0.0
deltaWasteFromDrilling     = 0.0



!---------------------
! mass into entry zone
!---------------------

if ((pumping) &
  .AND. ((wellV(numWellZones)*wellAreaInt(numWellZones) < stopPumpingExitVolRate) &
  .OR. (runTime < 5.0))) then


  if (firstWellZone == 1)then
    !dkr - limit pump pressure
    if(wellPres(1) > oldPres .and. wellPres(1) > maxPumpPressure) then
      oldPres     = wellPres(1)
      mudPumpRate = mudPumpRate*0.99

    elseif( oldPres > 0.0 .and. wellPres(1) < maxPumpPressure &
                          .and. wellPres(1) < oldPres) then
      oldPres     = wellPres(1)
      mudPumpRate = min(mudpumpRate*1.01, initialPumpRate)
    endif
    mudRho = wellRho(1)

  elseif(wellMudVol(firstWellZone) > 0.0)then
    mudRho = wellMudMass(firstWellZone)/wellMudVol(firstWellZone)

  else
    mudRho = initialMudDensity

  endif

    pumpedMudMass    = mudPumpRate*deltaTime*mudRho
    wellGRhoT  (firstWellZone) = wellGRhoT  (firstWellZone)+pumpedMudMass
    wellMudMass(firstWellZone) = wellMudMass(firstWellZone)+pumpedMudMass
    totalMudInWell   = totalMudInWell  +pumpedMudMass
else
  pumping = .FALSE.
end if


!--------------------------------------
! mass into cavity and well bottom zone
!--------------------------------------


! test case 5 bypasses normal mass loading logic at
! bottom of wellbore. Coupling with repository replaced
! with prescribed function.

if(validationTestCase.EQ.5) then

    Call TestCaseFiveMassLoading

else

!-----------------------------------------------------------------------
! 1. repository flow only (no material removal from drilling or fluidization)
!-----------------------------------------------------------------------
! RM2019 - changed to be a proper inequality comparison for REAL
  if ((.NOT.repositoryPenetrated ) .AND. (ABS(bitAboveRepository) > 0.0)) then


! cavity pressure

    ! flow from repository into cavity
    cavityGasMass = cavityGasmass+deltaGasFromWaste

    !************************************************
    ! Setting cavityPres = 0 for Chan validation test
    ! -D.Lord 12/16/2002
    !************************************************
    IF(validationTestCase == 1) THEN
      cavityPres = 0.0
    ELSE
      cavityPres = AtmosphericPressure*cavityGasMass*invGasBaseDensity &
                  *invInitialCavityGasVol
    ENDIF



! bleed from cavity through DRZ into well bottom zone
    if (bitAboveRepository > (dDZThickness+dRZThickness)) then
      lenTot  = dDZThickness +dRZThickness
      effPerm = 0.0
    else if (bitAboveRepository >= dDZThickness) then
      lenDam  = dDZThickness
      lenDRZ  = bitAboveRepository-dDZThickness
      lenTot  = lenDam+lenDRZ
      effPerm = (lenDam+lenDRZ)/(lenDam/dDZPerm+lenDRZ/dRZPerm)
    else
      lenTot = bitAboveRepository
      ! this for avoidance of near-zero divide at penetration time
      if (lenTot < 0.0001) lenTot = 0.0001
      effPerm = dDZPerm
    end if
    deltaGasIntoWell = deltaTime*effPerm*dRZArea*(cavityPres**2- &
    (wellPres(wellBottomIndex))**2)/(2.0d0*gasViscosity*GasConstant*ReposTemp*lenTot)

    !DKR avoids -wellGasMass
    IF( deltaGasIntoWell < 0.0 ) THEN
      deltaGasIntoWell = 0.0
!      write(6,*)' WARNING: corrected -wellGasMass '
    ENDIF

    cavityGasMass                = cavityGasMass   -deltaGasIntoWell
    wellGRhoT  (wellBottomIndex) = wellGRhoT  (wellBottomIndex)+deltaGasIntoWell
    wellGasMass(wellBottomIndex) = wellGasMass(wellBottomIndex)+deltaGasIntoWell
    totalGasInWell               = totalGasInWell  +deltaGasIntoWell
    totalGasInjected             = totalGasInjected+deltaGasIntoWell



    ! salt from drilling
    if ((drilling) &
      .AND. ((wellV(numWellZones)*wellAreaInt(numWellZones) < stopDrillingExitVolRate) &
      .OR. (runTime < 5.0))) then
      deltaSaltIntoWell = SaltDensity*drillingRate*deltaTime*bitArea
    else
      deltaSaltIntoWell = 0.0
    end if

    !dkr - add drilled salt volume to storage cell
    wellVol     (wellBottomIndex) = wellVol     (wellBottomIndex)+deltaSaltIntoWell/saltDensity
    wellGRhoT   (wellBottomIndex) = wellGRhoT   (wellBottomIndex)+deltaSaltIntoWell
    wellSaltMass(wellBottomIndex) = wellSaltMass(wellBottomIndex)+deltaSaltIntoWell
    totalSaltInWell               = totalSaltInWell+deltaSaltIntoWell




  else
    !----------------------------------------
    ! direct flow to well from repository
    !----------------------------------------
    cavityPres                   = wellPres   (wellBottomIndex)
    wellGRhoT  (wellBottomIndex) = wellGRhoT  (wellBottomIndex)+deltaGasFromWaste
    wellGasMass(wellBottomIndex) = wellGasMass(wellBottomIndex)+deltaGasFromWaste
    totalGasInWell               = totalGasInWell  +deltaGasFromWaste
    totalGasInjected             = totalGasInjected+deltaGasFromWaste
  end if




!-----------------------------------------------------------------------
! 2. drilling in repository
!-----------------------------------------------------------------------

  if (repositoryPenetrated  .AND. (-bitAboveRepository < repositoryThickness)) then
    equivDrilledCavityArea = initialCavityArea-bitAboveRepository*Pi*bitDiameter

    if (geometry == 'S') then
      equivDrilledCavityRadius = DSqrt(equivDrilledCavityArea/(2.0d0*Pi))

    else
      equivDrilledCavityRadius = equivDrilledCavityArea &
                                /(2.0d0*Pi*repositoryThickness)

    end if


    !dkr-geometry dependent cuttings volume
    if(equivDrilledCavityRadius >= reposRadiusH(firstUncutZone+1)) THEN
      cutVol = cutVol + reposVol(firstUncutZone)*(1.0-repositoryInitialPorosity)
      firstUncutZone = firstUncutZone+1
    endif

    !dkr - true cuttings volume
    cutTrueVol = (-bitAboveRepository*Pi*(0.5*bitDiameter)**2)*(1.0-repositoryInitialPorosity)



    ! test for possible drilling removal
    if (equivDrilledCavityRadius >= (reposRadiusH(firstIntactZone+1))) then
      tensileFailureStarted  (firstIntactZone) = .TRUE.
      tensileFailureCompleted(firstIntactZone) = .TRUE.
      fractionTensileFailed  (firstIntactZone) = 1.0d0
      shearFailed            (firstIntactZone) = .TRUE.
      fluidizationStarted    (firstIntactZone) = .TRUE.
      fluidizationCompleted  (firstIntactZone) = .TRUE.
      fractionFluidized      (firstIntactZone) = 1.0d0
      drillingfailure        (firstIntactZone) = -1.0d0

! Allow resetting of Lt region near wellbore after new zone
	  surfaceFailureAllowed = .True.

      deltaVol = reposVol(firstIntactZone)
      TotVol   = TotVol + deltaVol*(1.0d0-repositoryInitialPorosity)

      deltaVolFromDrilling = deltaVol
      deltaGasFromDrilling = deltaVol*repositoryInitialPorosity*gasBaseDensity* &
                             (reposPres(firstIntactZone)/AtmosphericPressure)
      deltaWasteFromDrilling = deltaVol*(1.0d0-repositoryInitialPorosity)*WasteDensity
      firstIntactZone = firstIntactZone+1

    end if

  end if



!-----------------------------------------------------------------------
! 3. fluidization
!-----------------------------------------------------------------------

  !calculate Ergun fluidization velocity (use only first intact zone velocity)
  curDensity = gasBaseDensity*(reposPres(firstIntactZone)/AtmosphericPressure)


  a = (1.75d0/(ShapeFactor*repositoryInitialPorosity**3)) &
     *(particleDiameter*curDensity*invGasViscosity)**2

  !DKR: added **2 to shape factor
  b = 150.0d0*((1.0-repositoryInitialPorosity) &
     /((ShapeFactor**2)*repositoryInitialPorosity**3)) &
     *(particleDiameter*curDensity*invGasViscosity)


  c = -(particleDiameter**3)*curDensity*(WasteDensity-curDensity)*Gravity &
      *invGasViscosity**2

  !DKR
  temp = b**2-4.0*a*c
  if(temp < 0.0) then
    call QAABORT ('SQRT(-x) fluidizationVelocity')  !apg was STOP
  endif

 fluidizationVelocity = MAX(1.0d-6, 0.5d0*(-b+DSqrt(temp))/a)


  if ((allowFluidization /= 'N' ) .AND. repositoryPenetrated ) then

   ! allow only those zones beginning at cavity that have completed failure
   ! to begin fluidization
   i = firstIntactZone
   do while (tensileFailureCompleted(i)  .and. i <= numReposZones)

      ! fluidization initialization
      if ((.NOT.fluidizationStarted(i) ) .AND. &
          (superficialVelocity     (i) >   fluidizationVelocity .OR. &
           allowFluidization          == 'A') ) then

        fluidizationStarted(i) = .TRUE.
        fluidStartTime     (i) = runTime
        fluidizationTime   (i) = reposRadius(i)/MAX(MinCharVel, superficialVelocity(i))

		! used by test case #4
		fluidizationSaveTime = runTime
!trz
        if(validationTestCase == 4)then
            CALL WriteToFluidizationTimeValidationFile(i) 
        endif 
!trz
      end if

      ! fluidization in progress - from wellbore outward
      if (fluidizationStarted(i)  .AND. .NOT.fluidizationCompleted(i) ) then


        ! dkr this is New - must stay above fluidization velocity
        if(superficialVelocity(i) > fluidizationVelocity .OR. &
            allowFluidization == 'A') then
          fractionFluidized(i) = fractionFluidized(i)+deltaTime/fluidizationTime(i)
        endif

        ! fluidization completed
        if (fractionFluidized(i) > 1.0 .AND. i == firstIntactZone) then
          fractionFluidized    (i) = 1.0
          fluidizationCompleted(i) = .TRUE.

!JFS3
          if (i == fluidizationWaitZone) then
		    surfaceFailureAllowed = .True.
          endif



          fluidStopTime(i) = runTime
!dkr
          deltaVol = reposVol(i)
          TotVol  = TotVol + deltaVol*(1.0d0-repositoryInitialPorosity)
          splVol  = splVol + deltaVol*(1.0d0-repositoryInitialPorosity)

          deltaVolFromFluidization = deltaVolFromFluidization+deltaVol
          deltaGasFromFluidization = deltaGasFromFluidization+deltaVol*&
             repositoryInitialPorosity*(reposPres(i)/AtmosphericPressure)*gasBaseDensity
          deltaWasteFromFluidization = deltaWasteFromFluidization &
                    +deltaVol*(1.0d0-repositoryInitialPorosity)*WasteDensity

          firstIntactZone = i+1

        elseif (fractionFluidized(i) > 1.0001) then !trz
          fractionFluidized (i) = 1.0001     !trz

        end if

      end if
      i = i + 1
    end do

     ! write to fluidization validation file
     if (ValidationTestCase == 4 .and. runtime >= fluidizationSaveTime) then 
        fluidizationSaveTime = fluidizationSaveTime + fluidizationSaveDelta
        Call WriteToFluidizationValidationFile
     endif


  end if





  !---------------------------------------------------------------------
  ! add to cavity
  !---------------------------------------------------------------------
  if (repositoryPenetrated) then
    deltaVol   = deltaVolFromFluidization  +deltaVolFromDrilling
    deltaGas   = deltaGasFromFluidization  +deltaGasFromDrilling
    deltaWaste = deltaWasteFromFluidization+deltaWasteFromDrilling
    totalGasFromWaste   = totalGasFromWaste   + deltaGas
    totalWasteFromRepos = totalWasteFromRepos + deltaWaste

    if (firstPenetration) then
!      deltaVol   = deltaVol  +initialCavityVol
!      deltaGas   = deltaGas  +cavityGasMass
!      deltaWaste = deltaWaste+cavityWasteMass
      firstPenetration = .FALSE.
    end if

    ! add to store
    volStore   = volStore  +deltaVol
    gasStore   = gasStore  +deltaGas
    wasteStore = wasteStore+deltaWaste

    !remove from store
    charVel = MAX(MinCharVel,superficialVelocity(firstIntactZone))
    cavityCharTime = reposRadius(firstIntactZone)/charVel
    storeLostRatio = deltaTime/cavityCharTime

    dVolStore   = storeLostRatio*volStore
    dGasStore   = storeLostRatio*gasStore
    dWasteStore = storeLostRatio*wasteStore

    wellVol      (wellBottomIndex) = wellVol      (wellBottomIndex)+dVolStore
    wellGasMass  (wellBottomIndex) = wellGasMass  (wellBottomIndex)+dGasStore
    wellWasteMass(wellBottomIndex) = wellWasteMass(wellBottomIndex)+dWasteStore
    wasteInjected =  wasteInjected + dwasteStore

    wellGRhoT(wellBottomIndex) = wellGRhoT(wellBottomIndex)+dGasStore+dwasteStore
    totalGasInWell   = totalGasInWell  +dGasStore
    totalGasInjected = totalGasInjected+dGasStore
    totalWasteInWell = totalWasteInWell+dWasteStore

    volStore   = volStore  - dVolStore
    gasStore   = gasStore  - dGasStore
    wasteStore = wasteStore- dWasteStore

  end if


end if	!end of test case 5

!----------------------
! mass out of exit zone
!----------------------

!dkr: modified to account for exit pipe
!ejectedTotalMass        = wellV(numWellZones)*deltaTime*wellRho(numWellZones)*pipeAnnulusArea
ejectedTotalMass        = wellV(numWellZones)*deltaTime*wellRho(numWellZones)*exitPipeArea
wellGRhoT(numWellZones) = wellGRhoT(numWellZones)-ejectedTotalMass


ejectedMassRatio  = ejectedTotalMass/wellMassOrig(numWellZones)

deltaMudEjected   = wellMudMass  (numWellZones)*ejectedMassRatio
deltaGasEjected   = wellGasMass  (numWellZones)*ejectedMassRatio
deltaWasteEjected = wellWasteMass(numWellZones)*ejectedMassRatio
deltaSaltEjected  = wellSaltMass (numWellZones)*ejectedMassRatio

wellMudMass  (numWellZones) = wellMudMass  (numWellZones)-deltaMudEjected
wellGasMass  (numWellZones) = wellGasMass  (numWellZones)-deltaGasEjected
wellWasteMass(numWellZones) = wellWasteMass(numWellZones)-deltaWasteEjected
wellSaltMass (numWellZones) = wellSaltMass (numWellZones)-deltaSaltEjected

mudMassEjected   = mudMassEjected  +deltaMudEjected
gasMassEjected   = gasMassEjected  +deltaGasEjected
wasteMassEjected = wasteMassEjected+deltaWasteEjected
saltMassEjected  = saltMassEjected +deltaSaltEjected

pitGain = pitGain + deltaMudEjected*invInitialMudDensity - mudPumpRate*deltaTime

totalMudInWell   = totalMudInWell  -deltaMudEjected
totalGasInWell   = totalGasInWell  -deltaGasEjected
totalWasteInWell = totalWasteInWell-deltaWasteEjected
totalSaltInWell  = totalSaltInWell -deltaSaltEjected

return
end

!-----------------------------------------------------------------------------------------

Subroutine WellMomentumSources

Use Globals
Implicit NONE

Integer i, wall, nvc  !apg Integer(4)
Real(8) gravTerm, presTerm, avgAreaInt, weightedViscosity, dissLength, dissTerm, &
  lam, turb, Re, lambda, pumpedMomentum, ejectedMomentum, wellSolidsFraction, accelTerm, &
  invDissLength, injectedGasMomentum, relRough, frictionFactor, &
  ffCole, dpdz, gasMF
Data nvc/0/


! gravity, pressure, and viscous loss
do i =firstWellZone, numWellZones

  !acceleration due to gravity
  gravTerm = wellVol(i)*deltaTime*wellRho(i)*Gravity

  ! inside pipe
  wall = 1

  if (i >numWellZones -numExitZones)Then
    ! horizontal exit pipe
    gravTerm = 0.0

  elseif (i >= (wellBottomIndex)) then
    ! up annulus
    gravTerm = -gravTerm
    wall = 2
  end if

  !acceleration due to pressure gradient
  if (i == 1) then
    ! inlet at pump
    dpdz       = (wellPresOrig(2)-wellPresOrig(1)) &
                 /(wellPos(2) - wellPos(1))

  else if (i == wellBottomIndex-1) then
    !above drill bit
    dpdz = (wellPresOrig(wellBottomIndex-1)-wellPresOrig(wellBottomIndex-2)) &
          /(wellPos(wellBottomIndex-1) - wellPos(wellBottomIndex-2))


  else if (i == wellBottomIndex) then
    !below drill bit
    !dkr - bottom hole cavity
    dpdz     = 0.0
    gravterm = 0.0


  else if (i == wellBottomIndex+1) then
    !above drill bit
    dpdz = (wellPresOrig(wellBottomIndex+2)-wellPresOrig(wellBottomIndex)) &
          /(wellPos(wellBottomIndex+2) - wellPos(wellBottomIndex))


  else if (i == numWellZones) then
    dpdz = (wellPresOrig(numwellzones+1)-wellPresOrig(numWellZones-1)) &
          /(wellPos(numwellzones+1) - wellPos(numwellzones-1))

  else
    ! interior
    dpdz =  (wellPresOrig(i+1)-wellPresOrig(i-1)) &
           /(wellPos(i+1) - wellPos(i-1))


  end if

!dkr correct  area centering
!  presTerm = -0.5d0*deltaTime*avgAreaInt*presDif
  presTerm   = -deltaTime*wellVol(i)*dpdz




!--------------------
!viscous dissipation
!--------------------

!find mass-weighted viscosity of mixture (this is a rough approximation)
!weightedViscosity = (mudViscosity*wellMudMass(i)+gasViscosity*wellGasMass(i)) &
!                   /(wellMudMass(i)+wellGasMass(i))

!**********************************
!DKR: gas viscosity by 95% gas mass fraction
! linear interpolation based on Mass fraction
  gasMF = wellGasMass(i)/(wellGasMass(i)+wellMudMass(i))
  weightedViscosity = mudViscosity + (MIN(gasMF,0.95)/0.95)*(gasViscosity-mudViscosity)
!**********************************

  !modify for solids based on Baree & Conway
  wellSolidsFraction = wellWasteMass(i)*invWasteDensity/wellVol(i)

  if (wellSolidsFraction < mudSolidsMax) then
    weightedViscosity = weightedViscosity &
                   *(1.0d0-wellSolidsFraction/mudSolidsMax)**mudSolidsViscosityExponent

    ! dissipation characteristic length
    dissLength    = DSqrt(wellVol(i)*invPi/WellZoneSize(i))
    invDissLength = 1.0d0/dissLength

    ! combine turbulent and laminar losses by estimation
    !dkr - changed
    turb    = 0.1d0*DSqrt(wallRoughness(wall)/hydraulicDia(i))

    Re    = DAbs(wellV(i))*wellRho(i)*hydraulicDia(i)/weightedViscosity


    !DKR FOR VALIDATION TEST CASE 2 and wellv=0
    IF(Re > 0.0) THEN

       ! Colebrook friciton factor model for pipes
       relRough = wallRoughness(wall)/hydraulicDia(i)
       Re       = DAbs(wellV(i))*wellRho(i)*hydraulicDia(i)/weightedViscosity
       ffCole   = FrictionFactor(relRough,Re)
       lambda   = 0.5*ffCole*wellV(i)**2/hydraulicDia(i)
       dissTerm = deltaTime*wellVol(i)*wellRho(i)*lambda
    ELSE
!      write(6,*)' WARNING: Re < 0.0'
      dissTerm = 0.0
    ENDIF

    ! viscous dissipation always impedes flow
    if (wellGRhoVStar(i) < 0.0) then
      dissTerm = -dissTerm
    end if

    accelTerm = wellGRhoVStar(i)+gravTerm+presTerm

    if (DAbs(dissTerm)>DAbs(accelTerm)) then
      nvc = nvc + 1
      wellGRhoVT(i) = 0.0 !viscous choking
      if(nvc < 100)then
        write(diagnosticFileID,*) ' WARNING: Viscous choking(2) at', i
        write(*,*) ' WARNING: Viscous choking(2) at', i
      endif
    else
      wellGRhoVT(i) = accelTerm-dissTerm
    end if

  else
    nvc = nvc + 1
    wellGRhoVT(i) = 0.0 !viscous choking
    if(nvc < 100)then
      write(diagnosticFileID,*) ' WARNING: Viscous choking(1) at', i
      write(*,*) ' WARNING: Viscous choking(1) at', i
    endif

    wellGRhoVT(i) = 0.0 !viscous choking
  end if

end do


! add mud pump
if(pumping  .AND. firstWellZone == 1)then
  pumpedMomentum = (mudPumpRate**2)*deltaTime*wellRho(1)/pipeArea
  wellGRhoVT(1)  = wellGRhoVT(1)+pumpedMomentum
end if



! add surface ejection
!dkr modified for exit pipe
!ejectedMomentum = (wellV(numWellZones))**2*deltaTime*wellRho(numWellZones)*pipeAnnulusArea
ejectedMomentum = (wellV(numWellZones))**2*deltaTime*wellRho(numWellZones)*exitPipeArea
wellGRhoVT(numWellZones) = wellGRhoVT(numWellZones)-ejectedMomentum


return
end

!-----------------------------------------------------------------------------------------

Subroutine WellMassDiffusion

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) massInFromRight, massOutToRight, massInFromLeft, massOutToLeft, &
  deltaMudMass, deltaGasMass, deltaWasteMass, deltaSaltMass, r81, r82, r83, &
  mdf, mass_old, mass_new, dzh, dmud(2000), dsalt(2000)

i = firstWellZone


do i = firstWellZone, numWellZones
    deltamudmass = 0.
    deltasaltmass = 0.

!    mdf = massDiffusionFactor*deltatime/wellZoneSize(i)**2
    mdf = massDiffusionFactor

    massInFromLeft  = wellRho(i-1)*wellAreaInt(i)  *wellZoneSize(i)*mdf
    massOutToLeft   = wellRho(i)  *wellAreaInt(i)  *wellZoneSize(i)*mdf
    massInFromRight = wellRho(i+1)*wellAreaInt(i+1)*wellZoneSize(i)*mdf
    massOutToRight  = wellRho(i)  *wellAreaInt(i+1)*wellZoneSize(i)*mdf


    r81 = 1.0d0/wellMassOrig(i-1)
    r82 = 1.0d0/wellMassOrig(i)
    r83 = 1.0d0/wellMassOrig(i+1)

  if(i == firstWellZone)then
    massInFromLeft  = 0.
    massOutToLeft   = 0.
    r81 = 1.0d0
  elseif(i == WellbottomIndex)then
    massInFromLeft  = 0.
    massOutToLeft   = 0.
    r81 = 1.0d0
    massInFromRight = 0.
    massOutToRight  = 0.
  elseif(i == WellbottomIndex+1)then
    massInFromLeft  = 0.
    massOutToLeft   = 0.
    r81 = 1.0d0
  elseif(i == numWellZones)then
    massInFromRight = 0.
    massOutToRight  = 0.
    r83 = 1.0d0
  endif


  deltaMudMass = massInFromLeft  *wellMudMassOrig(i-1)*r81 &
                +massInFromRight *wellMudMassOrig(i+1)*r83 &
  -(massOutToLeft+massOutToRight)*wellMudMassOrig(i)  *r82

  deltaSaltMass = massInFromLeft *wellSaltMassOrig(i-1)*r81 &
                 +massInFromRight*wellSaltMassOrig(i+1)*r83 &
  -(massOutToLeft+massOutToRight)*wellSaltMassOrig(i)  *r82

!Turns off salt diffusion
!  deltaSaltMass = 0.0


  mass_old = wellGRhoT(i)
  mass_new = wellGRhoT(i)+deltaMudMass+deltaSaltMass
  dsalt(i) = deltasaltMass
  dmud(i)  = deltamudmass

  wellRho(i) = (wellGRhoT(i)+deltaMudMass+deltaSaltMass)/wellVol(i)

  wellMudMass(i) = wellMudMass(i)+deltaMudMass
  wellSaltMass(i) = wellSaltMass(i)+deltaSaltMass

end do

return
end

!-----------------------------------------------------------------------------------------

Subroutine WellMomentumDiffusion

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) momentumInFromRight, momentumOutToRight, momentumInFromLeft, momentumOutToLeft


do i = firstWellZone, numWellZones

    momentumInFromLeft  = wellRhoV(i-1)*wellAreaInt(i)  *wellZoneSize(i)*momentumDiffusionFactor
    momentumOutToLeft   = wellRhoV(i)  *wellAreaInt(i)  *wellZoneSize(i)*momentumDiffusionFactor
    momentumInFromRight = wellRhoV(i+1)*wellAreaInt(i+1)*wellZoneSize(i)*momentumDiffusionFactor
    momentumOutToRight  = wellRhoV(i)  *wellAreaInt(i+1)*wellZoneSize(i)*momentumDiffusionFactor

  if(i == firstWellZone)then
    momentumInFromLeft  = 0.
    momentumOutToLeft   = 0.
  elseif(i == wellBottomIndex)then
    momentumInFromLeft  = 0.
    momentumOutToLeft   = 0.
!    momentumInFromLeft  = 0.
!    momentumInFromRight = 0.
  elseif(i == numwellzones)then
    momentumInFromRight = 0.
    momentumOutToRight  = 0.
  endif

  wellRhoV(i) = (wellGRhoVT(i)+momentumInFromRight+momentumInFromLeft &
    -momentumOutToLeft-momentumOutToRight)/wellVol(i)
end do

wellRhoV(wellBottomIndex) = 0.0
i=wellBottomIndex+1
wellRhoV(i) = wellGRhoVT(i)/wellVol(i)



!JFS limit to sonic at outlet
if ((wellRhoV(numWellZones)/wellRho(numWellZones)) > 1000.0) then
  wellRhoV(numWellZones) = 1000.0 * wellRho(numWellZones)
endif

return
end

!-----------------------------------------------------------------------------------------

Real(8) Function frictionFactor(relRough, Re)

Implicit None

Real(8) relRough, Re, ff1, ff2, x, f, aa, bb, cc, temp



if(Re >= 2100.) then
! turbulent Colebrook model
  aa = -2.0d0*log10(relRough/3.7d0 + 12.0d0/Re)
  bb = -2.0d0*log10(relRough/3.7d0 + 2.51d0*aa/Re)
  cc = -2.0d0*log10(relRough/3.7d0 + 2.51d0*bb/Re)
  frictionfactor= ( aa-((bb-aa)**2)/(cc-2.0d0*bb+aa) )**(-2)

else
! laminar
  frictionfactor = 64.0d0/Re

endif

return
end
!-------------------------------------------------------------------------------------------

Subroutine TestCaseFiveMassLoading

Use Globals
Implicit NONE

Real(8) gasLoadingFactor, wasteLoadingFactor, deltaWaste

if(validationSubcase <= 4) then
  gasLoadingFactor   = 0.0
  wasteLoadingFactor = 0.0

elseif(validationSubcase == 5) then
  gasLoadingFactor   = 0.25
  wasteLoadingFactor = 0.0

elseif(validationSubcase == 6) then
  gasLoadingFactor   = 2.5
  wasteLoadingFactor =  0.0

elseif(validationSubcase == 7) then
  gasLoadingFactor   =  2.5
  wasteLoadingFactor =  2.5

elseif(validationSubcase == 8) then
  wasteLoadingFactor =  0.0

! piecewise linear gas loading function
  if(RunTime < 25.0) then
    gasLoadingFactor   = 0.0
  elseif(RunTime < 150.0) then
    gasLoadingFactor   =  0.22*runTime-5.5
  else
    gasLoadingFactor   = -0.01*runTime+29.0
  endif

elseif(validationSubcase == 9) then

! piecewise linear gas loading function
  if(RunTime < 25.0) then
    gasLoadingFactor   = 0.0
  elseif(RunTime < 150.0) then
    gasLoadingFactor   =  0.22*runTime-5.5
  else
    gasLoadingFactor   = -0.01*runTime+29.0
  endif

! piecewise linear waste loading function
  if(RunTime < 28.0) then
    wasteLoadingFactor   = 0.0
  elseif(RunTime < 30.0) then
    wasteLoadingFactor   = 20.0*runTime-560.0
  elseif(RunTime < 40.0) then
    wasteLoadingFactor   = -2.9*runTime+127.0
  elseif(RunTime < 100.0) then
    wasteLoadingFactor   = 11.0
  elseif(RunTime < 150.0) then
    wasteLoadingFactor   = -0.22*runTime+33.0
  else
    wasteLoadingFactor   = 0.0
  endif

endif

deltaGasFromWaste = gasloadingFactor * (runtime - prevtime)

  wellGRhoT  (wellBottomIndex) = wellGRhoT  (wellBottomIndex) + deltaGasFromWaste
  wellGasMass(wellBottomIndex) = wellGasMass(wellBottomIndex) + deltaGasFromWaste
  totalGasInWell               = totalGasInWell   + deltaGasFromWaste
  totalGasInjected             = totalGasInjected + deltaGasFromWaste

deltaWaste = wasteLoadingFactor * (runtime - prevtime)

  wellGRhoT    (wellBottomIndex) = wellGRhoT    (wellBottomIndex) + deltaWaste
  wellWasteMass(wellBottomIndex) = wellWasteMass(wellBottomIndex) + deltaWaste
  totalWasteInWell    = totalWasteInWell    + deltaWaste
  totalWasteFromRepos = totalWasteFromRepos + deltaWaste
  wasteInjected       = wasteInjected       + deltaWaste

return
end
!-------------------------------------------------------------------------------------------
