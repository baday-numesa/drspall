!-----------------------------------------------------------------------------------------
!
! File MainCalc.f90, contains routines for the main calculation
!
!-----------------------------------------------------------------------------------------

Subroutine RunLoop
!note: this routine is used in both PC and VMS versions

Use Globals
!Use DFLogM

Implicit NONE
!Include 'resource.fd'

Integer i, i4, j  !apg Integer(4)
Logical l4  !apg Logical(4)
!Type(dialog) dlg

if(platform == 'PC')then
 Call Reset(RUNNING)
endif

stopRun   = .FALSE.
pauseRun  = .FALSE.
isRunning = .TRUE.

initialPumpRate = mudPumpRate

volStore   = 0.0
gasStore   = 0.0
wasteStore = 0.0
wasteInjected = 0.0

splVol  = 0.0
TotVol  = 0.0
cutTrueVol = 0.0

lastFailedZone  = 0


! March time.
do while (runTime < maxTime)

  if (.NOT.pauseRun ) then !dkr

    prevTime               = runTime
    prevWellBottomPressure = wellPres(wellBottomIndex)
    if (.NOT.repositoryPenetrated ) then !dkr
      prevCavityPressure = cavityPres
    else
      prevCavityPressure = wellPres(wellBottomIndex)
    end if
    prevPumpPressure         = wellPres(1)
    prevDrilledRadius        = equivDrilledCavityRadius
!dkr
    prevCavityRadius         = reposRadiusH(firstFailedZone)
    prevTensileRadius        = reposRadiusH(maxTensileFailedIndex+1)
    prevWasteBoundaryPoreVelocity = wasteBoundaryPoreVelocity
    prevWasteInWell          = totalWasteInWell
    prevWasteEjected         = wasteMassEjected
    prevMudEjectionVelocity  = wellV(numWellZones)
    prevFluidizationVelocity = fluidizationVelocity
    prevGasInjected          = totalGasInjected
    prevGasInWell            = totalGasInWell
    prevGasEjected           = gasMassEjected
    prevGasPosInWell         = curGasPosInWell
    prevWastePosInWell       = curWastePosInWell



    ! March time.
    Call CalculateTimestep

    runTime  = runTime+deltaTime
    runIndex = runIndex+1

    IF(validationTestCase == 1)Then
      ! Calculate flow in waste
      Call CalculateWasteFlow


    elseif(validationTestCase == 5)then

      ! Calculate wellbore flow.
      Call CalculateWellboreFlow

    else
                  CALL CPU_TIME(cpuT0)
      ! Calculate wellbore flow.
      Call CalculateWellboreFlow
                  CALL CPU_TIME(cpuT1)
                  cpuWBF = cpuWBF + (cpuT1-cpuT0)/60.
      ! Calculate flow in waste
      Call CalculateWasteFlow
                  CALL CPU_TIME(cpuT2)
                  cpuWSF = cpuWSF + (cpuT2-cpuT1)/60.
      ! Calculate stresses.
      Call CalculateWasteStresses
                  CALL CPU_TIME(cpuT3)
                  cpuWSS = cpuWSS + (cpuT3-cpuT2)/60.
    endif

    ! Plot results to screen
    curTime = runTime


    ! Time graphs

    ! pressure
    curPumpPressure = wellPres(1)
    curWellBottomPressure = wellPres(wellBottomIndex)
    if (.NOT. repositoryPenetrated ) then !dkr
      curCavityPressure = cavityPres
    else
      curCavityPressure = wellPres(wellBottomIndex)
    end if



    ! radius
    curDrilledRadius = equivDrilledCavityRadius
    curCavityRadius = reposRadiusH(firstFailedZone)
    if (maxTensileFailedIndex > 0) then
      curTensileRadius = reposRadiusH(maxTensileFailedIndex+1)
    else
      curTensileRadius = reposRadiusH(firstFailedZone)
    end if


    ! velocity
    curWasteBoundaryPoreVelocity = wasteBoundaryPoreVelocity
    curFluidizationVelocity      = fluidizationVelocity
    curMudEjectionVelocity       = wellV(numWellZones)


    ! waste mass
    curWasteInWell  = totalWasteInWell
    curWasteEjected = wasteMassEjected


    ! gas mass
    curGasInjected = totalGasInjected
    curGasInWell   = totalGasInWell
    curGasEjected  = gasMassEjected


    ! position
    i = numWellZones
    do while ((wellGasVol(i)/wellVol(i) .LT. 0.001) .AND. (i .GE. wellBottomIndex))
      i = i-1
    end do
    i = max(i,1)

    curGasPosInWell = -wellDepth + (wellpos(i) - wellpos(wellBottomIndex))

    i = numWellZones
! dkr change to mass fraction check
    do while ((wellWasteMass(i)/(wellRho(i)*wellVol(i)) .LT. 0.0001) .AND. (i .GE. wellBottomIndex))
      i = i-1
    end do
    i = max(i,1)
    curWastePosInWell = -wellDepth + (wellpos(i) - wellpos(wellBottomIndex))


    Call PlotControl

    ! Write results to files
    Call OutputControl

  else
!    Call ShowPauseDialog  ! APG Eliminate "unpause" dialog box

  end if

  ! Exit run loop if menu stop is used
  if (stopRun ) then  !dkr
    exit
  end if

end do

isRunning = .FALSE.
stopRun   = .FALSE.
pauseRun  = .FALSE.
Call RunShutDown

end


!-----------------------------------------------------------------------------------------


Subroutine CalculateTimestep

Use Globals
Implicit NONE

integer i  !apg Integer(4)
Real(8) allowedDeltaTime, maxDeltaTime, trialDeltaTime, maxDCoeff, cCoeff, &
  maxCCoeff, gasFac, soundSpeed, DT, &
  tensileFailureAllowedDeltaTime, wellboreAllowedDeltaTime, reposAllowedDeltaTime


trialDeltaTime   = 1.01*deltaTime
allowedDeltaTime = trialDeltaTime
oldStepControl   = stepControl
oldCellControl   = cellControl(oldStepControl)


         !------------------------
         !    repository flow
         !------------------------


!DKR explicit Courant timestep
!    implicit also, user specified factor*explicit

  gasFac = DSqrt(invGasViscosity)

!trz
  reposAllowedDeltaTime = maxTime !initialize with large time

! RM2019 - changed to be a proper inequality comparison for REAL
  if (ABS(cavityPres) > 0.0) then !skip for test case #1
    i = firstFailedZone
    dCoeff(i) = permeability(i)*invPorosity(i)*gasFac*DSqrt(psi(i))
    reposAllowedDeltaTime = reposDR(i)**2/DCoeff(i)
    cellControl(2)  = i
  endif
!trz

  do i = firstFailedZone+1, numReposZones
    dCoeff(i) = permeability(i)*invPorosity(i)*gasFac*DSqrt(psi(i))
    DT = reposDR(i)**2/DCoeff(i)
    if (DT < reposAllowedDeltaTime) then
      reposAllowedDeltaTime = DT
      cellControl(2)  = i
    end if
  end do

  reposAllowedDeltaTime = reposStabilityFactor*reposAllowedDeltaTime
  if (allowedDeltaTime > reposAllowedDeltaTime) then
    allowedDeltaTime = reposAllowedDeltaTime
    stepControl = 2
  end if




          !--------------------
          ! tensile failure
          !--------------------

if (runtime > (timeOfPenetration-1.0)) then
  ! Tensile-failure timestep limit uses firstIntactZone -- the inner edge of
  ! the intact region, where the active Lt batch progresses.
  tensileFailureAllowedDeltaTime = 0.1*tensileFailureTime(firstIntactZone)
  if (allowedDeltaTime > tensileFailureAllowedDeltaTime) then
    allowedDeltaTime = tensileFailureAllowedDeltaTime
    stepControl = 3
    cellControl(3) = firstIntactZone
  end if

end if



          !-----------------
          !  wellbore flow
          !-----------------


! uses Courant stability criteria
if(validationtestCase == 2)then
  soundspeed = 1100.0

else
  soundspeed = sqrt(1.0d0/(waterCompressibility*initialMudDensity))
!  soundspeed = 1600.0

endif

wellboreAllowedDeltaTime = wellZoneSize(firstWellZone)/soundSpeed
!dkr added for QA0110
cellControl(4) = firstWellZone
do i = firstWellZone+1, numWellZones
  if( i /= wellBottomIndex ) then
    DT = wellZoneSize(i)/soundSpeed
    if (DT < wellboreAllowedDeltaTime) then
      wellboreAllowedDeltaTime = DT
      cellControl(4) = i
    end if
  endif
end do


!DKR stability factor moveded to input
wellboreAllowedDeltaTime = WellStabilityFactor*wellboreAllowedDeltaTime
if (allowedDeltaTime > wellboreAllowedDeltaTime) then
  allowedDeltaTime = wellboreAllowedDeltaTime
  stepControl = 4
end if



          !------------------
          ! final timestep
          !------------------


if (trialDeltaTime > allowedDeltaTime) then
  deltaTime = allowedDeltaTime
else
  deltaTime = trialDeltaTime
  cellControl(1) = oldCellControl
  stepControl = 1
  stepControlName(1) = trim(stepControlName(oldstepcontrol))//'*'
end if



!---------------------------------------------------------
! timestep limit for plotting and startup accuracy
!---------------------------------------------------------


if (runTime > 10.0) then
  maxDeltaTime = 0.1
else if (runTime > 1.0) then
  maxDeltaTime = 0.01
else if (runTime > 0.1) then
  maxDeltaTime = 0.001
else if (runTime > 0.01) then
  maxDeltaTime = 0.0001
else if (runTime > 0.001) then
  maxDeltaTime = 0.0001
else
  maxDeltaTime = 0.00001
end if
if (deltaTime > maxDeltaTime) then
  deltaTime = maxDeltaTime
  Stepcontrol = 5
  cellcontrol(5) = 0
end if

displayedControl = TRIM(stepControlName(stepControl))
displayedCell    = cellcontrol(stepControl)

return
end

!-----------------------------------------------------------------------------------------


