Subroutine OutputControl
!----------------------------------------------------------------------
! Write results to files
!----------------------------------------------------------------------
Use Globals
Use CDBGlobals
Implicit None

Integer j, idum
integer, save :: kc = 0
integer, save :: ke = 0
Double Precision deltaWBV  !apg REAL
Logical summary, history

summary = .FALSE.
history = .FALSE.

displayedTime           = runTime
displayedPumpPressure   = wellPres(1)
displayedBitAbove       = bitAboveRepository
displayedCavityPressure = cavityPres
! 2026 BC-revision: cavity radius is the inner edge of the cavity-front cell
! (renamed from firstIntactZone -> firstFailedZone); physical meaning unchanged.
displayedCavityRadius   = reposRadius(firstFailedZone)-0.5*reposZoneSize
!displayedCavityRadius   = reposRadiusH(firstFailedZone)  !dkr V1.22
displayedWellBottomPressure        = wellPres(wellBottomIndex)
displayedWasteBoundaryPoreVelocity = wasteboundaryPoreVelocity
displayedFluidizationVelocity      = fluidizationVelocity
displayedEquivDrilledCavityRadius  = equivDrilledCavityRadius

if (maxTensileFailedIndex > 0) then
  displayedTensileFailedRadius = reposRadiusH(maxTensileFailedIndex+1)
else
  displayedTensileFailedRadius = displayedCavityRadius
end if
if (maxShearFailedIndex > 0) then
  displayedShearFailedRadius = reposRadiusH(maxShearFailedIndex+1)
else
  displayedShearFailedRadius = displayedCavityRadius
end if

! 2026 BC-revision: under the new scheme stress in the failed annulus is zero
! (no deviatoric capacity), so radEffStress(firstFailedZone) = 0 by
! construction. The physically meaningful displayed scalar is the radial
! effective stress at the failed/intact interface -- the inner edge of the
! intact region where the BC is applied.
displayedRadEffStress    = radEffStress(firstIntactZone)
displayedFirstFailedZone = firstFailedZone
displayedFirstIntactZone = firstIntactZone

displayedGasInjected  = totalGasInjected
displayedMudInWell    = totalMudInWell
displayedWasteInWell  = totalWasteInWell
displayedGasInWell    = totalGasInWell
displayedSaltInWell   = totalSaltInWell
displayedMudVelocity  = wellV(numWellZones)
displayedMudVolRate   = wellV(numWellZones)*exitPipeArea
displayedMudEjected   = mudMassEjected
displayedGasEjected   = gasMassEjected
displayedWasteEjected = wasteMassEjected
displayedSaltEjected  = saltMassEjected
displayedMudExitFraction   = wellMudVol(numWellZones)/wellVol(numWellZones)
displayedGasExitFraction   = wellGasVol(numWellZones)/wellVol(numWellZones)
displayedWasteExitFraction = wellWasteVol(numWellZones)/wellVol(numWellZones)
displayedSaltExitFraction  = wellSaltVol(numWellZones)/wellVol(numWellZones)
displayedGasPosInWell      = curGasPosInWell
displayedWastePosInWell    = curWastePosInWell
displayedControl = TRIM(stepControlName(stepControl))
displayedCell    = Float(cellcontrol(stepControl))
splvoleq         = (TotVol-cutVol)/(1.0-uncompactedWastePorosity)


!*********************
! spatial data control
!*********************

 WHOLE = .FALSE.
 if (runTime >= spatialSaveTime(spatialSaveIndex)) then

   CALL CPU_TIME(cpuEnd)
   cpuTime = cpuEnd/60. - cpuBegin

   WRITE(6,1001)   &
     runIndex, ' Whl Time:',curtime,displayedControl,displayedCell, &
     curCavityRadius,  splvol2, curWasteInWell, CurWastePosInWell, beddepth
     WRITE(6,1003)cpuWBF, cpuWSF, cpuWSS, cpuTime
   WRITE(diagnosticFileID,1001) &
     runIndex, ' Whl Time:',curtime,displayedControl,displayedCell, &
     curCavityRadius,  splvol2, curWasteInWell, CurWastePosInWell, beddepth
     WRITE(diagnosticFileID,1003)cpuWBF, cpuWSF, cpuWSS, cpuTime

   summary = .TRUE.


   !***********************************
   !  For test case #5 only
   if(validationtestCase == 5)Then
     Call WriteToWellboreValidationFile
   endif
   !***********************************




     WHOLE = .TRUE.
     CALL DBOTIME(IDBout,CurTime, WHOLE, IERRDB)
     numCDBsteps = numCDBsteps + 1
     CALL writeEleToCDB
     CALL writeHisToCDB
     history = .TRUE.
     CALL DBOSTEP( IDBout, idum,idum, IERRDB)


   j = spatialSaveIndex+1
   do while (runTime >= spatialSaveTime(j))
     j = j+1
   end do
   spatialSaveIndex = j

 endif

1001 Format(1x,i8,a,1pe10.3,1x,a,i4, &
   ' CavRad:',1pe10.3,' SplVol2:', 1pe10.3, /, 9x,&
   ' WasWel:',1pe10.3,' WasPos:', 1pe10.3,' BED:', 1pe10.3)

1003 Format(9x, ' cpuWBF:',1pe9.2,' cpuRPF:',1pe9.2,' cpuRPS:',1pe9.2, ' CPU:', 1pe9.2)



!*********************
! Time series data
!*********************

 if (runTime >= timeSaveTime(timeSaveIndex)) then

   CALL CPU_TIME(cpuEnd)
   cpuTime = cpuEnd/60. - cpuBegin

   IF( .NOT.summary)THEN
     WRITE(6,1001)   &
        runIndex, ' His Time:',curtime,displayedControl,displayedCell, &
        curCavityRadius,  splvol2, curWasteInWell, CurWastePosInWell, beddepth
     WRITE(6,1003)cpuWBF, cpuWSF, cpuWSS, cpuTime
     WRITE(diagnosticFileID,1001) &
        runIndex, ' His Time:',curtime,displayedControl,displayedCell, &
        curCavityRadius,  splvol2, curWasteInWell, CurWastePosInWell, beddepth
     WRITE(diagnosticFileID,1003)cpuWBF, cpuWSF, cpuWSS, cpuTime
     summary = .TRUE.
   ENDIF



   ! DKR: added for VMS CDB version
   IF(.NOT.WHOLE)THEN
     CALL DBOTIME(IDBout,CurTime, WHOLE, IERRDB)
     numCDBsteps = numCDBsteps + 1
     CALL WriteHisToCDB
     history = .TRUE.
     CALL DBOSTEP( IDBout, idum,idum, IERRDB)
   ENDIF

   j = timeSaveIndex+1
   do while (runTime >= timeSaveTime(j))
     j = j+1
   end do
   timeSaveIndex = j



! output history if waste boundary velocity changes by > 10%
  elseif(numCDBsteps < 2500 .AND. .NOT.WHOLE) then

    !ouput if boundary pore velocity changes by more than 20%
    ! RM2019 - changed to be a proper inequality comparison for REAL
    ! if(curWasteBoundaryPoreVelocity .ne. 0.0d0)then
    if (ABS(curWasteBoundaryPoreVelocity) > 0.0d0) then
       deltaWBV = (curWasteBoundaryPoreVelocity-prevWasteBoundaryPoreVelocity) &
              /curWasteBoundaryPoreVelocity
    elseif (abs(prevWasteBoundaryPoreVelocity) > 1.0d-6) then
      deltaWBV = 0.21
    else
      deltaWBV = 0.0
    endif


    if(abs(deltaWBV)  > 0.2)then
      CALL CPU_TIME(cpuEnd)
      cpuTime = cpuEnd/60. - cpuBegin

      if(.NOT.summary)then
        WRITE(6,1001)   &
          runIndex, ' His Time:',curtime,displayedControl,displayedCell, &
          curCavityRadius,  splvol2, curWasteInWell, CurWastePosInWell, beddepth
        WRITE(6,1003)cpuWBF, cpuWSF, cpuWSS, cpuTime
        WRITE(diagnosticFileID,1001) &
          runIndex, ' His Time:',curtime,displayedControl,displayedCell, &
          curCavityRadius,  splvol2, curWasteInWell, CurWastePosInWell, beddepth
        WRITE(diagnosticFileID,1003)cpuWBF, cpuWSF, cpuWSS, cpuTime
      endif


      CALL DBOTIME(IDBout,CurTime, WHOLE, IERRDB)
      numCDBsteps = numCDBsteps + 1
      CALL WriteHisToCDB
      CALL DBOSTEP( IDBout, idum,idum, IERRDB)
    endif

 endif




   !********************************************************
   !  For test case#1 only
   if (validationTestCase.EQ.1) then
      if (runTime >= chanSaveTime(chanSaveIndex) ) then

         Call WriteToChanValidationFile

         if (chanSaveIndex < 5) then
           chanSaveIndex = chanSaveIndex +1
         else
           Call CloseRunFiles
           stop  'Normal Termination Chan Test Case'
         end if
      end if

   end if


   !***************************************************
   !  test case #4
   if ((validationTestCase == 4).AND.(runIndex > 1)) then

     ! write to stress validation file
     if (runTime >= stressSaveTime) then
        stressSaveTime = stressSaveTime + stressSaveDelta
        Call WritetoStressValidationFile
     endif

! write to fluidization validation file - moved to wellborecalc after
! fluidization calculation
!     if (runtime >= fluidizationSaveTime) then
!        fluidizationSaveTime = fluidizationSaveTime + fluidizationSaveDelta
!        Call WriteToFluidizationValidationFile
!     endif

     ! write to coupling validation file
     if ((bitAboveRepository > -1.0).AND.(bitAboveRepository < 0.025)) then
       if (kc == 9999) then
         Call WriteToCouplingValidationFile
         kc = 0
       else
         kc = kc + 1
       endif
     endif

     ! write to expulsion validation file
     if (bitAboveRepository < 0.1) then
       if (ke == 9999) then
         Call WriteToExpulsionValidationFile
         ke = 0
       else
         ke = ke + 1
       endif
     endif

  endif
!***********************************************************
!endif

return
end
