!-----------------------------------------------------------------------------------------
!
! File FileWrite.f90, contains routines for the writing to save files
!
!-----------------------------------------------------------------------------------------


Subroutine CloseRunFiles
Use Globals


! test Case #1
IF(validationtestCase ==1)THEN
  close(chanValidationFileID)

! test case #4
ELSEIF(validationTestCase == 4)THEN
 Close(couplingValidationFileID)
 Close(stressValidationFileID)
 Close(fluidizationValidationFileID)
 Close(expulsionValidationFileID)
 Close(fluidizationTimeValidationFileID) !trz

! test case #5
ELSEIF(validationtestcase == 5)THEN
 Close(wellboreValidationFileID)

ENDIF

return
end
!-----------------------------------------------------------------------------------------



Subroutine writeParameters(FileID)
Use Globals
Implicit None

Integer FileID
Double Precision temp  !apg REAL


Write(FileID, '(///A)') '*******************************************************************************'
Write(FileID, '(A)')    'PARAMETERS USED IN THIS RUN'
Write(FileID, '(A)')    '*******************************************************************************'
Write(FileID, '(A)') ''
Write(FileID, '(A)') 'REPOSITORY'
Write(FileID, 10) 'Land Elevation                   (m): ', surfaceElevation
Write(FileID, 10) 'Repository top                   (m): ', repositoryTop
Write(FileID, 10) 'Total Thickness                  (m): ', repositoryThickness
Write(FileID, 20) 'DRZ Thickness                    (m): ', dRZThickness
Write(FileID, 60) 'DRZ Permeability               (m^2): ', dRZPerm
Write(FileID, 10) 'Outer Radius                     (m): ', repositoryOuterRadius
Write(FileID, 30) 'Initial Gas Pressure             (m): ', repositoryInitialPressure
!Write(FileID, 30) 'Far-Field Pore Pressure          (m): ', farFieldPorePressure
Write(FileID, 30) 'Far-Field In-Situ Stress         (m): ', farFieldStress
Write(FileID, '(A)') ''
Write(FileID, '(A)') 'WASTE'
Write(FileID, 10) 'Porosity                         (-): ', repositoryInitialPorosity
Write(FileID, 60) 'Permeability                   (m^2): ', repositoryInitialPerm
Write(FileID, 10) 'Forchheimer Beta                 (-): ', ForchBeta
Write(FileID, 10) 'Biot Beta                        (-): ', biotBeta
Write(FileID, 20) 'Poissons   Ratio                 (-): ', poissonsRatio
Write(FileID, 30) 'Cohesion                        (Pa): ', cohesion
Write(FileID, 10) 'Friction Angle                 (deg): ', frictionAngle
Write(FileID, 30) 'Tensile Strength                (Pa): ', tensileStrength
Write(FileID, 30) 'Failure Characteristic Length    (m): ', Lt
Write(FileID, 50) 'Particle Diameter                (m): ', particleDiameter
!Write(FileID, 40) 'Gas Density at STP          (kg/m^3): ', gasBaseDensity
Write(FileID, 60) 'Gas Viscosity                 (Pa-s): ', gasViscosity
Write(FileID, '(A)') ''
Write(FileID, '(A)') 'MUD'
Write(FileID, 10) 'Density                     (kg/m^3): ', initialMudDensity
Write(FileID, 60) 'Viscosity                     (Pa-s): ', mudViscosity
Write(FileID, 60) 'Wall Roughness Pipe              (m): ', wallRoughness(1)
Write(FileID, 60) 'Wall Roughness Annulus           (m): ', wallRoughness(2)
Write(FileID, 20) 'Max Solids Vol. Frac.         (Pa-s): ', mudSolidsMax
Write(FileID, 20) 'Solids Viscosity Exp.         (Pa-s): ', mudSolidsViscosityExponent
Write(FileID, '(A)') ''
Write(FileID, '(A)') 'WELLBORE/DRILLING'
Write(FileID, 50) 'Bit Diameter                   (m): ', bitDiameter
Write(FileID, 50) 'Pipe Diameter                  (m): ', pipeDiameter
Write(FileID, 50) 'Collar Diameter                (m): ', collarDiameter
Write(FileID, 50) 'Pipe Inside Diameter           (m): ', pipeInsideDiameter
Write(FileID, 60) 'Collar Length                  (m): ', collarLength
Write(FileID, 60) 'Exit Pipe Length               (m): ', exitPipeLength
Write(FileID, 60) 'Exit Pipe Diameter             (m): ', exitPipeDiameter
Write(FileID, 40) 'Drilling Rate                (m/s): ', drillingRate
Write(FileID, 20) 'Bit Above Respository          (m): ', initialBitAboveRepository
Write(FileID, 50) 'Mud Pump Rate              (m^3/s): ', mudPumpRate
Write(FileID, 50) 'Max Pump Pressure             (Pa): ', maxPumpPressure
Write(FileID, 20) 'DDZ Thickness                  (m): ', dDZThickness
Write(FileID, 60) 'DDZ Permeability             (m^2): ', dDZPerm
Write(FileID, 20) 'Stop Drill Exit Vol Rate   (m^3/s): ', stopDrillingExitVolRate
Write(FileID, 20) 'Stop Pump Exit Vol Rate    (m^3/s): ', stopPumpingExitVolRate
Write(FileID, 10) 'Stop Drilling Time             (s): ', stopDrillingTime
Write(FileID, '(A)') ''
Write(FileID, '(A)') 'COMPUTATIONAL'
!DKR
Write(FileID, 75) 'Spherical/Cylindrical        (S/C): ', geometry
Write(FileID, 75) 'Allow Fluidization         (Y/N/A): ', allowFluidization
Write(FileID, 50) 'Max Run Time                   (s): ', maxTime
Write(FileID, 50) 'Respository Cell Length        (m): ', initialReposZoneSize
Write(FileID, 52) 'Radius, Growth Rate          (m,-): ', ReposRadius1, growthRate
Write(FileID, 50) 'Wellbore Cell Length           (m): ', initialWellZoneSize
Write(FileID, 50) 'Wellbore Cell Growth Rate      (-): ', wellGrowthRate
Write(FileID, 65) 'First Wellbore Zone            (-): ', firstWellZone

Write(FileID, 60) 'Well Stability factor          (-): ', wellStabilityFactor
Write(FileID, 60) 'Repository Stability factor    (-): ', reposStabilityFactor
Write(FileID, 60) 'Mass Diffusion factor          (-): ', massDiffusionFactor
Write(FileID, 60) 'Momentum Diffusion factor      (-): ', momentumDiffusionFactor

if(validationTestCase >0)then
  Write(FileID, '(A)') ''
  Write(FileID, '(A)') 'VALIDATION'
  temp = dble(validationtestCase) + 0.1*dble(validationSubCase)  !apg real(
  Write(FileID, '(a,f3.1)') 'Validation Test Case           (-): ', temp
endif

Write(FileID, 60) 'Initial Cavity Radius          (-): ', inputCavityRadius
Write(FileID, 60) 'Minimum Characteristic Vel     (-): ', minCharVel
Write(FileID, 65) 'Minimum Number Zones/Lt        (-): ', minNumLt

Write(FileID, '(A)') ''
Write(FileID, '(A)') 'PARAMETERS'
Write(FileID, 60) 'Pi                             (-): ',Pi
Write(FileID, 60) 'Atmospheric Pressure          (Pa): ',AtmosphericPressure
Write(FileID, 60) 'gravity                    (m/s^2): ',gravity
Write(FileID, 60) 'Gas Constant              (J/kg K): ',GasConstant
Write(FileID, 60) 'Repository Temperature         (K): ',ReposTemp
Write(FileID, 60) 'Reference gas Density     (kg/m^3): ',gasBaseDensity
Write(FileID, 60) 'Water Compressibility       (1/Pa): ',WaterCompressibility
Write(FileID, 60) 'Waste Density             (kg/m^3): ',WasteDensity
Write(FileID, 60) 'Salt Density              (kg/m^3): ',SaltDensity
Write(FileID, 60) 'Shape Factor                   (-): ',ShapeFactor
Write(FileID, 60) 'Tensile Velocity             (m/s): ',TensileVelocity
Write(FileID, 60) 'Bit Nozzle Number              (-): ',BitNozzleNumber
Write(FileID, 60) 'Bit Nozzle Diameter            (m): ',BitNozzleDiameter
Write(FileID, 60) 'Choke Efficiency               (-): ',ChokeEfficiency
!dkr changed for QA0110
WRITE(FileID,'(//)')
Write(FileID, 40) 'Gas Density at STP          (kg/m^3): ', gasBaseDensity
WRITE(FileID,'(///)')

10 FORMAT (A, 1pe11.4)
20 FORMAT (A, 1pe11.4)
30 FORMAT (A, 1pe11.4)
40 FORMAT (A, 1pe11.4)
50 FORMAT (A, 1pe11.4)
52 FORMAT (A, f6.3, f7.3)
60 FORMAT (A, 1pe11.4)
65 FORMAT (A, I3)
70 FORMAT (' Perm from Porosity         (Y/N):', 5X, A1)
!DKR
75 FORMAT (A, A1)
95 FORMAT (A,1pe10.3,i4)

!
!75 FORMAT (' Numerical Method           (E/I):', 5X, A1)
!80 FORMAT (' Spherical/Cylindrical      (S/C):', 5X, A1)
!90 FORMAT (' Allow Fluidization         (Y/N):', 5X, A1)

return
end

!---------------------------------------------------------------------------------------


Subroutine writeTCFileHeaders
Use Globals
Implicit None




if(validationTestCase == 4)then


  ! coupling file
  Write(couplingValidationFileID, '(A)') 'Program DR_SPALL - WIPP PA 2003'
  Write(couplingValidationFileID, '(A)') 'ASCII Output file for Test Case #4'
  Write(couplingValidationFileID, '(A)') 'Verification of coupling between Repository and Wellbore'
  Write(couplingValidationFileID, '(A)') ''
  Write(couplingValidationFileID, '(A30, 1pE30.15)') 'Initial Repos Pressure (Pa)', RepositoryInitialPressure
  Write(couplingValidationFileID, '(A30, 1pE30.15)') 'Initial Gas in Repos (kg)  ', InitialGasInRepository
  Write(couplingValidationFileID, '(A)') ''
  Write(couplingValidationFileID, '(11A15)') 'Runtime', 'Bit Above', 'Repository', 'Cavity', &
              'Well Bottom', 'Total Gas', 'Total Gas', 'Gas Mass ', 'Gas Mass','        ', 'Mass Bal'
  Write(couplingValidationFileID, '(11A15)') '(sec)', 'Repository(m)', 'Penetrated(T/F)', 'Pressure(Pa)', &
              'Pressure (Pa)', 'In Well(kg)', 'Injected(kg)','In Repos(kg)', 'From Repos(kg)', 'Gas storage(kg)', &
              'Error(-)'

  ! stress file
  Write(stressValidationFileID, '(A)') 'Program DR_SPALL - WIPP PA 2003'
  Write(stressValidationFileID, '(A)') 'ASCII Output file for Test Case #4'
  Write(stressValidationFileID, '(A)') 'Verification of Stress Calculations'
  Write(stressValidationFileID, '(A)') ''
  Write(stressValidationFileID, 100) 'Tensile Str(Pa)', TensileStrength
  Write(stressValidationFileID, 100) 'PressureFF(Pa) ', farFieldPorePressure
  Write(stressValidationFileID, 100) 'SressFF  (Pa)  ' , farFieldStress
  Write(stressValidationFileID, 100) 'Forch Beta(-)  ' , forchbeta
  Write(stressValidationFileID, 100) 'Biot Beta (-)  ' , biotbeta
  Write(stressValidationFileID, 100) 'Poisson   (-)  ' , poissonsRatio
100 FORMAT(A15, 1pE15.7)

  ! fluidization file
  Write(fluidizationValidationFileID, '(A)') 'Program DR_SPALL - WIPP PA 2003'
  Write(fluidizationValidationFileID, '(A)') 'ASCII Output file for Test Case #4'
  Write(fluidizationValidationFileID, '(A)') ''
  Write(fluidizationValidationFileID, 200) 'Shape Factor (-)      ', ShapeFactor
  Write(fluidizationValidationFileID, 200) 'Particle Diam (m)     ', ParticleDiameter
  Write(fluidizationValidationFileID, 200) 'Gravity (kg*m/sec^2)  ', gravity
  Write(fluidizationValidationFileID, 200) 'Gas Viscosity (Pa*s)  ', GasViscosity
  Write(fluidizationValidationFileID, 200) 'Waste Density (kg/m^3)', WasteDensity
  Write(fluidizationValidationFileID, 200) 'Porosity (-)          ', repositoryInitialPorosity
  Write(fluidizationValidationFileID, '(A)') ''
200 FORMAT(A30, F15.7)
  ! expulsion file
  Write(expulsionValidationFileID, '(A)') 'Program DR_SPALL - WIPP PA 2003'
  Write(expulsionValidationFileID, '(A)') 'ASCII Output file for Test Case #4'
  Write(expulsionValidationFileID, '(A)') 'Verification of Mass Expelled from Wellbore'
  Write(expulsionValidationFileID, '(A)') ''
  write(expulsionValidationFileID, '(9A15)') 'Runtime', 'Repository',  'Zones', 'Mass Waste', &
                  'Waste in', 'Total Waste', 'Waste Mass', 'Waste Position', 'Mass Balance'
  write(expulsionValidationFileID, '(10A15)') '(sec)', 'Penetrated(T/F)', 'Removed(-)', 'Removed(kg)', &
                  'Store (kg)', 'In Well (kg)', 'Ejected (kg)', 'In Well (m)', 'Error (-)'
 
  ! fluidization time file !trz
  Write(fluidizationTimeValidationFileID, '(A)') 'Program DR_SPALL - WIPP PA 2003'
  Write(fluidizationTimeValidationFileID, '(A)') 'ASCII Output file for Test Case #4'
  Write(fluidizationTimeValidationFileID, '(A)') ''
  Write(fluidizationTimeValidationFileID, '(A)') 'Zone      Fluidization Time'

endif





if(validationTestCase == 5)then

  Write(wellboreValidationFileID, '(A)') 'Program DR_SPALL - WIPP PA 2003'
  Write(wellboreValidationFileID, '(A)') 'ASCII Output file for Test Case #5'
  Write(wellboreValidationFileID, '(A)') 'Wellbore Data'
  Write(wellboreValidationFileID, '(A)') ''
  Write(wellboreValidationFileID, 300) 'Mud Pump Rate (m^3/s)', mudPumpRate
  Write(wellboreValidationFileID, 300) 'Mud Density (kg/m^3)', initialMudDensity
  Write(wellboreValidationFileID, 300) 'Choke Efficiency (-)', chokeEfficiency
  Write(wellboreValidationFileID, 300) 'Pipe Roughness (m)', wallRoughness(1)
  Write(wellboreValidationFileID, 300) 'Wall Roughness (m)', wallRoughness(2)

endif
300 FORMAT(A25, 1pE15.7)

RETURN
END

!----------------------------------------------------------------------------------------

Subroutine WriteSummaryFile

Use Globals
Implicit None


! Data
Write(summaryFileID, 1000)'Stop Time          :', displayedTime
Write(summaryFileID, 1000)'Pump Pres.         :', displayedPumpPressure
Write(summaryFileID, 1000)'Bit Above          :', displayedBitAbove
Write(summaryFileID, 1000)'Bottom Pres.       :', displayedWellBottomPressure
Write(summaryFileID, 1000)'Bdy. Pres.         :', displayedCavityPressure
Write(summaryFileID, 1000)'Bdy. Pore Vel.     :', displayedWasteBoundaryPoreVelocity
Write(summaryFileID, 1000)'Fluidization V.    :', displayedFluidizationVelocity

Write(summaryFileID, 1000)'Drilled Rad.       :', displayedEquivDrilledCavityRadius
Write(summaryFileID, 1000)'Cavity Rad.        :', displayedCavityRadius
Write(summaryFileID, 1000)'Tensile F. Rad.    :', displayedTensileFailedRadius
Write(summaryFileID, 1000)'Shear F. Rad.      :', displayedShearFailedRadius
Write(summaryFileID, 1000)'Rad. Eff. Stress   :', displayedRadEffStress
! 2026 BC-revision: report both interface indices in the summary file.
Write(summaryFileID, 1000)'First Failed Zone  :', displayedFirstFailedZone
Write(summaryFileID, 1000)'First Intact Zone  :', displayedFirstIntactZone
Write(summaryFileID, 1000)'Gas Injected       :', displayedGasInjected
Write(summaryFileID, 1000)'Mud in Well        :', displayedMudInWell
Write(summaryFileID, 1000)'Gas in Well        :', displayedGasInWell
Write(summaryFileID, 1000)'Waste in Well      :', displayedWasteInWell
Write(summaryFileID, 1000)'Salt in Well       :', displayedSaltInWell

Write(summaryFileID, 1000)'Exit Velocity      :', displayedMudVelocity
Write(summaryFileID, 1000)'Mud Vol Ejec rate  :', displayedMudVolRate
Write(summaryFileID, 1000)'Mud Ejected        :', displayedMudEjected
Write(summaryFileID, 1000)'Gas Ejected        :', displayedGasEjected
Write(summaryFileID, 1000)'Waste Ejected      :', displayedWasteEjected
Write(summaryFileID, 1000)'Salt Ejected       :', displayedSaltEjected
Write(summaryFileID, 1000)'Mud Exit Frac.     :', displayedMudExitFraction
Write(summaryFileID, 1000)'Gas Exit Frac.     :', displayedGasExitFraction
Write(summaryFileID, 1000)'Waste Exit Frac.   :', displayedWasteExitFraction
Write(summaryFileID, 1000)'Salt Exit Frac.    :', displayedSaltExitFraction
Write(summaryFileID, 1000)'Gas Pos. In Well   :', displayedGasPosInWell
Write(summaryFileID, 1000)'Waste Pos. in Well :', displayedWastePosInWell

Write(summaryFileID, 1000)'Cuttings Mass      :', cutMass
Write(summaryFileID, 1000)'Total Mass         :', totMass
Write(summaryFileID, 1000)'Spall Mass         :', splmass
Write(summaryFileID, 1000)'Spall2 Mass        :', splMass2
Write(summaryFileID, 1000)'Cut   Vol (Equiv)  :', cutVoleq
Write(summaryFileID, 1000)'Total Vol (Equiv)  :', totvoleq
Write(summaryFileID, 1000)'Spall Vol (Equiv)  :', splvoleq
Write(summaryFileID, 1000)'Spall Vol2(Equiv)  :', splvol2


1000 FORMAT(1x,a,1x,1pe11.4)


return
end

!-------------------------------------------------------------------------------



Subroutine WriteToChanValidationFile

Use Globals
Implicit None

Integer i, int, imax
Real(8) zeta, dimensionlessPsi, tau

i   = 1
int = 1
tau = runTime/characteristicTime

zeta = 0.

do while (zeta.LT.4.0)

    zeta = (reposRadius(i)/initialCavityRadius-1.0)/sqrt(runTime/characteristicTime)
    imax = i
    i = i+1

end do

if (imax.gt.25) then
  int = imax/25
endif

!Write(chanValidationFileID, '(a5, i5)') 'imax: ', imax
!Write(chanValidationFileID, '(a5, i5)') 'int: ', int

Write(chanValidationFileID, '(1X)')
Write(chanValidationFileID, 50) 'Run time (s): ', runTime
Write(chanValidationFileID, 50) 'Tau(-): ', tau
Write(chanValidationFileID, '(5A20)') 'Radius (m)', 'PorePressure (Pa)', 'Zeta(-)', 'Psi(-)' !, 'RadEffStress (Pa)'

do i = 1, imax, int
    !zeta = (reposRadius(i)/(originalBitDiameter/2)-1.0)/sqrt(runTime/characteristicTime)
    !zeta = (reposRadius(i)/reposRadius(1)-1.0)/sqrt(runTime/characteristicTime)
    zeta = (reposRadius(i)/initialCavityRadius-1.0)/sqrt(runTime/characteristicTime)
    dimensionlessPsi = (reposPres(i)/farFieldPorePressure)**2

    Write(chanValidationFileID, 100) &
       reposRadius(i), reposPres(i), zeta, dimensionlessPsi !, radEffStress(i)
end do


50  FORMAT (A20, 1pE20.5)
100 FORMAT (1p6E20.7)

return
end

!-----------------------------------------------------------------------------------------------------
!test case #4
Subroutine WriteToStressValidationFile

Use Globals
Implicit None

Integer i

if (maxTensileFailedIndex <= 20 .or. &
   (maxTensileFailedIndex >  100 .and. maxTensileFailedIndex <250)) then


  write(stressValidationFileID, '(A15)') ''
  write(stressValidationFileID, 100) 'Runtime(sec)  =', runtime
  write(stressValidationFileID, 100) 'CavPres(Pa)   =', cavityPres
  write(stressValidationFileID, 100) 'CavRadius(m)  =', curCavityRadius
  write(stressValidationFileID, 100) 'DrilledRad(m) =', curDrilledRadius
  write(stressValidationFileID, 100) 'CavityVol(m^3)=', wellVol(wellbottomindex)
  write(stressValidationFileID, 100) 'Pff(Pa)       =', repospres(numReposZones)
  ! 2026 BC-revision: emit both interface indices for traceability.
  Write(stressValidationFileID, '(A17, I5)')'firstFailedZone=' , firstFailedZone
  Write(stressValidationFileID, '(A17, I5)')'firstIntactZone=' , firstIntactZone
  write(stressValidationFileID, '(A15, 8a21)') &
        'zone index', 'Radius(m)', 'PorePres(Pa)', &
        'ElastStr(Pa)', 'SeepStr(Pa)', 'EffStre(Pa)','Failed(T/F)', 'Fluidized(-)'

  do i = max(1,firstFailedZone-10), firstFailedZone+20
    write (stressValidationFileID, 200) &
      i , reposRadius(i), repospres(i), radElasticStress(i), radSeepageStress(i), &
      radEffStress(i), tensileFailureStarted(i), drillingfailure(i)*fractionFluidized(i)
  end do

endif

200 FORMAT(I15, 1p5E21.13, L15, 0pF15.4)
100 FORMAT(A17, 1pE21.13)


return
end
!--------------------------------------------------------------------------------------------------------
! test case #4
Subroutine WriteToFluidizationValidationFile

Use Globals
Implicit None

Integer i
Real(8) curGasDensity

 if (firstFailedZone <=70 .or. & !trz
   (firstFailedZone >100 .and. firstFailedZone < 150)) then 
  i = 0

  curGasDensity = gasBaseDensity*(reposPres(firstFailedZone)/AtmosphericPressure)

  write(fluidizationValidationFileID, '(A)') ''
  write(fluidizationValidationFileID, 100) 'Runtime (sec)              =', runtime
  write(fluidizationValidationFileID, 100) 'Cavity Pressure(Pa)        =', cavityPres
  write(fluidizationValidationFileID, 100) 'Cavity Radius(m)           =', curCavityRadius
  write(fluidizationValidationFileID, 100) 'Gas Density (kg/m^3)       =', curGasDensity
  write(fluidizationValidationFileID, 100) 'Fluidization Velocity(m)   =', fluidizationVelocity
  ! 2026 BC-revision: relabel for new nomenclature; emit both interface indices.
  write(fluidizationValidationFileID, 100) 'Superficial Gas Velocity(m) '
  write(fluidizationValidationFileID, 100) '        (First Failed Zone)=', superficialVelocity(firstFailedZone)
  write(fluidizationValidationFileID, 100) 'Waste In Well (kg)         =', totalWasteInWell
  Write(fluidizationValidationFileID, '(A30, I5)') 'firstFailedZone       ', firstFailedZone
  Write(fluidizationValidationFileID, '(A30, I5)') 'firstIntactZone       ', firstIntactZone
  write(fluidizationValidationFileID, '(A)') ''
  write(fluidizationValidationFileID, 150)  &
  'Cell', '          ',  'Failure',       'Fluidization', 'Fluidization',  'Fraction'
  write(fluidizationValidationFileID, 150) &
  'index', 'Radius(m)',  'Completed(T/F)','Start(T/F)',   'Complete(T/F)', 'Fluidized'

  do i = max(1,firstFailedZone-10), firstFailedZone+20
    write (fluidizationValidationFileID, 200) &
      i , reposRadius(i), tensileFailureCompleted(i), &
      fluidizationStarted(i), fluidizationCompleted(i), drillingFailure(i)*fractionFluidized(i)
   end do

endif

100 FORMAT (A30, 1pE21.13)
150 FORMAT (A8, A21, 6X, 4A15)  !apg V1.22 was (9A15)
200 FORMAT (I8, 1pE21.13, 6X, 3L15, 0pF15.4)  !apg V1.22


return
end

!-------------------------------------------------------------------------------------------------
! test case #4 !trz
Subroutine WriteToFluidizationTimeValidationFile(izone)
 
Use Globals
Implicit None
 
Integer izone

write (fluidizationTimeValidationFileID, 200) izone, fluidizationTime(izone)

200 FORMAT (I8, 1pE21.13) 

return
end

!-------------------------------------------------------------------------------------------------
! test case #4
Subroutine WriteToCouplingValidationFile

Use Globals
Implicit None
Real(8) currentGasInRepository, rPlus, rMinus, zonevolume, gasLostFromRepository, &
        massBalanceError, temp
Integer i

! include gas in cavity prior to bit penetration
!IF (repositoryPenetrated.EQ. .FALSE.) THEN
!   currentGasInRepository = (cavityPres*initialCavityGasVol)*(gasBaseDensity/AtmosphericPressure)
!ELSE
!   currentGasInRepository = gasStore
!ENDIF

! calculate the current mass of gas in repository from ideal gas eq. of state

! Initialize gas volume to zero
currentGasInRepository = 0.0


do i = firstFailedZone, numReposZones

    currentGasInRepository = currentGasInRepository + &
        (reposPres(i)*reposVol(i)*repositoryInitialPorosity) &
       *(gasBaseDensity/AtmosphericPressure)

end do

gasLostFromRepository = initialGasInRepository - currentGasInRepository

!massBalanceError = abs(totalGasInjected+gasstore - gasLostFromRepository) &
!                     /(gasstore + totalGasInjected)
massBalanceError = abs(initialGasInRepository  &
                    - (currentGasInRepository + gasStore + totalGasInjected)) &
                     /initialGasInRepository

write(couplingValidationFileID, 400) runTime, BitAboveRepository, repositoryPenetrated, &
                 cavityPres, wellPres(wellbottomIndex), totalGasInWell, totalGasInjected, &
                 currentGasInRepository, gasLostFromRepository, gasStore, massBalanceError


400 FORMAT(2F15.5, L15, 1p8E15.7)

return
end

!-----------------------------------------------------------------------------------------------------------
! test case #4
Subroutine WriteToExpulsionValidationFile

Use Globals
Implicit None
Real(8) massWasteRemoved, initialCavityWasteMass, error
Integer i, zonesRemoved

zonesremoved = firstFailedZone - 1

if(firstFailedZone >1)then
! mass=rho*volume*(1-porosity)
massWasteRemoved = wasteDensity*((2.0d0/3.0d0)*pi*(curCavityRadius**geomExponent)-initialCavityVol)* &
                        (1.0d0-repositoryInitialPorosity)
else
  massWasteRemoved = 0.0
endif

!initialCavityWasteMass = wasteDensity*initialCavityVol*(1.0d0-repositoryInitialPorosity)

if(massWasteRemoved > 0.0)then
  error =  ABS(massWasteRemoved-(wasteStore+wasteMassEjected+totalWasteInWell))/massWasteRemoved
else
  error = 0.0
endif


write(expulsionValidationFileID, 400) runTime, repositoryPenetrated, zonesRemoved, massWasteRemoved, &
                wasteStore, totalWasteInWell, wasteMassEjected, curWastePosInWell, error


400 FORMAT(F15.5, L15, I15, 1p4E15.7, 0pF15.1, 1pE15.4)

return
end
!------------------------------------------------------------------------------------
! test case #5
Subroutine WriteToWellboreValidationFile

Use Globals
Implicit None
Integer i

!write(wellboreValidationFileID, 100)
write(wellboreValidationFileID, '(A)') ''
write(wellboreValidationFileID, '(A15,F12.5)') 'Runtime (sec)=', runtime
write(wellboreValidationFileID, '(A15, 2A15)') 'Well Position (m)', 'Pressure (Pa)', 'Velocity (m/s)'



do i = 1, numwellzones, 5
    write(wellboreValidationFileID, 200) wellPos(i), wellPres(i), wellV(i)
end do

200 format (1pE15.3, 1p2E15.7)

return
end
!-------------------------------------------------------------------------------------
