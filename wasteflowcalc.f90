!--------------------------------------------------------------------------------
!
! File WasteFlowCalc.f90, contains routines for the waste porous flow calculation
!
!--------------------------------------------------------------------------------

Subroutine CalculateWasteFlow

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Integer j  !trz
Real(8) exitGasDen, exitGasArealFlux, boundaryPres, exitGasFlowRate, &
        deltaP, curGasDen, ForchNumber, temp, c1, c2, c3, c4, dr2
Real(8) permEnhanceFactor  !dkr

!if (repositoryPenetrated == .FALSE.) then
  boundaryPres = cavityPres
!else
!  boundaryPres = wellPres(wellBottomIndex)
!end if

!***********************************************
!JFS3 add velocity dependent k (Forchheimer)
permEnhanceFactor = 4.0  !dkr V1.22, was 4.0
do i = firstIntactZone, numReposZones !trz
  curGasDen       = gasBaseDensity*reposPres(i)*invAtmosphericPressure !trz
  ForchNumber     = ABS(ForchBeta)*curGasDen*abs(superficialVelocity(i))*invGasViscosity*invPorosity(i)
  permeability(i) = repositoryInitialPerm*(1.0+permEnhanceFactor*fractionFluidized(i))/(1.0+ForchNumber)
enddo
!permeability(firstIntactZone-1) = permeability(firstIntactZone) !trz
!***********************************************

! time-varying inner waste boundary condition
deltaP = reposPres(firstIntactZone) - boundaryPres
psi(firstIntactZone) = invGasViscosity*(boundaryPres)**2 !trz
!psi(firstIntactZone-1) = invGasViscosity*(boundaryPres)**2 !orig

! flow calculations
  Call CalculateWasteFlowImplicit


!*******************************************
!DKR check Forchhiemer assumption
maxForchRatio = 0.0
do i = firstIntactZone, numreposZones-1

  dr2 = reposRadius(i+1) - reposRadius(i-1)
  c1  = (psi(i+1)-psi(i-1))/dr2
  c2  = (permeability(i+1)-permeability(i-1))/dr2
  c3  = (psi(i+1)-psi(i)  )/reposDRH(i+1)
  c4  = (psi(i)  -psi(i-1))/reposDRH(i)

  temp =  (c3-c4)/reposDR(i)+(geomExponent-1)*c1/reposRadius(i)
! RM2019 - changed to be a proper inequality comparison for REAL
  if (ABS(temp) > 0.0)then
    ForchRatio(i) =  (geomExponent-1)*c1*c2/(permeability(i)*temp)
  else
    ForchRatio(i) = 0.0
  endif
  maxForchRatio = max(ABS(ForchRatio(i)),maxForchRatio)
enddo
!*******************************************



! finally, convert from pseudoPressure to pressure.
do i = firstIntactZone, numReposZones
  temp = gasViscosity*psi(i)  !apg V1.22 temp for negative check
  IF(temp < 0.0) then  !apg V1.22
      call QAABORT ('SQRT(-x) reposPres')  !apg V1.22
  ENDIF  !apg V1.22
  reposPres(i) = DSqrt(temp)  !apg was DSqrt(gasViscosity*psi(i))
end do
!dkr used for zone in cavity on CDB
reposPres(firstIntactZone) = boundaryPres !trz
!reposPres(0) = boundaryPres !orig

! Calculate pore velocity in waste, where positive is out of waste toward cavity
!dkr - added check on flow direction
if(reposPres(firstIntactZone+1) > boundaryPres) then
  poreVelocity(firstIntactZone) = &
      0.5*(reposPres(firstIntactZone+1)-(boundaryPres)) &
         *invGasViscosity*permeability(firstIntactZone)*invPorosity(firstIntactZone) &
         / reposDR(firstIntactZone)
else
  poreVelocity(firstIntactZone) = 0.0
  reposPres(firstIntactZone+1) = boundarypres
endif

do i = firstIntactZone+1, (numReposZones-1)
!dkr changed to account for variable zone size
  poreVelocity(i) = (reposPres(i+1)-reposPres(i-1)) &
      *invGasViscosity*permeability(i)*invPorosity(i)/(reposDRH(i)+reposDRH(i+1))
end do


! superficial or Darcy velocity ( m^3/s/m^2)
do i = firstIntactZone, (numReposZones-1)
  superficialVelocity(i) = porosity(i)*poreVelocity(i)
end do

!dkr: switched to superficial for consistancy with fluidization velocity
!wasteBoundaryPoreVelocity = poreVelocity(firstIntactZone)
wasteBoundaryPoreVelocity = superficialVelocity(firstIntactZone)


!dkr changed to improve centering => improved mass balance
!exitPoreVelocity = 2.0d0*(reposPres(firstIntactZone)-boundaryPres) &
exitPoreVelocity = (reposPres(firstIntactZone+1)-boundaryPres) & !trz
  *invGasViscosity*permeability(firstIntactZone)*invPorosity(firstIntactZone) &
  /reposDR(firstIntactZone)

!DKR
if(exitPoreVelocity < 0.0) Then
  exitporeVelocity = 0.0
!  reposPres(firstIntactZone) = boundaryPres
!  write(6,*) ' WARNING: exitPoreVelocity < 0 '
endif

! Boundary variables
!exitGasDen = 0.5d0*gasBaseDensity*(boundaryPres+reposPres(firstIntactZone)) &
exitGasDen = gasBaseDensity*(reposPres(firstIntactZone)) &
                *invAtmosphericPressure
exitGasArealFlux = exitGasDen*exitPoreVelocity*porosity(firstIntactZone)

! *******************************************************************************
!DKR - added exitPoreArea, moved exitporeVelocity to Globals
if (geometry == 'S') then
  exitGasFlowRate = exitGasArealFlux*2.0*Pi*reposRadiusH(firstIntactZone)**2
  exitPoreArea   = 2.0d0*Pi*(reposRadiusH(firstIntactZone))**2 &
                           *porosity(firstIntactZone)
else
  exitGasFlowRate = exitGasArealFlux*2.0d0*Pi*(reposRadiusH(firstIntactZone))*repositoryThickness
  exitPoreArea   = 2.0d0*Pi*reposRadiusH(firstIntactZone) &
                   *repositoryThickness*porosity(firstIntactZone)
end if


deltaGasFromWaste = deltaTime*exitGasFlowRate
totalGasFromWaste = totalGasFromWaste + deltaGasFromWaste

sumReposGasMass = 0.0
i = firstIntactZone
do while (i < numReposZones .and. reposRadius(i) < 30.0)
  reposGasMass(i) = gasBaseDensity*reposPres(i)*invAtmosphericPressure &
                   *(porosity(i)*reposVol(i))
  SumReposGasMass = sumReposGasMass + reposGasMass(i)
  i = i + 1
enddo

return
end


!--------------------------------------------------------------------------------

Subroutine CalculateWasteFlowImplicit

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
Real(8) dPrime, alpha1, alpha2, compressibility
Real(8) bet, Forchterm, Darcyterm1, DarcyTerm2

! Set up coefficients for tridiagonal inversion of pseudopressures.


! First cell coefficients
i = firstIntactZone + 1 !trz
compressibility = 1.0d0/reposPres(i)
dPrime = permeability(i)/(porosity(i)*gasViscosity*compressibility)
DarcyTerm1 = 1.0d0/reposDRH(i)   -0.5d0*(geomExponent-1)/reposRadius(i)
DarcyTerm2 = 1.0d0/reposDRH(i+1) +0.5d0*(geomExponent-1)/reposRadius(i)
IF(forchBeta > 0.0)THEN
  Forchterm = log(permeability(i+1)/permeability(i-1))  & !trz
            / (4.0*reposDR(i)) !trz
ELSE
  ForchTerm = 0.0
ENDIF
alpha1 = (DarcyTerm1 -ForchTerm) *dPrime*deltaTime/reposDR(i)
alpha2 = (DarcyTerm2 +Forchterm) *dPrime*deltaTime/reposDR(i)
bb(i) = 1.0d0 + alpha1 + alpha2
cc(i) = -alpha2
rr(i) = psi(i)+alpha1*psi(i-1)


! Interior cell coefficients
do i = (firstIntactZone+2), (numReposZones-1) !trz
  compressibility = 1.0d0/reposPres(i)
  dPrime = permeability(i)/(porosity(i)*gasViscosity*compressibility)
  DarcyTerm1 = 1.0d0/reposDRH(i)  -0.5d0*(geomExponent-1)/reposRadius(i)
  DarcyTerm2 = 1.0d0/reposDRH(i+1)+0.5d0*(geomExponent-1)/reposRadius(i)
IF(forchbeta > 0.0)THEN
  Forchterm = log(permeability(i+1)/permeability(i-1))  & !trz
            / (4.0*reposDR(i)) !trz
ELSE
  ForchTerm = 0.0
ENDIF
  alpha1 = (DarcyTerm1 -ForchTerm) *dPrime*deltaTime/reposDR(i)
  alpha2 = (DarcyTerm2 +ForchTerm) *dPrime*deltaTime/reposDR(i)
  aa(i) = -alpha1
  bb(i) =  1.0d0 + alpha1 + alpha2
  cc(i) = -alpha2
  rr(i) = psi(i)
end do


i = numReposZones
! Last cell coefficients
compressibility = 1.0d0/reposPres(i)
DarcyTerm1 = 1.0d0/reposDRH(i)  -0.5d0*(geomExponent-1)/reposRadius(i)
DarcyTerm2 = 1.0d0/reposDRH(i)  +0.5d0*(geomExponent-1)/reposRadius(i)
IF(forchbeta > 0.0)THEN
  Forchterm = log(permeability(i)/permeability(i-1))  & !trz
            / (4.0*reposDR(i)) !trz
ELSE
  ForchTerm = 0.0
ENDIF
dPrime = permeability(i)/(porosity(i)*gasViscosity*compressibility)
alpha1 = (DarcyTerm1 -Forchterm) *dPrime*deltaTime/reposDR(i)
alpha2 = (DarcyTerm2 +Forchterm) *dPrime*deltaTime/reposDR(i)
aa(i) = -alpha1 - alpha2 !trz
bb(i) = 1.0d0 + alpha1 + alpha2 !trz
rr(i) = psi(i)


! Perform inversion.
bet = bb(firstIntactZone+1) !trz
psi(firstIntactZone+1) = rr(firstIntactZone+1)/bet !trz
do i = (firstIntactZone+2), numReposZones !trz
  gam(i) = cc(i-1)/bet
  bet    = bb(i)-aa(i)*gam(i)
  psi(i) = (rr(i)-aa(i)*psi(i-1))/bet
end do


!dkr changed from contant pressure to zero gradient
!psi(numReposZones) = repositoryInitialPressure**2/gasViscosity
!New BC - Comment next line*******************
!psi(numReposZones) = reposPres(numReposZones-1)**2/gasViscosity
do i = (numReposZones-1), firstIntactZone+1, -1 !trz
  psi(i) = psi(i)-gam(i+1)*psi(i+1)
end do


return
end

!-----------------------------------------------------------------------------------------

