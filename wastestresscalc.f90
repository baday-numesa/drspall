!-----------------------------------------------------------------------------------------
!
! File WasteStressCalc.f90, contains routines for the stress calculation
!
! 2026 BC-revision (radial effective-stress inner boundary condition):
!   Tensile-failed (disaggregated) waste cannot transmit deviatoric stress;
!   it only transmits its local pore pressure isotropically. The intact-rock
!   stress problem must therefore be posed on the *intact* domain only --
!   i.e. starting at the failed/intact interface, not at the cavity wall --
!   and the inner radial-stress BC must be the local pore pressure at that
!   interface (taken from the gas-flow solution which spans the entire waste
!   domain), not the cavity pressure.
!
! Index conventions used here:
!   firstFailedZone : first non-fluidized cell (cavity wall). Was previously
!                     named firstIntactZone in the codebase.
!   firstIntactZone : first cell whose tensile failure has NOT completed --
!                     inner edge of the intact region; derived locally below
!                     by scanning outward from firstFailedZone.
!   When no tensile-failed annulus exists, firstIntactZone == firstFailedZone
!   and the new inner-BC pressure degenerates cleanly to cavityPres
!   (the flow solver overwrites reposPres(firstFailedZone) <- cavityPres).
!
!-----------------------------------------------------------------------------------------

Subroutine CalculateWasteStresses

Use Globals
Implicit NONE

Integer i  !apg Integer(4)
!RM2019 LtZone changed from REAL to INTEGER
integer LtZone
Real(8) temp1, temp2, temp3
Real(8) radTotStress, tanTotStress, mm1
Real(8) mu, S0, preFactor, integral1, integral2
Real(8) pPrime, pff, pPrimeR, Pnew, R1, R2, Rc, Pres1, Pres2, &
        aveRadStress, aveTanStress, aveShrStress, aveShrStrength, aveMeanStress
! 2026 BC-revision locals
Real(8) pInterface, rInterface, drFirst
Logical thinShell, internalFailureAllowed, Ltemp1, Ltemp2, Ltemp3
Data thinShell /.FALSE. /
Data internalFailureAllowed/.TRUE./

! Initialization (m-1)
mm1 = geomExponent-1

! Seepage Force Constant
preFactor = biotBeta*(1.0-2.0*poissonsRatio)/(1.0-poissonsRatio)

! Solid Failure Constants
! TODO(deg/rad): frictionAngle is stored in degrees (parameters.f90, summary
! file label, bounds 30-60); Fortran intrinsic Tan() expects radians. Out of
! scope for this commit -- flagged for a separate dedicated change.
mu = Tan(frictionAngle)
S0 = 0.5*cohesion/(mu+DSqrt(mu**2+1.0))


!------------------------------------------------------------
! 2026 BC-revision: derive firstIntactZone (NEW) and the
! inner-BC pressure/radius for the intact-rock stress problem.
!------------------------------------------------------------

! Scan outward from the cavity-front cell. firstIntactZone = first cell whose
! tensile failure has NOT completed. If firstFailedZone itself is intact, this
! immediately exits with firstIntactZone == firstFailedZone (graceful
! degeneration to old behavior in the no-failure regime).
firstIntactZone = firstFailedZone
do while (firstIntactZone <= numReposZones)
  if (.NOT. tensileFailureCompleted(firstIntactZone)) exit
  firstIntactZone = firstIntactZone + 1
end do

! Inner mechanical boundary is at the inner face of the first intact cell.
rInterface = reposRadiusH(firstIntactZone)

! Interpolate pore pressure to the failed/intact interface using the same
! reposFactor convention used elsewhere for half-grid pressures (see pL). When
! firstIntactZone == firstFailedZone the flow solver has already pinned
! reposPres(firstFailedZone) = cavityPres, so this gives pInterface = cavityPres
! and the new BC reduces to old-code behavior.
if (firstIntactZone == firstFailedZone) then
  pInterface = reposPres(firstFailedZone)   ! == cavityPres after flow solve
else
  pInterface = reposPres(firstIntactZone-1) &
             + reposFactor(firstIntactZone) &
              *(reposPres(firstIntactZone) - reposPres(firstIntactZone-1))
endif


!---------------------------------------------------------------
! Zero stress arrays in the failed-but-not-fluidized annulus so
! per-step radius plots show zero deviatoric stress where the
! waste is disaggregated (only pore pressure carried there).
!---------------------------------------------------------------

do i = firstFailedZone, firstIntactZone-1
  radElasticStress(i) = 0.0d0
  tanElasticStress(i) = 0.0d0
  radSeepageStress(i) = 0.0d0
  tanSeepageStress(i) = 0.0d0
  radEffStress    (i) = 0.0d0
  tanEffStress    (i) = 0.0d0
  shearStress     (i) = 0.0d0
  meanEffStress   (i) = 0.0d0
  shearStrength   (i) = 0.0d0
end do


!------------------------
!Calculate Elastic Stress
!------------------------

! 2026 BC-revision: inner BC moved from cavity wall (cavityPres at
! reposRadius(firstFailedZone)) to failed/intact interface (pInterface at
! rInterface = reposRadiusH(firstIntactZone)).
do i = firstIntactZone, numReposZones
  ! Elastic Stresses (Lame solution; inner BC = pInterface at rInterface)
  temp1 = (rInterface/reposRadius(i))**geomExponent
  radElasticStress(i) =  (pInterface-farfieldStress)*temp1    +farfieldStress
  tanElasticStress(i) = -(pInterface-farfieldStress)*temp1/mm1+farfieldStress
enddo


! interpolate pore pressure  accounting for zone size changes
! (pL is allocated and populated here for diagnostic / future use; not read by
!  the stress calculation itself.)
pL(firstIntactZone) = pInterface
pL(numReposZones+1) = reposPres(numreposZones)
do i = firstIntactZone+1,numreposZones
  pL(i) = reposPres(i-1)+reposFactor(i)*(reposPres(i)-reposPres(i-1))
enddo


!------------------------------------------------
! Seepage Stress, Total Stress and Effective Stress
!------------------------------------------------

! 2026 BC-revision: integration starts at the failed/intact interface
! rInterface with pressure pInterface. The first trapezoid is a special case
! spanning [rInterface, reposRadius(firstIntactZone)] -- a half-cell width --
! and uses pInterface at the lower limit. Subsequent trapezoids span between
! adjacent cell centers (width reposDRH(i)) as before.
!
! The (1 - fractionFluidized(i)) factor that was present here in the prior
! version has been removed: by construction every cell in the new integration
! range [firstIntactZone, numReposZones] has tensileFailureCompleted == FALSE
! and therefore fractionFluidized == 0 (fluidization can only initiate after
! tensile failure completes). Retaining the factor would falsely imply that
! partial fluidization could occur within the integration domain.

integral1 = 0.0
pff       = reposPres(numReposZones)

Do i = firstIntactZone, numReposZones

    pnew = reposPres(i)

    pPrime = pNew - pff

    if(i == firstIntactZone) then
      ! First trapezoid: from face (rInterface, pInterface) to cell center
      ! (reposRadius(i), reposPres(i)). Width = half cell = 0.5*reposDR(i).
      pPrimeR = 0.50d0*( (pInterface  -pff)*rInterface       **mm1 &
                        +(reposPres(i)-pff)*reposRadius(i)   **mm1 )
      drFirst = reposRadius(i) - rInterface       ! = 0.5*reposDR(i)
      integral1 = integral1 + pPrimeR*drFirst
    else
      ! Cell-center to cell-center trapezoid (unchanged form).
      pPrimeR = 0.50d0*( (reposPres(i-1)-pff)*reposRadius(i-1)**mm1 &
                        +(reposPres(i)  -pff)*reposRadius(i)  **mm1 )
      integral1 = integral1 + pPrimeR*reposDRH(i)
    endif

    temp1     = 1.0d0/reposRadius(i)**geomExponent

    !Seepage Stress
    radSeepageStress(i) = mm1*preFactor*integral1*temp1
    tanSeepageStress(i) = -preFactor*(integral1*temp1-pPrime)


    ! Total stresses
    radTotStress = radElasticStress(i) +radSeepageStress(i)
    tanTotStress = tanElasticStress(i) +tanSeepageStress(i)


    ! Effective stresses -
    ! DKR added BiotBeta, since it is included in seepage stress prefactor
    radEffStress(i) = (radTotStress-biotBeta*pNew)
    tanEffStress(i) = (tanTotStress-biotBeta*pNew)


	shearStress (i) = 0.5*DAbs(radEffStress(i)-tanEffStress(i))

enddo

! Failure only if repository is penetrated
if (repositoryPenetrated) then
    	

  ! check average of first N zones <= characteristic length, Lt
  i=firstIntactZone+1
  do while (reposRadiusH(i)-reposRadiusH(firstIntactZone) <= 0.9999*Lt)
   i = i+1
  enddo
  LtZone = i-1

  ! calculate average stress over characteristic failure length
  aveRadStress  = 0.0
  aveTanStress  = 0.0

  do i = firstIntactZone, LtZone
    aveRadStress  = aveRadStress  + radEffStress(i)
    aveTanStress  = aveTanStress  + tanEffStress(i)
  enddo

  temp1 = dble(LtZone-firstIntactZone+1)  !apg real(
  aveRadStress   = aveRadStress/temp1
  aveTanStress   = aveTanStress/temp1
  aveShrStress   = 0.5*DAbs(aveRadStress-aveTanStress)
  aveMeanStress  = (aveRadStress+mm1*aveTanStress)/geomExponent
  aveShrStrength = S0+aveMeanStress*Tan(frictionAngle)

  do i= firstIntactZone, numreposZones

    meanEffStress(i) = (radEffStress(i)+mm1*tanEffStress(i))/geomExponent
    shearStrength(i) = S0+meanEffStress(i)*Tan(frictionAngle)


    ! flag shear failed zones
    if ((i > LtZone .and. shearStress(i) > shearStrength(i)) .or. &
        (i <=LtZone .and. aveShrStress   > aveShrStrength  )) then
        shearFailed(i) = .TRUE.
    endif
  enddo


  do i= firstIntactZone, LtZone

    ! tensile failure initialization
    if ((.NOT.tensileFailureStarted(i)) .AND. &
        (runTime*TensileVelocity > reposRadius(i) - initialCavityRadius) ) then

      Ltemp1 = (surfaceFailureAllowed  .and. (-aveRadStress > tensileStrength))

      if ( Ltemp1 ) then

        tensileFailureStarted(i) = .TRUE.
        failStartTime(i) = runTime
		if (i == LtZone) then
		  surfaceFailureAllowed = .False.
		  fluidizationWaitZone = i
        endif

        ! test case #4 print out
		if (validationTestCase == 4) then
			stressSaveTime = runTime
			write (stressValidationFileID, '(A)') ''
			write (stressValidationFileID, '(A5, I5, A10, F10.5, A5)') 'zone', i,'failed at', runtime, 'sec'
		endif
      endif
    endif


    ! tensile failure in progress
    if ((tensileFailureStarted  (i)) .AND. &
        (.NOT.tensileFailureCompleted(i))) then

!JFS to not allow brief excursion over tensile to trigger ultimate failure
	  Ltemp1 = (-aveRadStress > tensileStrength)

	  if (Ltemp1) then
  	    fractionTensileFailed(i) = fractionTensileFailed(i) &
                                  +deltaTime/tensileFailureTime(i)
      endif

      if (fractionTensileFailed(i) > 1.0) then
        fractionTensileFailed  (i) = 1.0
        tensileFailureCompleted(i) = .TRUE.
        lastFailedZone = max(lastFailedZone, i)
      endif

    endif

  end do


  ! stress-related indices
  if (repositoryPenetrated) then
    do i = firstIntactZone, numReposZones
      if (tensileFailureCompleted(i)) then
        maxTensileFailedIndex = i
      endif
      if (shearFailed(i)) then
        maxShearFailedIndex = i

      endif
    end do
  endif

endif


return
end

!-----------------------------------------------------------------------------------------

