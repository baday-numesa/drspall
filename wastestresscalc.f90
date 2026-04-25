!-----------------------------------------------------------------------------------------
!
! File WasteStressCalc.f90, contains routines for the stress calculation
!
! 2026 BC-revision: the radial effective-stress inner boundary condition is
! applied at the failed/intact interface using the local pore pressure from
! the gas-flow solution. Stresses are posed on the intact domain only;
! tensile-failed (disaggregated) waste transmits pore pressure isotropically
! but carries no deviatoric stress.
!
! Index conventions:
!   firstFailedZone : inner edge of the failed/fluidizing annulus. Advances
!                     on fluidization completion.
!   firstIntactZone : inner edge of the intact waste; anchor for the stress
!                     problem. Advances in Lt-batch-sized jumps on batch
!                     completion.
!   batchEndZone    : outer cell of the currently progressing Lt batch; 0
!                     when no batch is active (surfaceFailureAllowed=.TRUE.).
!   When no failed annulus exists, firstIntactZone == firstFailedZone and
!   pInterface reduces to cavityPres (the flow solver pins
!   reposPres(firstFailedZone) <- cavityPres).
!
! Tensile-failure cycle:
!   1. No active batch: the first Lt of intact rock [firstIntactZone, LtZone]
!      is the candidate batch. If -aveRadStress over that window exceeds
!      tensileStrength, every cell in the window is flagged
!      (tensileFailureStarted) and batchEndZone is set to LtZone.
!   2. Active batch: fractionTensileFailed(i) accumulates per its
!      tensileFailureTime for each cell in [firstIntactZone, batchEndZone],
!      gated on the batch's current aveRadStress exceeding tensileStrength.
!   3. When the outer cell of the batch reaches fractionTensileFailed >= 1,
!      the whole batch flips tensileFailureCompleted simultaneously,
!      firstIntactZone jumps to batchEndZone+1, and surfaceFailureAllowed
!      returns to .TRUE.
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
! Inner-BC pressure/radius for the intact-rock stress problem.
! firstIntactZone is persistent state updated by the tensile-failure batch
! logic below; no per-call rescan is performed.
!------------------------------------------------------------

! Inner mechanical boundary is at the inner face of the first intact cell.
rInterface = reposRadiusH(firstIntactZone)

! Pore pressure at the failed/intact interface. When firstIntactZone ==
! firstFailedZone there is no failed annulus and the flow solver has pinned
! reposPres(firstFailedZone) = cavityPres, so pInterface = cavityPres.
if (firstIntactZone == firstFailedZone) then
  pInterface = reposPres(firstFailedZone)
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

! Lame solution on the intact domain with inner BC pInterface at rInterface.
do i = firstIntactZone, numReposZones
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

! Trapezoidal integration of (p' * r^(m-1)) across the intact domain starting
! at (rInterface, pInterface). The first trapezoid spans the half-cell
! [rInterface, reposRadius(firstIntactZone)]; subsequent trapezoids span
! between adjacent cell centers with width reposDRH(i).

integral1 = 0.0
pff       = reposPres(numReposZones)

Do i = firstIntactZone, numReposZones

    pnew = reposPres(i)

    pPrime = pNew - pff

    if(i == firstIntactZone) then
      pPrimeR = 0.50d0*( (pInterface  -pff)*rInterface       **mm1 &
                        +(reposPres(i)-pff)*reposRadius(i)   **mm1 )
      drFirst = reposRadius(i) - rInterface       ! = 0.5*reposDR(i)
      integral1 = integral1 + pPrimeR*drFirst
    else
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

! Failure only if repository is penetrated and intact waste remains
if (repositoryPenetrated .and. firstIntactZone <= numReposZones) then

  ! Lt window:
  !   surfaceFailureAllowed -> no active batch, candidate window is the first
  !                            Lt of intact rock starting at firstIntactZone.
  !   else                  -> active batch, window is [firstIntactZone,
  !                            batchEndZone] (firstIntactZone is pinned for
  !                            the lifetime of the batch, batchEndZone is
  !                            fixed at initiation).
  if (surfaceFailureAllowed) then
    i = firstIntactZone + 1
    do while (i <= numReposZones+1)
      if (reposRadiusH(i)-reposRadiusH(firstIntactZone) > 0.9999d0*Lt) exit
      i = i + 1
    enddo
    LtZone = min(i-1, numReposZones)
  else
    LtZone = batchEndZone
  endif

  ! Group-averaged stresses over the Lt window.
  aveRadStress = 0.0d0
  aveTanStress = 0.0d0
  do i = firstIntactZone, LtZone
    aveRadStress = aveRadStress + radEffStress(i)
    aveTanStress = aveTanStress + tanEffStress(i)
  enddo
  temp1 = dble(LtZone - firstIntactZone + 1)
  aveRadStress   = aveRadStress/temp1
  aveTanStress   = aveTanStress/temp1
  aveShrStress   = 0.5d0*DAbs(aveRadStress - aveTanStress)
  aveMeanStress  = (aveRadStress + mm1*aveTanStress)/geomExponent
  aveShrStrength = S0 + aveMeanStress*Tan(frictionAngle)

  ! Shear failure flag: Lt-averaged criterion inside the window, point-wise
  ! criterion outside.
  do i = firstIntactZone, numReposZones
    meanEffStress(i) = (radEffStress(i) + mm1*tanEffStress(i))/geomExponent
    shearStrength(i) = S0 + meanEffStress(i)*Tan(frictionAngle)

    if ((i >  LtZone .and. shearStress(i) > shearStrength(i)) .or. &
        (i <= LtZone .and. aveShrStress   > aveShrStrength  )) then
      shearFailed(i) = .TRUE.
    endif
  enddo

  ! Tensile failure batch state machine.
  if (surfaceFailureAllowed) then

    ! No active batch. The Lt window is the candidate group; if the group-
    ! averaged radial effective stress is tensile enough, flag the whole
    ! group for failure and open a batch with batchEndZone = LtZone.
    if (-aveRadStress > tensileStrength) then
      do i = firstIntactZone, LtZone
        if (.NOT. tensileFailureStarted(i)) then
          tensileFailureStarted(i) = .TRUE.
          failStartTime(i)         = runTime

          if (validationTestCase == 4) then
            stressSaveTime = runTime
            write (stressValidationFileID, '(A)') ''
            write (stressValidationFileID, '(A5, I5, A10, F10.5, A5)') &
              'zone', i, 'failed at', runtime, 'sec'
          endif
        endif
      enddo
      batchEndZone          = LtZone
      surfaceFailureAllowed = .FALSE.
    endif

  else

    ! Active batch. Each cell's fractionTensileFailed ticks at its own
    ! tensileFailureTime while the group's aveRadStress remains above
    ! tensileStrength (prevents brief tensile excursions from completing
    ! the batch).
    if (-aveRadStress > tensileStrength) then
      do i = firstIntactZone, batchEndZone
        if (tensileFailureStarted(i) .AND. .NOT. tensileFailureCompleted(i)) then
          fractionTensileFailed(i) = fractionTensileFailed(i) &
                                   + deltaTime/tensileFailureTime(i)
        endif
      enddo
    endif

    ! tensileFailureTime scales with radius, so the outer cell is the last
    ! of the group to reach fraction == 1. When it does, the whole batch
    ! flips to completed simultaneously, firstIntactZone jumps past the
    ! batch, and surfaceFailureAllowed re-opens for the next cycle.
    if (fractionTensileFailed(batchEndZone) >= 1.0d0) then
      do i = firstIntactZone, batchEndZone
        fractionTensileFailed  (i) = 1.0d0
        tensileFailureCompleted(i) = .TRUE.
        lastFailedZone = max(lastFailedZone, i)
      enddo
      firstIntactZone       = batchEndZone + 1
      batchEndZone          = 0
      surfaceFailureAllowed = .TRUE.
    endif

  endif

  ! but-not-yet-fluidized annulus is included in the max-extent tallies.
  do i = firstFailedZone, numReposZones
    if (tensileFailureCompleted(i)) maxTensileFailedIndex = i
    if (shearFailed            (i)) maxShearFailedIndex   = i
  enddo

endif


return
end

!-----------------------------------------------------------------------------------------

