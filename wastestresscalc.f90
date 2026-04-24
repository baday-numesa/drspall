!-----------------------------------------------------------------------------------------
!
! File WasteStressCalc.f90, contains routines for the stress calculation
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
Logical thinShell, internalFailureAllowed, Ltemp1, Ltemp2, Ltemp3
Data thinShell /.FALSE. /
Data internalFailureAllowed/.TRUE./

! Initialization (m-1)
mm1 = geomExponent-1

! Seepage Force Constant
preFactor = biotBeta*(1.0-2.0*poissonsRatio)/(1.0-poissonsRatio)

! Solid Failure Constants
mu = Tan(frictionAngle)
S0 = 0.5*cohesion/(mu+DSqrt(mu**2+1.0))



!------------------------
!Calculate Elastic Stress
!------------------------

do i = firstIntactZone, numReposZones
  ! Elastic Stresses
  temp1 = (reposRadius(firstIntactZone)/reposRadius(i))**geomExponent
  radElasticStress(i) =  (cavityPres-farfieldStress)*temp1    +farfieldStress
  tanElasticStress(i) = -(cavityPres-farfieldStress)*temp1/mm1+farfieldStress
enddo


! interpolate pore pressure  accounting for zone size changes
pL(firstIntactZone) = 0.5*(cavityPres+reposPres(firstIntactZone))
pL(numReposZones+1) = reposPres(numreposZones)
do i = firstIntactZone+1,numreposZones
  pL(i) = reposPres(i-1)+reposFactor(i)*(reposPres(i)-reposPres(i-1))
enddo


!------------------------------------------------
! Seepage Stress, Total Stress and Effective Stress
!------------------------------------------------

integral1 = 0.0
Do i = firstIntactZone, numReposZones

    pnew = reposPres(i)

    pff    = reposPres(numReposZones)
    pPrime = pNew - pff
    !DKR trapezoidal integration of (p' * r^(m-1))
    if(i == firstIntactZone) then
      pPrimeR = 0.50d0*( (cavityPres  -pff)*reposRadius(i-1)**mm1 &
                        +(reposPres(i)-pff)*reposRadius(i)**mm1 )
    else
      pPrimeR = 0.50d0*( (reposPres(i-1)-pff)*reposRadius(i-1)**mm1 &
                        +(reposPres(i)  -pff)*reposRadius(i)**mm1 )
    endif

!JFS3 this is a simple stress relaxation and could help to avoid numerically-generated cascading failure
!DKR3 moved from radeffstress calc

    integral1 = integral1+pPrimeR*reposDRH(i)*(1.0 - fractionFluidized(i))

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

