! @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
! @author David George
!
! @section DESCRIPTION
! Binds GeoClaw's Fortran routines to C.

! number of f-waves in the Problem (2 shock linearization = 2, augmented = 3)
#ifndef NUMBER_OF_FWAVES
#define NUMBER_OF_FWAVES 3
#endif

module c_bind_riemannsolvers
  implicit none

  contains

  ! C-interface to the augmented Riemann solver for the 1D shallow water equations.
  !   Remark: This routine is based on the respective Fortan interface (rpn2.f90)
  !           Wave-limiters are not implemented.
  subroutine c_bind_riemann_aug_JCP( i_maxIter,&
                                     i_variablesLeft,  i_variablesRight,&
                                     i_dryTol, i_g,&
                                     o_netUpdatesLeft, o_netUpdatesRight,&
                                     o_waveSpeeds &
#if AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
                                    ,o_eigenCoefficients &
#endif
                                     ) bind(c, name='c_bind_geoclaw_riemann_aug_JCP')
  !variable declaration
    !input
    integer                                   :: i_maxIter          !< maximum number of iterations over the Riemann problem.

    double precision, dimension(3)            :: i_variablesLeft,&  !< array containing variables of the left cell: $(h_l, hu_l, b_l)^T$
                                                 i_variablesRight   !< array containing variables of the right cell: $(h_r, hu_r, b_r)^T$

    double precision                          :: i_dryTol,&         !< dry tolerance
                                                 i_g                !< gravity

    !output
    double precision, dimension( 3 )                :: o_netUpdatesLeft,&   !< net-updates for the left cell
                                                       o_netUpdatesRight    !< net-updates for the right cell
    double precision, dimension( NUMBER_OF_FWAVES ) :: o_waveSpeeds         !< wave speeds/eigenvalues (negative/positive: waves travelling to the left/right cell)
#if AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
    double precision, dimension( NUMBER_OF_FWAVES ) :: o_eigenCoefficients  !< eigencoefficients computed in the riemann solution
#endif

    !local
    integer                                              :: l_waveNumber,&   !< wave number occuring in the Riemann solution (1: left wave, 2: (if defined) corrector wave, 3: right wave)
                                                            l_quantityNumber !< number of the quantity (h = 1, hu=2, hv=3)
    double precision, dimension(3, NUMBER_OF_FWAVES )    :: l_fWaves         !< f-waves
    double precision, dimension(3)                       :: l_wall           !< boundary conditions for the quantities
    double precision                                     :: l_hL,&           !< height, left cell
                                                            l_hR,&           !< height, right cell
                                                            l_huL,&          !< edge-normal momentum, left cell
                                                            l_huR,&          !< edge-normal momentum, right cell
                                                            l_hvL,&          !< edge-parallel momentum, left cell
                                                            l_hvR,&          !< edge-parallel momentum, right cell
                                                            l_bL,&           !< bathymetry, left cell
                                                            l_bR,&           !< bathymetry, right cell
                                                            l_uL,&           !< particle velocity, left cell in physical x-direction
                                                            l_uR,&           !< particle velocity, right cell in physical x-direction
                                                            l_vL,&           !< particle velocity, left cell in physical y-direction
                                                            l_vR,&           !< particle velocity, right cell in physical x-direction
                                                            l_hstar,&        !< height of the middle state
                                                            l_hstartest,&    !< maximum height h_l, h^*, h_r
                                                            l_uhat,&         !< Roe average of particle velocity
                                                            l_chat,&         !< wave speed relative to the fluid
                                                            l_phiL,&         !< normal momentum flux, left cell
                                                            l_phiR           !< normal momentum flux, right cell
    double precision                                    ::  l_sL,&           !< characteristic speed, left wave
                                                            l_sR,&           !< characteristic speed, right wave
                                                            l_sRoe1,&        !< Roe speed, left wave
                                                            l_sRoe2,&        !< Roe speed, right wave
                                                            l_sE1,&          !< Einfeldt speed, left wave
                                                            l_sE2,&          !< Einfeldt speed, right wave
                                                            l_s1m,&          !< Rarefaction speed, left side
                                                            l_s2m            !< Rarefaction speed, right side
    logical                                             ::  l_rare1,&        !< boolean whether the left wave is a rarefaction wave
                                                            l_rare2          !< boolean whether the right wave is a rarefaction wave

  ! copy cell values to local variables
  l_hL  = i_variablesLeft(1)
  l_huL = i_variablesLeft(2)
  l_bL  = i_variablesLeft(3)

  l_hR  = i_variablesRight(1)
  l_huR = i_variablesRight(2)
  l_bR  = i_variablesRight(3)

  ! set parallel momentums to zero (TODO: limiter)
  l_hvL = 0.d0
  l_hvR = 0.d0

  !******************************************************************
  !* necessary (changed) part of the GeoClaw subroutine rpn - start *
  !******************************************************************

  !reset wave speeds wave speed
  o_waveSpeeds = 0.d0;

  !reset net updates
  o_netUpdatesLeft  = 0.d0;
  o_netUpdatesRight = 0.d0;


  !Initialize Riemann problem for grid interface
  o_waveSpeeds = 0.d0
  l_fWaves     = 0.d0

  !zero (small) negative values if they exist
  if (l_hL.lt.0.d0) then
    l_hL = 0.d0; !rest of left variables will be set later
  endif

  if (l_hR.lt.0.d0) then
    l_hR = 0.d0; !rest of right variables will be set later
  endif

  !skip problem if in a completely dry area
  if (l_hL.le.i_dryTol.and.l_hR.le.i_dryTol) then
    return;
  endif

  !check for wet/dry boundary
  if (l_hR.gt.i_dryTol) then
    l_uR=l_huR/l_hR
    l_vR=l_hvR/l_hR
    l_phiR = 0.5d0*i_g*l_hR**2 + l_huR**2/l_hR
  else
    l_hR = 0.d0
    l_huR = 0.d0
    l_hvR = 0.d0
    l_uR = 0.d0
    l_vR = 0.d0
    l_phiR = 0.d0
  endif

  if (l_hL.gt.i_dryTol) then
    l_uL=l_huL/l_hL
    l_vL=l_hvL/l_hL
    l_phiL = 0.5d0*i_g*l_hL**2 + l_huL**2/l_hL
  else
    l_hL=0.d0
    l_huL=0.d0
    l_hvL=0.d0
    l_uL=0.d0
    l_vL=0.d0
    l_phiL = 0.d0
  endif

  !per default there is no wall
  l_wall(1) = 1.d0
  l_wall(2) = 1.d0
  l_wall(3) = 1.d0

  if (l_hR.le.i_dryTol) then
    call riemanntype(l_hL,l_hL,l_uL,-l_uL,l_hstar,l_s1m,l_s2m,l_rare1,l_rare2,1,i_dryTol,i_g)
    l_hstartest=max(l_hL,l_hstar)
    if (l_hstartest+l_bL.lt.l_bR) then !right state should become ghost values that mirror left for wall problem
      l_wall(2)=0.d0
      l_wall(3)=0.d0
      l_hR=l_hL
      l_huR=-l_huL
      l_bR=l_bL
      l_phiR=l_phiL
      l_uR=-l_uL
      l_vR=l_vL
    elseif (l_hL+l_bL.lt.l_bR) then
      l_bR=l_hL+l_bL
    endif
  elseif (l_hL.le.i_dryTol) then ! right surface is lower than left topo
    call riemanntype(l_hR,l_hR,-l_uR,l_uR,l_hstar,l_s1m,l_s2m,l_rare1,l_rare2,1,i_dryTol,i_g)
    l_hstartest=max(l_hR,l_hstar)
    if (l_hstartest+l_bR.lt.l_bL) then  !left state should become ghost values that mirror right
      l_wall(1)=0.d0
      l_wall(2)=0.d0
      l_hL=l_hR
      l_huL=-l_huR
      l_bL=l_bR
      l_phiL=l_phiR
      l_uL=-l_uR
      l_vL=l_vR
    elseif (l_hR+l_bR.lt.l_bL) then
      l_bL=l_hR+l_bR
    endif
  endif

  !determine wave speeds
  l_sL=l_uL-sqrt(i_g*l_hL) ! 1 wave speed of left state
  l_sR=l_uR+sqrt(i_g*l_hR) ! 2 wave speed of right state

  l_uhat=(sqrt(i_g*l_hL)*l_uL + sqrt(i_g*l_hR)*l_uR)/(sqrt(i_g*l_hR)+sqrt(i_g*l_hL)) ! Roe average
  l_chat=sqrt(i_g*0.5d0*(l_hR+l_hL)) ! Roe average
  l_sRoe1=l_uhat-l_chat ! Roe wave speed 1 wave
  l_sRoe2=l_uhat+l_chat ! Roe wave speed 2 wave

  l_sE1 = min(l_sL,l_sRoe1) ! Eindfeldt speed 1 wave
  l_sE2 = max(l_sR,l_sRoe2) ! Eindfeldt speed 2 wave

  !*******************
  !* call the solver *
  !*******************
  call riemann_aug_JCP(i_maxIter,3, NUMBER_OF_FWAVES, l_hL,l_hR,l_huL,l_huR,l_hvL,l_hvR,l_bL,l_bR,l_uL,l_uR,l_vL,l_vR,l_phiL,l_phiR,l_sE1,l_sE2,i_dryTol,i_g,o_waveSpeeds,&
#if AUGMENTED_RIEMANN_EIGEN_COEFFICIENTS
     o_eigenCoefficients, &
#endif
      l_fWaves)


  !eliminate ghost fluxes for wall
  do l_waveNumber=1,NUMBER_OF_FWAVES
    o_waveSpeeds(l_waveNumber)=o_waveSpeeds(l_waveNumber)*l_wall(l_waveNumber)
    do l_quantityNumber=1,3
      l_fWaves(l_quantityNumber,l_waveNumber)=l_fWaves(l_quantityNumber,l_waveNumber)*l_wall(l_waveNumber)
    enddo
  enddo

  !compute net updates
  do l_quantityNumber=1,3
    do  l_waveNumber=1,NUMBER_OF_FWAVES
      if (o_waveSpeeds(l_waveNumber).lt.0.d0) then
       o_netUpdatesLeft(l_quantityNumber)=o_netUpdatesLeft(l_quantityNumber) + l_fWaves(l_quantityNumber,l_waveNumber);
      elseif (o_waveSpeeds(l_waveNumber).gt.0.d0) then
       o_netUpdatesRight(l_quantityNumber)=o_netUpdatesRight(l_quantityNumber) + l_fWaves(l_quantityNumber,l_waveNumber);
      else
        o_netUpdatesLeft(l_quantityNumber)=o_netUpdatesLeft(l_quantityNumber) + .5d0*l_fWaves(l_quantityNumber,l_waveNumber);
        o_netUpdatesRight(l_quantityNumber)=o_netUpdatesRight(l_quantityNumber) + .5d0*l_fWaves(l_quantityNumber,l_waveNumber);
      endif
    enddo
  enddo

  !******************************************************
  !* necessary part of the GeoClaw subroutine rpn - end *
  !******************************************************

  !add second order correction terms
  !do l_waveNumber=1,NUMBER_OF_FWAVES
    ! TODO!
    !correctionL = SIGN(1.d0,waveSpeeds(l_waveNumber))
  !end do

  end subroutine c_bind_riemann_aug_JCP
end module
