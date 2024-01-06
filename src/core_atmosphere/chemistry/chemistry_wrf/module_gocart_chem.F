!=================================================================================================================
 module module_gocart_chem
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND
 use mpas_log


 implicit none
 private
 public:: gocart_chem


!TRACER INDEX:
!  1. DMS  = Dimethyl sulfide       = CH3SCH3
!  2. SO2  = Sulfur dioxide         = SO2
!  3. SO4  = Sulfate                = SO4
!  4. MSA  = Methane sulfonic acid  = CH3SO3H

 integer,parameter,public:: ndms = 1
 integer,parameter:: nso2 = 2
 integer,parameter:: nso4 = 3
 integer,parameter:: nmsa = 4

 real(kind=RKIND),parameter:: airmw  = 28.97
 real(kind=RKIND),parameter,public:: smw = 32.00
 real(kind=RKIND),parameter:: nh4_mfac = 1.375 ! increase sulf (output ond AOD only)
                                               ! to account for missing nh4
 real(kind=RKIND),parameter:: oc_mfac = 1.8    ! increase oc (output ond AOD nly)
                                               ! to account for carbon to organic mass


 contains


!=================================================================================================================
 subroutine gocart_chem(its,ite,kts,kte,grav,dt,xlon,xlat,area,backg_h2o2,backg_oh,backg_no3,qdms,qso2, &
                        qsulf,qmsa,qc,qi,cldfrac,rho,tt,zgrid)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte

 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in):: grav
 real(kind=RKIND),intent(in),dimension(its:ite):: xlon,xlat,area
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: cldfrac,rho,tt,qc,qi
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: backg_h2o2,backg_oh,backg_no3
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte+1):: zgrid


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qdms,qso2,qsulf,qmsa


!local variables and arrays:
 integer,parameter:: &
    imx = 1,  &
    lmx = 1,  &
    nmx = 4
 integer:: i,k

 real(kind=RKIND):: dz
 real(kind=RKIND),dimension(imx):: cossza
 real(kind=RKIND),dimension(imx,lmx):: airdens,airmass,cldf,temp
 real(kind=RKIND),dimension(imx,lmx):: oh,h2o2,no3
 real(kind=RKIND),dimension(imx,lmx):: chldms_oh,chldms_no3,chldms_x
 real(kind=RKIND),dimension(imx,lmx):: chlso2_oh,chlso2_aq
 real(kind=RKIND),dimension(imx,lmx):: chpso2,chpso4,chpmsa
 real(kind=RKIND),dimension(imx,lmx):: pmsa_dms,pso2_dms,pso4_so2
 real(kind=RKIND),dimension(imx,lmx,nmx):: tc

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_chem:')


!do i = its,ite
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r $r $r $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/xlon(i),xlat(i),area(i), &
!                          backg_h2o2(i,k),backg_oh(i,k),backg_no3(i,k),qdms(i,k),qso2(i,k),qsulf(i,k),qmsa(i,k)/))
!      call mpas_log_write('$i $i $r $r $r $r $r $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/qdms(i,k),qso2(i,k), &
!                           qsulf(i,k),qmsa(i,k),qc(i,k),qi(i,k),cldfrac(i,k),rho(i,k),tt(i,k),zgrid(i,k)/))
!      call mpas_log_write('$i $i $r $r $r',intArgs=(/i,k/),realArgs=(/backg_h2o2(i,k),backg_oh(i,k), &
!                          backg_no3(i,k)/))
!   enddo
!   call mpas_log_write(' ')
!enddo

 do i = its,ite
    do k = kts,kte

       dz = zgrid(i,k+1) - zgrid(i,k)

       !--- initialization of oh, h2o2, and no3:
       oh(imx,lmx)   = backg_oh(i,k)
       h2o2(imx,lmx) = backg_h2o2(i,k)
       no3(imx,lmx)  = backg_no3(i,k)

       !--- initialization of the air density, air mass, cloud fraction, and temperature:
       cldf(imx,lmx) = 1.0
       temp(imx,lmx) = tt(i,k)
       airdens(imx,lmx) = rho(i,k)
       airmass(imx,lmx) = airdens(imx,lmx)*dz*area(i)

       !--- assume that the cosine of the solar zenith angle equal 1:
       cossza(imx) = 1.0

       !--- initialization of chemical species (check the units for qdms,qso2,qsulf,qmsa): technically, I
       !    should loop over the length of imx and lmx.
       chldms_oh(imx,lmx)  = 0._RKIND
       chldms_no3(imx,lmx) = 0._RKIND
       chldms_x(imx,lmx)   = 0._RKIND
       chlso2_oh(imx,lmx)  = 0._RKIND
       chlso2_aq(imx,lmx)  = 0._RKIND

       chpso2(imx,lmx)     = 0._RKIND
       chpso4(imx,lmx)     = 0._RKIND
       chpmsa(imx,lmx)     = 0._RKIND

       pmsa_dms(imx,lmx)   = 0._RKIND
       pso2_dms(imx,lmx)   = 0._RKIND
       pso4_so2(imx,lmx)   = 0._RKIND

!      tc(imx,lmx,ndms) = qdms(i,k) * 1.e-6
!      tc(imx,lmx,nso2) = qso2(i,k) * 1.e-6
!      tc(imx,lmx,nso4) = qsulf(i,k)* 1.e-6
!      tc(imx,lmx,nmsa) = qmsa(i,k) * 1.e-6
       !--- all species are in mol/mol, so no need to multiply tc by 1.e-6:
       tc(imx,lmx,ndms) = qdms(i,k)
       tc(imx,lmx,nso2) = qso2(i,k)
       tc(imx,lmx,nso4) = qsulf(i,k)
       tc(imx,lmx,nmsa) = qmsa(i,k)


       !--- call to _dms, _so2, _so4, and _msa subroutines to compute production and loss of dms, so2,
       !    so4, and msa:
       call chem_dms(i,k,imx,lmx,nmx,dt,cossza,temp,airdens,airmass,oh,no3,tc,chldms_oh,chldms_no3, &
                     chldms_x,chpso2,chpmsa,pso2_dms,pmsa_dms)

       call chem_so2(i,k,imx,lmx,nmx,dt,cossza,temp,airdens,airmass,cldf,oh,h2o2,pso2_dms,tc, &
                     chlso2_oh,chlso2_aq,chpso4,pso4_so2)

       call chem_so4(i,k,imx,lmx,nmx,dt,pso4_so2,tc)

       call chem_msa(i,k,imx,lmx,nmx,dt,pmsa_dms,tc)

!      qdms(i,k)  = tc(imx,lmx,ndms) * 1.e6
!      qso2(i,k)  = tc(imx,lmx,nso2) * 1.e6
!      qsulf(i,k) = tc(imx,lmx,nso4) * 1.e6
!      qmsa(i,k)  = tc(imx,lmx,nmsa) * 1.e6
       !--- all species are in mol/mol, so no need to multiply tc by 1.e6:
       qdms(i,k)  = tc(imx,lmx,ndms)
       qso2(i,k)  = tc(imx,lmx,nso2)
       qsulf(i,k) = tc(imx,lmx,nso4)
       qmsa(i,k)  = tc(imx,lmx,nmsa)

    enddo
 enddo

 call mpas_log_write('--- end subroutine gocart_chem:')
 
 end subroutine gocart_chem

!=================================================================================================================
 subroutine chem_so2(i1,k1,imx,lmx,nmx,dt1,cossza,tmp,airden,airmas,cldf,oh,h2o2,pso2_dms,tc, &
                     chlso2_oh,chlso2_aq,chpso4,pso4_so2)
!=================================================================================================================

! ****************************************************************************
! *                                                                          *
! *  This is SO2 chemistry subroutine.                                       *
! *                                                                          *
! *  SO2 production:                                                         *
! *    DMS + OH, DMS + NO3 (saved in CHEM_DMS)                               *
! *                                                                          *
! *  SO2 loss:                                                               *
! *    SO2 + OH  -> SO4                                                      *
! *    SO2       -> drydep (NOT USED IN WRF/CHEM                             *
! *    SO2 + H2O2 or O3 (aq) -> SO4                                          *
! *                                                                          *
! *  SO2 = SO2_0 * exp(-bt)                                                  *
! *      + PSO2_DMS/bt * [1-exp(-bt)]                                        *
! *    where b is the sum of the reaction rate of SO2 + OH and the dry       *
! *    deposition rate of SO2, PSO2_DMS is SO2 production from DMS in        *
! *    MixingRatio/timestep.                                                 *
! *                                                                          *
! *  If there is cloud in the gridbox (fraction = fc), then the aqueous      *
! *  phase chemistry also takes place in cloud. The amount of SO2 oxidized   *
! *  by H2O2 in cloud is limited by the available H2O2; the rest may be      *
! *  oxidized due to additional chemistry, e.g, reaction with O3 or O2       *
! *  (catalyzed by trace metal).                                             *
! *                                                                          *
! ****************************************************************************

!input arguments:
 integer,intent(in):: i1,k1
 integer,intent(in):: imx,lmx,nmx

 real(kind=RKIND),intent(in):: dt1
 real(kind=RKIND),intent(in),dimension(imx,lmx):: cossza
 real(kind=RKIND),intent(in),dimension(imx,lmx):: airden,airmas,cldf,tmp
 real(kind=RKIND),intent(in),dimension(imx,lmx):: oh
 real(kind=RKIND),intent(in),dimension(imx,lmx):: pso2_dms


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: h2o2
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: chlso2_oh,chlso2_aq
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: chpso4
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: pso4_so2
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc


!local variables and arrays:
 integer:: i,l

 real(kind=RKIND),parameter:: ki = 1.5e-12
 real(kind=RKIND),parameter:: f = 1000./airmw*6.022d23*1.0d-6 !factor to convert kg of air/m3 to molecules/cm3.

 real(kind=RKIND):: l1,l2,ld,k0,k,kk,m
 real(kind=RKIND):: f1,fc,rk,rk1,rk2,rkt,so20,so2,so2_cd,tk

!-----------------------------------------------------------------------------------------------------------------

 do l = 1,lmx
    do i = 1,imx

       so20 = tc(i,l,nso2)

       ! RK1: SO2 + OH(g), in s-1
       tk = tmp(i,l)
       k0 = 3.0e-31 * (300.0/tk)**3.3
       m  = airden(i,l) * f
       kk = k0 * m / ki
       f1 = ( 1.0 + ( log10(kk) )**2 )**(-1)

       ! mozech: oh is in mol/mol
       ! convert to molecules/cm3
       rk1 = ( k0 * m / (1.0 + kk) ) * 0.6**f1 * oh(i,l)*airden(i,l) * f
       rk2 = 0._RKIND

       rk  = (rk1 + rk2)
       rkt =  rk * dt1

! ****************************************************************************
! *  Update SO2 concentration after gas phase chemistry and deposition.      *
! ****************************************************************************

       if(rk > 0.0) then
          so2_cd = so20 * exp(-rkt) + pso2_dms(i,l) * (1.0 - exp(-rkt)) / rkt
          l1     = (so20 - so2_cd + pso2_dms(i,l)) * rk1/rk
          if(l == 1) then
             ld = (so20 - so2_cd + pso2_dms(i,l)) * rk2/rk
          else
             ld    = 0.0
          endif
       else
          so2_cd = so20
          l1 = 0.0
       endif

! ****************************************************************************
! *  Update SO2 concentration after cloud chemistry.                         *
! *  SO2 chemical loss rate  = SO4 production rate (MixingRatio/timestep).   *
! ****************************************************************************

       ! Cloud chemistry (above 258K):
       fc = cldf(i,l)
       if(fc > 0.0 .and. so2_cd > 0.0 .and. tk > 258.0) then
          
          if(so2_cd > h2o2(i,l)) then
             fc = fc * (h2o2(i,l)/so2_cd)
             h2o2(i,l) = h2o2(i,l) * (1.0 - cldf(i,l))
          else
             h2o2(i,l) = h2o2(i,l) * (1.0 - cldf(i,l)*so2_cd/h2o2(i,l))
          endif
          so2 = so2_cd * (1.0 - fc)
          ! Aqueous phase SO2 loss rate (MixingRatio/timestep):
          l2  = so2_cd * fc
       else
          so2 = so2_cd
          l2 = 0.0
       endif

       so2 = max(so2, 1.0e-32)
       tc(i,l,nso2) = so2

! ****************************************************************************
! *  SO2 chemical loss rate  = SO4 production rate (MixingRatio/timestep).   *
! ****************************************************************************

       pso4_so2(i,l) = max(0.0,l1 + l2)

       !    ---------------------------------------------------------------
       !    DIAGNOSTICS:      SO2 gas-phase loss       (kgS/timestep)
       !                      SO2 aqueous-phase loss   (kgS/timestep)
       !                      SO2 dry deposition loss  (kgS/timestep)
       !                      SO4 production           (kgS/timestep)
       !    ---------------------------------------------------------------
       chlso2_oh(i,l) = chlso2_oh(i,l) + l1 * airmas(i,l) / airmw * smw
       chlso2_aq(i,l) = chlso2_aq(i,l) + l2 * airmas(i,l) / airmw * smw

       chpso4(i,l) = chpso4(i,l) + pso4_so2(i,l) * airmas(i,l) / airmw * smw

    enddo
 enddo

 end subroutine chem_so2

!=================================================================================================================
 subroutine chem_so4(ii,kk,imx,lmx,nmx,dt1,pso4_so2,tc)
!=================================================================================================================

!input arguments:
 integer,intent(in):: ii,kk
 integer,intent(in):: imx,lmx,nmx

 real(kind=RKIND),intent(in):: dt1
 real(kind=RKIND),intent(in),dimension(imx,lmx):: pso4_so2


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc


!local variables and arrays:
 integer:: i,l
 real(kind=RKIND):: so40,so4

!-----------------------------------------------------------------------------------------------------------------

 do l = 1,lmx
    do i = 1,imx

       so40 = tc(i,l,nso4)
       so4 = so40 + pso4_so2(i,l)
       so4 = max(so4,1.0e-32)

       tc(i,l,nso4) = so4

    enddo
 enddo

 end subroutine chem_so4

!=================================================================================================================
 subroutine chem_dms(ii,kk,imx,lmx,nmx,dt1,cossza,tmp,airden,airmas,oh,xno3,tc,chldms_oh,chldms_no3,chldms_x, &
                     chpso2,chpmsa,pso2_dms,pmsa_dms)
!=================================================================================================================

! ****************************************************************************
! *                                                                          *
! *  This is DMS chemistry subroutine.                                       *
! *                                                                          *
! *  R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel    *
! *         k1 = { 1.7e-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }  *
! *         a = 0.75, b = 0.25                                               *
! *                                                                          *
! *  R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel *
! *         k2 = 1.2e-11*exp(-260/T)                                         *
! *                                                                          *
! *     DMS_OH = DMS0 * exp(-(r1+r2)*DT1)                                    *
! *         where DMS0 is the DMS concentration at the beginning,            *
! *         r1 = k1*[OH], r2 = k2*[OH].                                      *
! *                                                                          *
! *  R3:    DMS + NO3 ->   SO2 + ...                                         *
! *         k3 = 1.9e-13*exp(500/T)                                          *
! *                                                                          *
! *     DMS = DMS_OH * exp(-r3*DT1)                                          *
! *         where r3 = k3*[NO3].                                             *
! *                                                                          *
! *  R4:    DMS + X   ->   SO2 + ...                                         *
! *         assume to be at the rate of DMS+OH and DMS+NO3 combined.         *
! *                                                                          *
! *  The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved    *
! *  for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They    *
! *  are in unit of MixingRatio/timestep.                                    *
! *                                                                          *
! ****************************************************************************


!input arguments:
 integer,intent(in):: ii,kk
 integer,intent(in):: imx,lmx,nmx

 real(kind=RKIND),intent(in):: dt1
 real(kind=RKIND),intent(in),dimension(imx):: cossza
 real(kind=RKIND),intent(in),dimension(imx,lmx):: airden,airmas,tmp
 real(kind=RKIND),intent(in),dimension(imx,lmx):: oh,xno3


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: chldms_oh,chldms_no3,chldms_x
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: chpso2,chpmsa
 real(kind=RKIND),intent(inout),dimension(imx,lmx):: pso2_dms,pmsa_dms
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc


!local variables and arrays:
 integer:: i,l

 real(kind=RKIND),parameter:: fx = 1.0
 real(kind=RKIND),parameter:: a = 0.75
 real(kind=RKIND),parameter:: b = 0.25
!From D4: only 0.8 efficiency, also some goes to DMSO and lost.
!So we assume 0.75 efficiency for DMS addtion channel to form
!products.
 real(kind=RKIND),parameter:: eff = 1.0
 real(kind=RKIND),parameter:: f = 1000./airmw*6.022d23*1.0d-6 !factor to convert kg of air/m3 to molecules/cm3.

 real(kind=RKIND):: dms0,dms_oh,dms,o2,tk
 real(kind=RKIND):: rk1,rk2,rk3
 real(kind=RKIND):: xoh,xn3,xx

!-----------------------------------------------------------------------------------------------------------------

 do l = 1,lmx
    do i = 1,imx

       tk = tmp(i,l)
       o2 = airden(i,l)*f*0.21
       dms0 = tc(i,l,ndms)

! ****************************************************************************
! *  (1) DMS + OH:  RK1 - addition channel;  RK2 - abstraction channel.      *
! ****************************************************************************

       rk1 = 0._RKIND
       rk2 = 0._RKIND
       rk3 = 0._RKIND

       if(oh(i,l) > 0.0) then
          !convert oh from  mol/mol to molecules/cm3:
          rk1 = (1.7e-42 * exp(7810.0/tk) * o2) &
              / (1.0 + 5.5e-31 * exp(7460.0/tk) * o2 ) * oh(i,l) * airden(i,l) * f
          rk2 = 1.2e-11 * exp(-260.0/tk) * oh(i,l) * airden(i,l) * f
       endif

! ****************************************************************************
! *  (2) DMS + NO3 (only happens at night):                                  *
! ****************************************************************************

       if(cossza(i) <= 0.0) then
          !convert xno3 from volume mixing ratio to molecules/cm3:
          rk3 = 1.9d-13 * exp(500.0/tk) * xno3(i,l) * airden(i,l) * f
       endif

! ****************************************************************************
! *  Update DMS concentrations after reaction with OH and NO3, and also      *
! *  account for DMS + X assuming at a rate as (DMS+OH)*Fx in the day and    *
! *  (DMS+NO3)*Fx at night:                                                  *
! *       DMS_OH       :  DMS concentration after reaction with OH           *
! *       DMS          :  DMS concentration after reaction with NO3          *
! *                           (min(DMS) = 1.0E-32)                           *
! ****************************************************************************

       dms_oh = dms0   * exp( -(rk1 + rk2) * fx * dt1 )
       dms    = dms_oh * exp( -(rk3) * fx * dt1 )
       dms    = max(dms,1.0e-32)

! ****************************************************************************
! *  Save SO2 and MSA production from DMS oxidation                          *
! *  (in MixingRatio/timestep):                                              *
! *                                                                          *
! *  SO2 is formed in DMS + OH addition (0.85) and abstraction (1.0)         *
! *      channels as well as DMS + NO3 reaction.  We also assume that        *
! *      SO2 yield from DMS + X is 1.0.                                      *
! *  MSA is formed in DMS + OH addition (0.15) channel.                      *
! ****************************************************************************

       if((rk1 + rk2) == 0.0) then
           pmsa_dms(i,l) = 0._RKIND
       else
!         pmsa_dms(i,l) = (dms0 - dms_oh) * b*rk1/((rk1+rk2)*fx)
          pmsa_dms(i,l) = max(0._RKIND,(dms0 - dms_oh) * b*rk1/((rk1+rk2) * fx) * eff)
       endif
!      pso2_dms(i,l) = (dms0 - dms - pmsa_dms(i,j,l)/eff) * eff
       pso2_dms(i,l) = max(0.0D0,dms0 - dms - pmsa_dms(i,l))

       !    ------------------------------------------------------------
       !    DIAGNOSTICS:      DMS loss       (kgS/timstep)
       !                      SO2 production (kgS/timestep)
       !                      MSA production (kgS/timestep)
       !    ------------------------------------------------------------
       xoh  = (dms0   - dms_oh) / fx  * airmas(i,l)/airmw*smw
       xn3  = (dms_oh - dms)    / fx  * airmas(i,l)/airmw*smw
       xx   = (dms0 - dms) * airmas(i,l)/airmw*smw - xoh - xn3

       chldms_oh (i,l) = chldms_oh(i,l) + xoh
       chldms_no3(i,l) = chldms_no3(i,l) + xn3
       chldms_x(i,l)   = chldms_x(i,l) + xx

       chpso2(i,l) = chpso2(i,l) + pso2_dms(i,l) * airmas(i,l) / airmw * smw
       chpmsa(i,l) = chpmsa(i,l) + pmsa_dms(i,l) * airmas(i,l) / airmw * smw

    enddo
 enddo

 end subroutine chem_dms

!=================================================================================================================
 subroutine chem_msa(ii,kk,imx,lmx,nmx,dt1,pmsa_dms,tc)
!=================================================================================================================

! ****************************************************************************
! *                                                                          *
! *  This is MSA chemistry subroutine.                                       *
! *                                                                          *
! *  The Only production is from DMS oxidation (save in CHEM_DMS), and the   *
! *  only loss is dry depsition here.  Wet deposition will be treated in     *
! *  WETDEP subroutine.                                                      *
! *                                                                          *
! *  MSA = MSA_0 * exp(-dt) + PMSA_DMS/kt * (1.-exp(-kt))                    *
! *    where k = dry deposition.                                             *
! *                                                                          *
! ****************************************************************************


!input arguments:
 integer,intent(in):: ii,kk
 integer,intent(in):: imx,lmx,nmx

 real(kind=RKIND),intent(in):: dt1
 real(kind=RKIND),intent(in),dimension(imx,lmx):: pmsa_dms

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc


!local variables and arrays:
 integer:: i,l

 real(kind=RKIND):: msa0,msa

!-----------------------------------------------------------------------------------------------------------------

 do l = 1,lmx                                                               
    do i = 1,imx
       msa0 = tc(i,l,nmsa)

       msa  = msa0 + pmsa_dms(i,l)
       msa  = max(msa,1.0e-32)

       tc(i,l,nmsa) = msa
    enddo
 enddo

 end subroutine chem_msa

!=================================================================================================================
 end module module_gocart_chem
!=================================================================================================================

