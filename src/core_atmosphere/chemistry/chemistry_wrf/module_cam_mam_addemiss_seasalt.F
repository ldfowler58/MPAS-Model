!=================================================================================================================
 module module_cam_mam_addemiss_seasalt
 use mpas_kind_types,only: RKIND
 use mpas_log


 implicit none
 private
 public:: cam_mam_mosaic_seasalt_emiss


 contains


!=================================================================================================================
 subroutine cam_mam_mosaic_seasalt_emiss(its,ite,kts,kte,ktep1,dt,grav,qseas1,qseas2,qseas3,u10,v10,xland, &
                                    rho,zgrid,emis_qseas1,emis_qseas2,emis_qseas3,emis_nseas1,emis_nseas2, &
                                    emis_nseas3)
!=================================================================================================================
!
! adds seasalt emissions for mosaic aerosol species
! (i.e., seasalt emissions tendencies over time dtstep are applied 
! to the aerosol mixing ratios)
!

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte,ktep1

 real(kind=RKIND),intent(in):: dt,grav
 real(kind=RKIND),intent(in),dimension(its:ite):: xland 
 real(kind=RKIND),intent(in),dimension(its:ite):: u10,v10
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho
 real(kind=RKIND),intent(in),dimension(its:ite,kts:ktep1):: zgrid


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: & 
    qseas1,& ! sea-salt mixing ratio (accumulation mode).
    qseas2,& ! sea-salt mixing ratio (aitken mode).
    qseas3   ! sea-salt mixing ratio (coarse mode).


!--- output arguments:
 real(kind=RKIND),intent(out),dimension(its:ite):: &
    emis_qseas1, &!
    emis_qseas2, &!
    emis_qseas3, &!
    emis_nseas1, &!
    emis_nseas2, &!
    emis_nseas3   !


!--- local variables and arrays:
 integer,parameter:: ntot_amode = 3
 integer,parameter:: modeptr_accum  = 1
 integer,parameter:: modeptr_aitken = 2
 integer,parameter:: modeptr_coarse = 3
 integer:: i,itype,k,kk,n

 real(kind=RKIND):: dum,dumdlo,dumdhi
 real(kind=RKIND),dimension(ntot_amode):: dlo_seas_emiss
 real(kind=RKIND),dimension(ntot_amode):: dhi_seas_emiss
 real(kind=RKIND),dimension(ntot_amode):: ssemfact_numb
 real(kind=RKIND),dimension(ntot_amode):: ssemfact_mass

 real(kind=RKIND):: dumoceanfrac,dumspd10,dz
 real(kind=RKIND):: factaa,factbb

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine cam_mam_mosaic_seasalt_emiss:')

 do i = its,ite
    emis_qseas1(i) = 0._RKIND
    emis_qseas2(i) = 0._RKIND
    emis_qseas3(i) = 0._RKIND
    emis_nseas1(i) = 0._RKIND
    emis_nseas2(i) = 0._RKIND
    emis_nseas3(i) = 0._RKIND
 enddo


!MODAL_AERO_3MODE (cut sizes are 0.02, 0.08, 1.0, 10.0 um)::
 dlo_seas_emiss(modeptr_aitken) = 0.02e-4
 dhi_seas_emiss(modeptr_aitken) = 0.08e-4
 dlo_seas_emiss(modeptr_accum ) = 0.08e-4
 dhi_seas_emiss(modeptr_accum ) = 1.00e-4
 dlo_seas_emiss(modeptr_coarse) = 1.00e-4
 dhi_seas_emiss(modeptr_coarse) = 10.0e-4


!compute emission factors for each size bin (limit emissions to dp > 0.1 micrometer):
 ssemfact_mass(:) = 0.0
 ssemfact_numb(:) = 0.0

 do itype = 1,ntot_amode
    dumdlo = max(dlo_seas_emiss(itype),0.1e-4)
    dumdhi = max(dhi_seas_emiss(itype),0.1e-4)
    if (dumdlo >= dumdhi) cycle
    call seasalt_emitfactors_1bin(1,dumdlo,dumdhi,ssemfact_numb(itype),dum,ssemfact_mass(itype))

    !convert mass emissions factor from (g/m2/s) to (kg/m2/s)
    ssemfact_mass(itype) = ssemfact_mass(itype)*1.0e-3
 enddo
!do itype = 1,ntot_amode
!   dumdlo = max(dlo_seas_emiss(itype),0.1e-4)
!   dumdhi = max(dhi_seas_emiss(itype),0.1e-4)
!   call mpas_log_write('$i $r $r $r $r',intArgs=(/itype/),realArgs=(/dumdlo,dumdhi,ssemfact_numb(itype), &
!                       ssemfact_mass(itype)/))
!enddo
!call mpas_log_write(' ')


!apply seasalt emissions:
 k = kts
 do i = its,ite
    if(xland(i) < 1.5) cycle
    dumoceanfrac = 1. !fraction of grid i,j that is salt water
    dumspd10 = dumoceanfrac*((u10(i)*u10(i) + v10(i)*v10(i))**(0.5*3.41))

!   factaa is (s*m2/kg-air)
!   factaa*ssemfact_mass(n) is (s*m2/kg-air)*(kg/m2/s) = kg/kg-air
!   factaa*ssemfact_numb(n) is (s*m2/kg-air)*( #/m2/s) =  #/kg-air
    dz = zgrid(i,k+1) - zgrid(i,k) 
    factaa = dt/(dz*rho(i,k))
    factbb = factaa*dumspd10
    qseas1(i,k) = qseas1(i,k) + factbb*ssemfact_mass(1) ! accumulation mode.
    qseas2(i,k) = qseas2(i,k) + factbb*ssemfact_mass(2) ! aitken mode.
    qseas3(i,k) = qseas3(i,k) + factbb*ssemfact_mass(3) ! coarse mode.

    emis_qseas1(i) = factbb*ssemfact_mass(1)/dt
    emis_qseas2(i) = factbb*ssemfact_mass(2)/dt
    emis_qseas3(i) = factbb*ssemfact_mass(3)/dt
    emis_nseas1(i) = factbb*ssemfact_numb(1)/dt
    emis_nseas2(i) = factbb*ssemfact_numb(2)/dt
    emis_nseas3(i) = factbb*ssemfact_numb(3)/dt
!   call mpas_log_write('$i $r $r $r $r $r $r',intArgs=(/i/),realArgs=(/emis_qseas1(i),emis_qseas2(i), &
!                       emis_qseas3(i),emis_nseas1(i),emis_nseas2(i),emis_nseas3(i)/))
 enddo

 call mpas_log_write('--- end subroutine cam_mam_mosaic_seasalt_emiss:')

 end subroutine cam_mam_mosaic_seasalt_emiss

!=================================================================================================================
 subroutine seasalt_emitfactors_1bin(ireduce_smallr_emit,dpdrylo_cm,dpdryhi_cm,emitfact_numb, &
                                     emitfact_surf,emitfact_mass)
!=================================================================================================================
!c
!c   computes seasalt emissions factors for a specifed 
!c   dry particle size range
!c      dpdrylo_cm  = lower dry diameter (cm)
!c      dpdryhi_cm  = upper dry diameter (cm)
!c
!c   number and mass emissions are then computed as
!c      number   emissions (#/m2/s)   == emitfact_numb * (spd10*3.41)
!c      dry-sfc  emissions (cm2/m2/s) == emitfact_surf * (spd10*3.41)
!c      dry-mass emissions (g/m2/s)   == emitfact_mass * (spd10*3.41)
!c
!c   where spd10 = 10 m windspeed in m/s
!c
!c   uses bubble emissions formula (eqn 5a) from 
!c      Gong et al. [JGR, 1997, p 3805-3818]
!c
!c   *** for rdry < rdry_star, this formula overpredicts emissions.
!c      A strictly ad hoc correction is applied to the formula,
!c      based on sea-salt size measurements of
!c      O'Dowd et al. [Atmos Environ, 1997, p 73-80]
!c
!c   *** the correction is only applied when ireduce_smallr_emit > 0
!c
        implicit none

!c   subr arguments
        integer ireduce_smallr_emit
        real dpdrylo_cm, dpdryhi_cm,                            &
                emitfact_numb, emitfact_surf, emitfact_mass

!c   local variables
        integer isub_bin, nsub_bin

        real alnrdrylo
        real drydens, drydens_43pi_em12, x_4pi_em8
        real dum, dumadjust, dumb, dumexpb
        real dumsum_na, dumsum_ma, dumsum_sa
        real drwet, dlnrdry
        real df0drwet, df0dlnrdry, df0dlnrdry_star
        real relhum
        real rdry, rdrylo, rdryhi, rdryaa, rdrybb
        real rdrylowermost, rdryuppermost, rdry_star
        real rwet, rwetaa, rwetbb
        real rdry_cm, rwet_cm
        real sigmag_star
        real xmdry, xsdry

        real pi
        parameter (pi = 3.1415936536)

!c   c1-c4 are constants for seasalt hygroscopic growth parameterization
!c   in Eqn 3 and Table 2 of Gong et al. [1997]
        real c1, c2, c3, c4, onethird
        parameter (c1 = 0.7674)
        parameter (c2 = 3.079)
        parameter (c3 = 2.573e-11)
        parameter (c4 = -1.424)
        parameter (onethird = 1.0/3.0)


!c   dry particle density (g/cm3)
        drydens = 2.165
!c   factor for radius (micrometers) to mass (g)
        drydens_43pi_em12 = drydens*(4.0/3.0)*pi*1.0e-12
!c   factor for radius (micrometers) to surface (cm2)
        x_4pi_em8 = 4.0*pi*1.0e-8
!c   bubble emissions formula assume 80% RH
        relhum = 0.80

!c   rdry_star = dry radius (micrometers) below which the
!c   dF0/dr emission formula is adjusted downwards
        rdry_star = 0.1
        if (ireduce_smallr_emit .le. 0) rdry_star = -1.0e20
!c   sigmag_star = geometric standard deviation used for
!c   rdry < rdry_star
        sigmag_star = 1.9

!c   initialize sums
        dumsum_na = 0.0
        dumsum_sa = 0.0
        dumsum_ma = 0.0

!c   rdrylowermost, rdryuppermost = lower and upper 
!c   dry radii (micrometers) for overall integration
        rdrylowermost = dpdrylo_cm*0.5e4
        rdryuppermost = dpdryhi_cm*0.5e4

!c
!c   "section 1"
!c   integrate over rdry > rdry_star, where the dF0/dr emissions 
!c   formula is applicable
!c   (when ireduce_smallr_emit <= 0, rdry_star = -1.0e20,
!c   and the entire integration is done here)
!c
        if (rdryuppermost .le. rdry_star) goto 2000

!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = max( rdrylowermost, rdry_star )
        rdryhi = rdryuppermost

        nsub_bin = 1000

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

!c   compute rdry, rwet (micrometers) at lowest size
        rdrybb = exp( alnrdrylo )
        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

        do 1900 isub_bin = 1, nsub_bin

!c   rdry, rwet at sub_bin lower boundary are those
!c   at upper boundary of previous sub_bin
        rdryaa = rdrybb
        rwetaa = rwetbb

!c   compute rdry, rwet (micrometers) at sub_bin upper boundary
        dum = alnrdrylo + isub_bin*dlnrdry
        rdrybb = exp( dum )

        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

!c   geometric mean rdry, rwet (micrometers) for sub_bin
        rdry = sqrt(rdryaa * rdrybb)
        rwet = sqrt(rwetaa * rwetbb)
        drwet = rwetbb - rwetaa

!c   xmdry is dry mass in g
        xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
        xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumb is "B" in Gong's Eqn 5a
!c   df0drwet is "dF0/dr" in Gong's Eqn 5a
        dumb = ( 0.380 - log10(rwet) ) / 0.650
        dumexpb = exp( -dumb*dumb)
        df0drwet = 1.373 * (rwet**(-3.0)) *                     &
                (1.0 + 0.057*(rwet**1.05)) *                    &
                (10.0**(1.19*dumexpb))

        dumsum_na = dumsum_na + drwet*df0drwet
        dumsum_ma = dumsum_ma + drwet*df0drwet*xmdry
        dumsum_sa = dumsum_sa + drwet*df0drwet*xsdry

1900    continue


!c
!c   "section 2"
!c   integrate over rdry < rdry_star, where the dF0/dr emissions 
!c   formula is just an extrapolation and predicts too many emissions
!c
!c   1.  compute dF0/dln(rdry) = (dF0/drwet)*(drwet/dlnrdry) 
!c      at rdry_star
!c   2.  for rdry < rdry_star, assume dF0/dln(rdry) is lognormal,
!c      with the same lognormal parameters observed in 
!c      O'Dowd et al. [1997]
!c
2000    if (rdrylowermost .ge. rdry_star) goto 3000

!c   compute dF0/dln(rdry) at rdry_star
        rdryaa = 0.99*rdry_star
        rdry_cm = rdryaa*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetaa = rwet_cm*1.0e4

        rdrybb = 1.01*rdry_star
        rdry_cm = rdrybb*1.0e-4
        rwet_cm = ( rdry_cm**3 + (c1*(rdry_cm**c2))/            &
                ( (c3*(rdry_cm**c4)) - log10(relhum) ) )**onethird
        rwetbb = rwet_cm*1.0e4

        rwet = 0.5*(rwetaa + rwetbb)
        dumb = ( 0.380 - log10(rwet) ) / 0.650
        dumexpb = exp( -dumb*dumb)
        df0drwet = 1.373 * (rwet**(-3.0)) *                     &
                (1.0 + 0.057*(rwet**1.05)) *                    &
                (10.0**(1.19*dumexpb))

        drwet = rwetbb - rwetaa
        dlnrdry = log( rdrybb/rdryaa )
        df0dlnrdry_star = df0drwet * (drwet/dlnrdry)


!c   rdrylo, rdryhi = lower and upper dry radii (micrometers) 
!c   for this part of the integration
        rdrylo = rdrylowermost
        rdryhi = min( rdryuppermost, rdry_star )

        nsub_bin = 1000

        alnrdrylo = log( rdrylo )
        dlnrdry = (log( rdryhi ) - alnrdrylo)/nsub_bin

        do 2900 isub_bin = 1, nsub_bin

!c   geometric mean rdry (micrometers) for sub_bin
        dum = alnrdrylo + (isub_bin-0.5)*dlnrdry
        rdry = exp( dum )

!c   xmdry is dry mass in g
        xmdry = drydens_43pi_em12 * (rdry**3.0)

!c   xsdry is dry surface in cm2
        xsdry = x_4pi_em8 * (rdry**2.0)

!c   dumadjust is adjustment factor to reduce dF0/dr
        dum = log( rdry/rdry_star ) / log( sigmag_star )
        dumadjust = exp( -0.5*dum*dum )

        df0dlnrdry = df0dlnrdry_star * dumadjust

        dumsum_na = dumsum_na + dlnrdry*df0dlnrdry
        dumsum_ma = dumsum_ma + dlnrdry*df0dlnrdry*xmdry
        dumsum_sa = dumsum_sa + dlnrdry*df0dlnrdry*xsdry

2900    continue


!c
!c  all done
!c
3000    emitfact_numb = dumsum_na
        emitfact_mass = dumsum_ma
        emitfact_surf = dumsum_sa

        return
        end subroutine seasalt_emitfactors_1bin

!=================================================================================================================
 end module module_cam_mam_addemiss_seasalt
!=================================================================================================================
