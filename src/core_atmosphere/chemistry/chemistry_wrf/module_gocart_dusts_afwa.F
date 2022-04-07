!=================================================================================================================
 module module_gocart_dusts_data_afwa
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND

 save

 integer,parameter :: ndust=5,ndcls=3,ndsrc=1,maxstypes=100
 integer,parameter :: ngsalt=9

!GAC-->
!20130219 - Drypoint no longer needed.  Use NOAA porosity for all schemes to
!allow use of AFWA dust scheme by all LSMs, not just NOAA and RUC.
!real(kind=RKIND),dimension(maxstypes):: drypoint
!real(kind=RKIND),dimension(maxstypes):: porosity

 real(kind=RKIND),dimension(19),parameter:: &
    porosity=(/0.339, 0.421, 0.434, 0.476, 0.476, 0.439, &
               0.404, 0.464, 0.465, 0.406, 0.468, 0.468, &
               0.439, 1.000, 0.200, 0.421, 0.468, 0.200, &
               0.339/)
!<--GAC

 integer,dimension(5),parameter:: ipoint(5)=(/3,2,2,2,2/)

 real(kind=RKIND),parameter:: dyn_visc = 1.5E-5
 real(kind=RKIND),dimension(5),parameter:: frac_s(5)=(/0.1,0.25,0.25,0.25,0.25/)

 real(kind=RKIND),dimension(ndust,12):: ch_dust

!real(kind=R8KIND),dimension(5),parameter:: den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
 real(kind=R8KIND),dimension(5)          :: den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
 real(kind=R8KIND),dimension(5),parameter:: reff_dust(5)=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/)

!saltation:
 integer,dimension(ngsalt),parameter:: spoint=(/1,2,2,2,2,2,3,3,3/)  ! 1 Clay, 2 Silt, 3 Sand

 real(kind=R8KIND),dimension (ngsalt),parameter:: &
      reff_salt=(/0.71D-6,1.37D-6,2.63D-6,5.00D-6,9.50D-6,18.1D-6,34.5D-6,65.5D-6,125.D-6/)
 real(kind=R8KIND),dimension (ngsalt),parameter:: den_salt=(/2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650./)
 real(kind=R8KIND),dimension (ngsalt),parameter:: frac_salt=(/1.,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333/)


!dust:
 real(kind=R8KIND),dimension(ndust),parameter:: distr_dust=(/1.074d-1,1.012d-1,2.078d-1,4.817d-1,1.019d-1/)
!real(kind=R8KIND),dimension(ndust),parameter:: lo_dust=(/0.1d-6,1.0d-6,1.8d-6,3.0d-6,6.0d-6/)
!real(kind=R8KIND),dimension(ndust),parameter:: up_dust=(/1.0d-6,1.8d-6,3.0d-6,6.0d-6,10.0d-6/)


!GOCART bins radii ranges
 real(kind=R8KIND),dimension(5),parameter:: ra_dust(5)=(/1.d-1,1.0d0,1.8d0,3.0d0,6.0d0/)
 real(kind=R8KIND),dimension(5),parameter:: rb_dust(5)=(/1.0d0,1.8d0,3.0d0,6.0d0,10.0d0/)

!=================================================================================================================
 end module module_gocart_dusts_data_afwa
!=================================================================================================================


!=================================================================================================================
 module module_gocart_dusts_afwa
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND
 use mpas_log
 use module_gocart_dusts_data_afwa


 implicit none
 private
 public:: gocart_dusts_emission_afwa


 contains


!=================================================================================================================
 subroutine gocart_dusts_emission_afwa(its,ite,kts,kte,ktep1,nerod,dt,grav,qdust1,qdust2,qdust3,qdust4,qdust5, &
                                       isltyp,xland,area,erod,clay,sand,zgrid,u10,v10,smois,snowh,ust,znt,rho, &
                                       emis_qdust1,emis_qdust2,emis_qdust3,emis_qdust4,emis_qdust5)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte,ktep1
 integer,intent(in):: nerod
 integer,intent(in),dimension(its:ite):: isltyp

 real(kind=RKIND),intent(in):: dt,grav
 real(kind=RKIND),intent(in),dimension(its:ite):: area,xland
 real(kind=RKIND),intent(in),dimension(its:ite):: u10,v10,smois,snowh,ust,znt
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho
 real(kind=RKIND),intent(in),dimension(its:ite,kts:ktep1):: zgrid

 real(kind=RKIND),intent(in),dimension(its:ite):: clay,sand
 real(kind=RKIND),intent(in),dimension(its:ite,nerod):: erod


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qdust1,qdust2,qdust3,qdust4,qdust5
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte),optional:: emis_qdust1,emis_qdust2, &
                                                                      emis_qdust3,emis_qdust4, &
                                                                      emis_qdust5


!local variables and arrays:
 integer,parameter:: &
    imx = 1,     &
    lmx = 1,     &
    nmx = ndust, &
    smx = ngsalt
 integer:: i,k,n
 integer:: smois_opt
 integer,dimension(imx):: ilwi

 real(kind=RKIND),parameter:: &
    alpha  = 1.0,& ! tuning parameter.
    gamma  = 1.0   ! tuning parameter.

 real(kind=R8KIND),parameter:: &
    smtune = 1.0,& !
    ustune = 1.0   ! tuning parameter.

 real(kind=R8KIND):: conver,converi
 real(kind=R8KIND):: ustart
 real(kind=R8KIND),dimension(3):: massfrac
 real(kind=R8KIND),dimension(imx):: airden,drylimit,dz,erodtot,gravsm,vegmask,volsm,ustar
 real(kind=R8KIND),dimension(imx,nmx):: bems
 real(kind=R8KIND),dimension(imx,lmx,nmx):: tc

 real(kind=RKIND):: tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_dusts_emission_afwa:')
 call mpas_log_write('--- nmx = $i',intArgs=(/nmx/))
 call mpas_log_write('--- smx = $i',intArgs=(/smx/))

 conver  = 1.d-9
 converi = 1.d9

!initialization:
 smois_opt = 0

 k = kts
 do i = its,ite
    ustart  = 0._R8KIND

    if(xland(i).lt.1.5)then
       ilwi(imx)  = 1
       airden(imx) = rho(i,kts)
       ustar(imx)  = ustune*ust(i)
       dz(imx) = zgrid(i,k+1)-zgrid(i,k)
       dz(imx) = max(0.,dz(imx))

       tc(imx,lmx,1) = qdust1(i,kts)
       tc(imx,lmx,2) = qdust2(i,kts)
       tc(imx,lmx,3) = qdust3(i,kts)
       tc(imx,lmx,4) = qdust4(i,kts)
       tc(imx,lmx,5) = qdust5(i,kts)

       !default value for total erodiblity (Ginoux):
       erodtot(:) = 0._R8KIND
       do n = 1,nerod
          erodtot(imx) = erodtot(imx) + erod(i,n)
       enddo

       !default value for the vegetation mask (Ginoux):
       if(erod(i,1) .eq. 0.) then
          vegmask(imx) = 0.0
       else
          vegmask(imx) = 1.0
       endif
       !remove vegetated areas (vegmask=0) from total erodibility:
       erodtot(imx) = vegmask(imx)*erodtot(imx)

       !mass fractions of clay, silt, and sand:
       massfrac(1) = clay(i)
       massfrac(2) = 1-(clay(i)+sand(i))
       massfrac(3) = sand(i)

       !don't allow roughness lengths greater than 20 cm to be lofted.
       !this kludge accounts for land use types like urban areas and
       !forests which would otherwise show up as high dust emitters.
       !this is a placeholder for a more widely accepted kludge
       !factor in the literature, which reduces lofting for rough areas.
       !forthcoming...
       if(znt(i) .gt. 0.2) then
          ilwi(imx)=0
       endif

       !do not allow areas with bedrock, lava, or land-ice to loft:
       if(isltyp(i).eq.15 .or. isltyp(i).eq.16. .or. isltyp(i).eq.18) then
          ilwi(imx) = 0
       endif

       !another hack to ensure dust does not loft from areas with snow
       !cover...because, well, that doesn't make sense:
       if(snowh(i) .gt. 0.01) then
          ilwi(imx) = 0
       endif

       !volumetric soil moisture can be tuned here with a dust_smtune
       !set in the namelist:
       volsm(imx) = max(smois(i)*smtune,0.)

       !calculate gravimetric soil moisture.
       gravsm(imx) = 100*volsm(imx)/((1.-porosity(isltyp(i)))*(2.65*(1-clay(i))+2.50*clay(i)))

       !use drylimit based on gravimetric soil moisture:
       drylimit(imx) = 14.0*clay(i)*clay(i)+17.0*clay(i)

       smois_opt = 0

       !call dust emission:
       call source_dust_afwa(i,imx,lmx,nmx,smx,dt,tc,ustar,massfrac,erodtot,ilwi,gravsm,volsm, &
                             airden,dz,bems,ustart,grav,drylimit,alpha,gamma,smois_opt)

       qdust1(i,kts) = real(tc(imx,lmx,1),RKIND)
       qdust2(i,kts) = real(tc(imx,lmx,2),RKIND)
       qdust3(i,kts) = real(tc(imx,lmx,3),RKIND)
       qdust4(i,kts) = real(tc(imx,lmx,4),RKIND)
       qdust5(i,kts) = real(tc(imx,lmx,5),RKIND)

       tem1 = real(bems(imx,1),RKIND)
       tem2 = real(bems(imx,2),RKIND)
       tem3 = real(bems(imx,3),RKIND)
       tem4 = real(bems(imx,4),RKIND)
       tem5 = real(bems(imx,5),RKIND)
    endif

 enddo

 call mpas_log_write('--- end subroutine gocart_dusts_emission_afwa:')

 end subroutine gocart_dusts_emission_afwa

!=================================================================================================================
 subroutine source_dust_afwa(ii,imx,lmx,nmx,smx,dt1,tc,ustar,massfrac,erod,ilwi,gravsm,volsm,airden, &
                             dz_lowest,bems,ustart,g0,drylimit,alpha,gamma,smois_opt)

 ! ****************************************************************************
 ! *  Evaluate the source of each dust particles size bin by soil emission  
 ! *
 ! *  Input:
 ! *         EROD      Fraction of erodible grid cell                (-)
 ! *         ILWI      Land/water flag                               (-)
 ! *         GRAVSM    Gravimetric soil moisture                     (g/g)
 ! *         VOLSM     Volumetric soil moisture                      (g/g)
 ! *         SOILM_OPT Soil moisture option (1:Use GRAVSM 2:VOLSM)   (-)
 ! *         DRYLIMIT  Upper GRAVSM (VOLSM) limit for air-dry soil   (g/g)
 ! *         ALPHA     Constant to fudge the total emission of dust  (1/m)
 ! *         GAMMA     Exponential tuning constant for erodibility   (-)
 ! *         AIRDEN    Density of air for each grid box              (kg/m3)
 ! *         USTAR     Friction velocity                             (m/s)
 ! *         MASSFRAC  Fraction of mass in each of 3 soil classes    (-)
 ! *         DT1       Time step                                     (s)
 ! *         NMX       Number of dust bins                           (-)
 ! *         SMX       Number of saltation bins                      (-)
 ! *         IMX       Number of I points                            (-)
 ! *         JMX       Number of J points                            (-)
 ! *         LMX       Number of L points                            (-)
 ! *         dz_lowest heigth of the lowest layer                    (m)
 ! *
 ! *  Data (see module_data_gocart_dust):
 ! *         SPOINT    Pointer to 3 soil classes                     (-)
 ! *         DEN_DUST  Dust density                                  (kg/m3)
 ! *         DEN_SALT  Saltation particle density                    (kg/m3)
 ! *         REFF_SALT Reference saltation particle diameter         (m)
 ! *         REFF_DUST Reference dust particle diameter              (m)
 ! *         LO_DUST   Lower diameter limits for dust bins           (m)
 ! *         UP_DUST   Upper diameter limits for dust bins           (m)
 ! *         FRAC_SALT Soil class mass fraction for saltation bins   (-)
 ! *
 ! *  Parameters:
 ! *         BETAMAX   Maximum sandblasting mass efficiency          (-)
 ! *         CMB       Constant of proportionality                   (-)
 ! *         MMD_DUST  Mass median diameter of dust                  (m)
 ! *         GSD_DUST  Geometric standard deviation of dust          (-)
 ! *         LAMBDA    Side crack propogation length                 (m)
 ! *         CV        Normalization constant                        (-)
 ! *         G0        Gravitational acceleration                    (m/s2)
 ! *         G         Gravitational acceleration in cgs             (cm/s2)
 ! *      
 ! *  Working:
 ! *         BETA      Sandblasting mass efficiency                  (-)
 ! *         U_TS0     "Dry" threshold friction velocity             (m/s)
 ! *         U_TS      Moisture-adjusted threshold friction velocity (m/s)
 ! *         RHOA      Density of air in cgs                         (g/cm3)
 ! *         DEN       Dust density in cgs                           (g/cm3)
 ! *         DIAM      Dust diameter in cgs                          (cm)
 ! *         DMASS     Saltation mass distribution                   (-)
 ! *         DSURFACE  Saltation surface area per unit mass          (m2/kg)
 ! *         DS_REL    Saltation surface area distribution           (-)
 ! *         SALT      Saltation flux                                (kg/m/s)
 ! *         DLNDP     Dust bin width                                (-)
 ! *         EMIT      Total vertical mass flux                      (kg/m2/s)
 ! *         EMIT_VOL  Total vertical volume flux                    (m/s)
 ! *         DSRC      Mass of emitted dust               (kg/timestep/m2)
 ! *
 ! *  Output:
 ! *         TC        Total concentration of dust                 (kg/kg)
 ! *         BEMS      Source of each dust type           (kg/timestep/m2)
 ! *         USTART    Threshold friction vel. (bin 7)               (m/s)
 ! *
 ! ****************************************************************************

!=================================================================================================================

!input arguments:
 integer,intent(in):: ii
 integer,intent(in):: nmx,imx,lmx,smx
 integer,intent(in):: smois_opt
 integer,intent(in),dimension(imx):: ilwi

 real(kind=RKIND),intent(in):: g0,dt1
 real(kind=R8KIND),intent(in),dimension(imx):: erod
 real(kind=R8KIND),intent(in),dimension(imx):: ustar
 real(kind=R8KIND),intent(in),dimension(imx):: gravsm
 real(kind=R8KIND),intent(in),dimension(imx):: drylimit
 real(kind=R8KIND),intent(in),dimension(imx):: dz_lowest
 real(kind=R8KIND),intent(in),dimension(imx):: volsm
 real(kind=R8KIND),intent(in),dimension(imx,lmx):: airden


!inout arguments:
 real(kind=R8KIND),intent(inout),dimension(imx,lmx,nmx):: tc


!output arguments:
 real(kind=R8KIND),intent(out):: ustart
 real(kind=R8KIND),intent(out),dimension(imx,nmx):: bems


!local variables and arrays:
 real(kind=R8KIND),dimension(smx):: den,diam
!real(kind=R8KIND),dimension(nmx):: dvol,dlndp
!real(kind=R8KIND),dimension(nmx):: distr_dust
 real(kind=R8KIND),dimension(smx):: dsurface,ds_rel
 real(kind=R8KIND),dimension(3)  :: massfrac
 real(kind=R8KIND):: u_ts0,u_ts,dsrc,dmass!, dvol_tot
 real(kind=R8KIND):: emit!, emit_vol
 real(kind=R8KIND):: salt, stotal
 real(kind=RKIND) :: rhoa,g
 integer   :: i, j, m, n, s

! Global tuning constant, alpha.  Sandblasting mass efficiency, beta.
! Beta maxes out for clay fractions above 0.2 => betamax.

  real(kind=RKIND),intent(in):: alpha
  real,parameter :: betamax=5.25e-4
  real(kind=R8KIND):: beta

! Experimental optional exponential tuning constant for erodibility.
! 0 < gamma < 1 -> more relative impact by low erodibility regions.
! 1 < gamma -> accentuates high erodibility regions.  Recommend this
! be set to 1 (default) unless looking for a way to increase spread
! within an ensemble framework.
  
  real(kind=RKIND),intent(in):: gamma

! Constant of proportionality from Marticorena et al, 1997 (unitless)
! Arguably more ~consistent~ fudge than alpha, which has many walnuts
! sprinkled throughout the literature. 

  real(kind=RKIND):: cmb=1.0    ! Marticorena et al,1997
! real(kind=RKIND):: cmb=2.61   ! White,1979

! Parameters used in Kok distribution function. Advise not to play with 
! these without the expressed written consent of someone who knows what
! they're doing. (See Kok, 2010 PNAS for details on this scheme).

! real(kind=R8KIND),parameter:: mmd_dust=3.4D-6  ! median mass diameter (m)
! real(kind=R8KIND),parameter:: gsd_dust=3.0     ! geometric standard deviation
! real(kind=R8KIND),parameter:: lambda=12.0D-6   ! crack propogation length (m)
! real(kind=R8KIND),parameter:: cv=12.62D-6      ! normalization constant


!-----------------------------------------------------------------------------------------------------------------

! Calculate saltation surface area distribution from sand, silt, and clay
! mass fractions and saltation bin fraction. This will later become a 
! modifier to the total saltation flux.  The reasoning here is that the 
! size and availability of saltators affects saltation efficiency. Based
! on Eqn. (32) in Marticorena & Bergametti, 1995 (hereon, MB95).

 do n = 1,smx
    dmass = massfrac(spoint(n))*frac_salt(n)
    dsurface(n) = 0.75*dmass/(den_salt(n)*reff_salt(n))  
 enddo
  
! The following equation yields relative surface area fraction.  It will only
! work if you are representing the "full range" of all three soil classes.
! For this reason alone, we have incorporated particle sizes that encompass
! the clay class, to account for its relative area over the basal
! surface, even though these smaller bins would be unlikely to play any large
! role in the actual saltation process.

 stotal = 0.
 do n = 1,smx
    stotal = stotal + dsurface(n)
 enddo
 do n = 1,smx
    ds_rel(n)=dsurface(n)/stotal
 enddo

! Calculate total dust emission due to saltation of sand sized particles.
! Begin by calculating DRY threshold friction velocity (u_ts0).  Next adjust
! u_ts0 for moisture to get threshold friction velocity (u_ts). Then
! calculate saltation flux (salt) where ustar has exceeded u_ts.  Finally, 
! calculate total dust emission (tot_emit), taking into account erodibility. 

 g = g0*1.0E2                         ! (cm s^-2)
 emit=0.0

 do n = 1,smx                         ! Loop over saltation bins
    den(n) = den_salt(n)*1.0D-3       ! (g cm^-3)
    diam(n) = 2.0*reff_salt(n)*1.0D2  ! (cm)

    do i = 1,imx
       rhoa = airden(i,1)*1.0D-3     ! (g cm^-3)
 
       ! Threshold friction velocity as a function of the dust density and
       ! diameter from Bagnold (1941) (m s^-1).

       u_ts0 = 0.13*1.0D-2*sqrt(den(n)*g*diam(n)/rhoa)* &
               sqrt(1.0+0.006/den(n)/g/(diam(n))**2.5)/ &
               sqrt(1.928*(1331.0*(diam(n))**1.56+0.38)**0.092-1.0) 

       ! Friction velocity threshold correction function based on physical
       ! properties related to moisture tension. Soil moisture greater than
       ! dry limit serves to increase threshold friction velocity (making
       ! it more difficult to loft dust). When soil moisture has not reached
       ! dry limit, treat as dry (no correction to threshold friction
       ! velocity). GAC

       ! Calculate threshold friction velocity. If volumetric (gravimetric)
       ! water content is less than the drylimit, then the threshold friction
       ! velocity (u_ts) will be equal to the dry threshold friction velocity
       ! (u_ts0). EDH

       if(smois_opt .eq. 1) then
          if(100.*volsm(i) > drylimit(i)) then
            u_ts = max(0.0d+0,u_ts0*(sqrt(1.0+1.21*((100.*volsm(i))-drylimit(i))**0.68)))
          else
            u_ts = u_ts0
          endif
       else
          if(gravsm(i) > drylimit(i)) then
            u_ts = MAX(0.0D+0,u_ts0*(sqrt(1.0+1.21*(gravsm(i)-drylimit(i))**0.68)))
          else
            u_ts = u_ts0
          endif
       endif

       ! Bin 7 threshold friction velocity for diagnostic dust lofting
       ! potential product

       if(n .eq. 7) then ! Saltation bin 7 is small sand
          ustart = u_ts
       endif
!      call mpas_log_write('uts: $i $i $r $r $r',intArgs=(/ii,n/), &
!                realArgs=(/real(u_ts0,RKIND),real(u_ts,RKIND),real(ustart,RKIND)/))
!      if(n.eq.smx) call mpas_log_write(' ')
 
       ! Saltation flux (kg m^-1 s^-1) from MB95
       ! ds_rel is the relative surface area distribution
       if(ustar(i) .gt. u_ts .and. erod(i) .gt. 0.0 .and. ilwi(i) == 1) then
          salt = cmb*ds_rel(n)*(airden(i,1)/g0)*(ustar(i)**3)* &
                 (1. + u_ts/ustar(i))*(1. - (u_ts**2)/(ustar(i)**2))
!         call mpas_log_write('--- salt: $i $i $r',intArgs=(/n,ii/),realArgs=(/real(salt,RKIND)/))
       else
          salt = 0.D0
       endif
 
       ! Calculate total vertical mass flux (note beta has units of m^-1)
       ! Beta acts to tone down dust in areas with so few dust-sized particles
       ! that the lofting efficiency decreases.  Otherwise, super sandy zones
       ! would be huge dust producers, which is generally not the case.
       ! Equation derived from wind-tunnel experiments (see MB95).
 
       beta=10**(13.6*massfrac(1)-6.0)  ! (unitless)
       if(beta .gt. betamax) then
          beta = betamax
       endif
       emit=emit+salt*(erod(i)**gamma)*alpha*beta    ! (kg m^-2 s^-1)
!      if(emit .gt. 0.) &
!         call mpas_log_write('--- emit: $i $i $r',intArgs=(/n,ii/),realArgs=(/real(emit,RKIND)/))
    enddo
 enddo ! End do over saltation bins

! Now that we have the total dust emission, distribute into dust bins using 
! lognormal distribution (Dr. Jasper Kok, 2010), and
! calculate total mass emitted over the grid box over the timestep. 
!
! In calculating the Kok distribution, we assume upper and lower limits to
! each bin. For reff_dust=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/) (default),
! lower limits were ASSUMED at lo_dust=(/0.1D-6,1.0D-6,1.8D-6,3.0D-6,6.0D-6/)
! upper limits were ASSUMED at up_dust=(/1.0D-6,1.8D-6,3.0D-6,6.0D-6,10.0D-6/)
! These may be changed within module_data_gocart_dust.F, but make sure it is
! consistent with reff_dust values.  These values were taken from the original
! GOCART bin configuration. We use them here to calculate dust bin width,
! dlndp. dvol is the volume distribution. GAC
!
! UPDATE: We bypass the calculation below and instead hardcode distr_dust for
! the five dust bins we are using here since this distribution is static and
! unnecessary to calculate at every time step. Keeping everything here to
! document the steps in obtaining distr_dust.  GAC 20140320
  
!dvol_tot=0.
!do n = 1,nmx  ! loop over all dust bins
!   dlndp(n) = log(up_dust(n)/lo_dust(n))
!   dvol(n)  = (2.0*reff_dust(n)/cv)*(1.+erf(log(2.0*reff_dust(n)/mmd_dust)/(sqrt(2.)*log(gsd_dust))))*&
!               exp(-(2.0*reff_dust(n)/lambda)**3.0)*dlndp(n)
!   dvol_tot = dvol_tot+dvol(n)
!
!   !convert mass flux to volume flux
!   emit_vol=emit/den_dust(n) ! (m s^-1)
!enddo
!do n = 1,nmx  ! loop over all dust bins
!   distr_dust(n)=dvol(n)/dvol_tot
!enddo
 
! now distribute total vertical emission into dust bins and update
! concentration:

 do n = 1,nmx  ! loop over all dust bins
    do i = 1,imx
       !calculate total mass emitted:
       dsrc = emit*distr_dust(n)*dt1  ! (kg m^-2) per dt1
       if(dsrc < 0.0) dsrc = 0.0

       ! update dust mixing ratio at first model level.
       tc(i,1,n) = tc(i,1,n) + dsrc/dz_lowest(i)/airden(i,1) ! (kg/kg)
       bems(i,n) = dsrc                                   ! diagnostic (kg/m2) per dt1
    enddo
 enddo

 end subroutine source_dust_afwa

!=================================================================================================================
 end module module_gocart_dusts_afwa
!=================================================================================================================

