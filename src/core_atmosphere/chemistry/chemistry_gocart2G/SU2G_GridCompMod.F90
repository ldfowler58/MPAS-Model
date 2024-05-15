!=================================================================================================================
 module SU2G_GridCompMod
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: Chem_Settling,SulfateChemDriver

 use SU2G_instance_SU,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                            fnum,rhFlag,pressure_lid_in_hPa,sigma,SO4_anthropogenic_fraction,      &
                            aircraft_fuel_emission_factor,aviation_vertical_layers
 use SU2G_StateSpecs,only: SU2G_State


 implicit none
 private


!--- constants (these parameters needs to be accessed from MPAS phys instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: airmw    = 28.97_RKIND
 real(kind=RKIND),parameter:: radTodeg = 180._RKIND/pi
 real(kind=RKIND),parameter:: Avogadro = 6.022e23


!--- relative position of sulfate tracers
 integer,parameter:: nDMS = 1, &
                     nSO2 = 2, &
                     nSO4 = 3, &
                     nMSA = 4


!--- molecular weights of sulfate species (grams):
 real(kind=RKIND),parameter:: fMassSulfur = 32._RKIND, &
                              fMassSO2    = 64._RKIND, &
                              fMassSO4    = 96._RKIND, &
                              fMassDMS    = 62._RKIND, &
                              fMassMSA    = 96._RKIND


!--- land/ocean/sea-ice mask (these parameters needs to be accessed from MPAS phys instead of redefined here):
 real(kind=RKIND),parameter:: OCEAN   = 2._RKIND, &
                              LAND    = 1._RKIND, &
                              SEA_ICE = 1._RKIND


!--- types needed to define SU2G:
 type:: ThreadWorkspace
    logical:: firstRun = .true.
    logical:: recycle_H2O2 = .false.

    integer:: nVolc = 0
    integer:: nPts = -1
    integer:: nymd_oxidants = -1 ! update the oxidant files?
    integer:: nymd_last = -1     ! previous nymd. updated daily.
    integer,dimension(:),allocatable:: pstart,pend
    integer,dimension(:),allocatable:: vStart,vEnd

    real(kind=RKIND),dimension(:),allocatable:: vLat,vLon,vSO2,vElev,vCloud
    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: SU2G_GridComp
!   logical:: diurnal_bb                              ! diurnal biomass burning
!   integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

    real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
    real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
    real(kind=RKIND):: fSO4anth                       ! fraction of anthropogenic emissions that are SO4
    real(kind=RKIND),dimension(:),allocatable:: sigma ! sigma of lognormal number distribution

!   !special handling for volcanic emissions:
!   character(len=strKIND):: volcano_srcfilen

!   !workspace for point emissions:
!   character(len=StrKIND):: point_emissions_srcfilen ! filename for pointwise emissions
!   logical:: doing_point_emissions = .false.
!   type(ThreadWorkspace),dimension(:),allocatable:: workspaces

    contains
       procedure:: load_GridComp => load_SU2G_GridComp
       procedure:: run2_GridComp => run2_SU2G_GridComp
 end type SU2G_GridComp

 type wrap_
    type(SU2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_SU2G_GridComp(self,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(SU2G_GridComp),intent(inout) :: self

!--- local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_SU2G_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization using parameters defined in SU2G_instance_SU:
 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 if(.not.allocated(self%sigma)) allocate(self%sigma(self%nbins))
 self%eAirCraftFuel = aircraft_fuel_emission_factor
 self%fSO4anth      = SO4_anthropogenic_fraction
 do n = 1,self%nbins
    self%sigma(n) = sigma(n)
 enddo
 do n = 1,4
    self%aviation_layers(n) = aviation_vertical_layers(n)
 enddo

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n),self%sigma(n)/))
 enddo


 call mpas_log_write('--- end subroutine load_SU2G_GridCOMP:')
!call mpas_log_write(' ')

 end subroutine load_SU2G_GridComp

!=================================================================================================================
 subroutine run2_SU2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(SU2G_GridComp),intent(inout):: self_params
 class(SU2G_State),intent(inout):: self

!--- local variables:
 logical:: correctionMaring
 integer:: istat
 integer:: ibin,i,j,k,n
 integer:: nymd,nhms

 real(kind=RKIND),dimension(:,:),allocatable:: drydepositionf
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qsu2G

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine run2_SU2G_GridComp:')


!--- initialization of klid:
!istat = 0
!call findKlid(self_params%klid,self_params%plid,self%ple,istat)


!--- SU2G oxidants:
!call SulfateUpdateOxidants( &
!--->         nymd_current   = nymd                          , &
!--->         nhms_current   = nhms                          , &
!--->         lonRad         = self%lons                     , &
!--->         latRad         = self%lats                     , &
!             rhoa           = self%airdens                  , &
!             km             = self%km                       , &
!             cdt            = self%cdt                      , &
!--->         nymd_last      = workspace%nymd_oxidants       , &
!--->         undefval       = MAPL_UNDEF                    , &
!--->         radToDeg       = real(MAPL_RADIANS_TO_DEGREES) , &
!--->         nAvogadro      = MAPL_AVOGAD/1000.             , &
!--->         pi             = pi                            , &
!--->         airMolWght     = airmw                         , &
!--->         oh_clim        = SU_OH                         , & !climatological OH source.
!--->         no3_clim       = SU_NO3                        , & !climatological NO3 source.
!--->         h2o2_clim      = SU_H2O2                       , & !climatological H2O2 source.
!--->         xoh            = xoh                           , &
!--->         xno3           = xno3                          , &
!--->         xh2o2          = xh2o2                         , &
!--->         recycle_h2o2   =  workspace%recycle_h2o2       , &
!--->         rc             = istat
!                         )
!subroutine SulfateUpdateOxidants( &
!             nymd_current, nhms_current, lonRad, latRad, &
!             rhoa, km, cdt, nymd_last, &
!             undefval, radToDeg, nAvogadro, pi, airMolWght, &
!             oh_clim, no3_clim, h2o2_clim, &
!             xoh, xno3, xh2o2, recycle_h2o2, rc)


!--- SU2G settling:
 if(.not.allocated(qsu2G)) allocate(qsu2G(its:ite,jts:jte,kts:kte,4))
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          qsu2G(i,j,k,1) = self%DMS(i,j,k)
          qsu2G(i,j,k,2) = self%SO2(i,j,k)
          qsu2G(i,j,k,3) = self%SO4(i,j,k)
          qsu2G(i,j,k,4) = self%MSA(i,j,k)
       enddo
    enddo
 enddo

 do ibin = 1, self_params%nbins
    !if radius == 0, then we're dealing with a gas which has no settling losses:
    if(self_params%radius(ibin) == 0.0) then
       if(associated(self%SUSD)) self%SUSD(:,:,ibin) = 0.0
       cycle
    endif

    call mpas_log_write('--- km        = $i',intArgs=(/self_params%km/))
    call mpas_log_write('--- klid      = $i',intArgs=(/self_params%klid/))
    call mpas_log_write('--- ibin      = $i',intArgs=(/ibin/))
    call mpas_log_write('--- rhFlag    = $i',intArgs=(/self_params%rhFlag/))
    call mpas_log_write('--- cdt       = $r',realArgs=(/self_params%cdt/))
    call mpas_log_write('--- grav      = $r',realArgs=(/grav/))
    call mpas_log_write('--- radiusInp = $r',realArgs=(/self_params%radius(ibin)/))
    call mpas_log_write('--- rhoInp    = $r',realArgs=(/self_params%rhop(ibin)/))
    call mpas_log_write(' ')

    call Chem_Settling( &
              km        = self_params%km                 , &
              klid      = self_params%klid               , &
              bin       = ibin                           , &
              flag      = self_params%rhFlag             , &
              cdt       = self_params%cdt                , &
              grav      = grav                           , &
              radiusInp = self_params%radius(ibin)*1.e-6 , &
              rhopInp   = self_params%rhop(ibin)         , &
              int_qa    = qsu2G(:,:,:,ibin)              , &
              tmpu      = self%t                         , &
              rhoa      = self%airdens                   , &
              rh        = self%rh2                       , &
              hghte     = self%zle                       , &
              ple       = self%ple                       , &
              delz      = self%delz                      , &
              delp      = self%delp                      , &
              fluxout   = self%SUSD                      , &
              rc        = istat                            &
                      )
 enddo
 if(allocated(qsu2G)) deallocate(qsu2G)


!--- SU2G chem driver:
 nymd = 0
 nhms = 0
 if(.not.allocated(drydepositionf)) allocate(drydepositionf(its:ite,jts:jte))

!call mpas_log_write('--- km         = $i',intArgs=(/self_params%km/))
!call mpas_log_write('--- klid       = $i',intArgs=(/self_params%klid/))
!call mpas_log_write('--- cdt        = $r',realArgs=(/self_params%cdt/))
!call mpas_log_write('--- pi         = $r',realArgs=(/pi/))
!call mpas_log_write('--- radToDeg   = $r',realArgs=(/radToDeg/))
!call mpas_log_write('--- von_karman = $r',realArgs=(/karman/))
!call mpas_log_write('--- airMolWght = $r',realArgs=(/airmw/))
!call mpas_log_write('--- nAvogadro  = $r',realArgs=(/Avogadro/1000./))
!call mpas_log_write('--- cpd        = $r',realArgs=(/cpd/))
!call mpas_log_write('--- grav       = $r',realArgs=(/grav/))
!call mpas_log_write('--- fMassMSA   = $r',realArgs=(/fMassMSA/))
!call mpas_log_write('--- fMassDMS   = $r',realArgs=(/fMassDMS/))
!call mpas_log_write('--- fMassSO2   = $r',realArgs=(/fMassSO2/))
!call mpas_log_write('--- fMassSO4   = $r',realArgs=(/fMassSO4/))
!call mpas_log_write('--- nymd       = $i',intArgs=(/nymd/))
!call mpas_log_write('--- nhms       = $i',intArgs=(/nhms/))
!call mpas_log_write('--- nDMS       = $i',intArgs=(/nDMS/))
!call mpas_log_write('--- nSO2       = $i',intArgs=(/nSO2/))
!call mpas_log_write('--- nSO4       = $i',intArgs=(/nSO4/))
!call mpas_log_write('--- nMSA       = $i',intArgs=(/nMSA/))
!call SulfateChemDriver( &
!             km             = self_params%km          , &
!             klid           = self_params%klid        , &
!             cdt            = self_params%cdt         , &
!             pi             = pi                      , &
!             radToDeg       = radToDeg                , &
!             von_karman     = karman                  , &
!             airMolWght     = airmw                   , &
!             nAvogadro      = Avogadro/1000.          , &
!             cpd            = cpd                     , &
!             grav           = grav                    , &
!             fMassMSA       = fMassMSA                , &
!             fMassDMS       = fMassDMS                , &
!             fMassSO2       = fMassSO2                , &
!             fMassSO4       = fMassSO4                , &
!             nymd           = nymd                    , &
!             nhms           = nhms                    , &
!             lonRad         = self%lons               , &
!             latRad         = self%lats               , &
!             dms            = self%DMS                , &
!             so2            = self%SO2                , &
!             so4            = self%SO4                , &
!             msa            = self%MSA                , &
!             nDMS           = nDMS                    , &
!             nSO2           = nSO2                    , &
!             nSO4           = nSO4                    , &
!             nMSA           = nMSA                    , &
!             xoh            = self%SU_OH              , &
!             xno3           = self%SU_NO3             , &
!             xh2o2          = self%SU_H2O2            , &
!             h2o2_init      = h2o2_init               , &
!             delp           = self%delp               , &
!             tmpu           = self%t                  , &
!             cloud          = self%fcld               , &
!             rhoa           = self%airdens            , &
!             hghte          = self%zle                , &
!             ustar          = self%ustar              , &
!             shflux         = self%sh                 , &
!             oro            = self%lwi                , &
!             pblh           = self%zpbl               , &
!             z0h            = self%z0h                , &
!             SU_dep         = self%SUDP               , &
!             SU_PSO2        = self%SUPSO2             , &
!             SU_PMSA        = self%SUPMSA             , &
!             SU_PSO4        = self%SUPSO4             , &
!             SU_PSO4g       = self%SUPSO4g            , &
!             SU_PSO4aq      = self%SUPSO4aq           , &
!             pso2           = self%PSO2               , &
!             pmsa           = self%PMSA               , &
!             pso4           = self%PSO4               , &
!             pso4g          = self%PSO4g              , &
!             pso4aq         = self%PSO4aq             , &
!             drydepositionfrequency = drydepositionf  , & ! 3d diagnostics
!             rc             = istat                     &
!                      )
 if(allocated(drydepositionf)) deallocate(drydepositionf)


!--- SU2G wet removal:
!kin = .true.
!call SU_Wet_Removal (self%km, self%nbins, self%klid, self%cdt, kin, MAPL_GRAV, MAPL_AIRMW, &
!                     delp, fMassSO4, fMassSO2, &
!                     h2o2_init, ple, airdens, cn_prcp, ncn_prcp, pfl_lsan, pfi_lsan, t, &
!                     nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, dummyMSA, &
!                     SUWT, SUPSO4, SUPSO4WT, PSO4, PSO4WET, __RC__ )

!subroutine SU_Wet_Removal (km, nbins, klid, cdt, kin, grav, airMolWght, delp, fMassSO4, fMassSO2, &
!                           h2o2_int, ple, rhoa, precc, precl, pfllsan, pfilsan, tmpu, &
!                           nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, MSA, &
!                           fluxout, pSO4_colflux, pSO4wet_colflux, &
!                           pso4, pso4wet, rc )


!--- SU2G diagnostics:
!Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
!call SU_Compute_Diags( &
!             km         = self_params%km                           , &
!             klid       = self_params%klid                         , &
!             rmed       = self%radius(nSO4)                        , &
!             sigma      = self%sigma(nSO4)                         , & 
!             rhop       = self%rhop(nSO4)                          , &
!             grav       = grav                                     , &
!             pi         = pi                                       , &
!             nSO4       = nSO4                                     , &
!             mie        = self%diag_Mie                            , &
!             wavelengths_profile = self%wavelengths_profile*1.0e-9 , &
!             wavelengths_vertint = self%wavelengths_vertint*1.0e-9 , &
!             tmpu       = self%t                                   , &
!             rhoa       = self%airdens                             , &
!             delp       = self%delp                                , &
!             ple        = self%ple                                 , &
!             tropp      = tropp                                    , &
!             rh2        = self%rh2                                 , &
!             u          = self%u                                   , &
!             v          = self%v                                   , &
!--->         DMS        = DMS                                      , &
!--->         SO2        = SO2                                      , &
!--->         SO4        = SO4                                      , &
!--->         dummyMSA   = dummyMSA                                 , &
!             dmssfcmass = self%DMSSMASS                            , &
!             dmscolmass = self%DMSCMASS                            , &
!             msasfcmass = self%MSASMASS                            , &
!             msacolmass = self%MSACMASS                            , &
!             so2sfcmass = self%SO2SMASS                            , &
!             so2colmass = self%SO2CMASS                            , &
!             so4sfcmass = self%SO4SMASS                            , &
!             so4colmass = self%SO4CMASS                            , &
!             exttau     = self%SUEXTTAU                            , &
!             stexttau   = self%SUSTEXTTAU                          , &
!             scatau     = self%SUSCATAU                            , &
!             stscatau   = self%SUSTSCATAU                          , &
!             so4mass    = self%SO4MASS                             , &
!             so4conc    = self%SUCONC                              , &
!             extcoef    = self%SUEXTCOEF                           , &
!             scacoef    = self%SUSCACOEF                           , &
!             bckcoef    = self%SUBCKCOEF                           , &
!             angstrom   = self%SUANGSTR                            , &
!             fluxu      = self%SUFLUXU                             , &
!             fluxv      = self%SUFLUXV                             , &
!             sarea      = self%SO4SAREA                            , &
!             snum       = self%SO4SNUM                             , &
!             rc         = istat
!                     )
!subroutine SU_Compute_Diags (km, klid, rmed, sigma, rhop, grav, pi, nSO4, mie, &
!                             wavelengths_profile, wavelengths_vertint, &
!                             tmpu, rhoa, delp, ple, tropp,rh, u, v, &
!                             DMS, SO2, SO4, MSA, &
!                             dmssfcmass, dmscolmass, &
!                             msasfcmass, msacolmass, &
!                             so2sfcmass, so2colmass, &
!                             so4sfcmass, so4colmass, &
!                             exttau, stexttau,scatau, stscatau,so4mass, so4conc, extcoef, &
!                             scacoef, bckcoef, angstrom, fluxu, fluxv, sarea, snum, rc )

!i1 = lbound(RH2, 1); i2 = ubound(RH2, 1) 
!j1 = lbound(RH2, 2); j2 = ubound(RH2, 2) 
!km = ubound(RH2, 3)                      
     
!allocate(RH20(i1:i2,j1:j2,km), __STAT__) 
!allocate(RH80(i1:i2,j1:j2,km), __STAT__)

!RH20(:,:,:) = 0.20
!call SU_Compute_Diags (km=self%km, klid=self%klid, rmed=self%radius(nSO4), sigma=self%sigma(nSO4),&
!                       rhop=self%rhop(nSO4), &
!                       grav=MAPL_GRAV, pi=MAPL_PI, nSO4=nSO4, mie=self%diag_Mie, &
!                       wavelengths_profile=self%wavelengths_profile*1.0e-9, &
!                       wavelengths_vertint=self%wavelengths_vertint*1.0e-9, &
!                       tmpu=t, rhoa=airdens, delp=delp, ple=ple,tropp=tropp, rh=rh20, u=u, v=v, &
!                       DMS=DMS, SO2=SO2, SO4=SO4, MSA=dummyMSA,extcoef=SUEXTCOEFRH20, &
!                       scacoef = SUSCACOEFRH20, __RC__)

!RH80(:,:,:) = 0.80
!call SU_Compute_Diags (km=self%km, klid=self%klid, rmed=self%radius(nSO4), sigma=self%sigma(nSO4),&
!                       rhop=self%rhop(nSO4), &
!                       grav=MAPL_GRAV, pi=MAPL_PI, nSO4=nSO4, mie=self%diag_Mie, &
!                       wavelengths_profile=self%wavelengths_profile*1.0e-9, &
!                       wavelengths_vertint=self%wavelengths_vertint*1.0e-9, &
!                       tmpu=t, rhoa=airdens, delp=delp, ple=ple,tropp=tropp, rh=rh80, u=u, v=v, &
!                       DMS=DMS, SO2=SO2, SO4=SO4, MSA=dummyMSA,extcoef=SUEXTCOEFRH80,&
!                       scacoef = SUSCACOEFRH80, __RC__)


 call mpas_log_write('--- end subroutine run2_SU2G_GridComp:')
 call mpas_log_write(' ')

 end subroutine run2_SU2G_GridComp

!=================================================================================================================
 end module SU2G_GridCompMod
!=================================================================================================================

