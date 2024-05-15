!=================================================================================================================
 module NI2G_StateSpecs
 use mpas_kind_types,only: RKIND

 implicit none
 public

!this module is the state variable specification file for nitrate parameters. it is the same as NI2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/NI2G_GridComp.

!schema_version: 2.0.0
!component: NI


 type NI2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer    :: LWI           => null() ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),pointer    :: TROPP         => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer    :: USTAR         => null() ! surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ZPBL          => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer    :: SH            => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer    :: Z0H           => null() !  surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer    :: CN_PRCP       => null() ! surface_conv._rain_flux_needed_by_land (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NCN_PRCP      => null() ! Non-convective precipitation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: AREA          => null() ! agrid_cell_area (^2)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer  :: AIRDENS       => null() ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: DELP          => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer  :: T             => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer  :: RH2           => null() ! Rel_Hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ZLE           => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer  :: PLE           => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer  :: PFL_LSAN      => null() ! 3D_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer  :: PFI_LSAN      => null() ! 3D_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer  :: U             => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: V             => null() ! northward_wind (m s-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_AG    => null() ! agriculture emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_BB    => null() ! biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_EN    => null() ! energy emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_IN    => null() ! industry emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_OC    => null() ! ocean emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_RE    => null() ! resedential emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: EMI_NH3_TR    => null() ! transport emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NITRATE_HNO3  => null() ! nitrate hno3 emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: DU            => null() ! Dust Mixing Ratio all bins (kg kg-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SS            => null() ! Sea Salt Mixing Ratio all bins (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: SO4           => null() ! Sulfate Mixing Ratio (kg kg-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: NH3MASS       => null() ! Ammonia Mass Mixing Ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NH4MASS       => null() ! Ammonium Aerosol Mass Mixing Ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIMASS        => null() ! Nitrate Mass Mixing Ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIMASS25      => null() ! Nitrate Mass Mixing Ratio [PM2.5] (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: HNO3CONC      => null() ! Nitric Acid Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NH3CONC       => null() ! Ammonia Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NH4CONC       => null() ! Ammonium Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NICONC        => null() ! Nitrate Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NICONC25      => null() ! Nitrate Mass Concentration [PM2.5] (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NIEXTCOEF     => null() ! Nitrate Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NIEXTCOEFRH20 => null() ! Nitrate Extinction Coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NIEXTCOEFRH80 => null() ! Nitrate Extinction Coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NISCACOEF     => null() ! Nitrate Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NISCACOEFRH20 => null() ! Nitrate Scattering Coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NISCACOEFRH80 => null() ! Nitrate Scattering Coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: NIBCKCOEF     => null() ! Nitrate Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer    :: NIPNO3AQ      => null() ! Nitrate Production from Aqueous Chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NIPNH4AQ      => null() ! Ammonium Production from Aqueous Chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NIPNH3AQ      => null() ! Ammonia Change from Aqueous Chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIHT          => null() ! Nitrate Production from Het Chem (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISD          => null() ! Nitrate Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIDP          => null() ! Nitrate Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIWT          => null() ! Nitrate Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISV          => null() ! Nitrate Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH3EM         => null() ! Ammonia Emission (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH3DP         => null() ! Ammonia Dry Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH3WT         => null() ! Ammonia Wet Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH3SV         => null() ! Ammonia Convective Scavenging (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH4SD         => null() ! Ammonium Settling (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH4DP         => null() ! Ammonium Dry Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH4WT         => null() ! Ammonium Wet Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NH4SV         => null() ! Ammonium Convective Scavenging (kg m-2 s-1)

 real(kind=RKIND),dimension(:,:),pointer    :: HNO3SMASS     => null() ! Nitric Acid Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NH3SMASS      => null() ! Ammonia Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NH4SMASS      => null() ! Ammonium Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NISMASS       => null() ! Nitrate Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NISMASS25     => null() ! Nitrate Surface Mass Concentration [PM2.5] (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: HNO3CMASS     => null() ! Nitric Acid Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NH3CMASS      => null() ! Ammonia Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NH4CMASS      => null() ! Ammonium Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: NICMASS       => null() ! Nitrate Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:),pointer    :: NICMASS25     => null() ! Nitrate Column Mass Density [PM2.5] (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIEXTTFM      => null() ! Nitrate Extinction AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISCATFM      => null() ! Nitrate Scattering AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIEXTT25      => null() ! Nitrate Extinction AOT - PM 2.5 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISCAT25      => null() ! Nitrate Scattering AOT - PM 2.5 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NIEXTTAU      => null() ! Nitrate Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISTEXTTAU    => null() ! Nitrate Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISCATAU      => null() ! Nitrate Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NISTSCATAU    => null() ! Nitrate Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer    :: NIANGSTR      => null() ! Nitrate Angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:),pointer    :: NIFLUXU       => null() ! Nitrate column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: NIFLUXV       => null() ! Nitrate column v-wind mass flux (kg m-1 s-1)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: NH3           => null() ! Ammonia (NH3, gas phase) (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NH4a          => null() ! Ammonium ion (NH4+, aerosol phase) (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NO3an1        => null() ! Nitrate size bin 001 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NO3an2        => null() ! Nitrate size bin 002 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: NO3an3        => null() ! Nitrate size bin 003 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: XHNO3         => null() ! buffer for NITRATE_HNO3 (kg m-2 s-1)


 contains
    procedure:: gocart2G_allocate   => NI2G_StateSpecsInit
    procedure:: gocart2G_deallocate => NI2G_StateSpecsFinalize

 end type NI2G_State


 contains


!=================================================================================================================
 subroutine NI2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(NI2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(.not.associated(self%lwi)          ) allocate(self%lwi(its:ite,jts:jte)                 )
 if(.not.associated(self%tropp)        ) allocate(self%tropp(its:ite,jts:jte)               )
 if(.not.associated(self%ustar)        ) allocate(self%tropp(its:ite,jts:jte)               )
 if(.not.associated(self%zpbl)         ) allocate(self%zpbl(its:ite,jts:jte)                )
 if(.not.associated(self%sh)           ) allocate(self%zpbl(its:ite,jts:jte)                )
 if(.not.associated(self%z0h)          ) allocate(self%zpbl(its:ite,jts:jte)                )
 if(.not.associated(self%cn_prcp)      ) allocate(self%cn_prcp(its:ite,jts:jte)             )
 if(.not.associated(self%ncn_prcp)     ) allocate(self%ncn_prcp(its:ite,jts:jte)            )
 if(.not.associated(self%area)         ) allocate(self%area(its:ite,jts:jte)                )
!.................................................................................................................
 if(.not.associated(self%airdens)      ) allocate(self%airdens(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%delp)         ) allocate(self%delp(its:ite,jts:jte,kts:kte)        )
 if(.not.associated(self%t)            ) allocate(self%t(its:ite,jts:jte,kts:kte)           )
 if(.not.associated(self%rh2)          ) allocate(self%rh2(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%zle)          ) allocate(self%zle(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%ple)          ) allocate(self%ple(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%pfl_lsan)     ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%pfi_lsan)     ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%u)            ) allocate(self%u(its:ite,jts:jte,kts:kte)           )
 if(.not.associated(self%v)            ) allocate(self%v(its:ite,jts:jte,kts:kte)           )
!.................................................................................................................
 if(.not.associated(self%emi_nh3_ag)   ) allocate(self%emi_nh3_ag(its:ite,jts:jte)          )
 if(.not.associated(self%emi_nh3_bb)   ) allocate(self%emi_nh3_bb(its:ite,jts:jte)          )
 if(.not.associated(self%emi_nh3_en)   ) allocate(self%emi_nh3_en(its:ite,jts:jte)          )
 if(.not.associated(self%emi_nh3_in)   ) allocate(self%emi_nh3_in(its:ite,jts:jte)          )
 if(.not.associated(self%emi_nh3_oc)   ) allocate(self%emi_nh3_oc(its:ite,jts:jte)          )
 if(.not.associated(self%emi_nh3_re)   ) allocate(self%emi_nh3_re(its:ite,jts:jte)          )
 if(.not.associated(self%emi_nh3_tr)   ) allocate(self%emi_nh3_tr(its:ite,jts:jte)          )
 if(.not.associated(self%nitrate_hno3) ) allocate(self%nitrate_hno3(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%du)           ) allocate(self%du(its:ite,jts:jte,kts:kte,5)        )
 if(.not.associated(self%ss)           ) allocate(self%ss(its:ite,jts:jte,kts:kte,5)        )
 if(.not.associated(self%so4)          ) allocate(self%so4(its:ite,jts:jte,kts:kte)         )

!category: EXPORT
 if(.not.associated(self%nh3mass)      ) allocate(self%nh3mass(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%nh4mass)      ) allocate(self%nh4mass(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%nimass)       ) allocate(self%nimass(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%nimass25)     ) allocate(self%nimass25(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%nh3conc)      ) allocate(self%nh3conc(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%nh4conc)      ) allocate(self%nh4conc(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%niconc)       ) allocate(self%niconc(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%niconc25)     ) allocate(self%niconc25(its:ite,jts:jte,kts:kte)    )
!if(.not.associated(self%niextcoef)    ) allocate(self%niextcoef(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)    )
!if(.not.associated(self%niextcoefrh20)) allocate(self%niextcoefrh20(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile))
!if(.not.associated(self%niextcoefrh80)) allocate(self%niextcoefrh80(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile))
!if(.not.associated(self%niscacoef)    ) allocate(self%niscacoef(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)    )
!if(.not.associated(self%niscacoefrh20)) allocate(self%niscacoefrh20(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile))
!if(.not.associated(self%niscacoefrh80)) allocate(self%niscacoefrh80(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile))
!if(.not.associated(self%nibckcoef)    ) allocate(self%nibckcoef(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)    )
!.................................................................................................................
 if(.not.associated(self%nipno3aq)     ) allocate(self%nipno3aq(its:ite,jts:jte)            )
 if(.not.associated(self%nipnh4aq)     ) allocate(self%nipnh4aq(its:ite,jts:jte)            )
 if(.not.associated(self%nipnh3aq)     ) allocate(self%nipnh3aq(its:ite,jts:jte)            )
 if(.not.associated(self%niht)         ) allocate(self%niht(its:ite,jts:jte,3)              )
 if(.not.associated(self%nisd)         ) allocate(self%nisd(its:ite,jts:jte,3)              )
 if(.not.associated(self%nidp)         ) allocate(self%nidp(its:ite,jts:jte,3)              )
 if(.not.associated(self%niwt)         ) allocate(self%niwt(its:ite,jts:jte,3)              )
 if(.not.associated(self%nisv)         ) allocate(self%nisv(its:ite,jts:jte,3)              )
 if(.not.associated(self%nh3em)        ) allocate(self%nh3em(its:ite,jts:jte)               )
 if(.not.associated(self%nh3dp)        ) allocate(self%nh3dp(its:ite,jts:jte)               )
 if(.not.associated(self%nh3wt)        ) allocate(self%nh3wt(its:ite,jts:jte)               )
 if(.not.associated(self%nh3sv)        ) allocate(self%nh3sv(its:ite,jts:jte)               )
 if(.not.associated(self%nh4sd)        ) allocate(self%nh4sd(its:ite,jts:jte)               )
 if(.not.associated(self%nh4dp)        ) allocate(self%nh4dp(its:ite,jts:jte)               )
 if(.not.associated(self%nh4dp)        ) allocate(self%nh4dp(its:ite,jts:jte)               )
 if(.not.associated(self%nh4sv)        ) allocate(self%nh4sv(its:ite,jts:jte)               ) 
 if(.not.associated(self%hno3smass)    ) allocate(self%hno3smass(its:ite,jts:jte)           )
 if(.not.associated(self%nh3smass)     ) allocate(self%nh3smass(its:ite,jts:jte)            )
 if(.not.associated(self%nh4smass)     ) allocate(self%nh4smass(its:ite,jts:jte)            )
 if(.not.associated(self%nismass)      ) allocate(self%nismass(its:ite,jts:jte)             )
 if(.not.associated(self%nismass25)    ) allocate(self%nismass25(its:ite,jts:jte)           )
 if(.not.associated(self%hno3cmass)    ) allocate(self%hno3cmass(its:ite,jts:jte)           )
 if(.not.associated(self%nh3cmass)     ) allocate(self%nh3cmass(its:ite,jts:jte)            )
 if(.not.associated(self%nh4cmass)     ) allocate(self%nh4cmass(its:ite,jts:jte)            )
 if(.not.associated(self%nicmass)      ) allocate(self%nicmass(its:ite,jts:jte)             )
 if(.not.associated(self%nicmass25)    ) allocate(self%nicmass25(its:ite,jts:jte)           )
!if(.not.associated(self%niexttfm)     ) allocate(self%niexttfm(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%niscatfm)     ) allocate(self%niscatfm(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%niextt25)     ) allocate(self%niextt25(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%niscat25)     ) allocate(self%niscat25(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%niexttau)     ) allocate(self%niexttau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%nistexttau)   ) allocate(self%nistexttau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%niscatau)     ) allocate(self%niscatau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%nistscatau)   ) allocate(self%nistscatau(its:ite,jts:jte,size(self%wavelengths_vertint)))
 if(.not.associated(self%niangstr)     ) allocate(self%niangstr(its:ite,jts:jte)            )
 if(.not.associated(self%nifluxu)      ) allocate(self%nifluxu(its:ite,jts:jte)             )
 if(.not.associated(self%nifluxv)      ) allocate(self%nifluxv(its:ite,jts:jte)             )

!category: INTERNAL
 if(.not.associated(self%nh3)          ) allocate(self%nh3(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%nh4a)         ) allocate(self%nh4a(its:ite,jts:jte,kts:kte)        )
 if(.not.associated(self%no3an1)       ) allocate(self%no3an1(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%no3an2)       ) allocate(self%no3an2(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%no3an3)       ) allocate(self%no3an3(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%xhno3)        ) allocate(self%xhno3(its:ite,jts:jte,kts:kte)       )

 end subroutine NI2G_StateSpecsInit

!=================================================================================================================
 subroutine NI2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(NI2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%lwi)          ) deallocate(self%lwi          )
 if(associated(self%tropp)        ) deallocate(self%tropp        )
 if(associated(self%ustar)        ) deallocate(self%tropp        )
 if(associated(self%zpbl)         ) deallocate(self%zpbl         )
 if(associated(self%sh)           ) deallocate(self%zpbl         )
 if(associated(self%z0h)          ) deallocate(self%z0h          )
 if(associated(self%cn_prcp)      ) deallocate(self%cn_prcp      )
 if(associated(self%ncn_prcp)     ) deallocate(self%ncn_prcp     )
 if(associated(self%area)         ) deallocate(self%area         )
!.................................................................................................................
 if(associated(self%airdens)      ) deallocate(self%airdens      )
 if(associated(self%delp)         ) deallocate(self%delp         )
 if(associated(self%t)            ) deallocate(self%t            )
 if(associated(self%rh2)          ) deallocate(self%rh2          )
 if(associated(self%zle)          ) deallocate(self%zle          )
 if(associated(self%ple)          ) deallocate(self%ple          )
 if(associated(self%pfl_lsan)     ) deallocate(self%pfl_lsan     )
 if(associated(self%pfi_lsan)     ) deallocate(self%pfi_lsan     )
 if(associated(self%u)            ) deallocate(self%u            )
 if(associated(self%v)            ) deallocate(self%v            )
!.................................................................................................................
 if(associated(self%emi_nh3_ag)   ) deallocate(self%emi_nh3_ag   )
 if(associated(self%emi_nh3_bb)   ) deallocate(self%emi_nh3_bb   )
 if(associated(self%emi_nh3_en)   ) deallocate(self%emi_nh3_en   )
 if(associated(self%emi_nh3_in)   ) deallocate(self%emi_nh3_in   )
 if(associated(self%emi_nh3_oc)   ) deallocate(self%emi_nh3_oc   )
 if(associated(self%emi_nh3_re)   ) deallocate(self%emi_nh3_re   )
 if(associated(self%emi_nh3_tr)   ) deallocate(self%emi_nh3_tr   )
 if(associated(self%nitrate_hno3) ) deallocate(self%nitrate_hno3 )
 if(associated(self%du)           ) deallocate(self%du           )
 if(associated(self%ss)           ) deallocate(self%ss           )
 if(associated(self%so4)          ) deallocate(self%so4          )

!category: EXPORT
 if(associated(self%nh3mass)      ) deallocate(self%nh3mass      )
 if(associated(self%nh4mass)      ) deallocate(self%nh4mass      )
 if(associated(self%nimass)       ) deallocate(self%nimass       )
 if(associated(self%nimass25)     ) deallocate(self%nimass25     )
 if(associated(self%hno3conc)     ) deallocate(self%hno3conc     )
 if(associated(self%nh3conc)      ) deallocate(self%nh3conc      )
 if(associated(self%nh4conc)      ) deallocate(self%nh4conc      )
 if(associated(self%niconc)       ) deallocate(self%niconc       )
 if(associated(self%niconc25)     ) deallocate(self%niconc25     )
!if(associated(self%niextcoef)    ) deallocate(self%niextcoef    )
!if(associated(self%niextcoefrh20)) deallocate(self%niextcoefrh20)
!if(associated(self%niextcoefrh80)) deallocate(self%niextcoefrh80)
!if(associated(self%niscacoef)    ) deallocate(self%niscaoef     )
!if(associated(self%niscacoefrh20)) deallocate(self%niscacoefrh20)
!if(associated(self%niscacoefrh80)) deallocate(self%niscacoefrh80)
!if(associated(self%nibckcoef)    ) deallocate(self%nibckcoef    )
!.................................................................................................................
 if(associated(self%nipno3aq)     ) deallocate(self%nipno3aq     )
 if(associated(self%nipnh4aq)     ) deallocate(self%nipnh4aq     )
 if(associated(self%nipnh3aq)     ) deallocate(self%nipnh3aq     )
 if(associated(self%niht)         ) deallocate(self%niht         )
 if(associated(self%nisd)         ) deallocate(self%nisd         )
 if(associated(self%nidp)         ) deallocate(self%nidp         )
 if(associated(self%niwt)         ) deallocate(self%niwt         )
 if(associated(self%nisv)         ) deallocate(self%nisv         )
 if(associated(self%nh3em)        ) deallocate(self%nh3em        )
 if(associated(self%nh3dp)        ) deallocate(self%nh3dp        )
 if(associated(self%nh3wt)        ) deallocate(self%nh3wt        )
 if(associated(self%nh3sv)        ) deallocate(self%nh3sv        )
 if(associated(self%nh4sd)        ) deallocate(self%nh4sd        )
 if(associated(self%nh4dp)        ) deallocate(self%nh4dp        )
 if(associated(self%nh4dp)        ) deallocate(self%nh4dp        )
 if(associated(self%nh4sv)        ) deallocate(self%nh4sv        ) 
 if(associated(self%hno3smass)    ) deallocate(self%hno3smass    )
 if(associated(self%nh3smass)     ) deallocate(self%nh3smass     )
 if(associated(self%nh4smass)     ) deallocate(self%nh4smass     )
 if(associated(self%nismass)      ) deallocate(self%nismass      )
 if(associated(self%nismass25)    ) deallocate(self%nismass25    )
 if(associated(self%hno3cmass)    ) deallocate(self%hno3cmass    )
 if(associated(self%nh3cmass)     ) deallocate(self%nh3cmass     )
 if(associated(self%nh4cmass)     ) deallocate(self%nh4cmass     )
 if(associated(self%nicmass)      ) deallocate(self%nicmass      )
 if(associated(self%nicmass25)    ) deallocate(self%nicmass25    )
!if(associated(self%niexttfm)     ) deallocate(self%niexttfm     )
!if(associated(self%niscatfm)     ) deallocate(self%niscatfm     )
!if(associated(self%niextt25)     ) deallocate(self%niextt25     )
!if(associated(self%niscat25)     ) deallocate(self%niscat25     )
!if(associated(self%niexttau)     ) deallocate(self%niexttau     )
!if(associated(self%nistexttau)   ) deallocate(self%nistexttau   )
!if(associated(self%niscatau)     ) deallocate(self%niscatau     )
!if(associated(self%nistscatau)   ) deallocate(self%nistscatau   )
 if(associated(self%niangstr)     ) deallocate(self%niangstr     )
 if(associated(self%nifluxu)      ) deallocate(self%nifluxu      )
 if(associated(self%nifluxv)      ) deallocate(self%nifluxv      )

!category: INTERNAL
 if(associated(self%nh3)          ) deallocate(self%nh3          )
 if(associated(self%nh4a)         ) deallocate(self%nh4a         )
 if(associated(self%no3an1)       ) deallocate(self%no3an1       )
 if(associated(self%no3an2)       ) deallocate(self%no3an2       )
 if(associated(self%no3an3)       ) deallocate(self%no3an3       )
 if(associated(self%xhno3)        ) deallocate(self%xhno3        )

 end subroutine NI2G_StateSpecsFinalize

!=================================================================================================================
 end module NI2G_StateSpecs
!=================================================================================================================
