!=================================================================================================================
 module DU2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use DU2G_instance_DU,only: nbins

 implicit none 
 public

!this module is the state variable specification file for dust parameters. it is the same as DU2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/DU2G_GridComp.

!schema_version: 2.0.0
!component: DU


 type DU2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer  :: DU_SRC         => null() ! erod - dust emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_Z0          => null() ! aerodynamic_surface_roughness_for_aeolian_processes (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_GVF         => null() ! GVF (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_SAND        => null() ! volume_fraction_of_sand_in_soil (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_SILT        => null() ! volume_fraction_of_silt_in_soil (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_CLAY        => null() ! volume_fraction_of_clay_in_soil (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_RDRAG       => null() ! drag_partition (m-1)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_SSM         => null() ! sediment_supply_map (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_UTHRES      => null() ! surface_dry_threshold_velocity (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: FRSNOW         => null() ! surface_snow_area_fraction (-)
 real(kind=RKIND),dimension(:,:),pointer  :: SLC            => null() ! liquid_water_content_of_soil_layer (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_TEXTURE     => null() ! soil_texture (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_VEG         => null() ! vegetation_type (-)
 real(kind=RKIND),dimension(:,:),pointer  :: FRLAKE         => null() ! fraction_of_lake (-)
 real(kind=RKIND),dimension(:,:),pointer  :: FRLAND         => null() ! fraction_of_land (-)
 real(kind=RKIND),dimension(:,:),pointer  :: ASNOW          => null() ! snow_covered_fraction_of_land (-)
 real(kind=RKIND),dimension(:,:),pointer  :: WET1           => null() ! surface_soil_wetness (-)
 real(kind=RKIND),dimension(:,:),pointer  :: LWI            => null() ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),pointer  :: TROPP          => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer  :: U10M           => null() ! 10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: V10M           => null() ! 10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: U10N           => null() ! equivalent_neutral_10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: V10N           => null() ! equivalent_neutral_10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: AREA           => null() ! agrid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer  :: USTAR          => null() ! equivalent_neutral_10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: CN_PRCP        => null() ! surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: NCN_PRCP       => null() ! Non-convective precipitation (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: ZPBL           => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer  :: SH             => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer  :: Z0H            => null() ! surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer  :: WCSF           => null() ! water_surface_layer (m3 m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: TSOIL1         => null() ! soil_temperatures_layer_1 (k)
 real(kind=RKIND),dimension(:,:),pointer  :: RHOS           => null() ! air_density_at_surface (kg m-3)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: AIRDENS        => null() ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer:: DELP           => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: RH2            => null() ! Rel_Hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: T              => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: ZLE            => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: PLE            => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: PFL_LSAN       => null() ! 3D_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: PFI_LSAN       => null() ! 3D_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: U              => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: V              => null() ! northward_wind (m s-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  ::DUMASS        => null() ! Dust Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  ::DUMASS25      => null() ! Dust Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  ::DUCONC        => null() ! Dust Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUEXTCOEF     => null() ! Dust Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUEXTCOEFRH20 => null() ! Dust Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUEXTCOEFRH80 => null() ! Dust Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUSCACOEF     => null() ! Dust Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUSCACOEFRH20 => null() ! Dust Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUSCACOEFRH80 => null() ! Dust Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::DUBCKCOEF     => null() ! Dust Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: DUSMASS        => null() ! Dust Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: DUCMASS        => null() ! Dust Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUEXTTAU       => null() ! Dust Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSTEXTTAU     => null() ! Dust Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSCATAU       => null() ! Dust Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSTSCATAU     => null() ! Dust Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DUSMASS25      => null() ! Dust Surface Mass Concentration - PM 2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: DUCMASS25      => null() ! Dust Column Mass Density - PM 2.5 (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUEXTT25       => null() ! Dust Extinction AOT - PM 2.5 (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSCAT25       => null() ! Dust Scattering AOT - PM 2.5 (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DUAERIDX       => null() ! Dust TOMS UV Aerosol Index (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DUFLUXU        => null() ! Dust column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: DUFLUXV        => null() ! Dust column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUEXTTFM       => null() ! Dust Extinction AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSCATFM       => null() ! Dust Scattering AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DUANGSTR       => null() ! Dust Angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUEM           => null() ! Dust Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSD           => null() ! Dust Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUDP           => null() ! Dust Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUWT           => null() ! Dust Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: DUSV           => null() ! Dust Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_UST         => null() ! aeolian_friction_velocity (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_UST_T       => null() ! aeolian_threshold_friction_velocity (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_UST_TS      => null() ! aeolian_threshold_friction_velocity_over_smooth_surface (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_DPC         => null() ! aeolian_drag_partition_correction (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_SMC         => null() ! aeolian_soil_moisture_correction (-)
 real(kind=RKIND),dimension(:,:),pointer  :: DU_EROD        => null() ! aeolian_erodibilitiy (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:,:),pointer:: DU           => null() ! Dust Mixing Ratio (Bin %d) (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => DU2G_StateSpecsInit
    procedure:: gocart2G_deallocate => DU2G_StateSpecsFinalize

 end type DU2G_State


 contains


!=================================================================================================================
 subroutine DU2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(DU2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(.not.associated(self%du_src)       ) allocate(self%du_src(its:ite,jts:jte)          )
 if(.not.associated(self%du_z0)        ) allocate(self%du_z0(its:ite,jts:jte)           )
 if(.not.associated(self%du_gvf)       ) allocate(self%du_gvf(its:ite,jts:jte)          )
 if(.not.associated(self%du_sand)      ) allocate(self%du_sand(its:ite,jts:jte)         )
 if(.not.associated(self%du_silt)      ) allocate(self%du_silt(its:ite,jts:jte)         )
 if(.not.associated(self%du_clay)      ) allocate(self%du_clay(its:ite,jts:jte)         )
 if(.not.associated(self%du_rdrag)     ) allocate(self%du_rdrag(its:ite,jts:jte)        )
 if(.not.associated(self%du_ssm)       ) allocate(self%du_ssm(its:ite,jts:jte)          )
 if(.not.associated(self%du_uthres)    ) allocate(self%du_uthres(its:ite,jts:jte)       )
 if(.not.associated(self%frsnow)       ) allocate(self%frsnow(its:ite,jts:jte)          )
 if(.not.associated(self%slc)          ) allocate(self%slc(its:ite,jts:jte)             )
 if(.not.associated(self%du_texture)   ) allocate(self%du_texture(its:ite,jts:jte)      )
 if(.not.associated(self%du_veg)       ) allocate(self%du_veg(its:ite,jts:jte)          )
 if(.not.associated(self%frlake)       ) allocate(self%frlake(its:ite,jts:jte)          )
 if(.not.associated(self%frland)       ) allocate(self%frland(its:ite,jts:jte)          )
 if(.not.associated(self%asnow)        ) allocate(self%asnow(its:ite,jts:jte)           )
 if(.not.associated(self%wet1)         ) allocate(self%wet1(its:ite,jts:jte)            )
 if(.not.associated(self%lwi)          ) allocate(self%lwi(its:ite,jts:jte)             )
 if(.not.associated(self%tropp)        ) allocate(self%tropp(its:ite,jts:jte)           )
 if(.not.associated(self%u10m)         ) allocate(self%u10m(its:ite,jts:jte)            )
 if(.not.associated(self%v10m)         ) allocate(self%v10m(its:ite,jts:jte)            )
 if(.not.associated(self%u10n)         ) allocate(self%u10n(its:ite,jts:jte)            )
 if(.not.associated(self%v10n)         ) allocate(self%v10n(its:ite,jts:jte)            )
 if(.not.associated(self%area)         ) allocate(self%area(its:ite,jts:jte)            )
 if(.not.associated(self%ustar)        ) allocate(self%ustar(its:ite,jts:jte)           )
 if(.not.associated(self%cn_prcp)      ) allocate(self%cn_prcp(its:ite,jts:jte)         )
 if(.not.associated(self%ncn_prcp)     ) allocate(self%ncn_prcp(its:ite,jts:jte)        )
 if(.not.associated(self%zpbl)         ) allocate(self%zpbl(its:ite,jts:jte)            )
 if(.not.associated(self%sh)           ) allocate(self%sh(its:ite,jts:jte)              )
 if(.not.associated(self%z0h)          ) allocate(self%z0h(its:ite,jts:jte)             )
 if(.not.associated(self%wcsf)         ) allocate(self%wcsf(its:ite,jts:jte)            )
 if(.not.associated(self%tsoil1)       ) allocate(self%tsoil1(its:ite,jts:jte)          )
 if(.not.associated(self%rhos)         ) allocate(self%rhos(its:ite,jts:jte)            )
!.................................................................................................................
 if(.not.associated(self%airdens)      ) allocate(self%airdens(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%delp)         ) allocate(self%delp(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%rh2)          ) allocate(self%rh2(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%t)            ) allocate(self%t(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%zle)          ) allocate(self%zle(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%ple)          ) allocate(self%ple(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%pfl_lsan)     ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%pfi_lsan)     ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%u)            ) allocate(self%u(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%v)            ) allocate(self%v(its:ite,jts:jte,kts:kte)       )

!category: EXPORT
 if(.not.associated(self%dumass)       ) allocate(self%dumass(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%dumass25)     ) allocate(self%dumass25(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%duconc)       ) allocate(self%duconc(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%duextcoef)    ) allocate(self%duextcoef(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(.not.associated(self%duextcoefrh20)) allocate(self%duextcoefrh20(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(.not.associated(self%duextcoefrh80)) allocate(self%duextcoefrh80(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(.not.associated(self%duscacoef)    ) allocate(self%dusacoef(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(.not.associated(self%duscacoefrh20)) allocate(self%dusacoefrh20(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(.not.associated(self%duscacoefrh80)) allocate(self%dusacoefrh80(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(.not.associated(self%dubckcoef)    ) allocate(self%dubckcoef(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!.................................................................................................................
 if(.not.associated(self%dusmass)      ) allocate(self%dusmass(its:ite,jts:jte)         )
 if(.not.associated(self%ducmass)      ) allocate(self%ducmass(its:ite,jts:jte)         )
!if(.not.associated(self%duexttau)     ) allocate(self%duexttau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%dustexttau)   ) allocate(self%dustexttau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%duscatau)     ) allocate(self%duscatau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%dustscatau)   ) allocate(self%dustscatau(its:ite,jts:jte,size(self%wavelengths_vertint))
 if(.not.associated(self%dusmass25)    ) allocate(self%dusmass25(its:ite,jts:jte)       )
 if(.not.associated(self%ducmass25)    ) allocate(self%ducmass25(its:ite,jts:jte)       )
!if(.not.associated(self%duextt25)     ) allocate(self%duextt25(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%duscat25)     ) allocate(self%duscat25(its:ite,jts:jte,size(self%wavelengths_vertint)))
 if(.not.associated(self%duaeridx)     ) allocate(self%duaeridx(its:ite,jts:jte)        )
 if(.not.associated(self%dufluxu)      ) allocate(self%dufluxu(its:ite,jts:jte)         )
 if(.not.associated(self%dufluxv)      ) allocate(self%dufluxv(its:ite,jts:jte)         )
!if(.not.associated(self%duexttfm)     ) allocate(self%duexttfm(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(.not.associated(self%duscatfm)     ) allocate(self%duscatfm(its:ite,jts:jte,size(self%wavelengths_vertint)))
 if(.not.associated(self%duangstr)     ) allocate(self%duangstr(its:ite,jts:jte)        )
 if(.not.associated(self%duem)         ) allocate(self%duem(its:ite,jts:jte,nbins)      )
 if(.not.associated(self%dusd)         ) allocate(self%dusd(its:ite,jts:jte,nbins)      )
 if(.not.associated(self%dudp)         ) allocate(self%dudp(its:ite,jts:jte,nbins)      )
 if(.not.associated(self%duwt)         ) allocate(self%duwt(its:ite,jts:jte,nbins)      )
 if(.not.associated(self%dusv)         ) allocate(self%dusv(its:ite,jts:jte,nbins)      )
 if(.not.associated(self%du_ust)       ) allocate(self%du_ust(its:ite,jts:jte)          )
 if(.not.associated(self%du_ust_t)     ) allocate(self%du_ust_t(its:ite,jts:jte)        )
 if(.not.associated(self%du_ust_ts)    ) allocate(self%du_ust_ts(its:ite,jts:jte)       )
 if(.not.associated(self%du_dpc)       ) allocate(self%du_dpc(its:ite,jts:jte)          )
 if(.not.associated(self%du_smc)       ) allocate(self%du_smc(its:ite,jts:jte)          )
 if(.not.associated(self%du_erod)      ) allocate(self%du_erod(its:ite,jts:jte)         )

!category: internal
 if(.not.associated(self%du)           ) allocate(self%du(its:ite,jts:jte,kts:kte,nbins))

 end subroutine DU2G_StateSpecsInit

!=================================================================================================================
 subroutine DU2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(DU2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%du_src)       ) deallocate(self%du_src    )
 if(associated(self%du_z0)        ) deallocate(self%du_z0     )
 if(associated(self%du_gvf)       ) deallocate(self%du_gvf    )
 if(associated(self%du_sand)      ) deallocate(self%du_sand   )
 if(associated(self%du_silt)      ) deallocate(self%du_silt   )
 if(associated(self%du_clay)      ) deallocate(self%du_clay   )
 if(associated(self%du_rdrag)     ) deallocate(self%du_rdrag  )
 if(associated(self%du_ssm)       ) deallocate(self%du_ssm    )
 if(associated(self%du_uthres)    ) deallocate(self%du_uthres )
 if(associated(self%frsnow)       ) deallocate(self%frsnow    )
 if(associated(self%slc)          ) deallocate(self%slc       )
 if(associated(self%du_texture)   ) deallocate(self%du_texture)
 if(associated(self%du_veg)       ) deallocate(self%du_veg    )
 if(associated(self%frlake)       ) deallocate(self%frlake    )
 if(associated(self%frland)       ) deallocate(self%frland    )
 if(associated(self%asnow)        ) deallocate(self%asnow     )
 if(associated(self%wet1)         ) deallocate(self%wet1      )
 if(associated(self%lwi)          ) deallocate(self%lwi       )
 if(associated(self%tropp)        ) deallocate(self%tropp     )
 if(associated(self%u10m)         ) deallocate(self%u10m      )
 if(associated(self%v10m)         ) deallocate(self%v10m      )
 if(associated(self%u10n)         ) deallocate(self%u10n      )
 if(associated(self%v10n)         ) deallocate(self%v10n      )
 if(associated(self%area)         ) deallocate(self%area      )
 if(associated(self%ustar)        ) deallocate(self%ustar     )
 if(associated(self%cn_prcp)      ) deallocate(self%cn_prcp   )
 if(associated(self%ncn_prcp)     ) deallocate(self%ncn_prcp  )
 if(associated(self%zpbl)         ) deallocate(self%zpbl      )
 if(associated(self%sh)           ) deallocate(self%sh        )
 if(associated(self%z0h)          ) deallocate(self%z0h       )
 if(associated(self%wcsf)         ) deallocate(self%wcsf      )
 if(associated(self%tsoil1)       ) deallocate(self%tsoil1    )
 if(associated(self%rhos)         ) deallocate(self%rhos      )
!.................................................................................................................
 if(associated(self%airdens)      ) deallocate(self%airdens   )
 if(associated(self%delp)         ) deallocate(self%delp      )
 if(associated(self%rh2)          ) deallocate(self%rh2       )
 if(associated(self%t)            ) deallocate(self%t         )
 if(associated(self%zle)          ) deallocate(self%zle       )
 if(associated(self%ple)          ) deallocate(self%ple       )
 if(associated(self%pfl_lsan)     ) deallocate(self%pfl_lsan  )
 if(associated(self%pfi_lsan)     ) deallocate(self%pfi_lsan  )
 if(associated(self%u)            ) deallocate(self%u         )
 if(associated(self%v)            ) deallocate(self%v         )

!category: EXPORT
 if(associated(self%dumass)       ) deallocate(self%dumass    )
 if(associated(self%dumass25)     ) deallocate(self%dumass25  )
 if(associated(self%duconc)       ) deallocate(self%duconc    )
!if(associated(self%duextcoef)    ) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(associated(self%duextcoefrh20)) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(associated(self%duextcoefrh80)) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(associated(self%duscacoef)    ) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(associated(self%duscacoefrh20)) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(associated(self%duscacoefrh80)) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!if(associated(self%::dubckcoef)  ) deallocate(self%(its:ite,jts:jte,kts:kte,size(self%wavelengths_profile)))
!.................................................................................................................
 if(associated(self%dusmass)      ) deallocate(self%dusmass   )
 if(associated(self%ducmass)      ) deallocate(self%ducmass   )
!if(associated(self%duexttau)     ) deallocate(self%duexttau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(associated(self%dustexttau)   ) deallocate(self%dustexttau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(associated(self%duscatau)     ) deallocate(self%duscatau(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(associated(self%dustscatau)   ) deallocate(self%dustscatau(its:ite,jts:jte,size(self%wavelengths_vertint))
 if(associated(self%dusmass25)    ) deallocate(self%dusmass25 )
 if(associated(self%ducmass25)    ) deallocate(self%ducmass25 )
!if(associated(self%duextt25)     ) deallocate(self%duextt25(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(associated(self%duscat25)     ) deallocate(self%duscat25(its:ite,jts:jte,size(self%wavelengths_vertint)))
 if(associated(self%duaeridx)     ) deallocate(self%duaeridx  )
 if(associated(self%dufluxu)      ) deallocate(self%dufluxu   )
 if(associated(self%dufluxv)      ) deallocate(self%dufluxv   )
!if(associated(self%duexttfm)     ) deallocate(self%duexttfm(its:ite,jts:jte,size(self%wavelengths_vertint)))
!if(associated(self%duscatfm)     ) deallocate(self%duscatfm(its:ite,jts:jte,size(self%wavelengths_vertint)))
 if(associated(self%duangstr)     ) deallocate(self%duangstr  )
 if(associated(self%duem)         ) deallocate(self%duem      )
 if(associated(self%dusd)         ) deallocate(self%dusd      )
 if(associated(self%dudp)         ) deallocate(self%dudp      )
 if(associated(self%duwt)         ) deallocate(self%duwt      )
 if(associated(self%dusv)         ) deallocate(self%dusv      )
 if(associated(self%du_ust)       ) deallocate(self%du_ust    )
 if(associated(self%du_ust_t)     ) deallocate(self%du_ust_t  )
 if(associated(self%du_ust_ts)    ) deallocate(self%du_ust_ts )
 if(associated(self%du_dpc)       ) deallocate(self%du_dpc    )
 if(associated(self%du_smc)       ) deallocate(self%du_smc    )
 if(associated(self%du_erod)      ) deallocate(self%du_erod   )

!category: internal
 if(associated(self%du)           ) deallocate(self%du        )

 end subroutine DU2G_StateSpecsFinalize

!=================================================================================================================
 end module DU2G_StateSpecs
!=================================================================================================================

