!==================================================================================================================
 module DU2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use DU2G_instance,only: nbins

 implicit none
 public

!this module is the state variable specification file for dust parameters. it is the same as DU2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/DU2G_GridComp.

!schema_version: 2.0.0
!component: DU


 type DU2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:,:),pointer:: du_src         => null() ! erod - dust emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_z0          => null() ! aerodynamic_surface_roughness_for_aeolian_processes (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_gvf         => null() ! gvf (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_sand        => null() ! volume_fraction_of_sand_in_soil (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_silt        => null() ! volume_fraction_of_silt_in_soil (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_clay        => null() ! volume_fraction_of_clay_in_soil (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_rdrag       => null() ! drag_partition (m-1)
 real(kind=RKIND),dimension(:,:),pointer  :: du_ssm         => null() ! sediment_supply_map (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_uthres      => null() ! surface_dry_threshold_velocity (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: frsnow         => null() ! surface_snow_area_fraction (-)
 real(kind=RKIND),dimension(:,:),pointer  :: slc            => null() ! liquid_water_content_of_soil_layer (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_texture     => null() ! soil_texture (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_veg         => null() ! vegetation_type (-)
 real(kind=RKIND),dimension(:,:),pointer  :: frlake         => null() ! fraction_of_lake (-)
 real(kind=RKIND),dimension(:,:),pointer  :: frland         => null() ! fraction_of_land (-)
 real(kind=RKIND),dimension(:,:),pointer  :: asnow          => null() ! snow_covered_fraction_of_land (-)
 real(kind=RKIND),dimension(:,:),pointer  :: wet1           => null() ! surface_soil_wetness (-)
 real(kind=RKIND),dimension(:,:),pointer  :: lwi            => null() ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),pointer  :: tropp          => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer  :: u10m           => null() ! 10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: v10m           => null() ! 10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: u10n           => null() ! equivalent_neutral_10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: v10n           => null() ! equivalent_neutral_10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: area           => null() ! agrid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer  :: ustar          => null() ! equivalent_neutral_10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: cn_prcp        => null() ! surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: ncn_prcp       => null() ! non-convective precipitation (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: zpbl           => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer  :: sh             => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer  :: z0h            => null() ! surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer  :: wcsf           => null() ! water_surface_layer (m3 m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: tsoil1         => null() ! soil_temperatures_layer_1 (k)
 real(kind=RKIND),dimension(:,:),pointer  :: rhos           => null() ! air_density_at_surface (kg m-3)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: airdens        => null() ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer:: delp           => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: delz           => null() ! geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: rh2            => null() ! rel_hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: t              => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: zle            => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: ple            => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfl_lsan       => null() ! 3d_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfi_lsan       => null() ! 3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: u              => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: v              => null() ! northward_wind (m s-1)

!category: export
 real(kind=RKIND),dimension(:,:,:),pointer  ::dumass        => null() ! dust mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  ::dumass25      => null() ! dust mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  ::duconc        => null() ! dust mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer::duextcoef     => null() ! dust extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::duextcoefrh20 => null() ! dust extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::duextcoefrh80 => null() ! dust extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::duscacoef     => null() ! dust scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::duscacoefrh20 => null() ! dust scattering coefficient - fixed Rh=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::duscacoefrh80 => null() ! dust scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer::dubckcoef     => null() ! dust backscatter coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: dusmass        => null() ! dust surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: ducmass        => null() ! dust column mass density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: duexttau       => null() ! dust extinction aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: dustexttau     => null() ! dust extinction aot stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: duscatau       => null() ! dust scattering aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: dustscatau     => null() ! dust scattering aot stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer  :: dusmass25      => null() ! dust surface mass concentration - pm 2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: ducmass25      => null() ! dust column mass density - pm 2.5 (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: duextt25       => null() ! dust extinction aot - pm 2.5 (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: duscat25       => null() ! dust scattering aot - pm 2.5 (-)
 real(kind=RKIND),dimension(:,:),pointer  :: duaeridx       => null() ! dust toms uv aerosol index (-)
 real(kind=RKIND),dimension(:,:),pointer  :: dufluxu        => null() ! dust column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: dufluxv        => null() ! dust column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: duexttfm       => null() ! dust extinction aot - pm 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: duscatfm       => null() ! dust scattering aot - pm 1.0 um (-)
 real(kind=RKIND),dimension(:,:),pointer  :: duangstr       => null() ! dust angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: duem           => null() ! dust emission (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: dusd           => null() ! dust sedimentation (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: dudp           => null() ! dust dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: duwt           => null() ! dust wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: dusv           => null() ! dust convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: du_ust         => null() ! aeolian_friction_velocity (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_ust_t       => null() ! aeolian_threshold_friction_velocity (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_ust_ts      => null() ! aeolian_threshold_friction_velocity_over_smooth_surface (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_dpc         => null() ! aeolian_drag_partition_correction (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_smc         => null() ! aeolian_soil_moisture_correction (-)
 real(kind=RKIND),dimension(:,:),pointer  :: du_erod        => null() ! aeolian_erodibilitiy (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:,:),pointer:: DU           => null() ! dust mixing Ratio (Bin %d) (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => DU2G_StateSpecsInit
    procedure:: gocart2G_deallocate => DU2G_StateSpecsFinalize

 end type DU2G_State


 contains


!==================================================================================================================
 subroutine DU2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte,nerod)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nerod

!--- inout arguments:
 class(DU2G_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!------------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
 if(.not.associated(self%du_src)       ) allocate(self%du_src(its:ite,jts:jte,nerod)    )
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
 if(.not.associated(self%delz)         ) allocate(self%delz(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%rh2)          ) allocate(self%rh2(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%t)            ) allocate(self%t(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%zle)          ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)   )
 if(.not.associated(self%ple)          ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)   )
 if(.not.associated(self%pfl_lsan)     ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%pfi_lsan)     ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%u)            ) allocate(self%u(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%v)            ) allocate(self%v(its:ite,jts:jte,kts:kte)       )

!category: EXPORT
 if(.not.associated(self%dumass)       ) allocate(self%dumass(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%dumass25)     ) allocate(self%dumass25(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%duconc)       ) allocate(self%duconc(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%duextcoef)    ) allocate(self%duextcoef(its:ite,jts:jte,kts:kte,nw_profile)     )
 if(.not.associated(self%duextcoefrh20)) allocate(self%duextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile) )
 if(.not.associated(self%duextcoefrh80)) allocate(self%duextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile) )
 if(.not.associated(self%duscacoef)    ) allocate(self%duscacoef(its:ite,jts:jte,kts:kte,nw_profile)     )
 if(.not.associated(self%duscacoefrh20)) allocate(self%duscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile) )
 if(.not.associated(self%duscacoefrh80)) allocate(self%duscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile) )
 if(.not.associated(self%dubckcoef)    ) allocate(self%dubckcoef(its:ite,jts:jte,kts:kte,nw_profile)     )  
!.................................................................................................................
 if(.not.associated(self%dusmass)      ) allocate(self%dusmass(its:ite,jts:jte)         )
 if(.not.associated(self%ducmass)      ) allocate(self%ducmass(its:ite,jts:jte)         )
 if(.not.associated(self%duexttau)     ) allocate(self%duexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%dustexttau)   ) allocate(self%dustexttau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%duscatau)     ) allocate(self%duscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%dustscatau)   ) allocate(self%dustscatau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%dusmass25)    ) allocate(self%dusmass25(its:ite,jts:jte)       )
 if(.not.associated(self%ducmass25)    ) allocate(self%ducmass25(its:ite,jts:jte)       )
 if(.not.associated(self%duextt25)     ) allocate(self%duextt25(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%duscat25)     ) allocate(self%duscat25(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%duaeridx)     ) allocate(self%duaeridx(its:ite,jts:jte)        )
 if(.not.associated(self%dufluxu)      ) allocate(self%dufluxu(its:ite,jts:jte)         )
 if(.not.associated(self%dufluxv)      ) allocate(self%dufluxv(its:ite,jts:jte)         )
 if(.not.associated(self%duexttfm)     ) allocate(self%duexttfm(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%duscatfm)     ) allocate(self%duscatfm(its:ite,jts:jte,nw_vertint))
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

!category: INTERNAL
 if(.not.associated(self%du)           ) allocate(self%du(its:ite,jts:jte,kts:kte,nbins))

 end subroutine DU2G_StateSpecsInit

!==================================================================================================================
 subroutine DU2G_StateSpecsFinalize(self)
!==================================================================================================================

!--- inout arguments:
 class(DU2G_State),intent(inout):: self

!------------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%du_src)       ) deallocate(self%du_src       )
 if(associated(self%du_z0)        ) deallocate(self%du_z0        )
 if(associated(self%du_gvf)       ) deallocate(self%du_gvf       )
 if(associated(self%du_sand)      ) deallocate(self%du_sand      )
 if(associated(self%du_silt)      ) deallocate(self%du_silt      )
 if(associated(self%du_clay)      ) deallocate(self%du_clay      )
 if(associated(self%du_rdrag)     ) deallocate(self%du_rdrag     )
 if(associated(self%du_ssm)       ) deallocate(self%du_ssm       )
 if(associated(self%du_uthres)    ) deallocate(self%du_uthres    )
 if(associated(self%frsnow)       ) deallocate(self%frsnow       )
 if(associated(self%slc)          ) deallocate(self%slc          )
 if(associated(self%du_texture)   ) deallocate(self%du_texture   )
 if(associated(self%du_veg)       ) deallocate(self%du_veg       )
 if(associated(self%frlake)       ) deallocate(self%frlake       )
 if(associated(self%frland)       ) deallocate(self%frland       )
 if(associated(self%asnow)        ) deallocate(self%asnow        )
 if(associated(self%wet1)         ) deallocate(self%wet1         )
 if(associated(self%lwi)          ) deallocate(self%lwi          )
 if(associated(self%tropp)        ) deallocate(self%tropp        )
 if(associated(self%u10m)         ) deallocate(self%u10m         )
 if(associated(self%v10m)         ) deallocate(self%v10m         )
 if(associated(self%u10n)         ) deallocate(self%u10n         )
 if(associated(self%v10n)         ) deallocate(self%v10n         )
 if(associated(self%area)         ) deallocate(self%area         )
 if(associated(self%ustar)        ) deallocate(self%ustar        )
 if(associated(self%cn_prcp)      ) deallocate(self%cn_prcp      )
 if(associated(self%ncn_prcp)     ) deallocate(self%ncn_prcp     )
 if(associated(self%zpbl)         ) deallocate(self%zpbl         )
 if(associated(self%sh)           ) deallocate(self%sh           )
 if(associated(self%z0h)          ) deallocate(self%z0h          )
 if(associated(self%wcsf)         ) deallocate(self%wcsf         )
 if(associated(self%tsoil1)       ) deallocate(self%tsoil1       )
 if(associated(self%rhos)         ) deallocate(self%rhos         )
!..................................................................................................................
 if(associated(self%airdens)      ) deallocate(self%airdens      )
 if(associated(self%delp)         ) deallocate(self%delp         )
 if(associated(self%delz)         ) deallocate(self%delz         )
 if(associated(self%rh2)          ) deallocate(self%rh2          )
 if(associated(self%t)            ) deallocate(self%t            )
 if(associated(self%zle)          ) deallocate(self%zle          )
 if(associated(self%ple)          ) deallocate(self%ple          )
 if(associated(self%pfl_lsan)     ) deallocate(self%pfl_lsan     )
 if(associated(self%pfi_lsan)     ) deallocate(self%pfi_lsan     )
 if(associated(self%u)            ) deallocate(self%u            )
 if(associated(self%v)            ) deallocate(self%v            )

!category: EXPORT
 if(associated(self%dumass)       ) deallocate(self%dumass       )
 if(associated(self%dumass25)     ) deallocate(self%dumass25     )
 if(associated(self%duconc)       ) deallocate(self%duconc       )
 if(associated(self%duextcoef)    ) deallocate(self%duextcoef    )
 if(associated(self%duextcoefrh20)) deallocate(self%duextcoefrh20)
 if(associated(self%duextcoefrh80)) deallocate(self%duextcoefrh80)
 if(associated(self%duscacoef)    ) deallocate(self%duscacoef    )
 if(associated(self%duscacoefrh20)) deallocate(self%duscacoefrh20)
 if(associated(self%duscacoefrh80)) deallocate(self%duscacoefrh80)
 if(associated(self%dubckcoef)    ) deallocate(self%dubckcoef    )
!..................................................................................................................
 if(associated(self%dusmass)      ) deallocate(self%dusmass      )
 if(associated(self%ducmass)      ) deallocate(self%ducmass      )
 if(associated(self%duexttau)     ) deallocate(self%duexttau     )
 if(associated(self%dustexttau)   ) deallocate(self%dustexttau   )
 if(associated(self%duscatau)     ) deallocate(self%duscatau     )
 if(associated(self%dustscatau)   ) deallocate(self%dustscatau   )
 if(associated(self%dusmass25)    ) deallocate(self%dusmass25    )
 if(associated(self%ducmass25)    ) deallocate(self%ducmass25    )
 if(associated(self%duextt25)     ) deallocate(self%duextt25     )
 if(associated(self%duscat25)     ) deallocate(self%duscat25     )
 if(associated(self%duaeridx)     ) deallocate(self%duaeridx     )
 if(associated(self%dufluxu)      ) deallocate(self%dufluxu      )
 if(associated(self%dufluxv)      ) deallocate(self%dufluxv      )
 if(associated(self%duexttfm)     ) deallocate(self%duexttfm     )
 if(associated(self%duscatfm)     ) deallocate(self%duscatfm     )
 if(associated(self%duangstr)     ) deallocate(self%duangstr     )
 if(associated(self%duem)         ) deallocate(self%duem         )
 if(associated(self%dusd)         ) deallocate(self%dusd         )
 if(associated(self%dudp)         ) deallocate(self%dudp         )
 if(associated(self%duwt)         ) deallocate(self%duwt         )
 if(associated(self%dusv)         ) deallocate(self%dusv         )
 if(associated(self%du_ust)       ) deallocate(self%du_ust       )
 if(associated(self%du_ust_t)     ) deallocate(self%du_ust_t     )
 if(associated(self%du_ust_ts)    ) deallocate(self%du_ust_ts    )
 if(associated(self%du_dpc)       ) deallocate(self%du_dpc       )
 if(associated(self%du_smc)       ) deallocate(self%du_smc       )
 if(associated(self%du_erod)      ) deallocate(self%du_erod      )

!category: internal
 if(associated(self%du)           ) deallocate(self%du           )

 end subroutine DU2G_StateSpecsFinalize

!==================================================================================================================
 end module DU2G_StateSpecs
!==================================================================================================================
