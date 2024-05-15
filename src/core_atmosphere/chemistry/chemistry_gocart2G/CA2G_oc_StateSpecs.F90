!=================================================================================================================
 module CA2G_oc_StateSpecs
 use mpas_kind_types,only: RKIND

 use CA2G_oc_instance_CA,only: nbins

 implicit none
 public

!this module is the state variable specification file for carbon parameters. it is the same as CA2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!schema_version: 2.0.0
!component: CA


 type CA2G_oc_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer  :: FROCEAN          => null() ! fraction_of_ocean (-)
 real(kind=RKIND),dimension(:,:),pointer  :: FRACI            => null() ! ice_covered_fraction_of_tile (-)
 real(kind=RKIND),dimension(:,:),pointer  :: LWI              => null() ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),pointer  :: TROPP            => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer  :: U10M             => null() ! 10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: V10M             => null() ! 10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: USTAR            => null() ! surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: FRLAKE           => null() ! fraction_of_lake (-)
 real(kind=RKIND),dimension(:,:),pointer  :: AREA             => null() ! agrid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer  :: ZPBL             => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer  :: SH               => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer  :: Z0H              => null() ! surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer  :: CN_PRCP          => null() ! surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: NCN_PRCP         => null() ! non-convective precipitation (kg/m^2/s)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: AIRDENS          => null() ! moist_air_density (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer:: DELP             => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: DELZ             => null() ! geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: T                => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: RH2              => null() ! Rel_Hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: ZLE              => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: PLE              => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: PFL_LSAN         => null() ! 3D_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: PFI_LSAN         => null() ! 3D_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: U                => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: V                => null() ! northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: OC_AIRCRAFT      => null() ! aircraft emissions (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSOA_ANTHRO_VOC  => null() ! SOA from Anthropogenic and biomass burning VOC (kg m-3 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSOA_BIOB_VOC    => null() ! SOA from Anthropogenic and biomass burning VOC (kg m-3 s-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: OC_BIOMASS       => null() ! biomass burning emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_ISOPRENE      => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_MTPA          => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_MTPO          => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_LIMO          => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_BIOFUEL       => null() ! biofuel emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_ANTEOC1       => null() ! anthropogenic BF emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_ANTEOC2       => null() ! anthropogenic FF emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_SHIP          => null() ! ship emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_AVIATION_LTO  => null() ! Landing/Take-off aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_AVIATION_CDS  => null() ! Climb/Descent aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: OC_AVIATION_CRS  => null() ! Cruise aircraft source species (-)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocMASS         => null() ! Organic carbon aerosol Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocCONC         => null() ! Organic carbon aerosol Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocEXTCOEF      => null() ! Organic carbon aerosol Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocEXTCOEFRH20  => null() ! Organic carbon aerosol Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocEXTCOEFRH80  => null() ! Organic carbon aerosol Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocSCACOEF      => null() ! Organic carbon aerosol Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocSCACOEFRH20  => null() ! Organic carbon aerosol Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocSCACOEFRH80  => null() ! Organic carbon aerosol Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocBCKCOEF      => null() ! Organic carbon aerosol Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: ocEM             => null() ! Organic carbon aerosol Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocSD             => null() ! Organic carbon aerosol Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocDP             => null() ! Organic carbon aerosol Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocWT             => null() ! Organic carbon aerosol Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocSV             => null() ! Organic carbon aerosol Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: ocEMAN             => null() ! Organic carbon aerosol Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: ocEMBB             => null() ! Organic carbon aerosol Biomass Burning Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocEMBF           => null() ! Organic carbon aerosol Biofuel Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocEMBG           => null() ! Organic carbon aerosol Biogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocHYPHIL         => null() ! Organic carbon aerosol Hydrophobic to Hydrophilic (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocPSOA           => null() ! Organic carbon aerosol SOA Production (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocSMASS          => null() ! Organic carbon aerosol Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: ocCMASS          => null() ! Organic carbon aerosol Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocEXTTAU         => null() ! Organic carbon aerosol Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocSTEXTTAU       => null() ! Organic carbon aerosol Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocSCATAU         => null() ! Organic carbon aerosol Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: ocSTSCATAU       => null() ! Organic carbon aerosol Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer  :: ocANGSTR         => null() ! Organic carbon aerosol Angstrom parameter [470-870 nm] (-)=> null() !
 real(kind=RKIND),dimension(:,:),pointer  :: ocFLUXU          => null() ! Organic carbon aerosol column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocFLUXV          => null() ! Organic carbon aerosol column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ocAERIDX         => null() ! Organic carbon aerosol TOMS UV Aerosol Index (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocPHOBIC       => null() ! Hydrophobic organic aerosol Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocPHILIC       => null() ! Hydrophilic organic aerosol Mixing Ratio (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => CA2G_oc_StateSpecsInit
    procedure:: gocart2G_deallocate => CA2G_oc_StateSpecsFinalize

 end type CA2G_oc_State


 contains


!=================================================================================================================
 subroutine CA2G_oc_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_oc_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(.not.associated(self%frocean)         ) allocate(self%frocean(its:ite,jts:jte)                )
 if(.not.associated(self%fraci)           ) allocate(self%fraci(its:ite,jts:jte)                  )
 if(.not.associated(self%frlake)          ) allocate(self%frlake(its:ite,jts:jte)                 )
 if(.not.associated(self%lwi)             ) allocate(self%lwi(its:ite,jts:jte)                    )
 if(.not.associated(self%tropp)           ) allocate(self%tropp(its:ite,jts:jte)                  )
 if(.not.associated(self%u10m)            ) allocate(self%u10m(its:ite,jts:jte)                   )
 if(.not.associated(self%v10m)            ) allocate(self%v10m(its:ite,jts:jte)                   )
 if(.not.associated(self%area)            ) allocate(self%area(its:ite,jts:jte)                   )
 if(.not.associated(self%zpbl)            ) allocate(self%zpbl(its:ite,jts:jte)                   )
 if(.not.associated(self%ustar)           ) allocate(self%ustar(its:ite,jts:jte)                  )
 if(.not.associated(self%sh)              ) allocate(self%sh(its:ite,jts:jte)                     )
 if(.not.associated(self%z0h)             ) allocate(self%z0h(its:ite,jts:jte)                    )
 if(.not.associated(self%cn_prcp)         ) allocate(self%cn_prcp(its:ite,jts:jte)                )
 if(.not.associated(self%ncn_prcp)        ) allocate(self%ncn_prcp(its:ite,jts:jte)               )
!........................................ .........................................................................
 if(.not.associated(self%airdens)         ) allocate(self%airdens(its:ite,jts:jte,kts:kte)        )
 if(.not.associated(self%delp)            ) allocate(self%delp(its:ite,jts:jte,kts:kte)           )
 if(.not.associated(self%delz)            ) allocate(self%delz(its:ite,jts:jte,kts:kte)           )
 if(.not.associated(self%t)               ) allocate(self%t(its:ite,jts:jte,kts:kte)              )
 if(.not.associated(self%rh2)             ) allocate(self%rh2(its:ite,jts:jte,kts:kte)            )
 if(.not.associated(self%zle)             ) allocate(self%zle(its:ite,jts:jte,kts:kte)            )
 if(.not.associated(self%ple)             ) allocate(self%ple(its:ite,jts:jte,kts:kte)            )
 if(.not.associated(self%pfl_lsan)        ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%pfi_lsan)        ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%u)               ) allocate(self%u(its:ite,jts:jte,kts:kte)              )
 if(.not.associated(self%v)               ) allocate(self%v(its:ite,jts:jte,kts:kte)              )
 if(.not.associated(self%oc_aircraft)     ) allocate(self%oc_aircraft(its:ite,jts:jte,kts:kte)    ) 
 if(.not.associated(self%psoa_anthro_voc) ) allocate(self%psoa_anthro_voc(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%psoa_biob_voc)   ) allocate(self%psoa_biob_voc(its:ite,jts:jte,kts:kte)  )
!.................................................................................................................
 if(.not.associated(self%oc_biomass)      ) allocate(self%oc_biomass(its:ite,jts:jte)             )
 if(.not.associated(self%oc_isoprene)     ) allocate(self%oc_isoprene(its:ite,jts:jte)            )
 if(.not.associated(self%oc_mtpa)         ) allocate(self%oc_mtpa(its:ite,jts:jte)                )
 if(.not.associated(self%oc_mtpo)         ) allocate(self%oc_mtpo(its:ite,jts:jte)                )
 if(.not.associated(self%oc_limo)         ) allocate(self%oc_limo(its:ite,jts:jte)                )
 if(.not.associated(self%oc_biofuel)      ) allocate(self%oc_biofuel(its:ite,jts:jte)             )
 if(.not.associated(self%oc_anteoc1)      ) allocate(self%oc_anteoc1(its:ite,jts:jte)             )
 if(.not.associated(self%oc_anteoc2)      ) allocate(self%oc_anteoc2(its:ite,jts:jte)             )
 if(.not.associated(self%oc_ship)         ) allocate(self%oc_ship(its:ite,jts:jte)                )
 if(.not.associated(self%oc_aviation_lto) ) allocate(self%oc_aviation_lto(its:ite,jts:jte)        )
 if(.not.associated(self%oc_aviation_cds) ) allocate(self%oc_aviation_cds(its:ite,jts:jte)        )
 if(.not.associated(self%oc_aviation_crs) ) allocate(self%oc_aviation_crs(its:ite,jts:jte)        )

!category: EXPORT
 if(.not.associated(self%ocmass)          ) allocate(self%ocmass(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%occonc)          ) allocate(self%occonc(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%ocextcoef)       ) allocate(self%ocextcoef(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%ocextcoefrh20)   ) allocate(self%ocextcoefrh20(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%ocextcoefrh80)   ) allocate(self%ocextcoefrh80(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%ocscacoef)       ) allocate(self%ocscacoef(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%ocscacoefrh20)   ) allocate(self%ocscacoefrh20(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%ocscacoefrh80)   ) allocate(self%ocscacoefrh80(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%ocbckcoef)       ) allocate(self%ocbckcoef(its:ite,jts:jte,kts:kte)      )
!.................................................................................................................
 if(.not.associated(self%ocem)            ) allocate(self%ocem(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocsd)            ) allocate(self%ocsd(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocdp)            ) allocate(self%ocdp(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocwt)            ) allocate(self%ocwt(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocsv)            ) allocate(self%ocsv(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%oceman)          ) allocate(self%oceman(its:ite,jts:jte)                 )
 if(.not.associated(self%ocembb)          ) allocate(self%ocembb(its:ite,jts:jte)                 )
 if(.not.associated(self%ocembf)          ) allocate(self%ocembf(its:ite,jts:jte)                 )
 if(.not.associated(self%ocembg)          ) allocate(self%ocembg(its:ite,jts:jte)                 )
 if(.not.associated(self%ochyphil)        ) allocate(self%ochyphil(its:ite,jts:jte)               )
 if(.not.associated(self%ocpsoa)          ) allocate(self%ocpsoa(its:ite,jts:jte)                 )
 if(.not.associated(self%ocsmass)         ) allocate(self%ocsmass(its:ite,jts:jte)                )
 if(.not.associated(self%occmass)         ) allocate(self%occmass(its:ite,jts:jte)                )
!if(.not.associated(self%ocexttau)        ) allocate(self%ocexttau(its:ite,jts:jte)               )
!if(.not.associated(self%ocstexttau)      ) allocate(self%ocstexttau(its:ite,jts:jte)             )
!if(.not.associated(self%ocscatau)        ) allocate(self%ocscatau(its:ite,jts:jte)               )
!if(.not.associated(self%ocstscatau)      ) allocate(self%ocstscatau(its:ite,jts:jte)             )
 if(.not.associated(self%ocangstr)        ) allocate(self%ocangstr(its:ite,jts:jte)               )
 if(.not.associated(self%ocfluxu)         ) allocate(self%ocfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%ocfluxu)         ) allocate(self%ocfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%ocaeridx)        ) allocate(self%ocaeridx(its:ite,jts:jte)               )

!category: INTERNAL
 if(.not.associated(self%ocPHOBIC)        ) allocate(self%ocPHOBIC(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%ocPHILIC)        ) allocate(self%ocPHILIC(its:ite,jts:jte,kts:kte)       )

 end subroutine CA2G_oc_StateSpecsInit

!=================================================================================================================
 subroutine CA2G_oc_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(CA2G_oc_State),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%frocean)         ) deallocate(self%frocean         )
 if(associated(self%fraci)           ) deallocate(self%fraci           )
 if(associated(self%frlake)          ) deallocate(self%frlake          )
 if(associated(self%lwi)             ) deallocate(self%lwi             )
 if(associated(self%tropp)           ) deallocate(self%tropp           )
 if(associated(self%u10m)            ) deallocate(self%u10m            )
 if(associated(self%v10m)            ) deallocate(self%v10m            )
 if(associated(self%area)            ) deallocate(self%area            )
 if(associated(self%zpbl)            ) deallocate(self%zpbl            )
 if(associated(self%ustar)           ) deallocate(self%ustar           )
 if(associated(self%sh)              ) deallocate(self%sh              )
 if(associated(self%z0h)             ) deallocate(self%z0h             )
 if(associated(self%cn_prcp)         ) deallocate(self%cn_prcp         )
 if(associated(self%ncn_prcp)        ) deallocate(self%ncn_prcp        )
!........................................ .........................................................................
 if(associated(self%airdens)         ) deallocate(self%airdens         )
 if(associated(self%delp)            ) deallocate(self%delp            )
 if(associated(self%delz)            ) deallocate(self%delz            )
 if(associated(self%t)               ) deallocate(self%t               )
 if(associated(self%rh2)             ) deallocate(self%rh2             )
 if(associated(self%zle)             ) deallocate(self%zle             )
 if(associated(self%ple)             ) deallocate(self%ple             )
 if(associated(self%pfl_lsan)        ) deallocate(self%pfl_lsan        )
 if(associated(self%pfi_lsan)        ) deallocate(self%pfi_lsan        )
 if(associated(self%u)               ) deallocate(self%u               )
 if(associated(self%v)               ) deallocate(self%v               )
 if(associated(self%oc_aircraft)     ) deallocate(self%oc_aircraft     ) 
 if(associated(self%psoa_anthro_voc) ) deallocate(self%psoa_anthro_voc )
 if(associated(self%psoa_biob_voc)   ) deallocate(self%psoa_biob_voc   )
!.................................................................................................................
 if(associated(self%oc_biomass)      ) deallocate(self%oc_biomass      )
 if(associated(self%oc_isoprene)     ) deallocate(self%oc_isoprene     )
 if(associated(self%oc_mtpa)         ) deallocate(self%oc_mtpa         )
 if(associated(self%oc_mtpo)         ) deallocate(self%oc_mtpo         )
 if(associated(self%oc_limo)         ) deallocate(self%oc_limo         )
 if(associated(self%oc_biofuel)      ) deallocate(self%oc_biofuel      )
 if(associated(self%oc_anteoc1)      ) deallocate(self%oc_anteoc1      )
 if(associated(self%oc_anteoc2)      ) deallocate(self%oc_anteoc2      )
 if(associated(self%oc_ship)         ) deallocate(self%oc_ship         )
 if(associated(self%oc_aviation_lto) ) deallocate(self%oc_aviation_lto )
 if(associated(self%oc_aviation_cds) ) deallocate(self%oc_aviation_cds )
 if(associated(self%oc_aviation_crs) ) deallocate(self%oc_aviation_crs )

!category: EXPORT
 if(associated(self%ocmass)          ) deallocate(self%ocmass          )
 if(associated(self%occonc)          ) deallocate(self%occonc          )
!if(associated(self%ocextcoef)       ) deallocate(self%ocextcoef       ) 
!if(associated(self%ocextcoefrh20)   ) deallocate(self%ocextcoefrh20   )
!if(associated(self%ocextcoefrh80)   ) deallocate(self%ocextcoefrh80   )
!if(associated(self%ocscacoef)       ) deallocate(self%ocscacoef       )
!if(associated(self%ocscacoefrh20)   ) deallocate(self%ocscacoefrh20   )
!if(associated(self%ocscacoefrh80)   ) deallocate(self%ocscacoefrh80   )
!if(associated(self%ocbckcoef)       ) deallocate(self%ocbckcoef       )
!.................................................................................................................
 if(associated(self%ocem)            ) deallocate(self%ocem            )
 if(associated(self%ocsd)            ) deallocate(self%ocsd            )
 if(associated(self%ocdp)            ) deallocate(self%ocdp            )
 if(associated(self%ocwt)            ) deallocate(self%ocwt            )
 if(associated(self%ocsv)            ) deallocate(self%ocsv            )
 if(associated(self%oceman)          ) deallocate(self%oceman          )
 if(associated(self%ocembb)          ) deallocate(self%ocembb          )
 if(associated(self%ocembf)          ) deallocate(self%ocembf          )
 if(associated(self%ocembg)          ) deallocate(self%ocembg          )
 if(associated(self%ochyphil)        ) deallocate(self%ochyphil        )
 if(associated(self%ocpsoa)          ) deallocate(self%ocpsoa          )
 if(associated(self%ocsmass)         ) deallocate(self%ocsmass         )
 if(associated(self%occmass)         ) deallocate(self%occmass         )
 if(associated(self%ocexttau)        ) deallocate(self%ocexttau        )
 if(associated(self%ocstexttau)      ) deallocate(self%ocstexttau      )
 if(associated(self%ocscatau)        ) deallocate(self%ocscatau        )
 if(associated(self%ocstscatau)      ) deallocate(self%ocstscatau      )
 if(associated(self%ocangstr)        ) deallocate(self%ocangstr        )
 if(associated(self%ocfluxu)         ) deallocate(self%ocfluxu         )
 if(associated(self%ocfluxu)         ) deallocate(self%ocfluxu         )
 if(associated(self%ocaeridx)        ) deallocate(self%ocaeridx        )

!category: INTERNAL
 if(associated(self%ocPHOBIC)        ) deallocate(self%ocPHOBIC        )
 if(associated(self%ocPHILIC)        ) deallocate(self%ocPHILIC        )

 end subroutine CA2G_oc_StateSpecsFinalize

!=================================================================================================================
 end module CA2G_oc_StateSpecs
!=================================================================================================================
