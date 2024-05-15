!=================================================================================================================
 module CA2G_br_StateSpecs
 use mpas_kind_types,only: RKIND

 use CA2G_br_instance_CA,only: nbins

 implicit none
 public

!this module is the state variable specification file for carbon parameters. it is the same as CA2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!schema_version: 2.0.0
!component: CA


 type CA2G_br_State

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
 real(kind=RKIND),dimension(:,:,:),pointer:: BRC_AIRCRAFT     => null() ! aircraft emissions (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSOA_ANTHRO_VOC  => null() ! SOA from Anthropogenic and biomass burning VOC (kg m-3 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSOA_BIOB_VOC    => null() ! SOA from Anthropogenic and biomass burning VOC (kg m-3 s-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_BIOMASS      => null() ! biomass burning emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_TERPENE      => null() ! terpene emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_BIOFUEL      => null() ! biofuel emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_ANTEBRC1     => null() ! anthropogenic BF emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_ANTEBRC2     => null() ! anthropogenic FF emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_SHIP         => null() ! ship emisisons (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_AVIATION_LTO => null() ! Landing/Take-off aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_AVIATION_CDS => null() ! Climb/Descent aircraft emissions(-)
 real(kind=RKIND),dimension(:,:),pointer  :: BRC_AVIATION_CRS => null() ! Cruise aircraft source species (-)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: brMASS         => null() ! Brown carbon aerosol Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brCONC         => null() ! Brown carbon aerosol Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brEXTCOEF      => null() ! Brown carbon aerosol Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brEXTCOEFRH20  => null() ! Brown carbon aerosol Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brEXTCOEFRH80  => null() ! Brown carbon aerosol Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brSCACOEF      => null() ! Brown carbon aerosol Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brSCACOEFRH20  => null() ! Brown carbon aerosol Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brSCACOEFRH80  => null() ! Brown carbon aerosol Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brBCKCOEF      => null() ! Brown carbon aerosol Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: brEM             => null() ! Brown carbon aerosol Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: brSD             => null() ! Brown carbon aerosol Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: brDP             => null() ! Brown carbon aerosol Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: brWT             => null() ! Brown carbon aerosol Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: brSV             => null() ! Brown carbon aerosol Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: brEMAN             => null() ! Brown carbon aerosol Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: brEMBB             => null() ! Brown carbon aerosol Biomass Burning Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brEMBF           => null() ! Brown carbon aerosol Biofuel Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brEMBG           => null() ! Brown carbon aerosol Biogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brHYPHIL         => null() ! Brown carbon aerosol Hydrophobic to Hydrophilic (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brPSOA           => null() ! Brown carbon aerosol SOA Production (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brSMASS          => null() ! Brown carbon aerosol Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: brCMASS          => null() ! Brown carbon aerosol Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: brEXTTAU         => null() ! Brown carbon aerosol Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: brSTEXTTAU       => null() ! Brown carbon aerosol Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: brSCATAU         => null() ! Brown carbon aerosol Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: brSTSCATAU       => null() ! Brown carbon aerosol Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer  :: brANGSTR         => null() ! Brown carbon aerosol Angstrom parameter [470-870 nm] (-)=> null() !
 real(kind=RKIND),dimension(:,:),pointer  :: brFLUXU          => null() ! Brown carbon aerosol column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brFLUXV          => null() ! Brown carbon aerosol column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: brAERIDX         => null() ! Brown carbon aerosol TOMS UV Aerosol Index (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: brPHOBIC       => null() ! Hydrophobic brown aerosol Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brPHILIC       => null() ! Hydrophilic brown aerosol Mixing Ratio (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => CA2G_br_StateSpecsInit
    procedure:: gocart2G_deallocate => CA2G_br_StateSpecsFinalize

 end type CA2G_br_State


 contains


!=================================================================================================================
 subroutine CA2G_br_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_br_State),intent(inout):: self

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
 if(.not.associated(self%brc_aircraft)    ) allocate(self%brc_aircraft(its:ite,jts:jte,kts:kte)   )
 if(.not.associated(self%psoa_anthro_voc) ) allocate(self%psoa_anthro_voc(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%psoa_biob_voc)   ) allocate(self%psoa_biob_voc(its:ite,jts:jte,kts:kte)  )
!.................................................................................................................
 if(.not.associated(self%brc_biomass)     ) allocate(self%brc_biomass(its:ite,jts:jte)            )
 if(.not.associated(self%brc_terpene)     ) allocate(self%brc_terpene(its:ite,jts:jte)            )
 if(.not.associated(self%brc_biofuel)     ) allocate(self%brc_biofuel(its:ite,jts:jte)            )
 if(.not.associated(self%brc_antebrc1)    ) allocate(self%brc_antebrc1(its:ite,jts:jte)           )
 if(.not.associated(self%brc_antebrc2)    ) allocate(self%brc_antebrc2(its:ite,jts:jte)           )
 if(.not.associated(self%brc_ship)        ) allocate(self%brc_ship(its:ite,jts:jte)               )
 if(.not.associated(self%brc_aviation_lto)) allocate(self%brc_aviation_lto(its:ite,jts:jte)       )
 if(.not.associated(self%brc_aviation_cds)) allocate(self%brc_aviation_cds(its:ite,jts:jte)       )
 if(.not.associated(self%brc_aviation_crs)) allocate(self%brc_aviation_crs(its:ite,jts:jte)       )

!category: EXPORT
 if(.not.associated(self%brmass)          ) allocate(self%brmass(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%brconc)          ) allocate(self%brconc(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%brextcoef)       ) allocate(self%brextcoef(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%brextcoefrh20)   ) allocate(self%brextcoefrh20(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%brextcoefrh80)   ) allocate(self%brextcoefrh80(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%brscacoef)       ) allocate(self%brscacoef(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%brscacoefrh20)   ) allocate(self%brscacoefrh20(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%brscacoefrh80)   ) allocate(self%brscacoefrh80(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%brbckcoef)       ) allocate(self%brbckcoef(its:ite,jts:jte,kts:kte)      )
!.................................................................................................................
 if(.not.associated(self%brem)            ) allocate(self%brem(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brsd)            ) allocate(self%brsd(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brdp)            ) allocate(self%brdp(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brwt)            ) allocate(self%brwt(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brsv)            ) allocate(self%brsv(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%breman)          ) allocate(self%breman(its:ite,jts:jte)                 )
 if(.not.associated(self%brembb)          ) allocate(self%brembb(its:ite,jts:jte)                 )
 if(.not.associated(self%brembf)          ) allocate(self%brembf(its:ite,jts:jte)                 )
 if(.not.associated(self%brembg)          ) allocate(self%brembg(its:ite,jts:jte)                 )
 if(.not.associated(self%brhyphil)        ) allocate(self%brhyphil(its:ite,jts:jte)               )
 if(.not.associated(self%brpsoa)          ) allocate(self%brpsoa(its:ite,jts:jte)                 )
 if(.not.associated(self%brsmass)         ) allocate(self%brsmass(its:ite,jts:jte)                )
 if(.not.associated(self%brcmass)         ) allocate(self%brcmass(its:ite,jts:jte)                )
!if(.not.associated(self%brexttau)        ) allocate(self%brexttau(its:ite,jts:jte)               )
!if(.not.associated(self%brstexttau)      ) allocate(self%brstexttau(its:ite,jts:jte)             )
!if(.not.associated(self%brscatau)        ) allocate(self%brscatau(its:ite,jts:jte)               )
!if(.not.associated(self%brstscatau)      ) allocate(self%brstscatau(its:ite,jts:jte)             )
 if(.not.associated(self%brangstr)        ) allocate(self%brangstr(its:ite,jts:jte)               )
 if(.not.associated(self%brfluxu)         ) allocate(self%brfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%brfluxv)         ) allocate(self%brfluxv(its:ite,jts:jte)                )
 if(.not.associated(self%braeridx)        ) allocate(self%braeridx(its:ite,jts:jte)               )

!category: INTERNAL
 if(.not.associated(self%brPHOBIC)        ) allocate(self%brPHOBIC(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%brPHILIC)        ) allocate(self%brPHILIC(its:ite,jts:jte,kts:kte)       )

 end subroutine CA2G_br_StateSpecsInit

!=================================================================================================================
 subroutine CA2G_br_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(CA2G_br_State),intent(inout) :: self

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
 if(associated(self%brc_aircraft)    ) deallocate(self%brc_aircraft    )
 if(associated(self%psoa_anthro_voc) ) deallocate(self%psoa_anthro_voc )
 if(associated(self%psoa_biob_voc)   ) deallocate(self%psoa_biob_voc   )
!.................................................................................................................
 if(associated(self%brc_biomass)     ) deallocate(self%brc_biomass     )
 if(associated(self%brc_terpene)     ) deallocate(self%brc_terpene     )
 if(associated(self%brc_biofuel)     ) deallocate(self%brc_biofuel     )
 if(associated(self%brc_antebrc1)    ) deallocate(self%brc_antebrc1    )
 if(associated(self%brc_antebrc2)    ) deallocate(self%brc_antebrc2    )
 if(associated(self%brc_ship)        ) deallocate(self%brc_ship        )
 if(associated(self%brc_aviation_lto)) deallocate(self%brc_aviation_lto)
 if(associated(self%brc_aviation_cds)) deallocate(self%brc_aviation_cds)
 if(associated(self%brc_aviation_crs)) deallocate(self%brc_aviation_crs)

!category: EXPORT
 if(associated(self%brmass)          ) deallocate(self%brmass          )
 if(associated(self%brconc)          ) deallocate(self%brconc          )
!if(associated(self%brextcoef)       ) deallocate(self%brextcoef       )
!if(associated(self%brextcoefrh20)   ) deallocate(self%brextcoefrh20   )
!if(associated(self%brextcoefrh80)   ) deallocate(self%brextcoefrh80   )
!if(associated(self%brscacoef)       ) deallocate(self%brscacoef       )
!if(associated(self%brscacoefrh20)   ) deallocate(self%brscacoefrh20   )
!if(associated(self%brscacoefrh80)   ) deallocate(self%brscacoefrh80   )
!if(associated(self%brbckcoef)       ) deallocate(self%brbckcoef       )
!.................................................................................................................
 if(associated(self%brem)            ) deallocate(self%brem            )
 if(associated(self%brsd)            ) deallocate(self%brsd            )
 if(associated(self%brdp)            ) deallocate(self%brdp            )
 if(associated(self%brwt)            ) deallocate(self%brwt            )
 if(associated(self%brsv)            ) deallocate(self%brsv            )
 if(associated(self%breman)          ) deallocate(self%breman          )
 if(associated(self%brembb)          ) deallocate(self%brembb          )
 if(associated(self%brembf)          ) deallocate(self%brembf          )
 if(associated(self%brembg)          ) deallocate(self%brembg          )
 if(associated(self%brhyphil)        ) deallocate(self%brhyphil        )
 if(associated(self%brpsoa)          ) deallocate(self%brpsoa          )
 if(associated(self%brsmass)         ) deallocate(self%brsmass         )
 if(associated(self%brcmass)         ) deallocate(self%brcmass         )
 if(associated(self%brexttau)        ) deallocate(self%brexttau        )
 if(associated(self%brstexttau)      ) deallocate(self%brstexttau      )
 if(associated(self%brscatau)        ) deallocate(self%brscatau        )
 if(associated(self%brstscatau)      ) deallocate(self%brstscatau      )
 if(associated(self%brangstr)        ) deallocate(self%brangstr        )
 if(associated(self%brfluxu)         ) deallocate(self%brfluxu         )
 if(associated(self%brfluxv)         ) deallocate(self%brfluxv         )
 if(associated(self%braeridx)        ) deallocate(self%braeridx        )

!category: INTERNAL
 if(associated(self%brPHOBIC)        ) deallocate(self%brPHOBIC        )
 if(associated(self%brPHILIC)        ) deallocate(self%brPHILIC        )

 end subroutine CA2G_br_StateSpecsFinalize

!=================================================================================================================
 end module CA2G_br_StateSpecs
!=================================================================================================================
