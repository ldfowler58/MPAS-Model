!=================================================================================================================
 module CA2G_bc_StateSpecs
 use mpas_kind_types,only: RKIND

 use CA2G_bc_instance_CA,only: nbins

 implicit none
 public

!this module is the state variable specification file for carbon parameters. it is the same as CA2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!schema_version: 2.0.0
!component: CA


 type CA2G_bc_State

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
 real(kind=RKIND),dimension(:,:,:),pointer:: BC_AIRCRAFT      => null() ! aircraft emissions (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSOA_ANTHRO_VOC  => null() ! SOA from Anthropogenic and biomass burning VOC (kg m-3 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSOA_BIOB_VOC    => null() ! SOA from Anthropogenic and biomass burning VOC (kg m-3 s-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: BC_BIOMASS       => null() ! biomass burning emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_BIOFUEL       => null() ! biofuel emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_ANTEBC1       => null() ! anthropogenic BF emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_ANTEBC2       => null() ! anthropogenic FF emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_SHIP          => null() ! ship emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_AVIATION_LTO  => null() ! Landing/Take-off aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_AVIATION_CDS  => null() ! Climb/Descent aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: BC_AVIATION_CRS  => null() ! Cruise aircraft source species (-)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcMASS         => null() ! Black carbon aerosol Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcCONC         => null() ! Black carbon aerosol Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcEXTCOEF      => null() ! Black carbon aerosol Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcEXTCOEFRH20  => null() ! Black carbon aerosol Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcEXTCOEFRH80  => null() ! Black carbon aerosol Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcSCACOEF      => null() ! Black carbon aerosol Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcSCACOEFRH20  => null() ! Black carbon aerosol Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcSCACOEFRH80  => null() ! Black carbon aerosol Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcBCKCOEF      => null() ! Black carbon aerosol Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: bcEM             => null() ! Black carbon aerosol Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcSD             => null() ! Black carbon aerosol Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcDP             => null() ! Black carbon aerosol Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcWT             => null() ! Black carbon aerosol Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcSV             => null() ! Black carbon aerosol Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcEMAN           => null() ! Black carbon aerosol Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcEMBB           => null() ! Black carbon aerosol Biomass Burning Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcEMBF           => null() ! Black carbon aerosol Biofuel Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcEMBG           => null() ! Black carbon aerosol Biogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcHYPHIL         => null() ! Black carbon aerosol Hydrophobic to Hydrophilic (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcPSOA           => null() ! Black carbon aerosol SOA Production (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcSMASS          => null() ! Black carbon aerosol Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: bcCMASS          => null() ! Black carbon aerosol Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcEXTTAU         => null() ! Black carbon aerosol Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcSTEXTTAU       => null() ! Black carbon aerosol Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcSCATAU         => null() ! Black carbon aerosol Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: bcSTSCATAU       => null() ! Black carbon aerosol Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bcANGSTR         => null() ! Black carbon aerosol Angstrom parameter [470-870 nm] (-)=> null() !
 real(kind=RKIND),dimension(:,:),pointer  :: bcFLUXU          => null() ! Black carbon aerosol column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcFLUXV          => null() ! Black carbon aerosol column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: bcAERIDX         => null() ! Black carbon aerosol TOMS UV Aerosol Index (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcPHOBIC       => null() ! Hydrophobic black aerosol Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcPHILIC       => null() ! Hydrophilic black aerosol Mixing Ratio (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => CA2G_bc_StateSpecsInit
    procedure:: gocart2G_deallocate => CA2G_bc_StateSpecsFinalize

 end type CA2G_bc_State


 contains


!=================================================================================================================
 subroutine CA2G_bc_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_bc_State),intent(inout):: self

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
 if(.not.associated(self%bc_aircraft)     ) allocate(self%bc_aircraft(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%psoa_anthro_voc) ) allocate(self%psoa_anthro_voc(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%psoa_biob_voc)   ) allocate(self%psoa_biob_voc(its:ite,jts:jte,kts:kte)  )
!.................................................................................................................
 if(.not.associated(self%bc_biomass)      ) allocate(self%bc_biomass(its:ite,jts:jte)             )
 if(.not.associated(self%bc_biofuel)      ) allocate(self%bc_biofuel(its:ite,jts:jte)             )
 if(.not.associated(self%bc_antebc1)      ) allocate(self%bc_antebc1(its:ite,jts:jte)             )
 if(.not.associated(self%bc_antebc2)      ) allocate(self%bc_antebc2(its:ite,jts:jte)             )
 if(.not.associated(self%bc_ship)         ) allocate(self%bc_ship(its:ite,jts:jte)                )
 if(.not.associated(self%bc_aviation_lto) ) allocate(self%bc_aviation_lto(its:ite,jts:jte)        )
 if(.not.associated(self%bc_aviation_cds) ) allocate(self%bc_aviation_cds(its:ite,jts:jte)        )
 if(.not.associated(self%bc_aviation_crs) ) allocate(self%bc_aviation_crs(its:ite,jts:jte)        )

!category: EXPORT
 if(.not.associated(self%bcmass)          ) allocate(self%bcmass(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%bcconc)          ) allocate(self%bcconc(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%bcextcoef)       ) allocate(self%bcextcoef(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%bcextcoefrh20)   ) allocate(self%bcextcoefrh20(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%bcextcoefrh80)   ) allocate(self%bcextcoefrh80(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%bcscacoef)       ) allocate(self%bcscacoef(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%bcscacoefrh20)   ) allocate(self%bcscacoefrh20(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%bcscacoefrh80)   ) allocate(self%bcscacoefrh80(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%bcbckcoef)       ) allocate(self%bcbckcoef(its:ite,jts:jte,kts:kte)      )
!.................................................................................................................
 if(.not.associated(self%bcem)            ) allocate(self%bcem(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcsd)            ) allocate(self%bcsd(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcdp)            ) allocate(self%bcdp(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcwt)            ) allocate(self%bcwt(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcsv)            ) allocate(self%bcsv(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bceman)          ) allocate(self%bceman(its:ite,jts:jte)                 )
 if(.not.associated(self%bcembb)          ) allocate(self%bcembb(its:ite,jts:jte)                 )
 if(.not.associated(self%bcembf)          ) allocate(self%bcembf(its:ite,jts:jte)                 )
 if(.not.associated(self%bcembg)          ) allocate(self%bcembg(its:ite,jts:jte)                 )
 if(.not.associated(self%bchyphil)        ) allocate(self%bchyphil(its:ite,jts:jte)               )
 if(.not.associated(self%bcpsoa)          ) allocate(self%bcpsoa(its:ite,jts:jte)                 )
 if(.not.associated(self%bcsmass)         ) allocate(self%bcsmass(its:ite,jts:jte)                )
 if(.not.associated(self%bccmass)         ) allocate(self%bccmass(its:ite,jts:jte)                )
!if(.not.associated(self%bcexttau)        ) allocate(self%bcexttau(its:ite,jts:jte)               )
!if(.not.associated(self%bcstexttau)      ) allocate(self%bcstexttau(its:ite,jts:jte)             )
!if(.not.associated(self%bcscatau)        ) allocate(self%bcscatau(its:ite,jts:jte)               )
!if(.not.associated(self%bcstscatau)      ) allocate(self%bcstscatau(its:ite,jts:jte)             )
 if(.not.associated(self%bcangstr)        ) allocate(self%bcangstr(its:ite,jts:jte)               )
 if(.not.associated(self%bcfluxu)         ) allocate(self%bcfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%bcfluxv)         ) allocate(self%bcfluxv(its:ite,jts:jte)                )
 if(.not.associated(self%bcaeridx)        ) allocate(self%bcaeridx(its:ite,jts:jte)               )

!category: INTERNAL
 if(.not.associated(self%bcPHOBIC)        ) allocate(self%bcPHOBIC(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%bcPHILIC)        ) allocate(self%bcPHILIC(its:ite,jts:jte,kts:kte)       )

 end subroutine CA2G_bc_StateSpecsInit

!=================================================================================================================
 subroutine CA2G_bc_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(CA2G_bc_State),intent(inout) :: self

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
 if(associated(self%t)               ) deallocate(self%t               )
 if(associated(self%rh2)             ) deallocate(self%rh2             )
 if(associated(self%zle)             ) deallocate(self%zle             )
 if(associated(self%ple)             ) deallocate(self%ple             )
 if(associated(self%pfl_lsan)        ) deallocate(self%pfl_lsan        )
 if(associated(self%pfi_lsan)        ) deallocate(self%pfi_lsan        )
 if(associated(self%u)               ) deallocate(self%u               )
 if(associated(self%v)               ) deallocate(self%v               )
 if(associated(self%bc_aircraft)     ) deallocate(self%bc_aircraft     )
 if(associated(self%psoa_anthro_voc) ) deallocate(self%psoa_anthro_voc )
 if(associated(self%psoa_biob_voc)   ) deallocate(self%psoa_biob_voc   )
!.................................................................................................................
 if(associated(self%bc_biomass)      ) deallocate(self%bc_biomass      )
 if(associated(self%bc_biofuel)      ) deallocate(self%bc_biofuel      )
 if(associated(self%bc_antebc1)      ) deallocate(self%bc_antebc1      )
 if(associated(self%bc_antebc2)      ) deallocate(self%bc_antebc2      )
 if(associated(self%bc_ship)         ) deallocate(self%bc_ship         )
 if(associated(self%bc_aviation_lto) ) deallocate(self%bc_aviation_lto )
 if(associated(self%bc_aviation_cds) ) deallocate(self%bc_aviation_cds )
 if(associated(self%bc_aviation_crs) ) deallocate(self%bc_aviation_crs )

!category: EXPORT
 if(associated(self%bcmass)          ) deallocate(self%bcmass          )
 if(associated(self%bcconc)          ) deallocate(self%bcconc          )
!if(associated(self%bcextcoef)       ) deallocate(self%bcextcoef       )
!if(associated(self%bcextcoefrh20)   ) deallocate(self%bcextcoefrh20   )
!if(associated(self%bcextcoefrh80)   ) deallocate(self%bcextcoefrh80   )
!if(associated(self%bcscacoef)       ) deallocate(self%bcscacoef       )
!if(associated(self%bcscacoefrh20)   ) deallocate(self%bcscacoefrh20   )
!if(associated(self%bcscacoefrh80)   ) deallocate(self%bcscacoefrh80   )
!if(associated(self%bcbckcoef)       ) deallocate(self%bcbckcoef       )
!.................................................................................................................
 if(associated(self%bcem)            ) deallocate(self%bcem            )
 if(associated(self%bcsd)            ) deallocate(self%bcsd            )
 if(associated(self%bcdp)            ) deallocate(self%bcdp            )
 if(associated(self%bcwt)            ) deallocate(self%bcwt            )
 if(associated(self%bcsv)            ) deallocate(self%bcsv            )
 if(associated(self%bceman)          ) deallocate(self%bceman          )
 if(associated(self%bcembb)          ) deallocate(self%bcembb          )
 if(associated(self%bcembf)          ) deallocate(self%bcembf          )
 if(associated(self%bcembg)          ) deallocate(self%bcembg          )
 if(associated(self%bchyphil)        ) deallocate(self%bchyphil        )
 if(associated(self%bcpsoa)          ) deallocate(self%bcpsoa          )
 if(associated(self%bcsmass)         ) deallocate(self%bcsmass         )
 if(associated(self%bccmass)         ) deallocate(self%bccmass         )
 if(associated(self%bcexttau)        ) deallocate(self%bcexttau        )
 if(associated(self%bcstexttau)      ) deallocate(self%bcstexttau      )
 if(associated(self%bcscatau)        ) deallocate(self%bcscatau        )
 if(associated(self%bcstscatau)      ) deallocate(self%bcstscatau      )
 if(associated(self%bcangstr)        ) deallocate(self%bcangstr        )
 if(associated(self%bcfluxu)         ) deallocate(self%bcfluxu         )
 if(associated(self%bcfluxv)         ) deallocate(self%bcfluxv         )
 if(associated(self%bcaeridx)        ) deallocate(self%bcaeridx        )

!category: INTERNAL
 if(associated(self%bcphobic)        ) deallocate(self%bcphobic        )
 if(associated(self%bcphilic)        ) deallocate(self%bcphilic        )

 end subroutine CA2G_bc_StateSpecsFinalize

!=================================================================================================================
 end module CA2G_bc_StateSpecs
!=================================================================================================================
