! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module CA2G_oc_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use CA2G_oc_instance,only: nbins

 implicit none
 public

!this module is the state variable specification file for carbon parameters. it is the same as CA2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!schema_version: 2.0.0
!component: CA


 type CA2G_oc_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer  :: lats             => null() ! latitude (radian)
 real(kind=RKIND),dimension(:,:),pointer  :: lons             => null() ! longitude (radian)
 real(kind=RKIND),dimension(:,:),pointer  :: frocean          => null() ! fraction_of_ocean (-)
 real(kind=RKIND),dimension(:,:),pointer  :: fraci            => null() ! ice_covered_fraction_of_tile (-)
 real(kind=RKIND),dimension(:,:),pointer  :: lwi              => null() ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),pointer  :: tropp            => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer  :: u10m             => null() ! 10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: v10m             => null() ! 10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ustar            => null() ! surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: frlake           => null() ! fraction_of_lake (-)
 real(kind=RKIND),dimension(:,:),pointer  :: area             => null() ! agrid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer  :: zpbl             => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer  :: sh               => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer  :: z0h              => null() ! surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer  :: cn_prcp          => null() ! surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: ncn_prcp         => null() ! non-convective precipitation (kg/m^2/s)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: airdens          => null() ! moist_air_density (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer:: delp             => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: delz             => null() ! geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: t                => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: rh2              => null() ! rel_hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: zle              => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: ple              => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfl_lsan         => null() ! 3d_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfi_lsan         => null() ! 3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: u                => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: v                => null() ! northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: psoa_anthro_voc  => null() ! soa from anthropogenic and biomass burning voc (kg m-3 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: psoa_biob_voc    => null() ! soa from anthropogenic and biomass burning voc (kg m-3 s-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: oc_biomass       => null() ! biomass burning emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_biofuel       => null() ! biofuel emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_anteoc1       => null() ! anthropogenic bf emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_anteoc2       => null() ! anthropogenic ff emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_ship          => null() ! ship emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_aviation_lto  => null() ! landing/take-off aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_aviation_cds  => null() ! climb/descent aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_aviation_crs  => null() ! cruise aircraft source species (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: oc_aircraft      => null() ! aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_isoprene      => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_mtpa          => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_mtpo          => null() ! source species (-)
 real(kind=RKIND),dimension(:,:),pointer  :: oc_limo          => null() ! source species (-)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocmass         => null() ! organic carbon aerosol mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: occonc         => null() ! organic carbon aerosol mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocextcoef      => null() ! organic carbon aerosol extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocextcoefrh20  => null() ! organic carbon aerosol extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocextcoefrh80  => null() ! organic carbon aerosol extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocscacoef      => null() ! organic carbon aerosol scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocscacoefrh20  => null() ! organic carbon aerosol scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocscacoefrh80  => null() ! organic carbon aerosol scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ocbckcoef      => null() ! organic carbon aerosol backscatter coefficient (m-1 sr-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocem           => null() ! organic carbon aerosol emission (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocsd           => null() ! organic carbon aerosol sedimentation (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocdp           => null() ! organic carbon aerosol dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocwt           => null() ! organic carbon aerosol wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocsv           => null() ! organic carbon aerosol convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: oceman         => null() ! organic carbon aerosol anthropogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocembb         => null() ! organic carbon aerosol biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocembf         => null() ! organic carbon aerosol biofuel emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocembg         => null() ! organic carbon aerosol biogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ochyphil       => null() ! organic carbon aerosol hydrophobic to hydrophilic (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocpsoa         => null() ! organic carbon aerosol soa production (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocsmass        => null() ! organic carbon aerosol surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: occmass        => null() ! organic carbon aerosol column mass density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocexttau       => null() ! organic carbon aerosol extinction aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocstexttau     => null() ! organic carbon aerosol extinction aot stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocscatau       => null() ! organic carbon aerosol scattering aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocstscatau     => null() ! organic carbon aerosol scattering aot stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer    :: ocangstr       => null() ! organic carbon aerosol angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:),pointer    :: ocfluxu        => null() ! organic carbon aerosol column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocfluxv        => null() ! organic carbon aerosol column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ocaeridx       => null() ! organic carbon aerosol toms uv aerosol index (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocphobic       => null() ! Hydrophobic organic carbon aerosol mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ocphilic       => null() ! Hydrophilic organic carbon aerosol mixing Ratio (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => CA2G_oc_StateSpecsInit
    procedure:: gocart2G_deallocate => CA2G_oc_StateSpecsFinalize

 end type CA2G_oc_State


 contains


!==================================================================================================================
 subroutine CA2G_oc_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_oc_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!------------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
 if(.not.associated(self%lats)           ) allocate(self%lats(its:ite,jts:jte)                   )
 if(.not.associated(self%lons)           ) allocate(self%lons(its:ite,jts:jte)                   )
 if(.not.associated(self%frocean)        ) allocate(self%frocean(its:ite,jts:jte)                )
 if(.not.associated(self%fraci)          ) allocate(self%fraci(its:ite,jts:jte)                  )
 if(.not.associated(self%frlake)         ) allocate(self%frlake(its:ite,jts:jte)                 )
 if(.not.associated(self%lwi)            ) allocate(self%lwi(its:ite,jts:jte)                    )
 if(.not.associated(self%tropp)          ) allocate(self%tropp(its:ite,jts:jte)                  )
 if(.not.associated(self%u10m)           ) allocate(self%u10m(its:ite,jts:jte)                   )
 if(.not.associated(self%v10m)           ) allocate(self%v10m(its:ite,jts:jte)                   )
 if(.not.associated(self%area)           ) allocate(self%area(its:ite,jts:jte)                   )
 if(.not.associated(self%zpbl)           ) allocate(self%zpbl(its:ite,jts:jte)                   )
 if(.not.associated(self%ustar)          ) allocate(self%ustar(its:ite,jts:jte)                  )
 if(.not.associated(self%sh)             ) allocate(self%sh(its:ite,jts:jte)                     )
 if(.not.associated(self%z0h)            ) allocate(self%z0h(its:ite,jts:jte)                    )
 if(.not.associated(self%cn_prcp)        ) allocate(self%cn_prcp(its:ite,jts:jte)                )
 if(.not.associated(self%ncn_prcp)       ) allocate(self%ncn_prcp(its:ite,jts:jte)               )
!........................................ .........................................................................
 if(.not.associated(self%airdens)        ) allocate(self%airdens(its:ite,jts:jte,kts:kte)        )
 if(.not.associated(self%delp)           ) allocate(self%delp(its:ite,jts:jte,kts:kte)           )
 if(.not.associated(self%delz)           ) allocate(self%delz(its:ite,jts:jte,kts:kte)           )
 if(.not.associated(self%t)              ) allocate(self%t(its:ite,jts:jte,kts:kte)              )
 if(.not.associated(self%rh2)            ) allocate(self%rh2(its:ite,jts:jte,kts:kte)            )
 if(.not.associated(self%u)              ) allocate(self%u(its:ite,jts:jte,kts:kte)              )
 if(.not.associated(self%v)              ) allocate(self%v(its:ite,jts:jte,kts:kte)              )
 if(.not.associated(self%pfl_lsan)       ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%pfi_lsan)       ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%psoa_anthro_voc)) allocate(self%psoa_anthro_voc(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%psoa_biob_voc)  ) allocate(self%psoa_biob_voc(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%zle)            ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)          )
 if(.not.associated(self%ple)            ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)          )
!..................................................................................................................
 if(.not.associated(self%oc_biomass)     ) allocate(self%oc_biomass(its:ite,jts:jte)             )
 if(.not.associated(self%oc_biofuel)     ) allocate(self%oc_biofuel(its:ite,jts:jte)             )
 if(.not.associated(self%oc_anteoc1)     ) allocate(self%oc_anteoc1(its:ite,jts:jte)             )
 if(.not.associated(self%oc_anteoc2)     ) allocate(self%oc_anteoc2(its:ite,jts:jte)             )
 if(.not.associated(self%oc_ship)        ) allocate(self%oc_ship(its:ite,jts:jte)                )
 if(.not.associated(self%oc_aviation_lto)) allocate(self%oc_aviation_lto(its:ite,jts:jte)        )
 if(.not.associated(self%oc_aviation_cds)) allocate(self%oc_aviation_cds(its:ite,jts:jte)        )
 if(.not.associated(self%oc_aviation_crs)) allocate(self%oc_aviation_crs(its:ite,jts:jte)        )
 if(.not.associated(self%oc_aircraft)    ) allocate(self%oc_aircraft(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%oc_isoprene)    ) allocate(self%oc_isoprene(its:ite,jts:jte)            )
 if(.not.associated(self%oc_mtpa)        ) allocate(self%oc_mtpa(its:ite,jts:jte)                )
 if(.not.associated(self%oc_mtpo)        ) allocate(self%oc_mtpo(its:ite,jts:jte)                )
 if(.not.associated(self%oc_limo)        ) allocate(self%oc_limo(its:ite,jts:jte)                )

!category: EXPORT
 if(.not.associated(self%ocmass)         ) allocate(self%ocmass(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%occonc)         ) allocate(self%occonc(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%ocextcoef)      ) allocate(self%ocextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%ocextcoefrh20)  ) allocate(self%ocextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ocextcoefrh80)  ) allocate(self%ocextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ocscacoef)      ) allocate(self%ocscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%ocscacoefrh20)  ) allocate(self%ocscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ocscacoefrh80)  ) allocate(self%ocscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ocbckcoef)      ) allocate(self%ocbckcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
!..................................................................................................................
 if(.not.associated(self%ocem)           ) allocate(self%ocem(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocsd)           ) allocate(self%ocsd(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocdp)           ) allocate(self%ocdp(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocwt)           ) allocate(self%ocwt(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%ocsv)           ) allocate(self%ocsv(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%oceman)         ) allocate(self%oceman(its:ite,jts:jte)                 )
 if(.not.associated(self%ocembb)         ) allocate(self%ocembb(its:ite,jts:jte)                 )
 if(.not.associated(self%ocembf)         ) allocate(self%ocembf(its:ite,jts:jte)                 )
 if(.not.associated(self%ocembg)         ) allocate(self%ocembg(its:ite,jts:jte)                 )
 if(.not.associated(self%ochyphil)       ) allocate(self%ochyphil(its:ite,jts:jte)               )
 if(.not.associated(self%ocpsoa)         ) allocate(self%ocpsoa(its:ite,jts:jte)                 )
 if(.not.associated(self%ocsmass)        ) allocate(self%ocsmass(its:ite,jts:jte)                )
 if(.not.associated(self%occmass)        ) allocate(self%occmass(its:ite,jts:jte)                )
 if(.not.associated(self%ocexttau)       ) allocate(self%ocexttau(its:ite,jts:jte,nw_vertint)    )
 if(.not.associated(self%ocstexttau)     ) allocate(self%ocstexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ocscatau)       ) allocate(self%ocscatau(its:ite,jts:jte,nw_vertint)    )
 if(.not.associated(self%ocstscatau)     ) allocate(self%ocstscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ocangstr)       ) allocate(self%ocangstr(its:ite,jts:jte)               )
 if(.not.associated(self%ocfluxu)        ) allocate(self%ocfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%ocfluxv)        ) allocate(self%ocfluxv(its:ite,jts:jte)                )
 if(.not.associated(self%ocaeridx)       ) allocate(self%ocaeridx(its:ite,jts:jte)               )

!category: INTERNAL
 if(.not.associated(self%ocphobic)       ) allocate(self%ocphobic(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%ocphilic)       ) allocate(self%ocphilic(its:ite,jts:jte,kts:kte)       )

 end subroutine CA2G_oc_StateSpecsInit

!==================================================================================================================
 subroutine CA2G_oc_StateSpecsFinalize(self)
!==================================================================================================================

!--- inout arguments:
 class(CA2G_oc_State),intent(inout) :: self

!------------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%lats)           ) deallocate(self%lats           )
 if(associated(self%lons)           ) deallocate(self%lons           )
 if(associated(self%frocean)        ) deallocate(self%frocean        )
 if(associated(self%fraci)          ) deallocate(self%fraci          )
 if(associated(self%frlake)         ) deallocate(self%frlake         )
 if(associated(self%lwi)            ) deallocate(self%lwi            )
 if(associated(self%tropp)          ) deallocate(self%tropp          )
 if(associated(self%u10m)           ) deallocate(self%u10m           )
 if(associated(self%v10m)           ) deallocate(self%v10m           )
 if(associated(self%area)           ) deallocate(self%area           )
 if(associated(self%zpbl)           ) deallocate(self%zpbl           )
 if(associated(self%ustar)          ) deallocate(self%ustar          )
 if(associated(self%sh)             ) deallocate(self%sh             )
 if(associated(self%z0h)            ) deallocate(self%z0h            )
 if(associated(self%cn_prcp)        ) deallocate(self%cn_prcp        )
 if(associated(self%ncn_prcp)       ) deallocate(self%ncn_prcp       )
!........................................ .........................................................................
 if(associated(self%airdens)        ) deallocate(self%airdens        )
 if(associated(self%delp)           ) deallocate(self%delp           )
 if(associated(self%delz)           ) deallocate(self%delz           )
 if(associated(self%t)              ) deallocate(self%t              )
 if(associated(self%rh2)            ) deallocate(self%rh2            )
 if(associated(self%u)              ) deallocate(self%u              )
 if(associated(self%v)              ) deallocate(self%v              )
 if(associated(self%pfl_lsan)       ) deallocate(self%pfl_lsan       )
 if(associated(self%pfi_lsan)       ) deallocate(self%pfi_lsan       )
 if(associated(self%psoa_anthro_voc)) deallocate(self%psoa_anthro_voc)
 if(associated(self%psoa_biob_voc)  ) deallocate(self%psoa_biob_voc  )
 if(associated(self%zle)            ) deallocate(self%zle            )
 if(associated(self%ple)            ) deallocate(self%ple            )
!..................................................................................................................
 if(associated(self%oc_biomass)     ) deallocate(self%oc_biomass     )
 if(associated(self%oc_biofuel)     ) deallocate(self%oc_biofuel     )
 if(associated(self%oc_anteoc1)     ) deallocate(self%oc_anteoc1     )
 if(associated(self%oc_anteoc2)     ) deallocate(self%oc_anteoc2     )
 if(associated(self%oc_ship)        ) deallocate(self%oc_ship        )
 if(associated(self%oc_aviation_lto)) deallocate(self%oc_aviation_lto)
 if(associated(self%oc_aviation_cds)) deallocate(self%oc_aviation_cds)
 if(associated(self%oc_aviation_crs)) deallocate(self%oc_aviation_crs)
 if(associated(self%oc_aircraft)    ) deallocate(self%oc_aircraft    )
 if(associated(self%oc_isoprene)    ) deallocate(self%oc_isoprene    )
 if(associated(self%oc_mtpa)        ) deallocate(self%oc_mtpa        )
 if(associated(self%oc_mtpo)        ) deallocate(self%oc_mtpo        )
 if(associated(self%oc_limo)        ) deallocate(self%oc_limo        )

!category: EXPORT
 if(associated(self%ocmass)         ) deallocate(self%ocmass         )
 if(associated(self%occonc)         ) deallocate(self%occonc         )
 if(associated(self%ocextcoef)      ) deallocate(self%ocextcoef      )
 if(associated(self%ocextcoefrh20)  ) deallocate(self%ocextcoefrh20  )
 if(associated(self%ocextcoefrh80)  ) deallocate(self%ocextcoefrh80  )
 if(associated(self%ocscacoef)      ) deallocate(self%ocscacoef      )
 if(associated(self%ocscacoefrh20)  ) deallocate(self%ocscacoefrh20  )
 if(associated(self%ocscacoefrh80)  ) deallocate(self%ocscacoefrh80  )
 if(associated(self%ocbckcoef)      ) deallocate(self%ocbckcoef      )
!..................................................................................................................
 if(associated(self%ocem)           ) deallocate(self%ocem           )
 if(associated(self%ocsd)           ) deallocate(self%ocsd           )
 if(associated(self%ocdp)           ) deallocate(self%ocdp           )
 if(associated(self%ocwt)           ) deallocate(self%ocwt           )
 if(associated(self%ocsv)           ) deallocate(self%ocsv           )
 if(associated(self%oceman)         ) deallocate(self%oceman         )
 if(associated(self%ocembb)         ) deallocate(self%ocembb         )
 if(associated(self%ocembf)         ) deallocate(self%ocembf         )
 if(associated(self%ocembg)         ) deallocate(self%ocembg         )
 if(associated(self%ochyphil)       ) deallocate(self%ochyphil       )
 if(associated(self%ocpsoa)         ) deallocate(self%ocpsoa         )
 if(associated(self%ocsmass)        ) deallocate(self%ocsmass        )
 if(associated(self%occmass)        ) deallocate(self%occmass        )
 if(associated(self%ocexttau)       ) deallocate(self%ocexttau       )
 if(associated(self%ocstexttau)     ) deallocate(self%ocstexttau     )
 if(associated(self%ocscatau)       ) deallocate(self%ocscatau       )
 if(associated(self%ocstscatau)     ) deallocate(self%ocstscatau     )
 if(associated(self%ocangstr)       ) deallocate(self%ocangstr       )
 if(associated(self%ocfluxu)        ) deallocate(self%ocfluxu        )
 if(associated(self%ocfluxv)        ) deallocate(self%ocfluxv        )
 if(associated(self%ocaeridx)       ) deallocate(self%ocaeridx       )

!category: INTERNAL
 if(associated(self%ocphobic)       ) deallocate(self%ocphobic       )
 if(associated(self%ocphilic)       ) deallocate(self%ocphilic       )

 end subroutine CA2G_oc_StateSpecsFinalize

!==================================================================================================================
 end module CA2G_oc_StateSpecs
!==================================================================================================================
