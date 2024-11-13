! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module CA2G_br_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use CA2G_br_instance,only: nbins

 implicit none
 public

!this module is the state variable specification file for carbon parameters. it is the same as CA2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!schema_version: 2.0.0
!component: CA


 type CA2G_br_State

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
 real(kind=RKIND),dimension(:,:),pointer  :: br_biomass       => null() ! biomass burning emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_biofuel       => null() ! biofuel emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_antebr1       => null() ! anthropogenic bf emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_antebr2       => null() ! anthropogenic ff emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_ship          => null() ! ship emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_aviation_lto  => null() ! landing/take-off aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_aviation_cds  => null() ! climb/descent aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_aviation_crs  => null() ! cruise aircraft source species (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: br_aircraft      => null() ! aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: br_terpene       => null() ! terpene emissions (-)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: brmass         => null() ! brown carbon aerosol mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brconc         => null() ! brown carbon aerosol mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brextcoef      => null() ! brown carbon aerosol extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brextcoefrh20  => null() ! brown carbon aerosol extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brextcoefrh80  => null() ! brown carbon aerosol extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brscacoef      => null() ! brown carbon aerosol scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brscacoefrh20  => null() ! brown carbon aerosol scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brscacoefrh80  => null() ! brown carbon aerosol scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: brbckcoef      => null() ! brown carbon aerosol backscatter coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer  :: brem           => null() ! brown carbon aerosol emission (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brsd           => null() ! brown carbon aerosol sedimentation (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brdp           => null() ! brown carbon aerosol dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brwt           => null() ! brown carbon aerosol wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brsv           => null() ! brown carbon aerosol convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: breman         => null() ! brown carbon aerosol anthropogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brembb         => null() ! brown carbon aerosol biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brembf         => null() ! brown carbon aerosol biofuel emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brembg         => null() ! brown carbon aerosol biogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brhyphil       => null() ! brown carbon aerosol hydrophobic to hydrophilic (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brpsoa         => null() ! brown carbon aerosol soa production (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brsmass        => null() ! brown carbon aerosol surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: brcmass        => null() ! brown carbon aerosol column mass density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brexttau       => null() ! brown carbon aerosol extinction aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brstexttau     => null() ! brown carbon aerosol extinction aot stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brscatau       => null() ! brown carbon aerosol scattering aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brstscatau     => null() ! brown carbon aerosol scattering aot stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer    :: brangstr       => null() ! brown carbon aerosol angstrom parameter [470-870 nm] (-)w
 real(kind=RKIND),dimension(:,:),pointer    :: brfluxu        => null() ! brown carbon aerosol column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: brfluxv        => null() ! brown carbon aerosol column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: braeridx       => null() ! brown carbon aerosol toms uv aerosol index (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: brphobic       => null() ! Hydrophobic brown carbon aerosol mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: brphilic       => null() ! Hydrophilic brown carbon aerosol mixing Ratio (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => CA2G_br_StateSpecsInit
    procedure:: gocart2G_deallocate => CA2G_br_StateSpecsFinalize

 end type CA2G_br_State


 contains


!==================================================================================================================
 subroutine CA2G_br_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_br_State),intent(inout):: self

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
 if(.not.associated(self%br_biomass)     ) allocate(self%br_biomass(its:ite,jts:jte)             )
 if(.not.associated(self%br_biofuel)     ) allocate(self%br_biofuel(its:ite,jts:jte)             )
 if(.not.associated(self%br_antebr1)     ) allocate(self%br_antebr1(its:ite,jts:jte)             )
 if(.not.associated(self%br_antebr2)     ) allocate(self%br_antebr2(its:ite,jts:jte)             )
 if(.not.associated(self%br_ship)        ) allocate(self%br_ship(its:ite,jts:jte)                )
 if(.not.associated(self%br_aviation_lto)) allocate(self%br_aviation_lto(its:ite,jts:jte)        )
 if(.not.associated(self%br_aviation_cds)) allocate(self%br_aviation_cds(its:ite,jts:jte)        )
 if(.not.associated(self%br_aviation_crs)) allocate(self%br_aviation_crs(its:ite,jts:jte)        )
 if(.not.associated(self%br_aircraft)    ) allocate(self%br_aircraft(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%br_terpene)     ) allocate(self%br_terpene(its:ite,jts:jte)             )

!category: EXPORT
 if(.not.associated(self%brmass)         ) allocate(self%brmass(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%brconc)         ) allocate(self%brconc(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%brextcoef)      ) allocate(self%brextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%brextcoefrh20)  ) allocate(self%brextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%brextcoefrh80)  ) allocate(self%brextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%brscacoef)      ) allocate(self%brscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%brscacoefrh20)  ) allocate(self%brscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%brscacoefrh80)  ) allocate(self%brscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%brbckcoef)      ) allocate(self%brbckcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
!..................................................................................................................
 if(.not.associated(self%brem)           ) allocate(self%brem(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brsd)           ) allocate(self%brsd(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brdp)           ) allocate(self%brdp(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brwt)           ) allocate(self%brwt(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%brsv)           ) allocate(self%brsv(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%breman)         ) allocate(self%breman(its:ite,jts:jte)                 )
 if(.not.associated(self%brembb)         ) allocate(self%brembb(its:ite,jts:jte)                 )
 if(.not.associated(self%brembf)         ) allocate(self%brembf(its:ite,jts:jte)                 )
 if(.not.associated(self%brembg)         ) allocate(self%brembg(its:ite,jts:jte)                 )
 if(.not.associated(self%brhyphil)       ) allocate(self%brhyphil(its:ite,jts:jte)               )
 if(.not.associated(self%brpsoa)         ) allocate(self%brpsoa(its:ite,jts:jte)                 )
 if(.not.associated(self%brsmass)        ) allocate(self%brsmass(its:ite,jts:jte)                )
 if(.not.associated(self%brcmass)        ) allocate(self%brcmass(its:ite,jts:jte)                )
 if(.not.associated(self%brexttau)       ) allocate(self%brexttau(its:ite,jts:jte,nw_vertint)    )
 if(.not.associated(self%brstexttau)     ) allocate(self%brstexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%brscatau)       ) allocate(self%brscatau(its:ite,jts:jte,nw_vertint)    )
 if(.not.associated(self%brstscatau)     ) allocate(self%brstscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%brangstr)       ) allocate(self%brangstr(its:ite,jts:jte)               )
 if(.not.associated(self%brfluxu)        ) allocate(self%brfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%brfluxv)        ) allocate(self%brfluxv(its:ite,jts:jte)                )
 if(.not.associated(self%braeridx)       ) allocate(self%braeridx(its:ite,jts:jte)               )

!category: INTERNAL
 if(.not.associated(self%brphobic)       ) allocate(self%brphobic(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%brphilic)       ) allocate(self%brphilic(its:ite,jts:jte,kts:kte)       )

 end subroutine CA2G_br_StateSpecsInit

!==================================================================================================================
 subroutine CA2G_br_StateSpecsFinalize(self)
!==================================================================================================================

!--- inout arguments:
 class(CA2G_br_State),intent(inout) :: self

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
!........................................ .........................................................................
 if(associated(self%br_biomass)     ) deallocate(self%br_biomass     )
 if(associated(self%br_biofuel)     ) deallocate(self%br_biofuel     )
 if(associated(self%br_antebr1)     ) deallocate(self%br_antebr1     )
 if(associated(self%br_antebr2)     ) deallocate(self%br_antebr2     )
 if(associated(self%br_ship)        ) deallocate(self%br_ship        )
 if(associated(self%br_aviation_lto)) deallocate(self%br_aviation_lto)
 if(associated(self%br_aviation_cds)) deallocate(self%br_aviation_cds)
 if(associated(self%br_aviation_crs)) deallocate(self%br_aviation_crs)
 if(associated(self%br_aircraft)    ) deallocate(self%br_aircraft    )
 if(associated(self%br_terpene)     ) deallocate(self%br_terpene     )

!category: EXPORT
 if(associated(self%brmass)         ) deallocate(self%brmass         )
 if(associated(self%brconc)         ) deallocate(self%brconc         )
 if(associated(self%brextcoef)      ) deallocate(self%brextcoef      )
 if(associated(self%brextcoefrh20)  ) deallocate(self%brextcoefrh20  )
 if(associated(self%brextcoefrh80)  ) deallocate(self%brextcoefrh80  )
 if(associated(self%brscacoef)      ) deallocate(self%brscacoef      )
 if(associated(self%brscacoefrh20)  ) deallocate(self%brscacoefrh20  )
 if(associated(self%brscacoefrh80)  ) deallocate(self%brscacoefrh80  )
 if(associated(self%brbckcoef)      ) deallocate(self%brbckcoef      )
!..................................................................................................................
 if(associated(self%brem)           ) deallocate(self%brem           )
 if(associated(self%brsd)           ) deallocate(self%brsd           )
 if(associated(self%brdp)           ) deallocate(self%brdp           )
 if(associated(self%brwt)           ) deallocate(self%brwt           )
 if(associated(self%brsv)           ) deallocate(self%brsv           )
 if(associated(self%breman)         ) deallocate(self%breman         )
 if(associated(self%brembb)         ) deallocate(self%brembb         )
 if(associated(self%brembf)         ) deallocate(self%brembf         )
 if(associated(self%brembg)         ) deallocate(self%brembg         )
 if(associated(self%brhyphil)       ) deallocate(self%brhyphil       )
 if(associated(self%brpsoa)         ) deallocate(self%brpsoa         )
 if(associated(self%brsmass)        ) deallocate(self%brsmass        )
 if(associated(self%brcmass)        ) deallocate(self%brcmass        )
 if(associated(self%brexttau)       ) deallocate(self%brexttau       )
 if(associated(self%brstexttau)     ) deallocate(self%brstexttau     )
 if(associated(self%brscatau)       ) deallocate(self%brscatau       )
 if(associated(self%brstscatau)     ) deallocate(self%brstscatau     )
 if(associated(self%brangstr)       ) deallocate(self%brangstr       )
 if(associated(self%brfluxu)        ) deallocate(self%brfluxu        )
 if(associated(self%brfluxv)        ) deallocate(self%brfluxv        )
 if(associated(self%braeridx)       ) deallocate(self%braeridx       )

!category: INTERNAL
 if(associated(self%brphobic)       ) deallocate(self%brphobic       )
 if(associated(self%brphilic)       ) deallocate(self%brphilic       )

 end subroutine CA2G_br_StateSpecsFinalize

!==================================================================================================================
 end module CA2G_br_StateSpecs
!==================================================================================================================
