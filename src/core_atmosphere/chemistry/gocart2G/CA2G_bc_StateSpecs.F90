! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module CA2G_bc_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use CA2G_bc_instance,only: nbins

 implicit none
 public

!this module is the state variable specification file for carbon parameters. it is the same as CA2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!schema_version: 2.0.0
!component: CA


 type CA2G_bc_State

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
!..................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: bc_biomass       => null() ! biomass burning emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_biofuel       => null() ! biofuel emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_antebc1       => null() ! anthropogenic bf emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_antebc2       => null() ! anthropogenic ff emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_ship          => null() ! ship emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_aviation_lto  => null() ! landing/take-off aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_aviation_cds  => null() ! climb/descent aircraft emissions (-)
 real(kind=RKIND),dimension(:,:),pointer  :: bc_aviation_crs  => null() ! cruise aircraft source species (-)
 real(kind=RKIND),dimension(:,:,:),pointer:: bc_aircraft      => null() ! aircraft emissions (-)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcmass         => null() ! black carbon aerosol mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcconc         => null() ! black carbon aerosol mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcextcoef      => null() ! black carbon aerosol extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcextcoefrh20  => null() ! black carbon aerosol extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcextcoefrh80  => null() ! black carbon aerosol extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcscacoef      => null() ! black carbon aerosol scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcscacoefrh20  => null() ! black carbon aerosol scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcscacoefrh80  => null() ! black carbon aerosol scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcbckcoef      => null() ! black carbon aerosol backscatter coefficient (m-1 sr-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcem           => null() ! black carbon aerosol emission (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcsd           => null() ! black carbon aerosol sedimentation (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcdp           => null() ! black carbon aerosol dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcwt           => null() ! black carbon aerosol wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcsv           => null() ! black carbon aerosol convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bceman         => null() ! black carbon aerosol anthropogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bcembb         => null() ! black carbon aerosol biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bcembf         => null() ! black carbon aerosol biofuel emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bcembg         => null() ! black carbon aerosol biogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bchyphil       => null() ! black carbon aerosol hydrophobic to hydrophilic (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bcsmass        => null() ! black carbon aerosol surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: bccmass        => null() ! black carbon aerosol column mass density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcexttau       => null() ! black carbon aerosol extinction aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcstexttau     => null() ! black carbon aerosol extinction aot stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcscatau       => null() ! black carbon aerosol scattering aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcstscatau     => null() ! black carbon aerosol scattering aot stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer    :: bcangstr       => null() ! black carbon aerosol angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:),pointer    :: bcfluxu        => null() ! black carbon aerosol column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bcfluxv        => null() ! black carbon aerosol column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: bcaeridx       => null() ! black carbon aerosol toms uv aerosol index (-)
 real(kind=RKIND),dimension(:,:),pointer    :: bcvdep         => null() ! dry deposition velocity (m s-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bctau_lw       => null() ! black carbon optical depth for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcasy_lw       => null() ! black carbon asymmetry factor for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcssa_lw       => null() ! black carbon single scattering albedo for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bctau_sw       => null() ! black carbon optical depth for shortwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcasy_sw       => null() ! black carbon asymmetry factor for shortwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: bcssa_sw       => null() ! black carbon single scattering albedo for shortwave RRTMG (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcphobic       => null() ! Hydrophobic black carbon aerosol mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: bcphilic       => null() ! Hydrophilic black carbon aerosol mixing Ratio (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => CA2G_bc_StateSpecsInit
    procedure:: gocart2G_deallocate => CA2G_bc_StateSpecsFinalize

 end type CA2G_bc_State


 contains


!==================================================================================================================
 subroutine CA2G_bc_StateSpecsInit(self,its,ite,jts,jte,kts,kte,nbndlw,nbndsw)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nbndlw,nbndsw

!--- inout arguments:
 class(CA2G_bc_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!------------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
!if(.not.associated(self%lats)           ) allocate(self%lats(its:ite,jts:jte)                   )
!if(.not.associated(self%lons)           ) allocate(self%lons(its:ite,jts:jte)                   )
!if(.not.associated(self%area)           ) allocate(self%area(its:ite,jts:jte)                   )
!if(.not.associated(self%frocean)        ) allocate(self%frocean(its:ite,jts:jte)                )
!if(.not.associated(self%fraci)          ) allocate(self%fraci(its:ite,jts:jte)                  )
!if(.not.associated(self%frlake)         ) allocate(self%frlake(its:ite,jts:jte)                 )
!if(.not.associated(self%lwi)            ) allocate(self%lwi(its:ite,jts:jte)                    )
!if(.not.associated(self%u10m)           ) allocate(self%u10m(its:ite,jts:jte)                   )
!if(.not.associated(self%v10m)           ) allocate(self%v10m(its:ite,jts:jte)                   )
!if(.not.associated(self%zpbl)           ) allocate(self%zpbl(its:ite,jts:jte)                   )
!if(.not.associated(self%ustar)          ) allocate(self%ustar(its:ite,jts:jte)                  )
!if(.not.associated(self%sh)             ) allocate(self%sh(its:ite,jts:jte)                     )
!if(.not.associated(self%z0h)            ) allocate(self%z0h(its:ite,jts:jte)                    )
!if(.not.associated(self%cn_prcp)        ) allocate(self%cn_prcp(its:ite,jts:jte)                )
!if(.not.associated(self%ncn_prcp)       ) allocate(self%ncn_prcp(its:ite,jts:jte)               )
!if(.not.associated(self%tropp)          ) allocate(self%tropp(its:ite,jts:jte)                  )
!........................................ .........................................................................
!if(.not.associated(self%airdens)        ) allocate(self%airdens(its:ite,jts:jte,kts:kte)        )
!if(.not.associated(self%delp)           ) allocate(self%delp(its:ite,jts:jte,kts:kte)           )
!if(.not.associated(self%delz)           ) allocate(self%delz(its:ite,jts:jte,kts:kte)           )
!if(.not.associated(self%t)              ) allocate(self%t(its:ite,jts:jte,kts:kte)              )
!if(.not.associated(self%rh2)            ) allocate(self%rh2(its:ite,jts:jte,kts:kte)            )
!if(.not.associated(self%pfl_lsan)       ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)       )
!if(.not.associated(self%pfi_lsan)       ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)       )
!if(.not.associated(self%u)              ) allocate(self%u(its:ite,jts:jte,kts:kte)              )
!if(.not.associated(self%v)              ) allocate(self%v(its:ite,jts:jte,kts:kte)              )
!if(.not.associated(self%zle)            ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)          )
!if(.not.associated(self%ple)            ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)          )
!..................................................................................................................
!if(.not.associated(self%bc_biomass)     ) allocate(self%bc_biomass(its:ite,jts:jte)             )
!if(.not.associated(self%bc_biofuel)     ) allocate(self%bc_biofuel(its:ite,jts:jte)             )
!if(.not.associated(self%bc_antebc1)     ) allocate(self%bc_antebc1(its:ite,jts:jte)             )
!if(.not.associated(self%bc_antebc2)     ) allocate(self%bc_antebc2(its:ite,jts:jte)             )
!if(.not.associated(self%bc_ship)        ) allocate(self%bc_ship(its:ite,jts:jte)                )
!if(.not.associated(self%bc_aviation_lto)) allocate(self%bc_aviation_lto(its:ite,jts:jte)        )
!if(.not.associated(self%bc_aviation_cds)) allocate(self%bc_aviation_cds(its:ite,jts:jte)        )
!if(.not.associated(self%bc_aviation_crs)) allocate(self%bc_aviation_crs(its:ite,jts:jte)        )
!if(.not.associated(self%bc_aircraft)    ) allocate(self%bc_aircraft(its:ite,jts:jte,kts:kte)    )

!category: EXPORT
 if(.not.associated(self%bcmass)         ) allocate(self%bcmass(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%bcconc)         ) allocate(self%bcconc(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%bcextcoef)      ) allocate(self%bcextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%bcextcoefrh20)  ) allocate(self%bcextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%bcextcoefrh80)  ) allocate(self%bcextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%bcscacoef)      ) allocate(self%bcscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%bcscacoefrh20)  ) allocate(self%bcscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%bcscacoefrh80)  ) allocate(self%bcscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%bcbckcoef)      ) allocate(self%bcbckcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
!..................................................................................................................
 if(.not.associated(self%bcem)           ) allocate(self%bcem(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcsd)           ) allocate(self%bcsd(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcdp)           ) allocate(self%bcdp(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcwt)           ) allocate(self%bcwt(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bcsv)           ) allocate(self%bcsv(its:ite,jts:jte,nbins)             )
 if(.not.associated(self%bceman)         ) allocate(self%bceman(its:ite,jts:jte)                 )
 if(.not.associated(self%bcembb)         ) allocate(self%bcembb(its:ite,jts:jte)                 )
 if(.not.associated(self%bcembf)         ) allocate(self%bcembf(its:ite,jts:jte)                 )
 if(.not.associated(self%bcembg)         ) allocate(self%bcembg(its:ite,jts:jte)                 )
 if(.not.associated(self%bchyphil)       ) allocate(self%bchyphil(its:ite,jts:jte)               )
 if(.not.associated(self%bcsmass)        ) allocate(self%bcsmass(its:ite,jts:jte)                )
 if(.not.associated(self%bccmass)        ) allocate(self%bccmass(its:ite,jts:jte)                )
 if(.not.associated(self%bcexttau)       ) allocate(self%bcexttau(its:ite,jts:jte,nw_vertint)    )
 if(.not.associated(self%bcstexttau)     ) allocate(self%bcstexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%bcscatau)       ) allocate(self%bcscatau(its:ite,jts:jte,nw_vertint)    )
 if(.not.associated(self%bcstscatau)     ) allocate(self%bcstscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%bcangstr)       ) allocate(self%bcangstr(its:ite,jts:jte)               )
 if(.not.associated(self%bcfluxu)        ) allocate(self%bcfluxu(its:ite,jts:jte)                )
 if(.not.associated(self%bcfluxv)        ) allocate(self%bcfluxv(its:ite,jts:jte)                )
 if(.not.associated(self%bcaeridx)       ) allocate(self%bcaeridx(its:ite,jts:jte)               )
 if(.not.associated(self%bcvdep)         ) allocate(self%bcvdep(its:ite,jts:jte)                 )
!..................................................................................................................
 if(.not.associated(self%bctau_lw)       ) allocate(self%bctau_lw(its:ite,jts:jte,kts:kte,nbndlw))
 if(.not.associated(self%bcssa_lw)       ) allocate(self%bcssa_lw(its:ite,jts:jte,kts:kte,nbndlw))
 if(.not.associated(self%bcasy_lw)       ) allocate(self%bcasy_lw(its:ite,jts:jte,kts:kte,nbndlw))
 if(.not.associated(self%bctau_sw)       ) allocate(self%bctau_sw(its:ite,jts:jte,kts:kte,nbndsw))
 if(.not.associated(self%bcssa_sw)       ) allocate(self%bcssa_sw(its:ite,jts:jte,kts:kte,nbndsw))
 if(.not.associated(self%bcasy_sw)       ) allocate(self%bcasy_sw(its:ite,jts:jte,kts:kte,nbndsw))

!category: INTERNAL
!if(.not.associated(self%bcphobic)       ) allocate(self%bcphobic(its:ite,jts:jte,kts:kte)       )
!if(.not.associated(self%bcphilic)       ) allocate(self%bcphilic(its:ite,jts:jte,kts:kte)       )


!--- initialization of diagnostics used in GOCART2G_GridComp:
 self%bcangstr(:,:) = 0._RKIND
 self%bcsmass(:,:)  = 0._RKIND

 self%bcexttau(:,:,:)   = 0._RKIND
 self%bcstexttau(:,:,:) = 0._RKIND
 self%bcscatau(:,:,:)   = 0._RKIND
 self%bcstscatau(:,:,:) = 0._RKIND

 self%bcextcoef(:,:,:,:)     = 0._RKIND
 self%bcextcoefrh20(:,:,:,:) = 0._RKIND
 self%bcextcoefrh80(:,:,:,:) = 0._RKIND
 self%bcscacoef(:,:,:,:)     = 0._RKIND
 self%bcscacoefrh20(:,:,:,:) = 0._RKIND
 self%bcscacoefrh80(:,:,:,:) = 0._RKIND
 self%bcbckcoef(:,:,:,:)     = 0._RKIND


!--- initialization of RRTMG longwave and shortwave optical properties:
 self%bctau_lw(:,:,:,:) = 0._RKIND
 self%bcssa_lw(:,:,:,:) = 0._RKIND
 self%bcasy_lw(:,:,:,:) = 0._RKIND
 self%bctau_sw(:,:,:,:) = 0._RKIND
 self%bcssa_sw(:,:,:,:) = 0._RKIND
 self%bcasy_sw(:,:,:,:) = 0._RKIND

 end subroutine CA2G_bc_StateSpecsInit

!==================================================================================================================
 subroutine CA2G_bc_StateSpecsFinalize(self)
!==================================================================================================================

!--- inout arguments:
 class(CA2G_bc_State),intent(inout) :: self

!------------------------------------------------------------------------------------------------------------------

!category: IMPORT
!if(associated(self%lats)           ) deallocate(self%lats           )
!if(associated(self%lons)           ) deallocate(self%lons           )
!if(associated(self%area)           ) deallocate(self%area           )
!if(associated(self%frocean)        ) deallocate(self%frocean        )
!if(associated(self%fraci)          ) deallocate(self%fraci          )
!if(associated(self%frlake)         ) deallocate(self%frlake         )
!if(associated(self%lwi)            ) deallocate(self%lwi            )
!if(associated(self%u10m)           ) deallocate(self%u10m           )
!if(associated(self%v10m)           ) deallocate(self%v10m           )
!if(associated(self%zpbl)           ) deallocate(self%zpbl           )
!if(associated(self%ustar)          ) deallocate(self%ustar          )
!if(associated(self%sh)             ) deallocate(self%sh             )
!if(associated(self%z0h)            ) deallocate(self%z0h            )
!if(associated(self%cn_prcp)        ) deallocate(self%cn_prcp        )
!if(associated(self%ncn_prcp)       ) deallocate(self%ncn_prcp       )
!if(associated(self%tropp)          ) deallocate(self%tropp          )
!........................................ .........................................................................
!if(associated(self%airdens)        ) deallocate(self%airdens        )
!if(associated(self%delp)           ) deallocate(self%delp           )
!if(associated(self%delz)           ) deallocate(self%delz           )
!if(associated(self%t)              ) deallocate(self%t              )
!if(associated(self%rh2)            ) deallocate(self%rh2            )
!if(associated(self%u)              ) deallocate(self%u              )
!if(associated(self%v)              ) deallocate(self%v              )
!if(associated(self%pfl_lsan)       ) deallocate(self%pfl_lsan       )
!if(associated(self%pfi_lsan)       ) deallocate(self%pfi_lsan       )
!if(associated(self%zle)            ) deallocate(self%zle            )
!if(associated(self%ple)            ) deallocate(self%ple            )
!..................................................................................................................
!if(associated(self%bc_biomass)     ) deallocate(self%bc_biomass     )
!if(associated(self%bc_biofuel)     ) deallocate(self%bc_biofuel     )
!if(associated(self%bc_antebc1)     ) deallocate(self%bc_antebc1     )
!if(associated(self%bc_antebc2)     ) deallocate(self%bc_antebc2     )
!if(associated(self%bc_ship)        ) deallocate(self%bc_ship        )
!if(associated(self%bc_aviation_lto)) deallocate(self%bc_aviation_lto)
!if(associated(self%bc_aviation_cds)) deallocate(self%bc_aviation_cds)
!if(associated(self%bc_aviation_crs)) deallocate(self%bc_aviation_crs)
!if(associated(self%bc_aircraft)    ) deallocate(self%bc_aircraft    )

!category: EXPORT
 if(associated(self%bcmass)         ) deallocate(self%bcmass         )
 if(associated(self%bcconc)         ) deallocate(self%bcconc         )
 if(associated(self%bcextcoef)      ) deallocate(self%bcextcoef      )
 if(associated(self%bcextcoefrh20)  ) deallocate(self%bcextcoefrh20  )
 if(associated(self%bcextcoefrh80)  ) deallocate(self%bcextcoefrh80  )
 if(associated(self%bcscacoef)      ) deallocate(self%bcscacoef      )
 if(associated(self%bcscacoefrh20)  ) deallocate(self%bcscacoefrh20  )
 if(associated(self%bcscacoefrh80)  ) deallocate(self%bcscacoefrh80  )
 if(associated(self%bcbckcoef)      ) deallocate(self%bcbckcoef      )
!..................................................................................................................
 if(associated(self%bcem)           ) deallocate(self%bcem           )
 if(associated(self%bcsd)           ) deallocate(self%bcsd           )
 if(associated(self%bcdp)           ) deallocate(self%bcdp           )
 if(associated(self%bcwt)           ) deallocate(self%bcwt           )
 if(associated(self%bcsv)           ) deallocate(self%bcsv           )
 if(associated(self%bceman)         ) deallocate(self%bceman         )
 if(associated(self%bcembb)         ) deallocate(self%bcembb         )
 if(associated(self%bcembf)         ) deallocate(self%bcembf         )
 if(associated(self%bcembg)         ) deallocate(self%bcembg         )
 if(associated(self%bchyphil)       ) deallocate(self%bchyphil       )
 if(associated(self%bcsmass)        ) deallocate(self%bcsmass        )
 if(associated(self%bccmass)        ) deallocate(self%bccmass        )
 if(associated(self%bcexttau)       ) deallocate(self%bcexttau       )
 if(associated(self%bcstexttau)     ) deallocate(self%bcstexttau     )
 if(associated(self%bcscatau)       ) deallocate(self%bcscatau       )
 if(associated(self%bcstscatau)     ) deallocate(self%bcstscatau     )
 if(associated(self%bcangstr)       ) deallocate(self%bcangstr       )
 if(associated(self%bcfluxu)        ) deallocate(self%bcfluxu        )
 if(associated(self%bcfluxv)        ) deallocate(self%bcfluxv        )
 if(associated(self%bcaeridx)       ) deallocate(self%bcaeridx       )
 if(associated(self%bcvdep)         ) deallocate(self%bcvdep         )
!..................................................................................................................
 if(associated(self%bctau_lw)       ) deallocate(self%bctau_lw       )
 if(associated(self%bcssa_lw)       ) deallocate(self%bcssa_lw       )
 if(associated(self%bcasy_lw)       ) deallocate(self%bcasy_lw       )
 if(associated(self%bctau_sw)       ) deallocate(self%bctau_sw       )
 if(associated(self%bcssa_sw)       ) deallocate(self%bcssa_sw       )
 if(associated(self%bcasy_sw)       ) deallocate(self%bcasy_sw       )

!category: INTERNAL
!if(associated(self%bcphobic)       ) deallocate(self%bcphobic       )
!if(associated(self%bcphilic)       ) deallocate(self%bcphilic       )

 end subroutine CA2G_bc_StateSpecsFinalize

!==================================================================================================================
 end module CA2G_bc_StateSpecs
!==================================================================================================================
