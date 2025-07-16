! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module SS2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use SS2G_instance,only: nbins

 implicit none
 public
 save

!this module is the state variable specification file for sea-salt parameters. it is the same as SS2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SG2G_GridComp.

!schema_version: 2.0.0
!component: SS


 type SS2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer  :: frocean           => null() ! fraction_of_ocean (1)
 real(kind=RKIND),dimension(:,:),pointer  :: fraci             => null() ! ice_covered_fraction_of_tile (1)
 real(kind=RKIND),dimension(:,:),pointer  :: lwi               => null() ! land-ocean-ice_mask (1)
 real(kind=RKIND),dimension(:,:),pointer  :: tropp             => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer  :: u10m              => null() ! !10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: v10m              => null() ! 10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ustar             => null() ! surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: ts                => null() ! surface skin temperature (K)
 real(kind=RKIND),dimension(:,:),pointer  :: dz                => null() ! surface_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer  :: frlake            => null() ! fraction_of_lake (1)
 real(kind=RKIND),dimension(:,:),pointer  :: area              => null() ! grid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer  :: zpbl              => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer  :: sh                => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer  :: z0h               => null() ! surface_roughness_for_heat(m)
 real(kind=RKIND),dimension(:,:),pointer  :: cn_prcp           => null() ! surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer  :: ncn_prcp          => null() ! non-convective precipitation (kg/m^2/s)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: airdens           => null() ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer:: delp              => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: delz              => null() ! geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: t                 => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: rh2               => null() ! rel_hum_after_moist (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: zle               => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: ple               => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfl_lsan          => null() ! 3d_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfi_lsan          => null() ! 3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: u                 => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: v                 => null() ! northward_wind (m s-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssmass          => null() ! sea salt mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssmass25        => null() ! sea salt mass mixing ratio - pm 2.5 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssconc          => null() ! sea salt mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssextcoef       => null() ! sea salt extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssextcoefrh20   => null() ! sea salt extinction coefficient - fixed Rh=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssextcoefrh80   => null() ! sea salt extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssscacoef       => null() ! sea salt scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssscacoefrh20   => null() ! sea salt scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssscacoefrh80   => null() ! sea salt scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssbckcoef       => null() ! sea salt backscatter coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssem            => null() ! sea salt emission (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: sssd            => null() ! sea salt sedimentation (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssdp            => null() ! sea salt dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: sswt            => null() ! sea salt wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: sssv            => null() ! sea salt convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: sssmass         => null() ! sea salt surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: sscmass         => null() ! sea salt column mass density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssexttau        => null() ! sea salt extinction aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssstexttau      => null() ! sea salt extinction aot stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssscatau        => null() ! sea salt scattering aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssstscatau      => null() ! sea salt scattering aot stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer    :: sssmass25       => null() ! sea salt surface mass concentration - pm 2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: sscmass25       => null() ! sea salt column mass density - pm 2.5 (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssextt25        => null() ! extinction aot - pm 2.5 (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssscat25        => null() ! sea salt scattering aot - pm 2.5 (-)
 real(kind=RKIND),dimension(:,:),pointer    :: ssaeridx        => null() ! sea salt toms uv aerosol index (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssexttfm        => null() ! sea salt extinction aot [550 nm] - pm 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ssscatfm        => null() ! sea salt scattering aot [550 nm] - pm 1.0 um (-)
 real(kind=RKIND),dimension(:,:),pointer    :: ssangstr        => null() ! sea salt angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:),pointer    :: ssfluxu         => null() ! sea salt column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ssfluxv         => null() ! sea salt column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ssvdep          => null() ! dry deposition velocity (m s-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:,:),pointer:: sstau_lw        => null() ! sea salt optical depth for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssasy_lw        => null() ! sea salt asymmetry factor for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssssa_lw        => null() ! sea salt single scattering albedo for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: sstau_sw        => null() ! sea salt optical depth for shortwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssasy_sw        => null() ! sea salt asymmetry factor for shortwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ssssa_sw        => null() ! sea salt single scattering albedo for shortwave RRTMG (-)


!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ss              => null() ! sea salt mixing ratio (bin %d) (kg kg-1)
 real(kind=RKIND),dimension(:,:),pointer    :: deep_lakes_mask => null() ! deep lakes mask


 contains
    procedure:: gocart2G_allocate   => SS2G_StateSpecsInit
    procedure:: gocart2G_deallocate => SS2G_StateSpecsFinalize

 end type SS2G_State


 contains


!=================================================================================================================
 subroutine SS2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte,nbndlw,nbndsw)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nbndlw,nbndsw

!--- inout arguments:
 class(SS2G_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!-----------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
!if(.not.associated(self%frocean)        ) allocate(self%frocean(its:ite,jts:jte)              )
!if(.not.associated(self%fraci)          ) allocate(self%fraci(its:ite,jts:jte)                )
!if(.not.associated(self%lwi)            ) allocate(self%lwi(its:ite,jts:jte)                  )
!if(.not.associated(self%tropp)          ) allocate(self%tropp(its:ite,jts:jte)                )
!if(.not.associated(self%u10m)           ) allocate(self%u10m(its:ite,jts:jte)                 )
!if(.not.associated(self%v10m)           ) allocate(self%v10m(its:ite,jts:jte)                 )
!if(.not.associated(self%ustar)          ) allocate(self%ustar(its:ite,jts:jte)                )
!if(.not.associated(self%ts)             ) allocate(self%ts(its:ite,jts:jte)                   )
!if(.not.associated(self%dz)             ) allocate(self%dz(its:ite,jts:jte)                   )
!if(.not.associated(self%frlake)         ) allocate(self%frlake(its:ite,jts:jte)               )
!if(.not.associated(self%area)           ) allocate(self%area(its:ite,jts:jte)                 )
!if(.not.associated(self%zpbl)           ) allocate(self%zpbl(its:ite,jts:jte)                 )
!if(.not.associated(self%sh)             ) allocate(self%sh(its:ite,jts:jte)                   )
!if(.not.associated(self%z0h)            ) allocate(self%z0h(its:ite,jts:jte)                  )
!if(.not.associated(self%cn_prcp)        ) allocate(self%cn_prcp(its:ite,jts:jte)              )
!if(.not.associated(self%ncn_prcp)       ) allocate(self%ncn_prcp(its:ite,jts:jte)             )
!.................................................................................................................
!if(.not.associated(self%airdens)        ) allocate(self%airdens(its:ite,jts:jte,kts:kte)      )
!if(.not.associated(self%delp)           ) allocate(self%delp(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%delz)           ) allocate(self%delz(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%t)              ) allocate(self%t(its:ite,jts:jte,kts:kte)            )
!if(.not.associated(self%rh2)            ) allocate(self%rh2(its:ite,jts:jte,kts:kte)          )
!if(.not.associated(self%zle)            ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)        )
!if(.not.associated(self%ple)            ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)        )
!if(.not.associated(self%pfl_lsan)       ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)     )
!if(.not.associated(self%pfi_lsan)       ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)     )
!if(.not.associated(self%u)              ) allocate(self%u(its:ite,jts:jte,kts:kte)            )
!if(.not.associated(self%v)              ) allocate(self%v(its:ite,jts:jte,kts:kte)            )

!category: EXPORT
 if(.not.associated(self%ssmass)         ) allocate(self%ssmass(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%ssmass25)       ) allocate(self%ssmass25(its:ite,jts:jte,kts:kte)     )
 if(.not.associated(self%ssconc)         ) allocate(self%ssconc(its:ite,jts:jte,kts:kte)       )
 if(.not.associated(self%ssextcoef)      ) allocate(self%ssextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%ssextcoefrh20)  ) allocate(self%ssextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ssextcoefrh80)  ) allocate(self%ssextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ssscacoef)      ) allocate(self%ssscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%ssscacoefrh20)  ) allocate(self%ssscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ssscacoefrh80)  ) allocate(self%ssscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%ssbckcoef)      ) allocate(self%ssbckcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
!.................................................................................................................
 if(.not.associated(self%ssem)           ) allocate(self%ssem(its:ite,jts:jte,nbins)           )
 if(.not.associated(self%sssd)           ) allocate(self%sssd(its:ite,jts:jte,nbins)           )
 if(.not.associated(self%ssdp)           ) allocate(self%ssdp(its:ite,jts:jte,nbins)           )
 if(.not.associated(self%sswt)           ) allocate(self%sswt(its:ite,jts:jte,nbins)           )
 if(.not.associated(self%sssv)           ) allocate(self%sssv(its:ite,jts:jte,nbins)           )
 if(.not.associated(self%sssmass)        ) allocate(self%sssmass(its:ite,jts:jte)              )
 if(.not.associated(self%sscmass)        ) allocate(self%sscmass(its:ite,jts:jte)              )
 if(.not.associated(self%ssexttau)       ) allocate(self%ssexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ssstexttau)     ) allocate(self%ssstexttau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%ssscatau)       ) allocate(self%ssscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ssstscatau)     ) allocate(self%ssstscatau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%sssmass25)      ) allocate(self%sssmass25(its:ite,jts:jte)            )
 if(.not.associated(self%sscmass25)      ) allocate(self%sscmass25(its:ite,jts:jte)            )
 if(.not.associated(self%ssextt25)       ) allocate(self%ssextt25(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ssscat25)       ) allocate(self%ssscat25(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ssaeridx)       ) allocate(self%ssaeridx(its:ite,jts:jte)             )
 if(.not.associated(self%ssexttfm)       ) allocate(self%ssexttfm(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ssscatfm)       ) allocate(self%ssscatfm(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%ssangstr)       ) allocate(self%ssangstr(its:ite,jts:jte)             )
 if(.not.associated(self%ssfluxu)        ) allocate(self%ssfluxu(its:ite,jts:jte)              )
 if(.not.associated(self%ssfluxv)        ) allocate(self%ssfluxv(its:ite,jts:jte)              )
 if(.not.associated(self%ssvdep)         ) allocate(self%ssvdep(its:ite,jts:jte)               )
!..................................................................................................................
 if(.not.associated(self%sstau_lw)       ) allocate(self%sstau_lw(its:ite,jts:jte,kts:kte,nbndlw))
 if(.not.associated(self%ssssa_lw)       ) allocate(self%ssssa_lw(its:ite,jts:jte,kts:kte,nbndlw))
 if(.not.associated(self%ssasy_lw)       ) allocate(self%ssasy_lw(its:ite,jts:jte,kts:kte,nbndlw))
 if(.not.associated(self%sstau_sw)       ) allocate(self%sstau_sw(its:ite,jts:jte,kts:kte,nbndsw))
 if(.not.associated(self%ssssa_sw)       ) allocate(self%ssssa_sw(its:ite,jts:jte,kts:kte,nbndsw))
 if(.not.associated(self%ssasy_sw)       ) allocate(self%ssasy_sw(its:ite,jts:jte,kts:kte,nbndsw))

!category: INTERNAL
 if(.not.associated(self%ss)             ) allocate(self%ss(its:ite,jts:jte,kts:kte,nbins)     )
!if(.not.associated(self%deep_lakes_mask)) allocate(self%deep_lakes_mask(its:ite,jts:jte)      )


!--- initialization of diagnostics used in GOCART2G_GridComp:
 self%ssangstr(:,:) = 0._RKIND
 self%sssmass(:,:)  = 0._RKIND

 self%ssexttau(:,:,:)   = 0._RKIND
 self%ssstexttau(:,:,:) = 0._RKIND
 self%ssscatau(:,:,:)   = 0._RKIND
 self%ssstscatau(:,:,:) = 0._RKIND

 self%ssextcoef(:,:,:,:)     = 0._RKIND
 self%ssextcoefrh20(:,:,:,:) = 0._RKIND
 self%ssextcoefrh80(:,:,:,:) = 0._RKIND
 self%ssscacoef(:,:,:,:)     = 0._RKIND
 self%ssscacoefrh20(:,:,:,:) = 0._RKIND
 self%ssscacoefrh80(:,:,:,:) = 0._RKIND
 self%ssbckcoef(:,:,:,:)     = 0._RKIND


!--- initialization of RRTMG longwave and shortwave optical properties:
 self%sstau_lw(:,:,:,:) = 0._RKIND
 self%ssssa_lw(:,:,:,:) = 0._RKIND
 self%ssasy_lw(:,:,:,:) = 0._RKIND
 self%sstau_sw(:,:,:,:) = 0._RKIND
 self%ssssa_sw(:,:,:,:) = 0._RKIND
 self%ssasy_sw(:,:,:,:) = 0._RKIND


 end subroutine SS2G_StateSpecsInit

!=================================================================================================================
 subroutine SS2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(SS2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
!if(associated(self%frocean)        ) deallocate(self%frocean      )
!if(associated(self%fraci)          ) deallocate(self%fraci        )
!if(associated(self%lwi)            ) deallocate(self%lwi          )
!if(associated(self%tropp)          ) deallocate(self%tropp        )
!if(associated(self%u10m)           ) deallocate(self%u10m         )
!if(associated(self%v10m)           ) deallocate(self%v10m         )
!if(associated(self%ustar)          ) deallocate(self%ustar        )
!if(associated(self%ts)             ) deallocate(self%ts           )
!if(associated(self%dz)             ) deallocate(self%dz           )
!if(associated(self%frlake)         ) deallocate(self%frlake       )
!if(associated(self%area)           ) deallocate(self%area         )
!if(associated(self%zpbl)           ) deallocate(self%zpbl         )
!if(associated(self%sh)             ) deallocate(self%sh           )
!if(associated(self%z0h)            ) deallocate(self%z0h          )
!if(associated(self%cn_prcp)        ) deallocate(self%cn_prcp      )
!if(associated(self%ncn_prcp)       ) deallocate(self%ncn_prcp     )
!.................................................................................................................
!if(associated(self%airdens)        ) deallocate(self%airdens      )
!if(associated(self%delp)           ) deallocate(self%delp         )
!if(associated(self%delz)           ) deallocate(self%delz         )
!if(associated(self%t)              ) deallocate(self%t            )
!if(associated(self%rh2)            ) deallocate(self%rh2          )
!if(associated(self%zle)            ) deallocate(self%zle          )
!if(associated(self%ple)            ) deallocate(self%ple          )
!if(associated(self%pfl_lsan)       ) deallocate(self%pfl_lsan     )
!if(associated(self%pfi_lsan)       ) deallocate(self%pfi_lsan     )
!if(associated(self%u)              ) deallocate(self%u            )
!if(associated(self%v)              ) deallocate(self%v            )

!category: EXPORT
 if(associated(self%ssmass)         ) deallocate(self%ssmass       )
 if(associated(self%ssmass25)       ) deallocate(self%ssmass25     )
 if(associated(self%ssconc)         ) deallocate(self%ssconc       )
 if(associated(self%ssextcoef)      ) deallocate(self%ssextcoef    )
 if(associated(self%ssextcoefrh20)  ) deallocate(self%ssextcoefrh20)
 if(associated(self%ssextcoefrh80)  ) deallocate(self%ssextcoefrh80)
 if(associated(self%ssscacoef)      ) deallocate(self%ssscacoef    )
 if(associated(self%ssscacoefrh20)  ) deallocate(self%ssscacoefrh20)
 if(associated(self%ssscacoefrh80)  ) deallocate(self%ssscacoefrh80)
 if(associated(self%ssbckcoef)      ) deallocate(self%ssbckcoef    )
!.................................................................................................................
 if(associated(self%ssem)           ) deallocate(self%ssem         )
 if(associated(self%sssd)           ) deallocate(self%sssd         )
 if(associated(self%ssdp)           ) deallocate(self%ssdp         )
 if(associated(self%sswt)           ) deallocate(self%sswt         )
 if(associated(self%sssv)           ) deallocate(self%sssv         )
 if(associated(self%sssmass)        ) deallocate(self%sssmass      )
 if(associated(self%sscmass)        ) deallocate(self%sscmass      )
 if(associated(self%ssexttau)       ) deallocate(self%ssexttau     )
 if(associated(self%ssstexttau)     ) deallocate(self%ssstexttau   )
 if(associated(self%ssscatau)       ) deallocate(self%ssscatau     )
 if(associated(self%ssstscatau)     ) deallocate(self%ssstscatau   )
 if(associated(self%sssmass25)      ) deallocate(self%sssmass25    )
 if(associated(self%sscmass25)      ) deallocate(self%sscmass25    )
 if(associated(self%ssextt25)       ) deallocate(self%ssextt25     )
 if(associated(self%ssscat25)       ) deallocate(self%ssscat25     )
 if(associated(self%ssaeridx)       ) deallocate(self%ssaeridx     )

 if(associated(self%ssexttfm)       ) deallocate(self%ssexttfm     )
 if(associated(self%ssscatfm)       ) deallocate(self%ssscatfm     )
 if(associated(self%ssangstr)       ) deallocate(self%ssangstr     )
 if(associated(self%ssfluxu)        ) deallocate(self%ssfluxu      )
 if(associated(self%ssfluxv)        ) deallocate(self%ssfluxv      )
 if(associated(self%ssvdep)         ) deallocate(self%ssvdep       )
!.................................................................................................................
 if(associated(self%sstau_lw)       ) deallocate(self%sstau_lw     )
 if(associated(self%ssssa_lw)       ) deallocate(self%ssssa_lw     )
 if(associated(self%ssasy_lw)       ) deallocate(self%ssasy_lw     )
 if(associated(self%sstau_sw)       ) deallocate(self%sstau_sw     )
 if(associated(self%ssssa_sw)       ) deallocate(self%ssssa_sw     )
 if(associated(self%ssasy_sw)       ) deallocate(self%ssasy_sw     )

!category: INTERNAL
 if(associated(self%ss)             ) deallocate(self%ss             )
!if(associated(self%deep_lakes_mask)) deallocate(self%deep_lakes_mask)

 end subroutine SS2G_StateSpecsFinalize

!=================================================================================================================
 end module SS2G_StateSpecs
!=================================================================================================================
