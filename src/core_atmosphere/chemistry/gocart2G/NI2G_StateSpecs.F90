! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module NI2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm

 implicit none
 public

!this module is the state variable specification file for nitrate parameters. it is the same as NI2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/NI2G_GridComp.

!schema_version: 2.0.0
!component: NI


 type NI2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer    :: lwi           => null() ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),pointer    :: tropp         => null() ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer    :: ustar         => null() ! surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: zpbl          => null() ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer    :: sh            => null() ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer    :: z0h           => null() ! surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer    :: cn_prcp       => null() ! surface_conv._rain_flux_needed_by_land (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: ncn_prcp      => null() ! non-convective precipitation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: area          => null() ! grid_cell_area (m^2)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer  :: airdens       => null() ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: delp          => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer  :: t             => null() ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer  :: rh2           => null() ! rel_hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: zle           => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer  :: ple           => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer  :: pfl_lsan      => null() ! 3d_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer  :: pfi_lsan      => null() ! 3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer  :: u             => null() ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: v             => null() ! northward_wind (m s-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_ag    => null() ! agriculture emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_bb    => null() ! biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_en    => null() ! energy emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_in    => null() ! industry emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_oc    => null() ! ocean emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_re    => null() ! residential emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: emi_nh3_tr    => null() ! transport emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nitrate_hno3  => null() ! nitrate hno3 emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: du            => null() ! dust mixing ratio all bins (kg kg-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: ss            => null() ! sea salt mixing ratio all bins (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: so4           => null() ! sulfate mixing ratio (kg kg-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer  :: nh3mass       => null() ! ammonia mass mixing ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nh4mass       => null() ! ammonium aerosol mass mixing ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nimass        => null() ! nitrate mass mixing ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nimass25      => null() ! nitrate mass mixing ratio [pm2.5] (kg/kg)
 real(kind=RKIND),dimension(:,:,:),pointer  :: hno3conc      => null() ! nitric acid mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nh3conc       => null() ! ammonia mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nh4conc       => null() ! ammonium mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niconc        => null() ! nitrate mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niconc25      => null() ! nitrate mass concentration [pm2.5] (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: niextcoef     => null() ! nitrate extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: niextcoefrh20 => null() ! nitrate extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: niextcoefrh80 => null() ! nitrate extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: niscacoef     => null() ! nitrate scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: niscacoefrh20 => null() ! nitrate scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: niscacoefrh80 => null() ! nitrate scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: nibckcoef     => null() ! nitrate backscatter coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer    :: nipno3aq      => null() ! nitrate production from aqueous chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nipnh4aq      => null() ! ammonium production from aqueous chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nipnh3aq      => null() ! ammonia change from aqueous chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niht          => null() ! nitrate production from het chem (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nisd          => null() ! nitrate sedimentation (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nidp          => null() ! nitrate dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niwt          => null() ! nitrate wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nisv          => null() ! nitrate convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh3em         => null() ! ammonia emission (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh3dp         => null() ! ammonia dry deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh3wt         => null() ! ammonia wet deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh3sv         => null() ! ammonia convective scavenging (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh4sd         => null() ! ammonium settling (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh4dp         => null() ! ammonium dry deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh4wt         => null() ! ammonium wet deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nh4sv         => null() ! ammonium convective scavenging (kg m-2 s-1)

 real(kind=RKIND),dimension(:,:),pointer    :: hno3smass     => null() ! nitric acid surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nh3smass      => null() ! ammonia surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nh4smass      => null() ! ammonium surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nismass       => null() ! nitrate surface mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nismass25     => null() ! nitrate surface mass concentration [pm2.5] (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: hno3cmass     => null() ! nitric acid column mass density (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nh3cmass      => null() ! ammonia column mass density (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nh4cmass      => null() ! ammonium column mass density (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer    :: nicmass       => null() ! nitrate column mass density (kg m-2)
 real(kind=RKIND),dimension(:,:),pointer    :: nicmass25     => null() ! nitrate column mass density [pm2.5] (kg m-2)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niexttfm      => null() ! nitrate extinction aot - pm 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niscatfm      => null() ! nitrate scattering aot - pm 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niextt25      => null() ! nitrate extinction aot - pm 2.5 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niscat25      => null() ! nitrate scattering aot - pm 2.5 um (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niexttau      => null() ! nitrate extinction aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nistexttau    => null() ! nitrate extinction aot stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: niscatau      => null() ! nitrate scattering aot (-)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nistscatau    => null() ! nitrate scattering aot stratosphere (-)
 real(kind=RKIND),dimension(:,:),pointer    :: niangstr      => null() ! nitrate angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:),pointer    :: nifluxu       => null() ! nitrate column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer    :: nifluxv       => null() ! nitrate column v-wind mass flux (kg m-1 s-1)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer  :: nh3           => null() ! ammonia (nh3, gas phase) (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: nh4a          => null() ! ammonium ion (nh4+, aerosol phase) (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: no3an1        => null() ! nitrate size bin 001 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: no3an2        => null() ! nitrate size bin 002 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: no3an3        => null() ! nitrate size bin 003 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer  :: xhno3         => null() ! buffer for nitrate_hno3 (kg m-2 s-1)


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

!--- local variables:
 integer:: nw_profile,nw_vertint

!-----------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
!if(.not.associated(self%lwi)         ) allocate(self%lwi(its:ite,jts:jte)                 )
!if(.not.associated(self%tropp)       ) allocate(self%tropp(its:ite,jts:jte)               )
!if(.not.associated(self%ustar)       ) allocate(self%tropp(its:ite,jts:jte)               )
!if(.not.associated(self%zpbl)        ) allocate(self%zpbl(its:ite,jts:jte)                )
!if(.not.associated(self%sh)          ) allocate(self%zpbl(its:ite,jts:jte)                )
!if(.not.associated(self%z0h)         ) allocate(self%zpbl(its:ite,jts:jte)                )
!if(.not.associated(self%cn_prcp)     ) allocate(self%cn_prcp(its:ite,jts:jte)             )
!if(.not.associated(self%ncn_prcp)    ) allocate(self%ncn_prcp(its:ite,jts:jte)            )
 if(.not.associated(self%area)        ) allocate(self%area(its:ite,jts:jte)                )
!.................................................................................................................
!if(.not.associated(self%airdens)     ) allocate(self%airdens(its:ite,jts:jte,kts:kte)     )
!if(.not.associated(self%delp)        ) allocate(self%delp(its:ite,jts:jte,kts:kte)        )
!if(.not.associated(self%t)           ) allocate(self%t(its:ite,jts:jte,kts:kte)           )
!if(.not.associated(self%rh2)         ) allocate(self%rh2(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%zle)         ) allocate(self%zle(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%ple)         ) allocate(self%ple(its:ite,jts:jte,kts:kte)         )
!if(.not.associated(self%pfl_lsan)    ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)    )
!if(.not.associated(self%pfi_lsan)    ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)    )
!if(.not.associated(self%u)           ) allocate(self%u(its:ite,jts:jte,kts:kte)           )
!if(.not.associated(self%v)           ) allocate(self%v(its:ite,jts:jte,kts:kte)           )
!.................................................................................................................
!if(.not.associated(self%emi_nh3_ag)  ) allocate(self%emi_nh3_ag(its:ite,jts:jte)          )
!if(.not.associated(self%emi_nh3_bb)  ) allocate(self%emi_nh3_bb(its:ite,jts:jte)          )
!if(.not.associated(self%emi_nh3_en)  ) allocate(self%emi_nh3_en(its:ite,jts:jte)          )
!if(.not.associated(self%emi_nh3_in)  ) allocate(self%emi_nh3_in(its:ite,jts:jte)          )
!if(.not.associated(self%emi_nh3_oc)  ) allocate(self%emi_nh3_oc(its:ite,jts:jte)          )
!if(.not.associated(self%emi_nh3_re)  ) allocate(self%emi_nh3_re(its:ite,jts:jte)          )
!if(.not.associated(self%emi_nh3_tr)  ) allocate(self%emi_nh3_tr(its:ite,jts:jte)          )
 if(.not.associated(self%nitrate_hno3)) allocate(self%nitrate_hno3(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%du)          ) allocate(self%du(its:ite,jts:jte,kts:kte,5)        )
 if(.not.associated(self%ss)          ) allocate(self%ss(its:ite,jts:jte,kts:kte,5)        )
 if(.not.associated(self%so4)         ) allocate(self%so4(its:ite,jts:jte,kts:kte)         )

!category: EXPORT
 if(.not.associated(self%nh3mass)      ) allocate(self%nh3mass(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%nh4mass)      ) allocate(self%nh4mass(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%nimass)       ) allocate(self%nimass(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%nimass25)     ) allocate(self%nimass25(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%hno3conc)     ) allocate(self%hno3conc(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%nh3conc)      ) allocate(self%nh3conc(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%nh4conc)      ) allocate(self%nh4conc(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%niconc)       ) allocate(self%niconc(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%niconc25)     ) allocate(self%niconc25(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%niextcoef)    ) allocate(self%niextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%niextcoefrh20)) allocate(self%niextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%niextcoefrh80)) allocate(self%niextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%niscacoef)    ) allocate(self%niscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%niscacoefrh20)) allocate(self%niscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%niscacoefrh80)) allocate(self%niscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%nibckcoef)    ) allocate(self%nibckcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
!.................................................................................................................
 if(.not.associated(self%nipno3aq)     ) allocate(self%nipno3aq(its:ite,jts:jte)        )
 if(.not.associated(self%nipnh4aq)     ) allocate(self%nipnh4aq(its:ite,jts:jte)        )
 if(.not.associated(self%nipnh3aq)     ) allocate(self%nipnh3aq(its:ite,jts:jte)        )
 if(.not.associated(self%niht)         ) allocate(self%niht(its:ite,jts:jte,3)          )
 if(.not.associated(self%nisd)         ) allocate(self%nisd(its:ite,jts:jte,3)          )
 if(.not.associated(self%nidp)         ) allocate(self%nidp(its:ite,jts:jte,3)          )
 if(.not.associated(self%niwt)         ) allocate(self%niwt(its:ite,jts:jte,3)          )
 if(.not.associated(self%nisv)         ) allocate(self%nisv(its:ite,jts:jte,3)          )
 if(.not.associated(self%nh3em)        ) allocate(self%nh3em(its:ite,jts:jte)           )
 if(.not.associated(self%nh3dp)        ) allocate(self%nh3dp(its:ite,jts:jte)           )
 if(.not.associated(self%nh3wt)        ) allocate(self%nh3wt(its:ite,jts:jte)           )
 if(.not.associated(self%nh3sv)        ) allocate(self%nh3sv(its:ite,jts:jte)           )
 if(.not.associated(self%nh4sd)        ) allocate(self%nh4sd(its:ite,jts:jte)           )
 if(.not.associated(self%nh4dp)        ) allocate(self%nh4dp(its:ite,jts:jte)           )
 if(.not.associated(self%nh4wt)        ) allocate(self%nh4wt(its:ite,jts:jte)           )
 if(.not.associated(self%nh4sv)        ) allocate(self%nh4sv(its:ite,jts:jte)           ) 
 if(.not.associated(self%hno3smass)    ) allocate(self%hno3smass(its:ite,jts:jte)       )
 if(.not.associated(self%nh3smass)     ) allocate(self%nh3smass(its:ite,jts:jte)        )
 if(.not.associated(self%nh4smass)     ) allocate(self%nh4smass(its:ite,jts:jte)        )
 if(.not.associated(self%nismass)      ) allocate(self%nismass(its:ite,jts:jte)         )
 if(.not.associated(self%nismass25)    ) allocate(self%nismass25(its:ite,jts:jte)       ) 
 if(.not.associated(self%hno3cmass)    ) allocate(self%hno3cmass(its:ite,jts:jte)       )
 if(.not.associated(self%nh3cmass)     ) allocate(self%nh3cmass(its:ite,jts:jte)        )
 if(.not.associated(self%nh4cmass)     ) allocate(self%nh4cmass(its:ite,jts:jte)        )
 if(.not.associated(self%nicmass)      ) allocate(self%nicmass(its:ite,jts:jte)         )
 if(.not.associated(self%nicmass25)    ) allocate(self%nicmass25(its:ite,jts:jte)       )
 if(.not.associated(self%niexttfm)     ) allocate(self%niexttfm(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%niscatfm)     ) allocate(self%niscatfm(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%niextt25)     ) allocate(self%niextt25(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%niscat25)     ) allocate(self%niscat25(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%niexttau)     ) allocate(self%niexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%nistexttau)   ) allocate(self%nistexttau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%niscatau)     ) allocate(self%niscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%nistscatau)   ) allocate(self%nistscatau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%niangstr)     ) allocate(self%niangstr(its:ite,jts:jte)        )
 if(.not.associated(self%nifluxu)      ) allocate(self%nifluxu(its:ite,jts:jte)         )
 if(.not.associated(self%nifluxv)      ) allocate(self%nifluxv(its:ite,jts:jte)         )

!category: INTERNAL
!if(.not.associated(self%nh3)   ) allocate(self%nh3(its:ite,jts:jte,kts:kte)   )
!if(.not.associated(self%nh4a)  ) allocate(self%nh4a(its:ite,jts:jte,kts:kte)  )
!if(.not.associated(self%no3an1)) allocate(self%no3an1(its:ite,jts:jte,kts:kte))
!if(.not.associated(self%no3an2)) allocate(self%no3an2(its:ite,jts:jte,kts:kte))
!if(.not.associated(self%no3an3)) allocate(self%no3an3(its:ite,jts:jte,kts:kte))
!if(.not.associated(self%xhno3) ) allocate(self%xhno3(its:ite,jts:jte,kts:kte) )

 end subroutine NI2G_StateSpecsInit

!=================================================================================================================
 subroutine NI2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(NI2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
!if(associated(self%lwi)         ) deallocate(self%lwi         )
!if(associated(self%tropp)       ) deallocate(self%tropp       )
!if(associated(self%ustar)       ) deallocate(self%tropp       )
!if(associated(self%zpbl)        ) deallocate(self%zpbl        )
!if(associated(self%sh)          ) deallocate(self%zpbl        )
!if(associated(self%z0h)         ) deallocate(self%z0h         )
!if(associated(self%cn_prcp)     ) deallocate(self%cn_prcp     )
!if(associated(self%ncn_prcp)    ) deallocate(self%ncn_prcp    )
 if(associated(self%area)        ) deallocate(self%area        )
!.................................................................................................................
!if(associated(self%airdens)     ) deallocate(self%airdens     )
!if(associated(self%delp)        ) deallocate(self%delp        )
!if(associated(self%t)           ) deallocate(self%t           )
!if(associated(self%rh2)         ) deallocate(self%rh2         )
!if(associated(self%zle)         ) deallocate(self%zle         )
!if(associated(self%ple)         ) deallocate(self%ple         )
!if(associated(self%pfl_lsan)    ) deallocate(self%pfl_lsan    )
!if(associated(self%pfi_lsan)    ) deallocate(self%pfi_lsan    )
!if(associated(self%u)           ) deallocate(self%u           )
!if(associated(self%v)           ) deallocate(self%v           )
!.................................................................................................................
!if(associated(self%emi_nh3_ag)  ) deallocate(self%emi_nh3_ag  )
!if(associated(self%emi_nh3_bb)  ) deallocate(self%emi_nh3_bb  )
!if(associated(self%emi_nh3_en)  ) deallocate(self%emi_nh3_en  )
!if(associated(self%emi_nh3_in)  ) deallocate(self%emi_nh3_in  )
!if(associated(self%emi_nh3_oc)  ) deallocate(self%emi_nh3_oc  )
!if(associated(self%emi_nh3_re)  ) deallocate(self%emi_nh3_re  )
!if(associated(self%emi_nh3_tr)  ) deallocate(self%emi_nh3_tr  )
 if(associated(self%nitrate_hno3)) deallocate(self%nitrate_hno3)
 if(associated(self%du)          ) deallocate(self%du          )
 if(associated(self%ss)          ) deallocate(self%ss          )
 if(associated(self%so4)         ) deallocate(self%so4         )

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
 if(associated(self%niextcoef)    ) deallocate(self%niextcoef    )
 if(associated(self%niextcoefrh20)) deallocate(self%niextcoefrh20)
 if(associated(self%niextcoefrh80)) deallocate(self%niextcoefrh80)
 if(associated(self%niscacoef)    ) deallocate(self%niscacoef    )
 if(associated(self%niscacoefrh20)) deallocate(self%niscacoefrh20)
 if(associated(self%niscacoefrh80)) deallocate(self%niscacoefrh80)
 if(associated(self%nibckcoef)    ) deallocate(self%nibckcoef    )
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
 if(associated(self%nh4wt)        ) deallocate(self%nh4wt        )
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
 if(associated(self%niexttfm)     ) deallocate(self%niexttfm     )
 if(associated(self%niscatfm)     ) deallocate(self%niscatfm     )
 if(associated(self%niextt25)     ) deallocate(self%niextt25     )
 if(associated(self%niscat25)     ) deallocate(self%niscat25     )
 if(associated(self%niexttau)     ) deallocate(self%niexttau     )
 if(associated(self%nistexttau)   ) deallocate(self%nistexttau   )
 if(associated(self%niscatau)     ) deallocate(self%niscatau     )
 if(associated(self%nistscatau)   ) deallocate(self%nistscatau   )
 if(associated(self%niangstr)     ) deallocate(self%niangstr     )
 if(associated(self%nifluxu)      ) deallocate(self%nifluxu      )
 if(associated(self%nifluxv)      ) deallocate(self%nifluxv      )

!category: INTERNAL
!if(associated(self%nh3)   ) deallocate(self%nh3   )
!if(associated(self%nh4a)  ) deallocate(self%nh4a  )
!if(associated(self%no3an1)) deallocate(self%no3an1)
!if(associated(self%no3an2)) deallocate(self%no3an2)
!if(associated(self%no3an3)) deallocate(self%no3an3)
!if(associated(self%xhno3) ) deallocate(self%xhno3 )

 end subroutine NI2G_StateSpecsFinalize

!=================================================================================================================
 end module NI2G_StateSpecs
!=================================================================================================================
