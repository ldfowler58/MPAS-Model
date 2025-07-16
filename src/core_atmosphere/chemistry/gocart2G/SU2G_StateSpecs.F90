! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module SU2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use SU2G_instance,only: nbins

 implicit none
 public

!this module is the state variable specification file for sulfur parameters. it is the same as SU2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SU2G_GridComp.

!schema_version: 2.0.0
!component: SU


 type SU2G_State
 
!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer:: lats              => null() !latitude (radian)
 real(kind=RKIND),dimension(:,:),pointer:: lons              => null() !longitude (radian)
 real(kind=RKIND),dimension(:,:),pointer:: frocean           => null() !fraction_of_ocean (1)
 real(kind=RKIND),dimension(:,:),pointer:: lwi               => null() !land-ocean-ice_mask (1)
 real(kind=RKIND),dimension(:,:),pointer:: tropp             => null() !tropopause_pressure_based_on_blended_estimate (pa)
 real(kind=RKIND),dimension(:,:),pointer:: u10m              => null() !10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer:: v10m              => null() !10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer:: area              => null() !grid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer:: zpbl              => null() !planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer:: ustar             => null() !surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer:: sh                => null() !sensible_heat_flux_from_turbulence (w m-2)
 real(kind=RKIND),dimension(:,:),pointer:: z0h               => null() !surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer:: cn_prcp           => null() !surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer:: ncn_prcp          => null() !non-convective precipitation (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer:: coszr             => null() !cosine of the solar zenith angle (-).
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: airdens         => null() !moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer:: delp            => null() !pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: delz            => null() !geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: t               => null() !air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: rh2             => null() !rel_hum_after_moist (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: zle             => null() !geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: ple             => null() !air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfl_lsan        => null() !3d_flux_of_liquid_nonconvective_precipitation (kg/m2s)
 real(kind=RKIND),dimension(:,:,:),pointer:: pfi_lsan        => null() !3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: u               => null() !eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: v               => null() !northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: fcld            => null() !cloud fraction for radiation (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso2_ocs        => null() !source species (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: su_no3          => null() !climatological no3 source (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: su_oh           => null() !climatological oh source (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: su_h2o2         => null() !climatological h2o2 source (1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer:: su_biomass        => null() !biomass burning emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_anthrol1       => null() !anthropogenic bf emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_anthrol2       => null() !anthropogenic ff emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_shipso2        => null() !so2 ship emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_shipso4        => null() !so4 ship emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_dmso           => null() !dms emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_aviation_lto   => null() !landing/take-off aircraft source species (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_aviation_cds   => null() !climb/descent aircraft source species (1)
 real(kind=RKIND),dimension(:,:),pointer:: su_aviation_crs   => null() !cruise aircraft source species (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: su_aircraft     => null() !fuel source species (1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer:: suem            => null() !sulfur emission (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: sudp            => null() !sulfate dry deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: susd            => null() !sulfate settling (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: suwt            => null() !sulfate wet deposition (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: susv            => null() !sulfate convective scavenging (bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: so4eman           => null() !so4 anthropogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: so2eman           => null() !so2 anthropogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: so2embb           => null() !so2 biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: so2emvn           => null() !so2 volcanic (non-explosive) emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: so2emve           => null() !so2 volcanic (explosive) emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso2            => null() !so2 prod from dms oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pmsa            => null() !msa prod from dms oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso4            => null() !so4 prod from all so2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso4g           => null() !so4 prod from gaseous so2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso4wet         => null() !so4 prod from wet so2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso4aq          => null() !so4 prod from aqueous so2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: supso2            => null() !so2 prod from dms oxidation [column] (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: supso4            => null() !so4 prod from all so2 oxidation [column] (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: supso4g           => null() !so4 prod from gaseous so2 oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: supso4aq          => null() !so4 prod from aqueous so2 oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: supso4wt          => null() !so4 prod from aqueous so2 oxidation [wet dep] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: supmsa            => null() !msa prod from dms oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: so2smass          => null() !so2 surface mass concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: so2cmass          => null() !so2 column mass density (kg m-2)     
 real(kind=RKIND),dimension(:,:),pointer:: so4smass          => null() !so4 surface mass concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: so4cmass          => null() !so4 column mass density (kg m-2)     
 real(kind=RKIND),dimension(:,:),pointer:: dmssmass          => null() !dms surface mass concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: dmscmass          => null() !dms column mass density (kg m-2)     
 real(kind=RKIND),dimension(:,:),pointer:: msasmass          => null() !msa surface mass concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: msacmass          => null() !msa column mass density (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: suconc        => null() !so4 aerosol mass concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suextcoef     => null() !so4 extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suextcoefrh20 => null() !so4 extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suextcoefrh80 => null() !so4 extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suscacoef     => null() !so4 scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suscacoefrh20 => null() !so4 scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suscacoefrh80 => null() !so4 scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: subckcoef     => null() !so4 backscatter coefficient (m-1 sr-1) 

 real(kind=RKIND),dimension(:,:),pointer:: suangstr          => null() !so4 angstrom parameter [470-870 nm] (1)
 real(kind=RKIND),dimension(:,:),pointer:: sufluxu           => null() !so4 column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: sufluxv           => null() !so4 column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: so4mass         => null() !so4 aerosol mass mixing ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: suexttau        => null() !so4 extinction aot (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: sustexttau      => null() !so4 extinction aot stratosphere (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: suscatau        => null() !so4 scattering aot (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: sustscatau      => null() !so4 scattering aot stratosphere (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: so4sarea        => null() !so4 surface area density (m2 m-3 )
 real(kind=RKIND),dimension(:,:,:),pointer:: so4snum         => null() !so4 number density (m-3)
 real(kind=RKIND),dimension(:,:),pointer  :: suvdep          => null() ! dry deposition velocity (m s-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:,:),pointer:: sutau_lw      => null() ! sulfate optical depth for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suasy_lw      => null() ! sulfate asymmetry factor for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: sussa_lw      => null() ! sulfate single scattering albedo for longwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: sutau_sw      => null() ! sulfate optical depth for shortwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: suasy_sw      => null() ! sulfate asymmetry factor for shortwave RRTMG (-)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: sussa_sw      => null() ! sulfate single scattering albedo for shortwave RRTMG (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer:: dms             => null() !dimethylsulphide (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: so2             => null() !sulphur dioxide (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: so4             => null() !sulphate aerosol (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: msa             => null() !methanesulphonic acid (kg kg-1) 
 real(kind=RKIND),dimension(:,:,:),pointer:: h2o2_init       => null() !private H2O2 (kg kg-1)


 contains
    procedure:: gocart2G_allocate   => SU2G_StateSpecsInit
    procedure:: gocart2G_deallocate => SU2G_StateSpecsFinalize

 end type SU2G_State


 contains


!=================================================================================================================
 subroutine SU2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte,nbndlw,nbndsw)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nbndlw,nbndsw

!--- inout arguments:
 class(SU2G_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!-----------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
!if(.not.associated(self%lats)           ) allocate(self%lats(its:ite,jts:jte)                            )
!if(.not.associated(self%lons)           ) allocate(self%lons(its:ite,jts:jte)                            )
!if(.not.associated(self%frocean)        ) allocate(self%frocean(its:ite,jts:jte)                         )
!if(.not.associated(self%lwi)            ) allocate(self%lwi(its:ite,jts:jte)                             )
!if(.not.associated(self%tropp)          ) allocate(self%tropp(its:ite,jts:jte)                           )
!if(.not.associated(self%u10m)           ) allocate(self%u10m(its:ite,jts:jte)                            )
!if(.not.associated(self%v10m)           ) allocate(self%v10m(its:ite,jts:jte)                            )
!if(.not.associated(self%area)           ) allocate(self%area(its:ite,jts:jte)                            )
!if(.not.associated(self%zpbl)           ) allocate(self%zpbl(its:ite,jts:jte)                            )
!if(.not.associated(self%ustar)          ) allocate(self%ustar(its:ite,jts:jte)                           )
!if(.not.associated(self%sh)             ) allocate(self%sh(its:ite,jts:jte)                              )
!if(.not.associated(self%z0h)            ) allocate(self%z0h(its:ite,jts:jte)                             )
!if(.not.associated(self%cn_prcp)        ) allocate(self%cn_prcp(its:ite,jts:jte)                         )
!if(.not.associated(self%ncn_prcp)       ) allocate(self%ncn_prcp(its:ite,jts:jte)                        )
!if(.not.associated(self%coszr)          ) allocate(self%coszr(its:ite,jts:jte)                           )
!.................................................................................................................
!if(.not.associated(self%airdens)        ) allocate(self%airdens(its:ite,jts:jte,kts:kte)                 )
!if(.not.associated(self%delp)           ) allocate(self%delp(its:ite,jts:jte,kts:kte)                    )
!if(.not.associated(self%delz)           ) allocate(self%delz(its:ite,jts:jte,kts:kte)                    )
!if(.not.associated(self%t)              ) allocate(self%t(its:ite,jts:jte,kts:kte)                       )
!if(.not.associated(self%rh2)            ) allocate(self%rh2(its:ite,jts:jte,kts:kte)                     )
!if(.not.associated(self%zle)            ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)                   )
!if(.not.associated(self%ple)            ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)                   )
!if(.not.associated(self%pfl_lsan)       ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)                )
!if(.not.associated(self%pfi_lsan)       ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)                )
!if(.not.associated(self%u)              ) allocate(self%u(its:ite,jts:jte,kts:kte)                       )
!if(.not.associated(self%v)              ) allocate(self%v(its:ite,jts:jte,kts:kte)                       )
!if(.not.associated(self%fcld)           ) allocate(self%fcld(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%pso2_ocs)       ) allocate(self%pso2_ocs(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%su_no3)         ) allocate(self%su_no3(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%su_oh)          ) allocate(self%su_oh(its:ite,jts:jte,kts:kte)                   )
 if(.not.associated(self%su_h2o2)        ) allocate(self%su_h2o2(its:ite,jts:jte,kts:kte)                 )
!.................................................................................................................
!if(.not.associated(self%su_biomass)     ) allocate(self%su_biomass(its:ite,jts:jte)                      )
!if(.not.associated(self%su_anthrol1)    ) allocate(self%su_anthrol1(its:ite,jts:jte)                     )
!if(.not.associated(self%su_anthrol2)    ) allocate(self%su_anthrol2(its:ite,jts:jte)                     )
!if(.not.associated(self%su_shipso2)     ) allocate(self%su_shipso2(its:ite,jts:jte)                      )
!if(.not.associated(self%su_shipso4)     ) allocate(self%su_shipso4(its:ite,jts:jte)                      )
!if(.not.associated(self%su_dmso)        ) allocate(self%su_dmso(its:ite,jts:jte)                         )
!if(.not.associated(self%su_aviation_lto)) allocate(self%su_aviation_lto(its:ite,jts:jte)                 )
!if(.not.associated(self%su_aviation_cds)) allocate(self%su_aviation_cds(its:ite,jts:jte)                 )
!if(.not.associated(self%su_aviation_crs)) allocate(self%su_aviation_crs(its:ite,jts:jte)                 )
!if(.not.associated(self%su_aircraft)    ) allocate(self%su_aircraft(its:ite,jts:jte,kts:kte)             )

!category: EXPORT
 if(.not.associated(self%suem)           ) allocate(self%suem(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%sudp)           ) allocate(self%sudp(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%susd)           ) allocate(self%susd(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%suwt)           ) allocate(self%suwt(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%susv)           ) allocate(self%susv(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%so4eman)        ) allocate(self%so4eman(its:ite,jts:jte)                         )
 if(.not.associated(self%so2eman)        ) allocate(self%so2eman(its:ite,jts:jte)                         )
 if(.not.associated(self%so2embb)        ) allocate(self%so2embb(its:ite,jts:jte)                         )
 if(.not.associated(self%so2emvn)        ) allocate(self%so2emvn(its:ite,jts:jte)                         )
 if(.not.associated(self%so2emve)        ) allocate(self%so2emve(its:ite,jts:jte)                         )
 if(.not.associated(self%pso2)           ) allocate(self%pso2(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%pmsa)           ) allocate(self%pmsa(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%pso4)           ) allocate(self%pso4(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%pso4g)          ) allocate(self%pso4g(its:ite,jts:jte,kts:kte)                   )
 if(.not.associated(self%pso4wet)        ) allocate(self%pso4wet(its:ite,jts:jte,kts:kte)                 )
 if(.not.associated(self%pso4aq)         ) allocate(self%pso4aq(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%supso2)         ) allocate(self%supso2(its:ite,jts:jte)                          )
 if(.not.associated(self%supso4)         ) allocate(self%supso4(its:ite,jts:jte)                          )
 if(.not.associated(self%supso4g)        ) allocate(self%supso4g(its:ite,jts:jte)                         )
 if(.not.associated(self%supso4aq)       ) allocate(self%supso4aq(its:ite,jts:jte)                        )
 if(.not.associated(self%supso4wt)       ) allocate(self%supso4wt(its:ite,jts:jte)                        )
 if(.not.associated(self%supmsa)         ) allocate(self%supmsa(its:ite,jts:jte)                          )
 if(.not.associated(self%so2smass)       ) allocate(self%so2smass(its:ite,jts:jte)                        )
 if(.not.associated(self%so2cmass)       ) allocate(self%so2cmass(its:ite,jts:jte)                        )
 if(.not.associated(self%so4smass)       ) allocate(self%so4smass(its:ite,jts:jte)                        )
 if(.not.associated(self%so4cmass)       ) allocate(self%so4cmass(its:ite,jts:jte)                        )
 if(.not.associated(self%dmssmass)       ) allocate(self%dmssmass(its:ite,jts:jte)                        )
 if(.not.associated(self%dmscmass)       ) allocate(self%dmscmass(its:ite,jts:jte)                        )
 if(.not.associated(self%msasmass)       ) allocate(self%msasmass(its:ite,jts:jte)                        )
 if(.not.associated(self%msacmass)       ) allocate(self%msacmass(its:ite,jts:jte)                        )
 if(.not.associated(self%suconc)         ) allocate(self%suconc(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%suextcoef)      ) allocate(self%suextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%suextcoefrh20)  ) allocate(self%suextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%suextcoefrh80)  ) allocate(self%suextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%suscacoef)      ) allocate(self%suscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%suscacoefrh20)  ) allocate(self%suscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%suscacoefrh80)  ) allocate(self%suscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile)) 
 if(.not.associated(self%subckcoef)      ) allocate(self%subckcoef(its:ite,jts:jte,kts:kte,nw_profile)    ) 
 if(.not.associated(self%suangstr)       ) allocate(self%suangstr(its:ite,jts:jte)                        )
 if(.not.associated(self%sufluxu)        ) allocate(self%sufluxu(its:ite,jts:jte)                         )
 if(.not.associated(self%sufluxu)        ) allocate(self%sufluxu(its:ite,jts:jte)                         )
 if(.not.associated(self%so4mass)        ) allocate(self%so4mass(its:ite,jts:jte,kts:kte)                 )
 if(.not.associated(self%suexttau)       ) allocate(self%suexttau(its:ite,jts:jte,nw_vertint)             )
 if(.not.associated(self%sustexttau)     ) allocate(self%sustexttau(its:ite,jts:jte,nw_vertint)           )
 if(.not.associated(self%suscatau)       ) allocate(self%suscatau(its:ite,jts:jte,nw_vertint)             )
 if(.not.associated(self%sustscatau)     ) allocate(self%sustscatau(its:ite,jts:jte,nw_vertint)           )
 if(.not.associated(self%so4sarea)       ) allocate(self%so4sarea(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%so4snum)        ) allocate(self%so4snum(its:ite,jts:jte,kts:kte)                 )
 if(.not.associated(self%suvdep)         ) allocate(self%suvdep(its:ite,jts:jte)                          )
!..................................................................................................................
 if(.not.associated(self%sutau_lw)       ) allocate(self%sutau_lw(its:ite,jts:jte,kts:kte,nbndlw)         )
 if(.not.associated(self%sussa_lw)       ) allocate(self%sussa_lw(its:ite,jts:jte,kts:kte,nbndlw)         )
 if(.not.associated(self%suasy_lw)       ) allocate(self%suasy_lw(its:ite,jts:jte,kts:kte,nbndlw)         )
 if(.not.associated(self%sutau_sw)       ) allocate(self%sutau_sw(its:ite,jts:jte,kts:kte,nbndsw)         )
 if(.not.associated(self%sussa_sw)       ) allocate(self%sussa_sw(its:ite,jts:jte,kts:kte,nbndsw)         )
 if(.not.associated(self%suasy_sw)       ) allocate(self%suasy_sw(its:ite,jts:jte,kts:kte,nbndsw)         )

!category: INTERNAL
!if(.not.associated(self%dms)            ) allocate(self%dms(its:ite,jts:jte,kts:kte)                     )
!if(.not.associated(self%so2)            ) allocate(self%so2(its:ite,jts:jte,kts:kte)                     )
!if(.not.associated(self%so4)            ) allocate(self%so4(its:ite,jts:jte,kts:kte)                     )
!if(.not.associated(self%msa)            ) allocate(self%msa(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%h2o2_init)      ) allocate(self%h2o2_init(its:ite,jts:jte,kts:kte)               )


!--- initialization of diagnostics used in GOCART2G_GridComp:
 self%suangstr(:,:) = 0._RKIND
 self%so4smass(:,:) = 0._RKIND

 self%suexttau(:,:,:)   = 0._RKIND
 self%sustexttau(:,:,:) = 0._RKIND
 self%suscatau(:,:,:)   = 0._RKIND
 self%sustscatau(:,:,:) = 0._RKIND

 self%suextcoef(:,:,:,:)     = 0._RKIND
 self%suextcoefrh20(:,:,:,:) = 0._RKIND
 self%suextcoefrh80(:,:,:,:) = 0._RKIND
 self%suscacoef(:,:,:,:)     = 0._RKIND
 self%suscacoefrh20(:,:,:,:) = 0._RKIND
 self%suscacoefrh80(:,:,:,:) = 0._RKIND
 self%subckcoef(:,:,:,:)     = 0._RKIND


!--- initialization of RRTMG longwave and shortwave optical properties:
 self%sutau_lw(:,:,:,:) = 0._RKIND
 self%sussa_lw(:,:,:,:) = 0._RKIND
 self%suasy_lw(:,:,:,:) = 0._RKIND
 self%sutau_sw(:,:,:,:) = 0._RKIND
 self%sussa_sw(:,:,:,:) = 0._RKIND
 self%suasy_sw(:,:,:,:) = 0._RKIND


 end subroutine SU2G_StateSpecsInit

!=================================================================================================================
 subroutine SU2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(SU2G_State),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
!if(associated(self%lats)           ) deallocate(self%lats           )
!if(associated(self%lons)           ) deallocate(self%lons           )
!if(associated(self%frocean)        ) deallocate(self%frocean        )
!if(associated(self%lwi)            ) deallocate(self%lwi            )
!if(associated(self%tropp)          ) deallocate(self%tropp          )
!if(associated(self%u10m)           ) deallocate(self%u10m           )
!if(associated(self%v10m)           ) deallocate(self%v10m           )
!if(associated(self%area)           ) deallocate(self%area           )
!if(associated(self%zpbl)           ) deallocate(self%zpbl           )
!if(associated(self%ustar)          ) deallocate(self%ustar          )
!if(associated(self%sh)             ) deallocate(self%sh             )
!if(associated(self%z0h)            ) deallocate(self%z0h            )
!if(associated(self%cn_prcp)        ) deallocate(self%cn_prcp        )
!if(associated(self%ncn_prcp)       ) deallocate(self%ncn_prcp       )
!if(associated(self%coszr)          ) deallocate(self%coszr          )
!.................................................................................................................
!if(associated(self%airdens)        ) deallocate(self%airdens        )
!if(associated(self%delp)           ) deallocate(self%delp           )
!if(associated(self%delz)           ) deallocate(self%delz           )
!if(associated(self%t)              ) deallocate(self%t              )
!if(associated(self%rh2)            ) deallocate(self%rh2            )
!if(associated(self%zle)            ) deallocate(self%zle            )
!if(associated(self%ple)            ) deallocate(self%ple            )
!if(associated(self%pfl_lsan)       ) deallocate(self%pfl_lsan       )
!if(associated(self%pfi_lsan)       ) deallocate(self%pfi_lsan       )
!if(associated(self%u)              ) deallocate(self%u              )
!if(associated(self%v)              ) deallocate(self%v              )
!if(associated(self%fcld)           ) deallocate(self%fcld           )
 if(associated(self%pso2_ocs)       ) deallocate(self%pso2_ocs       )
 if(associated(self%su_no3)         ) deallocate(self%su_no3         )
 if(associated(self%su_oh)          ) deallocate(self%su_oh          )
 if(associated(self%su_h2o2)        ) deallocate(self%su_h2o2        )
!.................................................................................................................
!if(associated(self%su_biomass)     ) deallocate(self%su_biomass     )
!if(associated(self%su_anthrol1)    ) deallocate(self%su_anthrol1    )
!if(associated(self%su_anthrol2)    ) deallocate(self%su_anthrol2    )
!if(associated(self%su_shipso2)     ) deallocate(self%su_shipso2     )
!if(associated(self%su_shipso4)     ) deallocate(self%su_shipso4     )
!if(associated(self%su_dmso)        ) deallocate(self%su_dmso        )
!if(associated(self%su_aviation_lto)) deallocate(self%su_aviation_lto)
!if(associated(self%su_aviation_cds)) deallocate(self%su_aviation_cds)
!if(associated(self%su_aviation_crs)) deallocate(self%su_aviation_crs)
!if(associated(self%su_aircraft)    ) deallocate(self%su_aircraft    )

!category: EXPORT
 if(associated(self%suem)           ) deallocate(self%suem           )
 if(associated(self%sudp)           ) deallocate(self%sudp           )
 if(associated(self%susd)           ) deallocate(self%susd           )
 if(associated(self%suwt)           ) deallocate(self%suwt           )
 if(associated(self%susv)           ) deallocate(self%susv           )
 if(associated(self%so4eman)        ) deallocate(self%so4eman        )
 if(associated(self%so2eman)        ) deallocate(self%so2eman        )
 if(associated(self%so2embb)        ) deallocate(self%so2embb        )
 if(associated(self%so2emvn)        ) deallocate(self%so2emvn        )
 if(associated(self%so2emve)        ) deallocate(self%so2emve        )
 if(associated(self%pso2)           ) deallocate(self%pso2           )
 if(associated(self%pmsa)           ) deallocate(self%pmsa           )
 if(associated(self%pso4)           ) deallocate(self%pso4           )
 if(associated(self%pso4g)          ) deallocate(self%pso4g          )
 if(associated(self%pso4wet)        ) deallocate(self%pso4wet        )
 if(associated(self%pso4aq)         ) deallocate(self%pso4aq         )
 if(associated(self%supso2)         ) deallocate(self%supso2         )
 if(associated(self%supso4)         ) deallocate(self%supso4         )
 if(associated(self%supso4g)        ) deallocate(self%supso4g        )
 if(associated(self%supso4aq)       ) deallocate(self%supso4aq       )
 if(associated(self%supso4wt)       ) deallocate(self%supso4wt       )
 if(associated(self%supmsa)         ) deallocate(self%supmsa         )
 if(associated(self%so2smass)       ) deallocate(self%so2smass       )
 if(associated(self%so2cmass)       ) deallocate(self%so2cmass       )
 if(associated(self%so4smass)       ) deallocate(self%so4smass       )
 if(associated(self%so4cmass)       ) deallocate(self%so4cmass       )
 if(associated(self%dmssmass)       ) deallocate(self%dmssmass       )
 if(associated(self%dmscmass)       ) deallocate(self%dmscmass       )
 if(associated(self%msasmass)       ) deallocate(self%msasmass       )
 if(associated(self%msacmass)       ) deallocate(self%msacmass       )
 if(associated(self%suconc)         ) deallocate(self%suconc         )

!if(associated(self%suextcoef)      ) deallocate(self%suextcoef      )
!if(associated(self%suextcoefrh20)  ) deallocate(self%suextcoefrh20  )
!if(associated(self%suextcoefrh80)  ) deallocate(self%suextcoefrh80  )
!if(associated(self%suscacoef)      ) deallocate(self%suscacoef      )
!if(associated(self%suscacoefrh20)  ) deallocate(self%suscacoefrh20  )
!if(associated(self%suscacoefrh80)  ) deallocate(self%suscacoefrh80  )
!if(associated(self%subckcoef)      ) deallocate(self%subckcoef      )
 if(associated(self%suangstr)       ) deallocate(self%suangstr       )
 if(associated(self%sufluxu)        ) deallocate(self%sufluxu        )
 if(associated(self%sufluxu)        ) deallocate(self%sufluxu        )
 if(associated(self%so4mass)        ) deallocate(self%so4mass        )
!if(associated(self%suexttau)       ) deallocate(self%suexttau       )
!if(associated(self%sustexttau)     ) deallocate(self%sustexttau     )
!if(associated(self%suscatau)       ) deallocate(self%suscatau       )
!if(associated(self%sustscatau)     ) deallocate(self%sustscatau     )
 if(associated(self%so4sarea)       ) deallocate(self%so4sarea       )
 if(associated(self%so4snum)        ) deallocate(self%so4snum        )
 if(associated(self%suvdep)         ) deallocate(self%suvdep         )
!..................................................................................................................
 if(associated(self%sutau_lw)       ) deallocate(self%sutau_lw       )
 if(associated(self%sussa_lw)       ) deallocate(self%sussa_lw       )
 if(associated(self%suasy_lw)       ) deallocate(self%suasy_lw       )
 if(associated(self%sutau_sw)       ) deallocate(self%sutau_sw       )
 if(associated(self%sussa_sw)       ) deallocate(self%sussa_sw       )
 if(associated(self%suasy_sw)       ) deallocate(self%suasy_sw       )

!category: INTERNAL
!if(associated(self%dms)            ) deallocate(self%dms            )
!if(associated(self%so2)            ) deallocate(self%so2            )
!if(associated(self%so4)            ) deallocate(self%so4            )
!if(associated(self%msa)            ) deallocate(self%msa            )
 if(associated(self%h2o2_init)      ) deallocate(self%h2o2_init      )

 end subroutine SU2G_StateSpecsFinalize

!=================================================================================================================
 end module SU2G_StateSpecs
!=================================================================================================================
