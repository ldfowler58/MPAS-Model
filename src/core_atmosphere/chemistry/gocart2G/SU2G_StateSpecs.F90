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
 real(kind=RKIND),dimension(:,:,:),pointer:: su_aircraft     => null() !fuel source species (1)
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
 subroutine SU2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(SU2G_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!-----------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
 if(.not.associated(self%lats)           ) allocate(self%lats(its:ite,jts:jte)                            )
 if(.not.associated(self%lons)           ) allocate(self%lons(its:ite,jts:jte)                            )
 if(.not.associated(self%frocean)        ) allocate(self%frocean(its:ite,jts:jte)                         )
 if(.not.associated(self%lwi)            ) allocate(self%lwi(its:ite,jts:jte)                             )
 if(.not.associated(self%tropp)          ) allocate(self%tropp(its:ite,jts:jte)                           )
 if(.not.associated(self%u10m)           ) allocate(self%u10m(its:ite,jts:jte)                            )
 if(.not.associated(self%v10m)           ) allocate(self%v10m(its:ite,jts:jte)                            )
 if(.not.associated(self%area)           ) allocate(self%area(its:ite,jts:jte)                            )
 if(.not.associated(self%zpbl)           ) allocate(self%zpbl(its:ite,jts:jte)                            )
 if(.not.associated(self%ustar)          ) allocate(self%ustar(its:ite,jts:jte)                           )
 if(.not.associated(self%sh)             ) allocate(self%sh(its:ite,jts:jte)                              )
 if(.not.associated(self%z0h)            ) allocate(self%z0h(its:ite,jts:jte)                             )
 if(.not.associated(self%cn_prcp)        ) allocate(self%cn_prcp(its:ite,jts:jte)                         )
 if(.not.associated(self%ncn_prcp)       ) allocate(self%ncn_prcp(its:ite,jts:jte)                        )
 if(.not.associated(self%coszr)          ) allocate(self%coszr(its:ite,jts:jte)                           )
!.................................................................................................................
 if(.not.associated(self%airdens)        ) allocate(self%airdens(its:ite,jts:jte,kts:kte)                 )
 if(.not.associated(self%delp)           ) allocate(self%delp(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%delz)           ) allocate(self%delz(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%t)              ) allocate(self%t(its:ite,jts:jte,kts:kte)                       )
 if(.not.associated(self%rh2)            ) allocate(self%rh2(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%zle)            ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)                   )
 if(.not.associated(self%ple)            ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)                   )
 if(.not.associated(self%pfl_lsan)       ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%pfi_lsan)       ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%u)              ) allocate(self%u(its:ite,jts:jte,kts:kte)                       )
 if(.not.associated(self%v)              ) allocate(self%v(its:ite,jts:jte,kts:kte)                       )
 if(.not.associated(self%fcld)           ) allocate(self%fcld(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%pso2_ocs)       ) allocate(self%pso2_ocs(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%su_aircraft)    ) allocate(self%su_aircraft(its:ite,jts:jte,kts:kte)             )
 if(.not.associated(self%su_no3)         ) allocate(self%su_no3(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%su_oh)          ) allocate(self%su_oh(its:ite,jts:jte,kts:kte)                   )
 if(.not.associated(self%su_h2o2)        ) allocate(self%su_h2o2(its:ite,jts:jte,kts:kte)                 )
!.................................................................................................................
 if(.not.associated(self%su_biomass)     ) allocate(self%su_biomass(its:ite,jts:jte)                      )
 if(.not.associated(self%su_anthrol1)    ) allocate(self%su_anthrol1(its:ite,jts:jte)                     )
 if(.not.associated(self%su_anthrol2)    ) allocate(self%su_anthrol2(its:ite,jts:jte)                     )
 if(.not.associated(self%su_shipso2)     ) allocate(self%su_shipso2(its:ite,jts:jte)                      )
 if(.not.associated(self%su_shipso4)     ) allocate(self%su_shipso4(its:ite,jts:jte)                      )
 if(.not.associated(self%su_dmso)        ) allocate(self%su_dmso(its:ite,jts:jte)                         )
 if(.not.associated(self%su_aviation_lto)) allocate(self%su_aviation_lto(its:ite,jts:jte)                 )
 if(.not.associated(self%su_aviation_cds)) allocate(self%su_aviation_cds(its:ite,jts:jte)                 )
 if(.not.associated(self%su_aviation_crs)) allocate(self%su_aviation_crs(its:ite,jts:jte)                 )

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
 if(.not.associated(self%so4smass)       ) allocate(self%so4smass(its:ite,jts:jte)                        )
 if(.not.associated(self%so4smass)       ) allocate(self%so4smass(its:ite,jts:jte)                        )
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

!category: INTERNAL
 if(.not.associated(self%dms)            ) allocate(self%dms(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%so2)            ) allocate(self%so2(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%so4)            ) allocate(self%so4(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%msa)            ) allocate(self%msa(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%h2o2_init)      ) allocate(self%h2o2_init(its:ite,jts:jte,kts:kte)               )

 end subroutine SU2G_StateSpecsInit

!=================================================================================================================
 subroutine SU2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(SU2G_State),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%lats)           ) deallocate(self%lats           )
 if(associated(self%lons)           ) deallocate(self%lons           )
 if(associated(self%frocean)        ) deallocate(self%frocean        )
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
 if(associated(self%coszr)          ) deallocate(self%coszr          )
!.................................................................................................................
 if(associated(self%airdens)        ) deallocate(self%airdens        )
 if(associated(self%delp)           ) deallocate(self%delp           )
 if(associated(self%delz)           ) deallocate(self%delz           )
 if(associated(self%t)              ) deallocate(self%t              )
 if(associated(self%rh2)            ) deallocate(self%rh2            )
 if(associated(self%zle)            ) deallocate(self%zle            )
 if(associated(self%ple)            ) deallocate(self%ple            )
 if(associated(self%pfl_lsan)       ) deallocate(self%pfl_lsan       )
 if(associated(self%pfi_lsan)       ) deallocate(self%pfi_lsan       )
 if(associated(self%u)              ) deallocate(self%u              )
 if(associated(self%v)              ) deallocate(self%v              )
 if(associated(self%fcld)           ) deallocate(self%fcld           )
 if(associated(self%PSO2_OCS)       ) deallocate(self%pSO2_OCS       )
 if(associated(self%SU_AIRCRAFT)    ) deallocate(self%SU_AIRCRAFT    )
 if(associated(self%SU_NO3)         ) deallocate(self%SU_NO3         )
 if(associated(self%SU_OH)          ) deallocate(self%SU_OH          )
 if(associated(self%SU_H2O2)        ) deallocate(self%SU_H2O2        )
!.................................................................................................................
 if(associated(self%SU_BIOMASS)     ) deallocate(self%SU_BIOMASS     )
 if(associated(self%SU_ANTHROL1)    ) deallocate(self%SU_ANTHROL1    )
 if(associated(self%SU_ANTHROL2)    ) deallocate(self%SU_ANTHROL2    )
 if(associated(self%SU_SHIPSO2)     ) deallocate(self%SU_SHIPSO2     )
 if(associated(self%SU_SHIPSO4)     ) deallocate(self%SU_SHIPSO4     )
 if(associated(self%SU_DMSO)        ) deallocate(self%SU_DMSO        )
 if(associated(self%SU_AVIATION_LTO)) deallocate(self%SU_AVIATION_LTO)
 if(associated(self%SU_AVIATION_CDS)) deallocate(self%SU_AVIATION_CDS)
 if(associated(self%SU_AVIATION_CRS)) deallocate(self%SU_AVIATION_CRS)

!category: EXPORT
 if(associated(self%SUEM)           ) deallocate(self%SUEM           )
 if(associated(self%SUDP)           ) deallocate(self%SUDP           )
 if(associated(self%SUSD)           ) deallocate(self%SUSD           )
 if(associated(self%SUWT)           ) deallocate(self%SUWT           )
 if(associated(self%SUSV)           ) deallocate(self%SUSV           )
 if(associated(self%SO4EMAN)        ) deallocate(self%SO4EMAN        )
 if(associated(self%SO2EMAN)        ) deallocate(self%SO2EMAN        )
 if(associated(self%SO2EMBB)        ) deallocate(self%SO2EMBB        )
 if(associated(self%SO2EMVN)        ) deallocate(self%SO2EMVN        )
 if(associated(self%SO2EMVE)        ) deallocate(self%SO2EMVE        )
 if(associated(self%PSO2)           ) deallocate(self%PSO2           )
 if(associated(self%PMSA)           ) deallocate(self%PMSA           )
 if(associated(self%PSO4)           ) deallocate(self%PSO4           )
 if(associated(self%PSO4G)          ) deallocate(self%PSO4G          )
 if(associated(self%PSO4WET)        ) deallocate(self%PSO4WET        )
 if(associated(self%PSO4AQ)         ) deallocate(self%PSO4AQ         )
 if(associated(self%SUPSO2)         ) deallocate(self%SUPSO2         )
 if(associated(self%SUPSO4)         ) deallocate(self%SUPSO4         )
 if(associated(self%SUPSO4G)        ) deallocate(self%SUPSO4G        )
 if(associated(self%SUPSO4AQ)       ) deallocate(self%SUPSO4AQ       )
 if(associated(self%SUPSO4WT)       ) deallocate(self%SUPSO4WT       )
 if(associated(self%SUPMSA)         ) deallocate(self%SUPMSA         )
 if(associated(self%SO2SMASS)       ) deallocate(self%SO2SMASS       )
 if(associated(self%SO2CMASS)       ) deallocate(self%SO2CMASS       )
 if(associated(self%SO4SMASS)       ) deallocate(self%SO4SMASS       )
 if(associated(self%SO4SMASS)       ) deallocate(self%SO4SMASS       )
 if(associated(self%SO4SMASS)       ) deallocate(self%SO4SMASS       )
 if(associated(self%DMSCMASS)       ) deallocate(self%DMSCMASS       )
 if(associated(self%MSASMASS)       ) deallocate(self%MSASMASS       )
 if(associated(self%MSACMASS)       ) deallocate(self%MSACMASS       )
 if(associated(self%SUCONC)         ) deallocate(self%SUCONC         )

!if(associated(self%SUEXTCOEF)      ) deallocate(self%SUEXTCOEF      )
!if(associated(self%SUEXTCOEFRH20)  ) deallocate(self%SUEXTCOEFRH20  )
!if(associated(self%SUEXTCOEFRH80)  ) deallocate(self%SUEXTCOEFRH80  )
!if(associated(self%SUSCACOEF)      ) deallocate(self%SUSCACOEF      )
!if(associated(self%SUSCACOEFRH20)  ) deallocate(self%SUSCACOEFRH20  )
!if(associated(self%SUSCACOEFRH80)  ) deallocate(self%SUSCACOEFRH80  ) 
!if(associated(self%SUBCKCOEF)      ) deallocate(self%SUBCKCOEF      ) 
 if(associated(self%SUANGSTR)       ) deallocate(self%SUANGSTR       )
 if(associated(self%SUFLUXU)        ) deallocate(self%SUFLUXU        )
 if(associated(self%SUFLUXU)        ) deallocate(self%SUFLUXU        )
 if(associated(self%SO4MASS)        ) deallocate(self%SO4MASS        )
!if(associated(self%SUEXTTAU)       ) deallocate(self%SUEXTTAU       )
!if(associated(self%SUSTEXTTAU)     ) deallocate(self%SUSTEXTTAU     )
!if(associated(self%SUSCATAU)       ) deallocate(self%SUSCATAU       )
!if(associated(self%SUSTSCATAU)     ) deallocate(self%SUSTSCATAU     )
 if(associated(self%SO4SAREA)       ) deallocate(self%SO4SAREA       )
 if(associated(self%SO4SNUM)        ) deallocate(self%SO4SNUM        )

!category: INTERNAL
 if(associated(self%DMS)            ) deallocate(self%DMS            )
 if(associated(self%SO2)            ) deallocate(self%SO2            )
 if(associated(self%SO4)            ) deallocate(self%SO4            )
 if(associated(self%MSA)            ) deallocate(self%MSA            )
 if(associated(self%H2O2_INIT)      ) deallocate(self%H2O2_INIT      )

 end subroutine SU2G_StateSpecsFinalize

!=================================================================================================================
 end module SU2G_StateSpecs
!=================================================================================================================
