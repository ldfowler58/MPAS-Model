!=================================================================================================================
 module SU2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_GridComp,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm
 use SU2G_instance_SU,only: nbins

 implicit none
 public

!this module is the state variable specification file for sulfur parameters. it is the same as SU2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SU2G_GridComp.

!schema_version: 2.0.0
!component: SU


 type SU2G_State
 
!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer:: LATS              => null() !latitude (radian)
 real(kind=RKIND),dimension(:,:),pointer:: LONS              => null() !longitude (radian)
 real(kind=RKIND),dimension(:,:),pointer:: FROCEAN           => null() !fraction_of_ocean (1)
 real(kind=RKIND),dimension(:,:),pointer:: LWI               => null() !land-ocean-ice_mask (1)
 real(kind=RKIND),dimension(:,:),pointer:: TROPP             => null() !tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),pointer:: U10M              => null() !10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer:: V10M              => null() !10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),pointer:: AREA              => null() !grid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),pointer:: ZPBL              => null() !planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),pointer:: USTAR             => null() !surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SH                => null() !sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),pointer:: Z0H               => null() !surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),pointer:: CN_PRCP           => null() !surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),pointer:: NCN_PRCP          => null() !non-convective precipitation (kg/m^2/s)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: AIRDENS         => null() !moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer:: DELP            => null() !pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: DELZ            => null() !geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: T               => null() !air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: RH2             => null() !rel_hum_after_moist (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: ZLE             => null() !geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: PLE             => null() !air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: PFL_LSAN        => null() !3d_flux_of_liquid_nonconvective_precipitation (kg/m2s)
 real(kind=RKIND),dimension(:,:,:),pointer:: PFI_LSAN        => null() !3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),pointer:: U               => null() !eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: V               => null() !northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: FCLD            => null() !cloud fraction for radiation (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSO2_OCS        => null() !source species (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SU_AIRCRAFT     => null() !fuel source species (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SU_NO3          => null() !climatological NO3 source (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SU_OH           => null() !climatological OH source (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SU_H2O2         => null() !climatological H2O2 source (1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer:: SU_BIOMASS        => null() !biomass burning emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_ANTHROL1       => null() !anthropogenic BF emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_ANTHROL2       => null() !anthropogenic FF emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_SHIPSO2        => null() !SO2 ship emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_SHIPSO4        => null() !SO4 ship emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_DMSO           => null() !DMS emissions (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_AVIATION_LTO   => null() !Landing/Take-off aircraft source species (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_AVIATION_CDS   => null() !Climb/Descent aircraft source species (1)
 real(kind=RKIND),dimension(:,:),pointer:: SU_AVIATION_CRS   => null() !Cruise aircraft source species (1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer:: SUEM            => null() !Sulfur Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUDP            => null() !Sulfate Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUSD            => null() !Sulfate Settling (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUWT            => null() !Sulfate Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUSV            => null() !Sulfate Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SO4EMAN           => null() !SO4 Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SO2EMAN           => null() !SO2 Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SO2EMBB           => null() !SO2 Biomass Burning Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SO2EMVN           => null() !SO2 Volcanic (non-explosive) Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SO2EMVE           => null() !SO2 Volcanic (explosive) Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSO2            => null() !SO2 Prod from DMS oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PMSA            => null() !MSA Prod from DMS oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSO4            => null() !SO4 Prod from all SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSO4G           => null() !SO4 Prod from gaseous SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSO4WET         => null() !SO4 Prod from wet SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: PSO4AQ          => null() !SO4 Prod from aqueous SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SUPSO2            => null() !SO2 Prod from DMS Oxidation [column] (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SUPSO4            => null() !SO4 Prod from All SO2 Oxidation [column] (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SUPSO4G           => null() !SO4 Prod from Gaseous SO2 Oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: SUPSO4AQ          => null() !SO4 Prod from Aqueous SO2 Oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: SUPSO4WT          => null() !SO4 Prod from Aqueous SO2 Oxidation (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: SUPMSA            => null() !MSA Prod from DMS Oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),pointer:: SO2SMASS          => null() !SO2 Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: SO2CMASS          => null() !SO2 Column Mass Density (kg m-2)     
 real(kind=RKIND),dimension(:,:),pointer:: SO4SMASS          => null() !SO4 Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: SO4CMASS          => null() !SO4 Column Mass Density (kg m-2)     
 real(kind=RKIND),dimension(:,:),pointer:: DMSSMASS          => null() !DMS Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: DMSCMASS          => null() !DMS Column Mass Density (kg m-2)     
 real(kind=RKIND),dimension(:,:),pointer:: MSASMASS          => null() !MSA Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),pointer:: MSACMASS          => null() !MSA Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer  :: SUCONC        => null() !SO4 Aerosol Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUEXTCOEF     => null() !SO4 Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUEXTCOEFRH20 => null() !SO4 Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUEXTCOEFRH80 => null() !SO4 Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUSCACOEF     => null() !SO4 Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUSCACOEFRH20 => null() !SO4 Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUSCACOEFRH80 => null() !SO4 Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: SUBCKCOEF     => null() !SO4 Backscatter Coefficient (m-1 sr-1) 

 real(kind=RKIND),dimension(:,:),pointer:: SUANGSTR          => null() !SO4 Angstrom parameter [470-870 nm] (1)
 real(kind=RKIND),dimension(:,:),pointer:: SUFLUXU           => null() !SO4 column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),pointer:: SUFLUXV           => null() !SO4 column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SO4MASS         => null() !SO4 Aerosol Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUEXTTAU        => null() !SO4 Extinction AOT (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: USTEXTTAU       => null() !SO4 Extinction AOT Stratosphere (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUSCATAU        => null() !SO4 Scattering AOT (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SUSTSCATAU      => null() !SO4 Scattering AOT Stratosphere (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SO4SAREA        => null() !SO4 Surface Area Density (m2 m-3 )
 real(kind=RKIND),dimension(:,:,:),pointer:: SO4SNUM         => null() !SO4 Number Density (m-3)

 real(kind=RKIND),dimension(:,:,:),pointer:: DMS             => null() !Dimethylsulphide (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SO2             => null() !Sulphur dioxide (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: SO4             => null() !Sulphate aerosol (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: MSA             => null() !Methanesulphonic acid (kg kg-1) 
 real(kind=RKIND),dimension(:,:,:),pointer:: H2O2_INIT       => null() !private H2O2 (kg kg-1)


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
 if(.not.associated(self%PSO2_OCS)       ) allocate(self%pSO2_OCS(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%SU_AIRCRAFT)    ) allocate(self%SU_AIRCRAFT(its:ite,jts:jte,kts:kte)             )
 if(.not.associated(self%SU_NO3)         ) allocate(self%SU_NO3(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%SU_OH)          ) allocate(self%SU_OH(its:ite,jts:jte,kts:kte)                   )
 if(.not.associated(self%SU_H2O2)        ) allocate(self%SU_H2O2(its:ite,jts:jte,kts:kte)                 )
!.................................................................................................................
 if(.not.associated(self%SU_BIOMASS)     ) allocate(self%SU_BIOMASS(its:ite,jts:jte)                      )
 if(.not.associated(self%SU_ANTHROL1)    ) allocate(self%SU_ANTHROL1(its:ite,jts:jte)                     )
 if(.not.associated(self%SU_ANTHROL2)    ) allocate(self%SU_ANTHROL2(its:ite,jts:jte)                     )
 if(.not.associated(self%SU_SHIPSO2)     ) allocate(self%SU_SHIPSO2(its:ite,jts:jte)                      )
 if(.not.associated(self%SU_SHIPSO4)     ) allocate(self%SU_SHIPSO4(its:ite,jts:jte)                      )
 if(.not.associated(self%SU_DMSO)        ) allocate(self%SU_DMSO(its:ite,jts:jte)                         )
 if(.not.associated(self%SU_AVIATION_LTO)) allocate(self%SU_AVIATION_LTO(its:ite,jts:jte)                 )
 if(.not.associated(self%SU_AVIATION_CDS)) allocate(self%SU_AVIATION_CDS(its:ite,jts:jte)                 )
 if(.not.associated(self%SU_AVIATION_CRS)) allocate(self%SU_AVIATION_CRS(its:ite,jts:jte)                 )

!category: EXPORT
 if(.not.associated(self%SUEM)           ) allocate(self%SUEM(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%SUDP)           ) allocate(self%SUDP(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%SUSD)           ) allocate(self%SUSD(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%SUWT)           ) allocate(self%SUWT(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%SUSV)           ) allocate(self%SUSV(its:ite,jts:jte,nbins)                      )
 if(.not.associated(self%SO4EMAN)        ) allocate(self%SO4EMAN(its:ite,jts:jte)                         )
 if(.not.associated(self%SO2EMAN)        ) allocate(self%SO2EMAN(its:ite,jts:jte)                         )
 if(.not.associated(self%SO2EMBB)        ) allocate(self%SO2EMBB(its:ite,jts:jte)                         )
 if(.not.associated(self%SO2EMVN)        ) allocate(self%SO2EMVN(its:ite,jts:jte)                         )
 if(.not.associated(self%SO2EMVE)        ) allocate(self%SO2EMVE(its:ite,jts:jte)                         )
 if(.not.associated(self%PSO2)           ) allocate(self%PSO2(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%PMSA)           ) allocate(self%PMSA(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%PSO4)           ) allocate(self%PSO4(its:ite,jts:jte,kts:kte)                    )
 if(.not.associated(self%PSO4G)          ) allocate(self%PSO4G(its:ite,jts:jte,kts:kte)                   )
 if(.not.associated(self%PSO4WET)        ) allocate(self%PSO4WET(its:ite,jts:jte,kts:kte)                 )
 if(.not.associated(self%PSO4AQ)         ) allocate(self%PSO4aq(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%SUPSO2)         ) allocate(self%SUPSO2(its:ite,jts:jte)                          )
 if(.not.associated(self%SUPSO4)         ) allocate(self%SUPSO4(its:ite,jts:jte)                          )
 if(.not.associated(self%SUPSO4G)        ) allocate(self%SUPSO4G(its:ite,jts:jte)                         )
 if(.not.associated(self%SUPSO4AQ)       ) allocate(self%SUPSO4AQ(its:ite,jts:jte)                        )
 if(.not.associated(self%SUPSO4WT)       ) allocate(self%SUPSO4WT(its:ite,jts:jte)                        )
 if(.not.associated(self%SUPMSA)         ) allocate(self%SUPMSA(its:ite,jts:jte)                          )
 if(.not.associated(self%SO2SMASS)       ) allocate(self%SO2SMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%SO2CMASS)       ) allocate(self%SO2CMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%SO4SMASS)       ) allocate(self%SO4SMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%SO4SMASS)       ) allocate(self%SO4SMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%SO4SMASS)       ) allocate(self%SO4SMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%DMSCMASS)       ) allocate(self%DMSCMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%MSASMASS)       ) allocate(self%MSASMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%MSACMASS)       ) allocate(self%MSACMASS(its:ite,jts:jte)                        )
 if(.not.associated(self%SUCONC)         ) allocate(self%SUCONC(its:ite,jts:jte,kts:kte)                  )
 if(.not.associated(self%SUEXTCOEF)      ) allocate(self%SUEXTCOEF(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%SUEXTCOEFRH20)  ) allocate(self%SUEXTCOEFRH20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%SUEXTCOEFRH80)  ) allocate(self%SUEXTCOEFRH80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%SUSCACOEF)      ) allocate(self%SUSCACOEF(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%SUSCACOEFRH20)  ) allocate(self%SUSCACOEFRH20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%SUSCACOEFRH80)  ) allocate(self%SUSCACOEFRH80(its:ite,jts:jte,kts:kte,nw_profile)) 
 if(.not.associated(self%SUBCKCOEF)      ) allocate(self%SUBCKCOEF(its:ite,jts:jte,kts:kte,nw_profile)    ) 
 if(.not.associated(self%SUANGSTR)       ) allocate(self%SUANGSTR(its:ite,jts:jte)                        )
 if(.not.associated(self%SUFLUXU)        ) allocate(self%SUFLUXU(its:ite,jts:jte)                         )
 if(.not.associated(self%SUFLUXU)        ) allocate(self%SUFLUXU(its:ite,jts:jte)                         )
 if(.not.associated(self%SO4MASS)        ) allocate(self%SO4MASS(its:ite,jts:jte,kts:kte)                 )
 if(.not.associated(self%SUEXTTAU)       ) allocate(self%SUEXTTAU(its:ite,jts:jte,nw_vertint)             )
 if(.not.associated(self%USTEXTTAU)      ) allocate(self%USTEXTTAU(its:ite,jts:jte,nw_vertint)            )
 if(.not.associated(self%SUSCATAU)       ) allocate(self%SUSCATAU(its:ite,jts:jte,nw_vertint)             )
 if(.not.associated(self%SUSTSCATAU)     ) allocate(self%SUSTSCATAU(its:ite,jts:jte,nw_vertint)           )
 if(.not.associated(self%SO4SAREA)       ) allocate(self%SO4SAREA(its:ite,jts:jte,kts:kte)                )
 if(.not.associated(self%SO4SNUM)        ) allocate(self%SO4SNUM(its:ite,jts:jte,kts:kte)                 )

!category: INTERNAL
 if(.not.associated(self%DMS)            ) allocate(self%DMS(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%SO2)            ) allocate(self%SO2(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%SO4)            ) allocate(self%SO4(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%MSA)            ) allocate(self%MSA(its:ite,jts:jte,kts:kte)                     )
 if(.not.associated(self%H2O2_INIT)      ) allocate(self%H2O2_INIT(its:ite,jts:jte,kts:kte)               )

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
!if(associated(self%USTEXTTAU)      ) deallocate(self%USTEXTTAU      )
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
