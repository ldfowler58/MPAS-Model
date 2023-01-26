! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_lsm_noahmp
 use mpas_kind_types
 use mpas_log
 use mpas_pool_routines
 use mpas_timer,only: mpas_timer_start, mpas_timer_stop

 use mpas_atmphys_constants,only: cp,degrad,R_d,rcp
 use mpas_atmphys_driver_lsm_shared,only: correct_tsk_over_seaice
 use mpas_atmphys_manager,only: curr_julday,year
 use mpas_atmphys_landuse,only: isurban
 use mpas_atmphys_lsm_noahmpinit,only: noahmp_init_formpas
 use mpas_atmphys_vars
 use mpas_atmphys_vars_noahmp


!wrf physics:
 use module_sf_noahmpdrv
 use module_sf_noah_seaice_drv
 use module_sf_sfcdiags_noahmp,only: sfcdiags_noahmp


 private
 public:: init_lsm_noahmp,       &
          allocate_lsm_noahmp,   &
          deallocate_lsm_noahmp, &
          driver_lsm_noahmp


!MPAS driver for the NOAHMP parameterization of land surface processes.
!Laura D. Fowler (laura@ucar.edu) / 2022-07-01

! subroutines in mpas_atmphys_driver_lsm_noahmp:
! ---------------------------------------
! allocate_lsm_noahmp  : allocate local arrays for the NOAHMP land surface parameterization.
! deallocate_lsm_noahmp: deallocate local arrays for the NOAHMP land surface parameterization.
! driver_lsm_noahmp    : main driver for the NOAHMP land surface parameterization.
! init_lsm_noahmp      : initialization of the NOAHMP land surface parameterization.
! lsm_noahmp_from_MPAS : initialize local arrays needed in the NOAHMP land surface parameterization.
! lsm_noahmp_to_MPAS   : copy local arrays to MPAS arrays.
!
! WRF physics called from driver_lsm:
! ------------------------ ----------
! * module_sf_noahmpdrv : NOAHMP land surface scheme.

 
 contains


!=================================================================================================================
 subroutine allocate_lsm_noahmp
!=================================================================================================================

!input arrays for soil layers properties:
 if(.not.allocated(dzs_p)         ) allocate(dzs_p(1:num_soils)                         )
 if(.not.allocated(smcrel_p)      ) allocate(smcrel_p(ims:ime,1:num_soils,jms:jme)      )
 if(.not.allocated(sh2o_p)        ) allocate(sh2o_p(ims:ime,1:num_soils,jms:jme)        )
 if(.not.allocated(smois_p)       ) allocate(smois_p(ims:ime,1:num_soils,jms:jme)       )
 if(.not.allocated(tslb_p)        ) allocate(tslb_p(ims:ime,1:num_soils,jms:jme)        )


!input arrays:
 if(.not.allocated(dx_p)          ) allocate(dx_p(ims:ime,jms:jme)                      )
 if(.not.allocated(area_p)        ) allocate(area_p(ims:ime,jms:jme)                    )
 if(.not.allocated(coszr_p)       ) allocate(coszr_p(ims:ime,jms:jme)                   )
 if(.not.allocated(xlat_p)        ) allocate(xlat_p(ims:ime,jms:jme)                    )
 if(.not.allocated(xlon_p)        ) allocate(xlon_p(ims:ime,jms:jme)                    )

 if(.not.allocated(isltyp_p)      ) allocate(isltyp_p(ims:ime,jms:jme)                  )
 if(.not.allocated(ivgtyp_p)      ) allocate(ivgtyp_p(ims:ime,jms:jme)                  )
 if(.not.allocated(shdmax_p)      ) allocate(shdmax_p(ims:ime,jms:jme)                  )
 if(.not.allocated(tmn_p)         ) allocate(tmn_p(ims:ime,jms:jme)                     )
 if(.not.allocated(vegfra_p)      ) allocate(vegfra_p(ims:ime,jms:jme)                  )
 if(.not.allocated(xice_p)        ) allocate(xice_p(ims:ime,jms:jme)                    )
 if(.not.allocated(xland_p)       ) allocate(xland_p(ims:ime,jms:jme)                   )

 if(.not.allocated(glw_p)         ) allocate(glw_p(ims:ime,jms:jme)                     )
 if(.not.allocated(swddir_p)      ) allocate(swddir_p(ims:ime,jms:jme)                  )
 if(.not.allocated(swddif_p)      ) allocate(swddif_p(ims:ime,jms:jme)                  )
 if(.not.allocated(swdown_p)      ) allocate(swdown_p(ims:ime,jms:jme)                  )
 if(.not.allocated(rainbl_p)      ) allocate(rainbl_p(ims:ime,jms:jme)                  )
 if(.not.allocated(sr_p)          ) allocate(sr_p(ims:ime,jms:jme)                      )


!arrays for seaice grid-points:
 if(.not.allocated(tsk_sea)       ) allocate(tsk_sea(ims:ime,jms:jme)                   )
 if(.not.allocated(tsk_ice)       ) allocate(tsk_ice(ims:ime,jms:jme)                   )
 if(.not.allocated(albsi_p)       ) allocate(albsi_p(ims:ime,jms:jme)                   )
 if(.not.allocated(icedepth_p)    ) allocate(icedepth_p(ims:ime,jms:jme)                )
 if(.not.allocated(snowsi_p)      ) allocate(snowsi_p(ims:ime,jms:jme)                  )

 if(.not.allocated(br_p)          ) allocate(br_p(ims:ime,jms:jme)                      )
 if(.not.allocated(chs_p)         ) allocate(chs_p(ims:ime,jms:jme)                     )
 if(.not.allocated(chs2_p)        ) allocate(chs2_p(ims:ime,jms:jme)                    )
 if(.not.allocated(cqs2_p)        ) allocate(cqs2_p(ims:ime,jms:jme)                    )
 if(.not.allocated(cpm_p)         ) allocate(cpm_p(ims:ime,jms:jme)                     )
 if(.not.allocated(qgh_p)         ) allocate(qgh_p(ims:ime,jms:jme)                     )
 if(.not.allocated(snoalb_p)      ) allocate(snoalb_p(ims:ime,jms:jme)                  )


!input arrays for NOAHMP only:
 if(.not.allocated(soilcl1_p)     ) allocate(soilcl1_p(ims:ime,jms:jme)                 )
 if(.not.allocated(soilcl2_p)     ) allocate(soilcl2_p(ims:ime,jms:jme)                 )
 if(.not.allocated(soilcl3_p)     ) allocate(soilcl3_p(ims:ime,jms:jme)                 )
 if(.not.allocated(soilcl4_p)     ) allocate(soilcl4_p(ims:ime,jms:jme)                 )
 if(.not.allocated(soilcomp_p)    ) allocate(soilcomp_p(ims:ime,1:num_soilcs,jms:jme)   )

 if(.not.allocated(cropcat_p)     ) allocate(cropcat_p(ims:ime,jms:jme)                 )
 if(.not.allocated(planting_p)    ) allocate(planting_p(ims:ime,jms:jme)                )
 if(.not.allocated(harvest_p)     ) allocate(harvest_p(ims:ime,jms:jme)                 )
 if(.not.allocated(planting_p)    ) allocate(planting_p(ims:ime,jms:jme)                )
 if(.not.allocated(season_gdd_p)  ) allocate(season_gdd_p(ims:ime,jms:jme)              )

 if(.not.allocated(irfract_p)     ) allocate(irfract_p(ims:ime,jms:jme)                 )
 if(.not.allocated(sifract_p)     ) allocate(sifract_p(ims:ime,jms:jme)                 )
 if(.not.allocated(mifract_p)     ) allocate(mifract_p(ims:ime,jms:jme)                 )
 if(.not.allocated(fifract_p)     ) allocate(fifract_p(ims:ime,jms:jme)                 )


!inout arrays for NOAHMP only:
 if(.not.allocated(tsnoxy_p)      ) allocate(tsnoxy_p(ims:ime,-2:0,jms:jme)             )
 if(.not.allocated(zsnsoxy_p)     ) allocate(zsnsoxy_p(ims:ime,-2:num_soils,jms:jme)    )
 if(.not.allocated(snicexy_p)     ) allocate(snicexy_p(ims:ime,-2:0,jms:jme)            )
 if(.not.allocated(snliqxy_p)     ) allocate(snliqxy_p(ims:ime,-2:0,jms:jme)            )
 if(.not.allocated(smoiseq_p)     ) allocate(smoiseq_p(ims:ime,1:num_soils,jms:jme)     )
 if(.not.allocated(grainxy_p)     ) allocate(grainxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(gddxy_p)       ) allocate(gddxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(pgsxy_p)       ) allocate(pgsxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qtdrain_p)     ) allocate(qtdrain_p(ims:ime,jms:jme)                 )
 if(.not.allocated(tdfraction_p)  ) allocate(tdfraction_p(ims:ime,jms:jme)              )
 if(.not.allocated(isnowxy_p)     ) allocate(isnowxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(tvxy_p)        ) allocate(tvxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(rs_p)          ) allocate(rs_p(ims:ime,jms:jme)                      )
 if(.not.allocated(tgxy_p)        ) allocate(tgxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(canicexy_p)    ) allocate(canicexy_p(ims:ime,jms:jme)                )
 if(.not.allocated(canliqxy_p)    ) allocate(canliqxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(eahxy_p)       ) allocate(eahxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(tahxy_p)       ) allocate(tahxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(cmxy_p)        ) allocate(cmxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(chxy_p)        ) allocate(chxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(fwetxy_p)      ) allocate(fwetxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(sneqvoxy_p)    ) allocate(sneqvoxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(alboldxy_p)    ) allocate(alboldxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(qsnowxy_p)     ) allocate(qsnowxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(qrainxy_p)     ) allocate(qrainxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(wslakexy_p)    ) allocate(wslakexy_p(ims:ime,jms:jme)                )
 if(.not.allocated(zwtxy_p)       ) allocate(zwtxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(waxy_p)        ) allocate(waxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(wtxy_p)        ) allocate(wtxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(lfmassxy_p)    ) allocate(lfmassxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(rtmassxy_p)    ) allocate(rtmassxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(stmassxy_p)    ) allocate(stmassxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(woodxy_p)      ) allocate(woodxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(stblcpxy_p)    ) allocate(stblcpxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(fastcpxy_p)    ) allocate(fastcpxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(xlaixy_p)      ) allocate(xlaixy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(xsaixy_p)      ) allocate(xsaixy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(taussxy_p)     ) allocate(taussxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(smcwtdxy_p)    ) allocate(smcwtdxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(deeprechxy_p)  ) allocate(deeprechxy_p(ims:ime,jms:jme)              )
 if(.not.allocated(rechxy_p)      ) allocate(rechxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(irnumsi_p)     ) allocate(irnumsi_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irnummi_p)     ) allocate(irnummi_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irnumfi_p)     ) allocate(irnumfi_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irwatsi_p)     ) allocate(irwatsi_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irwatmi_p)     ) allocate(irwatmi_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irwatfi_p)     ) allocate(irwatfi_p(ims:ime,jms:jme)                 )
 if(.not.allocated(ireloss_p)     ) allocate(ireloss_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irsivol_p)     ) allocate(irsivol_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irmivol_p)     ) allocate(irmivol_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irfivol_p)     ) allocate(irfivol_p(ims:ime,jms:jme)                 )
 if(.not.allocated(irrsplh_p)     ) allocate(irrsplh_p(ims:ime,jms:jme)                 )
 if(.not.allocated(gecros_state_p)) allocate(gecros_state_p(ims:ime,60,jms:jme)         )

 if(.not.allocated(acc_ssoilxy_p) ) allocate(acc_ssoilxy_p(ims:ime,jms:jme)             )
 if(.not.allocated(acc_qinsurxy_p)) allocate(acc_qinsurxy_p(ims:ime,jms:jme)            )
 if(.not.allocated(acc_qsevaxy_p) ) allocate(acc_qsevaxy_p(ims:ime,jms:jme)             )
 if(.not.allocated(acc_dwaterxy_p)) allocate(acc_dwaterxy_p(ims:ime,jms:jme)            )
 if(.not.allocated(acc_prcpxy_p)  ) allocate(acc_prcpxy_p(ims:ime,jms:jme)              )
 if(.not.allocated(acc_ecanxy_p)  ) allocate(acc_ecanxy_p(ims:ime,jms:jme)              )
 if(.not.allocated(acc_etranxy_p) ) allocate(acc_etranxy_p(ims:ime,jms:jme)             )
 if(.not.allocated(acc_edirxy_p)  ) allocate(acc_edirxy_p(ims:ime,jms:jme)              )
 if(.not.allocated(acc_etranixy_p)) allocate(acc_etranixy_p(ims:ime,1:num_soils,jms:jme))


!inout arrays shared with the NOAH lsm:
 if(.not.allocated(acsnom_p)      ) allocate(acsnom_p(ims:ime,jms:jme)                  )
 if(.not.allocated(acsnow_p)      ) allocate(acsnow_p(ims:ime,jms:jme)                  )
 if(.not.allocated(canwat_p)      ) allocate(canwat_p(ims:ime,jms:jme)                  )
 if(.not.allocated(grdflx_p)      ) allocate(grdflx_p(ims:ime,jms:jme)                  )
 if(.not.allocated(hfx_p)         ) allocate(hfx_p(ims:ime,jms:jme)                     )
 if(.not.allocated(lh_p)          ) allocate(lh_p(ims:ime,jms:jme)                      )
 if(.not.allocated(qfx_p)         ) allocate(qfx_p(ims:ime,jms:jme)                     )
 if(.not.allocated(qsfc_p)        ) allocate(qsfc_p(ims:ime,jms:jme)                    )
 if(.not.allocated(sfc_albedo_p)  ) allocate(sfc_albedo_p(ims:ime,jms:jme)              )
 if(.not.allocated(sfc_emiss_p)   ) allocate(sfc_emiss_p(ims:ime,jms:jme)               )
 if(.not.allocated(sfcrunoff_p)   ) allocate(sfcrunoff_p(ims:ime,jms:jme)               )
 if(.not.allocated(smstav_p)      ) allocate(smstav_p(ims:ime,jms:jme)                  )
 if(.not.allocated(smstot_p)      ) allocate(smstot_p(ims:ime,jms:jme)                  )
 if(.not.allocated(snow_p)        ) allocate(snow_p(ims:ime,jms:jme)                    )
 if(.not.allocated(snowc_p)       ) allocate(snowc_p(ims:ime,jms:jme)                   )
 if(.not.allocated(snowh_p)       ) allocate(snowh_p(ims:ime,jms:jme)                   )
 if(.not.allocated(tsk_p)         ) allocate(tsk_p(ims:ime,jms:jme)                     )
 if(.not.allocated(udrunoff_p)    ) allocate(udrunoff_p(ims:ime,jms:jme)                )
 if(.not.allocated(z0_p)          ) allocate(z0_p(ims:ime,jms:jme)                      )
 if(.not.allocated(znt_p)         ) allocate(znt_p(ims:ime,jms:jme)                     )
 if(.not.allocated(q2_p)          ) allocate(q2_p(ims:ime,jms:jme)                      )
 if(.not.allocated(t2m_p)         ) allocate(t2m_p(ims:ime,jms:jme)                     )
 if(.not.allocated(th2m_p)        ) allocate(th2m_p(ims:ime,jms:jme)                    )


 if(.not.allocated(sfcheadrt_p)   ) allocate(sfcheadrt_p(ims:ime,jms:jme)               )
 if(.not.allocated(infxsrt_p)     ) allocate(infxsrt_p(ims:ime,jms:jme)                 )
 if(.not.allocated(soldrain_p)    ) allocate(soldrain_p(ims:ime,jms:jme)                )
 if(.not.allocated(qtiledrain_p)  ) allocate(qtiledrain_p(ims:ime,jms:jme)              )
 if(.not.allocated(zwatble2d_p)   ) allocate(zwatble2d_p(ims:ime,jms:jme)               )


!output variables for NOAHMP only:
 if(.not.allocated(t2mvxy_p)      ) allocate(t2mvxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(t2mbxy_p)      ) allocate(t2mbxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(q2mvxy_p)      ) allocate(q2mvxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(q2mbxy_p)      ) allocate(q2mbxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(tradxy_p)      ) allocate(tradxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(neexy_p)       ) allocate(neexy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(gppxy_p)       ) allocate(gppxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(nppxy_p)       ) allocate(nppxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(fvegxy_p)      ) allocate(fvegxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(runsfxy_p)     ) allocate(runsfxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(runsbxy_p)     ) allocate(runsbxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(runsbxy_p)     ) allocate(runsbxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(ecanxy_p)      ) allocate(ecanxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(edirxy_p)      ) allocate(edirxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(etranxy_p)     ) allocate(etranxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(fsaxy_p)       ) allocate(fsaxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(firaxy_p)      ) allocate(firaxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(aparxy_p)      ) allocate(aparxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(psnxy_p)       ) allocate(psnxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(savxy_p)       ) allocate(savxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(sagxy_p)       ) allocate(sagxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(rssunxy_p)     ) allocate(rssunxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(rsshaxy_p)     ) allocate(rsshaxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(bgapxy_p)      ) allocate(bgapxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(wgapxy_p)      ) allocate(wgapxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(tgvxy_p)       ) allocate(tgvxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(tgbxy_p)       ) allocate(tgbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(chvxy_p)       ) allocate(chvxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(chbxy_p)       ) allocate(chbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(shgxy_p)       ) allocate(shgxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(shcxy_p)       ) allocate(shcxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(shbxy_p)       ) allocate(shbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(evgxy_p)       ) allocate(evgxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(evbxy_p)       ) allocate(evbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(ghvxy_p)       ) allocate(ghvxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(ghbxy_p)       ) allocate(ghbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(irgxy_p)       ) allocate(irgxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(ircxy_p)       ) allocate(ircxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(irbxy_p)       ) allocate(irbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(trxy_p)        ) allocate(trxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(evcxy_p)       ) allocate(evcxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(chleafxy_p)    ) allocate(chleafxy_p(ims:ime,jms:jme)                )
 if(.not.allocated(chucxy_p)      ) allocate(chucxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(chv2xy_p)      ) allocate(chv2xy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(chb2xy_p)      ) allocate(chb2xy_p(ims:ime,jms:jme)                  )


!additional output variables for NOAHP only:
 if(.not.allocated(pahxy_p)     ) allocate(pahxy_p(ims:ime,jms:jme)                     )
 if(.not.allocated(pahgxy_p)    ) allocate(pahgxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(pahbxy_p)    ) allocate(pahbxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(pahvxy_p)    ) allocate(pahvxy_p(ims:ime,jms:jme)                    )
 if(.not.allocated(qintsxy_p)   ) allocate(qintsxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qintrxy_p)   ) allocate(qintrxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qdripsxy_p)  ) allocate(qdripsxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qdriprxy_p)  ) allocate(qdriprxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qthrosxy_p)  ) allocate(qthrosxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qthrorxy_p)  ) allocate(qthrorxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qsnsubxy_p)  ) allocate(qsnsubxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qsnfroxy_p)  ) allocate(qsnfroxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qsubcxy_p)   ) allocate(qsubcxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qfrocxy_p)   ) allocate(qfrocxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qevacxy_p)   ) allocate(qevacxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qdewcxy_p)   ) allocate(qdewcxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qfrzcxy_p)   ) allocate(qfrzcxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(qmeltcxy_p)  ) allocate(qmeltcxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(qsnbotxy_p)  ) allocate(qsnbotxy_p(ims:ime,jms:jme)                  )
 if(.not.allocated(pondingxy_p) ) allocate(pondingxy_p(ims:ime,jms:jme)                 )
 if(.not.allocated(fpicexy_p)   ) allocate(fpicexy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(rainlsm_p)   ) allocate(rainlsm_p(ims:ime,jms:jme)                   )
 if(.not.allocated(snowlsm_p)   ) allocate(snowlsm_p(ims:ime,jms:jme)                   )
 if(.not.allocated(forctlsm_p)  ) allocate(forctlsm_p(ims:ime,jms:jme)                  )
 if(.not.allocated(forcqlsm_p)  ) allocate(forcqlsm_p(ims:ime,jms:jme)                  )
 if(.not.allocated(forcplsm_p)  ) allocate(forcplsm_p(ims:ime,jms:jme)                  )
 if(.not.allocated(forczlsm_p)  ) allocate(forczlsm_p(ims:ime,jms:jme)                  )
 if(.not.allocated(forcwlsm_p)  ) allocate(forcwlsm_p(ims:ime,jms:jme)                  )
 if(.not.allocated(eflxbxy_p)   ) allocate(eflxbxy_p(ims:ime,jms:jme)                   )
 if(.not.allocated(soilenergy_p)) allocate(soilenergy_p(ims:ime,jms:jme)                )
 if(.not.allocated(snowenergy_p)) allocate(snowenergy_p(ims:ime,jms:jme)                )
 if(.not.allocated(canhsxy_p)   ) allocate(canhsxy_p(ims:ime,jms:jme)                   )


 end subroutine allocate_lsm_noahmp

!=================================================================================================================
 subroutine deallocate_lsm_noahmp
!=================================================================================================================

!input arrays for soil layer properties:
 if(allocated(dzs_p)         ) deallocate(dzs_p         )
 if(allocated(smcrel_p)      ) deallocate(smcrel_p      ) 
 if(allocated(sh2o_p)        ) deallocate(sh2o_p        )
 if(allocated(smois_p)       ) deallocate(smois_p       )
 if(allocated(tslb_p)        ) deallocate(tslb_p        )


!input arrays:
 if(allocated(dx_p)          ) deallocate(dx_p          )
 if(allocated(area_p)        ) deallocate(area_p        )
 if(allocated(coszr_p)       ) deallocate(coszr_p       )
 if(allocated(xlat_p)        ) deallocate(xlat_p        )
 if(allocated(xlon_p)        ) deallocate(xlon_p        )

 if(allocated(isltyp_p)      ) deallocate(isltyp_p      )
 if(allocated(ivgtyp_p)      ) deallocate(ivgtyp_p      )
 if(allocated(shdmax_p)      ) deallocate(shdmax_p      )
 if(allocated(tmn_p)         ) deallocate(tmn_p         )
 if(allocated(vegfra_p)      ) deallocate(vegfra_p      )
 if(allocated(xice_p)        ) deallocate(xice_p        )
 if(allocated(xland_p)       ) deallocate(xland_p       )

 if(allocated(glw_p)         ) deallocate(glw_p         )
 if(allocated(swddir_p)      ) deallocate(swddir_p      )
 if(allocated(swddif_p)      ) deallocate(swddif_p      )
 if(allocated(swdown_p)      ) deallocate(swdown_p      )
 if(allocated(rainbl_p)      ) deallocate(rainbl_p      )
 if(allocated(sr_p)          ) deallocate(sr_p          )


!arrays for seaice grid-points:
 if(allocated(tsk_sea)       ) deallocate(tsk_sea       )
 if(allocated(tsk_ice)       ) deallocate(tsk_ice       )
 if(allocated(albsi_p)       ) deallocate(albsi_p       )
 if(allocated(icedepth_p)    ) deallocate(icedepth_p    )
 if(allocated(snowsi_p)      ) deallocate(snowsi_p      )
 if(allocated(chs_sea)       ) deallocate(chs_sea       )
 if(allocated(chs2_sea)      ) deallocate(chs2_sea      )
 if(allocated(cqs2_sea)      ) deallocate(cqs2_sea      )
 if(allocated(cpm_sea)       ) deallocate(cpm_sea       )
 if(allocated(hfx_sea)       ) deallocate(hfx_sea       )
 if(allocated(qfx_sea)       ) deallocate(qfx_sea       )
 if(allocated(qgh_sea)       ) deallocate(qgh_sea       )
 if(allocated(qsfc_sea)      ) deallocate(qsfc_sea      )
 if(allocated(lh_sea)        ) deallocate(lh_sea        )
 if(allocated(znt_sea)       ) deallocate(znt_sea       )

 if(allocated(br_p)          ) deallocate(br_p          )
 if(allocated(chs_p)         ) deallocate(chs_p         )
 if(allocated(chs2_p)        ) deallocate(chs2_p        )
 if(allocated(cqs2_p)        ) deallocate(cqs2_p        )
 if(allocated(cpm_p)         ) deallocate(cpm_p         )
 if(allocated(snoalb_p)      ) deallocate(snoalb_p      )


!input arrays for NOAHMP only:
 if(allocated(soilcl1_p)     ) deallocate(soilcl1_p     )
 if(allocated(soilcl2_p)     ) deallocate(soilcl2_p     )
 if(allocated(soilcl3_p)     ) deallocate(soilcl3_p     )
 if(allocated(soilcl4_p)     ) deallocate(soilcl4_p     )
 if(allocated(soilcomp_p)    ) deallocate(soilcomp_p    )

 if(allocated(cropcat_p)     ) deallocate(cropcat_p     )
 if(allocated(planting_p)    ) deallocate(planting_p    )
 if(allocated(harvest_p)     ) deallocate(harvest_p     )
 if(allocated(planting_p)    ) deallocate(planting_p    )
 if(allocated(season_gdd_p)  ) deallocate(season_gdd_p  )

 if(allocated(irfract_p)     ) deallocate(irfract_p     )
 if(allocated(sifract_p)     ) deallocate(sifract_p     )
 if(allocated(mifract_p)     ) deallocate(mifract_p     )
 if(allocated(fifract_p)     ) deallocate(fifract_p     )


!inout arrays for NOAHMP only:
 if(allocated(rs_p)          ) deallocate(rs_p          )
 if(allocated(tsnoxy_p)      ) deallocate(tsnoxy_p      )
 if(allocated(zsnsoxy_p)     ) deallocate(zsnsoxy_p     )
 if(allocated(snicexy_p)     ) deallocate(snicexy_p     )
 if(allocated(snliqxy_p)     ) deallocate(snliqxy_p     )
 if(allocated(smoiseq_p)     ) deallocate(smoiseq_p     )
 if(allocated(grainxy_p)     ) deallocate(grainxy_p     )
 if(allocated(gddxy_p)       ) deallocate(gddxy_p       )
 if(allocated(pgsxy_p)       ) deallocate(pgsxy_p       )
 if(allocated(qtdrain_p)     ) deallocate(qtdrain_p     )
 if(allocated(tdfraction_p)  ) deallocate(tdfraction_p  )
 if(allocated(isnowxy_p)     ) deallocate(isnowxy_p     )
 if(allocated(tvxy_p)        ) deallocate(tvxy_p        )
 if(allocated(rs_p)          ) deallocate(rs_p          )
 if(allocated(tgxy_p)        ) deallocate(tgxy_p        )
 if(allocated(canicexy_p)    ) deallocate(canicexy_p    )
 if(allocated(canliqxy_p)    ) deallocate(canliqxy_p    )
 if(allocated(eahxy_p)       ) deallocate(eahxy_p       )
 if(allocated(tahxy_p)       ) deallocate(tahxy_p       )
 if(allocated(cmxy_p)        ) deallocate(cmxy_p        )
 if(allocated(chxy_p)        ) deallocate(chxy_p        )
 if(allocated(fwetxy_p)      ) deallocate(fwetxy_p      )
 if(allocated(sneqvoxy_p)    ) deallocate(sneqvoxy_p    )
 if(allocated(alboldxy_p)    ) deallocate(alboldxy_p    )
 if(allocated(qsnowxy_p)     ) deallocate(qsnowxy_p     )
 if(allocated(qrainxy_p)     ) deallocate(qrainxy_p     )
 if(allocated(wslakexy_p)    ) deallocate(wslakexy_p    )
 if(allocated(zwtxy_p)       ) deallocate(zwtxy_p       )
 if(allocated(waxy_p)        ) deallocate(waxy_p        )
 if(allocated(wtxy_p)        ) deallocate(wtxy_p        )
 if(allocated(lfmassxy_p)    ) deallocate(lfmassxy_p    )
 if(allocated(rtmassxy_p)    ) deallocate(rtmassxy_p    )
 if(allocated(stmassxy_p)    ) deallocate(stmassxy_p    )
 if(allocated(woodxy_p)      ) deallocate(woodxy_p      )
 if(allocated(stblcpxy_p)    ) deallocate(stblcpxy_p    )
 if(allocated(fastcpxy_p)    ) deallocate(fastcpxy_p    )
 if(allocated(xlaixy_p)      ) deallocate(xlaixy_p      )
 if(allocated(xsaixy_p)      ) deallocate(xsaixy_p      )
 if(allocated(taussxy_p)     ) deallocate(taussxy_p     )
 if(allocated(smcwtdxy_p)    ) deallocate(smcwtdxy_p    )
 if(allocated(deeprechxy_p)  ) deallocate(deeprechxy_p  )
 if(allocated(rechxy_p)      ) deallocate(rechxy_p      )
 if(allocated(irnumsi_p)     ) deallocate(irnumsi_p     )
 if(allocated(irnummi_p)     ) deallocate(irnummi_p     ) 
 if(allocated(irnumfi_p)     ) deallocate(irnumfi_p     )
 if(allocated(irwatsi_p)     ) deallocate(irwatsi_p     )
 if(allocated(irwatmi_p)     ) deallocate(irwatmi_p     )
 if(allocated(irwatfi_p)     ) deallocate(irwatfi_p     )
 if(allocated(ireloss_p)     ) deallocate(ireloss_p     )
 if(allocated(irsivol_p)     ) deallocate(irsivol_p     )
 if(allocated(irmivol_p)     ) deallocate(irmivol_p     )
 if(allocated(irfivol_p)     ) deallocate(irfivol_p     )
 if(allocated(irrsplh_p)     ) deallocate(irrsplh_p     )
 if(allocated(gecros_state_p)) deallocate(gecros_state_p)

 if(allocated(acc_ssoilxy_p) ) deallocate(acc_ssoilxy_p )
 if(allocated(acc_qinsurxy_p)) deallocate(acc_qinsurxy_p)
 if(allocated(acc_qsevaxy_p) ) deallocate(acc_qsevaxy_p )
 if(allocated(acc_dwaterxy_p)) deallocate(acc_dwaterxy_p)
 if(allocated(acc_prcpxy_p)  ) deallocate(acc_prcpxy_p  )
 if(allocated(acc_ecanxy_p)  ) deallocate(acc_ecanxy_p  )
 if(allocated(acc_etranxy_p) ) deallocate(acc_etranxy_p )
 if(allocated(acc_edirxy_p)  ) deallocate(acc_edirxy_p  )
 if(allocated(acc_etranixy_p)) deallocate(acc_etranixy_p)


!inout arrays shared with the NOAH lsm:
 if(allocated(acsnom_p)      ) deallocate(acsnom_p      )
 if(allocated(acsnow_p)      ) deallocate(acsnow_p      )
 if(allocated(canwat_p)      ) deallocate(canwat_p      )
 if(allocated(grdflx_p)      ) deallocate(grdflx_p      )
 if(allocated(hfx_p)         ) deallocate(hfx_p         )
 if(allocated(lh_p)          ) deallocate(lh_p          )
 if(allocated(qfx_p)         ) deallocate(qfx_p         )
 if(allocated(qsfc_p)        ) deallocate(qsfc_p        )
 if(allocated(sfc_albedo_p)  ) deallocate(sfc_albedo_p  )
 if(allocated(sfc_emiss_p)   ) deallocate(sfc_emiss_p   )
 if(allocated(sfcrunoff_p)   ) deallocate(sfcrunoff_p   )
 if(allocated(smstav_p)      ) deallocate(smstav_p      )
 if(allocated(smstot_p)      ) deallocate(smstot_p      )
 if(allocated(snow_p)        ) deallocate(snow_p        )
 if(allocated(snowc_p)       ) deallocate(snowc_p       )
 if(allocated(snowh_p)       ) deallocate(snowh_p       )
 if(allocated(tsk_p)         ) deallocate(tsk_p         )
 if(allocated(udrunoff_p)    ) deallocate(udrunoff_p    )
 if(allocated(z0_p)          ) deallocate(z0_p          )
 if(allocated(znt_p)         ) deallocate(znt_p         )
 if(allocated(q2_p)          ) deallocate(q2_p          )
 if(allocated(t2m_p)         ) deallocate(t2m_p         )
 if(allocated(th2m_p)        ) deallocate(th2m_p        )

 if(allocated(sfcheadrt_p)   ) deallocate(sfcheadrt_p   )
 if(allocated(infxsrt_p)     ) deallocate(infxsrt_p     )
 if(allocated(soldrain_p)    ) deallocate(soldrain_p    )
 if(allocated(qtiledrain_p)  ) deallocate(qtiledrain_p  )
 if(allocated(zwatble2d_p)   ) deallocate(zwatble2d_p   )


!output variables for NOAHMP only:
 if(allocated(t2mvxy_p)      ) deallocate(t2mvxy_p      )
 if(allocated(t2mbxy_p)      ) deallocate(t2mbxy_p      )
 if(allocated(q2mvxy_p)      ) deallocate(q2mvxy_p      )
 if(allocated(q2mbxy_p)      ) deallocate(q2mbxy_p      )
 if(allocated(tradxy_p)      ) deallocate(tradxy_p      )
 if(allocated(neexy_p)       ) deallocate(neexy_p       )
 if(allocated(gppxy_p)       ) deallocate(gppxy_p       )
 if(allocated(nppxy_p)       ) deallocate(nppxy_p       )
 if(allocated(fvegxy_p)      ) deallocate(fvegxy_p      )
 if(allocated(runsfxy_p)     ) deallocate(runsfxy_p     )
 if(allocated(runsbxy_p)     ) deallocate(runsbxy_p     )
 if(allocated(runsbxy_p)     ) deallocate(runsbxy_p     )
 if(allocated(ecanxy_p)      ) deallocate(ecanxy_p      )
 if(allocated(edirxy_p)      ) deallocate(edirxy_p      )
 if(allocated(etranxy_p)     ) deallocate(etranxy_p     )
 if(allocated(fsaxy_p)       ) deallocate(fsaxy_p       )
 if(allocated(firaxy_p)      ) deallocate(firaxy_p      )
 if(allocated(aparxy_p)      ) deallocate(aparxy_p      )
 if(allocated(psnxy_p)       ) deallocate(psnxy_p       )
 if(allocated(savxy_p)       ) deallocate(savxy_p       )
 if(allocated(sagxy_p)       ) deallocate(sagxy_p       )
 if(allocated(rssunxy_p)     ) deallocate(rssunxy_p     )
 if(allocated(rsshaxy_p)     ) deallocate(rsshaxy_p     )
 if(allocated(bgapxy_p)      ) deallocate(bgapxy_p      )
 if(allocated(wgapxy_p)      ) deallocate(wgapxy_p      )
 if(allocated(tgvxy_p)       ) deallocate(tgvxy_p       )
 if(allocated(tgbxy_p)       ) deallocate(tgbxy_p       )
 if(allocated(chvxy_p)       ) deallocate(chvxy_p       )
 if(allocated(chbxy_p)       ) deallocate(chbxy_p       )
 if(allocated(shgxy_p)       ) deallocate(shgxy_p       )
 if(allocated(shcxy_p)       ) deallocate(shcxy_p       )
 if(allocated(shbxy_p)       ) deallocate(shbxy_p       )
 if(allocated(evgxy_p)       ) deallocate(evgxy_p       )
 if(allocated(evbxy_p)       ) deallocate(evbxy_p       )
 if(allocated(ghvxy_p)       ) deallocate(ghvxy_p       )
 if(allocated(ghbxy_p)       ) deallocate(ghbxy_p       )
 if(allocated(irgxy_p)       ) deallocate(irgxy_p       )
 if(allocated(ircxy_p)       ) deallocate(ircxy_p       )
 if(allocated(irbxy_p)       ) deallocate(irbxy_p       )
 if(allocated(trxy_p)        ) deallocate(trxy_p        )
 if(allocated(evcxy_p)       ) deallocate(evcxy_p       )
 if(allocated(chleafxy_p)    ) deallocate(chleafxy_p    )
 if(allocated(chucxy_p)      ) deallocate(chucxy_p      )
 if(allocated(chv2xy_p)      ) deallocate(chv2xy_p      )
 if(allocated(chb2xy_p)      ) deallocate(chb2xy_p      )


!additional output variables for NOAHMP only:
 if(allocated(pahxy_p)       ) deallocate(pahxy_p       )
 if(allocated(pahgxy_p)      ) deallocate(pahgxy_p      )
 if(allocated(pahbxy_p)      ) deallocate(pahbxy_p      )
 if(allocated(pahvxy_p)      ) deallocate(pahvxy_p      )
 if(allocated(qintsxy_p)     ) deallocate(qintsxy_p     )
 if(allocated(qintrxy_p)     ) deallocate(qintrxy_p     )
 if(allocated(qdripsxy_p)    ) deallocate(qdripsxy_p    )
 if(allocated(qdriprxy_p)    ) deallocate(qdriprxy_p    )
 if(allocated(qthrosxy_p)    ) deallocate(qthrosxy_p    )
 if(allocated(qthrorxy_p)    ) deallocate(qthrorxy_p    )
 if(allocated(qsnsubxy_p)    ) deallocate(qsnsubxy_p    )
 if(allocated(qsnfroxy_p)    ) deallocate(qsnfroxy_p    )
 if(allocated(qsubcxy_p)     ) deallocate(qsubcxy_p     )
 if(allocated(qfrocxy_p)     ) deallocate(qfrocxy_p     )
 if(allocated(qevacxy_p)     ) deallocate(qevacxy_p     )
 if(allocated(qdewcxy_p)     ) deallocate(qdewcxy_p     )
 if(allocated(qfrzcxy_p)     ) deallocate(qfrzcxy_p     )
 if(allocated(qmeltcxy_p)    ) deallocate(qmeltcxy_p    )
 if(allocated(qsnbotxy_p)    ) deallocate(qsnbotxy_p    )
 if(allocated(pondingxy_p)   ) deallocate(pondingxy_p   )
 if(allocated(fpicexy_p)     ) deallocate(fpicexy_p     )
 if(allocated(rainlsm_p)     ) deallocate(rainlsm_p     )
 if(allocated(snowlsm_p)     ) deallocate(snowlsm_p     )
 if(allocated(forctlsm_p)    ) deallocate(forctlsm_p    )
 if(allocated(forcqlsm_p)    ) deallocate(forcqlsm_p    )
 if(allocated(forcplsm_p)    ) deallocate(forcplsm_p    )
 if(allocated(forczlsm_p)    ) deallocate(forczlsm_p    )
 if(allocated(forcwlsm_p)    ) deallocate(forcwlsm_p    )
 if(allocated(eflxbxy_p)     ) deallocate(eflxbxy_p     )
 if(allocated(soilenergy_p)  ) deallocate(soilenergy_p  )
 if(allocated(snowenergy_p)  ) deallocate(snowenergy_p  )
 if(allocated(canhsxy_p)     ) deallocate(canhsxy_p     )


 end subroutine deallocate_lsm_noahmp

!=================================================================================================================
 subroutine lsm_noahmp_from_MPAS(configs,mesh,diag_physics,diag_physics_noahmp,sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite


!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: sfc_input


!local input pointers:
 logical,pointer:: fractional_seaice

 character(len=StrKIND),pointer:: microp_scheme,    &
                                  convection_scheme
 character(len=StrKIND),pointer:: mminlu

 integer,dimension(:),pointer:: isltyp,ivgtyp
 real(kind=RKIND),pointer:: len_disp
 real(kind=RKIND),dimension(:),pointer:: coszr
 real(kind=RKIND),dimension(:),pointer:: areaCell,latCell,lonCell,meshDensity
 real(kind=RKIND),dimension(:),pointer:: shdmax,tmn,vegfra,vegmax,xice,xland
 real(kind=RKIND),dimension(:),pointer:: glw,gsw,swddif,swddir
 real(kind=RKIND),dimension(:),pointer:: raincv,rainncv,sr
 real(kind=RKIND),dimension(:),pointer:: snoalb
 real(kind=RKIND),dimension(:,:),pointer:: dzs


!local input pointers for NOAHMP only:
 integer,dimension(:),pointer:: cropcat

 real(kind=RKIND),dimension(:),pointer  :: soilcl1,soilcl2,soilcl3,soilcl4
 real(kind=RKIND),dimension(:,:),pointer:: soilcomp


!local inout pointers shared with the NOAH lsm:
 real(kind=RKIND),dimension(:),pointer:: &
    skintemp,chs,chs2,cqs2,cpm,hfx,qfx,qgh,lai,lh,grdflx,smstav,smstot,sfcrunoff,udrunoff, &
    sfc_albedo,snow,snowc,snowh,canwat,acsnom,acsnow,sfc_emiss,qsfc,br,z0,znt,q2,t2m,th2m

 real(kind=RKIND),dimension(:,:),pointer:: &
    smois,sh2o,tslb 


!local inout NOAH-MP only pointers:
 integer,dimension(:),pointer:: &
    irnumsi,irnummi,irnumfi,isnowxy

 real(kind=RKIND),dimension(:),pointer:: &
    grainxy,gddxy,pgsxy

 real(kind=RKIND),dimension(:),pointer:: &
    qtdrain,td_fraction

 real(kind=RKIND),dimension(:),pointer:: &
    sfcheadrt,infxsrt,soldrain,qtiledrain,zwatble2d

 real(kind=RKIND),dimension(:),pointer:: &
    irwatsi,irwatmi,irwatfi,ireloss,irsivol,irmivol,irfivol,irrsplh

 real(kind=RKIND),dimension(:),pointer:: &
    tvxy,tgxy,canicexy,canliqxy,eahxy,tahxy,cmxy,chxy,fwetxy,sneqvoxy,alboldxy,qsnowxy,qrainxy,  &
    wslakexy,zwtxy,waxy,wtxy,lfmassxy,rtmassxy,stmassxy,woodxy,stblcpxy,fastcpxy,xsaixy,taussxy, &
    smcwtdxy,deeprechxy,rechxy,rs

 real(kind=RKIND),dimension(:),pointer:: &
    acc_ssoil,acc_qinsur,acc_qseva,acc_dwater,acc_prcp,acc_ecan,acc_etran,acc_edir

 real(kind=RKIND),dimension(:,:),pointer:: &
    smoiseq,snliqxy,snicexy,zsnsoxy,tsnoxy

 real(kind=RKIND),dimension(:,:),pointer:: &
    gecros_state

 real(kind=RKIND),dimension(:,:),pointer:: &
    acc_etrani

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine lsm_noahmp_from_MPAS:')


!local input pointers:
 call mpas_pool_get_config(configs,'config_frac_seaice'      ,fractional_seaice)
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,microp_scheme    )
 call mpas_pool_get_config(configs,'config_len_disp'         ,len_disp         )

 call mpas_pool_get_array(mesh,'areaCell'   ,areaCell   )
 call mpas_pool_get_array(mesh,'latCell'    ,latCell    )
 call mpas_pool_get_array(mesh,'lonCell'    ,lonCell    )
 call mpas_pool_get_array(mesh,'meshDensity',meshDensity)

 call mpas_pool_get_array(sfc_input,'mminlu',mminlu)
 call mpas_pool_get_array(sfc_input,'isltyp',isltyp)
 call mpas_pool_get_array(sfc_input,'ivgtyp',ivgtyp)
 call mpas_pool_get_array(sfc_input,'dzs'   ,dzs   )
 call mpas_pool_get_array(sfc_input,'shdmax',shdmax)
 call mpas_pool_get_array(sfc_input,'snoalb',snoalb)
 call mpas_pool_get_array(sfc_input,'tmn'   ,tmn   )
 call mpas_pool_get_array(sfc_input,'vegfra',vegfra)
 call mpas_pool_get_array(sfc_input,'xice'  ,xice  )
 call mpas_pool_get_array(sfc_input,'xland' ,xland )

 call mpas_pool_get_array(diag_physics,'coszr'  ,coszr  )
 call mpas_pool_get_array(diag_physics,'glw'    ,glw    )
 call mpas_pool_get_array(diag_physics,'gsw'    ,gsw    )
 call mpas_pool_get_array(diag_physics,'swddif' ,swddif )
 call mpas_pool_get_array(diag_physics,'swddir' ,swddir )
 call mpas_pool_get_array(diag_physics,'raincv' ,raincv )
 call mpas_pool_get_array(diag_physics,'rainncv',rainncv)
 call mpas_pool_get_array(diag_physics,'sr'     ,sr     )


!local input pointers for NOAHMP only:
 call mpas_pool_get_array(mesh,'soilcl1' ,soilcl1 )
 call mpas_pool_get_array(mesh,'soilcl2' ,soilcl2 )
 call mpas_pool_get_array(mesh,'soilcl3' ,soilcl3 )
 call mpas_pool_get_array(mesh,'soilcl4' ,soilcl4 )
 call mpas_pool_get_array(mesh,'soilcomp',soilcomp)
 
 call mpas_pool_get_array(diag_physics_noahmp,'cropcat',cropcat)


!inout pointers shared with the NOAH lsm:
 call mpas_pool_get_array(diag_physics,'chs'       ,chs       )
 call mpas_pool_get_array(diag_physics,'chs2'      ,chs2      )
 call mpas_pool_get_array(diag_physics,'cqs2'      ,cqs2      )
 call mpas_pool_get_array(diag_physics,'cpm'       ,cpm       )
 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx       )
 call mpas_pool_get_array(diag_physics,'qfx'       ,qfx       )
 call mpas_pool_get_array(diag_physics,'qgh'       ,qgh       )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh        )
 call mpas_pool_get_array(diag_physics,'grdflx'    ,grdflx    )
 call mpas_pool_get_array(diag_physics,'smstav'    ,smstav    )
 call mpas_pool_get_array(diag_physics,'smstot'    ,smstot    )
 call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
 call mpas_pool_get_array(diag_physics,'udrunoff'  ,udrunoff  )
 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'lai'       ,lai       )
 call mpas_pool_get_array(diag_physics,'canwat'    ,canwat    )
 call mpas_pool_get_array(diag_physics,'acsnom'    ,acsnom    )
 call mpas_pool_get_array(diag_physics,'acsnow'    ,acsnow    )
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
 call mpas_pool_get_array(diag_physics,'br'        ,br        )
 call mpas_pool_get_array(diag_physics,'z0'        ,z0        )
 call mpas_pool_get_array(diag_physics,'znt'       ,znt       )
 call mpas_pool_get_array(diag_physics,'q2'        ,q2        )
 call mpas_pool_get_array(diag_physics,'t2m'       ,t2m       )
 call mpas_pool_get_array(diag_physics,'th2m'      ,th2m      )

 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'snowc'   ,snowc   )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'sh2o'    ,sh2o    )
 call mpas_pool_get_array(sfc_input,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input,'tslb'    ,tslb    )


!local inout NOAH-MP only pointers:
 call mpas_pool_get_array(diag_physics_noahmp,'grainxy'     , grainxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'gddxy'       , gddxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'pgsxy'       , pgsxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'gecros_state',gecros_state)
 call mpas_pool_get_array(diag_physics_noahmp,'qtdrain'     ,qtdrain     )
 call mpas_pool_get_array(diag_physics_noahmp,'td_fraction' ,td_fraction )
 call mpas_pool_get_array(diag_physics_noahmp,'sfcheadrt'   ,sfcheadrt   )
 call mpas_pool_get_array(diag_physics_noahmp,'infxsrt'     ,infxsrt     )
 call mpas_pool_get_array(diag_physics_noahmp,'soldrain'    ,soldrain    )
 call mpas_pool_get_array(diag_physics_noahmp,'qtiledrain'  ,qtiledrain  )
 call mpas_pool_get_array(diag_physics_noahmp,'zwatble2d'   ,zwatble2d   )

 call mpas_pool_get_array(diag_physics_noahmp,'irnumsi'     ,irnumsi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irnummi'     ,irnummi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irnumfi'     ,irnumfi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatsi'     ,irwatsi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatmi'     ,irwatmi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatfi'     ,irwatfi     )
 call mpas_pool_get_array(diag_physics_noahmp,'ireloss'     ,ireloss     )
 call mpas_pool_get_array(diag_physics_noahmp,'irsivol'     ,irsivol     )
 call mpas_pool_get_array(diag_physics_noahmp,'irmivol'     ,irmivol     )
 call mpas_pool_get_array(diag_physics_noahmp,'irfivol'     ,irfivol     )
 call mpas_pool_get_array(diag_physics_noahmp,'irrsplh'     ,irrsplh     )
 call mpas_pool_get_array(diag_physics_noahmp,'isnowxy'     ,isnowxy     )

 call mpas_pool_get_array(diag_physics_noahmp,'rs'          ,rs          )
 call mpas_pool_get_array(diag_physics_noahmp,'tvxy'        ,tvxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'tgxy'        ,tgxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'canicexy'    ,canicexy    )
 call mpas_pool_get_array(diag_physics_noahmp,'canliqxy'    ,canliqxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'eahxy'       ,eahxy       )
 call mpas_pool_get_array(diag_physics_noahmp,'tahxy'       ,tahxy       )
 call mpas_pool_get_array(diag_physics_noahmp,'cmxy'        ,cmxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'chxy'        ,chxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'fwetxy'      ,fwetxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'sneqvoxy'    ,sneqvoxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'alboldxy'    ,alboldxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'qsnowxy'     ,qsnowxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'qrainxy'     ,qrainxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'wslakexy'    ,wslakexy    )
 call mpas_pool_get_array(diag_physics_noahmp,'zwtxy'       ,zwtxy       )
 call mpas_pool_get_array(diag_physics_noahmp,'waxy'        ,waxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'wtxy'        ,wtxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'tsnoxy'      ,tsnoxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'zsnsoxy'     ,zsnsoxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'snicexy'     ,snicexy     )
 call mpas_pool_get_array(diag_physics_noahmp,'snliqxy'     ,snliqxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'lfmassxy'    ,lfmassxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'rtmassxy'    ,rtmassxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'stmassxy'    ,stmassxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'woodxy'      ,woodxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'stblcpxy'    ,stblcpxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'fastcpxy'    ,fastcpxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'xsaixy'      ,xsaixy      )
 call mpas_pool_get_array(diag_physics_noahmp,'taussxy'     ,taussxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'smoiseq'     ,smoiseq     )
 call mpas_pool_get_array(diag_physics_noahmp,'smcwtdxy'    ,smcwtdxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'deeprechxy'  ,deeprechxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'rechxy'      ,rechxy      )

 call mpas_pool_get_array(diag_physics_noahmp,'acc_ssoil'   ,acc_ssoil   )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_qinsur'  ,acc_qinsur  )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_qseva'   ,acc_qseva   )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_dwater'  ,acc_dwater  )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_prcp'    ,acc_prcp    )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_ecan'    ,acc_ecan    )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_etran'   ,acc_etran   )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_edir'    ,acc_edir    )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_etrani'  ,acc_etrani  )

 llanduse = mminlu


!--- in Registry.xml, dzs is a function of nCells. In NOAH-MP, dzs is independent of cell locations:
 do n = 1,num_soils
    dzs_p(n) = dzs(n,its)
 enddo


!--- local input arrays:
 do j = jts,jte
    do i = its,ite
       isltyp_p(i,j) = isltyp(i)
       ivgtyp_p(i,j) = ivgtyp(i)
       shdmax_p(i,j) = shdmax(i)
       tmn_p(i,j)    = tmn(i)
       vegfra_p(i,j) = vegfra(i)
       xice_p(i,j)   = xice(i)
       xland_p(i,j)  = xland(i)

       coszr_p(i,j)  = coszr(i)
       dx_p(i,j)     = len_disp / meshDensity(i)**0.25
       area_p(i,j)   = areaCell(i)
       xlat_p(i,j)   = latCell(i)/degrad
       xlon_p(i,j)   = lonCell(i)/degrad

       glw_p(i,j)    = glw(i)
       swddif_p(i,j) = swddif(i)
       swddir_p(i,j) = swddif(i)
       swdown_p(i,j) = gsw(i) / (1._RKIND - sfc_albedo(i))
    enddo
 enddo

 do j = jts,jte
    do i = its,ite
       sr_p(i,j)     = 0._RKIND
       rainbl_p(i,j) = 0._RKIND
    enddo
 enddo
 if(microp_scheme .ne. 'off') then
    do j = jts,jte
       do i = its,ite
          sr_p(i,j) = sr(i)
          rainbl_p(i,j) = rainbl_p(i,j) + rainncv(i)
       enddo
    enddo
 endif
 if(convection_scheme .ne. 'off') then
    do j = jts,jte
       do i = its,ite
          rainbl_p(i,j) = rainbl_p(i,j) + raincv(i)
       enddo
    enddo
 endif


!--- local input arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       !soil composition:
       soilcl1_p(i,j)    = soilcl1(i)
       soilcl2_p(i,j)    = soilcl2(i)
       soilcl3_p(i,j)    = soilcl3(i)
       soilcl4_p(i,j)    = soilcl4(i)
       do n = 1,num_soilcs
          soilcomp_p(i,n,j) = soilcomp(n,i)
       enddo

       !crop model:
       cropcat_p(i,j)    = cropcat(i)
       planting_p(i,j)   = 0._RKIND
       harvest_p(i,j)    = 0._RKIND
       season_gdd_p(i,j) = 0._RKIND

       !irrigation model:
       irfract_p(i,j)    = 0._RKIND
       sifract_p(i,j)    = 0._RKIND
       mifract_p(i,j)    = 0._RKIND
       fifract_p(i,j)    = 0._RKIND
    enddo
 enddo


!--- inout arrays shared with NOAH lsm:
 do j = jts,jte
    do i = its,ite
       acsnom_p(i,j)     = acsnom(i)
       acsnow_p(i,j)     = acsnow(i)
       canwat_p(i,j)     = canwat(i)
       chs2_p(i,j)       = chs2(i)
       cqs2_p(i,j)       = cqs2(i)
       grdflx_p(i,j)     = grdflx(i)
       hfx_p(i,j)        = hfx(i)
       lh_p(i,j)         = lh(i)
       qfx_p(i,j)        = qfx(i)
       qsfc_p(i,j)       = qsfc(i)
       br_p(i,j)         = br(i)
       sfc_albedo_p(i,j) = sfc_albedo(i)
       sfc_emiss_p(i,j)  = sfc_emiss(i)
       sfcrunoff_p(i,j)  = sfcrunoff(i)
       smstav_p(i,j)     = smstav(i)
       smstot_p(i,j)     = smstot(i)
       snow_p(i,j)       = snow(i)
       snowc_p(i,j)      = snowc(i)
       snowh_p(i,j)      = snowh(i)
       tsk_p(i,j)        = skintemp(i)
       udrunoff_p(i,j)   = udrunoff(i)
       z0_p(i,j)         = z0(i)
       znt_p(i,j)        = znt(i)
       q2_p(i,j)         = q2(i)
       t2m_p(i,j)        = t2m(i)
       th2m_p(i,j)       = th2m(i)
!for WRF_HYDRO only:
       sfcheadrt_p(i,j)  = sfcheadrt(i)
       infxsrt_p(i,j)    = infxsrt(i)
       soldrain_p(i,j)   = soldrain(i)
       qtiledrain_p(i,j) = qtiledrain(i)
       zwatble2d_p(i,j)  = zwatble2d(i)
    enddo

    do n = 1,num_soils
       do i = its,ite
          sh2o_p(i,n,j)   = sh2o(n,i)
          smois_p(i,n,j)  = smois(n,i)
          tslb_p(i,n,j)   = tslb(n,i)
       enddo
    enddo
 enddo

 if(fractional_seaice) then
    do j = jts,jte
       do i = its,ite
          !modify the surface albedo and surface emissivity, and surface temperatures over sea-ice points:
          if(xice(i).ge.xice_threshold .and. xice(i).le.1._RKIND) then
             sfc_albedo_p(i,j) = (sfc_albedo(i) - 0.08_RKIND*(1._RKIND-xice(i))) / xice(i)
             sfc_emiss_p(i,j)  = (sfc_emiss(i) - 0.98_RKIND*(1._RKIND-xice(i))) / xice(i)
          else
             sfc_albedo_p(i,j) = sfc_albedo(i)
             sfc_emiss_p(i,j)  = sfc_emiss(i)
          endif

          !initialize tsk_sea and tsk_ice:
          tsk_sea(i,j) = tsk_p(i,j)
          tsk_ice(i,j) = tsk_p(i,j)

          !initialize extra variables (note that those variables are nor modified/used in NOAHMP):
          chs_p(i,j)    = chs(i)
          cpm_p(i,j)    = cpm(i)
          qgh_p(i,j)    = qgh(i)
          snoalb_p(i,j) = snoalb(i)
       enddo
    enddo

    !calculate sea-surface and sea-ice temperatures over sea-ice grid cells:
    call correct_tsk_over_seaice(ims,ime,jms,jme,its,ite,jts,jte,xice_threshold,xice_p, &
                                 tsk_p,tsk_sea,tsk_ice)

    do j = jts,jte
       do i = its,ite
          tsk_p(i,j) = tsk_ice(i,j)
       enddo
    enddo

    !initialize the surface albedo, the surface albedo over snow-covered seaice, and the
    !seaice thickness.
    do j = jts,jte
       do i = its,ite
          albsi_p(i,j)    = seaice_albedo_default
          icedepth_p(i,j) = seaice_thickness_default
          snowsi_p(i,j)   = seaice_snowdepth_min
       enddo
    enddo
 endif


!--- inout arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       irnumsi_p(i,j) = irnumsi(i)
       irnummi_p(i,j) = irnummi(i)
       irnumfi_p(i,j) = irnumfi(i)
       irwatsi_p(i,j) = irwatsi(i)
       irwatmi_p(i,j) = irwatmi(i)
       irwatfi_p(i,j) = irwatfi(i)
       ireloss_p(i,j) = ireloss(i)
       irsivol_p(i,j) = irsivol(i)
       irmivol_p(i,j) = irmivol(i)
       irfivol_p(i,j) = irfivol(i)
       irrsplh_p(i,j) = irrsplh(i)
    enddo

    do n = -2,0
       do i = its,ite
          tsnoxy_p(i,n,j)  = tsnoxy(n+3,i)
          snliqxy_p(i,n,j) = snliqxy(n+3,i)
          snicexy_p(i,n,j) = snicexy(n+3,i)
          zsnsoxy_p(i,n,j) = zsnsoxy(n+3,i)
       enddo
    enddo
    do n = 1,num_soils
       do i = its,ite
          smoiseq_p(i,n,j) = smoiseq(n,i)
          !finish the initialization of zsnsoxy_p:
          zsnsoxy_p(i,n,j) = zsnsoxy(n+num_soils-1,i)
       enddo
    enddo

    do n = 1,60
       do i = its,ite
          gecros_state_p(i,n,j) = gecros_state(n,i)
       enddo
    enddo

    do i = its,ite
       acc_ssoilxy_p(i,j)  = acc_ssoil(i)
       acc_qinsurxy_p(i,j) = acc_qinsur(i)
       acc_qsevaxy_p(i,j)  = acc_qseva(i)
       acc_dwaterxy_p(i,j) = acc_dwater(i)
       acc_prcpxy_p(i,j)   = acc_prcp(i)
       acc_ecanxy_p(i,j)   = acc_ecan(i)
       acc_etranxy_p(i,j)  = acc_etran(i)
       acc_edirxy_p(i,j)   = acc_edir(i)
    enddo
    do n = 1,num_soils
       do i = its,ite
          acc_etranixy_p(i,n,j) = acc_etrani(n,i)
       enddo
    enddo
 enddo


!--- inout arrays for NOAHMP ONLY:
 do j = jts,jte
    do i = its,ite
       grainxy_p(i,j)      = grainxy(i)
       gddxy_p(i,j)        = gddxy(i)
       pgsxy_p(i,j)        = pgsxy(i)
       qtdrain_p(i,j)      = qtdrain(i)
       tdfraction_p(i,j)   = td_fraction(i)
    enddo

    do i = its,ite
       xlaixy_p(i,j)       = lai(i)
       isnowxy_p(i,j)      = isnowxy(i)
       tvxy_p(i,j)         = tvxy(i)
       rs_p(i,j)           = rs(i)
       tgxy_p(i,j)         = tgxy(i)
       canicexy_p(i,j)     = canicexy(i)
       canliqxy_p(i,j)     = canliqxy(i)
       eahxy_p(i,j)        = eahxy(i)
       tahxy_p(i,j)        = tahxy(i)
       cmxy_p(i,j)         = cmxy(i)
       chxy_p(i,j)         = chxy(i)
       fwetxy_p(i,j)       = fwetxy(i)
       sneqvoxy_p(i,j)     = sneqvoxy(i)
       alboldxy_p(i,j)     = alboldxy(i)
       qsnowxy_p(i,j)      = qsnowxy(i)
       qrainxy_p(i,j)      = qrainxy(i)
       wslakexy_p(i,j)     = wslakexy(i)
       zwtxy_p(i,j)        = zwtxy(i)
       waxy_p(i,j)         = waxy(i)
       wtxy_p(i,j)         = wtxy(i)
       lfmassxy_p(i,j)     = lfmassxy(i)
       rtmassxy_p(i,j)     = rtmassxy(i)
       stmassxy_p(i,j)     = stmassxy(i)
       woodxy_p(i,j)       = woodxy(i)
       stblcpxy_p(i,j)     = stblcpxy(i)
       fastcpxy_p(i,j)     = fastcpxy(i)
       xsaixy_p(i,j)       = xsaixy(i)
       taussxy_p(i,j)      = taussxy(i)
       smcwtdxy_p(i,j)     = smcwtdxy(i)
       deeprechxy_p(i,j)   = deeprechxy(i)
       rechxy_p(i,j)       = rechxy(i)
    enddo
 enddo


!--- output arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       t2mvxy_p(i,j)     = 0._RKIND
       t2mbxy_p(i,j)     = 0._RKIND
       q2mvxy_p(i,j)     = 0._RKIND
       q2mbxy_p(i,j)     = 0._RKIND
       tradxy_p(i,j)     = 0._RKIND
       neexy_p(i,j)      = 0._RKIND
       gppxy_p(i,j)      = 0._RKIND
       nppxy_p(i,j)      = 0._RKIND
       fvegxy_p(i,j)     = 0._RKIND
       runsfxy_p(i,j)    = 0._RKIND
       runsbxy_p(i,j)    = 0._RKIND
       ecanxy_p(i,j)     = 0._RKIND
       edirxy_p(i,j)     = 0._RKIND
       etranxy_p(i,j)    = 0._RKIND
       fsaxy_p(i,j)      = 0._RKIND
       firaxy_p(i,j)     = 0._RKIND
       aparxy_p(i,j)     = 0._RKIND
       psnxy_p(i,j)      = 0._RKIND
       savxy_p(i,j)      = 0._RKIND
       sagxy_p(i,j)      = 0._RKIND
       rssunxy_p(i,j)    = 0._RKIND
       rsshaxy_p(i,j)    = 0._RKIND 
       bgapxy_p(i,j)     = 0._RKIND
       wgapxy_p(i,j)     = 0._RKIND
       tgvxy_p(i,j)      = 0._RKIND
       tgbxy_p(i,j)      = 0._RKIND
       chvxy_p(i,j)      = 0._RKIND
       chbxy_p(i,j)      = 0._RKIND
       shgxy_p(i,j)      = 0._RKIND
       shcxy_p(i,j)      = 0._RKIND
       shbxy_p(i,j)      = 0._RKIND
       evgxy_p(i,j)      = 0._RKIND
       evbxy_p(i,j)      = 0._RKIND
       ghvxy_p(i,j)      = 0._RKIND
       ghbxy_p(i,j)      = 0._RKIND
       irgxy_p(i,j)      = 0._RKIND
       ircxy_p(i,j)      = 0._RKIND
       irbxy_p(i,j)      = 0._RKIND
       trxy_p(i,j)       = 0._RKIND
       evcxy_p(i,j)      = 0._RKIND
       chleafxy_p(i,j)   = 0._RKIND
       chucxy_p(i,j)     = 0._RKIND
       chv2xy_p(i,j)     = 0._RKIND
       chb2xy_p(i,j)     = 0._RKIND
    enddo
 enddo


!--- additional output arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
    pahxy_p(i,j)         = 0._RKIND
    pahgxy_p(i,j)        = 0._RKIND
    pahbxy_p(i,j)        = 0._RKIND
    pahvxy_p(i,j)        = 0._RKIND
    qintsxy_p(i,j)       = 0._RKIND
    qintrxy_p(i,j)       = 0._RKIND
    qdripsxy_p(i,j)      = 0._RKIND
    qdriprxy_p(i,j)      = 0._RKIND
    qthrosxy_p(i,j)      = 0._RKIND
    qthrorxy_p(i,j)      = 0._RKIND
    qsnsubxy_p(i,j)      = 0._RKIND
    qsnfroxy_p(i,j)      = 0._RKIND
    qsubcxy_p(i,j)       = 0._RKIND
    qfrocxy_p(i,j)       = 0._RKIND
    qevacxy_p(i,j)       = 0._RKIND
    qdewcxy_p(i,j)       = 0._RKIND
    qfrzcxy_p(i,j)       = 0._RKIND
    qmeltcxy_p(i,j)      = 0._RKIND
    qsnbotxy_p(i,j)      = 0._RKIND
    pondingxy_p(i,j)     = 0._RKIND
    fpicexy_p(i,j)       = 0._RKIND
    rainlsm_p(i,j)       = 0._RKIND
    snowlsm_p(i,j)       = 0._RKIND
    forctlsm_p(i,j)      = 0._RKIND
    forcqlsm_p(i,j)      = 0._RKIND
    forcplsm_p(i,j)      = 0._RKIND
    forczlsm_p(i,j)      = 0._RKIND
    forcwlsm_p(i,j)      = 0._RKIND
    eflxbxy_p(i,j)       = 0._RKIND
    soilenergy_p(i,j)    = 0._RKIND
    snowenergy_p(i,j)    = 0._RKIND
    canhsxy_p(i,j)       = 0._RKIND
    enddo
 enddo


 call mpas_log_write('--- end subroutine lsm_noahmp_from_MPAS:')

 end subroutine lsm_noahmp_from_MPAS
 
!=================================================================================================================
 subroutine lsm_noahmp_to_MPAS(configs,mesh,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite


!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: output_noahmp
 type(mpas_pool_type),intent(inout):: sfc_input


!local input pointers:
 logical,pointer:: fractional_seaice


!local inout pointers shared with the NOAH lsm:
 real(kind=RKIND),dimension(:),pointer:: &
    skintemp,chs,chs2,cqs2,cpm,hfx,qfx,qgh,lai,lh,grdflx,smstav,smstot,sfcrunoff,udrunoff, &
    sfc_albedo,snow,snowc,snowh,canwat,acsnom,acsnow,sfc_emiss,qsfc,z0,znt,q2,t2m,th2m

 real(kind=RKIND),dimension(:,:),pointer:: &
    smois,sh2o,tslb


!local inout NOAH-MP only pointers:
 integer,dimension(:),pointer:: &
    irnumsi,irnummi,irnumfi,isnowxy

 real(kind=RKIND),dimension(:),pointer:: &
    grainxy,gddxy,pgsxy

 real(kind=RKIND),dimension(:),pointer:: &
    qtdrain,td_fraction

 real(kind=RKIND),dimension(:),pointer:: &
    sfcheadrt,infxsrt,soldrain,qtiledrain,zwatble2d

 real(kind=RKIND),dimension(:),pointer:: &
    irwatsi,irwatmi,irwatfi,ireloss,irsivol,irmivol,irfivol,irrsplh

 real(kind=RKIND),dimension(:),pointer:: &
    tvxy,tgxy,canicexy,canliqxy,eahxy,tahxy,cmxy,chxy,fwetxy,sneqvoxy,alboldxy,qsnowxy,qrainxy,  &
    wslakexy,zwtxy,waxy,wtxy,lfmassxy,rtmassxy,stmassxy,woodxy,stblcpxy,fastcpxy,xsaixy,taussxy, &
    smcwtdxy,deeprechxy,rechxy,rs

 real(kind=RKIND),dimension(:),pointer:: &
    acc_ssoil,acc_qinsur,acc_qseva,acc_dwater,acc_prcp,acc_ecan,acc_etran,acc_edir

 real(kind=RKIND),dimension(:,:),pointer:: &
    smoiseq,snliqxy,snicexy,zsnsoxy,tsnoxy

 real(kind=RKIND),dimension(:,:),pointer:: &
    gecros_state

 real(kind=RKIND),dimension(:,:),pointer:: &
    acc_etrani


!local out NOAHMP only pointers:
 real(kind=RKIND),dimension(:),pointer:: &
    q2mvxy,q2mbxy,tradxy,neexy,gppxy,nppxy,fvegxy,runsfxy,runsbxy,ecanxy,edirxy,etranxy,fsaxy, &
    firaxy,aparxy,psnxy,savxy,sagxy,rssunxy,rsshaxy,bgapxy,wgapxy,tgvxy,tgbxy,chvxy,chbxy,     &
    shgxy,shcxy,shbxy,evgxy,evbxy,ghvxy,ghbxy,irgxy,ircxy,irbxy,trxy,evcxy,chleafxy,chucxy,    &
    chv2xy,chb2xy

 real(kind=RKIND),dimension(:),pointer:: &
    t2mvxy,t2mbxy


!additional local out NOAHMP only pointers:
 real(kind=RKIND),dimension(:),pointer:: &
    pahxy,pahgxy,pahbxy,pahvxy,qintsxy,qintrxy,qdripsxy,qdriprxy,qthrosxy,qthrorxy,qsnsubxy,     &
    qsnfroxy,qsubcxy,qfrocxy,qevacxy,qdewcxy,qfrzcxy,qmeltcxy,qsnbotxy,pondingxy,fpicexy,        &
    rainlsm,snowlsm,forctlsm,forcqlsm,forcplsm,forczlsm,forcwlsm,eflxbxy,soilenergy,snowenergy,  &
    canhsxy

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine lsm_noahmp_to_MPAS:')

 call mpas_pool_get_config(configs,'config_frac_seaice',fractional_seaice)

!inout pointers shared with the NOAH lsm:
 call mpas_pool_get_array(diag_physics,'chs'       ,chs       )
 call mpas_pool_get_array(diag_physics,'chs2'      ,chs2      )
 call mpas_pool_get_array(diag_physics,'cqs2'      ,cqs2      )
 call mpas_pool_get_array(diag_physics,'cpm'       ,cpm       )
 call mpas_pool_get_array(diag_physics,'hfx'       ,hfx       )
 call mpas_pool_get_array(diag_physics,'qfx'       ,qfx       )
 call mpas_pool_get_array(diag_physics,'qgh'       ,qgh       )
 call mpas_pool_get_array(diag_physics,'lh'        ,lh        )
 call mpas_pool_get_array(diag_physics,'grdflx'    ,grdflx    )
 call mpas_pool_get_array(diag_physics,'smstav'    ,smstav    )
 call mpas_pool_get_array(diag_physics,'smstot'    ,smstot    )
 call mpas_pool_get_array(diag_physics,'sfcrunoff' ,sfcrunoff )
 call mpas_pool_get_array(diag_physics,'udrunoff'  ,udrunoff  )
 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'lai'       ,lai       )
 call mpas_pool_get_array(diag_physics,'canwat'    ,canwat    )
 call mpas_pool_get_array(diag_physics,'acsnom'    ,acsnom    )
 call mpas_pool_get_array(diag_physics,'acsnow'    ,acsnow    )
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'qsfc'      ,qsfc      )
 call mpas_pool_get_array(diag_physics,'z0'        ,z0        )
 call mpas_pool_get_array(diag_physics,'znt'       ,znt       )
 call mpas_pool_get_array(diag_physics,'q2'        ,q2        )
 call mpas_pool_get_array(diag_physics,'t2m'       ,t2m       )
 call mpas_pool_get_array(diag_physics,'th2m'      ,th2m      )

 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'snowc'   ,snowc   )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'sh2o'    ,sh2o    )
 call mpas_pool_get_array(sfc_input,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input,'tslb'    ,tslb    )


!local inout NOAHMP only pointers:
 call mpas_pool_get_array(diag_physics_noahmp,'grainxy'     , grainxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'gddxy'       , gddxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'pgsxy'       , pgsxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'gecros_state',gecros_state)
 call mpas_pool_get_array(diag_physics_noahmp,'qtdrain'     ,qtdrain     )
 call mpas_pool_get_array(diag_physics_noahmp,'td_fraction' ,td_fraction )
 call mpas_pool_get_array(diag_physics_noahmp,'sfcheadrt'   ,sfcheadrt   )
 call mpas_pool_get_array(diag_physics_noahmp,'infxsrt'     ,infxsrt     )
 call mpas_pool_get_array(diag_physics_noahmp,'soldrain'    ,soldrain    )
 call mpas_pool_get_array(diag_physics_noahmp,'qtiledrain'  ,qtiledrain  )
 call mpas_pool_get_array(diag_physics_noahmp,'zwatble2d'   ,zwatble2d   )

 call mpas_pool_get_array(diag_physics_noahmp,'irnumsi'     ,irnumsi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irnummi'     ,irnummi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irnumfi'     ,irnumfi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatsi'     ,irwatsi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatmi'     ,irwatmi     )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatfi'     ,irwatfi     )
 call mpas_pool_get_array(diag_physics_noahmp,'ireloss'     ,ireloss     )
 call mpas_pool_get_array(diag_physics_noahmp,'irsivol'     ,irsivol     )
 call mpas_pool_get_array(diag_physics_noahmp,'irmivol'     ,irmivol     )
 call mpas_pool_get_array(diag_physics_noahmp,'irfivol'     ,irfivol     )
 call mpas_pool_get_array(diag_physics_noahmp,'irrsplh'     ,irrsplh     )
 call mpas_pool_get_array(diag_physics_noahmp,'isnowxy'     ,isnowxy     )

 call mpas_pool_get_array(diag_physics_noahmp,'rs'          ,rs          )
 call mpas_pool_get_array(diag_physics_noahmp,'tvxy'        ,tvxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'tgxy'        ,tgxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'canicexy'    ,canicexy    )
 call mpas_pool_get_array(diag_physics_noahmp,'canliqxy'    ,canliqxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'eahxy'       ,eahxy       )
 call mpas_pool_get_array(diag_physics_noahmp,'tahxy'       ,tahxy       )
 call mpas_pool_get_array(diag_physics_noahmp,'cmxy'        ,cmxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'chxy'        ,chxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'fwetxy'      ,fwetxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'sneqvoxy'    ,sneqvoxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'alboldxy'    ,alboldxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'qsnowxy'     ,qsnowxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'qrainxy'     , qrainxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'wslakexy'    ,wslakexy    )
 call mpas_pool_get_array(diag_physics_noahmp,'zwtxy'       ,zwtxy       )
 call mpas_pool_get_array(diag_physics_noahmp,'waxy'        ,waxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'wtxy'        ,wtxy        )
 call mpas_pool_get_array(diag_physics_noahmp,'tsnoxy'      ,tsnoxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'zsnsoxy'     ,zsnsoxy     ) 
 call mpas_pool_get_array(diag_physics_noahmp,'snicexy'     ,snicexy     )
 call mpas_pool_get_array(diag_physics_noahmp,'snliqxy'     ,snliqxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'lfmassxy'    ,lfmassxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'rtmassxy'    ,rtmassxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'stmassxy'    ,stmassxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'woodxy'      ,woodxy      )
 call mpas_pool_get_array(diag_physics_noahmp,'stblcpxy'    ,stblcpxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'fastcpxy'    ,fastcpxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'xsaixy'      ,xsaixy      )
 call mpas_pool_get_array(diag_physics_noahmp,'taussxy'     ,taussxy     )
 call mpas_pool_get_array(diag_physics_noahmp,'smoiseq'     ,smoiseq     )
 call mpas_pool_get_array(diag_physics_noahmp,'smcwtdxy'    ,smcwtdxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'deeprechxy'  ,deeprechxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'rechxy'      ,rechxy      )

 call mpas_pool_get_array(diag_physics_noahmp,'acc_ssoil'   ,acc_ssoil   )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_qinsur'  ,acc_qinsur  )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_qseva'   ,acc_qseva   )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_dwater'  ,acc_dwater  )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_prcp'    ,acc_prcp    )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_ecan'    ,acc_ecan    )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_etran'   ,acc_etran   )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_edir'    ,acc_edir    )
 call mpas_pool_get_array(diag_physics_noahmp,'acc_etrani'  ,acc_etrani  )


!local out NOAHMP only pointers:
 call mpas_pool_get_array(output_noahmp,'t2mvxy'    ,t2mvxy    )
 call mpas_pool_get_array(output_noahmp,'t2mbxy'    ,t2mbxy    )
 call mpas_pool_get_array(output_noahmp,'q2mvxy'    ,q2mvxy    )
 call mpas_pool_get_array(output_noahmp,'q2mbxy'    ,q2mbxy    )
 call mpas_pool_get_array(output_noahmp,'tradxy'    ,tradxy    )
 call mpas_pool_get_array(output_noahmp,'neexy'     ,neexy     )
 call mpas_pool_get_array(output_noahmp,'gppxy'     ,gppxy     )
 call mpas_pool_get_array(output_noahmp,'nppxy'     ,nppxy     )
 call mpas_pool_get_array(output_noahmp,'fvegxy'    ,fvegxy    )
 call mpas_pool_get_array(output_noahmp,'runsfxy'   ,runsfxy   )
 call mpas_pool_get_array(output_noahmp,'runsbxy'   ,runsbxy   )
 call mpas_pool_get_array(output_noahmp,'ecanxy'    ,ecanxy    )
 call mpas_pool_get_array(output_noahmp,'edirxy'    ,edirxy    )
 call mpas_pool_get_array(output_noahmp,'etranxy'   ,etranxy   )
 call mpas_pool_get_array(output_noahmp,'fsaxy'     ,fsaxy     )
 call mpas_pool_get_array(output_noahmp,'firaxy'    ,firaxy    )
 call mpas_pool_get_array(output_noahmp,'aparxy'    ,aparxy    )
 call mpas_pool_get_array(output_noahmp,'psnxy'     ,psnxy     )
 call mpas_pool_get_array(output_noahmp,'savxy'     ,savxy     )
 call mpas_pool_get_array(output_noahmp,'sagxy'     ,sagxy     )
 call mpas_pool_get_array(output_noahmp,'rssunxy'   ,rssunxy   )
 call mpas_pool_get_array(output_noahmp,'rsshaxy'   ,rsshaxy   )
 call mpas_pool_get_array(output_noahmp,'bgapxy'    ,bgapxy    )
 call mpas_pool_get_array(output_noahmp,'wgapxy'    ,wgapxy    )
 call mpas_pool_get_array(output_noahmp,'tgvxy'     ,tgvxy     )
 call mpas_pool_get_array(output_noahmp,'tgbxy'     ,tgbxy     )
 call mpas_pool_get_array(output_noahmp,'chvxy'     ,chvxy     )
 call mpas_pool_get_array(output_noahmp,'chbxy'     ,chbxy     )
 call mpas_pool_get_array(output_noahmp,'shgxy'     ,shgxy     )
 call mpas_pool_get_array(output_noahmp,'shcxy'     ,shcxy     )
 call mpas_pool_get_array(output_noahmp,'shbxy'     ,shbxy     )
 call mpas_pool_get_array(output_noahmp,'evgxy'     ,evgxy     )
 call mpas_pool_get_array(output_noahmp,'evbxy'     ,evbxy     )
 call mpas_pool_get_array(output_noahmp,'ghvxy'     ,ghvxy     )
 call mpas_pool_get_array(output_noahmp,'ghbxy'     ,ghbxy     )
 call mpas_pool_get_array(output_noahmp,'irgxy'     ,irgxy     )
 call mpas_pool_get_array(output_noahmp,'ircxy'     ,ircxy     )
 call mpas_pool_get_array(output_noahmp,'irbxy'     ,irbxy     )
 call mpas_pool_get_array(output_noahmp,'trxy'      ,trxy      )
 call mpas_pool_get_array(output_noahmp,'evcxy'     ,evcxy     )
 call mpas_pool_get_array(output_noahmp,'chleafxy'  ,chleafxy  )
 call mpas_pool_get_array(output_noahmp,'chucxy'    ,chucxy    )
 call mpas_pool_get_array(output_noahmp,'chv2xy'    ,chv2xy    )
 call mpas_pool_get_array(output_noahmp,'chb2xy'    ,chb2xy    )


!additional local out NOAHMP only pointers:
 call mpas_pool_get_array(output_noahmp,'pahxy'     ,pahxy     )
 call mpas_pool_get_array(output_noahmp,'pahgxy'    ,pahgxy    )
 call mpas_pool_get_array(output_noahmp,'pahbxy'    ,pahbxy    )
 call mpas_pool_get_array(output_noahmp,'pahvxy'    ,pahvxy    )
 call mpas_pool_get_array(output_noahmp,'qintsxy'   ,qintsxy   )
 call mpas_pool_get_array(output_noahmp,'qintrxy'   ,qintrxy   )
 call mpas_pool_get_array(output_noahmp,'qdripsxy'  ,qdripsxy  )
 call mpas_pool_get_array(output_noahmp,'qdriprxy'  ,qdriprxy  )
 call mpas_pool_get_array(output_noahmp,'qthrosxy'  ,qthrosxy  )
 call mpas_pool_get_array(output_noahmp,'qthrorxy'  ,qthrorxy  )
 call mpas_pool_get_array(output_noahmp,'qsnsubxy'  ,qsnsubxy  )
 call mpas_pool_get_array(output_noahmp,'qsnfroxy'  ,qsnfroxy  )
 call mpas_pool_get_array(output_noahmp,'qsubcxy'   ,qsubcxy   )
 call mpas_pool_get_array(output_noahmp,'qfrocxy'   ,qfrocxy   )
 call mpas_pool_get_array(output_noahmp,'qevacxy'   ,qevacxy   )
 call mpas_pool_get_array(output_noahmp,'qdewcxy'   ,qdewcxy   )
 call mpas_pool_get_array(output_noahmp,'qfrzcxy'   ,qfrzcxy   )
 call mpas_pool_get_array(output_noahmp,'qmeltcxy'  ,qmeltcxy  )
 call mpas_pool_get_array(output_noahmp,'qsnbotxy'  ,qsnbotxy  )
 call mpas_pool_get_array(output_noahmp,'pondingxy' ,pondingxy )
 call mpas_pool_get_array(output_noahmp,'fpicexy'   ,fpicexy   )
 call mpas_pool_get_array(output_noahmp,'rainlsm'   ,rainlsm   )
 call mpas_pool_get_array(output_noahmp,'snowlsm'   ,snowlsm   )
 call mpas_pool_get_array(output_noahmp,'forctlsm'  ,forctlsm  )
 call mpas_pool_get_array(output_noahmp,'forcqlsm'  ,forcqlsm  )
 call mpas_pool_get_array(output_noahmp,'forcplsm'  ,forcplsm  )
 call mpas_pool_get_array(output_noahmp,'forczlsm'  ,forczlsm  )
 call mpas_pool_get_array(output_noahmp,'forcwlsm'  ,forcwlsm  )
 call mpas_pool_get_array(output_noahmp,'eflxbxy'   ,eflxbxy   )
 call mpas_pool_get_array(output_noahmp,'soilenergy',soilenergy)
 call mpas_pool_get_array(output_noahmp,'snowenergy',snowenergy)
 call mpas_pool_get_array(output_noahmp,'canhsxy'   ,canhsxy   )


!inout arrays shared with the NOAH lsm:
 do j = jts,jte
    do i = its,ite
       acsnom(i)     = acsnom_p(i,j)
       acsnow(i)     = acsnow_p(i,j)
       canwat(i)     = canwat_p(i,j)
       grdflx(i)     = grdflx_p(i,j)
       hfx(i)        = hfx_p(i,j)
       lh(i)         = lh_p(i,j)
       qfx(i)        = qfx_p(i,j)
       qsfc(i)       = qsfc_p(i,j)
       rs(i)         = rs_p(i,j)
       sfc_albedo(i) = sfc_albedo_p(i,j)
       sfc_emiss(i)  = sfc_emiss_p(i,j)
       sfcrunoff(i)  = sfcrunoff_p(i,j)
       smstav(i)     = smstav_p(i,j)
       smstot(i)     = smstot_p(i,j)
       snow(i)       = snow_p(i,j)
       snowc(i)      = snowc_p(i,j)
       snowh(i)      = snowh_p(i,j)
       skintemp(i)   = tsk_p(i,j)
       udrunoff(i)   = udrunoff_p(i,j)
       z0(i)         = z0_p(i,j)
       znt(i)        = znt_p(i,j)
       q2(i)         = q2_p(i,j)
       t2m(i)        = t2m_p(i,j)
       th2m(i)       = th2m_p(i,j)
!for WRF_HYDRO only:
       sfcheadrt(i)    = sfcheadrt_p(i,j)
       infxsrt(i)      = infxsrt_p(i,j)
       soldrain(i)     = soldrain_p(i,j)
       qtiledrain(i)   = qtiledrain_p(i,j)
       zwatble2d(i)    = zwatble2d_p(i,j)
    enddo

    do n = 1,num_soils
       do i = its,ite
          sh2o(n,i)  = sh2o_p(i,n,j)
          smois(n,i) = smois_p(i,n,j)
          tslb(n,i)  = tslb_p(i,n,j)
       enddo
    enddo
 enddo

 if(fractional_seaice) then
    do j = jts,jte
       do i = its,ite
          if(xice_p(i,j).ge.xice_threshold .and. xice_p(i,j).le.1._RKIND) then
             chs(i)        = xice_p(i,j)*chs_p(i,j)  + (1._RKIND-xice_p(i,j))*chs_sea(i,j)
             chs2(i)       = xice_p(i,j)*chs2_p(i,j) + (1._RKIND-xice_p(i,j))*chs2_sea(i,j)
             cqs2(i)       = xice_p(i,j)*cqs2_p(i,j) + (1._RKIND-xice_p(i,j))*cqs2_sea(i,j)
             cpm(i)        = xice_p(i,j)*cpm_p(i,j)  + (1._RKIND-xice_p(i,j))*cpm_sea(i,j)
             hfx(i)        = xice_p(i,j)*hfx_p(i,j)  + (1._RKIND-xice_p(i,j))*hfx_sea(i,j)
             lh(i)         = xice_p(i,j)*lh_p(i,j)   + (1._RKIND-xice_p(i,j))*lh_sea(i,j)
             qfx(i)        = xice_p(i,j)*qfx_p(i,j)  + (1._RKIND-xice_p(i,j))*qfx_sea(i,j)
             qgh(i)        = xice_p(i,j)*qgh_p(i,j)  + (1._RKIND-xice_p(i,j))*qgh_sea(i,j)
             qsfc(i)       = xice_p(i,j)*qsfc_p(i,j) + (1._RKIND-xice_p(i,j))*qsfc_sea(i,j)
             znt(i)        = xice_p(i,j)*znt_p(i,j)  + (1._RKIND-xice_p(i,j))*znt_sea(i,j)
             skintemp(i)   = xice_p(i,j)*tsk_p(i,j)  + (1._RKIND-xice_p(i,j))*tsk_sea(i,j)
             sfc_albedo(i) = xice_p(i,j)*sfc_albedo_p(i,j) + (1._RKIND-xice_p(i,j))*0.08_RKIND
             sfc_emiss(i)  = xice_p(i,j)*sfc_emiss_p(i,j)  + (1._RKIND-xice_p(i,j))*0.98_RKIND
          endif
       enddo
    enddo
 endif

!--- inout arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       irnumsi(i) = irnumsi_p(i,j)
       irnummi(i) = irnummi_p(i,j)
       irnumfi(i) = irnumfi_p(i,j)
       irwatsi(i) = irwatsi_p(i,j)
       irwatmi(i) = irwatmi_p(i,j)
       irwatfi(i) = irwatfi_p(i,j)
       ireloss(i) = ireloss_p(i,j)
       irsivol(i) = irsivol_p(i,j)
       irmivol(i) = irmivol_p(i,j)
       irfivol(i) = irfivol_p(i,j)
       irrsplh(i) = irrsplh_p(i,j)
    enddo

    do n = -2,0
       do i = its,ite
          tsnoxy(n+3,i)  = tsnoxy_p(i,n,j)
          snliqxy(n+3,i) = snliqxy_p(i,n,j)
          snicexy(n+3,i) = snicexy_p(i,n,j)
          zsnsoxy(n+3,i) = zsnsoxy_p(i,n,j)
       enddo
    enddo
    do n = 1,num_soils
       do i = its,ite
          smoiseq(n,i) = smoiseq_p(i,n,j)
          !finish the initialization of zsnsoxy_p:
          zsnsoxy(n+num_soils-1,i) = zsnsoxy_p(i,n,j)
       enddo
    enddo

    do n = 1,60
       do i = its,ite
          gecros_state(n,i) = gecros_state_p(i,n,j)
       enddo
    enddo

    do i = its,ite
       acc_ssoil(i)  = acc_ssoilxy_p(i,j)
       acc_qinsur(i) = acc_qinsurxy_p(i,j)
       acc_qseva(i)  = acc_qsevaxy_p(i,j)
       acc_dwater(i) = acc_dwaterxy_p(i,j)
       acc_prcp(i)   = acc_prcpxy_p(i,j)
       acc_ecan(i)   = acc_ecanxy_p(i,j)
       acc_etran(i)  = acc_etranxy_p(i,j)
       acc_edir(i)   = acc_edirxy_p(i,j)
    enddo
    do n = 1,num_soils
       do i = its,ite
          acc_etrani(n,i) = acc_etranixy_p(i,n,j)
       enddo
    enddo
 enddo


!--- inout arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       grainxy(i)      = grainxy_p(i,j)
       gddxy(i)        = gddxy_p(i,j)
       pgsxy(i)        = pgsxy_p(i,j)
       qtdrain(i)      = qtdrain_p(i,j)
       td_fraction(i)  = tdfraction_p(i,j)
    enddo

    do i = its,ite
       lai(i)          = xlaixy_p(i,j)
       isnowxy(i)      = isnowxy_p(i,j)
       rs(i)           = rs_p(i,j)
       tvxy(i)         = tvxy_p(i,j)
       tgxy(i)         = tgxy_p(i,j)
       canicexy(i)     = canicexy_p(i,j)
       canliqxy(i)     = canliqxy_p(i,j)
       eahxy(i)        = eahxy_p(i,j)
       tahxy(i)        = tahxy_p(i,j)
       cmxy(i)         = cmxy_p(i,j)
       chxy(i)         = chxy_p(i,j)
       fwetxy(i)       = fwetxy_p(i,j)
       sneqvoxy(i)     = sneqvoxy_p(i,j)
       alboldxy(i)     = alboldxy_p(i,j)
       qsnowxy(i)      = qsnowxy_p(i,j)
       qrainxy(i)      = qrainxy_p(i,j)
       wslakexy(i)     = wslakexy_p(i,j)
       zwtxy(i)        = zwtxy_p(i,j)
       waxy(i)         = waxy_p(i,j)
       wtxy(i)         = wtxy_p(i,j)
       lfmassxy(i)     = lfmassxy_p(i,j)
       rtmassxy(i)     = rtmassxy_p(i,j)
       stmassxy(i)     = stmassxy_p(i,j)
       woodxy(i)       = woodxy_p(i,j)
       stblcpxy(i)     = stblcpxy_p(i,j)
       fastcpxy(i)     = fastcpxy_p(i,j)
       xsaixy(i)       = xsaixy_p(i,j)
       taussxy(i)      = taussxy_p(i,j)
       smcwtdxy(i)     = smcwtdxy_p(i,j)
       deeprechxy(i)   = deeprechxy_p(i,j)
       rechxy(i)       = rechxy_p(i,j)
    enddo
 enddo


!--- output arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       t2mvxy(i)   = t2mvxy_p(i,j)
       t2mbxy(i)   = t2mbxy_p(i,j)
       q2mvxy(i)   = q2mvxy_p(i,j)
       q2mbxy(i)   = q2mbxy_p(i,j)
       tradxy(i)   = tradxy_p(i,j)
       neexy(i)    = neexy_p(i,j)
       gppxy(i)    = gppxy_p(i,j)
       nppxy(i)    = nppxy_p(i,j)
       fvegxy(i)   = fvegxy_p(i,j)
       runsfxy(i)  = runsfxy_p(i,j)
       runsbxy(i)  = runsbxy_p(i,j)
       ecanxy(i)   = ecanxy_p(i,j)
       edirxy(i)   = edirxy_p(i,j)
       etranxy(i)   = etranxy_p(i,j)
       fsaxy(i)    = fsaxy_p(i,j)
       firaxy(i)   = firaxy_p(i,j)
       aparxy(i)   = aparxy_p(i,j)
       psnxy(i)    = psnxy_p(i,j)
       savxy(i)    = savxy_p(i,j)
       sagxy(i)    = sagxy_p(i,j)
       rssunxy(i)  = rssunxy_p(i,j)
       rsshaxy(i)  = rsshaxy_p(i,j)
       bgapxy(i)   = bgapxy_p(i,j)
       wgapxy(i)   = wgapxy_p(i,j)
       tgvxy(i)    = tgvxy_p(i,j)
       tgbxy(i)    = tgbxy_p(i,j)
       chvxy(i)    = chvxy_p(i,j)
       chbxy(i)    = chbxy_p(i,j)
       shgxy(i)    = shgxy_p(i,j)
       shcxy(i)    = shcxy_p(i,j)
       shbxy(i)    = shbxy_p(i,j)
       evgxy(i)    = evgxy_p(i,j)
       evbxy(i)    = evbxy_p(i,j)
       ghvxy(i)    = ghvxy_p(i,j)
       ghbxy(i)    = ghbxy_p(i,j)
       irgxy(i)    = irgxy_p(i,j)
       ircxy(i)    = ircxy_p(i,j)
       irbxy(i)    = irbxy_p(i,j)
       trxy(i)     = trxy_p(i,j)
       evcxy(i)    = evcxy_p(i,j)
       chleafxy(i) = chleafxy_p(i,j)
       chucxy(i)   = chucxy_p(i,j)
       chv2xy(i)   = chv2xy_p(i,j)
       chb2xy(i)   = chb2xy_p(i,j)
    enddo
 enddo


!--- additional output arrays for NOAHMP only:
 do j = jts,jte
    do i = its,ite
       pahxy(i)      = pahxy_p(i,j)
       pahgxy(i)     = pahgxy_p(i,j)
       pahbxy(i)     = pahbxy_p(i,j)
       pahvxy(i)     = pahvxy_p(i,j)
       qintsxy(i)    = qintsxy_p(i,j)
       qintrxy(i)    = qintrxy_p(i,j)
       qdripsxy(i)   = qdripsxy_p(i,j)
       qdriprxy(i)   = qdriprxy_p(i,j)
       qthrosxy(i)   = qthrosxy_p(i,j)
       qthrorxy(i)   = qthrorxy_p(i,j)
       qsnsubxy(i)   = qsnsubxy_p(i,j)
       qsnfroxy(i)   = qsnfroxy_p(i,j)
       qsubcxy(i)    = qsubcxy_p(i,j)
       qfrocxy(i)    = qfrocxy_p(i,j)
       qevacxy(i)    = qevacxy_p(i,j)
       qdewcxy(i)    = qdewcxy_p(i,j)
       qfrzcxy(i)    = qfrzcxy_p(i,j)
       qmeltcxy(i)   = qmeltcxy_p(i,j)
       qsnbotxy(i)   = qsnbotxy_p(i,j)
       pondingxy(i)  = pondingxy_p(i,j)
       fpicexy(i)    = fpicexy_p(i,j)
       rainlsm(i)    = rainlsm_p(i,j)
       snowlsm(i)    = snowlsm_p(i,j)
       forctlsm(i)   = forctlsm_p(i,j)
       forcqlsm(i)   = forcqlsm_p(i,j)
       forcplsm(i)   = forcplsm_p(i,j)
       forczlsm(i)   = forczlsm_p(i,j)
       forcwlsm(i)   = forcwlsm_p(i,j)
       eflxbxy(i)    = eflxbxy_p(i,j)
       soilenergy(i) = soilenergy_p(i,j)
       snowenergy(i) = snowenergy_p(i,j)
       canhsxy(i)    = canhsxy_p(i,j)
    enddo
 enddo


 call mpas_log_write('--- end subroutine lsm_noahmp_to_MPAS:')

 end subroutine lsm_noahmp_to_MPAS
 
!=================================================================================================================
 subroutine init_lsm_noahmp(dminfo,mesh,configs,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)
!=================================================================================================================

!input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type):: mesh
 type(mpas_pool_type),intent(in):: configs

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: output_noahmp
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 character(len=StrKIND),pointer:: lsm_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_lsm_scheme',lsm_scheme)

 lsm_select: select case (trim(lsm_scheme))

    case ("sf_noahmp")
       call noahmp_init_forMPAS(mesh,configs,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)

    case default
 
 end select lsm_select

 end subroutine init_lsm_noahmp

!=================================================================================================================
 subroutine driver_lsm_noahmp(itimestep,configs,mesh,diag_physics,diag_physics_noahmp,output_noahmp, &
                              sfc_input,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs

 integer,intent(in):: its,ite
 integer,intent(in):: itimestep

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: output_noahmp
 type(mpas_pool_type),intent(inout):: sfc_input

 real(kind=RKIND),dimension(its:ite):: tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10

!local variables, arrays, and pointers:
 integer,parameter:: &
    lcz_1 = 0, lcz_2 = 0, lcz_3 = 0, lcz_4  = 0, lcz_5  = 0, lcz_6 = 0, &
    lcz_7 = 0, lcz_8 = 0, lcz_9 = 0, lcz_10 = 0, lcz_11 = 0

 integer,pointer:: isice,iswater

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_lsm_noahmp:')

 call mpas_pool_get_array(sfc_input,'isice'  ,isice  )
 call mpas_pool_get_array(sfc_input,'iswater',iswater)
 call mpas_log_write('--- isice   = $i',intArgs=(/isice/))
 call mpas_log_write('--- isurban = $i',intArgs=(/isurban/))
 call mpas_log_write('--- iswater = $i',intArgs=(/iswater/))


!copy MPAS arrays to local arrays:
 call lsm_noahmp_from_MPAS(configs,mesh,diag_physics,diag_physics_noahmp,sfc_input,its,ite)


 call noahmplsm( &
!--- input variables:
      itimestep   = itimestep    , nsoil        = num_soils     , yr           = year          , & !in
      julian      = curr_julday  , dt           = dt_pbl        , dzs          = dzs_p         , & !in
      isltyp      = isltyp_p     , ivgtyp       = ivgtyp_p      , tmn          = tmn_p         , & !in
      vegfra      = vegfra_p     , vegmax       = shdmax_p      , xland        = xland_p       , & !in
      xice        = xice_p       , coszin       = coszr_p       , xlat         = xlat_p        , & !in 
      xlong       = xlon_p       , dxin         = dx_p          , dz8w         = dz_p          , & !in
      t3d         = t_p          , qv3d         = qv_p          , u_phy        = u_p           , & !in
      v_phy       = v_p          , p8w3d        = pres2_hyd_p   , glw          = glw_p         , & !in
      swddif      = swddif_p     , swddir       = swddir_p      , swdown       = swdown_p      , & !in
      precip_in   = rainbl_p     , sr           = sr_p                                         , & !in
      iz0tlnd          = iz0tlnd          ,                                                      & !in
      xice_thres       = xice_threshold   ,                                                      & !in
      sf_urban_physics = sf_urban_physics ,                                                      & !in
      soiltstep        = soiltstep        ,                                                      & !in
      soilcl1     = soilcl1_p    , soilcl2      = soilcl2_p     , soilcl3      = soilcl3_p     , & !in NOAHMP
      soilcl4     = soilcl4_p    , soilcomp     = soilcomp_p    , cropcat      = cropcat_p     , & !in NOAHMP
      planting    = planting_p   , harvest      = harvest_p     , season_gdd   = season_gdd_p  , & !in NOAHMP
      irfract     = irfract_p    , sifract      = sifract_p     , mifract      = mifract_p     , & !in NOAHMP
      fifract     = fifract_p    ,                                                               & !in NOAHMP
      idveg       = idveg        , iopt_crs     = iopt_crs      , iopt_btr     = iopt_btr      , & !in NOAHMP
      iopt_run    = iopt_run     , iopt_sfc     = iopt_sfc      , iopt_frz     = iopt_frz      , & !in NOAHMP
      iopt_inf    = iopt_inf     , iopt_rad     = iopt_rad      , iopt_alb     = iopt_alb      , & !in NOAHMP
      iopt_snf    = iopt_snf     , iopt_tbot    = iopt_tbot     , iopt_stc     = iopt_stc      , & !in NOAHMP
      iopt_gla    = iopt_gla     , iopt_rsf     = iopt_rsf      , iopt_soil    = iopt_soil     , & !in NOAHMP
      iopt_pedo   = iopt_pedo    , iopt_crop    = iopt_crop     , iopt_irr     = iopt_irr      , & !in NOAHMP
      iopt_irrm   = iopt_irrm    , iopt_infdv   = iopt_infdv    , iopt_tdrn    = iopt_tdrn     , & !in NOAHMP
!--- inout variables shared with the NOAH lsm:
      tsk         = tsk_p        , hfx          = hfx_p         , qfx          = qfx_p         , & !in/out
      lh          = lh_p         , grdflx       = grdflx_p      , smstav       = smstav_p      , & !in/out
      smstot      = smstot_p     , sfcrunoff    = sfcrunoff_p   , udrunoff     = udrunoff_p    , & !in/out
      albedo      = sfc_albedo_p , snowc        = snowc_p       , smois        = smois_p       , & !in/out
      sh2o        = sh2o_p       , tslb         = tslb_p        , snow         = snow_p        , & !in/out
      snowh       = snowh_p      , canwat       = canwat_p      , acsnom       = acsnom_p      , & !in/out
      acsnow      = acsnow_p     , emiss        = sfc_emiss_p   , qsfc         = qsfc_p        , & !in/out
      z0          = z0_p         , znt          = znt_p         , llanduse     = llanduse      , & !in/out
!--- inout variables:
      isnowxy     = isnowxy_p    , tvxy         = tvxy_p        , tgxy         = tgxy_p        , & !in/out NOAHMP
      canicexy    = canicexy_p   , canliqxy     = canliqxy_p    , eahxy        = eahxy_p       , & !in/out NOAHMP
      tahxy       = tahxy_p      , cmxy         = cmxy_p        , chxy         = chxy_p        , & !in/out NOAHMP
      fwetxy      = fwetxy_p     , sneqvoxy     = sneqvoxy_p    , alboldxy     = alboldxy_p    , & !in/out NOAHMP
      qsnowxy     = qsnowxy_p    , qrainxy      = qrainxy_p     , wslakexy     = wslakexy_p    , & !in/out NOAHMP
      zwtxy       = zwtxy_p      , waxy         = waxy_p        , wtxy         = wtxy_p        , & !in/out NOAHMP
      tsnoxy      = tsnoxy_p     , zsnsoxy      = zsnsoxy_p     , snicexy      = snicexy_p     , & !in/out NOAHMP
      snliqxy     = snliqxy_p    , lfmassxy     = lfmassxy_p    , rtmassxy     = rtmassxy_p    , & !in/out NOAHMP
      stmassxy    = stmassxy_p   , woodxy       = woodxy_p      , stblcpxy     = stblcpxy_p    , & !in/out NOAHMP
      fastcpxy    = fastcpxy_p   , xlaixy       = xlaixy_p      , xsaixy       = xsaixy_p      , & !in/out NOAHMP
      taussxy     = taussxy_p    , smoiseq      = smoiseq_p     , smcwtdxy     = smcwtdxy_p    , & !in/out NOAHMP
      deeprechxy  = deeprechxy_p , rechxy       = rechxy_p      , grainxy      = grainxy_p     , & !in/out NOAHMP
      gddxy       = gddxy_p      , pgsxy        = pgsxy_p       , gecros_state = gecros_state_p, & !in/out NOAHMP
      qtdrain     = qtdrain_p    , td_fraction  = tdfraction_p  , irnumsi      = irnumsi_p     , & !in/out NOAHMP
      irnummi     = irnummi_p    , irnumfi      = irnumfi_p     , irwatsi      = irwatsi_p     , & !in/out NOAHMP
      irwatmi     = irwatmi_p    , irwatfi      = irwatfi_p     , ireloss      = ireloss_p     , & !in/out NOAHMP
      irsivol     = irsivol_p    , irmivol      = irmivol_p     , irfivol      = irfivol_p     , & !in/out NOAHMP
      irrsplh     = irrsplh_p    , acc_ssoilxy  = acc_ssoilxy_p , acc_qinsurxy = acc_qinsurxy_p, & !in/out NOAHMP
      acc_qsevaxy = acc_qsevaxy_p, acc_etranixy = acc_etranixy_p, acc_dwaterxy = acc_dwaterxy_p, & !in/out NOAHMP
      acc_prcpxy  = acc_prcpxy_p , acc_ecanxy   = acc_ecanxy_p  , acc_etranxy  = acc_etranxy_p , & !in/out NOAHMP
      acc_edirxy  = acc_edirxy_p ,                                                               & !in/out NOAHMP
!--- out variables:
      t2mvxy      = t2mvxy_p     , t2mbxy       = t2mbxy_p      , q2mvxy       = q2mvxy_p      , & !out NOAHMP
      q2mbxy      = q2mbxy_p     , tradxy       = tradxy_p      , neexy        = neexy_p       , & !out NOAHMP
      gppxy       = gppxy_p      , nppxy        = nppxy_p       , fvegxy       = fvegxy_p      , & !out NOAHMP
      runsfxy     = runsfxy_p    , runsbxy      = runsbxy_p     , ecanxy       = ecanxy_p      , & !out NOAHMP
      edirxy      = edirxy_p     , etranxy      = etranxy_p     , fsaxy        = fsaxy_p       , & !out NOAHMP
      firaxy      = firaxy_p     , aparxy       = aparxy_p      , psnxy        = psnxy_p       , & !out NOAHMP
      savxy       = savxy_p      , sagxy        = sagxy_p       , rssunxy      = rssunxy_p     , & !out NOAHMP
      rsshaxy     = rsshaxy_p    , bgapxy       = bgapxy_p      , wgapxy       = wgapxy_p      , & !out NOAHMP
      tgvxy       = tgvxy_p      , tgbxy        = tgbxy_p       , chvxy        = chvxy_p       , & !out NOAHMP
      chbxy       = chbxy_p      , shgxy        = shgxy_p       , shcxy        = shcxy_p       , & !out NOAHMP
      shbxy       = shbxy_p      , evgxy        = evgxy_p       , evbxy        = evbxy_p       , & !out NOAHMP
      ghvxy       = ghvxy_p      , ghbxy        = ghbxy_p       , irgxy        = irgxy_p       , & !out NOAHMP
      ircxy       = ircxy_p      , irbxy        = irbxy_p       , trxy         = trxy_p        , & !out NOAHMP
      evcxy       = evcxy_p      , chleafxy     = chleafxy_p    , chucxy       = chucxy_p      , & !out NOAHMP
      chv2xy      = chv2xy_p     , chb2xy       = chb2xy_p      , rs           = rs_p          , & !out NOAHMP
!--- additional out variables:
      pahxy       = pahxy_p      , pahgxy       = pahgxy_p      , pahbxy     = pahbxy_p        , & !out NOAHMP
      pahvxy      = pahvxy_p     , qintsxy      = qintsxy_p     , qintrxy    = qintrxy_p       , & !out NOAHMP
      qdripsxy    = qdripsxy_p   , qdriprxy     = qdriprxy_p    , qthrosxy   = qthrosxy_p      , & !out NOAHMP
      qthrorxy    = qthrorxy_p   , qsnsubxy     = qsnsubxy_p    , qsnfroxy   = qsnfroxy_p      , & !out NOAHMP
      qsubcxy     = qsubcxy_p    , qfrocxy      = qfrocxy_p     , qevacxy    = qevacxy_p       , & !out NOAHMP
      qdewcxy     = qdewcxy_p    , qfrzcxy      = qfrzcxy_p     , qmeltcxy   = qmeltcxy_p      , & !out NOAHMP
      qsnbotxy    = qsnbotxy_p   , pondingxy    = pondingxy_p   , fpicexy    = fpicexy_p       , & !out NOAHMP
      rainlsm     = rainlsm_p    , snowlsm      = snowlsm_p     , forctlsm   = forctlsm_p      , & !out NOAHMP
      forcqlsm    = forcqlsm_p   , forcplsm     = forcplsm_p    , forczlsm   = forczlsm_p      , & !out NOAHMP
      forcwlsm    = forcwlsm_p   , eflxbxy      = eflxbxy_p     , soilenergy = soilenergy_p    , & !out NOAHMP
      snowenergy  = snowenergy_p , canhsxy      = canhsxy_p     ,                                & !out NOAHMP
      ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,                    &
      ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,                    &
      its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte                      &
               )


 call seaice_noah( &
      dz8w   = dz_p         , p8w3d     = pres2_hyd_p , t3d      = t_p         , &
      qv3d   = qv_p         , xice      = xice_p      , snoalb2d = snoalb_p    , &
      glw    = glw_p        , swdown    = swdown_p    , rainbl   = rainbl_p    , &
      sr     = sr_p         , qgh       = qgh_p       , tsk      = tsk_p       , &
      hfx    = hfx_p        , qfx       = qfx_p       , lh       = lh_p        , &
      grdflx = grdflx_p     , qsfc      = qsfc_p      , emiss    = sfc_emiss_p , &
      albedo = sfc_albedo_p , rib       = br_p        , cqs2     = cqs2_p      , &
      chs    = chs_p        , chs2      = chs2_p      , z02d     = z0_p        , &
      znt    = znt_p        , tslb      = tslb_p      , snow     = snow_p      , &
      snowc  = snowc_p      , snowh2d   = snowh_p     , acsnow   = acsnow_p    , &
      acsnom = acsnom_p     , sfcrunoff = sfcrunoff_p , albsi    = albsi_p     , &
      snowsi = snowsi_p     , icedepth  = icedepth_p  , dt       = dt_pbl      , &
      frpcpn = frpcpn       ,                                                    &
      seaice_albedo_opt        = seaice_albedo_opt        ,                      &
      seaice_albedo_default    = seaice_albedo_default    ,                      &
      seaice_thickness_opt     = seaice_thickness_opt     ,                      &
      seaice_thickness_default = seaice_thickness_default ,                      &
      seaice_snowdepth_opt     = seaice_snowdepth_opt     ,                      &
      seaice_snowdepth_max     = seaice_snowdepth_max     ,                      &
      seaice_snowdepth_min     = seaice_snowdepth_min     ,                      &
      xice_threshold           = xice_threshold           ,                      &
      num_soil_layers          = num_soils                ,                      &
      sf_urban_physics         = sf_urban_physics         ,                      &
      ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,    &
      ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,    &
      its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte      &
                 )


 call sfcdiags_noahmp( &
      isice  = isice    , isurban = isurban  , iswater = iswater  , lcz_1   = lcz_1    , &
      lcz_2  = lcz_2    , lcz_3   = lcz_3    , lcz_4   = lcz_4    , lcz_5   = lcz_5    , &
      lcz_6  = lcz_6    , lcz_7   = lcz_7    , lcz_8   = lcz_8    , lcz_9   = lcz_9    , &
      lcz_10 = lcz_10   , lcz_11  = lcz_11   , cp      = cp       , R_d     = R_d      , &
      rcp    = rcp      , ivgtyp  = ivgtyp_p , chs2    = chs2_p   , cqs2    = cqs2_p   , &
      psfc   = psfc_p   , hfx     = hfx_p    , qfx     = qfx_p    , qsfc    = qsfc_p   , &
      tsk    = tsk_p    , xice    = xice_p   , q2mbxy  = q2mbxy_p , q2mvxy  = q2mvxy_p , &
      t2mbxy = t2mbxy_p , t2mvxy  = t2mvxy_p , fvegxy  = fvegxy_p , q2      = q2_p     , &
      t2     = t2m_p    , th2    = th2m_p    , xice_threshold = xice_threshold         , &
      ims = ims , ime = ime , jms = jms , jme = jme ,                                    &
      its = its , ite = ite , jts = jts , jte = jte)


!copy local arrays to MPAS grid:
 call lsm_noahmp_to_MPAS(configs,mesh,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input,its,ite)

 call mpas_log_write('--- end subroutine driver_lsm_noahmp:')
 call mpas_log_write(' ')


 end subroutine driver_lsm_noahmp

!=================================================================================================================
 end module mpas_atmphys_driver_lsm_noahmp
!=================================================================================================================
