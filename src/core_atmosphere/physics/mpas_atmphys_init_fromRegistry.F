!======================================================================================================================
 module mpas_atmphys_init_fromRegistry
!======================================================================================================================
 use mpas_log
 use mpas_pool_routines

 use pumas_vars

 implicit none
 private
 public:: pumas_init_fromRegistry


 contains


!======================================================================================================================
 subroutine pumas_init_fromRegistry(configs)
!======================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 real(kind=RKIND),pointer:: pumas_berg_eff_factor_sgl,pumas_dcs_sgl
 real(kind=RKIND),pointer:: pumas_ncnst_sgl,pumas_ninst_sgl,pumas_ngnst_sgl,pumas_nrnst_sgl,pumas_nsnst_sgl

!----------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine pumas_init_fromRegistry:')

!--- initializes logicals read from the namelisti physics_cam6_pumas:
 call mpas_pool_get_config(configs,'config_pumas_precip_method'    ,pumas_precip_method      )
 call mpas_pool_get_config(configs,'config_pumas_berg_eff_factor'  ,pumas_berg_eff_factor_sgl)
 call mpas_pool_get_config(configs,'config_pumas_dcs'              ,pumas_dcs_sgl            )

 call mpas_pool_get_config(configs,'config_pumas_do_graupel'       ,pumas_do_graupel         )
 call mpas_pool_get_config(configs,'config_pumas_do_hail'          ,pumas_do_hail            )
 call mpas_pool_get_config(configs,'config_pumas_do_uniform'       ,pumas_uniform            )
 call mpas_pool_get_config(configs,'config_pumas_do_cldice'        ,pumas_do_cldice          )
 call mpas_pool_get_config(configs,'config_pumas_hetfrz_classnuc'  ,pumas_hetfrz_classnuc    )
 call mpas_pool_get_config(configs,'config_pumas_rm_supersat'      ,pumas_rm_supersat        )
 call mpas_pool_get_config(configs,'config_pumas_do_sbphysics'     ,pumas_do_sb_physics      )

 call mpas_pool_get_config(configs,'config_pumas_evap_sed_off'     ,pumas_evap_sed_off       )
 call mpas_pool_get_config(configs,'config_pumas_icenuc_rh_off'    ,pumas_icenuc_rh_off      )
 call mpas_pool_get_config(configs,'config_pumas_icenuc_use_meyers',pumas_icenuc_use_meyers  )
 call mpas_pool_get_config(configs,'config_pumas_evap_scl_ifs'     ,pumas_evap_scl_ifs       )
 call mpas_pool_get_config(configs,'config_pumas_evap_rhthrsh_ifs' ,pumas_evap_rhthrsh_ifs   )
 call mpas_pool_get_config(configs,'config_pumas_rainfreeze_ifs'   ,pumas_rainfreeze_ifs     )
 call mpas_pool_get_config(configs,'config_pumas_ifs_sed'          ,pumas_ifs_sed            )
 call mpas_pool_get_config(configs,'config_pumas_precip_fall_corr' ,pumas_precip_fall_corr   )

 call mpas_pool_get_config(configs,'config_pumas_nccons',pumas_nccons   )
 call mpas_pool_get_config(configs,'config_pumas_nicons',pumas_nicons   )
 call mpas_pool_get_config(configs,'config_pumas_ngcons',pumas_ngcons   )
 call mpas_pool_get_config(configs,'config_pumas_nrcons',pumas_nrcons   )
 call mpas_pool_get_config(configs,'config_pumas_nscons',pumas_nscons   )
 call mpas_pool_get_config(configs,'config_pumas_ncnst' ,pumas_ncnst_sgl)
 call mpas_pool_get_config(configs,'config_pumas_ninst' ,pumas_ninst_sgl)
 call mpas_pool_get_config(configs,'config_pumas_ngnst' ,pumas_ngnst_sgl)
 call mpas_pool_get_config(configs,'config_pumas_nrnst' ,pumas_nrnst_sgl)
 call mpas_pool_get_config(configs,'config_pumas_nsnst' ,pumas_nsnst_sgl)

 pumas_berg_eff_factor = pumas_berg_eff_factor_sgl
 pumas_dcs             = pumas_dcs_sgl

 pumas_ncnst = pumas_ncnst_sgl
 pumas_ninst = pumas_ninst_sgl
 pumas_ngnst = pumas_ngnst_sgl
 pumas_nrnst = pumas_nrnst_sgl
 pumas_nsnst = pumas_nsnst_sgl

 call mpas_log_write('--- enter subroutine pumas_init_fromRegistry:')

 end subroutine pumas_init_fromRegistry

!======================================================================================================================
 end module mpas_atmphys_init_fromRegistry
!======================================================================================================================
