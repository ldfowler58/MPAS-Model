! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_lsm_noahmpinit
 use mpas_dmpar
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_log

 use mpas_atmphys_vars
 use mpas_atmphys_vars_noahmp
 use module_sf_noahmpdrv

 implicit none
 private
 public:: noahmp_init_forMPAS


!Initialization of the NOAHMP parameterization of land surface processes.
!Laura D. Fowler (laura@ucar.edu) / 2022-07-01.


 contains


!=================================================================================================================
 subroutine allocate_noahmp_init
!=================================================================================================================

 if(.not.allocated(xlat_p)    ) allocate(xlat_p(ims:ime,jms:jme)  )
 if(.not.allocated(isltyp_p)  ) allocate(isltyp_p(ims:ime,jms:jme))
 if(.not.allocated(ivgtyp_p)  ) allocate(ivgtyp_p(ims:ime,jms:jme))
 if(.not.allocated(canwat_p)  ) allocate(canwat_p(ims:ime,jms:jme))
 if(.not.allocated(lai_p)     ) allocate(lai_p(ims:ime,jms:jme)   )
 if(.not.allocated(snow_p)    ) allocate(snow_p(ims:ime,jms:jme)  )
 if(.not.allocated(snowh_p)   ) allocate(snowh_p(ims:ime,jms:jme) )
 if(.not.allocated(tmn_p)     ) allocate(tmn_p(ims:ime,jms:jme)   )
 if(.not.allocated(tsk_p)     ) allocate(tsk_p(ims:ime,jms:jme)   )
 if(.not.allocated(xice_p)    ) allocate(xice_p(ims:ime,jms:jme)  )
 if(.not.allocated(dzs_p)     ) allocate(dzs_p(1:num_soils)       )
 if(.not.allocated(sh2o_p)    ) allocate(sh2o_p(ims:ime,1:num_soils,jms:jme) )
 if(.not.allocated(smois_p)   ) allocate(smois_p(ims:ime,1:num_soils,jms:jme))
 if(.not.allocated(tslb_p)    ) allocate(tslb_p(ims:ime,1:num_soils,jms:jme) )

 if(.not.allocated(isnowxy_p) ) allocate(isnowxy_p(ims:ime,jms:jme) )
 if(.not.allocated(cropcat_p) ) allocate(cropcat_p(ims:ime,jms:jme) )
 if(.not.allocated(pgsxy_p)   ) allocate(pgsxy_p(ims:ime,jms:jme)   )
 if(.not.allocated(irnumsi_p) ) allocate(irnumsi_p(ims:ime,jms:jme) )
 if(.not.allocated(irnummi_p) ) allocate(irnummi_p(ims:ime,jms:jme) )
 if(.not.allocated(irnumfi_p) ) allocate(irnumfi_p(ims:ime,jms:jme) )
 if(.not.allocated(irwatsi_p) ) allocate(irwatsi_p(ims:ime,jms:jme) )
 if(.not.allocated(irwatmi_p) ) allocate(irwatmi_p(ims:ime,jms:jme) )
 if(.not.allocated(irwatfi_p) ) allocate(irwatfi_p(ims:ime,jms:jme) )
 if(.not.allocated(ireloss_p) ) allocate(ireloss_p(ims:ime,jms:jme) )
 if(.not.allocated(irsivol_p) ) allocate(irsivol_p(ims:ime,jms:jme) )
 if(.not.allocated(irmivol_p) ) allocate(irmivol_p(ims:ime,jms:jme) )
 if(.not.allocated(irfivol_p) ) allocate(irfivol_p(ims:ime,jms:jme) )
 if(.not.allocated(irrsplh_p) ) allocate(irrsplh_p(ims:ime,jms:jme) )

 if(.not.allocated(tvxy_p)    ) allocate(tvxy_p(ims:ime,jms:jme)    ) 
 if(.not.allocated(tgxy_p)    ) allocate(tgxy_p(ims:ime,jms:jme)    )
 if(.not.allocated(canicexy_p)) allocate(canicexy_p(ims:ime,jms:jme))
 if(.not.allocated(canliqxy_p)) allocate(canliqxy_p(ims:ime,jms:jme))
 if(.not.allocated(chstarxy_p)) allocate(chstarxy_p(ims:ime,jms:jme))
 if(.not.allocated(eahxy_p)   ) allocate(eahxy_p(ims:ime,jms:jme)   )
 if(.not.allocated(tahxy_p)   ) allocate(tahxy_p(ims:ime,jms:jme)   )
 if(.not.allocated(cmxy_p)    ) allocate(cmxy_p(ims:ime,jms:jme)    )
 if(.not.allocated(chxy_p)    ) allocate(chxy_p(ims:ime,jms:jme)    )
 if(.not.allocated(fwetxy_p)  ) allocate(fwetxy_p(ims:ime,jms:jme)  )
 if(.not.allocated(sneqvoxy_p)) allocate(sneqvoxy_p(ims:ime,jms:jme))
 if(.not.allocated(alboldxy_p)) allocate(alboldxy_p(ims:ime,jms:jme))
 if(.not.allocated(qrainxy_p) ) allocate(qrainxy_p(ims:ime,jms:jme) )
 if(.not.allocated(qsnowxy_p) ) allocate(qsnowxy_p(ims:ime,jms:jme) )
 if(.not.allocated(qtdrain_p) ) allocate(qtdrain_p(ims:ime,jms:jme) )
 if(.not.allocated(wslakexy_p)) allocate(wslakexy_p(ims:ime,jms:jme))
 if(.not.allocated(zwtxy_p)   ) allocate(zwtxy_p(ims:ime,jms:jme)   )
 if(.not.allocated(waxy_p)    ) allocate(waxy_p(ims:ime,jms:jme)    )
 if(.not.allocated(wtxy_p)    ) allocate(wtxy_p(ims:ime,jms:jme)    )
 if(.not.allocated(lfmassxy_p)) allocate(lfmassxy_p(ims:ime,jms:jme))
 if(.not.allocated(rtmassxy_p)) allocate(rtmassxy_p(ims:ime,jms:jme))
 if(.not.allocated(stmassxy_p)) allocate(stmassxy_p(ims:ime,jms:jme))
 if(.not.allocated(woodxy_p)  ) allocate(woodxy_p(ims:ime,jms:jme)  )
 if(.not.allocated(stblcpxy_p)) allocate(stblcpxy_p(ims:ime,jms:jme))
 if(.not.allocated(fastcpxy_p)) allocate(fastcpxy_p(ims:ime,jms:jme))
 if(.not.allocated(xsaixy_p)  ) allocate(xsaixy_p(ims:ime,jms:jme)  )
 if(.not.allocated(grainxy_p) ) allocate(grainxy_p(ims:ime,jms:jme) )
 if(.not.allocated(gddxy_p)   ) allocate(gddxy_p(ims:ime,jms:jme)   )
 if(.not.allocated(t2mvxy_p)  ) allocate(t2mvxy_p(ims:ime,jms:jme)  ) 
 if(.not.allocated(t2mbxy_p)  ) allocate(t2mbxy_p(ims:ime,jms:jme)  ) 
 if(.not.allocated(croptype_p)) allocate(croptype_p(ims:ime,5,jms:jme))
 if(.not.allocated(tsnoxy_p)  ) allocate(tsnoxy_p(ims:ime,-2:0,jms:jme)         )
 if(.not.allocated(zsnsoxy_p) ) allocate(zsnsoxy_p(ims:ime,-2:num_soils,jms:jme))
 if(.not.allocated(snicexy_p) ) allocate(snicexy_p(ims:ime,-2:0,jms:jme)        )
 if(.not.allocated(snliqxy_p) ) allocate(snliqxy_p(ims:ime,-2:0,jms:jme)        )

 end subroutine allocate_noahmp_init

!=================================================================================================================
 subroutine deallocate_noahmp_init
!=================================================================================================================

 if(allocated(xlat_p)    ) deallocate(xlat_p    )
 if(allocated(isltyp_p)  ) deallocate(isltyp_p  )
 if(allocated(ivgtyp_p)  ) deallocate(ivgtyp_p  )
 if(allocated(canwat_p)  ) deallocate(canwat_p  )
 if(allocated(lai_p)     ) deallocate(lai_p     )
 if(allocated(snow_p)    ) deallocate(snow_p    )
 if(allocated(snowh_p)   ) deallocate(snowh_p   )
 if(allocated(tmn_p)     ) deallocate(tmn_p     )
 if(allocated(tsk_p)     ) deallocate(tsk_p     )
 if(allocated(xice_p)    ) deallocate(xice_p    )
 if(allocated(dzs_p)     ) deallocate(dzs_p     )
 if(allocated(sh2o_p)    ) deallocate(sh2o_p    )
 if(allocated(smois_p)   ) deallocate(smois_p   )
 if(allocated(tslb_p)    ) deallocate(tslb_p    )

 if(allocated(isnowxy_p) ) deallocate(isnowxy_p )
 if(allocated(cropcat_p) ) deallocate(cropcat_p )
 if(allocated(pgsxy_p)   ) deallocate(pgsxy_p   )
 if(allocated(irnumsi_p) ) deallocate(irnumsi_p )
 if(allocated(irnummi_p) ) deallocate(irnummi_p )
 if(allocated(irnumfi_p) ) deallocate(irnumfi_p )
 if(allocated(irwatsi_p) ) deallocate(irwatsi_p )
 if(allocated(irwatmi_p) ) deallocate(irwatmi_p )
 if(allocated(irwatfi_p) ) deallocate(irwatfi_p )
 if(allocated(ireloss_p) ) deallocate(ireloss_p )
 if(allocated(irsivol_p) ) deallocate(irsivol_p )
 if(allocated(irmivol_p) ) deallocate(irmivol_p )
 if(allocated(irfivol_p) ) deallocate(irfivol_p )
 if(allocated(irrsplh_p) ) deallocate(irrsplh_p )

 if(allocated(tvxy_p)    ) deallocate(tvxy_p    ) 
 if(allocated(tgxy_p)    ) deallocate(tgxy_p    )
 if(allocated(canicexy_p)) deallocate(canicexy_p)
 if(allocated(canliqxy_p)) deallocate(canliqxy_p)
 if(allocated(chstarxy_p)) deallocate(chstarxy_p)
 if(allocated(eahxy_p)   ) deallocate(eahxy_p   )
 if(allocated(tahxy_p)   ) deallocate(tahxy_p   )
 if(allocated(cmxy_p)    ) deallocate(cmxy_p    )
 if(allocated(chxy_p)    ) deallocate(chxy_p    )
 if(allocated(fwetxy_p)  ) deallocate(fwetxy_p  )
 if(allocated(sneqvoxy_p)) deallocate(sneqvoxy_p)
 if(allocated(alboldxy_p)) deallocate(alboldxy_p)
 if(allocated(qrainxy_p) ) deallocate(qrainxy_p )
 if(allocated(qsnowxy_p) ) deallocate(qsnowxy_p )
 if(allocated(qtdrain_p) ) deallocate(qtdrain_p )
 if(allocated(wslakexy_p)) deallocate(wslakexy_p)
 if(allocated(zwtxy_p)   ) deallocate(zwtxy_p   )
 if(allocated(waxy_p)    ) deallocate(waxy_p    )
 if(allocated(wtxy_p)    ) deallocate(wtxy_p    )
 if(allocated(lfmassxy_p)) deallocate(lfmassxy_p)
 if(allocated(rtmassxy_p)) deallocate(rtmassxy_p)
 if(allocated(stmassxy_p)) deallocate(stmassxy_p)
 if(allocated(woodxy_p)  ) deallocate(woodxy_p  )
 if(allocated(stblcpxy_p)) deallocate(stblcpxy_p)
 if(allocated(fastcpxy_p)) deallocate(fastcpxy_p)
 if(allocated(xsaixy_p)  ) deallocate(xsaixy_p  )
 if(allocated(grainxy_p) ) deallocate(grainxy_p )
 if(allocated(gddxy_p)   ) deallocate(gddxy_p   )
 if(allocated(t2mvxy_p)  ) deallocate(t2mvxy_p  ) 
 if(allocated(t2mbxy_p)  ) deallocate(t2mbxy_p  ) 
 if(allocated(croptype_p)) deallocate(croptype_p)
 if(allocated(tsnoxy_p)  ) deallocate(tsnoxy_p  )
 if(allocated(zsnsoxy_p) ) deallocate(zsnsoxy_p )
 if(allocated(snicexy_p) ) deallocate(snicexy_p )
 if(allocated(snliqxy_p) ) deallocate(snliqxy_p )

 end subroutine deallocate_noahmp_init

!=================================================================================================================
 subroutine noahmp_init_forMPAS(mesh,configs,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: sfc_input

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: output_noahmp

!local pointers:
 logical,pointer:: config_do_restart 
 character(len=StrKIND),pointer:: mminlu

 integer,dimension(:),pointer:: isltyp,ivgtyp
 integer,dimension(:),pointer:: isnowxy,cropcat

 real(kind=RKIND),dimension(:),pointer:: latCell
 real(kind=RKIND),dimension(:),pointer:: canwat,lai,skintemp,snow,snowh,tmn,xice
 real(kind=RKIND),dimension(:),pointer:: alboldxy,canicexy,canliqxy,chxy,cmxy,eahxy,fastcpxy,fwetxy,gddxy,    &
                                         grainxy,lfmassxy,qrainxy,qsnowxy,qtdrain,rtmassxy,sneqvoxy,stblcpxy, &
                                         stmassxy,tahxy,tgxy,tvxy,xsaixy,waxy,woodxy,wslakexy,wtxy,zwtxy
 real(kind=RKIND),dimension(:),pointer::chstarxy,t2mbxy,t2mvxy

 real(kind=RKIND),dimension(:,:),pointer:: dzs,sh2o,smois,tslb
 real(kind=RKIND),dimension(:,:),pointer:: snicexy,snliqxy,tsnoxy,zsnsoxy
 real(kind=RKIND),dimension(:,:),pointer:: croptype

 real(kind=RKIND),dimension(:),pointer:: xland

!local variables:
 logical,parameter:: allowed_to_read = .false.
 logical,parameter:: fndsoilw = .true.
 logical,parameter:: fndsnowh = .true.
 integer:: i,j,n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine noahmp_init_forMPAS:')

!abort the initialization if iopt_run is equal to 5:
 if(iopt_run >= 5) &
    call physics_error_fatal("subroutine noahmpinit_forMPAS: iopt_run >= 5 is not available." // &
                             " reset iopt_run to a different value.")

!allocate local arrays:
 call allocate_noahmp_init()

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)

 call mpas_pool_get_array(mesh,'latCell',latCell)

 call mpas_pool_get_array(diag_physics,'canwat',canwat)
 call mpas_pool_get_array(diag_physics,'lai'   ,lai   )

 call mpas_pool_get_array(diag_physics_noahmp,'alboldxy',alboldxy)
 call mpas_pool_get_array(diag_physics_noahmp,'canicexy',canicexy)
 call mpas_pool_get_array(diag_physics_noahmp,'canliqxy',canliqxy)
 call mpas_pool_get_array(diag_physics_noahmp,'chxy'    ,chxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'cmxy'    ,cmxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'eahxy'   ,eahxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'fastcpxy',fastcpxy)
 call mpas_pool_get_array(diag_physics_noahmp,'fwetxy'  ,fwetxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'gddxy'   ,gddxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'grainxy' ,grainxy )
 call mpas_pool_get_array(diag_physics_noahmp,'isnowxy' ,isnowxy )
 call mpas_pool_get_array(diag_physics_noahmp,'lfmassxy',lfmassxy)
 call mpas_pool_get_array(diag_physics_noahmp,'qrainxy' ,qrainxy )
 call mpas_pool_get_array(diag_physics_noahmp,'qsnowxy' ,qsnowxy )
 call mpas_pool_get_array(diag_physics_noahmp,'qtdrain' ,qtdrain )
 call mpas_pool_get_array(diag_physics_noahmp,'rtmassxy',rtmassxy)
 call mpas_pool_get_array(diag_physics_noahmp,'sneqvoxy',sneqvoxy)
 call mpas_pool_get_array(diag_physics_noahmp,'snicexy' ,snicexy )
 call mpas_pool_get_array(diag_physics_noahmp,'snliqxy' ,snliqxy )
 call mpas_pool_get_array(diag_physics_noahmp,'stblcpxy',stblcpxy)
 call mpas_pool_get_array(diag_physics_noahmp,'stmassxy',stmassxy)
 call mpas_pool_get_array(diag_physics_noahmp,'tahxy'   ,tahxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'tgxy'    ,tgxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'tsnoxy'  ,tsnoxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'tvxy'    ,tvxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'xsaixy'  ,xsaixy  )
 call mpas_pool_get_array(diag_physics_noahmp,'waxy'    ,waxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'woodxy'  ,woodxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'wslakexy',wslakexy)
 call mpas_pool_get_array(diag_physics_noahmp,'wtxy'    ,wtxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'zsnsoxy' ,zsnsoxy )
 call mpas_pool_get_array(diag_physics_noahmp,'zwtxy'   ,zwtxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'cropcat' ,cropcat )
 call mpas_pool_get_array(diag_physics_noahmp,'croptype',croptype)

 call mpas_pool_get_array(output_noahmp,'chstarxy',chstarxy)
 call mpas_pool_get_array(output_noahmp,'t2mbxy'  ,t2mbxy  )
 call mpas_pool_get_array(output_noahmp,'t2mvxy'  ,t2mvxy  )

 call mpas_pool_get_array(sfc_input,'mminlu'  ,mminlu  )
 call mpas_pool_get_array(sfc_input,'isltyp'  ,isltyp  )
 call mpas_pool_get_array(sfc_input,'ivgtyp'  ,ivgtyp  )
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'tmn'     ,tmn     )
 call mpas_pool_get_array(sfc_input,'xice'    ,xice    )
 call mpas_pool_get_array(sfc_input,'dzs'     ,dzs     )
 call mpas_pool_get_array(sfc_input,'sh2o'    ,sh2o    )
 call mpas_pool_get_array(sfc_input,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input,'tslb'    ,tslb    )

 call mpas_pool_get_array(sfc_input,'xland'   ,xland   )

 do n = 1,num_soils
    dzs_p(n) = dzs(n,its)
 enddo

 do j = jts,jte
 do n = 1,5
 do i = its,ite
!   croptype_p(i,n,j) = croptype(n,i)
    croptype_p(i,n,j) = 0.001
 enddo
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
    sh2o_p(i,n,j)    = sh2o(n,i)
    smois_p(i,n,j)   = smois(n,i)
    tslb_p(i,n,j)    = tslb(n,i)
    !finish the initialization of zsnsoxy_p:
    zsnsoxy_p(i,n,j) = zsnsoxy(n+num_soils-1,i)
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    xlat_p(i,j)     = latCell(i)
    isltyp_p(i,j)   = isltyp(i) 
    ivgtyp_p(i,j)   = ivgtyp(i)
    cropcat_p(i,j)  = cropcat(i)
    canwat_p(i,j)   = canwat(i)
    snow_p(i,j)     = snow(i)
    snowh_p(i,j)    = snowh(i)
    tmn_p(i,j)      = tmn(i)
    tsk_p(i,j)      = skintemp(i)
    xice_p(i,j)     = xice(i)
    isnowxy_p(i,j)  = isnowxy(i)
    tvxy_p(i,j)     = tvxy(i)
    tgxy_p(i,j)     = tgxy(i)
    canicexy_p(i,j) = canicexy(i) 
    canliqxy_p(i,j) = canliqxy(i)
    eahxy_p(i,j)    = eahxy(i)
    tahxy_p(i,j)    = tahxy(i)
    cmxy_p(i,j)     = cmxy(i)
    chxy_p(i,j)     = chxy(i)
    fwetxy_p(i,j)   = fwetxy(i)
    sneqvoxy_p(i,j) = sneqvoxy(i)
    alboldxy_p(i,j) = alboldxy(i)
    qrainxy_p(i,j)  = qrainxy(i)
    qsnowxy_p(i,j)  = qsnowxy(i)
    qtdrain_p(i,j)  = qtdrain(i)
    wslakexy_p(i,j) = wslakexy(i)
    zwtxy_p(i,j)    = zwtxy(i)
    waxy_p(i,j)     = waxy(i)
    wtxy_p(i,j)     = wtxy(i)
    lfmassxy_p(i,j) = lfmassxy(i)
    rtmassxy_p(i,j) = rtmassxy(i)
    stmassxy_p(i,j) = stmassxy(i)
    woodxy_p(i,j)   = woodxy(i)
    grainxy_p(i,j)  = grainxy(i)
    gddxy_p(i,j)    = gddxy(i)
    stblcpxy_p(i,j) = stblcpxy(i)
    fastcpxy_p(i,j) = fastcpxy(i)
    xsaixy_p(i,j)   = xsaixy(i)
    lai_p(i,j)      = lai(i)
    t2mvxy_p(i,j)   = t2mvxy(i)
    t2mbxy_p(i,j)   = t2mbxy(i)
    chstarxy_p(i,j) = chstarxy(i)
 enddo
 enddo

!note: we do not include the optional groundwater-related variables:
 call noahmp_init( &
             mminlu    = mminlu       , snow      = snow_p     , snowh     = snowh_p    , &
             canwat    = canwat_p     , isltyp    = isltyp_p   , ivgtyp    = ivgtyp_p   , &
             tslb      = tslb_p       , smois     = smois_p    , sh2o      = sh2o_p     , &
             dzs       = dzs_p        , tsk       = tsk_p      , tmn       = tmn_p      , &
             xice      = xice_p       , fndsoilw  = fndsoilw   , fndsnowh  = fndsnowh   , &
             lai       = lai_p        , nsoil     = num_soils  , iopt_run  = iopt_run   , & 
             isnowxy   = isnowxy_p    , tvxy      = tvxy_p     , tgxy      = tgxy_p     , &  
             canicexy  = canicexy_p   , canliqxy  = canliqxy_p , eahxy     = eahxy_p    , & 
             tahxy     = tahxy_p      , cmxy      = cmxy_p     , chxy      = chxy_p     , & 
             fwetxy    = fwetxy_p     , sneqvoxy  = sneqvoxy_p , alboldxy  = alboldxy_p , & 
             qrainxy   = qrainxy_p    , qsnowxy   = qsnowxy_p  , wslakexy  = wslakexy_p , &
             zwtxy     = zwtxy_p      , waxy      = waxy_p     , wtxy      = wtxy_p     , &
             tsnoxy    = tsnoxy_p     , zsnsoxy   = zsnsoxy_p  , snicexy   = snicexy_p  , &
             snliqxy   = snliqxy_p    , lfmassxy  = lfmassxy_p , rtmassxy  = rtmassxy_p , &
             stmassxy  = stmassxy_p   , woodxy    = woodxy_p   , stblcpxy  = stblcpxy_p , &
             fastcpxy  = fastcpxy_p   , xsaixy    = xsaixy_p   , grainxy   = grainxy_p  , &
             gddxy     = gddxy_p      , t2mvxy    = t2mvxy_p   , t2mbxy    = t2mbxy_p   , &
             chstarxy  = chstarxy_p   , cropcat   = cropcat_p  , croptype  = croptype_p , &
             qtdrain   = qtdrain_p    ,                                                   &
             irnumsi   = irnumsi_p    , irnummi   = irnummi_p  , irnumfi   = irnumfi_p  , &
             irwatsi   = irwatsi_p    , irwatmi   = irwatmi_p  , irwatfi   = irwatfi_p  , &
             ireloss   = ireloss_p    , irsivol   = irsivol_p  , irmivol   = irmivol_p  , &
             irfivol   = irfivol_p    , irrsplh   = irrsplh_p  , iopt_crop = iopt_crop  , &
             iopt_irr  = iopt_irr     , iopt_irrm = iopt_irrm  , xlat      = xlat_p     , &
             restart   = config_do_restart                                              , &
             allowed_to_read  = allowed_to_read                                         , &
             sf_urban_physics = sf_urban_physics                                        , &
             ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde      , &
             ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme      , &
             its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                 )


 do j = jts,jte
 do n = -2,0
 do i = its,ite
    tsnoxy(n+3,i)  = tsnoxy_p(i,n,j)
    snicexy(n+3,i) = snicexy_p(i,n,j)
    snliqxy(n+3,i) = snliqxy_p(i,n,j)
    zsnsoxy(n+3,i) = zsnsoxy_p(i,n,j)
 enddo
 enddo
 do n = 1,num_soils
 do i = its,ite
    sh2o(n,i)  = sh2o_p(i,n,j)
    smois(n,i) = smois_p(i,n,j)
    tslb(n,i)  = tslb_p(i,n,j)
    !finish the initialization of zsnsoxy:
    zsnsoxy(n+num_soils-1,i) = zsnsoxy_p(i,n,j)
 enddo
 enddo
 enddo

 do j = jts,jte
 do i = its,ite
    cropcat(i)  = cropcat_p(i,j)
    canwat(i)   = canwat_p(i,j)
    snow(i)     = snow_p(i,j)
    snowh(i)    = snowh_p(i,j)
    tmn(i)      = tmn_p(i,j)
    isnowxy(i)  = isnowxy_p(i,j)
    tvxy(i)     = tvxy_p(i,j)
    tgxy(i)     = tgxy_p(i,j)
    canicexy(i) = canicexy_p(i,j)
    canliqxy(i) = canliqxy_p(i,j)
    eahxy(i)    = eahxy_p(i,j)
    tahxy(i)    = tahxy_p(i,j)
    cmxy(i)     = cmxy_p(i,j)
    chxy(i)     = chxy_p(i,j)
    fwetxy(i)   = fwetxy_p(i,j)
    sneqvoxy(i) = sneqvoxy_p(i,j)
    alboldxy(i) = alboldxy_p(i,j)
    qrainxy(i)  = qrainxy_p(i,j)
    qsnowxy(i)  = qsnowxy_p(i,j)
    qtdrain(i)  = qtdrain_p(i,j)
    wslakexy(i) = wslakexy_p(i,j)
    zwtxy(i)    = zwtxy_p(i,j)
    waxy(i)     = waxy_p(i,j)
    wtxy(i)     = wtxy_p(i,j)
    lfmassxy(i) = lfmassxy_p(i,j)
    rtmassxy(i) = rtmassxy_p(i,j)
    stmassxy(i) = stmassxy_p(i,j)
    woodxy(i)   = woodxy_p(i,j)
    grainxy(i)  = grainxy_p(i,j)
    gddxy(i)    = gddxy_p(i,j)
    stblcpxy(i) = stblcpxy_p(i,j)
    fastcpxy(i) = fastcpxy_p(i,j)
    xsaixy(i)   = xsaixy_p(i,j)
    lai(i)      = lai_p(i,j)
    t2mvxy(i)   = t2mvxy_p(i,j)
    t2mbxy(i)   = t2mbxy_p(i,j)
    chstarxy(i) = chstarxy_p(i,j)
 enddo
 enddo

!deallocate local arrays:
 call deallocate_noahmp_init()

 call mpas_log_write('--- exit subroutine noahmp_init_forMPAS:')
 call mpas_log_write(' ')

 end subroutine noahmp_init_forMPAS

!=================================================================================================================
 end module mpas_atmphys_lsm_noahmpinit
!=================================================================================================================
