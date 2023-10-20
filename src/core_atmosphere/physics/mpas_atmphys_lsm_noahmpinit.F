! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_lsm_noahmpinit
 use mpas_log
 use mpas_pool_routines

 use mpas_atmphys_vars,only: mpas_noahmp

 use NoahmpInitMainMod,only : NoahmpInitMain
 use NoahmpIOVarInitMod,only: NoahmpIOVarInitDefault
 use NoahmpIOVarType
 use NoahmpReadNamelistMod
 use NoahmpReadTableMod,only: NoahmpReadTable


 private
 public:: init_lsm_noahmp


 contains


!=================================================================================================================
 subroutine init_lsm_noahmp(configs,mesh,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: output_noahmp
 type(mpas_pool_type),intent(inout):: sfc_input

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_lsm_noahmp:')


!--- initialize dimensions:
 call noahmp_read_dimensions(mesh)

 
!--- initialize namelist options:
 call noahmp_read_namelist(configs)


!--- initialize noahmp:
 call noahmp_init(configs,mesh,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)


 call mpas_log_write('--- end subroutine init_lsm_noahmp:')
 call mpas_log_write(' ')

 end subroutine init_lsm_noahmp

!=================================================================================================================
 subroutine noahmp_read_dimensions(mesh)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: mesh

!--- local variables and pointers:
 integer,pointer:: nCellsSolve,nVertLevels
 integer,pointer:: nSoilLevels,nSnowLevels

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine noahmp_read_dimensions:')


 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nSoilLevels',nSoilLevels)
 call mpas_pool_get_dimension(mesh,'nSnowLevels',nSnowLevels)

 mpas_noahmp%ims = 1
 mpas_noahmp%ime = nCellsSolve
 mpas_noahmp%ids = 1
 mpas_noahmp%ide = nCellsSolve + 1
 mpas_noahmp%its = 1
 mpas_noahmp%ite = nCellsSolve

 mpas_noahmp%jms = 1
 mpas_noahmp%jme = 1
 mpas_noahmp%jds = 1
 mpas_noahmp%jde = 1
 mpas_noahmp%jts = 1
 mpas_noahmp%jte = 1

 mpas_noahmp%kms = 1
 mpas_noahmp%kme = nVertLevels + 1
 mpas_noahmp%kds = 1
 mpas_noahmp%kde = nVertLevels + 1
 mpas_noahmp%kts = 1
 mpas_noahmp%kte = nVertLevels

 mpas_noahmp%nsoil = nSoilLevels
 mpas_noahmp%nsnow = nSnowLevels

 call mpas_log_write('    ims = $i   ime = $i', intArgs=(/mpas_noahmp%ims,mpas_noahmp%ime/))
 call mpas_log_write('    jms = $i   jme = $i', intArgs=(/mpas_noahmp%jms,mpas_noahmp%jme/))
 call mpas_log_write('    kms = $i   kme = $i', intArgs=(/mpas_noahmp%kms,mpas_noahmp%kme/))
 call mpas_log_write('    ids = $i   ide = $i', intArgs=(/mpas_noahmp%ids,mpas_noahmp%ide/))
 call mpas_log_write('    jds = $i   jde = $i', intArgs=(/mpas_noahmp%jds,mpas_noahmp%jde/))
 call mpas_log_write('    kds = $i   kde = $i', intArgs=(/mpas_noahmp%kds,mpas_noahmp%kde/))
 call mpas_log_write('    its = $i   ite = $i', intArgs=(/mpas_noahmp%its,mpas_noahmp%ite/))
 call mpas_log_write('    jts = $i   jte = $i', intArgs=(/mpas_noahmp%jts,mpas_noahmp%jte/))
 call mpas_log_write('    kts = $i   kte = $i', intArgs=(/mpas_noahmp%kts,mpas_noahmp%kte/))


 call mpas_log_write('--- end subroutine noahmp_read_dimensions:')

 end subroutine noahmp_read_dimensions

!=================================================================================================================
 subroutine noahmp_read_namelist(configs)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs


!--- local variables and pointers:
 integer,pointer:: iopt_dveg  , iopt_crs  , iopt_btr  , iopt_runsrf , iopt_runsub , iopt_sfc  , iopt_frz  , &
                   iopt_inf   , iopt_rad  , iopt_alb  , iopt_snf    , iopt_tksno  , iopt_tbot , iopt_stc  , &
                   iopt_gla   , iopt_rsf  , iopt_soil , iopt_pedo   , iopt_crop   , iopt_irr  , iopt_irrm , &
                   iopt_infdv , iopt_tdrn

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine noahmp_read_namelist:')

 call mpas_pool_get_config(configs,'config_noahmp_iopt_dveg'  ,iopt_dveg  )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_crs'   ,iopt_crs   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_btr'   ,iopt_btr   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_runsrf',iopt_runsrf)
 call mpas_pool_get_config(configs,'config_noahmp_iopt_runsub',iopt_runsub)
 call mpas_pool_get_config(configs,'config_noahmp_iopt_sfc'   ,iopt_sfc   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_frz'   ,iopt_frz   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_inf'   ,iopt_inf   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_rad'   ,iopt_rad   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_alb'   ,iopt_alb   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_snf'   ,iopt_snf   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_tksno' ,iopt_tksno )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_tbot'  ,iopt_tbot  )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_stc'   ,iopt_stc   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_gla'   ,iopt_gla   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_rsf'   ,iopt_rsf   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_soil'  ,iopt_soil  )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_pedo'  ,iopt_pedo  )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_crop'  ,iopt_crop  )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_irr'   ,iopt_irr   )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_irrm'  ,iopt_irrm  )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_infdv' ,iopt_infdv )
 call mpas_pool_get_config(configs,'config_noahmp_iopt_tdrn'  ,iopt_tdrn  )

 mpas_noahmp%iopt_dveg   = iopt_dveg
 mpas_noahmp%iopt_crs    = iopt_crs
 mpas_noahmp%iopt_btr    = iopt_btr
 mpas_noahmp%iopt_runsrf = iopt_runsrf
 mpas_noahmp%iopt_runsub = iopt_runsub
 mpas_noahmp%iopt_sfc    = iopt_sfc
 mpas_noahmp%iopt_frz    = iopt_frz
 mpas_noahmp%iopt_inf    = iopt_inf
 mpas_noahmp%iopt_rad    = iopt_rad
 mpas_noahmp%iopt_alb    = iopt_alb
 mpas_noahmp%iopt_snf    = iopt_snf
 mpas_noahmp%iopt_tksno  = iopt_tksno
 mpas_noahmp%iopt_tbot   = iopt_tbot
 mpas_noahmp%iopt_stc    = iopt_stc
 mpas_noahmp%iopt_gla    = iopt_gla
 mpas_noahmp%iopt_rsf    = iopt_rsf
 mpas_noahmp%iopt_soil   = iopt_soil
 mpas_noahmp%iopt_pedo   = iopt_pedo
 mpas_noahmp%iopt_crop   = iopt_crop
 mpas_noahmp%iopt_irr    = iopt_irr
 mpas_noahmp%iopt_irrm   = iopt_irrm
 mpas_noahmp%iopt_infdv  = iopt_infdv
 mpas_noahmp%iopt_tdrn   = iopt_tdrn

 call mpas_log_write('--- end subroutine noahmp_read_namelist:')

 end subroutine noahmp_read_namelist

!=================================================================================================================
 subroutine noahmp_init(configs,mesh,diag_physics,diag_physics_noahmp,output_noahmp,sfc_input)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: diag_physics_noahmp
 type(mpas_pool_type),intent(inout):: output_noahmp
 type(mpas_pool_type),intent(inout):: sfc_input

!local pointers:
 logical,pointer:: urban_physics
 character(len=StrKIND),pointer:: mminlu

 integer,dimension(:),pointer:: isltyp,ivgtyp
 integer,dimension(:),pointer:: isnowxy,cropcat
 integer,dimension(:),pointer:: irnumsi,irnummi,irnumfi
 
 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:),pointer:: areaCell,latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: canwat,lai,skintemp,snow,snowh,tmn,xice
 real(kind=RKIND),dimension(:),pointer:: alboldxy,canicexy,canliqxy,chxy,cmxy,eahxy,fastcpxy,fwetxy,gddxy,    &
                                         grainxy,lfmassxy,qrainxy,qsnowxy,qtdrain,rtmassxy,sneqvoxy,stblcpxy, &
                                         stmassxy,tahxy,tgxy,tvxy,xsaixy,waxy,woodxy,wslakexy,wtxy,zwtxy
 real(kind=RKIND),dimension(:),pointer:: irwatsi,ireloss,irrsplh,irwatmi,irmivol,irwatfi,irfivol
 real(kind=RKIND),dimension(:),pointer:: t2mbxy,t2mvxy
 
 real(kind=RKIND),dimension(:,:),pointer:: dzs,sh2o,smois,tslb
 real(kind=RKIND),dimension(:,:),pointer:: snicexy,snliqxy,tsnoxy,zsnsoxy
 real(kind=RKIND),dimension(:,:),pointer:: croptype

!local variables and pointers:
 logical,parameter:: fndsnowh = .true.
 integer:: i,its,ite,ns,nsoil,nsnow

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine noahmp_init:')


!--- allocate Noahmp arrays:
 call mpas_log_write('--- enter subroutine NoahmpIOVarInitDefault:')
 call NoahmpIOVarInitDefault(mpas_noahmp)
 call mpas_log_write('--- end subroutine NoahmpIOVarInitDefault:')


!--- read NoahmpTable.TBL:
 call mpas_pool_get_array(sfc_input,'mminlu',mminlu)
 mpas_noahmp%llanduse = mminlu

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine NoahmpReadTable:')
 call NoahmpReadTable(mpas_noahmp)
 call mpas_log_write('--- end subroutine NoahmpReadTable:')


!--- initialize Noahmp:
 call mpas_pool_get_config(configs,'config_urban_physics',urban_physics)
 call mpas_pool_get_config(configs,'config_dt',dt)

 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 call mpas_pool_get_array(sfc_input,'dzs'     ,dzs     )
 call mpas_pool_get_array(sfc_input,'isltyp'  ,isltyp  )
 call mpas_pool_get_array(sfc_input,'ivgtyp'  ,ivgtyp  )
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'tmn'     ,tmn     )
 call mpas_pool_get_array(sfc_input,'xice'    ,xice    )
 call mpas_pool_get_array(sfc_input,'sh2o'    ,sh2o    )
 call mpas_pool_get_array(sfc_input,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input,'tslb'    ,tslb    )

 call mpas_pool_get_array(diag_physics,'canwat',canwat)
 call mpas_pool_get_array(diag_physics,'lai',lai)

 call mpas_pool_get_array(diag_physics_noahmp,'alboldxy',alboldxy)
 call mpas_pool_get_array(diag_physics_noahmp,'canicexy',canicexy)
 call mpas_pool_get_array(diag_physics_noahmp,'canliqxy',canliqxy)
 call mpas_pool_get_array(diag_physics_noahmp,'chxy'    ,chxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'cmxy'    ,cmxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'cropcat' ,cropcat )
 call mpas_pool_get_array(diag_physics_noahmp,'eahxy'   ,eahxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'fastcpxy',fastcpxy)
 call mpas_pool_get_array(diag_physics_noahmp,'fwetxy'  ,fwetxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'gddxy'   ,gddxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'grainxy' ,grainxy )
 call mpas_pool_get_array(diag_physics_noahmp,'lfmassxy',lfmassxy)
 call mpas_pool_get_array(diag_physics_noahmp,'qrainxy' ,qrainxy )
 call mpas_pool_get_array(diag_physics_noahmp,'qsnowxy' ,qsnowxy )
 call mpas_pool_get_array(diag_physics_noahmp,'qtdrain' ,qtdrain )
 call mpas_pool_get_array(diag_physics_noahmp,'rtmassxy',rtmassxy)
 call mpas_pool_get_array(diag_physics_noahmp,'sneqvoxy',sneqvoxy)
 call mpas_pool_get_array(diag_physics_noahmp,'stblcpxy',stblcpxy)
 call mpas_pool_get_array(diag_physics_noahmp,'stmassxy',stmassxy)
 call mpas_pool_get_array(diag_physics_noahmp,'tahxy'   ,tahxy   )
 call mpas_pool_get_array(diag_physics_noahmp,'tgxy'    ,tgxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'tvxy'    ,tvxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'waxy'    ,waxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'woodxy'  ,woodxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'wslakexy',wslakexy)
 call mpas_pool_get_array(diag_physics_noahmp,'wtxy'    ,wtxy    )
 call mpas_pool_get_array(diag_physics_noahmp,'xsaixy'  ,xsaixy  )
 call mpas_pool_get_array(diag_physics_noahmp,'zwtxy'   ,zwtxy   )

 call mpas_pool_get_array(diag_physics_noahmp,'croptype',croptype)

 call mpas_pool_get_array(diag_physics_noahmp,'irnumsi' ,irnumsi )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatsi' ,irwatsi )
 call mpas_pool_get_array(diag_physics_noahmp,'ireloss' ,ireloss )
 call mpas_pool_get_array(diag_physics_noahmp,'irrsplh' ,irrsplh )
 call mpas_pool_get_array(diag_physics_noahmp,'irnummi' ,irnummi )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatmi' ,irwatmi )
 call mpas_pool_get_array(diag_physics_noahmp,'irmivol' ,irmivol )
 call mpas_pool_get_array(diag_physics_noahmp,'irnumfi' ,irnumfi )
 call mpas_pool_get_array(diag_physics_noahmp,'irwatfi' ,irwatfi )
 call mpas_pool_get_array(diag_physics_noahmp,'irfivol', irfivol )

 call mpas_pool_get_array(diag_physics_noahmp,'isnowxy' ,isnowxy )
 call mpas_pool_get_array(diag_physics_noahmp,'snicexy' ,snicexy )
 call mpas_pool_get_array(diag_physics_noahmp,'snliqxy' ,snliqxy )
 call mpas_pool_get_array(diag_physics_noahmp,'tsnoxy'  ,tsnoxy  )
 call mpas_pool_get_array(diag_physics_noahmp,'zsnsoxy' ,zsnsoxy )

 call mpas_pool_get_array(output_noahmp,'t2mbxy',t2mbxy)
 call mpas_pool_get_array(output_noahmp,'t2mvxy',t2mvxy)
 
 its   = mpas_noahmp%its
 ite   = mpas_noahmp%ite
 nsoil = mpas_noahmp%nsoil
 nsnow = mpas_noahmp%nsnow

 mpas_noahmp%sf_urban_physics = 0
 if(urban_physics) mpas_noahmp%sf_urban_physics = 1

 mpas_noahmp%fndsnowh = fndsnowh
 mpas_noahmp%dtbl     = dt

 do i = its, ite
    mpas_noahmp%areaxy(i) = areaCell(i)
    mpas_noahmp%xlat(i)   = latCell(i)
    mpas_noahmp%xlong(i)  = lonCell(i)
    mpas_noahmp%isltyp(i) = isltyp(i)
    mpas_noahmp%ivgtyp(i) = ivgtyp(i)
    mpas_noahmp%snow(i)   = snow(i)
    mpas_noahmp%snowh(i)  = snowh(i)
    mpas_noahmp%tsk(i)    = skintemp(i)
    mpas_noahmp%tmn(i)    = tmn(i)
    mpas_noahmp%xice(i)   = xice(i)
 enddo

 do ns = 1, nsoil
    mpas_noahmp%dzs(ns) = dzs(ns,its)
    do i = its,ite
       mpas_noahmp%sh2o(i,ns)  = sh2o(ns,i)
       mpas_noahmp%smois(i,ns) = smois(ns,i)
       mpas_noahmp%tslb(i,ns)  = tslb(ns,i)
    enddo
 enddo


 call NoahmpInitMain(mpas_noahmp)


 do i = its, ite
    canwat(i) = mpas_noahmp%canwat(i)
    lai(i)    = mpas_noahmp%lai(i)
    snow(i)   = mpas_noahmp%snow(i)   ! in mm (check unit in noahmp driver).
    snowh(i)  = mpas_noahmp%snowh(i)  ! in m  (check unit in noahmp driver).
 enddo

 do ns = 1, nsoil
    do i = its,ite
       sh2o(ns,i)  = mpas_noahmp%sh2o(i,ns)
       smois(ns,i) = mpas_noahmp%smois(i,ns)
       tslb(ns,i)  = mpas_noahmp%tslb(i,ns)
    enddo
 enddo

 do i = its,ite
    alboldxy(i) = mpas_noahmp%alboldxy(i)
    canicexy(i) = mpas_noahmp%canicexy(i)
    canliqxy(i) = mpas_noahmp%canliqxy(i)
    chxy(i)     = mpas_noahmp%chxy(i)
    cmxy(i)     = mpas_noahmp%cmxy(i)
    cropcat(i)  = mpas_noahmp%cropcat(i)
    eahxy(i)    = mpas_noahmp%eahxy(i)
    fastcpxy(i) = mpas_noahmp%fastcpxy(i)
    fwetxy(i)   = mpas_noahmp%fwetxy(i)
    gddxy(i)    = mpas_noahmp%gddxy(i)
    grainxy(i)  = mpas_noahmp%grainxy(i)
    lfmassxy(i) = mpas_noahmp%lfmassxy(i)
    qrainxy(i)  = mpas_noahmp%qrainxy(i)
    qsnowxy(i)  = mpas_noahmp%qsnowxy(i)
    qtdrain(i)  = mpas_noahmp%qtdrain(i)
    rtmassxy(i) = mpas_noahmp%rtmassxy(i)
    sneqvoxy(i) = mpas_noahmp%sneqvoxy(i)
    stblcpxy(i) = mpas_noahmp%stblcpxy(i)
    stmassxy(i) = mpas_noahmp%stmassxy(i)
    tahxy(i)    = mpas_noahmp%tahxy(i)
    tgxy(i)     = mpas_noahmp%tgxy(i)
    tvxy(i)     = mpas_noahmp%tvxy(i)
    waxy(i)     = mpas_noahmp%waxy(i)
    woodxy(i)   = mpas_noahmp%woodxy(i)
    wslakexy(i) = mpas_noahmp%wslakexy(i)
    wtxy(i)     = mpas_noahmp%wtxy(i)
    xsaixy(i)   = mpas_noahmp%xsaixy(i)
    zwtxy(i)    = mpas_noahmp%zwtxy(i)

    t2mbxy(i)   = mpas_noahmp%t2mbxy(i)
    t2mvxy(i)   = mpas_noahmp%t2mvxy(i)
 enddo

 do i = its,ite
    isnowxy(i) = mpas_noahmp%isnowxy(i)
 enddo

 do i = its, ite
    irnumsi(i) = mpas_noahmp%irnumsi(i)
    irwatsi(i) = mpas_noahmp%irwatsi(i)
    ireloss(i) = mpas_noahmp%ireloss(i)
    irrsplh(i) = mpas_noahmp%irrsplh(i)
    irnummi(i) = mpas_noahmp%irnummi(i)
    irwatmi(i) = mpas_noahmp%irwatmi(i)
    irmivol(i) = mpas_noahmp%irmivol(i)
    irnumfi(i) = mpas_noahmp%irnumfi(i)
    irwatfi(i) = mpas_noahmp%irwatfi(i)
    irfivol(i) = mpas_noahmp%irfivol(i)
 enddo
 

 call mpas_log_write('--- end subroutine noahmp_init:')

 end subroutine noahmp_init

!=================================================================================================================
 end module mpas_atmphys_lsm_noahmpinit
!=================================================================================================================
