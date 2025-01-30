!==================================================================================================================
 module mpas_chemistry_gocart2G_driver
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_subpool

 use mpas_chemistry_gocart2G_diagnostics,only: CA2G_bc_diagnostics, &
                                               CA2G_br_diagnostics, &
                                               CA2G_oc_diagnostics, &
                                               DU2G_diagnostics,    &
                                               NI2G_diagnostics,    &
                                               SS2G_diagnostics,    &
                                               SU2G_diagnostics,    &
                                               GOCART2G_diagnostics
 use mpas_chemistry_gocart2G_interface
 use mpas_chemistry_gocart2G_manager,only: iyear,imonth,iday,ihour,iminute,isecond
 use mpas_chemistry_gocart2G_vars,only: mpas_gocart2G,           &
                                        CA2G_bc,CA2G_bc_params,  &
                                        CA2G_br,CA2G_br_params,  &
                                        CA2G_oc,CA2G_oc_params,  &
                                        DU2G,DU2G_params,        &
                                        NI2G,NI2G_params,        &
                                        SS2G,SS2G_params,        &
                                        SU2G,SU2G_params,        &
                                        GOCART2G,GOCART2G_params


 implicit none
 private
 public:: gocart2G_driver


 contains


!==================================================================================================================
 subroutine gocart2G_driver(domain,itimestep,xtime_s)
!==================================================================================================================

!inout arguments:
 type(domain_type),intent(inout):: domain
 integer,intent(in):: itimestep
 real(kind=RKIND),intent(in):: xtime_s

!local variables and pointers:
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: state
 type(mpas_pool_type),pointer:: diag
 type(mpas_pool_type),pointer:: diag_physics
 type(mpas_pool_type),pointer:: sfc_input
 type(mpas_pool_type),pointer:: CAMS_emissions
 type(mpas_pool_type),pointer:: gocart2G_backgrounds
 type(mpas_pool_type),pointer:: gocart2G_met

 type(mpas_pool_type),pointer:: CA2G_bc_diags
 type(mpas_pool_type),pointer:: CA2G_br_diags
 type(mpas_pool_type),pointer:: CA2G_oc_diags
 type(mpas_pool_type),pointer:: DU2G_diags
 type(mpas_pool_type),pointer:: NI2G_diags
 type(mpas_pool_type),pointer:: SS2G_diags
 type(mpas_pool_type),pointer:: SU2G_diags
 type(mpas_pool_type),pointer:: GOCART2G_diags

 type(mpas_pool_type),pointer:: CA2G_bc_aops
 type(mpas_pool_type),pointer:: CA2G_br_aops
 type(mpas_pool_type),pointer:: CA2G_oc_aops
 type(mpas_pool_type),pointer:: DU2G_aops
 type(mpas_pool_type),pointer:: NI2G_aops
 type(mpas_pool_type),pointer:: SS2G_aops
 type(mpas_pool_type),pointer:: SU2G_aops
 type(mpas_pool_type),pointer:: GOCART2G_aops

 type(block_type),pointer:: block

 logical,pointer:: do_CA2Gbc,do_CA2Gbr,do_CA2Goc
 logical,pointer:: do_NI2G,do_DU2G,do_SS2G,do_SU2G
 logical:: do_GOCART2G

 integer:: time_lev
 integer:: i,its,ite,j,jts,jte,k,kts,kte,n,nerod

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_driver:')

 time_lev = 1

 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_CA2Gbc',do_CA2Gbc)
 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_CA2Gbr',do_CA2Gbr)
 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_CA2Goc',do_CA2Goc)
 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_DU2G'  ,do_DU2G  )
 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_NI2G'  ,do_NI2G  )
 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_SS2G'  ,do_SS2G  )
 call mpas_pool_get_config(domain%configs,'config_gocart2G_do_SU2G'  ,do_SU2G  )


 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh                )
    call mpas_pool_get_subpool(block%structs,'state'               ,state               )
    call mpas_pool_get_subpool(block%structs,'diag'                ,diag                )
    call mpas_pool_get_subpool(block%structs,'diag_physics'        ,diag_physics        )
    call mpas_pool_get_subpool(block%structs,'sfc_input'           ,sfc_input           )
    call mpas_pool_get_subpool(block%structs,'CAMS_emissions'      ,CAMS_emissions      )
    call mpas_pool_get_subpool(block%structs,'gocart2G_backgrounds',gocart2G_backgrounds)
    call mpas_pool_get_subpool(block%structs,'gocart2G_met'        ,gocart2G_met        )

    call mpas_pool_get_subpool(block%structs,'CA2G_bc_diags' ,CA2G_bc_diags )
    call mpas_pool_get_subpool(block%structs,'CA2G_br_diags' ,CA2G_br_diags )
    call mpas_pool_get_subpool(block%structs,'CA2G_oc_diags' ,CA2G_oc_diags )
    call mpas_pool_get_subpool(block%structs,'DU2G_diags'    ,DU2G_diags    )
    call mpas_pool_get_subpool(block%structs,'NI2G_diags'    ,NI2G_diags    )
    call mpas_pool_get_subpool(block%structs,'SS2G_diags'    ,SS2G_diags    )
    call mpas_pool_get_subpool(block%structs,'SU2G_diags'    ,SU2G_diags    )
    call mpas_pool_get_subpool(block%structs,'GOCART2G_diags',GOCART2G_diags)

    call mpas_pool_get_subpool(block%structs,'CA2G_bc_aops' ,CA2G_bc_aops )
    call mpas_pool_get_subpool(block%structs,'CA2G_br_aops' ,CA2G_br_aops )
    call mpas_pool_get_subpool(block%structs,'CA2G_oc_aops' ,CA2G_oc_aops )
    call mpas_pool_get_subpool(block%structs,'DU2G_aops'    ,DU2G_aops    )
    call mpas_pool_get_subpool(block%structs,'NI2G_aops'    ,NI2G_aops    )
    call mpas_pool_get_subpool(block%structs,'SS2G_aops'    ,SS2G_aops    )
    call mpas_pool_get_subpool(block%structs,'SU2G_aops'    ,SU2G_aops    )
    call mpas_pool_get_subpool(block%structs,'GOCART2G_aops',GOCART2G_aops)


    !--- defines dimensions and allocate local arrays from MPAS needed to run GOCART2G:
    call mpas_gocart2G%gocart2G_dims(mesh)
    call mpas_gocart2G%gocart2G_allocate()

    its   = mpas_gocart2G%its
    ite   = mpas_gocart2G%ite
    jts   = mpas_gocart2G%jts
    jte   = mpas_gocart2G%jte
    kts   = mpas_gocart2G%kts
    kte   = mpas_gocart2G%kte
    nerod = mpas_gocart2G%nerod


    !--- fills local chemistry arrays with global chemistry arrays:
    call mpas_gocart2G%gocart2G_fromMPAS(CAMS_emissions,gocart2G_backgrounds,gocart2G_met,block%configs, &
                                         mesh,diag,state,diag_physics,sfc_input,time_lev)


    !--- CA2G_bc:
    if(do_CA2Gbc) then
       CA2G_bc_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       CA2G_bc%lats     => mpas_gocart2G%xlat      ; CA2G_bc%lons     => mpas_gocart2G%xlon
       CA2G_bc%area     => mpas_gocart2G%area      ; CA2G_bc%frocean  => mpas_gocart2G%frocean
       CA2G_bc%fraci    => mpas_gocart2G%frice     ; CA2G_bc%frlake   => mpas_gocart2G%frlake
       CA2G_bc%lwi      => mpas_gocart2G%lwi       ; CA2G_bc%u10m     => mpas_gocart2G%u10m
       CA2G_bc%v10m     => mpas_gocart2G%v10m      ; CA2G_bc%zpbl     => mpas_gocart2G%zpbl
       CA2G_bc%ustar    => mpas_gocart2G%ustar     ; CA2G_bc%sh       => mpas_gocart2G%sh
       CA2G_bc%z0h      => mpas_gocart2G%z0h       ; CA2G_bc%cn_prcp  => mpas_gocart2G%cn_prcp
       CA2G_bc%ncn_prcp => mpas_gocart2G%ncn_prcp  ; CA2G_bc%tropp    => mpas_gocart2G%backg_ptrop

       CA2G_bc%airdens  => mpas_gocart2G%airdens   ; CA2G_bc%delp     => mpas_gocart2G%delp
       CA2G_bc%delz     => mpas_gocart2G%delz      ; CA2G_bc%t        => mpas_gocart2G%t
       CA2G_bc%rh2      => mpas_gocart2G%rh2       ; CA2G_bc%zle      => mpas_gocart2G%zle
       CA2G_bc%ple      => mpas_gocart2G%ple       ; CA2G_bc%pfl_lsan => mpas_gocart2G%pfl_lsan
       CA2G_bc%pfi_lsan => mpas_gocart2G%pfi_lsan  ; CA2G_bc%u        => mpas_gocart2G%u
       CA2G_bc%v        => mpas_gocart2G%v

       !--- emissions:
!      CA2G_bc%bc_antebc1 => mpas_gocart2G%qbc1_em

       !--- chemistry fields:
       CA2G_bc%bcphobic => mpas_gocart2G%qbcphobic
       CA2G_bc%bcphilic => mpas_gocart2G%qbcphilic

       !--- gocart2G processes:
!      call CA2G_bc_params%emissions_GridComp(CA2G_bc,its,ite,jts,jte,kts,kte, &
!                                    iyear,imonth,iday,ihour,iminute,isecond)
       call CA2G_bc_params%processes_GridComp(CA2G_bc,its,ite,jts,jte,kts,kte)

       !--- global diagnostics:
       call CA2G_bc_diagnostics(mesh,CA2G_bc,CA2G_bc_diags,CA2G_bc_aops,its,ite,jts,jte,kts,kte)
    endif


    !--- CA2G_br:
    if(do_CA2Gbr) then
       CA2G_br_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       CA2G_br%lats     => mpas_gocart2G%xlat     ; CA2G_br%lons     => mpas_gocart2G%xlon
       CA2G_br%area     => mpas_gocart2G%area     ; CA2G_br%frocean  => mpas_gocart2G%frocean
       CA2G_br%fraci    => mpas_gocart2G%frice    ; CA2G_br%lwi      => mpas_gocart2G%lwi
       CA2G_br%u10m     => mpas_gocart2G%u10m     ; CA2G_br%v10m     => mpas_gocart2G%v10m
       CA2G_br%ustar    => mpas_gocart2G%ustar    ; CA2G_br%frlake   => mpas_gocart2G%frlake
       CA2G_br%zpbl     => mpas_gocart2G%zpbl     ; CA2G_br%sh       => mpas_gocart2G%sh
       CA2G_br%z0h      => mpas_gocart2G%z0h      ; CA2G_br%cn_prcp  => mpas_gocart2G%cn_prcp
       CA2G_br%ncn_prcp => mpas_gocart2G%ncn_prcp ; CA2G_br%tropp    => mpas_gocart2G%backg_ptrop

       CA2G_br%airdens  => mpas_gocart2G%airdens  ; CA2G_br%delp     => mpas_gocart2G%delp
       CA2G_br%delz     => mpas_gocart2G%delz     ; CA2G_br%t        => mpas_gocart2G%t
       CA2G_br%rh2      => mpas_gocart2G%rh2      ; CA2G_br%zle      => mpas_gocart2G%zle
       CA2G_br%ple      => mpas_gocart2G%ple      ; CA2G_br%pfl_lsan => mpas_gocart2G%pfl_lsan
       CA2G_br%pfi_lsan => mpas_gocart2G%pfi_lsan ; CA2G_br%u        => mpas_gocart2G%u
       CA2G_br%v        => mpas_gocart2G%v

       !--- chemistry fields:
       CA2G_br%brphobic => mpas_gocart2G%qbrphobic
       CA2G_br%brphilic => mpas_gocart2G%qbrphilic

       !--- gocart2G processes:
!      call CA2G_br_params%emissions_GridComp(CA2G_br,its,ite,jts,jte,kts,kte, &
!                                    iyear,imonth,iday,ihour,iminute,isecond)
       call CA2G_br_params%processes_GridComp(CA2G_br,its,ite,jts,jte,kts,kte)

       !--- global diagnostics:
       call CA2G_br_diagnostics(mesh,CA2G_br,CA2G_br_diags,CA2G_br_aops,its,ite,jts,jte,kts,kte)
    endif


    !--- CA2G_oc:
    if(do_CA2Goc) then
       CA2G_oc_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       CA2G_oc%lats     => mpas_gocart2G%xlat     ; CA2G_oc%lons     => mpas_gocart2G%xlon
       CA2G_oc%area     => mpas_gocart2G%area     ; CA2G_oc%frocean  => mpas_gocart2G%frocean
       CA2G_oc%fraci    => mpas_gocart2G%frice    ; CA2G_oc%frlake   => mpas_gocart2G%frlake
       CA2G_oc%lwi      => mpas_gocart2G%lwi      ; CA2G_oc%u10m     => mpas_gocart2G%u10m
       CA2G_oc%v10m     => mpas_gocart2G%v10m     ; CA2G_oc%zpbl     => mpas_gocart2G%zpbl
       CA2G_oc%ustar    => mpas_gocart2G%ustar    ; CA2G_oc%sh       => mpas_gocart2G%sh
       CA2G_oc%z0h      => mpas_gocart2G%z0h      ; CA2G_oc%cn_prcp  => mpas_gocart2G%cn_prcp
       CA2G_oc%ncn_prcp => mpas_gocart2G%ncn_prcp ; CA2G_oc%tropp    => mpas_gocart2G%backg_ptrop

       CA2G_oc%airdens  => mpas_gocart2G%airdens  ; CA2G_oc%delp     => mpas_gocart2G%delp
       CA2G_oc%delz     => mpas_gocart2G%delz     ; CA2G_oc%t        => mpas_gocart2G%t
       CA2G_oc%rh2      => mpas_gocart2G%rh2      ; CA2G_oc%zle      => mpas_gocart2G%zle
       CA2G_oc%ple      => mpas_gocart2G%ple      ; CA2G_oc%pfl_lsan => mpas_gocart2G%pfl_lsan
       CA2G_oc%pfi_lsan => mpas_gocart2G%pfi_lsan ; CA2G_oc%u        => mpas_gocart2G%u
       CA2G_oc%v        => mpas_gocart2G%v

       !--- emissions:
!      CA2G_oc%oc_ANTEOC1 => mpas_gocart2G%qoc1_em

       !--- chemistry fields:
       CA2G_oc%ocphobic => mpas_gocart2G%qocphobic
       CA2G_oc%ocphilic => mpas_gocart2G%qocphilic

       !--- gocart2G processes:
!      call CA2G_oc_params%emissions_GridComp(CA2G_oc,its,ite,jts,jte,kts,kte, &
!                                     iyear,imonth,iday,ihour,iminute,isecond)
       call CA2G_oc_params%processes_GridComp(CA2G_oc,its,ite,jts,jte,kts,kte)

       !--- global diagnostics:
       call CA2G_oc_diagnostics(mesh,CA2G_oc,CA2G_oc_diags,CA2G_oc_aops,its,ite,jts,jte,kts,kte)
    endif


    !--- DU2G:
    if(do_DU2G) then
       DU2G_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       DU2G%du_src     => mpas_gocart2G%erod     ; DU2G%frlake    => mpas_gocart2G%frlake
       DU2G%wet1       => mpas_gocart2G%wet1     ; DU2G%lwi       => mpas_gocart2G%lwi
       DU2G%area       => mpas_gocart2G%area     ; DU2G%ustar     => mpas_gocart2G%ustar
       DU2G%zpbl       => mpas_gocart2G%zpbl     ; DU2G%sh        => mpas_gocart2G%sh
       DU2G%z0h        => mpas_gocart2G%z0h      ; DU2G%u10m      => mpas_gocart2G%u10m
       DU2G%v10m       => mpas_gocart2G%v10m     ; DU2G%cn_prcp   => mpas_gocart2G%cn_prcp
       DU2G%ncn_prcp   => mpas_gocart2G%ncn_prcp ; DU2G%tropp     => mpas_gocart2G%backg_ptrop

       DU2G%airdens    => mpas_gocart2G%airdens  ; DU2G%delp     => mpas_gocart2G%delp
       DU2G%delz       => mpas_gocart2G%delz     ; DU2G%t        => mpas_gocart2G%t
       DU2G%rh2        => mpas_gocart2G%rh2      ; DU2G%zle      => mpas_gocart2G%zle
       DU2G%ple        => mpas_gocart2G%ple      ; DU2G%pfl_lsan => mpas_gocart2G%pfl_lsan
       DU2G%pfi_lsan   => mpas_gocart2G%pfi_lsan ; DU2G%u        => mpas_gocart2G%u
       DU2G%v          => mpas_gocart2G%v

       !--- chemistry fields:
       DU2G%du(:,:,:,1) = mpas_gocart2G%qdust1(:,:,:)
       DU2G%du(:,:,:,2) = mpas_gocart2G%qdust2(:,:,:)
       DU2G%du(:,:,:,3) = mpas_gocart2G%qdust3(:,:,:)
       DU2G%du(:,:,:,4) = mpas_gocart2G%qdust4(:,:,:)
       DU2G%du(:,:,:,5) = mpas_gocart2G%qdust5(:,:,:)

       !--- gocart2G processes:
       call DU2G_params%emissions_GridComp(DU2G,its,ite,jts,jte,kts,kte,nerod)
       call DU2G_params%processes_GridComp(DU2G,its,ite,jts,jte,kts,kte)

       !--- global diagnostics:
       call DU2G_diagnostics(mesh,DU2G,DU2G_diags,DU2G_aops,its,ite,jts,jte,kts,kte)

       mpas_gocart2G%qdust1(:,:,:) = DU2G%du(:,:,:,1)
       mpas_gocart2G%qdust2(:,:,:) = DU2G%du(:,:,:,2)
       mpas_gocart2G%qdust3(:,:,:) = DU2G%du(:,:,:,3)
       mpas_gocart2G%qdust4(:,:,:) = DU2G%du(:,:,:,4)
       mpas_gocart2G%qdust5(:,:,:) = DU2G%du(:,:,:,5)
    endif


    !--- NI2G:
    if(do_NI2G) then
       NI2G_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       NI2G%lwi      => mpas_gocart2G%lwi      ; NI2G%ustar    => mpas_gocart2G%ustar
       NI2G%zpbl     => mpas_gocart2G%zpbl     ; NI2G%sh       => mpas_gocart2G%sh
       NI2G%z0h      => mpas_gocart2G%z0h      ; NI2G%tropp    => mpas_gocart2G%backg_ptrop
       NI2G%cn_prcp  => mpas_gocart2G%cn_prcp  ; NI2G%ncn_prcp => mpas_gocart2G%ncn_prcp

       NI2G%airdens  => mpas_gocart2G%airdens  ; NI2G%delp     => mpas_gocart2G%delp
       NI2G%t        => mpas_gocart2G%t        ; NI2G%rh2      => mpas_gocart2G%rh2
       NI2G%u        => mpas_gocart2g%u        ; NI2G%v        => mpas_gocart2G%v
       NI2G%ple      => mpas_gocart2G%ple      ; NI2G%zle      => mpas_gocart2G%zle
       NI2G%pfl_lsan => mpas_gocart2G%pfl_lsan ; NI2G%pfi_lsan => mpas_gocart2G%pfi_lsan

       !--- emissions:
!       NI2G%emi_nh3_sum => mpas_gocart2G%qnh3_em

       !--- chemistry fields:
       NI2G%so4    => mpas_gocart2G%qso4  ; NI2G%nh3    => mpas_gocart2g%qnh3
       NI2G%nh4a   => mpas_gocart2G%qnh4a ; NI2G%no3an1 => mpas_gocart2G%qni1
       NI2G%no3an2 => mpas_gocart2g%qni2  ; NI2G%no3an3 => mpas_gocart2G%qni3

       NI2G%xhno3  => mpas_gocart2G%backg_hno3

       NI2G%du(:,:,:,1) = mpas_gocart2G%qdust1(:,:,:)
       NI2G%du(:,:,:,2) = mpas_gocart2G%qdust2(:,:,:)
       NI2G%du(:,:,:,3) = mpas_gocart2G%qdust3(:,:,:)
       NI2G%du(:,:,:,4) = mpas_gocart2G%qdust4(:,:,:)
       NI2G%du(:,:,:,5) = mpas_gocart2G%qdust5(:,:,:)

       NI2G%ss(:,:,:,1) = mpas_gocart2G%qseas1(:,:,:)
       NI2G%ss(:,:,:,2) = mpas_gocart2G%qseas2(:,:,:)
       NI2G%ss(:,:,:,3) = mpas_gocart2G%qseas3(:,:,:)
       NI2G%ss(:,:,:,4) = mpas_gocart2G%qseas4(:,:,:)
       NI2G%ss(:,:,:,5) = mpas_gocart2G%qseas5(:,:,:)

       !--- gocart2G processes:
!      call NI2G_params%emissions_GridComp(NI2G,its,ite,jts,jte,kts,kte)
       call NI2G_params%processes_GridComp(NI2G,its,ite,jts,jte,kts,kte)

       !--- global diagnostics:
       call NI2G_diagnostics(mesh,NI2G,NI2G_diags,NI2G_aops,its,ite,jts,jte,kts,kte)
    endif


    !--- SS2G:
    if(do_SS2G) then
       SS2G_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       SS2G%frocean         => mpas_gocart2G%frocean     ; SS2G%fraci    => mpas_gocart2G%frice
       SS2G%frlake          => mpas_gocart2G%frlake      ; SS2G%area     => mpas_gocart2G%area
       SS2G%lwi             => mpas_gocart2G%lwi         ; SS2G%u10m     => mpas_gocart2G%u10m
       SS2G%v10m            => mpas_gocart2G%v10m        ; SS2G%ustar    => mpas_gocart2G%ustar
       SS2g%ts              => mpas_gocart2G%ts          ; SS2G%zpbl     => mpas_gocart2G%zpbl
       SS2G%sh              => mpas_gocart2G%sh          ; SS2G%z0h      => mpas_gocart2G%z0h
       SS2G%cn_prcp         => mpas_gocart2G%cn_prcp     ; SS2G%ncn_prcp => mpas_gocart2G%ncn_prcp
       SS2G%deep_lakes_mask => mpas_gocart2G%dlakes_mask ; SS2G%tropp    => mpas_gocart2G%backg_ptrop

       SS2G%airdens  => mpas_gocart2G%airdens  ; SS2G%delp     => mpas_gocart2G%delp
       SS2G%delz     => mpas_gocart2G%delz     ; SS2G%t        => mpas_gocart2G%t
       SS2G%rh2      => mpas_gocart2G%rh2      ; SS2G%zle      => mpas_gocart2G%zle
       SS2G%ple      => mpas_gocart2G%ple      ; SS2G%pfl_lsan => mpas_gocart2G%pfl_lsan
       SS2G%pfi_lsan => mpas_gocart2G%pfi_lsan ; SS2G%u        => mpas_gocart2G%u
       SS2G%v        => mpas_gocart2G%v

       !--- chemistry fields:
       SS2G%ss(:,:,:,1) = mpas_gocart2G%qseas1(:,:,:)
       SS2G%ss(:,:,:,2) = mpas_gocart2G%qseas2(:,:,:)
       SS2G%ss(:,:,:,3) = mpas_gocart2G%qseas3(:,:,:)
       SS2G%ss(:,:,:,4) = mpas_gocart2G%qseas4(:,:,:)
       SS2G%ss(:,:,:,5) = mpas_gocart2G%qseas5(:,:,:)

       !--- gocart2G processes:
       call SS2G_params%emissions_GridComp(SS2G,its,ite,jts,jte,kts,kte)
       call SS2G_params%processes_GridComp(SS2G,its,ite,jts,jte,kts,kte)

       !--- global diagnostics:
       call SS2G_diagnostics(mesh,SS2G,SS2G_diags,SS2G_aops,its,ite,jts,jte,kts,kte)

       mpas_gocart2G%qseas1(:,:,:) = SS2G%ss(:,:,:,1)
       mpas_gocart2G%qseas2(:,:,:) = SS2G%ss(:,:,:,2)
       mpas_gocart2G%qseas3(:,:,:) = SS2G%ss(:,:,:,3)
       mpas_gocart2G%qseas4(:,:,:) = SS2G%ss(:,:,:,4)
       mpas_gocart2G%qseas5(:,:,:) = SS2G%ss(:,:,:,5)
    endif


    !--- SU2G:
    if(do_SU2G) then
       SU2G_params%cdt = mpas_gocart2G%dt

       !--- meteorological fields:
       SU2G%lats    => mpas_gocart2G%xlat        ; SU2G%lons     => mpas_gocart2G%xlon
       SU2G%area    => mpas_gocart2G%area        ; SU2G%coszr    => mpas_gocart2G%coszr
       SU2G%frocean => mpas_gocart2G%frocean     ; SU2G%lwi      => mpas_gocart2G%lwi
       SU2G%u10m    => mpas_gocart2G%u10m        ; SU2G%v10m     => mpas_gocart2G%v10m
       SU2G%zpbl    => mpas_gocart2G%zpbl        ; SU2G%ustar    => mpas_gocart2G%ustar
       SU2G%sh      => mpas_gocart2G%sh          ; SU2G%z0h      => mpas_gocart2G%z0h
       SU2G%cn_prcp => mpas_gocart2G%cn_prcp     ; SU2G%ncn_prcp => mpas_gocart2G%ncn_prcp
       SU2G%tropp   => mpas_gocart2G%backg_ptrop

       SU2G%airdens  => mpas_gocart2G%airdens  ; SU2G%delp     => mpas_gocart2G%delp
       SU2G%delz     => mpas_gocart2G%delz     ; SU2G%t        => mpas_gocart2G%t
       SU2G%rh2      => mpas_gocart2G%rh2      ; SU2G%zle      => mpas_gocart2G%zle
       SU2G%ple      => mpas_gocart2G%ple      ; SU2G%pfl_lsan => mpas_gocart2G%pfl_lsan
       SU2G%pfi_lsan => mpas_gocart2G%pfi_lsan ; SU2G%u        => mpas_gocart2G%u
       SU2G%v        => mpas_gocart2G%v        ; SU2G%fcld     => mpas_gocart2G%fcld

       !--- emissions:
!      SU2G%su_anthrol1 => mpas_gocart2G%qso2_em

       !--- chemistry fields and climatological background fields:
       SU2G%dms     => mpas_gocart2G%qdms      ; SU2G%so2     => mpas_gocart2G%qso2
       SU2G%so4     => mpas_gocart2G%qso4      ; SU2G%msa     => mpas_gocart2G%qmsa
       SU2G%su_dmso => mpas_gocart2g%backg_dms ; SU2G%su_no3  => mpas_gocart2G%backg_no3
       SU2G%su_oh   => mpas_gocart2G%backg_oh  ; SU2G%su_h2o2 => mpas_gocart2G%backg_h2o2

       !--- gocart2G processes:
!      call SU2G_params%emissions_GridComp(SU2G,its,ite,jts,jte,kts,kte,iyear,imonth,iday,ihour,iminute,isecond)
       call SU2G_params%processes_GridComp(SU2G,its,ite,jts,jte,kts,kte,iyear,imonth,iday,ihour,iminute,isecond)

       !--- global diagnostics:
       call SU2G_diagnostics(mesh,SU2G,SU2G_diags,SU2G_aops,its,ite,jts,jte,kts,kte)

    endif


    !--- GOCART2G global diagnostics:
    do_GOCART2G = .false.
    if(do_CA2Gbc .or. do_CA2Gbr .or. do_CA2Goc .or. do_DU2G .or. &
       do_NI2G .or. do_NI2G .or. do_SS2G .or. do_SU2G) do_GOCART2G = .true.
    if(do_GOCART2G) then
       call GOCART2G_params%processes_GridComp(GOCART2G,CA2G_bc,CA2G_br,CA2G_oc,DU2G,NI2G,SS2G,SU2G, &
                                               its,ite,jts,jte,kts,kte)
       call GOCART2G_diagnostics(mesh,GOCART2G,GOCART2G_diags,GOCART2G_aops,its,ite,jts,jte,kts,kte)
    endif


    !--- fills global chemistry arrays with local chemistry arrays:
    call mpas_gocart2G%gocart2G_toMPAS(state,time_lev)

    block => block % next
 end do


 call mpas_log_write('--- end subroutine gocart2G_driver:')
 call mpas_log_write(' ')

 end subroutine gocart2G_driver

!==================================================================================================================
 end module mpas_chemistry_gocart2G_driver
!==================================================================================================================

