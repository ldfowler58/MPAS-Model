!=================================================================================================================
 module mpas_chemistry_gocart2G_driver
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_subpool

 use mpas_chemistry_gocart2G_interface
 use mpas_chemistry_gocart2G_vars,only: mpas_gocart2G,          &
                                        CA2G_bc,CA2G_bc_params, &
                                        CA2G_br,CA2G_br_params, &
                                        CA2G_oc,CA2G_oc_params, &
                                        SU2G,SU2G_params


 implicit none
 private
 public:: gocart2G_driver


 contains


!=================================================================================================================
 subroutine gocart2G_driver(domain,itimestep,xtime_s)
!=================================================================================================================

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
 type(mpas_pool_type),pointer:: GOCART_backgrounds
 type(block_type),pointer:: block

 integer:: time_lev
 integer:: i,its,ite,j,jts,jte,k,kts,kte,ktep1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_driver:')

 time_lev = 1


 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'              ,mesh              )
    call mpas_pool_get_subpool(block%structs,'state'             ,state             )
    call mpas_pool_get_subpool(block%structs,'diag'              ,diag              )
    call mpas_pool_get_subpool(block%structs,'diag_physics'      ,diag_physics      )
    call mpas_pool_get_subpool(block%structs,'sfc_input'         ,sfc_input         )
    call mpas_pool_get_subpool(block%structs,'CAMS_emissions'    ,CAMS_emissions    )
    call mpas_pool_get_subpool(block%structs,'GOCART_backgrounds',GOCART_backgrounds)


    !--- defines dimensions and allocate local arrays from MPAS needed to run GOCART2G:
    call mpas_gocart2G%gocart2G_dims(mesh)
    call mpas_gocart2G%gocart2G_allocate()

    its = mpas_gocart2G%its
    ite = mpas_gocart2G%ite
    jts = mpas_gocart2G%jts
    jte = mpas_gocart2G%jte
    kts = mpas_gocart2G%kts
    kte = mpas_gocart2G%kte


    !--- fills local chemistry arrays with global chemistry arrays:
    call mpas_gocart2G%gocart2G_fromMPAS(CAMS_emissions,GOCART_backgrounds,block%configs,mesh,diag, &
                                         state,diag_physics,sfc_input,time_lev)

    !--- CA2G_bc:
    CA2G_bc_params%cdt = mpas_gocart2G%dt

    CA2G_bc%area     => mpas_gocart2G%area

    CA2G_bc%frocean  => mpas_gocart2G%frocean   ; CA2G_bc%fraci    => mpas_gocart2G%fraci
    CA2G_bc%lwi      => mpas_gocart2G%lwi       ; CA2G_bc%tropp    => mpas_gocart2G%tropp
    CA2G_bc%u10m     => mpas_gocart2G%u10m      ; CA2G_bc%v10m     => mpas_gocart2G%v10m
    CA2G_bc%ustar    => mpas_gocart2G%ustar     ; CA2G_bc%frlake   => mpas_gocart2G%frlake
    CA2G_bc%zpbl     => mpas_gocart2G%zpbl      ; CA2G_bc%sh       => mpas_gocart2G%sh
    CA2G_bc%z0h      => mpas_gocart2G%z0h       ; CA2G_bc%cn_prcp  => mpas_gocart2G%cn_prcp
    CA2G_bc%ncn_prcp => mpas_gocart2G%ncn_prcp

    CA2G_bc%airdens  => mpas_gocart2G%airdens   ; CA2G_bc%delp     => mpas_gocart2G%delp
    CA2G_bc%delz     => mpas_gocart2G%delz      ; CA2G_bc%t        => mpas_gocart2G%t
    CA2G_bc%rh2      => mpas_gocart2G%rh2       ; CA2G_bc%zle      => mpas_gocart2G%zle
    CA2G_bc%ple      => mpas_gocart2G%ple       ; CA2G_bc%pfl_lsan => mpas_gocart2G%pfl_lsan
    CA2G_bc%pfi_lsan => mpas_gocart2G%pfi_lsan  ; CA2G_bc%u        => mpas_gocart2G%u
    CA2G_bc%v        => mpas_gocart2G%v

    CA2G_bc%bcPHOBIC => mpas_gocart2G%qbcphobic ; CA2G_bc%bcPHILIC => mpas_gocart2G%qbcphilic

    call CA2G_bc_params%run2_GridComp(CA2G_bc,its,ite,jts,jte,kts,kte)


    !--- CA2G_br:
    CA2G_br_params%cdt = mpas_gocart2G%dt

    CA2G_br%area     => mpas_gocart2G%area

    CA2G_br%frocean  => mpas_gocart2G%frocean   ; CA2G_br%fraci    => mpas_gocart2G%fraci
    CA2G_br%lwi      => mpas_gocart2G%lwi       ; CA2G_br%tropp    => mpas_gocart2G%tropp
    CA2G_br%u10m     => mpas_gocart2G%u10m      ; CA2G_br%v10m     => mpas_gocart2G%v10m
    CA2G_br%ustar    => mpas_gocart2G%ustar     ; CA2G_br%frlake   => mpas_gocart2G%frlake
    CA2G_br%zpbl     => mpas_gocart2G%zpbl      ; CA2G_br%sh       => mpas_gocart2G%sh
    CA2G_br%z0h      => mpas_gocart2G%z0h       ; CA2G_br%cn_prcp  => mpas_gocart2G%cn_prcp
    CA2G_br%ncn_prcp => mpas_gocart2G%ncn_prcp

    CA2G_br%airdens  => mpas_gocart2G%airdens   ; CA2G_br%delp     => mpas_gocart2G%delp
    CA2G_br%delz     => mpas_gocart2G%delz      ; CA2G_br%t        => mpas_gocart2G%t
    CA2G_br%rh2      => mpas_gocart2G%rh2       ; CA2G_br%zle      => mpas_gocart2G%zle
    CA2G_br%ple      => mpas_gocart2G%ple       ; CA2G_br%pfl_lsan => mpas_gocart2G%pfl_lsan
    CA2G_br%pfi_lsan => mpas_gocart2G%pfi_lsan  ; CA2G_br%u        => mpas_gocart2G%u
    CA2G_br%v        => mpas_gocart2G%v

    CA2G_br%brPHOBIC => mpas_gocart2G%qbrphobic ; CA2G_br%brPHILIC => mpas_gocart2G%qbrphilic


    !--- CA2G_oc:
    CA2G_oc_params%cdt = mpas_gocart2G%dt

    CA2G_oc%area     => mpas_gocart2G%area

    CA2G_oc%frocean  => mpas_gocart2G%frocean   ; CA2G_oc%fraci    => mpas_gocart2G%fraci
    CA2G_oc%lwi      => mpas_gocart2G%lwi       ; CA2G_oc%tropp    => mpas_gocart2G%tropp
    CA2G_oc%u10m     => mpas_gocart2G%u10m      ; CA2G_oc%v10m     => mpas_gocart2G%v10m
    CA2G_oc%ustar    => mpas_gocart2G%ustar     ; CA2G_oc%frlake   => mpas_gocart2G%frlake
    CA2G_oc%zpbl     => mpas_gocart2G%zpbl      ; CA2G_oc%sh       => mpas_gocart2G%sh
    CA2G_oc%z0h      => mpas_gocart2G%z0h       ; CA2G_oc%cn_prcp  => mpas_gocart2G%cn_prcp
    CA2G_oc%ncn_prcp => mpas_gocart2G%ncn_prcp

    CA2G_oc%airdens  => mpas_gocart2G%airdens   ; CA2G_oc%delp     => mpas_gocart2G%delp
    CA2G_oc%delz     => mpas_gocart2G%delz      ; CA2G_oc%t        => mpas_gocart2G%t
    CA2G_oc%rh2      => mpas_gocart2G%rh2       ; CA2G_oc%zle      => mpas_gocart2G%zle
    CA2G_oc%ple      => mpas_gocart2G%ple       ; CA2G_oc%pfl_lsan => mpas_gocart2G%pfl_lsan
    CA2G_oc%pfi_lsan => mpas_gocart2G%pfi_lsan  ; CA2G_oc%u        => mpas_gocart2G%u
    CA2G_oc%v        => mpas_gocart2G%v

    CA2G_oc%ocPHOBIC => mpas_gocart2G%qocphobic ; CA2G_oc%ocPHILIC => mpas_gocart2G%qocphilic


    !--- SU2G:
    SU2G_params%cdt = mpas_gocart2G%dt

    SU2G%lats     => mpas_gocart2G%xlat       ; SU2G%lons     => mpas_gocart2G%xlon
    SU2G%area     => mpas_gocart2G%area

    SU2G%frocean  => mpas_gocart2G%frocean    ; SU2G%lwi      => mpas_gocart2G%lwi
    SU2G%tropp    => mpas_gocart2G%tropp      ; SU2G%u10m     => mpas_gocart2G%u10m
    SU2G%v10m     => mpas_gocart2G%v10m       ; SU2G%zpbl     => mpas_gocart2G%zpbl
    SU2G%ustar    => mpas_gocart2G%ustar      ; SU2G%sh       => mpas_gocart2G%sh
    SU2G%z0h      => mpas_gocart2G%z0h        ; SU2G%cn_prcp  => mpas_gocart2G%cn_prcp
    SU2G%ncn_prcp => mpas_gocart2G%ncn_prcp

    SU2G%airdens  => mpas_gocart2G%airdens    ; SU2G%delp     => mpas_gocart2G%delp
    SU2G%delz     => mpas_gocart2G%delz       ; SU2G%t        => mpas_gocart2G%t
    SU2G%rh2      => mpas_gocart2G%rh2        ; SU2G%zle      => mpas_gocart2G%zle
    SU2G%ple      => mpas_gocart2G%ple        ; SU2G%pfl_lsan => mpas_gocart2G%pfl_lsan
    SU2G%pfi_lsan => mpas_gocart2G%pfi_lsan   ; SU2G%u        => mpas_gocart2G%u
    SU2G%v        => mpas_gocart2G%v          ; SU2G%fcld     => mpas_gocart2G%fcld

    SU2G%DMS      => mpas_gocart2G%qdms       ; SU2G%SO2      => mpas_gocart2G%qso2
    SU2G%SO4      => mpas_gocart2G%qso4       ; SU2G%MSA      => mpas_gocart2G%qmsa 
    SU2G%SU_NO3   => mpas_gocart2G%backg_no3  ; SU2G%SU_OH    => mpas_gocart2G%backg_oh
    SU2G%SU_H2O2  => mpas_gocart2G%backg_h2o2

    call SU2G_params%run2_GridComp(SU2G,its,ite,jts,jte,kts,kte)

    block => block % next
 end do


 call mpas_log_write('--- end subroutine gocart2G_driver:')
 call mpas_log_write(' ')

 end subroutine gocart2G_driver

!=================================================================================================================
 end module mpas_chemistry_gocart2G_driver
!=================================================================================================================

