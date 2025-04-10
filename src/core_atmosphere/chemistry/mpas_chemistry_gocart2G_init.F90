! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_init
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use GOCART2G_instance
 use GOCART2G_MieMod_smiol
 use mpas_chemistry_gocart2G_vars


 implicit none
 private
 public:: init_gocart2G_chemistry


 contains


!==================================================================================================================
 subroutine init_gocart2G_chemistry(dminfo,configs,mesh,state)
!==================================================================================================================

!--- input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

!--- local variables:
 character(len=StrKIND):: fnm,message
 character(len=StrKIND),pointer:: fCA2G_bc,fCA2G_bc_RRTMG,fCA2G_br,fCA2G_br_RRTMG,fCA2G_oc,fCA2G_oc_RRTMG, &
                                  fDU2G,fDU2G_RRTMG,fNI2G,fNI2G_RRTMG,fSS2G,fSS2G_RRTMG,fSU2G,fSU2G_RRTMG

 logical:: l_exist
 logical,pointer:: do_CA2Gbc,do_CA2Gbr,do_CA2Goc
 logical,pointer:: do_NI2G,do_DU2G,do_SS2G,do_SU2G
 logical,pointer:: do_SOA2G
 logical:: do_GOCART2G

 logical,pointer:: to_MYNN

 integer:: its,ite,jts,jte,kts,kte,nerod
 integer:: ic,ch_size
 real(kind=RKIND),dimension(:),allocatable:: channels

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_gocart2G_chemistry:')

 call mpas_pool_get_config(configs,'config_gocart2G_opticsBC',fCA2G_bc)
 call mpas_pool_get_config(configs,'config_gocart2G_opticsBR',fCA2G_br)
 call mpas_pool_get_config(configs,'config_gocart2G_opticsOC',fCA2G_oc)
 call mpas_pool_get_config(configs,'config_gocart2G_opticsDU',fDU2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsNI',fNI2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsSS',fSS2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsSU',fSU2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsBC_RRTMG',fCA2G_bc_RRTMG)
 call mpas_pool_get_config(configs,'config_gocart2G_opticsBR_RRTMG',fCA2G_br_RRTMG)
 call mpas_pool_get_config(configs,'config_gocart2G_opticsOC_RRTMG',fCA2G_oc_RRTMG)
 call mpas_pool_get_config(configs,'config_gocart2G_opticsDU_RRTMG',fDU2G_RRTMG   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsNI_RRTMG',fNI2G_RRTMG   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsSS_RRTMG',fSS2G_RRTMG   )
 call mpas_pool_get_config(configs,'config_gocart2G_opticsSU_RRTMG',fSU2G_RRTMG   )

 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbc',do_CA2Gbc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbr',do_CA2Gbr)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Goc',do_CA2Goc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_DU2G'  ,do_DU2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_NI2G'  ,do_NI2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SS2G'  ,do_SS2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SU2G'  ,do_SU2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SOA2G' ,do_SOA2G )

 call mpas_pool_get_config(configs,'config_gocart2G_toMYNN'   ,to_MYNN  )


!--- reads input wavelengths from LUT:
 ch_size = size(aerosol_monochromatic_optics_wavelength_in_nm_from_LUT)
 allocate(channels(ch_size))

!call mpas_log_write('--- read input channels from LUT:')
 do ic = 1,ch_size
    channels(ic) = aerosol_monochromatic_optics_wavelength_in_nm_from_LUT(ic)
    channels(ic) = channels(ic)*1.e-9
!   call mpas_log_write('$i $r',intArgs=(/ic/),realArgs=(/channels(ic)/))
 enddo
!call mpas_log_write('--- end input channels from LUT.')


!--- initializes dimensions used to run the GOCART-2G chemistry:
 call mpas_gocart2G%gocart2G_dims(mesh)

 its   = mpas_gocart2G%its
 ite   = mpas_gocart2G%ite
 jts   = mpas_gocart2G%jts
 jte   = mpas_gocart2G%jte
 kts   = mpas_gocart2G%kts
 kte   = mpas_gocart2G%kte
 nerod = mpas_gocart2G%nerod


!--- CA2G_bc:
 call mpas_log_write(' ')
 call mpas_log_write('--- begin initialization of CA2G_bc:')

!initializes and allocates all parameters and arrays related to CA2G_bc:
 call CA2G_bc_params%load_GridComp(kts,kte)
 call CA2G_bc%gocart2G_allocate(its,ite,jts,jte,kts,kte)

 if(do_CA2Gbc) then
    !creates radiation Mie table for CA2G:
    l_exist = .false.
    fnm = trim(fCA2G_bc_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       CA2G_bc_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for CA2G:
    l_exist = .false.
    fnm = trim(fCA2G_bc)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       CA2G_bc_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of CA2G_bc.')
    call mpas_log_write(' ')
 endif


!--- CA2G_br:
 call mpas_log_write('--- begin initialization of CA2G_br:')

!initializes and allocates all parameters and arrays related to CA2G_br:
 call CA2G_br_params%load_GridComp(kts,kte)
 call CA2G_br%gocart2G_allocate(its,ite,jts,jte,kts,kte)

 if(do_CA2Gbr) then
    !creates radiation Mie table for CA2G_br:
    l_exist = .false.
    fnm = trim(fCA2G_br_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       CA2G_br_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for CA2G_br:
    l_exist = .false.
    fnm = trim(fCA2G_br)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       CA2G_br_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of CA2G_br.')
    call mpas_log_write(' ')
 endif


!--- CA2G_oc:
 call mpas_log_write('--- begin initialization of CA2G_oc:')

!initializes and allocates all parameters and arrays related to CA2G_oc:
 call CA2G_oc_params%load_GridComp(kts,kte)
 call CA2G_oc%gocart2G_allocate(its,ite,jts,jte,kts,kte)

 if(do_CA2Goc) then
    !creates radiation Mie table for CA2G_oc:
    l_exist = .false.
    fnm = trim(fCA2G_oc_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       CA2G_oc_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for CA2G_oc:
    l_exist = .false.
    fnm = trim(fCA2G_oc)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       CA2G_oc_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of CA2G_oc.')
    call mpas_log_write(' ')
 endif


!--- DU2G:
 call mpas_log_write('--- begin initialization of DU2G:')

!initializes and allocates all parameters and arrays related to DU2G:
 call DU2G_params%load_GridComp(kts,kte)
 call DU2G%gocart2G_allocate(its,ite,jts,jte,kts,kte,nerod)

 if(do_DU2G) then
    !create radiation Mie table for DU2G:
    l_exist = .false.
    fnm = trim(fDU2G_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       DU2G_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for DU2G:
    l_exist = .false.
    fnm = trim(fDU2G)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       DU2G_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of DU2G.')
    call mpas_log_write(' ')
 endif


!--- NI2G:
 call mpas_log_write('--- begin initialization of NI2G:')

!initializes and allocates all parameters and arrays related to NI2G:
 call DU2G_params%load_GridComp(kts,kte)
 call SS2G_params%load_GridComp(kts,kte)
 call NI2G_params%load_GridComp(DU2G_params,SS2G_params,kts,kte)
 call NI2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

 if(do_NI2G) then
    !create radiation Mie table for SU2G:
    l_exist = .false.
    fnm = trim(fNI2G_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       NI2G_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for SU2G:
    l_exist = .false.
    fnm = trim(fNI2G)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       NI2G_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of NI2G.')
    call mpas_log_write(' ')
 endif


!--- SS2G:
 call mpas_log_write('--- begin initialization of SS2G:')

!initializes and allocates all parameters and arrays related to SS2G:
 call SS2G_params%load_GridComp(kts,kte)
 call SS2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

 if(do_SS2G) then
    !create radiation Mie table for SU2G:
    l_exist = .false.
    fnm = trim(fSS2G_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       SS2G_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for SU2G:
    l_exist = .false.
    fnm = trim(fSS2G)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       SS2G_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of SS2G.')
    call mpas_log_write(' ')
 endif


!--- SU2G:
 call mpas_log_write('--- begin initialization of SU2G:')

!initializes and allocates all parameters and arrays related to SU2G:
 call SU2G_params%load_GridComp(kts,kte)
 call SU2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

 if(do_SU2G) then
    !create radiation Mie table for SU2G:
    l_exist = .false.
    fnm = trim(fSU2G_RRTMG)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       SU2G_params%rad_Mie = GOCART2G_Mie(dminfo,fnm)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    !create diagnostics Mie table for SU2G:
    l_exist = .false.
    fnm = trim(fSU2G)
    inquire(file=fnm,exist=l_exist)
    if(l_exist) then
       SU2G_params%diag_Mie = GOCART2G_Mie(dminfo,fnm,channels)
    else
       message = '--- file ''' //trim(fnm) //''' not in run directory.'
       call mpas_log_write(message,messageType=MPAS_LOG_CRIT)
    endif

    call mpas_log_write('--- end initialization of SU2G.')
    call mpas_log_write(' ')
 endif


!--- SOA2G:
 if(do_SOA2G .or. do_CA2Gbr .or. do_CA2Goc) then
    call mpas_log_write('--- begin initialization of SOA2G:')

    !initializes and allocates all parameters and arrasy related to SOA2G:
    call SOA2G_params%load_GridComp( )
    call SOA2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    call mpas_log_write('--- end initialization of SOA2G.')
 endif


!--- GOCART2G:
 do_GOCART2G = .false.
 if(do_CA2Gbc .or. do_CA2Gbr .or. do_CA2Goc .or. do_DU2G .or. &
    do_NI2G .or. do_NI2G .or. do_SS2G .or. do_SU2G) do_GOCART2G = .true.
 if(do_GOCART2G) then
    call GOCART2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)
 endif


!--- FEEBACKS TO PHYSICS:
 if(to_MYNN) then
    call mpas_chem_gocart2G%gocart2G_dims(mesh,state)
    call mpas_chem_gocart2G%gocart2G_allocate()
    call mpas_chem_gocart2G%gocart2G_tophysics_init(CA2G_bc_params,CA2G_br_params,CA2G_oc_params,DU2G_params, &
                                                    NI2G_params,SS2G_params,SU2G_params)
 endif


 call mpas_log_write('--- end subroutine init_gocart2G_chemistry.')
 call mpas_log_write(' ')

 end subroutine init_gocart2G_chemistry

!==================================================================================================================
 end module mpas_chemistry_gocart2G_init
!==================================================================================================================
