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
 subroutine init_gocart2G_chemistry(dminfo,configs,mesh)
!==================================================================================================================

!--- input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!--- local variables:
 character(len=StrKIND):: fnm

 logical,pointer:: do_CA2Gbc,do_CA2Gbr,do_CA2Goc
 logical,pointer:: do_NI2G,do_DU2G,do_SS2G,do_SU2G

 integer:: its,ite,jts,jte,kts,kte,nerod
 integer:: ic,ch_size
 real(kind=RKIND),dimension(:),allocatable:: channels

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_gocart2G_chemistry:')

 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbc',do_CA2Gbc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbr',do_CA2Gbr)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Goc',do_CA2Goc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_DU2G'  ,do_DU2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_NI2G'  ,do_NI2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SS2G'  ,do_SS2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SU2G'  ,do_SU2G  )


!--- reads input wavelengths from LUT:
 ch_size = size(aerosol_monochromatic_optics_wavelength_in_nm_from_LUT)
 allocate(channels(ch_size))

 call mpas_log_write(' ')
 call mpas_log_write('--- read input channels from LUT:')
 do ic = 1,ch_size
    channels(ic) = aerosol_monochromatic_optics_wavelength_in_nm_from_LUT(ic)
    channels(ic) = channels(ic)*1.e-9
    call mpas_log_write('$i $r',intArgs=(/ic/),realArgs=(/channels(ic)/))
 enddo
 call mpas_log_write('--- end input channels from LUT:')


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
 if(do_CA2Gbc) then
    call mpas_log_write('--- begin initialization of CA2G_bc:')

    !initializes and allocates all parameters and arrays related to CA2G_bc:
    call CA2G_bc_params%load_GridComp(kts,kte)
    call CA2G_bc%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    !creates radiation Mie table for CA2G:
    fnm = 'opticsBands_BC.v1_3.RRTMG.nc'
    CA2G_bc_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for CA2G:
    fnm = 'optics_BC.v1_3.nc'
    CA2G_bc_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of CA2G_bc:')
    call mpas_log_write(' ')
 endif


!--- CA2G_br:
 if(do_CA2Gbr) then
    call mpas_log_write('--- begin initialization of CA2G_br:')

    !initializes and allocates all parameters and arrays related to CA2G_br:
    call CA2G_br_params%load_GridComp(kts,kte)
    call CA2G_br%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    !creates radiation Mie table for CA2G_br:
    fnm = 'opticsBands_BRC.v1_5.RRTMG.nc'
    CA2G_br_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for CA2G_br:
    fnm = 'optics_BRC.v1_5.nc'
    CA2G_br_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of CA2G_br:')
    call mpas_log_write(' ')
 endif


!--- CA2G_oc:
 if(do_CA2Goc) then
    call mpas_log_write('--- begin initialization of CA2G_oc:')

    !initializes and allocates all parameters and arrays related to CA2G_oc:
    call CA2G_oc_params%load_GridComp(kts,kte)
    call CA2G_oc%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    !creates radiation Mie table for CA2G_oc:
    fnm = 'opticsBands_OC.v1_3.RRTMG.nc'
    CA2G_oc_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for CA2G_oc:
    fnm = 'optics_OC.v1_3.nc'
    CA2G_oc_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of CA2G_oc:')
    call mpas_log_write(' ')
 endif


!--- DU2G:
 if(do_DU2G) then
    call mpas_log_write('--- begin initialization of DU2G:')

    !initializes and allocates all parameters and arrays related to DU2G:
    call DU2G_params%load_GridComp(kts,kte)
    call DU2G%gocart2G_allocate(its,ite,jts,jte,kts,kte,nerod)

    !create radiation Mie table for DU2G:
    fnm = 'opticsBands_DU.v15_3.RRTMG.nc'
    DU2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for DU2G:
    !fnm = 'optics_DU.v15_3.nc'
    !DU2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of DU2G:')
    call mpas_log_write(' ')
 endif


!--- NI2G:
 if(do_NI2G) then
    call mpas_log_write('--- begin initialization of NI2G:')

    !initializes and allocates all parameters and arrays related to NI2G:
    call NI2G_params%load_GridComp(kts,kte)
    call NI2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    !create radiation Mie table for SU2G:
    fnm = 'opticsBands_NI.v2_5.RRTMG.nc'
    NI2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for SU2G:
    fnm = 'optics_NI.v2_5.nc'
    NI2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of NI2G:')
    call mpas_log_write(' ')
 endif


!--- SS2G:
 if(do_SS2G) then
    call mpas_log_write('--- begin initialization of SS2G:')

    !initializes and allocates all parameters and arrays related to SS2G:
    call SS2G_params%load_GridComp(kts,kte)
    call SS2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    !create radiation Mie table for SU2G:
    fnm = 'opticsBands_SS.v3_3.RRTMG.nc'
    SS2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for SU2G:
    !fnm = 'optics_SS.v3_3.nc'
    !SS2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of SS2G:')
    call mpas_log_write(' ')
 endif


!--- SU2G:
 if(do_SU2G) then
    call mpas_log_write('--- begin initialization of SU2G:')

    !initializes and allocates all parameters and arrays related to SU2G:
    call SU2G_params%load_GridComp(kts,kte)
    call SU2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

    !create radiation Mie table for SU2G:
    fnm = 'opticsBands_SU.v1_3.RRTMG.nc'
    SU2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

    !create diagnostics Mie table for SU2G:
    fnm = 'optics_SU.v1_3.nc'
    SU2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(fnm),channels)

    call mpas_log_write('--- end initialization of SU2G:')
    call mpas_log_write(' ')
 endif


 call mpas_log_write('--- end subroutine init_gocart2G_chemistry:')
 call mpas_log_write(' ')

 end subroutine init_gocart2G_chemistry

!==================================================================================================================
 end module mpas_chemistry_gocart2G_init
!==================================================================================================================
