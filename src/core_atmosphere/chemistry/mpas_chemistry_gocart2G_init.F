!=================================================================================================================
 module mpas_chemistry_gocart2G_init
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines

 use GOCART2G_MieMod_smiol
 use mpas_chemistry_gocart2G_vars


 implicit none
 private
 public:: init_chemistry_gocart2G


 contains


!=================================================================================================================
 subroutine init_chemistry_gocart2G(dminfo,mesh)
!=================================================================================================================

!--- input arguments:
 type(dm_info),intent(in):: dminfo
 type(mpas_pool_type),intent(in):: mesh

!--- local variables:
 character(len=StrKIND):: fnm
 integer:: its,ite,jts,jte,kts,kte

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_chemistry_gocart2G:')


!--- initializes dimensions used to run the GOCART-2G chemistry:
 call mpas_gocart2G%gocart2G_dims(mesh)

 its = mpas_gocart2G%its
 ite = mpas_gocart2G%ite
 jts = mpas_gocart2G%jts
 jte = mpas_gocart2G%jte
 kts = mpas_gocart2G%kts
 kte = mpas_gocart2G%kte


!--- CA2G (BC,BR,OC):
   !initializes and allocates all parameters and arrays related to CA2G:
   !call CA2G_bc_params%load_GridComp(kts,kte)
   !call CA2G_br_params%load_GridComp(kts,kte)
   !call CA2G_oc_params%load_GridComp(kts,kte)
   !call CA2G_bc%gocart2G_allocate(its,ite,jts,jte,kts,kte)
   !call CA2G_br%gocart2G_allocate(its,ite,jts,jte,kts,kte)
   !call CA2G_oc%gocart2G_allocate(its,ite,jts,jte,kts,kte)

   !creates radiation Mie table for CA2G:
   !CA2G_bc_params%rad_Mie = GOCART2G_Mie(dminfo,trim(file))
   !CA2G_br_params%rad_Mie = GOCART2G_Mie(dminfo,trim(file))
   !CA2G_oc_params%rad_Mie = GOCART2G_Mie(dminfo,trim(file))
   
   !create diagnostics Mie table for SU2G:
   !CA2G_bc_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)
   !CA2G_br_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)
   !CA2G_oc_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)


!--- DU2G:
   !initializes and allocates all parameters and arrays related to DU2G:
   !call DU2G_params%load_GridComp(kts,kte)
   !call DU2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

   !create radiation Mie table for DU2G:
   !DU2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(file))

   !create diagnostics Mie table for DU2G:
   !DU2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)


!--- NI2G:
   !initializes and allocates all parameters and arrays related to NI2G:
   !call NI2G_params%load_GridComp(kts,kte)
   !call NI2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

   !create radiation Mie table for SU2G:
   !NI2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(file))

   !create diagnostics Mie table for SU2G:
   !NI2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)


!--- SS2G:
   !initializes and allocates all parameters and arrays related to SS2G:
   !call SS2G_params%load_GridComp(kts,kte)
   !call SI2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

   !create radiation Mie table for SU2G:
   !NI2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(file))

   !create diagnostics Mie table for SU2G:
   !NI2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)


!--- SU2G:
   !initializes and allocates all parameters and arrays related to SU2G:
   call SU2G_params%load_GridComp(kts,kte)
   call SU2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)

   !create radiation Mie table for SU2G:
   fnm = 'opticsBands_SU.v1_3.RRTMG.nc'
   SU2G_params%rad_Mie = GOCART2G_Mie(dminfo,trim(fnm))

   !create diagnostics Mie table for SU2G:
   !SU2G_params%diag_Mie = GOCART2G_Mie(dminfo,trim(file_),channels_*1.e-9,nmom=nmom)


 call mpas_log_write('--- end subroutine init_chemistry_gocart2G:')
 call mpas_log_write(' ')

 end subroutine init_chemistry_gocart2G

!=================================================================================================================
 end module mpas_chemistry_gocart2G_init
!=================================================================================================================
