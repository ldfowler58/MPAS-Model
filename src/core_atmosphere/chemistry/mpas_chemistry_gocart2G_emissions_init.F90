
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_emissions_init
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_chemistry_gocart2G_emissions
 use mpas_chemistry_gocart2G_interface
 use mpas_chemistry_gocart2G_vars,only: mpas_emis_gocart2G,mpas_gocart2G


 contains


!==================================================================================================================
 subroutine init_gocart2G_emissions(mesh,anth_emissions,biob_emissions,BIOG_emissions)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: anth_emissions
 type(mpas_pool_type),intent(in):: biob_emissions
 type(mpas_pool_type),intent(in):: BIOG_emissions

!--- local variables:
 integer:: its,ite,jts,jte,kts,kte

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_gocart2G_emissions:')


!--- initializes dimensions needed to initialize emissions:
 call mpas_gocart2G%gocart2G_dims(mesh)

 its   = mpas_gocart2G%its
 ite   = mpas_gocart2G%ite
 jts   = mpas_gocart2G%jts
 jte   = mpas_gocart2G%jte
 kts   = mpas_gocart2G%kts
 kte   = mpas_gocart2G%kte


!--- allocates local arrays needed to run GOCART2G emissions:
 call mpas_emis_gocart2G%gocart2G_allocate(its,ite,jts,jte,kts,kte)


!--- initializes local arrays needed to run GOCART2G emissions:
 call mpas_emis_gocart2G%gocart2G_emissions(anth_emissions,biob_emissions,BIOG_emissions,its,ite,jts,jte,kts,kte)


 call mpas_log_write('--- end subroutine init_gocart2G_emissions.')

 end subroutine init_gocart2G_emissions

!==================================================================================================================
 end module mpas_chemistry_gocart2G_emissions_init
!==================================================================================================================
