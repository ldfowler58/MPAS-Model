! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_finalize
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config
 use mpas_chemistry_gocart2G_vars,only: CA2G_bc,CA2G_br,CA2G_oc,DU2G,NI2G,SS2G,SU2G


 implicit none
 private
 public:: gocart2G_finalize


 contains


!==================================================================================================================
 subroutine gocart2G_finalize(configs)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs

!--- local pointers:
 logical,pointer:: do_CA2Gbc,do_CA2Gbr,do_CA2Goc,do_DU2G,do_NI2G,do_SS2G,do_SU2G

!----------------------------------------------------------------------------------------------------------------- 

 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbc',do_CA2Gbc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbr',do_CA2Gbr)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Goc',do_CA2Goc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_DU2G' ,do_DU2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_do_NI2G' ,do_NI2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SS2G' ,do_SS2G   )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SU2G' ,do_SU2G   )

 if(do_CA2Gbc) call CA2G_bc%gocart2G_deallocate()
 if(do_CA2Gbr) call CA2G_br%gocart2G_deallocate()
 if(do_CA2Goc) call CA2G_oc%gocart2G_deallocate()
 if(do_DU2G) call DU2G%gocart2G_deallocate()
 if(do_NI2G) call NI2G%gocart2G_deallocate()
 if(do_SS2G) call SS2G%gocart2G_deallocate()
 if(do_SU2G) call SU2G%gocart2G_deallocate()

 end subroutine gocart2G_finalize

!==================================================================================================================
 end module mpas_chemistry_gocart2G_finalize
!==================================================================================================================
