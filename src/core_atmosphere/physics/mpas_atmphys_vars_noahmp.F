! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_vars_noahmp
 use mpas_kind_types

 use NoahmpIOVarType
 
 implicit none
 public
 save


!=================================================================================================================
!.. variables and arrays related to the Noahmp land-surface parameterization:
!=================================================================================================================

 type(NoahmpIO_type):: mpas_noahmp

!=================================================================================================================
 end module mpas_atmphys_vars_noahmp
!=================================================================================================================
