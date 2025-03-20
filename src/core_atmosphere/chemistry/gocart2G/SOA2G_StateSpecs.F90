! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module SOA2G_StateSpecs
 use mpas_kind_types,only: RKIND

 implicit none
 public

!this module is the state variable specification file for secondary organic aerosols.

!schema_version: 2.0.0
!component: SOA


 type SOA2G_State
 
!category: IMPORT
 real(kind=RKIND),dimension(:,:),pointer  :: zpbl          => null() ! planetary_boundary_layer_height (m)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:,:),pointer:: airdens       => null() ! moist_air_density (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer:: delp          => null() ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: delz          => null() ! geometric_layer_thickness (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: zle           => null() ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),pointer:: ple           => null() ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: soap_oh       => null() !climatological oh source (mole mole-1)
!..................................................................................................................
 real(kind=RKIND),dimension(:,:),pointer  :: soap_anthro   => null() ! SOAP anthropogenic emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: soap_biomass  => null() ! SOAP biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: soap_biofuel  => null() ! SOAP biofuel emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),pointer  :: soap_biogenic => null() ! SOAP biogenic emissions (kg m-2 s-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer:: soapa_prod    => null() ! SOAP anthropogenic production (kg m-3 s-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: soapbb_prod   => null() ! SOAP biomass burning production (kg m-3 s-1)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),pointer:: soap_a        => null() !
 real(kind=RKIND),dimension(:,:,:),pointer:: soap_bb       => null() ! 
 

 contains
    procedure:: gocart2G_allocate   => SOA2G_StateSpecsInit
    procedure:: gocart2G_deallocate => SOA2G_StateSpecsFinalize

 end type SOA2G_State


 contains


!=================================================================================================================
 subroutine SOA2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(SOA2G_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
!if(.not.associated(self%zpbl)         ) allocate(self%zpbl(its:ite,jts:jte)               )
!..................................................................................................................
!if(.not.associated(self%airdens)      ) allocate(self%airdens(its:ite,jts:jte,kts:kte)    )
!if(.not.associated(self%delp)         ) allocate(self%delp(its:ite,jts:jte,kts:kte)       )
!if(.not.associated(self%delz)         ) allocate(self%delz(its:ite,jts:jte,kts:kte)       )
!if(.not.associated(self%zle)          ) allocate(self%zle(its:ite,jts:jte,kts:kte+1)      )
!if(.not.associated(self%ple)          ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)      )
!if(.not.associated(self%soap_oh)      ) allocate(self%soap_oh(its:ite,jts:jte,kts:kte+1)  )
!..................................................................................................................
!if(.not.associated(self%soap_anthro)  ) allocate(self%soap_anthro(its:ite,jts:jte)        )
!if(.not.associated(self%soap_biomass) ) allocate(self%soap_biomass(its:ite,jts:jte)       )
!if(.not.associated(self%soap_biofuel) ) allocate(self%soap_biofuel(its:ite,jts:jte)       )
!if(.not.associated(self%soap_biogenic)) allocate(self%soap_biogenic(its:ite,jts:jte)      )

!category: EXPORT
 if(.not.associated(self%soapa_prod)   ) allocate(self%soapa_prod(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%soapbb_prod)  ) allocate(self%soapbb_prod(its:ite,jts:jte,kts:kte))

!category: INTERNAL
!if(.not.associated(self%soap_a)       ) allocate(self%soap_a(its:ite,jts:jte,kts_kte)     )
!if(.not.associated(self%soap_bb)      ) allocate(self%soap_bb(its:ite,jts:jte,kts_kte)    )

 end subroutine SOA2G_StateSpecsInit

!=================================================================================================================
 subroutine SOA2G_StateSpecsFinalize(self)
!=================================================================================================================

!--- inout arguments:
 class(SOA2G_State),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------

!category: IMPORT
!if(associated(self%zpbl)         ) deallocate(self%zpbl         )
!..................................................................................................................
!if(associated(self%airdens)      ) deallocate(self%airdens      )
!if(associated(self%delp)         ) deallocate(self%delp         )
!if(associated(self%delz)         ) deallocate(self%delz         )
!if(associated(self%zle)          ) deallocate(self%zle          )
!if(associated(self%ple)          ) deallocate(self%ple          )
!if(associated(self%soap_oh)      ) deallocate(self%soap_oh      )
!..................................................................................................................
!if(associated(self%soap_anthro)  ) deallocate(self%soap_anthro  )
!if(associated(self%soap_biomass) ) deallocate(self%soap_biomass )
!if(associated(self%soap_biofuel) ) deallocate(self%soap_biofuel )
!if(associated(self%soap_biogenic)) deallocate(self%soap_biogenic)

!category: EXPORT
 if(associated(self%soapa_prod)   ) deallocate(self%soapa_prod   )
 if(associated(self%soapbb_prod)  ) deallocate(self%soapbb_prod  )

!category: INTERNAL
!if(associated(self%soap_a)       ) deallocate(self%soap_a       )
!if(associated(self%soap_bb)      ) deallocate(self%soap_bb      )

 end subroutine SOA2G_StateSpecsFinalize

!=================================================================================================================
 end module SOA2G_StateSpecs
!=================================================================================================================
