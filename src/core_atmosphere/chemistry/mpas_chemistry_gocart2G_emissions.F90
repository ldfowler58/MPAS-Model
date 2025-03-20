! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_emissions
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array


 implicit none
 private


!--- parameters:
 real(kind=RKIND),parameter:: voc_AnthroFactor      = 0.069 ! (g/g CO)
 real(kind=RKIND),parameter:: voc_BiomassBurnFactor = 0.013 ! (g/g CO)
 real(kind=RKIND),parameter:: voc_BiogIsopFactor    = 0.015 ! (-)
 real(kind=RKIND),parameter:: voc_BiogMonxFactor    = 0.005 ! (-)


 type,public:: emis_gocart2G
    integer:: its,ite,jts,jte,kts,kte

    !--- anthropogenic emissions:
    real(kind=RKIND),dimension(:,:),pointer:: bc_antebc1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: bc_antebc2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: bc_ship         => null()
    real(kind=RKIND),dimension(:,:),pointer:: bc_aviation_lto => null()
    real(kind=RKIND),dimension(:,:),pointer:: bc_aviation_cds => null()
    real(kind=RKIND),dimension(:,:),pointer:: bc_aviation_crs => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: bc_aircraft   => null()

    real(kind=RKIND),dimension(:,:),pointer:: br_antebr1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_antebr2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_ship         => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_aviation_lto => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_aviation_cds => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_aviation_crs => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: br_aircraft   => null()

    real(kind=RKIND),dimension(:,:),pointer:: oc_anteoc1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_anteoc2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_ship         => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_aviation_lto => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_aviation_cds => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_aviation_crs => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: oc_aircraft   => null()

    real(kind=RKIND),dimension(:,:),pointer:: nh3_ag          => null()
    real(kind=RKIND),dimension(:,:),pointer:: nh3_en          => null()
    real(kind=RKIND),dimension(:,:),pointer:: nh3_in          => null()
    real(kind=RKIND),dimension(:,:),pointer:: nh3_oc          => null()
    real(kind=RKIND),dimension(:,:),pointer:: nh3_re          => null()
    real(kind=RKIND),dimension(:,:),pointer:: nh3_tr          => null()

    real(kind=RKIND),dimension(:,:),pointer:: soap_anthrop    => null()

    real(kind=RKIND),dimension(:,:),pointer:: su_anthrol1     => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_anthrol2     => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_shipso2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_shipso4      => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_dmso         => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_aviation_lto => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_aviation_cds => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_aviation_crs => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: su_aircraft   => null()

    !--- biomass burning emissions:
    real(kind=RKIND),dimension(:,:),pointer:: bc_biomass      => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_biomass      => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_biomass      => null()
    real(kind=RKIND),dimension(:,:),pointer:: su_biomass      => null()
    real(kind=RKIND),dimension(:,:),pointer:: nh3_bb          => null()

    real(kind=RKIND),dimension(:,:),pointer:: soap_biomass    => null()
    real(kind=RKIND),dimension(:,:),pointer:: soas_biomass    => null()

    !--- biofuel emissions:
    real(kind=RKIND),dimension(:,:),pointer:: bc_biofuel      => null()
    real(kind=RKIND),dimension(:,:),pointer:: br_biofuel      => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_biofuel      => null()

    real(kind=RKIND),dimension(:,:),pointer:: soap_biofuel    => null()
    real(kind=RKIND),dimension(:,:),pointer:: soas_biofuel    => null()

    !--- biogenic emissions:
    real(kind=RKIND),dimension(:,:),pointer:: br_terpene      => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_isoprene     => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_mtpa         => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_mtpo         => null()
    real(kind=RKIND),dimension(:,:),pointer:: oc_limo         => null()

    real(kind=RKIND),dimension(:,:),pointer:: soap_biogenic   => null()
    real(kind=RKIND),dimension(:,:),pointer:: soas_biogenic   => null()

    contains
       procedure:: gocart2G_allocate   => mpas_chemistry_gocart2G_emissions_allocate
       procedure:: gocart2G_deallocate => mpas_chemistry_gocart2G_emissions_deallocate
       procedure:: gocart2G_emissions  => mpas_chemistry_gocart2G_emissions_init
 end type


 contains


!==================================================================================================================
 subroutine mpas_chemistry_gocart2G_emissions_allocate(self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!inout arguments:
 class(emis_gocart2G),intent(inout):: self

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_emissions_allocate:')


!--- anthropogenic emissions:
 if(.not.associated(self%bc_antebc1)     ) allocate(self%bc_antebc1(its:ite,jts:jte)         )
 if(.not.associated(self%bc_antebc2)     ) allocate(self%bc_antebc2(its:ite,jts:jte)         )
 if(.not.associated(self%bc_ship)        ) allocate(self%bc_ship(its:ite,jts:jte)            )
 if(.not.associated(self%bc_aviation_lto)) allocate(self%bc_aviation_lto(its:ite,jts:jte)    )
 if(.not.associated(self%bc_aviation_cds)) allocate(self%bc_aviation_cds(its:ite,jts:jte)    )
 if(.not.associated(self%bc_aviation_crs)) allocate(self%bc_aviation_crs(its:ite,jts:jte)    )
 if(.not.associated(self%bc_aircraft)    ) allocate(self%bc_aircraft(its:ite,jts:jte,kts:kte))

 if(.not.associated(self%br_antebr1)     ) allocate(self%br_antebr1(its:ite,jts:jte)         )
 if(.not.associated(self%br_antebr2)     ) allocate(self%br_antebr2(its:ite,jts:jte)         )
 if(.not.associated(self%br_ship)        ) allocate(self%br_ship(its:ite,jts:jte)            )
 if(.not.associated(self%br_aviation_lto)) allocate(self%br_aviation_lto(its:ite,jts:jte)    )
 if(.not.associated(self%br_aviation_cds)) allocate(self%br_aviation_cds(its:ite,jts:jte)    )
 if(.not.associated(self%br_aviation_crs)) allocate(self%br_aviation_crs(its:ite,jts:jte)    )
 if(.not.associated(self%br_aircraft)    ) allocate(self%br_aircraft(its:ite,jts:jte,kts:kte))

 if(.not.associated(self%oc_anteoc1)     ) allocate(self%oc_anteoc1(its:ite,jts:jte)         )
 if(.not.associated(self%oc_anteoc2)     ) allocate(self%oc_anteoc2(its:ite,jts:jte)         )
 if(.not.associated(self%oc_ship)        ) allocate(self%oc_ship(its:ite,jts:jte)            )
 if(.not.associated(self%oc_aviation_lto)) allocate(self%oc_aviation_lto(its:ite,jts:jte)    )
 if(.not.associated(self%oc_aviation_cds)) allocate(self%oc_aviation_cds(its:ite,jts:jte)    )
 if(.not.associated(self%oc_aviation_crs)) allocate(self%oc_aviation_crs(its:ite,jts:jte)    )
 if(.not.associated(self%oc_aircraft)    ) allocate(self%oc_aircraft(its:ite,jts:jte,kts:kte))

 if(.not.associated(self%nh3_ag)         ) allocate(self%nh3_ag(its:ite,jts:jte)             )
 if(.not.associated(self%nh3_en)         ) allocate(self%nh3_en(its:ite,jts:jte)             )
 if(.not.associated(self%nh3_in)         ) allocate(self%nh3_in(its:ite,jts:jte)             )
 if(.not.associated(self%nh3_oc)         ) allocate(self%nh3_oc(its:ite,jts:jte)             )
 if(.not.associated(self%nh3_re)         ) allocate(self%nh3_re(its:ite,jts:jte)             )
 if(.not.associated(self%nh3_tr)         ) allocate(self%nh3_tr(its:ite,jts:jte)             )

 if(.not.associated(self%soap_anthrop)   ) allocate(self%soap_anthrop(its:ite,jts:jte)       )

 if(.not.associated(self%su_anthrol1)    ) allocate(self%su_anthrol1(its:ite,jts:jte)        )
 if(.not.associated(self%su_anthrol2)    ) allocate(self%su_anthrol2(its:ite,jts:jte)        )
 if(.not.associated(self%su_shipso2)     ) allocate(self%su_shipso2(its:ite,jts:jte)         )
 if(.not.associated(self%su_shipso4)     ) allocate(self%su_shipso4(its:ite,jts:jte)         )
 if(.not.associated(self%su_dmso)        ) allocate(self%su_dmso(its:ite,jts:jte)            )
 if(.not.associated(self%su_aviation_lto)) allocate(self%su_aviation_lto(its:ite,jts:jte)    )
 if(.not.associated(self%su_aviation_cds)) allocate(self%su_aviation_cds(its:ite,jts:jte)    )
 if(.not.associated(self%su_aviation_crs)) allocate(self%su_aviation_crs(its:ite,jts:jte)    )
 if(.not.associated(self%su_aircraft)    ) allocate(self%su_aircraft(its:ite,jts:jte,kts:kte))


!--- biomass burning emissions:
 if(.not.associated(self%bc_biomass)     ) allocate(self%bc_biomass(its:ite,jts:jte)         )
 if(.not.associated(self%br_biomass)     ) allocate(self%br_biomass(its:ite,jts:jte)         )
 if(.not.associated(self%oc_biomass)     ) allocate(self%oc_biomass(its:ite,jts:jte)         )
 if(.not.associated(self%su_biomass)     ) allocate(self%su_biomass(its:ite,jts:jte)         )
 if(.not.associated(self%nh3_bb)         ) allocate(self%nh3_bb(its:ite,jts:jte)             )

 if(.not.associated(self%soap_biomass)   ) allocate(self%soap_biomass(its:ite,jts:jte)       )
 if(.not.associated(self%soas_biomass)   ) allocate(self%soas_biomass(its:ite,jts:jte)       )

!--- biofuel emissions:
 if(.not.associated(self%bc_biofuel)     ) allocate(self%bc_biofuel(its:ite,jts:jte)         )
 if(.not.associated(self%br_biofuel)     ) allocate(self%br_biofuel(its:ite,jts:jte)         )
 if(.not.associated(self%oc_biofuel)     ) allocate(self%oc_biofuel(its:ite,jts:jte)         )

 if(.not.associated(self%soap_biofuel)   ) allocate(self%soap_biofuel(its:ite,jts:jte)       )
 if(.not.associated(self%soas_biofuel)   ) allocate(self%soas_biofuel(its:ite,jts:jte)       )


!--- biogenic emissions:
 if(.not.associated(self%br_terpene)     ) allocate(self%br_terpene(its:ite,jts:jte)         )
 if(.not.associated(self%oc_isoprene)    ) allocate(self%oc_isoprene(its:ite,jts:jte)        )
 if(.not.associated(self%oc_mtpa)        ) allocate(self%oc_mtpa(its:ite,jts:jte)            )
 if(.not.associated(self%oc_mtpo)        ) allocate(self%oc_mtpo(its:ite,jts:jte)            )
 if(.not.associated(self%oc_limo)        ) allocate(self%oc_limo(its:ite,jts:jte)            )

 if(.not.associated(self%soap_biogenic)  ) allocate(self%soap_biogenic(its:ite,jts:jte)      )
 if(.not.associated(self%soas_biogenic)  ) allocate(self%soas_biogenic(its:ite,jts:jte)      )


 call mpas_log_write('--- end subroutine mpas_chemistry_gocart2G_emissions_allocate.')

 end subroutine mpas_chemistry_gocart2G_emissions_allocate

!==================================================================================================================
 subroutine mpas_chemistry_gocart2G_emissions_deallocate(self)
!==================================================================================================================

!inout arguments:
 class(emis_gocart2G),intent(inout):: self

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_emissions_deallocate:')


!--- anthropogenic emissions:
 if(associated(self%bc_antebc1)     ) deallocate(self%bc_antebc1     )
 if(associated(self%bc_antebc2)     ) deallocate(self%bc_antebc2     )
 if(associated(self%bc_ship)        ) deallocate(self%bc_ship        )
 if(associated(self%bc_aviation_lto)) deallocate(self%bc_aviation_lto)
 if(associated(self%bc_aviation_cds)) deallocate(self%bc_aviation_cds)
 if(associated(self%bc_aviation_crs)) deallocate(self%bc_aviation_crs)
 if(associated(self%bc_aircraft)    ) deallocate(self%bc_aircraft    )

 if(associated(self%br_antebr1)     ) deallocate(self%br_antebr1     )
 if(associated(self%br_antebr2)     ) deallocate(self%br_antebr2     )
 if(associated(self%br_ship)        ) deallocate(self%br_ship        )
 if(associated(self%br_aviation_lto)) deallocate(self%br_aviation_lto)
 if(associated(self%br_aviation_cds)) deallocate(self%br_aviation_cds)
 if(associated(self%br_aviation_crs)) deallocate(self%br_aviation_crs)
 if(associated(self%br_aircraft)    ) deallocate(self%br_aircraft    )

 if(associated(self%oc_anteoc1)     ) deallocate(self%oc_anteoc1     )
 if(associated(self%oc_anteoc2)     ) deallocate(self%oc_anteoc2     )
 if(associated(self%oc_ship)        ) deallocate(self%oc_ship        )
 if(associated(self%oc_aviation_lto)) deallocate(self%oc_aviation_lto)
 if(associated(self%oc_aviation_cds)) deallocate(self%oc_aviation_cds)
 if(associated(self%oc_aviation_crs)) deallocate(self%oc_aviation_crs)
 if(associated(self%oc_aircraft)    ) deallocate(self%oc_aircraft    )

 if(associated(self%nh3_ag)         ) deallocate(self%nh3_ag         )
 if(associated(self%nh3_en)         ) deallocate(self%nh3_en         )
 if(associated(self%nh3_in)         ) deallocate(self%nh3_in         )
 if(associated(self%nh3_oc)         ) deallocate(self%nh3_oc         )
 if(associated(self%nh3_re)         ) deallocate(self%nh3_re         )
 if(associated(self%nh3_tr)         ) deallocate(self%nh3_tr         )

 if(associated(self%soap_anthrop)   ) deallocate(self%soap_anthrop   )

 if(associated(self%su_anthrol1)    ) deallocate(self%su_anthrol1    )
 if(associated(self%su_anthrol2)    ) deallocate(self%su_anthrol2    )
 if(associated(self%su_shipso2)     ) deallocate(self%su_shipso2     )
 if(associated(self%su_shipso4)     ) deallocate(self%su_shipso4     )
 if(associated(self%su_dmso)        ) deallocate(self%su_dmso        )
 if(associated(self%su_aviation_lto)) deallocate(self%su_aviation_lto)
 if(associated(self%su_aviation_cds)) deallocate(self%su_aviation_cds)
 if(associated(self%su_aviation_crs)) deallocate(self%su_aviation_crs)
 if(associated(self%su_aircraft)    ) deallocate(self%su_aircraft    )


!--- biomass burning emissions:
 if(associated(self%bc_biomass)     ) deallocate(self%bc_biomass     )
 if(associated(self%br_biomass)     ) deallocate(self%br_biomass     )
 if(associated(self%oc_biomass)     ) deallocate(self%oc_biomass     )
 if(associated(self%su_biomass)     ) deallocate(self%su_biomass     )
 if(associated(self%nh3_bb)         ) deallocate(self%nh3_bb         )

 if(associated(self%soap_biomass)   ) deallocate(self%soap_biomass   )
 if(associated(self%soas_biomass)   ) deallocate(self%soas_biomass   )

!--- biofuel emissions:
 if(associated(self%bc_biofuel)     ) deallocate(self%bc_biofuel     )
 if(associated(self%br_biofuel)     ) deallocate(self%br_biofuel     )
 if(associated(self%oc_biofuel)     ) deallocate(self%oc_biofuel     )

 if(associated(self%soap_biofuel)   ) deallocate(self%soap_biofuel   )
 if(associated(self%soas_biofuel)   ) deallocate(self%soas_biofuel   )


!--- biogenic emissions:
 if(associated(self%br_terpene)     ) deallocate(self%br_terpene     )
 if(associated(self%oc_isoprene)    ) deallocate(self%oc_isoprene    )
 if(associated(self%oc_mtpa)        ) deallocate(self%oc_mtpa        )
 if(associated(self%oc_mtpo)        ) deallocate(self%oc_mtpo        )
 if(associated(self%oc_limo)        ) deallocate(self%oc_limo        )

 if(associated(self%soap_biogenic)  ) deallocate(self%soap_biogenic  )
 if(associated(self%soas_biogenic)  ) deallocate(self%soas_biogenic  )


 call mpas_log_write('--- end subroutine mpas_chemistry_gocart2G_emissions_deallocate.')

 end subroutine mpas_chemistry_gocart2G_emissions_deallocate

!==================================================================================================================
 subroutine mpas_chemistry_gocart2G_emissions_init(self,anth_emissions,biob_emissions,biog_emissions, &
                                                   its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: anth_emissions
 type(mpas_pool_type),intent(in):: biob_emissions
 type(mpas_pool_type),intent(in):: biog_emissions
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(emis_gocart2G),intent(inout):: self

!--- local variables and pointers:
 real(kind=RKIND),dimension(:),pointer:: bc_anth_less100m,bc_anth_less500m,bc_anth_ship,bc_anth_aviation_lto, &
                                         bc_anth_aviation_cds,bc_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: bc_anth_aircraft
 real(kind=RKIND),dimension(:),pointer:: br_anth_less100m,br_anth_less500m,br_anth_ship,br_anth_aviation_lto, &
                                         br_anth_aviation_cds,br_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: br_anth_aircraft
 real(kind=RKIND),dimension(:),pointer:: oc_anth_less100m,oc_anth_less500m,oc_anth_ship,oc_anth_aviation_lto, &
                                         oc_anth_aviation_cds,oc_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: oc_anth_aircraft
 real(kind=RKIND),dimension(:),pointer:: nh3_anth_ag,nh3_anth_bb,nh3_anth_en,nh3_anth_in, &
                                         nh3_anth_oc,nh3_anth_re,nh3_anth_tr
 real(kind=RKIND),dimension(:),pointer:: su_anth_less100m,su_anth_less500m,su_anth_shipso2,su_anth_shipso4, &
                                         su_dmso_em,su_anth_aviation_lto,su_anth_aviation_cds,su_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: su_anth_aircraft
 real(kind=RKIND),dimension(:),pointer:: co_anth_em

 real(kind=RKIND),dimension(:),pointer:: bc_biob_em,br_biob_em,oc_biob_em,ni_biob_em,su_biob_em, &
                                         co_biob_em

 real(kind=RKIND),dimension(:),pointer:: iso_biog_em,mnt_biog_em,mnta_biog_em,mntb_biog_em

!--- local variables and arrays:
 integer:: i,j,k

 real(kind=RKIND),parameter:: Avogadro = 6.02214076e23 ! (per mole).
 real(kind=RKIND),parameter:: fMassBC  = 12.011        ! (grams per mole).
 real(kind=RKIND),parameter:: fMassBR  = 12.011        ! (grams per mole).
 real(kind=RKIND),parameter:: fMassOC  = 12.011        ! (grams per mole).
 real(kind=RKIND),parameter:: fMassNH3 = 17.031        ! (grams per mole).
 real(kind=RKIND),parameter:: fMassSO2 = 64.066        ! (grams per mole).

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_emissions_init:')


!--- anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'bc_anth_less100m'    ,bc_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'bc_anth_less500m'    ,bc_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'bc_anth_ship'        ,bc_anth_ship        )
 call mpas_pool_get_array(anth_emissions,'bc_anth_aviation_lto',bc_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'bc_anth_aviation_cds',bc_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'bc_anth_aviation_crs',bc_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'bc_anth_aircraft'    ,bc_anth_aircraft    )

 call mpas_pool_get_array(anth_emissions,'br_anth_less100m'    ,br_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'br_anth_less500m'    ,br_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'br_anth_ship'        ,br_anth_ship        )
 call mpas_pool_get_array(anth_emissions,'br_anth_aviation_lto',br_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'br_anth_aviation_cds',br_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'br_anth_aviation_crs',br_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'br_anth_aircraft'    ,br_anth_aircraft    )

 call mpas_pool_get_array(anth_emissions,'oc_anth_less100m'    ,oc_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'oc_anth_less500m'    ,oc_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'oc_anth_ship'        ,oc_anth_ship        )
 call mpas_pool_get_array(anth_emissions,'oc_anth_aviation_lto',oc_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'oc_anth_aviation_cds',oc_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'oc_anth_aviation_crs',oc_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'oc_anth_aircraft'    ,oc_anth_aircraft    )

 call mpas_pool_get_array(anth_emissions,'nh3_anth_ag'         ,nh3_anth_ag         )
 call mpas_pool_get_array(anth_emissions,'nh3_anth_bb'         ,nh3_anth_bb         )
 call mpas_pool_get_array(anth_emissions,'nh3_anth_en'         ,nh3_anth_en         )
 call mpas_pool_get_array(anth_emissions,'nh3_anth_in'         ,nh3_anth_in         )
 call mpas_pool_get_array(anth_emissions,'nh3_anth_oc'         ,nh3_anth_oc         )
 call mpas_pool_get_array(anth_emissions,'nh3_anth_re'         ,nh3_anth_re         )
 call mpas_pool_get_array(anth_emissions,'nh3_anth_tr'         ,nh3_anth_tr         )

 call mpas_pool_get_array(anth_emissions,'su_anth_less100m'    ,su_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'su_anth_less500m'    ,su_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'su_anth_shipso2'     ,su_anth_shipso2     )
 call mpas_pool_get_array(anth_emissions,'su_anth_shipso4'     ,su_anth_shipso4     )
 call mpas_pool_get_array(anth_emissions,'su_anth_aviation_lto',su_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'su_anth_aviation_cds',su_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'su_anth_aviation_crs',su_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'su_anth_aircraft'    ,su_anth_aircraft    )
 call mpas_pool_get_array(anth_emissions,'su_dmso_em'          ,su_dmso_em          )

 call mpas_pool_get_array(anth_emissions,'co_anth_em'          ,co_anth_em          )

 do j = jts,jte
    do i = its,ite
       !--- black carbon:
       self%bc_antebc1(i,j)      = bc_anth_less100m(i)
       self%bc_antebc2(i,j)      = bc_anth_less500m(i)
       self%bc_ship(i,j)         = bc_anth_ship(i)
       self%bc_aviation_lto(i,j) = bc_anth_aviation_lto(i)
       self%bc_aviation_cds(i,j) = bc_anth_aviation_cds(i)
       self%bc_aviation_crs(i,j) = bc_anth_aviation_crs(i)
       do k = kts,kte
          self%bc_aircraft(i,j,k) = bc_anth_aircraft(k,i)
       enddo

       !--- brown carbon:
       self%br_antebr1(i,j)      = br_anth_less100m(i)
       self%br_antebr2(i,j)      = br_anth_less500m(i)
       self%br_ship(i,j)         = br_anth_ship(i)
       self%br_aviation_lto(i,j) = br_anth_aviation_lto(i)
       self%br_aviation_cds(i,j) = br_anth_aviation_cds(i)
       self%br_aviation_crs(i,j) = br_anth_aviation_crs(i)
       do k = kts,kte
          self%br_aircraft(i,j,k) = br_anth_aircraft(k,i)
       enddo

       !--- organic carbon:
       self%oc_anteoc1(i,j)      = oc_anth_less100m(i)
       self%oc_anteoc2(i,j)      = oc_anth_less500m(i)
       self%oc_ship(i,j)         = oc_anth_ship(i)
       self%oc_aviation_lto(i,j) = oc_anth_aviation_lto(i)
       self%oc_aviation_cds(i,j) = oc_anth_aviation_cds(i)
       self%oc_aviation_crs(i,j) = oc_anth_aviation_crs(i)
       do k = kts,kte
          self%oc_aircraft(i,j,k) = oc_anth_aircraft(k,i)
       enddo

       !--- nitrate:
       self%nh3_ag(i,j)  = nh3_anth_ag(i)
       self%nh3_en(i,j)  = nh3_anth_en(i)
       self%nh3_in(i,j)  = nh3_anth_in(i)
       self%nh3_oc(i,j)  = nh3_anth_oc(i)
       self%nh3_re(i,j)  = nh3_anth_re(i)
       self%nh3_tr(i,j)  = nh3_anth_tr(i)

       !--- sulfate:
       self%su_anthrol1(i,j)     = su_anth_less100m(i)
       self%su_anthrol2(i,j)     = su_anth_less500m(i)
       self%su_shipso2(i,j)      = su_anth_shipso2(i)
       self%su_shipso4(i,j)      = su_anth_shipso4(i)
       self%su_aviation_lto(i,j) = su_anth_aviation_lto(i)
       self%su_aviation_cds(i,j) = su_anth_aviation_cds(i)
       self%su_aviation_crs(i,j) = su_anth_aviation_crs(i)
       self%su_dmso(i,j)         = su_dmso_em(i)
       do k = kts,kte
          self%su_aircraft(i,j,k) = su_anth_aircraft(k,i)
       enddo

       !--- secondary organic aerosols:
       self%soap_anthrop(i,j) = voc_AnthroFactor*co_anth_em(i)
    enddo
 enddo


!--- biomass burning emissions:
 call mpas_pool_get_array(biob_emissions,'bc_biob_em',bc_biob_em)
 call mpas_pool_get_array(biob_emissions,'br_biob_em',br_biob_em)
 call mpas_pool_get_array(biob_emissions,'oc_biob_em',oc_biob_em)
 call mpas_pool_get_array(biob_emissions,'ni_biob_em',ni_biob_em)
 call mpas_pool_get_array(biob_emissions,'su_biob_em',su_biob_em)
 call mpas_pool_get_array(biob_emissions,'co_biob_em',co_biob_em)

 do j = jts,jte
    do i = its,ite
       self%bc_biomass(i,j) = bc_biob_em(i)
       self%br_biomass(i,j) = br_biob_em(i)
       self%oc_biomass(i,j) = oc_biob_em(i)
       self%su_biomass(i,j) = su_biob_em(i)
       self%nh3_bb(i,j)     = ni_biob_em(i)

       !--- conversion from molecules/cm^2/s to kg/m^2/s:
       self%bc_biomass(i,j) = 10.*(fMassBC/Avogadro)*self%bc_biomass(i,j)
       self%br_biomass(i,j) = 10.*(fMassBR/Avogadro)*self%br_biomass(i,j)
       self%oc_biomass(i,j) = 10.*(fMassOC/Avogadro)*self%oc_biomass(i,j)
       self%su_biomass(i,j) = 10.*(fMassSO2/Avogadro)*self%su_biomass(i,j)
       self%nh3_bb(i,j)     = 10.*(fMassNH3/Avogadro)*self%nh3_bb(i,j)

       !--- secondary organic aerosols:
       self%soap_biomass(i,j) = voc_BiomassBurnFactor*co_biob_em(i)
       self%soap_biomass(i,j) = 10.*(fMassOC/Avogadro)*self%soap_biomass(i,j)
       self%soas_biomass(i,j) = 0._RKIND
    enddo
 enddo


!--- biofuel emissions:
 do j = jts,jte
    do i = its,ite
       self%bc_biofuel(i,j) = 0._RKIND
       self%br_biofuel(i,j) = 0._RKIND
       self%oc_biofuel(i,j) = 0._RKIND

       !--- secondary organic aerosols:
       self%soap_biofuel(i,j) = 0._RKIND
       self%soas_biofuel(i,j) = 0._RKIND
    enddo
 enddo


!--- biogenic emissions:
 call mpas_pool_get_array(biog_emissions,'iso_biog_em' ,iso_biog_em )
 call mpas_pool_get_array(biog_emissions,'mnt_biog_em' ,mnt_biog_em )
 call mpas_pool_get_array(biog_emissions,'mnta_biog_em',mnta_biog_em)
 call mpas_pool_get_array(biog_emissions,'mntb_biog_em',mntb_biog_em)

 do j = jts,jte
    do i = its,ite
       self%br_terpene(i,j)  = 0._RKIND
       self%oc_isoprene(i,j) = 0._RKIND
       self%oc_mtpa(i,j)     = 0._RKIND
       self%oc_mtpo(i,j)     = 0._RKIND
       self%oc_limo(i,j)     = 0._RKIND

       !--- secondary organic aerosols:
       self%soap_biogenic(i,j) = voc_BiogIsopFactor*iso_biog_em(i)
                               + voc_BiogMonxFactor*(mnt_biog_em(i)+mnta_biog_em(i)+mntb_biog_em(i))

       self%soas_biogenic(i,j) = voc_BiogIsopFactor*iso_biog_em(i)
                               + voc_BiogMonxFactor*(mnt_biog_em(i)+mnta_biog_em(i)+mntb_biog_em(i))
    enddo
 enddo


 call mpas_log_write('--- end subroutine mpas_chemistry_gocart2G_emissions_init.')

 end subroutine mpas_chemistry_gocart2G_emissions_init

!==================================================================================================================
 end module mpas_chemistry_gocart2G_emissions
!==================================================================================================================

