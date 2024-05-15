!=================================================================================================================
 module mpas_chemistry_gocart2G_interface
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension
 use mpas_atmphys_constants,only: R_d,R_v


 implicit none
 private


 type,public:: atm_gocart2G
    integer:: its,ite,jts,jte,kts,kte,ktep1


    !--- black,brown,and organic carbon mixing ratios (kg kg-1):
    real(kind=RKIND),dimension(:,:,:),pointer:: qbcphobic  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qbcphilic  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qbrphobic  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qbrphilic  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qocphobic  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qocphilic  => null()
    !--- dust mixing ratios (kg kg-1):
    real(kind=RKIND),dimension(:,:,:),pointer:: qdust1     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qdust2     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qdust3     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qdust4     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qdust5     => null()
    !--- nitrate mixing ratios (kg kg-1):
    real(kind=RKIND),dimension(:,:,:),pointer:: nitrate1   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nitrate2   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nitrate3   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nitrate4   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nitrate5   => null()
    !--- sea-salt mixing ratios (kg kg-1):
    real(kind=RKIND),dimension(:,:,:),pointer:: qseas1     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qseas2     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qseas3     => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qseas4     => null()
    !--- sulfate,sulfur dioxide,dimethyl sulfide,and methanesulfonic acid mixing ratios (kg kg-1):
    real(kind=RKIND),dimension(:,:,:),pointer:: qso4       => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qso2       => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qdms       => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: qmsa       => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: qsu2G    => null()


    !--- background DMS, OH, H2O2, and NO3:
    real(kind=RKIND),dimension(:,:),pointer  :: backg_dms  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: backg_oh   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: backg_h2o2 => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: backg_no3  => null()


    !--- CAMS surface emissions:
    real(kind=RKIND),dimension(:,:),pointer:: qbc1_em      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc1_em      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qnh3_em      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso2_em      => null()


    !--- miscellaneous variables:
    real(kind=RKIND):: dt


    !--- mesh fields:
    real(kind=RKIND),dimension(:,:),pointer:: xlat         => null()
    real(kind=RKIND),dimension(:,:),pointer:: xlon         => null()


    !--- surface fields:
    real(kind=RKIND),dimension(:,:),pointer:: frocean      => null()
    real(kind=RKIND),dimension(:,:),pointer:: fraci        => null()
    real(kind=RKIND),dimension(:,:),pointer:: frlake       => null()
    real(kind=RKIND),dimension(:,:),pointer:: lwi          => null()
    real(kind=RKIND),dimension(:,:),pointer:: tropp        => null()
    real(kind=RKIND),dimension(:,:),pointer:: u10m         => null()
    real(kind=RKIND),dimension(:,:),pointer:: v10m         => null()
    real(kind=RKIND),dimension(:,:),pointer:: area         => null()
    real(kind=RKIND),dimension(:,:),pointer:: zpbl         => null()
    real(kind=RKIND),dimension(:,:),pointer:: ustar        => null()
    real(kind=RKIND),dimension(:,:),pointer:: sh           => null()
    real(kind=RKIND),dimension(:,:),pointer:: z0h          => null()
    real(kind=RKIND),dimension(:,:),pointer:: cn_prcp      => null()
    real(kind=RKIND),dimension(:,:),pointer:: ncn_prcp     => null()


    !--- atmospheric fields:
    real(kind=RKIND),dimension(:,:,:),pointer:: airdens    => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: delp       => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: delz       => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: t          => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: rh2        => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: zle        => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: ple        => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: pfl_lsan   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: pfi_lsan   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: u          => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: v          => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: fcld       => null()


    contains
       procedure:: gocart2G_dims       => mpas_chemistry_gocart2G_dims
       procedure:: gocart2G_allocate   => mpas_chemistry_gocart2G_allocate
       procedure:: gocart2G_deallocate => mpas_chemistry_gocart2G_deallocate
       procedure:: gocart2G_fromMPAS   => mpas_chemistry_gocart2G_fromMPAS
       procedure:: gocart2G_toMPAS     => mpas_chemistry_gocart2G_toMPAS
 end type


 contains


!=================================================================================================================
 subroutine mpas_chemistry_gocart2G_dims(self,mesh)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 class(atm_gocart2G),intent(inout):: self

!local variables:
 integer,pointer:: nCellsSolve,nVertLevels

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_dims:')


 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 self%its = 1 ; self%ite = nCellsSolve
 self%jts = 1 ; self%jte = 1
 self%kts = 1 ; self%kte = nVertLevels ; self%ktep1 = nVertLevels+1
 call mpas_log_write('ITS = $i   ITE = $i',intArgs=(/self%its,self%ite/))
 call mpas_log_write('JTS = $i   JTE = $i',intArgs=(/self%jts,self%jte/))
 call mpas_log_write('KTS = $i   KTE = $i',intArgs=(/self%kts,self%kte/))


 call mpas_log_write('--- end subroutine mpas_chemistry_gocart2G_dims:')

 end subroutine mpas_chemistry_gocart2G_dims

!=================================================================================================================
 subroutine mpas_chemistry_gocart2G_allocate(self)
!=================================================================================================================

!inout arguments:
 class(atm_gocart2G),intent(inout):: self

!local variables and arrays:
 integer:: its,ite,jts,jte,kts,kte,ktep1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_allocate:')

 its   = self%its
 ite   = self%ite
 jts   = self%jts
 jte   = self%jte
 kts   = self%kts
 kte   = self%kte
 ktep1 = self%ktep1


!--- allocate gocart2G mixing ratios:
 if(.not.associated(self%qbcphobic) ) allocate(self%qbcphobic(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%qbcphilic) ) allocate(self%qbcphilic(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%qbrphobic) ) allocate(self%qbrphobic(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%qbrphilic) ) allocate(self%qbrphilic(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%qocphobic) ) allocate(self%qocphobic(its:ite,jts:jte,kts:kte) )
 if(.not.associated(self%qocphilic) ) allocate(self%qocphilic(its:ite,jts:jte,kts:kte) )

 if(.not.associated(self%qdust1)    ) allocate(self%qdust1(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qdust2)    ) allocate(self%qdust2(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qdust3)    ) allocate(self%qdust3(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qdust4)    ) allocate(self%qdust4(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qdust5)    ) allocate(self%qdust5(its:ite,jts:jte,kts:kte)    )

 if(.not.associated(self%nitrate1)  ) allocate(self%nitrate1(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%nitrate2)  ) allocate(self%nitrate2(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%nitrate3)  ) allocate(self%nitrate3(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%nitrate4)  ) allocate(self%nitrate4(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%nitrate5)  ) allocate(self%nitrate5(its:ite,jts:jte,kts:kte)  )

 if(.not.associated(self%qseas1)    ) allocate(self%qseas1(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qseas2)    ) allocate(self%qseas2(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qseas3)    ) allocate(self%qseas3(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%qseas4)    ) allocate(self%qseas4(its:ite,jts:jte,kts:kte)    )

 if(.not.associated(self%qso4)      ) allocate(self%qso4(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%qso2)      ) allocate(self%qso2(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%qdms)      ) allocate(self%qdms(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%qmsa)      ) allocate(self%qmsa(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%qsu2G)     ) allocate(self%qsu2G(its:ite,jts:jte,kts:kte,4)   )


!--- allocate surface emissions:
 if(.not.associated(self%qbc1_em)   ) allocate(self%qbc1_em(its:ite,jts:jte)           )
 if(.not.associated(self%qoc1_em)   ) allocate(self%qoc1_em(its:ite,jts:jte)           )
 if(.not.associated(self%qnh3_em)   ) allocate(self%qnh3_em(its:ite,jts:jte)           )
 if(.not.associated(self%qso2_em)   ) allocate(self%qso2_em(its:ite,jts:jte)           )


!--- allocate background fields:
 if(.not.associated(self%backg_dms) ) allocate(self%backg_dms(its:ite,jts:jte)         )
 if(.not.associated(self%backg_oh)  ) allocate(self%backg_oh(its:ite,jts:jte,kts:kte)  )
 if(.not.associated(self%backg_h2o2)) allocate(self%backg_h2o2(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%backg_no3) ) allocate(self%backg_no3(its:ite,jts:jte,kts:kte) )


!--- allocate mesh fields:
 if(.not.associated(self%xlat)      ) allocate(self%xlat(its:ite,jts:jte)              )
 if(.not.associated(self%xlon)      ) allocate(self%xlon(its:ite,jts:jte)              )


!--- allocate surface fields:
 if(.not.associated(self%frocean)   ) allocate(self%frocean(its:ite,jts:jte)           )
 if(.not.associated(self%fraci)     ) allocate(self%fraci(its:ite,jts:jte)             )
 if(.not.associated(self%frlake)    ) allocate(self%frlake(its:ite,jts:jte)            )
 if(.not.associated(self%lwi)       ) allocate(self%lwi(its:ite,jts:jte)               )
 if(.not.associated(self%tropp)     ) allocate(self%tropp(its:ite,jts:jte)             )
 if(.not.associated(self%u10m)      ) allocate(self%u10m(its:ite,jts:jte)              )
 if(.not.associated(self%v10m)      ) allocate(self%v10m(its:ite,jts:jte)              )
 if(.not.associated(self%area)      ) allocate(self%area(its:ite,jts:jte)              )
 if(.not.associated(self%zpbl)      ) allocate(self%zpbl(its:ite,jts:jte)              )
 if(.not.associated(self%ustar)     ) allocate(self%ustar(its:ite,jts:jte)             )
 if(.not.associated(self%sh)        ) allocate(self%sh(its:ite,jts:jte)                )
 if(.not.associated(self%z0h)       ) allocate(self%z0h(its:ite,jts:jte)               )
 if(.not.associated(self%cn_prcp)   ) allocate(self%cn_prcp(its:ite,jts:jte)           )
 if(.not.associated(self%ncn_prcp)  ) allocate(self%ncn_prcp(its:ite,jts:jte)          )


!--- allocate atmospheric fields:
 if(.not.associated(self%airdens)   ) allocate(self%airdens(its:ite,jts:jte,kts:kte)   ) 
 if(.not.associated(self%delp)      ) allocate(self%delp(its:ite,jts:jte,kts:kte)      ) 
 if(.not.associated(self%delz)      ) allocate(self%delz(its:ite,jts:jte,kts:kte)      ) 
 if(.not.associated(self%t)         ) allocate(self%t(its:ite,jts:jte,kts:kte)         ) 
 if(.not.associated(self%rh2)       ) allocate(self%rh2(its:ite,jts:jte,kts:kte)       ) 
 if(.not.associated(self%zle)       ) allocate(self%zle(its:ite,jts:jte,kts:ktep1)     ) 
 if(.not.associated(self%ple)       ) allocate(self%ple(its:ite,jts:jte,kts:ktep1)     ) 
 if(.not.associated(self%pfl_lsan)  ) allocate(self%pfl_lsan(its:ite,jts:jte,kts:kte)  ) 
 if(.not.associated(self%pfi_lsan)  ) allocate(self%pfi_lsan(its:ite,jts:jte,kts:kte)  ) 
 if(.not.associated(self%u)         ) allocate(self%u(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%v)         ) allocate(self%v(its:ite,jts:jte,kts:kte)         )
 if(.not.associated(self%fcld)      ) allocate(self%fcld(its:ite,jts:jte,kts:kte)      )

 end subroutine mpas_chemistry_gocart2G_allocate

!=================================================================================================================
 subroutine mpas_chemistry_gocart2G_deallocate(self)
!=================================================================================================================

!inout arguments:
 class(atm_gocart2G),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_deallocate:')

!--- deallocate gocart2G mixing ratios:
 if(associated(self%qbcphobic) ) deallocate(self%qbcphobic )
 if(associated(self%qbcphilic) ) deallocate(self%qbcphilic )
 if(associated(self%qbrphobic) ) deallocate(self%qbrphobic )
 if(associated(self%qbrphilic) ) deallocate(self%qbrphilic )
 if(associated(self%qocphobic) ) deallocate(self%qocphobic )
 if(associated(self%qocphilic) ) deallocate(self%qocphilic )

 if(associated(self%qdust1)    ) deallocate(self%qdust1    )   
 if(associated(self%qdust2)    ) deallocate(self%qdust2    )   
 if(associated(self%qdust3)    ) deallocate(self%qdust3    )   
 if(associated(self%qdust4)    ) deallocate(self%qdust4    )   
 if(associated(self%qdust5)    ) deallocate(self%qdust5    )   

 if(associated(self%nitrate1)  ) deallocate(self%nitrate1  )
 if(associated(self%nitrate2)  ) deallocate(self%nitrate2  )
 if(associated(self%nitrate3)  ) deallocate(self%nitrate3  )
 if(associated(self%nitrate4)  ) deallocate(self%nitrate4  )
 if(associated(self%nitrate5)  ) deallocate(self%nitrate5  )

 if(associated(self%qseas1)    ) deallocate(self%qseas1    )   
 if(associated(self%qseas2)    ) deallocate(self%qseas2    )   
 if(associated(self%qseas3)    ) deallocate(self%qseas3    )   
 if(associated(self%qseas4)    ) deallocate(self%qseas4    )   

 if(associated(self%qso4)      ) deallocate(self%qso4      )   
 if(associated(self%qso2)      ) deallocate(self%qso2      )
 if(associated(self%qdms)      ) deallocate(self%qdms      )   
 if(associated(self%qmsa)      ) deallocate(self%qmsa      )   
 if(associated(self%qsu2G)     ) deallocate(self%qsu2G     )


!--- allocate surface emissions:
 if(associated(self%qbc1_em)   ) deallocate(self%qbc1_em   )
 if(associated(self%qoc1_em)   ) deallocate(self%qoc1_em   )
 if(associated(self%qnh3_em)   ) deallocate(self%qnh3_em   )
 if(associated(self%qso2_em)   ) deallocate(self%qso2_em   )


!--- allocate background fields:
 if(associated(self%backg_dms) ) deallocate(self%backg_dms )
 if(associated(self%backg_oh)  ) deallocate(self%backg_oh  )
 if(associated(self%backg_h2o2)) deallocate(self%backg_h2o2)
 if(associated(self%backg_no3) ) deallocate(self%backg_no3 )


!--- deallocate mesh fields:
 if(associated(self%xlat)      ) deallocate(self%xlat      )
 if(associated(self%xlon)      ) deallocate(self%xlon      )


!--- deallocate surface fields:
 if(associated(self%frocean)   ) deallocate(self%frocean   )
 if(associated(self%fraci)     ) deallocate(self%fraci     )
 if(associated(self%frlake)    ) deallocate(self%frlake    )
 if(associated(self%lwi)       ) deallocate(self%lwi       )
 if(associated(self%tropp)     ) deallocate(self%tropp     )
 if(associated(self%u10m)      ) deallocate(self%u10m      )
 if(associated(self%v10m)      ) deallocate(self%v10m      )
 if(associated(self%area)      ) deallocate(self%area      )
 if(associated(self%zpbl)      ) deallocate(self%zpbl      )
 if(associated(self%ustar)     ) deallocate(self%ustar     )
 if(associated(self%sh)        ) deallocate(self%sh        )
 if(associated(self%z0h)       ) deallocate(self%z0h       )
 if(associated(self%cn_prcp)   ) deallocate(self%cn_prcp   )
 if(associated(self%ncn_prcp)  ) deallocate(self%ncn_prcp  )


!--- deallocate atmospheric fields:
 if(associated(self%airdens)   ) deallocate(self%airdens   )
 if(associated(self%delp)      ) deallocate(self%delp      )
 if(associated(self%delz)      ) deallocate(self%delz      )
 if(associated(self%t)         ) deallocate(self%t         )
 if(associated(self%rh2)       ) deallocate(self%rh2       )
 if(associated(self%zle)       ) deallocate(self%zle       )
 if(associated(self%ple)       ) deallocate(self%ple       )
 if(associated(self%pfl_lsan)  ) deallocate(self%pfl_lsan  )
 if(associated(self%pfi_lsan)  ) deallocate(self%pfi_lsan  )
 if(associated(self%u)         ) deallocate(self%u         )
 if(associated(self%v)         ) deallocate(self%v         )
 if(associated(self%fcld)      ) deallocate(self%fcld      )

 end subroutine mpas_chemistry_gocart2G_deallocate

!=================================================================================================================
 subroutine mpas_chemistry_gocart2G_fromMPAS(self,CAMS_emissions,GOCART_backgrounds,configs,mesh,diag, &
                                             state,diag_physics,sfc_input,time_lev)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: CAMS_emissions
 type(mpas_pool_type),intent(in):: GOCART_backgrounds

 integer,intent(in):: time_lev

!--- inout arguments:
 class(atm_gocart2G),intent(inout):: self

!--- local variables and pointers for gocart2G:
 integer,pointer:: num_chems
 integer,pointer:: index_qbc1,index_qbc2
 integer,pointer:: index_qoc1,index_qoc2
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4
 integer,pointer:: index_qsulf,index_qso2,index_qdms,index_qmsa
 integer:: i,its,ite,j,jts,jte,k,kts,kte,ktep1,kk

 real(kind=RKIND),dimension(:,:,:),pointer:: chems

!--- local variables and pointers for mesh, surface, and atmospheric fields:
 integer,pointer:: index_qv

 real(kind=RKIND),dimension(:),pointer:: bc1_em_anthro,oc1_em_anthro,nh3_em_anthro,so2_em_anthro
 real(kind=RKIND),dimension(:),pointer:: background_dms
 real(kind=RKIND),dimension(:,:),pointer:: background_h2o2,background_oh,background_no3

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:),pointer:: hfx,hpbl,raincv,rainncv,u10,v10,ust,z0h
 real(kind=RKIND),dimension(:),pointer:: areaCell,latCell,lonCell,xland
 real(kind=RKIND),dimension(:,:),pointer:: zgrid
 real(kind=RKIND),dimension(:,:),pointer:: airdens,qv,relhum
 real(kind=RKIND),dimension(:,:),pointer:: exner,pressure_b,pressure_p,theta,u,v
 real(kind=RKIND),dimension(:,:),pointer:: cldfrac
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND):: fzm,fzp,tem,z0,z1,z2
 real(kind=RKIND),dimension(:),allocatable:: pres,presl2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_fromMPAS:')

 its   = self%its
 ite   = self%ite
 jts   = self%jts
 jte   = self%jte
 kts   = self%kts
 kte   = self%kte
 ktep1 = self%ktep1


!--- initialization of gocart2G mixing ratios:
 call mpas_pool_get_dimension(state,'num_chems'   ,num_chems   )
 call mpas_pool_get_dimension(state,'index_qbc1'  ,index_qbc1  )
 call mpas_pool_get_dimension(state,'index_qbc2'  ,index_qbc2  )
 call mpas_pool_get_dimension(state,'index_qoc1'  ,index_qoc1  )
 call mpas_pool_get_dimension(state,'index_qoc2'  ,index_qoc2  )
 call mpas_pool_get_dimension(state,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(state,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(state,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(state,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(state,'index_qdust5',index_qdust5)
 call mpas_pool_get_dimension(state,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(state,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(state,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(state,'index_qseas4',index_qseas4)
 call mpas_pool_get_dimension(state,'index_qsulf' ,index_qsulf )
 call mpas_pool_get_dimension(state,'index_qso2'  ,index_qso2  )
 call mpas_pool_get_dimension(state,'index_qdms'  ,index_qdms  )
 call mpas_pool_get_dimension(state,'index_qmsa'  ,index_qmsa  )

 call mpas_log_write('--- num_chems    = $i',intArgs=(/num_chems/)   )
 call mpas_log_write('--- index_qbc1   = $i',intArgs=(/index_qbc1/)  )
 call mpas_log_write('--- index_qbc2   = $i',intArgs=(/index_qbc2/)  )
 call mpas_log_write('--- index_qoc1   = $i',intArgs=(/index_qoc1/)  )
 call mpas_log_write('--- index_qoc2   = $i',intArgs=(/index_qoc2/)  )
 call mpas_log_write('--- index_qdust1 = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2 = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3 = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4 = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5 = $i',intArgs=(/index_qdust5/))
 call mpas_log_write('--- index_qseas1 = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2 = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3 = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4 = $i',intArgs=(/index_qseas4/))
 call mpas_log_write('--- index_qsulf  = $i',intArgs=(/index_qsulf/) )
 call mpas_log_write('--- index_qso2   = $i',intArgs=(/index_qso2/)  )
 call mpas_log_write('--- index_qdms   = $i',intArgs=(/index_qdms/)  )
 call mpas_log_write('--- index_qmsa   = $i',intArgs=(/index_qmsa/)  )

 call mpas_pool_get_array(state,'chems',chems,time_lev)
 do k = kts,kte
    do j = jts,jte
       do i = its,ite
          self%qbcphobic(i,j,k) = chems(index_qbc1,k,i)
          self%qbcphilic(i,j,k) = chems(index_qbc2,k,i)
          self%qocphobic(i,j,k) = chems(index_qoc1,k,i)
          self%qocphilic(i,j,k) = chems(index_qoc2,k,i)
          self%qdust1(i,j,k)    = chems(index_qdust1,k,i)
          self%qdust2(i,j,k)    = chems(index_qdust2,k,i)
          self%qdust3(i,j,k)    = chems(index_qdust3,k,i)
          self%qdust4(i,j,k)    = chems(index_qdust4,k,i)
          self%qdust5(i,j,k)    = chems(index_qdust5,k,i)
          self%qseas1(i,j,k)    = chems(index_qseas1,k,i)
          self%qseas2(i,j,k)    = chems(index_qseas2,k,i)
          self%qseas3(i,j,k)    = chems(index_qseas3,k,i)
          self%qseas4(i,j,k)    = chems(index_qseas4,k,i)
          self%qso4(i,j,k)      = chems(index_qsulf,k,i)
          self%qso2(i,j,k)      = chems(index_qso2,k,i)
          self%qdms(i,j,k)      = chems(index_qdms,k,i)
          self%qmsa(i,j,k)      = chems(index_qmsa,k,i)

          self%qsu2G(i,j,k,1)   = chems(index_qdms,k,i)
          self%qsu2G(i,j,k,2)   = chems(index_qso2,k,i)
          self%qsu2G(i,j,k,3)   = chems(index_qsulf,k,i)
          self%qsu2G(i,j,k,4)   = chems(index_qmsa,k,i)

          !--- for now, set to 0._RKIND brown carbon and nitrate mixing ratios:
          self%qbrphobic(i,j,k) = 0._RKIND
          self%qbrphilic(i,j,k) = 0._RKIND
          self%nitrate1(i,j,k)  = 0._RKIND
          self%nitrate2(i,j,k)  = 0._RKIND
          self%nitrate3(i,j,k)  = 0._RKIND
          self%nitrate4(i,j,k)  = 0._RKIND
          self%nitrate5(i,j,k)  = 0._RKIND
       enddo
    enddo
 enddo


!--- initialization of emission fields:
 call mpas_pool_get_array(CAMS_emissions,'bc1_em_anthro',bc1_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'oc1_em_anthro',oc1_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'nh3_em_anthro',nh3_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'so2_em_anthro',so2_em_anthro)

 do j = jts,jte
    do i = its,ite
       self%qbc1_em(i,j) = bc1_em_anthro(i)
       self%qoc1_em(i,j) = oc1_em_anthro(i)
       self%qnh3_em(i,j) = nh3_em_anthro(i)
       self%qso2_em(i,j) = so2_em_anthro(i)
    enddo
 enddo


!--- initialization of background fields:
 call mpas_pool_get_array(GOCART_backgrounds,'background_dms' ,background_dms )
 call mpas_pool_get_array(GOCART_backgrounds,'background_h2o2',background_h2o2)
 call mpas_pool_get_array(GOCART_backgrounds,'background_oh'  ,background_oh  )
 call mpas_pool_get_array(GOCART_backgrounds,'background_no3' ,background_no3 )

 do j = jts,jte
    do i = its,ite
       self%backg_dms(i,j) = background_dms(i)

       do k = kts,kte
          self%backg_h2o2(i,j,k) = background_h2o2(k,i)
          self%backg_oh(i,j,k)   = background_oh(k,i)
          self%backg_no3(i,j,k)  = background_no3(k,i)
       enddo
    enddo
 enddo


!--- initialization of miscellaneous variables:
 call mpas_pool_get_config(configs,'config_dt',dt)
 self%dt = dt


!--- initialization of mesh fields:
 call mpas_pool_get_array(mesh,'areaCell',areaCell)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 do j = jts,jte
    do i = its,ite
       self%area(i,j) = areaCell(i)
       self%xlat(i,j) = latCell(i)
       self%xlon(i,j) = lonCell(i)
    enddo
 enddo


!--- initialization of surface fields:
 call mpas_pool_get_array(diag_physics,'u10'    ,u10    )
 call mpas_pool_get_array(diag_physics,'v10'    ,v10    )
 call mpas_pool_get_array(diag_physics,'ust'    ,ust    )
 call mpas_pool_get_array(diag_physics,'hfx'    ,hfx    )
 call mpas_pool_get_array(diag_physics,'hpbl'   ,hpbl   )
 call mpas_pool_get_array(diag_physics,'z0'     ,z0h    )
 call mpas_pool_get_array(diag_physics,'raincv' ,raincv )
 call mpas_pool_get_array(diag_physics,'rainncv',rainncv)

 call mpas_pool_get_array(sfc_input,'xland',xland)

 do j = jts,jte
    do i = its,ite
       self%frocean(i,j)  = 0._RKIND
       self%lwi(i,j)      = xland(i)
       self%tropp(i,j)    = 10000._RKIND
       self%u10m(i,j)     = u10(i)
       self%v10m(i,j)     = v10(i)
       self%zpbl(i,j)     = hpbl(i)
       self%ustar(i,j)    = ust(i)
       self%sh(i,j)       = hfx(i)
       self%z0h(i,j)      = z0h(i)
       self%cn_prcp(i,j)  = raincv(i)
       self%ncn_prcp(i,j) = rainncv(i)
    enddo 
 enddo


!--- initialization of meteorological fields:
 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(diag,'rho'          ,airdens   )
 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)
 call mpas_pool_get_array(diag,'relhum'       ,relhum    )
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,u)
 call mpas_pool_get_array(diag,'uReconstructMeridional',v)

 call mpas_pool_get_array(state,'theta_m',theta  ,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 qv => scalars(index_qv,:,:)

 call mpas_pool_get_array(diag_physics,'cldfrac',cldfrac)


 do k = kts,kte
    do j = jts,jte
       do i = its,ite
          self%delz(i,j,k)     = zgrid(k+1,i) - zgrid(k,i)
          self%airdens(i,j,k)  = airdens(k,i)*(1._RKIND+qv(k,i))
          self%t(i,j,k)        = theta(k,i)/(1._RKIND+R_v/R_d*max(0._RKIND,qv(k,i)))
          self%t(i,j,k)        = self%t(i,j,k)*exner(k,i)
          self%rh2(i,j,k)      = relhum(k,i)/100.
          self%u(i,j,k)        = u(k,i)
          self%v(i,j,k)        = v(k,i)

          self%fcld(i,j,k)     = cldfrac(k,i)
          self%pfl_lsan(i,j,k) = 0._RKIND
          self%pfi_lsan(i,j,k) = 0._RKIND
       enddo
    enddo
 enddo
 do k = kts,ktep1
    do j = jts,jte
       do i = its,ite
          self%zle(i,j,k) = zgrid(k,i)
       enddo
    enddo
 enddo


 if(.not.allocated(pres)) allocate(pres(kts:kte))
 if(.not.allocated(presl2)) allocate(presl2(kts:ktep1))
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          pres(k) = pressure_p(k,i) + pressure_b(k,i)
       enddo

       !pressure at the surface:
       k = kts
       z0  = zgrid(k,i)
       z1  = 0.5*(zgrid(k,i)+zgrid(k+1,i))
       z2  = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
       fzp = (z0-z2)/(z1-z2)
       fzm = 1.-fzp
       kk = ktep1+1-kts
       presl2(k) = fzp*pres(k) + fzm*pres(k+1)

       !pressure at interface between layers above the surface (k=kts) and below the model-top (k=kte+1):
       do k = kts+1,kte
          tem = 1./(zgrid(k+1,i)-zgrid(k-1,i))
          fzm = (zgrid(k,i)-zgrid(k-1,i))*tem
          fzp = (zgrid(k+1,i)-zgrid(k,i))*tem
          presl2(k) = fzm*pres(k) + fzp*pres(k-1)
       enddo
       !pressure at the model-top (k=kte+1):
       k = kte+1
       z0  = zgrid(k,i)
       z1  = 0.5*(zgrid(k,i)+zgrid(k-1,i))
       z2  = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
       fzm = (z0-z2)/(z1-z2)
       fzp = 1.-fzm
       !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
       presl2(k) = exp(fzm*log(pres(k-1)) + fzp*log(pres(k-2)))

       do k = kts,ktep1
          self%ple(i,j,k) = presl2(k)
       enddo
    enddo

    !pressure thickness:
    do i = its,ite
       do k = kts,kte
          self%delp(i,j,k) = self%ple(i,j,k)-self%ple(i,j,k+1)
       enddo
    enddo
 enddo
 if(allocated(pres)) deallocate(pres)
 if(allocated(presl2)) deallocate(presl2)


 call mpas_log_write('--- end subroutine mpas_chemistry_from_gocart2G:')

 end subroutine mpas_chemistry_gocart2G_fromMPAS

!=================================================================================================================
 subroutine mpas_chemistry_gocart2G_toMPAS(self,state,time_lev)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: state
 integer,intent(in):: time_lev

!--- inout arguments:
 class(atm_gocart2G),intent(inout):: self

!--- local variables and pointers:
 integer,pointer:: num_chems
 integer,pointer:: index_qbc1,index_qbc2
 integer,pointer:: index_qoc1,index_qoc2
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4
 integer,pointer:: index_qsulf,index_qdms,index_qmsa
 integer:: i,its,ite,j,jts,jte,k,kts,kte

 real(kind=RKIND),dimension(:,:,:),pointer:: chems

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart2G_toMPAS:')

 its = self%its
 ite = self%ite
 jts = self%jts
 jte = self%jte
 kts = self%kts
 kte = self%kte


!--- initialize gocart2G mixing ratios:
 call mpas_pool_get_dimension(state,'num_chems'   ,num_chems   )
 call mpas_pool_get_dimension(state,'index_qbc1'  ,index_qbc1  )
 call mpas_pool_get_dimension(state,'index_qbc2'  ,index_qbc2  )
 call mpas_pool_get_dimension(state,'index_qoc1'  ,index_qoc1  )
 call mpas_pool_get_dimension(state,'index_qoc2'  ,index_qoc2  )
 call mpas_pool_get_dimension(state,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(state,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(state,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(state,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(state,'index_qdust5',index_qdust5)
 call mpas_pool_get_dimension(state,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(state,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(state,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(state,'index_qseas4',index_qseas4)
 call mpas_pool_get_dimension(state,'index_qsulf' ,index_qsulf )
 call mpas_pool_get_dimension(state,'index_qdms'  ,index_qdms  )
 call mpas_pool_get_dimension(state,'index_qmsa'  ,index_qmsa  )


 call mpas_pool_get_array(state,'chems',chems,time_lev)
 do k = kts,kte
    do j = jts,jte
       do i = its,ite
          chems(index_qbc1,k,i)   = self%qbcphobic(i,j,k)
          chems(index_qbc2,k,i)   = self%qbcphilic(i,j,k)
          chems(index_qoc1,k,i)   = self%qbcphobic(i,j,k)
          chems(index_qoc2,k,i)   = self%qbcphilic(i,j,k)
          chems(index_qdust1,k,i) = self%qdust1(i,j,k)
          chems(index_qdust2,k,i) = self%qdust2(i,j,k)
          chems(index_qdust3,k,i) = self%qdust3(i,j,k)
          chems(index_qdust4,k,i) = self%qdust4(i,j,k)
          chems(index_qdust5,k,i) = self%qdust5(i,j,k)
          chems(index_qseas1,k,i) = self%qseas1(i,j,k)
          chems(index_qseas2,k,i) = self%qseas2(i,j,k)
          chems(index_qseas3,k,i) = self%qseas3(i,j,k)
          chems(index_qseas4,k,i) = self%qseas4(i,j,k)
          chems(index_qsulf,k,i)  = self%qso4(i,j,k)
          chems(index_qdms,k,i)   = self%qdms(i,j,k)
          chems(index_qmsa,k,i)   = self%qmsa(i,j,k)
       enddo
    enddo
 enddo

 call mpas_log_write('--- end subroutine mpas_chemistry_gocart2G_toMPAS:')

 end subroutine mpas_chemistry_gocart2G_toMPAS

!=================================================================================================================
 end module mpas_chemistry_gocart2G_interface
!=================================================================================================================

