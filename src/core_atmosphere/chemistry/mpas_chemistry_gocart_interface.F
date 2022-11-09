!=================================================================================================================
 module mpas_chemistry_gocart_interface
 use mpas_log
 use mpas_kind_types
 use mpas_constants,only: gravity
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension

 implicit none


 type,public:: chems_gocart
    integer:: its,ite,kts,kte,ktep1
    integer:: nerod

    real(kind=RKIND),dimension(:),pointer:: qbc1_em       => null()
    real(kind=RKIND),dimension(:),pointer:: qoc1_em       => null()
    real(kind=RKIND),dimension(:),pointer:: qco_em        => null()
    real(kind=RKIND),dimension(:),pointer:: qnh3_em       => null()
    real(kind=RKIND),dimension(:),pointer:: qso2_em       => null()

    real(kind=RKIND),dimension(:),pointer:: emis_qdust1   => null()
    real(kind=RKIND),dimension(:),pointer:: emis_qdust2   => null()
    real(kind=RKIND),dimension(:),pointer:: emis_qdust3   => null()
    real(kind=RKIND),dimension(:),pointer:: emis_qdust4   => null()
    real(kind=RKIND),dimension(:),pointer:: emis_qdust5   => null()

    real(kind=RKIND),dimension(:,:),pointer:: qbc1        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qbc2        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc1        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc2        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas3      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas4      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust3      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust4      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust5      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qco         => null()
    real(kind=RKIND),dimension(:,:),pointer:: qnh3        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso2        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso4a1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso4a2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso4a3      => null()

    real(kind=RKIND),dimension(:,:),pointer:: tend_qbc1   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qbc2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qoc1   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qoc2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas1 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas2 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas3 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas4 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust1 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust2 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust3 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust4 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust5 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qco    => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qnh3   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qso2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qso4a1 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qso4a2 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qso4a3 => null()


    !--- local "water-friendly" and "ice-friendly" Thompson aerosols:
    real(kind=RKIND),dimension(:),pointer:: nwfa2d_c      => null()
    real(kind=RKIND),dimension(:),pointer:: nifa2d_c      => null()

    real(kind=RKIND),dimension(:,:),pointer:: nwfa_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: nifa_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qwfa_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qifa_c      => null()

    real(kind=RKIND),dimension(:,:),pointer:: tend_nwfa_c => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_nifa_c => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qwfa_c => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qifa_c => null()


    !--- local variables from the MPAS dynamical core:
    integer,dimension(:),pointer           :: isltyp_c    => null()

    real(kind=RKIND):: dt_c
    real(kind=RKIND):: grav_c
    real(kind=RKIND),dimension(:),pointer  :: clay_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: sand_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: porosity_c  => null()
    real(kind=RKIND),dimension(:),pointer  :: area_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: xland_c     => null()
    real(kind=RKIND),dimension(:),pointer  :: u10_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: v10_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: smois_c     => null()
    real(kind=RKIND),dimension(:),pointer  :: snowh_c     => null()
    real(kind=RKIND),dimension(:),pointer  :: ust_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: znt_c       => null()
    real(kind=RKIND),dimension(:,:),pointer:: erod_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: zgrid_c     => null()
    real(kind=RKIND),dimension(:,:),pointer:: rho_c       => null()
    real(kind=RKIND),dimension(:,:),pointer:: pres_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: dpres_c     => null()
    real(kind=RKIND),dimension(:,:),pointer:: u_c         => null()
    real(kind=RKIND),dimension(:,:),pointer:: v_c         => null()


    contains
       procedure:: gocart_dims       => mpas_chemistry_gocart_dims
       procedure:: gocart_allocate   => mpas_chemistry_gocart_allocate
       procedure:: gocart_deallocate => mpas_chemistry_gocart_deallocate
       procedure:: from_gocart       => mpas_chemistry_from_gocart
       procedure:: to_gocart         => mpas_chemistry_to_gocart
 end type


 contains


!=================================================================================================================
 subroutine mpas_chemistry_gocart_dims(this,mesh)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 class(chems_gocart),intent(inout):: this

!local variables:
 integer,pointer:: nCellsSolve,nVertLevels,nDustErosion

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart_dims:')

 call mpas_pool_get_dimension(mesh,'nCellsSolve'    ,nCellsSolve )
 call mpas_pool_get_dimension(mesh,'nVertLevels'    ,nVertLevels )
 call mpas_pool_get_dimension(mesh,'nDustErosionDim',nDustErosion)

 this%its = 1 ; this%ite = nCellsSolve
 this%kts = 1 ; this%kte = nVertLevels ; this%ktep1 = nVertLevels+1
 call mpas_log_write('ITS = $i   ITE = $i',intArgs=(/this%its,this%ite/))
 call mpas_log_write('KTS = $i   KTE = $i',intArgs=(/this%kts,this%kte/))

 this%nerod = nDustErosion

 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_dims:')

 end subroutine mpas_chemistry_gocart_dims

!=================================================================================================================
 subroutine mpas_chemistry_gocart_allocate(this)
!=================================================================================================================

!inout arguments:
 class(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: its,ite,kts,kte,ktep1,nerod

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart_allocate:')

 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1
 nerod = this%nerod


!--- allocate surface emissions:
 if(.not.associated(this%qbc1_em)    ) allocate(this%qbc1_em(its:ite)            )
 if(.not.associated(this%qoc1_em)    ) allocate(this%qoc1_em(its:ite)            )
 if(.not.associated(this%qco_em)     ) allocate(this%qco_em(its:ite)             )
 if(.not.associated(this%qnh3_em)    ) allocate(this%qnh3_em(its:ite)            )
 if(.not.associated(this%qso2_em)    ) allocate(this%qso2_em(its:ite)            )


!--- allocate mixing ratios:
 if(.not.associated(this%qbc1)       ) allocate(this%qbc1(its:ite,kts:kte)       )
 if(.not.associated(this%qbc2)       ) allocate(this%qbc2(its:ite,kts:kte)       )
 if(.not.associated(this%qoc1)       ) allocate(this%qoc1(its:ite,kts:kte)       )
 if(.not.associated(this%qoc2)       ) allocate(this%qoc2(its:ite,kts:kte)       )
 if(.not.associated(this%qseas1)     ) allocate(this%qseas1(its:ite,kts:kte)     )
 if(.not.associated(this%qseas2)     ) allocate(this%qseas2(its:ite,kts:kte)     )
 if(.not.associated(this%qseas3)     ) allocate(this%qseas3(its:ite,kts:kte)     )
 if(.not.associated(this%qseas4)     ) allocate(this%qseas4(its:ite,kts:kte)     )
 if(.not.associated(this%qdust1)     ) allocate(this%qdust1(its:ite,kts:kte)     )
 if(.not.associated(this%qdust2)     ) allocate(this%qdust2(its:ite,kts:kte)     )
 if(.not.associated(this%qdust3)     ) allocate(this%qdust3(its:ite,kts:kte)     )
 if(.not.associated(this%qdust4)     ) allocate(this%qdust4(its:ite,kts:kte)     )
 if(.not.associated(this%qdust5)     ) allocate(this%qdust5(its:ite,kts:kte)     )

 if(.not.associated(this%qco)        ) allocate(this%qco(its:ite,kts:kte)        )
 if(.not.associated(this%qnh3)       ) allocate(this%qnh3(its:ite,kts:kte)       )
 if(.not.associated(this%qso2)       ) allocate(this%qso2(its:ite,kts:kte)       )
 if(.not.associated(this%qso4a1)     ) allocate(this%qso4a1(its:ite,kts:kte)     )
 if(.not.associated(this%qso4a2)     ) allocate(this%qso4a2(its:ite,kts:kte)     )
 if(.not.associated(this%qso4a3)     ) allocate(this%qso4a3(its:ite,kts:kte)     )


!--- allocate diagnostics:
 if(.not.associated(this%emis_qdust1)) allocate(this%emis_qdust1(its:ite)        )
 if(.not.associated(this%emis_qdust2)) allocate(this%emis_qdust2(its:ite)        )
 if(.not.associated(this%emis_qdust3)) allocate(this%emis_qdust3(its:ite)        )
 if(.not.associated(this%emis_qdust4)) allocate(this%emis_qdust4(its:ite)        )
 if(.not.associated(this%emis_qdust5)) allocate(this%emis_qdust5(its:ite)        )


!--- allocate tendencies:
 if(.not.associated(this%tend_qbc1)  ) allocate(this%tend_qbc1(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qbc2)  ) allocate(this%tend_qbc2(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qoc1)  ) allocate(this%tend_qoc1(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qoc2)  ) allocate(this%tend_qoc2(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qseas1)) allocate(this%tend_qseas1(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas2)) allocate(this%tend_qseas2(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas3)) allocate(this%tend_qseas3(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas4)) allocate(this%tend_qseas4(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust1)) allocate(this%tend_qdust1(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust2)) allocate(this%tend_qdust2(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust3)) allocate(this%tend_qdust3(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust4)) allocate(this%tend_qdust4(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust5)) allocate(this%tend_qdust5(its:ite,kts:kte))

 if(.not.associated(this%tend_qco)   ) allocate(this%tend_qco(its:ite,kts:kte)   )
 if(.not.associated(this%tend_qnh3)  ) allocate(this%tend_qnh3(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qso2)  ) allocate(this%tend_qso2(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qso4a1)) allocate(this%tend_qso4a1(its:ite,kts:kte))
 if(.not.associated(this%tend_qso4a2)) allocate(this%tend_qso4a2(its:ite,kts:kte))
 if(.not.associated(this%tend_qso4a3)) allocate(this%tend_qso4a3(its:ite,kts:kte))


!--- allocate local "water-friendly" and "ice-friendly" Thompson aerosols:
 if(.not.associated(this%nifa2d_c)   ) allocate(this%nifa2d_c(its:ite)           )
 if(.not.associated(this%nwfa2d_c)   ) allocate(this%nwfa2d_c(its:ite)           )

 if(.not.associated(this%nifa_c)     ) allocate(this%nifa_c(its:ite,kts:kte)     )
 if(.not.associated(this%nwfa_c)     ) allocate(this%nwfa_c(its:ite,kts:kte)     )
 if(.not.associated(this%qifa_c)     ) allocate(this%qifa_c(its:ite,kts:kte)     )
 if(.not.associated(this%qwfa_c)     ) allocate(this%qwfa_c(its:ite,kts:kte)     )

 if(.not.associated(this%tend_nifa_c)) allocate(this%tend_nifa_c(its:ite,kts:kte))
 if(.not.associated(this%tend_nwfa_c)) allocate(this%tend_nwfa_c(its:ite,kts:kte))
 if(.not.associated(this%tend_qifa_c)) allocate(this%tend_qifa_c(its:ite,kts:kte))
 if(.not.associated(this%tend_qwfa_c)) allocate(this%tend_qwfa_c(its:ite,kts:kte))


!--- allocate mesh and meteorological fields:
 if(.not.associated(this%isltyp_c)   ) allocate(this%isltyp_c(its:ite)           )
 if(.not.associated(this%clay_c)     ) allocate(this%clay_c(its:ite)             )
 if(.not.associated(this%sand_c)     ) allocate(this%sand_c(its:ite)             )
 if(.not.associated(this%porosity_c) ) allocate(this%porosity_c(its:ite)         )
 if(.not.associated(this%area_c)     ) allocate(this%area_c(its:ite)             )
 if(.not.associated(this%xland_c)    ) allocate(this%xland_c(its:ite)            )
 if(.not.associated(this%u10_c)      ) allocate(this%u10_c(its:ite)              )
 if(.not.associated(this%v10_c)      ) allocate(this%v10_c(its:ite)              )
 if(.not.associated(this%smois_c)    ) allocate(this%smois_c(its:ite)            )
 if(.not.associated(this%snowh_c)    ) allocate(this%snowh_c(its:ite)            )
 if(.not.associated(this%ust_c)      ) allocate(this%ust_c(its:ite)              )
 if(.not.associated(this%znt_c)      ) allocate(this%znt_c(its:ite)              )
 if(.not.associated(this%zgrid_c)    ) allocate(this%zgrid_c(its:ite,kts:ktep1)  )
 if(.not.associated(this%rho_c)      ) allocate(this%rho_c(its:ite,kts:kte)      )
 if(.not.associated(this%pres_c)     ) allocate(this%pres_c(its:ite,kts:kte)     )
 if(.not.associated(this%dpres_c)    ) allocate(this%dpres_c(its:ite,kts:kte)    )
 if(.not.associated(this%u_c)        ) allocate(this%u_c(its:ite,kts:kte)        )
 if(.not.associated(this%v_c)        ) allocate(this%v_c(its:ite,kts:kte)        )

 if(.not.associated(this%erod_c)     ) allocate(this%erod_c(its:ite,nerod)       )

 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_allocate:')

 end subroutine mpas_chemistry_gocart_allocate

!=================================================================================================================
 subroutine mpas_chemistry_gocart_deallocate(this)
!=================================================================================================================

!inout arguments:
 class(chems_gocart),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart_deallocate:')

!--- deallocate surface emissions:
 if(associated(this%qbc1_em)    ) deallocate(this%qbc1_em    )
 if(associated(this%qoc1_em)    ) deallocate(this%qoc1_em    )
 if(associated(this%qco_em)     ) deallocate(this%qco_em     )
 if(associated(this%qnh3_em)    ) deallocate(this%qnh3_em    )
 if(associated(this%qso2_em)    ) deallocate(this%qso2_em    )


!--- deallocate mixing ratios:
 if(associated(this%qbc1)       ) deallocate(this%qbc1       )
 if(associated(this%qbc2)       ) deallocate(this%qbc2       )
 if(associated(this%qoc1)       ) deallocate(this%qoc1       )
 if(associated(this%qoc2)       ) deallocate(this%qoc2       )
 if(associated(this%qseas1)     ) deallocate(this%qseas1     )
 if(associated(this%qseas2)     ) deallocate(this%qseas2     )
 if(associated(this%qseas3)     ) deallocate(this%qseas3     )
 if(associated(this%qseas4)     ) deallocate(this%qseas4     )
 if(associated(this%qdust1)     ) deallocate(this%qdust1     )
 if(associated(this%qdust2)     ) deallocate(this%qdust2     )
 if(associated(this%qdust3)     ) deallocate(this%qdust3     )
 if(associated(this%qdust4)     ) deallocate(this%qdust4     )
 if(associated(this%qdust5)     ) deallocate(this%qdust5     )

 if(associated(this%qco)        ) deallocate(this%qco        )
 if(associated(this%qnh3)       ) deallocate(this%qnh3       )
 if(associated(this%qso2)       ) deallocate(this%qso2       )
 if(associated(this%qso4a1)     ) deallocate(this%qso4a1     )
 if(associated(this%qso4a2)     ) deallocate(this%qso4a2     )
 if(associated(this%qso4a3)     ) deallocate(this%qso4a3     )


!--- deallocate diagnostics:
 if(associated(this%emis_qdust1)) deallocate(this%emis_qdust1)
 if(associated(this%emis_qdust2)) deallocate(this%emis_qdust2)
 if(associated(this%emis_qdust3)) deallocate(this%emis_qdust3)
 if(associated(this%emis_qdust4)) deallocate(this%emis_qdust4)
 if(associated(this%emis_qdust5)) deallocate(this%emis_qdust5)


!--- allocate tendencies:
 if(associated(this%tend_qbc1)  ) deallocate(this%tend_qbc1  )
 if(associated(this%tend_qbc2)  ) deallocate(this%tend_qbc2  )
 if(associated(this%tend_qoc1)  ) deallocate(this%tend_qoc1  )
 if(associated(this%tend_qoc2)  ) deallocate(this%tend_qoc2  )
 if(associated(this%tend_qseas1)) deallocate(this%tend_qseas1)
 if(associated(this%tend_qseas2)) deallocate(this%tend_qseas2)
 if(associated(this%tend_qseas3)) deallocate(this%tend_qseas3)
 if(associated(this%tend_qseas4)) deallocate(this%tend_qseas4)
 if(associated(this%tend_qdust1)) deallocate(this%tend_qdust1)
 if(associated(this%tend_qdust2)) deallocate(this%tend_qdust2)
 if(associated(this%tend_qdust3)) deallocate(this%tend_qdust3)
 if(associated(this%tend_qdust4)) deallocate(this%tend_qdust4)
 if(associated(this%tend_qdust5)) deallocate(this%tend_qdust5)

 if(associated(this%tend_qco)   ) deallocate(this%tend_qco   )
 if(associated(this%tend_qnh3)  ) deallocate(this%tend_qnh3  )
 if(associated(this%tend_qso2)  ) deallocate(this%tend_qso2  )
 if(associated(this%tend_qso4a1)) deallocate(this%tend_qso4a1)
 if(associated(this%tend_qso4a2)) deallocate(this%tend_qso4a2)
 if(associated(this%tend_qso4a3)) deallocate(this%tend_qso4a3)


!--- deallocate local "water-friendly" and "ice-friendly" Thompson aerosols:
 if(associated(this%nifa2d_c)   ) deallocate(this%nifa2d_c   )
 if(associated(this%nwfa2d_c)   ) deallocate(this%nwfa2d_c   )

 if(associated(this%nifa_c)     ) deallocate(this%nifa_c     )
 if(associated(this%nwfa_c)     ) deallocate(this%nwfa_c     )
 if(associated(this%qifa_c)     ) deallocate(this%qifa_c     )
 if(associated(this%qwfa_c)     ) deallocate(this%qwfa_c     )

 if(associated(this%tend_nifa_c)) deallocate(this%tend_nifa_c)
 if(associated(this%tend_nwfa_c)) deallocate(this%tend_nwfa_c)
 if(associated(this%tend_qifa_c)) deallocate(this%tend_qifa_c)
 if(associated(this%tend_qwfa_c)) deallocate(this%tend_qwfa_c)


!--- deallocate mesh and meteorological fields:
 if(associated(this%isltyp_c)   ) deallocate(this%isltyp_c   )
 if(associated(this%clay_c)     ) deallocate(this%clay_c     )
 if(associated(this%sand_c)     ) deallocate(this%sand_c     )
 if(associated(this%porosity_c) ) deallocate(this%porosity_c )
 if(associated(this%area_c)     ) deallocate(this%area_c     )
 if(associated(this%xland_c)    ) deallocate(this%xland_c    )
 if(associated(this%u10_c)      ) deallocate(this%u10_c      )
 if(associated(this%v10_c)      ) deallocate(this%v10_c      )
 if(associated(this%smois_c)    ) deallocate(this%smois_c    )
 if(associated(this%snowh_c)    ) deallocate(this%snowh_c    )
 if(associated(this%ust_c)      ) deallocate(this%ust_c      )
 if(associated(this%znt_c)      ) deallocate(this%znt_c      )
 if(associated(this%zgrid_c)    ) deallocate(this%zgrid_c    )
 if(associated(this%rho_c)      ) deallocate(this%rho_c      )
 if(associated(this%pres_c)     ) deallocate(this%pres_c     )
 if(associated(this%dpres_c)    ) deallocate(this%dpres_c    )
 if(associated(this%u_c)        ) deallocate(this%u_c        )
 if(associated(this%v_c)        ) deallocate(this%v_c        )

 if(associated(this%erod_c)     ) deallocate(this%erod_c     )

 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_deallocate:')

 end subroutine mpas_chemistry_gocart_deallocate

!=================================================================================================================
 subroutine mpas_chemistry_from_gocart(this,CAMS_emissions,configs,mesh,chems_input,state,time_lev,diag, &
                                       diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: chems_input
 type(mpas_pool_type),intent(in):: CAMS_emissions

 integer,intent(in):: time_lev

!inout arguments:
 class(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,n,nerod
 integer,pointer:: index_qbc1,index_qbc2,index_qoc1,index_qoc2,index_qseas1,index_qseas2,index_qseas3,  &
              index_qseas4,index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5,index_qnh3, &
              index_qso2,index_qso4a1,index_qso4a2,index_qso4a3

 real(kind=RKIND),dimension(:),pointer:: bc1_em_anthro,oc1_em_anthro,co_em_anthro, &
                                         nh3_em_anthro,so2_em_anthro
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2,qseas1,qseas2,qseas3,qseas4,  &
                                           qdust1,qdust2,qdust3,qdust4,qdust5,qnh3,qso2,     &
                                           qso4a1,qso4a2,qso4a3
 real(kind=RKIND),dimension(:,:,:),pointer:: chems_gocart

!local variables and arrays for meteorological fields:
 integer,pointer:: index_qv
 integer,pointer:: index_qco
 integer,dimension(:),pointer:: isltyp

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:),pointer:: area,clay,sand,xland
 real(kind=RKIND),dimension(:),pointer:: u10,v10
 real(kind=RKIND),dimension(:),pointer:: porosity,snowh,ust,znt
 real(kind=RKIND),dimension(:,:),pointer:: erod
 real(kind=RKIND),dimension(:,:),pointer:: zgrid,zz
 real(kind=RKIND),dimension(:,:),pointer:: pres_p,pres_b,qv,qco,rho_zz,uzonal,umerid
 real(kind=RKIND),dimension(:,:),pointer:: smois
 real(kind=RKIND),dimension(:,:,:),pointer:: chems
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND):: fzm,fzp,tem,z0,z1,z2
 real(kind=RKIND),dimension(:,:),pointer:: dpres,pres
 real(kind=RKIND),dimension(:,:),allocatable:: pres2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_from_gocart:')

 call mpas_pool_get_array(CAMS_emissions,'bc1_em_anthro',bc1_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'oc1_em_anthro',oc1_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'co_em_anthro' ,co_em_anthro )
 call mpas_pool_get_array(CAMS_emissions,'nh3_em_anthro',nh3_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'so2_em_anthro',so2_em_anthro)

 call mpas_pool_get_dimension(chems_input,'index_qbc1'  ,index_qbc1  )
 call mpas_pool_get_dimension(chems_input,'index_qbc2'  ,index_qbc2  )
 call mpas_pool_get_dimension(chems_input,'index_qoc1'  ,index_qoc1  )
 call mpas_pool_get_dimension(chems_input,'index_qoc2'  ,index_qoc2  )
 call mpas_pool_get_dimension(chems_input,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(chems_input,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(chems_input,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(chems_input,'index_qseas4',index_qseas4)
 call mpas_pool_get_dimension(chems_input,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(chems_input,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(chems_input,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(chems_input,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(chems_input,'index_qdust5',index_qdust5)
 call mpas_pool_get_dimension(chems_input,'index_qnh3'  ,index_qnh3  )
 call mpas_pool_get_dimension(chems_input,'index_qso2'  ,index_qso2  )
 call mpas_pool_get_dimension(chems_input,'index_qso4a1',index_qso4a1)
 call mpas_pool_get_dimension(chems_input,'index_qso4a2',index_qso4a2)
 call mpas_pool_get_dimension(chems_input,'index_qso4a3',index_qso4a3)

 call mpas_log_write('--- index_qbc1   = $i',intArgs=(/index_qbc1/)  )
 call mpas_log_write('--- index_qbc2   = $i',intArgs=(/index_qbc2/)  )
 call mpas_log_write('--- index_qoc1   = $i',intArgs=(/index_qoc1/)  )
 call mpas_log_write('--- index_qoc2   = $i',intArgs=(/index_qoc2/)  )
 call mpas_log_write('--- index_qseas1 = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2 = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3 = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4 = $i',intArgs=(/index_qseas4/))
 call mpas_log_write('--- index_qdust1 = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2 = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3 = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4 = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5 = $i',intArgs=(/index_qdust5/))
 call mpas_log_write('--- index_qnh3   = $i',intArgs=(/index_qnh3/)  )
 call mpas_log_write('--- index_qso2   = $i',intArgs=(/index_qso2/)  )
 call mpas_log_write('--- index_qso4a1 = $i',intArgs=(/index_qso4a1/))
 call mpas_log_write('--- index_qso4a2 = $i',intArgs=(/index_qso4a2/))
 call mpas_log_write('--- index_qso4a3 = $i',intArgs=(/index_qso4a3/))

 call mpas_pool_get_array(chems_input,'chems_gocart',chems_gocart)
 qbc1   => chems_gocart(index_qbc1,:,:)
 qbc2   => chems_gocart(index_qbc2,:,:)
 qoc1   => chems_gocart(index_qoc1,:,:)
 qoc2   => chems_gocart(index_qoc2,:,:)
 qseas1 => chems_gocart(index_qseas1,:,:)
 qseas2 => chems_gocart(index_qseas2,:,:)
 qseas3 => chems_gocart(index_qseas3,:,:)
 qseas4 => chems_gocart(index_qseas4,:,:)
 qdust1 => chems_gocart(index_qdust1,:,:)
 qdust2 => chems_gocart(index_qdust2,:,:)
 qdust3 => chems_gocart(index_qdust3,:,:)
 qdust4 => chems_gocart(index_qdust4,:,:)
 qdust5 => chems_gocart(index_qdust5,:,:)
 qnh3   => chems_gocart(index_qnh3,:,:)
 qso2   => chems_gocart(index_qso2,:,:)
 qso4a1 => chems_gocart(index_qso4a1,:,:)
 qso4a2 => chems_gocart(index_qso4a2,:,:)
 qso4a3 => chems_gocart(index_qso4a3,:,:)

 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 nerod = this%nerod

 do i = its,ite
    this%qbc1_em(i) = bc1_em_anthro(i)
    this%qoc1_em(i) = oc1_em_anthro(i)
    this%qco_em(i)  = co_em_anthro(i)
    this%qnh3_em(i) = nh3_em_anthro(i)
    this%qso2_em(i) = so2_em_anthro(i)

    do k = kts,kte
       this%qbc1(i,k)   = qbc1(k,i)
       this%qbc2(i,k)   = qbc2(k,i)
       this%qoc1(i,k)   = qoc1(k,i)
       this%qoc2(i,k)   = qoc2(k,i)
       this%qseas1(i,k) = qseas1(k,i)
       this%qseas2(i,k) = qseas2(k,i)
       this%qseas3(i,k) = qseas3(k,i)
       this%qseas4(i,k) = qseas4(k,i)
       this%qdust1(i,k) = qdust1(k,i)
       this%qdust2(i,k) = qdust2(k,i)
       this%qdust3(i,k) = qdust3(k,i)
       this%qdust4(i,k) = qdust4(k,i)
       this%qdust5(i,k) = qdust5(k,i)
       this%qnh3(i,k)   = qnh3(k,i)
       this%qso2(i,k)   = qso2(k,i)
       this%qso4a1(i,k) = qso4a1(k,i)
       this%qso4a2(i,k) = qso4a2(k,i)
       this%qso4a3(i,k) = qso4a3(k,i)
    enddo
 enddo


!--- initialization of diagnostics:
 do i = its,ite
    this%emis_qdust1(i) = 0._RKIND
    this%emis_qdust2(i) = 0._RKIND
    this%emis_qdust3(i) = 0._RKIND
    this%emis_qdust4(i) = 0._RKIND
    this%emis_qdust5(i) = 0._RKIND
 enddo


!--- initialization of tendencies before call to the GOCART model. Because mixing ratios
!    are only modified in the layer adjacent to the surface, all tendencies are set to
!    zero in the layers above k = kts:
 do i = its,ite
    k = kts
    this%tend_qbc1(i,k)   = this%qbc1(i,k)
    this%tend_qbc2(i,k)   = this%qbc2(i,k)
    this%tend_qoc1(i,k)   = this%qoc1(i,k)
    this%tend_qoc2(i,k)   = this%qoc2(i,k)
    this%tend_qseas1(i,k) = this%qseas1(i,k)
    this%tend_qseas2(i,k) = this%qseas2(i,k)
    this%tend_qseas3(i,k) = this%qseas3(i,k)
    this%tend_qseas4(i,k) = this%qseas4(i,k)
    this%tend_qdust1(i,k) = this%qdust1(i,k)
    this%tend_qdust2(i,k) = this%qdust2(i,k)
    this%tend_qdust3(i,k) = this%qdust3(i,k)
    this%tend_qdust4(i,k) = this%qdust4(i,k)
    this%tend_qdust5(i,k) = this%qdust5(i,k)
    this%tend_qnh3(i,k)   = this%qnh3(i,k)
    this%tend_qso2(i,k)   = this%qso2(i,k)
    this%tend_qso4a1(i,k) = this%qso4a1(i,k)
    this%tend_qso4a2(i,k) = this%qso4a2(i,k)
    this%tend_qso4a3(i,k) = this%qso4a3(i,k)

    do k = kts+1,kte
       this%tend_qbc1(i,k)   = 0._RKIND
       this%tend_qbc2(i,k)   = 0._RKIND
       this%tend_qoc1(i,k)   = 0._RKIND
       this%tend_qoc2(i,k)   = 0._RKIND
       this%tend_qseas1(i,k) = 0._RKIND
       this%tend_qseas2(i,k) = 0._RKIND
       this%tend_qseas3(i,k) = 0._RKIND
       this%tend_qseas4(i,k) = 0._RKIND
       this%tend_qdust1(i,k) = 0._RKIND
       this%tend_qdust2(i,k) = 0._RKIND
       this%tend_qdust3(i,k) = 0._RKIND
       this%tend_qdust4(i,k) = 0._RKIND
       this%tend_qdust5(i,k) = 0._RKIND
       this%tend_qnh3(i,k)   = 0._RKIND
       this%tend_qso2(i,k)   = 0._RKIND
       this%tend_qso4a1(i,k) = 0._RKIND
       this%tend_qso4a2(i,k) = 0._RKIND
       this%tend_qso4a3(i,k) = 0._RKIND
    enddo
 enddo


!--- initialization local "water-friendly" and "ice-friendly" Thompson aerosols:
 do i = its,ite
    do k = kts,kte
       this%tend_nifa_c(i,k) = 0._RKIND
       this%tend_nwfa_c(i,k) = 0._RKIND
       this%tend_qifa_c(i,k) = 0._RKIND
       this%tend_qwfa_c(i,k) = 0._RKIND
   enddo
 enddo


!--- initialization of meterological fields:
 call mpas_pool_get_config(configs,'config_dt',dt)

 call mpas_pool_get_array(mesh,'areaCell',area  )
 call mpas_pool_get_array(mesh,'zgrid'   ,zgrid )
 call mpas_pool_get_array(mesh,'zz'      ,zz    )
 call mpas_pool_get_array(mesh,'clayfrac_gocart',clay)
 call mpas_pool_get_array(mesh,'sandfrac_gocart',sand)
 call mpas_pool_get_array(mesh,'erod_gocart'    ,erod)

 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_dimension(state,'index_qco',index_qco)
 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 call mpas_pool_get_array(state,'chems'  ,chems  ,time_lev)
 qv  => scalars(index_qv,:,:)
 qco => chems(index_qco,:,:)

 call mpas_log_write(' ')
 call mpas_log_write('--- index_qv     = $i',intArgs=(/index_qv/))
 call mpas_log_write('--- index_qco    = $i',intArgs=(/index_qco/))
 call mpas_log_write(' ')

 call mpas_pool_get_array(sfc_input,'isltyp'  ,isltyp  )
 call mpas_pool_get_array(sfc_input,'porosity',porosity)
 call mpas_pool_get_array(sfc_input,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'xland'   ,xland   )

 call mpas_pool_get_array(diag,'pressure_base'         ,pres_b)
 call mpas_pool_get_array(diag,'pressure_p'            ,pres_p)
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,uzonal)
 call mpas_pool_get_array(diag,'uReconstructMeridional',umerid)

 call mpas_pool_get_array(diag_physics,'u10',u10)
 call mpas_pool_get_array(diag_physics,'v10',v10)
 call mpas_pool_get_array(diag_physics,'ust',ust)
 call mpas_pool_get_array(diag_physics,'znt',znt)

 this%dt_c   = dt
 this%grav_c = gravity

 do i = its,ite
    this%isltyp_c(i)   = isltyp(i)
    this%area_c(i)     = area(i)
    this%xland_c(i)    = xland(i)
    this%u10_c(i)      = u10(i)
    this%v10_c(i)      = v10(i)
    this%ust_c(i)      = ust(i)
    this%znt_c(i)      = znt(i)
    this%porosity_c(i) = porosity(i)
    this%smois_c(i)    = smois(kts,i)
    this%snowh_c(i)    = snowh(i)

    this%clay_c(i)   = clay(i)
    this%sand_c(i)   = sand(i)
    do n = 1,nerod
       this%erod_c(i,n) = erod(n,i)
    enddo

    do k = kts,kte
       this%qco(i,k)     = qco(k,i)
       this%zgrid_c(i,k) = zgrid(k,i)
       this%zgrid_c(i,k) = max(0._RKIND,this%zgrid_c(i,k))
       this%rho_c(i,k)   = zz(k,i)*rho_zz(k,i)*(1.+qv(k,i))
       this%pres_c(i,k)  = pres_b(k,i)+pres_p(k,i)
       this%u_c(i,k)     = uzonal(k,i)
       this%v_c(i,k)     = umerid(k,i)
    enddo
    k = kte+1
    this%zgrid_c(i,k) = zgrid(k,i)
 enddo

!--- interpolation of pressure thickness:
 pres  => this%pres_c
 dpres => this%dpres_c

 if(.not.allocated(pres2)) allocate(pres2(its:ite,kts:kte+1))
 do i = its,ite
    !pressure at the surface:
    k = kts
    z0  = zgrid(k,i)
    z1  = 0.5*(zgrid(k,i)+zgrid(k+1,i))
    z2  = 0.5*(zgrid(k+1,i)+zgrid(k+2,i))
    fzp = (z0-z2)/(z1-z2)
    fzm = 1.-fzp
    pres2(i,k) = fzp*pres(i,k) + fzm*pres(i,k+1)
    !pressure at interface between layers above the surface (k=kts) and below the model-top (k=kte+1):
    do k = kts+1,kte
       tem = 1./(zgrid(k+1,i)-zgrid(k-1,i))
       fzm = (zgrid(k,i)-zgrid(k-1,i))*tem
       fzp = (zgrid(k+1,i)-zgrid(k,i))*tem
       pres2(i,k) = fzm*pres(i,k) + fzp*pres(i,k-1)
    enddo
    !pressure at the model-top (k=kte+1):
    k = kte+1
    z0  = zgrid(k,i)
    z1  = 0.5*(zgrid(k,i)+zgrid(k-1,i))
    z2  = 0.5*(zgrid(k-1,i)+zgrid(k-2,i))
    fzm = (z0-z2)/(z1-z2)
    fzp = 1.-fzm
    !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
    pres2(i,k) = exp(fzm*log(pres(i,k-1)) + fzp*log(pres(i,k-2)))

    !pressure thickness:
    do k = kts,kte
       dpres(i,k) = pres2(i,k)-pres2(i,k+1)
    enddo
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r',intArgs=(/i,k/),realArgs=(/pres2(i,k),pres(i,k),dpres(i,k)/))
!   enddo
!   k = kte+1
!   call mpas_log_write('$i $i $r',intArgs=(/i,k/),realArgs=(/pres2(i,k)/))
!   call mpas_log_write(' ')
 enddo
 if(allocated(pres2)) deallocate(pres2)


 call mpas_log_write('--- end subroutine mpas_chemistry_from_gocart:')

 end subroutine mpas_chemistry_from_gocart

!=================================================================================================================
 subroutine mpas_chemistry_to_gocart(this,chems_input,state,time_lev,diag_chems,diag_physics)
!=================================================================================================================

!input arguments:
 integer,intent(in):: time_lev
 class(chems_gocart),intent(inout):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: chems_input
 type(mpas_pool_type),intent(inout):: diag_chems
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte

 real(kind=RKIND):: dt
 real(kind=RKIND),dimension(:),pointer:: nwfa2d_so4,nwfa2d_oc,nwfa2d_seas,nwfa2d
 real(kind=RKIND),dimension(:),pointer:: nifa2d_dust,nifa2d

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_to_gocart:')

 its = this%its
 ite = this%ite
 kts = this%kts
 kte = this%kte

 dt = this%dt_c

!--- finish to calculate the tendencies after call to the GOCART model. Because mixing ratios
!    are only modified in the layer adjacent to the surface, tendencies only need to be
!    updated in layer k = kts:
 k = kts
 do i = its,ite
    this%tend_qbc1(i,k)   = (this%qbc1(i,k)   - this%tend_qbc1(i,k)  ) / dt
    this%tend_qbc2(i,k)   = (this%qbc2(i,k)   - this%tend_qbc2(i,k)  ) / dt
    this%tend_qoc1(i,k)   = (this%qoc1(i,k)   - this%tend_qoc1(i,k)  ) / dt
    this%tend_qoc2(i,k)   = (this%qoc2(i,k)   - this%tend_qoc2(i,k)  ) / dt
    this%tend_qseas1(i,k) = (this%qseas1(i,k) - this%tend_qseas1(i,k)) / dt
    this%tend_qseas2(i,k) = (this%qseas2(i,k) - this%tend_qseas2(i,k)) / dt
    this%tend_qseas3(i,k) = (this%qseas3(i,k) - this%tend_qseas3(i,k)) / dt
    this%tend_qseas4(i,k) = (this%qseas4(i,k) - this%tend_qseas4(i,k)) / dt
    this%tend_qdust1(i,k) = (this%qdust1(i,k) - this%tend_qdust1(i,k)) / dt
    this%tend_qdust2(i,k) = (this%qdust2(i,k) - this%tend_qdust2(i,k)) / dt
    this%tend_qdust3(i,k) = (this%qdust3(i,k) - this%tend_qdust3(i,k)) / dt
    this%tend_qdust4(i,k) = (this%qdust4(i,k) - this%tend_qdust4(i,k)) / dt
    this%tend_qdust5(i,k) = (this%qdust5(i,k) - this%tend_qdust5(i,k)) / dt
    this%tend_qnh3(i,k)   = (this%qnh3(i,k)   - this%tend_qnh3(i,k)  ) / dt
    this%tend_qso2(i,k)   = (this%qso2(i,k)   - this%tend_qso2(i,k)  ) / dt
    this%tend_qso4a1(i,k) = (this%qso4a1(i,k) - this%tend_qso4a1(i,k)) / dt
    this%tend_qso4a2(i,k) = (this%qso4a2(i,k) - this%tend_qso4a2(i,k)) / dt
    this%tend_qso4a3(i,k) = (this%qso4a3(i,k) - this%tend_qso4a3(i,k)) / dt
 enddo
 call mpas_log_write(' ')


!--- update diagnostics of the surface emissions of nwfa and nifa:
!call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)
!call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)
!call mpas_pool_get_array(diag_chems,'nwfa2d_oc'  ,nwfa2d_oc)
!call mpas_pool_get_array(diag_chems,'nwfa2d_so4' ,nwfa2d_so4)
!call mpas_pool_get_array(diag_chems,'nwfa2d_seas',nwfa2d_seas)
!call mpas_pool_get_array(diag_chems,'nifa2d_dust',nifa2d_dust)

!k = kts
!do i = its,ite
!   nwfa2d(i) = 3*iconv_oc2*this%tend_qoc2(i,k)      &
!             + 0.4*iconv_so4*this%tend_qso4a1(i,k)  &
!             + 50*iconv_seas1*this%tend_qseas1(i,k) &
!             + 50*iconv_seas2*this%tend_qseas2(i,k)

!   nwfa2d_oc(i)   = 3*iconv_oc2*this%tend_qoc2(i,k)
!   nwfa2d_so4(i)  = 0.4*iconv_so4*this%tend_qso4a1(i,k)
!   nwfa2d_seas(i) = 50*iconv_seas1*this%tend_qseas1(i,k) &
!                  + 50*iconv_seas2*this%tend_qseas2(i,k)

!   nifa2d(i) = iconv_dust1*this%tend_qdust1(i,k) &
!             + iconv_dust2*this%tend_qdust2(i,k) &
!             + iconv_dust3*this%tend_qdust3(i,k) &
!             + iconv_dust4*this%tend_qdust4(i,k) &
!             + iconv_dust5*this%tend_qdust5(i,k)

!   nifa2d_dust(i) = nifa2d(i)
!enddo


 call compute_sfc_nbThompson(its,ite,kts,kte,this%tend_qbc2,this%tend_qoc2,this%tend_qso4a1, &
                        this%tend_qseas1,this%tend_qdust1,this%tend_qdust2,this%tend_qdust3, &
                        this%tend_qdust4,this%tend_qdust5,this%nwfa2d_c,this%nifa2d_c)

!---
 call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)
 call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)

 do i = its,ite
    nifa2d(i) = this%nifa2d_c(i)
    nwfa2d(i) = this%nwfa2d_c(i)
 enddo

 call mpas_log_write('--- end subroutine mpas_chemistry_to_gocart:')

 end subroutine mpas_chemistry_to_gocart

!=================================================================================================================
 subroutine compute_sfc_nbThompson(its,ite,kts,kte,tend_qbc2,tend_qoc2,tend_qso4a1,tend_qseas1,tend_qdust1, &
                                   tend_qdust2,tend_qdust3,tend_qdust4,tend_qdust5,tend_nwfa,tend_nifa)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte
 real(kind=RKIND),dimension(its:ite,kts:kte),intent(in):: tend_qbc2,tend_qoc2,tend_qso4a1,tend_qseas1
 real(kind=RKIND),dimension(its:ite,kts:kte),intent(in):: tend_qdust1,tend_qdust2,tend_qdust3,tend_qdust4, &
                                                          tend_qdust5

!inout arguments:
 real(kind=RKIND),dimension(its:ite),intent(inout):: tend_nwfa,tend_nifa

!local variables and arrays:
 integer:: i,k

!these parameters should be moved somewhere else:
 real(kind=RKIND),parameter:: pii = 3.141592653589793

 real(kind=RKIND),parameter:: dens_bc       = 1.7e03  ! density of black carbon (kg/m3).
 real(kind=RKIND),parameter:: dens_oc       = 1.0e03  ! density of organic carbon (kg/m3).
 real(kind=RKIND),parameter:: dens_so4      = 1.77e03 ! density of sulfates (kg/m3).
 real(kind=RKIND),parameter:: dens_ncl      = 1.9e03  ! density of seasalts (kg/m3).
 real(kind=RKIND),parameter:: dens_dst      = 2.6e03  ! density of dust (kg/m3).

 real(kind=RKIND),parameter:: dgnum_accum   = 1.1e-07 ! (m).
 real(kind=RKIND),parameter:: sigmag_accum  = 1.6     ! (-).
 real(kind=RKIND):: alnsg_accum,voltonumb_accum
 real(kind=RKIND):: sum_a1,ratio_a1

 real(kind=RKIND),parameter:: dgnum_aitken  = 2.6e-08 ! (m).
 real(kind=RKIND),parameter:: sigmag_aitken = 1.6     ! (-).
 real(kind=RKIND):: alnsg_aitken,voltonumb_aitken
 real(kind=RKIND):: sum_a2,ratio_a2

 real(kind=RKIND),parameter:: dgnum_coarse  = 9.0e-07 ! (m).
 real(kind=RKIND),parameter:: sigmag_coarse = 1.2     ! (-).
 real(kind=RKIND):: alnsg_coarse,voltonumb_coarse
 real(kind=RKIND):: sum_a3,ratio_a3
 real(kind=RKIND):: mr_to_nb_bc,mr_to_nb_oc,mr_to_nb_so4,mr_to_nb_ncl,mr_to_nb_dst

 real(kind=RKIND),dimension(its:ite):: tend_nbca1,tend_npoma1,tend_nso4a1,tend_nseasa1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine compute_sfc_nbThompson:')

!--- accumulation mode for nwfa:
 alnsg_accum = log(sigmag_accum)
 voltonumb_accum = 1./((pii/6.)*(dgnum_accum**3.)*exp(4.5*alnsg_accum**2.))

 mr_to_nb_bc  = voltonumb_accum/dens_bc
 mr_to_nb_oc  = voltonumb_accum/dens_oc
 mr_to_nb_so4 = voltonumb_accum/dens_so4
 mr_to_nb_ncl = voltonumb_accum/dens_ncl
 mr_to_nb_dst = voltonumb_accum/dens_dst

 call mpas_log_write('--- alnsg_accum      = $r',realArgs=(/alnsg_accum/))
 call mpas_log_write('--- voltonumb_accum  = $r',realArgs=(/voltonumb_accum/))
 call mpas_log_write('--- mr_to_nb_bc      = $r',realArgs=(/mr_to_nb_bc/))
 call mpas_log_write('--- mr_to_nb_dst     = $r',realArgs=(/mr_to_nb_dst/))
 call mpas_log_write('--- mr_to_nb_ncl     = $r',realArgs=(/mr_to_nb_ncl/))
 call mpas_log_write('--- mr_to_nb_oc      = $r',realArgs=(/mr_to_nb_oc/))
 call mpas_log_write('--- mr_to_nb_so4     = $r',realArgs=(/mr_to_nb_so4/))
 call mpas_log_write(' ')

 k = kts
 do i = its,ite
    tend_nbca1(i)   = mr_to_nb_bc*tend_qbc2(i,k)
    tend_npoma1(i)  = mr_to_nb_oc*tend_qoc2(i,k)
    tend_nso4a1(i)  = mr_to_nb_so4*tend_qso4a1(i,k)
    tend_nseasa1(i) = mr_to_nb_so4*tend_qseas1(i,k)

    tend_nwfa(i)   = tend_nbca1(i) + tend_npoma1(i) + tend_nso4a1(i) + tend_nseasa1(i)
!   call mpas_log_write('$i $r $r $r $r $r $r $r $r $r',intArgs=(/i/),realArgs=(/tend_qbc2(i,k),tend_nbca1(i), &
!              tend_qoc2(i,k),tend_npoma1(i),tend_qso4a1(i,k),tend_nso4a1(i),tend_qseas1(i,k),tend_nseasa1(i), &
!              tend_nwfa(i)/))
 enddo


!--- coarse mode for nifa:
 alnsg_coarse = log(sigmag_coarse)
 voltonumb_coarse = 1./((pii/6.)*(dgnum_coarse**3.)*exp(4.5*alnsg_coarse**2.))

 mr_to_nb_ncl = voltonumb_coarse/dens_ncl
 mr_to_nb_dst = voltonumb_coarse/dens_dst
 mr_to_nb_so4 = voltonumb_coarse/dens_so4

 call mpas_log_write(' ')
 call mpas_log_write('--- alnsg_coarse     = $r',realArgs=(/alnsg_coarse/))
 call mpas_log_write('--- voltonumb_coarse = $r',realArgs=(/voltonumb_coarse/))
 call mpas_log_write('--- mr_to_nb_so4     = $r',realArgs=(/mr_to_nb_so4/))
 call mpas_log_write('--- mr_to_nb_ncl     = $r',realArgs=(/mr_to_nb_ncl/))
 call mpas_log_write('--- mr_to_nb_dst     = $r',realArgs=(/mr_to_nb_dst/))
 call mpas_log_write(' ')

 k = kts
 do i = its,ite
    tend_nifa(i) = tend_qdust1(i,k) + tend_qdust2(i,k) + tend_qdust3(i,k) &
                 + tend_qdust4(i,k) + tend_qdust5(i,k)
    tend_nifa(i) = mr_to_nb_dst*tend_nifa(i)
!   call mpas_log_write('$i $r $r $r $r $r $r',intArgs=(/i/),realArgs=(/tend_qdust1(i,k),tend_qdust2(i,k), &
!             tend_qdust3(i,k),tend_qdust4(i,k),tend_qdust5(i,k),tend_nifa(i)/))
 enddo


 call mpas_log_write('--- end subroutine compute_sfc_nbThompson:')

 end subroutine compute_sfc_nbThompson

!=================================================================================================================
 end module mpas_chemistry_gocart_interface
!=================================================================================================================

