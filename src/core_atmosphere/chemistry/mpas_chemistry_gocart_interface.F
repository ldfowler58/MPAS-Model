!=================================================================================================================
 module mpas_chemistry_gocart_interface
 use mpas_log
 use mpas_kind_types
 use mpas_constants,only: gravity
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension

 use module_data_rrtmg_aeropt,only: nlwbands,nswbands


 implicit none
 private
 public:: mpas_chemistry_gocart_dims,       &
          mpas_chemistry_gocart_allocate,   &
          mpas_chemistry_gocart_deallocate, &
          mpas_chemistry_from_gocart,       &
          mpas_chemistry_to_gocart


 type,public:: chems_gocart
    integer:: its,ite,kts,kte,ktep1
    integer:: nerod

    real(kind=RKIND),dimension(:),pointer:: qbc1_em       => null()
    real(kind=RKIND),dimension(:),pointer:: qoc1_em       => null()
    real(kind=RKIND),dimension(:),pointer:: qco_em        => null()
    real(kind=RKIND),dimension(:),pointer:: qnh3_em       => null()
    real(kind=RKIND),dimension(:),pointer:: qso2_em       => null()

    real(kind=RKIND),dimension(:),pointer  :: backg_dms   => null()
    real(kind=RKIND),dimension(:,:),pointer:: backg_oh    => null()
    real(kind=RKIND),dimension(:,:),pointer:: backg_h2o2  => null()
    real(kind=RKIND),dimension(:,:),pointer:: backg_no3   => null()

    real(kind=RKIND),dimension(:),pointer:: emis_qseas1   => null() !GOCART seasalt emission (0.1-0.5 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qseas2   => null() !GOCART seasalt emission (0.5-1.5 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qseas3   => null() !GOCART seasalt emission (1.5-5.0 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qseas4   => null() !GOCART seasalt emission (5.0-10. microns).

    real(kind=RKIND),dimension(:),pointer:: emis_qdust1   => null() !GOCART dust emission (0.1-1.0 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qdust2   => null() !GOCART dust emission (1.0-1.8 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qdust3   => null() !GOCART dust emission (1.8-3.0 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qdust4   => null() !GOCART dust emission (3.0-6.0 microns).
    real(kind=RKIND),dimension(:),pointer:: emis_qdust5   => null() !GOCART dust emission (6.0-10. microns).

    real(kind=RKIND),dimension(:,:),pointer:: qbc1        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qbc2        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc1        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc2        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust3      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust4      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust5      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas1      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas2      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas3      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas4      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qsulf       => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdms        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qmsa        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qco         => null()
    real(kind=RKIND),dimension(:,:),pointer:: qh2o2       => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoh         => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso2        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qno3        => null()

    real(kind=RKIND),dimension(:,:),pointer:: tend_qbc1   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qbc2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qoc1   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qoc2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust1 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust2 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust3 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust4 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdust5 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas1 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas2 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas3 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qseas4 => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qsulf  => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qdms   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qmsa   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qco    => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qh2o2  => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qoh    => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qso2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: tend_qno3   => null()


    !--- diagnostics:
    real(kind=RKIND),dimension(:,:),pointer:: qpm10       => null()
    real(kind=RKIND),dimension(:,:),pointer:: qpm25       => null()


    !--- aerosol optical properties:
    real(kind=RKIND),dimension(:,:,:),pointer:: tauaerlw  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: extaerlw  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: extaersw  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: gaersw    => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: waersw    => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: bscoefsw  => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: sw12aer   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: sw13aer   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: sw14aer   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: sw15aer   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: sw16aer   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: sw17aer   => null()


    !--- local variables from the MPAS dynamical core:
    integer,dimension(:),pointer           :: isltyp_c    => null()

    real(kind=RKIND):: dt_c
    real(kind=RKIND):: grav_c
    real(kind=RKIND),dimension(:),pointer  :: clay_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: sand_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: porosity_c  => null()
    real(kind=RKIND),dimension(:),pointer  :: xlon_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: xlat_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: area_c      => null()
    real(kind=RKIND),dimension(:),pointer  :: xland_c     => null()
    real(kind=RKIND),dimension(:),pointer  :: u10_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: v10_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: smois_c     => null()
    real(kind=RKIND),dimension(:),pointer  :: snowh_c     => null()
    real(kind=RKIND),dimension(:),pointer  :: ust_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: znt_c       => null()
    real(kind=RKIND),dimension(:),pointer  :: tsk_c       => null()
    real(kind=RKIND),dimension(:,:),pointer:: erod_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: zgrid_c     => null()
    real(kind=RKIND),dimension(:,:),pointer:: rho_c       => null()
    real(kind=RKIND),dimension(:,:),pointer:: pres_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: dpres_c     => null()
    real(kind=RKIND),dimension(:,:),pointer:: u_c         => null()
    real(kind=RKIND),dimension(:,:),pointer:: v_c         => null()
    real(kind=RKIND),dimension(:,:),pointer:: tt_c        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qc_c        => null()
    real(kind=RKIND),dimension(:,:),pointer:: qi_c        => null()
    real(kind=RKIND),dimension(:,:),pointer:: cld_c       => null()
    real(kind=RKIND),dimension(:,:),pointer:: rh_c        => null()


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


!--- allocate background fields:
 if(.not.associated(this%backg_dms)  ) allocate(this%backg_dms(its:ite)          )
 if(.not.associated(this%backg_oh)   ) allocate(this%backg_oh(its:ite,kts:kte)   )
 if(.not.associated(this%backg_h2o2) ) allocate(this%backg_h2o2(its:ite,kts:kte) )
 if(.not.associated(this%backg_no3)  ) allocate(this%backg_no3(its:ite,kts:kte)  )


!--- allocate mixing ratios:
 if(.not.associated(this%qbc1)       ) allocate(this%qbc1(its:ite,kts:kte)       )
 if(.not.associated(this%qbc2)       ) allocate(this%qbc2(its:ite,kts:kte)       )
 if(.not.associated(this%qoc1)       ) allocate(this%qoc1(its:ite,kts:kte)       )
 if(.not.associated(this%qoc2)       ) allocate(this%qoc2(its:ite,kts:kte)       )
 if(.not.associated(this%qdust1)     ) allocate(this%qdust1(its:ite,kts:kte)     )
 if(.not.associated(this%qdust2)     ) allocate(this%qdust2(its:ite,kts:kte)     )
 if(.not.associated(this%qdust3)     ) allocate(this%qdust3(its:ite,kts:kte)     )
 if(.not.associated(this%qdust4)     ) allocate(this%qdust4(its:ite,kts:kte)     )
 if(.not.associated(this%qdust5)     ) allocate(this%qdust5(its:ite,kts:kte)     )
 if(.not.associated(this%qseas1)     ) allocate(this%qseas1(its:ite,kts:kte)     )
 if(.not.associated(this%qseas2)     ) allocate(this%qseas2(its:ite,kts:kte)     )
 if(.not.associated(this%qseas3)     ) allocate(this%qseas3(its:ite,kts:kte)     )
 if(.not.associated(this%qseas4)     ) allocate(this%qseas4(its:ite,kts:kte)     )
 if(.not.associated(this%qsulf)      ) allocate(this%qsulf(its:ite,kts:kte)      )
 if(.not.associated(this%qdms)       ) allocate(this%qdms(its:ite,kts:kte)       )
 if(.not.associated(this%qmsa)       ) allocate(this%qmsa(its:ite,kts:kte)       )
 if(.not.associated(this%qco)        ) allocate(this%qco(its:ite,kts:kte)        )
 if(.not.associated(this%qh2o2)      ) allocate(this%qh2o2(its:ite,kts:kte)      )
 if(.not.associated(this%qoh)        ) allocate(this%qoh(its:ite,kts:kte)        )
 if(.not.associated(this%qso2)       ) allocate(this%qso2(its:ite,kts:kte)       )
 if(.not.associated(this%qno3)       ) allocate(this%qno3(its:ite,kts:kte)       )


!--- allocate diagnostics:
 if(.not.associated(this%emis_qseas1)) allocate(this%emis_qseas1(its:ite)        )
 if(.not.associated(this%emis_qseas2)) allocate(this%emis_qseas2(its:ite)        )
 if(.not.associated(this%emis_qseas3)) allocate(this%emis_qseas3(its:ite)        )
 if(.not.associated(this%emis_qseas4)) allocate(this%emis_qseas4(its:ite)        )

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
 if(.not.associated(this%tend_qdust1)) allocate(this%tend_qdust1(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust2)) allocate(this%tend_qdust2(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust3)) allocate(this%tend_qdust3(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust4)) allocate(this%tend_qdust4(its:ite,kts:kte))
 if(.not.associated(this%tend_qdust5)) allocate(this%tend_qdust5(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas1)) allocate(this%tend_qseas1(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas2)) allocate(this%tend_qseas2(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas3)) allocate(this%tend_qseas3(its:ite,kts:kte))
 if(.not.associated(this%tend_qseas4)) allocate(this%tend_qseas4(its:ite,kts:kte))
 if(.not.associated(this%tend_qsulf) ) allocate(this%tend_qsulf(its:ite,kts:kte) )
 if(.not.associated(this%tend_qdms)  ) allocate(this%tend_qdms(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qmsa)  ) allocate(this%tend_qmsa(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qco)   ) allocate(this%tend_qco(its:ite,kts:kte)   )
 if(.not.associated(this%tend_qh2o2) ) allocate(this%tend_qh2o2(its:ite,kts:kte) )
 if(.not.associated(this%tend_qoh)   ) allocate(this%tend_qoh(its:ite,kts:kte)   )
 if(.not.associated(this%tend_qso2)  ) allocate(this%tend_qso2(its:ite,kts:kte)  )
 if(.not.associated(this%tend_qno3)  ) allocate(this%tend_qno3(its:ite,kts:kte)  )


!--- allocate diagnostics:
 if(.not.associated(this%qpm10)      ) allocate(this%qpm10(its:ite,kts:kte)      )
 if(.not.associated(this%qpm25)      ) allocate(this%qpm25(its:ite,kts:kte)      )


!--- allocate aerosol optical properties:
 if(.not.associated(this%tauaerlw)   ) allocate(this%tauaerlw(its:ite,kts:kte,nlwbands))
 if(.not.associated(this%extaerlw)   ) allocate(this%extaerlw(its:ite,kts:kte,nlwbands))
 if(.not.associated(this%tauaersw)   ) allocate(this%tauaersw(its:ite,kts:kte,nswbands))
 if(.not.associated(this%extaersw)   ) allocate(this%extaersw(its:ite,kts:kte,nswbands))
 if(.not.associated(this%gaersw  )   ) allocate(this%gaersw(its:ite,kts:kte,nswbands)  )
 if(.not.associated(this%waersw  )   ) allocate(this%waersw(its:ite,kts:kte,nswbands)  )
 if(.not.associated(this%bscoefsw)   ) allocate(this%bscoefsw(its:ite,kts:kte,nswbands))
 if(.not.associated(this%sw12aer )   ) allocate(this%sw12aer(its:ite,kts:kte,nswbands) )
 if(.not.associated(this%sw13aer )   ) allocate(this%sw13aer(its:ite,kts:kte,nswbands) )
 if(.not.associated(this%sw14aer )   ) allocate(this%sw14aer(its:ite,kts:kte,nswbands) )
 if(.not.associated(this%sw15aer )   ) allocate(this%sw15aer(its:ite,kts:kte,nswbands) )
 if(.not.associated(this%sw16aer )   ) allocate(this%sw16aer(its:ite,kts:kte,nswbands) )
 if(.not.associated(this%sw17aer )   ) allocate(this%sw17aer(its:ite,kts:kte,nswbands) )


!--- allocate mesh and meteorological fields:
 if(.not.associated(this%isltyp_c)   ) allocate(this%isltyp_c(its:ite)           )
 if(.not.associated(this%clay_c)     ) allocate(this%clay_c(its:ite)             )
 if(.not.associated(this%sand_c)     ) allocate(this%sand_c(its:ite)             )
 if(.not.associated(this%porosity_c) ) allocate(this%porosity_c(its:ite)         )
 if(.not.associated(this%xlon_c)     ) allocate(this%xlon_c(its:ite)             )
 if(.not.associated(this%xlat_c)     ) allocate(this%xlat_c(its:ite)             )
 if(.not.associated(this%area_c)     ) allocate(this%area_c(its:ite)             )
 if(.not.associated(this%xland_c)    ) allocate(this%xland_c(its:ite)            )
 if(.not.associated(this%u10_c)      ) allocate(this%u10_c(its:ite)              )
 if(.not.associated(this%v10_c)      ) allocate(this%v10_c(its:ite)              )
 if(.not.associated(this%smois_c)    ) allocate(this%smois_c(its:ite)            )
 if(.not.associated(this%snowh_c)    ) allocate(this%snowh_c(its:ite)            )
 if(.not.associated(this%ust_c)      ) allocate(this%ust_c(its:ite)              )
 if(.not.associated(this%znt_c)      ) allocate(this%znt_c(its:ite)              )
 if(.not.associated(this%tsk_c)      ) allocate(this%tsk_c(its:ite)              )
 if(.not.associated(this%zgrid_c)    ) allocate(this%zgrid_c(its:ite,kts:ktep1)  )
 if(.not.associated(this%rho_c)      ) allocate(this%rho_c(its:ite,kts:kte)      )
 if(.not.associated(this%pres_c)     ) allocate(this%pres_c(its:ite,kts:kte)     )
 if(.not.associated(this%dpres_c)    ) allocate(this%dpres_c(its:ite,kts:kte)    )
 if(.not.associated(this%u_c)        ) allocate(this%u_c(its:ite,kts:kte)        )
 if(.not.associated(this%v_c)        ) allocate(this%v_c(its:ite,kts:kte)        )
 if(.not.associated(this%tt_c)       ) allocate(this%tt_c(its:ite,kts:kte)       )
 if(.not.associated(this%qc_c)       ) allocate(this%qc_c(its:ite,kts:kte)       )
 if(.not.associated(this%qi_c)       ) allocate(this%qi_c(its:ite,kts:kte)       )
 if(.not.associated(this%cld_c)      ) allocate(this%cld_c(its:ite,kts:kte)      )
 if(.not.associated(this%rh_c)       ) allocate(this%rh_c(its:ite,kts:kte)       )

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


!--- allocate background fields:
 if(associated(this%backg_dms)  ) deallocate(this%backg_dms  )
 if(associated(this%backg_oh)   ) deallocate(this%backg_oh   )
 if(associated(this%backg_h2o2) ) deallocate(this%backg_h2o2 )
 if(associated(this%backg_no3)  ) deallocate(this%backg_no3  )


!--- deallocate mixing ratios:
 if(associated(this%qbc1)       ) deallocate(this%qbc1       )
 if(associated(this%qbc2)       ) deallocate(this%qbc2       )
 if(associated(this%qoc1)       ) deallocate(this%qoc1       )
 if(associated(this%qoc2)       ) deallocate(this%qoc2       )
 if(associated(this%qdust1)     ) deallocate(this%qdust1     )
 if(associated(this%qdust2)     ) deallocate(this%qdust2     )
 if(associated(this%qdust3)     ) deallocate(this%qdust3     )
 if(associated(this%qdust4)     ) deallocate(this%qdust4     )
 if(associated(this%qdust5)     ) deallocate(this%qdust5     )
 if(associated(this%qseas1)     ) deallocate(this%qseas1     )
 if(associated(this%qseas2)     ) deallocate(this%qseas2     )
 if(associated(this%qseas3)     ) deallocate(this%qseas3     )
 if(associated(this%qseas4)     ) deallocate(this%qseas4     )
 if(associated(this%qsulf)      ) deallocate(this%qsulf      )
 if(associated(this%qdms)       ) deallocate(this%qdms       )
 if(associated(this%qmsa)       ) deallocate(this%qmsa       )
 if(associated(this%qco)        ) deallocate(this%qco        )
 if(associated(this%qh2o2)      ) deallocate(this%qh2o2      )
 if(associated(this%qoh)        ) deallocate(this%qoh        )
 if(associated(this%qso2)       ) deallocate(this%qso2       )
 if(associated(this%qno3)       ) deallocate(this%qno3       )


!--- deallocate diagnostics:
 if(associated(this%emis_qseas1)) deallocate(this%emis_qseas1)
 if(associated(this%emis_qseas2)) deallocate(this%emis_qseas2)
 if(associated(this%emis_qseas3)) deallocate(this%emis_qseas3)
 if(associated(this%emis_qseas4)) deallocate(this%emis_qseas4)

 if(associated(this%emis_qdust1)) deallocate(this%emis_qdust1)
 if(associated(this%emis_qdust2)) deallocate(this%emis_qdust2)
 if(associated(this%emis_qdust3)) deallocate(this%emis_qdust3)
 if(associated(this%emis_qdust4)) deallocate(this%emis_qdust4)
 if(associated(this%emis_qdust5)) deallocate(this%emis_qdust5)


!--- deallocate tendencies:
 if(associated(this%tend_qbc1)  ) deallocate(this%tend_qbc1  )
 if(associated(this%tend_qbc2)  ) deallocate(this%tend_qbc2  )
 if(associated(this%tend_qoc1)  ) deallocate(this%tend_qoc1  )
 if(associated(this%tend_qoc2)  ) deallocate(this%tend_qoc2  )
 if(associated(this%tend_qdust1)) deallocate(this%tend_qdust1)
 if(associated(this%tend_qdust2)) deallocate(this%tend_qdust2)
 if(associated(this%tend_qdust3)) deallocate(this%tend_qdust3)
 if(associated(this%tend_qdust4)) deallocate(this%tend_qdust4)
 if(associated(this%tend_qdust5)) deallocate(this%tend_qdust5)
 if(associated(this%tend_qseas1)) deallocate(this%tend_qseas1)
 if(associated(this%tend_qseas2)) deallocate(this%tend_qseas2)
 if(associated(this%tend_qseas3)) deallocate(this%tend_qseas3)
 if(associated(this%tend_qseas4)) deallocate(this%tend_qseas4)
 if(associated(this%tend_qsulf) ) deallocate(this%tend_qsulf )
 if(associated(this%tend_qdms)  ) deallocate(this%tend_qdms  )
 if(associated(this%tend_qmsa)  ) deallocate(this%tend_qmsa  )
 if(associated(this%tend_qco)   ) deallocate(this%tend_qco   )
 if(associated(this%tend_qh2o2) ) deallocate(this%tend_qh2o2 )
 if(associated(this%tend_qoh)   ) deallocate(this%tend_qoh   )
 if(associated(this%tend_qso2)  ) deallocate(this%tend_qso2  )
 if(associated(this%tend_qno3)  ) deallocate(this%tend_qno3  )


!--- deallocate diagnostics:
 if(associated(this%qpm10)      ) deallocate(this%qpm10      )
 if(associated(this%qpm25)      ) deallocate(this%qpm25      )


!--- deallocate aerosol optical properties:
 if(associated(this%tauaerlw)   ) deallocate(this%tauaerlw   )
 if(associated(this%extaerlw)   ) deallocate(this%extaerlw   )
 if(associated(this%tauaersw)   ) deallocate(this%tauaersw   )
 if(associated(this%extaersw)   ) deallocate(this%extaersw   )
 if(associated(this%gaersw  )   ) deallocate(this%gaersw     )
 if(associated(this%waersw  )   ) deallocate(this%waersw     )
 if(associated(this%bscoefsw)   ) deallocate(this%bscoefsw   )
 if(associated(this%sw12aer )   ) deallocate(this%sw12aer    )
 if(associated(this%sw13aer )   ) deallocate(this%sw13aer    )
 if(associated(this%sw14aer )   ) deallocate(this%sw14aer    )
 if(associated(this%sw15aer )   ) deallocate(this%sw15aer    )
 if(associated(this%sw16aer )   ) deallocate(this%sw16aer    )
 if(associated(this%sw17aer )   ) deallocate(this%sw17aer    )


!--- deallocate mesh and meteorological fields:
 if(associated(this%isltyp_c)   ) deallocate(this%isltyp_c   )
 if(associated(this%clay_c)     ) deallocate(this%clay_c     )
 if(associated(this%sand_c)     ) deallocate(this%sand_c     )
 if(associated(this%porosity_c) ) deallocate(this%porosity_c )
 if(associated(this%xlon_c)     ) deallocate(this%xlon_c     )
 if(associated(this%xlat_c)     ) deallocate(this%xlat_c     )
 if(associated(this%area_c)     ) deallocate(this%area_c     )
 if(associated(this%xland_c)    ) deallocate(this%xland_c    )
 if(associated(this%u10_c)      ) deallocate(this%u10_c      )
 if(associated(this%v10_c)      ) deallocate(this%v10_c      )
 if(associated(this%smois_c)    ) deallocate(this%smois_c    )
 if(associated(this%snowh_c)    ) deallocate(this%snowh_c    )
 if(associated(this%ust_c)      ) deallocate(this%ust_c      )
 if(associated(this%znt_c)      ) deallocate(this%znt_c      )
 if(associated(this%tsk_c)      ) deallocate(this%tsk_c      )
 if(associated(this%zgrid_c)    ) deallocate(this%zgrid_c    )
 if(associated(this%rho_c)      ) deallocate(this%rho_c      )
 if(associated(this%pres_c)     ) deallocate(this%pres_c     )
 if(associated(this%dpres_c)    ) deallocate(this%dpres_c    )
 if(associated(this%u_c)        ) deallocate(this%u_c        )
 if(associated(this%v_c)        ) deallocate(this%v_c        )
 if(associated(this%tt_c)       ) deallocate(this%tt_c       )
 if(associated(this%qc_c)       ) deallocate(this%qc_c       )
 if(associated(this%qi_c)       ) deallocate(this%qi_c       )
 if(associated(this%cld_c)      ) deallocate(this%cld_c      )
 if(associated(this%rh_c)       ) deallocate(this%rh_c       )

 if(associated(this%erod_c)     ) deallocate(this%erod_c     )

 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_deallocate:')

 end subroutine mpas_chemistry_gocart_deallocate

!=================================================================================================================
 subroutine mpas_chemistry_from_gocart(this,CAMS_emissions,GOCART_backgrounds,configs,mesh,state,time_lev,diag, &
                                       diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: CAMS_emissions
 type(mpas_pool_type),intent(in):: GOCART_backgrounds

 integer,intent(in):: time_lev

!inout arguments:
 class(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,n,nerod

 integer,pointer:: num_chems
 integer,pointer:: index_qbc1,index_qbc2
 integer,pointer:: index_qoc1,index_qoc2
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4
 integer,pointer:: index_qsulf
 integer,pointer:: index_qdms,index_qmsa,index_qco,index_qh2o2,index_qoh,index_qso2,index_qno3

 real(kind=RKIND),dimension(:),pointer:: bc1_em_anthro,oc1_em_anthro,co_em_anthro, &
                                         nh3_em_anthro,so2_em_anthro
 real(kind=RKIND),dimension(:),pointer:: background_dms
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2
 real(kind=RKIND),dimension(:,:),pointer:: qoc1,qoc2
 real(kind=RKIND),dimension(:,:),pointer:: qseas1,qseas2,qseas3,qseas4
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2,qdust3,qdust4,qdust5
 real(kind=RKIND),dimension(:,:),pointer:: qsulf
 real(kind=RKIND),dimension(:,:),pointer:: qdms,qmsa,qco,qh2o2,qoh,qso2,qno3
 real(kind=RKIND),dimension(:,:),pointer:: background_h2o2,background_oh,background_no3
 real(kind=RKIND),dimension(:,:,:),pointer:: chems

!local variables and arrays for meteorological fields:
 integer,pointer:: index_qv,index_qc,index_qi
 integer,dimension(:),pointer:: isltyp

 real(kind=RKIND),parameter:: R_d = 287._RKIND
 real(kind=RKIND),parameter:: R_v = 461.6_RKIND

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:),pointer:: lat,lon,area,clay,sand,xland
 real(kind=RKIND),dimension(:),pointer:: u10,v10,skintemp
 real(kind=RKIND),dimension(:),pointer:: porosity,snowh,ust,znt
 real(kind=RKIND),dimension(:,:),pointer:: erod
 real(kind=RKIND),dimension(:,:),pointer:: zgrid,zz
 real(kind=RKIND),dimension(:,:),pointer:: cld,exner,pres_p,pres_b,qi,qc,qv,relhum,rho_zz,theta,uzonal,umerid
 real(kind=RKIND),dimension(:,:),pointer:: smois
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

 call mpas_pool_get_array(GOCART_backgrounds,'background_dms' ,background_dms )
 call mpas_pool_get_array(GOCART_backgrounds,'background_h2o2',background_h2o2)
 call mpas_pool_get_array(GOCART_backgrounds,'background_oh'  ,background_oh  )
 call mpas_pool_get_array(GOCART_backgrounds,'background_no3' ,background_no3 )

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
 call mpas_pool_get_dimension(state,'index_qco'   ,index_qco   )
 call mpas_pool_get_dimension(state,'index_qh2o2' ,index_qh2o2 )
 call mpas_pool_get_dimension(state,'index_qoh'   ,index_qoh   )
 call mpas_pool_get_dimension(state,'index_qso2'  ,index_qso2  )
 call mpas_pool_get_dimension(state,'index_qno3'  ,index_qno3  )

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
 call mpas_log_write('--- index_qdms   = $i',intArgs=(/index_qdms/)  )
 call mpas_log_write('--- index_qmsa   = $i',intArgs=(/index_qmsa/)  )
 call mpas_log_write('--- index_qco    = $i',intArgs=(/index_qco/)   )
 call mpas_log_write('--- index_qh2o2  = $i',intArgs=(/index_qh2o2/) )
 call mpas_log_write('--- index_qoh    = $i',intArgs=(/index_qoh/)   )
 call mpas_log_write('--- index_qso2   = $i',intArgs=(/index_qso2/)  )
 call mpas_log_write('--- index_qno3   = $i',intArgs=(/index_qno3/)  )

 call mpas_pool_get_array(state,'chems',chems,time_lev)
 qbc1   => chems(index_qbc1,:,:)
 qbc2   => chems(index_qbc2,:,:)
 qoc1   => chems(index_qoc1,:,:)
 qoc2   => chems(index_qoc2,:,:)
 qdust1 => chems(index_qdust1,:,:)
 qdust2 => chems(index_qdust2,:,:)
 qdust3 => chems(index_qdust3,:,:)
 qdust4 => chems(index_qdust4,:,:)
 qdust5 => chems(index_qdust5,:,:)
 qseas1 => chems(index_qseas1,:,:)
 qseas2 => chems(index_qseas2,:,:)
 qseas3 => chems(index_qseas3,:,:)
 qseas4 => chems(index_qseas4,:,:)
 qsulf  => chems(index_qsulf,:,:)
 qdms   => chems(index_qdms,:,:)
 qmsa   => chems(index_qmsa,:,:)
 qco    => chems(index_qco,:,:)
 qh2o2  => chems(index_qh2o2,:,:)
 qoh    => chems(index_qoh,:,:)
 qso2   => chems(index_qso2,:,:)
 qno3   => chems(index_qno3,:,:)

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

    this%backg_dms(i) = background_dms(i)

    do k = kts,kte
       this%qbc1(i,k)   = qbc1(k,i)
       this%qbc2(i,k)   = qbc2(k,i)
       this%qoc1(i,k)   = qoc1(k,i)
       this%qoc2(i,k)   = qoc2(k,i)
       this%qdust1(i,k) = qdust1(k,i)
       this%qdust2(i,k) = qdust2(k,i)
       this%qdust3(i,k) = qdust3(k,i)
       this%qdust4(i,k) = qdust4(k,i)
       this%qdust5(i,k) = qdust5(k,i)
       this%qseas1(i,k) = qseas1(k,i)
       this%qseas2(i,k) = qseas2(k,i)
       this%qseas3(i,k) = qseas3(k,i)
       this%qseas4(i,k) = qseas4(k,i)
       this%qsulf(i,k)  = qsulf(k,i)
       this%qdms(i,k)   = qdms(k,i)
       this%qmsa(i,k)   = qmsa(k,i)
       this%qco(i,k)    = qco(k,i)
       this%qh2o2(i,k)  = qh2o2(k,i)
       this%qoh(i,k)    = qoh(k,i)
       this%qso2(i,k)   = qso2(k,i)
       this%qno3(i,k)   = qno3(k,i)

       this%backg_h2o2(i,k) = background_h2o2(k,i)
       this%backg_oh(i,k)   = background_oh(k,i)
       this%backg_no3(i,k)  = background_no3(k,i)
    enddo
 enddo


!--- initialization of diagnostics:
 do i = its,ite
    this%emis_qseas1(i) = 0._RKIND
    this%emis_qseas2(i) = 0._RKIND
    this%emis_qseas3(i) = 0._RKIND
    this%emis_qseas4(i) = 0._RKIND

    this%emis_qdust1(i) = 0._RKIND
    this%emis_qdust2(i) = 0._RKIND
    this%emis_qdust3(i) = 0._RKIND
    this%emis_qdust4(i) = 0._RKIND
    this%emis_qdust5(i) = 0._RKIND
 enddo


!--- initialization of tendencies before call to the GOCART model:
 do i = its,ite
    do k = kts,kte
       this%tend_qbc1(i,k)   = this%qbc1(i,k)
       this%tend_qbc2(i,k)   = this%qbc2(i,k)
       this%tend_qoc1(i,k)   = this%qoc1(i,k)
       this%tend_qoc2(i,k)   = this%qoc2(i,k)
       this%tend_qdust1(i,k) = this%qdust1(i,k)
       this%tend_qdust2(i,k) = this%qdust2(i,k)
       this%tend_qdust3(i,k) = this%qdust3(i,k)
       this%tend_qdust4(i,k) = this%qdust4(i,k)
       this%tend_qdust5(i,k) = this%qdust5(i,k)
       this%tend_qseas1(i,k) = this%qseas1(i,k)
       this%tend_qseas2(i,k) = this%qseas2(i,k)
       this%tend_qseas3(i,k) = this%qseas3(i,k)
       this%tend_qseas4(i,k) = this%qseas4(i,k)
       this%tend_qsulf(i,k)  = this%qsulf(i,k)
       this%tend_qdms(i,k)   = this%qdms(i,k)
       this%tend_qmsa(i,k)   = this%qmsa(i,k)
       this%tend_qco(i,k)    = this%qco(i,k)
       this%tend_qh2o2(i,k)  = this%qh2o2(i,k)
       this%tend_qoh(i,k)    = this%qoh(i,k)
       this%tend_qso2(i,k)   = this%qso2(i,k)
       this%tend_qno3(i,k)   = this%qno3(i,k)
    enddo
 enddo


!do i = its,ite
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r $r $r $r $r $r $r $r $r $r $r',intArgs=(/i,k/), &
!           realArgs=(/this%tend_qbc1(i,k),this%tend_qbc2(i,k),this%tend_qoc1(i,k),this%tend_qoc2(i,k), &
!              this%tend_qdust1(i,k),this%tend_qdust2(i,k),this%tend_qdust3(i,k),this%tend_qdust4(i,k), &
!              this%tend_qdust5(i,k),this%tend_qseas1(i,k),this%tend_qseas2(i,k),this%tend_qseas3(i,k), &
!              this%tend_qseas4(i,k)/))
!   enddo
!   call mpas_log_write(' ')
!enddo


!--- initialization of meterological fields:
 call mpas_pool_get_config(configs,'config_dt',dt)

 call mpas_pool_get_array(mesh,'latCell' ,lat   )
 call mpas_pool_get_array(mesh,'lonCell' ,lon   )
 call mpas_pool_get_array(mesh,'areaCell',area  )
 call mpas_pool_get_array(mesh,'zgrid'   ,zgrid )
 call mpas_pool_get_array(mesh,'zz'      ,zz    )
 call mpas_pool_get_array(mesh,'clayfrac',clay  )
 call mpas_pool_get_array(mesh,'sandfrac',sand  )
 call mpas_pool_get_array(mesh,'erod'    ,erod  )

 call mpas_pool_get_dimension(state,'index_qc' ,index_qc )
 call mpas_pool_get_dimension(state,'index_qi' ,index_qi )
 call mpas_pool_get_dimension(state,'index_qv' ,index_qv )
 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta  ,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 call mpas_pool_get_array(state,'chems'  ,chems  ,time_lev)
 qc  => scalars(index_qc,:,:)
 qi  => scalars(index_qi,:,:)
 qv  => scalars(index_qv,:,:)

 call mpas_log_write(' ')
 call mpas_log_write('--- index_qv     = $i',intArgs=(/index_qv/))
 call mpas_log_write('--- index_qc     = $i',intArgs=(/index_qc/))
 call mpas_log_write('--- index_qi     = $i',intArgs=(/index_qi/))

 call mpas_pool_get_array(sfc_input,'isltyp'  ,isltyp  )
 call mpas_pool_get_array(sfc_input,'porosity',porosity)
 call mpas_pool_get_array(sfc_input,'smois'   ,smois   )
 call mpas_pool_get_array(sfc_input,'snowh'   ,snowh   )
 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'xland'   ,xland   )

 call mpas_pool_get_array(diag,'exner'                 ,exner )
 call mpas_pool_get_array(diag,'pressure_base'         ,pres_b)
 call mpas_pool_get_array(diag,'pressure_p'            ,pres_p)
 call mpas_pool_get_array(diag,'uReconstructZonal'     ,uzonal)
 call mpas_pool_get_array(diag,'uReconstructMeridional',umerid)
 call mpas_pool_get_array(diag,'relhum'                ,relhum)

 call mpas_pool_get_array(diag_physics,'cldfrac',cld   )
 call mpas_pool_get_array(diag_physics,'u10'    ,u10   )
 call mpas_pool_get_array(diag_physics,'v10'    ,v10   )
 call mpas_pool_get_array(diag_physics,'ust'    ,ust   )
 call mpas_pool_get_array(diag_physics,'znt'    ,znt   )

 this%dt_c   = dt
 this%grav_c = gravity

 do i = its,ite
    this%isltyp_c(i)   = isltyp(i)
    this%xlon_c(i)     = lon(i)
    this%xlat_c(i)     = lat(i)
    this%area_c(i)     = area(i)
    this%xland_c(i)    = xland(i)
    this%u10_c(i)      = u10(i)
    this%v10_c(i)      = v10(i)
    this%ust_c(i)      = ust(i)
    this%znt_c(i)      = znt(i)
    this%tsk_c(i)      = skintemp(i)
    this%porosity_c(i) = porosity(i)
    this%smois_c(i)    = smois(kts,i)
    this%snowh_c(i)    = snowh(i)

    this%clay_c(i)   = clay(i)
    this%sand_c(i)   = sand(i)
    do n = 1,nerod
       this%erod_c(i,n) = erod(n,i)
    enddo

    do k = kts,kte
       this%zgrid_c(i,k) = zgrid(k,i)
       this%rh_c(i,k)    = relhum(k,i)
       this%cld_c(i,k)   = cld(k,i)
       this%rho_c(i,k)   = zz(k,i)*rho_zz(k,i)*(1.+qv(k,i))
       this%pres_c(i,k)  = pres_b(k,i)+pres_p(k,i)
       this%u_c(i,k)     = uzonal(k,i)
       this%v_c(i,k)     = umerid(k,i)
       this%qc_c(i,k)    = qc(k,i)
       this%qi_c(i,k)    = qi(k,i)
       this%tt_c(i,k)    = theta(k,i) / (1._RKIND + R_v/R_d * max(0._RKIND,qv(k,i)))
       this%tt_c(i,k)    = this%tt_c(i,k)*exner(k,i)
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
 enddo
 if(allocated(pres2)) deallocate(pres2)

 call mpas_log_write('--- end subroutine mpas_chemistry_from_gocart:')

 end subroutine mpas_chemistry_from_gocart

!=================================================================================================================
 subroutine mpas_chemistry_to_gocart(this,state,time_lev,diag_chemistry,tend_chems)
!=================================================================================================================

!input arguments:
 integer,intent(in):: time_lev
 class(chems_gocart),intent(inout):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag_chemistry
 type(mpas_pool_type),intent(inout):: tend_chems

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,nl,ns

 integer,pointer:: index_qbc1,index_qbc2
 integer,pointer:: index_qoc1,index_qoc2
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4
 integer,pointer:: index_qsulf
 integer,pointer:: index_qdms,index_qmsa,index_qco,index_qh2o2,index_qoh,index_qso2,index_qno3

 real(kind=RKIND):: dt
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2
 real(kind=RKIND),dimension(:,:),pointer:: qoc1,qoc2
 real(kind=RKIND),dimension(:,:),pointer:: qseas1,qseas2,qseas3,qseas4
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2,qdust3,qdust4,qdust5
 real(kind=RKIND),dimension(:,:),pointer:: qsulf
 real(kind=RKIND),dimension(:,:),pointer:: qdms,qmsa,qco,qh2o2,qoh,qso2,qno3
 real(kind=RKIND),dimension(:,:,:),pointer:: chems

 real(kind=RKIND),dimension(:,:),pointer:: pm10,pm25
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaerlw,extaerlw
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw,extaersw,gaersw,waersw,bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer:: sw12aer,sw13aer,sw14aer,sw15aer,sw16aer,sw17aer

 integer,pointer:: num_gc_chems
 real(kind=RKIND),dimension(:,:,:),pointer:: gc_chems

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
 do i = its,ite
    do k = kts,kte
       this%tend_qbc1(i,k)   = (this%qbc1(i,k)   - this%tend_qbc1(i,k)  ) / dt
       this%tend_qbc2(i,k)   = (this%qbc2(i,k)   - this%tend_qbc2(i,k)  ) / dt
       this%tend_qoc1(i,k)   = (this%qoc1(i,k)   - this%tend_qoc1(i,k)  ) / dt
       this%tend_qoc2(i,k)   = (this%qoc2(i,k)   - this%tend_qoc2(i,k)  ) / dt
       this%tend_qdust1(i,k) = (this%qdust1(i,k) - this%tend_qdust1(i,k)) / dt
       this%tend_qdust2(i,k) = (this%qdust2(i,k) - this%tend_qdust2(i,k)) / dt
       this%tend_qdust3(i,k) = (this%qdust3(i,k) - this%tend_qdust3(i,k)) / dt
       this%tend_qdust4(i,k) = (this%qdust4(i,k) - this%tend_qdust4(i,k)) / dt
       this%tend_qdust5(i,k) = (this%qdust5(i,k) - this%tend_qdust5(i,k)) / dt
       this%tend_qseas1(i,k) = (this%qseas1(i,k) - this%tend_qseas1(i,k)) / dt
       this%tend_qseas2(i,k) = (this%qseas2(i,k) - this%tend_qseas2(i,k)) / dt
       this%tend_qseas3(i,k) = (this%qseas3(i,k) - this%tend_qseas3(i,k)) / dt
       this%tend_qseas4(i,k) = (this%qseas4(i,k) - this%tend_qseas4(i,k)) / dt

       this%tend_qsulf(i,k)  = (this%qsulf(i,k)  - this%tend_qsulf(i,k)) / dt
       this%tend_qdms(i,k)   = (this%qdms(i,k)   - this%tend_qdms(i,k))  / dt
       this%tend_qmsa(i,k)   = (this%qmsa(i,k)   - this%tend_qmsa(i,k))  / dt
       this%tend_qco(i,k)    = (this%qco(i,k)    - this%tend_qsulf(i,k)) / dt
       this%tend_qh2o2(i,k)  = (this%qh2o2(i,k)  - this%tend_qsulf(i,k)) / dt
       this%tend_qoh(i,k)    = (this%qoh(i,k)    - this%tend_qsulf(i,k)) / dt
       this%tend_qso2(i,k)   = (this%qso2(i,k)   - this%tend_qsulf(i,k)) / dt
       this%tend_qno3(i,k)   = (this%qno3(i,k)   - this%tend_qsulf(i,k)) / dt
    enddo
!   do k = kts,kts
!      if(this%tend_qdust1(i,k) .gt. 0._RKIND) then
!         call mpas_log_write('--- dust: $i $i $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/this%tend_qdust1(i,k), &
!                       this%tend_qdust2(i,k),this%tend_qdust3(i,k),this%tend_qdust4(i,k),this%tend_qdust5(i,k)/))
!      endif
!   enddo
!   if(this%tend_qseas1(i,kts) .gt. 0._RKIND) then
!      do k = kts,kts+5
!         call mpas_log_write('--- seas: $i $i $r $r $r $r',intArgs=(/i,k/),realArgs=(/this%tend_qseas1(i,k), &
!                             this%tend_qseas2(i,k),this%tend_qseas3(i,k),this%tend_qseas4(i,k)/))
!      enddo
!      call mpas_log_write(' ')
!   endif
 enddo


!--- update scalars:
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
 call mpas_pool_get_dimension(state,'index_qco'   ,index_qco   )
 call mpas_pool_get_dimension(state,'index_qh2o2' ,index_qh2o2 )
 call mpas_pool_get_dimension(state,'index_qoh'   ,index_qoh   )
 call mpas_pool_get_dimension(state,'index_qso2'  ,index_qso2  )
 call mpas_pool_get_dimension(state,'index_qno3'  ,index_qno3  )

 call mpas_pool_get_array(state,'chems',chems,time_lev)
 qbc1   => chems(index_qbc1,:,:)
 qbc2   => chems(index_qbc2,:,:)
 qoc1   => chems(index_qoc1,:,:)
 qoc2   => chems(index_qoc2,:,:)
 qdust1 => chems(index_qdust1,:,:)
 qdust2 => chems(index_qdust2,:,:)
 qdust3 => chems(index_qdust3,:,:)
 qdust4 => chems(index_qdust4,:,:)
 qdust5 => chems(index_qdust5,:,:)
 qseas1 => chems(index_qseas1,:,:)
 qseas2 => chems(index_qseas2,:,:)
 qseas3 => chems(index_qseas3,:,:)
 qseas4 => chems(index_qseas4,:,:)
 qsulf  => chems(index_qsulf,:,:)
 qdms   => chems(index_qdms,:,:)
 qmsa   => chems(index_qmsa,:,:)
 qco    => chems(index_qco,:,:)
 qh2o2  => chems(index_qh2o2,:,:)
 qoh    => chems(index_qoh,:,:)
 qso2   => chems(index_qso2,:,:)
 qno3   => chems(index_qno3,:,:)

 do i = its,ite
    do k = kts,kte
       qbc1(k,i)   = this%qbc1(i,k)
       qbc2(k,i)   = this%qbc2(i,k)
       qoc1(k,i)   = this%qoc1(i,k)
       qoc2(k,i)   = this%qoc2(i,k)
!      qdust1(k,i) = this%qdust1(i,k)
!      qdust2(k,i) = this%qdust2(i,k)
!      qdust3(k,i) = this%qdust3(i,k)
!      qdust4(k,i) = this%qdust4(i,k)
!      qdust5(k,i) = this%qdust5(i,k)
!      qseas1(k,i) = this%qseas1(i,k)
!      qseas2(k,i) = this%qseas2(i,k)
!      qseas3(k,i) = this%qseas3(i,k)
!      qseas4(k,i) = this%qseas4(i,k)
       qsulf(k,i)  = this%qsulf(i,k)
       qdms(k,i)   = this%qdms(i,k)
       qmsa(k,i)   = this%qmsa(i,k)
       qco(k,i)    = this%qco(i,k)
       qh2o2(k,i)  = this%qh2o2(i,k)
       qoh(k,i)    = this%qoh(i,k)
       qso2(k,i)   = this%qso2(i,k)
       qno3(k,i)   = this%qno3(i,k)
    enddo
 enddo

 call mpas_pool_get_array(tend_chems,'gc_chems',gc_chems)
 do i = its,ite
    do k = kts,kte
       gc_chems(index_qbc1,k,i)   = 0._RKIND
       gc_chems(index_qbc2,k,i)   = 0._RKIND
       gc_chems(index_qoc1,k,i)   = 0._RKIND
       gc_chems(index_qoc2,k,i)   = 0._RKIND
       gc_chems(index_qdust1,k,i) = this%tend_qdust1(i,k)
       gc_chems(index_qdust2,k,i) = this%tend_qdust2(i,k)
       gc_chems(index_qdust3,k,i) = this%tend_qdust3(i,k)
       gc_chems(index_qdust4,k,i) = this%tend_qdust4(i,k)
       gc_chems(index_qdust5,k,i) = this%tend_qdust5(i,k)
       gc_chems(index_qseas1,k,i) = this%tend_qseas1(i,k)
       gc_chems(index_qseas2,k,i) = this%tend_qseas2(i,k)
       gc_chems(index_qseas3,k,i) = this%tend_qseas3(i,k)
       gc_chems(index_qseas4,k,i) = this%tend_qseas4(i,k)
       gc_chems(index_qsulf,k,i)  = 0._RKIND
       gc_chems(index_qdms,k,i)   = 0._RKIND
       gc_chems(index_qmsa,k,i)   = 0._RKIND
       gc_chems(index_qco,k,i)    = 0._RKIND
       gc_chems(index_qh2o2,k,i)  = 0._RKIND
       gc_chems(index_qoh,k,i)    = 0._RKIND
       gc_chems(index_qso2,k,i)   = 0._RKIND
       gc_chems(index_qno3,k,i)   = 0._RKIND
    enddo
    if(this%tend_qdust1(i,kts) .gt. 0._RKIND) then
       do k = kts,kts+5
          call mpas_log_write('--- dust: $i $i $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/gc_chems(index_qdust1,k,i), &
                                    gc_chems(index_qdust2,k,i),gc_chems(index_qdust3,k,i),gc_chems(index_qdust4,k,i), &
                                    gc_chems(index_qdust5,k,i)/))
       enddo
       call mpas_log_write(' ')
    endif
 enddo


!--- compute diagnostics:
 call mpas_pool_get_array(diag_chemistry,'pm10',pm10)
 call mpas_pool_get_array(diag_chemistry,'pm25',pm25)

 call sum_pm_gocart(this) 
 do k = kts,kte
    do i = its,ite
       pm25(k,i) = this%qpm25(i,k)
       pm10(k,i) = this%qpm10(i,k)
    enddo
 enddo


!--- update longwave and shortwave aerosol optical properties:
 call mpas_pool_get_array(diag_chemistry,'tauaerlw',tauaerlw)
 call mpas_pool_get_array(diag_chemistry,'extaerlw',extaerlw)
 call mpas_pool_get_array(diag_chemistry,'tauaersw',tauaersw)
 call mpas_pool_get_array(diag_chemistry,'extaersw',extaersw)
 call mpas_pool_get_array(diag_chemistry,'gaersw'  ,gaersw  )
 call mpas_pool_get_array(diag_chemistry,'waersw'  ,waersw  )
 call mpas_pool_get_array(diag_chemistry,'bscoefsw',bscoefsw)
 call mpas_pool_get_array(diag_chemistry,'sw12aer' ,sw12aer )
 call mpas_pool_get_array(diag_chemistry,'sw13aer' ,sw13aer )
 call mpas_pool_get_array(diag_chemistry,'sw14aer' ,sw14aer )
 call mpas_pool_get_array(diag_chemistry,'sw15aer' ,sw15aer )
 call mpas_pool_get_array(diag_chemistry,'sw16aer' ,sw16aer )
 call mpas_pool_get_array(diag_chemistry,'sw17aer' ,sw17aer )

 do nl = 1,nlwbands
    do k = kts,kte
       do i = its,ite
          tauaerlw(nl,k,i) = this%tauaerlw(i,k,nl)
          extaerlw(nl,k,i) = this%extaerlw(i,k,nl)
       enddo
    enddo
 enddo

 do ns = 1,nswbands
    do k = kts,kte
       do i = its,ite
          tauaersw(ns,k,i) = this%tauaersw(i,k,ns)
          extaersw(ns,k,i) = this%extaersw(i,k,ns)
          gaersw(ns,k,i)   = this%gaersw(i,k,ns)
          waersw(ns,k,i)   = this%waersw(i,k,ns)
          bscoefsw(ns,k,i) = this%bscoefsw(i,k,ns)
          sw12aer(ns,k,i)  = this%sw12aer(i,k,ns)
          sw13aer(ns,k,i)  = this%sw13aer(i,k,ns)
          sw14aer(ns,k,i)  = this%sw14aer(i,k,ns)
          sw15aer(ns,k,i)  = this%sw15aer(i,k,ns)
          sw16aer(ns,k,i)  = this%sw16aer(i,k,ns)
          sw17aer(ns,k,i)  = this%sw17aer(i,k,ns)
       enddo
    enddo
 enddo

 call mpas_log_write('--- end subroutine mpas_chemistry_to_gocart:')

 end subroutine mpas_chemistry_to_gocart

!=================================================================================================================
 subroutine sum_pm_gocart(this)
!=================================================================================================================

!inout arguments:
 class(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte

 real(kind=RKIND),parameter:: d_2_5 = 0.380 ! ln(2.5/2)/ln(3.6/2)
 real(kind=RKIND),parameter:: s_2_5 = 0.834 ! ln(2.5/1)/ln(3/1)
 real(kind=RKIND),parameter:: d_10  = 0.737 ! ln(10/6)/ln(12/6)

 real(kind=RKIND),parameter:: airmw    = 28.97
 real(kind=RKIND),parameter:: mwso4    = 96.0576
 real(kind=RKIND),parameter:: nh4_mfac = 1.375
 real(kind=RKIND),parameter:: oc_mfac  = 1.8

 real(kind=RKIND):: sulf

 real(kind=RKIND),dimension(:,:),pointer:: qpm10,qpm25
 real(kind=RKIND),dimension(:,:),pointer:: qsulf
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2
 real(kind=RKIND),dimension(:,:),pointer:: qseas1,qseas2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine sum_pm_gocart:')

 its = this%its
 ite = this%ite
 kts = this%kts
 kte = this%kte

 qpm10  => this%qpm10
 qpm25  => this%qpm25
 qsulf  => this%qsulf
 qbc1   => this%qbc1
 qbc2   => this%qbc2
 qoc1   => this%qoc1
 qoc2   => this%qoc2
 qdust1 => this%qdust1
 qdust2 => this%qdust2
 qseas1 => this%qseas1
 qseas2 => this%qseas2

 do k = kts,kte
    do i = its,ite
       qpm25(i,k) = 0._RKIND
       qpm10(i,k) = 0._RKIND
    enddo
 enddo

 do k = kts,kte
    do i = its,ite
       sulf = qsulf(i,k)*mwso4/airmw
       qpm25(i,k) = qbc1(i,k) + qbc2(i,k) + qoc1(i,k) + qoc2(i,k) + qdust1(i,k) + d_2_5*qdust2(i,k) &
                  + qseas1(i,k) + s_2_5*qseas2(i,k) + nh4_mfac*sulf + (oc_mfac-1.)*(qoc1(i,k)+qoc2(i,k))
    enddo
 enddo

 call mpas_log_write('--- end subroutine sum_pm_gocart:')

 end subroutine sum_pm_gocart

!=================================================================================================================
 end module mpas_chemistry_gocart_interface
!=================================================================================================================

