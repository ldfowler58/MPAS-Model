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

    real(kind=RKIND),dimension(:),pointer:: qbc1_em    => null()
    real(kind=RKIND),dimension(:),pointer:: qoc1_em    => null()
    real(kind=RKIND),dimension(:),pointer:: qco_em     => null()
    real(kind=RKIND),dimension(:),pointer:: qnh3_em    => null()
    real(kind=RKIND),dimension(:),pointer:: qso2_em    => null()

    real(kind=RKIND),dimension(:,:),pointer:: qbc1     => null()
    real(kind=RKIND),dimension(:,:),pointer:: qbc2     => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc1     => null()
    real(kind=RKIND),dimension(:,:),pointer:: qoc2     => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas1   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas3   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qseas4   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust1   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust2   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust3   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust4   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qdust5   => null()
    real(kind=RKIND),dimension(:,:),pointer:: qco      => null()
    real(kind=RKIND),dimension(:,:),pointer:: qnh3     => null()
    real(kind=RKIND),dimension(:,:),pointer:: qso2     => null()

    !--- local "water-friendly" and "ice-friendly" Thompson aerosols:
    real(kind=RKIND),dimension(:,:),pointer:: nwfa_c   => null()
    real(kind=RKIND),dimension(:,:),pointer:: nifa_c   => null()

    !--- local variables from the MPAS dynamical core:
    integer,dimension(:),pointer           :: isltyp_c => null() 

    real(kind=RKIND):: dt_c
    real(kind=RKIND):: grav_c
    real(kind=RKIND),dimension(:),pointer  :: clay_c   => null()
    real(kind=RKIND),dimension(:),pointer  :: sand_c   => null()
    real(kind=RKIND),dimension(:),pointer  :: area_c   => null()
    real(kind=RKIND),dimension(:),pointer  :: xland_c  => null()
    real(kind=RKIND),dimension(:),pointer  :: u10_c    => null()
    real(kind=RKIND),dimension(:),pointer  :: v10_c    => null()
    real(kind=RKIND),dimension(:),pointer  :: smois_c  => null()
    real(kind=RKIND),dimension(:),pointer  :: snowh_c  => null()
    real(kind=RKIND),dimension(:),pointer  :: ust_c    => null()
    real(kind=RKIND),dimension(:),pointer  :: znt_c    => null()
    real(kind=RKIND),dimension(:,:),pointer:: erod_c   => null()
    real(kind=RKIND),dimension(:,:),pointer:: zgrid_c  => null()
    real(kind=RKIND),dimension(:,:),pointer:: rho_c    => null()
    real(kind=RKIND),dimension(:,:),pointer:: u_c      => null()
    real(kind=RKIND),dimension(:,:),pointer:: v_c      => null()


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
 if(.not.associated(this%qbc1_em)) allocate(this%qbc1_em(its:ite))
 if(.not.associated(this%qoc1_em)) allocate(this%qoc1_em(its:ite))
 if(.not.associated(this%qco_em) ) allocate(this%qco_em(its:ite) )
 if(.not.associated(this%qnh3_em)) allocate(this%qnh3_em(its:ite))
 if(.not.associated(this%qso2_em)) allocate(this%qso2_em(its:ite))


!--- allocate mixing ratios:
 if(.not.associated(this%qbc1)  ) allocate(this%qbc1(its:ite,kts:kte)  )
 if(.not.associated(this%qbc2)  ) allocate(this%qbc2(its:ite,kts:kte)  )
 if(.not.associated(this%qoc1)  ) allocate(this%qoc1(its:ite,kts:kte)  )
 if(.not.associated(this%qoc2)  ) allocate(this%qoc2(its:ite,kts:kte)  )
 if(.not.associated(this%qseas1)) allocate(this%qseas1(its:ite,kts:kte))
 if(.not.associated(this%qseas2)) allocate(this%qseas2(its:ite,kts:kte))
 if(.not.associated(this%qseas3)) allocate(this%qseas3(its:ite,kts:kte))
 if(.not.associated(this%qseas4)) allocate(this%qseas4(its:ite,kts:kte))
 if(.not.associated(this%qdust1)) allocate(this%qdust1(its:ite,kts:kte))
 if(.not.associated(this%qdust2)) allocate(this%qdust2(its:ite,kts:kte))
 if(.not.associated(this%qdust3)) allocate(this%qdust3(its:ite,kts:kte))
 if(.not.associated(this%qdust4)) allocate(this%qdust4(its:ite,kts:kte))
 if(.not.associated(this%qdust5)) allocate(this%qdust5(its:ite,kts:kte))

 if(.not.associated(this%qco)   ) allocate(this%qco(its:ite,kts:kte)   )
 if(.not.associated(this%qnh3)  ) allocate(this%qnh3(its:ite,kts:kte)  )
 if(.not.associated(this%qso2)  ) allocate(this%qso2(its:ite,kts:kte)  )


!--- allocate mesh and meteorological fields:
 if(.not.associated(this%isltyp_c)) allocate(this%isltyp_c(its:ite)         )
 if(.not.associated(this%clay_c)  ) allocate(this%clay_c(its:ite)           )
 if(.not.associated(this%sand_c)  ) allocate(this%sand_c(its:ite)           )
 if(.not.associated(this%area_c)  ) allocate(this%area_c(its:ite)           )
 if(.not.associated(this%xland_c) ) allocate(this%xland_c(its:ite)          )
 if(.not.associated(this%u10_c)   ) allocate(this%u10_c(its:ite)            )
 if(.not.associated(this%v10_c)   ) allocate(this%v10_c(its:ite)            )
 if(.not.associated(this%smois_c) ) allocate(this%smois_c(its:ite)          )
 if(.not.associated(this%snowh_c) ) allocate(this%snowh_c(its:ite)          )
 if(.not.associated(this%ust_c)   ) allocate(this%ust_c(its:ite)            )
 if(.not.associated(this%znt_c)   ) allocate(this%znt_c(its:ite)            )
 if(.not.associated(this%zgrid_c) ) allocate(this%zgrid_c(its:ite,kts:ktep1))
 if(.not.associated(this%rho_c)   ) allocate(this%rho_c(its:ite,kts:kte)    )
 if(.not.associated(this%u_c)     ) allocate(this%u_c(its:ite,kts:kte)      )
 if(.not.associated(this%v_c)     ) allocate(this%v_c(its:ite,kts:kte)      )

 if(.not.associated(this%erod_c) ) allocate(this%erod_c(its:ite,nerod)     )

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
 if(associated(this%qbc1_em)) deallocate(this%qbc1_em)
 if(associated(this%qoc1_em)) deallocate(this%qoc1_em)
 if(associated(this%qco_em) ) deallocate(this%qco_em )
 if(associated(this%qnh3_em)) deallocate(this%qnh3_em)
 if(associated(this%qso2_em)) deallocate(this%qso2_em)


!--- deallocate mixing ratios:
 if(associated(this%qbc1)  ) deallocate(this%qbc1  )
 if(associated(this%qbc2)  ) deallocate(this%qbc2  )
 if(associated(this%qoc1)  ) deallocate(this%qoc1  )
 if(associated(this%qoc2)  ) deallocate(this%qoc2  )
 if(associated(this%qseas1)) deallocate(this%qseas1)
 if(associated(this%qseas2)) deallocate(this%qseas2)
 if(associated(this%qseas3)) deallocate(this%qseas3)
 if(associated(this%qseas4)) deallocate(this%qseas4)
 if(associated(this%qdust1)) deallocate(this%qdust1)
 if(associated(this%qdust2)) deallocate(this%qdust2)
 if(associated(this%qdust3)) deallocate(this%qdust3)
 if(associated(this%qdust4)) deallocate(this%qdust4)
 if(associated(this%qdust5)) deallocate(this%qdust5)

 if(associated(this%qco)   ) deallocate(this%qco   )
 if(associated(this%qnh3)  ) deallocate(this%qnh3  )
 if(associated(this%qso2)  ) deallocate(this%qso2  )


!--- deallocate mesh and meteorological fields:
 if(associated(this%isltyp_c)) deallocate(this%isltyp_c)
 if(associated(this%clay_c)  ) deallocate(this%clay_c  )
 if(associated(this%sand_c)  ) deallocate(this%sand_c  )
 if(associated(this%area_c)  ) deallocate(this%area_c  )
 if(associated(this%xland_c) ) deallocate(this%xland_c )
 if(associated(this%u10_c)   ) deallocate(this%u10_c   )
 if(associated(this%v10_c)   ) deallocate(this%v10_c   )
 if(associated(this%smois_c) ) deallocate(this%smois_c )
 if(associated(this%snowh_c) ) deallocate(this%snowh_c )
 if(associated(this%ust_c)   ) deallocate(this%ust_c   )
 if(associated(this%znt_c)   ) deallocate(this%znt_c   )
 if(associated(this%zgrid_c) ) deallocate(this%zgrid_c )
 if(associated(this%rho_c)   ) deallocate(this%rho_c   )
 if(associated(this%u_c)     ) deallocate(this%u_c     )
 if(associated(this%v_c)     ) deallocate(this%v_c     )

 if(associated(this%erod_c) ) deallocate(this%erod_c   )

 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_deallocate:')

 end subroutine mpas_chemistry_gocart_deallocate

!=================================================================================================================
 subroutine mpas_chemistry_from_gocart(this,CAMS_emissions,configs,mesh,state,time_lev,diag,diag_physics, &
                                       sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input
 type(mpas_pool_type),intent(in):: CAMS_emissions

 integer,intent(in):: time_lev

!inout arguments:
 class(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,n,nerod
 integer,pointer:: index_qbc1,index_qbc2,index_qoc1,index_qoc2,index_qseas1,index_qseas2,index_qseas3, &
              index_qseas4,index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5,index_qco, &
              index_qnh3,index_qso2

 real(kind=RKIND),dimension(:),pointer:: bc1_em_anthro,oc1_em_anthro,co_em_anthro, &
                                         nh3_em_anthro,so2_em_anthro
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2,qseas1,qseas2,qseas3,qseas4, &
                                           qdust1,qdust2,qdust3,qdust4,qdust5,qco,qnh3,qso2
 real(kind=RKIND),dimension(:,:,:),pointer:: chems

!local variables and arrays for meteorological fields:
 integer,pointer:: index_qv
 integer,dimension(:),pointer:: isltyp

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:),pointer:: area,clay,sand,xland
 real(kind=RKIND),dimension(:),pointer:: u10,v10
 real(kind=RKIND),dimension(:),pointer:: snowh,ust,znt
 real(kind=RKIND),dimension(:,:),pointer:: erod
 real(kind=RKIND),dimension(:,:),pointer:: zgrid,zz
 real(kind=RKIND),dimension(:,:),pointer:: qv,rho_zz,uzonal,umerid
 real(kind=RKIND),dimension(:,:),pointer:: smois
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_from_gocart:')

 call mpas_pool_get_array(CAMS_emissions,'bc1_em_anthro',bc1_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'oc1_em_anthro',oc1_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'co_em_anthro' ,co_em_anthro )
 call mpas_pool_get_array(CAMS_emissions,'nh3_em_anthro',nh3_em_anthro)
 call mpas_pool_get_array(CAMS_emissions,'so2_em_anthro',so2_em_anthro)

 call mpas_pool_get_dimension(state,'index_qbc1'  ,index_qbc1  )
 call mpas_pool_get_dimension(state,'index_qbc2'  ,index_qbc2  )
 call mpas_pool_get_dimension(state,'index_qoc1'  ,index_qoc1  )
 call mpas_pool_get_dimension(state,'index_qoc2'  ,index_qoc2  )
 call mpas_pool_get_dimension(state,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(state,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(state,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(state,'index_qseas4',index_qseas4)
 call mpas_pool_get_dimension(state,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(state,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(state,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(state,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(state,'index_qdust5',index_qdust5)
 call mpas_pool_get_dimension(state,'index_qco'   ,index_qco   )
 call mpas_pool_get_dimension(state,'index_qnh3'  ,index_qnh3  )
 call mpas_pool_get_dimension(state,'index_qso2'  ,index_qso2  )
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
 call mpas_log_write('--- index_qco    = $i',intArgs=(/index_qco/)   )
 call mpas_log_write('--- index_qnh3   = $i',intArgs=(/index_qnh3/)  )
 call mpas_log_write('--- index_qso2   = $i',intArgs=(/index_qso2/)  )

 call mpas_pool_get_array(state,'chems',chems,time_lev)
 qbc1   => chems(index_qbc1,:,:)
 qbc2   => chems(index_qbc2,:,:)
 qoc1   => chems(index_qoc1,:,:)
 qoc2   => chems(index_qoc2,:,:)
 qseas1 => chems(index_qseas1,:,:)
 qseas2 => chems(index_qseas2,:,:)
 qseas3 => chems(index_qseas3,:,:)
 qseas4 => chems(index_qseas4,:,:)
 qdust1 => chems(index_qdust1,:,:)
 qdust2 => chems(index_qdust2,:,:)
 qdust3 => chems(index_qdust3,:,:)
 qdust4 => chems(index_qdust4,:,:)
 qdust5 => chems(index_qdust5,:,:)
 qco    => chems(index_qco,:,:)
 qnh3   => chems(index_qnh3,:,:)
 qso2   => chems(index_qso2,:,:)

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
       this%qco(i,k)    = qco(k,i)
       this%qnh3(i,k)   = qnh3(k,i)
       this%qso2(i,k)   = qso2(k,i)
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

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_array(state,'rho_zz' ,rho_zz ,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)

 call mpas_pool_get_array(sfc_input,'isltyp',isltyp)
 call mpas_pool_get_array(sfc_input,'smois' ,smois )
 call mpas_pool_get_array(sfc_input,'snowh' ,snowh )
 call mpas_pool_get_array(sfc_input,'xland' ,xland )

 call mpas_pool_get_array(diag,'uReconstructZonal'     ,uzonal)
 call mpas_pool_get_array(diag,'uReconstructMeridional',umerid)

 call mpas_pool_get_array(diag_physics,'u10',u10)
 call mpas_pool_get_array(diag_physics,'v10',v10)
 call mpas_pool_get_array(diag_physics,'ust',ust)
 call mpas_pool_get_array(diag_physics,'znt',znt)


 this%dt_c   = dt
 this%grav_c = gravity

 do i = its,ite
    this%isltyp_c(i) = isltyp(i)
    this%area_c(i)   = area(i)
    this%xland_c(i)  = xland(i)
    this%u10_c(i)    = u10(i)
    this%v10_c(i)    = v10(i)
    this%ust_c(i)    = ust(i)
    this%znt_c(i)    = znt(i)
    this%smois_c(i)  = smois(kts,i)
    this%snowh_c(i)  = snowh(i)

    this%clay_c(i)   = clay(i)
    this%sand_c(i)   = sand(i)
    do n = 1,nerod
       this%erod_c(i,n) = erod(n,i)
    enddo

    do k = kts,kte
       this%zgrid_c(i,k) = zgrid(k,i)
       this%rho_c(i,k)   = zz(k,i)*rho_zz(k,i)*(1.+qv(k,i))
       this%u_c(i,k)     = uzonal(k,i)
       this%v_c(i,k)     = umerid(k,i)
    enddo
    k = kte+1
    this%zgrid_c(i,k) = zgrid(k,i)
 enddo

 call mpas_log_write('--- end subroutine mpas_chemistry_from_gocart:')

 end subroutine mpas_chemistry_from_gocart

!=================================================================================================================
 subroutine mpas_chemistry_to_gocart(this,state,time_lev)
!=================================================================================================================

!input arguments:
 integer,intent(in):: time_lev
 class(chems_gocart),intent(inout):: this

!inout arguments:
 type(mpas_pool_type),intent(inout):: state

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte
 integer,pointer:: index_qbc1,index_qbc2,index_qoc1,index_qoc2,index_qseas1,index_qseas2,index_qseas3, &
              index_qseas4,index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5,index_qco, &
              index_qnh3,index_qso2

 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2,qseas1,qseas2,qseas3,qseas4, &
                                           qdust1,qdust2,qdust3,qdust4,qdust5,qco,qnh3,qso2
 real(kind=RKIND),dimension(:,:,:),pointer:: chems

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_to_gocart:')

 call mpas_pool_get_dimension(state,'index_qbc1'  ,index_qbc1  )
 call mpas_pool_get_dimension(state,'index_qbc2'  ,index_qbc2  )
 call mpas_pool_get_dimension(state,'index_qoc1'  ,index_qoc1  )
 call mpas_pool_get_dimension(state,'index_qoc2'  ,index_qoc2  )
 call mpas_pool_get_dimension(state,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(state,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(state,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(state,'index_qseas4',index_qseas4)
 call mpas_pool_get_dimension(state,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(state,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(state,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(state,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(state,'index_qdust5',index_qdust5)
 call mpas_pool_get_dimension(state,'index_qco'   ,index_qco   )
 call mpas_pool_get_dimension(state,'index_qnh3'  ,index_qnh3  )
 call mpas_pool_get_dimension(state,'index_qso2'  ,index_qso2  )
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
 call mpas_log_write('--- index_qco    = $i',intArgs=(/index_qco/)   )
 call mpas_log_write('--- index_qnh3   = $i',intArgs=(/index_qnh3/)  )
 call mpas_log_write('--- index_qso2   = $i',intArgs=(/index_qso2/)  )

 call mpas_pool_get_array(state,'chems',chems,time_lev)
 qbc1   => chems(index_qbc1,:,:)
 qbc2   => chems(index_qbc2,:,:)
 qoc1   => chems(index_qoc1,:,:)
 qoc2   => chems(index_qoc2,:,:)
 qseas1 => chems(index_qseas1,:,:)
 qseas2 => chems(index_qseas2,:,:)
 qseas3 => chems(index_qseas3,:,:)
 qseas4 => chems(index_qseas4,:,:)
 qdust1 => chems(index_qdust1,:,:)
 qdust2 => chems(index_qdust2,:,:)
 qdust3 => chems(index_qdust3,:,:)
 qdust4 => chems(index_qdust4,:,:)
 qdust5 => chems(index_qdust5,:,:)
 qco    => chems(index_qco,:,:)
 qnh3   => chems(index_qnh3,:,:)
 qso2   => chems(index_qso2,:,:)

 its = this%its
 ite = this%ite
 kts = this%kts
 kte = this%kte

 do i = its,ite
    do k = kts,kte
       qbc1(k,i)   = this%qbc1(i,k)
       qbc2(k,i)   = this%qbc2(i,k)
       qoc1(k,i)   = this%qoc1(i,k)
       qoc2(k,i)   = this%qoc2(i,k)
       qseas1(k,i) = this%qseas1(i,k)
       qseas2(k,i) = this%qseas2(i,k)
       qseas3(k,k) = this%qseas3(i,k)
       qseas4(k,i) = this%qseas4(i,k)
       qdust1(k,i) = this%qdust1(i,k)
       qdust2(k,i) = this%qdust2(i,k)
       qdust3(k,i) = this%qdust3(i,k)
       qdust4(k,i) = this%qdust4(i,k)
       qdust5(k,i) = this%qdust5(i,k)
       qco(k,i)    = this%qco(i,k)
       qnh3(k,i)   = this%qnh3(i,k)
       qso2(k,k)   = this%qso2(i,k)
    enddo
 enddo

 call mpas_log_write('--- end subroutine mpas_chemistry_to_gocart:')

 end subroutine mpas_chemistry_to_gocart

!=================================================================================================================
 end module mpas_chemistry_gocart_interface
!=================================================================================================================

