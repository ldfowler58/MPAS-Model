!=================================================================================================================
 module mpas_chemistry_gocart_driver
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_subpool

 use mpas_chemistry_gocart_interface
 use mpas_chemistry_gocart_driver_aerosols
 use mpas_chemistry_gocart_driver_chem
 use mpas_chemistry_gocart_driver_co
 use mpas_chemistry_gocart_driver_dms
 use mpas_chemistry_gocart_driver_dust
 use mpas_chemistry_gocart_driver_opt
 use mpas_chemistry_gocart_driver_seasalt
 use mpas_chemistry_gocart_driver_sulfates


 implicit none
 private
 public:: chemistry_driver


!MPAS top gocart chemistry driver.
!Laura D. Fowler (laura@ucar.edu) / (2022-02-08).
!
!subroutine chemistry_driver is the top chemistry driver that calls the chemistry drivers for the treatment of
!gocart aerosols,gocart seasalt,and gocart dust.


 contains


!=================================================================================================================
 subroutine chemistry_driver(domain)
!=================================================================================================================

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables and pointers:
 type(mpas_pool_type),pointer:: configs
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: state
 type(mpas_pool_type),pointer:: diag
 type(mpas_pool_type),pointer:: diag_physics
 type(mpas_pool_type),pointer:: sfc_input
 type(mpas_pool_type),pointer:: CAMS_emissions
 type(mpas_pool_type),pointer:: GOCART_backgrounds
 type(mpas_pool_type),pointer:: diag_chemistry
 type(mpas_pool_type),pointer:: tend_chems
 type(block_type),pointer:: block

 type(chems_gocart):: mpas_chems

 integer:: time_lev

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart_driver:')

 time_lev = 1


 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'              ,mesh              )
    call mpas_pool_get_subpool(block%structs,'state'             ,state             )
    call mpas_pool_get_subpool(block%structs,'diag'              ,diag              )
    call mpas_pool_get_subpool(block%structs,'diag_physics'      ,diag_physics      )
    call mpas_pool_get_subpool(block%structs,'sfc_input'         ,sfc_input         )
    call mpas_pool_get_subpool(block%structs,'CAMS_emissions'    ,CAMS_emissions    )
    call mpas_pool_get_subpool(block%structs,'GOCART_backgrounds',GOCART_backgrounds)
    call mpas_pool_get_subpool(block%structs,'diag_chemistry'    ,diag_chemistry    )
    call mpas_pool_get_subpool(block%structs,'tend_chems'        ,tend_chems        )


    !--- defines dimensions used to run the GOCART chemistry:
    call mpas_chems%gocart_dims(mesh)


    !--- allocates local arrays needed to run the GOCART chemistry:
    call mpas_chems%gocart_allocate()


    !--- fills local chemistry arrays with global chemistry arrays:
    call mpas_chems%from_gocart(CAMS_emissions,GOCART_backgrounds,block%configs,mesh,state,time_lev,diag, &
                                diag_physics,sfc_input)


    !--- gocart chemistry for anthropogenic aerosols:
    call driver_gocart_aerosols(mpas_chems)


    !--- gocart chemistry for dust:
    call driver_gocart_dust(mpas_chems)


    !--- gocart chemistry for seasalt:
    call driver_gocart_seasalt(mpas_chems)


    !--- gocart_chemistry for sulfates:
    !call driver_gocart_sulfates(mpas_chems)


    !--- gocart chemistry for carbon monoxide:
    call driver_gocart_co(mpas_chems)


    !--- gocart chemistry for chemical species:
    call driver_gocart_chem(mpas_chems)


    !--- gocart chemistry for dms emissions:
    call driver_gocart_dms(mpas_chems)


    !--- gocart aerosol optical depth:
    call driver_gocart_opt(mpas_chems)


    !--- fills global chemistry arrays with local chemistry arrays:
    call mpas_chems%to_gocart(state,time_lev,diag_chemistry,tend_chems)


    !--- deallocates local arrays used to run the GOCART chemistry:
    call mpas_chems%gocart_deallocate()

    block => block % next
 end do


 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_driver:')
 call mpas_log_write(' ')

 end subroutine chemistry_driver

!=================================================================================================================
 end module mpas_chemistry_gocart_driver
!=================================================================================================================

