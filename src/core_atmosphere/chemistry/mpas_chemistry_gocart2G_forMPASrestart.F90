! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_forMPASrestart
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config, &
                              mpas_pool_get_dimension,mpas_pool_get_subpool

 use mpas_chemistry_gocart2G_vars


 implicit none
 private
 public:: gocart2G_forMPASrestart


 contains


!==================================================================================================================
 subroutine gocart2G_forMPASrestart(domain)
!==================================================================================================================

!--- input arguments:
 type(domain_type),intent(in):: domain

!--- local variables and pointers:
 type(block_type),pointer:: block
 type(mpas_pool_type),pointer:: mesh,        &
                                state,       &
                                diag_physics

 logical,pointer:: do_restart

 integer,pointer:: gocart2G_start,gocart2G_end
 integer:: i,ic,ig,its,ite,k,kts,kte
 integer:: kdvel,ndvel
 integer:: time_lev

 real(kind=RKIND),dimension(:,:,:),pointer:: drydepv
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!------------------------------------------------------------------------------------------------------------------
 call mpas_pool_get_config(domain%configs,'config_do_restart',do_restart)
 if(.not. do_restart) return


 time_lev = 1

 block => domain % blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'         ,mesh       )
    call mpas_pool_get_subpool(block%structs,'state'        ,state      )
    call mpas_pool_get_subpool(block%structs,'diag_physics',diag_physics)

    call mpas_pool_get_array(diag_physics,'bl_drydepv',drydepv)

    call mpas_pool_get_dimension(state,'gocart2G_start',gocart2G_start)
    call mpas_pool_get_dimension(state,'gocart2G_end'  ,gocart2G_end  )
    call mpas_pool_get_array(state,'scalars',scalars,time_lev)

    call mpas_chem_gocart2G%gocart2G_dims(mesh,state)
    call mpas_chem_gocart2G%gocart2G_allocate()

    its = mpas_chem_gocart2G%its
    ite = mpas_chem_gocart2G%ite
    kts = mpas_chem_gocart2G%kts
    kte = mpas_chem_gocart2G%kte
    ndvel = mpas_chem_gocart2G%ndepvel
    kdvel = mpas_chem_gocart2G%kdepvel

    ic = 0
    do ig = gocart2G_start,gocart2G_end
       ic = ic+1
       do i = its,ite
          do k = kts,kte
             mpas_chem_gocart2G%chem_mr(i,k,ic) = scalars(ig,k,i)
          enddo
          do k = kts,kdvel
             mpas_chem_gocart2G%drydepv(i,k,ic) = drydepv(ic,k,i)
          enddo
       enddo
    enddo

    block => block%next
 end do


 end subroutine gocart2G_forMPASrestart

!==================================================================================================================
 end module mpas_chemistry_gocart2G_forMPASrestart
!==================================================================================================================

