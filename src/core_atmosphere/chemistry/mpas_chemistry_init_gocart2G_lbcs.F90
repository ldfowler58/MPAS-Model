! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_init_gocart2G_lbcs
 use mpas_dmpar
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_stream_manager
 use mpas_timekeeping,only: sub_t_t
 use mpas_derived_types,only     : mpas_time_type,mpas_timeinterval_type
 use mpas_timekeeping,only       : mpas_get_clock_time,mpas_set_time,mpas_get_timeInterval
 use init_atm_read_met,only      : met_data,read_met_init,read_met_close,read_next_met_field
 use init_atm_hinterp,only       : interp_sequence,FOUR_POINT,SEARCH,W_AVERAGE4,W_AVERAGE16
 use init_atm_llxy,only          : latlon_to_ij,map_init,map_set,proj_info,PROJ_LATLON,PROJ_GAUSS,DEG_PER_RAD

 use mpas_chemistry_init_gocart2G_interp,only: init_hinterp_gocart2G,pressure_interp


 implicit none
 private
 public:: init_gocart2G_aerosols_lbcs


!--- this module contains all the subroutines needed to build lateral boundary conditions for the gocart2G aerosol
!    species.
!    Laura D. Fowler (laura@ucar.edu) / 2025-09-09.


 contains


!==================================================================================================================
 subroutine init_gocart2G_aerosols_lbcs(timestamp,timestart,configs,mesh,fg,diag,lbc_state)
!==================================================================================================================

!input arguments:
 character(len=StrKIND),intent(in):: timestart,timestamp
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg
 type(mpas_pool_type),intent(inout):: lbc_state

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_gocart2G_aerosols_lbcs:')

 call init_hinterp_gocart2G(configs,mesh,fg)
 call init_vinterp_gocart2G_lbcs(mesh,fg,diag,lbc_state)

 call mpas_log_write('--- end subroutine init_gocart2G_aerosols_lbcs.')
 call mpas_log_write(' ')

 end subroutine init_gocart2G_aerosols_lbcs

!==================================================================================================================
 subroutine init_vinterp_gocart2G_lbcs(mesh,fg,diag,lbc_state)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: fg
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: lbc_state

!local variables and arrays:
 integer:: k,iCell,n,nn
 integer,pointer:: nCells,nAerLevels,nVertLevels
 integer,pointer:: num_scalars,gocart2G_start,gocart2G_end
 integer,pointer:: num_scalars_fg,gocart2G_fg_start,gocart2G_fg_end
 integer,pointer:: index_qnh3,index_qnh4a,index_qso4
 integer,pointer:: index_qsoapa,index_qsoapbb,index_qsoapbg

 real(kind=RKIND),dimension(:,:),pointer:: pgoc,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars_fg
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND):: fMassAir,fMassNH3,fMassNH4a,fMassSO4,fmult
 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_vinterp_gocart2G:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nAerLevels' ,nAerLevels )

 call mpas_pool_get_array(diag,'pressure',pressure)

 call mpas_pool_get_dimension(fg,'num_scalars_fg',num_scalars_fg   )
 call mpas_pool_get_dimension(fg,'gocart2G_start',gocart2G_fg_start)
 call mpas_pool_get_dimension(fg,'gocart2G_end'  ,gocart2G_fg_end  )
 call mpas_log_write('--- num_scalars_fg     = $i',intArgs=(/num_scalars_fg/)   )
 call mpas_log_write('--- gocart2G_fg_start  = $i',intArgs=(/gocart2G_fg_start/))
 call mpas_log_write('--- gocart2G_fg_end    = $i',intArgs=(/gocart2G_fg_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(lbc_state,'num_lbc_scalars',num_scalars)
 call mpas_pool_get_dimension(lbc_state,'gocart2G_start',gocart2G_start)
 call mpas_pool_get_dimension(lbc_state,'gocart2G_end'  ,gocart2G_end  )
 call mpas_log_write('--- num_scalars    = $i',intArgs=(/num_scalars/)   )
 call mpas_log_write('--- gocart2G_start = $i',intArgs=(/gocart2G_start/))
 call mpas_log_write('--- gocart2G_end   = $i',intArgs=(/gocart2G_end/)  )

 call mpas_pool_get_array(fg,'pgoc',pgoc)
 call mpas_pool_get_array(fg,'scalars_fg',scalars_fg)
 call mpas_pool_get_array(lbc_state,'lbc_scalars',scalars)
 do nn = gocart2G_start,gocart2G_end
    scalars(nn,:,:) = 0._RKIND
 enddo


 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nAerLevels))
 n = 0
 do nn = gocart2G_start,gocart2G_end-2
    n = n+1
    do iCell = 1,nCells
       sorted_arr(1,1:nAerLevels) = 0._RKIND
       sorted_arr(2,1:nAerLevels) = 0._RKIND
       do k = 1,nAerLevels
          sorted_arr(1,k) = pgoc(k,iCell)
          sorted_arr(2,k) = scalars_fg(n,k,iCell)
       enddo
       do k = nVertLevels,1,-1
          target_p = pressure(k,iCell)
          scalars(nn,k,iCell) = pressure_interp(iCell,k,target_p,nAerLevels,sorted_arr(:,1:nAerLevels))
          if(target_p.gt.sorted_arr(1,1)) scalars(nn,k,iCell) = scalars(nn,k+1,iCell)
       enddo
    enddo
 enddo
 if(allocated(sorted_arr)) deallocate(sorted_arr)


!--- convert qnh3 from mole per mole to kg per kg. initialize qnh4a as a function of qso4 (we assume
!    that the number of moles of NH4a is equal to the number of moles of SO4, and then convert mole
!    per mole to kg per kg:
 fMassNH3  = 62._RKIND   ! as defined in NI2G_GridCompMod.F90
 fMassNH4a = 18._RKIND   ! as defined in NI2G_GridCompMod.F90
 fMassAir  = 28.97_RKIND ! as defined in NI2G_GridCompMod.F90
 fMassSO4  = 96._RKIND   ! as defined in SU2G_GridCompMod.F90

 call mpas_pool_get_dimension(lbc_state,'index_lbc_qnh3' ,index_qnh3 )
 call mpas_pool_get_dimension(lbc_state,'index_lbc_qnh4a',index_qnh4a)
 call mpas_pool_get_dimension(lbc_state,'index_lbc_qso4' ,index_qso4 )

 fmult = fMassNH3/fMassAir
 scalars(index_qnh3,:,:) = fmult*scalars(index_qnh3,:,:)

 fmult = fMassNH4a/fMassSO4
 scalars(index_qnh4a,:,:) = fmult*scalars(index_qso4,:,:)


!--- initializes precursor and simple secondary organic aerosols mixing ratios:
 call mpas_pool_get_dimension(lbc_state,'index_lbc_qsoapa' ,index_qsoapa )
 call mpas_pool_get_dimension(lbc_state,'index_lbc_qsoapbb',index_qsoapbb)
 call mpas_pool_get_dimension(lbc_state,'index_lbc_qsoapbg',index_qsoapbg)
!call mpas_log_write('--- index_qsoapa  = $i',intArgs=(/index_qsoapa/) )
 call mpas_log_write('--- index_qsoapbb = $i',intArgs=(/index_qsoapbb/))
 call mpas_log_write('--- index_qsoapbg = $i',intArgs=(/index_qsoapbg/))
 scalars(index_qsoapa,:,:)  = 0._RKIND
 scalars(index_qsoapbb,:,:) = 0._RKIND
 scalars(index_qsoapbg,:,:) = 0._RKIND


 call mpas_log_write('--- end subroutine init_vinterp_gocart2G.')

 end subroutine init_vinterp_gocart2G_lbcs

!==================================================================================================================
 end module mpas_chemistry_init_gocart2G_lbcs
!==================================================================================================================
