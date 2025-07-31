! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_init_gocart2G_static
!==================================================================================================================
 use mpas_kind_types
 use mpas_log
 use mpas_derived_types,only  : mpas_pool_type,MPAS_LOG_CRIT
 use mpas_pool_routines,only  : mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension
 use mpas_kd_tree, only       : mpas_kd_type,mpas_kd_construct
 use mpas_geotile_manager,only: mpas_geotile_mgr_type
 use mpas_init_atm_static,only: init_atm_map_static_data,max_cell_diameter


 implicit none
 private
 public:: init_gocart2G_static


 abstract interface
    function interp_criteria_function(iCell)
        integer, intent(in) :: iCell
        logical :: interp_criteria_function
    end function interp_criteria_function
 end interface

 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: nhs
 integer(kind=I8KIND),dimension(:),pointer:: soil_int
 integer(kind=I8KIND),dimension(:,:),pointer:: erod_int

 real(kind=RKIND):: soil_msgval = 0._RKIND
 real(kind=RKIND):: erod_msgval = 9.99e9
 real(kind=RKIND):: max_kdtree_distance2


 contains


!=================================================================================================================
 subroutine init_gocart2G_static(configs,mesh)
!=================================================================================================================

!input arguments:
 type (mpas_pool_type),intent(in):: configs

!inout arguments:
 type (mpas_pool_type),intent(inout):: mesh

!local variables:
 type(mpas_kd_type),pointer:: tree
 type (mpas_kd_type),dimension(:),pointer:: kd_points

 character(len=StrKIND),pointer:: config_geog_data_path
 character(len=StrKIND)  :: fname
 character(len=StrKIND+1):: geog_data_path !same as config_geog_data_path but guaranteed to have a trailing slash
 character(len=StrKIND+1):: geog_sub_path  !subdirectory names in config_geog_data_path, with trailing slash

 integer,pointer:: nCells
 integer,pointer:: supersample_fac_30s
 integer,dimension(:),pointer:: nEdgesOnCell
 integer,dimension(:,:),pointer:: verticesOnCell
 integer:: i

 real(kind=RKIND),pointer:: sphere_radius
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latVertex,lonVertex
 real(kind=RKIND),dimension(:),pointer:: xCell,yCell,zCell
 real(kind=RKIND):: max_diameter

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine init_gocart2G_static:')

 call mpas_pool_get_config(configs,'config_30s_supersample_factor',supersample_fac_30s)
 call mpas_pool_get_config(mesh,'sphere_radius',sphere_radius)

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'nEdgesOnCell'  ,nEdgesOnCell  )
 call mpas_pool_get_array(mesh,'verticesOnCell',verticesOnCell)
 call mpas_pool_get_array(mesh,'latCell'       ,latCell       )
 call mpas_pool_get_array(mesh,'lonCell'       ,lonCell       )
 call mpas_pool_get_array(mesh,'latVertex'     ,latVertex     )
 call mpas_pool_get_array(mesh,'lonVertex'     ,lonVertex     )
 call mpas_pool_get_array(mesh,'xCell'         ,xCell         )
 call mpas_pool_get_array(mesh,'yCell'         ,yCell         )
 call mpas_pool_get_array(mesh,'zCell'         ,zCell         )
 call mpas_pool_get_array(mesh,'landmask'      ,landmask      )


!
!Set max squared distance for k-d tree search to twice the squared cell diameter
!The factor of two is simply a safety factor to account for possible inaccuracies
!in the distance function used in the k-d tree
!
 max_diameter = max_cell_diameter(nCells,nEdgesOnCell,verticesOnCell,latCell,lonCell, &
                                  latVertex,lonVertex,sphere_radius)
 max_kdtree_distance2 = 2.0_RKIND * max_diameter**2


!
!Initialize the KD-Tree
!
 allocate(kd_points(nCells))
 do i = 1, nCells
    allocate(kd_points(i)%point(3))
    kd_points(i)%point = (/xCell(i),yCell(i),zCell(i)/)
    kd_points(i)%id = i ! Cell ID
 enddo
 tree => null()
 tree => mpas_kd_construct(kd_points,3)
 if(.not. associated(tree)) then
    call mpas_log_write('Error creating the KD-Tree for static interpolation', messageType=MPAS_LOG_CRIT)
 endif


 call mpas_pool_get_config(configs, 'config_geog_data_path',config_geog_data_path)
 write(geog_data_path,'(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if(geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 endif


!
! Interpolate EROD
!
 geog_sub_path = 'erod/'

 call mpas_log_write('--- start interpolate EROD')
 call interp_erod(mesh,tree,trim(geog_data_path)//trim(geog_sub_path),supersample_fac=supersample_fac_30s)
 call mpas_log_write('--- end interpolate EROD')


!
! Interpolate CLAYFRAC
!
 geog_sub_path = 'clayfrac_5m/'

 call mpas_log_write(' ')
 call mpas_log_write('--- start interpolate CLAYFRAC')
 call interp_soilfrac('clayfrac',mesh,tree,trim(geog_data_path)//trim(geog_sub_path), &
                      supersample_fac=supersample_fac_30s)
 call mpas_log_write('--- end interpolate CLAYFRAC')


!
! Interpolate SANDFRAC
!
 geog_sub_path = 'sandfrac_5m/'

 call mpas_log_write(' ')
 call mpas_log_write('--- start interpolate SANDFRAC')
 call interp_soilfrac('sandfrac',mesh,tree,trim(geog_data_path)//trim(geog_sub_path), &
                      supersample_fac=supersample_fac_30s)
 call mpas_log_write('--- end interpolate SANDFRAC')


 call mpas_log_write('--- end subroutine init_gocart2G_static.')

 end subroutine init_gocart2G_static

!==================================================================================================================
 subroutine interp_erod(mesh,kdtree,geog_data_path,supersample_fac)
!==================================================================================================================

!--- input arguments:
 type(mpas_kd_type),intent(in),pointer:: kdtree
 character(len=*),intent(in):: geog_data_path
 integer,intent(in),optional:: supersample_fac

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: mesh

!--- local variables:
 type(mpas_geotile_mgr_type):: mgr

 integer,pointer:: nCells
 integer,pointer:: nDustErosionDim
 integer:: iCell,ierr

 real(kind=RKIND),pointer:: scalefactor
 real(kind=RKIND),pointer:: missing_value
 real(kind=RKIND),dimension(:,:),pointer:: erod
 real(kind=RKIND):: erod_msgval

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine interp_erod:')
!call mpas_log_write('--- geog_data_path = '//geog_data_path)


 ierr = mgr%init(trim(geog_data_path))
 if(ierr /= 0) then
    call mpas_log_write('Error occurred initializing interpolation for '//trim(geog_data_path), &
                        messageType=MPAS_LOG_CRIT)
    return !program execution should not reach this statement since the preceding message is a critical error
 endif


 call mpas_pool_get_dimension(mesh,'nCells'         ,nCells         )
 call mpas_pool_get_dimension(mesh,'nDustErosionDim',nDustErosionDim)
 call mpas_pool_get_array(mesh,'erod',erod)

 call mpas_pool_get_config(mgr%pool,'scale_factor' ,scalefactor  )
 call mpas_pool_get_config(mgr%pool,'missing_value',missing_value)
 erod_msgval = missing_value

 allocate(nhs(nCells))
 allocate(erod_int(nDustErosionDim,nCells))


!
!store tile values as a I8KIND integer temporarily to avoid floating
!point round off differences and to have +/- 9.22x10^18 range of representative
!values. For example, a 120 km mesh with a 1 meter data set with 6 decimal of
!precision will allow for values of 180x10^12.
!
 nhs(:) = 0
 erod_int(:,:) = 0_I8KIND
 erod(:,:) = 0._RKIND

 call init_atm_map_static_data(mesh,mgr,kdtree,max_kdtree_distance2,continuous_interp_criteria, &
                               erod_interp_accumulation,supersample_fac=supersample_fac)


 do iCell = 1,nCells
    if(nhs(iCell) > 0) then
       erod(:,iCell) = real(real(erod_int(:,iCell),kind=R8KIND),kind=RKIND)
       erod(:,iCell) = erod(:,iCell) / nhs(iCell)
    endif
 enddo
 erod(:,:) = scalefactor*erod(:,:)


 deallocate(nhs)
 deallocate(erod_int)

!call mpas_log_write('--- end subroutine interp_erod.')

 end subroutine interp_erod

!==================================================================================================================
 subroutine interp_soilfrac(fieldname,mesh,kdtree,geog_data_path,supersample_fac)
!==================================================================================================================

!--- input arguments:
 type(mpas_kd_type),intent(in),pointer:: kdtree
 character(len=*),intent(in):: fieldname
 character(len=*),intent(in):: geog_data_path
 integer,intent(in),optional:: supersample_fac

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: mesh

!--- local variables:
 type(mpas_geotile_mgr_type):: mgr

 integer,pointer:: nCells
 integer:: iCell,ierr

 real(kind=RKIND),pointer:: scalefactor
 real(kind=RKIND),pointer:: missing_value
 real(kind=RKIND),dimension(:),pointer :: soil
 real(kind=RKIND):: soil_msgval

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- enter subroutine interp_soilfrac:')
!call mpas_log_write('--- geog_data_path = '//geog_data_path)


 ierr = mgr%init(trim(geog_data_path))
 if(ierr /= 0) then
    call mpas_log_write('Error occurred initializing interpolation for '//trim(geog_data_path), &
                        messageType=MPAS_LOG_CRIT)
    return !program execution should not reach this statement since the preceding message is a critical error
 endif


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_array(mesh,trim(fieldname),soil)

 call mpas_pool_get_config(mgr%pool,'scale_factor' ,scalefactor  )
 call mpas_pool_get_config(mgr%pool,'missing_value',missing_value)
 erod_msgval = missing_value

 allocate(nhs(nCells))
 allocate(soil_int(nCells))

!
!store tile values as a I8KIND integer temporarily to avoid floating
!point round off differences and to have +/- 9.22x10^18 range of representative
!values. For example, a 120 km mesh with a 1 meter data set with 6 decimal of
!precision will allow for values of 180x10^12.
!
 nhs(:) = 0
 soil_int(:) = 0_I8KIND
 soil(:) = 0._RKIND

 call init_atm_map_static_data(mesh,mgr,kdtree,max_kdtree_distance2,continuous_interp_criteria, &
                               soil_interp_accumulation,supersample_fac=supersample_fac)


 do iCell = 1,nCells
    if(nhs(iCell) > 0) then
       soil(iCell) = real(real(soil_int(iCell),kind=R8KIND),kind=RKIND)
       soil(iCell) = soil(iCell) / nhs(iCell)
    endif
 enddo
 soil(:) = scalefactor*soil(:)


 deallocate(nhs)
 deallocate(soil_int)

!call mpas_log_write('--- end subroutine interp_soilfrac.')

 end subroutine interp_soilfrac

!==================================================================================================================
 subroutine erod_interp_accumulation(iCell,pixel)
!==================================================================================================================

!input arguments:
 integer(kind=I8KIND),intent(in),dimension(:):: pixel
 integer,intent(in):: iCell

!------------------------------------------------------------------------------------------------------------------

 if(landmask(iCell) == 0) return

 if(pixel(1) /= erod_msgval) then
    erod_int(:,iCell) = erod_int(:,iCell) + int(pixel(:),kind=I8KIND)
    nhs(iCell) = nhs(iCell) + 1
 endif

 end subroutine erod_interp_accumulation

!==================================================================================================================
 subroutine soil_interp_accumulation(iCell,pixel)
!==================================================================================================================

!input arguments:
 integer(kind=I8KIND),intent(in),dimension(:):: pixel
 integer,intent(in):: iCell

!------------------------------------------------------------------------------------------------------------------

 if(landmask(iCell) == 0) return
 
 if(pixel(1) /= soil_msgval) then
    soil_int(iCell) = soil_int(iCell) + int(pixel(1),kind=I8KIND)
    nhs(iCell) = nhs(iCell) + 1
 endif

 end subroutine soil_interp_accumulation

!==================================================================================================================
 function continuous_interp_criteria(iCell)
!==================================================================================================================
 integer,intent(in) :: iCell
 logical:: continuous_interp_criteria

!------------------------------------------------------------------------------------------------------------------

 continuous_interp_criteria = .false.

 if(nhs(iCell) == 0) then
    continuous_interp_criteria = .true.
 endif

 end function continuous_interp_criteria

!==================================================================================================================
 end module  mpas_chemistry_init_gocart2G_static
!==================================================================================================================
