!=================================================================================================================
 module mpas_init_atm_chems_static
!=================================================================================================================
 use mpas_kind_types
 use mpas_log
 use mpas_c_interfacing,only  : mpas_f_to_c_string
 use mpas_derived_types,only  : mpas_pool_type
 use mpas_pool_routines,only  : mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension
 use mpas_init_atm_static,only: init_atm_check_read_error
!use mpas_init_atm_static,only: init_atm_check_read_error,read_geogrid

 use init_atm_hinterp,only    : interp_sequence,FOUR_POINT,SEARCH,W_AVERAGE4,W_AVERAGE16
 use init_atm_llxy,only       : latlon_to_ij,map_set,proj_info,DEG_PER_RAD,PROJ_LATLON

 use iso_c_binding,only: c_char,c_int,c_float,c_loc,c_ptr

 implicit none
 private
 public:: init_atm_chems_static


!additional interpolation of static fields to run the AFWA GOCART dust model: subroutine init_atm_chems_static
!interpolate the variables clayfrac,sandfrac,and erod from WPS_GEOG to MPAS meshes.
!Laura D. Fowler (laura@ucar.edu) / 2022-01-28.


 interface
    subroutine read_geogrid(fname,rarray,nx,ny,nz,isigned,endian, &
                            scalefactor,wordsize,status) bind(C)
       use iso_c_binding, only : c_char, c_int, c_float, c_ptr
       character (c_char), dimension(*), intent(in) :: fname
       type (c_ptr), value :: rarray
       integer (c_int), intent(in), value :: nx
       integer (c_int), intent(in), value :: ny
       integer (c_int), intent(in), value :: nz
       integer (c_int), intent(in), value :: isigned
       integer (c_int), intent(in), value :: endian
       real (c_float), intent(in), value :: scalefactor
       integer (c_int), intent(in), value :: wordsize
       integer (c_int), intent(inout) :: status
    end subroutine read_geogrid
 end interface


 contains


!=================================================================================================================
 subroutine init_atm_chems_static(configs,mesh)
!=================================================================================================================

!input arguments:
 type (mpas_pool_type),intent(in):: configs

!inout arguments:
 type (mpas_pool_type),intent(inout):: mesh

!local variables:
 type(proj_info):: proj

 character(len=StrKIND)  :: fname
 character(len=StrKIND+1):: geog_data_path
 character(len=StrKIND),pointer :: config_geog_data_path

 integer:: i,iCell,k
 integer,dimension(5):: interp_list
 integer,pointer:: nCells
 integer,dimension(:),pointer:: landmask

 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:,:),allocatable  :: clayfrac_r,sandfrac_r
 real(kind=RKIND),dimension(:,:,:),allocatable:: erod_r
 real(kind=RKIND),dimension(:),pointer  :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer  :: clayfrac,sandfrac
 real(kind=RKIND),dimension(:,:),pointer:: erod

 character(kind=c_char),dimension(StrKIND+1):: c_fname
 integer(c_int):: nx,ny,nz
 integer(c_int):: endian,isigned,istatus,wordsize
 real(kind=c_float):: missing_value,scalefactor
 real(kind=c_float),dimension(:,:,:),pointer,contiguous :: rarray
 type(c_ptr):: rarray_ptr

 real(kind=RKIND):: max_x,min_x,max_y,min_y
 real(kind=RKIND):: max_lat,min_lat,max_lon,min_lon

!--------------------------------------------------------------------------------------------------
 call mpas_log_write('')
 call mpas_log_write('--- enter subroutine init_atm_chems_static:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 call mpas_pool_get_config(configs, 'config_geog_data_path',config_geog_data_path)
 write(geog_data_path,'(a)') config_geog_data_path
 i = len_trim(geog_data_path)
 if(geog_data_path(i:i) /= '/') then
    geog_data_path(i+1:i+1) = '/'
 endif

!
! Interpolate EROD:
!
 call mpas_pool_get_array(mesh,'erod',erod)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 nx = 1440
 ny = 720
 nz = 3
 isigned  = 0
 endian   = 0
 wordsize = 2
 scalefactor   = 0.0001
 missing_value = 9.99E9 

 allocate(rarray(nx,ny,nz))
 allocate(erod_r(1:1440,1:720,3))

 call map_set(PROJ_LATLON, proj,      &
              latinc = 0.25_RKIND,    &
              loninc = 0.25_RKIND,    &
              knowni = 1.0_RKIND,     &
              knownj = 1.0_RKIND,     &
              lat1   = -89.875_RKIND, &
              lon1   = -179.875_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'erod/',1,'-',1440,'.',1,'-',0720
 call mpas_log_write(trim(fname))
 call mpas_f_to_c_string(fname,c_fname)

 rarray_ptr = c_loc(rarray)
 call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian,scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 erod_r(1:1440,1:720,1:3) = rarray(1:1440,1:720,1:3)

 do k = 1,nz
    call mpas_log_write('--- max erod_r = $r',realArgs=(/maxval(erod_r(1:1140,1:720,k))/))
    call mpas_log_write('--- min erod_r = $r',realArgs=(/minval(erod_r(1:1140,1:720,k))/))
 enddo

 erod(:,:) = 0._RKIND
 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj,lat,lon,x,y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj,lat,lon,x,y)
       else if (x >= 1440.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj,lat,lon,x,y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 719.0) y = 719.0
       do k = 1,3
          erod(k,iCell) = interp_sequence(x,y,k,erod_r,1,1440,1,720, &
                                      1,3,missing_value,interp_list,1)
       enddo
    else
       erod(:,iCell) = 0.0
    endif
 enddo

 do k = 1,nz
    call mpas_log_write('--- min erod   = $r',realArgs=(/minval(erod(k,1:nCells))/))
    call mpas_log_write('--- max erod   = $r',realArgs=(/maxval(erod(k,1:nCells))/))
 enddo

 deallocate(rarray)
 deallocate(erod_r)
 call mpas_log_write('--- end interpolate EROD.')
 call mpas_log_write(' ')


!
!--- Interpolate CLAYFRAC:
!
 call mpas_pool_get_array(mesh,'clayfrac',clayfrac)

 interp_list(1) = FOUR_POINT
 interp_list(2) = W_AVERAGE4
 interp_list(3) = W_AVERAGE16
 interp_list(4) = SEARCH
 interp_list(5) = 0

 nx = 4500
 ny = 2250
 nz = 1
 isigned  = 0
 endian   = 0
 wordsize = 2
 scalefactor   = 0.0001
 missing_value = 0.

 allocate(rarray(nx,ny,nz))
 allocate(clayfrac_r(1:4500,1:2250))
 clayfrac = 0._RKIND

 call map_set(PROJ_LATLON, proj,      &
              latinc = 0.08_RKIND,    &
              loninc = 0.08_RKIND,    &
              knowni = 1.0_RKIND,     &
              knownj = 1.0_RKIND,     &
              lat1   = -90.0_RKIND,   &
              lon1   = -180.0_RKIND) 

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'clayfrac_5m/',1,'-',4500,'.',1,'-',2250
 call mpas_log_write(trim(fname))
 call mpas_f_to_c_string(fname,c_fname)

 rarray_ptr = c_loc(rarray)
 call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian,scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 clayfrac_r(1:4500,1:2250) = rarray(1:4500,1:2250,1)

 call mpas_log_write('--- max clayfrac_r = $r',realArgs=(/maxval(clayfrac_r(1:4500,1:2250))/))
 call mpas_log_write('--- min clayfrac_r = $r',realArgs=(/minval(clayfrac_r(1:4500,1:2250))/))

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj,lat,lon,x,y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj,lat,lon,x,y)
       else if (x >= 4500.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj,lat,lon,x,y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 2249.0) y = 2249.0
       clayfrac(iCell) = interp_sequence(x,y,1,clayfrac_r,1,4500,1,2250, &
                                         1,1,missing_value,interp_list,1)
    else
       clayfrac(iCell) = 0.0
    endif
 enddo

 call mpas_log_write('--- max clayfrac  = $r',realArgs=(/maxval(clayfrac(1:nCells))/))
 call mpas_log_write('--- min clayfrac  = $r',realArgs=(/minval(clayfrac(1:nCells))/))

 deallocate(rarray)
 deallocate(clayfrac_r)
 call mpas_log_write('--- end interpolate CLAYFRAC.')
 call mpas_log_write(' ')


!
!--- Interpolate SANDFRAC:
!
 call mpas_pool_get_array(mesh,'sandfrac',sandfrac)

 nx = 4500
 ny = 2250
 nz = 1
 isigned  = 0
 endian   = 0
 wordsize = 2
 scalefactor   = 0.0001
 missing_value = 0.

 allocate(rarray(nx,ny,nz))
 allocate(sandfrac_r(1:4500,1:2250))
 sandfrac = 0._RKIND

 call map_set(PROJ_LATLON, proj,      &
              latinc = 0.08_RKIND,    &
              loninc = 0.08_RKIND,    &
              knowni = 1.0_RKIND,     &
              knownj = 1.0_RKIND,     &
              lat1   = -90.0_RKIND,   &
              lon1   = -180.0_RKIND)

 write(fname,'(a,i5.5,a1,i5.5,a1,i5.5,a1,i5.5)') trim(geog_data_path)// &
       'sandfrac_5m/',1,'-',4500,'.',1,'-',2250
 call mpas_log_write(trim(fname))
 call mpas_f_to_c_string(fname,c_fname)

 rarray_ptr = c_loc(rarray)
 call read_geogrid(c_fname,rarray_ptr,nx,ny,nz,isigned,endian,scalefactor,wordsize,istatus)
 call init_atm_check_read_error(istatus,fname)
 sandfrac_r(1:4500,1:2250) = rarray(1:4500,1:2250,1)

 call mpas_log_write('--- max sandfrac_r = $r',realArgs=(/maxval(sandfrac_r(1:4500,1:2250))/))
 call mpas_log_write('--- min sandfrac_r = $r',realArgs=(/minval(sandfrac_r(1:4500,1:2250))/))

 do iCell = 1,nCells
    if(landmask(iCell) == 1) then
       lat = latCell(iCell)*DEG_PER_RAD
       lon = lonCell(iCell)*DEG_PER_RAD
       call latlon_to_ij(proj,lat,lon,x,y)
       if(x < 0.5) then
          lon = lon + 360.0
          call latlon_to_ij(proj,lat,lon,x,y)
       else if (x >= 4500.5) then
          lon = lon - 360.0
          call latlon_to_ij(proj,lat,lon,x,y)
       end if
       if (y < 1.0) y = 1.0
       if (y > 2249.0) y = 2249.0
       sandfrac(iCell) = interp_sequence(x,y,1,sandfrac_r,1,4500,1,2250, &
                                         1,1,missing_value,interp_list,1)
    else
       sandfrac(iCell) = 0.0
    endif
 enddo

 call mpas_log_write('--- max sandfrac  = $r',realArgs=(/maxval(sandfrac(1:nCells))/))
 call mpas_log_write('--- min sandfrac  = $r',realArgs=(/minval(sandfrac(1:nCells))/))

 deallocate(rarray)
 deallocate(sandfrac_r)
 call mpas_log_write('--- end interpolate SANDFRAC.')


 call mpas_log_write('--- end subroutine init_atm_chems_static:')

 end subroutine init_atm_chems_static

!=================================================================================================================
 end module mpas_init_atm_chems_static
!=================================================================================================================
