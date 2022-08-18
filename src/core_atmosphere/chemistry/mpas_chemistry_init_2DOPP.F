!=================================================================================================================
 module mpas_init_chem2DOPP
 use mpas_kind_types
 use mpas_log,only          : mpas_log_write
 use mpas_io_units,only     : mpas_new_unit
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_dimension,mpas_pool_get_array

 use mpas_atmphys_constants,only: degrad
 use mpas_atmphys_utilities,only: physics_error_fatal

 implicit none
 private
 public:: init_chem2DOPP


 integer,public:: latsozp,levozp,oz_coeff,timeoz
 real(kind=RKIND),dimension(:),allocatable,public:: oz_lat,oz_pres,oz_time
 real(kind=RKIND),dimension(:,:,:,:),allocatable,public:: ozplin


 contains


!=================================================================================================================
 subroutine init_chem2DOPP(mesh,chemopp_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh


!inout arguments:
 type(mpas_pool_type),intent(inout):: chemopp_input


!local variables:
 logical,pointer:: o3climatology,o3prognostic
 integer:: i,iCell,i1,i2,k,m,n
 integer,pointer:: nCells,nMonths,nO3Coeffs,nO3Levels

 real(kind=RKIND):: dlat,dlatCell,lat,lon
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: O3pres
 real(kind=RKIND),dimension(:,:,:,:),pointer:: O3data

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_chem2DOPP:')

 call mpas_pool_get_dimension(mesh,'nCells'    ,nCells   )
 call mpas_pool_get_dimension(mesh,'nMonths'   ,nMonths  )
 call mpas_pool_get_dimension(mesh,'nOPPLevels',nO3Levels)
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs',nO3Coeffs)

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)

 call mpas_pool_get_array(chemopp_input,'OPPpres',O3pres)
 call mpas_pool_get_array(chemopp_input,'OPPdata',O3data)


!--- read monthly-mean CHEM2-OPP linearized coefficients:
 call read_o3data()


!--- make sure that the dimensions in the input O3CHem2D-OPP file are the same as the ones set in Registry.xml:
 if(nO3Levels /= levozp) then
    call mpas_log_write('--- nO3Levels = $i whereas levozp = $i',intArgs=(/nO3Levels,levozp/))
    call physics_error_fatal('--- nO3Levels is different than levozp from input data file:')
 elseif(nO3Coeffs /= oz_coeff) then
    call mpas_log_write('--- nO3Coeffs = $i whereas oz_coeff = $i',intArgs=(/nO3Coeffs,oz_coeff/))
 elseif(nMonths /= timeoz) then
    call mpas_log_write('--- nMonths = $i whereas timeiz = $i',intArgs=(/nMonths,timeoz/))
 endif


!--- interpolation of input CHEM2D-OPP coefficients to the MPAS grid: this loop is the same as in subroutine
!    init_o3climatology (see module mpas_atmphys_o3climatology, lines 149-184):
 do k = 1,nO3Levels
    O3pres(k) = oz_pres(k)
 enddo

 do iCell = 1,nCells
    lat = latCell(iCell)/degrad
    lon = lonCell(iCell)/degrad
    if(lat .gt. oz_lat(latsozp)) then
     i1 = latsozp
     i2 = latsozp
    elseif(lat .lt. oz_lat(1)) then
       i1 = 1
       i2 = 1
    else
       do i = 1, latsozp
          if(lat.ge.oz_lat(i) .and. lat.lt.oz_lat(i+1)) exit
       enddo
       i1 = i
       i2 = i+1
    endif

    do m = 1,nMonths
       do n = 1,nO3Coeffs
          do k = 1,nO3Levels
             dlat     = oz_lat(i2)-oz_lat(i1)
             dlatCell = lat-oz_lat(i1)

             if(dlat == 0.) then
                O3data(m,n,k,iCell) = ozplin(i1,k,n,m)
             else
                O3data(m,n,k,iCell) = ozplin(i1,k,n,m) &
                     + (ozplin(i2,k,n,m)-ozplin(i1,k,n,m))*dlatCell/dlat
             endif
          enddo
       enddo
    enddo
 enddo

!do m = 1,nMonths
!   do iCell = 1,nCells
!      do k = 1,nO3Levels
!         if(m.eq.1 .or. m.eq.nMonths) &
!            call mpas_log_write('$i $i $i $r $r $r $r $r $r $r',intArgs=(/m,iCell,k/),realArgs=(/O3pres(k), &
!                               (O3data(m,n,k,iCell),n=1,nO3Coeffs)/))
!      enddo
!      call mpas_log_write(' ')
!   enddo
!   call mpas_log_write(' ')
!enddo

 call mpas_log_write('--- end subroutine init_chem2DOPP:')

 end subroutine init_chem2DOPP

!=================================================================================================================
 subroutine read_o3data()
!=================================================================================================================

!local variables and arrays:
 integer:: i,k,ll,n,o3unit
 integer:: latsozc,levozc,timeozc

 real(kind=4),dimension(:),allocatable:: oz_lat4,oz_pres4,oz_time4,tempin

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine read_o3data:')

!--- define unit:
 call mpas_new_unit(o3unit)
 if(o3unit < 0) &
    call physics_error_fatal('GFSprogO3 for MPAS: All file units are taken. Change maxUnits in mpas_io_units.F')


!--- prognostic ozone input data file:
 open(unit=o3unit,file='global_o3prdlos.f77',form='unformatted',convert='big_endian')


!--- read in indices:
 read (o3unit) oz_coeff,latsozp,levozp,timeoz
 call mpas_log_write('reading in o3data from global_o3prdlos.f77: ')
 call mpas_log_write('oz_coeff = $i',intArgs=(/oz_coeff/))
 call mpas_log_write('latsozp  = $i',intArgs=(/latsozp/))
 call mpas_log_write('levozp   = $i',intArgs=(/levozp/))
 call mpas_log_write('timeoz   = $i',intArgs=(/timeoz/))


!--- read in data
!---   oz_lat   -  latitude of data        (-90 to 90)
!---   oz_pres  -  vertical pressure level (mb)
!---   oz_time  -  time coordinate         (days)
!---

 if(.not.allocated(oz_lat) ) allocate(oz_lat(latsozp)    )
 if(.not.allocated(oz_pres)) allocate(oz_pres(levozp)    )
 if(.not.allocated(oz_time)) allocate(oz_time(timeoz+1)  )

 if(.not.allocated(oz_lat4) ) allocate(oz_lat4(latsozp)  )
 if(.not.allocated(oz_pres4)) allocate(oz_pres4(levozp)  )
 if(.not.allocated(oz_time4)) allocate(oz_time4(timeoz+1))
 if(.not.allocated(tempin)  ) allocate(tempin(latsozp)   )

 rewind(o3unit)
 read(o3unit) oz_coeff,latsozp,levozp,timeoz,oz_lat4,oz_pres4,oz_time4

!--- convert pressure levels from mb to ln(Pa):
 oz_pres(:) = oz_pres4(:)
 oz_pres(:) = log(100.0*oz_pres(:))
 oz_lat(:)  = oz_lat4(:)
 oz_time(:) = oz_time4(:)

!call mpas_log_write(' ')
!call mpas_log_write('--- pressure:')
!do k = 1,levozp
!   call mpas_log_write('$i $r $r',intArgs=(/k/),realArgs=(/oz_pres4(k),oz_pres(k)/))
!enddo
!call mpas_log_write(' ')
!call mpas_log_write('--- latitude:')
!do k = 1,latsozp
!   call mpas_log_write('$i $r',intArgs=(/k/),realArgs=(/oz_lat(k)/))
!enddo
!call mpas_log_write(' ')
!call mpas_log_write('--- oz_time:')
!do k = 1,timeoz+1
!   call mpas_log_write('$i $r',intArgs=(/k/),realArgs=(/oz_time(k)/))
!enddo

!--- read in ozplin which is in order of (latitudes,ozone levels,coeff number,time) where latitudes are on
!    a uniform gaussian grid:
 if(.not.allocated(ozplin)) allocate(ozplin(latsozp,levozp,oz_coeff,timeoz))

 do i = 1,timeoz
    do n = 1,oz_coeff
       do k = 1,levozp
          read(o3unit) tempin
          ozplin(:,k,n,i) = tempin(:)
       enddo
    enddo
 enddo

 if(allocated(oz_lat4) ) deallocate(oz_lat4 )
 if(allocated(oz_pres4)) deallocate(oz_pres4)
 if(allocated(oz_time4)) deallocate(oz_time4)
 if(allocated(tempin)  ) deallocate(tempin  )


 close(o3unit)
 call mpas_log_write('--- end subroutine read_o3data:')

 end subroutine read_o3data

!=================================================================================================================
 end module mpas_init_chem2DOPP
!=================================================================================================================
