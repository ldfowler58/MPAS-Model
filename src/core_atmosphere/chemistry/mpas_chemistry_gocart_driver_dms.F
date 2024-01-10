!=================================================================================================================
 module mpas_chemistry_gocart_driver_dms
 use mpas_kind_types
 use mpas_log

 use mpas_chemistry_gocart_interface
 use module_gocart_dms_emission


 implicit none
 private
 public:: driver_gocart_dms


 contains


!=================================================================================================================
 subroutine driver_gocart_dms(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: its,ite,kts,kte
 integer:: i,k

 real(kind=RKIND):: dt

 real(kind=RKIND),dimension(:),pointer:: backg_dms
 real(kind=RKIND),dimension(:),pointer:: area,xland
 real(kind=RKIND),dimension(:),pointer:: tsk,u10,v10
 real(kind=RKIND),dimension(:,:),pointer:: rho,tt,u,v,zgrid
 real(kind=RKIND),dimension(:,:),pointer:: qdms

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_dms:')

!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte

 dt    = this%dt_c

 backg_dms => this%backg_dms

 qdms  => this%qdms

 area  => this%area_c
 xland => this%xland_c
 tsk   => this%tsk_c
 u10   => this%u10_c
 v10   => this%v10_c
 rho   => this%rho_c
 tt    => this%tt_c
 u     => this%u_c
 v     => this%v_c
 zgrid => this%zgrid_c

!do i = its,ite
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r $r $r $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/xland(i), &
!                          area(i),tsk(i),u10(i),v10(i),rho(i,k),tt(i,k),u(i,k),v(i,k),zgrid(i,k)/))
!   enddo
!   call mpas_log_write(' ')
!enddo

 call gocart_dms_emission(its,ite,kts,kte,dt,backg_dms,area,xland,tsk,u10,v10,rho,tt,u,v,zgrid,qdms)

 call mpas_log_write('--- end subroutine driver_gocart_dms:')

 end subroutine driver_gocart_dms

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_dms
!=================================================================================================================

