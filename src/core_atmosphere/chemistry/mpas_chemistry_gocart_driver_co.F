!=================================================================================================================
 module mpas_chemistry_gocart_driver_co
 use mpas_kind_types
 use mpas_log

 use mpas_chemistry_gocart_interface

 implicit none
 private
 public:: driver_gocart_co


 contains


!=================================================================================================================
 subroutine driver_gocart_co(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,ktep1

 real(kind=RKIND):: dt,grav
 real(kind=RKIND),dimension(:),pointer:: area
 real(kind=RKIND),dimension(:,:),pointer:: dpres,rho,zgrid

 real(kind=RKIND),dimension(:),pointer  :: qco_em
 real(kind=RKIND),dimension(:,:),pointer:: qco

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_co:')

!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1

 dt    = this%dt_c
 grav  = this%grav_c

 qco_em => this%qco_em
 qco    => this%qco

 area   => this%area_c
 dpres  => this%dpres_c
 rho    => this%rho_c
 zgrid  => this%zgrid_c


!--- adds CAMS surface emissions:
 call gocart_add_CAMS_emission_co(its,ite,kts,kte,dt,grav,qco_em,rho,dpres,zgrid,qco)


 call mpas_log_write('--- end subroutine driver_gocart_co:')

 end subroutine driver_gocart_co

!=================================================================================================================
 subroutine gocart_add_CAMS_emission_co(its,ite,kts,kte,dt,grav,qco_em,rho,dpres,zgrid,qco)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte

 real(kind=RKIND),intent(in):: dt,grav
 real(kind=RKIND),intent(in),dimension(its:ite):: qco_em
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: dpres,rho
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte+1):: zgrid

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qco

!local variables and arrays:
 integer:: i,k
 real(kind=RKIND):: dz,mult

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_add_CAMS_emissions:')

 k = kts
 do i = its,ite
    dz   = zgrid(i,k+1)-zgrid(i,k)
    dz = max(0.,dz)
!   mult = rho(i,k)*(zgrid(i,k+1)-zgrid(i,k))
    mult = grav*dt/dpres(i,k)

    qco(i,k) = qco(i,k) + mult*qco_em(i)
 enddo

 call mpas_log_write('--- end subroutine gocart_add_CAMS_emissions:')

 end subroutine gocart_add_CAMS_emission_co

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_co
!=================================================================================================================
