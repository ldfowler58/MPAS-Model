!=================================================================================================================
 module mpas_chemistry_gocart_driver_sulfates
 use mpas_kind_types

 use mpas_chemistry_gocart_interface

 implicit none
 private
 public:: driver_gocart_sulfates


 contains


!=================================================================================================================
 subroutine driver_gocart_sulfates(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,ktep1

 real(kind=RKIND):: dt,grav
 real(kind=RKIND),dimension(:,:),pointer:: rho,dpres,zgrid

 real(kind=RKIND),dimension(:),pointer  :: qso2_em
 real(kind=RKIND),dimension(:,:),pointer:: qso4a1,qso4a2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_sulfates:')


!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1

 dt    = this%dt_c
 grav  = this%grav_c

 qso2_em => this%qso2_em
 qso4a1  => this%qso4a1
 qso4a2  => this%qso4a2

 rho     => this%rho_c
 dpres   => this%dpres_c
 zgrid   => this%zgrid_c


!--- adds CAMS surface emissions:
 call gocart_add_CAMS_emission(its,ite,kts,kte,dt,grav,qso2_em,rho,dpres,zgrid,qso4a1,qso4a2)


 call mpas_log_write('--- end subroutine driver_gocart_sulfates:')

 end subroutine driver_gocart_sulfates

!=================================================================================================================
 subroutine gocart_add_CAMS_emission(its,ite,kts,kte,dt,grav,qso2_em,rho,dpres,zgrid,qso4a1,qso4a2)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte

 real(kind=RKIND),intent(in):: dt,grav
 real(kind=RKIND),intent(in),dimension(its:ite):: qso2_em
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: dpres,rho
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte+1):: zgrid

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qso4a1,qso4a2

!local variables and arrays:
 integer:: i,k
 real(kind=RKIND):: dz,mult

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart_add_CAMS_emissions:')

 k = kts
 do i = its,ite
    dz   = zgrid(i,k+1)-zgrid(i,k)
    dz = max(0.,dz)
!   mult = rho(i,k)*(zgrid(i,k+1)-zgrid(i,k))
    mult = grav*dt/dpres(i,k)

!   qso4a1(i,k) = qso4a1(i,k) + 0.5*mult*qso2_em(i)
!   qso4a2(i,k) = qso4a2(i,k) + 0.5*mult*qso2_em(i)
!--- increase source of so2 to 80% for accumulation mode:
    qso4a1(i,k) = qso4a1(i,k) + 0.8*mult*qso2_em(i)
    qso4a2(i,k) = qso4a2(i,k) + 0.2*mult*qso2_em(i)
 enddo

 call mpas_log_write('--- end subroutine gocart_add_CAMS_emissions:')

 end subroutine gocart_add_CAMS_emission

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_sulfates
!=================================================================================================================

