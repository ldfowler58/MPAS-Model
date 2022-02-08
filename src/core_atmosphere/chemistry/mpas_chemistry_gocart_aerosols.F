!=================================================================================================================
 module mpas_chemistry_gocart_aerosols

 use mpas_chemistry_gocart_interface

 implicit none
 private
 public:: mpas_chemistry_gocart_aerosols_driver


 contains


!=================================================================================================================
 subroutine mpas_chemistry_gocart_aerosols_driver(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,ktep1

 real(kind=RKIND):: dt
 real(kind=RKIND),dimension(:,:),pointer:: rho,zgrid

 real(kind=RKIND),dimension(:),pointer  :: qbc1_em,qoc1_em
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chemistry_gocart_aerosols_driver:')


!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1

 dt = this%dt_c

 rho     => this%rho_c
 zgrid   => this%zgrid_c

 qbc1_em => this%qbc1_em
 qoc1_em => this%qoc1_em
 qbc1    => this%qbc1
 qbc2    => this%qbc2
 qoc1    => this%qoc1
 qoc2    => this%qoc2


!--- adds CAMS surface emissions:
 call gocart_add_CAMS_emission(its,ite,kts,kte,dt,qbc1_em,qoc1_em,rho,zgrid,qbc1,qbc2,qoc1,qoc2)


!--- converts a fraction of hydrophobic to hydrophilic black and organic carbon mixing ratios:
!call gocart_carbon_hydrophobic_to_hydrophilic(


 call mpas_log_write('--- end subroutine mpas_chemistry_gocart_aerosols_driver:')

 end subroutine mpas_chemistry_gocart_aerosols_driver

!=================================================================================================================
 subroutine gocart_add_CAMS_emission(its,ite,kts,kte,dt,qbc1_em,qoc1_em,rho,zgrid,qbc1,qbc2,qoc1,qoc2)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte

 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in),dimension(its:ite):: qbc1_em,qoc1_em
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte+1):: zgrid

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qbc1,qbc2,qoc1,qoc2

!local variables and arrays:
 integer:: i,k
 real(kind=RKIND):: mult

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart_add_CAMS_emissions:')

 k = kts
 do i = its,ite
    mult = rho(i,k)*zgrid(i,k+1)-zgrid(i,k)
    mult = dt/mult

    !--- black carbon:
    qbc1(i,k) = qbc1(i,k) + 1.e09*mult*qbc1_em(i)
    qbc2(i,k) = qbc2(i,k) + 1.e09*mult*qbc1_em(i)

    !--- organic carbon:
    qoc1(i,k) = qoc1(i,k) + 1.e09*mult*qoc1_em(i)
    qoc2(i,k) = qoc2(i,k) + 1.e09*mult*qoc1_em(i)
 enddo

 call mpas_log_write('--- end subroutine gocart_add_CAMS_emissions:')

 end subroutine gocart_add_CAMS_emission

!=================================================================================================================
 end module mpas_chemistry_gocart_aerosols
!=================================================================================================================

