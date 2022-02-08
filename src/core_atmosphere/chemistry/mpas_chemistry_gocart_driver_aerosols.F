!=================================================================================================================
 module mpas_chemistry_gocart_driver_aerosols

 use mpas_chemistry_gocart_interface
 use module_gocart_aerosols


 implicit none
 private
 public:: driver_gocart_aerosols


!MPAS gocart chemistry driver for the treatment of aerosols.
!Laura D. Fowler (laura@ucar.edu) / (2022-02-08).


 contains


!=================================================================================================================
 subroutine driver_gocart_aerosols(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,ktep1

 real(kind=RKIND):: dt
 real(kind=RKIND),dimension(:),pointer:: area
 real(kind=RKIND),dimension(:,:),pointer:: rho,zgrid

 real(kind=RKIND),dimension(:),pointer  :: qbc1_em,qoc1_em
 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_aerosols:')


!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1

 dt = this%dt_c

 area    => this%area_c
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
 call gocart_carbon_hydrophobic_to_hydrophilic(its,ite,kts,kte,ktep1,dt,qbc1,qbc2,qoc1,qoc2,rho,area,zgrid)


 call mpas_log_write('--- end subroutine driver_gocart_aerosols:')

 end subroutine driver_gocart_aerosols

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
 call mpas_log_write('--- enter subroutine gocart_add_CAMS_emissions:')

 k = kts
 do i = its,ite
    mult = rho(i,k)*zgrid(i,k+1)-zgrid(i,k)
    mult = dt/mult

    !--- black carbon (following Chin et al. (2002), we assume that 80% (20%) of black carbon is hydrophobic
    !    (hydrophilic):
    qbc1(i,k) = qbc1(i,k) + 0.8*1.e09*mult*qbc1_em(i)
    qbc2(i,k) = qbc2(i,k) + 0.2*1.e09*mult*qbc1_em(i)

    !--- organic carbon (following Chin et al. (2002), we assume that 50% (50%) of organic carbon is hydrophobic
    !    (hydrophilic):
    qoc1(i,k) = qoc1(i,k) + 0.5*1.e09*mult*qoc1_em(i)
    qoc2(i,k) = qoc2(i,k) + 0.5*1.e09*mult*qoc1_em(i)
 enddo

 call mpas_log_write('--- end subroutine gocart_add_CAMS_emissions:')

 end subroutine gocart_add_CAMS_emission

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_aerosols
!=================================================================================================================

