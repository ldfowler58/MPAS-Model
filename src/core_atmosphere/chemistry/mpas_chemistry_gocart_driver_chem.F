!=================================================================================================================
 module mpas_chemistry_gocart_driver_chem
 use mpas_kind_types
 use mpas_log

 use mpas_chemistry_gocart_interface
 use module_gocart_chem


 implicit none
 private
 public:: driver_gocart_chem


!MPAS gocart chemistry driver for the treatment of chemical species.
!Laura D. Fowler (laura@ucar.edu) / (2022-02-08).


 contains


!=================================================================================================================
 subroutine driver_gocart_chem(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,ktep1

 real(kind=RKIND):: dt,grav
 real(kind=RKIND),dimension(:),pointer:: xlon,xlat,area
 real(kind=RKIND),dimension(:,:),pointer:: cld,dpres,qc,qi,rho,tt,zgrid

 real(kind=RKIND),dimension(:,:),pointer:: backg_h2o2,backg_oh,backg_no3
 real(kind=RKIND),dimension(:,:),pointer:: qdms,qmsa,qso2,qsulf

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_chem:')

!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1

 dt    = this%dt_c
 grav  = this%grav_c

 backg_h2o2 => this%backg_h2o2
 backg_oh   => this%backg_oh
 backg_no3  => this%backg_no3

 qdms    => this%qdms
 qmsa    => this%qmsa
 qso2    => this%qso2
 qsulf   => this%qsulf

 xlon    => this%xlon_c
 xlat    => this%xlat_c
 area    => this%area_c
 cld     => this%cld_c
 dpres   => this%dpres_c
 qc      => this%qc_c
 qi      => this%qi_c
 rho     => this%rho_c
 tt      => this%tt_c
 zgrid   => this%zgrid_c


 call gocart_chem(its,ite,kts,kte,grav,dt,xlon,xlat,area,backg_h2o2,backg_oh,backg_no3,qdms,qso2, &
                  qsulf,qmsa,qc,qi,cld,rho,tt,zgrid)


 call mpas_log_write('--- end subroutine driver_gocart_chem:')

 end subroutine driver_gocart_chem

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_chem
!=================================================================================================================

