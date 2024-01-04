!=================================================================================================================
 module mpas_chemistry_gocart_driver_opt
 use mpas_kind_types
 use mpas_log

 use mpas_chemistry_gocart_interface
 use module_gocart_opt


 implicit none
 private
 public:: driver_gocart_opt


 contains


!=================================================================================================================
 subroutine driver_gocart_opt(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 integer:: i,its,ite,k,kts,kte

 real(kind=RKIND),dimension(:,:),pointer:: qbc1,qbc2,qoc1,qoc2
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2,qdust3,qdust4,qdust5
 real(kind=RKIND),dimension(:,:),pointer:: qseas1,qseas2,qseas3,qseas4
 real(kind=RKIND),dimension(:,:),pointer:: qsulf,qmsa

 real(kind=RKIND),dimension(:,:),pointer:: relhum,rho,zgrid

 real(kind=RKIND),dimension(:,:,:),pointer:: tauaersw,extaersw,gaersw,waersw,bscoefsw
 real(kind=RKIND),dimension(:,:,:),pointer:: l2aer,l3aer,l4aer,l5aer,l6aer,l7aer
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaerlw,extaerlw

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_opt:')

!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte

 qbc1   => this%qbc1
 qbc2   => this%qbc2
 qoc1   => this%qoc1
 qoc2   => this%qoc2
 qdust1 => this%qdust1
 qdust2 => this%qdust2
 qdust3 => this%qdust3
 qdust4 => this%qdust4
 qdust5 => this%qdust5
 qseas1 => this%qseas1
 qseas2 => this%qseas2
 qseas3 => this%qseas3
 qseas4 => this%qseas4
 qsulf  => this%qsulf
 qmsa   => this%qmsa

 tauaerlw => this%tauaerlw
 extaerlw => this%extaerlw
 tauaersw => this%tauaersw
 extaersw => this%extaersw
 gaersw   => this%gaersw
 waersw   => this%waersw
 bscoefsw => this%bscoefsw
 l2aer    => this%sw12aer
 l3aer    => this%sw13aer
 l4aer    => this%sw14aer
 l5aer    => this%sw15aer
 l6aer    => this%sw16aer
 l7aer    => this%sw17aer

 relhum => this%rh_c
 rho    => this%rho_c
 zgrid  => this%zgrid_c


 call gocart_optical_averaging(its,ite,kts,kte,rho,relhum,zgrid,qbc1,qbc2,qoc1,qoc2,qdust1,qdust2,qdust3, &
                           qdust4,qdust5,qseas1,qseas2,qseas3,qseas4,qmsa,qsulf,tauaersw,extaersw,gaersw, &
                           waersw,bscoefsw,l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,tauaerlw,extaerlw)


 call mpas_log_write('--- end subroutine driver_gocart_opt:')

 end subroutine driver_gocart_opt

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_opt
!=================================================================================================================

