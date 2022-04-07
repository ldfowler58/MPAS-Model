!=================================================================================================================
 module mpas_chemistry_gocart_driver_seasalts

 use mpas_chemistry_gocart_interface
 use module_gocart_seasalts


 implicit none
 private
 public:: driver_gocart_seasalts


!MPAS gocart chemistry driver for the treatment of seasalts.
!Laura D. Fowler (laura@ucar.edu) / (2022-02-08).


 contains


!=================================================================================================================
 subroutine driver_gocart_seasalts(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this


!local variables and arrays:
 integer:: i,its,ite,k,kts,kte,ktep1

 real(kind=RKIND):: dt,grav
 real(kind=RKIND),dimension(:),pointer:: area,xland,u10,v10
 real(kind=RKIND),dimension(:,:),pointer:: dpres,rho,u,v,zgrid

 real(kind=RKIND),dimension(:,:),pointer:: qseas1,qseas2,qseas3,qseas4

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_seasalts:')


!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1

 dt    = this%dt_c
 grav  = this%grav_c

 qseas1 => this%qseas1
 qseas2 => this%qseas2
 qseas3 => this%qseas3
 qseas4 => this%qseas4

 area   => this%area_c
 xland  => this%xland_c
 u10    => this%u10_c
 v10    => this%v10_c
 dpres  => this%dpres_c
 rho    => this%rho_c
 u      => this%u_c
 v      => this%v_c
 zgrid  => this%zgrid_c


!--- calls seasalts surface emission:
 call gocart_seasalts_emission(its,ite,kts,kte,ktep1,dt,grav,qseas1,qseas2,qseas3,qseas4,xland,area, &
                               dpres,zgrid,u10,v10,rho,u,v)


 call mpas_log_write('--- end subroutine driver_gocart_seasalts:')

 end subroutine driver_gocart_seasalts

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_seasalts
!=================================================================================================================

