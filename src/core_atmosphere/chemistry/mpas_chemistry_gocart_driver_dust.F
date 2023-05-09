!=================================================================================================================
 module mpas_chemistry_gocart_driver_dust

 use mpas_chemistry_gocart_interface
 use module_gocart_dust
 use module_gocart_dust_afwa


 implicit none
 private
 public:: driver_gocart_dust


!MPAS gocart chemistry driver for the treatment of dust.
!laura D. Fowler (laura@ucar.edu) / (2022-02-08).


 contains


!=================================================================================================================
 subroutine driver_gocart_dust(this)
!=================================================================================================================

!inout arguments:
 type(chems_gocart),intent(inout):: this

!local variables and arrays:
 logical:: do_AFWA
 integer:: i,its,ite,k,kts,kte,ktep1,n
 integer:: nerod
 integer,dimension(:),pointer:: isltyp

 real(kind=RKIND):: dt,grav
 real(kind=RKIND),dimension(:),pointer:: area,xland,clay,sand,porosity,smois,snowh,ust,u10,v10,znt
 real(kind=RKIND),dimension(:,:),pointer:: erod,dpres,zgrid
 real(kind=RKIND),dimension(:,:),pointer:: rho,u,v

 real(kind=RKIND),dimension(:),pointer  :: emis_qdust1,emis_qdust2,emis_qdust3,emis_qdust4,emis_qdust5
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2,qdust3,qdust4,qdust5

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_gocart_dust:')

 do_AFWA = .false.

!--- initialization of local variables:
 its   = this%its
 ite   = this%ite
 kts   = this%kts
 kte   = this%kte
 ktep1 = this%ktep1
 nerod = this%nerod

 dt    = this%dt_c
 grav  = this%grav_c

 qdust1      => this%qdust1
 qdust2      => this%qdust2
 qdust3      => this%qdust3
 qdust4      => this%qdust4
 qdust5      => this%qdust5

 emis_qdust1 => this%emis_qdust1
 emis_qdust2 => this%emis_qdust2
 emis_qdust3 => this%emis_qdust3
 emis_qdust4 => this%emis_qdust4
 emis_qdust5 => this%emis_qdust5

 isltyp      => this%isltyp_c
 area        => this%area_c
 xland       => this%xland_c
 smois       => this%smois_c
 snowh       => this%snowh_c
 porosity    => this%porosity_c
 u10         => this%u10_c
 v10         => this%v10_c
 ust         => this%ust_c
 znt         => this%znt_c
 rho         => this%rho_c
 u           => this%u_c
 v           => this%v_c
 dpres       => this%dpres_c
 zgrid       => this%zgrid_c

 erod        => this%erod_c
 clay        => this%clay_c
 sand        => this%sand_c


!--- calls dust surface emission:
 if(do_AFWA) then
    call gocart_dust_emission_afwa(its,ite,kts,kte,ktep1,nerod,dt,grav,qdust1,qdust2,qdust3,qdust4,qdust5, &
                                    isltyp,xland,area,erod,clay,sand,zgrid,u10,v10,smois,snowh,ust,znt,rho)
 else
    call gocart_dust_emission(its,ite,kts,kte,ktep1,nerod,dt,grav,qdust1,qdust2,qdust3,qdust4,qdust5, &
                               isltyp,xland,porosity,smois,area,erod,dpres,zgrid,rho,u,v,emis_qdust1,  &
                               emis_qdust2,emis_qdust3,emis_qdust4,emis_qdust5)
 endif

 call mpas_log_write('--- end subroutine driver_gocart_dust:')

 end subroutine driver_gocart_dust

!=================================================================================================================
 end module mpas_chemistry_gocart_driver_dust
!=================================================================================================================

