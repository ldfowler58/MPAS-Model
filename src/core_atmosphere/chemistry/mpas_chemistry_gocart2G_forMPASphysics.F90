! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_forMPASphysics
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types,only: mpas_pool_type,MPAS_LOG_CRIT
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_dimension

 use CA2G_bc_GridCompMod
 use CA2G_br_GridCompMod
 use CA2G_oc_GridCompMod
 use DU2G_GridCompMod
 use NI2G_GridCompMod
 use SS2G_GridCompMod
 use SU2G_GridCompMod

 use CA2G_bc_StateSpecs
 use CA2G_br_StateSpecs
 use CA2G_oc_StateSpecs
 use DU2G_StateSpecs
 use NI2G_StateSpecs
 use SS2G_StateSpecs
 use SU2G_StateSpecs
 use SOA2G_StateSpecs


 implicit none
 private


 type,public:: chem_gocart2G
    integer:: its,ite,jts,jte,kts,kte,ktep1
    integer:: kdepvel
    integer:: ndepvel
    integer:: nchem

    logical,dimension(:),pointer:: fthom_w
    logical,dimension(:),pointer:: fthom_i

    real(kind=RKIND),dimension(:),pointer    :: fdens     => null()
    real(kind=RKIND),dimension(:),pointer    :: fscav     => null()
    real(kind=RKIND),dimension(:),pointer    :: fnum      => null()
    real(kind=RKIND),dimension(:),pointer    :: fmeanr    => null()
    real(kind=RKIND),dimension(:),pointer    :: fhygros   => null()
    real(kind=RKIND),dimension(:,:),pointer  :: chem_nwfa => null()
    real(kind=RKIND),dimension(:,:),pointer  :: chem_ka   => null()
    real(kind=RKIND),dimension(:,:),pointer  :: chem_ra   => null()
    real(kind=RKIND),dimension(:,:),pointer  :: chem_nifa => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: drydepv   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: chem_mr   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: chem_nc   => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: chem_nr   => null()

    real(kind=RKIND),dimension(:,:,:),pointer:: chemblten => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: chemcuten => null()


    contains
       procedure:: gocart2G_dims       => gocart2G_forMPASphys_dims
       procedure:: gocart2G_allocate   => gocart2G_forMPASphys_allocate
       procedure:: gocart2G_deallocate => gocart2G_forMPASphys_deallocate
       procedure:: gocart2G_forMPASdyn
       procedure:: gocart2G_forMPASphys
       procedure:: gocart2G_forMPASphys_init
       procedure:: gocart2G_forMPASphys_mr
       procedure:: gocart2G_forMPASthom
       procedure:: gocart2G_fromMPASthom
 end type


 contains


!==================================================================================================================
 subroutine gocart2G_forMPASphys_dims(self,mesh,state)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self

!--- local variables and arrays:
 integer,pointer:: nCellsSolve,nVertLevels,kDepLevels
 integer,pointer:: num_scalars
 integer,pointer:: moist_start,moist_end
 integer,pointer:: number_start,number_end
 integer,pointer:: gocart2G_start,gocart2G_end
 integer:: num_moist,num_number,num_chem,num_gocart2G

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASphys_dims:')


 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 self%its = 1 ; self%ite = nCellsSolve
 self%jts = 1 ; self%jte = 1
 self%kts = 1 ; self%kte = nVertLevels ; self%ktep1 = nVertLevels+1
 call mpas_log_write('ITS = $i   ITE = $i',intArgs=(/self%its,self%ite/))
 call mpas_log_write('JTS = $i   JTE = $i',intArgs=(/self%jts,self%jte/))
 call mpas_log_write('KTS = $i   KTE = $i',intArgs=(/self%kts,self%kte/))


!--- computes the number of aerosol species:
!    num_scalars = number of moist arrays (array_group = "moist"), plus number of number concentration arrays
!    (array_group = "number"),  plus number of gocart2G aerosols (array_group = "gocart2G") if GOCART2G = true.

 call mpas_pool_get_dimension(state,'num_scalars' ,num_scalars )
 call mpas_pool_get_dimension(state,'moist_start' ,moist_start )
 call mpas_pool_get_dimension(state,'moist_end'   ,moist_end   )
 call mpas_pool_get_dimension(state,'number_start',number_start)
 call mpas_pool_get_dimension(state,'number_end'  ,number_end  )

 num_moist  = moist_end - moist_start + 1
 num_number = number_end - number_start + 1
 num_gocart2G = num_scalars - num_moist - num_number
 call mpas_log_write(' ')
 call mpas_log_write('--- num_scalars  = $i',intArgs=(/num_scalars/) )
 call mpas_log_write('--- num_moist    = $i',intArgs=(/num_moist/)   )
 call mpas_log_write('--- num_number   = $i',intArgs=(/num_number/)  )
 call mpas_log_write('--- num_gocart2G = $i',intArgs=(/num_gocart2G/))

 if(num_gocart2G > 0) then
    call mpas_pool_get_dimension(state,'gocart2G_start',gocart2G_start)
    call mpas_pool_get_dimension(state,'gocart2G_end'  ,gocart2G_end  )
    num_chem = gocart2G_end - gocart2G_start + 1
    if(num_chem /= num_gocart2G) &
       call mpas_log_write('--- gocart2G_tophysics: error in calculation of num_chem',messageType=MPAS_LOG_CRIT)
 else
    num_chem = 0
 endif

 self%nchem   = num_chem   ! number of aerosol species.
 self%ndepvel = self%nchem ! number of aerosol species undergoing dry deposition (as self%nchem).
 self%kdepvel = 1          ! vertical dimension of dry deposition vertical velocity from gocart2G.
 call mpas_log_write(' ')
 call mpas_log_write('--- nchem        = $i',intArgs=(/self%nchem/))
 call mpas_log_write('--- kdepvel      = $i',intArgs=(/self%kdepvel/))
 call mpas_log_write('--- ndepvel      = $i',intArgs=(/self%ndepvel/))


 call mpas_log_write('--- end subroutine gocart2G_forMPASphys_dims.')

 end subroutine gocart2G_forMPASphys_dims

!==================================================================================================================
 subroutine gocart2G_forMPASphys_allocate(self)
!==================================================================================================================

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self

!--- local variables and arrays:
 integer:: its,ite,kts,kte
 integer:: nchem,kdvel,ndvel

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASphys_allocate:')


 its   = self%its
 ite   = self%ite
 kts   = self%kts
 kte   = self%kte
 nchem = self%nchem
 ndvel = self%ndepvel
 kdvel = self%kdepvel

 if(.not.associated(self%fdens)  ) allocate(self%fdens(nchem)  )
 if(.not.associated(self%fhygros)) allocate(self%fhygros(nchem))
 if(.not.associated(self%fmeanr) ) allocate(self%fmeanr(nchem) )
 if(.not.associated(self%fnum)   ) allocate(self%fnum(nchem)   )
 if(.not.associated(self%fscav)  ) allocate(self%fscav(nchem)  )
 if(.not.associated(self%fthom_w)) allocate(self%fthom_w(nchem))
 if(.not.associated(self%fthom_i)) allocate(self%fthom_i(nchem))

 if(.not.associated(self%drydepv)  ) allocate(self%drydepv(its:ite,kdvel,ndvel)  )
 if(.not.associated(self%chem_nwfa)) allocate(self%chem_nwfa(its:ite,kts:kte)    )
 if(.not.associated(self%chem_ka)  ) allocate(self%chem_ka(its:ite,kts:kte)      )
 if(.not.associated(self%chem_ra)  ) allocate(self%chem_ra(its:ite,kts:kte)      )
 if(.not.associated(self%chem_nifa)) allocate(self%chem_nifa(its:ite,kts:kte)    )
 if(.not.associated(self%chem_mr)  ) allocate(self%chem_mr(its:ite,kts:kte,nchem))
 if(.not.associated(self%chem_nc)  ) allocate(self%chem_nc(its:ite,kts:kte,nchem))
 if(.not.associated(self%chem_nr)  ) allocate(self%chem_nr(its:ite,kts:kte,nchem))

 if(.not.associated(self%chemblten)) allocate(self%chemblten(its:ite,kts:kte,nchem))
 if(.not.associated(self%chemcuten)) allocate(self%chemcuten(its:ite,kts:kte,nchem))


 call mpas_log_write('--- end subroutine gocart2G_forMPASphys_allocate.')

 end subroutine gocart2G_forMPASphys_allocate

!==================================================================================================================
 subroutine gocart2G_forMPASphys_deallocate(self)
!==================================================================================================================

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASphys_deallocate:')


 if(associated(self%fdens)  ) deallocate(self%fdens  )
 if(associated(self%fhygros)) deallocate(self%fhygros)
 if(associated(self%fmeanr) ) deallocate(self%fmeanr )
 if(associated(self%fnum)   ) deallocate(self%fnum   )
 if(associated(self%fscav)  ) deallocate(self%fscav  )
 if(associated(self%fthom_w)) deallocate(self%fthom_w)
 if(associated(self%fthom_i)) deallocate(self%fthom_i)

 if(associated(self%drydepv)  ) deallocate(self%drydepv  )
 if(associated(self%chem_nwfa)) deallocate(self%chem_nwfa)
 if(associated(self%chem_ka)  ) deallocate(self%chem_ka  )
 if(associated(self%chem_ra)  ) deallocate(self%chem_ra  )
 if(associated(self%chem_nifa)) deallocate(self%chem_nifa)
 if(associated(self%chem_mr)  ) deallocate(self%chem_mr  )
 if(associated(self%chem_nc)  ) deallocate(self%chem_nc  )
 if(associated(self%chem_nr)  ) deallocate(self%chem_nr  )

 if(associated(self%chemblten)) deallocate(self%chemblten)
 if(associated(self%chemcuten)) deallocate(self%chemcuten)


 call mpas_log_write('--- end subroutine gocart2G_forMPASphys_deallocate.')

 end subroutine gocart2G_forMPASphys_deallocate

!==================================================================================================================
 subroutine gocart2G_forMPASphys_init(self,CA2G_bc_params,CA2G_br_params,CA2G_oc_params,DU2G_params,NI2G_params, &
                         SS2G_params,SU2G_params)
!==================================================================================================================

!--- input arguments:
 type(CA2G_bc_GridComp),intent(in):: CA2G_bc_params
 type(CA2G_br_GridComp),intent(in):: CA2G_br_params
 type(CA2G_oc_GridComp),intent(in):: CA2G_oc_params
 type(DU2G_GridComp),intent(in):: DU2G_params
 type(NI2G_GridComp),intent(in):: NI2G_params
 type(SS2G_GridComp),intent(in):: SS2G_params
 type(SU2G_GridComp),intent(in):: SU2G_params

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self

!--- local variables and pointers:
 integer:: n,nn

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASphys_init:')

 if(associated(self%fdens)  ) self%fdens(:)   = 0._RKIND
 if(associated(self%fhygros)) self%fhygros(:) = 0._RKIND
 if(associated(self%fmeanr) ) self%fmeanr(:)  = 0._RKIND
 if(associated(self%fnum)   ) self%fnum(:)    = 0._RKIND
 if(associated(self%fscav)  ) self%fscav(:)   = 0._RKIND
 if(associated(self%fthom_w)) self%fthom_w(:) = .false.
 if(associated(self%fthom_i)) self%fthom_i(:) = .false.


 n = 0
!--- black carbon:
 n = n+1
 self%fdens(n)   = CA2G_bc_params%rhop(1)  ! hydrophobic
 self%fhygros(n) = 0._RKIND                ! hydrophobic
 self%fmeanr(n)  = CA2G_bc_params%rmed(1)  ! hydrophobic
 self%fnum(n)    = CA2G_bc_params%fnum(1)  ! hydrophobic
 self%fscav(n)   = CA2G_bc_params%fscav(1) ! hydrophobic
 self%fthom_w(n) = .false.                 ! hydrophobic
 n = n+ 1
 self%fdens(n)   = CA2G_bc_params%rhop(2)  ! hydrophilic
 self%fhygros(n) = 0.10                    ! hydrophilic (as hygro_ca_aer in WRF/chem)
 self%fmeanr(n)  = CA2G_bc_params%rmed(2)  ! hydrophilic
 self%fnum(n)    = CA2G_bc_params%fnum(2)  ! hydrophilic
 self%fscav(n)   = CA2G_bc_params%fscav(2) ! hydrophilic
 self%fthom_w(n) = .true.                  ! hydrophilic

!--- brown carbon:
 n = n+1
 self%fdens(n)   = CA2G_br_params%rhop(1)  ! hydrophobic
 self%fhygros(n) = 0._RKIND                ! hydrophobic
 self%fmeanr(n)  = CA2G_br_params%rmed(1)  ! hydrophobic
 self%fnum(n)    = CA2G_br_params%fnum(1)  ! hydrophobic
 self%fscav(n)   = CA2G_br_params%fscav(1) ! hydrophobic
 self%fthom_w(n) = .false.                 ! hydrophobic
 n = n+1
 self%fdens(n)   = CA2G_br_params%rhop(2)  ! hydrophilic
 self%fhygros(n) = 0.10                    ! hydrophilic (as hygro_ca_aer in WRF/chem)
 self%fmeanr(n)  = CA2G_br_params%rmed(2)  ! hydrophilic
 self%fnum(n)    = CA2G_br_params%fnum(2)  ! hydrophilic
 self%fscav(n)   = CA2G_br_params%fscav(2) ! hydrophilic
 self%fthom_w(n) = .true.                  ! hydrophilic

!--- organic carbon:
 n = n+1
 self%fdens(n)   = CA2G_oc_params%rhop(1)  ! hydrophobic
 self%fhygros(n) = 0._RKIND                ! hydrophobic
 self%fmeanr(n)  = CA2G_oc_params%rmed(1)  ! hydrophobic
 self%fnum(n)    = CA2G_oc_params%fnum(1)  ! hydrophobic
 self%fscav(n)   = CA2G_oc_params%fscav(1) ! hydrophobic
 self%fthom_w(n) = .false.                 ! hydrophobic
 n = n+1
 self%fdens(n)   = CA2G_oc_params%rhop(2)  ! hydrophilic
 self%fhygros(n) = 0.10                    ! hydrophilic (as hygro_ca_aer in WRF/chem)
 self%fmeanr(n)  = CA2G_oc_params%rmed(2)  ! hydrophilic
 self%fnum(n)    = CA2G_oc_params%fnum(2)  ! hydrophilic
 self%fscav(n)   = CA2G_oc_params%fscav(2) ! hydrophilic
 self%fthom_w(n) = .true.                  ! hydrophilic

!--- mineral dust:
 n = n+1
 self%fdens(n)   = DU2G_params%rhop(1)     ! dust bin 1
 self%fhygros(n) = 0.1                     ! dust bin 1 (as hygro_dust_aer in WRF/chem)
 self%fmeanr(n)  = DU2G_params%rmed(1)     ! dust bin 1
 self%fnum(n)    = DU2G_params%fnum(1)     ! dust bin 1
 self%fscav(n)   = DU2G_params%fscav(1)    ! dust bin 1
 self%fthom_w(n) = .true.                  ! dust bin 1
 n = n+1
 self%fdens(n)   = DU2G_params%rhop(2)     ! dust bin 2
 self%fhygros(n) = 0.1                     ! dust bin 2 (as hygro_dust_aer in WRF/chem)
 self%fmeanr(n)  = DU2G_params%rmed(2)     ! dust bin 2
 self%fnum(n)    = DU2G_params%fnum(2)     ! dust bin 2
 self%fscav(n)   = DU2G_params%fscav(2)    ! dust bin 2
 self%fthom_w(n) = .true.                  ! dust bin 2
 n = n+1
 self%fdens(n)   = DU2G_params%rhop(3)     ! dust bin 3
 self%fhygros(n) = 0.1                     ! dust bin 3 (as hygro_dust_aer in WRF/chem)
 self%fmeanr(n)  = DU2G_params%rmed(3)     ! dust bin 3
 self%fnum(n)    = DU2G_params%fnum(3)     ! dust bin 3
 self%fscav(n)   = DU2G_params%fscav(3)    ! dust bin 3
 self%fthom_w(n) = .true.                  ! dust bin 3
 n = n+1
 self%fdens(n)   = DU2G_params%rhop(4)     ! dust bin 4
 self%fhygros(n) = 0.1                     ! dust bin 4 (as hygro_dust_aer in WRF/chem)
 self%fmeanr(n)  = DU2G_params%rmed(4)     ! dust bin 4
 self%fnum(n)    = DU2G_params%fnum(4)     ! dust bin 4
 self%fscav(n)   = DU2G_params%fscav(4)    ! dust bin 4
 self%fthom_w(n) = .true.                  ! dust bin 4
 n = n+1
 self%fdens(n)   = DU2G_params%rhop(5)     ! dust bin 5
 self%fhygros(n) = 0.1                     ! dust bin 5
 self%fmeanr(n)  = -1._RKIND               ! dust bin 5
 self%fnum(n)    = DU2G_params%fnum(5)     ! dust bin 5
 self%fscav(n)   = DU2G_params%fscav(5)    ! dust bin 5
 self%fthom_i(n) = .true.                  ! dust bin 5

!--- nitrate:
 n = n+1
 self%fdens(n)   = NI2G_params%rhop(3)     ! no3an1
 self%fhygros(n) = 0.5                     ! no3an1 (as hygro_no3_aer in WRF/chem)
 self%fmeanr(n)  = NI2G_params%rmed(3)     ! no3an1
 self%fnum(n)    = NI2G_params%fnum(3)     ! no3an1
 self%fscav(n)   = NI2G_params%fscav(3)    ! no3an1
 self%fthom_w(n) = .true.                  ! no3an1
 n = n+1
 self%fdens(n)   = NI2G_params%rhop(4)     ! no3an2
 self%fhygros(n) = 0.5                     ! no3an2 (as hygro_no3_aer in WRF/chem)
 self%fmeanr(n)  = NI2G_params%rmed(4)     ! no3an2
 self%fnum(n)    = NI2G_params%fnum(4)     ! no3an2
 self%fscav(n)   = NI2G_params%fscav(4)    ! no3an2
 self%fthom_w(n) = .true.                  ! no3an2
 n = n+1
 self%fdens(n)   = NI2G_params%rhop(5)     ! no3an3
 self%fhygros(n) = 0.5                     ! no3an3 (as hygro_no3_aer in WRF/chem)
 self%fmeanr(n)  = NI2G_params%rmed(5)     ! no3an3
 self%fnum(n)    = NI2G_params%fnum(5)     ! no3an3
 self%fscav(n)   = NI2G_params%fscav(5)    ! no3an3
 self%fthom_w(n) = .true.                  ! no3an3

!--- sulfate:
 n = n+1
 self%fdens(n)   = SU2G_params%rhop(2)     ! so2
 self%fhygros(n) = 0._RKIND                ! so2
 self%fmeanr(n)  = SU2G_params%rmed(2)     ! so2
 self%fnum(n)    = SU2G_params%fnum(2)     ! so2
 self%fscav(n)   = SU2G_params%fscav(2)    ! so2
 self%fthom_w(n) = .false.                 ! so2
 n = n+1
 self%fdens(n)   = SU2G_params%rhop(2)     ! volcanic so2
 self%fhygros(n) = 0._RKIND                ! volcanic so2
 self%fmeanr(n)  = SU2G_params%rmed(2)     ! volcanic so2
 self%fnum(n)    = SU2G_params%fnum(2)     ! volcanic so2
 self%fscav(n)   = SU2G_params%fscav(2)    ! volcanic so2
 self%fthom_w(n) = .false.                 ! volcanic so2
 n = n+1
 self%fdens(n)   = SU2G_params%rhop(3)     ! so4
 self%fhygros(n) = 0.5                     ! so4 (as hygro_so4_aer in WRF/chem)
 self%fmeanr(n)  = SU2G_params%rmed(3)     ! so4
 self%fnum(n)    = SU2G_params%fnum(3)     ! so4
 self%fscav(n)   = SU2G_params%fscav(3)    ! so4
 self%fthom_w(n) = .true.                  ! so4
 n = n+1
 self%fdens(n)   = SU2G_params%rhop(3)     ! volcanic so4
 self%fhygros(n) = 0.5                     ! volcanic so4 (as hygro_so4_aer in WRF/chem)
 self%fmeanr(n)  = SU2G_params%rmed(3)     ! volcanic so4
 self%fnum(n)    = SU2G_params%fnum(3)     ! volcanic so4
 self%fscav(n)   = SU2G_params%fscav(3)    ! volcanic so4
 self%fthom_w(n) = .true.                  ! volcanic so4

!--- sea salt:
 n = n+1
 self%fdens(n)   = SS2G_params%rhop(1)     ! sea salt bin 1
 self%fhygros(n) = 1.16                    ! sea salt bin 1 (as hygro_seas_aer in WRF/chem)
 self%fmeanr(n)  = SS2G_params%rmed(1)     ! sea salt bin 1
 self%fnum(n)    = SS2G_params%fnum(1)     ! sea salt bin 1
 self%fscav(n)   = SS2G_params%fscav(1)    ! sea salt bin 1
 self%fthom_w(n) = .true.                  ! sea salt bin 1
 n = n+1
 self%fdens(n)   = SS2G_params%rhop(2)     ! sea salt bin 2
 self%fhygros(n) = 1.16                    ! sea salt bin 2 (as hygro_seas_aer in WRF/chem)
 self%fmeanr(n)  = SS2G_params%rmed(2)     ! sea salt bin 2
 self%fnum(n)    = SS2G_params%fnum(2)     ! sea salt bin 2
 self%fscav(n)   = SS2G_params%fscav(2)    ! sea salt bin 2
 self%fthom_w(n) = .true.                  ! sea salt bin 2
 n = n+1
 self%fdens(n)   = SS2G_params%rhop(3)     ! sea salt bin 3
 self%fhygros(n) = 1.16                    ! sea salt bin 3 (as hygro_seas_aer in WRF/chem)
 self%fmeanr(n)  = SS2G_params%rmed(3)     ! sea salt bin 3
 self%fnum(n)    = SS2G_params%fnum(3)     ! sea salt bin 3
 self%fscav(n)   = SS2G_params%fscav(3)    ! sea salt bin 3
 self%fthom_w(n) = .true.                  ! sea salt bin 3
 n = n+1
 self%fdens(n)   = SS2G_params%rhop(4)     ! sea salt bin 4
 self%fhygros(n) = 1.16                    ! sea salt bin 4 (as hygro_seas_aer in WRF/chem)
 self%fmeanr(n)  = SS2G_params%rmed(4)     ! sea salt bin 4
 self%fnum(n)    = SS2G_params%fnum(4)     ! sea salt bin 4
 self%fscav(n)   = SS2G_params%fscav(4)    ! sea salt bin 4
 self%fthom_w(n) = .true.                  ! sea salt bin 4
 n = n+1
 self%fdens(n)   = SS2G_params%rhop(5)     ! sea salt bin 5
 self%fhygros(n) = 1.16                    ! sea salt bin 5 (as hygro_seas_aer in WRF/chem)
 self%fmeanr(n)  = SS2G_params%rmed(5)     ! sea salt bin 5
 self%fnum(n)    = SS2G_params%fnum(5)     ! sea salt bin 5
 self%fscav(n)   = SS2G_params%fscav(5)    ! sea salt bin 5
 self%fthom_w(n) = .true.                  ! sea salt bin 5

!--- dms and msa needed for sulfate:
 n = n+1
 self%fdens(n)   = SU2G_params%rhop(1)     ! dms
 self%fhygros(n) = 0._RKIND                ! dms
 self%fmeanr(n)  = SU2G_params%rmed(1)     ! dms
 self%fnum(n)    = SU2G_params%fnum(1)     ! dms
 self%fscav(n)   = SU2G_params%fscav(1)    ! dms
 self%fthom_w(n) = .false.                 ! dms
 n = n+1
 self%fdens(n)   = SU2G_params%rhop(4)     ! msa
 self%fhygros(n) = 0._RKIND                ! msa
 self%fmeanr(n)  = SU2G_params%rmed(4)     ! msa
 self%fnum(n)    = SU2G_params%fnum(4)     ! msa
 self%fscav(n)   = SU2G_params%fscav(4)    ! msa
 self%fthom_w(n) = .false.                 ! msa

!--- ammonia and ammonium ion needed for nitrate:
 n = n+1
 self%fdens(n)   = NI2G_params%rhop(1)     ! nh3
 self%fhygros(n) = 0._RKIND                ! nh3
 self%fmeanr(n)  = NI2G_params%rmed(1)     ! nh3
 self%fnum(n)    = NI2G_params%fnum(1)     ! nh3
 self%fscav(n)   = NI2G_params%fscav(1)    ! nh3
 self%fthom_w(n) = .false.                 ! nh3
 n = n+1
 self%fdens(n)   = NI2G_params%rhop(2)     ! nh4a
 self%fhygros(n) = 0._RKIND                ! nh4a
 self%fmeanr(n)  = NI2G_params%rmed(2)     ! nh4a
 self%fnum(n)    = NI2G_params%fnum(2)     ! nh4a
 self%fscav(n)   = NI2G_params%fscav(2)    ! nh4a
 self%fthom_w(n) = .false.                 ! nh4a

!--- secondary organic aerosols: here, we set the conversion factor between mass mixing ratios and
!    wet scavenging coefficients to that of organic carbon:
 n = n+1
 self%fdens(n)   = CA2G_oc_params%rhop(1)  ! soa (anthropogenic)
 self%fhygros(n) = 0._RKIND                ! soa (anthropogenic)
 self%fmeanr(n)  = CA2G_oc_params%rmed(1)  ! soa (anthropogenic)
 self%fnum(n)    = CA2G_oc_params%fnum(1)  ! soa (anthropogenic)
 self%fscav(n)   = CA2G_oc_params%fscav(1) ! soa (anthropogenic)
 self%fthom_w(n) = .false.                 ! soa (anthropogenic)
 n = n+1
 self%fdens(n)   = CA2G_oc_params%rhop(1)  ! soa (biomass burning)
 self%fhygros(n) = 0._RKIND                ! soa (biomass burning)
 self%fmeanr(n)  = CA2G_oc_params%rmed(1)  ! soa (biomass burning)
 self%fnum(n)    = CA2G_oc_params%fnum(1)  ! soa (biomass burning)
 self%fscav(n)   = CA2G_oc_params%fscav(1) ! soa (biomass burning)
 self%fthom_w(n) = .false.                 ! soa (biomass burning)
 n = n+1
 self%fdens(n)   = CA2G_oc_params%rhop(1)  ! soa (biogenic)
 self%fhygros(n) = 0._RKIND                ! soa (biogenic)
 self%fmeanr(n)  = CA2G_oc_params%rmed(1)  ! soa (biogenic)
 self%fnum(n)    = CA2G_oc_params%fnum(1)  ! soa (biogenic)
 self%fscav(n)   = CA2G_oc_params%fscav(1) ! soa (biogenic)
 self%fthom_w(n) = .false.                 ! soa (biogenic)


!--- initialization of arrays needed in physics: on the first timestep of a cold start run, physics is called
!    first and need to access the arrays below depending on the namelist.atmosphere options.
 self%drydepv(:,:,:) = 0._RKIND
 self%chem_mr(:,:,:) = 0._RKIND
 self%chem_nc(:,:,:) = 0._RKIND
 self%chem_nr(:,:,:) = 0._RKIND

 self%chem_ka(:,:)   = 0._RKIND
 self%chem_ra(:,:)   = 0._RKIND
 self%chem_nwfa(:,:) = 0._RKIND
 self%chem_nifa(:,:) = 0._RKIND


 call mpas_log_write('--- nchem = $i',intArgs=(/n/))
 do nn = 1,n
    call mpas_log_write('$i $l $l $r $r $r $r $r',intArgs=(/nn/),logicArgs=(/self%fthom_w(nn), &
                    self%fthom_i(nn)/),realArgs=(/self%fnum(nn),self%fscav(nn),self%fdens(nn), &
                    self%fhygros(nn),self%fmeanr(nn)/))
 enddo


 call mpas_log_write('--- end subroutine gocart2G_forMPASphys_init.')

 end subroutine gocart2G_forMPASphys_init

!==================================================================================================================
 subroutine gocart2G_forMPASphys_mr(self,state)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: state

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self

!--- local variables and pointers:
 integer,pointer:: gocart2G_start,gocart2G_end
 integer:: its,ite,kts,kte
 integer:: i,ic,ig,k

 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASphys_mr:')


 its = self%its
 ite = self%ite
 kts = self%kts
 kte = self%kte


 call mpas_pool_get_dimension(state,'gocart2G_start',gocart2G_start)
 call mpas_pool_get_dimension(state,'gocart2G_end'  ,gocart2G_end  )
 call mpas_log_write('--- gocart2G_start = $i',intArgs=(/gocart2G_start/))
 call mpas_log_write('--- gocart2G_end   = $i',intArgs=(/gocart2G_end/)  )
 call mpas_log_write('--- nchem          = $i',intArgs=(/self%nchem/)    )


 call mpas_pool_get_array(state,'scalars',scalars,1)
 ic = 0
 do ig = gocart2G_start,gocart2G_end
    ic = ic+1
    do i = its,ite
       do k = kts,kte
          self%chem_mr(i,k,ic) = scalars(ig,k,i)
       enddo
    enddo
 enddo
 call mpas_log_write('--- nchem          = $i',intArgs=(/ic/))


 call mpas_log_write('--- end subroutine gocart2G_forMPASphys_mr.')

 end subroutine gocart2G_forMPASphys_mr

!==================================================================================================================
 subroutine gocart2G_forMPASphys(self,diag_physics,CA2G_bc,CA2G_br,CA2G_oc,DU2G,NI2G,SS2G,SU2G,SOA2G)
!==================================================================================================================

!--- input arguments:
 type(CA2G_bc_State),intent(in):: CA2G_bc
 type(CA2G_br_State),intent(in):: CA2G_br
 type(CA2G_oc_State),intent(in):: CA2G_oc
 type(DU2G_State),intent(in)   :: DU2G
 type(NI2G_State),intent(in)   :: NI2G
 type(SS2G_State),intent(in)   :: SS2G
 type(SU2G_State),intent(in)   :: SU2G
 type(SOA2G_State),intent(in)  :: SOA2G

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self
 type(mpas_pool_type),intent(inout):: diag_physics

!--- local variables:
 integer:: its,ite,jts,jte,kts,kte
 integer:: i,j,k,kk,kdvel,n,nn

 real(kind=RKIND),dimension(:,:,:),pointer:: drydepv

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASphys:')


 its = self%its
 ite = self%ite
 jts = self%jts
 jte = self%jte
 kts = self%kts
 kte = self%kte
 kdvel = self%kdepvel


!--- initialize aerosol species for physics:
 if(associated(self%drydepv)) self%drydepv(:,:,:) = 0._RKIND
 if(associated(self%chem_mr)) self%chem_mr(:,:,:) = 0._RKIND
 if(associated(self%chem_nc)) self%chem_nc(:,:,:) = 0._RKIND
 if(associated(self%chem_nr)) self%chem_nr(:,:,:) = 0._RKIND

 do j = jts,jte
    do k = 1,kdvel
       do i = its,ite

          n = 0
          !--- black carbon:
          n = n + 1
          if(associated(CA2G_bc%bcphobic)) self%drydepv(i,k,n) = CA2G_bc%bcvdep(i,j)
          n = n + 1
          if(associated(CA2G_bc%bcphilic)) self%drydepv(i,k,n) = CA2G_bc%bcvdep(i,j)

          !--- brown carbon:
          n = n + 1
          if(associated(CA2G_br%brphobic)) self%drydepv(i,k,n) = CA2G_br%brvdep(i,j)
          n = n + 1
          if(associated(CA2G_br%brphilic)) self%drydepv(i,k,n) = CA2G_br%brvdep(i,j)

          !--- organic carbon:
          n = n + 1
          if(associated(CA2G_oc%ocphobic)) self%drydepv(i,k,n) = CA2G_oc%ocvdep(i,j)
          n = n + 1
          if(associated(CA2G_oc%ocphilic)) self%drydepv(i,k,n) = CA2G_oc%ocvdep(i,j)

          !--- mineral dust:
          n = n+1
          if(associated(DU2G%du)) self%drydepv(i,k,n) = DU2G%duvdep(i,j)
          n = n+1
          if(associated(DU2G%du)) self%drydepv(i,k,n) = DU2G%duvdep(i,j)
          n = n+1
          if(associated(DU2G%du)) self%drydepv(i,k,n) = DU2G%duvdep(i,j)
          n = n+1
          if(associated(DU2G%du)) self%drydepv(i,k,n) = DU2G%duvdep(i,j)
          n = n+1
          if(associated(DU2G%du)) self%drydepv(i,k,n) = DU2G%duvdep(i,j)

          !--- nitrate:
          n = n+1
          if(associated(NI2G%no3an1)) self%drydepv(i,k,n) = NI2G%nivdep(i,j)
          n = n+1
          if(associated(NI2G%no3an2)) self%drydepv(i,k,n) = NI2G%nivdep(i,j)
          n = n+1
          if(associated(NI2G%no3an3)) self%drydepv(i,k,n) = NI2G%nivdep(i,j)

          !--- sulfate:
          n = n+1
          if(associated(SU2G%so2)) self%drydepv(i,k,n) = SU2G%suvdep(i,j)
          n = n+1
          self%drydepv(i,k,n) = SU2G%suvdep(i,j) ! volcanic so2.
          n = n+1
          if(associated(SU2G%so4)) self%drydepv(i,k,n) = SU2G%suvdep(i,j)
          n = n+1
          self%drydepv(i,k,n) = SU2G%suvdep(i,j) ! volcanic so4.

          !--- sea salt:
          n = n+1
          if(associated(SS2G%ss)) self%drydepv(i,k,n) = SS2G%ssvdep(i,j)
          n = n+1
          if(associated(SS2G%ss)) self%drydepv(i,k,n) = SS2G%ssvdep(i,j)
          n = n+1
          if(associated(SS2G%ss)) self%drydepv(i,k,n) = SS2G%ssvdep(i,j)
          n = n+1
          if(associated(SS2G%ss)) self%drydepv(i,k,n) = SS2G%ssvdep(i,j)
          n = n+1
          if(associated(SS2G%ss)) self%drydepv(i,k,n) = SS2G%ssvdep(i,j)

          !--- dms and msa needed for sulfate:
          n = n+1
          if(associated(SU2G%dms)) self%drydepv(i,k,n) = SU2G%suvdep(i,j)
          n = n+1
          if(associated(SU2G%msa)) self%drydepv(i,k,n) = SU2G%suvdep(i,j)

          !--- ammonia and ammonium ion needed for nitrate:
          n = n+1
          if(associated(NI2G%nh3)) self%drydepv(i,k,n)  = NI2G%nivdep(i,j)
          n = n+1
          if(associated(NI2G%nh4a)) self%drydepv(i,k,n) = NI2G%nivdep(i,j)

          !--- secondary organic aerosols: here, we set the dry deposition vertical velocity
          !    to that of organic carbon:
          n = n+1
          if(associated(SOA2G%soap_a))  self%drydepv(i,k,n) = CA2G_oc%ocvdep(i,j)
          n = n+1
          if(associated(SOA2G%soap_bb)) self%drydepv(i,k,n) = CA2G_oc%ocvdep(i,j)
          n = n+1
          if(associated(SOA2G%soap_bg)) self%drydepv(i,k,n) = CA2G_oc%ocvdep(i,j)

       enddo
    enddo
 enddo
 call mpas_log_write('--- n = $i',intArgs=(/n/))


 do k = kts,kte
    kk = kte+1-k
    do j = jts,jte
       do i = its,ite
          n = 0
          !--- black carbon:
          n = n+1
          if(associated(CA2G_bc%bcphobic)) then
             self%chem_mr(i,kk,n) = CA2G_bc%bcphobic(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(CA2G_bc%bcphilic)) then
             self%chem_mr(i,kk,n) = CA2G_bc%bcphilic(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- brown carbon:
          n = n+1
          if(associated(CA2G_br%brphobic)) then
             self%chem_mr(i,kk,n) = CA2G_br%brphobic(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(CA2G_br%brphilic)) then
             self%chem_mr(i,kk,n) = CA2G_br%brphilic(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- organic carbon:
          n = n+1
          if(associated(CA2G_oc%ocphobic)) then
             self%chem_mr(i,kk,n) = CA2G_oc%ocphobic(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(CA2G_oc%ocphilic)) then
             self%chem_mr(i,kk,n) = CA2G_oc%ocphilic(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- mineral dust:
          n = n+1
          if(associated(DU2G%du)) then
             self%chem_mr(i,kk,n) = DU2G%du(i,j,k,1)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(DU2G%du)) then
             self%chem_mr(i,kk,n) = DU2G%du(i,j,k,2)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(DU2G%du)) then
             self%chem_mr(i,kk,n) = DU2G%du(i,j,k,3)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(DU2G%du)) then
             self%chem_mr(i,kk,n) = DU2G%du(i,j,k,4)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(DU2G%du)) then
             self%chem_mr(i,kk,n) = DU2G%du(i,j,k,5)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- nitrate:
          n = n+1
          if(associated(NI2G%no3an1)) then
             self%chem_mr(i,kk,n) = NI2G%no3an1(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(NI2G%no3an2)) then
             self%chem_mr(i,kk,n) = NI2G%no3an2(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(NI2G%no3an3)) then
             self%chem_mr(i,kk,n) = NI2G%no3an3(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- sulfate:
          n = n+1
          if(associated(SU2G%so2)) then
             self%chem_mr(i,kk,n) = SU2G%so2(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          self%chem_mr(i,kk,n)    = 0._RKIND ! volcanic so2.
          self%chem_nc(i,kk,n)    = 0._RKIND ! volcanic so2.
          n = n+1
          if(associated(SU2G%so4)) then
             self%chem_mr(i,kk,n) = SU2G%so4(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          self%chem_mr(i,kk,n)    = 0._RKIND ! volcanic so4.
          self%chem_nc(i,kk,n)    = 0._RKIND ! volcanic so4.

          !--- sea salt:
          n = n+1
          if(associated(SS2G%ss)) then
             self%chem_mr(i,kk,n) = SS2G%ss(i,j,k,1)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SS2G%ss)) then
             self%chem_mr(i,kk,n) = SS2G%ss(i,j,k,2)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SS2G%ss)) then
             self%chem_mr(i,kk,n) = SS2G%ss(i,j,k,3)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SS2G%ss)) then
             self%chem_mr(i,kk,n) = SS2G%ss(i,j,k,4)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SS2G%ss)) then
             self%chem_mr(i,kk,n) = SS2G%ss(i,j,k,5)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- dms and msa needed for sulfate:
          n = n+1
          if(associated(SU2G%dms)) then
             self%chem_mr(i,kk,n) = SU2G%dms(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SU2G%msa)) then
             self%chem_mr(i,kk,n) = SU2G%msa(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- ammonia and ammonium ion needed for nitrate:
          n = n+1
          if(associated(NI2G%nh3)) then
             self%chem_mr(i,kk,n) = NI2G%nh3(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(NI2G%nh4a)) then
             self%chem_mr(i,kk,n) = NI2G%nh4a(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif

          !--- secondary organic aerosols:
          n = n+1
          if(associated(SOA2G%soap_a)) then
             self%chem_mr(i,kk,n) = SOA2G%soap_a(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SOA2G%soap_bb)) then
             self%chem_mr(i,kk,n) = SOA2G%soap_bb(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
          n = n+1
          if(associated(SOA2G%soap_bg)) then
             self%chem_mr(i,kk,n) = SOA2G%soap_bg(i,j,k)
             self%chem_nc(i,kk,n) = self%fnum(n)*self%chem_mr(i,kk,n)
          endif
       enddo
    enddo
 enddo
 call mpas_log_write('--- n = $i',intArgs=(/n/))


!--- save dry deposition vertical velocities to diag_physics:
 call mpas_pool_get_array(diag_physics,'bl_drydepv',drydepv)

 do n = 1,self%ndepvel
    do k = kts,kdvel
       do i = its,ite
          drydepv(n,k,i) = self%drydepv(i,k,n)
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine gocart2G_forMPASphys.')

 end subroutine gocart2G_forMPASphys

!==================================================================================================================
 subroutine gocart2G_forMPASthom(self,diag,state,time_lev)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: state
 integer,intent(in):: time_lev

!--- inout arguments:
 class(chem_gocart2G),intent(inout):: self

!--- local variables and arrays:
 integer,pointer:: gocart2G_start,gocart2G_end
 integer:: its,ite,kts,kte
 integer:: i,ic,ig,k

 real(kind=RKIND):: sum_hygro,sum_meanr,sum_vol
 real(kind=RKIND),dimension(:,:),pointer:: rho
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASthom:')


!--- initializes the mass mixing ratios of gocart2G mass mixing ratios and number concentration
!    after call to dynamical solver:
 its = self%its
 ite = self%ite
 kts = self%kts
 kte = self%kte

 call mpas_pool_get_array(diag,'rho',rho)

 call mpas_pool_get_dimension(state,'gocart2G_start',gocart2G_start)
 call mpas_pool_get_dimension(state,'gocart2G_end'  ,gocart2G_end  )
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)


 ic = 0
 do ig = gocart2G_start,gocart2G_end
    ic = ic+1
    do i = its,ite
       do k = kts,kte
          self%chem_mr(i,k,ic) = scalars(ig,k,i)
       enddo
    enddo
    if(self%fthom_w(ic) .or. self%fthom_i(ic)) self%chem_nc(:,:,ic) = self%fnum(ic)*self%chem_mr(:,:,ic)
 enddo


!--- computes the gocart2G water-friendly and ice-friendly aerosols for use in the Thompson cloud microphysics
!    scheme (see Thompson and Eidhammer, 2014). Water-friendly aerosols (nwfa) include hydrophilic black carbon,
!    brown carbon, organic carbon, sulfate, sea-salt (all five bins), and mineral dust (first four bins). Ice-
!    friendly aerosols (nifa) include the last bin of mineral dust:
 self%chem_nwfa(:,:) = 0._RKIND
 self%chem_nifa(:,:) = 0._RKIND
 do ic = 1,self%nchem
    if(self%fthom_w(ic)) self%chem_nwfa(:,:) = self%chem_nwfa(:,:) + self%chem_nc(:,:,ic)
    if(self%fthom_i(ic)) self%chem_nifa(:,:) = self%chem_nifa(:,:) + self%chem_nc(:,:,ic)
 enddo
 do ic = 1,self%nchem
    if(self%fthom_w(ic)) then
       where(self%chem_nwfa(:,:).gt.0._RKIND) self%chem_nr(:,:,ic) = self%chem_nc(:,:,ic) / self%chem_nwfa(:,:)
    endif
    if(self%fthom_i(ic)) then
       where(self%chem_nifa(:,:).gt.0._RKIND) self%chem_nr(:,:,ic) = self%chem_nc(:,:,ic) / self%chem_nifa(:,:)
    endif
 enddo


!--- computes the mean hygroscopicity and mean radius of the gocart2G water-friendly aerosols for use in the
!    Thompson cloud microphysics scheme (see Thompson and Eidhammer, 2014).
 do i = its,ite
    do k = kts,kte
       sum_hygro = 0._RKIND
       sum_meanr = 0._RKIND
       sum_vol   = 0._RKIND
       do ic = 1,self%nchem
          if(self%fthom_w(ic)) then
             sum_hygro = sum_hygro + self%fhygros(ic)*self%chem_mr(i,k,ic)/self%fdens(ic)
             sum_meanr = sum_meanr + self%fmeanr(ic)*self%chem_mr(i,k,ic)/self%fdens(ic)
             sum_vol   = sum_vol + self%chem_mr(i,k,ic)/self%fdens(ic)
          endif
       enddo
       if(sum_vol .gt. 0._RKIND) then
          self%chem_ka(i,k) = sum_hygro/sum_vol
          self%chem_ra(i,k) = sum_meanr/sum_vol
       endif
    enddo
 enddo


 call mpas_log_write('--- end subroutine gocart2G_forMPASthom.')

 end subroutine gocart2G_forMPASthom

!==================================================================================================================
 subroutine gocart2G_fromMPASthom(self,state,time_lev)
!==================================================================================================================

!--- input arguments:
 class(chem_gocart2G),intent(in):: self
 integer,intent(in):: time_lev

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: state

!--- local variables and arrays:
 integer,pointer:: gocart2G_start,gocart2G_end
 integer:: its,ite,kts,kte
 integer:: i,ic,ig,k

 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_fromMPASthom:')


!---
 its = self%its
 ite = self%ite
 kts = self%kts
 kte = self%kte

 do ic = 1,self%nchem
    if(self%fthom_w(ic)) then
       self%chem_nc(:,:,ic) = self%chem_nwfa(:,:)*self%chem_nr(:,:,ic)
       self%chem_mr(:,:,ic) = self%chem_nc(:,:,ic)/self%fnum(ic)
    elseif(self%fthom_i(ic)) then
       self%chem_nc(:,:,ic) = self%chem_nifa(:,:)*self%chem_nr(:,:,ic)
       self%chem_mr(:,:,ic) = self%chem_nc(:,:,ic)/self%fnum(ic)
    endif
 enddo
!do i = its,ite
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r $r $r $r $r $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/ &
!                          (self%chem_mr(i,k,ic),ic=1,12)/))
!   enddo
!   call mpas_log_write(' ')
!enddo


 call mpas_pool_get_dimension(state,'gocart2G_start',gocart2G_start)
 call mpas_pool_get_dimension(state,'gocart2G_end'  ,gocart2G_end  )
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 ic = 0
 do ig = gocart2G_start,gocart2G_end
    ic = ic+1
    do i = its,ite
       do k = kts,kte
          scalars(ig,k,i) = self%chem_mr(i,k,ic)
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine gocart2G_fromMPASthom.')

 end subroutine gocart2G_fromMPASthom

!==================================================================================================================
 subroutine gocart2G_forMPASdyn(self,tend_physics)
!==================================================================================================================

!--- input arguments:
 class(chem_gocart2G),intent(in):: self

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: tend_physics

!local variables and pointers:
 character(len=StrKIND):: message

 integer,pointer:: bl_start,bl_end
 integer,pointer:: cu_start,cu_end
 integer:: bl,cu,i,ic,its,ite,k,kts,kte

 real(kind=RKIND),dimension(:,:,:),pointer:: bl_chemistry
 real(kind=RKIND),dimension(:,:,:),pointer:: cu_chemistry

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_forMPASdyn:')

 its = self%its
 ite = self%ite
 kts = self%kts
 kte = self%kte


 call mpas_pool_get_dimension(tend_physics,'gocart2G_bl_start',bl_start)
 call mpas_pool_get_dimension(tend_physics,'gocart2G_bl_end'  ,bl_end  )
 call mpas_log_write('--- gocart2G_bl_start = $i',intArgs=(/bl_start/))
 call mpas_log_write('--- gocart2G_bl_end   = $i',intArgs=(/bl_end/)  )

 call mpas_pool_get_array(tend_physics,'bl_chemistry',bl_chemistry)
 ic = 0
 do bl = bl_start,bl_end
    ic = ic+1
    do k = kts,kte
       do i = its,ite
          bl_chemistry(bl,k,i) = self%chemblten(i,k,ic)
       enddo
    enddo
 enddo
 call mpas_log_write('--- ic                = $i',intArgs=(/ic/))
 if(ic /= bl_end-bl_start+1) then
    message = '--- subroutine gocart2G_todynamics: bl_end-bl_start different than nb of gocart2G aerosol species'
    call mpas_log_write(message,messageType=MPAS_LOG_crit)
 endif


 call mpas_pool_get_dimension(tend_physics,'gocart2G_cu_start',cu_start)
 call mpas_pool_get_dimension(tend_physics,'gocart2G_cu_end'  ,cu_end  )
 call mpas_log_write('--- gocart2G_cu_start = $i',intArgs=(/cu_start/))
 call mpas_log_write('--- gocart2G_cu_end   = $i',intArgs=(/cu_end/)  )

 call mpas_pool_get_array(tend_physics,'cu_chemistry',cu_chemistry)
 ic = 0
 do cu = cu_start,cu_end
    ic = ic+1
    do k = kts,kte
       do i = its,ite
          cu_chemistry(cu,k,i) = self%chemcuten(i,k,ic)
       enddo
    enddo
 enddo
 call mpas_log_write('--- ic                = $i',intArgs=(/ic/))
 if(ic /= cu_end-cu_start+1) then
    message = '--- subroutine gocart2G_todynamics: cu_end-cu_start different than nb of gocart2G aerosol species'
    call mpas_log_write(message,messageType=MPAS_LOG_crit)
 endif

 call mpas_log_write('--- end subroutine gocart2G_forMPASdyn:')

 end subroutine gocart2G_forMPASdyn

!==================================================================================================================
 end module mpas_chemistry_gocart2G_forMPASphysics
!==================================================================================================================
