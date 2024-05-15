!=================================================================================================================
 module CA2G_bc_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: phobicTophilic,Chem_Settling,DryDeposition

 use CA2G_bc_instance_CA,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                               fnum,rhFlag,pressure_lid_in_hPa,sigma,hydrophobic_fraction,            &
                               aircraft_fuel_emission_factor,aviation_vertical_layers
 use CA2G_bc_StateSpecs,only: CA2G_bc_State


 implicit none
 private


!--- constants (these parameters needs to be accessed from MPAS phys instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND


!--- types needed to define CA2G_bc:
 type :: ThreadWorkspace
    integer:: nPts = -1
    integer,dimension(:),allocatable:: pstart,pend

    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: CA2G_bc_GridComp
    logical:: diurnal_bb                              ! diurnal biomass burning
    integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

    real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
    real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
    real(kind=RKIND):: fMonoterpenes = 0.0            ! fraction of monoterpene emissions -> aerosol
    real(kind=RKIND):: fIsoprene = 0.0                ! fraction of isoprene emissions -> aerosol
    real(kind=RKIND):: fHydrophobic                   ! initially hydrophobic portion
    real(kind=RKIND):: ratPOM = 1.0                   ! ratio of POM to OC mass
    real(kind=RKIND),dimension(:),allocatable:: sigma ! sigma of lognormal number distribution

    !workspace for point emissions:
    logical:: doing_point_emissions = .false.
    character(len=255):: point_emissions_srcfilen  ! filename for pointwise emissions
    integer:: nPts = -1

    contains
       procedure:: load_GridComp =>load_CA2G_bc_GridComp
       procedure:: run2_GridComp => run2_CA2G_bc_GridComp
 end type CA2G_bc_GridComp

 type wrap_
    type(CA2G_bc_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_CA2G_bc_GridComp(self,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(CA2G_bc_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_CA2G_bc_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization using parameters defined in CA2G_bc_instance_CA:
 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 if(.not.allocated(self%sigma)) allocate(self%sigma(self%nbins))
 self%eAirCraftFuel = aircraft_fuel_emission_factor
 do n = 1,self%nbins
    self%sigma(n) = sigma(n)
 enddo
 do n = 1,4
    self%aviation_layers(n) = aviation_vertical_layers(n)
 enddo

 self%fHydrophobic = hydrophobic_fraction

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum,sigma:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n),self%sigma(n)/))
 enddo

 call mpas_log_write('--- end subroutine load_CA2G_bc_GridCOMP:')
!call mpas_log_write(' ')

 end subroutine load_CA2G_bc_GridComp

!=================================================================================================================
 subroutine run2_CA2G_bc_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!=================================================================================================================
!--- input arguments:

 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_bc_GridComp),intent(inout):: self_params
 class(CA2G_bc_State),intent(inout):: self


!--- local variables:
 integer:: i,j,k,ibin
 integer:: istat

 real(kind=RKIND),dimension(:,:),allocatable:: drydepf,dqa
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qca2G

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine run2_CA2G_bc_GridComp:')

!Ad Hoc transfer of hydrophobic to hydrophilic aerosols
!Following Chin's parameterization, the rate constant is
!k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)

 istat = 0
 call phobicTophilic( &
           aerosol_phobic        = self%bcPHOBIC   , &
           aerosol_philic        = self%bcPHILIC   , &
           aerosol_toHydrophilic = self%bcHYPHIL   , &
           km                    = self_params%km  , &
           cdt                   = self_params%cdt , &
           grav                  = grav            , &
           delp                  = self%delp       , &
           rc = istat                                &
                    )
!do j = jts,jte
!   do i = its,ite
!      do k = kts,kte
!         call mpas_log_write('$i $i $i $r $r $r $r',intArgs=(/j,i,k/),realArgs=(/self%delp(i,j,k), &
!                   self%bcPHOBIC(i,j,k),self%bcPHILIC(i,j,k),self%bcHYPHIL(i,j)/))
!      enddo
!      call mpas_log_write(' ')
!   enddo
!enddo


!--- CA2G_bc settling:
!if(.not.allocated(qca2G)) allocate(qca2G(its:ite,jts:jte,kts:kte,self_params%nbins))
!do j = jts,jte
!   do i = its,ite
!      do k = kts,kte
!         qca2G(i,j,k,1) = self%bcPHOBIC(i,j,k)
!         qca2G(i,j,k,2) = self%bcPHILIC(i,j,k)
!      enddo
!   enddo
!enddo

!do ibin = 1, self_params%nbins
!   !if radius == 0, then we're dealing with a gas which has no settling losses:
!   if(self_params%radius(ibin) == 0.0) then
!      if(associated(self%bcSD)) self%bcSD(:,:,ibin) = 0.0
!      cycle
!   endif

!   call mpas_log_write('--- km        = $i',intArgs=(/self_params%km/))
!   call mpas_log_write('--- klid      = $i',intArgs=(/self_params%klid/))
!   call mpas_log_write('--- ibin      = $i',intArgs=(/ibin/))
!   call mpas_log_write('--- rhFlag    = $i',intArgs=(/self_params%rhFlag/))
!   call mpas_log_write('--- cdt       = $r',realArgs=(/self_params%cdt/))
!   call mpas_log_write('--- grav      = $r',realArgs=(/grav/))
!   call mpas_log_write('--- radiusInp = $r',realArgs=(/self_params%radius(ibin)/))
!   call mpas_log_write('--- rhoInp    = $r',realArgs=(/self_params%rhop(ibin)/))
!   call mpas_log_write(' ')
!   istat = 0
!   call Chem_Settling( &
!             km        = self_params%km                 , &
!             klid      = self_params%klid               , &
!             bin       = ibin                           , &
!             flag      = self_params%rhFlag             , &
!             cdt       = self_params%cdt                , &
!             grav      = grav                           , &
!             radiusInp = self_params%radius(ibin)*1.e-6 , &
!             rhopInp   = self_params%rhop(ibin)         , &
!             int_qa    = qca2G(:,:,:,ibin)              , &
!             tmpu      = self%t                         , &
!             rhoa      = self%airdens                   , &
!             rh        = self%rh2                       , &
!             hghte     = self%zle                       , &
!             ple       = self%ple                       , &
!             delz      = self%delz                      , &
!             delp      = self%delp                      , &
!             fluxout   = self%bcSD                      , &
!             rc        = istat                            &
!                     )
!enddo
!if(allocated(qca2G)) deallocate(qca2G)


!--- CA2G_bc deposition:
!if(.not.allocated(dqa)    ) allocate(dqa(its:ite,jts:jte)    )
!if(.not.allocated(drydepf)) allocate(drydepf(its:ite,jts:jte))
!drydepf = 0.

!istat = 0
!call DryDeposition( &
!          km         = self_params%km , &
!          tmpu       = self%t         , &
!          rhoa       = self%airdens   , &
!          hghte      = self%zle       , &
!          oro        = self%lwi       , &
!          ustar      = self%ustar     , &
!          pblh       = self%zpbl      , &
!          shflux     = self%sh        , &
!          von_karman = karman         , &
!          cpd        = cpd            , &
!          grav       = grav           , &
!          z0h        = self%z0h       , &
!          drydepf    = drydepf        , &
!          rc         = istat            &
!                  )

!do n = 1, self%nbins
!   dqa = 0.
!   dqa = max(0.0, int_ptr(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
!   int_ptr(:,:,self%km) = int_ptr(:,:,self%km) - dqa
!   if (associated(DP)) then
!      DP(:,:,n) = dqa * delp(:,:,self%km) / MAPL_GRAV / self%cdt
!   end if
!end do
!if(allocated(dqa)    ) deallocate(dqa    )
!if(allocated(drydepf)) deallocate(drydepf)


 call mpas_log_write('--- end subroutine run2_CA2G_bc_GridComp:')

 end subroutine run2_CA2G_bc_GridComp

!=================================================================================================================
 end module CA2G_bc_GridCompMod
!=================================================================================================================
