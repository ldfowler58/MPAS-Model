!=================================================================================================================
 module CA2G_oc_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: phobicTophilic

 use CA2G_oc_instance_CA,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                               fnum,rhFlag,pressure_lid_in_hPa
 use CA2G_oc_StateSpecs,only: CA2G_oc_State


 implicit none
 private
!public:: load_CA2G_oc_GridComp


 type :: ThreadWorkspace
    integer:: nPts = -1
    integer,dimension(:),allocatable:: pstart,pend

    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: CA2G_oc_GridComp
!   logical:: diurnal_bb                              ! diurnal biomass burning
!   integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

!   real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
!   real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
!   real(kind=RKIND):: fMonoterpenes = 0.0            ! fraction of monoterpene emissions -> aerosol
!   real(kind=RKIND):: fIsoprene = 0.0                ! fraction of isoprene emissions -> aerosol
!   real(kind=RKIND):: fHydrophobic                   ! initially hydrophobic portion
!   real(kind=RKIND):: ratPOM = 1.0                   ! ratio of POM to OC mass

!   !workspace for point emissions:
!      logical:: doing_point_emissions = .false.
!      character(len=255):: point_emissions_srcfilen  ! filename for pointwise emissions

    contains
       procedure:: load_GridComp => load_CA2G_oc_GridComp
       procedure:: run2_GridComp => run2_CA2G_oc_GridComp
 end type CA2G_oc_GridComp

 type wrap_
    type(CA2G_oc_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_CA2G_oc_GridComp(self,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(CA2G_oc_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_CA2G_oc_GridComp:')

!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization using parameters defined in CA2G_bc_instance_CA:
 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo

 call mpas_log_write('--- end subroutine load_CA2G_oc_GridCOMP:')
!call mpas_log_write(' ')

 end subroutine load_CA2G_oc_GridComp

!=================================================================================================================
 subroutine run2_CA2G_oc_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!=================================================================================================================
!--- input arguments:

 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_oc_GridComp),intent(inout):: self_params
 class(CA2G_oc_State),intent(inout):: self

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine run2_CA2G_oc_GridComp:')

!Ad Hoc transfer of hydrophobic to hydrophilic aerosols
!Following Chin's parameterization, the rate constant is
!k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)
!call phobicTophilic( &
!   intPtr_phobic, &
!   intPtr_philic, &
!   HYPHIL       , &
!   self%km      , &
!   self%cdt     , &
!   MAPL_GRAV    , &
!   delp         , &
!   __RC__)


 call mpas_log_write('--- end subroutine run2_CA2G_oc_GridComp:')

 end subroutine run2_CA2G_oc_GridComp

!=================================================================================================================
 end module CA2G_oc_GridCompMod
!=================================================================================================================
