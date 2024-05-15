!=================================================================================================================
 module SS2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process

 use SS2G_instance_SS,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                            fnum,rhFlag,pressure_lid_in_hPa
 use SS2G_StateSpecs,only: SS2G_State


 implicit none
!private
!public:: load_SS2G_GridComp


 integer,parameter:: NHRES = 6
 type,extends(GA_Environment),public:: SS2G_GridComp
!      logical:: hoppelFlag     ! apply the Hoppel correction to emissions (Fan and Toon, 2011)
!      logical:: weibullFlag    ! apply the Weibull distribution to wind speed for emissions (Fan and Toon, 2011)

!      integer:: emission_scheme
!      integer:: sstEmisFlag    ! choice of SST correction to emissions:
!                               ! 0 - none; 1 - Jaegle et al. 2011; 2 - GEOS5

!      real(kind=RKIND):: emission_scale                             ! global scaling factor
!      real(kind=RKIND):: emission_scale_res(NHRES)                  ! global scaling factor
!      real(kind=RKIND),dimension(:),allocatable:: rlow              ! particle effective radius lower bound [um]
!      real(kind=RKIND),dimension(:),allocatable:: rup               ! particle effective radius upper bound [um]
!      real(kind=RKIND),dimension(:),allocatable:: rmed              ! number median radius [um]
!      real(kind=RKIND),dimension(:,:),allocatable:: deep_lakes_mask ! mask for deep lakes

    contains
       procedure:: load_GridComp => load_SS2G_GridComp
 end type SS2G_GridComp

 type wrap_
    type(SS2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_SS2G_GridComp(self)
!=================================================================================================================

!--- inout arguments:
 class(SS2G_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_SS2G_GridComp:')

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo

 call mpas_log_write('--- end subroutine load_SS2G_GridCOMP:')
!call mpas_log_write(' ')

 end subroutine load_SS2G_GridComp

!=================================================================================================================
 endmodule SS2G_GridCompMod
!=================================================================================================================
