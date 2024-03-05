!=================================================================================================================
 module GA_EnvironmentMod
 use mpas_kind_types,only: RKIND

 use GOCART2G_MieMod

 implicit none
 private
 public :: GA_Environment


 type :: GA_Environment
    type(GOCART2G_Mie):: rad_Mie,diag_Mie

!   logical:: scav_byColdCloud                                      ! new flag example

    integer:: rhFlag
    integer:: nbins
    integer:: km                                                    ! vertical grid dimension
    integer:: instance                                              ! data or computational instance
    integer:: klid                                                  ! vertical index of pressure lid

    real(kind=RKIND):: cdt                                          ! chemistry timestep (secs)
    real(kind=RKIND):: plid                                         ! pressure lid [hPa]
    real(kind=RKIND),dimension(:),allocatable:: radius              ! particle effective radius [um]
    real(kind=RKIND),dimension(:),allocatable:: rhop                ! soil class density [kg m-3]
    real(kind=RKIND),dimension(:),allocatable:: fscav               ! scavenging efficiency
    real(kind=RKIND),dimension(:),allocatable:: molwght             ! molecular weight
    real(kind=RKIND),dimension(:),allocatable:: fnum                ! number of particles per kg mass
    real(kind=RKIND),dimension(:),allocatable:: wavelengths_profile ! wavelengths for profile aop [nm]
    real(kind=RKIND),dimension(:),allocatable:: wavelengths_vertint ! wavelengths for vertically integrated aop [nm]

    contains
       procedure:: load_from_config
    end type GA_Environment


 contains


!=================================================================================================================
 subroutine load_from_config(self,nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                             rhFlag,pressure_lid_in_hPa)
!=================================================================================================================

!--- inut arguments:
 integer,intent(in):: nbins

 real(kind=RKIND),intent(in):: pressure_lid_in_hPa
 real(kind=RKIND),intent(in):: rhFlag

 real(kind=RKIND),intent(in),dimension(nbins):: particle_radius_microns
 real(kind=RKIND),intent(in),dimension(nbins):: particle_density
 real(kind=RKIND),intent(in),dimension(nbins):: molecular_weight
 real(kind=RKIND),intent(in),dimension(nbins):: fscav
 real(kind=RKIND),intent(in),dimension(nbins):: fnum


!--- inout arguments:
 class(GA_Environment),intent(inout) :: self


!--- local variables:
 integer:: n
 integer:: n_wavelengths_profile
 integer:: n_wavelengths_vertint

!-----------------------------------------------------------------------------------------------------------------

 self%nbins = nbins

!--- allocate arrays contained in GA_Environment:
 if(.not.allocated(self%radius) ) allocate(self%radius(nbins) )
 if(.not.allocated(self%rhop)   ) allocate(self%rhop(nbins)   )
 if(.not.allocated(self%fscav)  ) allocate(self%fscav(nbins)  )
 if(.not.allocated(self%molwght)) allocate(self%molwght(nbins))
 if(.not.allocated(self%fnum)   ) allocate(self%fnum(nbins)   )

!if(.not.allocated(self%wavelengths_profile)) allocate(self%wavelengths_profile(n_wavelengths_profile))
!if(.not.allocated(self%wavelengths_vertint)) allocate(self%wavelengths_vertint(n_wavelengths_profile))


!--- initialize arrays in GA_Environment:
 do n = 1, nbins
    self%radius(n)  = particle_radius_microns(n)
    self%rhop(n)    = particle_density(n)
    self%fscav(n)   = fscav(n)
    self%molwght(n) = molecular_weight(n)
    self%fnum(n)    = fnum(n)
 enddo

 self%rhFlag = rhFlag
 self%plid   = pressure_lid_in_hPa


!--- DO NOT KNOW YET HOW TO INITIALIZE THOSE ARRAYS:
!n_wavelengths_profile = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_profile_aop_in_nm:', __RC__)
!n_wavelengths_vertint = ESMF_ConfigGetLen (universal_cfg, label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)
!do n = 1, n_wavelengths_profile
!   self%wavelengths_profile(n) =
!   self%wavelengths_vertint(n) = 
!enddo
!call ESMF_ConfigGetAttribute(universal_cfg,self%wavelengths_profile, &
!                             label='wavelengths_for_profile_aop_in_nm:', __RC__)
!call ESMF_ConfigGetAttribute(universal_cfg,self%wavelengths_vertint, &
!                             label='wavelengths_for_vertically_integrated_aop_in_nm:', __RC__)

 end subroutine load_from_config

!=================================================================================================================
 end module GA_EnvironmentMod
!=================================================================================================================
