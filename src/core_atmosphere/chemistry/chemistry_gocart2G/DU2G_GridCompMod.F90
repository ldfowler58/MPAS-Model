!=================================================================================================================
 module DU2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process

 use DU2G_instance_DU,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                            fnum,rhFlag,pressure_lid_in_hPa
 use DU2G_StateSpecs,only: DU2G_State


 implicit none
!private
!public:: load_DU2G_GridComp


 integer,parameter:: NHRES = 6
 type,extends(GA_Environment),public:: DU2G_GridComp
!      character(len=:),allocatable :: emission_scheme  ! emission scheme selector

!      logical:: maringFlag=.false.                     ! maring settling velocity correction
!      integer:: day_save = -1
!      integer:: clayFlag                               ! clay and silt term in K14

!      real(kind=RKIND):: f_swc                         ! soil mosture scaling factor
!      real(kind=RKIND):: f_scl                         ! clay content scaling factor
!      real(kind=RKIND):: uts_gamma                     ! threshold friction velocity parameter 'gamma'
!      real(kind=RKIND):: alpha                         ! FENGSHA scaling factor
!      real(kind=RKIND):: gamma                         ! FENGSHA tuning exponent
!      real(kind=RKIND):: kvhmax                        ! FENGSHA max. vertical/horizontal mass flux ratio [1]
!      real(kind=RKIND):: Ch_DU                         ! dust emission tuning coefficient [kg s2 m-5].
!      real(kind=RKIND),dimension(NHRES):: Ch_DU_res    ! resolutions used for Ch_DU

!      real(kind=RKIND),dimension(:),allocatable:: rlow ! particle effective radius lower bound [um]
!      real(kind=RKIND),dimension(:),allocatable:: rup  ! particle effective radius upper bound [um]
!      real(kind=RKIND),dimension(:),allocatable:: sfra ! fraction of total source
!      real(kind=RKIND),dimension(:),allocatable:: sdist! FENGSHA aerosol fractional size distribution [1]

!      !workspace for point emissions:
!      logical           :: doing_point_emissions = .false.
!      character(len=255):: point_emissions_srcfilen    ! filename for pointwise emissions

    contains
       procedure:: load_GridComp => load_DU2G_GridComp
 end type DU2G_GridComp


 contains


!=================================================================================================================
 subroutine load_DU2G_GridComp(self)
!=================================================================================================================

!--- inout arguments:
 class(DU2G_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_DU2G_GridComp:')

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo

 call mpas_log_write('--- end subroutine load_DU2G_GridCOMP:')
!call mpas_log_write(' ')

 end subroutine load_DU2G_GridComp

!=================================================================================================================
 endmodule DU2G_GridCompMod
!=================================================================================================================
