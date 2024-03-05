!=================================================================================================================
 module DU2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use GA_EnvironmentMod


 integer,parameter:: NHRES = 6
 type,extends(GA_Environment),public:: DU2G_GridComp
       character(len=:),allocatable :: emission_scheme  ! emission scheme selector

       logical:: maringFlag=.false.                     ! maring settling velocity correction
       integer:: day_save = -1
       integer:: clayFlag                               ! clay and silt term in K14

       real(kind=RKIND):: f_swc                         ! soil mosture scaling factor
       real(kind=RKIND):: f_scl                         ! clay content scaling factor
       real(kind=RKIND):: uts_gamma                     ! threshold friction velocity parameter 'gamma'
       real(kind=RKIND):: alpha                         ! FENGSHA scaling factor
       real(kind=RKIND):: gamma                         ! FENGSHA tuning exponent
       real(kind=RKIND):: kvhmax                        ! FENGSHA max. vertical/horizontal mass flux ratio [1]
       real(kind=RKIND):: Ch_DU                         ! dust emission tuning coefficient [kg s2 m-5].
       real(kind=RKIND),dimension(NHRES):: Ch_DU_res    ! resolutions used for Ch_DU

       real(kind=RKIND),dimension(:),allocatable:: rlow ! particle effective radius lower bound [um]
       real(kind=RKIND),dimension(:),allocatable:: rup  ! particle effective radius upper bound [um]
       real(kind=RKIND),dimension(:),allocatable:: sfra ! fraction of total source
       real(kind=RKIND),dimension(:),allocatable:: sdist! FENGSHA aerosol fractional size distribution [1]

       !workspace for point emissions:
       logical           :: doing_point_emissions = .false.
       character(len=255):: point_emissions_srcfilen    ! filename for pointwise emissions

    contains
       procedure:: load_forDU2G
 end type DU2G_GridComp


 contains


!=================================================================================================================
 subroutine load_forDU2G(self)
!=================================================================================================================

!--- inout arguments:
 class(DU2G_GridComp),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------


 end subroutine load_forDU2G

!=================================================================================================================
 endmodule DU2G_GridCompMod
!=================================================================================================================
