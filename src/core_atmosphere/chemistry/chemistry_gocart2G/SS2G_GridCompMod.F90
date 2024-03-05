!=================================================================================================================
 module SS2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use GA_EnvironmentMod

 integer,parameter:: NHRES = 6
 type,extends(GA_Environment),public:: SS2G_GridComp
       logical:: hoppelFlag     ! apply the Hoppel correction to emissions (Fan and Toon, 2011)
       logical:: weibullFlag    ! apply the Weibull distribution to wind speed for emissions (Fan and Toon, 2011)

       integer:: emission_scheme
       integer:: sstEmisFlag    ! choice of SST correction to emissions:
                                ! 0 - none; 1 - Jaegle et al. 2011; 2 - GEOS5

       real(kind=RKIND):: emission_scale                ! global scaling factor
       real(kind=RKIND):: emission_scale_res(NHRES)     ! global scaling factor
       real(kind=RKIND),dimension(:),allocatable:: rlow ! particle effective radius lower bound [um]
       real(kind=RKIND),dimension(:),allocatable:: rup  ! particle effective radius upper bound [um]
       real(kind=RKIND),dimension(:),allocatable:: rmed ! number median radius [um]
!      real(kind=RKIND),dimension(:,:),allocatable:: deep_lakes_mask

    contains
       procedure:: load_forSS2G
 end type SS2G_GridComp

 type wrap_
    type(SS2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_forSS2G(self)
!=================================================================================================================

!--- inout arguments:
 class(SS2G_GridComp),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------


 end subroutine load_forSS2G

!=================================================================================================================
 endmodule SS2G_GridCompMod
!=================================================================================================================
