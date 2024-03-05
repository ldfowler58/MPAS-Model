!=================================================================================================================
 module SU2G_GridCompMod
 use mpas_kind_types,only: RKIND,StrKIND
 use GA_EnvironmentMod


 type:: ThreadWorkspace
    logical:: firstRun = .true.
    logical:: recycle_H2O2 = .false.

    integer:: nVolc = 0
    integer:: nPts = -1
    integer:: nymd_oxidants = -1 ! update the oxidant files?
    integer:: nymd_last = -1     ! previous nymd. updated daily.
    integer,dimension(:),allocatable:: pstart,pend
    integer,dimension(:),allocatable:: vStart,vEnd

    real(kind=RKIND),dimension(:),allocatable:: vLat,vLon,vSO2,vElev,vCloud
    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: SU2G_GridComp
    logical:: diurnal_bb                              ! diurnal biomass burning
    integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

    real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
    real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
    real(kind=RKIND):: fSO4anth                       ! fraction of anthropogenic emissions that are SO4
    real(kind=RKIND),dimension(:),allocatable:: sigma ! sigma of lognormal number distribution

    !special handling for volcanic emissions:
    character(len=strKIND):: volcano_srcfilen

    !workspace for point emissions:
    character(len=StrKIND):: point_emissions_srcfilen ! filename for pointwise emissions
    logical:: doing_point_emissions = .false.
    type(ThreadWorkspace),dimension(:),allocatable:: workspaces

    contains
       procedure:: load_forSU2G
 end type SU2G_GridComp

 type wrap_
    type(SU2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_forSU2G(self)
!=================================================================================================================

!--- inout arguments:
 class(SU2G_GridComp),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------


 end subroutine load_forSU2G

!=================================================================================================================
 end module SU2G_GridCompMod
!=================================================================================================================

