!=================================================================================================================
 module CA2G_bc_GridCompMod
 use mpas_kind_types,only: RKIND
 use GA_EnvironmentMod


 type,extends(GA_Environment),public:: CA2G_bc_GridComp
    logical:: diurnal_bb                              ! diurnal biomass burning
    integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

    real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
    real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
    real(kind=RKIND):: fMonoterpenes = 0.0            ! fraction of monoterpene emissions -> aerosol
    real(kind=RKIND):: fIsoprene = 0.0                ! fraction of isoprene emissions -> aerosol
    real(kind=RKIND):: fHydrophobic                   ! initially hydrophobic portion
    real(kind=RKIND):: ratPOM = 1.0                   ! ratio of POM to OC mass

    !workspace for point emissions:
       logical:: doing_point_emissions = .false.
       character(len=255):: point_emissions_srcfilen  ! filename for pointwise emissions

    contains
       procedure:: load_forCA2G_bc
 end type CA2G_bc_GridComp


 contains


!=================================================================================================================
 subroutine load_forCA2G_bc(self)
!=================================================================================================================

!--- inout arguments:
 class(CA2G_bc_GridComp),intent(inout) :: self

!-----------------------------------------------------------------------------------------------------------------


 end subroutine load_forCA2G_bc

!=================================================================================================================
 end module CA2G_bc_GridCompMod
!=================================================================================================================
