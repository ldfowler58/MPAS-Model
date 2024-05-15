!=================================================================================================================
 module GOCART2G_GridComp
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the resource file for gocart2G parameters. it is the same as GOCART2G_GridComp.rc in the GOCART-2G
!directory ./GOCART-2G/ESMF/GOCART2G_GridComp. it only contains the definition related to wavelengths and not to
!"active" or "passive" instances.

!=================================================================================================================

!set optics parameters:
 real(kind=RKIND),dimension(4),parameter:: &
    aerosol_monochromatic_optics_wavelength_in_nm_from_LUT = (/470.,550.,670.,870./)

 real(kind=RKIND),dimension(1),parameter:: wavelengths_for_profile_aop_in_nm = 550.
 real(kind=RKIND),dimension(1),parameter:: wavelengths_for_vertically_integrated_aop_in_nm = 550.

!=================================================================================================================
 end module GOCART2G_GridComp
!=================================================================================================================

