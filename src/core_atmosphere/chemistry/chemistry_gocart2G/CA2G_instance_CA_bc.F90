!=================================================================================================================
 module CA2G_instance_CA_bc
 use mpas_kind_types,only: RKIND,StrKIND

 public
 save

!this module is the resource file for black carbon parameters. it is the same as CA2G_instance_CA.bc.rc in the
!GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!=================================================================================================================

!the location of those two files needs to be added if used:
!aerosol_radBands_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_BC.v1_3.RRTMG.nc
!aerosol_monochromatic_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_BC.v1_3.nc

!aircraft emission factor: convert input unit to kg C:
 real(kind=RKIND),parameter:: aircraft_fuel_emission_factor = 1.0000

!heights [m] of LTO, CDS and CRS aviation emissions layers:
 real(kind=RKIND),dimension(4),parameter:: aviation_vertical_layers = (/0.0,100.0,9.0e3,10.0e3/)

!number of bins:
 integer,parameter:: nbins = 2

!initially hydrophobic portion:
 real(kind=RKIND),parameter:: hydrophobic_fraction = 0.8

!scavenging efficiency per bin [km-1] (NOT USED UNLESS RAS IS CALLED):
 real(kind=RKIND),dimension(nbins),parameter:: fscav = (/0.0,0.4/)

!dry particle density [kg m-3]
 real(kind=RKIND),dimension(nbins),parameter:: particle_density = (/1800,1800/)

!molecular weight of species [kg mole-1]:
 real(kind=RKIND),dimension(nbins),parameter:: molecular_weight = (/0.18,0.18/)

!number of particles per kg mass:
 real(kind=RKIND),dimension(nbins),parameter:: fnum = (/1.50e19,1.50e19/)

!number median radius [um]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_microns =  (/0.35,0.35/)

 real(kind=RKIND),parameter:: rhFlag = 0

!sigma of lognormal number distribution:
 real(kind=RKIND),dimension(nbins),parameter:: sigma = (/2.0,2.0/)

 real(kind=RKIND),parameter:: pressure_lid_in_hPa = 0.01

!character(kind=StrKIND),parameter:: point_emissions_srcfilen = /dev/null

!=================================================================================================================
 end module CA2G_instance_CA_bc
!=================================================================================================================
