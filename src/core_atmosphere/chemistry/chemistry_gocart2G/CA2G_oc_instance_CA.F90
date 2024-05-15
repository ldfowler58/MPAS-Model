!=================================================================================================================
 module CA2G_oc_instance_CA
 use mpas_kind_types,only: RKIND,StrKIND

 public
 save

!this module is the resource file for organic carbon parameters. it is the same as CA2G_instance_CA.oc.rc in the
!GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/CA2G_GridComp.

!=================================================================================================================

!the location of those two files needs to be added if used:
!aerosol_radBands_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_OC.v1_3.RRTMG.nc
!aerosol_monochromatic_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_OC.v1_3.nc

!aircraft emission factor: convert input unit to kg C:
 real(kind=RKIND),parameter:: aircraft_fuel_emission_factor = 1.0000

!heights [m] of LTO, CDS and CRS aviation emissions layers:
 real(kind=RKIND),dimension(4),parameter:: aviation_vertical_layers = (/0.0,100.0,9.0e3,10.0e3/)

!number of bins:
 integer,parameter:: nbins = 2

!fraction of biogenic VOCs emissions for SOA production:
 real(kind=RKIND),parameter:: monoterpenes_emission_fraction = 0.05
 real(kind=RKIND),parameter:: isoprene_emission_fraction = 0.03

!ratio of POM/OC -> convert source masses from carbon to POM:
 real(kind=RKIND),parameter:: pom_ca_ratio = 1.8

!particle radius:
 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_microns = (/0.35,0.35/)

 real(kind=RKIND),parameter:: rhFlag = 0

!initially hydrophobic portion:
 real(kind=RKIND),parameter:: hydrophobic_fraction = 0.5

!scavenging efficiency per bin [km-1] (NOT USED UNLESS RAS IS CALLED):
 real(kind=RKIND),dimension(nbins),parameter:: fscav = (/0.0,0.4/)

!dry particle density [kg m-3]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_density = (/1800,1800/)

!molecular weight of species [kg mole-1]:
 real(kind=RKIND),dimension(nbins),parameter:: molecular_weight = (/0.18,0.18/)

!number of particles per kg mass:
 real(kind=RKIND),dimension(nbins),parameter:: fnum = (/9.76e17,9.76e17/)

!sigma of lognormal number distribution:
 real(kind=RKIND),dimension(nbins),parameter:: sigma = (/2.20,2.20/)

 real(kind=RKIND),parameter:: pressure_lid_in_hPa = 0.01

!character(kind=StrKIND),parameter:: point_emissions_srcfilen = /dev/null

!=================================================================================================================
 end module CA2G_oc_instance_CA
!=================================================================================================================
