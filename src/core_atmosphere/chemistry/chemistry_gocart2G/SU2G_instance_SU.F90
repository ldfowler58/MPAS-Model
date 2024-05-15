!=================================================================================================================
 module SU2G_instance_SU
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the resource file for sulfur parameters. it is the same as SU2G_instance_SU.rc in the GOCART-2G
!directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SU2G_GridComp.

!=================================================================================================================

!the location of those two files needs to be added if used:
!aerosol_radBands_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_SU.v1_3.RRTMG.nc
!aerosol_monochromatic_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_SU.v1_3.nc

!the location of this file needs to be added if needed:
!Volcanic pointwise sources:
!volcano_srcfilen: ExtData/chemistry/CARN/v202106/sfc/so2_volcanic_emissions_CARN_v202106.degassing_only.rc

!number of bins:
 integer,parameter:: nbins = 4

!heights [m] of LTO, CDS and CRS aviation emissions layers:
 real(kind=RKIND),dimension(4),parameter:: aviation_vertical_layers = (/0.0,100.0,9.0e3,10.0e3/)

!fraction of anthropogenic emissions that are SO4:
 real(kind=RKIND),parameter:: so4_anthropogenic_fraction = 0.03

!aircraft emission factor: convert input unit to kg SO2:
 real(kind=RKIND),parameter:: aircraft_fuel_emission_factor = 1.0000

!scavenging efficiency per bin [km-1] (NOT USED UNLESS RAS IS CALLED):
 real(kind=RKIND),dimension(nbins),parameter:: fscav = (/0.0,0.0,0.4,0.0/)

!dry particle radius [um], used for settling:
 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_microns =  (/0.0,0.0,0.35,0.0/)

!type of settling to use (see Chem_SettlingMod):
 real(kind=RKIND),parameter:: rhFlag = 4

!dry particle density [kg m-3]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_density =(/-1.,-1.,1700.,-1./)

 real(kind=RKIND),parameter:: pressure_lid_in_hPa = 0.01

!molecular weight of species [kg mole-1]:
 real(kind=RKIND),dimension(nbins),parameter:: molecular_weight = (/-1.,-1.,0.132,-1./)

!number of particles per kg mass:
 real(kind=RKIND),dimension(nbins),parameter:: fnum = (/-1.,-1.,9.01e16,-1./)

!number median radius [um]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_number = (/-1.,-1.,0.0695,-1./)

!sigma of lognormal number distribution:
 real(kind=RKIND),dimension(nbins),parameter:: sigma = (/-1.,-1.,2.03,-1./)

!=================================================================================================================
 end module SU2G_instance_SU
!=================================================================================================================

