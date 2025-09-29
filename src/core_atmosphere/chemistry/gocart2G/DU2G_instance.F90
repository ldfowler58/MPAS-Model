! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module DU2G_instance
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the resource file for dust parameters. it is the same as DU2G_instance_DU.rc in the GOCART-2G
!directory ./GOCART-2G/ESMF/GOCART2G_GridComp/DU2G_GridComp.

!==================================================================================================================

!the location of those two files needs to be added if used:
!aerosol_radBands_optics_file:      ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_DU.v15_3.RRTMG.nc
!aerosol_monochromatic_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_DU.v15_3.nc

!number of bins:
 integer,parameter:: nbins = 5

 real(kind=RKIND),dimension(nbins),parameter:: particle_radius_microns = (/0.73,1.4,2.4,4.5,8.0/)

 real(kind=RKIND),dimension(nbins),parameter:: radius_lower = (/0.1,1.0,1.8,3.0,6.0/)

 real(kind=RKIND),dimension(nbins),parameter:: radius_upper = (/1.0,1.8,3.0,6.0,10.0/)

!units [kg/m-3]:
 real(kind=RKIND),dimension(nbins),parameter:: particle_density = (/2500.,2650.,2650.,2650.,2650./)

!Ginoux emission scheme:
!----------------------
 character(len=*),parameter:: emission_scheme = 'ginoux'

!Ch_DU: 0.2 0.2 0.07 0.07 0.07 0.056 #original values for (a,b,c,d,e,f)
!real(kind=RKIND),dimension(nbins),parameter:: source_fraction = (/.0435465,0.106903,0.220117,0.484606,0.144828/)

 real(kind=RKIND),dimension(6),parameter:: Ch_DU = (/0.3,0.3,0.11,0.11,0.11,0.088/)
 real(kind=RKIND),dimension(nbins),parameter:: source_fraction = (/0.1,0.25,0.25,0.25,0.25/)

!LDF: added parameters from WRF model (see ./phys/module_data_gocart_dust):
 integer,dimension(nbins),parameter:: ipoint = (/3,2,2,2,2/)

!K14 emission scheme:
!-------------------
!emission_scheme: k14
!source_fraction:  0.043 0.106 0.219 0.485 0.144
!Ch_DU: 0.02 0.02 0.02 0.0161 0.015 0.015 # resolution dependent tuning constant for emissions (a,b,c,d,e,f)

!threshold friction velocity parameter 'gamma':
 real(kind=RKIND),parameter:: uts_gamma = 1.65e-4

!formulation of the clay and silt factor in K14 that modulates the strength of the dust emissions:
 real(kind=RKIND),parameter:: clayFlag = 1 ! 0 - original K14, 1 - I&K2017, 2 - I&K2017

!soil mosture scaling factor:
 real(kind=RKIND),parameter:: soil_moisture_factor = 0.8

!clay fraction scaling factor:
 real(kind=RKIND),parameter:: soil_clay_factor = 1.0

!scavenging efficiency per bin [km-1]:
 real(kind=RKIND),dimension(nbins),parameter:: fscav = (/0.2,0.2,0.2,0.2,0.2/)

!molecular weight of species [kg mole-1]:
 real(kind=RKIND),dimension(nbins),parameter:: molecular_weight = (/0.1,0.1,0.1,0.1,0.1/)

!number of particles per kg mass:
 real(kind=RKIND),dimension(nbins),parameter:: fnum = (/2.45e14,3.28e13,6.52e12,9.89e11,1.76e11/)

!number median radius [um]:
real(kind=RKIND),dimension(nbins),parameter:: particle_radius_number = (/0.123,0.421,0.722,1.354,1.354/)

!sigma of lognormal number distribution:
 real(kind=RKIND),dimension(nbins),parameter:: sigma = (/2.0,2.0,2.0,2.0,2.0/)

 integer,parameter:: rhFlag = 0

!maring settling velocity correction:
 logical,parameter:: maringFlag = .true.

 real(kind=RKIND),parameter:: pressure_lid_in_hPa = 0.01

!FENGSHA settings:
 real(kind=RKIND),parameter:: alpha = 0.3
 real(kind=RKIND),parameter:: gamma = 1.3
 real(kind=RKIND),parameter:: vertical_to_horizontal_flux_ratio_limit = 2.e-04

!==================================================================================================================
 end module DU2G_instance
!==================================================================================================================
