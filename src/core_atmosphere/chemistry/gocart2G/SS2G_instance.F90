! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module SS2G_instance
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the resource file for sea salt parameters. it is the same as SS2G_instance_SS.rc in the GOCART-2G
!directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SS2G_GridComp.

!==================================================================================================================

!the location of those two files needs to be added if used:
!aerosol_radBands_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/opticsBands_SS.v3_3.RRTMG.nc
!aerosol_monochromatic_optics_file: ExtData/chemistry/AerosolOptics/v0.0.0/x/optics_SS.v3_3.nc

!number of bins:
 integer,parameter:: nbins = 5

real(kind=RKIND),dimension(nbins),parameter:: particle_radius_microns = (/0.079,0.316,1.119,2.818,7.772/)

real(kind=RKIND),dimension(nbins),parameter:: radius_lower = (/0.03,0.1,0.5,1.5,5.0/)

real(kind=RKIND),dimension(nbins),parameter:: radius_upper = (/0.1,0.5,1.5,5.0,10.0/)

real(kind=RKIND),dimension(nbins),parameter:: particle_density = (/2200.,2200.,2200.,2200.,2200./)

!scavenging efficiency per bin [km-1]:
real(kind=RKIND),dimension(nbins),parameter:: fscav = (/0.4,0.4,0.4,0.4,0.4/)

!emissions methods and scaling:
logical,parameter:: hoppelFlag  =  .false.    ! apply Hoppel correction (see Fan and Toon 2011)
logical,parameter:: weibullFlag =  .false.    ! apply Weibull distribution (see Fan and Toon 2011)
integer,parameter:: emission_scheme = 3       ! 1 for Gong 2003, 2 for ...
integer,parameter:: sstEmisFlag     = 2       ! apply a correction to emissions based on SST (see code)

real(kind=RKIND),dimension(6),parameter:: emission_scale_res = (/0.613,0.613,0.613,0.429,0.429,0.429/)

!method of apply relative humidity to particle radius:
integer,parameter:: rhFlag = 2                ! RH swelling of Seasalt (1 for Fitzgerald 1975,
                                              ! 2 for Gerber 1985 method)

!molecular weight of species [kg mole-1]:
real(kind=RKIND),dimension(nbins),parameter:: molecular_weight = (/0.058,0.058,0.058,0.058,0.058/)

!number of particles per kg mass:
real(kind=RKIND),dimension(nbins),parameter:: fnum = (/3.017e17,1.085e16,1.207e14,9.391e12,2.922e11/)

!number median radius [um]:
real(kind=RKIND),dimension(nbins),parameter:: particle_radius_number = (/0.066,0.176,0.885,2.061,6.901/)

real(kind=RKIND),parameter:: pressure_lid_in_hPa = 0.01

!==================================================================================================================
 end module SS2G_instance
!==================================================================================================================
