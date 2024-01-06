!=================================================================================================================
 module module_data_gocart_dust
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND

 implicit none
 save

 integer,parameter:: ndust=5,ndcls=3,ndsrc=1,maxstypes=100

 real(kind=RKIND),dimension(maxstypes):: porosity
 
 real(kind=RKIND),parameter:: ch_dust  = 1.0E-9 !empirical constant used in dust emission calculation. default is
                                                !1 ug s2 m-5 == 1.e-09 kg s2 m-5.

 real(kind=RKIND),parameter:: dyn_visc = 1.5E-5
!dust density (clay = 2500 kg m-3; Silt = 2650 kg m-3):
!real(kind=RKIND),dimension(ndust),parameter:: den_dust=(/2500.,2650.,2650.,2650.,2650./)
 real(kind=RKIND),dimension(ndust):: den_dust=(/2500.,2650.,2650.,2650.,2650./)

!dust effective radii (m):
 real(kind=RKIND),dimension(ndust),parameter:: reff_dust=(/0.73E-6,1.4E-6,2.4E-6,4.5E-6,8.0E-6/)

!USDA  defines dust particles with radius below 1 um as clay (class 3) and those with
!radius between 1 um and 25 um as silt (class 2). Array pointing to the class of particles
 integer,dimension(ndust),parameter:: ipoint=(/3,2,2,2,2/)

! fraction of each class available for erosion. fraction of clay is assumed to
! be 10% of the total mass of the emitted silt. 
 real(kind=RKIND),dimension(ndust),parameter:: frac_s=(/0.1,0.25,0.25,0.25,0.25/)

 real(kind=RKIND),dimension(ndust),parameter:: distr_dust=(/1.074e-1,1.012e-1,2.078e-1,4.817e-1,1.019e-1/)
!real(kind=RKIND),dimension(ndust),parameter:: lo_dust=(/0.1e-6,1.0e-6,1.8e-6,3.0e-6,6.0e-6/)
!real(kind=RKIND),dimension(ndust),parameter:: up_dust=(/1.0e-6,1.8e-6,3.0e-6,6.0e-6,10.0e-6/)

!bins radii ranges:
 real(kind=RKIND),dimension(5),parameter:: ra_dust(5)=(/1.e-1,1.0e0,1.8e0,3.0e0,6.0e0/)
 real(kind=RKIND),dimension(5),parameter:: rb_dust(5)=(/1.0e0,1.8e0,3.0e0,6.0e0,10.0e0/)

!=================================================================================================================
 end module module_data_gocart_dust
!=================================================================================================================


