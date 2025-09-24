!==================================================================================================================
 module mp_thompson

 use mpas_log
 use mpas_kind_types
 use mpas_atmphys_functions, only: gammp,wgamma,rslf,rsif
 use mpas_atmphys_utilities
 use mpas_io_units, only : mpas_new_unit, mpas_release_unit
 use mp_radar

 implicit none


 logical,parameter,private:: iiwarm = .false.
 logical,public:: is_aerosol_aware = .false.

 logical,parameter,private:: dustyIce = .true.
 logical,parameter,private:: homogIce = .true.

 integer,parameter,private:: ifdry = 0
 real(kind=RKIND),parameter,private:: t_0 = 273.15
 real(kind=RKIND),parameter,private:: pi = 3.1415926536

!..Densities of rain, snow, graupel, and cloud ice.
 real(kind=RKIND),parameter,private:: rho_w = 1000.0
 real(kind=RKIND),parameter,private:: rho_s = 100.0
 real(kind=RKIND),parameter,private:: rho_g = 500.0
 real(kind=RKIND),parameter,private:: rho_i = 890.0

!..Prescribed number of cloud droplets.  Set according to known data or
!.. roughly 100 per cc (100.E6 m^-3) for Maritime cases and
!.. 300 per cc (300.E6 m^-3) for Continental.  Gamma shape parameter,
!.. mu_c, calculated based on Nt_c is important in autoconversion
!.. scheme.  In 2-moment cloud water, Nt_c represents a maximum of
!.. droplet concentration and nu_c is also variable depending on local
!.. droplet number concentration.
!.. MPAS: Nt_c is initialized to 100.E6 over oceans and 300.E6 over land as
!   a function of landmask in subroutine init_thompson_clouddroplets_forMPAS.
!real(kind=RKIND),parameter,private:: Nt_c = 100.E6
 real(kind=RKIND),parameter,private:: Nt_c_max = 1999.E6
 real(kind=RKIND),private:: Nt_c

!..Declaration of constants for assumed CCN/IN aerosols when none in
!.. the input data.  Look inside the init routine for modifications
!.. due to surface land-sea points or vegetation characteristics.
!.. MPAS: naIN0, naIN1, naCCN0, and naCCN1 are used in init_thompson_aerosols_forMPAS
!.. for initialization of nwfa. and nifa.
 real(kind=RKIND),parameter, PUBLIC:: naIN0 = 1.5E6
 real(kind=RKIND),parameter, PUBLIC:: naIN1 = 0.5E6
 real(kind=RKIND),parameter, PUBLIC:: naCCN0 = 300.0E6
 real(kind=RKIND),parameter, PUBLIC:: naCCN1 = 50.0E6

!..Generalized gamma distributions for rain, graupel and cloud ice.
!.. N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
 real(kind=RKIND),parameter,private:: mu_r = 0.0
 real(kind=RKIND),parameter,private:: mu_g = 0.0
 real(kind=RKIND),parameter,private:: mu_i = 0.0
 real(kind=RKIND),private:: mu_c

!..Sum of two gamma distrib for snow (Field et al. 2005).
!.. N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
!..    + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
!.. M2 and M3 are the (bm_s)th and (bm_s+1)th moments respectively
!.. calculated as function of ice water content and temperature.
 real(kind=RKIND),parameter,private:: mu_s = 0.6357
 real(kind=RKIND),parameter,private:: Kap0 = 490.6
 real(kind=RKIND),parameter,private:: Kap1 = 17.46
 real(kind=RKIND),parameter,private:: Lam0 = 20.78
 real(kind=RKIND),parameter,private:: Lam1 = 3.29

!..Y-intercept parameter for graupel is not constant and depends on
!.. mixing ratio.  Also, when mu_g is non-zero, these become equiv
!.. y-intercept for an exponential distrib and proper values are
!.. computed based on same mixing ratio and total number concentration.
 real(kind=RKIND),parameter,private:: gonv_min = 1.E4
 real(kind=RKIND),parameter,private:: gonv_max = 3.E6

!..Mass power law relations:  mass = am*D**bm
!.. Snow from Field et al. (2005), others assume spherical form.
 real(kind=RKIND),parameter,private:: am_r = PI*rho_w/6.0
 real(kind=RKIND),parameter,private:: bm_r = 3.0
 real(kind=RKIND),parameter,private:: am_s = 0.069
 real(kind=RKIND),parameter,private:: bm_s = 2.0
 real(kind=RKIND),parameter,private:: am_g = PI*rho_g/6.0
 real(kind=RKIND),parameter,private:: bm_g = 3.0
 real(kind=RKIND),parameter,private:: am_i = PI*rho_i/6.0
 real(kind=RKIND),parameter,private:: bm_i = 3.0

!..Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
!.. Rain from Ferrier (1994), ice, snow, and graupel from
!.. Thompson et al (2008). Coefficient fv is zero for graupel/ice.
 real(kind=RKIND),parameter,private:: av_r = 4854.0
 real(kind=RKIND),parameter,private:: bv_r = 1.0
 real(kind=RKIND),parameter,private:: fv_r = 195.0
 real(kind=RKIND),parameter,private:: av_s = 40.0
 real(kind=RKIND),parameter,private:: bv_s = 0.55
 real(kind=RKIND),parameter,private:: fv_s = 100.0
 real(kind=RKIND),parameter,private:: av_g = 442.0
 real(kind=RKIND),parameter,private:: bv_g = 0.89
 real(kind=RKIND),parameter,private:: av_i = 1847.5
 real(kind=RKIND),parameter,private:: bv_i = 1.0
 real(kind=RKIND),parameter,private:: av_c = 0.316946E8
 real(kind=RKIND),parameter,private:: bv_c = 2.0

!..Capacitance of sphere and plates/aggregates: D**3, D**2
 real(kind=RKIND),parameter,private:: C_cube = 0.5
 real(kind=RKIND),parameter,private:: C_sqrd = 0.15

!..Collection efficiencies.  Rain/snow/graupel collection of cloud
!.. droplets use variables (Ef_rw, Ef_sw, Ef_gw respectively) and
!.. get computed elsewhere because they are dependent on stokes
!.. number.
 real(kind=RKIND),parameter,private:: Ef_si = 0.05
 real(kind=RKIND),parameter,private:: Ef_rs = 0.95
 real(kind=RKIND),parameter,private:: Ef_rg = 0.75
 real(kind=RKIND),parameter,private:: Ef_ri = 0.95

!..Minimum microphys values
!.. R1 value, 1.E-12, cannot be set lower because of numerical
!.. problems with Paul Field's moments and should not be set larger
!.. because of truncation problems in snow/ice growth.
 real(kind=RKIND),parameter,private:: R1 = 1.E-12
 real(kind=RKIND),parameter,private:: R2 = 1.E-6
 real(kind=RKIND),parameter,private:: eps = 1.E-15

!..Constants in Cooper curve relation for cloud ice number.
 real(kind=RKIND),parameter,private:: TNO = 5.0
 real(kind=RKIND),parameter,private:: ATO = 0.304

!..Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
 real(kind=RKIND),parameter,private:: rho_not = 101325.0/(287.05*298.0)

!..Schmidt number
 real(kind=RKIND),parameter,private:: Sc = 0.632
 real(kind=RKIND),private:: Sc3

!..Homogeneous freezing temperature
 real(kind=RKIND),parameter,private:: HGFR = 235.16

!..Water vapor and air gas constants at constant pressure
 real(kind=RKIND),parameter,private:: Rv = 461.5
 real(kind=RKIND),parameter,private:: oRv = 1./Rv
 real(kind=RKIND),parameter,private:: R = 287.04
 real(kind=RKIND),parameter,private:: Cp = 1004.0
 real(kind=RKIND),parameter,private:: R_uni = 8.314                    !J (mol K)-1

 real(kind=R8KIND),parameter,private:: k_b = 1.38065E-23               !Boltzmann constant [J/K]
 real(kind=R8KIND),parameter,private:: M_w = 18.01528E-3               !molecular mass of water [kg/mol]
 real(kind=R8KIND),parameter,private:: M_a = 28.96E-3                  !molecular mass of air [kg/mol]
 real(kind=R8KIND),parameter,private:: N_avo = 6.022E23                !Avogadro number [1/mol]
 real(kind=R8KIND),parameter,private:: ma_w = M_w / N_avo              !mass of water molecule [kg]
 real(kind=RKIND),parameter,private:: ar_volume = 4./3.*PI*(2.5e-6)**3 !assume radius of 0.025 micrometer, 2.5e-6 cm

!..Enthalpy of sublimation, vaporization, and fusion at 0C.
 real(kind=RKIND),parameter,private:: lsub = 2.834E6
 real(kind=RKIND),parameter,private:: lvap0 = 2.5E6
 real(kind=RKIND),parameter,private:: lfus = lsub - lvap0
 real(kind=RKIND),parameter,private:: olfus = 1./lfus

!..Ice initiates with this mass (kg), corresponding diameter calc.
!..Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
 real(kind=RKIND),parameter,private:: xm0i = 1.E-12
 real(kind=RKIND),parameter,private:: D0c = 1.E-6
 real(kind=RKIND),parameter,private:: D0r = 50.E-6
 real(kind=RKIND),parameter,private:: D0s = 200.E-6
 real(kind=RKIND),parameter,private:: D0g = 250.E-6
 real(kind=RKIND),private:: D0i, xm0s, xm0g

!..Lookup table dimensions
 integer,parameter,private:: nbins = 100
 integer,parameter,private:: nbc = nbins
 integer,parameter,private:: nbi = nbins
 integer,parameter,private:: nbr = nbins
 integer,parameter,private:: nbs = nbins
 integer,parameter,private:: nbg = nbins
 integer,parameter,private:: ntb_c = 37
 integer,parameter,private:: ntb_i = 64
 integer,parameter,private:: ntb_r = 37
 integer,parameter,private:: ntb_s = 28
 integer,parameter,private:: ntb_g = 28
 integer,parameter,private:: ntb_g1 = 28
 integer,parameter,private:: ntb_r1 = 37
 integer,parameter,private:: ntb_i1 = 55
 integer,parameter,private:: ntb_t = 9
 integer,private:: nic1, nic2, nii2, nii3, nir2, nir3, nis2, nig2, nig3
 integer,parameter,public:: ntb_arc = 7
 integer,parameter,public:: ntb_arw = 9
 integer,parameter,public:: ntb_art = 7
 integer,parameter,public:: ntb_arr = 5
 integer,parameter,public:: ntb_ark = 4
 integer,parameter,public:: ntb_IN = 55
 integer,private:: niIN2

 real(kind=R8KIND),dimension(nbins+1):: xDx
 real(kind=R8KIND),dimension(nbc):: Dc, dtc
 real(kind=R8KIND),dimension(nbi):: Di, dti
 real(kind=R8KIND),dimension(nbr):: Dr, dtr
 real(kind=R8KIND),dimension(nbs):: Ds, dts
 real(kind=R8KIND),dimension(nbg):: Dg, dtg
 real(kind=R8KIND),dimension(nbc):: t_Nc

!..Lookup tables for cloud water content (kg/m**3).
 real(kind=RKIND),dimension(ntb_c),parameter,private:: &
 r_c = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
         1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
         1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
         1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
         1.e-2/)

!..Lookup tables for cloud ice content (kg/m**3).
 real(kind=RKIND),dimension(ntb_i),parameter,private:: &
 r_i = (/1.e-10,2.e-10,3.e-10,4.e-10, &
         5.e-10,6.e-10,7.e-10,8.e-10,9.e-10, &
         1.e-9,2.e-9,3.e-9,4.e-9,5.e-9,6.e-9,7.e-9,8.e-9,9.e-9, &
         1.e-8,2.e-8,3.e-8,4.e-8,5.e-8,6.e-8,7.e-8,8.e-8,9.e-8, &
         1.e-7,2.e-7,3.e-7,4.e-7,5.e-7,6.e-7,7.e-7,8.e-7,9.e-7, &
         1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
         1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
         1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
         1.e-3/)

!..Lookup tables for rain content (kg/m**3).
 real(kind=RKIND),dimension(ntb_r),parameter,private:: &
 r_r = (/1.e-6,2.e-6,3.e-6,4.e-6,5.e-6,6.e-6,7.e-6,8.e-6,9.e-6, &
         1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
         1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
         1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
         1.e-2/)

!..Lookup tables for graupel content (kg/m**3).
 real(kind=RKIND),dimension(ntb_g),parameter,private:: &
 r_g = (/1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
         1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
         1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
         1.e-2/)

!..Lookup tables for snow content (kg/m**3).
 real(kind=RKIND),dimension(ntb_s),parameter,private:: &
 r_s = (/1.e-5,2.e-5,3.e-5,4.e-5,5.e-5,6.e-5,7.e-5,8.e-5,9.e-5, &
         1.e-4,2.e-4,3.e-4,4.e-4,5.e-4,6.e-4,7.e-4,8.e-4,9.e-4, &
         1.e-3,2.e-3,3.e-3,4.e-3,5.e-3,6.e-3,7.e-3,8.e-3,9.e-3, &
         1.e-2/)

!..Lookup tables for rain y-intercept parameter (/m**4).
 real(kind=RKIND),dimension(ntb_r1),parameter,private:: &
 N0r_exp = (/1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, &
             1.e7,2.e7,3.e7,4.e7,5.e7,6.e7,7.e7,8.e7,9.e7, &
             1.e8,2.e8,3.e8,4.e8,5.e8,6.e8,7.e8,8.e8,9.e8, &
             1.e9,2.e9,3.e9,4.e9,5.e9,6.e9,7.e9,8.e9,9.e9, &
             1.e10/)

!..Lookup tables for graupel y-intercept parameter (/m**4).
 real(kind=RKIND),dimension(ntb_g1),parameter,private:: &
 N0g_exp = (/1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
             1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
             1.e6,2.e6,3.e6,4.e6,5.e6,6.e6,7.e6,8.e6,9.e6, &
             1.e7/)

!..Lookup tables for ice number concentration (/m**3).
 real(kind=RKIND),dimension(ntb_i1),parameter,private:: &
 Nt_i = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
          1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
          1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
          1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
          1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
          1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
          1.e6/)

!..Aerosol table parameter: Number of available aerosols, vertical
!.. velocity, temperature, aerosol mean radius, and hygroscopicity.
 real(kind=RKIND),dimension(ntb_arc),parameter,private:: &
    ta_Na = (/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/)
 real(kind=RKIND),dimension(ntb_arw),parameter,private:: &
    ta_Ww = (/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/)
 real(kind=RKIND),dimension(ntb_art),parameter,private:: &
    ta_Tk = (/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/)
 real(kind=RKIND),dimension(ntb_arr),parameter,private:: &
    ta_Ra = (/0.01, 0.02, 0.04, 0.08, 0.16/)
 real(kind=RKIND),dimension(ntb_ark),parameter,private:: &
    ta_Ka = (/0.2, 0.4, 0.6, 0.8/)

!..Lookup tables for IN concentration (/m**3) from 0.001 to 1000/Liter.
 real(kind=RKIND),dimension(ntb_IN),parameter,private:: &
    Nt_IN = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
              1.e1,2.e1,3.e1,4.e1,5.e1,6.e1,7.e1,8.e1,9.e1, &
              1.e2,2.e2,3.e2,4.e2,5.e2,6.e2,7.e2,8.e2,9.e2, &
              1.e3,2.e3,3.e3,4.e3,5.e3,6.e3,7.e3,8.e3,9.e3, &
              1.e4,2.e4,3.e4,4.e4,5.e4,6.e4,7.e4,8.e4,9.e4, &
              1.e5,2.e5,3.e5,4.e5,5.e5,6.e5,7.e5,8.e5,9.e5, &
              1.e6/)

!..For snow moments conversions (from Field et al. 2005)
 real(kind=RKIND),dimension(10),parameter,private:: &
    sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
            0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
 real(kind=RKIND),dimension(10),parameter,private:: &
    sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
            0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

!..Temperatures (5 C interval 0 to -40) used in lookup tables.
 real(kind=RKIND),dimension(ntb_t),parameter,private:: &
    Tc = (/-0.01, -5., -10., -15., -20., -25., -30., -35., -40./)

!..Lookup tables for various accretion/collection terms.
!.. ntb_x refers to the number of elements for rain, snow, graupel,
!.. and temperature array indices.  Variables beginning with t-p/c/m/n
!.. represent lookup tables.  Save compile-time memory by making
!.. allocatable (2009Jun12, J. Michalakes).
!integer,parameter,private:: R8SIZE = 8
!integer,parameter,private:: R4SIZE = 4
 integer,parameter,private:: R8SIZE = R8KIND
 integer,parameter,private:: R4SIZE = R4KIND
 real(kind=R8SIZE),allocatable,dimension(:,:,:,:):: &
    tcg_racg,tmr_racg,tcr_gacr,tmg_gacr, &
    tnr_racg,tnr_gacr
 real(kind=R8SIZE),allocatable,dimension(:,:,:,:):: &
    tcs_racs1,tmr_racs1,tcs_racs2, tmr_racs2, &
    tcr_sacr1,tms_sacr1,tcr_sacr2, tms_sacr2, &
    tnr_racs1,tnr_racs2,tnr_sacr1, tnr_sacr2
 real(kind=R8SIZE),allocatable,dimension(:,:,:,:):: &
    tpi_qcfz,tni_qcfz
 real(kind=R8SIZE),allocatable,dimension(:,:,:,:):: &
    tpi_qrfz,tpg_qrfz,tni_qrfz,tnr_qrfz
 real(kind=R8SIZE),allocatable,dimension(:,:):: &
    tps_iaus,tni_iaus,tpi_ide
 real(kind=R8SIZE),allocatable,dimension(:,:):: t_Efrw
 real(kind=R8SIZE),allocatable,dimension(:,:):: t_Efsw
 real(kind=R8SIZE),allocatable,dimension(:,:,:):: tnr_rev
 real(kind=R8SIZE),allocatable,dimension(:,:,:):: &
    tpc_wev,tnc_wev
 real(kind=R4SIZE),allocatable,dimension(:,:,:,:,:):: tnccn_act

!..Variables holding a bunch of exponents and gamma values (cloud water,
!.. cloud ice, rain, snow, then graupel).
 real(kind=RKIND),dimension(5,15),private:: cce, ccg
 real(kind=RKIND),dimension(15),private::  ocg1, ocg2
 real(kind=RKIND),dimension(7),private:: cie, cig
 real(kind=RKIND),private:: oig1, oig2, obmi
 real(kind=RKIND),dimension(13),private:: cre, crg
 real(kind=RKIND),private:: ore1, org1, org2, org3, obmr
 real(kind=RKIND),dimension(18),private:: cse, csg
 real(kind=RKIND),private:: oams, obms, ocms
 real(kind=RKIND),dimension(12),private:: cge, cgg
 real(kind=RKIND),private:: oge1, ogg1, ogg2, ogg3, oamg, obmg, ocmg

!..Declaration of precomputed constants in various rate eqns.
 real(kind=RKIND):: t1_qr_qc, t1_qr_qi, t2_qr_qi, t1_qg_qc, t1_qs_qc, t1_qs_qi
 real(kind=RKIND):: t1_qr_ev, t2_qr_ev
 real(kind=RKIND):: t1_qs_sd, t2_qs_sd, t1_qg_sd, t2_qg_sd
 real(kind=RKIND):: t1_qs_me, t2_qs_me, t1_qg_me, t2_qg_me


 contains


!==================================================================================================================
 subroutine thompson_init(l_mp_tables)
 implicit none
!==================================================================================================================

!input arguments:
 logical,intent(in):: l_mp_tables
 
 integer,parameter:: open_OK = 0
 integer:: i,j,k,l,m,n
 integer:: istat
 logical:: micro_init
 integer:: mp_unit

!..Allocate space for lookup tables (J. Michalakes 2009Jun08).
 micro_init = .false.

!==================================================================================================================

 if(.not. allocated(tcg_racg) ) then
    allocate(tcg_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
    micro_init = .true.
 endif

 if(.not.allocated(tmr_racg)) allocate(tmr_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
 if(.not.allocated(tcr_gacr)) allocate(tcr_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
 if(.not.allocated(tmg_gacr)) allocate(tmg_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))
 if(.not.allocated(tnr_racg)) allocate(tnr_racg(ntb_g1,ntb_g,ntb_r1,ntb_r))
 if(.not.allocated(tnr_gacr)) allocate(tnr_gacr(ntb_g1,ntb_g,ntb_r1,ntb_r))

 if(.not.allocated(tcs_racs1)) allocate(tcs_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tmr_racs1)) allocate(tmr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tcs_racs2)) allocate(tcs_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tmr_racs2)) allocate(tmr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tcr_sacr1)) allocate(tcr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tms_sacr1)) allocate(tms_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tcr_sacr2)) allocate(tcr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tms_sacr2)) allocate(tms_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tnr_racs1)) allocate(tnr_racs1(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tnr_racs2)) allocate(tnr_racs2(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tnr_sacr1)) allocate(tnr_sacr1(ntb_s,ntb_t,ntb_r1,ntb_r))
 if(.not.allocated(tnr_sacr2)) allocate(tnr_sacr2(ntb_s,ntb_t,ntb_r1,ntb_r))

 if(.not.allocated(tpi_qcfz)) allocate(tpi_qcfz(ntb_c,nbc,45,ntb_IN))
 if(.not.allocated(tni_qcfz)) allocate(tni_qcfz(ntb_c,nbc,45,ntb_IN))

 if(.not.allocated(tpi_qrfz)) allocate(tpi_qrfz(ntb_r,ntb_r1,45,ntb_IN))
 if(.not.allocated(tpg_qrfz)) allocate(tpg_qrfz(ntb_r,ntb_r1,45,ntb_IN))
 if(.not.allocated(tni_qrfz)) allocate(tni_qrfz(ntb_r,ntb_r1,45,ntb_IN))
 if(.not.allocated(tnr_qrfz)) allocate(tnr_qrfz(ntb_r,ntb_r1,45,ntb_IN))

 if(.not.allocated(tps_iaus)) allocate(tps_iaus(ntb_i,ntb_i1))
 if(.not.allocated(tni_iaus)) allocate(tni_iaus(ntb_i,ntb_i1))
 if(.not.allocated(tpi_ide)) allocate(tpi_ide(ntb_i,ntb_i1))

 if(.not.allocated(t_Efrw)) allocate(t_Efrw(nbr,nbc))
 if(.not.allocated(t_Efsw)) allocate(t_Efsw(nbs,nbc))

 if(.not.allocated(tnr_rev)) allocate(tnr_rev(nbr, ntb_r1, ntb_r))
 if(.not.allocated(tpc_wev)) allocate(tpc_wev(nbc,ntb_c,nbc))
 if(.not.allocated(tnc_wev)) allocate(tnc_wev(nbc,ntb_c,nbc))
 
 if(.not.allocated(tnccn_act)) &
    allocate(tnccn_act(ntb_arc,ntb_arw,ntb_art,ntb_arr,ntb_ark))


 if(micro_init) then

!..From Martin et al. (1994), assign gamma shape parameter mu for cloud
!.. drops according to general dispersion characteristics (disp=~0.25
!.. for Maritime and 0.45 for Continental).
!.. disp=SQRT((mu+2)/(mu+1) - 1) so mu varies from 15 for Maritime
!.. to 2 for really dirty air.  This not used in 2-moment cloud water
!.. scheme and nu_c used instead and varies from 2 to 15 (integer-only).
!   mu_c = MIN(15., (1000.E6/Nt_c + 2.))

!..Schmidt number to one-third used numerous times.
    Sc3 = Sc**(1./3.)

!..Compute min ice diam from mass, min snow/graupel mass from diam.
    D0i = (xm0i/am_i)**(1./bm_i)
    xm0s = am_s * D0s**bm_s
    xm0g = am_g * D0g**bm_g

!..These constants various exponents and gamma() assoc with cloud,
!.. rain, snow, and graupel.
    do n = 1, 15
       cce(1,n) = n + 1.
       cce(2,n) = bm_r + n + 1.
       cce(3,n) = bm_r + n + 4.
       cce(4,n) = n + bv_c + 1.
       cce(5,n) = bm_r + n + bv_c + 1.
       ccg(1,n) = wgamma(cce(1,n))
       ccg(2,n) = wgamma(cce(2,n))
       ccg(3,n) = wgamma(cce(3,n))
       ccg(4,n) = wgamma(cce(4,n))
       ccg(5,n) = wgamma(cce(5,n))
       ocg1(n) = 1./ccg(1,n)
       ocg2(n) = 1./ccg(2,n)
    enddo

    cie(1) = mu_i + 1.
    cie(2) = bm_i + mu_i + 1.
    cie(3) = bm_i + mu_i + bv_i + 1.
    cie(4) = mu_i + bv_i + 1.
    cie(5) = mu_i + 2.
    cie(6) = bm_i*0.5 + mu_i + bv_i + 1.
    cie(7) = bm_i*0.5 + mu_i + 1.
    cig(1) = wgamma(cie(1))
    cig(2) = wgamma(cie(2))
    cig(3) = wgamma(cie(3))
    cig(4) = wgamma(cie(4))
    cig(5) = wgamma(cie(5))
    cig(6) = wgamma(cie(6))
    cig(7) = wgamma(cie(7))
    oig1 = 1./cig(1)
    oig2 = 1./cig(2)
    obmi = 1./bm_i

    cre(1) = bm_r + 1.
    cre(2) = mu_r + 1.
    cre(3) = bm_r + mu_r + 1.
    cre(4) = bm_r*2. + mu_r + 1.
    cre(5) = mu_r + bv_r + 1.
    cre(6) = bm_r + mu_r + bv_r + 1.
    cre(7) = bm_r*0.5 + mu_r + bv_r + 1.
    cre(8) = bm_r + mu_r + bv_r + 3.
    cre(9) = mu_r + bv_r + 3.
    cre(10) = mu_r + 2.
    cre(11) = 0.5*(bv_r + 5. + 2.*mu_r)
    cre(12) = bm_r*0.5 + mu_r + 1.
    cre(13) = bm_r*2. + mu_r + bv_r + 1.
    do n = 1, 13
       crg(n) = wgamma(cre(n))
    enddo
    obmr = 1./bm_r
    ore1 = 1./cre(1)
    org1 = 1./crg(1)
    org2 = 1./crg(2)
    org3 = 1./crg(3)

    cse(1) = bm_s + 1.
    cse(2) = bm_s + 2.
    cse(3) = bm_s*2.
    cse(4) = bm_s + bv_s + 1.
    cse(5) = bm_s*2. + bv_s + 1.
    cse(6) = bm_s*2. + 1.
    cse(7) = bm_s + mu_s + 1.
    cse(8) = bm_s + mu_s + 2.
    cse(9) = bm_s + mu_s + 3.
    cse(10) = bm_s + mu_s + bv_s + 1.
    cse(11) = bm_s*2. + mu_s + bv_s + 1.
    cse(12) = bm_s*2. + mu_s + 1.
    cse(13) = bv_s + 2.
    cse(14) = bm_s + bv_s
    cse(15) = mu_s + 1.
    cse(16) = 1.0 + (1.0 + bv_s)/2.
    cse(17) = cse(16) + mu_s + 1.
    cse(18) = bv_s + mu_s + 3.
    do n = 1, 18
       csg(n) = wgamma(cse(n))
    enddo
    oams = 1./am_s
    obms = 1./bm_s
    ocms = oams**obms

    cge(1) = bm_g + 1.
    cge(2) = mu_g + 1.
    cge(3) = bm_g + mu_g + 1.
    cge(4) = bm_g*2. + mu_g + 1.
    cge(5) = bm_g*2. + mu_g + bv_g + 1.
    cge(6) = bm_g + mu_g + bv_g + 1.
    cge(7) = bm_g + mu_g + bv_g + 2.
    cge(8) = bm_g + mu_g + bv_g + 3.
    cge(9) = mu_g + bv_g + 3.
    cge(10) = mu_g + 2.
    cge(11) = 0.5*(bv_g + 5. + 2.*mu_g)
    cge(12) = 0.5*(bv_g + 5.) + mu_g
    do n = 1, 12
       cgg(n) = wgamma(cge(n))
    enddo
    oamg = 1./am_g
    obmg = 1./bm_g
    ocmg = oamg**obmg
    oge1 = 1./cge(1)
    ogg1 = 1./cgg(1)
    ogg2 = 1./cgg(2)
    ogg3 = 1./cgg(3)

!+---+-----------------------------------------------------------------+
!..Simplify various rate eqns the best we can now.
!+---+-----------------------------------------------------------------+

!..Rain collecting cloud water and cloud ice
    t1_qr_qc = PI*.25*av_r * crg(9)
    t1_qr_qi = PI*.25*av_r * crg(9)
    t2_qr_qi = PI*.25*am_r*av_r * crg(8)

!..Graupel collecting cloud water
    t1_qg_qc = PI*.25*av_g * cgg(9)

!..Snow collecting cloud water
    t1_qs_qc = PI*.25*av_s

!..Snow collecting cloud ice
    t1_qs_qi = PI*.25*av_s

!..Evaporation of rain; ignore depositional growth of rain.
    t1_qr_ev = 0.78 * crg(10)
    t2_qr_ev = 0.308*Sc3*SQRT(av_r) * crg(11)

!..Sublimation/depositional growth of snow
    t1_qs_sd = 0.86
    t2_qs_sd = 0.28*Sc3*SQRT(av_s)

!..Melting of snow
    t1_qs_me = PI*4.*C_sqrd*olfus * 0.86
    t2_qs_me = PI*4.*C_sqrd*olfus * 0.28*Sc3*SQRT(av_s)

!..Sublimation/depositional growth of graupel
    t1_qg_sd = 0.86 * cgg(10)
    t2_qg_sd = 0.28*Sc3*SQRT(av_g) * cgg(11)

!..Melting of graupel
    t1_qg_me = PI*4.*C_cube*olfus * 0.86 * cgg(10)
    t2_qg_me = PI*4.*C_cube*olfus * 0.28*Sc3*SQRT(av_g) * cgg(11)

!..Constants for helping find lookup table indexes.
    nic2 = nint(alog10(r_c(1)))
    nii2 = nint(alog10(r_i(1)))
    nii3 = nint(alog10(Nt_i(1)))
    nir2 = nint(alog10(r_r(1)))
    nir3 = nint(alog10(N0r_exp(1)))
    nis2 = nint(alog10(r_s(1)))
    nig2 = nint(alog10(r_g(1)))
    nig3 = nint(alog10(N0g_exp(1)))
    niIN2 = nint(alog10(Nt_IN(1)))

!..Create bins of cloud water (from min diameter up to 100 microns).
    Dc(1) = D0c*1.0d0
    dtc(1) = D0c*1.0d0
    do n = 2, nbc
       Dc(n) = Dc(n-1) + 1.0D-6
       dtc(n) = (Dc(n) - Dc(n-1))
    enddo

!..Create bins of cloud ice (from min diameter up to 5x min snow size).
    xDx(1) = D0i*1.0d0
    xDx(nbi+1) = 5.0d0*D0s
    do n = 2, nbi
       xDx(n) = dexp(real(n-1,kind=R8KIND)/real(nbi,kind=R8KIND) &
              * dlog(xDx(nbi+1)/xDx(1)) + dlog(xDx(1)))
    enddo
    do n = 1, nbi
       Di(n) = dsqrt(xDx(n)*xDx(n+1))
       dti(n) = xDx(n+1) - xDx(n)
    enddo

!..Create bins of rain (from min diameter up to 5 mm).
    xDx(1) = D0r*1.0d0
    xDx(nbr+1) = 0.005d0
    do n = 2, nbr
       xDx(n) = dexp(real(n-1,kind=R8KIND)/real(nbr,kind=R8KIND) &
              * dlog(xDx(nbr+1)/xDx(1)) + dlog(xDx(1)))
    enddo
    do n = 1, nbr
       Dr(n) = dsqrt(xDx(n)*xDx(n+1))
       dtr(n) = xDx(n+1) - xDx(n)
    enddo

!..Create bins of snow (from min diameter up to 2 cm).
    xDx(1) = D0s*1.0d0
    xDx(nbs+1) = 0.02d0
    do n = 2, nbs
       xDx(n) = dexp(real(n-1,kind=R8KIND)/real(nbs,kind=R8KIND) &
              * dlog(xDx(nbs+1)/xDx(1)) + dlog(xDx(1)))
    enddo
    do n = 1, nbs
       Ds(n) = dsqrt(xDx(n)*xDx(n+1))
       dts(n) = xDx(n+1) - xDx(n)
    enddo

!..Create bins of graupel (from min diameter up to 5 cm).
    xDx(1) = D0g*1.0d0
    xDx(nbg+1) = 0.05d0
    do n = 2, nbg
       xDx(n) = dexp(real(n-1,kind=R8KIND)/real(nbg,kind=R8KIND) &
              * dlog(xDx(nbg+1)/xDx(1)) + dlog(xDx(1)))
    enddo
    do n = 1, nbg
       Dg(n) = dsqrt(xDx(n)*xDx(n+1))
       dtg(n) = xDx(n+1) - xDx(n)
    enddo

!..Create bins of cloud droplet number concentration (1 to 3000 per cc).
    xDx(1) = 1.0d0
    xDx(nbc+1) = 3000.0d0
    do n = 2, nbc
       xDx(n) = dexp(real(n-1,kind=R8KIND)/real(nbc,KIND=R8KIND) &
              * dlog(xDx(nbc+1)/xDx(1)) + dlog(xDx(1)))
    enddo
    do n = 1, nbc
       t_Nc(n) = dsqrt(xDx(n)*xDx(n+1)) * 1.D6
    enddo
    nic1 = dlog(t_Nc(nbc)/t_Nc(1))

!+---+-----------------------------------------------------------------+
!..Create lookup tables for most costly calculations.
!+---+-----------------------------------------------------------------+

    do m = 1,ntb_r
       do k = 1,ntb_r1
          do j = 1,ntb_g
             do i = 1,ntb_g1
                tcg_racg(i,j,k,m) = 0._R8KIND
                tmr_racg(i,j,k,m) = 0._R8KIND
                tcr_gacr(i,j,k,m) = 0._R8KIND
                tmg_gacr(i,j,k,m) = 0._R8KIND
                tnr_racg(i,j,k,m) = 0._R8KIND
                tnr_gacr(i,j,k,m) = 0._R8KIND
             enddo
          enddo
       enddo
    enddo

    do m = 1,ntb_r
       do k = 1,ntb_r1
          do j = 1,ntb_t
             do i = 1,ntb_s
                tcs_racs1(i,j,k,m) = 0._R8KIND
                tmr_racs1(i,j,k,m) = 0._R8KIND
                tcs_racs2(i,j,k,m) = 0._R8KIND
                tmr_racs2(i,j,k,m) = 0._R8KIND
                tcr_sacr1(i,j,k,m) = 0._R8KIND
                tms_sacr1(i,j,k,m) = 0._R8KIND
                tcr_sacr2(i,j,k,m) = 0._R8KIND
                tms_sacr2(i,j,k,m) = 0._R8KIND
                tnr_racs1(i,j,k,m) = 0._R8KIND
                tnr_racs2(i,j,k,m) = 0._R8KIND
                tnr_sacr1(i,j,k,m) = 0._R8KIND
                tnr_sacr2(i,j,k,m) = 0._R8KIND
             enddo
          enddo
       enddo
    enddo

    do m = 1,ntb_IN
       do k = 1,45
          do j = 1,ntb_r1
             do i = 1,ntb_r
                tpi_qrfz(i,j,k,m) = 0._R8KIND
                tni_qrfz(i,j,k,m) = 0._R8KIND
                tpg_qrfz(i,j,k,m) = 0._R8KIND
                tnr_qrfz(i,j,k,m) = 0._R8KIND
             enddo
          enddo
          do j = 1,nbc
             do i = 1,ntb_c
                tpi_qcfz(i,j,k,m) = 0._R8KIND
                tni_qcfz(i,j,k,m) = 0._R8KIND
             enddo
          enddo
       enddo
    enddo

    do j = 1,ntb_i1
       do i = 1,ntb_i
          tps_iaus(i,j) = 0._R8KIND
          tni_iaus(i,j) = 0._R8KIND
          tpi_ide(i,j) = 0._R8KIND
       enddo
    enddo

    do j = 1,nbc
       do i = 1, nbr
          t_Efrw(i,j) = 0._RKIND
       enddo
       do i = 1,nbs
          t_Efsw(i,j) = 0._RKIND
       enddo
    enddo

    do k = 1,ntb_r
       do j = 1,ntb_r1
          do i = 1,nbr
             tnr_rev(i,j,k) = 0._R8KIND
          enddo
       enddo
    enddo

    do k = 1,nbc
       do j = 1,ntb_c
          do i = 1,nbc
             tpc_wev(i,j,k) = 0._R8KIND
             tnc_wev(i,j,k) = 0._R8KIND
          enddo
       enddo
    enddo

    do m = 1,ntb_ark
       do l = 1,ntb_arr
          do k = 1,ntb_art
             do j = 1,ntb_arw
                do i = 1,ntb_arc
                   tnccn_act(i,j,k,l,m) = 1._RKIND
                enddo
             enddo
          enddo
       enddo
    enddo

!..Check that the look-up tables are available.
    if(.not. l_mp_tables) return

!..Collision efficiency between rain/snow and cloud water.
!   call physics_message('--- creating qc collision eff tables')
    call table_Efrw
    call table_Efsw

!..Drop evaporation.
!   call physics_message('--- creating rain evap table')
    call table_dropEvap

!..Rain collecting graupel & graupel collecting rain.
    call mpas_new_unit(mp_unit, unformatted = .true.)
    open(unit=mp_unit,file='MP_THOMPSON_QRacrQG_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
         iostat = istat)
    if(istat /= open_OK) &
       call physics_error_fatal('subroutine thompson_init: ' // &
                                'failure opening MP_THOMPSON_QRacrQG.DBL')
    read(mp_unit) tcg_racg
    read(mp_unit) tmr_racg
    read(mp_unit) tcr_gacr
    read(mp_unit) tmg_gacr
    read(mp_unit) tnr_racg
    read(mp_unit) tnr_gacr
    close(unit=mp_unit)
!   call mpas_log_write('--- end read MP_THOMPSON_QRacrQG.DBL')
!   call mpas_log_write('max tcg_racg =',realArgs=(/maxval(tcg_racg)/))
!   call mpas_log_write('min tcg_racg =',realArgs=(/minval(tcg_racg)/))
!   call mpas_log_write('max tmr_racg =',realArgs=(/maxval(tmr_racg)/))
!   call mpas_log_write('min tmr_racg =',realArgs=(/minval(tmr_racg)/))
!   call mpas_log_write('max tcr_gacr =',realArgs=(/maxval(tcr_gacr)/))
!   call mpas_log_write('min tcr_gacr =',realArgs=(/minval(tcr_gacr)/))
!   call mpas_log_write('max tmg_gacr =',realArgs=(/maxval(tmg_gacr)/))
!   call mpas_log_write('min tmg_gacr =',realArgs=(/minval(tmg_gacr)/))
!   call mpas_log_write('max tnr_racg =',realArgs=(/maxval(tnr_racg)/))
!   call mpas_log_write('min tnr_racg =',realArgs=(/minval(tnr_racg)/))
!   call mpas_log_write('max tnr_gacr =',realArgs=(/maxval(tnr_gacr)/))
!   call mpas_log_write('min tnr_gacr =',realArgs=(/minval(tnr_gacr)/))

!..Rain collecting snow & snow collecting rain.
    open(unit=mp_unit,file='MP_THOMPSON_QRacrQS_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
         iostat=istat)
    if(istat /= open_OK) &
       call physics_error_fatal('subroutine thompson_init: ' // &
                                'failure opening MP_THOMPSON_QRacrQS.DBL')
    read(mp_unit) tcs_racs1
    read(mp_unit) tmr_racs1
    read(mp_unit) tcs_racs2
    read(mp_unit) tmr_racs2
    read(mp_unit) tcr_sacr1
    read(mp_unit) tms_sacr1
    read(mp_unit) tcr_sacr2
    read(mp_unit) tms_sacr2
    read(mp_unit) tnr_racs1
    read(mp_unit) tnr_racs2
    read(mp_unit) tnr_sacr1
    read(mp_unit) tnr_sacr2
    close(unit=mp_unit)
!   call mpas_log_write('--- end read MP_THOMPSON_QRacrQS.DBL')
!   call mpas_log_write('max tcs_racs1 =',realArgs=(/maxval(tcs_racs1)/))
!   call mpas_log_write('min tcs_racs1 =',realArgs=(/minval(tcs_racs1)/))
!   call mpas_log_write('max tmr_racs1 =',realArgs=(/maxval(tmr_racs1)/))
!   call mpas_log_write('min tmr_racs1 =',realArgs=(/minval(tmr_racs1)/))
!   call mpas_log_write('max tcs_racs2 =',realArgs=(/maxval(tcs_racs2)/))
!   call mpas_log_write('min tcs_racs2 =',realArgs=(/minval(tcs_racs2)/))
!   call mpas_log_write('max tmr_racs2 =',realArgs=(/maxval(tmr_racs2)/))
!   call mpas_log_write('min tmr_racs2 =',realArgs=(/minval(tmr_racs2)/))
!   call mpas_log_write('max tcr_sacr1 =',realArgs=(/maxval(tcr_sacr1)/))
!   call mpas_log_write('min tcr_sacr1 =',realArgs=(/minval(tcr_sacr1)/))
!   call mpas_log_write('max tms_sacr1 =',realArgs=(/maxval(tms_sacr1)/))
!   call mpas_log_write('min tms_sacr1 =',realArgs=(/minval(tms_sacr1)/))
!   call mpas_log_write('max tcr_sacr2 =',realArgs=(/maxval(tcr_sacr2)/))
!   call mpas_log_write('min tcr_sacr2 =',realArgs=(/minval(tcr_sacr2)/))
!   call mpas_log_write('max tms_sacr2 =',realArgs=(/maxval(tms_sacr2)/))
!   call mpas_log_write('min tms_sacr2 =',realArgs=(/minval(tms_sacr2)/))
!   call mpas_log_write('max tnr_racs1 =',realArgs=(/maxval(tnr_racs1)/))
!   call mpas_log_write('min tnr_racs1 =',realArgs=(/minval(tnr_racs1)/))
!   call mpas_log_write('max tnr_racs2 =',realArgs=(/maxval(tnr_racs2)/))
!   call mpas_log_write('min tnr_racs2 =',realArgs=(/minval(tnr_racs2)/))
!   call mpas_log_write('max tnr_sacr1 =',realArgs=(/maxval(tnr_sacr1)/))
!   call mpas_log_write('min tnr_sacr1 =',realArgs=(/minval(tnr_sacr1)/))
!   call mpas_log_write('max tnr_sacr2 =',realArgs=(/maxval(tnr_sacr2)/))
!   call mpas_log_write('min tnr_sacr2 =',realArgs=(/minval(tnr_sacr2)/))

!..Cloud water and rain freezing (Bigg, 1953).
    open(unit=mp_unit,file='MP_THOMPSON_freezeH2O_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
         iostat=istat)
    if(istat /= open_OK) &
       call physics_error_fatal('subroutine thompson_init: ' // &
                                'failure opening MP_THOMPSON_freezeH2O.DBL')
    read(mp_unit) tpi_qrfz
    read(mp_unit) tni_qrfz
    read(mp_unit) tpg_qrfz
    read(mp_unit) tnr_qrfz
    read(mp_unit) tpi_qcfz
    read(mp_unit) tni_qcfz
    close(unit=mp_unit)
!   call mpas_log_write('--- end read MP_THOMPSON_freezeH2O.DBL:')
!   call mpas_log_write('max tpi_qrfz =',realArgs=(/maxval(tpi_qrfz)/))
!   call mpas_log_write('min tpi_qrfz =',realArgs=(/minval(tpi_qrfz)/))
!   call mpas_log_write('max tni_qrfz =',realArgs=(/maxval(tni_qrfz)/))
!   call mpas_log_write('min tni_qrfz =',realArgs=(/minval(tni_qrfz)/))
!   call mpas_log_write('max tpg_qrfz =',realArgs=(/maxval(tpg_qrfz)/))
!   call mpas_log_write('min tpg_qrfz =',realArgs=(/minval(tpg_qrfz)/))
!   call mpas_log_write('max tnr_qrfz =',realArgs=(/maxval(tnr_qrfz)/))
!   call mpas_log_write('min tnr_qrfz =',realArgs=(/minval(tnr_qrfz)/))
!   call mpas_log_write('max tpi_qcfz =',realArgs=(/maxval(tpi_qcfz)/))
!   call mpas_log_write('min tpi_qcfz =',realArgs=(/minval(tpi_qcfz)/))
!   call mpas_log_write('max tni_qcfz =',realArgs=(/maxval(tni_qcfz)/))
!   call mpas_log_write('min tni_qcfz =',realArgs=(/minval(tni_qcfz)/))

!..Conversion of some ice mass into snow category.
    open(unit=mp_unit,file='MP_THOMPSON_QIautQS_DATA.DBL',form='UNFORMATTED',status='OLD',action='READ', &
         iostat=istat)
    if(istat /= open_OK) &
       call physics_error_fatal('subroutine thompson_init: ' // &
                                'failure opening MP_THOMPSON_QIautQS.DBL')
    read(mp_unit) tpi_ide
    read(mp_unit) tps_iaus
    read(mp_unit) tni_iaus
    close(unit=mp_unit)
    call mpas_release_unit(mp_unit)
!   call mpas_log_write('--- end read MP_THOMPSON_QIautQS.DBL')
!   call mpas_log_write('max tps_iaus =',realArgs=(/maxval(tps_iaus)/))
!   call mpas_log_write('min tps_iaus =',realArgs=(/minval(tps_iaus)/))
!   call mpas_log_write('max tni_iaus =',realArgs=(/maxval(tni_iaus)/))
!   call mpas_log_write('min tni_iaus =',realArgs=(/minval(tni_iaus)/))

!..Initialize various constants for computing radar reflectivity.
    xam_r = am_r
    xbm_r = bm_r
    xmu_r = mu_r
    xam_s = am_s
    xbm_s = bm_s
    xmu_s = mu_s
    xam_g = am_g
    xbm_g = bm_g
    xmu_g = mu_g
    call radar_init

 endif ! end "if(micro_init)"

 end subroutine thompson_init

!==================================================================================================================
 subroutine mp_thompson_run(qv,qc,qr,qi,qs,qg,ni,nr,nc,nwfa,nifa,nwfa2d,nifa2d, &
                            th,pi,p,w,dz,dt_in,itimestep,rainnc,rainncv,snownc, &
                            snowncv,graupelnc,graupelncv,sr,rainprod,evapprod,  &
                            re_cloud,re_ice,re_snow,has_reqc,has_reqi,has_reqs, &
                            ntc,muc,pflrain,pflice,pflsnow,pflgraul,            &
                            ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
 implicit none
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte
 integer,intent(in):: has_reqc,has_reqi,has_reqs
 integer,intent(in):: itimestep

 real(kind=RKIND),intent(in):: dt_in
 real(kind=RKIND),intent(in),dimension(its:ite):: ntc,muc
 real(kind=RKIND),intent(in),dimension(its:ite),optional:: nwfa2d,nifa2d
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: pi,p,w,dz

!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite):: rainnc,rainncv,sr
 real(kind=RKIND),intent(inout),dimension(its:ite),optional:: snownc,snowncv,graupelnc,graupelncv

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: th,qv,qc,qr,qi,qs,qg,ni,nr
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: re_cloud,re_ice,re_snow
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte),optional:: rainprod,evapprod
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte),optional:: pflrain,pflice,pflsnow,pflgraul
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte),optional:: nc,nwfa,nifa

!--- local variables and arrays:
 integer:: i,j,k

 real(kind=RKIND):: dt
 real(kind=RKIND):: pptrain,pptsnow,pptgraul,pptice
 real(kind=RKIND),dimension(kts:kte):: qv1d,qc1d,qi1d,qr1d,qs1d,qg1d
 real(kind=RKIND),dimension(kts:kte):: ni1d,nr1d,nc1d,nwfa1d,nifa1d
 real(kind=RKIND),dimension(kts:kte):: t1d,p1d,w1d,dz1d,rho
 real(kind=RKIND),dimension(kts:kte):: re_qc1d,re_qi1d,re_qs1d
 real(kind=RKIND),dimension(kts:kte):: rainprod1d, evapprod1d
 real(kind=RKIND),dimension(kts:kte):: ppfr1d,ppfi1d,ppfs1d,ppfg1d

!------------------------------------------------------------------------------------------------------------------

 dt = dt_in
   
!--- BEGIN CLOUD MICROPHYSICS SCHEME:
 do i = its,ite

    pptrain  = 0._RKIND
    pptsnow  = 0._RKIND
    pptgraul = 0._RKIND
    pptice   = 0._RKIND

    rainncv(i)    = 0._RKIND
    snowncv(i)    = 0._RKIND
    graupelncv(i) = 0._RKIND
    sr(i)         = 0._RKIND

    Nt_c = ntc(i)
    mu_c = muc(i)

    !--- initialization of 1D column arrays:
    do k = kts, kte
       t1d(k)  = th(i,k)*pi(i,k)
       p1d(k)  = p(i,k)
       w1d(k)  = w(i,k)
       dz1d(k) = dz(i,k)
       qv1d(k) = qv(i,k)
       qc1d(k) = qc(i,k)
       qi1d(k) = qi(i,k)
       qr1d(k) = qr(i,k)
       qs1d(k) = qs(i,k)
       qg1d(k) = qg(i,k)
       ni1d(k) = ni(i,k)
       nr1d(k) = nr(i,k)
       rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))

       rainprod1d(k) = 0._RKIND
       evapprod1d(k) = 0._RKIND
       ppfr1d(k)     = 0._RKIND
       ppfi1d(k)     = 0._RKIND
       ppfs1d(k)     = 0._RKIND
       ppfg1d(k)     = 0._RKIND
    enddo

    if(is_aerosol_aware) then
       do k = kts,kte
          nc1d(k)   = nc(i,k)
          nwfa1d(k) = nwfa(i,k)
          nifa1d(k) = nifa(i,k)
       enddo
    else
       do k = kts, kte
          nc1d(k)   = Nt_c/rho(k)
          nwfa1d(k) = 11.1E6/rho(k)
          nifa1d(k) = naIN1*0.01/rho(k)
       enddo
    endif


    !--- call to 1D Thompson cloud microphysics scheme:
    call thompson_1d( &
            qv1d   , qc1d    , qi1d    , qr1d     , qs1d   , qg1d       , ni1d       , &
            nr1d   , nc1d    , nwfa1d  , nifa1d   , t1d    , p1d        , w1d        , &
            dz1d   , pptrain , pptsnow , pptgraul , pptice , rainprod1d , evapprod1d , &
            ppfr1d , ppfi1d  , ppfs1d  , ppfg1d   , kts    , kte        , dt           &
                    )


    rainncv(i)    = pptrain + pptsnow + pptgraul + pptice
    rainnc(i)     = rainnc(i) + pptrain + pptsnow + pptgraul + pptice
    snowncv(i)    = pptsnow + pptice
    snownc(i)     = snownc(i) + pptsnow + pptice
    graupelncv(i) = pptgraul
    graupelnc(i)  = graupelnc(i) + pptgraul
    sr(i)         = (pptsnow + pptgraul + pptice)/(rainncv(i)+1.e-12)

    !--- update 2D arrays after call to cloud microphysics:
    do k = kts,kte
       rainprod(i,k) = rainprod1d(k)
       evapprod(i,k) = evapprod1d(k)
       pflrain(i,k)  = ppfr1d(k)
       pflice(i,k)   = ppfi1d(k)
       pflsnow(i,k)  = ppfs1d(k)
       pflgraul(i,k) = ppfg1d(k)
    enddo

    do k = kts, kte
       th(i,k) = t1d(k)/pi(i,k)
       qv(i,k) = qv1d(k)
       qc(i,k) = qc1d(k)
       qi(i,k) = qi1d(k)
       qr(i,k) = qr1d(k)
       qs(i,k) = qs1d(k)
       qg(i,k) = qg1d(k)
       ni(i,k) = ni1d(k)
       nr(i,k) = nr1d(k)
    enddo

    if(is_aerosol_aware) then
       if(present(nwfa2d) .and. present(nifa2d)) then
          nwfa1d(kts) = nwfa1d(kts) + nwfa2d(i)*dt_in
          nifa1d(kts) = nifa1d(kts) + nifa2d(i)*dt_in
       endif

       do k = kts,kte
          nc(i,k)   = nc1d(k)
          nwfa(i,k) = nwfa1d(k)
          nifa(i,k) = nifa1d(k)
       enddo
    endif

    !--- calculate effective radius of cloud water, cloud ice, and snow:
    if(has_reqc.ne.0 .and. has_reqi.ne.0 .and. has_reqs.ne.0) then
       do k = kts,kte
          re_qc1d(k) = 2.49E-6
          re_qi1d(k) = 4.99E-6
          re_qs1d(k) = 9.99E-6
       enddo
       call calc_effectRad (t1d, p1d , qv1d    , qc1d    , nc1d , qi1d , ni1d , qs1d ,  &
                            re_qc1d  , re_qi1d , re_qs1d , kts  , kte)
       do k = kts, kte
          re_cloud(i,k) = MAX(2.49E-6, MIN(re_qc1d(k),50.E-6))
          re_ice(i,k)   = MAX(4.99E-6, MIN(re_qi1d(k),125.E-6))
          re_snow(i,k)  = MAX(9.99E-6, MIN(re_qs1d(k),999.E-6))
       enddo
    endif

 enddo

 end subroutine mp_thompson_run

!+---+-----------------------------------------------------------------+
!ctrlL
!+---+-----------------------------------------------------------------+
!+---+-----------------------------------------------------------------+
!.. This subroutine computes the moisture tendencies of water vapor,
!.. cloud droplets, rain, cloud ice (pristine), snow, and graupel.
!.. Previously this code was based on Reisner et al (1998), but few of
!.. those pieces remain.  A complete description is now found in
!.. Thompson et al. (2004, 2008).
!+---+-----------------------------------------------------------------+
!
      subroutine thompson_1d(qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d,   &
                         nr1d, nc1d, nwfa1d, nifa1d, t1d, p1d, w1d, dzq, &
                         pptrain, pptsnow, pptgraul, pptice,             &
                         rainprod, evapprod,                             &
                         ppfr,ppfi,ppfs,ppfg,                            &
                         kts, kte, dt)

      implicit none

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: &
                          qv1d, qc1d, qi1d, qr1d, qs1d, qg1d, ni1d, &
                          nr1d, nc1d, nwfa1d, nifa1d, t1d
      REAL, DIMENSION(kts:kte), INTENT(IN):: p1d, w1d, dzq
      REAL, INTENT(INOUT):: pptrain, pptsnow, pptgraul, pptice
      REAL, INTENT(IN):: dt
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: &
                          rainprod, evapprod
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: &
                          ppfr,ppfi,ppfs,ppfg

!..Local variables
      REAL, DIMENSION(kts:kte):: tten, qvten, qcten, qiten, &
           qrten, qsten, qgten, niten, nrten, ncten, nwfaten, nifaten

      DOUBLE PRECISION, DIMENSION(kts:kte):: prw_vcd

      DOUBLE PRECISION, DIMENSION(kts:kte):: pnc_wcd, pnc_wau, pnc_rcw, &
           pnc_scw, pnc_gcw

      DOUBLE PRECISION, DIMENSION(kts:kte):: pna_rca, pna_sca, pna_gca, &
           pnd_rcd, pnd_scd, pnd_gcd

      DOUBLE PRECISION, DIMENSION(kts:kte):: prr_wau, prr_rcw, prr_rcs, &
           prr_rcg, prr_sml, prr_gml, &
           prr_rci, prv_rev,          &
           pnr_wau, pnr_rcs, pnr_rcg, &
           pnr_rci, pnr_sml, pnr_gml, &
           pnr_rev, pnr_rcr, pnr_rfz

      DOUBLE PRECISION, DIMENSION(kts:kte):: pri_inu, pni_inu, pri_ihm, &
           pni_ihm, pri_wfz, pni_wfz, &
           pri_rfz, pni_rfz, pri_ide, &
           pni_ide, pri_rci, pni_rci, &
           pni_sci, pni_iau, pri_iha, pni_iha

      DOUBLE PRECISION, DIMENSION(kts:kte):: prs_iau, prs_sci, prs_rcs, &
           prs_scw, prs_sde, prs_ihm, &
           prs_ide

      DOUBLE PRECISION, DIMENSION(kts:kte):: prg_scw, prg_rfz, prg_gde, &
           prg_gcw, prg_rci, prg_rcs, &
           prg_rcg, prg_ihm

      DOUBLE PRECISION, PARAMETER:: zeroD0 = 0.0d0

      REAL, DIMENSION(kts:kte):: temp, pres, qv
      REAL, DIMENSION(kts:kte):: rc, ri, rr, rs, rg, ni, nr, nc, nwfa, nifa
      REAL, DIMENSION(kts:kte):: rho, rhof, rhof2
      REAL, DIMENSION(kts:kte):: qvs, qvsi, delQvs
      REAL, DIMENSION(kts:kte):: satw, sati, ssatw, ssati
      REAL, DIMENSION(kts:kte):: diffu, visco, vsc2, &
           tcond, lvap, ocp, lvt2

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
      REAL, DIMENSION(kts:kte):: mvd_r, mvd_c
      REAL, DIMENSION(kts:kte):: smob, smo2, smo1, smo0, &
           smoc, smod, smoe, smof

      REAL, DIMENSION(kts:kte):: sed_r, sed_s, sed_g, sed_i, sed_n,sed_c

      REAL:: rgvm, delta_tp, orho, lfus2
      REAL, DIMENSION(5):: onstep
      DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamc, lamr, lamg
      DOUBLE PRECISION:: lami, ilami, ilamc
      REAL:: xDc, Dc_b, Dc_g, xDi, xDr, xDs, xDg, Ds_m, Dg_m
      DOUBLE PRECISION:: Dr_star, Dc_star
      REAL:: zeta1, zeta, taud, tau
      REAL:: stoke_r, stoke_s, stoke_g, stoke_i
      REAL:: vti, vtr, vts, vtg, vtc
      REAL, DIMENSION(kts:kte+1):: vtik, vtnik, vtrk, vtnrk, vtsk, vtgk,  &
           vtck, vtnck
      REAL, DIMENSION(kts:kte):: vts_boost
      REAL:: Mrat, ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts, C_snow
      REAL:: a_, b_, loga_, A1, A2, tf
      REAL:: tempc, tc0, r_mvd1, r_mvd2, xkrat
      REAL:: xnc, xri, xni, xmi, oxmi, xrc, xrr, xnr
      REAL:: xsat, rate_max, sump, ratio
      REAL:: clap, fcd, dfcd
      REAL:: otemp, rvs, rvs_p, rvs_pp, gamsc, alphsc, t1_evap, t1_subl
      REAL:: r_frac, g_frac
      REAL:: Ef_rw, Ef_sw, Ef_gw, Ef_rr
      REAL:: Ef_ra, Ef_sa, Ef_ga
      REAL:: dtsave, odts, odt, odzq, hgt_agl, SR
      REAL:: xslw1, ygra1, zans1, eva_factor
      INTEGER:: i, k, k2, n, nn, nstep, k_0, kbot, IT, iexfrq
      INTEGER, DIMENSION(5):: ksed1
      INTEGER:: nir, nis, nig, nii, nic, niin
      INTEGER:: idx_tc, idx_t, idx_s, idx_g1, idx_g, idx_r1, idx_r,     &
                idx_i1, idx_i, idx_c, idx, idx_d, idx_n, idx_in

      LOGICAL:: melti, no_micro
      LOGICAL, DIMENSION(kts:kte):: L_qc, L_qi, L_qr, L_qs, L_qg
      LOGICAL:: debug_flag
      CHARACTER*256:: mp_debug
      INTEGER:: nu_c

!     modifications proposed by Ted Mansell for MPAS.
!     Laura D. Fowler (laura@ucar.edu) / 2017-03-27.
!     real, parameter:: mvd_r_breakup = 1.e-3
!...  end modifications.
      LOGICAL, DIMENSION(kts:kte):: L_nifa,L_nwfa
      REAL:: tem
!+---+

      debug_flag = .false.
!     if (ii.eq.901 .and. jj.eq.379) debug_flag = .true.
      if(debug_flag) then
!       write(mp_debug, *) 'DEBUG INFO, mp_thompson at (i,j) ', ii, ', ', jj
!       CALL wrf_debug(550, mp_debug)
      endif

      no_micro = .true.
      dtsave = dt
      odt = 1./dt
      odts = 1./dtsave
      iexfrq = 1

!+---+-----------------------------------------------------------------+
!.. Source/sink terms.  First 2 chars: "pr" represents source/sink of
!.. mass while "pn" represents source/sink of number.  Next char is one
!.. of "v" for water vapor, "r" for rain, "i" for cloud ice, "w" for
!.. cloud water, "s" for snow, and "g" for graupel.  Next chars
!.. represent processes: "de" for sublimation/deposition, "ev" for
!.. evaporation, "fz" for freezing, "ml" for melting, "au" for
!.. autoconversion, "nu" for ice nucleation, "hm" for Hallet/Mossop
!.. secondary ice production, and "c" for collection followed by the
!.. character for the species being collected.  ALL of these terms are
!.. positive (except for deposition/sublimation terms which can switch
!.. signs based on super/subsaturation) and are treated as negatives
!.. where necessary in the tendency equations.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         tten(k) = 0.
         qvten(k) = 0.
         qcten(k) = 0.
         qiten(k) = 0.
         qrten(k) = 0.
         qsten(k) = 0.
         qgten(k) = 0.
         niten(k) = 0.
         nrten(k) = 0.
         ncten(k) = 0.
         nwfaten(k) = 0.
         nifaten(k) = 0.

         prw_vcd(k) = 0.

         pnc_wcd(k) = 0.
         pnc_wau(k) = 0.
         pnc_rcw(k) = 0.
         pnc_scw(k) = 0.
         pnc_gcw(k) = 0.

         prv_rev(k) = 0.
         prr_wau(k) = 0.
         prr_rcw(k) = 0.
         prr_rcs(k) = 0.
         prr_rcg(k) = 0.
         prr_sml(k) = 0.
         prr_gml(k) = 0.
         prr_rci(k) = 0.
         pnr_wau(k) = 0.
         pnr_rcs(k) = 0.
         pnr_rcg(k) = 0.
         pnr_rci(k) = 0.
         pnr_sml(k) = 0.
         pnr_gml(k) = 0.
         pnr_rev(k) = 0.
         pnr_rcr(k) = 0.
         pnr_rfz(k) = 0.

         pri_inu(k) = 0.
         pni_inu(k) = 0.
         pri_ihm(k) = 0.
         pni_ihm(k) = 0.
         pri_wfz(k) = 0.
         pni_wfz(k) = 0.
         pri_rfz(k) = 0.
         pni_rfz(k) = 0.
         pri_ide(k) = 0.
         pni_ide(k) = 0.
         pri_rci(k) = 0.
         pni_rci(k) = 0.
         pni_sci(k) = 0.
         pni_iau(k) = 0.
         pri_iha(k) = 0.
         pni_iha(k) = 0.

         prs_iau(k) = 0.
         prs_sci(k) = 0.
         prs_rcs(k) = 0.
         prs_scw(k) = 0.
         prs_sde(k) = 0.
         prs_ihm(k) = 0.
         prs_ide(k) = 0.

         prg_scw(k) = 0.
         prg_rfz(k) = 0.
         prg_gde(k) = 0.
         prg_gcw(k) = 0.
         prg_rci(k) = 0.
         prg_rcs(k) = 0.
         prg_rcg(k) = 0.
         prg_ihm(k) = 0.

         pna_rca(k) = 0.
         pna_sca(k) = 0.
         pna_gca(k) = 0.

         pnd_rcd(k) = 0.
         pnd_scd(k) = 0.
         pnd_gcd(k) = 0.
      enddo
      do k = kts, kte
         rainprod(k) = 0.
         evapprod(k) = 0.
         ppfr(k)     = 0.
         ppfi(k)     = 0.
         ppfs(k)     = 0.
         ppfg(k)     = 0.
      enddo
!.. initialize the logicals L_nifa and L_nwfa used to detect instances of the cloud
!.. ice and cloud liquid water mixing ratios being greater than R1 but their number
!.. concentration being less than 2. and R2:
      if(is_aerosol_aware) then
         do k = kts, kte
            L_nifa(k) = .false.
            L_nwfa(k) = .false.
         enddo
      endif

!..Bug fix (2016Jun15), prevent use of uninitialized value(s) of snow moments.
      do k = kts, kte
         smo0(k) = 0.
         smo1(k) = 0.
         smo2(k) = 0.
         smob(k) = 0.
         smoc(k) = 0.
         smod(k) = 0.
         smoe(k) = 0.
         smof(k) = 0.
         mvd_r(k) = 0.
         mvd_c(k) = 0.
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         nwfa(k) = MAX(11.1E6, MIN(9999.E6, nwfa1d(k)*rho(k)))
         nifa(k) = MAX(naIN1*0.01, MIN(9999.E6, nifa1d(k)*rho(k)))

         if (qc1d(k) .gt. R1) then
            no_micro = .false.
            rc(k) = qc1d(k)*rho(k)
            nc(k) = MAX(2., MIN(nc1d(k)*rho(k), Nt_c_max))
            L_qc(k) = .true.
!.. set L_nwfa to true when the cloud liquid water number concentration is less than 2.:
            if(is_aerosol_aware .and. nc(k) .le. 2.) L_nwfa(k) = .true.
            nu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
            lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
            xDc = (bm_r + nu_c + 1.) / lamc
            if (xDc.lt. D0c) then
             lamc = cce(2,nu_c)/D0c
            elseif (xDc.gt. D0r*2.) then
             lamc = cce(2,nu_c)/(D0r*2.)
            endif
            nc(k) = MIN( DBLE(Nt_c_max), ccg(1,nu_c)*ocg2(nu_c)*rc(k)   &
                  / am_r*lamc**bm_r)
            if (.NOT. is_aerosol_aware) nc(k) = Nt_c
         else
            qc1d(k) = 0.0
            nc1d(k) = 0.0
            rc(k) = R1
            nc(k) = 2.
            L_qc(k) = .false.
         endif

         if (qi1d(k) .gt. R1) then
            no_micro = .false.
            ri(k) = qi1d(k)*rho(k)
            ni(k) = MAX(R2, ni1d(k)*rho(k))
            L_qi(k) = .true.
!.. set L_nifa to true when the cloud ice number concentration is less than R2:
            if(is_aerosol_aware .and. ni(k) .le. R2) L_nifa(k) = .true.
            if (ni(k).le. R2) then
               lami = cie(2)/5.E-6
               ni(k) = MIN(9999.D3, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
            endif
!           L_qi(k) = .true.
            lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
            ilami = 1./lami
            xDi = (bm_i + mu_i + 1.) * ilami
            if (xDi.lt. 5.E-6) then
             lami = cie(2)/5.E-6
             ni(k) = MIN(9999.D3, cig(1)*oig2*ri(k)/am_i*lami**bm_i)
            elseif (xDi.gt. 300.E-6) then
             lami = cie(2)/300.E-6
             ni(k) = cig(1)*oig2*ri(k)/am_i*lami**bm_i
            endif
         else
            qi1d(k) = 0.0
            ni1d(k) = 0.0
            ri(k) = R1
            ni(k) = R2
            L_qi(k) = .false.
         endif

         if (qr1d(k) .gt. R1) then
            no_micro = .false.
            rr(k) = qr1d(k)*rho(k)
            nr(k) = MAX(R2, nr1d(k)*rho(k))
            if (nr(k).le. R2) then
               mvd_r(k) = 1.0E-3
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            endif
            L_qr(k) = .true.
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
            if (mvd_r(k) .gt. 2.5E-3) then
               mvd_r(k) = 2.5E-3
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            elseif (mvd_r(k) .lt. D0r*0.75) then
               mvd_r(k) = D0r*0.75
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            endif
         else
            qr1d(k) = 0.0
            nr1d(k) = 0.0
            rr(k) = R1
            nr(k) = R2
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R1) then
            no_micro = .false.
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            qs1d(k) = 0.0
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R1) then
            no_micro = .false.
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            qg1d(k) = 0.0
            rg(k) = R1
            L_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!     if (debug_flag) then
!      write(mp_debug,*) 'DEBUG-VERBOSE at (i,j) ', ii, ', ', jj
!      CALL wrf_debug(550, mp_debug)
!      do k = kts, kte
!        write(mp_debug, '(a,i3,f8.2,1x,f7.2,1x, 11(1x,e13.6))')        &
!    &              'VERBOSE: ', k, pres(k)*0.01, temp(k)-273.15, qv(k), rc(k), rr(k), ri(k), rs(k), rg(k), nc(k), nr(k), ni(k), nwfa(k), nifa(k)
!        CALL wrf_debug(550, mp_debug)
!      enddo
!     endif
!+---+-----------------------------------------------------------------+

!+---+-----------------------------------------------------------------+
!..Derive various thermodynamic variables frequently used.
!.. Saturation vapor pressure (mixing ratio) over liquid/ice comes from
!.. Flatau et al. 1992; enthalpy (latent heat) of vaporization from
!.. Bohren & Albrecht 1998; others from Pruppacher & Klett 1978.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         tempc = temp(k) - 273.15
         rhof(k) = SQRT(RHO_NOT/rho(k))
         rhof2(k) = SQRT(rhof(k))
         qvs(k) = rslf(pres(k), temp(k))
         delQvs(k) = MAX(0.0, rslf(pres(k), 273.15)-qv(k))
         if (tempc .le. 0.0) then
          qvsi(k) = rsif(pres(k), temp(k))
         else
          qvsi(k) = qvs(k)
         endif
         satw(k) = qv(k)/qvs(k)
         sati(k) = qv(k)/qvsi(k)
         ssatw(k) = satw(k) - 1.
         ssati(k) = sati(k) - 1.
         if (abs(ssatw(k)).lt. eps) ssatw(k) = 0.0
         if (abs(ssati(k)).lt. eps) ssati(k) = 0.0
         if (no_micro .and. ssati(k).gt. 0.0) no_micro = .false.
         diffu(k) = 2.11E-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
         if (tempc .ge. 0.0) then
            visco(k) = (1.718+0.0049*tempc)*1.0E-5
         else
            visco(k) = (1.718+0.0049*tempc-1.2E-5*tempc*tempc)*1.0E-5
         endif
         ocp(k) = 1./(Cp*(1.+0.887*qv(k)))
         vsc2(k) = SQRT(rho(k)/visco(k))
         lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
         tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
      enddo

!+---+-----------------------------------------------------------------+
!..If no existing hydrometeor species and no chance to initiate ice or
!.. condense cloud water, just exit quickly!
!+---+-----------------------------------------------------------------+

      if (no_micro) return

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         if (.not. L_qs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
               + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
               + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
               + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
               + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
               + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
               + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
               + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
               + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif

!..Calculate 0th moment.  Represents snow number concentration.
         loga_ = sa(1) + sa(2)*tc0 + sa(5)*tc0*tc0 + sa(9)*tc0*tc0*tc0
         a_ = 10.0**loga_
         b_ = sb(1) + sb(2)*tc0 + sb(5)*tc0*tc0 + sb(9)*tc0*tc0*tc0
         smo0(k) = a_ * smo2(k)**b_

!..Calculate 1st moment.  Useful for depositional growth and melting.
         loga_ = sa(1) + sa(2)*tc0 + sa(3) &
               + sa(4)*tc0 + sa(5)*tc0*tc0 &
               + sa(6) + sa(7)*tc0*tc0 &
               + sa(8)*tc0 + sa(9)*tc0*tc0*tc0 &
               + sa(10)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3) + sb(4)*tc0 &
              + sb(5)*tc0*tc0 + sb(6) &
              + sb(7)*tc0*tc0 + sb(8)*tc0 &
              + sb(9)*tc0*tc0*tc0 + sb(10)
         smo1(k) = a_ * smo2(k)**b_

!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
              + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
              + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_

!..Calculate bv_s+2 (th) moment.  Useful for riming.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(13) &
               + sa(4)*tc0*cse(13) + sa(5)*tc0*tc0 &
               + sa(6)*cse(13)*cse(13) + sa(7)*tc0*tc0*cse(13) &
               + sa(8)*tc0*cse(13)*cse(13) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(13)*cse(13)*cse(13)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(13) + sb(4)*tc0*cse(13) &
              + sb(5)*tc0*tc0 + sb(6)*cse(13)*cse(13) &
              + sb(7)*tc0*tc0*cse(13) + sb(8)*tc0*cse(13)*cse(13) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(13)*cse(13)*cse(13)
         smoe(k) = a_ * smo2(k)**b_

!..Calculate 1+(bv_s+1)/2 (th) moment.  Useful for depositional growth.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(16) &
               + sa(4)*tc0*cse(16) + sa(5)*tc0*tc0 &
               + sa(6)*cse(16)*cse(16) + sa(7)*tc0*tc0*cse(16) &
               + sa(8)*tc0*cse(16)*cse(16) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(16)*cse(16)*cse(16)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(16) + sb(4)*tc0*cse(16) &
              + sb(5)*tc0*tc0 + sb(6)*cse(16)*cse(16) &
              + sb(7)*tc0*tc0*cse(16) + sb(8)*tc0*cse(16)*cse(16) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(16)*cse(16)*cse(16)
         smof(k) = a_ * smo2(k)**b_

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+
      N0_min = gonv_max
      k_0 = kts
      do k = kte, kts, -1
         if (temp(k).ge.270.65) k_0 = MAX(k_0, k)
      enddo
      do k = kte, kts, -1
         if (k.gt.k_0 .and. L_qr(k) .and. mvd_r(k).gt.100.E-6) then
            xslw1 = 4.01 + alog10(mvd_r(k))
         else
            xslw1 = 0.01
         endif
         ygra1 = 4.31 + alog10(max(5.E-5, rg(k)))
         zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         N0_min = MIN(N0_exp, N0_min)
         N0_exp = N0_min
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo

      endif

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for rain.
!+---+-----------------------------------------------------------------+
      do k = kte, kts, -1
         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
         ilamr(k) = 1./lamr
         mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
         N0_r(k) = nr(k)*org2*lamr**cre(2)
      enddo

!+---+-----------------------------------------------------------------+
!..Compute warm-rain process terms (except evap done later).
!+---+-----------------------------------------------------------------+

      do k = kts, kte

!..Rain self-collection follows Seifert, 1994 and drop break-up
!.. follows Verlinde and Cotton, 1993.                                        RAIN2M
         if (L_qr(k) .and. mvd_r(k).gt. D0r) then
!-GT      Ef_rr = 1.0
!-GT      if (mvd_r(k) .gt. 1500.0E-6) then
             Ef_rr = 1.0 - EXP(2300.0*(mvd_r(k)-1950.0E-6))
!-GT      endif
          pnr_rcr(k) = Ef_rr * 2.0*nr(k)*rr(k)
         endif

         mvd_c(k) = D0c
         if (L_qc(k)) then
          nu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
          xDc = MAX(D0c*1.E6, ((rc(k)/(am_r*nc(k)))**obmr) * 1.E6)
          lamc = (nc(k)*am_r* ccg(2,nu_c) * ocg1(nu_c) / rc(k))**obmr
          mvd_c(k) = (3.0+nu_c+0.672) / lamc
         endif

!..Autoconversion follows Berry & Reinhardt (1974) with characteristic
!.. diameters correctly computed from gamma distrib of cloud droplets.
         if (rc(k).gt. 0.01e-3) then
          Dc_g = ((ccg(3,nu_c)*ocg2(nu_c))**obmr / lamc) * 1.E6
          Dc_b = (xDc*xDc*xDc*Dc_g*Dc_g*Dc_g - xDc*xDc*xDc*xDc*xDc*xDc) &
                 **(1./6.)
          zeta1 = 0.5*((6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4) &
                     + abs(6.25E-6*xDc*Dc_b*Dc_b*Dc_b - 0.4))
          zeta = 0.027*rc(k)*zeta1
          taud = 0.5*((0.5*Dc_b - 7.5) + abs(0.5*Dc_b - 7.5)) + R1
          tau  = 3.72/(rc(k)*taud)
          prr_wau(k) = zeta/tau
          prr_wau(k) = MIN(DBLE(rc(k)*odts), prr_wau(k))
          pnr_wau(k) = prr_wau(k) / (am_r*nu_c*200.*D0r*D0r*D0r)         ! RAIN2M
          pnc_wau(k) = MIN(DBLE(nc(k)*odts), prr_wau(k)                 &
                     / (am_r*mvd_c(k)*mvd_c(k)*mvd_c(k)))                   ! Qc2M
         endif

!..Rain collecting cloud water.  In CE, assume Dc<<Dr and vtc=~0.
         if (L_qr(k) .and. mvd_r(k).gt. D0r .and. mvd_c(k).gt. D0c) then
          lamr = 1./ilamr(k)
          idx = 1 + INT(nbr*DLOG(mvd_r(k)/Dr(1))/DLOG(Dr(nbr)/Dr(1)))
          idx = MIN(idx, nbr)
          Ef_rw = t_Efrw(idx, INT(mvd_c(k)*1.E6))
          prr_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*rc(k)*N0_r(k) &
                         *((lamr+fv_r)**(-cre(9)))
          prr_rcw(k) = MIN(DBLE(rc(k)*odts), prr_rcw(k))
          pnc_rcw(k) = rhof(k)*t1_qr_qc*Ef_rw*nc(k)*N0_r(k)             &
                         *((lamr+fv_r)**(-cre(9)))                          ! Qc2M
          pnc_rcw(k) = MIN(DBLE(nc(k)*odts), pnc_rcw(k))
         endif

!..Rain collecting aerosols, wet scavenging.
         if (L_qr(k) .and. mvd_r(k).gt. D0r) then
          Ef_ra = Eff_aero(mvd_r(k),0.04E-6,visco(k),rho(k),temp(k),'r')
          lamr = 1./ilamr(k)
          pna_rca(k) = rhof(k)*t1_qr_qc*Ef_ra*nwfa(k)*N0_r(k)           &
                         *((lamr+fv_r)**(-cre(9)))
          pna_rca(k) = MIN(DBLE(nwfa(k)*odts), pna_rca(k))

          Ef_ra = Eff_aero(mvd_r(k),0.8E-6,visco(k),rho(k),temp(k),'r')
          pnd_rcd(k) = rhof(k)*t1_qr_qc*Ef_ra*nifa(k)*N0_r(k)           &
                         *((lamr+fv_r)**(-cre(9)))
          pnd_rcd(k) = MIN(DBLE(nifa(k)*odts), pnd_rcd(k))
         endif

      enddo

!+---+-----------------------------------------------------------------+
!..Compute all frozen hydrometeor species' process terms.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      !..vts_boost is the factor applied to snow terminal
      !..fallspeed due to riming of snow
      do k = kts, kte
         vts_boost(k) = 1.0
         xDs = 0.0
         if (L_qs(k)) xDs = smoc(k) / smob(k)

!..Temperature lookup table indexes.
         tempc = temp(k) - 273.15
         idx_tc = MAX(1, MIN(NINT(-tempc), 45) )
         idx_t = INT( (tempc-2.5)/5. ) - 1
         idx_t = MAX(1, -idx_t)
         idx_t = MIN(idx_t, ntb_t)
         IT = MAX(1, MIN(NINT(-tempc), 31) )

!..Cloud water lookup table index.
         if (rc(k).gt. r_c(1)) then
          nic = NINT(ALOG10(rc(k)))
          do nn = nic-1, nic+1
             n = nn
             if ( (rc(k)/10.**nn).ge.1.0 .and. &
                  (rc(k)/10.**nn).lt.10.0) goto 141
          enddo
 141      continue
          idx_c = INT(rc(k)/10.**n) + 10*(n-nic2) - (n-nic2)
          idx_c = MAX(1, MIN(idx_c, ntb_c))
         else
          idx_c = 1
         endif

!..Cloud droplet number lookup table index.
         idx_n = NINT(1.0 + FLOAT(nbc) * DLOG(nc(k)/t_Nc(1)) / nic1)
         idx_n = MAX(1, MIN(idx_n, nbc))

!..Cloud ice lookup table indexes.
         if (ri(k).gt. r_i(1)) then
          nii = NINT(ALOG10(ri(k)))
          do nn = nii-1, nii+1
             n = nn
             if ( (ri(k)/10.**nn).ge.1.0 .and. &
                  (ri(k)/10.**nn).lt.10.0) goto 142
          enddo
 142      continue
          idx_i = INT(ri(k)/10.**n) + 10*(n-nii2) - (n-nii2)
          idx_i = MAX(1, MIN(idx_i, ntb_i))
         else
          idx_i = 1
         endif

         if (ni(k).gt. Nt_i(1)) then
          nii = NINT(ALOG10(ni(k)))
          do nn = nii-1, nii+1
             n = nn
             if ( (ni(k)/10.**nn).ge.1.0 .and. &
                  (ni(k)/10.**nn).lt.10.0) goto 143
          enddo
 143      continue
          idx_i1 = INT(ni(k)/10.**n) + 10*(n-nii3) - (n-nii3)
          idx_i1 = MAX(1, MIN(idx_i1, ntb_i1))
         else
          idx_i1 = 1
         endif

!..Rain lookup table indexes.
         if (rr(k).gt. r_r(1)) then
          nir = NINT(ALOG10(rr(k)))
          do nn = nir-1, nir+1
             n = nn
             if ( (rr(k)/10.**nn).ge.1.0 .and. &
                  (rr(k)/10.**nn).lt.10.0) goto 144
          enddo
 144      continue
          idx_r = INT(rr(k)/10.**n) + 10*(n-nir2) - (n-nir2)
          idx_r = MAX(1, MIN(idx_r, ntb_r))

          lamr = 1./ilamr(k)
          lam_exp = lamr * (crg(3)*org2*org1)**bm_r
          N0_exp = org1*rr(k)/am_r * lam_exp**cre(1)
          nir = NINT(DLOG10(N0_exp))
          do nn = nir-1, nir+1
             n = nn
             if ( (N0_exp/10.**nn).ge.1.0 .and. &
                  (N0_exp/10.**nn).lt.10.0) goto 145
          enddo
 145      continue
          idx_r1 = INT(N0_exp/10.**n) + 10*(n-nir3) - (n-nir3)
          idx_r1 = MAX(1, MIN(idx_r1, ntb_r1))
         else
          idx_r = 1
          idx_r1 = ntb_r1
         endif

!..Snow lookup table index.
         if (rs(k).gt. r_s(1)) then
          nis = NINT(ALOG10(rs(k)))
          do nn = nis-1, nis+1
             n = nn
             if ( (rs(k)/10.**nn).ge.1.0 .and. &
                  (rs(k)/10.**nn).lt.10.0) goto 146
          enddo
 146      continue
          idx_s = INT(rs(k)/10.**n) + 10*(n-nis2) - (n-nis2)
          idx_s = MAX(1, MIN(idx_s, ntb_s))
         else
          idx_s = 1
         endif

!..Graupel lookup table index.
         if (rg(k).gt. r_g(1)) then
          nig = NINT(ALOG10(rg(k)))
          do nn = nig-1, nig+1
             n = nn
             if ( (rg(k)/10.**nn).ge.1.0 .and. &
                  (rg(k)/10.**nn).lt.10.0) goto 147
          enddo
 147      continue
          idx_g = INT(rg(k)/10.**n) + 10*(n-nig2) - (n-nig2)
          idx_g = MAX(1, MIN(idx_g, ntb_g))

          lamg = 1./ilamg(k)
          lam_exp = lamg * (cgg(3)*ogg2*ogg1)**bm_g
          N0_exp = ogg1*rg(k)/am_g * lam_exp**cge(1)
          nig = NINT(DLOG10(N0_exp))
          do nn = nig-1, nig+1
             n = nn
             if ( (N0_exp/10.**nn).ge.1.0 .and. &
                  (N0_exp/10.**nn).lt.10.0) goto 148
          enddo
 148      continue
          idx_g1 = INT(N0_exp/10.**n) + 10*(n-nig3) - (n-nig3)
          idx_g1 = MAX(1, MIN(idx_g1, ntb_g1))
         else
          idx_g = 1
          idx_g1 = ntb_g1
         endif

!..Deposition/sublimation prefactor (from Srivastava & Coen 1992).
         otemp = 1./temp(k)
         rvs = rho(k)*qvsi(k)
         rvs_p = rvs*otemp*(lsub*otemp*oRv - 1.)
         rvs_pp = rvs * ( otemp*(lsub*otemp*oRv - 1.) &
                         *otemp*(lsub*otemp*oRv - 1.) &
                         + (-2.*lsub*otemp*otemp*otemp*oRv) &
                         + otemp*otemp)
         gamsc = lsub*diffu(k)/tcond(k) * rvs_p
         alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                    * rvs_pp/rvs_p * rvs/rvs_p
         alphsc = MAX(1.E-9, alphsc)
         xsat = ssati(k)
         if (abs(xsat).lt. 1.E-9) xsat=0.
         t1_subl = 4.*PI*( 1.0 - alphsc*xsat &
                + 2.*alphsc*alphsc*xsat*xsat &
                - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                / (1.+gamsc)

!..Snow collecting cloud water.  In CE, assume Dc<<Ds and vtc=~0.
         if (L_qc(k) .and. mvd_c(k).gt. D0c) then
          if (xDs .gt. D0s) then
           idx = 1 + INT(nbs*DLOG(xDs/Ds(1))/DLOG(Ds(nbs)/Ds(1)))
           idx = MIN(idx, nbs)
           Ef_sw = t_Efsw(idx, INT(mvd_c(k)*1.E6))
           prs_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*rc(k)*smoe(k)
           prs_scw(k) = MIN(DBLE(rc(k)*odts), prs_scw(k))
           pnc_scw(k) = rhof(k)*t1_qs_qc*Ef_sw*nc(k)*smoe(k)                ! Qc2M
           pnc_scw(k) = MIN(DBLE(nc(k)*odts), pnc_scw(k))
          endif

!..Graupel collecting cloud water.  In CE, assume Dc<<Dg and vtc=~0.
          if (rg(k).ge. r_g(1) .and. mvd_c(k).gt. D0c) then
           xDg = (bm_g + mu_g + 1.) * ilamg(k)
           vtg = rhof(k)*av_g*cgg(6)*ogg3 * ilamg(k)**bv_g
           stoke_g = mvd_c(k)*mvd_c(k)*vtg*rho_w/(9.*visco(k)*xDg)
           if (xDg.gt. D0g) then
            if (stoke_g.ge.0.4 .and. stoke_g.le.10.) then
             Ef_gw = 0.55*ALOG10(2.51*stoke_g)
            elseif (stoke_g.lt.0.4) then
             Ef_gw = 0.0
            elseif (stoke_g.gt.10) then
             Ef_gw = 0.77
            endif
            prg_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*rc(k)*N0_g(k) &
                          *ilamg(k)**cge(9)
            pnc_gcw(k) = rhof(k)*t1_qg_qc*Ef_gw*nc(k)*N0_g(k)           &
                          *ilamg(k)**cge(9)                                 ! Qc2M
            pnc_gcw(k) = MIN(DBLE(nc(k)*odts), pnc_gcw(k))
           endif
          endif
         endif

!..Snow and graupel collecting aerosols, wet scavenging.
         if (rs(k) .gt. r_s(1)) then
          Ef_sa = Eff_aero(xDs,0.04E-6,visco(k),rho(k),temp(k),'s')
          pna_sca(k) = rhof(k)*t1_qs_qc*Ef_sa*nwfa(k)*smoe(k)
          pna_sca(k) = MIN(DBLE(nwfa(k)*odts), pna_sca(k))

          Ef_sa = Eff_aero(xDs,0.8E-6,visco(k),rho(k),temp(k),'s')
          pnd_scd(k) = rhof(k)*t1_qs_qc*Ef_sa*nifa(k)*smoe(k)
          pnd_scd(k) = MIN(DBLE(nifa(k)*odts), pnd_scd(k))
         endif
         if (rg(k) .gt. r_g(1)) then
          xDg = (bm_g + mu_g + 1.) * ilamg(k)
          Ef_ga = Eff_aero(xDg,0.04E-6,visco(k),rho(k),temp(k),'g')
          pna_gca(k) = rhof(k)*t1_qg_qc*Ef_ga*nwfa(k)*N0_g(k)           &
                        *ilamg(k)**cge(9)
          pna_gca(k) = MIN(DBLE(nwfa(k)*odts), pna_gca(k))

          Ef_ga = Eff_aero(xDg,0.8E-6,visco(k),rho(k),temp(k),'g')
          pnd_gcd(k) = rhof(k)*t1_qg_qc*Ef_ga*nifa(k)*N0_g(k)           &
                        *ilamg(k)**cge(9)
          pnd_gcd(k) = MIN(DBLE(nifa(k)*odts), pnd_gcd(k))
         endif

!..Rain collecting snow.  Cannot assume Wisner (1972) approximation
!.. or Mizuno (1990) approach so we solve the CE explicitly and store
!.. results in lookup table.
         if (rr(k).ge. r_r(1)) then
          if (rs(k).ge. r_s(1)) then
           if (temp(k).lt.T_0) then
            prr_rcs(k) = -(tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                           + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                           + tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                           + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r))
            prs_rcs(k) = tmr_racs2(idx_s,idx_t,idx_r1,idx_r) &
                         + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r) &
                         - tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                         - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
            prg_rcs(k) = tmr_racs1(idx_s,idx_t,idx_r1,idx_r) &
                         + tcr_sacr1(idx_s,idx_t,idx_r1,idx_r) &
                         + tcs_racs1(idx_s,idx_t,idx_r1,idx_r) &
                         + tms_sacr1(idx_s,idx_t,idx_r1,idx_r)
            prr_rcs(k) = MAX(DBLE(-rr(k)*odts), prr_rcs(k))
            prs_rcs(k) = MAX(DBLE(-rs(k)*odts), prs_rcs(k))
            prg_rcs(k) = MIN(DBLE((rr(k)+rs(k))*odts), prg_rcs(k))
            pnr_rcs(k) = tnr_racs1(idx_s,idx_t,idx_r1,idx_r)            &   ! RAIN2M
                         + tnr_racs2(idx_s,idx_t,idx_r1,idx_r)          &
                         + tnr_sacr1(idx_s,idx_t,idx_r1,idx_r)          &
                         + tnr_sacr2(idx_s,idx_t,idx_r1,idx_r)
            pnr_rcs(k) = MIN(DBLE(nr(k)*odts), pnr_rcs(k))
           else
            prs_rcs(k) = -tcs_racs1(idx_s,idx_t,idx_r1,idx_r)           &
                         - tms_sacr1(idx_s,idx_t,idx_r1,idx_r)          &
                         + tmr_racs2(idx_s,idx_t,idx_r1,idx_r)          &
                         + tcr_sacr2(idx_s,idx_t,idx_r1,idx_r)
            prs_rcs(k) = MAX(DBLE(-rs(k)*odts), prs_rcs(k))
            prr_rcs(k) = -prs_rcs(k)
           endif
          endif

!..Rain collecting graupel.  Cannot assume Wisner (1972) approximation
!.. or Mizuno (1990) approach so we solve the CE explicitly and store
!.. results in lookup table.
          if (rg(k).ge. r_g(1)) then
           if (temp(k).lt.T_0) then
            prg_rcg(k) = tmr_racg(idx_g1,idx_g,idx_r1,idx_r) &
                         + tcr_gacr(idx_g1,idx_g,idx_r1,idx_r)
            prg_rcg(k) = MIN(DBLE(rr(k)*odts), prg_rcg(k))
            prr_rcg(k) = -prg_rcg(k)
            pnr_rcg(k) = tnr_racg(idx_g1,idx_g,idx_r1,idx_r)            &   ! RAIN2M
                         + tnr_gacr(idx_g1,idx_g,idx_r1,idx_r)
            pnr_rcg(k) = MIN(DBLE(nr(k)*odts), pnr_rcg(k))
           else
            prr_rcg(k) = tcg_racg(idx_g1,idx_g,idx_r1,idx_r)
            prr_rcg(k) = MIN(DBLE(rg(k)*odts), prr_rcg(k))
            prg_rcg(k) = -prr_rcg(k)
!..Put in explicit drop break-up due to collisions.
!...  modifications to drop break-up due to collision of rain collecting graupel, ad-hoc
!     change proposed by Ted Mansell for MPAS: basically reduces rain breakup in graupel
!     collisions when the mean volume diameter of rain is less than 1mm.
!     Laura D. Fowler (laura@ucar.edu) / 2017-03-27.
!           if ( mvd_r(k) > mvd_r_breakup ) then
              pnr_rcg(k) = -5.0*tnr_gacr(idx_g1,idx_g,idx_r1,idx_r)       ! RAIN2M
!           else
!             pnr_rcg(k) = -3.0*tnr_gacr(idx_g1,idx_g,idx_r1,idx_r)       ! RAIN2M
!           endif
           endif
          endif
         endif

!+---+-----------------------------------------------------------------+
!..Next IF block handles only those processes below 0C.
!+---+-----------------------------------------------------------------+

         if (temp(k).lt.T_0) then

          vts_boost(k) = 1.0
          rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999

!+---+---------------- BEGIN NEW ICE NUCLEATION -----------------------+
!..Freezing of supercooled water (rain or cloud) is influenced by dust
!.. but still using Bigg 1953 with a temperature adjustment of a few
!.. degrees depending on dust concentration.  A default value by way
!.. of idx_IN is 1.0 per Liter of air is used when dustyIce flag is
!.. false.  Next, a combination of deposition/condensation freezing
!.. using DeMott et al (2010) dust nucleation when water saturated or
!.. Phillips et al (2008) when below water saturation; else, without
!.. dustyIce flag, use the previous Cooper (1986) temperature-dependent
!.. value.  Lastly, allow homogeneous freezing of deliquesced aerosols
!.. following Koop et al. (2001, Nature).
!.. Implemented by T. Eidhammer and G. Thompson 2012Dec18
!+---+-----------------------------------------------------------------+

          if (dustyIce .AND. is_aerosol_aware) then
           xni = iceDeMott(tempc,qvs(k),qvs(k),qvsi(k),rho(k),nifa(k))
          else
           xni = 1.0 *1000.                                               ! Default is 1.0 per Liter
          endif

!..Ice nuclei lookup table index.
          if (xni.gt. Nt_IN(1)) then
           niin = NINT(ALOG10(xni))
           do nn = niin-1, niin+1
              n = nn
              if ( (xni/10.**nn).ge.1.0 .and. &
                   (xni/10.**nn).lt.10.0) goto 149
           enddo
 149       continue
           idx_IN = INT(xni/10.**n) + 10*(n-niin2) - (n-niin2)
           idx_IN = MAX(1, MIN(idx_IN, ntb_IN))
          else
           idx_IN = 1
          endif

!..Freezing of water drops into graupel/cloud ice (Bigg 1953).
          if (rr(k).gt. r_r(1)) then
           prg_rfz(k) = tpg_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
           pri_rfz(k) = tpi_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
           pni_rfz(k) = tni_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts
           pnr_rfz(k) = tnr_qrfz(idx_r,idx_r1,idx_tc,idx_IN)*odts          ! RAIN2M
           pnr_rfz(k) = MIN(DBLE(nr(k)*odts), pnr_rfz(k))
          elseif (rr(k).gt. R1 .and. temp(k).lt.HGFR) then
           pri_rfz(k) = rr(k)*odts
           pni_rfz(k) = nr(k)*odts                                         ! RAIN2M
          endif

          if (rc(k).gt. r_c(1)) then
           pri_wfz(k) = tpi_qcfz(idx_c,idx_n,idx_tc,idx_IN)*odts
           pri_wfz(k) = MIN(DBLE(rc(k)*odts), pri_wfz(k))
           pni_wfz(k) = tni_qcfz(idx_c,idx_n,idx_tc,idx_IN)*odts
           pni_wfz(k) = MIN(DBLE(nc(k)*odts), pri_wfz(k)/(2.*xm0i),     &
                                pni_wfz(k))
          elseif (rc(k).gt. R1 .and. temp(k).lt.HGFR) then
           pri_wfz(k) = rc(k)*odts
           pni_wfz(k) = nc(k)*odts
          endif

!..Deposition nucleation of dust/mineral from DeMott et al (2010)
!.. we may need to relax the temperature and ssati constraints.
          if ( (ssati(k).ge. 0.25) .or. (ssatw(k).gt. eps &
                                .and. temp(k).lt.253.15) ) then
           if (dustyIce .AND. is_aerosol_aware) then
            xnc = iceDeMott(tempc,qv(k),qvs(k),qvsi(k),rho(k),nifa(k))
           else
            xnc = MIN(250.E3, TNO*EXP(ATO*(T_0-temp(k))))
           endif
           xni = ni(k) + (pni_rfz(k)+pni_wfz(k))*dtsave
           pni_inu(k) = 0.5*(xnc-xni + abs(xnc-xni))*odts
           pri_inu(k) = MIN(DBLE(rate_max), xm0i*pni_inu(k))
           pni_inu(k) = pri_inu(k)/xm0i
          endif

!..Freezing of aqueous aerosols based on Koop et al (2001, Nature)
          xni = smo0(k)+ni(k) + (pni_rfz(k)+pni_wfz(k)+pni_inu(k))*dtsave
          if (is_aerosol_aware .AND. homogIce .AND. (xni.le.999.E3)     &
     &                .AND.(temp(k).lt.238).AND.(ssati(k).ge.0.4) ) then
            xnc = iceKoop(temp(k),qv(k),qvs(k),nwfa(k), dtsave)
            pni_iha(k) = xnc*odts
            pri_iha(k) = MIN(DBLE(rate_max), xm0i*0.1*pni_iha(k))
            pni_iha(k) = pri_iha(k)/(xm0i*0.1)
          endif
!+---+------------------ END NEW ICE NUCLEATION -----------------------+


!..Deposition/sublimation of cloud ice (Srivastava & Coen 1992).
          if (L_qi(k)) then
           lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
           ilami = 1./lami
           xDi = MAX(DBLE(D0i), (bm_i + mu_i + 1.) * ilami)
           xmi = am_i*xDi**bm_i
           oxmi = 1./xmi
           pri_ide(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                  *oig1*cig(5)*ni(k)*ilami

           if (pri_ide(k) .lt. 0.0) then
            pri_ide(k) = MAX(DBLE(-ri(k)*odts), pri_ide(k), DBLE(rate_max))
            pni_ide(k) = pri_ide(k)*oxmi
            pni_ide(k) = MAX(DBLE(-ni(k)*odts), pni_ide(k))
           else
            pri_ide(k) = MIN(pri_ide(k), DBLE(rate_max))
            prs_ide(k) = (1.0D0-tpi_ide(idx_i,idx_i1))*pri_ide(k)
            pri_ide(k) = tpi_ide(idx_i,idx_i1)*pri_ide(k)
           endif

!..Some cloud ice needs to move into the snow category.  Use lookup
!.. table that resulted from explicit bin representation of distrib.
           if ( (idx_i.eq. ntb_i) .or. (xDi.gt. 5.0*D0s) ) then
            prs_iau(k) = ri(k)*.99*odts
            pni_iau(k) = ni(k)*.95*odts
           elseif (xDi.lt. 0.1*D0s) then
            prs_iau(k) = 0.
            pni_iau(k) = 0.
           else
            prs_iau(k) = tps_iaus(idx_i,idx_i1)*odts
            prs_iau(k) = MIN(DBLE(ri(k)*.99*odts), prs_iau(k))
            pni_iau(k) = tni_iaus(idx_i,idx_i1)*odts
            pni_iau(k) = MIN(DBLE(ni(k)*.95*odts), pni_iau(k))
           endif
          endif

!..Deposition/sublimation of snow/graupel follows Srivastava & Coen
!.. (1992).
          if (L_qs(k)) then
           C_snow = C_sqrd + (tempc+1.5)*(C_cube-C_sqrd)/(-30.+1.5)
           C_snow = MAX(C_sqrd, MIN(C_snow, C_cube))
           prs_sde(k) = C_snow*t1_subl*diffu(k)*ssati(k)*rvs &
                        * (t1_qs_sd*smo1(k) &
                         + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
           if (prs_sde(k).lt. 0.) then
            prs_sde(k) = MAX(DBLE(-rs(k)*odts), prs_sde(k), DBLE(rate_max))
           else
            prs_sde(k) = MIN(prs_sde(k), DBLE(rate_max))
           endif
          endif

          if (L_qg(k) .and. ssati(k).lt. -eps) then
           prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
               * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10) &
               + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11))
           if (prg_gde(k).lt. 0.) then
            prg_gde(k) = MAX(DBLE(-rg(k)*odts), prg_gde(k), DBLE(rate_max))
           else
            prg_gde(k) = MIN(prg_gde(k), DBLE(rate_max))
           endif
          endif

!..Snow collecting cloud ice.  In CE, assume Di<<Ds and vti=~0.
          if (L_qi(k)) then
           lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
           ilami = 1./lami
           xDi = MAX(DBLE(D0i), (bm_i + mu_i + 1.) * ilami)
           xmi = am_i*xDi**bm_i
           oxmi = 1./xmi
           if (rs(k).ge. r_s(1)) then
            prs_sci(k) = t1_qs_qi*rhof(k)*Ef_si*ri(k)*smoe(k)
            pni_sci(k) = prs_sci(k) * oxmi
           endif

!..Rain collecting cloud ice.  In CE, assume Di<<Dr and vti=~0.
           if (rr(k).ge. r_r(1) .and. mvd_r(k).gt. 4.*xDi) then
            lamr = 1./ilamr(k)
            pri_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ri(k)*N0_r(k) &
                           *((lamr+fv_r)**(-cre(9)))
            pnr_rci(k) = rhof(k)*t1_qr_qi*Ef_ri*ni(k)*N0_r(k)           &   ! RAIN2M
                           *((lamr+fv_r)**(-cre(9)))
            pni_rci(k) = pri_rci(k) * oxmi
            prr_rci(k) = rhof(k)*t2_qr_qi*Ef_ri*ni(k)*N0_r(k) &
                           *((lamr+fv_r)**(-cre(8)))
            prr_rci(k) = MIN(DBLE(rr(k)*odts), prr_rci(k))
            prg_rci(k) = pri_rci(k) + prr_rci(k)
           endif
          endif

!..Ice multiplication from rime-splinters (Hallet & Mossop 1974).
          if (prg_gcw(k).gt. eps .and. tempc.gt.-8.0) then
           tf = 0.
           if (tempc.ge.-5.0 .and. tempc.lt.-3.0) then
            tf = 0.5*(-3.0 - tempc)
           elseif (tempc.gt.-8.0 .and. tempc.lt.-5.0) then
            tf = 0.33333333*(8.0 + tempc)
           endif
           pni_ihm(k) = 3.5E8*tf*prg_gcw(k)
           pri_ihm(k) = xm0i*pni_ihm(k)
           prs_ihm(k) = prs_scw(k)/(prs_scw(k)+prg_gcw(k)) &
                          * pri_ihm(k)
           prg_ihm(k) = prg_gcw(k)/(prs_scw(k)+prg_gcw(k)) &
                          * pri_ihm(k)
          endif

!..A portion of rimed snow converts to graupel but some remains snow.
!.. Interp from 15 to 95% as riming factor increases from 2.0 to 30.0
!.. 0.028 came from (.95-.15)/(30.-2.).  This remains ad-hoc and should
!.. be revisited.
          if (prs_scw(k).gt.2.0*prs_sde(k) .and. &
                         prs_sde(k).gt.eps) then
           r_frac = MIN(30.0D0, prs_scw(k)/prs_sde(k))
           g_frac = MIN(0.95, 0.15 + (r_frac-2.)*.028)
           vts_boost(k) = MIN(1.5, 1.1 + (r_frac-2.)*.014)
           prg_scw(k) = g_frac*prs_scw(k)
           prs_scw(k) = (1. - g_frac)*prs_scw(k)
          endif

         else

!..Melt snow and graupel and enhance from collisions with liquid.
!.. We also need to sublimate snow and graupel if subsaturated.
          if (L_qs(k)) then
           prr_sml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                      * (t1_qs_me*smo1(k) + t2_qs_me*rhof2(k)*vsc2(k)*smof(k))
           if (prr_sml(k) .gt. 0.) then
              prr_sml(k) = prr_sml(k) + 4218.*olfus*tempc               &
                                      * (prr_rcs(k)+prs_scw(k))
           endif
           prr_sml(k) = MIN(DBLE(rs(k)*odts), MAX(0.D0, prr_sml(k)))
           pnr_sml(k) = smo0(k)/rs(k)*prr_sml(k) * 10.0**(-0.25*tempc)      ! RAIN2M
           pnr_sml(k) = MIN(DBLE(smo0(k)*odts), pnr_sml(k))

           if (ssati(k).lt. 0.) then
            prs_sde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                         * (t1_qs_sd*smo1(k) &
                          + t2_qs_sd*rhof2(k)*vsc2(k)*smof(k))
            prs_sde(k) = MAX(DBLE(-rs(k)*odts), prs_sde(k))
           endif
          endif

          if (L_qg(k)) then
           prr_gml(k) = (tempc*tcond(k)-lvap0*diffu(k)*delQvs(k))       &
                      * N0_g(k)*(t1_qg_me*ilamg(k)**cge(10)             &
                      + t2_qg_me*rhof2(k)*vsc2(k)*ilamg(k)**cge(11))
!-GT       prr_gml(k) = prr_gml(k) + 4218.*olfus*tempc &
!-GT                               * (prr_rcg(k)+prg_gcw(k))
           prr_gml(k) = MIN(DBLE(rg(k)*odts), MAX(0.D0, prr_gml(k)))
           pnr_gml(k) = N0_g(k)*cgg(2)*ilamg(k)**cge(2) / rg(k)         &   ! RAIN2M
                      * prr_gml(k) * 10.0**(-0.5*tempc)

           if (ssati(k).lt. 0.) then
            prg_gde(k) = C_cube*t1_subl*diffu(k)*ssati(k)*rvs &
                * N0_g(k) * (t1_qg_sd*ilamg(k)**cge(10) &
                + t2_qg_sd*vsc2(k)*rhof2(k)*ilamg(k)**cge(11))
            prg_gde(k) = MAX(DBLE(-rg(k)*odts), prg_gde(k))
           endif
          endif

!.. This change will be required if users run adaptive time step that
!.. results in delta-t that is generally too long to allow cloud water
!.. collection by snow/graupel above melting temperature.
!.. Credit to Bjorn-Egil Nygaard for discovering.
          if (dt .gt. 120.) then
             prr_rcw(k)=prr_rcw(k)+prs_scw(k)+prg_gcw(k)
             prs_scw(k)=0.
             prg_gcw(k)=0.
          endif

         endif

      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Ensure we do not deplete more hydrometeor species than exists.
!+---+-----------------------------------------------------------------+
      do k = kts, kte

!..If ice supersaturated, ensure sum of depos growth terms does not
!.. deplete more vapor than possibly exists.  If subsaturated, limit
!.. sum of sublimation terms such that vapor does not reproduce ice
!.. supersat again.
         sump = pri_inu(k) + pri_ide(k) + prs_ide(k) &
              + prs_sde(k) + prg_gde(k) + pri_iha(k)
         rate_max = (qv(k)-qvsi(k))*rho(k)*odts*0.999
         if ( (sump.gt. eps .and. sump.gt. rate_max) .or. &
              (sump.lt. -eps .and. sump.lt. rate_max) ) then
          ratio = rate_max/sump
          pri_inu(k) = pri_inu(k) * ratio
          pri_ide(k) = pri_ide(k) * ratio
          pni_ide(k) = pni_ide(k) * ratio
          prs_ide(k) = prs_ide(k) * ratio
          prs_sde(k) = prs_sde(k) * ratio
          prg_gde(k) = prg_gde(k) * ratio
          pri_iha(k) = pri_iha(k) * ratio
         endif

!..Cloud water conservation.
         sump = -prr_wau(k) - pri_wfz(k) - prr_rcw(k) &
                - prs_scw(k) - prg_scw(k) - prg_gcw(k)
         rate_max = -rc(k)*odts
         if (sump.lt. rate_max .and. L_qc(k)) then
          ratio = rate_max/sump
          prr_wau(k) = prr_wau(k) * ratio
          pri_wfz(k) = pri_wfz(k) * ratio
          prr_rcw(k) = prr_rcw(k) * ratio
          prs_scw(k) = prs_scw(k) * ratio
          prg_scw(k) = prg_scw(k) * ratio
          prg_gcw(k) = prg_gcw(k) * ratio
         endif

!..Cloud ice conservation.
         sump = pri_ide(k) - prs_iau(k) - prs_sci(k) &
                - pri_rci(k)
         rate_max = -ri(k)*odts
         if (sump.lt. rate_max .and. L_qi(k)) then
          ratio = rate_max/sump
          pri_ide(k) = pri_ide(k) * ratio
          prs_iau(k) = prs_iau(k) * ratio
          prs_sci(k) = prs_sci(k) * ratio
          pri_rci(k) = pri_rci(k) * ratio
         endif

!..Rain conservation.
         sump = -prg_rfz(k) - pri_rfz(k) - prr_rci(k) &
                + prr_rcs(k) + prr_rcg(k)
         rate_max = -rr(k)*odts
         if (sump.lt. rate_max .and. L_qr(k)) then
          ratio = rate_max/sump
          prg_rfz(k) = prg_rfz(k) * ratio
          pri_rfz(k) = pri_rfz(k) * ratio
          prr_rci(k) = prr_rci(k) * ratio
          prr_rcs(k) = prr_rcs(k) * ratio
          prr_rcg(k) = prr_rcg(k) * ratio
         endif

!..Snow conservation.
         sump = prs_sde(k) - prs_ihm(k) - prr_sml(k) &
                + prs_rcs(k)
         rate_max = -rs(k)*odts
         if (sump.lt. rate_max .and. L_qs(k)) then
          ratio = rate_max/sump
          prs_sde(k) = prs_sde(k) * ratio
          prs_ihm(k) = prs_ihm(k) * ratio
          prr_sml(k) = prr_sml(k) * ratio
          prs_rcs(k) = prs_rcs(k) * ratio
         endif

!..Graupel conservation.
         sump = prg_gde(k) - prg_ihm(k) - prr_gml(k) &
              + prg_rcg(k)
         rate_max = -rg(k)*odts
         if (sump.lt. rate_max .and. L_qg(k)) then
          ratio = rate_max/sump
          prg_gde(k) = prg_gde(k) * ratio
          prg_ihm(k) = prg_ihm(k) * ratio
          prr_gml(k) = prr_gml(k) * ratio
          prg_rcg(k) = prg_rcg(k) * ratio
         endif

!..Re-enforce proper mass conservation for subsequent elements in case
!.. any of the above terms were altered.  Thanks P. Blossey. 2009Sep28
         pri_ihm(k) = prs_ihm(k) + prg_ihm(k)
         ratio = MIN( ABS(prr_rcg(k)), ABS(prg_rcg(k)) )
         prr_rcg(k) = ratio * SIGN(1.0, SNGL(prr_rcg(k)))
         prg_rcg(k) = -prr_rcg(k)
         if (temp(k).gt.T_0) then
            ratio = MIN( ABS(prr_rcs(k)), ABS(prs_rcs(k)) )
            prr_rcs(k) = ratio * SIGN(1.0, SNGL(prr_rcs(k)))
            prs_rcs(k) = -prr_rcs(k)
         endif

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate tendencies of all species but constrain the number of ice
!.. to reasonable values.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         orho = 1./rho(k)
         lfus2 = lsub - lvap(k)

!..Aerosol number tendency
         if (is_aerosol_aware) then
            nwfaten(k) = nwfaten(k) - (pna_rca(k) + pna_sca(k)          &
                       + pna_gca(k) + pni_iha(k)) * orho
            nifaten(k) = nifaten(k) - (pnd_rcd(k) + pnd_scd(k)          &
                       + pnd_gcd(k)) * orho
            if (dustyIce) then
               nifaten(k) = nifaten(k) - pni_inu(k)*orho
            else
               nifaten(k) = 0.
            endif
         endif

!..Water vapor tendency
         qvten(k) = qvten(k) + (-pri_inu(k) - pri_iha(k) - pri_ide(k)   &
                      - prs_ide(k) - prs_sde(k) - prg_gde(k)) &
                      * orho

!..Cloud water tendency
         qcten(k) = qcten(k) + (-prr_wau(k) - pri_wfz(k) &
                      - prr_rcw(k) - prs_scw(k) - prg_scw(k) &
                      - prg_gcw(k)) &
                      * orho

!..Cloud water number tendency
         ncten(k) = ncten(k) + (-pnc_wau(k) - pnc_rcw(k) &
                      - pni_wfz(k) - pnc_scw(k) - pnc_gcw(k)) &
                      * orho

!..Cloud water mass/number balance; keep mass-wt mean size between
!.. 1 and 50 microns.  Also no more than Nt_c_max drops total.
         xrc=MAX(R1, (qc1d(k) + qcten(k)*dtsave)*rho(k))
         xnc=MAX(2., (nc1d(k) + ncten(k)*dtsave)*rho(k))
         if (xrc .gt. R1) then
          nu_c = MIN(15, NINT(1000.E6/xnc) + 2)
          lamc = (xnc*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
          xDc = (bm_r + nu_c + 1.) / lamc
          if (xDc.lt. D0c) then
           lamc = cce(2,nu_c)/D0c
           xnc = ccg(1,nu_c)*ocg2(nu_c)*xrc/am_r*lamc**bm_r
           ncten(k) = (xnc-nc1d(k)*rho(k))*odts*orho
          elseif (xDc.gt. D0r*2.) then
           lamc = cce(2,nu_c)/(D0r*2.)
           xnc = ccg(1,nu_c)*ocg2(nu_c)*xrc/am_r*lamc**bm_r
           ncten(k) = (xnc-nc1d(k)*rho(k))*odts*orho
          endif
         else
          ncten(k) = -nc1d(k)*odts
         endif
         xnc=MAX(0.,(nc1d(k) + ncten(k)*dtsave)*rho(k))
         if (xnc.gt.Nt_c_max) &
                ncten(k) = (Nt_c_max-nc1d(k)*rho(k))*odts*orho

!..Cloud ice mixing ratio tendency
         qiten(k) = qiten(k) + (pri_inu(k) + pri_iha(k) + pri_ihm(k)    &
                      + pri_wfz(k) + pri_rfz(k) + pri_ide(k) &
                      - prs_iau(k) - prs_sci(k) - pri_rci(k)) &
                      * orho

!..Cloud ice number tendency.
         niten(k) = niten(k) + (pni_inu(k) + pni_iha(k) + pni_ihm(k)    &
                      + pni_wfz(k) + pni_rfz(k) + pni_ide(k) &
                      - pni_iau(k) - pni_sci(k) - pni_rci(k)) &
                      * orho

!..Cloud ice mass/number balance; keep mass-wt mean size between
!.. 5 and 300 microns.  Also no more than 500 xtals per liter.
         xri=MAX(R1,(qi1d(k) + qiten(k)*dtsave)*rho(k))
         xni=MAX(R2,(ni1d(k) + niten(k)*dtsave)*rho(k))
         if (xri.gt. R1) then
           lami = (am_i*cig(2)*oig1*xni/xri)**obmi
           ilami = 1./lami
           xDi = (bm_i + mu_i + 1.) * ilami
           if (xDi.lt. 5.E-6) then
            lami = cie(2)/5.E-6
            xni = MIN(9999.D3, cig(1)*oig2*xri/am_i*lami**bm_i)
            niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
           elseif (xDi.gt. 300.E-6) then
            lami = cie(2)/300.E-6
            xni = cig(1)*oig2*xri/am_i*lami**bm_i
            niten(k) = (xni-ni1d(k)*rho(k))*odts*orho
           endif
         else
          niten(k) = -ni1d(k)*odts
         endif
         xni=MAX(0.,(ni1d(k) + niten(k)*dtsave)*rho(k))
         if (xni.gt.9999.E3) &
                niten(k) = (9999.E3-ni1d(k)*rho(k))*odts*orho

!..Rain tendency
         qrten(k) = qrten(k) + (prr_wau(k) + prr_rcw(k) &
                      + prr_sml(k) + prr_gml(k) + prr_rcs(k) &
                      + prr_rcg(k) - prg_rfz(k) &
                      - pri_rfz(k) - prr_rci(k)) &
                      * orho

!..Rain number tendency
         nrten(k) = nrten(k) + (pnr_wau(k) + pnr_sml(k) + pnr_gml(k)    &
                      - (pnr_rfz(k) + pnr_rcr(k) + pnr_rcg(k)           &
                      + pnr_rcs(k) + pnr_rci(k) + pni_rfz(k)) )         &
                      * orho

!..Rain mass/number balance; keep median volume diameter between
!.. 37 microns (D0r*0.75) and 2.5 mm.
         xrr=MAX(R1,(qr1d(k) + qrten(k)*dtsave)*rho(k))
         xnr=MAX(R2,(nr1d(k) + nrten(k)*dtsave)*rho(k))
         if (xrr.gt. R1) then
           lamr = (am_r*crg(3)*org2*xnr/xrr)**obmr
           mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
           if (mvd_r(k) .gt. 2.5E-3) then
              mvd_r(k) = 2.5E-3
              lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
              xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
              nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
           elseif (mvd_r(k) .lt. D0r*0.75) then
              mvd_r(k) = D0r*0.75
              lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
              xnr = crg(2)*org3*xrr*lamr**bm_r / am_r
              nrten(k) = (xnr-nr1d(k)*rho(k))*odts*orho
           endif
         else
           qrten(k) = -qr1d(k)*odts
           nrten(k) = -nr1d(k)*odts
         endif

!..Snow tendency
         qsten(k) = qsten(k) + (prs_iau(k) + prs_sde(k) &
                      + prs_sci(k) + prs_scw(k) + prs_rcs(k) &
                      + prs_ide(k) - prs_ihm(k) - prr_sml(k)) &
                      * orho

!..Graupel tendency
         qgten(k) = qgten(k) + (prg_scw(k) + prg_rfz(k) &
                      + prg_gde(k) + prg_rcg(k) + prg_gcw(k) &
                      + prg_rci(k) + prg_rcs(k) - prg_ihm(k) &
                      - prr_gml(k)) &
                      * orho

!..Temperature tendency
         if (temp(k).lt.T_0) then
          tten(k) = tten(k) &
                    + ( lsub*ocp(k)*(pri_inu(k) + pri_ide(k) &
                                     + prs_ide(k) + prs_sde(k) &
                                     + prg_gde(k) + pri_iha(k)) &
                     + lfus2*ocp(k)*(pri_wfz(k) + pri_rfz(k) &
                                     + prg_rfz(k) + prs_scw(k) &
                                     + prg_scw(k) + prg_gcw(k) &
                                     + prg_rcs(k) + prs_rcs(k) &
                                     + prr_rci(k) + prg_rcg(k)) &
                       )*orho * (1-IFDRY)
         else
          tten(k) = tten(k) &
                    + ( lfus*ocp(k)*(-prr_sml(k) - prr_gml(k) &
                                     - prr_rcg(k) - prr_rcs(k)) &
                      + lsub*ocp(k)*(prs_sde(k) + prg_gde(k)) &
                       )*orho * (1-IFDRY)
         endif

      enddo

!+---+-----------------------------------------------------------------+
!..Update variables for TAU+1 before condensation & sedimention.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k) + DT*tten(k)
         otemp = 1./temp(k)
         tempc = temp(k) - 273.15
         qv(k) = MAX(1.E-10, qv1d(k) + DT*qvten(k))
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rhof(k) = SQRT(RHO_NOT/rho(k))
         rhof2(k) = SQRT(rhof(k))
         qvs(k) = rslf(pres(k), temp(k))
         ssatw(k) = qv(k)/qvs(k) - 1.
         if (abs(ssatw(k)).lt. eps) ssatw(k) = 0.0
         diffu(k) = 2.11E-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
         if (tempc .ge. 0.0) then
            visco(k) = (1.718+0.0049*tempc)*1.0E-5
         else
            visco(k) = (1.718+0.0049*tempc-1.2E-5*tempc*tempc)*1.0E-5
         endif
         vsc2(k) = SQRT(rho(k)/visco(k))
         lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
         tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
         ocp(k) = 1./(Cp*(1.+0.887*qv(k)))
         lvt2(k)=lvap(k)*lvap(k)*ocp(k)*oRv*otemp*otemp

         nwfa(k) = MAX(11.1E6, (nwfa1d(k) + nwfaten(k)*DT)*rho(k))
      enddo

      do k = kts, kte
         if ((qc1d(k) + qcten(k)*DT) .gt. R1) then
            rc(k) = (qc1d(k) + qcten(k)*DT)*rho(k)
            nc(k) = MAX(2., MIN((nc1d(k)+ncten(k)*DT)*rho(k), Nt_c_max))
            if (.NOT. is_aerosol_aware) nc(k) = Nt_c
            L_qc(k) = .true.
         else
            rc(k) = R1
            nc(k) = 2.
            L_qc(k) = .false.
         endif

         if ((qi1d(k) + qiten(k)*DT) .gt. R1) then
            ri(k) = (qi1d(k) + qiten(k)*DT)*rho(k)
            ni(k) = MAX(R2, (ni1d(k) + niten(k)*DT)*rho(k))
            L_qi(k) = .true. 
         else
            ri(k) = R1
            ni(k) = R2
            L_qi(k) = .false.
         endif

         if ((qr1d(k) + qrten(k)*DT) .gt. R1) then
            rr(k) = (qr1d(k) + qrten(k)*DT)*rho(k)
            nr(k) = MAX(R2, (nr1d(k) + nrten(k)*DT)*rho(k))
            L_qr(k) = .true.
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
            if (mvd_r(k) .gt. 2.5E-3) then
               mvd_r(k) = 2.5E-3
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            elseif (mvd_r(k) .lt. D0r*0.75) then
               mvd_r(k) = D0r*0.75
               lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
               nr(k) = crg(2)*org3*rr(k)*lamr**bm_r / am_r
            endif
         else
            rr(k) = R1
            nr(k) = R2
            L_qr(k) = .false.
         endif
               
         if ((qs1d(k) + qsten(k)*DT) .gt. R1) then
            rs(k) = (qs1d(k) + qsten(k)*DT)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif

         if ((qg1d(k) + qgten(k)*DT) .gt. R1) then
            rg(k) = (qg1d(k) + qgten(k)*DT)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..With tendency-updated mixing ratios, recalculate snow moments and
!.. intercepts/slopes of graupel and rain.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         smo2(k) = 0.
         smob(k) = 0.
         smoc(k) = 0.
         smod(k) = 0.
      enddo
      do k = kts, kte
         if (.not. L_qs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
               + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
               + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
               + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
               + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
               + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
               + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
               + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
               + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif

!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
               + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
               + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
               + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
              + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
              + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_

!..Calculate bm_s+bv_s (th) moment.  Useful for sedimentation.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(14) &
               + sa(4)*tc0*cse(14) + sa(5)*tc0*tc0 &
               + sa(6)*cse(14)*cse(14) + sa(7)*tc0*tc0*cse(14) &
               + sa(8)*tc0*cse(14)*cse(14) + sa(9)*tc0*tc0*tc0 &
               + sa(10)*cse(14)*cse(14)*cse(14)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(14) + sb(4)*tc0*cse(14) &
              + sb(5)*tc0*tc0 + sb(6)*cse(14)*cse(14) &
              + sb(7)*tc0*tc0*cse(14) + sb(8)*tc0*cse(14)*cse(14) &
              + sb(9)*tc0*tc0*tc0 + sb(10)*cse(14)*cse(14)*cse(14)
         smod(k) = a_ * smo2(k)**b_
      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+
      N0_min = gonv_max
      k_0 = kts
      do k = kte, kts, -1
         if (temp(k).ge.270.65) k_0 = MAX(k_0, k)
      enddo
      do k = kte, kts, -1
         if (k.gt.k_0 .and. L_qr(k) .and. mvd_r(k).gt.100.E-6) then
            xslw1 = 4.01 + alog10(mvd_r(k))
         else
            xslw1 = 0.01
         endif
         ygra1 = 4.31 + alog10(max(5.E-5, rg(k)))
         zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         N0_min = MIN(N0_exp, N0_min)
         N0_exp = N0_min
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo

      endif

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for rain.
!+---+-----------------------------------------------------------------+
      do k = kte, kts, -1
         lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
         ilamr(k) = 1./lamr
         mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
         N0_r(k) = nr(k)*org2*lamr**cre(2)
      enddo

!+---+-----------------------------------------------------------------+
!..Cloud water condensation and evaporation.  Nucleate cloud droplets
!.. using explicit CCN aerosols with hygroscopicity like sulfates using
!.. parcel model lookup table results provided by T. Eidhammer.  Evap
!.. drops using calculation of max drop size capable of evaporating in
!.. single timestep and explicit number of drops smaller than Dc_star
!.. from lookup table.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         orho = 1./rho(k)
         if ( (ssatw(k).gt. eps) .or. (ssatw(k).lt. -eps .and. &
                   L_qc(k)) ) then
          clap = (qv(k)-qvs(k))/(1. + lvt2(k)*qvs(k))
          do n = 1, 3
             fcd = qvs(k)* EXP(lvt2(k)*clap) - qv(k) + clap
             dfcd = qvs(k)*lvt2(k)* EXP(lvt2(k)*clap) + 1.
             clap = clap - fcd/dfcd
          enddo
          xrc = rc(k) + clap*rho(k)
          xnc = 0.
          if (xrc.gt. R1) then
           prw_vcd(k) = clap*odt
!+---+-----------------------------------------------------------------+ !  DROPLET NUCLEATION
           if (clap .gt. eps) then
            if (is_aerosol_aware) then
               xnc = MAX(2., activ_ncloud(temp(k), w1d(k), nwfa(k)))
            else
               xnc = Nt_c
            endif
            pnc_wcd(k) = 0.5*(xnc-nc(k) + abs(xnc-nc(k)))*odts*orho

!+---+-----------------------------------------------------------------+ !  EVAPORATION
           elseif (clap .lt. -eps .AND. ssatw(k).lt.-1.E-6 .AND. is_aerosol_aware) then
            tempc = temp(k) - 273.15
            otemp = 1./temp(k)
            rvs = rho(k)*qvs(k)
            rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.)
            rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.) &
                            *otemp*(lvap(k)*otemp*oRv - 1.) &
                            + (-2.*lvap(k)*otemp*otemp*otemp*oRv) &
                            + otemp*otemp)
            gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
            alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                       * rvs_pp/rvs_p * rvs/rvs_p
            alphsc = MAX(1.E-9, alphsc)
            xsat = ssatw(k)
            if (abs(xsat).lt. 1.E-9) xsat=0.
            t1_evap = 2.*PI*( 1.0 - alphsc*xsat  &
                   + 2.*alphsc*alphsc*xsat*xsat  &
                   - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                   / (1.+gamsc)

            Dc_star = DSQRT(-2.D0*DT * t1_evap/(2.*PI) &
                    * 4.*diffu(k)*ssatw(k)*rvs/rho_w)
            idx_d = MAX(1, MIN(INT(1.E6*Dc_star), nbc))

            idx_n = NINT(1.0 + FLOAT(nbc) * DLOG(nc(k)/t_Nc(1)) / nic1)
            idx_n = MAX(1, MIN(idx_n, nbc))

!..Cloud water lookup table index.
            if (rc(k).gt. r_c(1)) then
             nic = NINT(ALOG10(rc(k)))
             do nn = nic-1, nic+1
                n = nn
                if ( (rc(k)/10.**nn).ge.1.0 .and. &
                     (rc(k)/10.**nn).lt.10.0) goto 159
             enddo
 159         continue
             idx_c = INT(rc(k)/10.**n) + 10*(n-nic2) - (n-nic2)
             idx_c = MAX(1, MIN(idx_c, ntb_c))
            else
             idx_c = 1
            endif

           !prw_vcd(k) = MAX(DBLE(-rc(k)*orho*odt),                     &
           !           -tpc_wev(idx_d, idx_c, idx_n)*orho*odt)
            prw_vcd(k) = MAX(DBLE(-rc(k)*0.99*orho*odt), prw_vcd(k))
            pnc_wcd(k) = MAX(DBLE(-nc(k)*0.99*orho*odt),                &
                       -tnc_wev(idx_d, idx_c, idx_n)*orho*odt)

           endif
           if(is_aerosol_aware .and. L_nwfa(k)) L_nwfa(k) = .false.
          else
           prw_vcd(k) = -rc(k)*orho*odt
           pnc_wcd(k) = -nc(k)*orho*odt
          endif

!+---+-----------------------------------------------------------------+

          qvten(k) = qvten(k) - prw_vcd(k)
          qcten(k) = qcten(k) + prw_vcd(k)
          ncten(k) = ncten(k) + pnc_wcd(k)
          nwfaten(k) = nwfaten(k) - pnc_wcd(k)
          tten(k) = tten(k) + lvap(k)*ocp(k)*prw_vcd(k)*(1-IFDRY)
          rc(k) = MAX(R1, (qc1d(k) + DT*qcten(k))*rho(k))
          if (rc(k).eq.R1) L_qc(k) = .false.
          nc(k) = MAX(2., MIN((nc1d(k)+ncten(k)*DT)*rho(k), Nt_c_max))
          if (.NOT. is_aerosol_aware) nc(k) = Nt_c
          qv(k) = MAX(1.E-10, qv1d(k) + DT*qvten(k))
          temp(k) = t1d(k) + DT*tten(k)
          rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
          qvs(k) = rslf(pres(k), temp(k))
          ssatw(k) = qv(k)/qvs(k) - 1.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!.. If still subsaturated, allow rain to evaporate, following
!.. Srivastava & Coen (1992).
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         if ( (ssatw(k).lt. -eps) .and. L_qr(k) &
                     .and. (.not.(prw_vcd(k).gt. 0.)) ) then
          tempc = temp(k) - 273.15
          otemp = 1./temp(k)
          orho = 1./rho(k)
          rhof(k) = SQRT(RHO_NOT*orho)
          rhof2(k) = SQRT(rhof(k))
          diffu(k) = 2.11E-5*(temp(k)/273.15)**1.94 * (101325./pres(k))
          if (tempc .ge. 0.0) then
             visco(k) = (1.718+0.0049*tempc)*1.0E-5
          else
             visco(k) = (1.718+0.0049*tempc-1.2E-5*tempc*tempc)*1.0E-5
          endif
          vsc2(k) = SQRT(rho(k)/visco(k))
          lvap(k) = lvap0 + (2106.0 - 4218.0)*tempc
          tcond(k) = (5.69 + 0.0168*tempc)*1.0E-5 * 418.936
          ocp(k) = 1./(Cp*(1.+0.887*qv(k)))

          rvs = rho(k)*qvs(k)
          rvs_p = rvs*otemp*(lvap(k)*otemp*oRv - 1.)
          rvs_pp = rvs * ( otemp*(lvap(k)*otemp*oRv - 1.) &
                          *otemp*(lvap(k)*otemp*oRv - 1.) &
                          + (-2.*lvap(k)*otemp*otemp*otemp*oRv) &
                          + otemp*otemp)
          gamsc = lvap(k)*diffu(k)/tcond(k) * rvs_p
          alphsc = 0.5*(gamsc/(1.+gamsc))*(gamsc/(1.+gamsc)) &
                     * rvs_pp/rvs_p * rvs/rvs_p
          alphsc = MAX(1.E-9, alphsc)
          xsat   = MIN(-1.E-9, ssatw(k))
          t1_evap = 2.*PI*( 1.0 - alphsc*xsat  &
                 + 2.*alphsc*alphsc*xsat*xsat  &
                 - 5.*alphsc*alphsc*alphsc*xsat*xsat*xsat ) &
                 / (1.+gamsc)

          lamr = 1./ilamr(k)
!..Rapidly eliminate near zero values when low humidity (<95%)
          if (qv(k)/qvs(k) .lt. 0.95 .AND. rr(k)*orho.le.1.E-8) then
          prv_rev(k) = rr(k)*orho*odts
          else
          prv_rev(k) = t1_evap*diffu(k)*(-ssatw(k))*N0_r(k)*rvs &
              * (t1_qr_ev*ilamr(k)**cre(10) &
              + t2_qr_ev*vsc2(k)*rhof2(k)*((lamr+0.5*fv_r)**(-cre(11))))
          rate_max = MIN((rr(k)*orho*odts), (qvs(k)-qv(k))*odts)
          prv_rev(k) = MIN(DBLE(rate_max), prv_rev(k)*orho)

!..TEST: G. Thompson  10 May 2013
!..Reduce the rain evaporation in same places as melting graupel occurs.
!..Rationale: falling and simultaneous melting graupel in subsaturated
!..regions will not melt as fast because particle temperature stays
!..at 0C.  Also not much shedding of the water from the graupel so
!..likely that the water-coated graupel evaporating much slower than
!..if the water was immediately shed off.
          IF (prr_gml(k).gt.0.0) THEN
             eva_factor = MIN(1.0, 0.01+(0.99-0.01)*(tempc/20.0))
             prv_rev(k) = prv_rev(k)*eva_factor
          ENDIF
          endif

          pnr_rev(k) = MIN(DBLE(nr(k)*0.99*orho*odts),                  &   ! RAIN2M
                       prv_rev(k) * nr(k)/rr(k))

          qrten(k) = qrten(k) - prv_rev(k)
          qvten(k) = qvten(k) + prv_rev(k)
          nrten(k) = nrten(k) - pnr_rev(k)
          nwfaten(k) = nwfaten(k) + pnr_rev(k)
          tten(k) = tten(k) - lvap(k)*ocp(k)*prv_rev(k)*(1-IFDRY)

          rr(k) = MAX(R1, (qr1d(k) + DT*qrten(k))*rho(k))
          qv(k) = MAX(1.E-10, qv1d(k) + DT*qvten(k))
          nr(k) = MAX(R2, (nr1d(k) + DT*nrten(k))*rho(k))
          temp(k) = t1d(k) + DT*tten(k)
          rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         endif
      enddo
      do k = kts, kte
         evapprod(k) = prv_rev(k) - (min(zeroD0,prs_sde(k)) + &
                                     min(zeroD0,prg_gde(k)))
         rainprod(k) = prr_wau(k) + prr_rcw(k) + prs_scw(k) + &
                                    prg_scw(k) + prs_iau(k) + &
                                    prg_gcw(k) + prs_sci(k) + &
                                    pri_rci(k)
      enddo

!+---+-----------------------------------------------------------------+
!..Find max terminal fallspeed (distribution mass-weighted mean
!.. velocity) and use it to determine if we need to split the timestep
!.. (var nstep>1).  Either way, only bother to do sedimentation below
!.. 1st level that contains any sedimenting particles (k=ksed1 on down).
!.. New in v3.0+ is computing separate for rain, ice, snow, and
!.. graupel species thus making code faster with credit to J. Schmidt.
!+---+-----------------------------------------------------------------+
      nstep = 0
      onstep(:) = 1.0
      ksed1(:) = 1
      do k = kte+1, kts, -1
         vtrk(k) = 0.
         vtnrk(k) = 0.
         vtik(k) = 0.
         vtnik(k) = 0.
         vtsk(k) = 0.
         vtgk(k) = 0.
         vtck(k) = 0.
         vtnck(k) = 0.
      enddo

      if (ANY(L_qr .eqv. .true.)) then
      do k = kte, kts, -1
         vtr = 0.
         rhof(k) = SQRT(RHO_NOT/rho(k))

         if (rr(k).gt. R1) then
          lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
          vtr = rhof(k)*av_r*crg(6)*org3 * lamr**cre(3)                 &
                      *((lamr+fv_r)**(-cre(6)))
          vtrk(k) = vtr
! First below is technically correct:
!         vtr = rhof(k)*av_r*crg(5)*org2 * lamr**cre(2)                 &
!                     *((lamr+fv_r)**(-cre(5)))
! Test: make number fall faster (but still slower than mass)
! Goal: less prominent size sorting
          vtr = rhof(k)*av_r*crg(7)/crg(12) * lamr**cre(12)             &
                      *((lamr+fv_r)**(-cre(7)))
          vtnrk(k) = vtr
         else
          vtrk(k) = vtrk(k+1)
          vtnrk(k) = vtnrk(k+1)
         endif

         if (MAX(vtrk(k),vtnrk(k)) .gt. 1.E-3) then
            ksed1(1) = MAX(ksed1(1), k)
            delta_tp = dzq(k)/(MAX(vtrk(k),vtnrk(k)))
            nstep = MAX(nstep, INT(DT/delta_tp + 1.))
         endif
      enddo
      if (ksed1(1) .eq. kte) ksed1(1) = kte-1
      if (nstep .gt. 0) onstep(1) = 1./REAL(nstep)
      endif

!+---+-----------------------------------------------------------------+

      if (ANY(L_qc .eqv. .true.)) then
      hgt_agl = 0.
      do k = kts, kte-1
         if (rc(k) .gt. R2) ksed1(5) = k
         hgt_agl = hgt_agl + dzq(k)
         if (hgt_agl .gt. 500.0) goto 151
      enddo
 151  continue

      do k = ksed1(5), kts, -1
         vtc = 0.
         if (rc(k) .gt. R1 .and. w1d(k) .lt. 1.E-1) then
          nu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
          lamc = (nc(k)*am_r*ccg(2,nu_c)*ocg1(nu_c)/rc(k))**obmr
          ilamc = 1./lamc
          vtc = rhof(k)*av_c*ccg(5,nu_c)*ocg2(nu_c) * ilamc**bv_c
          vtck(k) = vtc
          vtc = rhof(k)*av_c*ccg(4,nu_c)*ocg1(nu_c) * ilamc**bv_c
          vtnck(k) = vtc
         endif
      enddo
      endif

!+---+-----------------------------------------------------------------+

      if (.not. iiwarm) then

       if (ANY(L_qi .eqv. .true.)) then
       nstep = 0
       do k = kte, kts, -1
          vti = 0.

          if (ri(k).gt. R1) then
           lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
           ilami = 1./lami
           vti = rhof(k)*av_i*cig(3)*oig2 * ilami**bv_i
           vtik(k) = vti
! First below is technically correct:
!          vti = rhof(k)*av_i*cig(4)*oig1 * ilami**bv_i
! Goal: less prominent size sorting
           vti = rhof(k)*av_i*cig(6)/cig(7) * ilami**bv_i
           vtnik(k) = vti
          else
           vtik(k) = vtik(k+1)
           vtnik(k) = vtnik(k+1)
          endif

          if (vtik(k) .gt. 1.E-3) then
             ksed1(2) = MAX(ksed1(2), k)
             delta_tp = dzq(k)/vtik(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.))
          endif
       enddo
       if (ksed1(2) .eq. kte) ksed1(2) = kte-1
       if (nstep .gt. 0) onstep(2) = 1./REAL(nstep)
       endif

!+---+-----------------------------------------------------------------+

       if (ANY(L_qs .eqv. .true.)) then
       nstep = 0
       do k = kte, kts, -1
          vts = 0.

          if (rs(k).gt. R1) then
           xDs = smoc(k) / smob(k)
           Mrat = 1./xDs
           ils1 = 1./(Mrat*Lam0 + fv_s)
           ils2 = 1./(Mrat*Lam1 + fv_s)
           t1_vts = Kap0*csg(4)*ils1**cse(4)
           t2_vts = Kap1*Mrat**mu_s*csg(10)*ils2**cse(10)
           ils1 = 1./(Mrat*Lam0)
           ils2 = 1./(Mrat*Lam1)
           t3_vts = Kap0*csg(1)*ils1**cse(1)
           t4_vts = Kap1*Mrat**mu_s*csg(7)*ils2**cse(7)
           vts = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
           if (temp(k).gt. (T_0+0.1)) then
            SR = rs(k)/(rs(k)+rr(k))
            vtsk(k) = vts*SR + (1.-SR)*vtrk(k)
           else
            vtsk(k) = vts*vts_boost(k)
           endif
          else
            vtsk(k) = vtsk(k+1)
          endif

          if (vtsk(k) .gt. 1.E-3) then
             ksed1(3) = MAX(ksed1(3), k)
             delta_tp = dzq(k)/vtsk(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.))
          endif
       enddo
       if (ksed1(3) .eq. kte) ksed1(3) = kte-1
       if (nstep .gt. 0) onstep(3) = 1./REAL(nstep)
       endif

!+---+-----------------------------------------------------------------+

       if (ANY(L_qg .eqv. .true.)) then
       nstep = 0
       do k = kte, kts, -1
          vtg = 0.

          if (rg(k).gt. R1) then
           vtg = rhof(k)*av_g*cgg(6)*ogg3 * ilamg(k)**bv_g
           if (temp(k).gt. T_0) then
            vtgk(k) = MAX(vtg, vtrk(k))
           else
            vtgk(k) = vtg
           endif
          else
            vtgk(k) = vtgk(k+1)
          endif

          if (vtgk(k) .gt. 1.E-3) then
             ksed1(4) = MAX(ksed1(4), k)
             delta_tp = dzq(k)/vtgk(k)
             nstep = MAX(nstep, INT(DT/delta_tp + 1.))
          endif
       enddo
       if (ksed1(4) .eq. kte) ksed1(4) = kte-1
       if (nstep .gt. 0) onstep(4) = 1./REAL(nstep)
       endif
      endif

!+---+-----------------------------------------------------------------+
!..Sedimentation of mixing ratio is the integral of v(D)*m(D)*N(D)*dD,
!.. whereas neglect m(D) term for number concentration.  Therefore,
!.. cloud ice has proper differential sedimentation.
!+---+-----------------------------------------------------------------+

      if (ANY(L_qr .eqv. .true.)) then
      !--- initializes sedimentation flux of rain:
      do k = kts,kte
         ppfr(k) = qrten(k)
      enddo
      nstep = NINT(1./onstep(1))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_r(k) = vtrk(k)*rr(k)
            sed_n(k) = vtnrk(k)*nr(k)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qrten(k) = qrten(k) - sed_r(k)*odzq*onstep(1)*orho
         nrten(k) = nrten(k) - sed_n(k)*odzq*onstep(1)*orho
         rr(k) = MAX(R1, rr(k) - sed_r(k)*odzq*DT*onstep(1))
         nr(k) = MAX(R2, nr(k) - sed_n(k)*odzq*DT*onstep(1))
         do k = ksed1(1), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qrten(k) = qrten(k) + (sed_r(k+1)-sed_r(k))                 &
                                               *odzq*onstep(1)*orho
            nrten(k) = nrten(k) + (sed_n(k+1)-sed_n(k))                 &
                                               *odzq*onstep(1)*orho
            rr(k) = MAX(R1, rr(k) + (sed_r(k+1)-sed_r(k)) &
                                           *odzq*DT*onstep(1))
            nr(k) = MAX(R2, nr(k) + (sed_n(k+1)-sed_n(k)) &
                                           *odzq*DT*onstep(1))
         enddo

         if (rr(kts).gt.R1*1000.) &
         pptrain = pptrain + sed_r(kts)*DT*onstep(1)
      enddo
      !--- computes sedimentation flux of rain:
      do k = kts,kte
         ppfr(k) = (qrten(k)-ppfr(k))/(odzq*orho)
      enddo
      endif

!+---+-----------------------------------------------------------------+

      if (ANY(L_qc .eqv. .true.)) then
      do k = kte, kts, -1
         sed_c(k) = vtck(k)*rc(k)
         sed_n(k) = vtnck(k)*nc(k)
      enddo
      do k = ksed1(5), kts, -1
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qcten(k) = qcten(k) + (sed_c(k+1)-sed_c(k)) *odzq*orho
         ncten(k) = ncten(k) + (sed_n(k+1)-sed_n(k)) *odzq*orho
         rc(k) = MAX(R1, rc(k) + (sed_c(k+1)-sed_c(k)) *odzq*DT)
         nc(k) = MAX(10., nc(k) + (sed_n(k+1)-sed_n(k)) *odzq*DT)
      enddo
      endif

!+---+-----------------------------------------------------------------+

      if (ANY(L_qi .eqv. .true.)) then
      !--- initializes sedimentation flux of ice:
      do k = kts,kte
         ppfi(k) = qiten(k)
      enddo
      nstep = NINT(1./onstep(2))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_i(k) = vtik(k)*ri(k)
            sed_n(k) = vtnik(k)*ni(k)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qiten(k) = qiten(k) - sed_i(k)*odzq*onstep(2)*orho
         niten(k) = niten(k) - sed_n(k)*odzq*onstep(2)*orho
         ri(k) = MAX(R1, ri(k) - sed_i(k)*odzq*DT*onstep(2))
         ni(k) = MAX(R2, ni(k) - sed_n(k)*odzq*DT*onstep(2))
         do k = ksed1(2), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qiten(k) = qiten(k) + (sed_i(k+1)-sed_i(k))                 &
                                               *odzq*onstep(2)*orho
            niten(k) = niten(k) + (sed_n(k+1)-sed_n(k))                 &
                                               *odzq*onstep(2)*orho
            ri(k) = MAX(R1, ri(k) + (sed_i(k+1)-sed_i(k)) &
                                           *odzq*DT*onstep(2))
            ni(k) = MAX(R2, ni(k) + (sed_n(k+1)-sed_n(k)) &
                                           *odzq*DT*onstep(2))
         enddo

         if (ri(kts).gt.R1*1000.) &
         pptice = pptice + sed_i(kts)*DT*onstep(2)
      enddo
      !--- computes sedimentation flux of ice:
      do k = kts,kte
         ppfi(k) = (qiten(k)-ppfi(k))/(odzq*orho)
      enddo
      endif

!+---+-----------------------------------------------------------------+

      if (ANY(L_qs .eqv. .true.)) then
      !--- initializes sedimentation flux of snow:
      do k = kts,kte
         ppfs(k) = qsten(k)
      enddo
      nstep = NINT(1./onstep(3))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_s(k) = vtsk(k)*rs(k)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qsten(k) = qsten(k) - sed_s(k)*odzq*onstep(3)*orho
         rs(k) = MAX(R1, rs(k) - sed_s(k)*odzq*DT*onstep(3))
         do k = ksed1(3), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qsten(k) = qsten(k) + (sed_s(k+1)-sed_s(k))                 &
                                               *odzq*onstep(3)*orho
            rs(k) = MAX(R1, rs(k) + (sed_s(k+1)-sed_s(k)) &
                                           *odzq*DT*onstep(3))
         enddo

         if (rs(kts).gt.R1*1000.) &
         pptsnow = pptsnow + sed_s(kts)*DT*onstep(3)
      enddo
      !--- computes sedimentation flux of snow:
      do k = kts,kte
         ppfs(k) = (qsten(k)-ppfs(k))/(odzq*orho)
      enddo
      endif

!+---+-----------------------------------------------------------------+

      if (ANY(L_qg .eqv. .true.)) then
      !--- initializes sedimentation flux of graupel:
      do k = kts,kte
         ppfg(k) = qgten(k)
      enddo
      nstep = NINT(1./onstep(4))
      do n = 1, nstep
         do k = kte, kts, -1
            sed_g(k) = vtgk(k)*rg(k)
         enddo
         k = kte
         odzq = 1./dzq(k)
         orho = 1./rho(k)
         qgten(k) = qgten(k) - sed_g(k)*odzq*onstep(4)*orho
         rg(k) = MAX(R1, rg(k) - sed_g(k)*odzq*DT*onstep(4))
         do k = ksed1(4), kts, -1
            odzq = 1./dzq(k)
            orho = 1./rho(k)
            qgten(k) = qgten(k) + (sed_g(k+1)-sed_g(k))                 &
                                               *odzq*onstep(4)*orho
            rg(k) = MAX(R1, rg(k) + (sed_g(k+1)-sed_g(k)) &
                                           *odzq*DT*onstep(4))
         enddo

         if (rg(kts).gt.R1*1000.) &
         pptgraul = pptgraul + sed_g(kts)*DT*onstep(4)
      enddo
      !--- computes sedimentation flux of graupel:
      do k = kts,kte
         ppfg(k) = (qgten(k)-ppfg(k))/(odzq*orho)
      enddo
      endif

!+---+-----------------------------------------------------------------+
!.. Instantly melt any cloud ice into cloud water if above 0C and
!.. instantly freeze any cloud water found below HGFR.
!+---+-----------------------------------------------------------------+
      if (.not. iiwarm) then
      do k = kts, kte
         xri = MAX(0.0, qi1d(k) + qiten(k)*DT)
         if ( (temp(k).gt. T_0) .and. (xri.gt. 0.0) ) then
          qcten(k) = qcten(k) + xri*odt
          ncten(k) = ncten(k) + ni1d(k)*odt
          qiten(k) = qiten(k) - xri*odt
          niten(k) = -ni1d(k)*odt
          tten(k) = tten(k) - lfus*ocp(k)*xri*odt*(1-IFDRY)
         endif

         xrc = MAX(0.0, qc1d(k) + qcten(k)*DT)
         if ( (temp(k).lt. HGFR) .and. (xrc.gt. 0.0) ) then
          lfus2 = lsub - lvap(k)
          xnc = nc1d(k) + ncten(k)*DT
          qiten(k) = qiten(k) + xrc*odt
          niten(k) = niten(k) + xnc*odt
          qcten(k) = qcten(k) - xrc*odt
          ncten(k) = ncten(k) - xnc*odt
          tten(k) = tten(k) + lfus2*ocp(k)*xrc*odt*(1-IFDRY)
         endif
      enddo
      endif

!+---+-----------------------------------------------------------------+
!.. All tendencies computed, apply and pass back final values to parent.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         t1d(k)  = t1d(k) + tten(k)*DT
         qv1d(k) = MAX(1.E-10, qv1d(k) + qvten(k)*DT)
         qc1d(k) = qc1d(k) + qcten(k)*DT
         nc1d(k) = MAX(2./rho(k), MIN(nc1d(k) + ncten(k)*DT, Nt_c_max))
         nwfa1d(k) = MAX(11.1E6, MIN(9999.E6,                           &
                       (nwfa1d(k)+nwfaten(k)*DT)))
         nifa1d(k) = MAX(naIN1*0.01, MIN(9999.E6,                       &
                       (nifa1d(k)+nifaten(k)*DT)))

         if (qc1d(k) .le. R1) then
           qc1d(k) = 0.0
           nc1d(k) = 0.0
         else
           nu_c = MIN(15, NINT(1000.E6/(nc1d(k)*rho(k))) + 2)
           lamc = (am_r*ccg(2,nu_c)*ocg1(nu_c)*nc1d(k)/qc1d(k))**obmr
           xDc = (bm_r + nu_c + 1.) / lamc
           if (xDc.lt. D0c) then
            lamc = cce(2,nu_c)/D0c
           elseif (xDc.gt. D0r*2.) then
            lamc = cce(2,nu_c)/(D0r*2.)
           endif
           nc1d(k) = MIN(ccg(1,nu_c)*ocg2(nu_c)*qc1d(k)/am_r*lamc**bm_r,&
                         DBLE(Nt_c_max)/rho(k))
         endif

         qi1d(k) = qi1d(k) + qiten(k)*DT
         ni1d(k) = MAX(R2/rho(k), ni1d(k) + niten(k)*DT)
         if (qi1d(k) .le. R1) then
           qi1d(k) = 0.0
           ni1d(k) = 0.0
         else
           lami = (am_i*cig(2)*oig1*ni1d(k)/qi1d(k))**obmi
           ilami = 1./lami
           xDi = (bm_i + mu_i + 1.) * ilami
           if (xDi.lt. 5.E-6) then
            lami = cie(2)/5.E-6
           elseif (xDi.gt. 300.E-6) then
            lami = cie(2)/300.E-6
           endif
           ni1d(k) = MIN(cig(1)*oig2*qi1d(k)/am_i*lami**bm_i,           &
                         9999.D3/rho(k))
         endif
         qr1d(k) = qr1d(k) + qrten(k)*DT
         nr1d(k) = MAX(R2/rho(k), nr1d(k) + nrten(k)*DT)
         if (qr1d(k) .le. R1) then
           qr1d(k) = 0.0
           nr1d(k) = 0.0
         else
           lamr = (am_r*crg(3)*org2*nr1d(k)/qr1d(k))**obmr
           mvd_r(k) = (3.0 + mu_r + 0.672) / lamr
           if (mvd_r(k) .gt. 2.5E-3) then
              mvd_r(k) = 2.5E-3
           elseif (mvd_r(k) .lt. D0r*0.75) then
              mvd_r(k) = D0r*0.75
           endif
           lamr = (3.0 + mu_r + 0.672) / mvd_r(k)
           nr1d(k) = crg(2)*org3*qr1d(k)*lamr**bm_r / am_r
         endif
         qs1d(k) = qs1d(k) + qsten(k)*DT
         if (qs1d(k) .le. R1) qs1d(k) = 0.0
         qg1d(k) = qg1d(k) + qgten(k)*DT
         if (qg1d(k) .le. R1) qg1d(k) = 0.0
      enddo

      end subroutine thompson_1d

!==================================================================================================================
!..Creation of the lookup tables and support functions found below here.
!..Rain collecting graupel (and inverse). Explicit CE integration.

 subroutine qr_acr_qg
 implicit none
!==================================================================================================================

!... local variables:
 integer:: i,j,k,m,n,n2
 integer:: km,km_s,km_e
 real(kind=R8KIND),dimension(nbg):: vg, N_g
 real(kind=R8KIND),dimension(nbr):: vr, N_r
 real(kind=R8KIND):: N0_r,N0_g,lam_exp,lamg,lamr
 real(kind=R8KIND):: massg,massr,dvg,dvr,t1,t2,z1,z2,y1,y2

!------------------------------------------------------------------------------------------------------------------

 do n2 = 1,nbr
!   vr(n2) = av_r*Dr(n2)**bv_r * dexp(-fv_r*Dr(n2))
    vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2) &
           + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                    &
           - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
 enddo
 do n = 1,nbg
    vg(n) = av_g*Dg(n)**bv_g
 enddo

 km_s = 0
 km_e = ntb_r*ntb_r1 - 1

 do km = km_s,km_e
    m = km / ntb_r1 + 1
    k = mod( km , ntb_r1 ) + 1

    lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
    lamr = lam_exp * (crg(3)*org2*org1)**obmr
    N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
    do n2 = 1,nbr
       N_r(n2) = N0_r*Dr(n2)**mu_r *dexp(-lamr*Dr(n2))*dtr(n2)
    enddo

    do j = 1,ntb_g
       do i = 1,ntb_g1
          lam_exp = (N0g_exp(i)*am_g*cgg(1)/r_g(j))**oge1
          lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
          N0_g = N0g_exp(i)/(cgg(2)*lam_exp) * lamg**cge(2)
          do n = 1,nbg
             N_g(n) = N0_g*Dg(n)**mu_g * dexp(-lamg*Dg(n))*dtg(n)
          enddo

          t1 = 0._R8KIND
          t2 = 0._R8KIND
          z1 = 0._R8KIND
          z2 = 0._R8KIND
          y1 = 0._R8KIND
          y2 = 0._R8KIND
          do n2 = 1,nbr
             massr = am_r * Dr(n2)**bm_r
             do n = 1,nbg
                massg = am_g * Dg(n)**bm_g

                dvg = 0.5d0*((vr(n2) - vg(n)) + dabs(vr(n2)-vg(n)))
                dvr = 0.5d0*((vg(n) - vr(n2)) + dabs(vg(n)-vr(n2)))

                t1 = t1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                   * dvg*massg * N_g(n)* N_r(n2)
                z1 = z1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                   * dvg*massr * N_g(n)* N_r(n2)
                y1 = y1+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                   * dvg * N_g(n)* N_r(n2)

                t2 = t2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                   * dvr*massr * N_g(n)* N_r(n2)
                y2 = y2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                   * dvr * N_g(n)* N_r(n2)
                z2 = z2+ PI*.25*Ef_rg*(Dg(n)+Dr(n2))*(Dg(n)+Dr(n2)) &
                   * dvr*massg * N_g(n)* N_r(n2)
             enddo
 97          continue
          enddo
          tcg_racg(i,j,k,m) = t1
          tmr_racg(i,j,k,m) = dmin1(z1,r_r(m)*1.0d0)
          tcr_gacr(i,j,k,m) = t2
          tmg_gacr(i,j,k,m) = dmin1(z2,r_g(j)*1.0d0)
          tnr_racg(i,j,k,m) = y1
          tnr_gacr(i,j,k,m) = y2
       enddo
    enddo
 enddo

 end subroutine qr_acr_qg

!==================================================================================================================
!..Creation of the lookup tables and support functions found below here.
!..Rain collecting snow (and inverse). Explicit CE integration.

 subroutine qr_acr_qs
 implicit none
!==================================================================================================================

!... local variables:
 integer:: i,j,k,m,n,n2
 integer:: km,km_s,km_e
 real(kind=R8KIND),dimension(nbr):: vr,D1,N_r
 real(kind=R8KIND),dimension(nbs):: vs,N_s
 real(kind=R8KIND):: loga_,a_,b_,second,M0,M2,M3,Mrat,oM3
 real(kind=R8KIND):: N0_r,lam_exp,lamr,slam1,slam2
 real(kind=R8KIND):: dvs,dvr,masss,massr
 real(kind=R8KIND):: t1,t2,t3,t4,z1,z2,z3,z4
 real(kind=R8KIND):: y1,y2,y3,y4

!------------------------------------------------------------------------------------------------------------------

 do n2 = 1,nbr
!   vr(n2) = av_r*Dr(n2)**bv_r * dexp(-fv_r*Dr(n2))
    vr(n2) = -0.1021 + 4.932E3*Dr(n2) - 0.9551E6*Dr(n2)*Dr(n2) &
           + 0.07934E9*Dr(n2)*Dr(n2)*Dr(n2)                    &
           - 0.002362E12*Dr(n2)*Dr(n2)*Dr(n2)*Dr(n2)
    D1(n2) = (vr(n2)/av_s)**(1./bv_s)
 enddo
 do n = 1,nbs
    vs(n) = 1.5*av_s*Ds(n)**bv_s * dexp(-fv_s*Ds(n))
 enddo

 km_s = 0
 km_e = ntb_r*ntb_r1 - 1

 do km = km_s,km_e
    m = km / ntb_r1 + 1
    k = mod( km , ntb_r1 ) + 1

    lam_exp = (N0r_exp(k)*am_r*crg(1)/r_r(m))**ore1
    lamr = lam_exp * (crg(3)*org2*org1)**obmr
    N0_r = N0r_exp(k)/(crg(2)*lam_exp) * lamr**cre(2)
    do n2 = 1, nbr
       N_r(n2) = N0_r*Dr(n2)**mu_r * dexp(-lamr*Dr(n2))*dtr(n2)
    enddo

    do j = 1,ntb_t
       do i = 1,ntb_s

!..From the bm_s moment, compute plus one moment.  If we are not
!..using bm_s=2, then we must transform to the pure 2nd moment
!..(variable called "second") and then to the bm_s+1 moment.

          M2 = r_s(i)*oams *1.0d0
          if(bm_s.gt.2.0-1.E-3 .and. bm_s.lt.2.0+1.E-3) then
             loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*bm_s &
                   + sa(4)*Tc(j)*bm_s + sa(5)*Tc(j)*Tc(j) &
                   + sa(6)*bm_s*bm_s + sa(7)*Tc(j)*Tc(j)*bm_s &
                   + sa(8)*Tc(j)*bm_s*bm_s + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                   + sa(10)*bm_s*bm_s*bm_s
             a_ = 10.0**loga_
             b_ = sb(1) + sb(2)*Tc(j) + sb(3)*bm_s &
                + sb(4)*Tc(j)*bm_s + sb(5)*Tc(j)*Tc(j) &
                + sb(6)*bm_s*bm_s + sb(7)*Tc(j)*Tc(j)*bm_s &
                + sb(8)*Tc(j)*bm_s*bm_s + sb(9)*Tc(j)*Tc(j)*Tc(j) &
                + sb(10)*bm_s*bm_s*bm_s
             second = (M2/a_)**(1./b_)
          else
             second = M2
          endif

          loga_ = sa(1) + sa(2)*Tc(j) + sa(3)*cse(1) &
                + sa(4)*Tc(j)*cse(1) + sa(5)*Tc(j)*Tc(j) &
                + sa(6)*cse(1)*cse(1) + sa(7)*Tc(j)*Tc(j)*cse(1) &
                + sa(8)*Tc(j)*cse(1)*cse(1) + sa(9)*Tc(j)*Tc(j)*Tc(j) &
                + sa(10)*cse(1)*cse(1)*cse(1)
          a_ = 10.0**loga_
          b_ = sb(1)+sb(2)*Tc(j)+sb(3)*cse(1) + sb(4)*Tc(j)*cse(1) &
                + sb(5)*Tc(j)*Tc(j) + sb(6)*cse(1)*cse(1) &
                + sb(7)*Tc(j)*Tc(j)*cse(1) + sb(8)*Tc(j)*cse(1)*cse(1) &
                + sb(9)*Tc(j)*Tc(j)*Tc(j)+sb(10)*cse(1)*cse(1)*cse(1)
          M3 = a_ * second**b_

          oM3 = 1./M3
          Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
          M0   = (M2*oM3)**mu_s
          slam1 = M2 * oM3 * Lam0
          slam2 = M2 * oM3 * Lam1

          do n = 1,nbs
             N_s(n) = Mrat*(Kap0*dexp(-slam1*Ds(n)) &
                    + Kap1*M0*Ds(n)**mu_s * DEXP(-slam2*Ds(n)))*dts(n)
          enddo

          t1 = 0._R8KIND
          t2 = 0._R8KIND
          t3 = 0._R8KIND
          t4 = 0._R8KIND
          z1 = 0._R8KIND
          z2 = 0._R8KIND
          z3 = 0._R8KIND
          z4 = 0._R8KIND
          y1 = 0._R8KIND
          y2 = 0._R8KIND
          y3 = 0._R8KIND
          y4 = 0._R8KIND
          do n2 = 1,nbr
             massr = am_r * Dr(n2)**bm_r
             do n = 1,nbs
                masss = am_s * Ds(n)**bm_s
      
                dvs = 0.5d0*((vr(n2) - vs(n)) + dabs(vr(n2)-vs(n)))
                dvr = 0.5d0*((vs(n) - vr(n2)) + dabs(vs(n)-vr(n2)))

                if(massr .gt. 1.5*masss) then
                   t1 = t1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvs*masss * N_s(n)* N_r(n2)
                   z1 = z1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvs*massr * N_s(n)* N_r(n2)
                   y1 = y1+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvs * N_s(n)* N_r(n2)
                else
                   t3 = t3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvs*masss * N_s(n)* N_r(n2)
                   z3 = z3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvs*massr * N_s(n)* N_r(n2)
                   y3 = y3+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvs * N_s(n)* N_r(n2)
                endif

                if(massr .gt. 1.5*masss) then
                   t2 = t2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvr*massr * N_s(n)* N_r(n2)
                   y2 = y2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvr * N_s(n)* N_r(n2)
                   z2 = z2+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvr*masss * N_s(n)* N_r(n2)
                else
                   t4 = t4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvr*massr * N_s(n)* N_r(n2)
                   y4 = y4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvr * N_s(n)* N_r(n2)
                   z4 = z4+ PI*.25*Ef_rs*(Ds(n)+Dr(n2))*(Ds(n)+Dr(n2)) &
                      * dvr*masss * N_s(n)* N_r(n2)
                endif

             enddo
          enddo
          tcs_racs1(i,j,k,m) = t1
          tmr_racs1(i,j,k,m) = dmin1(z1,r_r(m)*1.0d0)
          tcs_racs2(i,j,k,m) = t3
          tmr_racs2(i,j,k,m) = z3
          tcr_sacr1(i,j,k,m) = t2
          tms_sacr1(i,j,k,m) = z2
          tcr_sacr2(i,j,k,m) = t4
          tms_sacr2(i,j,k,m) = z4
          tnr_racs1(i,j,k,m) = y1
          tnr_racs2(i,j,k,m) = y3
          tnr_sacr1(i,j,k,m) = y2
          tnr_sacr2(i,j,k,m) = y4
       enddo
    enddo
 enddo

 end subroutine qr_acr_qs

!==================================================================================================================
!..This is a literal adaptation of Bigg (1954) probability of drops of a particular volume freezing.
!..Given this probability, simply freeze the proportion of drops summing their masses.

 subroutine freezeH2O
 implicit none
!==================================================================================================================

!... local variables
 integer:: i,j,k,m,n,n2
 integer:: km,km_s,km_e
 integer:: nu_c
 real(kind=R8KIND):: N_r, N_c
 real(kind=R8KIND),dimension(nbr):: massr
 real(kind=R8KIND),dimension(nbc):: massc
 real(kind=R8KIND):: sum1,sum2,sumn1,sumn2,prob,vol,Texp,orho_w,lam_exp,lamr,N0_r,lamc,N0_c,y
 real(kind=RKIND):: T_adjust

!------------------------------------------------------------------------------------------------------------------

 orho_w = 1./rho_w

 do n2 = 1,nbr
    massr(n2) = am_r*Dr(n2)**bm_r
 enddo
 do n = 1,nbc
    massc(n) = am_r*Dc(n)**bm_r
 enddo

 km_s = 0
 km_e = ntb_IN*45 - 1

!..Freeze water (smallest drops become cloud ice, otherwise graupel).
 do km = km_s,km_e
    m = km / 45 + 1
    k = mod(km ,45) + 1
    T_adjust = max(-3.0,min(3.0 - alog10(Nt_IN(m)),3.0))
    Texp = dexp(real(k,kind=R8SIZE) - T_adjust*1.0D0 ) - 1._R8KIND
    do j = 1,ntb_r1
       do i = 1,ntb_r
          lam_exp = (N0r_exp(j)*am_r*crg(1)/r_r(i))**ore1
          lamr = lam_exp * (crg(3)*org2*org1)**obmr
          N0_r = N0r_exp(j)/(crg(2)*lam_exp) * lamr**cre(2)
          sum1 = 0._R8KIND
          sum2 = 0._R8KIND
          sumn1 = 0._R8KIND
          sumn2 = 0._R8KIND
          do n2 = nbr,1,-1
             N_r = N0_r*Dr(n2)**mu_r*dexp(-lamr*Dr(n2))*dtr(n2)
             vol = massr(n2)*orho_w
             prob = max(0._R8KIND, 1._R8KIND - dexp(-120._R8KIND*vol*5.2D-4 * Texp))
             if(massr(n2) .lt. xm0g) then
                sumn1 = sumn1 + prob*N_r
                sum1 = sum1 + prob*N_r*massr(n2)
             else
                sumn2 = sumn2 + prob*N_r
                sum2 = sum2 + prob*N_r*massr(n2)
             endif
             if((sum1+sum2).ge.r_r(i)) exit
          enddo
          tpi_qrfz(i,j,k,m) = sum1
          tni_qrfz(i,j,k,m) = sumn1
          tpg_qrfz(i,j,k,m) = sum2
          tnr_qrfz(i,j,k,m) = sumn2
       enddo
    enddo

    do j = 1,nbc
       nu_c = min(15,nint(1000.E6/t_Nc(j)) + 2)
       do i = 1,ntb_c
          lamc = (t_Nc(j)*am_r* ccg(2,nu_c) * ocg1(nu_c) / r_c(i))**obmr
          N0_c = t_Nc(j)*ocg1(nu_c) * lamc**cce(1,nu_c)
          sum1 = 0._R8KIND
          sumn2 = 0._R8KIND
          do n = nbc,1,-1
             vol = massc(n)*orho_w
             prob = max(0._R8KIND, 1._R8KIND - dexp(-120._R8KIND*vol*5.2D-4 * Texp))
             N_c = N0_c*Dc(n)**nu_c*exp(-lamc*Dc(n))*dtc(n)
             sumn2 = min(t_Nc(j), sumn2 + prob*N_c)
             sum1 = sum1 + prob*N_c*massc(n)
             if(sum1 .ge. r_c(i)) exit
          enddo
          tpi_qcfz(i,j,k,m) = sum1
          tni_qcfz(i,j,k,m) = sumn2
       enddo
    enddo
 enddo

 end subroutine freezeH2O

!==================================================================================================================
!..Cloud ice converting to snow since portion greater than min snow
!.. size.  Given cloud ice content (kg/m**3), number concentration
!.. (#/m**3) and gamma shape parameter, mu_i, break the distrib into
!.. bins and figure out the mass/number of ice with sizes larger than
!.. D0s.  Also, compute incomplete gamma function for the integration
!.. of ice depositional growth from diameter=0 to D0s.  Amount of
!.. ice depositional growth is this portion of distrib while larger
!.. diameters contribute to snow growth (as in Harrington et al. 1995).
!+---+-----------------------------------------------------------------+

      subroutine qi_aut_qs

      implicit none

!..Local variables
      INTEGER:: i, j, n2
      DOUBLE PRECISION, DIMENSION(nbi):: N_i
      DOUBLE PRECISION:: N0_i, lami, Di_mean, t1, t2
      REAL:: xlimit_intg

!+---+

      do j = 1, ntb_i1
         do i = 1, ntb_i
            lami = (am_i*cig(2)*oig1*Nt_i(j)/r_i(i))**obmi
            Di_mean = (bm_i + mu_i + 1.) / lami
            N0_i = Nt_i(j)*oig1 * lami**cie(1)
            t1 = 0.0d0
            t2 = 0.0d0
            if (SNGL(Di_mean) .gt. 5.*D0s) then
             t1 = r_i(i)
             t2 = Nt_i(j)
             tpi_ide(i,j) = 0.0D0
            elseif (SNGL(Di_mean) .lt. D0i) then
             t1 = 0.0D0
             t2 = 0.0D0
             tpi_ide(i,j) = 1.0D0
            else
             xlimit_intg = lami*D0s
             tpi_ide(i,j) = GAMMP(mu_i+2.0, xlimit_intg) * 1.0D0
             do n2 = 1, nbi
               N_i(n2) = N0_i*Di(n2)**mu_i * DEXP(-lami*Di(n2))*dti(n2)
               if (Di(n2).ge.D0s) then
                  t1 = t1 + N_i(n2) * am_i*Di(n2)**bm_i
                  t2 = t2 + N_i(n2)
               endif
             enddo
            endif
            tps_iaus(i,j) = t1
            tni_iaus(i,j) = t2
         enddo
      enddo

      end subroutine qi_aut_qs

!==================================================================================================================
!..Variable collision efficiency for rain collecting cloud water using method of Beard and Grover, 1974 if
!  a/A less than 0.25; otherwise uses polynomials to get close match of Pruppacher & Klett Fig 14-9.

 subroutine table_Efrw
 implicit none
!==================================================================================================================

!... local variables:
 integer:: i,j
 real(kind=R8KIND):: vtr,stokes,reynolds,Ef_rw
 real(kind=R8KIND):: p,yc0,F,G,H,z,K0,X

!------------------------------------------------------------------------------------------------------------------

 do j = 1,nbc
    do i = 1,nbr
       Ef_rw = 0.0
       p = Dc(j)/Dr(i)
       if(Dr(i).lt.50.E-6 .or. Dc(j).lt.3.E-6) then
          t_Efrw(i,j) = 0.0
       elseif (p.gt.0.25) then
          X = Dc(j)*1.D6
          if(Dr(i) .lt. 75.e-6) then
             Ef_rw = 0.026794*X - 0.20604
          elseif(Dr(i) .lt. 125.e-6) then
             Ef_rw = -0.00066842*X*X + 0.061542*X - 0.37089
          elseif(Dr(i) .lt. 175.e-6) then
             Ef_rw = 4.091e-06*X*X*X*X - 0.00030908*X*X*X               &
                   + 0.0066237*X*X - 0.0013687*X - 0.073022
          elseif(Dr(i) .lt. 250.e-6) then
             Ef_rw = 9.6719e-5*X*X*X - 0.0068901*X*X + 0.17305*X        &
                   - 0.65988
          elseif(Dr(i) .lt. 350.e-6) then
             Ef_rw = 9.0488e-5*X*X*X - 0.006585*X*X + 0.16606*X         &
                   - 0.56125
          else
             Ef_rw = 0.00010721*X*X*X - 0.0072962*X*X + 0.1704*X        &
                   - 0.46929
          endif
       else
          vtr = -0.1021 + 4.932E3*Dr(i) - 0.9551E6*Dr(i)*Dr(i) &
              + 0.07934E9*Dr(i)*Dr(i)*Dr(i) &
              - 0.002362E12*Dr(i)*Dr(i)*Dr(i)*Dr(i)
          stokes = Dc(j)*Dc(j)*vtr*rho_w/(9.*1.718E-5*Dr(i))
          reynolds = 9.*stokes/(p*p*rho_w)

          F = dlog(reynolds)
          G = -0.1007D0 - 0.358D0*F + 0.0261D0*F*F
          K0 = dexp(G)
          z = dlog(stokes/(K0+1.D-15))
          H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
          yc0 = 2.0D0/PI * atan(H)
          Ef_rw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))
       endif

       t_Efrw(i,j) = max(0._RKIND,min(sngl(Ef_rw),0.95))
    enddo
 enddo

 end subroutine table_Efrw

!==================================================================================================================
!..Variable collision efficiency for snow collecting cloud water using method of Wang and Ji, 2000 except
!  equate melted snow diameter to their "effective collision cross-section."

 subroutine table_Efsw
 implicit none
!==================================================================================================================

!... local variables:
 integer:: i,j
 real(kind=R8KIND):: Ds_m,vts,vtc,stokes,reynolds,Ef_sw
 real(kind=R8KIND):: p,yc0,F,G,H,z,K0

!------------------------------------------------------------------------------------------------------------------

 do j = 1,nbc
    vtc = 1.19D4 * (1.0D4*Dc(j)*Dc(j)*0.25D0)
    do i = 1,nbs
       vts = av_s*Ds(i)**bv_s * dexp(-fv_s*Ds(i)) - vtc
       Ds_m = (am_s*Ds(i)**bm_s / am_r)**obmr
       p = Dc(j)/Ds_m
       if(p.gt.0.25 .or. Ds(i).lt.D0s .or. Dc(j).lt.6.E-6 .or. vts.lt.1.E-3) then
          t_Efsw(i,j) = 0.0
       else
          stokes = Dc(j)*Dc(j)*vts*rho_w/(9.*1.718E-5*Ds_m)
          reynolds = 9.*stokes/(p*p*rho_w)

          F = dlog(reynolds)
          G = -0.1007D0 - 0.358D0*F + 0.0261D0*F*F
          K0 = dexp(G)
          z = dlog(stokes/(K0+1.D-15))
          H = 0.1465D0 + 1.302D0*z - 0.607D0*z*z + 0.293D0*z*z*z
          yc0 = 2.0D0/PI * atan(H)
          Ef_sw = (yc0+p)*(yc0+p) / ((1.+p)*(1.+p))

          t_Efsw(i,j) = max(0._RKIND,min(sngl(Ef_sw),0.95))
       endif
    enddo
 enddo

 end subroutine table_Efsw

!==================================================================================================================
!..Function to compute collision efficiency of collector species (rain,snow, graupel) of aerosols.  Follows Wang
!  et al, 2010, ACP, which follows Slinn (1983).

 real function Eff_aero(D, Da, visc,rhoa,Temp,species)
 implicit none
!==================================================================================================================

!... local variables:
 character(LEN=1):: species
 real:: D, Da, visc, rhoa, Temp
 real:: aval, Cc, diff, Re, Sc, St, St2, vt, Eff
 real, parameter:: boltzman = 1.3806503E-23
 real, parameter:: meanPath = 0.0256E-6

!------------------------------------------------------------------------------------------------------------------

 vt = 1.
 if(species .eq. 'r') then
    vt = -0.1021 + 4.932E3*D - 0.9551E6*D*D      &
       + 0.07934E9*D*D*D - 0.002362E12*D*D*D*D
 elseif(species .eq. 's') then
    vt = av_s*D**bv_s
 elseif(species .eq. 'g') then
    vt = av_g*D**bv_g
 endif

 Cc   = 1. + 2.*meanPath/Da *(1.257+0.4*exp(-0.55*Da/meanPath))
 diff = boltzman*Temp*Cc/(3.*PI*visc*Da)

 Re   = 0.5*rhoa*D*vt/visc
 Sc   = visc/(rhoa*diff)

 St   = Da*Da*vt*1000./(9.*visc*D)
 aval = 1.+log(1.+Re)
 St2  = (1.2 + 1./12.*aval)/(1.+aval)

 Eff = 4./(Re*Sc) * (1. + 0.4*sqrt(Re)*Sc**0.3333  &
     + 0.16*sqrt(Re)*sqrt(Sc))                     &
     + 4.*Da/D * (0.02 + Da/D*(1.+2.*sqrt(Re)))

 if(St.gt.St2) Eff = Eff  + ( (St-St2)/(St-St2+0.666667))**1.5
 Eff_aero = max(1.E-5,min(Eff,1.0))

 end function Eff_aero

!==================================================================================================================
!..Integrate rain size distribution from zero to D-star to compute the number of drops smaller than D-star
!  that evaporate in a single timestep.  Drops larger than D-star dont evaporate entirely so do not affect
!  number concentration.

 subroutine table_dropEvap
 implicit none
!==================================================================================================================

!... local variables:
 integer:: i,j,k,n
 integer:: nu_c

 real(kind=R8KIND):: summ,summ2,lamc,N0_c
 real(kind=R8KIND),dimension(nbc):: N_c,massc
!real(kind=R8KIND:: Nt_r, N0, lam_exp, lam
!real(kind=RKIND):: xlimit_intg

!------------------------------------------------------------------------------------------------------------------

 do n = 1,nbc
    massc(n) = am_r*Dc(n)**bm_r
 enddo

 do k = 1,nbc
    nu_c = min(15,nint(1000.E6/t_Nc(k)) + 2)
    do j = 1,ntb_c
       lamc = (t_Nc(k)*am_r* ccg(2,nu_c)*ocg1(nu_c) / r_c(j))**obmr
       N0_c = t_Nc(k)*ocg1(nu_c) * lamc**cce(1,nu_c)
       do i = 1,nbc
!-GT      tnc_wev(i,j,k) = gammp(nu_c+1.,sngl(Dc(i)*lamc))*t_Nc(k)
          N_c(i) = N0_c* Dc(i)**nu_c*exp(-lamc*Dc(i))*dtc(i)
          summ = 0._R8KIND
          summ2 = 0._R8KIND
          do n = 1,i
             summ = summ + massc(n)*N_c(n)
             summ2 = summ2 + N_c(n)
          enddo
          tpc_wev(i,j,k) = summ
          tnc_wev(i,j,k) = summ2
       enddo
    enddo
 enddo

 end subroutine table_dropEvap

!==================================================================================================================
!..Retrieve fraction of CCN that gets activated given the model temp, vertical velocity, and available CCN
!  concentration.  The lookup table (read from external file) has CCN concentration varying the quickest, 
!  then updraft, then temperature, then mean aerosol radius, and finally hygroscopicity, kappa.
!..TO_DO ITEM:  For radiation cooling producing fog, in which case the updraft velocity could easily be negative,
!  we could use the temp and its tendency to diagnose a pretend postive updraft velocity.

 real function activ_ncloud(Tt,Ww,NCCN)
 implicit none
!==================================================================================================================

!... input arguments:
 real(kind=RKIND),intent(in):: Tt,Ww,NCCN

!... local variables:
 integer:: i,j,k,l,m,n
 real(kind=RKIND):: n_local,w_local
 real(kind=RKIND):: A,B,C,D,t,u,x1,x2,y1,y2,nx,wy,fraction

!------------------------------------------------------------------------------------------------------------------

!ta_Na = (/10.0, 31.6, 100.0, 316.0, 1000.0, 3160.0, 10000.0/)  ntb_arc
!ta_Ww = (/0.01, 0.0316, 0.1, 0.316, 1.0, 3.16, 10.0, 31.6, 100.0/)  ntb_arw
!ta_Tk = (/243.15, 253.15, 263.15, 273.15, 283.15, 293.15, 303.15/)  ntb_art
!ta_Ra = (/0.01, 0.02, 0.04, 0.08, 0.16/)  ntb_arr
!ta_Ka = (/0.2, 0.4, 0.6, 0.8/)  ntb_ark


 n_local = NCCN * 1.E-6
 w_local = Ww

 if(n_local .ge. ta_Na(ntb_arc)) then
    n_local = ta_Na(ntb_arc) - 1.0
 elseif(n_local .le. ta_Na(1)) then
    n_local = ta_Na(1) + 1.0
 endif
 do n = 2,ntb_arc
    if(n_local.ge.ta_Na(n-1) .and. n_local.lt.ta_Na(n)) goto 8003
 enddo
 8003 continue
 i = n
 x1 = log(ta_Na(i-1))
 x2 = log(ta_Na(i))


 if(w_local .ge. ta_Ww(ntb_arw)) then
    w_local = ta_Ww(ntb_arw) - 1.0
 elseif(w_local .le. ta_Ww(1)) then
    w_local = ta_Ww(1) + 0.001
 endif
 do n = 2,ntb_arw
    if(w_local.ge.ta_Ww(n-1) .and. w_local.lt.ta_Ww(n)) goto 8005
 enddo
 8005 continue
 j = n
 y1 = log(ta_Ww(j-1))
 y2 = log(ta_Ww(j))

 k = max(1,min(nint( (Tt - ta_Tk(1))*0.1) + 1, ntb_art))


!..The next two values are indexes of mean aerosol radius and hygroscopicity.  Currently these are constant
!  but a future version should implement other variables to allow more freedom such as at least simple separation
!  of tiny size sulfates from larger sea salts:
 l = 3
 m = 2

 A = tnccn_act(i-1,j-1,k,l,m)
 B = tnccn_act(i,j-1,k,l,m)
 C = tnccn_act(i,j,k,l,m)
 D = tnccn_act(i-1,j,k,l,m)
 nx = log(n_local)
 wy = log(w_local)

 t = (nx-x1)/(x2-x1)
 u = (wy-y1)/(y2-y1)
!t = (n_local-ta(Na(i-1))/(ta_Na(i)-ta_Na(i-1))
!u = (w_local-ta_Ww(j-1))/(ta_Ww(j)-ta_Ww(j-1))

 fraction = (1.0-t)*(1.0-u)*A + t*(1.0-u)*B + t*u*C + (1.0-t)*u*D

!if(NCCN*fraction .gt. 0.75*Nt_c_max) then
!   write(*,*) ' DEBUG-GT ', n_local, w_local, Tt, i, j, k
!endif

 activ_ncloud = NCCN*fraction

 end function activ_ncloud

!==================================================================================================================
 real function iceDeMott(tempc,qv,qvs,qvsi,rho,nifa)
 implicit none
!==================================================================================================================

!... input arguments:
 real(kind=RKIND),intent(in):: tempc,qv,qvs,qvsi,rho,nifa

!... local variables:
 real(kind=RKIND):: satw, sati, siw, p_x, si0x, dtt, dsi, dsw, dab, fc, hx
 real(kind=RKIND):: ntilde, n_in, nmax, nhat, mux, xni, nifa_cc
 real(kind=RKIND),parameter:: p_c1    = 1000.
 real(kind=RKIND),parameter:: p_rho_c = 0.76
 real(kind=RKIND),parameter:: p_alpha = 1.0
 real(kind=RKIND),parameter:: p_gam   = 2.
 real(kind=RKIND),parameter:: delT    = 5.
 real(kind=RKIND),parameter:: T0x     = -40.
 real(kind=RKIND),parameter:: Sw0x    = 0.97
 real(kind=RKIND),parameter:: delSi   = 0.1
 real(kind=RKIND),parameter:: hdm     = 0.15
 real(kind=RKIND),parameter:: p_psi   = 0.058707*p_gam/p_rho_c
 real(kind=RKIND),parameter:: aap     = 1.
 real(kind=RKIND),parameter:: bbp     = 0.
 real(kind=RKIND),parameter:: y1p     = -35.
 real(kind=RKIND),parameter:: y2p     = -25.
 real(kind=RKIND),parameter:: rho_not0 = 101325./(287.05*273.15)

!------------------------------------------------------------------------------------------------------------------

 xni = 0._RKIND

!satw = qv/qvs
!sati = qv/qvsi
!siw = qvs/qvsi
!p_x = -1.0261+(3.1656e-3*tempc)+(5.3938e-4*(tempc*tempc)) &
!    +  (8.2584e-6*(tempc*tempc*tempc))
!si0x = 1.+(10.**p_x)
!if(sati.ge.si0x .and. satw.lt.0.985) then
!   dtt = delta_p (tempc, T0x, T0x+delT, 1., hdm)
!   dsi = delta_p (sati, Si0x, Si0x+delSi, 0., 1.)
!   dsw = delta_p (satw, Sw0x, 1., 0., 1.)
!   fc = dtt*dsi*0.5
!   hx = min(fc+((1.-fc)*dsw), 1.)
!   ntilde = p_c1*p_gam*((exp(12.96*(sati-1.1)))**0.3) / p_rho_c
!   if(tempc .le. y1p) then
!      n_in = ntilde
!   elseif(tempc .ge. y2p) then
!      n_in = p_psi*p_c1*exp(12.96*(sati-1.)-0.639)
!   else
!      if(tempc .le. -30.) then
!         nmax = p_c1*p_gam*(exp(12.96*(siw-1.1)))**0.3/p_rho_c
!      else
!         nmax = p_psi*p_c1*exp(12.96*(siw-1.)-0.639)
!      endif
!      ntilde = min(ntilde, nmax)
!      nhat = min(p_psi*p_c1*exp(12.96*(sati-1.)-0.639), nmax)
!      dab = delta_p(tempc,y1p,y2p,aap,bbp)
!      n_in = MIN(nhat*(ntilde/nhat)**dab, nmax)
!   endif
!   mux = hx*p_alpha*n_in*rho
!   xni = mux*((6700.*nifa)-200.)/((6700.*5.E5)-200.)
!elseif (satw.ge.0.985 .and. tempc.gt.HGFR-273.15) then
    nifa_cc = max(0.5, nifa*RHO_NOT0*1.E-6/rho)
!   xni = 3.*nifa_cc**(1.25)*exp((0.46*(-tempc))-11.6)              !  [DeMott, 2015]
    xni = (5.94e-5*(-tempc)**3.33)                                & !  [DeMott, 2010]
        * (nifa_cc**((-0.0264*(tempc))+0.0033))
    xni = xni*rho/RHO_NOT0 * 1000.
!endif

 iceDeMott = max(0._RKIND,xni)

 end function iceDeMott

!==================================================================================================================
!..Newer research since Koop et al (2001) suggests that the freezing rate should be lower than original paper,
!  so J_rate is reduced by two orders of magnitude.

 real function iceKoop(temp,qv,qvs,naero,dt)
 implicit none
!==================================================================================================================

!... input arguments:
 real(kind=RKIND),intent(in):: temp,qv,qvs,naero,dt

!... local variables:
 real(kind=RKIND):: mu_diff,a_w_i,delta_aw,log_J_rate,J_rate,prob_h,satw
 real(kind=RKIND):: xni

!------------------------------------------------------------------------------------------------------------------

 xni = 0._RKIND
 satw = qv/qvs
 mu_diff = 210368.0 + (131.438*temp) - (3.32373E6/temp) - (41729.1*alog(temp))
 a_w_i   = exp(mu_diff/(R_uni*temp))
 delta_aw   = satw - a_w_i
 log_J_rate = -906.7 + (8502.0*delta_aw)           &
            - (26924.0*delta_aw*delta_aw)          &
            + (29180.0*delta_aw*delta_aw*delta_aw)
 log_J_rate = min(20.0, log_J_rate)
 J_rate     = 10.**log_J_rate                            ! cm-3 s-1
 prob_h     = min(1.-exp(-J_rate*ar_volume*dt),1._RKIND)
 if(prob_h .gt. 0._RKIND) then
    xni = min(prob_h*naero,1000.E3)
 endif

 iceKoop = max(0._RKIND,xni)

 end function iceKoop

!==================================================================================================================
!.. Helper routine for Phillips et al (2008) ice nucleation.  Trude

 real function delta_p(yy,y1,y2,aa,bb)
 implicit none
!==================================================================================================================

!... input arguments:
 real(kind=RKIND),intent(in):: yy,y1,y2,aa,bb

!... local variables:
 real(kind=RKIND):: dab,A,B,a0,a1,a2,a3

!------------------------------------------------------------------------------------------------------------------

 A  = 6.*(aa-bb)/((y2-y1)*(y2-y1)*(y2-y1))
 B  = aa+(A*y1*y1*y1/6.)-(A*y1*y1*y2*0.5)
 a0 = B
 a1 = A*y1*y2
 a2 = -A*(y1+y2)*0.5
 a3 = A/3.

 if(yy.le.y1) then 
    dab = aa
 elseif(yy.ge.y2) then
    dab = bb
 else
    dab = a0+(a1*yy)+(a2*yy*yy)+(a3*yy*yy*yy)
 endif

 if(dab.lt.aa) then 
    dab = aa
 endif
 if(dab.gt.bb) then 
    dab = bb
 endif

 delta_p = dab

 end function delta_p 

!==================================================================================================================
!..Compute _radiation_ effective radii of cloud water, ice, and snow.
!.. These are entirely consistent with microphysics assumptions, not
!.. constant or otherwise ad hoc as is internal to most radiation
!.. schemes.  Since only the smallest snowflakes should impact
!.. radiation, compute from first portion of complicated Field number
!.. distribution, not the second part, which is the larger sizes.
!+---+-----------------------------------------------------------------+

      subroutine calc_effectRad (t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d,   &
     &                re_qc1d, re_qi1d, re_qs1d, kts, kte)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
     &                    t1d, p1d, qv1d, qc1d, nc1d, qi1d, ni1d, qs1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: re_qc1d, re_qi1d, re_qs1d
!..Local variables
      INTEGER:: k
      REAL, DIMENSION(kts:kte):: rho, rc, nc, ri, ni, rs
      REAL:: smo2, smob, smoc
      REAL:: tc0, loga_, a_, b_
      DOUBLE PRECISION:: lamc, lami
      LOGICAL:: has_qc, has_qi, has_qs
      INTEGER:: inu_c
      real, dimension(15), parameter:: g_ratio = (/24,60,120,210,336,   &
     &                504,720,990,1320,1716,2184,2730,3360,4080,4896/)

      has_qc = .false.
      has_qi = .false.
      has_qs = .false.

      do k = kts, kte
         rho(k) = 0.622*p1d(k)/(R*t1d(k)*(qv1d(k)+0.622))
         rc(k) = MAX(R1, qc1d(k)*rho(k))
         nc(k) = MAX(2., MIN(nc1d(k)*rho(k), Nt_c_max))
         if (.NOT. is_aerosol_aware) nc(k) = Nt_c
         if (rc(k).gt.R1 .and. nc(k).gt.R2) has_qc = .true.
         ri(k) = MAX(R1, qi1d(k)*rho(k))
         ni(k) = MAX(R2, ni1d(k)*rho(k))
         if (ri(k).gt.R1 .and. ni(k).gt.R2) has_qi = .true.
         rs(k) = MAX(R1, qs1d(k)*rho(k))
         if (rs(k).gt.R1) has_qs = .true.
      enddo

      if (has_qc) then
      do k = kts, kte
         re_qc1d(k) = 2.49E-6
         if (rc(k).le.R1 .or. nc(k).le.R2) CYCLE
         if (nc(k).lt.100) then
            inu_c = 15
         elseif (nc(k).gt.1.E10) then
            inu_c = 2
         else
            inu_c = MIN(15, NINT(1000.E6/nc(k)) + 2)
         endif
         lamc = (nc(k)*am_r*g_ratio(inu_c)/rc(k))**obmr
         re_qc1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+inu_c)/lamc), 50.E-6))
      enddo
      endif

      if (has_qi) then
      do k = kts, kte
         re_qi1d(k) = 2.49E-6
         if (ri(k).le.R1 .or. ni(k).le.R2) CYCLE
         lami = (am_i*cig(2)*oig1*ni(k)/ri(k))**obmi
         re_qi1d(k) = MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_i)/lami), 125.E-6))
      enddo
      endif

      if (has_qs) then
      do k = kts, kte
         re_qs1d(k) = 4.99E-6
         if (rs(k).le.R1) CYCLE
         tc0 = MIN(-0.1, t1d(k)-273.15)
         smob = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2 = smob
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
     &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
     &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
     &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
     &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
     &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
     &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
     &         + sb(10)*bm_s*bm_s*bm_s
            smo2 = (smob/a_)**(1./b_)
         endif
!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
     &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
     &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
     &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc = a_ * smo2**b_
         re_qs1d(k) = MAX(5.01E-6, MIN(0.5*(smoc/smob), 999.E-6))
      enddo
      endif

      end subroutine calc_effectRad

!+---+-----------------------------------------------------------------+
!..Compute radar reflectivity assuming 10 cm wavelength radar and using
!.. Rayleigh approximation.  Only complication is melted snow/graupel
!.. which we treat as water-coated ice spheres and use Uli Blahak's
!.. library of routines.  The meltwater fraction is simply the amount
!.. of frozen species remaining from what initially existed at the
!.. melting level interface.
!+---+-----------------------------------------------------------------+

      subroutine calc_refl10cm (qv1d, qc1d, qr1d, nr1d, qs1d, qg1d,     &
                          t1d, p1d, dBZ, kts, kte, ii, jj)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                          qv1d, qc1d, qr1d, nr1d, qs1d, qg1d, t1d, p1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ
!     REAL, DIMENSION(kts:kte), INTENT(INOUT):: vt_dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho, rhof
      REAL, DIMENSION(kts:kte):: rc, rr, nr, rs, rg

      DOUBLE PRECISION, DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g
      REAL, DIMENSION(kts:kte):: mvd_r
      REAL, DIMENSION(kts:kte):: smob, smo2, smoc, smoz
      REAL:: oM3, M0, Mrat, slam1, slam2, xDs
      REAL:: ils1, ils2, t1_vts, t2_vts, t3_vts, t4_vts
      REAL:: vtr_dbz_wt, vts_dbz_wt, vtg_dbz_wt

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

      DOUBLE PRECISION:: N0_exp, N0_min, lam_exp, lamr, lamg
      REAL:: a_, b_, loga_, tc0
      DOUBLE PRECISION:: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0, kbot, n
      LOGICAL:: melti
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

      DOUBLE PRECISION:: cback, x, eta, f_d
      REAL:: xslw1, ygra1, zans1

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         rhof(k) = SQRT(RHO_NOT/rho(k))
         rc(k) = MAX(R1, qc1d(k)*rho(k))
         if (qr1d(k) .gt. R1) then
            rr(k) = qr1d(k)*rho(k)
            nr(k) = MAX(R2, nr1d(k)*rho(k))
            lamr = (am_r*crg(3)*org2*nr(k)/rr(k))**obmr
            ilamr(k) = 1./lamr
            N0_r(k) = nr(k)*org2*lamr**cre(2)
            mvd_r(k) = (3.0 + mu_r + 0.672) * ilamr(k)
            L_qr(k) = .true.
         else
            rr(k) = R1
            nr(k) = R1
            mvd_r(k) = 50.E-6
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R2) then
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R2) then
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif
      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         smo2(k) = 0.
         smob(k) = 0.
         smoc(k) = 0.
         smoz(k) = 0.
      enddo
      if (ANY(L_qs .eqv. .true.)) then
      do k = kts, kte
         if (.not. L_qs(k)) CYCLE
         tc0 = MIN(-0.1, temp(k)-273.15)
         smob(k) = rs(k)*oams

!..All other moments based on reference, 2nd moment.  If bm_s.ne.2,
!.. then we must compute actual 2nd moment and use as reference.
         if (bm_s.gt.(2.0-1.e-3) .and. bm_s.lt.(2.0+1.e-3)) then
            smo2(k) = smob(k)
         else
            loga_ = sa(1) + sa(2)*tc0 + sa(3)*bm_s &
     &         + sa(4)*tc0*bm_s + sa(5)*tc0*tc0 &
     &         + sa(6)*bm_s*bm_s + sa(7)*tc0*tc0*bm_s &
     &         + sa(8)*tc0*bm_s*bm_s + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*bm_s*bm_s*bm_s
            a_ = 10.0**loga_
            b_ = sb(1) + sb(2)*tc0 + sb(3)*bm_s &
     &         + sb(4)*tc0*bm_s + sb(5)*tc0*tc0 &
     &         + sb(6)*bm_s*bm_s + sb(7)*tc0*tc0*bm_s &
     &         + sb(8)*tc0*bm_s*bm_s + sb(9)*tc0*tc0*tc0 &
     &         + sb(10)*bm_s*bm_s*bm_s
            smo2(k) = (smob(k)/a_)**(1./b_)
         endif

!..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(1) &
     &         + sa(4)*tc0*cse(1) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(1)*cse(1) + sa(7)*tc0*tc0*cse(1) &
     &         + sa(8)*tc0*cse(1)*cse(1) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(1)*cse(1)*cse(1)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(1) + sb(4)*tc0*cse(1) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(1)*cse(1) &
     &        + sb(7)*tc0*tc0*cse(1) + sb(8)*tc0*cse(1)*cse(1) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(1)*cse(1)*cse(1)
         smoc(k) = a_ * smo2(k)**b_

!..Calculate bm_s*2 (th) moment.  Useful for reflectivity.
         loga_ = sa(1) + sa(2)*tc0 + sa(3)*cse(3) &
     &         + sa(4)*tc0*cse(3) + sa(5)*tc0*tc0 &
     &         + sa(6)*cse(3)*cse(3) + sa(7)*tc0*tc0*cse(3) &
     &         + sa(8)*tc0*cse(3)*cse(3) + sa(9)*tc0*tc0*tc0 &
     &         + sa(10)*cse(3)*cse(3)*cse(3)
         a_ = 10.0**loga_
         b_ = sb(1)+ sb(2)*tc0 + sb(3)*cse(3) + sb(4)*tc0*cse(3) &
     &        + sb(5)*tc0*tc0 + sb(6)*cse(3)*cse(3) &
     &        + sb(7)*tc0*tc0*cse(3) + sb(8)*tc0*cse(3)*cse(3) &
     &        + sb(9)*tc0*tc0*tc0 + sb(10)*cse(3)*cse(3)*cse(3)
         smoz(k) = a_ * smo2(k)**b_
      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

      if (ANY(L_qg .eqv. .true.)) then
      N0_min = gonv_max
      k_0 = kts
      do k = kte, kts, -1
         if (temp(k).ge.270.65) k_0 = MAX(k_0, k)
      enddo
      do k = kte, kts, -1
         if (k.gt.k_0 .and. L_qr(k) .and. mvd_r(k).gt.100.E-6) then
            xslw1 = 4.01 + alog10(mvd_r(k))
         else
            xslw1 = 0.01
         endif
         ygra1 = 4.31 + alog10(max(5.E-5, rg(k)))
         zans1 = 3.1 + (100./(300.*xslw1*ygra1/(10./xslw1+1.+0.25*ygra1)+30.+10.*ygra1))
         N0_exp = 10.**(zans1)
         N0_exp = MAX(DBLE(gonv_min), MIN(N0_exp, DBLE(gonv_max)))
         N0_min = MIN(N0_exp, N0_min)
         N0_exp = N0_min
         lam_exp = (N0_exp*am_g*cgg(1)/rg(k))**oge1
         lamg = lam_exp * (cgg(3)*ogg2*ogg1)**obmg
         ilamg(k) = 1./lamg
         N0_g(k) = N0_exp/(cgg(2)*lam_exp) * lamg**cge(2)
      enddo
      endif

!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
     &                            .and. (L_qs(k+1).or.L_qg(k+1)) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*crg(4)*ilamr(k)**cre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
     &                           * (am_s/900.0)*(am_s/900.0)*smoz(k)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
     &                              * (am_g/900.0)*(am_g/900.0)         &
     &                              * N0_g(k)*cgg(4)*ilamg(k)**cge(4)
      enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (.not. iiwarm .and. melti .and. k_0.ge.2) then
       do k = k_0-1, kts, -1

!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.05d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           oM3 = 1./smoc(k)
           M0 = (smob(k)*oM3)
           Mrat = smob(k)*M0*M0*M0
           slam1 = M0 * Lam0
           slam2 = M0 * Lam1
           do n = 1, nrbins
              x = am_s * xxDs(n)**bm_s
              call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
     &              fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
     &              CBACK, mixingrulestring_s, matrixstring_s,          &
     &              inclusionstring_s, hoststring_s,                    &
     &              hostmatrixstring_s, hostinclusionstring_s)
              f_d = Mrat*(Kap0*DEXP(-slam1*xxDs(n))                     &
     &              + Kap1*(M0*xxDs(n))**mu_s * DEXP(-slam2*xxDs(n)))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

!..Reflectivity contributed by melting graupel

          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.05d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = am_g * xxDg(n)**bm_g
              call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
     &              fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
     &              CBACK, mixingrulestring_g, matrixstring_g,          &
     &              inclusionstring_g, hoststring_g,                    &
     &              hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**mu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo


!..Reflectivity-weighted terminal velocity (snow, rain, graupel, mix).
!     do k = kte, kts, -1
!        vt_dBZ(k) = 1.E-3
!        if (rs(k).gt.R2) then
!         Mrat = smob(k) / smoc(k)
!         ils1 = 1./(Mrat*Lam0 + fv_s)
!         ils2 = 1./(Mrat*Lam1 + fv_s)
!         t1_vts = Kap0*csg(5)*ils1**cse(5)
!         t2_vts = Kap1*Mrat**mu_s*csg(11)*ils2**cse(11)
!         ils1 = 1./(Mrat*Lam0)
!         ils2 = 1./(Mrat*Lam1)
!         t3_vts = Kap0*csg(6)*ils1**cse(6)
!         t4_vts = Kap1*Mrat**mu_s*csg(12)*ils2**cse(12)
!         vts_dbz_wt = rhof(k)*av_s * (t1_vts+t2_vts)/(t3_vts+t4_vts)
!         if (temp(k).ge.273.15 .and. temp(k).lt.275.15) then
!            vts_dbz_wt = vts_dbz_wt*1.5
!         elseif (temp(k).ge.275.15) then
!            vts_dbz_wt = vts_dbz_wt*2.0
!         endif
!        else
!         vts_dbz_wt = 1.E-3
!        endif

!        if (rr(k).gt.R1) then
!         lamr = 1./ilamr(k)
!         vtr_dbz_wt = rhof(k)*av_r*crg(13)*(lamr+fv_r)**(-cre(13))      &
!    &               / (crg(4)*lamr**(-cre(4)))
!        else
!         vtr_dbz_wt = 1.E-3
!        endif

!        if (rg(k).gt.R2) then
!         lamg = 1./ilamg(k)
!         vtg_dbz_wt = rhof(k)*av_g*cgg(5)*lamg**(-cge(5))               &
!    &               / (cgg(4)*lamg**(-cge(4)))
!        else
!         vtg_dbz_wt = 1.E-3
!        endif

!        vt_dBZ(k) = (vts_dbz_wt*ze_snow(k) + vtr_dbz_wt*ze_rain(k)      &
!    &                + vtg_dbz_wt*ze_graupel(k))                        &
!    &                / (ze_rain(k)+ze_snow(k)+ze_graupel(k))
!     enddo

      end subroutine calc_refl10cm
!
!+---+-----------------------------------------------------------------+

!+---+-----------------------------------------------------------------+
END MODULE mp_thompson
!+---+-----------------------------------------------------------------+
