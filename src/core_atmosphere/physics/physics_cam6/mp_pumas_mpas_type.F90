!======================================================================================================================
 module mp_pumas_type
!======================================================================================================================
 use mpas_kind_types

 type,public:: mp_pumas
    integer:: mgncol
    integer:: nlev

    !--- variables used in the initialization of PUMAS cloud microphysics scheme:
    logical,pointer:: &
       micro_mg_do_graupel,        &
       micro_mg_do_hail,           &
       microp_uniform,             &
       do_cldice,                  &
       use_hetfrz_classnuc,        &
       remove_supersat,            &
       do_sb_physics

    character(len=16):: &
       micro_mg_precip_frac_method = 'max_overlap' ! type of precipitation fraction method
    real(kind=R8KIND):: &
       micro_mg_berg_eff_factor    = 1.0_R8KIND    ! berg efficiency factor

    !IFS-like switches:
    logical,pointer:: &
       micro_mg_evap_sed_off,      &
       micro_mg_icenuc_rh_off,     &
       micro_mg_icenuc_use_meyers, &
       micro_mg_evap_scl_ifs,      &
       micro_mg_evap_rhthrsh_ifs,  &
       micro_mg_rainfreeze_ifs,    &
       micro_mg_ifs_sed,           &
       micro_mg_precip_fall_corr

    !namelist variables for option to specify constant cloud droplet/ice number:
    logical,pointer:: &
       micro_mg_nccons,            &
       micro_mg_nicons,            &
       micro_mg_ngcons,            &
       micro_mg_nrcons,            &
       micro_mg_nscons

    real(kind=R8KIND):: micro_mg_dcs = -1._R8KIND

    !parameters for specified ice and droplet number concentration
    !note: these are local in-cloud values, not grid-mean
    real(kind=R8KIND):: micro_mg_ncnst = 50.e6_R8KIND    !constant liquid droplet num concentration (m-3)
    real(kind=R8KIND):: micro_mg_ninst = 0.05e6_R8KIND   !ice num concentration when nicons=.true. (m-3)
    real(kind=R8KIND):: micro_mg_nrnst = 0.2e6_R8KIND    !rain  num concentration when nrcons=.true. (m-3)
    real(kind=R8KIND):: micro_mg_nsnst = 0.005e6_R8KIND  !snow num concentration when nscons=.true. (m-3)
    real(kind=R8KIND):: micro_mg_ngnst = 0.0005e6_R8KIND !graupel/hail num concentration when ngcons=.true. (m-3)

    !temporary:
    real(kind=R8KIND):: rhmini = 10._R8KIND


    !--- inputs:
    real(kind=R8KIND),pointer:: delta

    real(kind=R8KIND),dimension(:,:),pointer:: &
       pthick      => null(), &
       press       => null(), &
       temp        => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       q           => null(), &
       qcn         => null(), &
       qrn         => null(), &
       qin         => null(), &
       qsn         => null(), &
       qgn         => null(), &
       ncn         => null(), &
       nrn         => null(), &
       nin         => null(), &
       nsn         => null(), &
       ngn         => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       relvar      => null(), &
       cldn        => null(), &
       liqcldf     => null(), &
       icecldf     => null(), &
       qsatfac     => null(), &
       accre_enhan => null(), &
       naai        => null(), &
       npccn       => null()

    real(kind=R8KIND),dimension(:,:,:),pointer:: &
       nacon       => null(), &
       rndst       => null()

    !--- outputs:
    real(kind=R8KIND),dimension(:,:),pointer:: &
       tlat       => null(), &
       qvlat      => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       effc       => null(), &
       effc_fn    => null(), &
       effi       => null(), &
       sadice     => null(), &
       sadsnow    => null(), &
       prect      => null(), &
       preci      => null(), &
       nevapr     => null(), &
       evapsnow   => null(), &
       am_evp_st  => null(), &
       prain      => null(), &
       prodsnow   => null(), &
       cmeout     => null(), &
       deffi      => null(), &
       pgamrad    => null(), &
       lamcrad    => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       qctend     => null(), &
       qrtend     => null(), &
       qitend     => null(), &
       qstend     => null(), &
       qgtend     => null(), &
       nctend     => null(), &
       nrtend     => null(), &
       nitend     => null(), &
       nstend     => null(), &
       ngtend     => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       qcsedten   => null(), &
       qrsedten   => null(), &
       qisedten   => null(), &
       qssedten   => null(), &
       qgsedten   => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       qrout      => null(), &
       qsout      => null(), &
       qgout      => null(), &
       nrout      => null(), &
       nsout      => null(), &
       ngout      => null(), &
       drout      => null(), &
       dsout      => null(), &
       dgout      => null(), &
       qrout2     => null(), &
       qsout2     => null(), &
       qgout2     => null(), &
       nrout2     => null(), &
       nsout2     => null(), &
       ngout2     => null(), &
       drout2     => null(), &
       dsout2     => null(), &
       dgout2     => null(), &
       umr        => null(), &
       ums        => null(), &
       umg        => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       pratot     => null(), & ! accretion of cloud by rain
       prctot     => null(), & ! autoconversion of cloud to rain
       mnuccctot  => null(), & ! mixing ratio tend due to immersion freezing
       mnuccttot  => null(), & ! mixing ratio tend due to contact freezing
       msacwitot  => null(), & ! mixing ratio tend due to H-M splintering
       psacwstot  => null(), & ! collection of cloud water by snow
       bergstot   => null(), & ! bergeron process on snow
       bergtot    => null(), & ! bergeron process on cloud ice
       melttot    => null(), & ! melting of cloud ice
       meltstot   => null(), & ! melting of snow
       meltgtot   => null(), & ! melting of graupel
       mnudeptot  => null(), & ! deposition nucleation to ice
       homotot    => null(), & ! homogeneous freezing cloud water
       qcrestot   => null(), & ! residual cloud condensation due to removal of excess supersat
       prcitot    => null(), & ! autoconversion of cloud ice to snow
       praitot    => null(), & ! accretion of cloud ice by snow
       qirestot   => null(), & ! residual ice deposition due to removal of excess supersat
       mnuccrtot  => null(), & ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
       mnuccritot => null(), & ! mixing ratio tendency due to heterogeneous freezing of rain to ice (1/s)
       pracstot   => null(), & ! mixing ratio tendency due to accretion of rain by snow (1/s)
       meltsdttot => null(), & ! latent heating rate due to melting of snow  (W/kg)
       frzrdttot  => null(), & ! latent heating rate due to homogeneous freezing of rain (W/kg)
       mnuccdtot  => null(), & ! mass tendency from ice nucleation
       pracgtot   => null(), & ! change in q collection rain by graupel  (precipf)
       psacwgtot  => null(), & ! change in q collection droplets by graupel (lcldm)
       pgsacwtot  => null(), & ! conversion q to graupel due to collection droplets by snow  (lcldm)
       pgracstot  => null(), & ! conversion q to graupel due to collection rain by snow (precipf)
       prdgtot    => null(), & ! dep of graupel (precipf)
       qmultgtot  => null(), & ! change q due to ice mult droplets/graupel  (lcldm)
       qmultrgtot => null(), & ! change q due to ice mult rain/graupel (precipf)
       psacrtot   => null(), & ! conversion due to coll of snow by rain (precipf)
       npracgtot  => null(), & ! change n collection rain by graupel  (precipf)
       nscngtot   => null(), & ! change n conversion to graupel due to collection droplets by snow (lcldm)
       ngracstot  => null(), & ! change n conversion to graupel due to collection rain by snow (precipf)
       nmultgtot  => null(), & ! ice mult due to acc droplets by graupel  (lcldm)
       nmultrgtot => null(), & ! ice mult due to acc rain by graupel  (precipf)
       npsacwgtot => null(), & ! change n collection droplets by graupel (lcldm?)
       refl       => null(), & ! analytic radar reflectivity
       arefl      => null(), & ! average reflectivity will zero points outside valid range
       areflz     => null(), & ! average reflectivity in z.
       frefl      => null(), & ! fractional occurrence of radar reflectivity
       csrfl      => null(), & ! cloudsat reflectivity
       acsrfl     => null(), & ! cloudsat average
       fcsrfl     => null(), & ! cloudsat fractional occurrence of radar reflectivity
       rercld     => null(), & ! effective radius calculation for rain + cloud
       ncai       => null(), & ! output number conc of ice nuclei available (1/m3)
       ncal       => null(), & ! output number conc of CCN (1/m3)
       nfice      => null(), & ! fractional occurrence of ice
       qcrat      => null(), & ! limiter for qc process rates (1=no limit --> 0. no qc)
       prer_evap  => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       freqr      => null(), &
       freqs      => null(), &
       freqg      => null(), &
       reff_rain  => null(), &
       reff_snow  => null(), &
       reff_grau  => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       qcsevap    => null(), &
       qisevap    => null(), &
       qvres      => null(), &
       cmeitot    => null(), &
       vtrmc      => null(), &
       vtrmi      => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       qcsinksum_rate1ord => null()

    real(kind=R8KIND),dimension(:,:),pointer:: &
       lflx       => null(), &
       rflx       => null(), &
       iflx       => null(), &
       sflx       => null(), &
       gflx       => null()

    !--- tendencies calculated by external schemes that can replace MG's native process tendencies:
    !used with CARMA cirrus microphysics (or similar external microphysics model):
    real(kind=R8KIND),dimension(:,:),pointer:: &
       tnd_qsnow => null(), &
       tnd_nsnow => null(), &
       re_ice    => null()

    !from external ice nucleation:
    real(kind=R8KIND),dimension(:,:),pointer:: &
       frzimm    => null(), &
       frzcnt    => null(), &
       frzdep    => null()


    contains
!      procedure:: allocate_mem   => mp_pumas_allocate
!      procedure:: deallocate_mem => mp_pumas_deallocate
!      procedure:: MPAS_to_pumas  => MPAS_to_PUMAS
!      procedure:: mp_pumas_allocate
!      procedure:: mp_pumas_deallocate
!      procedure:: MPAS_to_pumas
 end type

!======================================================================================================================
 end module mp_pumas_type
!======================================================================================================================
