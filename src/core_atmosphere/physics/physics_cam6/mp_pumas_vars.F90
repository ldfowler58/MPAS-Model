!======================================================================================================================
 module pumas_vars
!======================================================================================================================
 use mpas_kind_types
 use mpas_log

 implicit none
 public
 save


!--- local logicals and constants that point to logicals and constants originally defined in Registry_physics_cam6.xml.
!    these logicals and constants can be modified in namelist.atmosphere:
 character(len=StrKIND),pointer:: &
    pumas_precip_method

 logical,pointer:: &
    pumas_do_graupel,        &
    pumas_do_hail,           &
    pumas_uniform,           &
    pumas_do_cldice,         &
    pumas_hetfrz_classnuc,   &
    pumas_rm_supersat,       &
    pumas_do_sb_physics

logical,pointer:: &
    pumas_evap_sed_off,      &
    pumas_icenuc_rh_off,     &
    pumas_icenuc_use_meyers, &
    pumas_evap_scl_ifs,      &
    pumas_evap_rhthrsh_ifs,  &
    pumas_rainfreeze_ifs,    &
    pumas_ifs_sed,           &
    pumas_precip_fall_corr

 logical,pointer:: &
    pumas_nccons,            &
    pumas_nicons,            &
    pumas_ngcons,            &
    pumas_nrcons,            &
    pumas_nscons

 real(kind=R8KIND):: &
    pumas_ncnst,             &
    pumas_ninst,             &
    pumas_ngnst,             &
    pumas_nrnst,             &
    pumas_nsnst

!real(kind=R8KIND),pointer:: &
 real(kind=R8KIND):: &
    pumas_berg_eff_factor,   &
    pumas_dcs


!--- physical constants needed in PUMAS and initialized as functions of the physics constants used in MPAS physics,
!    except for rhmini:
 integer:: &
    mgncol,                  &
    nlev

 real(kind=R8KIND),parameter:: &
    rhmini = 10._R8KIND

 real(kind=R8KIND):: &
    cpair,                   &
    gravit,                  &
    latice,                  &
    latvap,                  &
    rair,                    &
    rh2o,                    &
    tmelt


!---
 type,public:: mp_pumas
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
       procedure:: allocate_mem   => pumas_allocate
       procedure:: deallocate_mem => pumas_deallocate
 end type


 contains


!======================================================================================================================
 subroutine pumas_allocate(this)
!======================================================================================================================

!inout arguments:
 class(mp_pumas),intent(inout):: this

!local variables and pointers:
 integer,pointer:: mgncol,nlev

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mp_pumas_allocate:')


!--- inputs:
 if(.not.associated(this%q)  ) allocate(this%q(mgncol,nlev))
 if(.not.associated(this%qcn)) allocate(this%qcn(mgncol,nlev))
 if(.not.associated(this%qrn)) allocate(this%qrn(mgncol,nlev))
 if(.not.associated(this%qin)) allocate(this%qin(mgncol,nlev))
 if(.not.associated(this%qsn)) allocate(this%qsn(mgncol,nlev))
 if(.not.associated(this%qgn)) allocate(this%qgn(mgncol,nlev))
 if(.not.associated(this%ncn)) allocate(this%ncn(mgncol,nlev))
 if(.not.associated(this%nrn)) allocate(this%nrn(mgncol,nlev))
 if(.not.associated(this%nin)) allocate(this%nin(mgncol,nlev))
 if(.not.associated(this%nsn)) allocate(this%nsn(mgncol,nlev))
 if(.not.associated(this%ngn)) allocate(this%ngn(mgncol,nlev))

 if(.not.associated(this%relvar)     ) allocate(this%relvar(mgncol,nlev)     )
 if(.not.associated(this%cldn)       ) allocate(this%cldn(mgncol,nlev)       )
 if(.not.associated(this%liqcldf)    ) allocate(this%liqcldf(mgncol,nlev)    )
 if(.not.associated(this%icecldf)    ) allocate(this%icecldf(mgncol,nlev)    )
 if(.not.associated(this%qsatfac)    ) allocate(this%qsatfac(mgncol,nlev)    )
 if(.not.associated(this%accre_enhan)) allocate(this%accre_enhan(mgncol,nlev))
 if(.not.associated(this%naai)       ) allocate(this%naai(mgncol,nlev)       )
 if(.not.associated(this%npccn)      ) allocate(this%npccn(mgncol,nlev)      )

 if(.not.associated(this%rndst)      ) allocate(this%rndst(mgncol,nlev,1)    )
 if(.not.associated(this%nacon)      ) allocate(this%nacon(mgncol,nlev,1)    )

!--- outputs:
 if(.not.associated(this%tlat)      ) allocate(this%tlat(mgncol,nlev)      )
 if(.not.associated(this%qvlat)     ) allocate(this%qvlat(mgncol,nlev)     )

 if(.not.associated(this%effc)      ) allocate(this%effc(mgncol,nlev)      )
 if(.not.associated(this%effc_fn)   ) allocate(this%effc_fn(mgncol,nlev)   )
 if(.not.associated(this%effi)      ) allocate(this%effi(mgncol,nlev)      )
 if(.not.associated(this%sadice)    ) allocate(this%sadice(mgncol,nlev)    )
 if(.not.associated(this%sadsnow)   ) allocate(this%sadsnow(mgncol,nlev)   )
 if(.not.associated(this%prect)     ) allocate(this%prect(mgncol,nlev)     )
 if(.not.associated(this%preci)     ) allocate(this%preci(mgncol,nlev)     )
 if(.not.associated(this%nevapr)    ) allocate(this%nevapr(mgncol,nlev)    )
 if(.not.associated(this%evapsnow)  ) allocate(this%evapsnow(mgncol,nlev)  )
 if(.not.associated(this%am_evp_st) ) allocate(this%am_evp_st(mgncol,nlev) )
 if(.not.associated(this%prain)     ) allocate(this%prain(mgncol,nlev)     )
 if(.not.associated(this%prodsnow)  ) allocate(this%prodsnow(mgncol,nlev)  )
 if(.not.associated(this%cmeout)    ) allocate(this%cmeout(mgncol,nlev)    )
 if(.not.associated(this%deffi)     ) allocate(this%deffi(mgncol,nlev)     )
 if(.not.associated(this%pgamrad)   ) allocate(this%pgamrad(mgncol,nlev)   )
 if(.not.associated(this%lamcrad)   ) allocate(this%lamcrad(mgncol,nlev)   )

 if(.not.associated(this%qctend)    ) allocate(this%qctend(mgncol,nlev)    )
 if(.not.associated(this%qrtend)    ) allocate(this%qrtend(mgncol,nlev)    )
 if(.not.associated(this%qitend)    ) allocate(this%qitend(mgncol,nlev)    )
 if(.not.associated(this%qstend)    ) allocate(this%qstend(mgncol,nlev)    )
 if(.not.associated(this%qgtend)    ) allocate(this%qgtend(mgncol,nlev)    )
 if(.not.associated(this%nctend)    ) allocate(this%nctend(mgncol,nlev)    )
 if(.not.associated(this%nrtend)    ) allocate(this%nrtend(mgncol,nlev)    )
 if(.not.associated(this%nitend)    ) allocate(this%nitend(mgncol,nlev)    )
 if(.not.associated(this%nstend)    ) allocate(this%nstend(mgncol,nlev)    )
 if(.not.associated(this%ngtend)    ) allocate(this%ngtend(mgncol,nlev)    )

 if(.not.associated(this%qcsedten)  ) allocate(this%qcsedten(mgncol,nlev)  )
 if(.not.associated(this%qrsedten)  ) allocate(this%qrsedten(mgncol,nlev)  )
 if(.not.associated(this%qisedten)  ) allocate(this%qisedten(mgncol,nlev)  )
 if(.not.associated(this%qssedten)  ) allocate(this%qssedten(mgncol,nlev)  )
 if(.not.associated(this%qgsedten)  ) allocate(this%qgsedten(mgncol,nlev)  )

 if(.not.associated(this%qrout)     ) allocate(this%qrout(mgncol,nlev)     )
 if(.not.associated(this%qsout)     ) allocate(this%qsout(mgncol,nlev)     )
 if(.not.associated(this%qgout)     ) allocate(this%qgout(mgncol,nlev)     )
 if(.not.associated(this%nrout)     ) allocate(this%nrout(mgncol,nlev)     )
 if(.not.associated(this%nsout)     ) allocate(this%nsout(mgncol,nlev)     )
 if(.not.associated(this%ngout)     ) allocate(this%ngout(mgncol,nlev)     )
 if(.not.associated(this%drout)     ) allocate(this%drout(mgncol,nlev)     )
 if(.not.associated(this%dsout)     ) allocate(this%dsout(mgncol,nlev)     )
 if(.not.associated(this%dgout)     ) allocate(this%dgout(mgncol,nlev)     )
 if(.not.associated(this%qrout2)    ) allocate(this%qrout2(mgncol,nlev)    )
 if(.not.associated(this%qsout2)    ) allocate(this%qsout2(mgncol,nlev)    )
 if(.not.associated(this%qgout2)    ) allocate(this%qgout2(mgncol,nlev)    )
 if(.not.associated(this%nrout2)    ) allocate(this%nrout2(mgncol,nlev)    )
 if(.not.associated(this%nsout2)    ) allocate(this%nsout2(mgncol,nlev)    )
 if(.not.associated(this%ngout2)    ) allocate(this%ngout2(mgncol,nlev)    )
 if(.not.associated(this%drout2)    ) allocate(this%drout2(mgncol,nlev)    )
 if(.not.associated(this%dsout2)    ) allocate(this%dsout2(mgncol,nlev)    )
 if(.not.associated(this%dgout2)    ) allocate(this%dgout2(mgncol,nlev)    )

 if(.not.associated(this%pratot)    ) allocate(this%pratot(mgncol,nlev)    )
 if(.not.associated(this%prctot)    ) allocate(this%prctot(mgncol,nlev)    )
 if(.not.associated(this%mnuccctot) ) allocate(this%mnuccctot(mgncol,nlev) )
 if(.not.associated(this%mnuccctot) ) allocate(this%mnuccctot(mgncol,nlev) )
 if(.not.associated(this%msacwitot) ) allocate(this%msacwitot(mgncol,nlev) )
 if(.not.associated(this%psacwstot) ) allocate(this%psacwstot(mgncol,nlev) )
 if(.not.associated(this%bergstot)  ) allocate(this%bergstot(mgncol,nlev)  )
 if(.not.associated(this%bergstot)  ) allocate(this%bergstot(mgncol,nlev)  )
 if(.not.associated(this%melttot)   ) allocate(this%melttot(mgncol,nlev)   )
 if(.not.associated(this%meltstot)  ) allocate(this%meltstot(mgncol,nlev)  )
 if(.not.associated(this%meltgtot)  ) allocate(this%meltgtot(mgncol,nlev)  )
 if(.not.associated(this%mnudeptot) ) allocate(this%mnudeptot(mgncol,nlev) )
 if(.not.associated(this%homotot)   ) allocate(this%homotot(mgncol,nlev)   )
 if(.not.associated(this%qcrestot)  ) allocate(this%qcrestot(mgncol,nlev)  )
 if(.not.associated(this%prcitot)   ) allocate(this%prcitot(mgncol,nlev)   )
 if(.not.associated(this%praitot)   ) allocate(this%praitot(mgncol,nlev)   )
 if(.not.associated(this%qirestot)  ) allocate(this%qirestot(mgncol,nlev)  )
 if(.not.associated(this%mnuccrtot )) allocate(this%mnuccrtot (mgncol,nlev))
 if(.not.associated(this%mnuccritot)) allocate(this%mnuccritot(mgncol,nlev))
 if(.not.associated(this%pracstot)  ) allocate(this%pracstot(mgncol,nlev)  )
 if(.not.associated(this%meltsdttot)) allocate(this%meltsdttot(mgncol,nlev))
 if(.not.associated(this%frzrdttot) ) allocate(this%frzrdttot(mgncol,nlev) )
 if(.not.associated(this%mnuccdtot) ) allocate(this%mnuccdtot(mgncol,nlev) )
 if(.not.associated(this%pracgtot)  ) allocate(this%pracgtot(mgncol,nlev)  )
 if(.not.associated(this%psacwgtot) ) allocate(this%psacwgtot(mgncol,nlev) )
 if(.not.associated(this%pgsacwtot) ) allocate(this%pgsacwtot(mgncol,nlev) )
 if(.not.associated(this%pgracstot) ) allocate(this%pgracstot(mgncol,nlev) )
 if(.not.associated(this%prdgtot)   ) allocate(this%prdgtot(mgncol,nlev)   )
 if(.not.associated(this%qmultgtot) ) allocate(this%qmultgtot(mgncol,nlev) )
 if(.not.associated(this%qmultrgtot)) allocate(this%qmultrgtot(mgncol,nlev))
 if(.not.associated(this%psacrtot)  ) allocate(this%psacrtot(mgncol,nlev)  )
 if(.not.associated(this%npracgtot) ) allocate(this%npracgtot(mgncol,nlev) )
 if(.not.associated(this%nscngtot)  ) allocate(this%nscngtot(mgncol,nlev)  )
 if(.not.associated(this%ngracstot) ) allocate(this%ngracstot(mgncol,nlev) )
 if(.not.associated(this%nmultgtot) ) allocate(this%nmultgtot(mgncol,nlev) )
 if(.not.associated(this%nmultrgtot)) allocate(this%nmultrgtot(mgncol,nlev))
 if(.not.associated(this%npsacwgtot)) allocate(this%npsacwgtot(mgncol,nlev))
 if(.not.associated(this%refl)      ) allocate(this%refl(mgncol,nlev)      )
 if(.not.associated(this%arefl)     ) allocate(this%arefl(mgncol,nlev)     )
 if(.not.associated(this%areflz)    ) allocate(this%areflz(mgncol,nlev)    )
 if(.not.associated(this%frefl)     ) allocate(this%frefl(mgncol,nlev)     )
 if(.not.associated(this%csrfl)     ) allocate(this%csrfl(mgncol,nlev)     )
 if(.not.associated(this%acsrfl)    ) allocate(this%acsrfl(mgncol,nlev)    )
 if(.not.associated(this%fcsrfl)    ) allocate(this%fcsrfl(mgncol,nlev)    )
 if(.not.associated(this%rercld)    ) allocate(this%rercld(mgncol,nlev)    )
 if(.not.associated(this%ncai)      ) allocate(this%ncai(mgncol,nlev)      )
 if(.not.associated(this%ncal)      ) allocate(this%ncal(mgncol,nlev)      )
 if(.not.associated(this%nfice)     ) allocate(this%nfice(mgncol,nlev)     )
 if(.not.associated(this%qcrat)     ) allocate(this%qcrat(mgncol,nlev)     )
 if(.not.associated(this%prer_evap) ) allocate(this%prer_evap(mgncol,nlev) )

 if(.not.associated(this%freqr)     ) allocate(this%freqr(mgncol,nlev)     )
 if(.not.associated(this%freqs)     ) allocate(this%freqs(mgncol,nlev)     )
 if(.not.associated(this%freqg)     ) allocate(this%freqg(mgncol,nlev)     )
 if(.not.associated(this%reff_rain) ) allocate(this%reff_rain(mgncol,nlev) )
 if(.not.associated(this%reff_snow) ) allocate(this%reff_snow(mgncol,nlev) )
 if(.not.associated(this%reff_grau) ) allocate(this%reff_grau(mgncol,nlev) )
 if(.not.associated(this%umr)       ) allocate(this%umr(mgncol,nlev)       )
 if(.not.associated(this%ums)       ) allocate(this%ums(mgncol,nlev)       )
 if(.not.associated(this%umg)       ) allocate(this%umg(mgncol,nlev)       )
 if(.not.associated(this%qcsevap)   ) allocate(this%qcsevap(mgncol,nlev)   )
 if(.not.associated(this%qisevap)   ) allocate(this%qisevap(mgncol,nlev)   )
 if(.not.associated(this%qvres)     ) allocate(this%qvres(mgncol,nlev)     )
 if(.not.associated(this%cmeitot)   ) allocate(this%cmeitot(mgncol,nlev)   )
 if(.not.associated(this%vtrmc)     ) allocate(this%vtrmc(mgncol,nlev)     )
 if(.not.associated(this%vtrmi)     ) allocate(this%vtrmi(mgncol,nlev)     )

 if(.not.associated(this%qcsinksum_rate1ord)) allocate(this%qcsinksum_rate1ord(mgncol,nlev))

 if(.not.associated(this%lflx)      ) allocate(this%lflx(mgncol,nlev+1)    )
 if(.not.associated(this%rflx)      ) allocate(this%rflx(mgncol,nlev+1)    )
 if(.not.associated(this%iflx)      ) allocate(this%iflx(mgncol,nlev+1)    )
 if(.not.associated(this%sflx)      ) allocate(this%sflx(mgncol,nlev+1)    )
 if(.not.associated(this%gflx)      ) allocate(this%gflx(mgncol,nlev+1)    )

 call mpas_log_write('--- end subroutine mp_pumas_allocate:')

 end subroutine pumas_allocate

!======================================================================================================================
 subroutine pumas_deallocate(this)
!======================================================================================================================

!inout arguments:
 class(mp_pumas),intent(inout):: this

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mp_pumas_deallocate:')

!--- inputs:
 if(associated(this%q)      ) deallocate(this%q  )
 if(associated(this%qcn)    ) deallocate(this%qcn)
 if(associated(this%qrn)    ) deallocate(this%qrn)
 if(associated(this%qin)    ) deallocate(this%qin)
 if(associated(this%qsn)    ) deallocate(this%qsn)
 if(associated(this%qgn)    ) deallocate(this%qgn)
 if(associated(this%ncn)    ) deallocate(this%ncn)
 if(associated(this%nrn)    ) deallocate(this%nrn)
 if(associated(this%nin)    ) deallocate(this%nin)
 if(associated(this%nsn)    ) deallocate(this%nsn)
 if(associated(this%ngn)    ) deallocate(this%ngn)

 if(associated(this%relvar)     ) deallocate(this%relvar     )
 if(associated(this%cldn)       ) deallocate(this%cldn       )
 if(associated(this%liqcldf)    ) deallocate(this%liqcldf    )
 if(associated(this%icecldf)    ) deallocate(this%icecldf    )
 if(associated(this%qsatfac)    ) deallocate(this%qsatfac    )
 if(associated(this%accre_enhan)) deallocate(this%accre_enhan)
 if(associated(this%naai)       ) deallocate(this%naai       )
 if(associated(this%npccn)      ) deallocate(this%npccn      )

 if(associated(this%rndst)      ) deallocate(this%rndst      )
 if(associated(this%nacon)      ) deallocate(this%nacon      )


!--- outputs:
 if(associated(this%tlat)      ) deallocate(this%tlat      )
 if(associated(this%qvlat)     ) deallocate(this%qvlat     )

 if(associated(this%effc)      ) deallocate(this%effc      )
 if(associated(this%effc_fn)   ) deallocate(this%effc_fn   )
 if(associated(this%effi)      ) deallocate(this%effi      )
 if(associated(this%sadice)    ) deallocate(this%sadice    )
 if(associated(this%sadsnow)   ) deallocate(this%sadsnow   )
 if(associated(this%prect)     ) deallocate(this%prect     )
 if(associated(this%preci)     ) deallocate(this%preci     )
 if(associated(this%nevapr)    ) deallocate(this%nevapr    )
 if(associated(this%evapsnow)  ) deallocate(this%evapsnow  )
 if(associated(this%am_evp_st) ) deallocate(this%am_evp_st )
 if(associated(this%prain)     ) deallocate(this%prain     )
 if(associated(this%prodsnow)  ) deallocate(this%prodsnow  )
 if(associated(this%cmeout)    ) deallocate(this%cmeout    )
 if(associated(this%deffi)     ) deallocate(this%deffi     )
 if(associated(this%pgamrad)   ) deallocate(this%pgamrad   )
 if(associated(this%lamcrad)   ) deallocate(this%lamcrad   )

 if(associated(this%qctend)    ) deallocate(this%qctend    )
 if(associated(this%qrtend)    ) deallocate(this%qrtend    )
 if(associated(this%qitend)    ) deallocate(this%qitend    )
 if(associated(this%qstend)    ) deallocate(this%qstend    )
 if(associated(this%qgtend)    ) deallocate(this%qgtend    )
 if(associated(this%nctend)    ) deallocate(this%nctend    )
 if(associated(this%nrtend)    ) deallocate(this%nrtend    )
 if(associated(this%nitend)    ) deallocate(this%nitend    )
 if(associated(this%nstend)    ) deallocate(this%nstend    )
 if(associated(this%ngtend)    ) deallocate(this%ngtend    )

 if(associated(this%qcsedten)  ) deallocate(this%qcsedten  )
 if(associated(this%qrsedten)  ) deallocate(this%qrsedten  )
 if(associated(this%qisedten)  ) deallocate(this%qisedten  )
 if(associated(this%qssedten)  ) deallocate(this%qssedten  )
 if(associated(this%qgsedten)  ) deallocate(this%qgsedten  )

 if(associated(this%qrout )    ) deallocate(this%qrout     )
 if(associated(this%qsout )    ) deallocate(this%qsout     )
 if(associated(this%qgout )    ) deallocate(this%qgout     )
 if(associated(this%nrout )    ) deallocate(this%nrout     )
 if(associated(this%nsout )    ) deallocate(this%nsout     )
 if(associated(this%ngout )    ) deallocate(this%ngout     )
 if(associated(this%drout )    ) deallocate(this%drout     )
 if(associated(this%dsout )    ) deallocate(this%dsout     )
 if(associated(this%dgout )    ) deallocate(this%dgout     )
 if(associated(this%qrout2)    ) deallocate(this%qrout2    )
 if(associated(this%qsout2)    ) deallocate(this%qsout2    )
 if(associated(this%qgout2)    ) deallocate(this%qgout2    )
 if(associated(this%nrout2)    ) deallocate(this%nrout2    )
 if(associated(this%nsout2)    ) deallocate(this%nsout2    )
 if(associated(this%ngout2)    ) deallocate(this%ngout2    )
 if(associated(this%drout2)    ) deallocate(this%drout2    )
 if(associated(this%dsout2)    ) deallocate(this%dsout2    )
 if(associated(this%dgout2)    ) deallocate(this%dgout2    )

 if(associated(this%pratot)    ) deallocate(this%pratot    )
 if(associated(this%prctot)    ) deallocate(this%prctot    )
 if(associated(this%mnuccctot) ) deallocate(this%mnuccctot )
 if(associated(this%mnuccctot) ) deallocate(this%mnuccctot )
 if(associated(this%msacwitot) ) deallocate(this%msacwitot )
 if(associated(this%psacwstot) ) deallocate(this%psacwstot )
 if(associated(this%bergstot)  ) deallocate(this%bergstot  )
 if(associated(this%bergstot)  ) deallocate(this%bergstot  )
 if(associated(this%melttot)   ) deallocate(this%melttot   )
 if(associated(this%meltstot)  ) deallocate(this%meltstot  )
 if(associated(this%meltgtot)  ) deallocate(this%meltgtot  )
 if(associated(this%mnudeptot) ) deallocate(this%mnudeptot )
 if(associated(this%homotot)   ) deallocate(this%homotot   )
 if(associated(this%qcrestot)  ) deallocate(this%qcrestot  )
 if(associated(this%prcitot)   ) deallocate(this%prcitot   )
 if(associated(this%praitot)   ) deallocate(this%praitot   )
 if(associated(this%qirestot)  ) deallocate(this%qirestot  )
 if(associated(this%mnuccrtot )) deallocate(this%mnuccrtot )
 if(associated(this%mnuccritot)) deallocate(this%mnuccritot)
 if(associated(this%pracstot)  ) deallocate(this%pracstot  )
 if(associated(this%meltsdttot)) deallocate(this%meltsdttot)
 if(associated(this%frzrdttot) ) deallocate(this%frzrdttot )
 if(associated(this%mnuccdtot) ) deallocate(this%mnuccdtot )
 if(associated(this%pracgtot)  ) deallocate(this%pracgtot  )
 if(associated(this%psacwgtot) ) deallocate(this%psacwgtot )
 if(associated(this%pgsacwtot) ) deallocate(this%pgsacwtot )
 if(associated(this%pgracstot) ) deallocate(this%pgracstot )
 if(associated(this%prdgtot)   ) deallocate(this%prdgtot   )
 if(associated(this%qmultgtot) ) deallocate(this%qmultgtot )
 if(associated(this%qmultrgtot)) deallocate(this%qmultrgtot)
 if(associated(this%psacrtot)  ) deallocate(this%psacrtot  )
 if(associated(this%npracgtot) ) deallocate(this%npracgtot )
 if(associated(this%nscngtot)  ) deallocate(this%nscngtot  )
 if(associated(this%ngracstot) ) deallocate(this%ngracstot )
 if(associated(this%nmultgtot) ) deallocate(this%nmultgtot )
 if(associated(this%nmultrgtot)) deallocate(this%nmultrgtot)
 if(associated(this%npsacwgtot)) deallocate(this%npsacwgtot)
 if(associated(this%refl)      ) deallocate(this%refl      )
 if(associated(this%arefl)     ) deallocate(this%arefl     )
 if(associated(this%areflz)    ) deallocate(this%areflz    )
 if(associated(this%frefl)     ) deallocate(this%frefl     )
 if(associated(this%csrfl)     ) deallocate(this%csrfl     )
 if(associated(this%acsrfl)    ) deallocate(this%acsrfl    )
 if(associated(this%fcsrfl)    ) deallocate(this%fcsrfl    )
 if(associated(this%rercld)    ) deallocate(this%rercld    )
 if(associated(this%ncai)      ) deallocate(this%ncai      )
 if(associated(this%ncal)      ) deallocate(this%ncal      )
 if(associated(this%nfice)     ) deallocate(this%nfice     )
 if(associated(this%qcrat)     ) deallocate(this%qcrat     )
 if(associated(this%prer_evap) ) deallocate(this%prer_evap )

 if(associated(this%freqr)     ) deallocate(this%freqr     )
 if(associated(this%freqs)     ) deallocate(this%freqs     )
 if(associated(this%freqg)     ) deallocate(this%freqg     )
 if(associated(this%reff_rain) ) deallocate(this%reff_rain )
 if(associated(this%reff_snow) ) deallocate(this%reff_snow )
 if(associated(this%reff_grau) ) deallocate(this%reff_grau )
 if(associated(this%umr      ) ) deallocate(this%umr       )
 if(associated(this%ums      ) ) deallocate(this%ums       )
 if(associated(this%umg      ) ) deallocate(this%umg       )
 if(associated(this%qcsevap)   ) deallocate(this%qcsevap   )
 if(associated(this%qisevap)   ) deallocate(this%qisevap   )
 if(associated(this%qvres)     ) deallocate(this%qvres     )
 if(associated(this%cmeitot)   ) deallocate(this%cmeitot   )
 if(associated(this%vtrmc)     ) deallocate(this%vtrmc     )
 if(associated(this%vtrmi)     ) deallocate(this%vtrmi     )

 if(associated(this%qcsinksum_rate1ord)) deallocate(this%qcsinksum_rate1ord)

 if(associated(this%lflx)      ) deallocate(this%lflx      )
 if(associated(this%rflx)      ) deallocate(this%rflx      )
 if(associated(this%iflx)      ) deallocate(this%iflx      )
 if(associated(this%sflx)      ) deallocate(this%sflx      )
 if(associated(this%gflx)      ) deallocate(this%gflx      )

 call mpas_log_write('--- end subroutine mp_pumas_deallocate:')

 end subroutine pumas_deallocate

!======================================================================================================================
 end module pumas_vars
!======================================================================================================================

