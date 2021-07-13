!======================================================================================================================
 module mp_pumas_mpas
!======================================================================================================================

 use shr_kind_mod, only: &
    r8 => shr_kind_r8

 use shr_const_mod, only: &
    gravit => shr_const_g      , &
    rair   => shr_const_rdair  , &
    rh2o   => shr_const_rwv    , &
    cpair  => shr_const_cpdair , &
    tmelt  => shr_const_tkfrz  , &
    latvap => shr_const_latvap , &
    latice => shr_const_latice

 use micro_mg3_0, only: &
    micro_mg_init


 implicit none
 private
 public:: mp_pumas_mpas_init


!--- variables that will need to be moved to Registry.xml:
!MG3 dense precipitating ice. Note, only 1 can be true, or both false.
 logical,parameter:: micro_mg_do_graupel = .false. ! .true.  = configure with graupel
                                                   ! .false. = no graupel (hail possible)
 logical,parameter:: micro_mg_do_hail    = .false. ! .true.  = configure with hail
                                                   ! .false. = no hail (graupel possible)
 logical,parameter:: microp_uniform      = .false. ! .true.  = configure uniform for sub-columns
                                                   ! .false. = use w/o sub-columns (standard)
 logical,parameter:: do_cldice           = .true.  ! .true   = do all processes (standard)
                                                   ! .false. = skip all processes affecting cloud ice
 logical,parameter:: use_hetfrz_classnuc = .true.  ! .true.  = use heterogeneous freezing
                                                   ! .false. = skip heterogeneous freezing
 logical,parameter:: remove_supersat     = .true.  ! .true.  = remove supersaturation after sedimentation loop
                                                   ! .false. = skip removing supersaturation after sedimentation loop
 logical,parameter:: do_sb_physics       = .false. ! .true.  = do SB autoconversion and accretion physics
                                                   ! .false. = default autoconversion and accretion physics

 character(len=16),parameter:: micro_mg_precip_frac_method = 'max_overlap' ! type of precipitation fraction method
 real(r8),parameter:: micro_mg_berg_eff_factor             = 1.0_r8        ! berg efficiency factor

!IFS-like switches:
 logical,parameter:: micro_mg_evap_sed_off      = .false. ! .true. = Turn off evaporation/sublimation based on
                                                          !          cloud fraction for sedimenting condensate
 logical,parameter:: micro_mg_icenuc_rh_off     = .false. ! .true. = remove RH conditional from ice nucleation
 logical,parameter:: micro_mg_icenuc_use_meyers = .false. ! .true. = use Meyers Ice Nucleation
 logical,parameter:: micro_mg_evap_scl_ifs      = .false. ! .true. = scale evaporation as IFS does (*0.3)
 logical,parameter:: micro_mg_evap_rhthrsh_ifs  = .false. ! .true. = evap RH threhold following ifs
 logical,parameter:: micro_mg_rainfreeze_ifs    = .false. ! .true. = rain freezing temp following ifs
 logical,parameter:: micro_mg_ifs_sed           = .false. ! .true. = snow sedimentation = 1m/s following ifs
 logical,parameter:: micro_mg_precip_fall_corr  = .false. ! .true. = ensure rain fall speed non-zero if rain above


 real(r8):: micro_mg_dcs = -1._r8

!namelist variables for option to specify constant cloud droplet/ice number
 logical :: micro_mg_nccons = .false.                     ! set .true. to specify constant cloud droplet number
 logical :: micro_mg_nicons = .false.                     ! set .true. to specify constant cloud ice number
 logical :: micro_mg_ngcons = .false.                     ! set .true. to specify constant graupel/hail number
 logical :: micro_mg_nrcons = .false.                     ! set .true. to specify constant rain number
 logical :: micro_mg_nscons = .false.                     ! set .true. to specify constant snow number

!parameters for specified ice and droplet number concentration
!note: these are local in-cloud values, not grid-mean
 real(r8) :: micro_mg_ncnst = 50.e6_r8                    ! constant liquid droplet num concentration (m-3)
 real(r8) :: micro_mg_ninst = 0.05e6_r8                   ! ice num concentration when nicons=.true. (m-3)
 real(r8) :: micro_mg_nrnst = 0.2e6_r8                    ! rain  num concentration when nrcons=.true. (m-3)
 real(r8) :: micro_mg_nsnst = 0.005e6_r8                  ! snow num concentration when nscons=.true. (m-3)
 real(r8) :: micro_mg_ngnst = 0.0005e6_r8                 ! graupel/hail num concentration when ngcons=.true. (m-3)

!temporary:
 real(r8):: rhmini = 10.


 contains


!======================================================================================================================
 subroutine mp_pumas_mpas_init()
!======================================================================================================================
 
 character(128):: errstring

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mp_pumas_init:')

 errstring = ' '

 call micro_mg_init( &
            kind      = r8        , gravit    = gravit    , rair      = rair   , &
            rh2o      = rh2o      , cpair     = cpair     , tmelt_in  = tmelt  , &
            latvap    = latvap    , latice    = latice    , rhmini_in = rhmini , &
            nccons_in = micro_mg_nccons    , nicons_in = micro_mg_nicons       , &
            ncnst_in  = micro_mg_ncnst     , ninst_in  = micro_mg_ninst        , &
            ngcons_in = micro_mg_ngcons    , ngnst_in  = micro_mg_ngnst        , &
            nrcons_in = micro_mg_nrcons    , nrnst_in  = micro_mg_nrnst        , &
            nscons_in = micro_mg_nscons    , nsnst_in  = micro_mg_nsnst        , &
            micro_mg_dcs                   = micro_mg_dcs                      , &
            micro_mg_do_hail_in            = micro_mg_do_hail                  , &
            micro_mg_do_graupel_in         = micro_mg_do_graupel               , &
            microp_uniform_in              = microp_uniform                    , &
            do_cldice_in                   = do_cldice                         , &
            use_hetfrz_classnuc_in         = use_hetfrz_classnuc               , &
            micro_mg_precip_frac_method_in = micro_mg_precip_frac_method       , &
            micro_mg_berg_eff_factor_in    = micro_mg_berg_eff_factor          , &
            remove_supersat_in             = remove_supersat                   , &
            do_sb_physics_in               = do_sb_physics                     , &
            micro_mg_evap_sed_off_in       = micro_mg_evap_sed_off             , &
            micro_mg_icenuc_rh_off_in      = micro_mg_icenuc_rh_off            , &
            micro_mg_icenuc_use_meyers_in  = micro_mg_icenuc_use_meyers        , &
            micro_mg_evap_scl_ifs_in       = micro_mg_evap_scl_ifs             , &
            micro_mg_evap_rhthrsh_ifs_in   = micro_mg_evap_rhthrsh_ifs         , &
            micro_mg_rainfreeze_ifs_in     = micro_mg_rainfreeze_ifs           , &
            micro_mg_ifs_sed_in            = micro_mg_ifs_sed                  , &
            micro_mg_precip_fall_corr      = micro_mg_ifs_sed                  , &
            errstring                      = errstring                           &
                   )

 call mpas_log_write('--- end subroutine mp_pumas_init:')

 end subroutine mp_pumas_mpas_init

!======================================================================================================================
 subroutine mp_pumas_mpas_run()
!======================================================================================================================

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mp_pumas_run:')

 end subroutine mp_pumas_mpas_run

!======================================================================================================================
 subroutine mp_pumas_mpas_final()
!======================================================================================================================


 end subroutine mp_pumas_mpas_final

!======================================================================================================================
 end module mp_pumas_mpas
!======================================================================================================================
