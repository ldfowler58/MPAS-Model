!======================================================================================================================
 module mpas_atmphys_driver_microphysics_pumas
!======================================================================================================================
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_constants,only: cp,gravity,R_d,R_v,svpt0,xlf,xlv
 use mpas_atmphys_init_fromregistry,only: pumas_init_fromRegistry

 use micro_mg3_0
 use pumas_mpas
 use pumas_vars

 implicit none
 private
 public:: init_microphysics_pumas


 contains


!======================================================================================================================
 subroutine init_microphysics_pumas(configs,mesh)
!======================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 real(kind=R8KIND):: &
    pumas_ncnst_sgl,pumas_ninst_sgl,pumas_ngnst_dbl,pumas_nrnst_dbl,pumas_nsnst_dbl

!local variables and arrays:
 character(len=StrKIND):: errmsg
 integer:: errflg

 integer,pointer:: nCells,nVertLevels

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_microphysics_pumas:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)


!--- initializes the local variables (logicals,constants) needed to run the PUMAS cloud microphysics scheme. these
!    variables are initially available from Registry_physics_cam6.xml and their default values can be modified in
!    namelist.atmosphere. in a CCPP framework, this call should be moved up in the calling sequence.
 call pumas_init_fromRegistry(configs)


!--- initializes pumas constants:
 call pumas_mpas_init(nCells,nVertLevels,cp,gravity,R_d,R_v,svpt0,xlf,xlv,errmsg,errflg)


!--- initializes pumas cloud microphysics:
 call mpas_log_write('--- enter subroutine micro_mg_init:')
 call micro_mg_init( &
            kind   = R8KIND , gravit = gravit , rair      = rair   , &
            rh2o   = rh2o   , cpair  = cpair  , tmelt_in  = tmelt  , &
            latvap = latvap , latice = latice , rhmini_in = rhmini , &
            nccons_in = pumas_nccons       , ncnst_in  = pumas_ncnst  , & 
            nicons_in = pumas_nicons       , ninst_in  = pumas_ninst  , &
            ngcons_in = pumas_ngcons       , ngnst_in  = pumas_ngnst  , &
            nrcons_in = pumas_nrcons       , nrnst_in  = pumas_nrnst  , &
            nscons_in = pumas_nscons       , nsnst_in  = pumas_nsnst  , &
            micro_mg_dcs                   = pumas_dcs                , &
            micro_mg_do_hail_in            = pumas_do_hail            , &
            micro_mg_do_graupel_in         = pumas_do_graupel         , &
            microp_uniform_in              = pumas_uniform            , &
            do_cldice_in                   = pumas_do_cldice          , &
            use_hetfrz_classnuc_in         = pumas_hetfrz_classnuc    , &
            micro_mg_precip_frac_method_in = pumas_precip_method      , &
            micro_mg_berg_eff_factor_in    = pumas_berg_eff_factor    , &
            remove_supersat_in             = pumas_rm_supersat        , &
            do_sb_physics_in               = pumas_do_sb_physics      , &
            micro_mg_evap_sed_off_in       = pumas_evap_sed_off       , &
            micro_mg_icenuc_rh_off_in      = pumas_icenuc_rh_off      , &
            micro_mg_icenuc_use_meyers_in  = pumas_icenuc_use_meyers  , &
            micro_mg_evap_scl_ifs_in       = pumas_evap_scl_ifs       , &
            micro_mg_evap_rhthrsh_ifs_in   = pumas_evap_rhthrsh_ifs   , &
            micro_mg_rainfreeze_ifs_in     = pumas_rainfreeze_ifs     , &
            micro_mg_ifs_sed_in            = pumas_ifs_sed            , &
            micro_mg_precip_fall_corr      = pumas_precip_fall_corr   , &
            errstring                      = errmsg                     &
                   )
 call mpas_log_write('--- end subroutine micro_mg_init:')

 call mpas_log_write('--- end subroutine init_microphysics_pumas')

 end subroutine init_microphysics_pumas

!======================================================================================================================
 subroutine driver_microphysics_pumas()
!======================================================================================================================

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_microphysics_pumas:')

!call microphysics_pumas_from_MPAS()

!--- call to PUMAS parameterization of cloud microphysics:
!call micro_mg_tend( &
!              mgncol     = mgncol          , nlev        = nlev             , &
!              deltatin   = delta           , t           = temp             , &
!              p          = press           , pdel        = pthick           , &
!              q          = q               , qcn         = qcn              , &
!              qin        = qin             , qrn         = qrn              , &
!              qsn        = qsn             , qgr         = qgn              , &
!              ncn        = ncn             , nin         = nin              , &
!              nrn        = nrn             , nsn         = nsn              , &
!              ngr        = ngn             , relvar      = this%relvar      , &
!              cldn       = this%cldn       , liqcldf     = this%liqcldf     , &
!              icecldf    = this%icecldf    , qsatfac     = this%qsatfac     , &
!              tlat       = this%tlat       , qvlat       = this%qvlat       , &
!              qctend     = this%qctend     , qrtend      = this%qrtend      , &
!              qitend     = this%qitend     , qstend      = this%qstend      , &
!              qgtend     = this%qgtend     , nctend      = this%nctend      , &
!              nrtend     = this%nrtend     , nitend      = this%nitend      , &
!              nstend     = this%nstend     , ngtend      = this%ngtend      , &
!              qcsedten   = this%qcsedten   , qrsedten    = this%qrsedten    , &
!              qisedten   = this%qisedten   , qssedten    = this%qssedten    , &
!              qgsedten   = this%qgsedten   , qrout       = this%qrout       , &
!              qsout      = this%qsout      , qgout       = this%qgout       , &
!              nrout      = this%nrout      , nsout       = this%nsout       , &
!              ngout      = this%ngout      ,                                  &
!              dsout      = this%dsout      , dgout       = this%dgout       , &
!              qrout2     = this%qrout2     , qsout2      = this%qsout2      , &
!              qgout2     = this%qgout2     , nrout2      = this%nrout2      , &
!              nsout2     = this%nsout2     , ngout2      = this%ngout2      , &
!              drout2     = this%drout2     , dsout2      = this%dsout2      , &
!              dgout2     = this%dgout2     , freqr       = this%freqr       , &
!              freqs      = this%freqs      , freqg       = this%freqg       , &
!              reff_rain  = this%reff_rain  , reff_snow   = this%reff_snow   , &
!              reff_grau  = this%reff_grau  , umr         = this%umr         , &
!              ums        = this%ums        , umg         = this%umg         , &
!              pratot     = this%pratot     , prctot      = this%prctot      , &
!              mnuccctot  = this%mnuccctot  , mnuccttot   = this%mnuccttot   , &
!              msacwitot  = this%msacwitot  , psacwstot   = this%psacwstot   , &
!              bergstot   = this%bergstot   , bergtot     = this%bergtot     , &
!              melttot    = this%melttot    , meltstot    = this%meltstot    , &
!              meltgtot   = this%meltgtot   , mnudeptot   = this%mnudeptot   , &
!              homotot    = this%homotot    , qcrestot    = this%qcrestot    , &
!              prcitot    = this%prcitot    , praitot     = this%praitot     , &
!              qirestot   = this%qirestot   , mnuccrtot   = this%mnuccrtot   , &
!              mnuccritot = this%mnuccritot , pracstot    = this%pracstot    , &
!              meltsdttot = this%meltsdttot , frzrdttot   = this%frzrdttot   , &
!              mnuccdtot  = this%mnuccdtot  , pracgtot    = this%pracgtot    , &
!              psacwgtot  = this%psacwgtot  , pgsacwtot   = this%pgsacwtot   , &
!              pgracstot  = this%pgracstot  , prdgtot     = this%prdgtot     , &
!              qmultgtot  = this%qmultgtot  , qmultrgtot  = this%qmultrgtot  , &
!              psacrtot   = this%psacrtot   , npracgtot   = this%npracgtot   , &
!              nscngtot   = this%nscngtot   , ngracstot   = this%ngracstot   , &
!              nmultgtot  = this%nmultgtot  , nmultrgtot  = this%nmultrgtot  , &
!              npsacwgtot = this%npsacwgtot , refl        = this%refl        , &
!              arefl      = this%arefl      , areflz      = this%areflz      , &
!              frefl      = this%frefl      , csrfl       = this%csrfl       , &
!              acsrfl     = this%acsrfl     , fcsrfl      = this%fcsrfl      , &
!              rercld     = this%rercld     , ncai        = this%ncai        , &
!              ncal       = this%ncal       , nfice       = this%nfice       , &
!              qcrat      = this%qcrat      , effc        = this%effc        , &
!              effc_fn    = this%effc_fn    , effi        = this%effi        , &
!              sadice     = this%sadice     , sadsnow     = this%sadsnow     , &
!              prect      = this%prect      , preci       = this%preci       , &
!              nevapr     = this%nevapr     , evapsnow    = this%evapsnow    , &
!              am_evp_st  = this%am_evp_st  , prain       = this%prain       , &
!              prodsnow   = this%prodsnow   , cmeout      = this%cmeout      , &
!              deffi      = this%deffi      , pgamrad     = this%pgamrad     , &
!              lamcrad    = this%lamcrad    , accre_enhan = this%accre_enhan , &
!              lflx       = this%lflx       , rflx        = this%rflx        , &
!              iflx       = this%iflx       , sflx        = this%sflx        , &
!              gflx       = this%gflx       , naai        = this%naai        , &
!              npccn      = this%npccn      , nacon       = this%nacon       , &
!              rndst      = this%rndst      , qcsevap     = this%qcsevap     , &
!              qisevap    = this%qisevap    , qvres       = this%qvres       , &
!              cmeitot    = this%cmeitot    , vtrmc       = this%vtrmc       , &
!              vtrmi      = this%vtrmi      , tnd_qsnow   = this%tnd_qsnow   , &
!              tnd_nsnow  = this%tnd_nsnow  , re_ice      = this%re_ice      , &
!              frzimm     = this%frzimm     , frzcnt      = this%frzcnt      , &
!              frzdep     = this%frzdep     , prer_evap   = this%prer_evap   , &
!              qcsinksum_rate1ord = this%qcsinksum_rate1ord                  , &
!              errstring  = errstring                                          &
!                  )

!call microphysics_pumas_to_MPAS()

 call mpas_log_write('--- end subroutine driver_microphysics_pumas:')

 end subroutine driver_microphysics_pumas

!======================================================================================================================
 subroutine microphysics_pumas_from_MPAS(this,configs,mesh,state,time_lev,diag)
!======================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag
 integer,intent(in):: time_lev

!inout arguments:
 class(mp_pumas),intent(inout):: this


!local variables and arrays:
 character(len=128):: errstring
 integer:: mgncol,nlev
 integer:: i,k

!MPAS pointers:
 integer,pointer:: index_qv,index_qc,index_qr,index_qi,index_qs,index_qg
 integer,pointer:: index_nc,index_nr,index_ni,index_ns,index_ng

 real(kind=RKIND),pointer:: dt_dyn
 real(kind=RKIND),dimension(:,:),pointer:: exner,pressure_b,pressure_p,theta
 real(kind=RKIND),dimension(:,:),pointer:: qv,qc,qr,qi,qs,qg,nc,nr,ni,ns,ng
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

!PUMAS pointers:
 real(kind=R8KIND),pointer:: delta
 real(kind=R8KIND),dimension(:,:),pointer:: press,pthick,temp
 real(kind=R8KIND),dimension(:,:),pointer:: q,qcn,qrn,qin,qsn,qgn
 real(kind=R8KIND),dimension(:,:),pointer:: ncn,nrn,nin,nsn,ngn

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine driver_microphysics_pumas:')

 call mpas_pool_get_config(configs,'config_dt',dt_dyn)


!--- MPAS state and diagnostics variables:
 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_dimension(state,'index_qc',index_qc)
 call mpas_pool_get_dimension(state,'index_qr',index_qr)
 call mpas_pool_get_dimension(state,'index_qi',index_qi)
 call mpas_pool_get_dimension(state,'index_qs',index_qs)
 call mpas_pool_get_dimension(state,'index_qg',index_qg)

 call mpas_pool_get_dimension(state,'index_nc',index_nc)
 call mpas_pool_get_dimension(state,'index_nr',index_nr)
 call mpas_pool_get_dimension(state,'index_ni',index_ni)
 call mpas_pool_get_dimension(state,'index_ns',index_ns)
 call mpas_pool_get_dimension(state,'index_ng',index_ng)

 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 qv => scalars(index_qv,:,:)
 qc => scalars(index_qc,:,:)
 qr => scalars(index_qr,:,:)
 qi => scalars(index_qi,:,:)
 qs => scalars(index_qs,:,:)
 qg => scalars(index_qg,:,:)

 nc => scalars(index_nc,:,:)
 nr => scalars(index_nr,:,:)
 ni => scalars(index_ni,:,:)
 ns => scalars(index_ns,:,:)
 ng => scalars(index_ng,:,:)


 call mpas_pool_get_array(state,'theta_m',theta,time_lev)

 call mpas_pool_get_array(diag,'exner'        ,exner     )
 call mpas_pool_get_array(diag,'pressure_base',pressure_b)
 call mpas_pool_get_array(diag,'pressure_p'   ,pressure_p)


!--- initialization of PUMAS variables:
!mgncol = this%mgncol
!nlev   = this%nlev

 q   => this%q
 qcn => this%qcn
 qrn => this%qrn
 qin => this%qin
 qsn => this%qsn
 qgn => this%qgn
 ncn => this%ncn
 nrn => this%nrn
 nin => this%nin
 nsn => this%nsn
 ngn => this%ngn

 press  => this%press
 pthick => this%pthick
 temp   => this%temp

 delta => this%delta
 delta = dt_dyn

 do k = 1,nlev
    do i = 1,mgncol
       !--- input mixing ratios:
       q(i,k)   = qv(k,i)
       qcn(i,k) = qc(k,i)
       qrn(i,k) = qr(k,i)
       qin(i,k) = qi(k,i)
       qsn(i,k) = qs(k,i)
       qgn(i,k) = qg(k,i)

       !--- input number concentrations:
       ncn(i,k) = nc(k,i)
       nrn(i,k) = nr(k,i)
       nin(i,k) = ni(k,i)
       nsn(i,k) = ns(k,i)
       ngn(i,k) = ng(k,i)
    enddo
 enddo

 do k = 1,nlev
    do i = 1,mgncol
       temp(i,k)  = theta(k,i)/(1.+ R_v/R_d*max(0._R8KIND,qv(k,i)))
       temp(i,k)  = temp(i,k)*exner(k,i)
       press(i,k) = pressure_p(k,i) + pressure_b(k,i)
    enddo
 enddo

 end subroutine microphysics_pumas_from_MPAS

!======================================================================================================================
 subroutine microphysics_pumas_to_MPAS()
!======================================================================================================================

 end subroutine microphysics_pumas_to_MPAS

!======================================================================================================================
 end module mpas_atmphys_driver_microphysics_pumas
!======================================================================================================================

