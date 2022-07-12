!=================================================================================================================
 module mpas_atmphys_vars_noahmp
 use mpas_kind_types

 implicit none
 public
 save


!=================================================================================================================
!.. variables and arrays related to land-surface parameterization, only used in the NOAH MP scheme: additional
!.. description of the iopt_* options are in README.namelist in the WRF sourcecode.
!=================================================================================================================

!.. default options as in README.namelist from WRF version 3.9.1.1:
 integer,parameter:: idveg     = 4 ! Noah-MP dynamic vegetation option.
 integer,parameter:: iopt_crs  = 1 ! Noah-MP stomatal resistance option.
 integer,parameter:: iopt_sfc  = 1 ! Noah-MP surface layer drag coefficient calculation.
 integer,parameter:: iopt_btr  = 1 ! Noah-MP soil moisture factor for stomatal resistance.
 integer,parameter:: iopt_run  = 3 ! Noah-MP runoff and groundwater option.
 integer,parameter:: iopt_frz  = 1 ! Noah-MP supercooled liquid water option.
 integer,parameter:: iopt_inf  = 1 ! Noah-MP soil permeability option.
 integer,parameter:: iopt_rad  = 3 ! Noah-MP radiative transfer option.
 integer,parameter:: iopt_alb  = 2 ! Noah-MP ground surface albedo option.
 integer,parameter:: iopt_snf  = 1 ! Noah-MP precipitation partitioning between snow and rain.
 integer,parameter:: iopt_tbot = 2 ! Noah-MP Soil Temperature Lower Boundary Condition.
 integer,parameter:: iopt_stc  = 1 ! Noah-MP Snow/Soil temperature time scheme.
 integer,parameter:: iopt_gla  = 1 ! Noah-MP glacier treatment option (new in WRF 3.8)
 integer,parameter:: iopt_rsf  = 1 ! Noah-MP surface evaporation resistance option (new in WRF 3.8).
 integer,parameter:: iopt_soil = 1 ! Noah-MP options for defining soil properties.
 integer,parameter:: iopt_pedo = 1 ! Noah-MP options for pedotransfer functions (used when OPT_SOIL = 3).
 integer,parameter:: iopt_crop = 0 ! Noah-MP options for crop model.
 integer,parameter:: iopt_irr  = 0 !
 integer,parameter:: iopt_irrm = 0 !


!inputs:
 character(len=StrKIND):: &
    llanduse           !

 integer:: &
    num_soilcs         !number of soil composition textures.

 integer,dimension(:,:),allocatable:: &
    cropcat_p          !vegetation/soil characteristics (crop model).

 real(kind=RKIND),parameter:: &
    soiltstep = 0      !soil timestep(s), default=0 ->same as main timestep

 real(kind=RKIND),dimension(:,:),allocatable:: &
    soilcl1_p,        &!
    soilcl2_p,        &!
    soilcl3_p,        &!
    soilcl4_p          !
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    soilcomp_p         !

 real(kind=RKIND),dimension(:,:),allocatable:: &
    irfract_p,        &!irrigation fraction (irrigation model).
    sifract_p,        &!sprinkler irrigation fraction (irrigation model).
    mifract_p,        &!micro irrigation fraction (irrigation model).
    fifract_p          !flood irrigation fraction (irrigation model).


!inouts:
 integer,dimension(:,:),allocatable:: &
    isnowxy_p,        &!
    pgsxy_p,          &!
    irnumsi_p,        &!
    irnummi_p,        &!
    irnumfi_p          !

 real(kind=RKIND),dimension(:,:),allocatable:: &
    irwatsi_p,        &!
    irwatmi_p,        &!
    irwatfi_p,        &!
    ireloss_p,        &!
    irsivol_p,        &!
    irmivol_p,        &!
    irfivol_p,        &!
    irrsplh_p          !

 real(kind=RKIND),dimension(:,:),allocatable:: &
    rs_p,             &!surface resistance (s m-1)
    tvxy_p,           &!vegetation leaf temperature
    tgxy_p,           &!bulk ground surface temperature
    canicexy_p,       &!canopy-intercepted ice (mm)
    canliqxy_p,       &!canopy-intercepted liquid water (mm)
    chstarxy_p,       &!dummy exchange coefficient
    eahxy_p,          &!canopy air vapor pressure (pa)
    tahxy_p,          &!canopy air temperature (k)
    cmxy_p,           &!bulk momentum drag coefficient
    chxy_p,           &!bulk sensible heat exchange coefficient
    fwetxy_p,         &!wetted or snowed fraction of the canopy (-)
    sneqvoxy_p,       &!snow mass at last time step(mm h2o)
    alboldxy_p,       &!snow albedo at last time step (-)
    qrainxy_p,        &!rainfall on the ground (mm/s)
    qsnowxy_p,        &!snowfall on the ground (mm/s)
    wslakexy_p,       &!lake water storage (mm)
    zwtxy_p,          &!water table depth [m)
    waxy_p,           &!water in the "aquifer" [mm]
    wtxy_p,           &!groundwater storage [mm]
    lfmassxy_p,       &!leaf mass [g/m2]
    rtmassxy_p,       &!mass of fine roots [g/m2]
    stmassxy_p,       &!stem mass [g/m2]
    woodxy_p,         &!mass of wood (incl. woody roots) [g/m2]
    grainxy_p,        &!mass of grain XING [g/m2]
    planting_p,       &!planting date
    harvest_p,        &!harvest date
    season_gdd_p,     &!growing season GDD
    gddxy_p,          &!growing degree days XING (based on 10C)
    stblcpxy_p,       &!stable carbon in deep soil [g/m2]
    fastcpxy_p,       &!short-lived carbon, shallow soil [g/m2]
    xlaixy_p,         &!leaf area index
    xsaixy_p,         &!stem area index
    taussxy_p,        &!snow age factor
    smcwtdxy_p,       &!soil moisture content in the layer to the water table when deep
    deeprechxy_p,     &!recharge to the water table when deep
    rechxy_p,         &!recharge to the water table (diagnostic)
    qtdrain_p,        &!tile drainage (mm)
    tdfraction_p,     &!
    acc_ssoilxy_p,    &!
    acc_qinsurxy_p,   &!
    acc_qsevaxy_p,    &!
    acc_dwaterxy_p,   &!
    acc_prcpxy_p,     &!
    acc_ecanxy_p,     &!
    acc_etranxy_p,    &!
    acc_edirxy_p       !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    acc_etranixy_p     !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    tsnoxy_p,         &!snow temperature [K]
    zsnsoxy_p,        &!snow layer depth [m]
    snicexy_p,        &!snow layer ice [mm]
    snliqxy_p,        &!snow layer liquid water [mm]
    smoiseq_p          !eq volumetric soil moisture [m3/m3]
 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    croptype_p         !croptype


!outputs:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    t2mvxy_p,         &!2m temperature of vegetation part
    t2mbxy_p,         &!2m temperature of bare ground part
    q2mvxy_p,         &!2m mixing ratio of vegetation part
    q2mbxy_p,         &!2m mixing ratio of bare ground part
    tradxy_p,         &!surface radiative temperature (K)
    neexy_p,          &!net ecosys exchange (g/m2/s CO2)
    gppxy_p,          &!gross primary assimilation [g/m2/s C]
    nppxy_p,          &!net primary productivity [g/m2/s C]
    fvegxy_p,         &!noah-mp vegetation fraction [-]
    runsfxy_p,        &!surface runoff [mm/s]
    runsbxy_p,        &!subsurface runoff [mm/s]
    ecanxy_p,         &!evaporation of intercepted water (mm/s)
    edirxy_p,         &!soil surface evaporation rate (mm/s]
    etranxy_p,        &!transpiration rate (mm/s)
    fsaxy_p,          &!total absorbed solar radiation (W/m2)
    firaxy_p,         &!total net longwave rad (W/m2) [+ to atm]
    aparxy_p,         &!photosyn active energy by canopy (w/m2)
    psnxy_p,          &!total photosynthesis (umol co2/m2/s) [+]
    savxy_p,          &!solar rad absorbed by veg. (w/m2)
    sagxy_p,          &!solar rad absorbed by ground (w/m2)
    rssunxy_p,        &!sunlit leaf stomatal resistance (s/m)
    rsshaxy_p,        &!shaded leaf stomatal resistance (s/m)
    bgapxy_p,         &!between gap fraction
    wgapxy_p,         &!within gap fraction
    tgvxy_p,          &!under canopy ground temperature[k]
    tgbxy_p,          &!bare ground temperature [k]
    chvxy_p,          &!sensible heat exchange coefficient vegetated
    chbxy_p,          &!sensible heat exchange coefficient bare-ground
    shgxy_p,          &!veg ground sen. heat [w/m2]   [+ to atm]
    shcxy_p,          &!canopy sen. heat [w/m2]   [+ to atm]
    shbxy_p,          &!bare sensible heat [w/m2]     [+ to atm]
    evgxy_p,          &!veg ground evap. heat [w/m2]  [+ to atm]
    evbxy_p,          &!bare soil evaporation [w/m2]  [+ to atm]
    ghvxy_p,          &!veg ground heat flux [w/m2]  [+ to soil]
    ghbxy_p,          &!bare ground heat flux [w/m2] [+ to soil]
    irgxy_p,          &!veg ground net lw rad. [w/m2] [+ to atm]
    ircxy_p,          &!canopy net lw rad. [w/m2] [+ to atm]
    irbxy_p,          &!bare net longwave rad. [w/m2] [+ to atm]
    trxy_p,           &!transpiration [w/m2]  [+ to atm]
    evcxy_p,          &!canopy evaporation heat [w/m2]  [+ to atm]
    chleafxy_p,       &!leaf exchange coefficient
    chucxy_p,         &!under canopy exchange coefficient
    chv2xy_p,         &!veg 2m exchange coefficient
    chb2xy_p           !bare 2m exchange coefficient


!additional output variables:
 real(kind=RKIND),dimension(:,:),allocatable:: &
    pahxy_p,         &!
    pahgxy_p,        &!
    pahbxy_p,        &!
    pahvxy_p,        &!
    qintsxy_p,       &!
    qintrxy_p,       &!
    qdripsxy_p,      &!
    qdriprxy_p,      &!
    qthrosxy_p,      &!
    qthrorxy_p,      &!
    qsnsubxy_p,      &!
    qsnfroxy_p,      &!
    qsubcxy_p,       &!
    qfrocxy_p,       &!
    qevacxy_p,       &!
    qdewcxy_p,       &!
    qfrzcxy_p,       &!
    qmeltcxy_p,      &!
    qsnbotxy_p,      &!
    pondingxy_p,     &!
    fpicexy_p,       &!
    rainlsm_p,       &!
    snowlsm_p,       &!
    forctlsm_p,      &!
    forcqlsm_p,      &!
    forcplsm_p,      &!
    forczlsm_p,      &!
    forcwlsm_p,      &!
    eflxbxy_p,       &!
    soilenergy_p,    &!
    snowenergy_p,    &!
    canhsxy_p         !

 real(kind=RKIND),dimension(:,:,:),allocatable:: &
    gecros_state_p    !


!for WRF_HYDRO only:
 real(kind=RKIND),dimension(:,:),allocatable:: &
!   sfcheadrt_p,     &!surface water depth [mm]
!   infxsrt_p,       &!time step infiltration excess [mm]
!   soldrain_p,      &!soil column drainage [mm]
    qtiledrain_p,    &!tile drainage [mm]
    zwatble2d_p       !water table depth for tile drainage [m]

!=================================================================================================================
 end module mpas_atmphys_vars_noahmp
!=================================================================================================================

