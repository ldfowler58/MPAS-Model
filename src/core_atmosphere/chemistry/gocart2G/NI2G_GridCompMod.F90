! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module NI2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: MPAS_LOG_CRIT
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: NIthermo,NIheterogenousChem,Chem_SettlingSimpleOrig,Chem_SettlingSimple, &
                            DryDeposition,WetRemovalGOCART2G,Aero_Compute_Diags
 use DU2G_GridCompMod
 use SS2G_GridCompMod

 use NI2G_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                         fnum,rhFlag,pressure_lid_in_hPa
 use NI2G_StateSpecs,only: NI2G_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: Avogadro = 6.02214076e23
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: runiv    = 8314.462618_RKIND
 real(kind=RKIND),parameter:: undefval = 1.0e15



!--- land/ocean/sea-ice mask (these parameters needs to be accessed from MPAS phys instead of redefined here):
 real(kind=RKIND),parameter:: OCEAN   = 2._RKIND, &
                              LAND    = 1._RKIND, &
                              SEA_ICE = 1._RKIND


!--- relative position of nitrate tracers
 integer,parameter:: nNH3    = 1, &
                     nNH4a   = 2, &
                     nNO3an1 = 3, &
                     nNO3an2 = 4, &
                     nNO3an3 = 5


!--- molecular weights of nitrate species (grams) and dry air:
 real(kind=RKIND),parameter:: fMassAir  = 28.97_RKIND, &
                              fMassHNO3 = 63._RKIND,   &
                              fMassNO3  = 62._RKIND


!--- types needed to define NI2G:
 type:: ThreadWorkspace
    logical:: first = .true.
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: NI2G_GridComp
    logical:: recycle_HNO3 = .false.
    real(kind=RKIND),dimension(:),allocatable:: rmedDU,rmedSS ! DU and SS radius
    real(kind=RKIND),dimension(:),allocatable:: fnumDU,fnumSS ! DU and SS particles per kg mass

    type(ThreadWorkspace),dimension(:),allocatable:: workspaces

    contains
       procedure:: emissions_GridComp => emissions_NI2G_GridComp
       procedure:: load_GridComp      => load_NI2G_GridComp
       procedure:: processes_GridComp => processes_NI2G_GridComp
       procedure:: rrtmg_GridComp     => rrtmg_NI2G_GridComp
 end type NI2G_GridComp

 type wrap_
    type(NI2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!==================================================================================================================
 subroutine load_NI2G_GridComp(self,du2G_params,ss2G_params,kts,kte)
!==================================================================================================================

!--- input arguments:
 type(DU2G_GridComp),intent(in):: du2G_params
 type(SS2G_GridComp),intent(in):: ss2G_params
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(NI2G_GridComp),intent(inout) :: self

!local variables:
 integer:: n,du_nbins,ss_nbins
 real(kind=RKIND),dimension(:),allocatable:: du_fnum,ss_fnum
 real(kind=RKIND),dimension(:),allocatable:: du_radius,ss_radius

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine load_NI2G_GridComp:')

 du_nbins = du2G_params%nbins
 if(.not.allocated(du_fnum)) allocate(du_fnum(du_nbins))
 if(.not.allocated(du_radius)) allocate(du_radius(du_nbins))
 du_fnum = du2G_params%fnum
 du_radius = du2G_params%radius

 ss_nbins = ss2G_params%nbins
 if(.not.allocated(ss_fnum)) allocate(ss_fnum(ss_nbins))
 if(.not.allocated(ss_radius)) allocate(ss_radius(ss_nbins))
 ss_fnum = ss2G_params%fnum
 ss_radius = ss2G_params%radius


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization of variables in GA_Environment:
 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)


!--- initialization of all other variables in NI2_GridComp:
 if(.not.allocated(self%fnumDU)) allocate(self%fnumDU(du_nbins))
 if(.not.allocated(self%rmedDU)) allocate(self%rmedDU(du_nbins))
 if(.not.allocated(self%fnumSS)) allocate(self%fnumSS(ss_nbins))
 if(.not.allocated(self%rmedSS)) allocate(self%rmedSS(ss_nbins))
 do n = 1,du_nbins
    self%fnumDU(n) = du_fnum(n)
    self%rmedDU(n) = du_radius(n)
 enddo
 do n = 1,ss_nbins
    self%fnumSS(n) = ss_fnum(n)
    self%rmedSS(n) = ss_radius(n)
 enddo


 if(.not.allocated(du_fnum)) deallocate(du_fnum)
 if(.not.allocated(ss_fnum)) deallocate(ss_fnum)
 if(.not.allocated(du_radius)) deallocate(du_radius)
 if(.not.allocated(du_radius)) deallocate(ss_radius)

!call mpas_log_write('--- end subroutine load_NI2G_GridCOMP.')

 end subroutine load_NI2G_GridComp

!==================================================================================================================
 subroutine emissions_NI2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 class(NI2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(NI2G_State),intent(inout):: self

!--- local variables:
 integer:: i,j

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_NI2G_GridComp:')


 if(associated(self%nh3em)) then
    self%nh3em = 0.
    if(associated(self%emi_nh3_bb) ) self%nh3em = self%nh3em + self%emi_nh3_bb
    if(associated(self%emi_nh3_ag) ) self%nh3em = self%nh3em + self%emi_nh3_ag
    if(associated(self%emi_nh3_en) ) self%nh3em = self%nh3em + self%emi_nh3_en
    if(associated(self%emi_nh3_tr) ) self%nh3em = self%nh3em + self%emi_nh3_tr
    if(associated(self%emi_nh3_re) ) self%nh3em = self%nh3em + self%emi_nh3_re
    if(associated(self%emi_nh3_in) ) self%nh3em = self%nh3em + self%emi_nh3_in
    if(associated(self%emi_nh3_oc) ) self%nh3em = self%nh3em + self%emi_nh3_oc
 endif

 if(associated(self%emi_nh3_bb)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_bb
 if(associated(self%emi_nh3_ag)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_ag
 if(associated(self%emi_nh3_en)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_en
 if(associated(self%emi_nh3_in)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_in
 if(associated(self%emi_nh3_re)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_re
 if(associated(self%emi_nh3_tr)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_tr
 if(associated(self%emi_nh3_oc)) self%nh3(:,:,self_params%km) = &
    self%nh3(:,:,self_params%km) + self_params%cdt*grav/self%delp(:,:,self_params%km)*self%emi_nh3_oc


 call mpas_log_write('--- end subroutine emissions_NI2G_GridComp.')

 end subroutine emissions_NI2G_GridComp

!==================================================================================================================
 subroutine processes_NI2G_GridComp(self_params,self,to_MYNN,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 logical,intent(in):: to_MYNN
 integer,intent(in):: its,ite,jts,jte,kts,kte
 class(NI2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(NI2G_State),intent(inout):: self

!--- local variables:
 logical:: KIN

 integer:: i,i1,i2,j,j1,j2,k,km,n
 integer:: istat
 integer:: rhFlag_l
 integer:: nw_profile,nw_vertint

 real(kind=RKIND):: fwet
 real(kind=RKIND),dimension(:,:,:),allocatable,target:: fluxoutWT
 real(kind=RKIND),dimension(:,:),pointer:: flux_ptr
 real(kind=RKIND),dimension(:,:,:),pointer:: fluxWT_ptr
 real(kind=RKIND),dimension(:,:),allocatable:: drydepf,dqa
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: aerosol

 real(kind=RKIND),allocatable,dimension(:,:,:),target:: rh20,rh80

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine processes_NI2G_GridComp:')

 nw_profile = size(self_params%wavelengths_profile)
 nw_vertint = size(self_params%wavelengths_vertint)


!---
!call mpas_log_write('--- enter subroutine NIthermo:')
 if(associated(self%nipno3aq)) self%nipno3aq(:,:) = 0._RKIND
 if(associated(self%nipnh4aq)) self%nipnh4aq(:,:) = 0._RKIND
 if(associated(self%nipnh3aq)) self%nipnh3aq(:,:) = 0._RKIND
 istat = 0
 call NIthermo( &
           km        = self_params%km , klid      = self_params%klid , cdt       = self_params%cdt  , &
           grav      = grav           , delp      = self%delp        , rhoa      = self%airdens     , &
           tmpu      = self%t         , rh        = self%rh2         , fMassHNO3 = fmassHNO3        , &
           fMassAir  = fMassAir       , no3an1    = self%no3an1      , nh3       = self%nh3         , &
           nh4a      = self%nh4a      , xhno3     = self%xhno3       , so4       = self%so4         , &
           ni_pno3aq = self%nipno3aq  , ni_pnh4aq = self%nipnh4aq    , ni_pnh3aq = self%nipnh3aq    , &
           rc        = istat                                                                          &
              )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine NIthermo.',messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine NIthermo.')
 endif


!---
!call mpas_log_write('--- enter subroutine NIheterogenous:')
 if(associated(self%hno3conc) ) self%hno3conc  = 0._RKIND
 if(associated(self%hno3smass)) self%hno3smass = 0._RKIND
 if(associated(self%hno3cmass)) self%hno3cmass = 0._RKIND
 if(associated(self%niht     )) self%niht      = 0._RKIND
 istat = 0
 call NIheterogenousChem( &
           km           = self_params%km           , klid         = self_params%klid         , &
           cdt          = self_params%cdt          , rmedDU       = self_params%rmedDU*1.e-6 , &
           rmedSS       = self_params%rmedSS*1.e-6 , fnumDU       = self_params%fnumDU       , &
           fnumSS       = self_params%fnumSS       , grav         = grav                     , &
           avogad       = Avogadro                 , undef        = undefval                 , &
           pi           = pi                       , runiv        = runiv/1000.              , &
           airmw        = fMassAir                 , fMassHNO3    = fMassHNO3                , &
           fMassNO3     = fMassNO3                 , delp         = self%delp                , &
           rhoa         = self%airdens             , tmpu         = self%t                   , &
           relhum       = self%rh2                 , nno3an1      = self%no3an1              , &
           nno3an2      = self%no3an2              , nno3an3      = self%no3an3              , &
           xhno3        = self%xhno3               , du           = self%du                  , &
           ss           = self%ss                  , hno3_conc    = self%hno3conc            , &
           hno3_sfcmass = self%hno3smass           , hno3_colmass = self%hno3cmass           , &
           ni_phet      = self%niht                , rc           = istat                      &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine NIheterogenousChem.',messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine NIheterogenousChem.')
 endif


!--- NI2G settling: because different bins have different swelling coefficients, chemical settling
!    must treated as a function of bins.

!- ammonium ion settles like ammonium sulfate.
!call mpas_log_write('--- enter subroutine Chem_SettlingSimple NH4a:')
 if(associated(self%nh4sd)) self%nh4sd  = 0._RKIND
 istat = 0
 rhflag_l = 3
!call mpas_log_write('--- rhop(nNH4a)   = $r',realArgs=(/self_params%rhop(nNH4a)/))
!call mpas_log_write('--- radius(nNH4a) = $r',realArgs=(/self_params%radius(nNH4a)/))
 call Chem_SettlingSimple( &
           km      = self_params%km          , klid      = self_params%klid                , &
           flag    = rhFlag_l                , cdt       = self_params%cdt                 , &
           grav    = grav                    , radiusInp = self_params%radius(nNH4a)*1.e-6 , &
           rhopInp = self_params%rhop(nNH4a) , tmpu      = self%t                          , &
           rhoa    = self%airdens            , rh        = self%rh2                        , &
           hghte   = self%zle                , delp      = self%delp                       , &
           int_qa  = self%nh4a               , fluxout   = self%nh4sd                      , &
           rc      = istat                                                                   &
                         )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Chem_SettlingSimple NH4a.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Chem_SettlingSimple NH4a.')
 endif


!- nitrate bin 1 settles like ammonium sulfate:
!call mpas_log_write('--- enter subroutine Chem_SettlingSimple NO3AN1:')
 if(associated(self%nisd)) self%nisd(:,:,1) = 0._RKIND
 istat = 0
 rhflag_l = 3
!call mpas_log_write('--- rhop(nNO3AN1)   = $r',realArgs=(/self_params%rhop(nNO3AN1)/))
!call mpas_log_write('--- radius(nNO3AN1) = $r',realArgs=(/self_params%radius(nNO3AN1)/))
 nullify(flux_ptr)
 if(associated(self%nisd)) flux_ptr => self%nisd(:,:,1)
 call Chem_SettlingSimple( &
           km      = self_params%km            , klid      = self_params%klid                  , &
           flag    = rhFlag_l                  , cdt       = self_params%cdt                   , &
           grav    = grav                      , radiusInp = self_params%radius(nNO3an1)*1.e-6 , &
           rhopInp = self_params%rhop(nNO3an1) , tmpu      = self%t                            , &
           rhoa    = self%airdens              , rh        = self%rh2                          , &
           hghte   = self%zle                  , delp      = self%delp                         , &
           int_qa  = self%no3an1               , fluxout   = flux_ptr                          , &
           rc      = istat                                                                       &
                         )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Chem_SettlingSimple NO3AN1.', &
           messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Chem_SettlingSimple NO3AN1.')
 endif


!- nitrate bin 2 settles like sea salt:
!call mpas_log_write('--- enter subroutine Chem_SettlingSimple NO3AN2:')
 if(associated(self%nisd)) self%nisd(:,:,2) = 0._RKIND
 istat = 0
 rhflag_l = 2
!call mpas_log_write('--- rhop(nNO3AN2)   = $r',realArgs=(/self_params%rhop(nNO3AN2)/))
!call mpas_log_write('--- radius(nNO3AN2) = $r',realArgs=(/self_params%radius(nNO3AN2)/))
 nullify(flux_ptr)
 if(associated(self%nisd)) flux_ptr => self%nisd(:,:,2)
 call Chem_SettlingSimple( &
           km      = self_params%km            , klid      = self_params%klid                  , &
           flag    = rhFlag_l                  , cdt       = self_params%cdt                   , &
           grav    = grav                      , radiusInp = self_params%radius(nNO3an2)*1.e-6 , &
           rhopInp = self_params%rhop(nNO3an2) , tmpu      = self%t                            , &
           rhoa    = self%airdens              , rh        = self%rh2                          , &
           hghte   = self%zle                  , delp      = self%delp                         , &
           int_qa  = self%no3an2               , fluxout   = flux_ptr                          , &
           rc      = istat                                                                       &
                         )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Chem_SettlingSimple NO3AN2.', &
           messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Chem_SettlingSimple NO3AN2.')
 endif


!- nitrate bin 3 settles like dust:
!call mpas_log_write('--- enter subroutine Chem_SettlingSimple NO3AN3:')
 if(associated(self%nisd)) self%nisd(:,:,3) = 0._RKIND
 istat = 0
 rhflag_l = 0
!call mpas_log_write('--- rhop(nNO3AN3)   = $r',realArgs=(/self_params%rhop(nNO3AN3)/))
!call mpas_log_write('--- radius(nNO3AN3) = $r',realArgs=(/self_params%radius(nNO3AN3)/))
 nullify(flux_ptr)
 if(associated(self%nisd)) flux_ptr => self%nisd(:,:,3)
 call Chem_SettlingSimple( &
           km      = self_params%km            , klid      = self_params%klid                  , &
           flag    = rhFlag_l                  , cdt       = self_params%cdt                   , &
           grav    = grav                      , radiusInp = self_params%radius(nNO3an3)*1.e-6 , &
           rhopInp = self_params%rhop(nNO3an3) , tmpu      = self%t                            , &
           rhoa    = self%airdens              , rh        = self%rh2                          , &
           hghte   = self%zle                  , delp      = self%delp                         , &
           int_qa  = self%no3an3               , fluxout   = flux_ptr                          , &
           rc      = istat                                                                       &
                         )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Chem_SettlingSimple NO3AN3.', &
           messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Chem_SettlingSimple NO3AN3.')
 endif


!--- NI2G dry deposition:
!call mpas_log_write('--- enter subroutine DryDeposition:')
 if(associated(self%nh3dp)) self%nh3dp(:,:) = 0._RKIND
 if(associated(self%nh3dp)) self%nh4dp(:,:) = 0._RKIND
 if(associated(self%nh3dp)) self%nidp(:,:,:)  = 0._RKIND
 if(associated(self%nivdep)) self%nivdep(:,:) = 0._RKIND
 if(.not.allocated(dqa)    ) allocate(dqa(its:ite,jts:jte)    )
 if(.not.allocated(drydepf)) allocate(drydepf(its:ite,jts:jte))
 drydepf = 0._RKIND
 istat   = 0
 call DryDeposition( &
              km         = self_params%km , &
              tmpu       = self%t         , &
              rhoa       = self%airdens   , &
              hghte      = self%zle       , &
              oro        = self%lwi       , &
              ustar      = self%ustar     , &
              pblh       = self%zpbl      , &
              shflux     = self%sh        , &
              von_karman = karman         , &
              cpd        = cpd            , &
              grav       = grav           , &
              z0h        = self%z0h       , &
              drydepf    = drydepf        , &
              rc         = istat            &
                   )
!--- if gocart2G does not interact with the MYNN PBL parameterization, then we update the mixing ratios
!    due to dry deposition in the model layer adjacent to the surface, otherwise we simply skip this step:
 if(associated(self%nivdep)) self%nivdep(:,:) = drydepf(:,:)*self%delz(:,:,self_params%km)
 if(.not. to_MYNN) then
    !- nh3:
    dqa = 0._RKIND
    do i = 1,ubound(self%lwi,1)
       do j = 1,ubound(self%lwi,2)
          if(abs(self%lwi(i,j) - OCEAN) < 0.5) then
             dqa(i,j) = max(0._RKIND,self%nh3(i,j,self_params%km)*(1.-exp(-10.*drydepf(i,j)*self_params%cdt)))
          else
             dqa(i,j) = max(0._RKIND,self%nh3(i,j,self_params%km)*(1.-exp(-3.*drydepf(i,j)*self_params%cdt)))
          endif
       enddo
    enddo
    self%nh3(:,:,self_params%km) = self%nh3(:,:,self_params%km) - dqa
    if(associated(self%nh3dp)) self%nh3dp = dqa(:,:)*self%delp(:,:,self_params%km)/grav/self_params%cdt

    !- nh4a:
    dqa = 0._RKIND
    dqa = max(0._RKIND,self%nh4a(:,:,self_params%km)*(1.-exp(-drydepf*self_params%cdt)))
    self%nh4a(:,:,self_params%km) = self%nh4a(:,:,self_params%km) - dqa
    if(associated(self%nh4dp)) self%nh4dp = dqa(:,:)*self%delp(:,:,self_params%km)/grav/self_params%cdt

    !- no3anx:
    dqa = 0._RKIND
    dqa = max(0._RKIND,self%no3an1(:,:,self_params%km)*(1.-exp(-drydepf*self_params%cdt)))
    self%no3an1(:,:,self_params%km) = self%no3an1(:,:,self_params%km) - dqa
    if(associated(self%nidp)) self%nidp(:,:,1) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt

    dqa = 0._RKIND
    dqa = max(0._RKIND,self%no3an2(:,:,self_params%km)*(1.-exp(-drydepf*self_params%cdt)))
    self%no3an2(:,:,self_params%km) = self%no3an2(:,:,self_params%km) - dqa
    if(associated(self%nidp)) self%nidp(:,:,2) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt

    dqa = 0._RKIND
    dqa = max(0._RKIND,self%no3an3(:,:,self_params%km)*(1.-exp(-drydepf*self_params%cdt)))
    self%no3an3(:,:,self_params%km) = self%no3an3(:,:,self_params%km) - dqa
    if(associated(self%nidp)) self%nidp(:,:,3) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt
 endif
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine DryDeposition.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine DryDeposition.')
 endif


!--- NI2G large-scale wet removal:
!call mpas_log_write('--- enter subroutine WetRemovalGOCART2G NH3:')
 if(associated(self%nh3wt) .or. associated(self%nh4wt)) allocate(fluxoutWT(ubound(self%t,1),ubound(self%t,2),1))

 istat = 0
 KIN  = .false.
 fwet = 1._RKIND
 nullify(fluxWT_ptr)
 if(associated(self%nh3wt)) self%nh3wt = 0._RKIND
 if(associated(self%nh3wt)) fluxWT_ptr => fluxoutWT
 call WetRemovalGOCART2G( &
           km        = self_params%km    , klid    = self_params%klid , n1      = self_params%nbins , &
           n2        = self_params%nbins , bin_ind = 1                , cdt     = self_params%cdt   , &
           aero_type = 'NH3'             , kin     = KIN              , grav    = grav              , &
           fwet      = fwet              , aerosol = self%nh3         , ple     = self%ple          , &
           tmpu      = self%t            , rhoa    = self%airdens     , pfllsan = self%pfl_lsan     , &
           pfilsan   = self%pfi_lsan     , precc   = self%cn_prcp     , precl   = self%ncn_prcp     , &
           fluxout   = fluxWT_ptr        , rc      = istat                                            &
                        )
 if (associated(self%nh3wt)) self%nh3wt = fluxWT_ptr(:,:,1)
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine WetRemovalGOCART2G NH3.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine WetRemovalGOCART2G NH3.')
 endif


!call mpas_log_write('--- enter subroutine WetRemovalGOCART2G NH4a:')
 istat = 0
 KIN  = .true.
 fwet = 1._RKIND
 nullify(fluxWT_ptr)
 if(associated(self%nh4wt)) self%nh4wt = 0._RKIND
 if(associated(self%nh4wt)) fluxWT_ptr => fluxoutWT
 call WetRemovalGOCART2G( &
           km        = self_params%km    , klid    = self_params%klid , n1      = self_params%nbins , &
           n2        = self_params%nbins , bin_ind = 1                , cdt     = self_params%cdt   , &
           aero_type = 'NH4a'            , kin     = KIN              , grav    = grav              , &
           fwet      = fwet              , aerosol = self%nh4a        , ple     = self%ple          , &
           tmpu      = self%t            , rhoa    = self%airdens     , pfllsan = self%pfl_lsan     , &
           pfilsan   = self%pfi_lsan     , precc   = self%cn_prcp     , precl   = self%ncn_prcp     , &
           fluxout   = fluxWT_ptr        , rc      = istat                                            &
                        )
 if (associated(self%nh4wt)) self%nh4wt = fluxWT_ptr(:,:,1)
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine WetRemovalGOCART2G NH4a.', &
                        messageType=MPAS_LOG_CRIT)
 else
    if(allocated(fluxoutWT)) deallocate(fluxoutWT)
!   call mpas_log_write('--- end subroutine WetRemovalGOCART2G NH4a.')
 endif


!call mpas_log_write('--- enter subroutine WetRemovalGOCART2G NO3AN1:')
 istat = 0
 KIN  = .true.
 fwet = 1._RKIND
 if(associated(self%niwt)) self%niwt = 0._RKIND
 call WetRemovalGOCART2G( &
           km        = self_params%km    , klid    = self_params%klid , n1      = self_params%nbins , &
           n2        = self_params%nbins , bin_ind = 1                , cdt     = self_params%cdt   , &
           aero_type = 'nitrate'         , kin     = KIN              , grav    = grav              , &
           fwet      = fwet              , aerosol = self%no3an1      , ple     = self%ple          , &
           tmpu      = self%t            , rhoa    = self%airdens     , pfllsan = self%pfl_lsan     , &
           pfilsan   = self%pfi_lsan     , precc   = self%cn_prcp     , precl   = self%ncn_prcp     , &
           fluxout   = self%niwt         , rc      = istat                                            &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine WetRemovalGOCART2G NO3AN1.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine WetRemovalGOCART2G NO3AN1:')
 endif


!call mpas_log_write('--- enter subroutine WetRemovalGOCART2G NO3AN2:')
 istat = 0
 KIN  = .true.
 fwet = 1._RKIND
 call WetRemovalGOCART2G( &
           km        = self_params%km    , klid    = self_params%klid , n1      = self_params%nbins , &
           n2        = self_params%nbins , bin_ind = 2                , cdt     = self_params%cdt   , &
           aero_type = 'nitrate'         , kin     = KIN              , grav    = grav              , &
           fwet      = fwet              , aerosol = self%no3an2      , ple     = self%ple          , &
           tmpu      = self%t            , rhoa    = self%airdens     , pfllsan = self%pfl_lsan     , &
           pfilsan   = self%pfi_lsan     , precc   = self%cn_prcp     , precl   = self%ncn_prcp     , &
           fluxout   = self%niwt         , rc      = istat                                            &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine WetRemovalGOCART2G NO3AN2.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine WetRemovalGOCART2G NO3AN2.')
 endif


!call mpas_log_write('--- enter subroutine WetRemovalGOCART2G NO3AN3:')
 istat = 0
 KIN  = .true.
 fwet = 0.3_RKIND
 call WetRemovalGOCART2G( &
           km        = self_params%km    , klid    = self_params%klid , n1      = self_params%nbins , &
           n2        = self_params%nbins , bin_ind = 3                , cdt     = self_params%cdt   , &
           aero_type = 'nitrate'         , kin     = KIN              , grav    = grav              , &
           fwet      = fwet              , aerosol = self%no3an3      , ple     = self%ple          , &
           tmpu      = self%t            , rhoa    = self%airdens     , pfllsan = self%pfl_lsan     , &
           pfilsan   = self%pfi_lsan     , precc   = self%cn_prcp     , precl   = self%ncn_prcp     , &
           fluxout   = self%niwt         , rc      = istat                                            &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine WetRemovalGOCART2G NO3AN3.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine WetRemovalGOCART2G NO3AN3.')
 endif


!--- NI2G diagnostics:
!Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
!call mpas_log_write('--- nw_profile = $i',intArgs=(/nw_profile/))
!call mpas_log_write('--- nw_vertint = $i',intArgs=(/nw_vertint/))

 if(.not.allocated(aerosol)) allocate(aerosol(ubound(self%nh4a,1),ubound(self%nh4a,2),ubound(self%nh4a,3),3))
 aerosol(:,:,:,:) = 0._RKIND

!call mpas_log_write('--- enter subroutine Aero_Compute_Diags NH4a:')
 aerosol(:,:,:,1) = self%nh4a(:,:,:)
 if(associated(self%nh4smass)) self%nh4smass(:,:)  = 0._RKIND
 if(associated(self%nh4cmass)) self%nh4cmass(:,:)  = 0._RKIND
 if(associated(self%nh4mass) ) self%nh4mass(:,:,:) = 0._RKIND
 if(associated(self%nh4conc) ) self%nh4conc(:,:,:) = 0._RKIND
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 1                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = aerosol                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = self%rh2                               , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           sfcmass             = self%nh4smass                          , &
           colmass             = self%nh4cmass                          , &
           mass                = self%nh4mass                           , &
           conc                = self%nh4conc                           , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Aero_Compute_Diags NH4a.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags NH4a.')
 endif


!call mpas_log_write('--- enter subroutine Aero_Compute_Diags NH3:')
 aerosol(:,:,:,1) = self%nh3(:,:,:)
 if(associated(self%nh3smass)) self%nh3smass(:,:)  = 0._RKIND
 if(associated(self%nh3cmass)) self%nh3cmass(:,:)  = 0._RKIND
 if(associated(self%nh3mass) ) self%nh3mass(:,:,:) = 0._RKIND
 if(associated(self%nh3conc) ) self%nh3conc(:,:,:) = 0._RKIND
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 1                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = aerosol                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = self%rh2                               , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           sfcmass             = self%nh3smass                          , &
           colmass             = self%nh3cmass                          , &
           mass                = self%nh3mass                           , &
           conc                = self%nh3conc                           , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Aero_Compute_Diags NH3.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags NH3.')
 endif


!call mpas_log_write('--- enter subroutine Aero_Compute_Diags NO3AN1:')
 aerosol(:,:,:,1) = self%no3an1(:,:,:)
 if(associated(self%nismass25)) self%nismass25(:,:)  = 0._RKIND
 if(associated(self%nicmass25)) self%nicmass25(:,:)  = 0._RKIND
 if(associated(self%nimass25) ) self%nimass25(:,:,:) = 0._RKIND
 if(associated(self%niconc25) ) self%niconc25(:,:,:) = 0._RKIND
 if(associated(self%niextt25) ) self%niextt25(:,:,:) = 0._RKIND
 if(associated(self%niscat25) ) self%niscat25(:,:,:) = 0._RKIND
 if(associated(self%niexttfm) ) self%niexttfm(:,:,:) = 0._RKIND
 if(associated(self%niscatfm) ) self%niscatfm(:,:,:) = 0._RKIND
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 1                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = aerosol                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = self%rh2                               , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           sfcmass             = self%nismass25                         , &
           colmass             = self%nicmass25                         , &
           mass                = self%nimass25                          , &
           conc                = self%niconc25                          , &
           exttau25            = self%niextt25                          , &
           scatau25            = self%niscat25                          , &
           exttaufm            = self%niexttfm                          , &
           scataufm            = self%niscatfm                          , &
           NO3nFlag            = .true.                                 , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Aero_Compute_Diags NO3AN1.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags NO3AN1.')
 endif


!call mpas_log_write('--- enter subroutine Aero_Compute_Diags NO3AN1 NO3AN2 NO3AN3:')
 aerosol(:,:,:,1) = self%no3an1(:,:,:)
 aerosol(:,:,:,2) = self%no3an2(:,:,:)
 aerosol(:,:,:,3) = self%no3an3(:,:,:)
 if(associated(self%nismass)   ) self%nismass(:,:)       = 0._RKIND
 if(associated(self%nicmass)   ) self%nicmass(:,:)       = 0._RKIND
 if(associated(self%nifluxu)   ) self%nifluxu(:,:)       = 0._RKIND
 if(associated(self%nifluxv)   ) self%nifluxv(:,:)       = 0._RKIND
 if(associated(self%niangstr)  ) self%niangstr(:,:)      = 0._RKIND
 if(associated(self%nimass)    ) self%nimass(:,:,:)      = 0._RKIND
 if(associated(self%niconc)    ) self%niconc(:,:,:)      = 0._RKIND
 if(associated(self%niexttau)  ) self%niexttau(:,:,:)    = 0._RKIND
 if(associated(self%nistexttau)) self%nistexttau(:,:,:)  = 0._RKIND
 if(associated(self%niscatau)  ) self%niscatau(:,:,:)    = 0._RKIND
 if(associated(self%nistscatau)) self%nistscatau(:,:,:)  = 0._RKIND
 if(associated(self%niextcoef) ) self%niextcoef(:,:,:,:) = 0._RKIND
 if(associated(self%niscacoef) ) self%niscacoef(:,:,:,:) = 0._RKIND
 if(associated(self%nibckcoef) ) self%nibckcoef(:,:,:,:) = 0._RKIND
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 3                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = aerosol                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = self%rh2                               , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           sfcmass             = self%nismass                           , &
           colmass             = self%nicmass                           , &
           mass                = self%nimass                            , &
           conc                = self%niconc                            , &
           exttau              = self%niexttau                          , &
           scatau              = self%niscatau                          , &
!          stexttau            = self%nistexttau                        , &
!          stscatau            = self%nistscatau                        , &
           fluxu               = self%nifluxu                           , &
           fluxv               = self%nifluxv                           , &
           extcoef             = self%niextcoef                         , &
           scacoef             = self%niscacoef                         , &
           bckcoef             = self%nibckcoef                         , &
           angstrom            = self%niangstr                          , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Aero_Compute_Diags NO3AN1 NO3AN2 NO3AN3.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags NO3AN1 NO3AN2 NO3AN3.')
 endif


 i1 = lbound(self%rh2,1); i2 = ubound(self%rh2,1)
 j1 = lbound(self%rh2,2); j2 = ubound(self%rh2,2)
 km = ubound(self%rh2,3)

!call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH20:')
 if(.not.associated(self%niextcoefrh20)) self%niextcoefrh20(:,:,:,:) = 0._RKIND
 if(.not.associated(self%niscacoefrh20)) self%niscacoefrh20(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh20)) allocate(rh20(i1:i2,j1:j2,km))
 rh20(:,:,:) = 0.20
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 3                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = aerosol                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh20                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%niextcoefrh20                     , &
           scacoef             = self%niscacoefrh20                     , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Aero_Compute_Diags RH20.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags RH20.')
 endif


!call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH80:')
 if(.not.associated(self%niextcoefrh80)) self%niextcoefrh80(:,:,:,:) = 0._RKIND
 if(.not.associated(self%niscacoefrh80)) self%niscacoefrh80(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh80)) allocate(rh80(i1:i2,j1:j2,km))
 rh80(:,:,:) = 0.80
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 3                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = aerosol                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh80                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%niextcoefrh80                     , &
           scacoef             = self%niscacoefrh80                     , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- NI2G_GridComp: error in subroutine Aero_Compute_Diags RH80.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags RH80.')
 endif
 if(allocated(rh20)) deallocate(rh20)
 if(allocated(rh80)) deallocate(rh80)
 if(allocated(aerosol)) deallocate(aerosol)


 call mpas_log_write('--- end subroutine processes_NI2G_GridCOMP.')

 end subroutine processes_NI2G_GridComp

!==================================================================================================================
 subroutine rrtmg_NI2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte,nbndlw,nbndsw)
!==================================================================================================================
!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nbndlw,nbndsw
 class(NI2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(NI2G_State),intent(inout):: self

!--- local variables:
 integer:: nbands,km,nbins
 integer:: istat
 integer:: i,k,kk,j,n,nl,ns

 real(kind=RKIND),dimension(:,:,:),allocatable :: asy_s,ext_s,ssa_s
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qni3G

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_NI2G_GridComp:')


!--- the calculation of the RRTMG optical properties only includes the three kinds of nitrate out of the five
!    prognostic variables. therefore, we explicitly set nbins to three instead of self_params%nbins:
 km     = self_params%km
 nbins  = 3
 nbands = self_params%rad_Mie%nch
!call mpas_log_write('--- km     = $i',intArgs=(/km/))
!call mpas_log_write('--- nbands = $i',intArgs=(/nbands/))


 if(associated(self%nitau_sw)) self%nitau_sw(:,:,:,:) = 0._RKIND
 if(associated(self%niasy_sw)) self%niasy_sw(:,:,:,:) = 0._RKIND
 if(associated(self%nissa_sw)) self%nissa_sw(:,:,:,:) = 0._RKIND
 if(associated(self%nitau_lw)) self%nitau_lw(:,:,:,:) = 0._RKIND
 if(associated(self%niasy_lw)) self%niasy_lw(:,:,:,:) = 0._RKIND
 if(associated(self%nissa_lw)) self%nissa_lw(:,:,:,:) = 0._RKIND

 if(.not.allocated(asy_s)) allocate(asy_s(its:ite,jts:jte,kts:kte))
 if(.not.allocated(ext_s)) allocate(ext_s(its:ite,jts:jte,kts:kte))
 if(.not.allocated(ssa_s)) allocate(ssa_s(its:ite,jts:jte,kts:kte))
 if(.not.allocated(qni3G)) allocate(qni3G(its:ite,jts:jte,kts:kte,nbins))

 asy_s(:,:,:) = 0._RKIND
 ext_s(:,:,:) = 0._RKIND
 ssa_s(:,:,:) = 0._RKIND

 qni3G(:,:,:,1) = self%no3an1(:,:,:)*self%delp(:,:,:)/grav
 qni3G(:,:,:,2) = self%no3an2(:,:,:)*self%delp(:,:,:)/grav
 qni3G(:,:,:,3) = self%no3an3(:,:,:)*self%delp(:,:,:)/grav


 nl = 0
 ns = 0
 do n = 1, nbands
    istat = 0
    asy_s(:,:,:) = 0._RKIND
    ext_s(:,:,:) = 0._RKIND
    ssa_s(:,:,:) = 0._RKIND
    call mie_(self_params%rad_Mie,its,ite,jts,jte,kts,kte,nbins,n,qni3G,self%rh2,ext_s,ssa_s,asy_s,istat)
    if(istat /=0) then
       call mpas_log_write('--- NI2G_GridComp: error in subroutine rrtmg_NI2G_GridComp.', &
                           messageType=MPAS_LOG_CRIT)
    else
       if(n .le. nbndsw) then
          ns = ns+1
          self%nitau_sw(:,:,:,ns) = real(ext_s(:,:,:),kind=RKIND)
          self%niasy_sw(:,:,:,ns) = real(asy_s(:,:,:),kind=RKIND)
          self%nissa_sw(:,:,:,ns) = real(ssa_s(:,:,:),kind=RKIND)
       elseif(n .gt. nbndsw) then
          nl = nl+1
          self%nitau_lw(:,:,:,nl) = real(ext_s(:,:,:),kind=RKIND)
          self%niasy_lw(:,:,:,nl) = real(asy_s(:,:,:),kind=RKIND)
          self%nissa_lw(:,:,:,nl) = real(ssa_s(:,:,:),kind=RKIND)
       endif
    endif
 enddo


 if(allocated(asy_s)) deallocate(asy_s)
 if(allocated(ext_s)) deallocate(ext_s)
 if(allocated(ssa_s)) deallocate(ssa_s)
 if(allocated(qni3G)) deallocate(qni3G)

 call mpas_log_write('--- end subroutine rrtmg_NI2G_GridComp.')

 return


 contains


    subroutine mie_(mie,its,ite,jts,jte,kts,kte,nbins,ichannel,q,rh,bext_s,bssa_s,basy_s,rc)

    type(GOCART2G_Mie),intent(in):: mie ! Mie table
    integer,intent(in):: its,ite,jts,jte,kts,kte
    integer,intent(in):: nbins          ! number of bins
    integer,intent(in):: ichannel       ! channel
    integer,intent(out) :: rc

    real(kind=RKIND),intent(in),dimension(:,:,:):: rh  ! relative humidity
    real(kind=RKIND),intent(in),dimension(:,:,:,:):: q ! aerosol
    real(kind=RKIND),intent(out):: bext_s(size(ext_s,1),size(ext_s,2),size(ext_s,3)) ! dimensionless.
    real(kind=RKIND),intent(out):: bssa_s(size(ext_s,1),size(ext_s,2),size(ext_s,3)) ! dimensionless.
    real(kind=RKIND),intent(out):: basy_s(size(ext_s,1),size(ext_s,2),size(ext_s,3)) ! dimensionless.

    !local variables and arrays:
    integer:: l,n
    real(kind=RKIND):: bext(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! extinction
    real(kind=RKIND):: bssa (size(ext_s,1),size(ext_s,2),size(ext_s,3)) ! single scattering
    real(kind=RKIND):: gasy(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! asymmetry parameter


    bext_s(:,:,:) = 0._RKIND
    bssa_s(:,:,:) = 0._RKIND
    basy_s(:,:,:) = 0._RKIND

    do l = 1, nbins
       call mie%Query(ichannel,l,q(:,:,:,l),rh,tau=bext,gasym=gasy,ssa=bssa,rc=rc)
       bext_s = bext_s + bext           ! extinction
       bssa_s = bssa_s + bssa*bext      ! scattering extinction
       basy_s = basy_s + gasy*bssa*bext ! asymetry parameter multiplied by scatering extiction
    enddo


    end subroutine mie_

 end subroutine rrtmg_NI2G_GridComp

!==================================================================================================================
 end module NI2G_GridCompMod
!==================================================================================================================
