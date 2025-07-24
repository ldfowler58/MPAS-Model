! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module SS2G_GridCompMod
 use mpas_kind_types,only: RKIND,R8KIND
 use mpas_derived_types,only: MPAS_LOG_CRIT
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: jeagleSSTcorrection, &
                            weibullDistribution, &
                            SeasaltEmission,     &
                            hoppelCorrection,    &
                            Chem_Settling,       &
                            DryDeposition,       &
                            WetRemovalGOCART2G,  &
                            Aero_Compute_Diags

 use SS2G_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight,           &
                         fnum,rhFlag,pressure_lid_in_hPa,emission_scheme,sstEmisFlag,hoppelFlag,          &
                         weibullFlag,radius_lower,radius_upper,particle_radius_number,emission_scale_res
 use SS2G_StateSpecs,only: SS2G_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 integer,parameter:: NHRES = 6

 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: airmw    = 28.97_RKIND
 real(kind=RKIND),parameter:: radTodeg = 180._RKIND/pi
 real(kind=RKIND),parameter:: Avogadro = 6.02214076e23
 real(kind=RKIND),parameter:: undefval = 1.0e15

 real(kind=RKIND),parameter:: OCEAN   = 2._RKIND
 real(kind=RKIND),parameter:: LAND    = 1._RKIND
 real(kind=RKIND),parameter:: SEA_ICE = 1._RKIND


!--- types needed to define SS2G:
 type,extends(GA_Environment),public:: SS2G_GridComp
    logical:: hoppelFlag     ! apply the Hoppel correction to emissions (Fan and Toon, 2011)
    logical:: weibullFlag    ! apply the Weibull distribution to wind speed for emissions (Fan and Toon, 2011)

    integer:: emission_scheme
    integer:: sstEmisFlag    ! choice of SST correction to emissions:
                             ! 0 - none; 1 - Jaegle et al. 2011; 2 - GEOS5

    real(kind=RKIND):: emission_scale                             ! global scaling tuning coefficient
    real(kind=RKIND),dimension(NHRES):: emission_scale_res        ! global scaling tuning resolution

    real(kind=RKIND),dimension(:),allocatable:: rlow              ! particle effective radius lower bound [um]
    real(kind=RKIND),dimension(:),allocatable:: rup               ! particle effective radius upper bound [um]
    real(kind=RKIND),dimension(:),allocatable:: rmed              ! number median radius [um]
!   real(kind=RKIND),dimension(:,:),allocatable:: deep_lakes_mask ! mask for deep lakes

    contains
       procedure:: emissions_GridComp => emissions_SS2G_GridComp
       procedure:: load_GridComp      => load_SS2G_GridComp
       procedure:: processes_GridComp => processes_SS2G_GridComp
       procedure:: rrtmg_GridComp     => rrtmg_SS2G_GridComp
 end type SS2G_GridComp

 type wrap_
    type(SS2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!==================================================================================================================
 subroutine load_SS2G_GridComp(self,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(SS2G_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine load_SS2G_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization of variables in GA_Environment:
 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

!call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
!call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
!do n = 1,self%nbins
!   call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
!                       self%fscav(n),self%molwght(n),self%fnum(n)/))
!enddo


!--- initialization of all other variables in SS2_GridComp:
 self%emission_scheme = emission_scheme
 self%sstEmisFlag     = sstEmisFlag
 self%hoppelFlag      = hoppelFlag
 self%weibullFlag     = weibullFlag

 if(.not.allocated(self%rlow) ) allocate(self%rlow(self%nbins))
 if(.not.allocated(self%rup)  ) allocate(self%rup(self%nbins) )
 if(.not.allocated(self%rmed) ) allocate(self%rmed(self%nbins))

 do n = 1,self%nbins
    self%rlow(n) = radius_lower(n)
    self%rup(n)  = radius_upper(n)
    self%rmed(n) = particle_radius_number(n)
 enddo


!--- initialization of sea-salt emission tuning coefficient (emission_scale) and resolutions used for
!    emission_scale (emission_scale_res). dims(1) and dims(2) are the number of grid-points in the longitude
!    and latitude directions. In MPAS,  we will need to adjust subroutine Chem_UtilResVal since MPAS uses an
!    unstructured grid. For now, we set self%Ch_DU to 1:
 do n = 1, NHRES
    self%emission_scale_res(n) = emission_scale_res(n)
 enddo
!self%emission_scale = Chem_UtilResVal(dims(1),dims(2),self%emission_scale_res(:),stat)
 self%emission_scale = 1._RKIND


!call mpas_log_write('--- end subroutine load_SS2G_GridCOMP.')

 end subroutine load_SS2G_GridComp

!==================================================================================================================
 subroutine emissions_SS2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!==================================================================================================================
!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 class(SS2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(SS2G_State),intent(inout):: self

!--- local variables and arrays:
 integer:: i,j,ibin
 integer:: istat

 real(kind=RKIND),dimension(:,:),allocatable:: fgridefficiency
 real(kind=RKIND),dimension(:,:),allocatable:: fsstemis
 real(kind=RKIND),dimension(:,:),allocatable:: fhoppel,memissions,nemissions,dqa
 real(kind=R8KIND),dimension(:,:),allocatable:: gweibull

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_SS2G_GridComp:')


 if(.not.allocated(fgridefficiency)) allocate(fgridefficiency(its:ite,jts:jte))
 if(.not.allocated(fsstemis)       ) allocate(fsstemis(its:ite,jts:jte)       )
 if(.not.allocated(gweibull)       ) allocate(gweibull(its:ite,jts:jte)       )
 if(.not.allocated(fhoppel)        ) allocate(fhoppel(its:ite,jts:jte)        )
 if(.not.allocated(memissions)     ) allocate(memissions(its:ite,jts:jte)     )
 if(.not.allocated(nemissions)     ) allocate(nemissions(its:ite,jts:jte)     )
 if(.not.allocated(dqa)            ) allocate(dqa(its:ite,jts:jte)            )


!--- grid box efficiency to emission (fraction of sea water)
 fgridefficiency = min(max(0.,(self%frocean-self%fraci)*self%deep_lakes_mask),1.)


!--- apply SST correction (Jaegle et al., 2011):
!call mpas_log_write('--- enter subroutine jeagleSSTcorrection:')
 istat = 0
 call jeagleSSTcorrection( &
    sstEmisFlag = self_params%sstEmisFlag, &
    fsstemis    = fsstemis,                &
    ts          = self%ts,                 &
    rc          = istat                    &
                         )
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine jeagleSSTcorrection.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine jeagleSSTcorrection.')
 endif


!--- apply Weibull distribution to emissions wind speeds (Fan and Toon, 2011):
!call mpas_log_write('--- enter subroutine weibullDistribution:')
 istat = 0
 call weibullDistribution( &
    gweibull    = gweibull,                &
    weibullFlag = self_params%weibullFlag, &
    u10m        = self%u10m,               &
    v10m        = self%v10m,               &
    rc          = istat                    &
                         )
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine weibullDistribution.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine weibullDistribution.')
 endif


!--- seasalt emission and possibly apply Hoppel correction (Fan and Toon, 2011):
!call mpas_log_write('--- enter subroutine SeasaltEmission:')
 fhoppel(:,:) = 1._RKIND
 do ibin = 1, self_params%nbins
    if(associated(self%SSEM)) self%SSEM(:,:,ibin) = 0._RKIND
    memissions(:,:) = 0._RKIND
    nemissions(:,:) = 0._RKIND
    dqa(:,:)        = 0._RKIND

    istat = 0
    call SeasaltEmission( &
       rLow       = self_params%rlow(ibin),      &
       rUP        = self_params%rup(ibin),       &
       method     = self_params%emission_scheme, &
       u10m       = self%u10m,                   &
       v10m       = self%v10m,                   &
       ustar      = self%ustar,                  &
       pi         = pi,                          &
       memissions = memissions,                  &
       nemissions = nemissions,                  &
       rc         = istat                        &
                        )

    !--- Hoppel correction:
    if(self_params%hoppelFlag) then
       call hoppelCorrection( &
          radius  = self_params%radius(ibin)*1.e-6,   &
          rhop    = self_params%rhop(ibin),           &
          rhFlag  = self_params%rhFlag,               &
          rh      = self%rh2(:,:,self_params%km),     &
          dz      = self%delz(:,:,self_params%km),    &
          ustar   = self%ustar,                       &
          airdens = self%airdens(:,:,self_params%km), &
          t       = self%t(:,:,self_params%km),       &
          grav    = grav,                             &
          karman  = karman,                           &
          fhoppel = fhoppel,                          &
          rc      = istat                             &
                            )
    endif

    memissions = self_params%emission_scale*fgridefficiency*fsstemis*fhoppel*gweibull*memissions
    dqa = memissions*self_params%cdt*grav/self%delp(:,:,self_params%km)

    self%SS(:,:,self_params%km,ibin) = self%SS(:,:,self_params%km,ibin) + dqa
    if(associated(self%SSEM)) self%SSEM(:,:,ibin) = memissions
 enddo
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine SeasaltEmission.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SeasaltEmission.')
 endif


 if(allocated(fgridefficiency)) deallocate(fgridefficiency)
 if(allocated(fsstemis)       ) deallocate(fsstemis       )
 if(allocated(gweibull)       ) deallocate(gweibull       )
 if(allocated(fhoppel)        ) deallocate(fhoppel        )
 if(allocated(memissions)     ) deallocate(memissions     )
 if(allocated(nemissions)     ) deallocate(nemissions     )
 if(allocated(dqa)            ) deallocate(dqa            )


 call mpas_log_write('--- end subroutine emissions_SS2G_GridComp.')

 end subroutine emissions_SS2G_GridComp

!==================================================================================================================
 subroutine processes_SS2G_GridComp(self_params,self,to_MYNN,its,ite,jts,jte,kts,kte)
!==================================================================================================================
!--- input arguments:
 logical,intent(in):: to_MYNN
 integer,intent(in):: its,ite,jts,jte,kts,kte
 class(SS2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(SS2G_State),intent(inout):: self

!--- local variables and arrays:
 logical:: KIN

 integer:: i,i1,i2,j,j1,j2,k,km,ibin
 integer:: istat
 integer:: n_profile,n_vertint

 real(kind=RKIND):: fwet
!real(kind=RKIND),dimension(:,:,:),pointer:: rh20,rh80
 real(kind=RKIND),dimension(:,:),allocatable:: drydepf,dqa

 real(kind=RKIND),allocatable,dimension(:,:,:),target:: rh20,rh80

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_SS2G_GridComp:')

!--- SS2G settling:
!call mpas_log_write('--- enter subroutine Chem_Settling:')
 do ibin = 1, self_params%nbins
    if(associated(self%sssd)) self%sssd(:,:,ibin) = 0._RKIND
    istat = 0
    call Chem_Settling( &
              km        = self_params%km                 , &
              klid      = self_params%klid               , &
              bin       = ibin                           , &
              flag      = self_params%rhFlag             , &
              cdt       = self_params%cdt                , &
              grav      = grav                           , &
              radiusInp = self_params%radius(ibin)*1.e-6 , &
              rhopInp   = self_params%rhop(ibin)         , &
              int_qa    = self%SS(:,:,:,ibin)            , &
              tmpu      = self%t                         , &
              rhoa      = self%airdens                   , &
              rh        = self%rh2                       , &
              hghte     = self%zle                       , &
              delp      = self%delp                      , &
              fluxout   = self%sssd                      , &
              rc        = istat                            &
                      )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine Chem_Settling.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Chem_Settling.')
 endif


!--- SS2G dry deposition:
!call mpas_log_write('--- enter subroutine DryDeposition:')
 if(associated(self%ssvdep)) self%ssvdep(:,:) = 0._RKIND
 if(.not.allocated(dqa)    ) allocate(dqa(its:ite,jts:jte)    )
 if(.not.allocated(drydepf)) allocate(drydepf(its:ite,jts:jte))
 drydepf(:,:) = 0._RKIND
 istat = 0
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

!increase drydeposition velocity over land: 
 where(abs(self%lwi - LAND) < 0.5)
    drydepf = 5._RKIND*drydepf
 end where

!--- if gocart2G does not interact with the MYNN PBL parameterization, then we update the mixing ratios
!    due to dry deposition in the model layer adjacent to the surface, otherwise we simply skip this step:
 if(associated(self%ssvdep)) self%ssvdep(:,:) = drydepf(:,:)*self%delz(:,:,self_params%km)
 if(.not. to_MYNN) then
    do ibin = 1, self_params%nbins
       if(associated(self%ssdp)) self%ssdp(:,:,ibin) = 0._RKIND
       dqa = 0._RKIND
       dqa = max(0._RKIND,self%ss(:,:,self_params%km,ibin)*(1.-exp(-drydepf*self_params%cdt)))
       self%ss(:,:,self_params%km,ibin) = self%ss(:,:,self_params%km,ibin) - dqa
       if(associated(self%ssdp)) then
          self%ssdp(:,:,ibin) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt
       end if
    enddo
 endif
 if(istat /=0) then
    call mpas_log_write('--- SS2G_bc_GridComp: error in subroutine DryDeposition.', &
                        messageType=MPAS_LOG_CRIT)
 else
    if(allocated(dqa)    ) deallocate(dqa    )
    if(allocated(drydepf)) deallocate(drydepf)
!   call mpas_log_write('--- end subroutine DryDeposition.')
 endif


!--- SS2G large-scale wet removal:
!call mpas_log_write('--- enter subroutine WetRemovalGOCART2G:')
 do ibin = 1, self_params%nbins
    if(associated(self%sswt)) self%sswt(:,:,ibin) = 0._RKIND
    KIN   = .true.
    fwet  = 0.8_RKIND
    istat = 0
    call WetRemovalGOCART2G( &
              km        = self_params%km      , &
              klid      = self_params%klid    , &
              n1        = self_params%nbins   , &
              n2        = self_params%nbins   , &
              bin_ind   = ibin                , &
              cdt       = self_params%cdt     , &
              aero_type = 'DUST'              , &
              kin       = KIN                 , &
              grav      = grav                , &
              fwet      = fwet                , &
              aerosol   = self%SS(:,:,:,ibin) , &
              ple       = self%ple            , &
              tmpu      = self%t              , &
              rhoa      = self%airdens        , &
              pfllsan   = self%pfl_lsan       , &
              pfilsan   = self%pfi_lsan       , &
              precc     = self%cn_prcp        , &
              precl     = self%ncn_prcp       , &
              fluxout   = self%sswt           , &
              rc        = istat                 &
                        )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine WetRemovalGOCART2Gi.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine WetRemovalGOCART2G.')
 endif


!--- SS2G diagnostics:
 n_profile = size(self_params%wavelengths_profile)
 n_vertint = size(self_params%wavelengths_vertint)
!call mpas_log_write('--- enter subroutine Aero_Compute_Diags:')
!call mpas_log_write('--- nbins     = $i',intArgs=(/nbins/))
!call mpas_log_write('--- n_profile = $i',intArgs=(/n_profile/))
!call mpas_log_write('--- n_vertint = $i',intArgs=(/n_vertint/))
 if(associated(self%sssmass)   ) self%sssmass(:,:)       = 0._RKIND
 if(associated(self%sscmass)   ) self%sscmass(:,:)       = 0._RKIND
 if(associated(self%ssmass )   ) self%ssmass(:,:,:)      = 0._RKIND
 if(associated(self%ssexttau)  ) self%ssexttau(:,:,:)    = 0._RKIND
 if(associated(self%ssstexttau)) self%ssstexttau(:,:,:)  = 0._RKIND
 if(associated(self%ssscatau)  ) self%ssscatau(:,:,:)    = 0._RKIND
 if(associated(self%ssstscatau)) self%ssstscatau(:,:,:)  = 0._RKIND
 if(associated(self%ssfluxu)   ) self%ssfluxu(:,:)       = 0._RKIND
 if(associated(self%ssfluxv)   ) self%ssfluxv(:,:)       = 0._RKIND
 if(associated(self%ssconc)    ) self%ssconc(:,:,:)      = 0._RKIND
 if(associated(self%ssextcoef) ) self%ssextcoef(:,:,:,:) = 0._RKIND
 if(associated(self%ssscacoef) ) self%ssscacoef(:,:,:,:) = 0._RKIND
 if(associated(self%ssbckcoef) ) self%ssbckcoef(:,:,:,:) = 0._RKIND
 if(associated(self%ssangstr)  ) self%ssangstr(:,:)      = 0._RKIND
 if(associated(self%ssaeridx)  ) self%ssaeridx(:,:)      = 0._RKIND
 if(associated(self%sssmass25) ) self%sssmass25(:,:)     = 0._RKIND
 if(associated(self%sscmass25) ) self%sscmass25(:,:)     = 0._RKIND
 if(associated(self%ssmass25)  ) self%ssmass25(:,:,:)    = 0._RKIND
 if(associated(self%ssextt25)  ) self%ssextt25(:,:,:)    = 0._RKIND
 if(associated(self%ssscat25)  ) self%ssscat25(:,:,:)    = 0._RKIND

 istat = 0
 call Aero_Compute_Diags( &
              mie                 = self_params%diag_Mie                   , &
              km                  = self_params%km                         , &
              klid                = self_params%klid                       , &
              rlow                = self_params%rlow                       , &
              rup                 = self_params%rup                        , &
              nbegin              = 1                                      , &
              nbins               = self_params%nbins                      , &
              wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
              wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
              aerosol             = self%ss                                , &
              grav                = grav                                   , &
              tmpu                = self%t                                 , &
              rhoa                = self%airdens                           , &
              rh                  = self%rh2                               , &
              u                   = self%u                                 , &
              v                   = self%v                                 , &
              delp                = self%delp                              , &
              ple                 = self%ple                               , &
              tropp               = self%tropp                             , &
              sfcmass             = self%sssmass                           , &
              colmass             = self%sscmass                           , &
              mass                = self%ssmass                            , &
              exttau              = self%ssexttau                          , &
              scatau              = self%ssscatau                          , &
!             stexttau            = self%ssstexttau                        , &
!             stscatau            = self%ssstscatau                        , &
              fluxu               = self%ssfluxu                           , &
              fluxv               = self%ssfluxv                           , &
              conc                = self%ssconc                            , &
              extcoef             = self%ssextcoef                         , &
              scacoef             = self%ssscacoef                         , &
              bckcoef             = self%ssbckcoef                         , &
              angstrom            = self%ssangstr                          , &
              aerindx             = self%ssaeridx                          , &
              sfcmass25           = self%sssmass25                         , &
              colmass25           = self%sscmass25                         , &
              mass25              = self%ssmass25                          , &
              exttau25            = self%ssextt25                          , &
              scatau25            = self%ssscat25                          , &
              NO3nFlag            = .false.                                , &
              rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine Aero_Compute_Diags.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags.')
 endif


 i1 = lbound(self%rh2,1); i2 = ubound(self%rh2,1)
 j1 = lbound(self%rh2,2); j2 = ubound(self%rh2,2)
 km = ubound(self%rh2,3)


!call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH20:')
 if(associated(self%ssextcoefrh20)) self%ssextcoefrh20(:,:,:,:) = 0._RKIND
 if(associated(self%ssscacoefrh20)) self%ssscacoefrh20(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh20)) allocate(rh20(i1:i2,j1:j2,km))
 rh20(:,:,:) = 0.20
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = self_params%nbins                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = self%ss                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh20                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%ssextcoefrh20                     , &
           scacoef             = self%ssscacoefrh20                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- SS2G_GridComp: error in subroutine Aero_Compute_Diags RH20.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags RH20.')
 endif


!call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH80:')
 if(associated(self%ssextcoefrh80)) self%ssextcoefrh80(:,:,:,:) = 0._RKIND
 if(associated(self%ssscacoefrh80)) self%ssscacoefrh80(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh80)) allocate(rh80(i1:i2,j1:j2,km))
 rh80(:,:,:) = 0.80
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = self_params%nbins                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = self%ss                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh80                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%ssextcoefrh80                     , &
           scacoef             = self%ssscacoefrh80                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine Aero_Compute_Diags RH80.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine Aero_Compute_Diags RH80.')
 endif
 if(allocated(rh20)) deallocate(rh20)
 if(allocated(rh80)) deallocate(rh80)


 call mpas_log_write('--- end subroutine processes_SS2G_GridComp.')

 end subroutine processes_SS2G_GridComp

!==================================================================================================================
 subroutine rrtmg_SS2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte,nbndlw,nbndsw)
!==================================================================================================================
!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nbndlw,nbndsw
 class(SS2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(SS2G_State),intent(inout):: self

!--- local variables:
 integer:: nbands,km,nbins
 integer:: istat
 integer:: i,k,kk,j,n,nl,ns

 real(kind=RKIND),dimension(:,:,:),allocatable :: asy_s,ext_s,ssa_s
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qss5G

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_SS2G_GridComp:')

 km     = self_params%km
 nbins  = self_params%nbins
 nbands = self_params%rad_Mie%nch
!call mpas_log_write('--- km     = $i',intArgs=(/km/))
!call mpas_log_write('--- nbands = $i',intArgs=(/nbands/))

 if(associated(self%sstau_sw)) self%sstau_sw(:,:,:,:) = 0._RKIND
 if(associated(self%ssasy_sw)) self%ssasy_sw(:,:,:,:) = 0._RKIND
 if(associated(self%ssssa_sw)) self%ssssa_sw(:,:,:,:) = 0._RKIND
 if(associated(self%sstau_lw)) self%sstau_lw(:,:,:,:) = 0._RKIND
 if(associated(self%ssasy_lw)) self%ssasy_lw(:,:,:,:) = 0._RKIND
 if(associated(self%ssssa_lw)) self%ssssa_lw(:,:,:,:) = 0._RKIND

 if(.not.allocated(asy_s)) allocate(asy_s(its:ite,jts:jte,kts:kte))
 if(.not.allocated(ext_s)) allocate(ext_s(its:ite,jts:jte,kts:kte))
 if(.not.allocated(ssa_s)) allocate(ssa_s(its:ite,jts:jte,kts:kte))
 if(.not.allocated(qss5G)) allocate(qss5G(its:ite,jts:jte,kts:kte,nbins))

 asy_s(:,:,:) = 0._RKIND
 ext_s(:,:,:) = 0._RKIND
 ssa_s(:,:,:) = 0._RKIND

 do n = 1,nbins
    qss5G(:,:,:,n) = self%ss(:,:,:,n)*self%delp(:,:,:)/grav
 enddo


 nl = 0
 ns = 0
 do n = 1, nbands
    istat = 0
    asy_s(:,:,:) = 0._RKIND
    ext_s(:,:,:) = 0._RKIND
    ssa_s(:,:,:) = 0._RKIND
    call mie_(self_params%rad_Mie,its,ite,jts,jte,kts,kte,nbins,n,qss5G,self%rh2,ext_s,ssa_s,asy_s,istat)
    if(istat /=0) then
       call mpas_log_write('--- SS2G_GridComp: error in subroutine rrtmg_SS2G_GridComp.', &
                           messageType=MPAS_LOG_CRIT)
    else
       if(n .le. nbndsw) then
          ns = ns+1
          self%sstau_sw(:,:,:,ns) = real(ext_s(:,:,:),kind=RKIND)
          self%ssasy_sw(:,:,:,ns) = real(asy_s(:,:,:),kind=RKIND)
          self%ssssa_sw(:,:,:,ns) = real(ssa_s(:,:,:),kind=RKIND)
       elseif(n .gt. nbndsw) then
          nl = nl+1
          self%sstau_lw(:,:,:,nl) = real(ext_s(:,:,:),kind=RKIND)
          self%ssasy_lw(:,:,:,nl) = real(asy_s(:,:,:),kind=RKIND)
          self%ssssa_lw(:,:,:,nl) = real(ssa_s(:,:,:),kind=RKIND)
       endif
    endif
 enddo


 if(allocated(asy_s)) deallocate(asy_s)
 if(allocated(ext_s)) deallocate(ext_s)
 if(allocated(ssa_s)) deallocate(ssa_s)
 if(allocated(qss5G)) deallocate(qss5G)

 call mpas_log_write('--- end subroutine rrtmg_SS2G_GridComp.')

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

 end subroutine rrtmg_SS2G_GridComp

!==================================================================================================================
 end module SS2G_GridCompMod
!==================================================================================================================
