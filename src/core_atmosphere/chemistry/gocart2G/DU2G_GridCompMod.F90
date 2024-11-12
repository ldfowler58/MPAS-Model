!==================================================================================================================
 module DU2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: MPAS_LOG_CRIT
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: DustAerosolDistributionKok,Chem_UtilResVal,Chem_Settling,DryDeposition, &
                            WetRemovalGOCART2G,Aero_Compute_Diags,UpdateAerosolState

 use DU2G_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight,          &
                         fnum,rhFlag,pressure_lid_in_hPa,emission_scheme,clayFlag,soil_moisture_factor,  &
                         soil_clay_factor,uts_gamma,alpha,gamma,vertical_to_horizontal_flux_ratio_limit, &
                         radius_lower,radius_upper,source_fraction,ipoint,Ch_DU
 use DU2G_StateSpecs,only: DU2G_State


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

 real(kind=RKIND),parameter:: OCEAN   = 0._RKIND
 real(kind=RKIND),parameter:: LAND    = 1._RKIND
 real(kind=RKIND),parameter:: SEA_ICE = 2._RKIND

!--- types needed to define DU2G:
 type:: ThreadWorkspace
    integer:: day_save = -1
    integer:: nPts = -1
    integer,dimension(:),allocatable:: pstart,pend
    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: DU2G_GridComp
    character(len=255):: point_emissions_srcfilen    ! filename for pointwise emissions
    character(len=:),allocatable:: emission_scheme   ! emission scheme selector

    logical:: doing_point_emissions = .false.
    logical:: maringFlag=.false.                     ! maring settling velocity correction

    integer:: clayFlag                               ! clay and silt term in K14
    integer,dimension(:),allocatable:: ipoint        ! added parameter from WRF model.

    real(kind=RKIND):: f_swc                         ! soil mosture scaling factor
    real(kind=RKIND):: f_scl                         ! clay content scaling factor
    real(kind=RKIND):: uts_gamma                     ! threshold friction velocity parameter 'gamma'
    real(kind=RKIND):: alpha                         ! FENGSHA scaling factor
    real(kind=RKIND):: gamma                         ! FENGSHA tuning exponent
    real(kind=RKIND):: kvhmax                        ! FENGSHA max. vertical/horizontal mass flux ratio [1]
    real(kind=RKIND):: Ch_DU                         ! dust emission tuning coefficient [kg s2 m-5].
    real(kind=RKIND),dimension(NHRES):: Ch_DU_res    ! resolutions used for Ch_DU

    real(kind=RKIND),dimension(:),allocatable:: rlow ! particle effective radius lower bound [um]
    real(kind=RKIND),dimension(:),allocatable:: rup  ! particle effective radius upper bound [um]
    real(kind=RKIND),dimension(:),allocatable:: sfra ! fraction of total source
    real(kind=RKIND),dimension(:),allocatable:: sdist! FENGSHA aerosol fractional size distribution [1]

    !workspace for point emissions:
    type(ThreadWorkspace),allocatable :: workspaces(:)

    contains
       procedure:: emissions_GridComp => emissions_DU2G_GridComp
       procedure:: load_GridComp      => load_DU2G_GridComp
       procedure:: processes_GridComp => processes_DU2G_GridComp
 end type DU2G_GridComp

 type wrap_
    type(DU2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!==================================================================================================================
 subroutine load_DU2G_GridComp(self,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(DU2G_GridComp),intent(inout) :: self

!local variables:
 integer:: n,stat
 integer,dimension(3):: dims

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine load_DU2G_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization of variables in GA_Environment:
 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo


!--- initialization of all other variables in DU2_GridComp:
 self%emission_scheme = trim(emission_scheme)

 self%clayFlag  = clayFlag
 self%f_swc     = soil_moisture_factor
 self%f_scl     = soil_clay_factor
 self%uts_gamma = uts_gamma
 self%alpha     = alpha
 self%gamma     = gamma
 self%kvhmax    = vertical_to_horizontal_flux_ratio_limit

 if(.not.allocated(self%rlow)  ) allocate(self%rlow(self%nbins)  )
 if(.not.allocated(self%rup)   ) allocate(self%rup(self%nbins)   )
 if(.not.allocated(self%sfra)  ) allocate(self%sfra(self%nbins)  )
 if(.not.allocated(self%sdist) ) allocate(self%sdist(self%nbins) )
 if(.not.allocated(self%ipoint)) allocate(self%ipoint(self%nbins))

 do n = 1,self%nbins
    self%rlow(n)   = radius_lower(n)
    self%rup(n)    = radius_upper(n)
    self%sfra(n)   = source_fraction(n)
    self%ipoint(n) = ipoint(n)
 enddo
 call DustAerosolDistributionKok(self%radius,self%rup,self%rlow,self%sdist)


!--- initialization of dust emission tuning coefficient (Ch_DU) and resolutions used for Ch_DU (Ch_DU_res). 
!    dims(1) and dims(2) are the number of grid-points in the longitude and latitude directions. In MPAS,
!    we will need to adjust subroutine Chem_UtilResVal since MPAS uses an unstructured grid. For now, we
!    set self%Ch_DU to the same value as the one used in WRF-chem:
!do n = 1, NHRES
!   self%Ch_DU_res(n) = Ch_DU(n)
!enddo
!self%Ch_DU = Chem_UtilResVal(dims(1),dims(2),self%Ch_DU_res(:),stat)
 self%Ch_DU = 0.8e-09


 call mpas_log_write('--- end subroutine load_DU2G_GridCOMP:')

 end subroutine load_DU2G_GridComp

!==================================================================================================================
 subroutine emissions_DU2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte,nerod)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: nerod
 class(DU2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(DU2G_State),intent(inout):: self

!--- local variables and arrays:
 type(ThreadWorkspace),pointer:: emisPoints

 integer:: i,j,k,n
 integer:: istat
 integer:: nPts

 real(kind=RKIND),dimension(:,:),allocatable:: ustar_,ustar_t_,ustar_ts_
 real(kind=RKIND),dimension(:,:),allocatable:: R_,H_w_,z_
 real(kind=RKIND),dimension(:,:),allocatable:: f_erod_

 real(kind=RKIND),dimension(:,:,:),allocatable:: emissions_point
 real(kind=RKIND),dimension(:,:,:),allocatable:: emissions_surface
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: emissions

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_DU2G_GridComp: nerod = $i',intArgs=(/nerod/))


!--- total emissions:
 if(.not.allocated(emissions)) allocate(emissions(its:ite,jts:jte,kts:kte,self_params%nbins))
 emissions(:,:,:,:) = 0._RKIND


!--- point emissions: for now, the sourcecode does not support point emissions. therefore, we simply set the total
!    number of point emissions to 0 and allocate emissions_point accordingly.
 nPts = 0
 if(.not.allocated(emissions_point)) then
    if(nPts > 0) then
       allocate(emissions_point(its:ite,jts:jte,kts:kte))
       emissions_point(:,:,:) = 0._RKIND
    else
       allocate(emissions_point(0,0,0))
    endif
 endif


!--- surface emissions: for now, the sourcecode supports only the option ginoux (Ginoux et al. 2001). additional
!    options (K14, Fengsha) will be tested after the flow of the sourcecode is fully implemented.
 if(.not.allocated(emissions_surface)) allocate(emissions_surface(its:ite,jts:jte,self_params%nbins))

 select case(self_params%emission_scheme)
    case('ginoux')
       call DustEmissionGOCART2G_revised(               &
          radius    = self_params%radius*1.e-6,         &
          rhop      = self_params%rhop,                 &
          sfra      = self_params%sfra,                 &
          Ch_DU     = self_params%Ch_Du,                &
          ipoint    = self_params%ipoint,               &
          grav      = grav,                             &
          du_src    = self%du_src,                      &
          oro       = self%lwi,                         &
          frlake    = self%frlake,                      &
          u10m      = self%u10m,                        &
          v10m      = self%v10m,                        &
          gwet      = self%wet1,                        &
          dz        = self%delz(:,:,self_params%km),    &
          rhoa      = self%airdens(:,:,self_params%km), &
          u         = self%u(:,:,self_params%km),       &
          v         = self%v(:,:,self_params%km),       &
          emissions = emissions_surface                 &
                                        )

    case default
 end select


!--- updates mineral dust mixing ratios after calculation of surface emissions:
 call UpdateAerosolState( &
    nbins             = self_params%nbins, &
    km                = self_params%km,    &
    cdt               = self_params%cdt,   &
    sfrac             = self_params%sfra,  &
    nPts              = nPts,              &
    grav              = grav,              &
    emissions         = emissions,         &
    emissions_surface = emissions_surface, &
    emissions_point   = emissions_point,   &
    delp              = self%delp,         &
    aero              = self%du,           &
    rc                = istat              &
                        )


!--- saves surface emissions to diagnostic DUEM:
 do n = 1,self_params%nbins
    self%duem(:,:,n) = emissions(:,:,self_params%km,n)
 enddo
  

 call mpas_log_write('--- end subroutine emissions_DU2G_GridComp:')

 end subroutine emissions_DU2G_GridComp

!==================================================================================================================
 subroutine processes_DU2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 class(DU2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(DU2G_State),intent(inout):: self

!--- local variables:
 logical:: KIN

 integer:: i,i1,i2,j,j1,j2,k,km,ibin
 integer:: istat
 integer:: n_profile,n_vertint

 real(kind=RKIND):: fwet
 real(kind=RKIND),dimension(:,:,:),pointer:: rh20,rh80
 real(kind=RKIND),dimension(:,:),allocatable:: drydepf,dqa

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_DU2G_GridComp: nbins = $i',intArgs=(/self_params%nbins/))

!--- DU2G settling:
 call mpas_log_write('--- enter subroutine Chem_Settling:')
 do ibin = 1, self_params%nbins
    if(associated(self%dusd)) self%dusd(:,:,ibin) = 0._RKIND
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
              int_qa    = self%du(:,:,:,ibin)            , &
              tmpu      = self%t                         , &
              rhoa      = self%airdens                   , &
              rh        = self%rh2                       , &
              hghte     = self%zle                       , &
              delp      = self%delp                      , &
              fluxout   = self%dusd                      , &
              rc        = istat                            &
                      )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- DU2G_GridComp: error in subroutine Chem_Settling', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Chem_Settling:')
 endif


!--- DU2G dry deposition:
 call mpas_log_write('--- enter subroutine DryDeposition:')
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
 do ibin = 1, self_params%nbins
    if(associated(self%dudp)) self%dudp(:,:,ibin) = 0._RKIND
    dqa = 0._RKIND
    dqa = max(0._RKIND,self%du(:,:,self_params%km,ibin)*(1.-exp(-drydepf*self_params%cdt)))
    self%du(:,:,self_params%km,ibin) = self%du(:,:,self_params%km,ibin) - dqa
    if(associated(self%dudp)) then
       self%dudp(:,:,ibin) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt
    end if
 enddo
 if(istat /=0) then
    call mpas_log_write('--- DU2G_bc_GridComp: error in subroutine DryDeposition', &
                        messageType=MPAS_LOG_CRIT)
 else
    if(allocated(dqa)    ) deallocate(dqa    )
    if(allocated(drydepf)) deallocate(drydepf)
    call mpas_log_write('--- end subroutine DryDeposition:')
 endif


!--- DU2G large-scale wet removal:
 call mpas_log_write('--- enter subroutine WetRemovalGOCART2G:')
 do ibin = 1, self_params%nbins
    if(associated(self%duwt)) self%duwt(:,:,ibin) = 0._RKIND
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
              aerosol   = self%du(:,:,:,ibin) , &
              ple       = self%ple            , &
              tmpu      = self%t              , &
              rhoa      = self%airdens        , &
              pfllsan   = self%pfl_lsan       , &
              pfilsan   = self%pfi_lsan       , &
              precc     = self%cn_prcp        , &
              precl     = self%ncn_prcp       , &
              fluxout   = self%duwt           , &
              rc        = istat                 &
                        )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- DU2G_GridComp: error in subroutine WetRemovalGOCART2G', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine WetRemovalGOCART2G:')
 endif


!--- DU2G  diagnostics:
 n_profile = size(self_params%wavelengths_profile)
 n_vertint = size(self_params%wavelengths_vertint)
 call mpas_log_write('--- enter subroutine Aero_Compute_Diags:')
 call mpas_log_write('--- nbins     = $i',intArgs=(/nbins/))
 call mpas_log_write('--- n_profile = $i',intArgs=(/n_profile/))
 call mpas_log_write('--- n_vertint = $i',intArgs=(/n_vertint/))
 if(associated(self%dusmass)   ) self%dusmass(:,:)       = 0._RKIND
 if(associated(self%ducmass)   ) self%ducmass(:,:)       = 0._RKIND
 if(associated(self%dumass )   ) self%dumass(:,:,:)      = 0._RKIND
 if(associated(self%duexttau)  ) self%duexttau(:,:,:)    = 0._RKIND
 if(associated(self%dustexttau)) self%dustexttau(:,:,:)  = 0._RKIND
 if(associated(self%duscatau)  ) self%duscatau(:,:,:)    = 0._RKIND
 if(associated(self%dustscatau)) self%dustscatau(:,:,:)  = 0._RKIND
 if(associated(self%dufluxu)   ) self%dufluxu(:,:)       = 0._RKIND
 if(associated(self%dufluxv)   ) self%dufluxv(:,:)       = 0._RKIND
 if(associated(self%duconc)    ) self%duconc(:,:,:)      = 0._RKIND
 if(associated(self%duextcoef) ) self%duextcoef(:,:,:,:) = 0._RKIND
 if(associated(self%duscacoef) ) self%duscacoef(:,:,:,:) = 0._RKIND
 if(associated(self%dubckcoef) ) self%dubckcoef(:,:,:,:) = 0._RKIND
 if(associated(self%duangstr)  ) self%duangstr(:,:)      = 0._RKIND
 if(associated(self%duaeridx)  ) self%duaeridx(:,:)      = 0._RKIND
 istat = 0
 call Aero_Compute_Diags( &
              mie                 = self_params%diag_Mie                   , &
              km                  = self_params%km                         , &
              klid                = self_params%klid                       , &
              nbegin              = 1                                      , &
              nbins               = 2                                      , &
              wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
              wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
              aerosol             = self%du                                , &
              grav                = grav                                   , &
              tmpu                = self%t                                 , &
              rhoa                = self%airdens                           , &
              rh                  = self%rh2                               , &
              u                   = self%u                                 , &
              v                   = self%v                                 , &
              delp                = self%delp                              , &
              ple                 = self%ple                               , &
              tropp               = self%tropp                             , &
              sfcmass             = self%dusmass                           , &
              colmass             = self%ducmass                           , &
              mass                = self%dumass                            , &
              exttau              = self%duexttau                          , &
              scatau              = self%duscatau                          , &
!             stexttau            = self%dustexttau                        , &
!             stscatau            = self%dustscatau                        , &
              fluxu               = self%dufluxu                           , &
              fluxv               = self%dufluxv                           , &
              conc                = self%duconc                            , &
              extcoef             = self%duextcoef                         , &
              scacoef             = self%duscacoef                         , &
              bckcoef             = self%dubckcoef                         , &
              angstrom            = self%duangstr                          , &
              aerindx             = self%duaeridx                          , &
              NO3nFlag            = .false.                                , &
              rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- DU2G_GridComp: error in subroutine Aero_Compute_Diags', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags:')
 endif


 i1 = lbound(self%rh2,1); i2 = ubound(self%rh2,1)
 j1 = lbound(self%rh2,2); j2 = ubound(self%rh2,2)
 km = ubound(self%rh2,3)

 call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH20:')
 if(.not.associated(rh20)) allocate(rh20(i1:i2,j1:j2,km))
 if(associated(self%duextcoefrh20)) self%duextcoefrh20(:,:,:,:) = 0._RKIND
 if(associated(self%duscacoefrh20)) self%duscacoefrh20(:,:,:,:) = 0._RKIND
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
           aerosol             = self%du                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh20                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%duextcoefrh20                     , &
           scacoef             = self%duscacoefrh20                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- DU2G_GridComp: error in subroutine Aero_Compute_Diags RH20', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags RH20:')
 endif


 call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH80:')
 if(.not.associated(rh80)) allocate(rh80(i1:i2,j1:j2,km))
 if(associated(self%duextcoefrh80)) self%duextcoefrh80(:,:,:,:) = 0._RKIND
 if(associated(self%duscacoefrh80)) self%duscacoefrh80(:,:,:,:) = 0._RKIND
 rh80(:,:,:) = 0.80
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 2                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = self%du                                , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh80                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%duextcoefrh80                     , &
           scacoef             = self%duscacoefrh80                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- DU2G_GridComp: error in subroutine Aero_Compute_Diags RH80', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags RH80:')
 endif
 if(associated(rh20)) deallocate(rh20)
 if(associated(rh80)) deallocate(rh80)


 call mpas_log_write('--- end subroutine processes_DU2G_GridComp:')

 end subroutine processes_DU2G_GridComp

!==================================================================================================================
 subroutine DustEmissionGOCART2G_revised(radius,rhop,sfra,grav,oro,frlake,dz,gwet,u10m,v10m, &
                                         rhoa,u,v,Ch_DU,ipoint,du_src,emissions)
!==================================================================================================================

!--- input arguments:
 integer,intent(in),dimension(:):: ipoint

 real(kind=RKIND),intent(in):: grav
 real(kind=RKIND),intent(in):: Ch_DU
 real(kind=RKIND),intent(in),dimension(:):: radius
 real(kind=RKIND),intent(in),dimension(:):: rhop
 real(kind=RKIND),intent(in),dimension(:):: sfra
 real(kind=RKIND),intent(in),dimension(:,:):: frlake
 real(kind=RKIND),intent(in),dimension(:,:):: dz,gwet,oro
 real(kind=RKIND),intent(in),dimension(:,:):: u10m,v10m
 real(kind=RKIND),intent(in),dimension(:,:):: rhoa,u,v
 real(kind=RKIND),intent(in),dimension(:,:,:):: du_src

!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(:,:,:):: emissions

!--- local variables and arrays:
 integer:: i,i1,i2,j,j1,j2,m,n
 integer:: nbins
 integer,dimension(2):: dims

 real(kind=RKIND),parameter:: LAND = 1._RKIND
 real(kind=RKIND):: u_thresh0,u_thresh
 real(kind=RKIND):: air_dens,soil_dens,diam,gwet1,w10m

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine DustEmissionGOCART2G_revised:')

!--- get dimensions:
 nbins = size(radius)

 dims = shape(u10m)
 i1 = 1 ; i2 = dims(1)
 j1 = 1 ; j2 = dims(2)

 do n = 1,nbins
    m = ipoint(n)
    diam = 2*radius(n)
    soil_dens = rhop(n)
!   call mpas_log_write('$i $i $i $i $r $r',intArgs=(/j2,i2,n,m/),realArgs=(/diam,soil_dens/))

    emissions(:,:,n) = 0._RKIND

    do j = j1,j2
       do i = i1,i2

          if(oro(i,j) /= LAND) cycle ! only compute emissions over land points:

          !--- compute wind speed:
          if(dz(i,j).lt.12.) then
             w10m = sqrt(u(i,j)*u(i,j) + v(i,j)*v(i,j))
          else
             w10m = sqrt(u10m(i,j)*u10m(i,j) + v10m(i,j)*v10m(i,j))
          endif

          !--- compute the threshold velocity of wind erosion for dry soil in each dust bin,
          !    following Marticorena et al. (1997):
          air_dens  = rhoa(i,j)

          u_thresh0 = 0.13 * sqrt(soil_dens*grav*diam/air_dens) &
                    * sqrt(1.+6.e-7/(soil_dens*grav*diam**2.5)) &
                    / sqrt(1.928*(1331.*(100.*diam)**1.56+0.38)**0.092 - 1.)
 
          !--- adjust threshold velocity as function of soil moisture, following Ginoux et al. (2001), and
          !    compute emissions:
          gwet1 = gwet(i,j)

          if(gwet1 .lt. 0.5) then
             u_thresh = max(0.,u_thresh0*(1.2+0.2*alog10(max(1.e-3,gwet1))))
             if(w10m .gt. u_thresh) then
                emissions(i,j,n) = (1.-frlake(i,j))*w10m**2*(w10m-u_thresh)
                emissions(i,j,n) = Ch_DU*du_src(i,j,m)*emissions(i,j,n)
!               call mpas_log_write('--- land: $i $i $r $r $r $r $r $r',intArgs=(/n,i/), &
!                         realArgs=(/w10m,air_dens,u_thresh0,u_thresh,gwet1,emissions(i,j,n)/))
             endif
          endif

       enddo
    enddo
!   call mpas_log_write(' ')
 enddo

 call mpas_log_write('--- end subroutine DustEmissionGOCART2G_revised:')

 end subroutine DustEmissionGOCART2G_revised

!==================================================================================================================
 end module DU2G_GridCompMod
!==================================================================================================================
