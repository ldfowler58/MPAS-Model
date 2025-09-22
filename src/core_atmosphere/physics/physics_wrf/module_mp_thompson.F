 module module_mp_thompson
 use mpas_log
 use mpas_kind_types,only: RKIND

 use mp_thompson,only: mp_thompson_run

 implicit none
 private
 public:: mp_thompson_driver


 contains


!==================================================================================================================
 subroutine mp_thompson_driver(qv,qc,qr,qi,qs,qg,ni,nr,nc,nwfa,nifa,nwfa2d,nifa2d,th,pii,p,w,dz,dt, &
                         ntc,muc,has_reqc,has_reqi,has_reqs,re_cloud,re_ice,re_snow,rainnc,rainncv, &
                         snownc,snowncv,graupelnc,graupelncv,sr,rainprod,evapprod,pflrain,pflice,   &
                         pflsnow,pflgraul,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte
 integer,intent(in):: has_reqc,has_reqi,has_reqs

 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in),dimension(ims:ime,jms:jme):: ntc,muc
 real(kind=RKIND),intent(inout),dimension(ims:ime,jms:jme),optional:: nifa2d,nwfa2d
 real(kind=RKIND),intent(in),dimension(ims:ime,kms:kme,jms:jme):: pii,p,w,dz

!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(ims:ime,jms:jme):: rainnc
 real(kind=RKIND),intent(inout),dimension(ims:ime,jms:jme),optional:: snownc,graupelnc

 real(kind=RKIND),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: th,qv,qc,qr,qi,qs,qg,ni,nr
 real(kind=RKIND),intent(inout),dimension(ims:ime,kms:kme,jms:jme):: re_cloud,re_ice,re_snow
 real(kind=RKIND),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: nc,nifa,nwfa

!--- output arguments:
 real(kind=RKIND),intent(out),dimension(ims:ime,jms:jme):: rainncv,sr
 real(kind=RKIND),intent(out),dimension(ims:ime,jms:jme),optional:: snowncv,graupelncv
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: rainprod,evapprod
 real(kind=RKIND),intent(out),dimension(ims:ime,kms:kme,jms:jme),optional:: pflrain,pflice,pflsnow,pflgraul

!--- local variables and arrays:
 integer:: i,j,k
 integer:: itimestep

 real(kind=RKIND),dimension(its:ite):: muc_hv,ntc_hv
 real(kind=RKIND),dimension(its:ite):: nifa2d_hv,nwfa2d_hv
 real(kind=RKIND),dimension(its:ite):: rainnc_hv,graupelnc_hv,snownc_hv,sr_hv
 real(kind=RKIND),dimension(its:ite):: rainncv_hv,graupelncv_hv,snowncv_hv
 real(kind=RKIND),dimension(its:ite,kts:kte):: th_hv,dz_hv,pi_hv,p_hv,w_hv
 real(kind=RKIND),dimension(its:ite,kts:kte):: qv_hv,qc_hv,qr_hv,qi_hv,qs_hv,qg_hv,ni_hv,nr_hv
 real(kind=RKIND),dimension(its:ite,kts:kte):: nc_hv,nwfa_hv,nifa_hv
 real(kind=RKIND),dimension(its:ite,kts:kte):: rainprod_hv,evapprod_hv
 real(kind=RKIND),dimension(its:ite,kts:kte):: pflrain_hv,pflice_hv,pflsnow_hv,pflgraul_hv
 real(kind=RKIND),dimension(its:ite,kts:kte):: recloud_hv,reice_hv,resnow_hv

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mp_thompson_driver:')

 itimestep = 0

 do j = jts,jte

    do i = its,ite
       !--- input arguments:
       muc_hv(i) = muc(i,j)
       ntc_hv(i) = ntc(i,j)

       do k = kts,kte
          pi_hv(i,k) = pii(i,k,j)
          p_hv(i,k)  = p(i,k,j)
          dz_hv(i,k) = dz(i,k,j)
          w_hv(i,k)  = w(i,k,j)
       enddo

       !--- inout arguments:
       rainnc_hv(i) = rainnc(i,j)

       do k = kts,kte
          th_hv(i,k) = th(i,k,j)
          qv_hv(i,k) = qv(i,k,j)
          qc_hv(i,k) = qc(i,k,j)
          qr_hv(i,k) = qr(i,k,j)
          qi_hv(i,k) = qi(i,k,j)
          qs_hv(i,k) = qs(i,k,j)
          qg_hv(i,k) = qg(i,k,j)
          ni_hv(i,k) = ni(i,k,j)
          nr_hv(i,k) = nr(i,k,j)

          recloud_hv(i,k) = re_cloud(i,k,j)
          reice_hv(i,k)   = re_ice(i,k,j)
          resnow_hv(i,k)  = re_snow(i,k,j)
       enddo
    enddo

    if(present(nc) .and. present(nwfa) .and. present(nifa)) then
       do i = its,ite
          do k = kts,kte
             nc_hv(i,k)   = nc(i,k,j)
             nwfa_hv(i,k) = nwfa(i,k,j)
             nifa_hv(i,k) = nifa(i,k,j)
          enddo
       enddo
    endif

    if(present(nwfa2d) .and. present(nifa2d)) then
       do i = its,ite
          nifa2d_hv(i) = nifa2d(i,j)
          nwfa2d_hv(i) = nwfa2d(i,j)
       enddo
    else
       do i = its,ite
          nifa2d_hv(i) = 0._RKIND
          nwfa2d_hv(i) = 0._RKIND
       enddo
    endif

    if(present(snownc)) then
       do i = its,ite
          snownc_hv(i) = snownc(i,j)
       enddo
    endif
    if(present(graupelnc)) then
       do i = its,ite
          graupelnc_hv(i) = graupelnc(i,j)
       enddo
    endif


    !--- BEGIN CLOUD MICROPHYSICS SCHEME:
    do i = its,10
       do k = kts,kte
          call mpas_log_write('$i $i $r',intArgs=(/i,k/),realArgs=(/p_hv(i,k)/))
       enddo
    enddo
    call mpas_log_write(' ')
    call mp_thompson_run( &
            qv         = qv_hv         , qc       = qc_hv     , qr       = qr_hv       , qi        = qi_hv        , &
            qs         = qs_hv         , qg       = qg_hv     , ni       = ni_hv       , nr        = nr_hv        , &
            nc         = nc_hv         , nwfa     = nwfa_hv   , nifa     = nifa_hv     , nwfa2d    = nwfa2d_hv    , &
            nifa2d     = nifa2d_hv     , th       = th_hv     , pi       = pi_hv       , p         = p_hv         , &
            w          = w_hv          , dz       = dz_hv     , dt_in    = dt          , rainnc    = rainnc_hv    , &
            rainncv    = rainncv_hv    , snownc   = snownc_hv , snowncv  = snowncv_hv  , graupelnc = graupelnc_hv , &
            graupelncv = graupelncv_hv , sr       = sr_hv     , rainprod = rainprod_hv , evapprod  = evapprod_hv  , &
            re_cloud   = recloud_hv    , re_ice   = reice_hv  , re_snow  = resnow_hv   , has_reqc  = has_reqc     , &
            has_reqi   = has_reqi      , has_reqs = has_reqs  , ntc      = ntc_hv      , muc       = muc_hv       , &
            pflrain    = pflrain_hv    , pflice   = pflice_hv , pflsnow  = pflsnow_hv  , pflgraul  = pflgraul_hv  , &
            itimestep  = itimestep     ,                                                                            &
            ims        = ims           , ime      = ime       , jms      = jms         , jme       = jme          , &
            kms        = kms           , kme      = kme       , its      = its         , ite       = ite          , &
            jts        = jts           , jte      = jte       , kts      = kts         , kte       = kte            &
                        )
    !--- END CLOUD MICROPHYSICS SCHEME:


    do i = its,ite
       rainnc(i,j)  = rainnc_hv(i)
       rainncv(i,j) = rainncv_hv(i)
       sr(i,j)      = sr_hv(i)

       do k = kts,kte
          th(i,k,j) = th_hv(i,k)
          qv(i,k,j) = qv_hv(i,k)
          qc(i,k,j) = qc_hv(i,k)
          qr(i,k,j) = qr_hv(i,k)
          qi(i,k,j) = qi_hv(i,k)
          qs(i,k,j) = qs_hv(i,k)
          qg(i,k,j) = qg_hv(i,k)
          ni(i,k,j) = ni_hv(i,k)
          nr(i,k,j) = nr_hv(i,k)
       enddo
    enddo

    if(present(nc) .and. present(nwfa) .and. present(nifa)) then
       do i = its,ite
          do k = kts,kte
             nc(i,k,j)   = nc_hv(i,k)
             nwfa(i,k,j) = nwfa_hv(i,k)
             nifa(i,k,j) = nifa_hv(i,k)
          enddo
       enddo
    endif

    if(present(snowncv)) then
       do i = its,ite
          snowncv(i,j) = snowncv_hv(i)
       enddo
    endif
    if(present(snownc)) then
       do i = its,ite
          snownc(i,j) = snownc_hv(i)
       enddo
    endif

    if(present(graupelncv)) then
       do i = its,ite
          graupelncv(i,j) = graupelncv_hv(i)
       enddo
    endif
    if(present(graupelnc)) then
       do i = its,ite
          graupelnc(i,j) = graupelnc_hv(i)
       enddo
    endif

    if(present(rainprod) .and. present(evapprod)) then
       do i = its,ite
          do k = kts,kte
             rainprod(i,k,j) = rainprod_hv(i,k)
             evapprod(i,k,j) = evapprod_hv(i,k)
          enddo
       enddo
    endif

    if(present(pflrain) .and. present(pflice) .and. present(pflsnow) .and. present(pflgraul)) then
       do i = its,ite
          do k = kts,kte
             pflrain(i,k,j)  = pflrain_hv(i,k)
             pflice(i,k,j)   = pflice_hv(i,k)
             pflsnow(i,k,j)  = pflsnow_hv(i,k)
             pflgraul(i,k,j) = pflgraul_hv(i,k)
          enddo
       enddo
    endif

 enddo

 call mpas_log_write('--- end subroutine mp_thompson_driver.')

 end subroutine mp_thompson_driver

!==================================================================================================================
 end module module_mp_thompson
!==================================================================================================================
