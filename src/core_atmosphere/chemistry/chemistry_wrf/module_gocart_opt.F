!=================================================================================================================
 module module_gocart_opt
 use mpas_kind_types
 use mpas_log

 use module_data_gocart_chem,only: nh4_mfac,oc_mfac
 use module_data_sorgam,only: dginia,dginic,dginin,sginia,sginic,sginin,hygro_cl_aer,hygro_dust_aer, &
                              hygro_nh4_aer,hygro_na_aer,hygro_oc_aer,hygro_oin_aer,hygro_so4_aer 
 use module_data_gocart_dust,only: ndust,den_dust,ra_dust,rb_dust
 use module_data_gocart_seas,only: ra_seas,rb_seas
 use module_data_mosaic_asect,only: hygro_msa_aer
 use module_data_rrtmg_aeropt
 use module_gocart_mie,only: mieaer

! n = 0
! mass_seas = 0.0
! do m =p_seas_1,  p_seas_3 ! loop over seasalt size bins less than 10 um diam
! n = n+1
! mass_seas=mass_seas+seasfrc_goc8bin(n,isize)*chem(i,k,j,m)     
! end do
! n = 0
! mass_soil = 0.0
! do m =p_dust_1,  p_dust_1+ndust-2 ! loop over dust size bins less than 10 um diam
! n = n+1
! mass_soil=mass_soil+dustfrc_goc8bin(n,isize)*chem(i,k,j,m)
! end do


 implicit none
 private
 public:: gocart_optical_averaging


 integer,parameter:: p_bc1    = 1
 integer,parameter:: p_bc2    = 2
 integer,parameter:: p_oc1    = 3
 integer,parameter:: p_oc2    = 4
 integer,parameter:: p_dust_1 = 5
 integer,parameter:: p_dust_2 = 6
 integer,parameter:: p_dust_3 = 7
 integer,parameter:: p_dust_4 = 8
 integer,parameter:: p_dust_5 = 9
 integer,parameter:: p_seas_1 = 10
 integer,parameter:: p_seas_2 = 11
 integer,parameter:: p_seas_3 = 12
 integer,parameter:: p_seas_4 = 13
 integer,parameter:: p_sulf   = 14
 integer,parameter:: p_msa    = 15
 integer,parameter:: p_p25    = 16


 contains


!=================================================================================================================
 subroutine gocart_optical_averaging(its,ite,kts,kte,rho,relhum,zgrid,qbc1,qbc2,qoc1,qoc2,qdust1,qdust2,qdust3, &
                          qdust4,qdust5,qseas1,qseas2,qseas3,qseas4,qmsa,qsulf,tauaersw,extaersw,gaersw,waersw, &
                          bscoefsw,l2aer,l3aer,l4aer,l5aer,l6aer,l7aer,tauaerlw,extaerlw)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: qbc1,qbc2,qoc1,qoc2,qdust1,qdust2,qdust3,qdust4, &
                                                          qdust5,qseas1,qseas2,qseas3,qseas4,qmsa,qsulf
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho,relhum
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte+1):: zgrid


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte,nswbands):: tauaersw,extaersw,gaersw,waersw,bscoefsw
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte,nswbands):: l2aer,l3aer,l4aer,l5aer,l6aer,l7aer
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte,nlwbands):: tauaerlw,extaerlw


!--- local variables and arrays:
 integer,parameter:: nspint = nswbands !number of spectral interval shortwave bands.
 integer,parameter:: nbin_o = 8
 integer,parameter:: uoc    = 0        !flag for UoC dust emissions.
 integer:: i,isize,k,nl,n,ns

 real(kind=RKIND),dimension(kts:kte):: dz_col
 real(kind=RKIND),dimension(its:ite,kts:kte):: alt,dz

 real(kind=RKIND),dimension(its:ite,kts:kte,1:nbin_o):: radius_core,radius_wet,number_bin

 complex,dimension(1:nbin_o,kts:kte):: swrefindx_col1,swrefindx_core_col1,swrefindx_shell_col1
 complex,dimension(1:nbin_o,kts:kte,1:nspint) :: swrefindx_col,swrefindx_core_col,swrefindx_shell_col
 complex,dimension(its:ite,kts:kte,1:nbin_o,1:nspint):: swrefindx,swrefindx_core,swrefindx_shell

 complex,dimension(1:nbin_o,kts:kte):: lwrefindx_col1,lwrefindx_core_col1,lwrefindx_shell_col1
 complex,dimension(1:nbin_o,kts:kte,1:nlwbands):: lwrefindx_col,lwrefindx_core_col,lwrefindx_shell_col
 complex,dimension(its:ite,kts:kte,1:nbin_o,1:nlwbands):: lwrefindx,lwrefindx_core,lwrefindx_shell

 real(kind=RKIND),dimension(1:nbin_o,kts:kte):: radius_wet_col,number_bin_col,radius_core_col
 real(kind=RKIND),dimension(nspint,kts:kte)  :: l2,l3,l4,l5,l6,l7
 real(kind=RKIND),dimension(nspint,kts:kte)  :: swsizeaer,swextaer,swwaer,swgaer,swtauaer,swbscoef
 real(kind=RKIND),dimension(nlwbands,kts:kte):: lwtauaer,lwextaer

 real(kind=R8KIND):: curr_secs

 integer,parameter:: num_chem = p_p25-p_bc1+1
 real(kind=RKIND),dimension(its:ite,kts:kte,num_chem):: chem

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter module_gocart_opt:')

 curr_secs = 0._R8KIND

!--- initialization of inverse of density and layer thickness:
 do k = kts,kte
    do i = its,ite
       alt(i,k) = 1._RKIND/rho(i,k)
       dz(i,k)  = zgrid(i,k+1) - zgrid(i,k)
    enddo
 enddo


!--- initialization of the chem array for use in optical_prep_gocart:
 do k = kts,kte
    do i = its,ite
       chem(i,k,p_bc1)    = qbc1(i,k)
       chem(i,k,p_bc2)    = qbc2(i,k)
       chem(i,k,p_oc1)    = qoc1(i,k)
       chem(i,k,p_oc2)    = qoc2(i,k)
!      chem(i,k,p_dust_1) = qdust1(i,k)
!      chem(i,k,p_dust_2) = qdust2(i,k)
!      chem(i,k,p_dust_3) = qdust3(i,k)
!      chem(i,k,p_dust_4) = qdust4(i,k)
!      chem(i,k,p_dust_5) = qdust5(i,k)
!      chem(i,k,p_seas_1) = qseas1(i,k)
!      chem(i,k,p_seas_2) = qseas2(i,k)
!      chem(i,k,p_seas_3) = qseas3(i,k)
!      chem(i,k,p_seas_4) = qseas4(i,k)
       chem(i,k,p_dust_1) = 0._RKIND
       chem(i,k,p_dust_2) = 0._RKIND
       chem(i,k,p_dust_3) = 0._RKIND
       chem(i,k,p_dust_4) = 0._RKIND
       chem(i,k,p_dust_5) = 0._RKIND
       chem(i,k,p_seas_1) = 0._RKIND
       chem(i,k,p_seas_2) = 0._RKIND
       chem(i,k,p_seas_3) = 0._RKIND
       chem(i,k,p_seas_4) = 0._RKIND
       chem(i,k,p_sulf)   = qsulf(i,k)
       chem(i,k,p_msa)    = qmsa(i,k)
       chem(i,k,p_p25)    = 0._RKIND
    enddo
 enddo


!--- conversion from kg/kg to ug/kg:
 do k = kts,kte
    do i = its,ite
       chem(i,k,p_bc1)    = chem(i,k,p_bc1)*1.e09
       chem(i,k,p_bc2)    = chem(i,k,p_bc2)*1.e09
       chem(i,k,p_oc1)    = chem(i,k,p_oc1)*1.e09
       chem(i,k,p_oc2)    = chem(i,k,p_oc2)*1.e09
       chem(i,k,p_dust_1) = chem(i,k,p_dust_1)*1.e09
       chem(i,k,p_dust_2) = chem(i,k,p_dust_2)*1.e09
       chem(i,k,p_dust_3) = chem(i,k,p_dust_3)*1.e09
       chem(i,k,p_dust_4) = chem(i,k,p_dust_4)*1.e09
       chem(i,k,p_dust_5) = chem(i,k,p_dust_5)*1.e09
       chem(i,k,p_seas_1) = chem(i,k,p_seas_1)*1.e09
       chem(i,k,p_seas_2) = chem(i,k,p_seas_2)*1.e09
       chem(i,k,p_seas_3) = chem(i,k,p_seas_3)*1.e09
       chem(i,k,p_seas_4) = chem(i,k,p_seas_4)*1.e09
       chem(i,k,p_p25)    = chem(i,k,p_p25)*1.e09
    enddo
 enddo


!--- call optical_prep_gocart:
 call optical_prep_gocart(nbin_o,num_chem,chem,alt,relhum,radius_core,radius_wet,number_bin, &
                          swrefindx,swrefindx_core,swrefindx_shell,lwrefindx,lwrefindx_core, &
                          lwrefindx_shell,uoc,its,ite,kts,kte)
!do n = 1,nbin_o
!   do i = its,ite
!      do k = kts,kte
!         call mpas_log_write('$i $i $i $r $r $r',intArgs=(/n,i,k/),realArgs=(/radius_core(i,k,n), &
!                             radius_wet(i,k,n),number_bin(i,k,n)/))
!      enddo
!      call mpas_log_write(' ')
!   enddo
!enddo
!write(*,*) swrefindx
!write(*,*) swrefindx_core
!write(*,*) swrefindx_shell
!write(*,*) lwrefindx
!write(*,*) lwrefindx_core
!write(*,*) lwrefindx_shell


 do i = its,ite

    do k = kts,kte
       do isize = 1,nbin_o
          number_bin_col(isize,k)        = number_bin(i,k,isize)
          radius_wet_col(isize,k)        = radius_wet(i,k,isize)
          swrefindx_col(isize,k,:)       = swrefindx(i,k,isize,:)
          lwrefindx_col(isize,k,:)       = lwrefindx(i,k,isize,:)
          radius_core_col(isize,k)       = radius_core(i,k,isize)
          swrefindx_core_col(isize,k,:)  = swrefindx_core(i,k,isize,:)
          swrefindx_shell_col(isize,k,:) = swrefindx_shell(i,k,isize,:)
          lwrefindx_core_col(isize,k,:)  = lwrefindx_core(i,k,isize,:)
          lwrefindx_shell_col(isize,k,:) = lwrefindx_shell(i,k,isize,:)

          radius_core_col(isize,k)       = 0.0                          !set to 0 if option_method <= 2 in WRF/chem.
          swrefindx_core_col(isize,k,:)  = cmplx(0.0,0.0)               !set to 0 if option_method <= 2 in WRF/chem.

          !--- at 600 nm:
          swrefindx_col1(isize,k)        = swrefindx(i,k,isize,3)
          swrefindx_core_col1(isize,k)   = swrefindx_core(i,k,isize,3)
          swrefindx_shell_col1(isize,k)  = swrefindx_shell(i,k,isize,3)
          swrefindx_core_col1(isize,k)   = cmplx(0.0,0.0)              !set to 0 if option_method <= 2 in WRF/chem.

          !--- following JB (see lines 345-362 of WRF/chem/module_optical_averaging.F):
          swrefindx_col(isize,k,:) = swrefindx_shell(i,k,isize,:) ! JCB
          swrefindx_col1(isize,k)  = swrefindx_shell(i,k,isize,3) ! JCB

          if(radius_wet_col(isize,k) < 1e-20) then
               radius_core_col(isize,k)=0.0
          else if(radius_core_col(isize,k)/radius_wet_col(isize,k)**3.le.0.0001) then
               radius_core_col(isize,k)=0.0
          end if
       enddo

       lwtauaer(:,k) = 1.e-20
       swtauaer(:,k) = 1.e-20

       dz_col(k) = dz(i,k)
    enddo


    !--- call to Mie scattering:
    call mpas_log_write('--- enter subroutine miaer - i = $i',intArgs=(/i/))
    call mieaer(1,i,1,nbin_o,number_bin_col,radius_wet_col,swrefindx_col,lwrefindx_col,  &
                dz,curr_secs,kts,kte,swsizeaer,swextaer,swwaer,swgaer,swtauaer,lwextaer, &
                lwtauaer,l2,l3,l4,l5,l6,l7,swbscoef)


    !--- shortwave optical properties:
    do k = kts,kte
       do ns = 1,nswbands
          tauaersw(i,k,ns) = amax1(swtauaer(ns,k),1.e-20)
          extaersw(i,k,ns) = amax1(swextaer(ns,k),1.e-20)
          gaersw(i,k,ns)   = amax1(amin1(swgaer(ns,k),1.0-1.e-8),1.e-20)
          waersw(i,k,ns)   = amax1(amin1(swwaer(ns,k),1.0-1.e-8),1.e-20)
          bscoefsw(i,k,ns) = amax1(swbscoef(ns,k),1.e-20)
       enddo
    enddo
    do k = kts,kte
       do ns = 1,nswbands
          l2aer(i,k,ns) = l2(ns,k)
          l3aer(i,k,ns) = l3(ns,k)
          l4aer(i,k,ns) = l4(ns,k)
          l5aer(i,k,ns) = l5(ns,k)
          l6aer(i,k,ns) = l6(ns,k)
          l7aer(i,k,ns) = l7(ns,k)
       enddo
    enddo
!   if(i.eq.100) then
!      do ns = 1,nswbands
!         do k = kts,kte
!            call mpas_log_write('$i $i $i $r $r $r $r $r',intArgs=(/i,ns,k/),realArgs=(/tauaersw(i,k,ns), &
!                                extaersw(i,k,ns),gaersw(i,k,ns),waersw(i,k,ns),bscoefsw(i,k,ns)/))
!         enddo
!         call mpas_log_write(' ')
!      enddo
!   endif


    !--- longwave optical properties:
    do k = kts,kte
       do nl = 1,nlwbands
          tauaerlw(i,k,nl) = amax1(lwtauaer(nl,k),1.e-20)
          extaerlw(i,k,nl) = amax1(lwextaer(nl,k),1.e-20)
       enddo
    enddo

 enddo

 call mpas_log_write('--- end module_gocart_opt:')

 end subroutine gocart_optical_averaging

!=================================================================================================================
 subroutine optical_prep_gocart(nbin_o,num_chem,chem,alt,     &
                    relhum,radius_core,radius_wet,number_bin, &
                    swrefindx,swrefindx_core,swrefindx_shell, &
                    lwrefindx,lwrefindx_core,lwrefindx_shell, &
                    uoc,its,ite,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte
 integer,intent(in):: nbin_o
 integer,intent(in):: num_chem

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: alt,relhum
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte,num_chem):: chem


!--- output arguments:
 real(kind=RKIND),intent(out),dimension(its:ite,kts:kte,1:nbin_o):: & 
    radius_wet,number_bin,radius_core

 complex,intent(out),dimension(its:ite,kts:kte,1:nbin_o,nswbands):: &
    swrefindx,swrefindx_core,swrefindx_shell
 complex,intent(out),dimension(its:ite,kts:kte,1:nbin_o,nlwbands):: &
    lwrefindx,lwrefindx_core,lwrefindx_shell


!--- local variables and arrays:
 integer:: i,j,k,l,m,n,isize,itype,iphase
 integer:: ns
 integer:: iflag

 complex:: ref_index_lvcite  , ref_index_nh4hso4 ,                   &
           ref_index_nh4msa  , ref_index_nh4no3  , ref_index_nh4cl , &
           ref_index_nano3   , ref_index_na2so4  ,                   &
           ref_index_na3hso4 , ref_index_nahso4  , ref_index_namsa , &
           ref_index_caso4   , ref_index_camsa2  , ref_index_cano3 , &
           ref_index_cacl2   , ref_index_caco3   , ref_index_h2so4 , &
           ref_index_hhso4   , ref_index_hno3    , ref_index_hcl   , &
           ref_index_msa     , ref_index_bc      ,                   &
           ref_index_oin     , ref_index_aro1    , ref_index_aro2  , &
           ref_index_alk1    , ref_index_ole1    , ref_index_api1  , &
           ref_index_api2    , ref_index_lim1    , ref_index_lim2  , &
           ri_dum            , ri_ave_a

 complex,dimension(nswbands):: & ! now only 5 aerosols have wave-dependent refr 
    swref_index_oc,swref_index_dust,swref_index_nh4so4,swref_index_nacl,swref_index_h2o
 complex,dimension(nlwbands):: & ! now only 5 aerosols have wave-dependent refr 
    lwref_index_oc,lwref_index_dust,lwref_index_nh4so4,lwref_index_nacl,lwref_index_h2o

 real(kind=RKIND),parameter:: FRAC2Aitken=0.25 ! Fraction of modal mass in Aitken mode
                                               ! (applied globally to each species).


 real(kind=RKIND)::  dens_so4  , dens_no3  , dens_cl   , dens_msa  , dens_co3 ,  &
                     dens_nh4  , dens_na   , dens_ca   , dens_oin  , dens_oc  ,  &
                     dens_bc   , dens_aro1 , dens_aro2 , dens_alk1 , dens_ole1,  &
                     dens_api1 , dens_api2 , dens_lim1 , dens_lim2 , dens_h2o ,  &
                     dens_dust
 real(kind=RKIND)::  mass_so4  , mass_no3  , mass_cl   , mass_msa  , mass_co3 ,  &
                     mass_nh4  , mass_na   , mass_ca   , mass_oin  , mass_oc  ,  &
                     mass_bc   , mass_aro1 , mass_aro2 , mass_alk1 , mass_ole1,  &
                     mass_api1 , mass_api2 , mass_lim1 , mass_lim2 , mass_h2o,   &
                     mass_dust
 real(kind=RKIND)::  mass_so4i , mass_no3i , mass_cli  , mass_msai , mass_co3i,  &
                     mass_nh4i , mass_nai  , mass_cai  , mass_oini , mass_oci ,  &
                     mass_bci  , mass_aro1i, mass_aro2i, mass_alk1i, mass_ole1i, &
                     mass_ba1i , mass_ba2i,  mass_ba3i , mass_ba4i , mass_pai,   &
                     mass_h2oi , mass_dusti
 real(kind=RKIND)::  mass_so4j , mass_no3j , mass_clj  , mass_msaj , mass_co3j,  &
                     mass_nh4j , mass_naj  , mass_caj  , mass_oinj , mass_ocj ,  &
                     mass_bcj  , mass_aro1j, mass_aro2j, mass_alk1j, mass_ole1j, &
                     mass_ba1j , mass_ba2j,  mass_ba3j , mass_ba4j , mass_paj,   &
                     mass_h2oj , mass_dustj
 real(kind=RKIND)::  mass_antha, mass_seas, mass_soil

 real(kind=RKIND):: vol_so4   , vol_no3   , vol_cl    , vol_msa   , vol_co3  ,  &
                    vol_nh4   , vol_na    , vol_ca    , vol_oin   , vol_oc   ,  &
                    vol_bc    , vol_aro1  , vol_aro2  , vol_alk1  , vol_ole1 ,  &
                    vol_api1  , vol_api2  , vol_lim1  , vol_lim2  , vol_h2o  ,  & 
                    vol_dust

 real(kind=RKIND):: conv1a, conv1b, conv1sulf
 real(kind=RKIND):: mass_dry_a, mass_wet_a, vol_dry_a , vol_wet_a , vol_shell,  &
                    dp_dry_a  , dp_wet_a  , num_a     , dp_bc_a
 real(kind=RKIND):: ifac, jfac, cfac

 real(kind=RKIND):: dgnum_um,drydens,duma,dlo_um,dhi_um,dgmin,sixpi,ss1,ss2,ss3,dtemp
 real(kind=RKIND),dimension(1:nbin_o):: xnum_secti,xnum_sectj,xnum_sectc
 real(kind=RKIND),dimension(1:nbin_o):: xmas_secti,xmas_sectj,xmas_sectc
 real(kind=RKIND),dimension(1:nbin_o):: xdia_um,xdia_cm

!--- SAM variables needed to convert GOCART sectional dust and seasalt to MOZAIC sections:
 integer,parameter:: nbin_omoz = 8
 integer:: istop
 integer:: uoc

 integer,save:: kcall
 data  kcall / 0 /

 real(kind=RKIND),dimension(4,nbin_omoz),save    :: seasfrc_goc8bin ! GOCART seasalt size distibution
                                                                    ! mass fracs in MOSAIC 8-bins
 real(kind=RKIND),dimension(ndust,nbin_omoz),save:: dustfrc_goc8bin ! GOCART dust size distibution
                                                                    ! mass fracs in MOSAIC 8-bins
 real(kind=RKIND):: dgnum, dhi, dlo, xlo, xhi, dxbin, relh_frc
 real(kind=RKIND):: mass_bc1 ,mass_bc2 ,vol_bc2  , mass_bc1j, &
                    mass_bc2j,mass_bc1i,mass_bc2i, vol_soil
 real(kind=RKIND),dimension(nbin_o):: dlo_sectm,dhi_sectm

 real(kind=R8KIND):: dlogoc,dhigoc

!-----------------------------------------------------------------------------------------------------------------

  if (uoc == 1) then     ! mklose
     den_dust(1) = 2650. ! change dust density in first bin for UoC dust emission schemes
  endif

!
!  real  sginin,sginia,sginic from module_data_sorgam.F
! 
! Mass from modal distribution is divided into individual sections before
! being passed back into the Mie routine.
! * currently use the same size bins as 8 default MOSAIC size bins
! * dlo_um and dhi_um define the lower and upper bounds of individual sections
!   used to compute optical properties
! * sigmas for 3 modes taken from module_sorgam_data.F
! * these parameters are needed by sect02 that is called later
! * sginin=1.7, sginia=2.0, sginic=2.5
!
      sixpi=6.0/3.14159265359
      dlo_um=0.0390625
      dhi_um=10.0
      drydens=1.8
      iflag=2
      duma=1.0
      dgmin=1.0e-07 ! in (cm)
      dtemp=dlo_um
      do isize=1,nbin_o
        xdia_um(isize)=(dtemp+dtemp*2.0)/2.0
        dtemp=dtemp*2.0
      enddo
        if (kcall .eq. 0) then
! 7/21/09 SAM calculate sectional contributions from GOCART seasalt and dust
        dlo = dlo_um*1.0e-6
        dhi = dhi_um*1.0e-6
        xlo = log( dlo )
        xhi = log( dhi )
        dxbin = (xhi - xlo)/nbin_o
        do n = 1, nbin_o
            dlo_sectm(n) = exp( xlo + dxbin*(n-1) )
            dhi_sectm(n) = exp( xlo + dxbin*n )
        end do
!  real, save :: seasfrc_goc8bin(4,nbin_o)   ! GOCART seasalt size distibution - mass fracs in MOSAIC 8-bins
!  real, save :: dustfrc_goc8bin(ndust,nbin_o)   ! GOCART dust size distibution - mass fracs in MOSAIC 8-bins
!  USE module_data_gocart_seas
! real*8, DIMENSION (4), PARAMETER :: ra(4)=(/1.d-1,5.d-1,1.5d0,5.0d0/)
! real*8, DIMENSION (4), PARAMETER :: rb(4)=(/5.d-1,1.5d0,5.d0,1.d1/)
! real*8, DIMENSION (4), PARAMETER :: den_seas(4)=(/2.2d3,2.2d3,2.2d3,2.2d3/)
! real*8, DIMENSION (4), PARAMETER :: reff_seas(4)=(/0.30D-6,1.00D-6,3.25D-6,7.50D-6/)
!  USE module_data_gocart_dust, only:  ndust, reff_dust, den_dust
!  real*8, DIMENSION (5), PARAMETER :: den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
!  real*8, DIMENSION (5), PARAMETER :: reff_dust(5)=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/)
! Seasalt bin mass fractions
        seasfrc_goc8bin=0.
!      WRITE(*,*)'Seasalt mass fractions'
!      WRITE(*,*)'         ','       ',(dlo_sectm(n),n=1,nbin_o)
!      WRITE(*,*)'         ','       ',(dhi_sectm(n),n=1,nbin_o)
       do m =1, 4  ! loop over seasalt size bins
       dlogoc = ra_seas(m)*2.E-6  ! low diameter limit (m)
       dhigoc = rb_seas(m)*2.E-6  ! hi diameter limit (m)
        do n = 1, nbin_o
        seasfrc_goc8bin(m,n)=max(DBLE(0.),min(DBLE(log(dhi_sectm(n))),log(dhigoc))- &
                             max(log(dlogoc),DBLE(log(dlo_sectm(n)))) )/(log(dhigoc)-log(dlogoc))

       end do
!      WRITE(*,*)m,dlogoc,dhigoc,(seasfrc_goc8bin(m,n),n=1,nbin_o)
       end do
! Dust bin mass fractions
!      WRITE(*,*)'Dust mass fractions'
!      WRITE(*,*)'         ','       ',(dlo_sectm(n),n=1,nbin_o)
!      WRITE(*,*)'         ','       ',(dhi_sectm(n),n=1,nbin_o)
        dustfrc_goc8bin=0.
       do m =1, ndust  ! loop over dust size bins
        dlogoc = ra_dust(m)*2.E-6  ! low diameter limit (m)
        dhigoc = rb_dust(m)*2.E-6  ! hi diameter limit (m)
        do n = 1, nbin_o
        dustfrc_goc8bin(m,n)=max(DBLE(0.),min(DBLE(log(dhi_sectm(n))),log(dhigoc))- &
                             max(log(dlogoc),DBLE(log(dlo_sectm(n)))) )/(log(dhigoc)-log(dlogoc))

       end do
!      WRITE(*,*)m,dlogoc,dhigoc,(dustfrc_goc8bin(m,n),n=1,nbin_o)
       end do
        kcall=kcall+1
!       ISTOP=1
!       IF(ISTOP.EQ.1)THEN
!       STOP
!       ENDIF
        endif
!
! Define refractive indicies
! * assume na and cl are the same as nacl
! * assume so4, no3, and nh4 are the same as nh4no3
! * assume ca and co3 are the same as caco3
! * assume msa is just msa
! Further work:
! * to be more precise, need to compute electrolytes to apportion
!   so4, no3, nh4, na, cl, msa, ca, co3 among various componds
!   as was done previously in module_mosaic_therm.F
!
      do ns = 1, nswbands
      swref_index_nh4so4(ns) = cmplx(refrsw_sulf(ns),refisw_sulf(ns))
      swref_index_oc(ns) = cmplx(refrsw_oc(ns),refisw_oc(ns))
      swref_index_dust(ns) = cmplx(refrsw_dust(ns),refisw_dust(ns))
      swref_index_nacl(ns) = cmplx(refrsw_seas(ns),refisw_seas(ns))
      swref_index_h2o(ns) = cmplx(refrwsw(ns),refiwsw(ns))
      enddo
      do ns = 1, nlwbands
      lwref_index_nh4so4(ns) = cmplx(refrlw_sulf(ns),refilw_sulf(ns))
      lwref_index_oc(ns) = cmplx(refrlw_oc(ns),refilw_oc(ns))
      lwref_index_dust(ns) = cmplx(refrlw_dust(ns),refilw_dust(ns))
      lwref_index_nacl(ns) = cmplx(refrlw_seas(ns),refilw_seas(ns))
      lwref_index_h2o(ns) = cmplx(refrwlw(ns),refiwlw(ns))
      enddo
!     ref_index_nh4so4 = cmplx(1.52,0.)
      ref_index_lvcite = cmplx(1.50,0.)
      ref_index_nh4hso4= cmplx(1.47,0.)
      ref_index_nh4msa = cmplx(1.50,0.)     ! assumed
      ref_index_nh4no3 = cmplx(1.50,0.)
      ref_index_nh4cl  = cmplx(1.50,0.)
!     ref_index_nacl   = cmplx(1.45,0.)
      ref_index_nano3  = cmplx(1.50,0.)
      ref_index_na2so4 = cmplx(1.50,0.)
      ref_index_na3hso4= cmplx(1.50,0.)
      ref_index_nahso4 = cmplx(1.50,0.)
      ref_index_namsa  = cmplx(1.50,0.)     ! assumed
      ref_index_caso4  = cmplx(1.56,0.006)
      ref_index_camsa2 = cmplx(1.56,0.006)  ! assumed
      ref_index_cano3  = cmplx(1.56,0.006)
      ref_index_cacl2  = cmplx(1.52,0.006)
      ref_index_caco3  = cmplx(1.68,0.006)
      ref_index_h2so4  = cmplx(1.43,0.)
      ref_index_hhso4  = cmplx(1.43,0.)
      ref_index_hno3   = cmplx(1.50,0.)
      ref_index_hcl    = cmplx(1.50,0.)
      ref_index_msa    = cmplx(1.43,0.)     ! assumed
!     ref_index_oc     = cmplx(1.45,0.)  ! JCB, Feb. 20, 2008: no complex part?
! JCB, Feb. 20, 2008:  set the refractive index of BC equal to the
! midpoint of ranges given in Bond and Bergstrom, Light absorption by
! carboneceous particles: an investigative review 2006, Aerosol Sci.
! and Tech., 40:27-67.
!     ref_index_bc     = cmplx(1.82,0.74) old value
      ref_index_bc     = cmplx(1.85,0.71)
      ref_index_oin    = cmplx(1.55,0.006)  ! JCB, Feb. 20, 2008:  "other inorganics"
      ref_index_aro1   = cmplx(1.45,0.)
      ref_index_aro2   = cmplx(1.45,0.)
      ref_index_alk1   = cmplx(1.45,0.)
      ref_index_ole1   = cmplx(1.45,0.)
      ref_index_api1   = cmplx(1.45,0.)
      ref_index_api2   = cmplx(1.45,0.)
      ref_index_lim1   = cmplx(1.45,0.)
      ref_index_lim2   = cmplx(1.45,0.)
!     ref_index_h2o    = cmplx(1.33,0.)
!
! densities in g/cc
!
      dens_so4   = 1.8        ! used
      dens_no3   = 1.8        ! used
      dens_cl    = 2.2        ! used
      dens_msa   = 1.8        ! used
      dens_co3   = 2.6        ! used
      dens_nh4   = 1.8        ! used
      dens_na    = 2.2        ! used
      dens_ca    = 2.6        ! used
      dens_oin   = 2.6        ! used
      dens_dust  = 2.6        ! used
      dens_oc    = 1.0        ! used
! JCB, Feb. 20, 2008:  the density of BC is updated to reflect values
! published by Bond and Bergstrom, Light absorption by carboneceous
! particles: an investigative review 2006, Aerosol Sci. and Tech., 40:27-67.
!     dens_bc    = 1.7        ! used, old value
      dens_bc    =  1.8       ! midpoint of Bond and Bergstrom value
      dens_aro1  = 1.0
      dens_aro2  = 1.0
      dens_alk1  = 1.0
      dens_ole1  = 1.0
      dens_api1  = 1.0
      dens_api2  = 1.0
      dens_lim1  = 1.0
      dens_lim2  = 1.0
      dens_h2o   = 1.0
!
      swrefindx=0.0
      lwrefindx=0.0
      radius_wet=0.0
      number_bin=0.0
      radius_core=0.0
      swrefindx_core=0.0
      swrefindx_shell=0.0
      lwrefindx_core=0.0
      lwrefindx_shell=0.0
!
! units:
! * mass     - g/cc(air)
! * number   - #/cc(air)
! * volume   - cc(air)/cc(air)
! * diameter - cm
!
      do k = kts, kte
      do i = its, ite
        mass_so4i = 0.0
        mass_so4j = 0.0
        mass_no3i = 0.0
        mass_no3j = 0.0
        mass_nh4i = 0.0
        mass_nh4j = 0.0
        mass_oini = 0.0
        mass_oinj = 0.0
        mass_dusti = 0.0
        mass_dustj = 0.0
        mass_aro1i = 0.0
        mass_aro1j = 0.0
        mass_aro2i = 0.0
        mass_aro2j = 0.0
        mass_alk1i = 0.0
        mass_alk1j = 0.0
        mass_ole1i = 0.0
        mass_ole1j = 0.0
        mass_ba1i = 0.0
        mass_ba1j = 0.0
        mass_ba2i = 0.0
        mass_ba2j = 0.0
        mass_ba3i = 0.0
        mass_ba3j = 0.0
        mass_ba4i = 0.0
        mass_ba4j = 0.0
        mass_pai = 0.0
        mass_paj = 0.0
        mass_oci = 0.0
        mass_ocj = 0.0
        mass_bci = 0.0
        mass_bcj = 0.0
        mass_bc1i = 0.0
        mass_bc1j = 0.0
        mass_bc2i = 0.0
        mass_bc2j = 0.0
        mass_cai = 0.0
        mass_caj = 0.0
        mass_co3i = 0.0
        mass_co3j = 0.0
        mass_nai = 0.0
        mass_naj = 0.0
        mass_cli = 0.0
        mass_clj = 0.0
        mass_msai = 0.0
        mass_msaj = 0.0
        mass_nai = 0.0
        mass_naj = 0.0
        mass_cli = 0.0
        mass_clj = 0.0
        mass_h2oi = 0.0
        mass_h2oj = 0.0
        mass_antha = 0.0
        mass_seas = 0.0
        mass_soil = 0.0
        mass_cl = 0.0
        mass_na = 0.0
        mass_msa = 0.0

! convert ug / kg dry air to g / cc  air
        conv1a = (1.0/alt(i,k)) * 1.0e-12
! convert # / kg dry air to # / cc  air
        conv1b = (1.0/alt(i,k)) * 1.0e-6
! convert ppmv sulfate (and coincidentally MSA) to g / cc  air
        conv1sulf = (1.0/alt(i,k)) * 1.0e-9 * 96./28.97

! Accumulation mode...
! SAM 7/18/09 - Put fraction of GOCART sulfate, organic, black carbon masses into modal accumulation mode
        mass_oinj = (1.-FRAC2Aitken)*chem(i,k,p_p25)*conv1a
        mass_so4j= (1.-FRAC2Aitken)*chem(i,k,p_sulf)*conv1sulf
        mass_nh4j= (1.-FRAC2Aitken)*chem(i,k,p_sulf)*conv1sulf*(nh4_mfac-1.)
        mass_aro1j= (1.-FRAC2Aitken)*chem(i,k,p_oc1)*conv1a*oc_mfac
        mass_aro2j= (1.-FRAC2Aitken)*chem(i,k,p_oc2)*conv1a*oc_mfac
        mass_bc1j= (1.-FRAC2Aitken)*chem(i,k,p_bc1)*conv1a
        mass_bc2j= (1.-FRAC2Aitken)*chem(i,k,p_bc2)*conv1a
        mass_bcj= mass_bc1j + mass_bc2j
        if( p_msa .gt. 1) mass_msaj= (1.-FRAC2Aitken)*chem(i,k,p_msa)*conv1sulf


! Aitken mode...
! SAM 7/18/09 - Put fraction of GOCART sulfate, organic, black carbon masses into modal Aitken mode
        mass_oini = FRAC2Aitken*chem(i,k,p_p25)*conv1a
        mass_so4i= FRAC2Aitken*chem(i,k,p_sulf)*conv1sulf
        mass_nh4i= FRAC2Aitken*chem(i,k,p_sulf)*conv1sulf*(nh4_mfac-1.)
        mass_aro1i= FRAC2Aitken*chem(i,k,p_oc1)*conv1a*oc_mfac
        mass_aro2i= FRAC2Aitken*chem(i,k,p_oc2)*conv1a*oc_mfac
        mass_bc1i= FRAC2Aitken*chem(i,k,p_bc1)*conv1a
        mass_bc2i= FRAC2Aitken*chem(i,k,p_bc2)*conv1a
        mass_bci= mass_bc1i + mass_bc2i
        if( p_msa .gt. 1) mass_msai= FRAC2Aitken*chem(i,k,p_msa)*conv1sulf

!
! Now divide mass into sections which is done by sect02:
! * xmas_secti is for aiken mode
! * xmas_sectj is for accumulation mode
! * xmas_sectc is for coarse mode
! * sect02 expects input in um
! * pass in generic mass of 1.0 just to get a percentage distribution of mass among bins
!
!!      ss1=log(sginin)
!!      ss2=exp(ss1*ss1*36.0/8.0) 
!!      ss3=(sixpi*vol_ai/(num_ai*ss2))**0.3333333
!!      dgnum_um=amax1(dgmin,ss3)*1.0e+04
        dgnum_um=dginin*1.E6
        call sect02(dgnum_um,sginin,drydens,iflag,duma,nbin_o,dlo_um,dhi_um, &
          xnum_secti,xmas_secti)
!!      ss1=log(sginia)
!!      ss2=exp(ss1*ss1*36.0/8.0) 
!!      ss3=(sixpi*vol_aj/(num_aj*ss2))**0.3333333
!!      dgnum_um=amax1(dgmin,ss3)*1.0e+04
        dgnum_um=dginia*1.E6
        call sect02(dgnum_um,sginia,drydens,iflag,duma,nbin_o,dlo_um,dhi_um, &
          xnum_sectj,xmas_sectj)
!!      ss1=log(sginic)
!!      ss2=exp(ss1*ss1*36.0/8.0) 
!!      ss3=(sixpi*vol_ac/(num_ac*ss2))**0.3333333
        dgnum_um=dginic*1.E6
        call sect02(dgnum_um,sginic,drydens,iflag,duma,nbin_o,dlo_um,dhi_um, &
          xnum_sectc,xmas_sectc)

        do isize = 1, nbin_o
          xdia_cm(isize)=xdia_um(isize)*1.0e-04
          mass_so4 = mass_so4i*xmas_secti(isize) + mass_so4j*xmas_sectj(isize)
          mass_no3 = mass_no3i*xmas_secti(isize) + mass_no3j*xmas_sectj(isize)
          mass_nh4 = mass_nh4i*xmas_secti(isize) + mass_nh4j*xmas_sectj(isize)
          mass_oin = mass_oini*xmas_secti(isize) + mass_oinj*xmas_sectj(isize) + &
                     mass_soil*xmas_sectc(isize) + mass_antha*xmas_sectc(isize)
          if( p_msa .gt. 1) then
             mass_msa = mass_msai*xmas_secti(isize) + mass_msaj*xmas_sectj(isize)
          endif
! GOCART OC mass_aero1 is hydrophobic, mass_aero2 is hydrophylic
          mass_aro1 = mass_aro1j*xmas_sectj(isize) + mass_aro1i*xmas_secti(isize)
          mass_aro2 = mass_aro2j*xmas_sectj(isize) + mass_aro2i*xmas_secti(isize)
          mass_oc = mass_aro1 + mass_aro2
! GOCART BC mass_bc1 is hydrophobic, mass_bc2 is hydrophylic
          mass_bc1  = mass_bc1i*xmas_secti(isize) + mass_bc1j*xmas_sectj(isize)
          mass_bc2  = mass_bc2i*xmas_secti(isize) + mass_bc2j*xmas_sectj(isize)
          mass_bc = mass_bc1 + mass_bc2
! Add in seasalt and dust from GOCART sectional distributions
       n = 0
       mass_seas = 0.0
       do m =p_seas_1,  p_seas_3 ! loop over seasalt size bins less than 10 um diam
       n = n+1
        mass_seas=mass_seas+seasfrc_goc8bin(n,isize)*chem(i,k,m)
       end do
       n = 0
       mass_soil = 0.0
       do m =p_dust_1,  p_dust_1+ndust-2 ! loop over dust size bins less than 10 um diam
       n = n+1
        mass_soil=mass_soil+dustfrc_goc8bin(n,isize)*chem(i,k,m)
       end do
       mass_cl=mass_seas*conv1a*35.4530/58.4428
       mass_na=mass_seas*conv1a*22.9898/58.4428
       mass_soil=mass_soil*conv1a
!         mass_h2o = 0.0 ! testing purposes only
          vol_so4 = mass_so4 / dens_so4
          vol_no3 = mass_no3 / dens_no3
          vol_nh4 = mass_nh4 / dens_nh4
          vol_oin = mass_oin / dens_oin
          vol_oc  = mass_oc  / dens_oc
          vol_aro2 = mass_aro2 / dens_oc
          vol_bc  = mass_bc  / dens_bc
          vol_bc2 = mass_bc2 / dens_bc
          vol_na  = mass_na  / dens_na
          vol_cl  = mass_cl  / dens_cl
          vol_soil  = mass_soil  / dens_dust
          vol_msa  = mass_msa  / dens_msa
!         vol_h2o = mass_h2o / dens_h2o
! 7/23/09 SAM calculate vol_h2o from kappas in Petters and Kreidenweis ACP, 2007, vol. 7, 1961-1971.
!  Their kappas are the hygroscopicities used in Abdul-Razzak and Ghan, 2004, JGR, V105, p. 6837-6844. 
!  These kappas are defined in module_data_sorgam and module_data_mosaic_asect.
!  Note that hygroscopicities are at 298K and specific surface tension - further refinement could
!  include temperature dependence in Petters and Kreidenweis
! Also, for hygroscopic BC part, assume kappa of OC (how can BC be hydrophylic?)
               relh_frc=amin1(.9,relhum(i,k)) !0.8   ! Put in fractional relative humidity, max of .9, here
                  vol_h2o = vol_so4*hygro_so4_aer + vol_aro2*hygro_oc_aer + &
                            vol_nh4*hygro_nh4_aer                         + &
                  vol_cl*hygro_cl_aer + vol_na*hygro_na_aer + vol_msa*hygro_msa_aer + &
                  vol_oin*hygro_oin_aer + vol_bc2*hygro_oc_aer  + vol_soil*hygro_dust_aer
                  vol_h2o = relh_frc*vol_h2o/(1.-relh_frc)
                  mass_h2o = vol_h2o*dens_h2o
          mass_dry_a = mass_so4 + mass_no3 + mass_nh4 + mass_oin + &
                       mass_oc  + mass_bc  + mass_na  + mass_cl  + &
                       mass_soil
          mass_wet_a = mass_dry_a + mass_h2o 
          vol_dry_a  = vol_so4  + vol_no3  + vol_nh4  + vol_oin  + &
                       vol_oc   + vol_bc   + vol_na   + vol_cl  + &
                       vol_soil
          vol_wet_a  = vol_dry_a + vol_h2o
          vol_shell  = vol_wet_a - vol_bc
          num_a      = vol_wet_a / (0.52359877*xdia_cm(isize)*xdia_cm(isize)*xdia_cm(isize))

          !shortwave 
          do ns=1,nswbands
          ri_dum     = (0.0,0.0)
          ri_dum     = (swref_index_nh4so4(ns) * mass_so4 / dens_so4) +  &
                       (ref_index_nh4no3 * mass_no3 / dens_no3) +  &
                       (ref_index_nh4no3 * mass_nh4 / dens_nh4) +  &
                       (ref_index_oin    * mass_oin / dens_oin) +  &
                       (swref_index_dust(ns) * mass_soil / dens_dust) + &
                       (swref_index_oc(ns)     * mass_oc  / dens_oc) +   &
                       (ref_index_bc     * mass_bc  / dens_bc) +   &
                       (swref_index_nacl(ns)   * mass_na  / dens_na) +   &
                       (swref_index_nacl(ns)   * mass_cl  / dens_cl) +   &
                       (ref_index_msa    * mass_msa / dens_msa) +  &
                       (swref_index_h2o(ns) * mass_h2o / dens_h2o) 
!
! for some reason MADE/SORGAM occasionally produces zero aerosols so
! need to add a check here to avoid divide by zero
!
          IF(num_a .lt. 1.0e-20 .or. vol_wet_a .lt. 1.0e-20) then
            dp_dry_a   = xdia_cm(isize)
            dp_wet_a   = xdia_cm(isize)
            dp_bc_a    = xdia_cm(isize)
            ri_ave_a   = 0.0
            ri_dum     = 0.0
          else
            dp_dry_a   = (1.90985*vol_dry_a/num_a)**0.3333333
            dp_wet_a   = (1.90985*vol_wet_a/num_a)**0.3333333
            dp_bc_a    = (1.90985*vol_bc/num_a)**0.3333333
            ri_ave_a   = ri_dum/vol_wet_a
            ri_dum   = (swref_index_nh4so4(ns) * mass_so4 / dens_so4) +  &
                       (ref_index_nh4no3 * mass_no3 / dens_no3) +  &
                       (ref_index_nh4no3 * mass_nh4 / dens_nh4) +  &
                       (ref_index_msa    * mass_msa / dens_msa) +  &
                       (ref_index_oin    * mass_oin / dens_oin) +  &
                       (swref_index_dust(ns)    * mass_soil / dens_dust) + &
                       (swref_index_oc(ns)     * mass_oc  / dens_oc) +   &
                       (swref_index_nacl(ns)   * mass_na  / dens_na) +   &
                       (swref_index_nacl(ns)   * mass_cl  / dens_cl) +   &
                       (swref_index_h2o(ns)    * mass_h2o / dens_h2o) 
          endif
          if(dp_wet_a/2.0 .lt. dlo_um*1.0e-4/2.0) then
            swrefindx(i,k,isize,ns)    = (1.5,0.0)
            radius_wet(i,k,isize) =dlo_um*1.0e-4/2.0
            number_bin(i,k,isize) =num_a
            radius_core(i,k,isize) =0.0
            swrefindx_core(i,k,isize,ns) = ref_index_bc
            swrefindx_shell(i,k,isize,ns) = ref_index_oin
          elseif(num_a .lt. 1.e-20 .or. vol_shell .lt. 1.0e-20) then
            swrefindx(i,k,isize,ns)    = (1.5,0.0)
            radius_wet(i,k,isize) =dlo_um*1.0e-4/2.0
            number_bin(i,k,isize) =num_a
            radius_core(i,k,isize) =0.0
            swrefindx_core(i,k,isize,ns) = ref_index_bc
            swrefindx_shell(i,k,isize,ns) = ref_index_oin
          else
            swrefindx(i,k,isize,ns)    =ri_ave_a
            radius_wet(i,k,isize) =dp_wet_a/2.0
            number_bin(i,k,isize) =num_a
            radius_core(i,k,isize) =dp_bc_a/2.0
            swrefindx_core(i,k,isize,ns) =ref_index_bc
            swrefindx_shell(i,k,isize,ns) =ri_dum/vol_shell
          endif
        enddo  ! ns shortwave

          !longwave 
          do ns=1,nlwbands
          ri_dum     = (0.0,0.0)
          ri_dum     = (lwref_index_nh4so4(ns) * mass_so4 / dens_so4) +  &
                       (ref_index_nh4no3 * mass_no3 / dens_no3) +  &
                       (ref_index_nh4no3 * mass_nh4 / dens_nh4) +  &
                       (ref_index_oin    * mass_oin / dens_oin) +  &
                       (lwref_index_dust(ns) * mass_soil / dens_dust) + &
                       (lwref_index_oc(ns)     * mass_oc  / dens_oc) +   &
                       (ref_index_bc     * mass_bc  / dens_bc) +   &
                       (lwref_index_nacl(ns)   * mass_na  / dens_na) +   &
                       (lwref_index_nacl(ns)   * mass_cl  / dens_cl) +   &
                       (ref_index_msa    * mass_msa / dens_msa) +  &
                       (lwref_index_h2o(ns) * mass_h2o / dens_h2o) 
!
! for some reason MADE/SORGAM occasionally produces zero aerosols so
! need to add a check here to avoid divide by zero
!
          IF(num_a .lt. 1.0e-20 .or. vol_wet_a .lt. 1.0e-20) then
            dp_dry_a   = xdia_cm(isize) 
            dp_wet_a   = xdia_cm(isize) 
            dp_bc_a    = xdia_cm(isize)
            ri_ave_a   = 0.0
            ri_dum     = 0.0
          else
            dp_dry_a   = (1.90985*vol_dry_a/num_a)**0.3333333
            dp_wet_a   = (1.90985*vol_wet_a/num_a)**0.3333333
            dp_bc_a    = (1.90985*vol_bc/num_a)**0.3333333
            ri_ave_a   = ri_dum/vol_wet_a
            ri_dum   = (lwref_index_nh4so4(ns) * mass_so4 / dens_so4) +  &
                       (ref_index_nh4no3 * mass_no3 / dens_no3) +  &
                       (ref_index_nh4no3 * mass_nh4 / dens_nh4) +  &
                       (ref_index_oin    * mass_oin / dens_oin) +  &
                       (lwref_index_dust(ns)    * mass_oin / dens_dust)  +  &
                       (lwref_index_oc(ns)     * mass_oc  / dens_oc) +   &
                       (lwref_index_nacl(ns)   * mass_na  / dens_na) +   &
                       (lwref_index_nacl(ns)   * mass_cl  / dens_cl) +   &
                       (lwref_index_h2o(ns)    * mass_h2o / dens_h2o)
          endif
          if(dp_wet_a/2.0 .lt. dlo_um*1.0e-4/2.0) then
            lwrefindx(i,k,isize,ns)    = (1.5,0.0)
            radius_wet(i,k,isize) =dlo_um*1.0e-4/2.0
            number_bin(i,k,isize) =num_a
            radius_core(i,k,isize) =0.0
            lwrefindx_core(i,k,isize,ns) = ref_index_bc
            lwrefindx_shell(i,k,isize,ns) = ref_index_oin
          elseif(num_a .lt. 1.e-20 .or. vol_shell .lt. 1.0e-20) then
            lwrefindx(i,k,isize,ns)    = (1.5,0.0)
            radius_wet(i,k,isize) =dlo_um*1.0e-4/2.0
            number_bin(i,k,isize) =num_a
            radius_core(i,k,isize) =0.0
            lwrefindx_core(i,k,isize,ns) = ref_index_bc
            lwrefindx_shell(i,k,isize,ns) = ref_index_oin
          else
            lwrefindx(i,k,isize,ns)    =ri_ave_a
            radius_wet(i,k,isize) =dp_wet_a/2.0
            number_bin(i,k,isize) =num_a
            radius_core(i,k,isize) =dp_bc_a/2.0
            lwrefindx_core(i,k,isize,ns) =ref_index_bc
            lwrefindx_shell(i,k,isize,ns) =ri_dum/vol_shell
          endif
        enddo  ! ns longwave
!  refr=real(refindx(i,k,j,isize))
        enddo  !isize
      enddo  !i
      enddo  !k

 end subroutine optical_prep_gocart

!=================================================================================================================
!
      subroutine sect02(dgnum_um,sigmag,drydens,iflag,duma,nbin,dlo_um,dhi_um, &
        xnum_sect,xmas_sect)
!
!   user specifies a single log-normal mode and a set of section boundaries
!   prog calculates mass and number for each section
!
        implicit none
        REAL, DIMENSION(nbin), INTENT(OUT) :: xnum_sect, xmas_sect
        integer iflag, n, nbin
        real &
          dgnum, dgnum_um, dhi, dhi_um, dlo, dlo_um,   &
          drydens, dstar, duma, dumfrac, dx,   &
          sigmag, sumnum, summas,   &
          sx, sxroot2, thi, tlo, vtot,   &
          x0, x3, xhi, xlo, xmtot, xntot, xvtot
        real dlo_sect(nbin), dhi_sect(nbin)
!       real erfc_num_recipes
        real pi
        parameter (pi = 3.1415926536)
!
        if (iflag .le. 1) then
            xntot = duma
        else
            xmtot = duma
            xntot = duma   !czhao
        end if
!   compute total volume and number for mode
!       dgnum = dgnum_um*1.0e-4
!       sx = log( sigmag )
!       x0 = log( dgnum )
!       x3 = x0 + 3.*sx*sx
!       dstar = dgnum * exp(1.5*sx*sx)
!       if (iflag .le. 1) then
!           xvtot = xntot*(pi/6.0)*dstar*dstar*dstar
!           xmtot = xvtot*drydens*1.0e12
!       else
!           xvtot = xmtot/(drydens*1.0e12)
!           xntot = xvtot/((pi/6.0)*dstar*dstar*dstar)
!       end if
!   compute section boundaries
        dlo = dlo_um*1.0e-4
        dhi = dhi_um*1.0e-4
        xlo = log( dlo )
        xhi = log( dhi )
        dx = (xhi - xlo)/nbin
        do n = 1, nbin
            dlo_sect(n) = exp( xlo + dx*(n-1) )
            dhi_sect(n) = exp( xlo + dx*n )
        end do
!   compute modal "working" parameters including total num/vol/mass
        dgnum = dgnum_um*1.0e-4
        sx = log( sigmag )
        x0 = log( dgnum )
        x3 = x0 + 3.*sx*sx
        dstar = dgnum * exp(1.5*sx*sx)
        if (iflag .le. 1) then
            xvtot = xntot*(pi/6.0)*dstar*dstar*dstar
            xmtot = xvtot*drydens*1.0e12
        else
!czhao      xvtot = xmtot/(drydens*1.0e12)
!czhao      xntot = xvtot/((pi/6.0)*dstar*dstar*dstar)
        end if
!   compute number and mass for each section
        sxroot2 = sx * sqrt( 2.0 )
        sumnum = 0.
        summas = 0.
!       write(22,*)
!       write(22,*) 'dgnum_um, sigmag = ', dgnum_um, sigmag
!       write(22,*) 'drydens =', drydens
!       write(22,*) 'ntot (#/cm3), mtot (ug/m3) = ', xntot, xmtot
!        write(22,9220)
!9220    format( /   &
!        '  n   dlo(um)   dhi(um)       number         mass' / )
!9225    format(   i3, 2f10.6, 2(1pe13.4) )
!9230    format( / 'sum over all sections  ', 2(1pe13.4) )
!9231    format(   'modal totals           ', 2(1pe13.4) )
        do n = 1, nbin
            xlo = log( dlo_sect(n) )
            xhi = log( dhi_sect(n) )
            tlo = (xlo - x0)/sxroot2
            thi = (xhi - x0)/sxroot2
            if (tlo .le. 0.) then
                dumfrac = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
            else
                dumfrac = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
            end if
            xnum_sect(n) = xntot*dumfrac
            tlo = (xlo - x3)/sxroot2
            thi = (xhi - x3)/sxroot2
            if (tlo .le. 0.) then
                dumfrac = 0.5*( erfc_num_recipes(-thi) - erfc_num_recipes(-tlo) )
            else
                dumfrac = 0.5*( erfc_num_recipes(tlo) - erfc_num_recipes(thi) )
            end if
            xmas_sect(n) = xmtot*dumfrac
            sumnum = sumnum + xnum_sect(n)
            summas = summas + xmas_sect(n)
!           write(22,9225) n, 1.e4*dlo_sect(n), 1.e4*dhi_sect(n),   &
!               xnum_sect(n), xmas_sect(n)
        end do
!       write(22,9230) sumnum, summas
!       write(22,9231) xntot, xmtot

      end subroutine  sect02

!-----------------------------------------------------------------------
        real function erfc_num_recipes( x )
!
!   from press et al, numerical recipes, 1990, page 164
!
        implicit none
        real x
        double precision erfc_dbl, dum, t, z
        z = abs(x)
        t = 1.0/(1.0 + 0.5*z)
!       erfc_num_recipes =
!     &   t*exp( -z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
!     &   t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
!     &                                    t*(-1.13520398 +
!     &   t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))
        dum =  ( -z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 +   &
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +   &
                                           t*(-1.13520398 +   &
          t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))
        erfc_dbl = t * exp(dum)
        if (x .lt. 0.0) erfc_dbl = 2.0d0 - erfc_dbl
        erfc_num_recipes = erfc_dbl
        return

        end function erfc_num_recipes

!=================================================================================================================
 end module module_gocart_opt
!=================================================================================================================
