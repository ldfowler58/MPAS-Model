!
!.............................................................................
!
!    Multi-Scale Kain-Fritsch (MSKF) scheme for MPAS
!              Laura Fowler, NCAR 
!                     & 
!              Kiran Alapaty, EPA  
!NOTES:
!   (1) Based off wrf4.2 version of MSKF , June 2020
!   (2) WITHOUT CONVECTIVE MICROPHYSICS & AEROSOL EFFECTS codes,June 2020
!
!..............................................................................
!
!CKAY=Kiran Alapaty, EPA
!CGM = Chris Marciano, NCSU
!PSH = Patrick Hawbecker, NCSU
!TWG = Tim Glotfelty, ORISE/EPA
!JTR = Jacob Radford, NCSU
!
!MULTI-SCALE KAIN-FRITSCH scheme features & publications:
!  (1) With diagnosed deep and shallow KF cloud fraction using
!       CAM3-CAM5 methodology, along with captured liquid and ice condensates.
!       and linking with the RRTMG & Other radiation schemes
! (2) Scale-dependent Dynamic adjustment timescale for KF clouds (both shallow
!     and deep)
! (3) Scale-dependent LCL-based entrainment Methodology that avoids 2-km cloud
!     radius method
! (4) Scale-dependent Fallout Rate
! (5) Scale-dependent Stabilization Capacity
! (6) Elimination of "double counting" when environment is saturated
! (7) Inclusion of Convective Momentum Transport (CMT):Zhang&McFarlane,JGR,100,1995
!
! Alapaty et al., 2012: Introducing subgrid-scale cloud feedbacks to radiation
!    for regional meteorological and climate modeling. GRL, V39, I24.
!
! Alapaty et al., 2013:  The Kain-Fritsch Scheme: Science Updates and revisiting
!     gray-scale issues from the NWP and regional climate perspectives. 2013 WRF
!     workshop: URL:
!     http://www.mmm.ucar.edu/wrf/users/workshops/WS2013/ppts/9.2.pdf
!
! Herwehe et al., 2014: Increasing the credibility of regional climate
!     simulations by introducing subgrid-scale cloud-radiation interactions. JGR, 119,
!     5317-5330, doi:10.1002/2014JD021504.
!
! Zheng et al., 2016: Improving High-Resolution Weather Forecasts using the
!     Weather Research and Forecasting (WRF) Model with an Updated Kain-Fritsch
!     Scheme. Mon. Wea. Rev., 144, 833-860
!
! He, J., and K. Alapaty, 2018:  Precipitation partitioning in multiscale
!     atmospheric simulations: Impacts of stability restoration methods. Journal of
!     Geophysical Research: Atmospheres, 123. https://doi.org/10.1029/2018JD028710
!
! Glotfelty, T., K. Alapaty, J. He, P. Hawbecker, X. Song, and G. Zhang, 2019:
!     The Weather Research and Forecasting Model with Aerosol Cloud Interactions
!     (WRF-ACI): Development, Evaluation, and Initial Application. Mon. Wea.
!     Rev., 147, 1491-1511
!................................................................

!=================================================================================================================
 module module_cu_mskf

#if defined(mpas)
 use mpas_kind_types,only: r8=>R8KIND
 use mpas_atmphys_utilities,only: physics_message,physics_error_fatal
#define FATAL_ERROR(M) call physics_error_fatal( M )
#define WRITE_MESSAGE(M) call physics_message( M )
#else
 use module_wrf_error
 use shr_kind_mod,only: r8=>shr_kind_r8
#define FATAL_ERROR(M) call wrf_error_fatal( M )
#define WRITE_MESSAGE(M) call wrf_message( M )
#endif

 use module_cu_mskf_cmp

!
!--------------------------------------------------------------------
!lookup table variables:
 integer, parameter :: kfnt=250,kfnp=220
 real,dimension(kfnt,kfnp),private,save :: ttab,qstab
 real,dimension(kfnp),private,save      :: the0k
 real,dimension(200),private,save       :: alu
 real,private,save                      :: rdpr,rdthk,plutop
!note: KF Lookup table is used by subroutines KF_eta_PARA,TPMIX2,TPMIX2DD,ENVIRTHT
!end of Lookup table variables.


 contains


!=================================================================================================================
   SUBROUTINE MSKF_CPS(                                      &
              ids,ide, jds,jde, kds,kde                      &
             ,ims,ime, jms,jme, kms,kme                      &
             ,its,ite, jts,jte, kts,kte                      &
             ,DT,KTAU,DX,CUDT,ADAPT_STEP_FLAG                &
             ,rho,RAINCV,PRATEC,NCA                          &
             ,U,V,TH,T,W,dz8w,Pcps,pi                        &
             ,W0AVG,XLV0,XLV1,XLS0,XLS1,CP,R,G,EP1           &
             ,EP2,SVP1,SVP2,SVP3,SVPT0                       &
             ,STEPCU,CU_ACT_FLAG,warm_rain,CUTOP,CUBOT       &
             ,QV                                             &
            ! optionals
             ,F_QV    ,F_QC    ,F_QR    ,F_QI    ,F_QS       &
             ,RTHCUTEN,RQVCUTEN,RQCCUTEN,RQRCUTEN            &
             ,RQICUTEN,RQSCUTEN, RQVFTEN                     &
!ckay
             ,cldfra_dp_KF,cldfra_sh_KF,w_up                 &
             ,qc_KF,qi_KF,qr_KF,qs_KF                        & ! TWG
             ,nc_KF,ni_KF,nr_KF,ns_KF                        & ! TWG
             ,ccn_KF,ainc_frac                               & ! TWG
!kf_edrates
             ,UDR_KF,DDR_KF                                  &
             ,UER_KF,DER_KF                                  &
             ,TIMEC_KF,KF_EDRATES                            & 
             ,ZOL,HFX,UST,PBLH                               &   !ckay
             ,aerocu,no_src_types_cu,aercu_fct,aercu_opt     & !PSH/TWG
             ,EFCS,EFIS,EFSS                                 &
             ,RUCUTEN,RVCUTEN,XLAND)                                 !JTR
!
!-------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 Add to avoid memory issues
!-------------------------------------------------------------
   INTEGER,      INTENT(IN   ) ::                            &
                                  ids,ide, jds,jde, kds,kde, &
                                  ims,ime, jms,jme, kms,kme, &
                                  its,ite, jts,jte, kts,kte

   INTEGER,      INTENT(IN   ) :: STEPCU
   LOGICAL,      INTENT(IN   ) :: warm_rain

   REAL,         INTENT(IN   ) :: XLV0,XLV1,XLS0,XLS1
   REAL,         INTENT(IN   ) :: CP,R,G,EP1,EP2
   REAL,         INTENT(IN   ) :: SVP1,SVP2,SVP3,SVPT0

   INTEGER,      INTENT(IN   ) :: KTAU           

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(IN   ) ::                                   &
                                                          U, &
                                                          V, &
                                                          W, &
                                                         TH, &
                                                          T, &
                                                         QV, &
                                                       dz8w, &
                                                       Pcps, &
                                                        rho, &
                                                         pi

  INTEGER,      INTENT(IN   ) :: no_src_types_cu !PSH/TWG
  INTEGER,      INTENT(IN   ) :: aercu_opt       !PSH/TWG
  REAL,         INTENT(IN   ) :: aercu_fct       !PSH/TWG
  REAL,  DIMENSION( ims:ime, kms:kme, jms:jme, no_src_types_cu), OPTIONAL, &
          INTENT(INOUT) ::                                   aerocu !PSH/TWG
!
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(INOUT) ::                                   &
                                                      W0AVG   

   REAL,  INTENT(IN   ) :: DT, DX
   REAL,  INTENT(IN   ) :: CUDT
   LOGICAL,OPTIONAL,INTENT(IN   ) :: ADAPT_STEP_FLAG
!
   REAL, DIMENSION( ims:ime , jms:jme ),                     &
          INTENT(INOUT) ::                           RAINCV

   REAL,    DIMENSION( ims:ime , jms:jme ),                  &
          INTENT(INOUT) ::                           PRATEC

   REAL,    DIMENSION( ims:ime , jms:jme ),                  &
            INTENT(INOUT) ::                            NCA

   REAL, DIMENSION( ims:ime , jms:jme ),                     &
          INTENT(OUT) ::                              CUBOT, &
                                                      CUTOP    

   LOGICAL, DIMENSION( ims:ime , jms:jme ),                  &
          INTENT(INOUT) :: CU_ACT_FLAG
!
! Optional arguments
!
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),           &
         INTENT(INOUT) ::                                    &
                                                   RTHCUTEN, &
                                                   RQVCUTEN, &
                                                   RQCCUTEN, &
                                                   RQRCUTEN, &
                                                   RQICUTEN, &
                                                   RQSCUTEN, &
                                                   RQVFTEN,  &
                                                   RUCUTEN,  & !JTR
                                                   RVCUTEN
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, OPTIONAL ::                                      &
                                                   F_QV      &
                                                  ,F_QC      &
                                                  ,F_QR      &
                                                  ,F_QI      &
                                                  ,F_QS

!ckay
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(INOUT) ::                                   &
                                               cldfra_dp_KF, &
                                               cldfra_sh_KF, &
                                                      qc_KF, &
                                                      qi_KF

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(INOUT) ::                                   &
                                                      qr_KF, & ! TWG
                                                      qs_KF, & ! TWG
                                                      nc_KF, & ! TWG
                                                      ni_KF, & ! TWG
                                                      nr_KF, & ! TWG
                                                      ns_KF, & ! TWG
                                                     ccn_KF, & ! TWG
                                                       EFCS, & ! TWG
                                                       EFIS, & ! TWG
                                                       EFSS

   REAL, DIMENSION( ims:ime , jms:jme ),                     & !TWG
          INTENT(INOUT) ::                           ainc_frac

!kf_edrates
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         , &
          INTENT(INOUT) ::                         &
                                                     UDR_KF, &
                                                     DDR_KF, &
                                                     UER_KF, &
                                                     DER_KF

   REAL,  DIMENSION( ims:ime , jms:jme )                   , &
          INTENT(INOUT) ::                         &
                                                   TIMEC_KF

   INTEGER, INTENT(IN) ::              KF_EDRATES

!ckay
   REAL, DIMENSION( ims:ime, jms:jme )                     , &
         INTENT(   IN) ::                               ZOL, &
                                                        HFX, &
                                                        UST, &
                                                       PBLH, &
                                                       XLAND
!ckaywup
   REAL, DIMENSION( ims:ime,  kms:kme , jms:jme )          , &
           INTENT(INOUT) ::                             w_up

! LOCAL VARS

   LOGICAL :: flag_qr, flag_qi, flag_qs

   REAL, DIMENSION( kts:kte ) ::                             &
                                                        U1D, &
                                                        V1D, &
                                                        T1D, &
                                                       DZ1D, &
                                                       QV1D, &
                                                        P1D, &
                                                      RHO1D, &
                                                  tpart_v1D, &
                                                  tpart_h1D, &
                                                    W0AVG1D

   REAL, DIMENSION( kts:kte )::                              &
                                                       DQDT, &
                                                      DQIDT, &
                                                      DQCDT, &
                                                      DQRDT, &
                                                      DQSDT, &
                                                       DTDT

  INTEGER :: ii,jj,kk

  REAL, DIMENSION (kts:kte)  :: z0

   REAL    ::         TST,tv,PRS,RHOE,W0,SCR1,DXSQ,tmp
   integer :: ibegh,iendh,jbegh,jendh
   integer :: istart,iend,jstart,jend
   INTEGER :: i,j,k,NTST,ISP
   REAL    :: W0AVGfctr, W0fctr, W0den
   
!JTR 06/26/19: Added tendency variables and CMT flag
   REAL, DIMENSION( kts:kte ) :: DUDT, DVDT
   LOGICAL :: cmt_opt_flag
!local arrays for CCPP formatting
   REAL :: ZOL1,HFX1,UST1,PBLH1,NCA1,AINC_FRAC1,XLAND1,DYNTAU
   REAL :: RAINCV1,PRATEC1,CUTOP1,CUBOT1
   REAL, DIMENSION (kts:kte,NO_SRC_TYPES_CU)  :: aerocu2d
   REAL, DIMENSION (kts:kte)  :: dp_cf1d,sh_cf1d,w_up1d,QLIQ,QICE
   REAL, DIMENSION (kts:kte)  :: qrain,qsnow,nLIQ,nice,nrain,nsnow,ccn
   REAL, DIMENSION (kts:kte)  :: EFFCH,EFFIH,EFFSH,UDR,DDR,UER,DER   

!JTR: CMT on by default
   cmt_opt_flag = .TRUE.
!
   DXSQ=DX*DX

!----------------------
   NTST=STEPCU
   TST=float(NTST*2)
   flag_qr = .FALSE.
   flag_qi = .FALSE.
   flag_qs = .FALSE.
   IF ( PRESENT(F_QR) ) flag_qr = F_QR
   IF ( PRESENT(F_QI) ) flag_qi = F_QI
   IF ( PRESENT(F_QS) ) flag_qs = F_QS
!
   if (ADAPT_STEP_FLAG) then
      W0AVGfctr = 2 * MAX(CUDT*60,dt) - dt
      W0fctr = dt
      W0den = 2 * MAX(CUDT*60,dt)
   else
      W0AVGfctr = (TST-1.)
      W0fctr = 1.
      W0den = TST
   endif

  DO J = jts,jte
      DO K=kts,kte
         DO I= its,ite
!            SCR1=-5.0E-4*G*rho(I,K,J)*(w(I,K,J)+w(I,K+1,J))
!            TV=T(I,K,J)*(1.+EP1*QV(I,K,J))
!            RHOE=Pcps(I,K,J)/(R*TV)
!            W0=-101.9368*SCR1/RHOE
            W0=0.5*(w(I,K,J)+w(I,K+1,J))

!           Old:            
!
!            W0AVG(I,K,J)=(W0AVG(I,K,J)*(TST-1.)+W0)/TST            
!
!           New, to support adaptive time step:
!
            W0AVG(I,K,J) = ( W0AVG(I,K,J) * W0AVGfctr + W0 * W0fctr ) / W0den
!ckaywup
!           w(I,K,J)=w(I,K,J)+w_up(i,K,j)  

         ENDDO
      ENDDO
   ENDDO
!
     DO J = jts,jte
     DO I= its,ite
        CU_ACT_FLAG(i,j) = .true.
     ENDDO
     ENDDO

     DO J = jts,jte
       DO I=its,ite
          

         IF ( NCA(I,J) .ge. 0.5*DT ) then
            CU_ACT_FLAG(i,j) = .false.
         ELSE

            DO k=kts,kte
               DQDT(k)=0.
               DQIDT(k)=0.
               DQCDT(k)=0.
               DQRDT(k)=0.
               DQSDT(k)=0.
               DTDT(k)=0.
               DUDT(k)=0.
               DVDT(k)=0.
!ckay
               cldfra_dp_KF(I,k,J)=0.
               cldfra_sh_KF(I,k,J)=0.
               qc_KF(I,k,J)=0.
               qi_KF(I,k,J)=0.
             IF (aercu_opt.gt.0) THEN
               qr_KF(I,k,J)=0.
               qs_KF(I,k,J)=0.
               nc_KF(I,k,J)=0.
               ni_KF(I,k,J)=0.
               nr_KF(I,k,J)=0.
               ns_KF(I,k,J)=0.
               ccn_KF(I,k,J)=0.
               EFSS(I,k,J)=10.01
               EFCS(I,k,J)=2.51
               EFIS(I,k,J)=5.01
              END IF
               w_up(I,k,J)=0.
            ENDDO
             IF (aercu_opt.gt.0) THEN
               ainc_frac(I,J) = 0. ! TWG
             END IF
            IF (KF_EDRATES == 1) THEN
               DO k=kts,kte
                  UDR_KF(I,k,J)=0.
                  DDR_KF(I,k,J)=0.
                  UER_KF(I,k,J)=0.
                  DER_KF(I,k,J)=0.
               ENDDO
               TIMEC_KF(I,J)=0.
            ENDIF
            RAINCV(I,J)=0.
            CUTOP(I,J)=KTS
            CUBOT(I,J)=KTE+1
            PRATEC(I,J)=0.
!
! assign vars from 3D to 1D

            DO K=kts,kte
               U1D(K) =U(I,K,J)
               V1D(K) =V(I,K,J)
               T1D(K) =T(I,K,J)
               RHO1D(K) =rho(I,K,J)
               QV1D(K)=QV(I,K,J)
               P1D(K) =Pcps(I,K,J)
               W0AVG1D(K) =W0AVG(I,K,J)
               DZ1D(k)=dz8w(I,K,J)
            ENDDO
!
!dkay
            IF (aercu_opt.gt.0) THEN
                call mskf_mphyi ()

               DO K=KTS,KTE
                qr_KF(I,K,J)=0.
                qs_KF(I,K,J)=0.
                nc_KF(I,K,J)=0.
                ni_KF(I,K,J)=0.
                nr_KF(I,K,J)=0.
                ns_KF(I,K,J)=0.
                ccn_KF(I,K,J)=0.
                EFCS(I,K,J)=2.51
                EFIS(I,K,J)=5.01
                EFSS(I,K,J)=10.01

                 DO ISP=1,NO_SRC_TYPES_CU
                  aerocu2d(k,ISP) = aerocu(i,k,j,isp)
                 END DO
               END DO
            END IF
!ckay
! prepare inputs to CCPP version of the mskf_eta_para and fill in all 
! one-D (spatial) from 2D & 3D fields and other dimensions 
!
! first, fill in single value variables
            ZOL1       = ZOL(I,J)
            HFX1       = HFX(I,J)
            UST1       = UST(I,J)
            PBLH1      = PBLH(I,J)
            NCA1       = NCA(I,J) 
            AINC_FRAC1 = AINC_FRAC(I,J)
            XLAND1     = XLAND(I,J)
            RAINCV1    = 0.0
            PRATEC1    = 0.0
            CUTOP1     = KTS
            CUBOT1     = KTE+1 

            DO K=KTS,KTE
            IF (KF_EDRATES == 1) THEN
              UDR(K)=0.0
              DDR(K)=0.0
              UER(K)=0.0
              DER(K)=0.0
              DYNTAU = 0.0
            END IF
              dp_cf1d(K)=0.0
              sh_cf1d(K)=0.0
              w_up1d(K)=0.0
              QLIQ(K)=0.0
              QICE(K)=0.0
              DUDT(K)=0.0
              DVDT(K)=0.0
            END DO
              
            CALL MSKF_eta_PARA(                       &
                 U1D,V1D,T1D,QV1D,P1D,DZ1D,W0AVG1D, &
                 DT,DX,DXSQ,RHO1D,                  &
                 XLV0,XLV1,XLS0,XLS1,CP,R,G,        &
                 EP2,SVP1,SVP2,SVP3,SVPT0,          &
                 DQDT,DQIDT,DQCDT,DQRDT,DQSDT,DTDT, &
                 RAINCV1,PRATEC1,NCA1,              &
                 flag_QI,flag_QS,warm_rain,         &
                 CUTOP1,CUBOT1,CUDT,                  &
                                   kts,kte,         &
!ckay
                 dp_cf1d,sh_cf1d,w_up1d,    &
                 QLIQ,qice,qrain,qsnow,           &
                 nLIQ,nice,nrain,nsnow,ccn,       & !TWG
                 ainc_frac1,                         & !TWG
!kf_edrates
                 UDR,DDR,                     &
                 UER,DER,                     &
                 DYNTAU,KF_EDRATES,               &                 
                 ZOL1,HFX1,UST1,PBLH1,              &
                 aerocu2d,no_src_types_cu,aercu_fct,  &
                 aercu_opt,EFFCH,EFFIH,EFFSH,          & !PSH/TWG
                 DUDT, DVDT, cmt_opt_flag,XLAND1)      !JTR
!
!ckay fill in needed arrays back to original arrays/variables
                 NCA(I,J)       = NCA1
                 AINC_FRAC(I,J) = AINC_FRAC1
                 RAINCV(i,j)    = RAINCV1
                 PRATEC(I,J)    = PRATEC1
                 CUTOP(I,J)     = CUTOP1
                 CUBOT(I,J)     = CUBOT1

                  DO K=KTS,KTE
                    cldfra_dp_KF(I,K,J)=dp_cf1d(K)
                    cldfra_sh_KF(I,K,J)=sh_cf1d(K)
                    w_up(I,K,J)        =w_up1d(K)
                    qc_KF(I,K,J)       =QLIQ(K)
                    qi_KF(I,K,J)       =QICE(K)

                    IF (KF_EDRATES == 1) THEN
                       TIMEC_KF(I,J)  = DynTau
                         UDR_KF(I,K,J)= UDR(K)
                         DDR_KF(I,K,J)= DDR(K)
                         UER_KF(I,K,J)= UER(K)
                         DER_KF(I,K,J)= DER(K)
                    ENDIF

                    IF(aercu_opt.GT.0) then
                      qr_KF(I,K,J)=QRAIN(K)
                      qs_KF(I,K,J)=QSNOW(K)
                      nc_KF(I,K,J)=NLIQ(K)
                      ni_KF(I,K,J)=NICE(K)
                      nr_KF(I,K,J)=NRAIN(K)
                      ns_KF(I,K,J)=NSNOW(K)
                      ccn_KF(I,K,J)=CCN(K)
                      EFCS(I,K,J)=EFFCH(K)
                      EFIS(I,K,J)=EFFIH(K)
                      EFSS(I,K,J)=EFFSH(K)
                     END IF
                  END DO
            RAINCV1    = 0.0
            PRATEC1    = 0.0
            CUTOP1     = KTS
            CUBOT1     = KTE+1

!JTR: Pass 1D tendency arrays to 3D arrays              
              IF(cmt_opt_flag) THEN
                 DO K=kts,kte
                    RUCUTEN(I,K,J) = DUDT(K)
                    RVCUTEN(I,K,J) = DVDT(K)
                 ENDDO
              ENDIF
              DO K=kts,kte
                 RTHCUTEN(I,K,J)=DTDT(K)/pi(I,K,J)
                 RQVCUTEN(I,K,J)=DQDT(K)
              ENDDO

              IF( F_QR )THEN
                DO K=kts,kte
                   RQRCUTEN(I,K,J)=DQRDT(K)
                   RQCCUTEN(I,K,J)=DQCDT(K)
                ENDDO
              ELSE
! This is the case for Eta microphysics without 3d rain field
                DO K=kts,kte
                   RQRCUTEN(I,K,J)=0.
                   RQCCUTEN(I,K,J)=DQRDT(K)+DQCDT(K)
                ENDDO
              ENDIF

!......     QSTEN STORES GRAUPEL TENDENCY IF IT EXISTS, OTHERISE SNOW (V2)

              IF ( F_QI ) THEN
                DO K=kts,kte
                   RQICUTEN(I,K,J)=DQIDT(K)
                ENDDO
              ENDIF

              IF ( F_QS ) THEN
                DO K=kts,kte
                   RQSCUTEN(I,K,J)=DQSDT(K)
                ENDDO
              ENDIF
!
         ENDIF  ! for cu acct flag
       ENDDO     ! i-loop
     ENDDO       ! j-loop
!
   END SUBROUTINE MSKF_CPS
! ****************************************************************************
! insert here the ccpp mskf_etap_para codes here 

   SUBROUTINE MSKF_eta_PARA (                              &
                      U0,V0,T0,QV0,P0,DZQ,W0AVG1D,         &
                      DT,DX,DXSQ,rhoe,                     &
                      XLV0,XLV1,XLS0,XLS1,CP,R,G,          &
                      EP2,SVP1,SVP2,SVP3,SVPT0,            &
                      DQDT,DQIDT,DQCDT,DQRDT,DQSDT,DTDT,   &
                      RAINCV1,PRATEC1,NCA1,                   &
                      F_QI,F_QS,warm_rain,                 &
                      CUTOP1,CUBOT1,CUDT,                    &
                                        kts,kte,           &
!ckay
                      dp_cf1d,sh_cf1d,W_up1d,      &
                      QLIQ,QICE,QRAIN,QSNOW,             & !TWG
                      NLIQ,NICE,NRAIN,NSNOW,CCN,         & !TWG
                      ainc_frac1,                           & !TWG
!kf_edrates
                      UDR,DDR,                       &
                      UER,DER,                       &
                      DynTau,KF_EDRATES,                 &
                      ZOL1,HFX1,UST1,PBLH1,                &
                      aerocu2d,no_src_types_cu,aercu_fct,    &
                      aercu_opt,EFFCH,EFFIH,EFFSH,         & !PSH/TWG
                      DUDT,DVDT,cmt_opt_flag,XLAND1)          !JTR

!-----------------------------------------------------------
!***** The KF scheme that is currently used in experimental runs of EMCs 
!***** Eta model....jsk 8/00
!
      IMPLICIT NONE
!     SAVE  !TWG 2017 Add to avoid memory issues
!-----------------------------------------------------------
      INTEGER, INTENT(IN   ) :: kts,kte   

      LOGICAL, INTENT(IN   ) :: F_QI, F_QS

      LOGICAL, INTENT(IN   ) :: warm_rain
!
      REAL, DIMENSION( kts:kte ),                          &
            INTENT(IN   ) ::                           U0, &
                                                       V0, &
                                                       T0, &
                                                      QV0, &
                                                       P0, &
                                                     rhoe, &
                                                      DZQ, &
                                                  W0AVG1D
!
      REAL,  INTENT(IN   ) :: DT,DX,DXSQ
!

      REAL,  INTENT(IN   ) :: XLV0,XLV1,XLS0,XLS1,CP,R,G
      REAL,  INTENT(IN   ) :: EP2,SVP1,SVP2,SVP3,SVPT0

      INTEGER, INTENT(IN   ) :: no_src_types_cu  !PSH/TWG
      REAL,    INTENT(IN   ) :: aercu_fct        !PSH/TWG
      INTEGER, INTENT(IN   ) :: aercu_opt        !PSH/TWG
      REAL,  DIMENSION( kts:kte, no_src_types_cu), OPTIONAL, &
      INTENT(INOUT) ::                             aerocu2d !PSH/TWG


!ckay
      REAL,                                                &
            INTENT(   IN) ::                          ZOL1, &
                                                      HFX1, &
                                                      UST1, &
                                                     PBLH1, &
                                                     XLAND1
!
      REAL, DIMENSION( kts:kte ), INTENT(INOUT) ::         &
                                                     DQDT, &
                                                    DQIDT, &
                                                    DQCDT, &
                                                    DQRDT, &
                                                    DQSDT, &
                                                     DTDT

      REAL,                                                &
            INTENT(INOUT) ::                          NCA1

      REAL,                                                & !TWG
            INTENT(INOUT) ::                      ainc_frac1

      REAL, DIMENSION(           kts:kte           ),      &
            INTENT(INOUT) ::                 dp_cf1d, &
                                             sh_cf1d


      REAL, INTENT(INOUT) ::     DynTau,RAINCV1,PRATEC1

      INTEGER, INTENT(IN) ::         KF_EDRATES

      REAL, INTENT(OUT) ::                          CUBOT1, &
                                                    CUTOP1
      REAL,  INTENT(IN   ) :: CUDT

!ckaywup
  REAL, DIMENSION(           kts:kte           )         , &
        INTENT(  OUT) ::                             W_up1d                       
!
!...DEFINE LOCAL VARIABLES...
!
      REAL, DIMENSION( kts:kte ) ::                        &
            Q0,Z0,TV0,TU,TVU,QU,TZ,TVD,                    &
            QD,QES,THTES,TG,TVG,QG,WU,WD,W0,EMS,EMSD,      &
            UMF,UER,UDR,DMF,DER,DDR,UMF2,UER2,             &
            UDR2,DMF2,DER2,DDR2,DZA,THTA0,THETEE,          &
            THTAU,THETEU,THTAD,THETED,QLIQ,QICE,           &
!TWG 06/14/16
            QRAIN,QSNOW,NLIQ,NICE,NRAIN,NSNOW,CCN,         &
            EFFCH,EFFIH,EFFSH, &

!EBD TWG
            QLQOUT,QICOUT,PPTLIQ,PPTICE,DETLQ,DETIC,       &
            DETLQ2,DETIC2,RATIO,RATIO2


      REAL, DIMENSION( kts:kte ) ::                        &
            DOMGDP,EXN,TVQU,DP,RH,EQFRC,WSPD,              &
            QDT,FXM,THTAG,THPA,THFXOUT,                    &
            THFXIN,QPA,QFXOUT,QFXIN,QLPA,QLFXIN,           &
            QLFXOUT,QIPA,QIFXIN,QIFXOUT,QRPA,              &
            QRFXIN,QRFXOUT,QSPA,QSFXIN,QSFXOUT,            &
            QL0,QLG,QI0,QIG,QR0,QRG,QS0,QSG


      REAL, DIMENSION( kts:kte+1 ) :: OMG
      REAL, DIMENSION( kts:kte ) :: RAINFB,SNOWFB
      REAL, DIMENSION( kts:kte ) ::                        &
            CLDHGT,QSD,DILFRC,DDILFRC,TKE,TGU,QGU,THTEEG

! LOCAL VARS

      REAL    :: P00,T00,RLF,RHIC,RHBC,PIE,         &
                 TTFRZ,TBFRZ,C5,RATE
      REAL    :: GDRY,ROCP,ALIQ,BLIQ,                      &
                 CLIQ,DLIQ
      REAL    :: FBFRC,P300,DPTHMX,THMIX,QMIX,ZMIX,PMIX,   &
                 ROCPQ,TMIX,EMIX,TLOG,TDPT,TLCL,TVLCL,     &
                 CPORQ,PLCL,ES,DLP,TENV,QENV,TVEN,TVBAR,   &
                 ZLCL,WKL,WABS,TRPPT,WSIGNE,DTLCL,GDT,WLCL,&
                 TVAVG,QESE,WTW,RHOLCL,AU0,VMFLCL,UPOLD,   &
                 UPNEW,ABE,WKLCL,TTEMP,FRC1,   &
                 QNEWIC,RL,R1,QNWFRZ,EFFQ,BE,BOTERM,ENTERM,&
                 DZZ,UDLBE,REI,EE2,UD2,TTMP,F1,F2,         &
                 THTTMP,QTMP,TMPLIQ,TMPICE,TU95,TU10,EE1,  &
                 UD1,DPTT,QNEWLQ,DUMFDP,EE,TSAT,           &
                 THTA,VCONV,TIMEC,SHSIGN,VWS,PEF, &
                 CBH,RCBH,PEFCBH,PEFF,PEFF2,TDER,THTMIN,   &
                 DTMLTD,QS,TADVEC,DPDD,FRC,DPT,RDD,A1,     &
                 DSSDT,DTMP,T1RH,QSRH,PPTFLX,CPR,CNDTNF,   &
                 UPDINC,AINCM2,DEVDMF,PPR,RCED,DPPTDF,     &
                 DMFLFS,DMFLFS2,RCED2,DDINC,AINCMX,AINCM1, &
                 AINC,TDER2,PPTFL2,FABE,STAB,DTT,DTT1,     &
                 DTIME,TMA,TMB,TMM,BCOEFF,ACOEFF,QVDIFF,   &
                 TOPOMG,CPM,DQ,ABEG,DABE,DFDA,FRC2,DR,     &
                 UDFRC,TUC,QGS,RH0,RHG,QINIT,QFNL,ERR2,    &
                 RELERR,RLC,RLS,RNC,FABEOLD,AINCOLD,UEFRC, &
                 DDFRC,TDC,DEFRC,RHBAR,DMFFRC,DPMIN,DILBE
   REAL    ::    ASTRT,TP,VALUE,AINTRP,TKEMAX,QFRZ,&
                 QSS,PPTMLT,DTMELT,RHH,EVAC,BINC
!
      INTEGER :: INDLU,NU,NUCHM,NNN,KLFS
   REAL    :: CHMIN,PM15,CHMAX,DTRH,RAD,DPPP
   REAL    :: TVDIFF,DTTOT,ABSOMG,ABSOMGTC,FRDP
!ckay
   REAL    :: xcldfra,UMF_new,DMF_new,FXM_new
   REAL    :: sourceht, Scale_Fac, TOKIOKA, RATE_kay
   REAL    :: capeDX, tempKay
   REAL    :: SCLvel, ZLCL_KAY, zz_kay

!ckaywup
   REAL    :: envEsat, envQsat, envRH, envRHavg, denSplume
   REAL    :: updil, Drag, WST, thetav
!TWG Mar 2017
   REAL, PARAMETER :: P1_HU10 = 7.6725
   REAL, PARAMETER :: P2_HU10 = 1.0118
   REAL, PARAMETER :: P3_HU10 = 0.1422
   REAL, PARAMETER :: P4_HU10 = 0.0106
   REAL, PARAMETER :: P5_HU10 = 3.39E-4
   REAL, PARAMETER :: P6_HU10 = 3.95E-6
   REAL    :: SF_HU10, TC_HU10
!END TWG
!dkay for dccmp
   real :: a1kay
   LOGICAL :: DCCMP
   REAL :: eps1u, alatent, Qsu
   LOGICAL :: onetime
   Data onetime/.true./
   integer, parameter :: naer_cu = 10
   integer, parameter :: pcols = 1
   REAL(r8) muu(pcols, KTS:KTE)
   REAL(r8) su(pcols, KTS:KTE)
   REAL(r8) quu(pcols, KTS:KTE)
   REAL(r8) duu(pcols, KTS:KTE)
   REAL(r8) euu(pcols, KTS:KTE)
   REAL(r8) cmel(pcols, KTS:KTE)
   REAL(r8) cmei(pcols, KTS:KTE)
   REAL(r8) zfu(pcols, KTS:KTE+1)
   REAL(r8) zf_wrf(0:KTE)
   REAL(r8) pru(pcols, KTS:KTE)
   REAL(r8) tee(pcols, KTS:KTE)
   REAL(r8) qee(pcols, KTS:KTE)
   REAL(r8) qsatzm(pcols, KTS:KTE)
   REAL(r8) gamhat(pcols, KTS:KTE)
   REAL(r8) aer_mmr(pcols, KTS:KTE,naer_cu)
   REAL(r8) Aqnewic(KTS:KTE)
   REAL(r8) Aqnewlq(KTS:KTE)
   REAL(r8) wu_mskf_act(KTS:KTE)
   REAL(r8) qc_mskf_act(KTS:KTE)
   REAL(r8) qi_mskf_act(KTS:KTE)
   REAL(r8) effc(pcols, KTS:KTE)
   REAL(r8) effi(pcols, KTS:KTE)
   REAL(r8) effs(pcols, KTS:KTE)
   real(r8) QSATu(KTS:KTE), oldQU(KTS:KTE),oldTU(KTS:KTE)      ! rate of freezing
   REAL(r8) EPSI0(pcols)
   REAL(r8) dLfmzmp(pcols,KTS:KTE),dIfmzmp(pcols,KTS:KTE)
!junk
   REAL(r8) oldpptliq(KTS:KTE)
   REAL(r8) oldpptice(KTS:KTE)

   REAL(r8) wump(pcols, KTS:KTE)
   real(r8) zmqliq(pcols,KTS:KTE)       ! cloud water mixing ratio (kg/kg)
   real(r8) zmqice(pcols,KTS:KTE)       ! cloud ice mixing ratio (kg/kg)
   real(r8) zmqrain(pcols,KTS:KTE)      ! rain mixing ratio (kg/kg) !TWG
   real(r8) zmqsnow(pcols,KTS:KTE)      ! snow mixing ratio (kg/kg) !TWG
   real(r8) ncmp(pcols,KTS:KTE)       ! cloud water number conc (1/kg)
   real(r8) nimp(pcols,KTS:KTE)       ! cloud ice number conc (1/kg)
   real(r8) nrmp(pcols,KTS:KTE)       ! rain number conc (1/kg) !TWG
   real(r8) nsmp(pcols,KTS:KTE)       ! snow number conc (1/kg) !TWG
   real(r8) zmccn(pcols,KTS:KTE)       ! ccn conc (1/kg) !TWG
   real(r8) rprd(pcols,KTS:KTE)     ! rate of production of precip at that layer
   real(r8) sprd(pcols,KTS:KTE)     ! rate of production of snow at that layer
   real(r8) frz(pcols,KTS:KTE)      ! rate of freezing

   REAL(r8) grav, Rdry , DTZMP, CPIN, psh_fac

   Integer KQ, JK, JBB(1), JTT(1), JLCL(1), msg1, il2g , JZM, KA
   Integer NLEVZM, NLEVZMP1, KKAY, Miter, Itest, KC
!dkay

      INTEGER :: KX,K,KL
!
      INTEGER :: NCHECK
      INTEGER, DIMENSION (kts:kte) :: KCHECK

      INTEGER :: ISTOP,ML,L5,KMIX,LOW,                     &
                 LC,MXLAYR,LLFC,NLAYRS,NK,                 &
                 KPBL,KLCL,LCL,LET,IFLAG,                  &
                 NK1,LTOP,NJ,LTOP1,                        &
                 LTOPM1,LVF,KSTART,KMIN,LFS,               &
                 ND,NIC,LDB,LDT,ND1,NDK,                   &
                 NM,LMAX,NCOUNT,NOITR,                     &
                 NSTEP,NTC,NCHM,ISHALL,NSHALL
      LOGICAL :: IPRNT
      REAL :: u00,qslcl,rhlcl,dqssdt    !jfb
      CHARACTER*1024 message

!JTR 06/18/2019: Variables needed for CMT
      REAL, DIMENSION (kts:kte), INTENT(INOUT) :: DUDT, DVDT
      REAL, DIMENSION (kts:kte) :: TDN,TUP
      REAL, DIMENSION (kts:kte) :: stat_energy
      REAL(r8), DIMENSION (pcols,kts:kte) :: DUDTnew, DVDTnew,   &
                                             DPDX, DPDY, MC,     &
                                             DERconvF, UERconvF, &
                                             UMFconvF, DMFconvF, &
                                             DPconvF, U0F, V0F,  &
                                             Z0F, SUU, SDD,      &
                                             SHAT, QHAT, QDN,QUP
      Logical :: CMTprint
      Data CMTprint/.false./
      INTEGER :: JDD(1), IL1G, ILG
      REAL(r8) :: DSUBCLD(1), VMFLCLconv(1), DTnew
      LOGICAL :: cmt_opt_flag

!
      DATA P00,T00/1.E5,273.16/
      DATA RLF/3.339E5/
      DATA RHIC,RHBC/1.,0.90/
      DATA PIE,TTFRZ,TBFRZ,C5/3.141592654,268.16,248.16,1.0723E-3/
      DATA RATE/0.03/   ! wrf default
!      DATA RATE/0.01/  ! value used in NRCM
!      DATA RATE/0.001/  ! effectively turn off autoconversion
!-----------------------------------------------------------
   IF (aercu_opt.gt.0) THEN
     DCCMP = .TRUE.
   ELSE
     DCCMP = .FALSE.
   END IF

      IPRNT=.FALSE.
      GDRY=-G/CP
      ROCP=R/CP
      NSHALL = 0
      KL=kte
      KX=kte
!
!     ALIQ = 613.3
!     BLIQ = 17.502
!     CLIQ = 4780.8
!     DLIQ = 32.19
      ALIQ = SVP1*1000.
      BLIQ = SVP2
      CLIQ = SVP2*SVPT0
      DLIQ = SVP3
!
      IF(DX.GE.24.999E3) THEN
         Scale_Fac = 1.0
         capeDX = 0.1
      ELSE
         Scale_Fac = 1.0 + (log(25.E3/DX))
         capeDX = 0.1 *SQRT(Scale_Fac)
      END IF
!
!****************************************************************************
!                                                      ! PPT FB MODS
!...OPTION TO FEED CONVECTIVELY GENERATED RAINWATER    ! PPT FB MODS
!...INTO GRID-RESOLVED RAINWATER (OR SNOW/GRAUPEL)     ! PPT FB MODS
!...FIELD.  "FBFRC" IS THE FRACTION OF AVAILABLE       ! PPT FB MODS
!...PRECIPITATION TO BE FED BACK (0.0 - 1.0)...        ! PPT FB MODS
      FBFRC=0.0                                        ! PPT FB MODS
!...mods to allow shallow convection...
      NCHM = 0
      ISHALL = 0
      DPMIN = 5.E3
!...
      P300=P0(1)-30000.
!
!...PRESSURE PERTURBATION TERM IS ONLY DEFINED AT MID-POINT OF
!...VERTICAL LAYERS...SINCE TOTAL PRESSURE IS NEEDED AT THE TOP AND
!...BOTTOM OF LAYERS BELOW, DO AN INTERPOLATION...
!
!...INPUT A VERTICAL SOUNDING ... NOTE THAT MODEL LAYERS ARE NUMBERED
!...FROM BOTTOM-UP IN THE KF SCHEME...
!
      ML=0 
!SUE  tmprpsb=1./PSB(I,J)
!SUE  CELL=PTOP*tmprpsb
!
      DO K=1,KX
!
!  Saturation vapor pressure (ES) is calculated following Buck (1981)
!...IF Q0 IS ABOVE SATURATION VALUE, REDUCE IT TO SATURATION LEVEL...
!
         ES=ALIQ*EXP((BLIQ*T0(K)-CLIQ)/(T0(K)-DLIQ))
         QES(K)=0.622*ES/(P0(K)-ES)
         Q0(K)=AMIN1(QES(K),QV0(K))
         Q0(K)=AMAX1(0.000001,Q0(K))
         QL0(K)=0.
         QI0(K)=0.
         QR0(K)=0.
         QS0(K)=0.
         RH(K) = Q0(K)/QES(K)
         DILFRC(K) = 1.
         TV0(K)=T0(K)*(1.+0.608*Q0(K))
!        RHOE(K)=P0(K)/(R*TV0(K))
!   DP IS THE PRESSURE INTERVAL BETWEEN FULL SIGMA LEVELS...
         DP(K)=rhoe(k)*g*DZQ(k)
! IF Turbulent Kinetic Energy (TKE) is available from turbulent mixing scheme
! use it for shallow convection...For now, assume it is not available....
!         TKE(K) = Q2(I,J,NK)
         TKE(K) = 0.
         CLDHGT(K) = 0.
!        IF(P0(K).GE.500E2)L5=K
         IF(P0(K).GE.0.5*P0(1))L5=K
         IF(P0(K).GE.P300)LLFC=K
      ENDDO
!
!...DZQ IS DZ BETWEEN SIGMA SURFACES, DZA IS DZ BETWEEN MODEL HALF LEVEL
        Z0(1)=.5*DZQ(1)
!cdir novector
        DO K=2,KL
          Z0(K)=Z0(K-1)+.5*(DZQ(K)+DZQ(K-1))
          DZA(K-1)=Z0(K)-Z0(K-1)
        ENDDO   
        DZA(KL)=0.
!
!
!  To save time, specify a pressure interval to move up in sequential
!  check of different ~50 mb deep groups of adjacent model layers in
!  the process of identifying updraft source layer (USL).  Note that 
!  this search is terminated as soon as a buoyant parcel is found and 
!  this parcel can produce a cloud greater than specifed minimum depth
!  (CHMIN)...For now, set interval at 15 mb...
!
       NCHECK = 1
       KCHECK(NCHECK)=1
       PM15 = P0(1)-15.E2
       DO K=2,LLFC
         IF(P0(K).LT.PM15)THEN
           NCHECK = NCHECK+1
           KCHECK(NCHECK) = K
           PM15 = PM15-15.E2
         ENDIF
       ENDDO
!
       NU=0
       NUCHM=0
usl:   DO
           NU = NU+1
           IF(NU.GT.NCHECK)THEN 
             IF(ISHALL.EQ.1)THEN
               CHMAX = 0.
               NCHM = 0
               DO NK = 1,NCHECK
                 NNN=KCHECK(NK)
                 IF(CLDHGT(NNN).GT.CHMAX)THEN
                   NCHM = NNN
                   NUCHM = NK
                   CHMAX = CLDHGT(NNN)
                 ENDIF
               ENDDO
               NU = NUCHM-1
               FBFRC=1.
               CYCLE usl
             ELSE
               RETURN
             ENDIF
           ENDIF      
           KMIX = KCHECK(NU)
           LOW=KMIX
!...
           LC = LOW
!
!...ASSUME THAT IN ORDER TO SUPPORT A DEEP UPDRAFT YOU NEED A LAYER OF
!...UNSTABLE AIR AT LEAST 50 mb DEEP...TO APPROXIMATE THIS, ISOLATE A
!...GROUP OF ADJACENT INDIVIDUAL MODEL LAYERS, WITH THE BASE AT LEVEL
!...LC, SUCH THAT THE COMBINED DEPTH OF THESE LAYERS IS AT LEAST 50 mb..
!   
           NLAYRS=0
           DPTHMX=0.
           NK=LC-1
           IF ( NK+1 .LT. KTS ) THEN
             WRITE(message,*)'WOULD GO OFF BOTTOM: MSKF_PARA I,J,NK',NK
             CALL wrf_message (TRIM(message)) 
           ELSE
             DO 
               NK=NK+1   
               IF ( NK .GT. KTE ) THEN
                 WRITE(message,*)'WOULD GO OFF TOP: MSKF_PARA I,J,DPTHMX,DPMIN',DPTHMX,DPMIN
                 CALL wrf_message (TRIM(message))
                 EXIT
               ENDIF
               DPTHMX=DPTHMX+DP(NK)
               NLAYRS=NLAYRS+1
               IF(DPTHMX.GT.DPMIN)THEN
                 EXIT 
               ENDIF
             END DO    
           ENDIF
           IF(DPTHMX.LT.DPMIN)THEN 
             RETURN
           ENDIF
           KPBL=LC+NLAYRS-1   
!
!...********************************************************
!...for computational simplicity without much loss in accuracy,
!...mix temperature instead of theta for evaluating convective
!...initiation (triggering) potential...
!          THMIX=0.
           TMIX=0.
           QMIX=0.
           ZMIX=0.
           PMIX=0.
!
!...FIND THE THERMODYNAMIC CHARACTERISTICS OF THE LAYER BY
!...MASS-WEIGHTING THE CHARACTERISTICS OF THE INDIVIDUAL MODEL
!...LAYERS...
!
!cdir novector
           DO NK=LC,KPBL
             TMIX=TMIX+DP(NK)*T0(NK)
             QMIX=QMIX+DP(NK)*Q0(NK)
             ZMIX=ZMIX+DP(NK)*Z0(NK)
             PMIX=PMIX+DP(NK)*P0(NK)
           ENDDO   
!         THMIX=THMIX/DPTHMX
          TMIX=TMIX/DPTHMX
          QMIX=QMIX/DPTHMX
          ZMIX=ZMIX/DPTHMX
          PMIX=PMIX/DPTHMX
          EMIX=QMIX*PMIX/(0.622+QMIX)
!
!...FIND THE TEMPERATURE OF THE MIXTURE AT ITS LCL...
!
!        TLOG=ALOG(EMIX/ALIQ)
! ...calculate dewpoint using lookup table...
!
          astrt=1.e-3
          ainc=0.075
          a1=emix/aliq
          tp=(a1-astrt)/ainc
          indlu=int(tp)+1
          value=(indlu-1)*ainc+astrt
          aintrp=(a1-value)/ainc
          tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
          TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
          TLCL=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(TMIX-T00))*(TMIX-TDPT)
          TLCL=AMIN1(TLCL,TMIX)
          TVLCL=TLCL*(1.+0.608*QMIX)
          ZLCL = ZMIX+(TLCL-TMIX)/GDRY
     !     NK = LC-1
     !     DO 
     !       NK = NK+1
     !       KLCL=NK
     !       IF(ZLCL.LE.Z0(NK) .or. NK.GT.KL)THEN
     !         EXIT
     !       ENDIF 
     !     ENDDO   
     !     IF(NK.GT.KL)THEN
     !       RETURN  
     !     ENDIF

       DO NK = LC, KL
         KLCL = NK
         IF ( ZLCL.LE.Z0(NK) )  EXIT
     END DO
     IF ( ZLCL.GT.Z0(KL) )  RETURN

          K=KLCL-1
! calculate DLP using Z instead of log(P)
          DLP=(ZLCL-Z0(K))/(Z0(KLCL)-Z0(K))
!     
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
!     
          TENV=T0(K)+(T0(KLCL)-T0(K))*DLP
          QENV=Q0(K)+(Q0(KLCL)-Q0(K))*DLP
          TVEN=TENV*(1.+0.608*QENV)
!     
! ww: this needs to be initialized
          DTRH = 0.

! Bechtold 2001 trigger with my Beta parameter
           DTLCL = W0AVG1D(KLCL)/Scale_Fac
           if(DTLCL.lt.0.0) then
              tempKay = -1.0
              DTLCL = tempKay * DTLCL
              DTLCL = (DTLCL)**0.3333
           else
              tempKay = 1.0
              DTLCL = tempKay * DTLCL
              DTLCL = (DTLCL)**0.3333
           end if

           DTLCL = 6.0 * tempKay * DTLCL 
!
! old trigger
! Stick with the old trigger for now... CGM July 2015
!
          IF(ZLCL.LT.2.E3)THEN        ! Kain (2004) Eq. 2
            WKLCL=0.02*ZLCL/2.E3
          ELSE
            WKLCL=0.02                ! units of m/s
          ENDIF
!TWG.ckay c
         if(DX.GE.25.E3) then
           WKL=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)*DX/25.E3-WKLCL
         else
           WKL=(W0AVG1D(K)+(W0AVG1D(KLCL)-W0AVG1D(K))*DLP)-WKLCL
         end if 
!TWG ckay, Modified WKL
          IF(WKL.LT.0.0001)THEN
            DTLCL=0.
          ELSE
            DTLCL=4.64*WKL**0.33      ! Kain (2004) Eq. 1
          ENDIF


!         IF(ISHALL.EQ.1)IPRNT=.TRUE.
!         IPRNT=.TRUE.
!         IF(TLCL+DTLCL.GT.TENV)GOTO 45

           IF(TLCL+DTLCL.LT.TENV)THEN
!
! Parcel not buoyant, CYCLE back to start of trigger and evaluate next potential
! USL...
!
            CYCLE usl
!
          ELSE                            ! Parcel is buoyant, determine updraft
!     
!...CONVECTIVE TRIGGERING CRITERIA HAS BEEN SATISFIED...COMPUTE
!...EQUIVALENT POTENTIAL TEMPERATURE
!...(THETEU) AND VERTICAL VELOCITY OF THE RISING PARCEL AT THE LCL...
!     
            CALL ENVIRTHT(PMIX,TMIX,QMIX,THETEU(K),ALIQ,BLIQ,CLIQ,DLIQ)
!
!...modify calculation of initial parcel vertical velocity...jsk 11/26/97
!
            DTTOT = DTLCL+DTRH
            IF(DTTOT.GT.1.E-4)THEN
              GDT=2.*G*DTTOT*500./TVEN     ! Kain (2004) Eq. 3  (sort of)
              WLCL=1.+0.5*SQRT(GDT)
              WLCL = AMIN1(WLCL,3.)
            ELSE
              WLCL=1.
            ENDIF
            PLCL=P0(K)+(P0(KLCL)-P0(K))*DLP
            WTW=WLCL*WLCL
!
            TVLCL=TLCL*(1.+0.608*QMIX)
            RHOLCL=PLCL/(R*TVLCL)
!        
            LCL=KLCL
            LET=LCL
!ckay
! new formulation based on the LCL replacing the cloud radius concept
!introduce LCL instead of RAD based on WKL here
            RAD = ZLCL
!ckay Dec20
            sourceht = Z0(KPBL)
            RAD = amax1(sourceht, RAD)

            RAD = AMIN1(4000.,RAD)  ! max trap
            RAD = AMAX1(500.,RAD)  ! min trap
 
!     
!*******************************************************************
!                                                                  *
!                 COMPUTE UPDRAFT PROPERTIES                       *
!                                                                  *
!*******************************************************************
!     
!     
!...
!...ESTIMATE INITIAL UPDRAFT MASS FLUX (UMF(K))...
!     
            WU(K)=WLCL
            AU0=0.01*DXSQ
            UMF(K)=RHOLCL*AU0

            VMFLCL=UMF(K)
            UPOLD=VMFLCL
            UPNEW=UPOLD
!     
!...RATIO2 IS THE DEGREE OF GLACIATION IN THE CLOUD (0 TO 1),
!...UER IS THE ENVIR ENTRAINMENT RATE, ABE IS AVAILABLE
!...BUOYANT ENERGY, TRPPT IS THE TOTAL RATE OF PRECIPITATION
!...PRODUCTION...
!     
            RATIO2(K)=0.
            UER(K)=0.
            ABE=0.
            TRPPT=0.
            TU(K)=TLCL
            TVU(K)=TVLCL
            QU(K)=QMIX
            EQFRC(K)=1.
            QLIQ(K)=0.
            QICE(K)=0.
         IF (aercu_opt .GT. 0) THEN
            QRAIN(K)=0.
            QSNOW(K)=0.
            NLIQ(K)=0.
            NICE(K)=0.
            NRAIN(K)=0.
            NSNOW(K)=0.
            CCN(K)=0.
            EFFCH(K) = 2.5
            EFFIH(K) = 4.99
            EFFSH(K) = 9.99
         END IF
            QLQOUT(K)=0.
            QICOUT(K)=0.
            DETLQ(K)=0.
            DETIC(K)=0.
            PPTLIQ(K)=0.
            PPTICE(K)=0.
            IFLAG=0
!     
!...TTEMP IS USED DURING CALCULATION OF THE LINEAR GLACIATION
!...PROCESS; IT IS INITIALLY SET TO THE TEMPERATURE AT WHICH
!...FREEZING IS SPECIFIED TO BEGIN.  WITHIN THE GLACIATION
!...INTERVAL, IT IS SET EQUAL TO THE UPDRAFT TEMP AT THE
!...PREVIOUS MODEL LEVEL...
!     
            TTEMP=TTFRZ
!     
!...ENTER THE LOOP FOR UPDRAFT CALCULATIONS...CALCULATE UPDRAFT TEMP,
!...MIXING RATIO, VERTICAL MASS FLUX, LATERAL DETRAINMENT OF MASS AND
!...MOISTURE, PRECIPITATION RATES AT EACH MODEL LEVEL...
!     
!    **1 variables indicate the bottom of a model layer
!    **2 variables indicate the top of a model layer
!     
            EE1=1.
            UD1=0.
            REI = 0.
            DILBE = 0.

!dkay
      IF (aercu_opt.gt.0) THEN
                zf_wrf(0) = 0.0  ! ground
                DO KQ=KTS,KTE
                 zf_wrf(KQ) = zf_wrf(KQ-1)+DZQ(KQ)
                 Aqnewlq(kq) = 0.0
                 Aqnewic(kq) = 0.0
                 rprd(1,kq) = 0.0
                 wump(1,kq) =0.0
                 ncmp(1,kq) =0.0
                 nimp(1,kq) =0.0
                 sprd(1,kq) =0.0
                 frz(1,kq) =0.0
                 jk = kq
                 muu(1,JK) = 0.0
                 duu(1,JK) =0.0
                 EUU(1,JK) =0.0
                 cmel(1,JK) =0.0
                 cmei(1,JK) =0.0
                 oldTU(kq) = t0(kq)
                 oldQU(kq) = Q0(kq)
                End do
              Miter = 0
      END IF

updraft:    DO NK=K,KL-1
              NK1=NK+1
              RATIO2(NK1)=RATIO2(NK)
              FRC1=0.
              TU(NK1)=T0(NK1)
              THETEU(NK1)=THETEU(NK)
              QU(NK1)=QU(NK)
!dkay
     IF (aercu_opt.gt.0) THEN
              oldQU(NK) = QU(NK)
              oldTU(NK) = TU(NK)
     END IF
!dkay
              QLIQ(NK1)=QLIQ(NK)
              QICE(NK1)=QICE(NK)
              call tpmix2(p0(nk1),theteu(nk1),tu(nk1),qu(nk1),qliq(nk1),        &
                     qice(nk1),qnewlq,qnewic,XLV1,XLV0,QSu)

!dkay QSu has been added to the tpmix2

!dkay
! saturation value of Q of updraft for use with gamma hat in DCCMP routine
!   IF (aercu_opt.gt.0) THEN
!        QSATu(NK) = QSu/(1.+QSu)   ! saturated specific hum
!         Aqnewlq(NK) = qnewlq
!         Aqnewic(NK) = qnewic

!        Aqnewlq(NK) = qnewlq + Qliq(nk )   This is to be removed
!        Aqnewic(NK) = qnewic + Qice(nk )

!dkaydec26
     !    if(TU(NK).le.273.) then
     !     Aqnewlq(NK) = 0.0
     !     Aqnewic(NK) = qnewlq + qnewic
      !   else
       !   Aqnewlq(NK) = qnewlq + qnewic
        !  Aqnewic(NK) = 0.0
        ! end if
!   END IF
!
!
!...CHECK TO SEE IF UPDRAFT TEMP IS ABOVE THE TEMPERATURE AT WHICH
!/dec26...GLACIATION IS ASSUMED TO INITIATE; IF IT IS, CALCULATE THE
!...FRACTION OF REMAINING LIQUID WATER TO FREEZE...TTFRZ IS THE
!...TEMP AT WHICH FREEZING BEGINS, TBFRZ THE TEMP BELOW WHICH ALL
!...LIQUID WATER IS FROZEN AT EACH LEVEL...
!
              IF(TU(NK1).LE.TTFRZ)THEN
                IF(TU(NK1).GT.TBFRZ)THEN
                  IF(TTEMP.GT.TTFRZ)TTEMP=TTFRZ
                  FRC1=(TTEMP-TU(NK1))/(TTEMP-TBFRZ)
                ELSE
                  FRC1=1.
                  IFLAG=1
                ENDIF
                TTEMP=TU(NK1)
                !print*,'OLD FRC1',FRC1

!TWG Mar 2017
!Refine FRC1 to match super cooled water path estimate from Hu et al. [2010]
           IF (aercu_opt.gt.0) THEN
            IF(TU(NK1).GT.TBFRZ)THEN
             TC_HU10 = TU(NK1)-273.15
             SF_HU10 = -1.0*(P1_HU10+(P2_HU10*TC_HU10)+(P3_HU10*(TC_HU10**2))+ &
             (P4_HU10*(TC_HU10**3))+(P5_HU10*(TC_HU10**4))+(P6_HU10*(TC_HU10**5)))
             FRC1 = 1.0 - (1.0/(1.0 + EXP(SF_HU10)))
            ELSE
             FRC1=1.
             IFLAG=1
            ENDIF
           END IF
!END TWG 
!
!  DETERMINE THE EFFECTS OF LIQUID WATER FREEZING WHEN TEMPERATURE
!...IS BELOW TTFRZ...
!
                QFRZ = (QLIQ(NK1)+QNEWLQ)*FRC1
                QNEWIC=QNEWIC+QNEWLQ*FRC1
                QNEWLQ=QNEWLQ-QNEWLQ*FRC1
                QICE(NK1) = QICE(NK1)+QLIQ(NK1)*FRC1
                QLIQ(NK1) = QLIQ(NK1)-QLIQ(NK1)*FRC1
                CALL DTFRZNEW(TU(NK1),P0(NK1),THETEU(NK1),QU(NK1),QFRZ,         &
                          QICE(NK1),ALIQ,BLIQ,CLIQ,DLIQ)
              ENDIF
              TVU(NK1)=TU(NK1)*(1.+0.608*QU(NK1))

   IF (aercu_opt.gt.0) THEN
        QSATu(NK) = QSu/(1.+QSu)   ! saturated specific hum
         Aqnewlq(NK) = qnewlq
         Aqnewic(NK) = qnewic
   END IF

!
!  CALCULATE UPDRAFT VERTICAL VELOCITY AND PRECIPITATION FALLOUT...
!
              IF(NK.EQ.K)THEN
                BE=(TVLCL+TVU(NK1))/(TVEN+TV0(NK1))-1.
                BOTERM=2.*(Z0(NK1)-ZLCL)*G*BE/1.5
                DZZ=Z0(NK1)-ZLCL
              ELSE
                BE=(TVU(NK)+TVU(NK1))/(TV0(NK)+TV0(NK1))-1.
                BOTERM=2.*DZA(NK)*G*BE/1.5
                DZZ=DZA(NK)
              ENDIF
              ENTERM=2.*REI*WTW/UPOLD
!
! ckay
! using corrected RATE_kay for Test simulation #2... CGM July 2015
!
              IF(DX.GE.24.999E3) then
                RATE_kay = RATE
              else
                RATE_kay = RATE / Scale_Fac 
               end if
              CALL CONDLOAD(QLIQ(NK1),QICE(NK1),WTW,DZZ,BOTERM,ENTERM,      &
                   RATE_kay,QNEWLQ,QNEWIC,QLQOUT(NK1),QICOUT(NK1),G)

!
!...IF VERT VELOCITY IS LESS THAN ZERO, EXIT THE UPDRAFT LOOP AND,
!...IF CLOUD IS TALL ENOUGH, FINALIZE UPDRAFT CALCULATIONS...
!
              IF(WTW.LT.1.E-3)THEN
                EXIT
              ELSE
                WU(NK1)=SQRT(WTW)
              ENDIF
!...Calculate value of THETA-E in environment to entrain into updraft...
!
              CALL ENVIRTHT(P0(NK1),T0(NK1),Q0(NK1),THETEE(NK1),ALIQ,BLIQ,CLIQ,DLIQ)
!
!...REI IS THE RATE OF ENVIRONMENTAL INFLOW...
! New formulation for entrainment
!ckay introduce DX dependcy for the TOKIOKA Parameter =0.03
!ckay Kim et al 2011; Kang et al 2009; Lin et al 2013;  GCM findings

              TOKIOKA = 0.03
              TOKIOKA = TOKIOKA * Scale_Fac
              REI=VMFLCL*DP(NK1)*TOKIOKA/RAD         
!ckay
              TVQU(NK1)=TU(NK1)*(1.+0.608*QU(NK1)-QLIQ(NK1)-QICE(NK1))
              IF(NK.EQ.K)THEN
                DILBE=((TVLCL+TVQU(NK1))/(TVEN+TV0(NK1))-1.)*DZZ
              ELSE
                DILBE=((TVQU(NK)+TVQU(NK1))/(TV0(NK)+TV0(NK1))-1.)*DZZ
              ENDIF
              IF(DILBE.GT.0.)ABE=ABE+DILBE*G
!
!...IF CLOUD PARCELS ARE VIRTUALLY COLDER THAN THE ENVIRONMENT, MINIMAL 
!...ENTRAINMENT (0.5*REI) IS IMPOSED...
!
              IF(TVQU(NK1).LE.TV0(NK1))THEN    ! Entrain/Detrain IF BLOCK
                EE2=0.5       ! Kain (2004)  Eq. 4
                UD2=1.
                EQFRC(NK1)=0.
              ELSE
                LET=NK1
                TTMP=TVQU(NK1)
!
!...DETERMINE THE CRITICAL MIXED FRACTION OF UPDRAFT AND ENVIRONMENTAL AIR...
!
                F1=0.95
                F2=1.-F1
                THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)
                QTMP=F1*Q0(NK1)+F2*QU(NK1)
                TMPLIQ=F2*QLIQ(NK1)
                TMPICE=F2*QICE(NK1)
                call tpmix2(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,        &
                           qnewlq,qnewic,XLV1,XLV0,QSu)
                TU95=TTMP*(1.+0.608*QTMP-TMPLIQ-TMPICE)
                IF(TU95.GT.TV0(NK1))THEN
                  EE2=1.
                  UD2=0.
                  EQFRC(NK1)=1.0
                ELSE
                  F1=0.10
                  F2=1.-F1
                  THTTMP=F1*THETEE(NK1)+F2*THETEU(NK1)
                  QTMP=F1*Q0(NK1)+F2*QU(NK1)
                  TMPLIQ=F2*QLIQ(NK1)
                  TMPICE=F2*QICE(NK1)
                  call tpmix2(p0(nk1),thttmp,ttmp,qtmp,tmpliq,tmpice,        &
                               qnewlq,qnewic,XLV1,XLV0,QSu)
                  TU10=TTMP*(1.+0.608*QTMP-TMPLIQ-TMPICE)
                  TVDIFF = ABS(TU10-TVQU(NK1))
                  IF(TVDIFF.LT.1.e-3)THEN
                    EE2=1.
                    UD2=0.
                    EQFRC(NK1)=1.0
                  ELSE
                    EQFRC(NK1)=(TV0(NK1)-TVQU(NK1))*F1/(TU10-TVQU(NK1))
                    EQFRC(NK1)=AMAX1(0.0,EQFRC(NK1))
                    EQFRC(NK1)=AMIN1(1.0,EQFRC(NK1))
                    IF(EQFRC(NK1).EQ.1)THEN
                      EE2=1.
                      UD2=0.
                    ELSEIF(EQFRC(NK1).EQ.0.)THEN
                      EE2=0.
                      UD2=1.
                    ELSE
!
!...SUBROUTINE PROF5 INTEGRATES OVER THE GAUSSIAN DIST TO DETERMINE THE
!   FRACTIONAL ENTRAINMENT AND DETRAINMENT RATES...
!
                      CALL PROF5(EQFRC(NK1),EE2,UD2)
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF                            ! End of Entrain/Detrain IF BLOCK
!
!
!...NET ENTRAINMENT AND DETRAINMENT RATES ARE GIVEN BY THE AVERAGE FRACTIONAL
!   VALUES IN THE LAYER...
!
              EE2 = AMAX1(EE2,0.5)
              UD2 = 1.5*UD2
              UER(NK1)=0.5*REI*(EE1+EE2)
              UDR(NK1)=0.5*REI*(UD1+UD2)
!
!...IF THE CALCULATED UPDRAFT DETRAINMENT RATE IS GREATER THAN THE TOTAL
!   UPDRAFT MASS FLUX, ALL CLOUD MASS DETRAINS, EXIT UPDRAFT CALCULATIONS...
!
              IF(UMF(NK)-UDR(NK1).LT.10.)THEN
!
!...IF THE CALCULATED DETRAINED MASS FLUX IS GREATER THAN THE TOTAL UPD MASS
!   FLUX, IMPOSE TOTAL DETRAINMENT OF UPDRAFT MASS AT THE PREVIOUS MODEL LVL..
!   First, correct ABE calculation if needed...
!
                IF(DILBE.GT.0.)THEN
                  ABE=ABE-DILBE*G
                ENDIF
                LET=NK
!               WRITE(98,1015)P0(NK1)/100.
                EXIT 
              ELSE
                EE1=EE2
                UD1=UD2
                UPOLD=UMF(NK)-UDR(NK1)
                UPNEW=UPOLD+UER(NK1)
                UMF(NK1)=UPNEW
                DILFRC(NK1) = UPNEW/UPOLD
!
!...DETLQ AND DETIC ARE THE RATES OF DETRAINMENT OF LIQUID AND
!...ICE IN THE DETRAINING UPDRAFT MASS...
!
                DETLQ(NK1)=QLIQ(NK1)*UDR(NK1)
                DETIC(NK1)=QICE(NK1)*UDR(NK1)
                QDT(NK1)=QU(NK1)
                QU(NK1)=(UPOLD*QU(NK1)+UER(NK1)*Q0(NK1))/UPNEW
                THETEU(NK1)=(THETEU(NK1)*UPOLD+THETEE(NK1)*UER(NK1))/UPNEW
                QLIQ(NK1)=QLIQ(NK1)*UPOLD/UPNEW
                QICE(NK1)=QICE(NK1)*UPOLD/UPNEW
!
!...PPTLIQ IS THE RATE OF GENERATION (FALLOUT) OF
!...LIQUID PRECIP AT A GIVEN MODEL LVL, PPTICE THE SAME FOR ICE,
!...TRPPT IS THE TOTAL RATE OF PRODUCTION OF PRECIP UP TO THE
!...CURRENT MODEL LEVEL...
!
                PPTLIQ(NK1)=QLQOUT(NK1)*UMF(NK)
                PPTICE(NK1)=QICOUT(NK1)*UMF(NK)
!
                TRPPT=TRPPT+PPTLIQ(NK1)+PPTICE(NK1)
                IF(NK1.LE.KPBL)UER(NK1)=UER(NK1)+VMFLCL*DP(NK1)/DPTHMX
              ENDIF

!dkay
              IF (aercu_opt.gt.0) THEN
                  eps1u = 0.622
                  alatent = 2.54E6
                  KQ = NK
                  JK = KTE-KQ+1
                   muu(1,JK) = UMF(KQ)/VMFLCL     ! normalized updraft mass flux
                   duu(1,JK) = UDR(KQ)/DZQ(KQ)/VMFLCL  ! fractional detrainment rate in units of per meter
                   EUU(1,JK) = UER(KQ)/DZQ(KQ)/VMFLCL  ! normalized entrainment rate in unts of per meter
                   cmel(1,JK) = muu(1,JK)*AQNEWLQ(KQ)/DZQ(KQ)
                   cmei(1,JK) = muu(1,JK)*AQNEWIC(KQ)/DZQ(KQ)
                  gamhat(1,JK) = QSATu(KQ)*(1.+QSATu(KQ)/eps1u)   &
                 *eps1u*alatent/(R*oldTU(KQ)**2)*alatent/CP
                   wu_mskf_act(JK) = WU(KQ)  ! kf updraft velocity incloud
                   qc_mskf_act(JK) = AQNEWLQ(KQ)
                   qi_mskf_act(JK)=AQNEWIC(KQ)
              END IF
!end dkay


!
            END DO updraft

!dkay
          IF (aercu_opt.gt.0) THEN
            Zfu(1,KTE+1) = 0.0

                CPin   = CP

            EPSI0(1) = 2.0E-4
            DO KQ=KTS,KTE
              JK = KTE-KQ+1
                   zfu(1,JK) = zf_wrf(KQ)
                   su(1,JK) = oldTU(KQ)*(1.0+0.622*QSATu(KQ)) + (G*zf_wrf(KQ))/CP !TWG updraft temperature calulation
                   quu(1,JK) = oldQU(KQ)/(1.+oldQU(KQ)) ! specific humidity of updraft
                   pru(1,JK) = P0(KQ)/100.0    ! in millibars
                   TEE(1,JK) = T0(KQ)    ! ccc
                   QEE(1,JK) = Q0(KQ)/(1.+Q0(KQ))           ! specific humidity of environment

                   qee(1,JK) = oldQU(KQ)/(1.+oldQU(KQ)) ! specific humidity of updraft
                   QSATZM(1,JK) = QSATu(KQ)
!
!psh: Now, using aerosol concs from CESM
!
                    denSplume = P0(KQ)/(R*oldTU(KQ))
                    psh_fac = 1.0E-09/denSplume        ! convert ug/m3 to kg/kg

                    aer_mmr(1,JK, 1) = aercu_fct*aerocu2d(KQ,6)*psh_fac
                    aer_mmr(1,JK, 2) = aercu_fct*aerocu2d(KQ,5)*psh_fac
                    aer_mmr(1,JK, 3) = aercu_fct*1.44*aerocu2d(KQ,1)*psh_fac
                    aer_mmr(1,JK, 4) = aercu_fct*1.44*aerocu2d(KQ,2)*psh_fac
                    aer_mmr(1,JK, 5) = aercu_fct*1.44*aerocu2d(KQ,3)*psh_fac
                    aer_mmr(1,JK, 6) = aercu_fct*1.44*aerocu2d(KQ,4)*psh_fac
                    aer_mmr(1,JK, 7) = aercu_fct*1.54*aerocu2d(KQ,9)*psh_fac
                    aer_mmr(1,JK, 8) = aercu_fct*1.37*aerocu2d(KQ,7)*psh_fac
                    aer_mmr(1,JK, 9) = aercu_fct*1.25*aerocu2d(KQ,10)*psh_fac
                    aer_mmr(1,JK,10) = aercu_fct*1.37*aerocu2d(KQ,8)*psh_fac

!psh
                  gamhat(1,JK) = QSATu(KQ)*(1.+QSATu(KQ)/eps1u)   &
                 *eps1u*alatent/(R*oldTU(KQ)**2)*alatent/CP

              END DO
                  JTT(1) = KX-NK+1
                  JBB(1) = KX-K+1  ! updraft base level   =====>>> flipped for CAM5 indexing
                  if(jtt(1).gt.jbb(1))  then
                   JTT(1) = JBB(1)
                  end if
                  JLCL(1) = JBB(1) - 1
                  msg1 = 0
                  il2g = 1
                  grav  = G
                  Rdry = R
                  DTzmp = DT
!                 print *,'jtt,jbb=', JTT(1), JBB(1)

!dkay: call the new DCCMP scheme here
              NLEVZM = KTE-KTS+1      !  this is equal to pver in zm_mp
              NLEVZMP1 = NLEVZM + 1   !  pverp

                if(jtt(1).eq.1) then
                   print *,' cloud bottom is on ground!'
!                  print*,'I ',I,' J ',J
                   FATAL_ERROR('MSKF Cloud Bottom IS ON THE GROUND, diags')
                 end if
                if(jbb(1).eq.KTE) then
                   print *,' cloud top went through the roof!'
                   print *,'JTT, jbb, jlcl=',JTT(1),JBB(1),JLCL(1)
                   FATAL_ERROR('MSKF CLOUD TOP WENT OVER MODEL TOP, diags')
                 end if
            if(DCCMP) then

!           do kq=KTE,1,-1
!             print *,'wrf dz=',dzq(kq),(KTE-KQ+1)
!           end do

          call mskf_mphy(su,quu,muu,duu,cmel,cmei,zfu, pru,tee,qee,epsi0, &
               jbb,jtt,jlcl, msg1,il2g, grav, cpin, rdry,zmqliq,zmqice,zmqrain,zmqsnow,&
               rprd,wump, euu, ncmp,nimp,nrmp,nsmp,zmccn,sprd, frz, aer_mmr, dtzmp, &
               NLEVZM,NLEVZMP1,gamhat,qsatzm,wu_mskf_act,qc_mskf_act,qi_mskf_act,effc,effi,effs)
              end if

               Itest = 0
                if(Itest.eq.1) then

               write(121,*) 'k,nk, kq,jk,su,quuE3,muu,duu,cmel,zfu,pru,tee,&
              &qeeE3,zmqliqE4,zmqiceE4,rprd,wump,euu,ncmp,nimp,sprd,frz'
               do kq=K,NK
               JK = KTE-KQ+1
               write (121,2021) k,nk,kq,jk,su(1,jk),quu(1,jk)*1000,muu(1,jk),duu(1,jk),cmel(1,jk)
               write (121,2022) zfu(1,jk),pru(1,jk),tee(1,jk),qee(1,jk)*1000,zmqliq(1,jk)*1.e3,zmqice(1,jk)*1.e3
               write (121,2022) rprd(1,jk),wump(1,jk),euu(1,jk),ncmp(1,jk),nimp(1,jk),sprd(1,jk)
               write (121,2023) frz(1,jk)

2021          format(4I3,6(1x,E13.6))
2022          format(6(1x,e13.6))
2023          format(2(1x,e13.6))
               end do
                end if ! itest
            

               if(DCCMP) then
              do kq=KTS,KTE
               QLIQ(KQ) = 0.0
               QICE(KQ) = 0.0
               QRAIN(KQ) = 0.0
               QSNOW(KQ) = 0.0
               NLIQ(KQ) = 0.0
               NICE(KQ) = 0.0
               NRAIN(KQ) = 0.0
               NSNOW(KQ) = 0.0
               CCN(KQ) = 0.0
               EFFCH(KQ) = 2.51
               EFFIH(KQ) = 4.99
               EFFSH(KQ) = 9.99
               PPTLIQ(KQ)=0.0  ! nov23
               PPTICE(KQ)=0.0  ! nov23
               QLQOUT(KQ)=0.0  ! nov23
               QICOUT(KQ)=0.0  ! nov23
               DETLQ(KQ)=0.0  ! dec26
               DETIC(KQ)=0.0  ! dec26
              end do
              TRPPT = 0.0
              DO KQ=KTS, KTE

               JK = KX-KQ+1
!               print *,'kf qliq=', QLIQ(KQ)
                 QLIQ(KQ) = amax1(0.0,zmqliq(1,JK))
                 QICE(KQ) = amax1(0.0,zmqice(1,JK))
!TWG 06/14/16
                 QRAIN(KQ) = amax1(0.0,zmqrain(1,JK))
                 QSNOW(KQ) = amax1(0.0,zmqsnow(1,JK))
                 NLIQ(KQ) = amax1(0.0,ncmp(1,JK))
                 NICE(KQ) = amax1(0.0,nimp(1,JK))
                 NRAIN(KQ) = amax1(0.0,nrmp(1,JK))
                 NSNOW(KQ) = amax1(0.0,nsmp(1,JK))
                 CCN(KQ) = amax1(0.0,zmccn(1,JK))
                 EFFCH(KQ) = MAX(2.49, MIN(effc(1,JK), 50.))
                 EFFIH(KQ) = MAX(4.99, MIN(effi(1,JK), 125.))
                 EFFSH(KQ) = MAX(9.99, MIN(effs(1,JK), 999.))
! END TWG        
                  DETLQ(KQ)= QLIQ(KQ)*UDR(KQ)
                  DETIC(KQ)= QICE(KQ)*UDR(KQ)
!               print *,'zm qliq=', QLIQ(KQ)
                   densPlume = PPTLIQ(KQ)
!nov23
                  if(rprd(1,JK).lt.0.0) rprd(1,JK) = 0.0
                  if(sprd(1,JK).lt.0.0) sprd(1,JK) = 0.0

                   QLQOUT(KQ)=rprd(1,JK)*dzq(KQ)
                  QICOUT(KQ)=sprd(1,JK)*dzq(KQ)
                  PPTLIQ(KQ)=QLQOUT(KQ)*VMFLCL   ! check this out
                  PPTICE(KQ)=QICOUT(KQ)*VMFLCL   !   ditto

                 TRPPT=TRPPT+PPTLIQ(KQ)+PPTICE(KQ)
!                  if(densPlume.gt.0.0) then
!               print *,'zm pptliq=', &
!        PPTLIQ(KQ),'kf pptliq=',oldPPTLIQ(kq),'KQ=',KQ
!               end if
!                 print *,'zmQliqout=',kq,densPlume,VMFLCL,pptliq(kq)
!                if((i.ge.60.and.i.le.65).and.(j.ge.60.and.j.le.65)) then
!                     print *, 'KF & MP qliq=', Aqnewlq(nk),QLIQ(NK)
!                     print *, 'mu & du=', muu(1,JK), duu(1,JK)
!                      print *,'i,j,k=',I,J,KQ
!                   end if

              END DO

              end if  ! dccmp
!dkay
     
2999         CONTINUE
     END IF


!
!...CHECK CLOUD DEPTH...IF CLOUD IS TALL ENOUGH, ESTIMATE THE EQUILIBRIUM
!   TEMPERATURE LEVEL (LET) AND ADJUST MASS FLUX PROFILE AT CLOUD TOP SO
!   THAT MASS FLUX DECREASES TO ZERO AS A LINEAR FUNCTION OF PRESSURE BETWEEN
!   THE LET AND CLOUD TOP...
!     
!...LTOP IS THE MODEL LEVEL JUST BELOW THE LEVEL AT WHICH VERTICAL VELOCITY
!   FIRST BECOMES NEGATIVE...
!     
            LTOP=NK
            CLDHGT(LC)=Z0(LTOP)-ZLCL 
!
!...Instead of using the same minimum cloud height (for deep convection)
!...everywhere, try specifying minimum cloud depth as a function of TLCL...
!
!   Kain (2004)  Eq. 7
!
            IF(TLCL.GT.293.)THEN
              CHMIN = 4.E3
            ELSEIF(TLCL.LE.293. .and. TLCL.GE.273)THEN
              CHMIN = 2.E3 + 100.*(TLCL-273.)
            ELSEIF(TLCL.LT.273.)THEN
              CHMIN = 2.E3
            ENDIF
!ckay
!ckay: if mean env RH with respect to water/ice is over 99% then dont allow KF
!ckay: added saturation w.r.to ice june 10, 2015
! to avoid double counting
          envRHavg = 0.0
          DO NK=K-1,LTOP+1
           if(T0(NK).LE.273.16) then
             envEsat = 6.112*exp(21.87*(T0(NK)-273.16)/(T0(NK)-7.66))
           else
             envEsat = 6.112*exp(17.67*(T0(NK)-273.16)/(243.5+T0(NK)-273.16))
           end if
           envEsat = envEsat * 100.0  ! to hPa
           envQsat = 0.622*envEsat/(P0(NK)-envEsat)
           envRH  = Q0(NK)/envQsat
            if(NK.GT.K.and.envRH.LT.0.99)  then
              envRHavg = 0.0
              goto 2020
            end if
            envRHavg = envRHavg + envRH
          END DO
!ckay ; get vertically averaged envRHavg
         envRHavg = envRHavg / float(LTOP-K+1+2)
2020     continue

!     
!...If cloud top height is less than the specified minimum for deep 
!...convection, save value to consider this level as source for 
!...shallow convection, go back up to check next level...
!     
!...Try specifying minimum cloud depth as a function of TLCL...
!
!
!...DO NOT ALLOW ANY CLOUD FROM THIS LAYER IF:
!
!...            1.) if there is no CAPE, or 
!...            2.) cloud top is at model level just above LCL, or
!...            3.) cloud top is within updraft source layer, or
!...            4.) cloud-top detrainment layer begins within 
!...                updraft source layer.
!...ckay        5.) if the environment is supersaturated i.e., RH > 100%
!...ckay            For now, with respect to water
!

            IF(LTOP.LE.KLCL .or. LTOP.LE.KPBL .or. LET+1.LE.KPBL &
              .or. envRHavg.ge.1.01)THEN  ! No Convection Allowed
!ckay

              CLDHGT(LC)=0.
              DO NK=K,LTOP
                UMF(NK)=0.
                UDR(NK)=0.
                UER(NK)=0.
                DETLQ(NK)=0.
                DETIC(NK)=0.
                PPTLIQ(NK)=0.
                PPTICE(NK)=0.
!ckay
                dp_cf1d(NK)=0.
                sh_cf1d(NK)=0.
                w_up1d(NK)=0.
              ENDDO
!        
            ELSEIF(CLDHGT(LC).GT.CHMIN .and. ABE.GT.1)THEN      ! Deep Convection allowed
              ISHALL=0
!ckay
              DO NK=K,LTOP
                sh_cf1d(NK)=0.
              ENDDO
              EXIT usl
            ELSE
!
!...TO DISALLOW SHALLOW CONVECTION, COMMENT OUT NEXT LINE !!!!!!!!
              ISHALL = 1
!ckay
              DO NK=K,LTOP
                dp_cf1d(NK)=0.
                w_up1d(NK)=0.
              ENDDO
              IF(NU.EQ.NUCHM)THEN
                EXIT usl               ! Shallow Convection from this layer
              ELSE
! Remember this layer (by virtue of non-zero CLDHGT) as potential shallow-cloud layer
                DO NK=K,LTOP
                  UMF(NK)=0.
                  UDR(NK)=0.
                  UER(NK)=0.
                  DETLQ(NK)=0.
                  DETIC(NK)=0.
                  PPTLIQ(NK)=0.
                  PPTICE(NK)=0.
!ckay
                  dp_cf1d(NK)=0.
                  sh_cf1d(NK)=0.
                  w_up1d(NK)=0.
                ENDDO
              ENDIF
            ENDIF
          ENDIF  ! for trigger
        END DO usl
    IF(ISHALL.EQ.1)THEN
      KSTART=MAX0(KPBL,KLCL)
      LET=KSTART
    endif
!     
!...IF THE LET AND LTOP ARE THE SAME, DETRAIN ALL OF THE UPDRAFT MASS FL
!   THIS LEVEL...
!     
    IF(LET.EQ.LTOP)THEN
      UDR(LTOP)=UMF(LTOP)+UDR(LTOP)-UER(LTOP)
      DETLQ(LTOP)=QLIQ(LTOP)*UDR(LTOP)*UPNEW/UPOLD
      DETIC(LTOP)=QICE(LTOP)*UDR(LTOP)*UPNEW/UPOLD
      UER(LTOP)=0.
      UMF(LTOP)=0.
    ELSE 
!     
!   BEGIN TOTAL DETRAINMENT AT THE LEVEL ABOVE THE LET...
!     
      DPTT=0.
      DO NJ=LET+1,LTOP
        DPTT=DPTT+DP(NJ)
      ENDDO
      DUMFDP=UMF(LET)/DPTT
!     
!...ADJUST MASS FLUX PROFILES, DETRAINMENT RATES, AND PRECIPITATION FALL
!   RATES TO REFLECT THE LINEAR DECREASE IN MASS FLX BETWEEN THE LET AND
!     
      DO NK=LET+1,LTOP
!
!...entrainment is allowed at every level except for LTOP, so disallow
!...entrainment at LTOP and adjust entrainment rates between LET and LTOP
!...so the the dilution factor due to entrainment is not changed but 
!...the actual entrainment rate will change due due forced total 
!...detrainment in this layer...
!
        IF(NK.EQ.LTOP)THEN
          UDR(NK) = UMF(NK-1)
          UER(NK) = 0.
          DETLQ(NK) = UDR(NK)*QLIQ(NK)*DILFRC(NK)
          DETIC(NK) = UDR(NK)*QICE(NK)*DILFRC(NK)
        ELSE
          UMF(NK)=UMF(NK-1)-DP(NK)*DUMFDP
          UER(NK)=UMF(NK)*(1.-1./DILFRC(NK))
          UDR(NK)=UMF(NK-1)-UMF(NK)+UER(NK)
          DETLQ(NK)=UDR(NK)*QLIQ(NK)*DILFRC(NK)
          DETIC(NK)=UDR(NK)*QICE(NK)*DILFRC(NK)
        ENDIF
        IF(NK.GE.LET+2)THEN
          TRPPT=TRPPT-PPTLIQ(NK)-PPTICE(NK)
          PPTLIQ(NK)=UMF(NK-1)*QLQOUT(NK)
          PPTICE(NK)=UMF(NK-1)*QICOUT(NK)
          TRPPT=TRPPT+PPTLIQ(NK)+PPTICE(NK)
        ENDIF
      ENDDO
    ENDIF
!     
! Initialize some arrays below cloud base and above cloud top...
!
    DO NK=1,LTOP
      IF(T0(NK).GT.T00)ML=NK
    ENDDO
    DO NK=1,K
      IF(NK.GE.LC)THEN
        IF(NK.EQ.LC)THEN
          UMF(NK)=VMFLCL*DP(NK)/DPTHMX
          UER(NK)=VMFLCL*DP(NK)/DPTHMX
        ELSEIF(NK.LE.KPBL)THEN
          UER(NK)=VMFLCL*DP(NK)/DPTHMX
          UMF(NK)=UMF(NK-1)+UER(NK)
        ELSE
          UMF(NK)=VMFLCL
          UER(NK)=0.
        ENDIF
        TU(NK)=TMIX+(Z0(NK)-ZMIX)*GDRY
        QU(NK)=QMIX
        WU(NK)=WLCL
      ELSE
        TU(NK)=0.
        QU(NK)=0.
        UMF(NK)=0.
        WU(NK)=0.
        UER(NK)=0.
!ckay
        dp_cf1d(NK)=0.
        sh_cf1d(NK)=0.
        W_up1d(NK)=0.
      ENDIF
      UDR(NK)=0.
      QDT(NK)=0.
      QLIQ(NK)=0.
      QICE(NK)=0.
      QLQOUT(NK)=0.
      QICOUT(NK)=0.
      PPTLIQ(NK)=0.
      PPTICE(NK)=0.
      DETLQ(NK)=0.
      DETIC(NK)=0.
      RATIO2(NK)=0.
      CALL ENVIRTHT(P0(NK),T0(NK),Q0(NK),THETEE(NK),ALIQ,BLIQ,CLIQ,DLIQ)
      EQFRC(NK)=1.0
    ENDDO
!     
      LTOP1=LTOP+1
      LTOPM1=LTOP-1
!     
!...DEFINE VARIABLES ABOVE CLOUD TOP...
!     
      DO NK=LTOP1,KX
        UMF(NK)=0.
        UDR(NK)=0.
        UER(NK)=0.
        QDT(NK)=0.
        QLIQ(NK)=0.
        QICE(NK)=0.
        QLQOUT(NK)=0.
        QICOUT(NK)=0.
        DETLQ(NK)=0.
        DETIC(NK)=0.
        PPTLIQ(NK)=0.
        PPTICE(NK)=0.
        IF(NK.GT.LTOP1)THEN
          TU(NK)=0.
          QU(NK)=0.
          WU(NK)=0.
!ckay
          dp_cf1d(NK)=0.
          sh_cf1d(NK)=0.
          w_up1d(NK)=0.
        ENDIF
        THTA0(NK)=0.
        THTAU(NK)=0.
        EMS(NK)=0.
        EMSD(NK)=0.
        TG(NK)=T0(NK)
        QG(NK)=Q0(NK)
        QLG(NK)=0.
        QIG(NK)=0.
        QRG(NK)=0.
        QSG(NK)=0.
        OMG(NK)=0.
      ENDDO
        OMG(KX+1)=0.
        DO NK=1,LTOP
          EMS(NK)=DP(NK)*DXSQ/G
          EMSD(NK)=1./EMS(NK)
!     
!...INITIALIZE SOME VARIABLES TO BE USED LATER IN THE VERT ADVECTION SCHEME
!     
          EXN(NK)=(P00/P0(NK))**(0.2854*(1.-0.28*QDT(NK)))
          THTAU(NK)=TU(NK)*EXN(NK)
          EXN(NK)=(P00/P0(NK))**(0.2854*(1.-0.28*Q0(NK)))
          THTA0(NK)=T0(NK)*EXN(NK)
          DDILFRC(NK) = 1./DILFRC(NK)
          OMG(NK)=0.
        ENDDO
!     IF (XTIME.LT.10.)THEN
!      WRITE(98,1025)KLCL,ZLCL,DTLCL,LTOP,P0(LTOP),IFLAG,
!    * TMIX-T00,PMIX,QMIX,ABE
!      WRITE(98,1030)P0(LET)/100.,P0(LTOP)/100.,VMFLCL,PLCL/100.,
!    * WLCL,CLDHGT
!     ENDIF
!     
!...COMPUTE CONVECTIVE TIME SCALE(TIMEC). THE MEAN WIND AT THE LCL
!...AND MIDTROPOSPHERE IS USED.
!     
        WSPD(KLCL)=SQRT(U0(KLCL)*U0(KLCL)+V0(KLCL)*V0(KLCL))
        WSPD(L5)=SQRT(U0(L5)*U0(L5)+V0(L5)*V0(L5))
        WSPD(LTOP)=SQRT(U0(LTOP)*U0(LTOP)+V0(LTOP)*V0(LTOP))
        VCONV=.5*(WSPD(KLCL)+WSPD(L5))
!...for ETA model, DX is a function of location...
        TIMEC=DX/VCONV
        TADVEC=TIMEC
!
!ckay
!new dynTau based on subcloud layer scales : note Z0(KPBL)=altitude of source layer

         TIMEC = Amax1(CHMIN,CLDHGT(LC))
         TIMEC = TIMEC*Scale_Fac

!ckay SCLvel = SubCloudLayerVELOCITY = Wsb
!ckay estimate WSTAR to allow most of the PBL schemes here...
         densPlume = P0(1)/R/T0(1)
         WST = HFX1/densPlume/CP  ! hfx in kinematic units  CKAY
         WST = amax1(0.,WST)  ! +ve hfx only
         WST = G*PBLH1*WST
         thetav = (1.E5/P0(1))**(R/CP)
         thetav = T0(1) * thetav
         eps1u = 0.622
         thetav = thetav*(1.+Q0(1)*eps1u)
         WST = WST/thetav
         WST = WST**0.3333

         SCLvel = WST**3
         ZLCL_KAY = amax1(ZLCL,Z0(KPBL))
         SCLvel = SCLvel/PBLH1
         SCLvel = SCLvel*ZLCL_kay   
         SCLvel = SCLvel**0.333   ! Wsb=SubCloudLayerVelocity for ConvectivePBL
         if(ZOL1.le.0.0) then
           FRC2=3.8*UST1*UST1
           FRC2 = FRC2 + 0.22*SCLvel*SCLvel
           zz_kay = -1.0*ZOL1
           ZLCL_KAY = zz_kay**(2./3.)
           ZLCL_KAY = ZLCL_KAY * (1.9*UST1*UST1)
           FRC2 = FRC2 + ZLCL_KAY
         else
           FRC2=3.8*UST1*UST1
         end if 

          FRC2 = SQRT(FRC2)  
          SCLvel =  FRC2  ! Wsb=new subcloud layer velocity scale for all conditions
         IF(SCLvel.lt.0.1) SCLvel = 0.1
         if(ABE.le.0.0) ABE = 1.0 
         TIMEC = TIMEC/((0.03*SCLvel*ABE)**0.3333)

!ckay: this dynTau is good for the Deep as well as Shallow Cu clouds
        TIMEC = AMAX1(TADVEC, TIMEC)
 
         NIC=NINT(TIMEC/DT)
         TIMEC=FLOAT(NIC)*DT
         TIMEC=MIN(TIMEC,86400.)  !JRJ Ramboll: cap convective time scale at 24 hrs
!     
!...COMPUTE WIND SHEAR AND PRECIPITATION EFFICIENCY.
!     
        IF(WSPD(LTOP).GT.WSPD(KLCL))THEN
          SHSIGN=1.
        ELSE
          SHSIGN=-1.
        ENDIF
        VWS=(U0(LTOP)-U0(KLCL))*(U0(LTOP)-U0(KLCL))+(V0(LTOP)-V0(KLCL))*   &
            (V0(LTOP)-V0(KLCL))
        VWS=1.E3*SHSIGN*SQRT(VWS)/(Z0(LTOP)-Z0(LCL))
        PEF=1.591+VWS*(-.639+VWS*(9.53E-2-VWS*4.96E-3))
        PEF=AMAX1(PEF,.2)
        PEF=AMIN1(PEF,.9)
!     
!...PRECIPITATION EFFICIENCY IS A FUNCTION OF THE HEIGHT OF CLOUD BASE.
!     
        CBH=(ZLCL-Z0(1))*3.281E-3
        IF(CBH.LT.3.)THEN
          RCBH=.02
        ELSE
          RCBH=.96729352+CBH*(-.70034167+CBH*(.162179896+CBH*(-            &
               1.2569798E-2+CBH*(4.2772E-4-CBH*5.44E-6))))
        ENDIF
        IF(CBH.GT.25)RCBH=2.4
        PEFCBH=1./(1.+RCBH)
        PEFCBH=AMIN1(PEFCBH,.9)
!     
!... MEAN PEF. IS USED TO COMPUTE RAINFALL.
!     
        PEFF=.5*(PEF+PEFCBH)
        PEFF2 = PEFF                                ! JSK MODS
       IF(IPRNT)THEN  
!         WRITE(98,1035)PEF,PEFCBH,LC,LET,WKL,VWS
         WRITE(message,1035)PEF,PEFCBH,LC,LET,WKL,VWS
         CALL wrf_message( message )
!       flush(98)   
       endif     
!        WRITE(98,1035)PEF,PEFCBH,LC,LET,WKL,VWS
!*****************************************************************
!                                                                *
!                  COMPUTE DOWNDRAFT PROPERTIES                  *
!                                                                *
!*****************************************************************
!     
!     
       TDER=0.
 devap:IF(ISHALL.EQ.1)THEN
         LFS = 1
       ELSE
!
!...start downdraft about 150 mb above cloud base...
!
!        KSTART=MAX0(KPBL,KLCL)
!        KSTART=KPBL                                  ! Changed 7/23/99
         KSTART=KPBL+1                                ! Changed 7/23/99
         KLFS = LET-1
         DO NK = KSTART+1,KL
           DPPP = P0(KSTART)-P0(NK)
!          IF(DPPP.GT.200.E2)THEN
           IF(DPPP.GT.150.E2)THEN
             KLFS = NK
             EXIT 
           ENDIF
         ENDDO
         KLFS = MIN0(KLFS,LET-1)
         LFS = KLFS
!
!...if LFS is not at least 50 mb above cloud base (implying that the 
!...level of equil temp, LET, is just above cloud base) do not allow a
!...downdraft...
!
        IF((P0(KSTART)-P0(LFS)).GT.50.E2)THEN
          THETED(LFS) = THETEE(LFS)
          QD(LFS) = Q0(LFS)
!
!...call tpmix2dd to find wet-bulb temp, qv...
!
          call tpmix2dd(p0(lfs),theted(lfs),tz(lfs),qss)
          THTAD(LFS)=TZ(LFS)*(P00/P0(LFS))**(0.2854*(1.-0.28*QSS))
!     
!...TAKE A FIRST GUESS AT THE INITIAL DOWNDRAFT MASS FLUX...
!     
          TVD(LFS)=TZ(LFS)*(1.+0.608*QSS)
          RDD=P0(LFS)/(R*TVD(LFS))
          A1=(1.-PEFF)*AU0
          DMF(LFS)=-A1*RDD
          DER(LFS)=DMF(LFS)
          DDR(LFS)=0.
          RHBAR = RH(LFS)*DP(LFS)
          DPTT = DP(LFS)
          DO ND = LFS-1,KSTART,-1
            ND1 = ND+1
            DER(ND)=DER(LFS)*EMS(ND)/EMS(LFS)
            DDR(ND)=0.
            DMF(ND)=DMF(ND1)+DER(ND)
            THETED(ND)=(THETED(ND1)*DMF(ND1)+THETEE(ND)*DER(ND))/DMF(ND)
            QD(ND)=(QD(ND1)*DMF(ND1)+Q0(ND)*DER(ND))/DMF(ND)    
            DPTT = DPTT+DP(ND)
            RHBAR = RHBAR+RH(ND)*DP(ND)
          ENDDO
          RHBAR = RHBAR/DPTT
          DMFFRC = 2.*(1.-RHBAR)    ! Kain (2004) eq. 11
          DPDD = 0.
!...Calculate melting effect
!... first, compute total frozen precipitation generated...
!
          pptmlt = 0.
          DO NK = KLCL,LTOP
            PPTMLT = PPTMLT+PPTICE(NK)
          ENDDO
          if(lc.lt.ml)then
!...For now, calculate melting effect as if DMF = -UMF at KLCL, i.e., as
!...if DMFFRC=1.  Otherwise, for small DMFFRC, DTMELT gets too large!
!...12/14/98 jsk...
            DTMELT = RLF*PPTMLT/(CP*UMF(KLCL))
          else
            DTMELT = 0.
          endif
          LDT = MIN0(LFS-1,KSTART-1)
!
          call tpmix2dd(p0(kstart),theted(kstart),tz(kstart),qss)
!
          tz(kstart) = tz(kstart)-dtmelt
          ES=ALIQ*EXP((BLIQ*TZ(KSTART)-CLIQ)/(TZ(KSTART)-DLIQ))
          QSS=0.622*ES/(P0(KSTART)-ES)
          THETED(KSTART)=TZ(KSTART)*(1.E5/P0(KSTART))**(0.2854*(1.-0.28*QSS))*    &
                EXP((3374.6525/TZ(KSTART)-2.5403)*QSS*(1.+0.81*QSS))
!....  
          LDT = MIN0(LFS-1,KSTART-1)
          DO ND = LDT,1,-1
            DPDD = DPDD+DP(ND)
            THETED(ND) = THETED(KSTART)
            QD(ND)     = QD(KSTART)       
!
!...call tpmix2dd to find wet bulb temp, saturation mixing ratio...
!
            call tpmix2dd(p0(nd),theted(nd),tz(nd),qss)
            qsd(nd) = qss
!
!...specify RH decrease of 20%/km in downdraft...
!
            RHH = 1.-0.2/1000.*(Z0(KSTART)-Z0(ND))
!
!...adjust downdraft TEMP, Q to specified RH:
!
            IF(RHH.LT.1.)THEN
              DSSDT=(CLIQ-BLIQ*DLIQ)/((TZ(ND)-DLIQ)*(TZ(ND)-DLIQ))
              RL=XLV0-XLV1*TZ(ND)
              DTMP=RL*QSS*(1.-RHH)/(CP+RL*RHH*QSS*DSSDT)
              T1RH=TZ(ND)+DTMP
              ES=RHH*ALIQ*EXP((BLIQ*T1RH-CLIQ)/(T1RH-DLIQ))
              QSRH=0.622*ES/(P0(ND)-ES)
!
!...CHECK TO SEE IF MIXING RATIO AT SPECIFIED RH IS LESS THAN ACTUAL
!...MIXING RATIO...IF SO, ADJUST TO GIVE ZERO EVAPORATION...
!
              IF(QSRH.LT.QD(ND))THEN
                QSRH=QD(ND)
                T1RH=TZ(ND)+(QSS-QSRH)*RL/CP
              ENDIF
              TZ(ND)=T1RH
              QSS=QSRH
              QSD(ND) = QSS
            ENDIF         
            TVD(nd) = tz(nd)*(1.+0.608*qsd(nd))
            IF(TVD(ND).GT.TV0(ND).OR.ND.EQ.1)THEN
              LDB=ND
              EXIT
            ENDIF
          ENDDO
          IF((P0(LDB)-P0(LFS)) .gt. 50.E2)THEN   ! minimum Downdraft depth! 
            DO ND=LDT,LDB,-1
              ND1 = ND+1
              DDR(ND) = -DMF(KSTART)*DP(ND)/DPDD
              DER(ND) = 0.
              DMF(ND) = DMF(ND1)+DDR(ND)
              TDER=TDER+(QSD(nd)-QD(ND))*DDR(ND)
              QD(ND)=QSD(nd)
              THTAD(ND)=TZ(ND)*(P00/P0(ND))**(0.2854*(1.-0.28*QD(ND)))
            ENDDO
          ENDIF
        ENDIF
      ENDIF devap

!...IF DOWNDRAFT DOES NOT EVAPORATE ANY WATER FOR SPECIFIED RELATIVE
!...HUMIDITY, NO DOWNDRAFT IS ALLOWED...
!
d_mf:   IF(TDER.LT.1.)THEN
!           WRITE(98,3004)I,J 
!3004       FORMAT(' ','No Downdraft!;  I=',I3,2X,'J=',I3,'ISHALL =',I2)
          PPTFLX=TRPPT
          CPR=TRPPT
          TDER=0.
          CNDTNF=0.
          UPDINC=1.
          LDB=LFS
          DO NDK=1,LTOP
            DMF(NDK)=0.
            DER(NDK)=0.
            DDR(NDK)=0.
            THTAD(NDK)=0.
            WD(NDK)=0.
            TZ(NDK)=0.
            QD(NDK)=0.
          ENDDO
          AINCM2=100.
        ELSE 
          DDINC = -DMFFRC*UMF(KLCL)/DMF(KSTART)
          UPDINC=1.
          IF(TDER*DDINC.GT.TRPPT)THEN
            DDINC = TRPPT/TDER
          ENDIF
          TDER = TDER*DDINC
          DO NK=LDB,LFS
            DMF(NK)=DMF(NK)*DDINC
            DER(NK)=DER(NK)*DDINC
            DDR(NK)=DDR(NK)*DDINC
          ENDDO
         CPR=TRPPT
         PPTFLX = TRPPT-TDER
         PEFF=PPTFLX/TRPPT
         IF(IPRNT)THEN
!           write(98,*)'PRECIP EFFICIENCY =',PEFF
           write(message,*)'PRECIP EFFICIENCY =',PEFF
           CALL wrf_message(message)
!          flush(98)   
         ENDIF
!
!
!...ADJUST UPDRAFT MASS FLUX, MASS DETRAINMENT RATE, AND LIQUID WATER AN
!   DETRAINMENT RATES TO BE CONSISTENT WITH THE TRANSFER OF THE ESTIMATE
!   FROM THE UPDRAFT TO THE DOWNDRAFT AT THE LFS...
!     
!         DO NK=LC,LFS
!           UMF(NK)=UMF(NK)*UPDINC
!           UDR(NK)=UDR(NK)*UPDINC
!           UER(NK)=UER(NK)*UPDINC
!           PPTLIQ(NK)=PPTLIQ(NK)*UPDINC
!           PPTICE(NK)=PPTICE(NK)*UPDINC
!           DETLQ(NK)=DETLQ(NK)*UPDINC
!           DETIC(NK)=DETIC(NK)*UPDINC
!         ENDDO
!
  
!...ZERO OUT THE ARRAYS FOR DOWNDRAFT DATA AT LEVELS ABOVE AND BELOW THE
!...DOWNDRAFT...
!     
         IF(LDB.GT.1)THEN
           DO NK=1,LDB-1
             DMF(NK)=0.
             DER(NK)=0.
             DDR(NK)=0.
             WD(NK)=0.
             TZ(NK)=0.
             QD(NK)=0.
             THTAD(NK)=0.
           ENDDO
         ENDIF
         DO NK=LFS+1,KX
           DMF(NK)=0.
           DER(NK)=0.
           DDR(NK)=0.
           WD(NK)=0.
           TZ(NK)=0.
           QD(NK)=0.
           THTAD(NK)=0.
         ENDDO
         DO NK=LDT+1,LFS-1
           TZ(NK)=0.
           QD(NK)=0.
           THTAD(NK)=0.
         ENDDO
       ENDIF d_mf
!
!...SET LIMITS ON THE UPDRAFT AND DOWNDRAFT MASS FLUXES SO THAT THE INFLOW
!   INTO CONVECTIVE DRAFTS FROM A GIVEN LAYER IS NO MORE THAN IS AVAILABLE
!   IN THAT LAYER INITIALLY...
!     
       AINCMX=1000.
       LMAX=MAX0(KLCL,LFS)
       DO NK=LC,LMAX
         IF((UER(NK)-DER(NK)).GT.1.e-3)THEN
           AINCM1=EMS(NK)/((UER(NK)-DER(NK))*TIMEC)
           AINCMX=AMIN1(AINCMX,AINCM1)
         ENDIF
       ENDDO
       AINC=1.
       IF(AINCMX.LT.AINC)AINC=AINCMX
!     
!...SAVE THE RELEVENT VARIABLES FOR A UNIT UPDRAFT AND DOWNDRAFT...THEY WILL 
!...BE ITERATIVELY ADJUSTED BY THE FACTOR AINC TO SATISFY THE STABILIZATION
!...CLOSURE...
!     
       TDER2=TDER
       PPTFL2=PPTFLX
       DO NK=1,LTOP
         DETLQ2(NK)=DETLQ(NK)
         DETIC2(NK)=DETIC(NK)
         UDR2(NK)=UDR(NK)
         UER2(NK)=UER(NK)
         DDR2(NK)=DDR(NK)
         DER2(NK)=DER(NK)
         UMF2(NK)=UMF(NK)
         DMF2(NK)=DMF(NK)
       ENDDO
       FABE=1.
       STAB=0.95
       NOITR=0
       ISTOP=0
!
        IF(ISHALL.EQ.1)THEN                              ! First for shallow convection
!
! No iteration for shallow convection; if turbulent kinetic energy (TKE) is available
! from a turbulence parameterization, scale cloud-base updraft mass flux as a function
! of TKE, but for now, just specify shallow-cloud mass flux using TKEMAX = 5...
!
!...find the maximum TKE value between LC and KLCL...
!         TKEMAX = 0.
          TKEMAX = 5.
!          DO 173 K = LC,KLCL
!            NK = KX-K+1
!            TKEMAX = AMAX1(TKEMAX,Q2(I,J,NK))
! 173      CONTINUE
!          TKEMAX = AMIN1(TKEMAX,10.)
!          TKEMAX = AMAX1(TKEMAX,5.)
!c         TKEMAX = 10.
!c...3_24_99...DPMIN was changed for shallow convection so that it is the
!c...          the same as for deep convection (5.E3).  Since this doubles
!c...          (roughly) the value of DPTHMX, add a factor of 0.5 to calcu-
!c...          lation of EVAC...
!c         EVAC  = TKEMAX*0.1
          EVAC  = 0.5*TKEMAX*0.1
!         AINC = 0.1*DPTHMX*DXIJ*DXIJ/(VMFLCL*G*TIMEC)
!          AINC = EVAC*DPTHMX*DX(I,J)*DX(I,J)/(VMFLCL*G*TIMEC)
          AINC = EVAC*DPTHMX*DXSQ/(VMFLCL*G*TIMEC)
          TDER=TDER2*AINC
          PPTFLX=PPTFL2*AINC
          DO NK=1,LTOP
            UMF(NK)=UMF2(NK)*AINC
            DMF(NK)=DMF2(NK)*AINC
            DETLQ(NK)=DETLQ2(NK)*AINC
            DETIC(NK)=DETIC2(NK)*AINC
            UDR(NK)=UDR2(NK)*AINC
            UER(NK)=UER2(NK)*AINC
            DER(NK)=DER2(NK)*AINC
            DDR(NK)=DDR2(NK)*AINC
          ENDDO
        ENDIF                                           ! Otherwise for deep convection
! use iterative procedure to find mass fluxes...
iter:     DO NCOUNT=1,10
!     
!*****************************************************************
!                                                                *
!           COMPUTE PROPERTIES FOR COMPENSATIONAL SUBSIDENCE     *
!                                                                *
!*****************************************************************
!     
!...DETERMINE OMEGA VALUE NECESSARY AT TOP AND BOTTOM OF EACH LAYER TO
!...SATISFY MASS CONTINUITY...
!     
            DTT=TIMEC
            DO NK=1,LTOP
              DOMGDP(NK)=-(UER(NK)-DER(NK)-UDR(NK)-DDR(NK))*EMSD(NK)
              IF(NK.GT.1)THEN
                OMG(NK)=OMG(NK-1)-DP(NK-1)*DOMGDP(NK-1)
                ABSOMG = ABS(OMG(NK))
                ABSOMGTC = ABSOMG*TIMEC
                FRDP = 0.75*DP(NK-1)
                IF(ABSOMGTC.GT.FRDP)THEN
                  DTT1 = FRDP/ABSOMG
                  DTT=AMIN1(DTT,DTT1)
                ENDIF
              ENDIF
            ENDDO
            DO NK=1,LTOP
              THPA(NK)=THTA0(NK)
              QPA(NK)=Q0(NK)
              NSTEP=NINT(TIMEC/DTT+1)
              DTIME=TIMEC/FLOAT(NSTEP)
              FXM(NK)=OMG(NK)*DXSQ/G
            ENDDO
!     
!...DO AN UPSTREAM/FORWARD-IN-TIME ADVECTION OF THETA, QV...
!     
        DO NTC=1,NSTEP
!     
!...ASSIGN THETA AND Q VALUES AT THE TOP AND BOTTOM OF EACH LAYER BASED ON
!...SIGN OF OMEGA...
!     
            DO  NK=1,LTOP
              THFXIN(NK)=0.
              THFXOUT(NK)=0.
              QFXIN(NK)=0.
              QFXOUT(NK)=0.
            ENDDO
            DO NK=2,LTOP
              IF(OMG(NK).LE.0.)THEN
                THFXIN(NK)=-FXM(NK)*THPA(NK-1)
                QFXIN(NK)=-FXM(NK)*QPA(NK-1)
                THFXOUT(NK-1)=THFXOUT(NK-1)+THFXIN(NK)
                QFXOUT(NK-1)=QFXOUT(NK-1)+QFXIN(NK)
              ELSE
                THFXOUT(NK)=FXM(NK)*THPA(NK)
                QFXOUT(NK)=FXM(NK)*QPA(NK)
                THFXIN(NK-1)=THFXIN(NK-1)+THFXOUT(NK)
                QFXIN(NK-1)=QFXIN(NK-1)+QFXOUT(NK)
              ENDIF
            ENDDO
!     
!...UPDATE THE THETA AND QV VALUES AT EACH LEVEL...
!     
            DO NK=1,LTOP
              THPA(NK)=THPA(NK)+(THFXIN(NK)+UDR(NK)*THTAU(NK)+DDR(NK)*      &
                       THTAD(NK)-THFXOUT(NK)-(UER(NK)-DER(NK))*THTA0(NK))*  &
                       DTIME*EMSD(NK)
              QPA(NK)=QPA(NK)+(QFXIN(NK)+UDR(NK)*QDT(NK)+DDR(NK)*QD(NK)-    &
                      QFXOUT(NK)-(UER(NK)-DER(NK))*Q0(NK))*DTIME*EMSD(NK)
            ENDDO   
          ENDDO   
          DO NK=1,LTOP
            THTAG(NK)=THPA(NK)
            QG(NK)=QPA(NK)
          ENDDO
!     
!...CHECK TO SEE IF MIXING RATIO DIPS BELOW ZERO ANYWHERE;  IF SO, BORROW
!...MOISTURE FROM ADJACENT LAYERS TO BRING IT BACK UP ABOVE ZERO...
!     
        DO NK=1,LTOP
          IF(QG(NK).LT.0.)THEN
            IF(NK.EQ.1)THEN                             ! JSK MODS
!              PRINT *,' PROBLEM WITH KF SCHEME:  ' ! JSK MODS
!              PRINT *,'QG = 0 AT THE SURFACE!!!!!!!'    ! JSK MODS
               FATAL_ERROR('QG, QG(NK).LT.0')
            ENDIF
            NK1=NK+1
            IF(NK.EQ.LTOP)THEN
              NK1=KLCL
            ENDIF
            TMA=QG(NK1)*EMS(NK1)
            TMB=QG(NK-1)*EMS(NK-1)
            TMM=(QG(NK)-1.E-9)*EMS(NK  )
            BCOEFF=-TMM/((TMA*TMA)/TMB+TMB)
            ACOEFF=BCOEFF*TMA/TMB
            TMB=TMB*(1.-BCOEFF)
            TMA=TMA*(1.-ACOEFF)
            IF(NK.EQ.LTOP)THEN
              QVDIFF=(QG(NK1)-TMA*EMSD(NK1))*100./QG(NK1)
!              IF(ABS(QVDIFF).GT.1.)THEN
!             PRINT *,'!!!WARNING!!! CLOUD BASE WATER VAPOR CHANGES BY ',     &
!                      QVDIFF,                                                &
!                     '% WHEN MOISTURE IS BORROWED TO PREVENT NEGATIVE ',     &
!                     'VALUES IN KAIN-FRITSCH'
!              ENDIF
            ENDIF
            QG(NK)=1.E-9
            QG(NK1)=TMA*EMSD(NK1)
            QG(NK-1)=TMB*EMSD(NK-1)
          ENDIF
        ENDDO
        TOPOMG=(UDR(LTOP)-UER(LTOP))*DP(LTOP)*EMSD(LTOP)
        IF(ABS(TOPOMG-OMG(LTOP)).GT.1.E-3)THEN
!       WRITE(99,*)'ERROR:  MASS DOES NOT BALANCE IN KF SCHEME;            &
!      TOPOMG, OMG =',TOPOMG,OMG(LTOP)
!      TOPOMG, OMG =',TOPOMG,OMG(LTOP)
          ISTOP=1
          IPRNT=.TRUE.
          EXIT iter
        ENDIF
!     
!...CONVERT THETA TO T...
!     
        DO NK=1,LTOP
          EXN(NK)=(P00/P0(NK))**(0.2854*(1.-0.28*QG(NK)))
          TG(NK)=THTAG(NK)/EXN(NK)
          TVG(NK)=TG(NK)*(1.+0.608*QG(NK))
        ENDDO
        IF(ISHALL.EQ.1)THEN
          EXIT iter
        ENDIF
!     
!*******************************************************************
!                                                                  *
!     COMPUTE NEW CLOUD AND CHANGE IN AVAILABLE BUOYANT ENERGY.    *
!                                                                  *
!*******************************************************************
!     
!...THE FOLLOWING COMPUTATIONS ARE SIMILAR TO THAT FOR UPDRAFT
!     
!        THMIX=0.
          TMIX=0.
          QMIX=0.
!
!...FIND THE THERMODYNAMIC CHARACTERISTICS OF THE LAYER BY
!...MASS-WEIGHTING THE CHARACTERISTICS OF THE INDIVIDUAL MODEL
!...LAYERS...
!
          DO NK=LC,KPBL
            TMIX=TMIX+DP(NK)*TG(NK)
            QMIX=QMIX+DP(NK)*QG(NK)  
          ENDDO
          TMIX=TMIX/DPTHMX
          QMIX=QMIX/DPTHMX
          ES=ALIQ*EXP((TMIX*BLIQ-CLIQ)/(TMIX-DLIQ))
          QSS=0.622*ES/(PMIX-ES)
!     
!...REMOVE SUPERSATURATION FOR DIAGNOSTIC PURPOSES, IF NECESSARY...
!     
          IF(QMIX.GT.QSS)THEN
            RL=XLV0-XLV1*TMIX
            CPM=CP*(1.+0.887*QMIX)
            DSSDT=QSS*(CLIQ-BLIQ*DLIQ)/((TMIX-DLIQ)*(TMIX-DLIQ))
            DQ=(QMIX-QSS)/(1.+RL*DSSDT/CPM)
            TMIX=TMIX+RL/CP*DQ
            QMIX=QMIX-DQ
            TLCL=TMIX
          ELSE
            QMIX=AMAX1(QMIX,0.)
            EMIX=QMIX*PMIX/(0.622+QMIX)
            astrt=1.e-3
            binc=0.075
            a1=emix/aliq
            tp=(a1-astrt)/binc
            indlu=int(tp)+1
            value=(indlu-1)*binc+astrt
            aintrp=(a1-value)/binc
            tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
            TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
            TLCL=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(TMIX-T00))*(TMIX-TDPT)
            TLCL=AMIN1(TLCL,TMIX)
          ENDIF
          TVLCL=TLCL*(1.+0.608*QMIX)
          ZLCL = ZMIX+(TLCL-TMIX)/GDRY
          DO NK = LC,KL
            KLCL=NK
            IF(ZLCL.LE.Z0(NK))THEN
              EXIT 
            ENDIF
          ENDDO
          K=KLCL-1
          DLP=(ZLCL-Z0(K))/(Z0(KLCL)-Z0(K))
!     
!...ESTIMATE ENVIRONMENTAL TEMPERATURE AND MIXING RATIO AT THE LCL...
!     
          TENV=TG(K)+(TG(KLCL)-TG(K))*DLP
          QENV=QG(K)+(QG(KLCL)-QG(K))*DLP
          TVEN=TENV*(1.+0.608*QENV)
          PLCL=P0(K)+(P0(KLCL)-P0(K))*DLP
          THETEU(K)=TMIX*(1.E5/PMIX)**(0.2854*(1.-0.28*QMIX))*             &
                  EXP((3374.6525/TLCL-2.5403)*QMIX*(1.+0.81*QMIX))
!     
!...COMPUTE ADJUSTED ABE(ABEG).
!     
          ABEG=0.
          DO NK=K,LTOPM1
            NK1=NK+1

!new ckay adding FRZ effect 01-30-2015
         IF (aercu_opt.GT.0.0) THEN
            JK = KX-NK+1
            a1kay =  FRZ(1,JK)*DZA(NK)*3.337E5/CP
            a1kay = a1kay * ((1.E5/P0(NK))**ROCP)
            THETEU(NK)  =  a1kay + THETEU(NK)
         END IF
!ckay freezing effect included in ThetaU for cape calculation

            THETEU(NK1) = THETEU(NK)
!
            call tpmix2dd(p0(nk1),theteu(nk1),tgu(nk1),qgu(nk1))
!
            TVQU(NK1)=TGU(NK1)*(1.+0.608*QGU(NK1)-QLIQ(NK1)-QICE(NK1))
            IF(NK.EQ.K)THEN
              DZZ=Z0(KLCL)-ZLCL
              DILBE=((TVLCL+TVQU(NK1))/(TVEN+TVG(NK1))-1.)*DZZ
            ELSE
              DZZ=DZA(NK)
              DILBE=((TVQU(NK)+TVQU(NK1))/(TVG(NK)+TVG(NK1))-1.)*DZZ
            ENDIF
            IF(DILBE.GT.0.)ABEG=ABEG+DILBE*G
!
!...DILUTE BY ENTRAINMENT BY THE RATE AS ORIGINAL UPDRAFT...
!
            CALL ENVIRTHT(P0(NK1),TG(NK1),QG(NK1),THTEEG(NK1),ALIQ,BLIQ,CLIQ,DLIQ)
            THETEU(NK1)=THETEU(NK1)*DDILFRC(NK1)+THTEEG(NK1)*(1.-DDILFRC(NK1))
          ENDDO
!     
!...ASSUME AT LEAST 90% OF CAPE (ABE) IS REMOVED BY CONVECTION DURING
!...THE PERIOD TIMEC...
!     
          IF(NOITR.EQ.1)THEN
!         write(98,*)' '
!         write(98,*)'TAU, I, J, =',NTSD,I,J
!         WRITE(98,1060)FABE
!          GOTO 265
          EXIT iter
          ENDIF
          DABE=AMAX1(ABE-ABEG,capeDX*ABE)

          FABE=ABEG/ABE
          IF(FABE.GT.1. .and. ISHALL.EQ.0)THEN
!          WRITE(98,*)'UPDRAFT/DOWNDRAFT COUPLET INCREASES CAPE AT THIS
!     *GRID POINT; NO CONVECTION ALLOWED!'
            RETURN  
          ENDIF
          IF(NCOUNT.NE.1)THEN
            IF(ABS(AINC-AINCOLD).LT.0.0001)THEN
              NOITR=1
              AINC=AINCOLD
              CYCLE iter
            ENDIF
            DFDA=(FABE-FABEOLD)/(AINC-AINCOLD)
            IF(DFDA.GT.0.)THEN
              NOITR=1
              AINC=AINCOLD
              CYCLE iter
            ENDIF
          ENDIF
          AINCOLD=AINC
          FABEOLD=FABE
          IF(AINC/AINCMX.GT.0.999.AND.FABE.GT.1.05-STAB)THEN
!           write(98,*)' '
!           write(98,*)'TAU, I, J, =',NTSD,I,J
!           WRITE(98,1055)FABE
!            GOTO 265
            EXIT
          ENDIF
          IF((FABE.LE.1.05-STAB.AND.FABE.GE.0.95-STAB) .or. NCOUNT.EQ.10)THEN
            EXIT iter
          ELSE
            IF(NCOUNT.GT.10)THEN
!             write(98,*)' '
!             write(98,*)'TAU, I, J, =',NTSD,I,J
!             WRITE(98,1060)FABE
!             GOTO 265
              EXIT
            ENDIF
!     
!...IF MORE THAN 10% OF THE ORIGINAL CAPE REMAINS, INCREASE THE CONVECTIVE
!...MASS FLUX BY THE FACTOR AINC:
!     
            IF(FABE.EQ.0.)THEN
              AINC=AINC*0.5
            ELSE
              IF(DABE.LT.1.e-4)THEN
                NOITR=1
                AINC=AINCOLD
                CYCLE iter
              ELSE
                AINC=AINC*STAB*ABE/DABE
              ENDIF
            ENDIF
!           AINC=AMIN1(AINCMX,AINC)
            AINC=AMIN1(AINCMX,AINC)
!...IF AINC BECOMES VERY SMALL, EFFECTS OF CONVECTION ! JSK MODS
!...WILL BE MINIMAL SO JUST IGNORE IT...              ! JSK MODS
            IF(AINC.LT.0.05)then
              RETURN                          ! JSK MODS
            ENDIF
!            AINC=AMAX1(AINC,0.05)                        ! JSK MODS
            TDER=TDER2*AINC
            PPTFLX=PPTFL2*AINC
!           IF (XTIME.LT.10.)THEN
!           WRITE(98,1080)LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,
!          *              FABEOLD,AINCOLD 
!           ENDIF
            DO NK=1,LTOP
              UMF(NK)=UMF2(NK)*AINC
              DMF(NK)=DMF2(NK)*AINC
              DETLQ(NK)=DETLQ2(NK)*AINC
              DETIC(NK)=DETIC2(NK)*AINC
              UDR(NK)=UDR2(NK)*AINC
              UER(NK)=UER2(NK)*AINC
              DER(NK)=DER2(NK)*AINC
              DDR(NK)=DDR2(NK)*AINC
            ENDDO
!     
!...GO BACK UP FOR ANOTHER ITERATION...
!     
          ENDIF
        ENDDO iter

!ckay
! get the cloud fraction for layer NK+1=NK1
            updil = (100.-AINC)
            updil = updil/100.
           IF (aercu_opt .GT. 0) THEN
            ainc_frac1 = 1.0-updil !TWG
           END IF
            updil = updil*dxsq  
            Drag = 0.5   

        IF(ISHALL.EQ.1) THEN
          DO NK=KLCL, LTOP
            UMF_new = UMF(NK)/updil
            denSplume = P0(NK)/(R*TU(NK))
            xcldfra = 0.07*alog(1.+(500.*UMF_new))
            xcldfra = amax1(0.01,xcldfra)
            sh_cf1d(NK) = amin1(0.2,xcldfra)
!ckaywup
            DMF_new=DMF(NK)/updil
            FXM_new=FXM(NK)/dxsq
!           w_up1d(NK) = (UMF_new+DMF_new-FXM_new)/denSplume
!           w_up1d(NK) = w_up1d(NK)*Drag*DT/TIMEC
           w_up1d(NK) = (UMF_new/denSplume)*Drag*DT/TIMEC
          ENDDO
        ELSE 
          DO NK=KLCL, LTOP
! ww: moved the next line up
            UMF_new = UMF(NK)/updil
            denSplume = P0(NK)/(R*TU(NK))
            xcldfra = 0.14*alog(1.+(500.*UMF_new))
            xcldfra = amax1(0.01,xcldfra)
            dp_cf1d(NK) = amin1(0.6,xcldfra)
!new added downdraft impact
            DMF_new = DMF(NK)/updil
            FXM_new = FXM(NK)/dxsq
!           w_up1d(NK) = (UMF_new+DMF_new-FXM_new)/denSplume
!           w_up1d(NK) = w_up1d(NK)*Drag*DT/TIMEC
           w_up1d(NK) = (UMF_new/denSplume)*Drag*DT/TIMEC

          ENDDO
        ENDIF
!ckaywup
          envRHavg = 0.0
          DO NK=KLCL-1,LTOP1
           envEsat = 6.112*exp(17.67*(T0(NK)-273.16)/(243.5+T0(NK)-273.16))
           envEsat = envEsat * 100.0  ! to hPa
           envQsat = 0.622*envEsat/(P0(NK)-envEsat)
           envRH  = Q0(NK)/envQsat
            envRHavg = envRHavg + envRH
           if(envRH.gt.1.01) then
             w_up1d(NK) = 0.0
           end if
          END DO

!     
!...COMPUTE HYDROMETEOR TENDENCIES AS IS DONE FOR T, QV...
!     
!...FRC2 IS THE FRACTION OF TOTAL CONDENSATE      !  PPT FB MODS
!...GENERATED THAT GOES INTO PRECIPITIATION       !  PPT FB MODS
!
!  Redistribute hydormeteors according to the final mass-flux values:
!
        IF(CPR.GT.0.)THEN 
          FRC2=PPTFLX/(CPR*AINC)                    !  PPT FB MODS
        ELSE
           FRC2=0.
        ENDIF
        DO NK=1,LTOP
          QLPA(NK)=QL0(NK)
          QIPA(NK)=QI0(NK)
          QRPA(NK)=QR0(NK)
          QSPA(NK)=QS0(NK)
          RAINFB(NK)=PPTLIQ(NK)*AINC*FBFRC*FRC2   !  PPT FB MODS
          SNOWFB(NK)=PPTICE(NK)*AINC*FBFRC*FRC2   !  PPT FB MODS
        ENDDO
        DO NTC=1,NSTEP
!     
!...ASSIGN HYDROMETEORS CONCENTRATIONS AT THE TOP AND BOTTOM OF EACH LAYER
!...BASED ON THE SIGN OF OMEGA...
!     
          DO NK=1,LTOP
            QLFXIN(NK)=0.
            QLFXOUT(NK)=0.
            QIFXIN(NK)=0.
            QIFXOUT(NK)=0.
            QRFXIN(NK)=0.
            QRFXOUT(NK)=0.
            QSFXIN(NK)=0.
            QSFXOUT(NK)=0.
          ENDDO   
          DO NK=2,LTOP
            IF(OMG(NK).LE.0.)THEN
              QLFXIN(NK)=-FXM(NK)*QLPA(NK-1)
              QIFXIN(NK)=-FXM(NK)*QIPA(NK-1)
              QRFXIN(NK)=-FXM(NK)*QRPA(NK-1)
              QSFXIN(NK)=-FXM(NK)*QSPA(NK-1)
              QLFXOUT(NK-1)=QLFXOUT(NK-1)+QLFXIN(NK)
              QIFXOUT(NK-1)=QIFXOUT(NK-1)+QIFXIN(NK)
              QRFXOUT(NK-1)=QRFXOUT(NK-1)+QRFXIN(NK)
              QSFXOUT(NK-1)=QSFXOUT(NK-1)+QSFXIN(NK)
            ELSE
              QLFXOUT(NK)=FXM(NK)*QLPA(NK)
              QIFXOUT(NK)=FXM(NK)*QIPA(NK)
              QRFXOUT(NK)=FXM(NK)*QRPA(NK)
              QSFXOUT(NK)=FXM(NK)*QSPA(NK)
              QLFXIN(NK-1)=QLFXIN(NK-1)+QLFXOUT(NK)
              QIFXIN(NK-1)=QIFXIN(NK-1)+QIFXOUT(NK)
              QRFXIN(NK-1)=QRFXIN(NK-1)+QRFXOUT(NK)
              QSFXIN(NK-1)=QSFXIN(NK-1)+QSFXOUT(NK)
            ENDIF
          ENDDO   
!     
!...UPDATE THE HYDROMETEOR CONCENTRATION VALUES AT EACH LEVEL...
!     
          DO NK=1,LTOP
            QLPA(NK)=QLPA(NK)+(QLFXIN(NK)+DETLQ(NK)-QLFXOUT(NK))*DTIME*EMSD(NK)
            QIPA(NK)=QIPA(NK)+(QIFXIN(NK)+DETIC(NK)-QIFXOUT(NK))*DTIME*EMSD(NK)
            QRPA(NK)=QRPA(NK)+(QRFXIN(NK)-QRFXOUT(NK)+RAINFB(NK))*DTIME*EMSD(NK)         !  PPT FB MODS
            QSPA(NK)=QSPA(NK)+(QSFXIN(NK)-QSFXOUT(NK)+SNOWFB(NK))*DTIME*EMSD(NK)         !  PPT FB MODS
          ENDDO     
        ENDDO
        DO NK=1,LTOP
          QLG(NK)=QLPA(NK)
          QIG(NK)=QIPA(NK)
          QRG(NK)=QRPA(NK)
          QSG(NK)=QSPA(NK)
        ENDDO   

!kf_edrates
!Save convective timescale (TIMEC) as 2D variable
        IF (KF_EDRATES == 1) THEN
           DynTau=TIMEC
        ENDIF

!
!...CLEAN THINGS UP, CALCULATE CONVECTIVE FEEDBACK TENDENCIES FOR THIS
!...GRID POINT...
!     
!     IF (XTIME.LT.10.)THEN
!     WRITE(98,1080)LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,FABE,AINC 
!     ENDIF
       IF(IPRNT)THEN  
!         WRITE(98,1080)LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,FABE,AINC
         WRITE(message,1080)LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,FABE,AINC
         CALL wrf_message(message)
!        flush(98)   
       endif  
!     
!...SEND FINAL PARAMETERIZED VALUES TO OUTPUT FILES...
!     
!297   IF(IPRNT)then 
       IF(IPRNT)then 
!    if(I.eq.16 .and. J.eq.41)then
!      IF(ISTOP.EQ.1)THEN
!        write(98,*)
!        write(98,*)'At t(h), I, J =',float(NTSD)*72./3600.,I,J
         write(message,*)'P(LC), DTP, WKL, WKLCL =',p0(LC)/100.,       &
                     TLCL+DTLCL+dtrh-TENV,WKL,WKLCL
         call wrf_message(message)
         write(message,*)'TLCL, DTLCL, DTRH, TENV =',TLCL,DTLCL,       &
                      DTRH,TENV   
         call wrf_message(message)
         WRITE(message,1025)KLCL,ZLCL,DTLCL,LTOP,P0(LTOP),IFLAG,       &
         TMIX-T00,PMIX,QMIX,ABE
         call wrf_message(message)
         WRITE(message,1030)P0(LET)/100.,P0(LTOP)/100.,VMFLCL,PLCL/100.,  &
         WLCL,CLDHGT(LC)
         call wrf_message(message)
         WRITE(message,1035)PEF,PEFCBH,LC,LET,WKL,VWS 
         call wrf_message(message)
         write(message,*)'PRECIP EFFICIENCY =',PEFF 
         call wrf_message(message)
      WRITE(message,1080)LFS,LDB,LDT,TIMEC,TADVEC,NSTEP,NCOUNT,FABE,AINC
         call wrf_message(message)
!      ENDIF
!!!!! HERE !!!!!!!
           WRITE(message,1070)'  P  ','   DP ',' DT K/D ',' DR K/D ','   OMG  ',        &
          ' DOMGDP ','   UMF  ','   UER  ','   UDR  ','   DMF  ','   DER  '        &
          ,'   DDR  ','   EMS  ','    W0  ','  DETLQ ',' DETIC '
         call wrf_message(message)
           write(message,*)'just before DO 300...'
         call wrf_message(message)
!          flush(98)
           DO NK=1,LTOP
             K=LTOP-NK+1
             DTT=(TG(K)-T0(K))*86400./TIMEC
             RL=XLV0-XLV1*TG(K)
             DR=-(QG(K)-Q0(K))*RL*86400./(TIMEC*CP)
             UDFRC=UDR(K)*TIMEC*EMSD(K)
             UEFRC=UER(K)*TIMEC*EMSD(K)
             DDFRC=DDR(K)*TIMEC*EMSD(K)
             DEFRC=-DER(K)*TIMEC*EMSD(K)
             WRITE(message,1075)P0(K)/100.,DP(K)/100.,DTT,DR,OMG(K),DOMGDP(K)*1.E4,       &
             UMF(K)/1.E6,UEFRC,UDFRC,DMF(K)/1.E6,DEFRC,DDFRC,EMS(K)/1.E11,           &
             W0AVG1D(K)*1.E2,DETLQ(K)*TIMEC*EMSD(K)*1.E3,DETIC(K)*                   &
             TIMEC*EMSD(K)*1.E3
         call wrf_message(message)
           ENDDO
           WRITE(message,1085)'K','P','Z','T0','TG','DT','TU','TD','Q0','QG',             &
                  'DQ','QU','QD','QLG','QIG','QRG','QSG','RH0','RHG'
         call wrf_message(message)
           DO NK=1,KL
             K=KX-NK+1
             DTT=TG(K)-T0(K)
             TUC=TU(K)-T00
             IF(K.LT.LC.OR.K.GT.LTOP)TUC=0.
             TDC=TZ(K)-T00
             IF((K.LT.LDB.OR.K.GT.LDT).AND.K.NE.LFS)TDC=0.
             IF(T0(K).LT.T00)THEN
               ES=ALIQ*EXP((BLIQ*TG(K)-CLIQ)/(TG(K)-DLIQ))
             ELSE
               ES=ALIQ*EXP((BLIQ*TG(K)-CLIQ)/(TG(K)-DLIQ))
             ENDIF  
             QGS=ES*0.622/(P0(K)-ES)
             RH0=Q0(K)/QES(K)
             RHG=QG(K)/QGS
             WRITE(message,1090)K,P0(K)/100.,Z0(K),T0(K)-T00,TG(K)-T00,DTT,TUC,            &
             TDC,Q0(K)*1000.,QG(K)*1000.,(QG(K)-Q0(K))*1000.,QU(K)*                   &
             1000.,QD(K)*1000.,QLG(K)*1000.,QIG(K)*1000.,QRG(K)*1000.,                &
             QSG(K)*1000.,RH0,RHG
         call wrf_message(message)
           ENDDO
!     
!...IF CALCULATIONS ABOVE SHOW AN ERROR IN THE MASS BUDGET, PRINT OUT A
!...TO BE USED LATER FOR DIAGNOSTIC PURPOSES, THEN ABORT RUN...
!     
!         IF(ISTOP.EQ.1 .or. ISHALL.EQ.1)THEN

!         IF(ISHALL.NE.1)THEN
!            write(98,4421)i,j,iyr,imo,idy,ihr,imn
!           write(98)i,j,iyr,imo,idy,ihr,imn,kl
! 4421       format(7i4)
!            write(98,4422)kl
! 4422       format(i6) 
            DO 310 NK = 1,KL
              k = kl - nk + 1
!             write(98,4455) p0(k)/100.,t0(k)-273.16,q0(k)*1000.,       &
!                      u0(k),v0(k),W0AVG1D(K),dp(k),tke(k)
!             write(98) p0,t0,q0,u0,v0,w0,dp,tke
!           WRITE(98,1115)Z0(K),P0(K)/100.,T0(K)-273.16,Q0(K)*1000.,
!    *               U0(K),V0(K),DP(K)/100.,W0AVG(I,J,K)
 310        CONTINUE
            IF(ISTOP.EQ.1)THEN
              FATAL_ERROR('KAIN-FRITSCH, istop=1, diags')
            ENDIF
!         ENDIF
  4455  format(8f11.3) 
       ENDIF
        CNDTNF=(1.-EQFRC(LFS))*(QLIQ(LFS)+QICE(LFS))*DMF(LFS)
        PRATEC1=PPTFLX*(1.-FBFRC)/DXSQ
        RAINCV1=DT*PRATEC1     !  PPT FB MODS
!        RAINCV(I,J)=.1*.5*DT*PPTFLX/DXSQ               !  PPT FB MODS
!         RNC=0.1*TIMEC*PPTFLX/DXSQ
        RNC=RAINCV1*NIC
!       IF(ISHALL.EQ.0.AND.IPRNT)write (98,909)I,J,RNC

!     WRITE(98,1095)CPR*AINC,TDER+PPTFLX+CNDTNF
!     
!  EVALUATE MOISTURE BUDGET...
!     

        QINIT=0.
        QFNL=0.
        DPT=0.
        DO 315 NK=1,LTOP
          DPT=DPT+DP(NK)
          QINIT=QINIT+Q0(NK)*EMS(NK)
          QFNL=QFNL+QG(NK)*EMS(NK)
          QFNL=QFNL+(QLG(NK)+QIG(NK)+QRG(NK)+QSG(NK))*EMS(NK)
  315   CONTINUE
        QFNL=QFNL+PPTFLX*TIMEC*(1.-FBFRC)       !  PPT FB MODS
!        QFNL=QFNL+PPTFLX*TIMEC                 !  PPT FB MODS
        ERR2=(QFNL-QINIT)*100./QINIT
!       IF(IPRNT)WRITE(98,1110)QINIT,QFNL,ERR2
      IF(ABS(ERR2).GT.0.05 .AND. ISTOP.EQ.0)THEN 
!       write(99,*)'!!!!!!!! MOISTURE BUDGET ERROR IN KFPARA !!!'
!       WRITE(99,1110)QINIT,QFNL,ERR2
        IPRNT=.FALSE.
        ISTOP=1
!            write(98,4422)kl
 4422       format(i6)
            DO 311 NK = 1,KL
              k = kl - nk + 1
!             write(99,4455) p0(k)/100.,t0(k)-273.16,q0(k)*1000.,       &
!                      u0(k),v0(k),W0AVG1D(K),dp(k)
!             write(98) p0,t0,q0,u0,v0,w0,dp,tke
!           WRITE(98,1115)P0(K)/100.,T0(K)-273.16,Q0(K)*1000.,          &
!                    U0(K),V0(K),W0AVG1D(K),dp(k)/100.,tke(k)
!            WRITE(98,4456)P0(K)/100.,T0(K)-273.16,Q0(K)*1000.,          &
!                     U0(K),V0(K),W0AVG1D(K),dp(k)/100.,tke(k)
 311        CONTINUE
!           flush(98)

!        GOTO 297
!         STOP 'QVERR'
      ENDIF
 1115 FORMAT (2X,F7.2,2X,F5.1,2X,F6.3,2(2X,F5.1),2X,F7.2,2X,F7.4)
 4456  format(8f12.3)
        IF(PPTFLX.GT.0.)THEN
          RELERR=ERR2*QINIT/(PPTFLX*TIMEC)
        ELSE
          RELERR=0.
        ENDIF
     IF(IPRNT)THEN
!        WRITE(98,1120)RELERR
!        WRITE(98,*)'TDER, CPR, TRPPT =',              &
!          TDER,CPR*AINC,TRPPT*AINC
     ENDIF
!     
!...FEEDBACK TO RESOLVABLE SCALE TENDENCIES.
!     
!...IF THE ADVECTIVE TIME PERIOD (TADVEC) IS LESS THAN SPECIFIED MINIMUM
!...TIMEC, ALLOW FEEDBACK TO OCCUR ONLY DURING TADVEC...
!     
        IF(TADVEC.LT.TIMEC)NIC=NINT(TADVEC/DT)
        NCA1 = REAL(NIC)*DT
        IF(ISHALL.EQ.1)THEN
!          TIMEC = 2400.
           NCA1 = CUDT*60.
           NSHALL = NSHALL+1
        ENDIF

        DO K=1,KX
!         IF(IMOIST(INEST).NE.2)THEN
!
!...IF HYDROMETEORS ARE NOT ALLOWED, THEY MUST BE EVAPORATED OR SUBLIMATED
!...AND FED BACK AS VAPOR, ALONG WITH ASSOCIATED CHANGES IN TEMPERATURE.
!...NOTE:  THIS WILL INTRODUCE CHANGES IN THE CONVECTIVE TEMPERATURE AND
!...WATER VAPOR FEEDBACK TENDENCIES AND MAY LEAD TO SUPERSATURATED VALUE
!...OF QG...
!
!           RLC=XLV0-XLV1*TG(K)
!           RLS=XLS0-XLS1*TG(K)
!           CPM=CP*(1.+0.887*QG(K))
!           TG(K)=TG(K)-(RLC*(QLG(K)+QRG(K))+RLS*(QIG(K)+QSG(K)))/CPM
!           QG(K)=QG(K)+(QLG(K)+QRG(K)+QIG(K)+QSG(K))
!           DQLDT(I,J,NK)=0.
!           DQIDT(I,J,NK)=0.
!           DQRDT(I,J,NK)=0.
!           DQSDT(I,J,NK)=0.
!         ELSE
!
!...IF ICE PHASE IS NOT ALLOWED, MELT ALL FROZEN HYDROMETEORS...
!
          IF(warm_rain)THEN

            CPM=CP*(1.+0.887*QG(K))
            TG(K)=TG(K)-(QIG(K)+QSG(K))*RLF/CPM
            DQCDT(K)=(QLG(K)+QIG(K)-QL0(K)-QI0(K))/TIMEC
            DQIDT(K)=0.
            DQRDT(K)=(QRG(K)+QSG(K)-QR0(K)-QS0(K))/TIMEC
            DQSDT(K)=0.
          ELSEIF(.NOT. F_QS)THEN
!
!...IF ICE PHASE IS ALLOWED, BUT MIXED PHASE IS NOT, MELT FROZEN HYDROMETEORS
!...BELOW THE MELTING LEVEL, FREEZE LIQUID WATER ABOVE THE MELTING LEVEL
!
            CPM=CP*(1.+0.887*QG(K))
            IF(K.LE.ML)THEN
              TG(K)=TG(K)-(QIG(K)+QSG(K))*RLF/CPM
            ELSEIF(K.GT.ML)THEN
              TG(K)=TG(K)+(QLG(K)+QRG(K))*RLF/CPM
            ENDIF
            DQCDT(K)=(QLG(K)+QIG(K)-QL0(K)-QI0(K))/TIMEC
            DQIDT(K)=0.
            DQRDT(K)=(QRG(K)+QSG(K)-QR0(K)-QS0(K))/TIMEC
            DQSDT(K)=0.
          ELSEIF(F_QS) THEN
!
!...IF MIXED PHASE HYDROMETEORS ARE ALLOWED, FEED BACK CONVECTIVE TENDENCIES
!...OF HYDROMETEORS DIRECTLY...
!
            DQCDT(K)=(QLG(K)-QL0(K))/TIMEC
            DQSDT(K)=(QSG(K)-QS0(K))/TIMEC
            DQRDT(K)=(QRG(K)-QR0(K))/TIMEC
            IF (F_QI) THEN
               DQIDT(K)=(QIG(K)-QI0(K))/TIMEC
            ELSE
               DQSDT(K)=DQSDT(K)+(QIG(K)-QI0(K))/TIMEC
            ENDIF
          ELSE
!              PRINT *,'THIS COMBINATION OF IMOIST, IEXICE, IICE NOT ALLOWED!'
               FATAL_ERROR('KAIN-FRITSCH, THIS MICROPHYSICS CHOICE IS NOT ALLOWED')
          ENDIF
          DTDT(K)=(TG(K)-T0(K))/TIMEC
          DQDT(K)=(QG(K)-Q0(K))/TIMEC
        ENDDO

!JTR Begin CMT
        IF(cmt_opt_flag) THEN

              DO KQ=KTS,KTE
                 JK = KTE-KQ+1
                 DPconvF(1,JK)=0.0
                 DERconvF(1,JK)=0.0
                 UERconvF(1,JK)=0.0
                 MC(1,JK)=0.0
                 UMFconvF(1,JK)=0.0
                 PRU(1,JK)=0.0
                 QDN(1,JK)=0.0
                 QUP(1,JK)=0.0
                 QHAT(1,JK)=0.0
                 SDD(1,JK)=0.0
                 SUU(1,JK)=0.0
                 SHAT(1,JK)=0.0
                 TEE(1,JK)=0.0
                 U0F(1,JK)=0.0
                 V0F(1,JK)=0.0
                 Z0F(1,JK)=0.0
                 ZFU(1,JK)=0.0
              END DO

              zf_wrf(0) = 0.0  ! ground
              grav = 9.8
              cpin = CP
              Rdry = R
              Zfu(1,KTE+1) = 0.0
                 DTnew = DT
                 DSUBCLD(1) = ZLCL
                 VMFLCLconv(1) = ((VMFLCL/DXSQ)*G)/100.
                 IL1G = 1
                 IL2G = 1
                 ILG = 1
                 MSG1 = 0
!ckay
                 JBB(1) = KTE-KLCL+1  ! updraft base level   =====>>> flipped for CAM5 indexing
                 JDD(1) = KTE-LFS+1
                 JTT(1) = KTE-LTOP+1
!ckay
                  if(JDD(1).LT.JTT(1).or.JDD(1).GT.JBB(1)) then 
                      JDD(1)=JBB(1)-1  ! for cases no downdraft
                   end if
                  if(jtt(1).gt.jbb(1))  then
                   JTT(1) = JBB(1)
                  end if     

!JTR Begin CMT Variable Prep
! CKAY fill in grid-scale values for below cloud and above cloud portions
! Ckay and then fill in-cloud updraft and downdraft properties
               DO KQ=KTS,KTE
                 TDN(KQ) = T0(KQ)
                 TUP(KQ) = T0(KQ)
               end do


               !JTRnew: Added conditionals in case up/downdraft temps are 0.0
               DO KQ= KLCL,LTOP
                 if(TZ(KQ).NE.0.0) then
                        TDN(KQ) = TZ(KQ)
                 endif
                 if(TU(KQ).NE.0.0) then
                        TUP(KQ) = TU(KQ)
                 endif
               end do
               
               !JTRnew: Pulled this out of the main loop so the entire column
               !gets defined before use in the main loop
               DO KQ=KTS,KTE
                  zf_wrf(KQ) = zf_wrf(KQ-1)+DZQ(KQ)
                  stat_energy(KQ) = CP*(T0(KQ)*(1.0+0.622*Q0(KQ))) + G*Zf_wrf(KQ)
               ENDDO

               DO KQ=KTS,KTE
                 JK = KTE-KQ+1
                 DUDTnew(1,JK) = 0.0
                 DVDTnew(1,JK) = 0.0
                 DPDX(1,JK) = 0.0
                 DPDY(1,JK) = 0.0

                 zf_wrf(KQ) = zf_wrf(KQ-1)+DZQ(KQ)
                 zfu(1,JK) = zf_wrf(KQ) 

                 IF(KQ.EQ.KTE) THEN
                     qhat(1,JK) = Q0(KQ)
                     shat(1,JK) = stat_energy(KQ)
                 ELSE
                     qhat(1,JK) = 0.5*(Q0(KQ)/(1.+Q0(KQ)) + Q0(KQ+1)/(1.+Q0(KQ+1))) 
                     shat(1,JK) = 0.5*(stat_energy(KQ) + stat_energy(KQ+1))
                 ENDIF
                  
!ckay fill in clear and cloudy layers properly
               IF(KQ.LT.KLCL .OR. KQ.GT.LTOP) then  ! subcloud layer or above-cloud layer
                   QUP(1,JK) = Q0(KQ)/(1.+Q0(KQ))
                   QDN(1,JK) = Q0(KQ)/(1.+Q0(KQ))
                   SUU(1,JK) = CP*(T0(KQ)*(1.0+0.622*Q0(KQ))) + G*Zf_wrf(KQ)
                   SDD(1,JK) = CP*(T0(KQ)*(1.0+0.622*Q0(KQ))) + G*Zf_wrf(KQ)
               Else  
                 QUP(1,JK) = QU(KQ)/(1.+QU(KQ))
                 QDN(1,JK) = QD(KQ)/(1.+QD(KQ))

                 !JTRnew: Replaced TU and TZ with TUP and TDN
                 !Ckay replaced Qu and Qd with QUP & QDN
                 SUU(1,JK) = CP*(TUP(KQ)*(1.0+0.622*QUP(1,JK))) + G*Zf_wrf(KQ)
                 SDD(1,JK) = CP*(TDN(KQ)*(1.0+0.622*QDN(1,JK))) + G*Zf_wrf(KQ)
               End if ! for subcloud and above-cloud layers
                 DERconvF(1,JK) = (((DER(KQ)/DXSQ)*G)/100.)/DZQ(KQ)
                 UERconvF(1,JK) = (((UER(KQ)/DXSQ)*G)/100.)/DZQ(KQ)
                 DMFconvF(1,JK) = ((DMF(KQ)/DXSQ)*G)/100.
                 !JTR Updraft mass flux from kg/s to hpa/s
                 UMFconvF(1,JK) = ((UMF(KQ)/DXSQ)*G)/100.
                 MC(1,JK) = DMFconvF(1,JK) + UMFconvF(1,JK)
                 DPconvF(1,JK) = DP(KQ)/100.
                 PRU(1,JK) = P0(KQ)/100.0    ! in millibars or hPa
                 U0F(1,JK) = U0(KQ)
                 V0F(1,JK) = V0(KQ)
                 Z0F(1,JK) = Z0(KQ)
                 TEE(1,JK) = T0(KQ)
               END DO ! for k loop

!JTR End CMT Variable Prep
              if(CMTprint) then
                print *,'DP',minval(DPconvF),maxval(DPconvF)
                print *,'DER',minval(DERconvF),maxval(DERconvF)
                print *,'UER',minval(UERconvF),maxval(UERconvF)
                print *,'MC',minval(MC),maxval(MC)
                print *,'DMF',minval(DMFconvF),maxval(DMFconvF)
                print *,'UMF',minval(UMFconvF),maxval(UMFconvF)
                print *,'VMFLCL',VMFLCLconv(ILG)
                print *,'PRU',minval(PRU),maxval(PRU)
                print *,'QDD',maxval(QDn)
                print *,'QUU',minval(QUp),maxval(QUp)
                print *,'QHAT',minval(QHAT),maxval(QHAT)
                print *,'SDD',minval(SDD),maxval(SDD)
                print *,'SUU',minval(SUU),maxval(SUU)
                print *,'SHAT',minval(SHAT),maxval(SHAT)
                print *,'TEE',minval(TEE),maxval(TEE)
                print *,'U0F',minval(U0F),maxval(U0F)
                print *,'V0F',minval(V0F),maxval(V0F)
                print *,'Z0F',minval(Z0F),maxval(Z0F)
                print *,'ZFU',minval(ZFU),maxval(ZFU)
                print *,'DSUBCLD',DSUBCLD(ILG)
                print *,'JBB',JBB(ILG)
                print *,'JDD',JDD(ILG)
                print *,'JTT',JTT(ILG)
                print *,'msg1',msg1
                print *,'DT',DTnew
                print *,'grav',grav
                print *,'cpin',CPIN
                print *,'rdry',Rdry
                print *,'KTE',KTE
                print *,'IL1G',IL1G
                print *,'IL2G',IL2G
                print *,'ILG',ILG
              end if ! CMTpring
           
                CALL MSKF_CMT(DUDTnew,DVDTnew,dpdx,dpdy,     &
                       DPconvF,DERconvF,UERconvF,          &
                       MC,DMFconvF,                        &
                       UMFconvF,VMFLCLconv,PRU,       &
                       QDN,QUP,QHAT,SDD,SUU,SHAT,TEE,U0F,  &
                       V0F,Z0F,ZFU,DSUBCLD,JBB,  &
                       JDD,JTT,msg1,DTnew,       &
                       grav,CPIN,Rdry,KTE,IL1G,IL2G,ILG)
           !Invert tendency arrays
           DO KQ=kts,kte
              JK = KTE-KQ+1
              DUDT(KQ) =  DUDTnew(1,JK)
              DVDT(KQ) =  DVDTnew(1,JK)
           ENDDO

              if(CMTprint) then
           print *,'max/min dudt=',maxval(DUDT), minval(DUDT)
           print *,'max/min dVdt=',maxval(DVDT), minval(DVDT)
              end if ! for CMTprint
        ENDIF  ! for cmt flag
!JTR End CMT
        
        PRATEC1=PPTFLX*(1.-FBFRC)/DXSQ
        RAINCV1=DT*PRATEC1
!        RAINCV(I,J)=.1*.5*DT*PPTFLX/DXSQ               !  PPT FB MODS
!         RNC=0.1*TIMEC*PPTFLX/DXSQ
        RNC=RAINCV1*NIC
 909     FORMAT('AT I, J =',i3,1x,i3,' CONVECTIVE RAINFALL =',F8.4,' mm')
!      write (98,909)I,J,RNC
!      write (6,909)I,J,RNC
!      WRITE(98,*)'at NTSD =',NTSD,',No. of KF points activated =',
!     *            NCCNT
!      flush(98)
1000  FORMAT(' ',10A8)
1005  FORMAT(' ',F6.0,2X,F6.4,2X,F7.3,1X,F6.4,2X,4(F6.3,2X),2(F7.3,1X))
1010  FORMAT(' ',' VERTICAL VELOCITY IS NEGATIVE AT ',F4.0,' MB')
1015   FORMAT(' ','ALL REMAINING MASS DETRAINS BELOW ',F4.0,' MB')
1025   FORMAT(5X,' KLCL=',I2,' ZLCL=',F7.1,'M',                         &
        ' DTLCL=',F5.2,' LTOP=',I2,' P0(LTOP)=',-2PF5.1,'MB FRZ LV=',   &
        I2,' TMIX=',0PF4.1,1X,'PMIX=',-2PF6.1,' QMIX=',3PF5.1,          &
        ' CAPE=',0PF7.1)
1030   FORMAT(' ',' P0(LET) = ',F6.1,' P0(LTOP) = ',F6.1,' VMFLCL =',   &
      E12.3,' PLCL =',F6.1,' WLCL =',F6.3,' CLDHGT =',                  &
      F8.1)
1035  FORMAT(1X,'PEF(WS)=',F4.2,'(CB)=',F4.2,'LC,LET=',2I3,'WKL='       &
      ,F6.3,'VWS=',F5.2)
!1055  FORMAT('*** DEGREE OF STABILIZATION =',F5.3,                  &
!      ', NO MORE MASS FLUX IS ALLOWED!')
!1060     FORMAT(' ITERATION DOES NOT CONVERGE TO GIVE THE SPECIFIED    &
!      &DEGREE OF STABILIZATION!  FABE= ',F6.4) 
 1070 FORMAT (16A8) 
 1075 FORMAT (F8.2,3(F8.2),2(F8.3),F8.2,2F8.3,F8.2,6F8.3) 
 1080 FORMAT(2X,'LFS,LDB,LDT =',3I3,' TIMEC, TADVEC, NSTEP=',           &
              2(1X,F5.0),I3,'NCOUNT, FABE, AINC=',I2,1X,F5.3,F6.2) 
 1085 FORMAT (A3,16A7,2A8) 
 1090 FORMAT (I3,F7.2,F7.0,10F7.2,4F7.3,2F8.3) 
 1095 FORMAT(' ','  PPT PRODUCTION RATE= ',F10.0,' TOTAL EVAP+PPT= ',F10.0)
1105   FORMAT(' ','NET LATENT HEAT RELEASE =',E12.5,' ACTUAL HEATING =',&
       E12.5,' J/KG-S, DIFFERENCE = ',F9.3,'%')
1110   FORMAT(' ','INITIAL WATER =',E12.5,' FINAL WATER =',E12.5,       &
       ' TOTAL WATER CHANGE =',F8.2,'%')
! 1115 FORMAT (2X,F6.0,2X,F7.2,2X,F5.1,2X,F6.3,2(2X,F5.1),2X,F7.2,2X,F7.4)
1120   FORMAT(' ','MOISTURE ERROR AS FUNCTION OF TOTAL PPT =',F9.3,'%')
!
!-----------------------------------------------------------------------
!--------------SAVE CLOUD TOP AND BOTTOM FOR RADIATION------------------
!-----------------------------------------------------------------------
!
      CUTOP1=REAL(LTOP)
      CUBOT1=REAL(LCL)
!
!-----------------------------------------------------------------------
   END SUBROUTINE  MSKF_eta_PARA
!********************************************************************
!-----------------------------------------------------------
! ***********************************************************************
!dkay
!dkay: added QSu as output to get saturated Q of updraft

   SUBROUTINE TPMIX2(p,thes,tu,qu,qliq,qice,qnewlq,qnewic,XLV1,XLV0,Qsu)
!
! Lookup table variables:
!     INTEGER, PARAMETER :: (KFNT=250,KFNP=220)
!     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
!     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
!     REAL, SAVE, DIMENSION(1:200) :: ALU
!     REAL, SAVE :: RDPR,RDTHK,PLUTOP
! End of Lookup table variables:
!-----------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memory issues
!-----------------------------------------------------------------------
   REAL,         INTENT(IN   )   :: P,THES,XLV1,XLV0
   REAL,         INTENT(OUT  )   :: QNEWLQ,QNEWIC,QSu
   REAL,         INTENT(INOUT)   :: TU,QU,QLIQ,QICE
   REAL    ::    TP,QQ,BTH,TTH,PP,T00,T10,T01,T11,Q00,Q10,Q01,Q11,          &
                 TEMP,QS,QNEW,DQ,QTOT,RLL,CPP
   INTEGER ::    IPTB,ITHTB
!-----------------------------------------------------------------------

!c******** LOOKUP TABLE VARIABLES... ****************************
!      parameter(kfnt=250,kfnp=220)
!c
!      COMMON/KFLUT/ ttab(kfnt,kfnp),qstab(kfnt,kfnp),the0k(kfnp),
!     *              alu(200),rdpr,rdthk,plutop 
!C*************************************************************** 
!c
!c***********************************************************************
!c     scaling pressure and tt table index                         
!c***********************************************************************
!c
      tp=(p-plutop)*rdpr
      qq=tp-aint(tp)
      iptb=int(tp)+1

!
!***********************************************************************
!              base and scaling factor for the                           
!***********************************************************************
!
!  scaling the and tt table index                                        
      bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
      tth=(thes-bth)*rdthk
      pp   =tth-aint(tth)
      ithtb=int(tth)+1
       IF(IPTB.GE.220 .OR. IPTB.LE.1 .OR. ITHTB.GE.250 .OR. ITHTB.LE.1)THEN
!         write(98,*)'**** OUT OF BOUNDS *********'
!        flush(98)
       ENDIF
!
      t00=ttab(ithtb  ,iptb  )
      t10=ttab(ithtb+1,iptb  )
      t01=ttab(ithtb  ,iptb+1)
      t11=ttab(ithtb+1,iptb+1)
!
      q00=qstab(ithtb  ,iptb  )
      q10=qstab(ithtb+1,iptb  )
      q01=qstab(ithtb  ,iptb+1)
      q11=qstab(ithtb+1,iptb+1)
!
!***********************************************************************
!              parcel temperature                                        
!***********************************************************************
!
      temp=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
!
      qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)

!dkay
       QSu = qs
!

!
      DQ=QS-QU
      IF(DQ.LE.0.)THEN
        QNEW=QU-QS
        QU=QS
      ELSE 
!
!   IF THE PARCEL IS SUBSATURATED, TEMPERATURE AND MIXING RATIO MUST BE
!   ADJUSTED...IF LIQUID WATER IS PRESENT, IT IS ALLOWED TO EVAPORATE
! 
        QNEW=0.
        QTOT=QLIQ+QICE
!
!   IF THERE IS ENOUGH LIQUID OR ICE TO SATURATE THE PARCEL, TEMP STAYS AT ITS
!   WET BULB VALUE, VAPOR MIXING RATIO IS AT SATURATED LEVEL, AND THE MIXING
!   RATIOS OF LIQUID AND ICE ARE ADJUSTED TO MAKE UP THE ORIGINAL SATURATION
!   DEFICIT... OTHERWISE, ANY AVAILABLE LIQ OR ICE VAPORIZES AND APPROPRIATE
!   ADJUSTMENTS TO PARCEL TEMP; VAPOR, LIQUID, AND ICE MIXING RATIOS ARE MADE.
!
!...subsaturated values only occur in calculations involving various mixtures of
!...updraft and environmental air for estimation of entrainment and detrainment.
!...For these purposes, assume that reasonable estimates can be given using 
!...liquid water saturation calculations only - i.e., ignore the effect of the
!...ice phase in this process only...will not affect conservative properties...
!
        IF(QTOT.GE.DQ)THEN
          qliq=qliq-dq*qliq/(qtot+1.e-10)
          qice=qice-dq*qice/(qtot+1.e-10)
          QU=QS
        ELSE
          RLL=XLV0-XLV1*TEMP
          CPP=1004.5*(1.+0.89*QU)
          IF(QTOT.LT.1.E-10)THEN
!
!...IF NO LIQUID WATER OR ICE IS AVAILABLE, TEMPERATURE IS GIVEN BY:
            TEMP=TEMP+RLL*(DQ/(1.+DQ))/CPP
          ELSE
!
!...IF SOME LIQ WATER/ICE IS AVAILABLE, BUT NOT ENOUGH TO ACHIEVE SATURATION,
!   THE TEMPERATURE IS GIVEN BY:
!
            TEMP=TEMP+RLL*((DQ-QTOT)/(1+DQ-QTOT))/CPP
            QU=QU+QTOT
            QTOT=0.
            QLIQ=0.
            QICE=0.
          ENDIF
        ENDIF
      ENDIF
      TU=TEMP
      qnewlq=qnew
      qnewic=0.
!
   END SUBROUTINE TPMIX2
!******************************************************************************
      SUBROUTINE DTFRZNEW(TU,P,THTEU,QU,QFRZ,QICE,ALIQ,BLIQ,CLIQ,DLIQ)
!-----------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 Add to avoid memory issues
!-----------------------------------------------------------------------
   REAL,         INTENT(IN   )   :: P,QFRZ,ALIQ,BLIQ,CLIQ,DLIQ
   REAL,         INTENT(INOUT)   :: TU,THTEU,QU,QICE
   REAL    ::    RLC,RLS,RLF,CPP,A,DTFRZ,ES,QS,DQEVAP,PII
!-----------------------------------------------------------------------
!
!...ALLOW THE FREEZING OF LIQUID WATER IN THE UPDRAFT TO PROCEED AS AN 
!...APPROXIMATELY LINEAR FUNCTION OF TEMPERATURE IN THE TEMPERATURE RANGE 
!...TTFRZ TO TBFRZ...
!...FOR COLDER TEMPERATURES, FREEZE ALL LIQUID WATER...
!...THERMODYNAMIC PROPERTIES ARE STILL CALCULATED WITH RESPECT TO LIQUID WATER
!...TO ALLOW THE USE OF LOOKUP TABLE TO EXTRACT TMP FROM THETAE...
!
      RLC=2.5E6-2369.276*(TU-273.16)
      RLS=2833922.-259.532*(TU-273.16)
      RLF=RLS-RLC
      CPP=1004.5*(1.+0.89*QU)
!
!  A = D(es)/DT IS THAT CALCULATED FROM BUCK (1981) EMPERICAL FORMULAS
!  FOR SATURATION VAPOR PRESSURE...
!
      A=(CLIQ-BLIQ*DLIQ)/((TU-DLIQ)*(TU-DLIQ))
      DTFRZ = RLF*QFRZ/(CPP+RLS*QU*A)
      TU = TU+DTFRZ
      
      ES = ALIQ*EXP((BLIQ*TU-CLIQ)/(TU-DLIQ))
      QS = ES*0.622/(P-ES)
!
!...FREEZING WARMS THE AIR AND IT BECOMES UNSATURATED...ASSUME THAT SOME OF THE 
!...LIQUID WATER THAT IS AVAILABLE FOR FREEZING EVAPORATES TO MAINTAIN SATURA-
!...TION...SINCE THIS WATER HAS ALREADY BEEN TRANSFERRED TO THE ICE CATEGORY,
!...SUBTRACT IT FROM ICE CONCENTRATION, THEN SET UPDRAFT MIXING RATIO AT THE NEW
!...TEMPERATURE TO THE SATURATION VALUE...
!
      DQEVAP = QS-QU
      QICE = QICE-DQEVAP
      QU = QU+DQEVAP
      PII=(1.E5/P)**(0.2854*(1.-0.28*QU))
      THTEU=TU*PII*EXP((3374.6525/TU-2.5403)*QU*(1.+0.81*QU))
!
   END SUBROUTINE DTFRZNEW
! --------------------------------------------------------------------------------

      SUBROUTINE CONDLOAD(QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,RATE,QNEWLQ,           &
                          QNEWIC,QLQOUT,QICOUT,G)

!-----------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memory issues
!-----------------------------------------------------------------------
!  9/18/88...THIS PRECIPITATION FALLOUT SCHEME IS BASED ON THE SCHEME US
!  BY OGURA AND CHO (1973).  LIQUID WATER FALLOUT FROM A PARCEL IS CAL-
!  CULATED USING THE EQUATION DQ=-RATE*Q*DT, BUT TO SIMULATE A QUASI-
!  CONTINUOUS PROCESS, AND TO ELIMINATE A DEPENDENCY ON VERTICAL
!  RESOLUTION THIS IS EXPRESSED AS Q=Q*EXP(-RATE*DZ).

      REAL, INTENT(IN   )   :: G
      REAL, INTENT(IN   )   :: DZ,BOTERM,ENTERM,RATE
      REAL, INTENT(INOUT)   :: QLQOUT,QICOUT,WTW,QLIQ,QICE,QNEWLQ,QNEWIC
      REAL :: QTOT,QNEW,QEST,G1,WAVG,CONV,RATIO3,OLDQ,RATIO4,DQ,PPTDRG
!
      QTOT=QLIQ+QICE                                                    
      QNEW=QNEWLQ+QNEWIC                                                
!                                                                       
!  ESTIMATE THE VERTICAL VELOCITY SO THAT AN AVERAGE VERTICAL VELOCITY 
!  BE CALCULATED TO ESTIMATE THE TIME REQUIRED FOR ASCENT BETWEEN MODEL 
!  LEVELS...                                                            
!                                                                       
      QEST=0.5*(QTOT+QNEW)                                              
      G1=WTW+BOTERM-ENTERM-2.*G*DZ*QEST/1.5                             
      IF(G1.LT.0.0)G1=0.                                                
      WAVG=0.5*(SQRT(WTW)+SQRT(G1))                                      
      CONV=RATE*DZ/WAVG               ! KF90  Eq. 9
!                                                                       
!  RATIO3 IS THE FRACTION OF LIQUID WATER IN FRESH CONDENSATE, RATIO4 IS
!  THE FRACTION OF LIQUID WATER IN THE TOTAL AMOUNT OF CONDENSATE INVOLV
!  IN THE PRECIPITATION PROCESS - NOTE THAT ONLY 60% OF THE FRESH CONDEN
!  SATE IS IS ALLOWED TO PARTICIPATE IN THE CONVERSION PROCESS...       
!                                                                       
      RATIO3=QNEWLQ/(QNEW+1.E-8)                                       
!     OLDQ=QTOT                                                         
      QTOT=QTOT+0.6*QNEW                                                
      OLDQ=QTOT                                                         
      RATIO4=(0.6*QNEWLQ+QLIQ)/(QTOT+1.E-8)                            
      QTOT=QTOT*EXP(-CONV)            ! KF90  Eq. 9                                              
!                                                                       
!  DETERMINE THE AMOUNT OF PRECIPITATION THAT FALLS OUT OF THE UPDRAFT  
!  PARCEL AT THIS LEVEL...                                              
!                                                                       
      DQ=OLDQ-QTOT                                                      
      QLQOUT=RATIO4*DQ                                                  
      QICOUT=(1.-RATIO4)*DQ                                             
!                                                                       
!  ESTIMATE THE MEAN LOAD OF CONDENSATE ON THE UPDRAFT IN THE LAYER, CAL
!  LATE VERTICAL VELOCITY                                               
!                                                                       
      PPTDRG=0.5*(OLDQ+QTOT-0.2*QNEW)                                   
      WTW=WTW+BOTERM-ENTERM-2.*G*DZ*PPTDRG/1.5                          
      IF(ABS(WTW).LT.1.E-4)WTW=1.E-4
!                                                                       
!  DETERMINE THE NEW LIQUID WATER AND ICE CONCENTRATIONS INCLUDING LOSSE
!  DUE TO PRECIPITATION AND GAINS FROM CONDENSATION...                  
!                                                                       
      QLIQ=RATIO4*QTOT+RATIO3*0.4*QNEW                                  
      QICE=(1.-RATIO4)*QTOT+(1.-RATIO3)*0.4*QNEW                        
      QNEWLQ=0.                                                         
      QNEWIC=0.                                                         

   END SUBROUTINE CONDLOAD

! ----------------------------------------------------------------------
   SUBROUTINE PROF5(EQ,EE,UD)                                        
!
!***********************************************************************
!*****    GAUSSIAN TYPE MIXING PROFILE....******************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  THIS SUBROUTINE INTEGRATES THE AREA UNDER THE CURVE IN THE GAUSSIAN  
!  DISTRIBUTION...THE NUMERICAL APPROXIMATION TO THE INTEGRAL IS TAKEN FROM
!  "HANDBOOK OF MATHEMATICAL FUNCTIONS WITH FORMULAS, GRAPHS AND MATHEMATICS TABLES"
!  ED. BY ABRAMOWITZ AND STEGUN, NATL BUREAU OF STANDARDS APPLIED
!  MATHEMATICS SERIES.  JUNE, 1964., MAY, 1968.                         
!                                     JACK KAIN                         
!                                     7/6/89                            
!  Solves for KF90 Eq. 2
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!-----------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memory issues
!-----------------------------------------------------------------------
   REAL,         INTENT(IN   )   :: EQ
   REAL,         INTENT(INOUT)   :: EE,UD
   REAL ::       SQRT2P,A1,A2,A3,P,SIGMA,FE,X,Y,EY,E45,T1,T2,C1,C2

      DATA SQRT2P,A1,A2,A3,P,SIGMA,FE/2.506628,0.4361836,-0.1201676,       &
           0.9372980,0.33267,0.166666667,0.202765151/                        
      X=(EQ-0.5)/SIGMA                                                  
      Y=6.*EQ-3.                                                        
      EY=EXP(Y*Y/(-2))                                                  
      E45=EXP(-4.5)                                                     
      T2=1./(1.+P*ABS(Y))                                               
      T1=0.500498                                                       
      C1=A1*T1+A2*T1*T1+A3*T1*T1*T1                                     
      C2=A1*T2+A2*T2*T2+A3*T2*T2*T2                                     
      IF(Y.GE.0.)THEN                                                   
        EE=SIGMA*(0.5*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY))-E45*EQ*EQ/2.
        UD=SIGMA*(0.5*(EY*C2-E45*C1)+SIGMA*(E45-EY))-E45*(0.5+EQ*EQ/2.-    &
           EQ)                                                          
      ELSE                                                              
        EE=SIGMA*(0.5*(EY*C2-E45*C1)+SIGMA*(E45-EY))-E45*EQ*EQ/2.       
        UD=SIGMA*(0.5*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY))-E45*(0.5+EQ*   &
           EQ/2.-EQ)                                                    
      ENDIF                                                             
      EE=EE/FE                                                          
      UD=UD/FE                                                          

   END SUBROUTINE PROF5

! ------------------------------------------------------------------------
   SUBROUTINE TPMIX2DD(p,thes,ts,qs)
!
! Lookup table variables:
!     INTEGER, PARAMETER :: (KFNT=250,KFNP=220)
!     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
!     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
!     REAL, SAVE, DIMENSION(1:200) :: ALU
!     REAL, SAVE :: RDPR,RDTHK,PLUTOP
! End of Lookup table variables:
!-----------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memory issues
!-----------------------------------------------------------------------
   REAL,         INTENT(IN   )   :: P,THES
   REAL,         INTENT(INOUT)   :: TS,QS
   REAL    ::    TP,QQ,BTH,TTH,PP,T00,T10,T01,T11,Q00,Q10,Q01,Q11
   INTEGER ::    IPTB,ITHTB
   CHARACTER*256 :: MESS
!-----------------------------------------------------------------------

!
!******** LOOKUP TABLE VARIABLES (F77 format)... ****************************
!     parameter(kfnt=250,kfnp=220)
!
!     COMMON/KFLUT/ ttab(kfnt,kfnp),qstab(kfnt,kfnp),the0k(kfnp),        &
!                   alu(200),rdpr,rdthk,plutop 
!*************************************************************** 
!
!***********************************************************************
!     scaling pressure and tt table index                         
!***********************************************************************
!
      tp=(p-plutop)*rdpr
      qq=tp-aint(tp)
      iptb=int(tp)+1
!
!***********************************************************************
!              base and scaling factor for the                           
!***********************************************************************
!
!  scaling the and tt table index                                        
      bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
      tth=(thes-bth)*rdthk
      pp   =tth-aint(tth)
      ithtb=int(tth)+1
!
      t00=ttab(ithtb  ,iptb  )
      t10=ttab(ithtb+1,iptb  )
      t01=ttab(ithtb  ,iptb+1)
      t11=ttab(ithtb+1,iptb+1)
!
      q00=qstab(ithtb  ,iptb  )
      q10=qstab(ithtb+1,iptb  )
      q01=qstab(ithtb  ,iptb+1)
      q11=qstab(ithtb+1,iptb+1)
!
!***********************************************************************
!              parcel temperature and saturation mixing ratio                                        
!***********************************************************************
!
      ts=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
!
      qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)
!
   END SUBROUTINE TPMIX2DD

! -----------------------------------------------------------------------
  SUBROUTINE ENVIRTHT(P1,T1,Q1,THT1,ALIQ,BLIQ,CLIQ,DLIQ)                       
!
!-----------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memory issues
!-----------------------------------------------------------------------
   REAL,         INTENT(IN   )   :: P1,T1,Q1,ALIQ,BLIQ,CLIQ,DLIQ
   REAL,         INTENT(INOUT)   :: THT1
   REAL    ::    EE,TLOG,ASTRT,AINC,A1,TP,VALUE,AINTRP,TDPT,TSAT,THT,      &
                 T00,P00,C1,C2,C3,C4,C5
   INTEGER ::    INDLU
!-----------------------------------------------------------------------
      DATA T00,P00,C1,C2,C3,C4,C5/273.16,1.E5,3374.6525,2.5403,3114.834,   &
           0.278296,1.0723E-3/                                          
!                                                                       
!  CALCULATE ENVIRONMENTAL EQUIVALENT POTENTIAL TEMPERATURE...          
!                                                                       
! NOTE: Calculations for mixed/ice phase no longer used...jsk 8/00
!        For example, KF90 Eq. 10 no longer used
!
      EE=Q1*P1/(0.622+Q1)                                             
!     TLOG=ALOG(EE/ALIQ)                                              
! ...calculate LOG term using lookup table...
!
      astrt=1.e-3
      ainc=0.075
      a1=ee/aliq
      tp=(a1-astrt)/ainc
      indlu=int(tp)+1
      value=(indlu-1)*ainc+astrt
      aintrp=(a1-value)/ainc
      tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
!
      TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)                               
      TSAT=TDPT-(.212+1.571E-3*(TDPT-T00)-4.36E-4*(T1-T00))*(T1-TDPT) 
      THT=T1*(P00/P1)**(0.2854*(1.-0.28*Q1))                          
      THT1=THT*EXP((C1/TSAT-C2)*Q1*(1.+0.81*Q1))                      
!
  END SUBROUTINE ENVIRTHT                                                              
! ***********************************************************************
!====================================================================
   SUBROUTINE mskf_init(RTHCUTEN,RQVCUTEN,RQCCUTEN,RQRCUTEN,        &
                     RQICUTEN,RQSCUTEN,NCA,W0AVG,P_QI,P_QS,         &
                     SVP1,SVP2,SVP3,SVPT0,                          &
                     P_FIRST_SCALAR,restart,allowed_to_read,        &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte,                  &
                     RUCUTEN, RVCUTEN                               ) !JTR
!--------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memeory issues
!--------------------------------------------------------------------
   LOGICAL , INTENT(IN)           ::  restart,allowed_to_read
   INTEGER , INTENT(IN)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte
   INTEGER , INTENT(IN)           ::  P_QI,P_QS,P_FIRST_SCALAR


   REAL,     DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) ::       &
                                                          RTHCUTEN, &
                                                          RQVCUTEN, &
                                                          RQCCUTEN, &
                                                          RQRCUTEN, &
                                                          RQICUTEN, &
                                                          RQSCUTEN, &
                                                           RUCUTEN, &
                                                           RVCUTEN

   REAL ,   DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: W0AVG

   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT):: NCA

   INTEGER :: i, j, k, itf, jtf, ktf
   REAL, INTENT(IN)    :: SVP1,SVP2,SVP3,SVPT0

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   IF(.not.restart)THEN

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         RTHCUTEN(i,k,j)=0.
         RQVCUTEN(i,k,j)=0.
         RQCCUTEN(i,k,j)=0.
         RQRCUTEN(i,k,j)=0.
         !JTR Momentum tendencies
         RUCUTEN(i,k,j)=0.
         RVCUTEN(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO

      IF (P_QI .ge. P_FIRST_SCALAR) THEN
         DO j=jts,jtf
         DO k=kts,ktf
         DO i=its,itf
            RQICUTEN(i,k,j)=0.
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QS .ge. P_FIRST_SCALAR) THEN
         DO j=jts,jtf
         DO k=kts,ktf
         DO i=its,itf
            RQSCUTEN(i,k,j)=0.
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      DO j=jts,jtf
      DO i=its,itf
         NCA(i,j)=-100.
      ENDDO
      ENDDO

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         W0AVG(i,k,j)=0.
      ENDDO
      ENDDO
      ENDDO

   endif
 
   CALL MSKF_LUTAB(SVP1,SVP2,SVP3,SVPT0)


   END SUBROUTINE mskf_init

!-------------------------------------------------------

      subroutine mskf_lutab(SVP1,SVP2,SVP3,SVPT0)
!
!  This subroutine is a lookup table.
!  Given a series of series of saturation equivalent potential 
!  temperatures, the temperature is calculated.
!
!--------------------------------------------------------------------
   IMPLICIT NONE
!  SAVE !TWG 2017 add to avoid memory issues
!--------------------------------------------------------------------
! Lookup table variables
!     INTEGER, SAVE, PARAMETER :: KFNT=250,KFNP=220
!     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
!     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
!     REAL, SAVE, DIMENSION(1:200) :: ALU
!     REAL, SAVE :: RDPR,RDTHK,PLUTOP
! End of Lookup table variables

     INTEGER :: KP,IT,ITCNT,I
     REAL :: DTH,TMIN,TOLER,PBOT,DPR,                               &
             TEMP,P,ES,QS,PI,THES,TGUES,THGUES,F0,T1,T0,THGS,F1,DT, &
             ASTRT,AINC,A1,THTGS
!    REAL    :: ALIQ,BLIQ,CLIQ,DLIQ,SVP1,SVP2,SVP3,SVPT0
     REAL    :: ALIQ,BLIQ,CLIQ,DLIQ
     REAL, INTENT(IN)    :: SVP1,SVP2,SVP3,SVPT0
!
! equivalent potential temperature increment
      data dth/1./
! minimum starting temp 
      data tmin/150./
! tolerance for accuracy of temperature 
      data toler/0.001/
! top pressure (pascals)
      plutop=5000.0
! bottom pressure (pascals)
      pbot=110000.0

      ALIQ = SVP1*1000.
      BLIQ = SVP2
      CLIQ = SVP2*SVPT0
      DLIQ = SVP3

!
! compute parameters
!
! 1._over_(sat. equiv. theta increment)
      rdthk=1./dth
! pressure increment
!
      DPR=(PBOT-PLUTOP)/REAL(KFNP-1)
!      dpr=(pbot-plutop)/REAL(kfnp-1)
! 1._over_(pressure increment)
      rdpr=1./dpr
! compute the spread of thes
!     thespd=dth*(kfnt-1)
!
! calculate the starting sat. equiv. theta
!
      temp=tmin 
      p=plutop-dpr
      do kp=1,kfnp
        p=p+dpr
        es=aliq*exp((bliq*temp-cliq)/(temp-dliq))
        qs=0.622*es/(p-es)
        pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
        the0k(kp)=temp*pi*exp((3374.6525/temp-2.5403)*qs*        &
               (1.+0.81*qs))
      enddo   
!
! compute temperatures for each sat. equiv. potential temp.
!
      p=plutop-dpr
      do kp=1,kfnp
        thes=the0k(kp)-dth
        p=p+dpr
        do it=1,kfnt
! define sat. equiv. pot. temp.
          thes=thes+dth
! iterate to find temperature
! find initial guess
          if(it.eq.1) then
            tgues=tmin
          else
            tgues=ttab(it-1,kp)
          endif
          es=aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
          qs=0.622*es/(p-es)
          pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
          thgues=tgues*pi*exp((3374.6525/tgues-2.5403)*qs*      &
               (1.+0.81*qs))
          f0=thgues-thes
          t1=tgues-0.5*f0
          t0=tgues
          itcnt=0
! iteration loop
          do itcnt=1,11
            es=aliq*exp((bliq*t1-cliq)/(t1-dliq))
            qs=0.622*es/(p-es)
            pi=(1.e5/p)**(0.2854*(1.-0.28*qs))
            thtgs=t1*pi*exp((3374.6525/t1-2.5403)*qs*(1.+0.81*qs))
            f1=thtgs-thes
            if(abs(f1).lt.toler)then
              exit
            endif
!           itcnt=itcnt+1
            dt=f1*(t1-t0)/(f1-f0)
            t0=t1
            f0=f1
            t1=t1-dt
          enddo 
          ttab(it,kp)=t1 
          qstab(it,kp)=qs
        enddo
      enddo   
!
! lookup table for tlog(emix/aliq)
!
! set up intial values for lookup tables
!
       astrt=1.e-3
       ainc=0.075
!
       a1=astrt-ainc
       do i=1,200
         a1=a1+ainc
         alu(i)=alog(a1)
       enddo   
!
   END SUBROUTINE MSKF_LUTAB

!JTR 06/18/2019: Inserted momentum transport subroutine
    SUBROUTINE MSKF_CMT(DUDT,DVDT,dpdx,dpdy, &
    DP,ED,EU,MC,MD,MU,MB, &
    P,QD,QU,QHAT,SD,SU,SHAT,T,U,V,Z,ZF, &
    DSUBCLD,JB,JD,JT, &
    MSG,DT,GRAV,CPRES,RGAS,ILEV,IL1G,IL2G,ILG)

!     * JULY 17/92. - GUANG JUN ZHANG, M.LAZARE.

!     * PERFORMS MOMENTUM MIXING DUE TO CUMULUS PARAMETRIZATION.

! ongxl 20060901---------------------
    use shr_kind_mod, only: r8=>shr_kind_r8
    implicit none
    integer, parameter :: NBF = 10
! ongxl      PARAMETER(NBF=10)
! ongxl      PARAMETER(JLG=128, JLEV=18)
! ongxl 20060901---------------------
         
    real(r8) DUDT(ILG,ILEV),  DVDT(ILG,ILEV), &
    dpdx(ILG,ILEV),  dpdy(ILG,ILEV), &
    ALPHA(ILG,ILEV), DP(ILG,ILEV), &
    ED(ILG,ILEV),    EU(ILG,ILEV),    MB(ILG), &
    MC(ILG,ILEV),    MD(ILG,ILEV),    MU(ILG,ILEV), &
    P(ILG,ILEV),     QD(ILG,ILEV),    QU(ILG,ILEV), &
    QHAT(ILG,ILEV),  SD(ILG,ILEV),    SU(ILG,ILEV), &
    SHAT(ILG,ILEV),  T(ILG,ILEV),     U(ILG,ILEV), &
    V(ILG,ILEV),     Z(ILG,ILEV),     ZF(ILG,ILEV+1)    !songxl
! ongxl     7     V(ILG,ILEV),     Z(ILG,ILEV),     ZF(ILG,ILEV)

    real(r8) DSUBCLD(ILG)
! ongxl 20060901-----------------
    real(r8) DT, GRAV, CPRES, RGAS
    integer :: MSG, ILEV, IL1G, IL2G, ILG
! ongxl 20060901-----------------

    INTEGER ::   JB(ILG),    JD(ILG),         JT(ILG)

!     * INTERNAL WORK FIELDS.

    real(r8) AC(ILG,ILEV),    AD(ILG,ILEV),    AU(ILG,ILEV), &
    ACFL(ILG,ILEV),  ADFL(ILG,ILEV),  AUFL(ILG,ILEV), &
    B1(ILG,ILEV),    B1FL(ILG,ILEV),  BD(ILG,ILEV), &
    BU(ILG,ILEV),    D0(ILG,ILEV),    D0HAT(ILG,ILEV), &
    DELPX(ILG,ILEV), DELPY(ILG,ILEV), DZ(ILG,ILEV), &
    DZF(ILG,ILEV),   EC(ILG,ILEV),    E(ILG,ILEV), &
    EHAT(ILG,ILEV),  RHO(ILG,ILEV),   RHOHAT(ILG,ILEV), &
    UHAT(ILG,ILEV),  VHAT(ILG,ILEV),  UC(ILG,ILEV), &
    VC(ILG,ILEV),    W0(ILG,ILEV),    W1(ILG,ILEV), &
    W0FL(ILG,ILEV),  W1FL(ILG,ILEV)

    real(r8) COSA(ILG),    SINA(ILG),    CCX(ILG),     CCY(ILG), &
    UMN(ILG),     VMN(ILG),     DEP(ILG)

    real(r8) FX(ILG,ILEV,NBF), DX(ILG,ILEV,NBF), &
    AA(ILG,ILEV,NBF), BB(ILG,ILEV,NBF),  CC(ILG,ILEV,NBF), &
    FY(ILG,ILEV,NBF), DY(ILG,ILEV,NBF)

! ongxl      COMMON/CONST/ALFA1(NBF),ALFA2(NBF),BSJ0(NBF),BSJ1(NBF),
! ongxl     1             BSJ0CHI,FACTOR
! ongxl      COMMON/TAU/TAU(NBF)
! ongxl 20060901-----------------
    real(r8) ALFA1(NBF),ALFA2(NBF),BSJ0(NBF),BSJ1(NBF),BSJ0CHI,FACTOR,TAU(NBF)
! ongxl 20060901-----------------
    DATA ALFA1/4.33675e-05_r8, 7.88429e-05_r8, 1.13065e-04_r8, 1.45836e-04_r8, &
    &            1.76823e-04_r8, 2.05660e-04_r8, 2.32101e-04_r8, 2.55829e-04_r8, &
    &            2.76672e-04_r8, 2.94347e-04_r8/
    DATA ALFA2/-8.64731e-05_r8,-1.56092e-04_r8,-2.21315e-04_r8,-2.80995e-04_r8, &
    -3.33799e-04_r8,-3.78455e-04_r8,-4.14020e-04_r8,-4.39600e-04_r8, &
    -4.54652e-04_r8,-4.58761e-04_r8/
    DATA BSJ0/-4.02759e-01_r8, 3.00116e-01_r8, -2.49705e-01_r8, 2.18359e-01_r8, &
    -1.96465e-01_r8, 1.80062e-01_r8, -1.67183e-01_r8, 1.56722e-01_r8, &
    -1.48011e-01_r8, 1.40605e-01_r8/
    DATA BSJ1/1.14193e-01_r8, 2.05841e-01_r8, 2.91209e-01_r8, 3.68619e-01_r8, &
    &           4.36201e-01_r8, 4.92233e-01_r8, 5.35486e-01_r8, 5.64886e-01_r8, &
    &           5.79879e-01_r8, 5.80195e-01_r8/
    DATA BSJ0CHI/-0.402759_r8/
    DATA FACTOR/0.179503_r8/
    DATA TAU/ 3.83170e+00_r8, 7.01560e+00_r8, 1.01735e+01_r8, 1.33237e+01_r8, &
    &           1.64705e+01_r8, 1.96120e+01_r8, 2.27560e+01_r8, 2.58980e+01_r8, &
    &           2.90480e+01_r8, 3.21926e+01_r8/
! ongxl 20060901------------------------
    integer :: J, IL, N ,jj
    real(r8) RC, RMAX, CHI, SHAPE, RNU, WINDMAG
! ongxl 20060901-----------------------

!----------------------------------------------------------------------
!***********************************************************************
! CCC INITIALIZE RELEVANT INTERFACIAL AND MIDLAYER VARIABLES     CCCCC
!***********************************************************************
    RC=3000._r8         ! cloud radius (m)
    RMAX=50000._r8      ! distance where perturb. vanishes (m)
    CHI=3.8317_r8/RC
    SHAPE=0.8_r8
    DO 5 J=MSG+1,ILEV
        DO 5 IL=IL1G,IL2G
            BU(IL,J)=0._r8
            ALPHA(IL,J)=0.5_r8
            BD(IL,J)=0._r8
            DELPX(IL,J)=0._r8
            DELPY(IL,J)=0._r8
            EC(IL,J)=EU(IL,J)+ED(IL,J)     ! unit: 1/s
            IF(T(IL,J) > 0._r8)        THEN
                RHO(IL,J)=100._r8*P(IL,J)/(RGAS*T(IL,J))
            ENDIF
    5 END DO

    DO 10 N=1,NBF
        DO 10 J=MSG+1,ILEV
            DO 10 IL=IL1G,IL2G
                FX(IL,J,N)=0._r8
                FY(IL,J,N)=0._r8
    10 END DO

    DO 15 J=MSG+2,ILEV
        DO 15 IL=IL1G,IL2G
            BU(IL,J)=(   SU(IL,J)-SHAT(IL,J)+0.608_r8*( QU(IL,J) &
            *(SU(IL,J)-GRAV/CPRES*ZF(IL,J))-QHAT(IL,J) &
            *(SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) )  ) &
            /( (SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) &
            *(1._r8+.608_r8*QHAT(IL,J)) )*GRAV
            BD(IL,J)=(   SD(IL,J)-SHAT(IL,J)+0.608_r8*( QD(IL,J) &
            *(SD(IL,J)-GRAV/CPRES*ZF(IL,J))-QHAT(IL,J) &
            *(SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) )  ) &
            /( (SHAT(IL,J)-GRAV/CPRES*ZF(IL,J)) &
            *(1._r8+.608_r8*QHAT(IL,J)) )*GRAV
    15 END DO

    DO 20 IL=IL1G,IL2G
        UMN(IL)=0._r8
        VMN(IL)=0._r8
    !       DEP(IL)=0.
        COSA(IL)=0._r8
        SINA(IL)=0._r8
        CCX(IL)=0._r8
        CCY(IL)=0._r8
    20 END DO

    DO 25 J=1,ILEV
        DO 25 IL=IL1G,IL2G
            IF(J >= JT(IL) .AND. J < JB(IL))     THEN
                UMN(IL)=UMN(IL)+U(IL,J)*P(IL,J)/T(IL,J)* &
                (ZF(IL,J)-ZF(IL,J+1))
                VMN(IL)=VMN(IL)+V(IL,J)*P(IL,J)/T(IL,J)* &
                (ZF(IL,J)-ZF(IL,J+1))
            !         DEP(IL)=DEP(IL)+P(IL,J)/T(IL,J)*(ZF(IL,J)-ZF(IL,J+1))
            ENDIF
    25 END DO

    DO 30 IL=IL1G,IL2G
        WINDMAG=SQRT(UMN(IL)**2+VMN(IL)**2)
    !       IF(DEP(IL).NE.0. .AND. WINDMAG.NE.0.)   THEN
        IF(WINDMAG /= 0._r8)                       THEN
            COSA(IL)=UMN(IL)/WINDMAG
            SINA(IL)=VMN(IL)/WINDMAG
        !         CCX(IL)=0.*UMN(IL)/DEP(IL)
        !         CCY(IL)=0.*VMN(IL)/DEP(IL)
        ENDIF
    30 END DO

    DO 35 J=MSG+1,ILEV
        DO 35 IL=IL1G,IL2G
            IF(J > MSG+1)                          THEN
                RHOHAT(IL,J)=ALPHA(IL,J)*RHO(IL,J-1)+(1._r8-ALPHA(IL,J))*RHO(IL,J)
                DZF(IL,J)=Z(IL,J-1)-Z(IL,J)
                UHAT(IL,J)=ALPHA(IL,J)*U(IL,J-1)+(1._r8-ALPHA(IL,J))*U(IL,J)
                VHAT(IL,J)=ALPHA(IL,J)*V(IL,J-1)+(1._r8-ALPHA(IL,J))*V(IL,J)
            ELSE
                RHOHAT(IL,J)=RHO(IL,J)
                DZF(IL,J)=ZF(IL,J)-Z(IL,J)
                UHAT(IL,J)=U(IL,J)
                VHAT(IL,J)=V(IL,J)
            ENDIF
            UC(IL,J)=UHAT(IL,J)
            VC(IL,J)=VHAT(IL,J)
    35 END DO

    DO 40 J=MSG+1,ILEV
        DO 40 IL=IL1G,IL2G
            IF(J < ILEV)                       THEN
                DZ(IL,J)=ZF(IL,J)-ZF(IL,J+1)
            ELSE
            !          DZ(IL,ILEV)=ZF(IL,ILEV)
                DZ(IL,ILEV)=DP(IL,ILEV)*100._r8/(RHO(IL,ILEV)*GRAV) !m
            ENDIF
    40 END DO
!     ******************************************************************
!     AU, AD,AC ARE ACTUAL FRACTIONAL CLOUD AREA TIMES GRAV
!     ******************************************************************
    DO 75 J=MSG+1,ILEV
        DO 75 IL=IL1G,IL2G
            IF(J >= JT(IL))                     THEN
                AU(IL,J)=MU(IL,JB(IL))*100._r8/GRAV/RHOHAT(IL,JB(IL))
            ELSE
                AU(IL,J)=0._r8
            ENDIF
            IF(J >= JD(IL))                     THEN
                AD(IL,J)=-MD(IL,JD(IL))*100._r8/GRAV/RHOHAT(IL,JD(IL))
            ELSE
                AD(IL,J)=0._r8
            ENDIF
            AC(IL,J)=AU(IL,J)+AD(IL,J)
            IF(AC(IL,J) > 0._r8)                  THEN
                W0(IL,J)=MC(IL,J)*100._r8/GRAV/(RHOHAT(IL,J)*AC(IL,J))
            ELSE
                W0(IL,J)=0._r8
            ENDIF
    75 END DO

    DO 80 J=MSG+1,ILEV
        DO 80 IL=IL1G,IL2G
            IF(J < ILEV)                       THEN
                ACFL(IL,J)=ALPHA(IL,J)*AC(IL,J)+(1._r8-ALPHA(IL,J))*AC(IL,J+1)
                AUFL(IL,J)=ALPHA(IL,J)*AU(IL,J)+(1._r8-ALPHA(IL,J))*AU(IL,J+1)
                ADFL(IL,J)=ALPHA(IL,J)*AD(IL,J)+(1._r8-ALPHA(IL,J))*AD(IL,J+1)
                W0FL(IL,J)=ALPHA(IL,J)*W0(IL,J)+(1._r8-ALPHA(IL,J))*W0(IL,J+1)
            ELSE
                ACFL(IL,J)=ALPHA(IL,J)*AC(IL,J)
                AUFL(IL,J)=ALPHA(IL,J)*AU(IL,J)
                ADFL(IL,J)=ALPHA(IL,J)*AD(IL,J)
                W0FL(IL,J)=ALPHA(IL,J)*W0(IL,J)
            ENDIF
    80 END DO

    DO 250 J=MSG+1,ILEV
        DO 250 IL=IL1G,IL2G
            IF(J < JB(IL) .AND. J >= JT(IL) .AND. ACFL(IL,J) > 0._r8) THEN
                D0(IL,J)=( MC(IL,J+1)-MC(IL,J) )*100._r8/GRAV &
                /(RHO(IL,J)*DZ(IL,J)*ACFL(IL,J))
            ELSE
                IF(J == JB(IL) .AND. ACFL(IL,J) > 0._r8 )               THEN
                    D0(IL,J)= -MC(IL,J)*100._r8/GRAV &
                    /(RHO(IL,J)*DZ(IL,J)*ACFL(IL,J))
                ELSE
                    D0(IL,J)=0._r8
                ENDIF
            ENDIF
            IF(J < JB(IL) .AND. J >= JT(IL) .AND. ACFL(IL,J) > 0._r8 &
             .AND. ADFL(IL,J) > 0._r8 .AND. AUFL(IL,J) > 0._r8)        THEN
                E(IL,J)=1._r8/CHI**2*ADFL(IL,J)/(ACFL(IL,J)*RHO(IL,J)) &
                *( 1._r8/(AUFL(IL,J)*DZ(IL,J)) &
                *(MU(IL,J)-MU(IL,J+1))*100._r8/GRAV - &
                &                 1._r8/(ADFL(IL,J)*DZ(IL,J)) &
                *(MD(IL,J)-MD(IL,J+1))*100._r8/GRAV ) &
                *SHAPE/FACTOR

            ELSE
                E(IL,J)=0._r8
            ENDIF
    250 END DO

    DO 275 J=MSG+1,ILEV
        DO 275 IL=IL1G,IL2G
            IF(J > MSG+1)                                           THEN
                D0HAT(IL,J)=ALPHA(IL,J)*D0(IL,J-1)+(1._r8-ALPHA(IL,J))*D0(IL,J)
                EHAT(IL,J)=ALPHA(IL,J)*E(IL,J-1)+(1._r8-ALPHA(IL,J))*E(IL,J)
            ELSE
                D0HAT(IL,J)=(1._r8-ALPHA(IL,J))*D0(IL,J)
                EHAT(IL,J)=(1._r8-ALPHA(IL,J))*E(IL,J)
            ENDIF
    275 END DO

!     **********************************************
!     CALCULATE FIRST HARMONICS IN THERMODYNAMICS
!     **********************************************
    DO 325 J=MSG+1,ILEV
        DO 325 IL=IL1G,IL2G
            IF( J >= JD(IL) .AND. J <= JB(IL) .AND. AU(IL,J) > 0._r8 &
             .AND. AD(IL,J) > 0._r8 .AND. AC(IL,J) > 0._r8 )          THEN
                W1(IL,J)=SHAPE/FACTOR*(MU(IL,J)/AU(IL,J)-MD(IL,J)/AD(IL,J)) &
                *100._r8/GRAV*AD(IL,J)/(AC(IL,J)*RHOHAT(IL,J))
                B1(IL,J)=SHAPE/FACTOR*(BU(IL,J)-BD(IL,J))*AD(IL,J)/AC(IL,J)
            ELSE
                W1(IL,J)=0._r8
                B1(IL,J)=0._r8
            ENDIF
    325 END DO

    DO 350 J=MSG+1,ILEV
        DO 350 IL=IL1G,IL2G
            IF(J < ILEV)                                      THEN
                W1FL(IL,J)=ALPHA(IL,J)*W1(IL,J)+(1._r8-ALPHA(IL,J))*W1(IL,J+1)
                B1FL(IL,J)=ALPHA(IL,J)*B1(IL,J)+(1._r8-ALPHA(IL,J))*B1(IL,J+1)
            ELSE
                W1FL(IL,J)=ALPHA(IL,J)*W1(IL,J)
                B1FL(IL,J)=ALPHA(IL,J)*B1(IL,J)
            ENDIF
    350 END DO

    DO 500 N=1,NBF
        DO 500 J=MSG+1,ILEV
            DO 500 IL=IL1G,IL2G
                IF(J < JB(IL) .AND. J >= JT(IL) .AND. MC(IL,J) > 0._r8 &
                 .AND. MC(IL,J+1) > 0._r8)                          THEN
                    FX(IL,J,N)=2._r8*RHO(IL,J)/BSJ0(N)**2 &
                    *( D0(IL,J)*E(IL,J)*CHI**2*COSA(IL)*ALFA1(N) &
                    +2._r8*W0FL(IL,J)*RC/RMAX**2*BSJ1(N) &
                    *( (EHAT(IL,J)-EHAT(IL,J+1))/DZ(IL,J)*CHI*COSA(IL) &
                    *BSJ0CHI+(UHAT(IL,J)-UHAT(IL,J+1))/DZ(IL,J) ) &
                    -(D0HAT(IL,J)-D0HAT(IL,J+1))/DZ(IL,J)*W1FL(IL,J)*COSA(IL) &
                    *(ALFA2(N)-ALFA1(N))    -2._r8*( (W0(IL,J)-W0(IL,J+1)) &
                    *(W1(IL,J)-W1(IL,J+1))/DZ(IL,J)**2 &
                    -W0FL(IL,J)*W1FL(IL,J) &
                    *( LOG(RHO(IL,J-1))/(DZ(IL,J)*DZF(IL,J)) &
                    -(1._r8/DZF(IL,J)+1._r8/DZF(IL,J+1))/DZ(IL,J) &
                    *LOG(RHO(IL,J)) &
                    +LOG(RHO(IL,J+1))/(DZ(IL,J)*DZF(IL,J+1)) ) ) &
                    *COSA(IL)*ALFA1(N) &
                    +(RHOHAT(IL,J)*B1(IL,J)-RHOHAT(IL,J+1)*B1(IL,J+1)) &
                    /(DZ(IL,J)*RHO(IL,J))*COSA(IL)*ALFA1(N) )
                    FY(IL,J,N)=2._r8*RHO(IL,J)/BSJ0(N)**2 &
                    *( D0(IL,J)*E(IL,J)*CHI**2*SINA(IL)*ALFA1(N) &
                    +2._r8*W0FL(IL,J)*RC/RMAX**2*BSJ1(N) &
                    *( (EHAT(IL,J)-EHAT(IL,J+1))/DZ(IL,J)*CHI*SINA(IL) &
                    *BSJ0CHI+(VHAT(IL,J)-VHAT(IL,J+1))/DZ(IL,J) ) &
                    -(D0HAT(IL,J)-D0HAT(IL,J+1))/DZ(IL,J)*W1FL(IL,J)*SINA(IL) &
                    *(ALFA2(N)-ALFA1(N))    -2._r8*( (W0(IL,J)-W0(IL,J+1)) &
                    *(W1(IL,J)-W1(IL,J+1))/DZ(IL,J)**2 &
                    -W0FL(IL,J)*W1FL(IL,J) &
                    *( LOG(RHO(IL,J-1))/(DZ(IL,J)*DZF(IL,J)) &
                    -(1._r8/DZF(IL,J)+1._r8/DZF(IL,J+1))/DZ(IL,J) &
                    *LOG(RHO(IL,J)) &
                    +LOG(RHO(IL,J+1))/(DZ(IL,J)*DZF(IL,J+1)) ) ) &
                    *SINA(IL)*ALFA1(N) &
                    +(RHOHAT(IL,J)*B1(IL,J)-RHOHAT(IL,J+1)*B1(IL,J+1)) &
                    /(DZ(IL,J)*RHO(IL,J))*SINA(IL)*ALFA1(N) )
                ENDIF
    500 END DO

    DO 525 N=1,NBF
        DO 525 IL=IL1G,IL2G
            AA(IL,MSG+1,N)=0._r8
            BB(IL,MSG+1,N)=1._r8/ DZF(IL,MSG+2)
            CC(IL,MSG+1,N)=-1._r8/ DZF(IL,MSG+2)
              
            AA(IL,ILEV,N)=1._r8/ DZF(IL,ILEV)
            BB(IL,ILEV,N)=-1._r8/ DZF(IL,ILEV)
            CC(IL,ILEV,N)=0._r8
              
            DX(IL,MSG+1,N)=B1(IL,MSG+2)*ALFA1(N)*2._r8 &
            *COSA(IL)/BSJ0(N)**2
            DX(IL,ILEV,N)=B1(IL,ILEV)*ALFA1(N)*2._r8 &
            *COSA(IL)/BSJ0(N)**2
              
            DY(IL,MSG+1,N)=B1(IL,MSG+2)*ALFA1(N)*2._r8 &
            *SINA(IL)/BSJ0(N)**2
            DY(IL,ILEV,N)=B1(IL,ILEV)*ALFA1(N)*2._r8 &
            *SINA(IL)/BSJ0(N)**2
    525 END DO

    DO 550 N=1,NBF
        DO 550 J=MSG+2,ILEV-1
            DO 550 IL=IL1G,IL2G
                AA(IL,J,N)=1._r8/( DZ(IL,J)*DZF(IL,J) )
                BB(IL,J,N)=-( 1._r8/DZ(IL,J)*(1._r8/DZF(IL,J)+1._r8/DZF(IL,J+1)) &
                +(TAU(N)/RMAX)**2 )
                CC(IL,J,N)=1._r8/( DZ(IL,J)*DZF(IL,J+1) )
                  
                DX(IL,J,N)=FX(IL,J,N)
                DY(IL,J,N)=FY(IL,J,N)
    550 END DO


    DO 575 N=1,NBF
        DO 575 IL=IL1G,IL2G
            CC(IL,MSG+1,N)=CC(IL,MSG+1,N)/BB(IL,MSG+1,N)
            DX(IL,MSG+1,N)=DX(IL,MSG+1,N)/BB(IL,MSG+1,N)
            DY(IL,MSG+1,N)=DY(IL,MSG+1,N)/BB(IL,MSG+1,N)
    575 END DO

    DO 600 N=1,NBF
        DO 600 J=MSG+2,ILEV
            DO 600 IL=IL1G,IL2G
                CC(IL,J,N)=CC(IL,J,N)/(BB(IL,J,N)-AA(IL,J,N)*CC(IL,J-1,N))
                DX(IL,J,N)=(DX(IL,J,N)-AA(IL,J,N)*DX(IL,J-1,N)) &
                /(BB(IL,J,N)-AA(IL,J,N)*CC(IL,J-1,N))
                DY(IL,J,N)=(DY(IL,J,N)-AA(IL,J,N)*DY(IL,J-1,N)) &
                /(BB(IL,J,N)-AA(IL,J,N)*CC(IL,J-1,N))
    600 END DO

    DO 650 N=1,NBF
        DO 650 J=ILEV-1,MSG+1,-1
            DO 650 IL=IL1G,IL2G
                DX(IL,J,N)=DX(IL,J,N)-CC(IL,J,N)*DX(IL,J+1,N)
                DY(IL,J,N)=DY(IL,J,N)-CC(IL,J,N)*DY(IL,J+1,N)
    650 END DO

    DO 700 N=1,NBF
        DO 700 J=MSG+1,ILEV
            DO 700 IL=IL1G,IL2G
                DELPX(IL,J)=DELPX(IL,J)+DX(IL,J,N)*BSJ1(N)   !kg/m/s**2
                DELPY(IL,J)=DELPY(IL,J)+DY(IL,J,N)*BSJ1(N)   !kg/m/s**2
    700 END DO

    DO 850 J=MSG+1,ILEV
        DO 850 IL=IL1G,IL2G
            DELPX(IL,J)=DELPX(IL,J)/RC             !kg/m**2/s**2
            DELPY(IL,J)=DELPY(IL,J)/RC             !kg/m**2/s**2
        ! to get cloud-scale pressure gradient by multiplying cloud fraction
            DELPX(IL,J)=ACFL(IL,J)*DELPX(IL,J)/RHO(IL,J)  !m/s/s
            DELPY(IL,J)=ACFL(IL,J)*DELPY(IL,J)/RHO(IL,J)  !m/s/s
            dpdx(IL,J)=DELPX(IL,J)
            dpdy(IL,J)=DELPY(IL,J)
    850 END DO

!     ************************************
!     CALCULATE THE CLOUD MEAN WIND
!     ************************************

    DO 875 J=ILEV-1,MSG+1,-1
        DO 875 IL=IL1G,IL2G
            IF(MC(IL,J) > 0._r8 .AND. MC(IL,J+1) > 0._r8 &
             .AND. J > JT(IL) .AND. J < JB(IL))      THEN
                UC(IL,J)=UC(IL,J+1) + RHO(IL,J)*DZ(IL,J) &
                /((MC(IL,J)+MC(IL,J+1))*0.5_r8*100._r8/GRAV) &
                *( EC(IL,J)*(U(IL,J)-UC(IL,J+1))-DELPX(IL,J) )
                VC(IL,J)=VC(IL,J+1) + RHO(IL,J)*DZ(IL,J) &
                /((MC(IL,J)+MC(IL,J+1))*0.5_r8*100._r8/GRAV) &
                *( EC(IL,J)*(V(IL,J)-VC(IL,J+1))-DELPY(IL,J) )
            ENDIF
    875 END DO

!     RNU=1.
    RNU=0._r8
    DO 950 J=MSG+1,ILEV
        DO 950 IL=IL1G,IL2G
            IF( J >= JT(IL) .AND. J <= JB(IL) )                 THEN
                UHAT(IL,J)=UHAT(IL,J)+RNU*ALPHA(IL,J)*DT*DUDT(IL,J-1)
                VHAT(IL,J)=VHAT(IL,J)+RNU*ALPHA(IL,J)*DT*DVDT(IL,J-1)
                IF(J == JT(IL))                                     THEN
                    DUDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J))* &
                    MC(IL,J+1)*(UC(IL,J+1)-UHAT(IL,J+1))/DP(IL,J)
                    DVDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J))* &
                    MC(IL,J+1)*(VC(IL,J+1)-VHAT(IL,J+1))/DP(IL,J)
                ELSE IF(J < JB(IL))                                THEN
                    DUDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J)) &
                    *(MC(IL,J+1)*(UC(IL,J+1)-UHAT(IL,J+1)) &
                    -MC(IL,J)*(UC(IL,J)-UHAT(IL,J)))/DP(IL,J)
                    DVDT(IL,J)=1._r8/(1._r8+RNU*ALPHA(IL,J+1)*MC(IL,J+1)*DT/DP(IL,J)) &
                    *(MC(IL,J+1)*(VC(IL,J+1)-VHAT(IL,J+1)) &
                    -MC(IL,J)*(VC(IL,J)-VHAT(IL,J)))/DP(IL,J)
                ELSE
                    DUDT(IL,J)=-1._r8/DSUBCLD(IL)*MC(IL,J)*(UC(IL,J)-UHAT(IL,J))
                    DVDT(IL,J)=-1._r8/DSUBCLD(IL)*MC(IL,J)*(VC(IL,J)-VHAT(IL,J))
                ENDIF
            ENDIF
        ! dudt and dvdt (m/s/s)

            if (abs(dudt(il,j)) > 5.0e-2_r8 .OR. abs(dvdt(il,j)) > 5.0e-2_r8) then
                print*,'moment',il,j,dt,jt(il),jb(il) &
                ,dudt(il,j),dvdt(il,j),dp(il,j) &
                ,MC(IL,J+1),UC(IL,J+1),VC(IL,J+1),uhat(il,j+1),vhat(il,j+1) &
                ,DSUBCLD(IL),MC(IL,J),UC(IL,J),VC(IL,J),UHAT(IL,J),VHAT(IL,J)
                print*,'mb,msg,ilev=',mb(il),msg,ilev
                print*,'uc=',(uc(il,jj),jj=msg+1,ilev)
                print*,'vc=',(vc(il,jj),jj=msg+1,ilev)
                print*,'mc=',(mc(il,jj),jj=msg+1,ilev)
                print*,'mu=',(mu(il,jj),jj=msg+1,ilev)
                print*,'md=',(md(il,jj),jj=msg+1,ilev)
                print*,'u=',(u(il,jj),jj=msg+1,ilev)
                print*,'v=',(v(il,jj),jj=msg+1,ilev)
                print*,'ec=',(ec(il,jj),jj=msg+1,ilev)
                print*,'delpx=',(delpx(il,jj),jj=msg+1,ilev)
                print*,'delpy=',(delpy(il,jj),jj=msg+1,ilev)
                print*,'RHO=',(RHO(il,jj),jj=msg+1,ilev)
                print*,'dz=',(dz(il,jj),jj=msg+1,ilev)
            endif

        !     if (abs(dudt(il,j)).gt.1.0e-2.or.abs(dvdt(il,j)).gt.1.0e-2) then
        !      do i9=IL1G,IL2G
        !      do j9=MSG+1,ILEV
        !      print* ,'bad',i9,j9,dt,jt(i9),jb(i9)
        !    $       ,dudt(i9,j9),dvdt(i9,j9),dp(i9,j9),dz(i9,j9)
        !    $       ,mc(i9,j9+1),uc(i9,j9+1),vc(i9,j9+1),ec(i9,j9+1)
        !    $       ,du(i9,j9+1)
        !    $       ,uhat(i9,j9+1),vhat(i9,j9+1),delpx(i9,j9),delpx(i9,j9+1)
        !    $       ,delpy(i9,j9),delpy(i9,j9+1)
        !      enddo
        !      enddo
        !     endif
    950 END DO


    RETURN
    END SUBROUTINE MSKF_CMT

!=================================================================================================================
 end module module_cu_mskf
!=================================================================================================================
