!===================================================================================================

module module_cu_ntiedtke

  use mpas_kind_types,only: kind_phys => RKIND

  use mpas_atmphys_constants, only: g=>gravity, rd=>r_d, rv=>r_v

  use cu_ntiedtke_mpas, only: cu_ntiedtke_mpas_timestep_init, &
                              cu_ntiedtke_mpas_timestep_final
  use cu_ntiedtke, only: cu_ntiedtke_run

  implicit none
  private
  public :: ntiedtke

contains
!-----------------------------------------------------------------------
     subroutine ntiedtke(                                      &
                dt,itimestep,stepcu                            &
               ,raincv,pratec,qfx,hfx                          &
               ,u3d,v3d,w3d,t3d,qv3d,qc3d,qi3d,pi3d,rho3d      &
               ,qvften,thften                                  &
               ,dz8w,pcps,p8w,xland,cu_act_flag,dx             &
               ,ids,ide, jds,jde, kds,kde                      &
               ,ims,ime, jms,jme, kms,kme                      &
               ,its,ite, jts,jte, kts,kte                      &
               ,rthcuten,rqvcuten,rqccuten,rqicuten            &
               ,rucuten, rvcuten                               &
               ,f_qv    ,f_qc    ,f_qr    ,f_qi    ,f_qs       &
               ,errmsg,errflg                                  &
                                                               )
!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!-- rho3d       3d air density (kg/m^3)
!-- p8w         3d hydrostatic pressure at full levels (pa)
!-- pcps        3d hydrostatic pressure at half levels (pa)
!-- pi3d        3d exner function (dimensionless)
!-- qvften      3d total advective + PBL moisture tendency (kg kg-1 s-1)
!-- thften      3d total advective + PBL + radiative temperature tendency (k s-1)
!-- rthcuten    theta tendency due to 
!               cumulus scheme precipitation (k/s)
!-- rucuten     u wind tendency due to 
!               cumulus scheme precipitation (k/s)
!-- rvcuten     v wind tendency due to 
!               cumulus scheme precipitation (k/s)
!-- rqvcuten    qv tendency due to 
!               cumulus scheme precipitation (kg/kg/s)
!-- rqccuten    qc tendency due to 
!               cumulus scheme precipitation (kg/kg/s)
!-- rqicuten    qi tendency due to 
!               cumulus scheme precipitation (kg/kg/s)
!-- rainc       accumulated total cumulus scheme precipitation (mm)
!-- raincv      cumulus scheme precipitation (mm)
!-- pratec      precipitiation rate from cumulus scheme (mm/s)
!-- dz8w        dz between full levels (m)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- hfx         upward heat flux at the surface (w/m^2) 
!-- dt          time step (s)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------

      integer, intent(in) ::            ids,ide, jds,jde, kds,kde,      &
                                        ims,ime, jms,jme, kms,kme,      &
                                        its,ite, jts,jte, kts,kte,      &
                                        itimestep,                      &
                                        stepcu

      real(kind=kind_phys), intent(in) ::                               &
                                        dt
      real(kind=kind_phys), dimension(ims:ime, jms:jme), intent(in) ::  &
                                        dx

      real(kind=kind_phys), dimension(ims:ime, jms:jme), intent(in) ::  &
                                        xland

      real(kind=kind_phys), dimension(ims:ime, jms:jme), intent(inout) :: &
                                        raincv, pratec

      logical, dimension(ims:ime,jms:jme), intent(out) ::               &
                                        cu_act_flag

      real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(in) ::      &
                                        dz8w,                           &
                                        pcps,                           &
                                        p8w,                            &
                                        pi3d,                           &
                                        qc3d,                           &
                                        qvften,                         &
                                        thften,                         &
                                        qi3d,                           &
                                        qv3d,                           &
                                        rho3d,                          &
                                        t3d,                            &
                                        u3d,                            &
                                        v3d,                            &
                                        w3d

      real(kind=kind_phys), dimension(ims:ime, jms:jme), intent(in) ::  &
                                        qfx,                            &
                                        hfx                            

      character(len=*), intent(out) ::                                  &
                                        errmsg
      integer,          intent(out) ::                                  &
                                        errflg

      real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme),       &
               intent(out) ::                                           &
                                        rqccuten,                       &
                                        rqicuten,                       &
                                        rqvcuten,                       &
                                        rthcuten,                       &
                                        rucuten,                        &
                                        rvcuten
!
! flags relating to the tendency arrays declared above
! determine at run-time whether a particular tracer is in
! use or not.
!
     logical ::                                              &
                                                   f_qv      &
                                                  ,f_qc      &
                                                  ,f_qr      &
                                                  ,f_qi      &
                                                  ,f_qs

!--------------------------- local vars ------------------------------
      real(kind=kind_phys)      ::                                      &
                                        delt,                           &
                                        rdelt

      real(kind=kind_phys) , dimension(its:ite) ::                      &
                                        evap,                           &
                                        heatflux,                       &
                                        dx2d

      real(kind=kind_phys) , dimension(its:ite, kts:kte+1) ::           &
                                        dzi,                            &
                                         w

      real(kind=kind_phys) , dimension(its:ite, kts:kte) ::             &
                                        rho

      integer ::                                                        &
                                        i,                              &
                                        im,                             &
                                        j,                              &
                                        k,                              &
                                        km,                             &
                                        kp,                             &
                                        kx,                             &
                                        kx1

!-------other local variables----
      integer                      :: zz, pp
!
!------- new local variables ----
!
      real(kind=kind_phys) , dimension(its:ite, kts:kte) ::             &
                                        pap,                            &
                                        pqv,                            &
                                        pqc,                            &
                                        pqi,                            &
                                       pqvf,                            &
                                        ptf,                            &
                                         pt,                            &
                                         pu,                            &
                                         pv

      real(kind=kind_phys) , dimension(its:ite, kts:kte) ::             &
                                       pomg,                            &
                                       pgeo,                            &
                                         pq,                            &
                                      pqsat

      real(kind=kind_phys) , dimension(its:ite, kts:kte+1) ::           &
                                       paph,                            &
                                      pgeoh

      real(kind=kind_phys) , dimension(its:ite, kts:kte) ::             &
                                        ztt,                            &
                                        zqq

      real(kind=kind_phys) , dimension(its:ite) ::                      &
                                     zprecc
!
! more local arrays for cu which are not output out of the routine
!
      real(kind=kind_phys) , dimension(its:ite, kts:kte) ::             &
                                       pcte ,                           &
                                       pvom ,                           &
                                       pvol

      logical  , dimension(its:ite) ::                                  &
                                      cu_act_flag_local

      real(kind=kind_phys) , dimension(its:ite) ::                      &
                                      prsfc,                            &
                                      pssfc,                            &
                                      zrain,                            &
                                      xland_local

      real(kind=kind_phys) , dimension(its:ite, kts:kte) ::             &
                                         ztu ,                          &
                                         zqu ,                          &
                                         zlu ,                          &
                                       zlude ,                          &
                                        zmfu ,                          &
                                        zmfd

      integer  , dimension(its:ite) ::                                  &
                                       icbot,                           &
                                       ictop,                           &
                                       ktype,                           &
                                      slimsk

      logical  , dimension(its:ite) :: locum

      real(kind=kind_phys)     :: zew,zqs,zcor

      im=ite-its+1
      kx=kte-kts+1
      kx1=kx+1

!  begin first j loop

   do j=jts,jte

      do k=kts,kte+1
        do i=its,ite
          w(i,k)=w3d(i,k,j)
          dzi(i,k)=dz8w(i,k,j)
        enddo
      enddo
 
      do k=kts,kte
        do i=its,ite
          rho(i,k)=rho3d(i,k,j)
        enddo
      enddo

      do i=its,ite
        dx2d(i) = dx(i,j)
        xland_local(i)=xland(i,j)
      enddo

      pp = 0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          pu(i,zz)=u3d(i,k,j)
          pv(i,zz)=v3d(i,k,j)
          pt(i,zz)=t3d(i,k,j)
          pqv(i,zz)=qv3d(i,k,j)
          pqvf(i,zz)=qvften(i,k,j)
          ptf(i,zz)=thften(i,k,j)
          pqc(i,zz)=qc3d(i,k,j)
          pqi(i,zz)=qi3d(i,k,j)
          pap(i,zz)=pcps(i,k,j)
        enddo
        pp = pp + 1
      enddo

      pp = 0
      do k=kts,kte+1
        zz = kte+1-pp
        do i=its,ite
          paph(i,zz) = p8w(i,k,j) 
        enddo
        pp = pp + 1
      enddo

      do i=its,ite
        evap(i) = qfx(i,j)
        heatflux(i) = hfx(i,j)
      enddo

!  enddo ! first j loop

   call cu_ntiedtke_mpas_timestep_init(                         &
                 its,ite, kts,kte, kx1, itimestep               &
                ,w,pt,pqv,rho,dzi,pqvf,ptf,pap                  &
                ,ztt,zqq                                        &
                ,pomg,pq,pqsat,pgeo,pgeoh                       &
                ,zrain,locum,prsfc,pssfc                        &
                ,pcte,pvom,pvol                                 &
                ,g,rv,rd                                        &
                ,cu_act_flag_local,xland_local,slimsk           &
                ,dt,stepcu,delt,rdelt                           &
                ,errmsg,errflg                                  &
                                                                )
!
!     call 'cu_ntiedtke_run'(master-routine for cumulus parameterization)
!          (originally cumastrn)
!
   call cu_ntiedtke_run                                                          &
     &    (klon=im,       klev=kx,        klevp1=kx1,  klevm1=kx-1, pten=pt,     &
     &     pqen=pq,       puen=pu,        pven=pv,     pverv=pomg,  pqsen=pqsat, &
     &     pqhfl=evap,    ztmst=delt,     pap=pap,     paph=paph,   pgeo=pgeo,   &
     &     ptte=ptf,      pqte=pqvf,      pvom=pvom,   pvol=pvol,                &
     &     prsfc=prsfc,   pssfc=pssfc,    ldcum=locum,                           &
     &     ktype=ktype,   kcbot=icbot,    kctop=ictop, ptu=ztu,     pqu=zqu,     &
     &     plu=zlu,       plude=zlude,    pmfu=zmfu,   pmfd=zmfd,   prain=zrain, &
     &     pcte=pcte,     phhfl=heatflux, lndj=slimsk, zgeoh=pgeoh, dx=dx2d,     &  
     &     errmsg=errmsg, errflg=errflg                                          &
          )

   call cu_ntiedtke_mpas_timestep_final(                        &
                 its,ite, kts,kte                               &
                ,delt,stepcu,dt,rdelt                           &    
                ,pcte,ptf,pqvf,prsfc,pssfc                      &
                ,pvom,pvol,ztt,zqq                              &
                ,pqc,pqi,pt,pq,pqv,pu,pv                        &
                ,zprecc                                         &
                ,errmsg,errflg                                  &
                                                                )
!
!     compute tendencies for model variables
!
!  do j=jts,jte

      do i=its,ite
         raincv(i,j)=zprecc(i)/stepcu
         pratec(i,j)=zprecc(i)/(stepcu * dt)
      enddo

      pp = 0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          rthcuten(i,k,j)=(pt(i,zz)-t3d(i,k,j))/pi3d(i,k,j)*rdelt
          rqvcuten(i,k,j)=(pqv(i,zz)-qv3d(i,k,j))*rdelt
          rucuten(i,k,j) =(pu(i,zz)-u3d(i,k,j))*rdelt
          rvcuten(i,k,j) =(pv(i,zz)-v3d(i,k,j))*rdelt
        enddo
        pp = pp + 1
      enddo

        if ( f_qc ) then
          pp = 0
          do k=kts,kte
            zz = kte-pp
            do i=its,ite
              rqccuten(i,k,j)=(pqc(i,zz)-qc3d(i,k,j))*rdelt
            enddo
            pp = pp + 1
          enddo
        endif

        if ( f_qi ) then
          pp = 0
          do k=kts,kte
            zz = kte-pp
            do i=its,ite
              rqicuten(i,k,j)=(pqi(i,zz)-qi3d(i,k,j))*rdelt
            enddo
            pp = pp + 1
          enddo
        endif

      do i=its,ite
         cu_act_flag(i,j)=cu_act_flag_local(i)
      enddo

   enddo
!  end of second j loop

   end subroutine ntiedtke

!====================================================================

end module module_cu_ntiedtke
