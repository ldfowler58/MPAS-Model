

module cu_ntiedtke_mpas

USE ccpp_kinds, ONLY: kind_phys
 
  use cu_ntiedtke, only: foeewm, foealfa

  implicit none

  PUBLIC :: cu_ntiedtke_mpas_run
  PUBLIC :: cu_ntiedtke_mpas_init
  PUBLIC :: cu_ntiedtke_mpas_final
  PUBLIC :: cu_ntiedtke_mpas_timestep_init
  PUBLIC :: cu_ntiedtke_mpas_timestep_final

contains

!-----------------------------------------------------------------------
 subroutine cu_ntiedtke_mpas_run(errmsg,errflg)
!-----------------------------------------------------------------------

implicit none
!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out)::errflg

 errmsg = 'cu_ntiedtke_mpas_run OK'
 errflg = 0

 end subroutine cu_ntiedtke_mpas_run

!-----------------------------------------------------------------------
     subroutine cu_ntiedtke_mpas_timestep_init(                 &
                 its,ite, kts,kte, kx1, itimestep               &
                ,w,pt,pqv,rho,dzi,pqvf,ptf,pap                  &
                ,ztt,zqq                                        &
                ,pomg,pq,pqsat,pgeo,pgeoh                       &
                ,pu,pv,pqc,pqi,paph                             &
                ,zrain,locum,prsfc,pssfc                        &
                ,pcte,pvom,pvol                                 &
                ,grav,rv,rd                                     &
                ,cu_act_flag,xland,slimsk                       &
                ,dt,stepcu,delt,rdelt                           &
                ,f_qc,f_qi                                      &
                ,t3d_l,qv3d_l,qc3d_l,qi3d_l,u3d_l,v3d_l         &
                ,qvften_l,thften_l,p8w_l,pcps_l                 &
                ,errmsg,errflg                                  &
                                                                )
!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-- itimestep   time step count
!-- pt          2d temperature (k)
!-- pqv         2d water vapor mixing ratio (kg/kg)
!-- w           2d vertical velocity (m s-1) (bottom up)
!-- rho         2d air density (kg/m^3) (bottom up)
!-- dzi         2d dz between full levels (m) (bottom up)
!-- pqvf        2d total advective + PBL moisture tendency (kg kg-1 s-1)
!-- ptf         2d total advective + PBL + radiative temperature tendency (k s-1)
!-- zi          height at interface level (full level)
!-- zl          height at middle of the layer (half level)
!-- dot         vertical pressure velocity at layer center (bottom up)
!-- pomg        vertical pressure velocity at layer center (top down)
!-- zrain       surface rain field out of cu
!-- locum       logical for convective grid point out of cu
!-- prsfc       liquid surface rainfall out of cu
!-- pssfc       frozen surface rainfall out of cu
!-- pcte        total surface rainfall
!-- pvom        u component of wind out of cu
!-- pvol        v component of wind out of cu
!-- pq          specific humidity
!-- pqsat       saturation specific humidity
!-- pgeo        geopotential at layer center
!-- pgeoh       geopotential at interface level
!-------------------------------------------------------------------

      integer, intent(in) ::            its,ite, kts,kte, kx1,            &
                                        itimestep

      real(kind=kind_phys), intent(in) ::            grav,rv,rd

      real(kind=kind_phys), dimension(its:ite, kts:kte), intent(in) :: &
                                        rho

      real(kind=kind_phys), dimension(its:ite, kts:kx1), intent(in) :: &
                                        dzi,                           &
                                        w

      logical, dimension(its:ite), intent(out)         ::              &
                                        cu_act_flag

      real(kind=kind_phys), dimension(its:ite, kts:kte), intent(out) ::&
                                        pu,                            &
                                        pv,                            &
                                        pt,                            &
                                        pqv,                           &
                                        pqc,                           &
                                        pqi,                           &
                                        pap,                           &
                                        pqvf,                          &
                                        ptf,                           &
                                        pomg,                          &
                                        pgeo

      real(kind=kind_phys), dimension(its:ite, kts:kx1), intent(out) ::&
                                        pgeoh,                         &
                                        paph

      real(kind=kind_phys), dimension(its:ite, kts:kte), intent(out) ::&
                                        pq,                            &
                                        pqsat,                         &
                                        ztt,                           &
                                        zqq,                           &
                                        pcte,                          &
                                        pvom,                          &
                                        pvol

      real(kind=kind_phys), dimension(its:ite), intent(out) ::         &
                                        prsfc,                         &
                                        pssfc,                         &
                                        zrain

      real(kind=kind_phys), dimension(its:ite), intent(in) ::          &
                                        xland

      integer  , dimension(its:ite), intent(out) ::                    &
                                        slimsk

      logical  , dimension(its:ite), intent(out) ::                    &
                                        locum

      real(kind=kind_phys), intent(in) ::                              &
                                        dt

      real(kind=kind_phys), intent(out) ::                             &
                                        delt,                          &
                                        rdelt
      integer, intent(in) ::                                           &
                                        stepcu

      logical, intent(in) ::                                           &
                                        f_qc,                          &
                                        f_qi

      real(kind=kind_phys), dimension(its:ite, kts:kte), intent(in) :: &
                                        t3d_l,                         &
                                        qv3d_l,                        &
                                        qc3d_l,                        &
                                        qi3d_l,                        &
                                        u3d_l,                         &
                                        v3d_l,                         &
                                        qvften_l,                      &
                                        thften_l,                      &
                                        pcps_l

      real(kind=kind_phys), dimension(its:ite, kts:kx1), intent(in) :: &
                                        p8w_l

      character(len=*),intent(out) ::                                  &
                                        errmsg

      integer,intent(out) ::                                           &
                                        errflg

!--------------------------- local vars ------------------------------

      real(kind=kind_phys), dimension(its:ite, kts:kx1) ::             &
                                        zi

      real(kind=kind_phys), dimension(its:ite, kts:kte) ::             &
                                        dot,                           &
                                        zl
!
      real(kind=kind_phys)     :: zew,zqs,zcor,tt,vtmpc1
      integer  :: i,k,kp,pp,zz
!
!     initializations for output
!
      delt=dt*stepcu
      rdelt = 1./delt

      do i=its,ite
         cu_act_flag(i)=.true.
         slimsk(i)=int(abs(xland(i)-2.))
      enddo

      vtmpc1 = rv/rd-1.0

      call upside_down(p8w_l    ,paph,its,ite,kts,kx1)
      call upside_down(u3d_l    ,pu  ,its,ite,kts,kte)
      call upside_down(v3d_l    ,pv  ,its,ite,kts,kte)
      call upside_down(t3d_l    ,pt  ,its,ite,kts,kte)
      call upside_down(qv3d_l   ,pqv ,its,ite,kts,kte)
      call upside_down(qvften_l ,pqvf,its,ite,kts,kte)
      call upside_down(thften_l ,ptf ,its,ite,kts,kte)
      if ( f_qc ) then
         call upside_down(qc3d_l   ,pqc ,its,ite,kts,kte)
      end if
      if ( f_qi ) then
         call upside_down(qi3d_l   ,pqi ,its,ite,kts,kte)
      end if
      call upside_down(pcps_l   ,pap ,its,ite,kts,kte)

!     compute height at full (zi) and half (zl) levels

      do i=its,ite
        zi(i,kts)=0.0
      enddo
 
      do k=kts,kte
        do i=its,ite
          zi(i,k+1)=zi(i,k)+dzi(i,k)
        enddo
      enddo
 
      do k=kts,kte
        do i=its,ite
          zl(i,k)=0.5*(zi(i,k)+zi(i,k+1))
        enddo
      enddo

      do k=kts,kte
        kp=k+1
        do i=its,ite
          dot(i,k)=-0.5*grav*rho(i,k)*(w(i,k)+w(i,kp))
        enddo
      enddo

      pp = 0
      do k=kts,kte
        zz = kte-pp
        do i=its,ite
          pomg(i,zz)=dot(i,k)
          pgeo(i,zz)=zl(i,k)*grav
        enddo
        pp = pp + 1
      enddo

      pp = 0
      do k=kts,kx1
        zz = kx1-pp
        do i=its,ite
          pgeoh(i,zz) = zi(i,k)*grav
        enddo
        pp = pp + 1
      enddo

! DAVE this should be moved to an init routine
      if(itimestep == 1) then
        do k=kts,kte
          do i=its,ite
             ptf(i,k)=0.
             pqvf(i,k)=0.
          enddo
        enddo
      endif

! need to save the input values
      do k=kts,kte
        do i=its,ite
           ztt(i,k)=ptf(i,k)
           zqq(i,k)=pqvf(i,k)
        enddo
      enddo
!
!  initializing arrays
!
      do i=its,ite
        zrain(i)=0.0
        locum(i)=.false.
        prsfc(i)=0.0
        pssfc(i)=0.0
      end do

      do k=kts,kte
        do i=its,ite
          pcte(i,k)=0.0
          pvom(i,k)=0.0
          pvol(i,k)=0.0
        end do
      end do
!
!     convert model variables for cu scheme
!
      do k=kts,kte
        do i=its,ite
          pq(i,k)=pqv(i,k)/(1.0+pqv(i,k))
          tt=pt(i,k)
          zew  = foeewm(tt)
          zqs  = zew/pap(i,k)
          zqs  = min(0.5,zqs)
          zcor = 1./(1.-vtmpc1*zqs)
          pqsat(i,k)=zqs*zcor
        end do
      end do

      errmsg = 'cu_ntiedtke_mpas_timestep_init OK'
      errflg = 0

end subroutine cu_ntiedtke_mpas_timestep_init

!-----------------------------------------------------------------------
     subroutine cu_ntiedtke_mpas_timestep_final(                &
                 its,ite, kts,kte                               &
                ,delt,stepcu,dt,rdelt                           &
                ,pcte,ptf,pqvf,prsfc,pssfc                      &
                ,pvom,pvol,ztt,zqq                              &
                ,pqc,pqi,pt,pq,pqv,pu,pv                        &
                ,f_qc,f_qi                                      &
                ,raincv_l,pratec_l                              &
                ,t3d_l,qv3d_l,qc3d_l,qi3d_l,u3d_l,v3d_l,pi3d_l  &
                ,rthcuten_l,rqvcuten_l,rqccuten_l,rqicuten_l    &
                ,rucuten_l,rvcuten_l                            &
                ,errmsg,errflg                                  &
                                                                )
!-------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------
!-- pu          2d u-velocity interpolated to theta points (m/s)
!-- pv          2d v-velocity interpolated to theta points (m/s)
!-- pt          2d temperature (k)
!-- pq          2d specific humidity (kg/kg)
!-- pqv         2d water vapor mixing ratio (kg/kg)
!-- pqc         2d cloud mixing ratio (kg/kg)
!-- pqi         2d ice mixing ratio (kg/kg)
!-- pcte        2d hydrometeos 
!-- pqvf        2d total advective + PBL + cumulus moisture tendency (kg kg-1 s-1)
!-- ptf         2d total advective + PBL + radiative + cumulus temperature tendency (k s-1)
!-- pvom        2d u change due to cumulus
!-- pvol        2d v change due to cumulus
!-- ztt         2d total advective + PBL  moisture tendency before cumulus (kg kg-1 s-1)
!-- zqq         2d total advective + PBL  + radiation temperature tendenc before cumulus y (k s-1)
!-- zprecc      total cumulus scheme precipitation (mm)
!-- delt        1 over time step (dt)
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------

      integer, intent(in) ::            its,ite, kts,kte

      real(kind=kind_phys), intent(in) ::                             &
                                       delt,                          &
                                         dt,                          &
                                      rdelt

      integer ,             intent(in) ::                             &
                                     stepcu

      real(kind=kind_phys), dimension(its:ite, kts:kte) ::            &
                                        pqv,                          &
                                        pqc,                          &
                                        pqi,                          &
                                       pqvf,                          &
                                        ptf,                          &
                                         pt,                          &
                                         pq,                          &
                                         pu,                          &
                                         pv
 
      real(kind=kind_phys), dimension(its:ite, kts:kte) ::            &
                                        ztt,                          &
                                        zqq,                          &
                                       pcte,                          &
                                       pvom,                          &
                                       pvol

      real(kind=kind_phys), dimension(its:ite) ::                     &
                                     zprecc

      real(kind=kind_phys), dimension(its:ite) ::                     &
                                      prsfc,                          &
                                      pssfc

      logical, intent(in) ::                                          &
                                        f_qc,                         &
                                        f_qi

      real(kind=kind_phys), dimension(its:ite), intent(out) ::        &
                                   raincv_l,                          &
                                   pratec_l

      real(kind=kind_phys), dimension(its:ite, kts:kte), intent(in) ::&   
                                       t3d_l ,                        &   
                                      qv3d_l ,                        &   
                                      qc3d_l ,                        &   
                                      qi3d_l ,                        &   
                                       u3d_l ,                        &   
                                       v3d_l ,                        &   
                                      pi3d_l

      real(kind=kind_phys), dimension(its:ite, kts:kte), intent(out) :: &
                                  rthcuten_l ,                          &
                                  rqvcuten_l ,                          &
                                  rqccuten_l ,                          &
                                  rqicuten_l ,                          &
                                   rucuten_l ,                          &
                                   rvcuten_l

      real(kind=kind_phys)     :: fliq,fice

      integer  :: i,k,kk

      character(len=*),intent(out):: errmsg

      integer,intent(out)::errflg

!
!     compute cloud water and cloud ice detrained from convection
!
!DAVE, this should have qc qi IF test
      if ( f_qc ) then 
        do k=kts,kte
          do i=its,ite
          if(pcte(i,k).gt.0.) then
            fliq=foealfa(pt(i,k))
            pqc(i,k)=pqc(i,k)+fliq*pcte(i,k)*delt
          endif
          end do
        end do
      end if

      if ( f_qi ) then 
        do k=kts,kte
          do i=its,ite
          if(pcte(i,k).gt.0.) then
            fliq=foealfa(pt(i,k))
            fice=1.0-fliq
            pqi(i,k)=pqi(i,k)+fice*pcte(i,k)*delt
          endif
          end do
        end do
      end if
!
!     compute temperature and water vapor mixing ratio after convection
!
      do k=kts,kte
        do i=its,ite
          pt(i,k)=pt(i,k)+(ptf(i,k)-ztt(i,k))*delt
          pq(i,k)=pq(i,k)+(pqvf(i,k)-zqq(i,k))*delt
          pqv(i,k)=pq(i,k)/(1.0-pq(i,k))
        end do
      end do
!
!     compute surface rainfall
!
      do i=its,ite
        zprecc(i)=amax1(0.0,(prsfc(i)+pssfc(i))*delt)
      end do
!
!     compute u and v winds after convection
!
      do k=kts,kte
        do i=its,ite
          pu(i,k)=pu(i,k)+pvom(i,k)*delt
          pv(i,k)=pv(i,k)+pvol(i,k)*delt
        end do
      end do
!
!     precipitation rate
!
      do i=its,ite
         raincv_l(i)=zprecc(i)/stepcu
         pratec_l(i)=zprecc(i)/(stepcu * dt) 
      enddo
!
!     tendencies due to CU scheme
!
      do k=kts,kte
        kk=kte+1-k
        do i=its,ite
          rthcuten_l(i,k)=(pt(i,kk)-t3d_l(i,k))/pi3d_l(i,k)*rdelt
          rqvcuten_l(i,k)=(pqv(i,kk)-qv3d_l(i,k))*rdelt
          rucuten_l(i,k) =(pu(i,kk)-u3d_l(i,k))*rdelt
          rvcuten_l(i,k) =(pv(i,kk)-v3d_l(i,k))*rdelt
        end do
        if ( f_qc ) then
          do i=its,ite
            rqccuten_l(i,k)=(pqc(i,kk)-qc3d_l(i,k))*rdelt
          end do
        end if
        if ( f_qi ) then
          do i=its,ite
            rqicuten_l(i,k)=(pqi(i,kk)-qi3d_l(i,k))*rdelt
          end do
        end if
      end do

      errmsg = 'cu_ntiedtke_mpas_timestep_final OK'
      errflg = 0

   end subroutine cu_ntiedtke_mpas_timestep_final

  subroutine cu_ntiedtke_mpas_init (errmsg, errflg)

    implicit none
    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine cu_ntiedtke_mpas_init

  subroutine cu_ntiedtke_mpas_final (errmsg, errflg)
  
    implicit none
    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine cu_ntiedtke_mpas_final

  subroutine upside_down (in,out,is,ie,ks,ke)
    implicit none
    integer, intent(in) :: is, ie, ks, ke
    real(kind=kind_phys), intent(in) , dimension(is:ie,ks:ke) :: in
    real(kind=kind_phys), intent(out), dimension(is:ie,ks:ke) :: out
    integer :: i, k, kk
    do k = ks, ke
      do i = is, ie
        kk= ke+1-k
        out(i,kk) = in(i,k)      
      end do
    end do 
  end subroutine upside_down

end module cu_ntiedtke_mpas
