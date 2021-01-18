!-----------------------------------------------------------------------

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
                ,zrain,locum,prsfc,pssfc                        &
                ,pcte,pvom,pvol                                 &
                ,grav,rv,rd                                     &
                ,cu_act_flag,xland,slimsk                       &
                ,dt,stepcu,delt,rdelt                           &
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

      real(kind=kind_phys), dimension(its:ite, kts:kte) ::                          &
                                        pt,                            &
                                        pqv,                           &
                                        pap,                           &
                                        pqvf,                          &
                                        ptf,                           &
                                        pomg,                          &
                                        pgeo

      real(kind=kind_phys), dimension(its:ite, kts:kx1) ::             &
                                        pgeoh

      real(kind=kind_phys), dimension(its:ite, kts:kte) ::             &
                                        pq,                            &
                                        pqsat,                         &
                                        ztt,                           &
                                        zqq,                           &
                                        pcte,                          &
                                        pvom,                          &
                                        pvol

      real(kind=kind_phys), dimension(its:ite) ::                      &
                                        prsfc,                         &
                                        pssfc,                         &
                                        zrain

      real(kind=kind_phys), dimension(its:ite), intent(in) ::          &
                                        xland

      integer  , dimension(its:ite), intent(out) ::                    &
                                        slimsk

      logical  , dimension(its:ite) ::                                 &
                                        locum

      real(kind=kind_phys), intent(in) ::                              &
                                        dt

      real(kind=kind_phys), intent(out) ::                             &
                                        delt,                          &
                                        rdelt
      integer, intent(in) ::                                           &
                                        stepcu

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
                ,zprecc                                         &
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

      real(kind=kind_phys)     :: fliq,fice
      integer  :: i,k

      character(len=*),intent(out):: errmsg

      integer,intent(out)::errflg

!
!     compute cloud water and cloud ice detrained from convection
!
      do k=kts,kte
        do i=its,ite
        if(pcte(i,k).gt.0.) then
          fliq=foealfa(pt(i,k))
          fice=1.0-fliq
          pqc(i,k)=pqc(i,k)+fliq*pcte(i,k)*delt
          pqi(i,k)=pqi(i,k)+fice*pcte(i,k)*delt
        endif
        end do
      end do
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

      errmsg = 'cu_ntiedtke_mpas_timestep_final OK'
      errflg = 0

   end subroutine cu_ntiedtke_mpas_timestep_final

  subroutine cu_ntiedtke_mpas_init (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine cu_ntiedtke_mpas_init

  subroutine cu_ntiedtke_mpas_final (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

  end subroutine cu_ntiedtke_mpas_final

end module cu_ntiedtke_mpas
