!=================================================================================================================
 module module_gocart_dust_data
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND
 use mpas_log


 implicit none
 save

 integer,parameter:: ndust=5,ndcls=3,ndsrc=1,maxstypes=100

 real(kind=RKIND),dimension(maxstypes):: porosity
 
 real(kind=RKIND),parameter:: ch_dust  = 1.0E-9 !empirical constant used in dust emission calculation. default is
                                                !1 ug s2 m-5 == 1.e-09 kg s2 m-5.

 real(kind=RKIND),parameter:: dyn_visc = 1.5E-5
!dust density (clay = 2500 kg m-3; Silt = 2650 kg m-3):
 real(kind=RKIND),dimension(ndust),parameter:: den_dust=(/2500.,2650.,2650.,2650.,2650./)

!dust effective radii (m):
 real(kind=RKIND),dimension(ndust),parameter:: reff_dust=(/0.73E-6,1.4E-6,2.4E-6,4.5E-6,8.0E-6/)

!USDA  defines dust particles with radius below 1 um as clay (class 3) and those with
!radius between 1 um and 25 um as silt (class 2). Array pointing to the class of particles
 integer,dimension(ndust),parameter:: ipoint=(/3,2,2,2,2/)

! fraction of each class available for erosion. fraction of clay is assumed to
! be 10% of the total mass of the emitted silt. 
 real(kind=RKIND),dimension(ndust),parameter:: frac_s=(/0.1,0.25,0.25,0.25,0.25/)

!=================================================================================================================
 end module module_gocart_dust_data
!=================================================================================================================


!=================================================================================================================
 module module_gocart_dust
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND
 use module_gocart_dust_data


 implicit none
 private
 public:: gocart_dust_emission


 contains


!=================================================================================================================
 subroutine gocart_dust_emission(its,ite,kts,kte,ktep1,nerod,dt,grav,qdust1,qdust2,qdust3,qdust4,qdust5, &
                      isltyp,xland,porosity,smois,area,erod,dpres,zgrid,rho,u,v,emis_qdust1,emis_qdust2, &
                      emis_qdust3,emis_qdust4,emis_qdust5)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte,ktep1
 integer,intent(in):: nerod
 integer,intent(in),dimension(its:ite):: isltyp

 real(kind=RKIND),intent(in):: dt,grav
 real(kind=RKIND),intent(in),dimension(its:ite):: area,porosity,smois,xland
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: dpres,zgrid
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho,u,v
 real(kind=RKIND),intent(in),dimension(its:ite,nerod):: erod


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qdust1,qdust2,qdust3,qdust4,qdust5


!output arguments:
 real(kind=RKIND),intent(out),dimension(its:ite),optional:: emis_qdust1,emis_qdust2,emis_qdust3, &
                                                            emis_qdust4,emis_qdust5


!local variables and arrays:
 integer,parameter:: &
    imx = 1,    &
    lmx = 1,    &
    nmx = ndust
 integer:: i,k,n
 integer,dimension(imx):: ilwi

 real(kind=RKIND):: dz,gravi
 real(kind=RKIND),dimension(imx):: dxy,gwet,w10m
 real(kind=RKIND),dimension(imx,nmx):: bems
 real(kind=RKIND),dimension(imx,lmx):: airden,airmas
 real(kind=RKIND),dimension(imx,lmx,nmx):: tc
 real(kind=RKIND),dimension(imx,ndcls,ndsrc):: erodin

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_dust_emission:')

 gravi = 1./grav

 k = kts
 do i = its,ite
    if(xland(i).lt.1.5)then
       ilwi(imx) = 1
       dxy(imx)  = area(i)
       dz = zgrid(i,k+1)-zgrid(i,k)
       w10m(imx) = sqrt(u(i,k)*u(i,k)+v(i,k)*v(i,k))
       gwet(imx) = smois(i)/porosity(i)
       airden(imx,lmx) = rho(i,k)
       airmas(imx,lmx) = rho(i,k)*dz*area(i)
!      airmas(imx,lmx) = gravi*dpres(i,k)*area(i)

       tc(imx,lmx,1) = qdust1(i,k)
       tc(imx,lmx,2) = qdust2(i,k)
       tc(imx,lmx,3) = qdust3(i,k)
       tc(imx,lmx,4) = qdust4(i,k)
       tc(imx,lmx,5) = qdust5(i,k)

       do n = 1,ndust
          bems(imx,n) = 0._RKIND
       enddo

       do n = 1,nerod
          erodin(imx,1,1) = erod(i,1)
          erodin(imx,2,1) = erod(i,2)
          erodin(imx,3,1) = erod(i,3)
       enddo

       call source_dust(imx,lmx,nmx,dt,tc,erodin,ilwi,dxy,w10m,gwet,airden,airmas,bems,grav)

       qdust1(i,k) = tc(imx,lmx,1)
       qdust2(i,k) = tc(imx,lmx,2)
       qdust3(i,k) = tc(imx,lmx,3)
       qdust4(i,k) = tc(imx,lmx,4)
       qdust5(i,k) = tc(imx,lmx,5)

       if(present(emis_qdust1) .and. present(emis_qdust2) .and. present(emis_qdust3) .and. &
          present(emis_qdust4) .and. present(emis_qdust5)) then
          emis_qdust1(i) = bems(imx,1)/(dt*airmas(imx,lmx))
          emis_qdust2(i) = bems(imx,2)/(dt*airmas(imx,lmx))
          emis_qdust3(i) = bems(imx,3)/(dt*airmas(imx,lmx))
          emis_qdust4(i) = bems(imx,4)/(dt*airmas(imx,lmx))
          emis_qdust5(i) = bems(imx,5)/(dt*airmas(imx,lmx))
       endif
    endif
 enddo

 call mpas_log_write('--- end subroutine gocart_dust_emission:')

 end subroutine gocart_dust_emission

!=================================================================================================================
 subroutine source_dust(imx,lmx,nmx,dt1,tc,erod,ilwi,dxy,w10m,gwet,airden,airmas,bems,g0)

! ****************************************************************************
! *  Evaluate the source of each dust particles size classes  (kg/m3)        
! *  by soil emission.
! *  Input:
! *         EROD      Fraction of erodible grid cell                (-)
! *                   for 1: Sand, 2: Silt, 3: Clay
! *         DUSTDEN   Dust density                                  (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         AIRVOL    Volume occupy by each grid boxes              (m3)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *         u_tresh   Threshold velocity for particule uplifting    (m/s)
! *         CH_dust   Constant to fudge the total emission of dust  (s2/m2)
! *      
! *  Output:
! *         DSRC      Source of each dust type           (kg/timestep/cell) 
! *
! *  Working:
! *         SRC       Potential source                   (kg/m/timestep/cell)
! *
! ****************************************************************************

!=================================================================================================================

!input arguments:
 integer,intent(in):: imx,lmx,nmx
 integer,intent(in),dimension(imx):: ilwi

 real(kind=RKIND),intent(in):: dt1,g0
 real(kind=RKIND),intent(in),dimension(imx):: dxy
 real(kind=RKIND),intent(in),dimension(imx):: gwet,w10m
 real(kind=RKIND),intent(in),dimension(imx,lmx):: airden,airmas
 real(kind=RKIND),intent(in),dimension(imx,ndcls,ndsrc):: erod

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc


!output arguments:
 real(kind=RKIND),intent(out),dimension(imx,nmx):: bems

!local variables and arrays:
 integer:: i,k,m,n

 real(kind=RKIND),parameter:: c0 = .2_RKIND
 real(kind=RKIND):: cw,dsrc,g,rhoa,srce,u_ts0,u_ts
 real(kind=RKIND),dimension(nmx):: diam,den

!-----------------------------------------------------------------------------------------------------------------

 do n = 1,nmx
    !Threshold velocity as a function of the dust density and the diameter
    !from Bagnold (1941):
    den(n)  = den_dust(n)*1.e-3_RKIND          ! conversion from kg/m3 to g/cm3.
    diam(n) = 2._RKIND*reff_dust(n)*1.e2_RKIND ! conversion from m to cm.
    g = g0*1.E2_RKIND

    !Pointer to the 3 classes considered in the source data files:
    m = ipoint(n)

    do k = 1,ndsrc
       !No flux if wet soil:
       do i = 1,imx
          rhoa = airden(i,1)*1.0e-3_RKIND

          u_ts0 = 0.13_RKIND*1.e-2_RKIND*sqrt(den(n)*g*diam(n)/rhoa)  &
                * sqrt(1._RKIND + .006_RKIND/den(n)/g/(diam(n))**2.5) &
                / sqrt(1.928_RKIND*(1331._RKIND*(diam(n))**1.56 + .38_RKIND)**.092 - 1._RKIND)

          !Fraction of emerged surfaces (subtract lakes, coastal ocean,..):
          !cw = 1.0 - water(i,j)

          !Case of surface dry enough to erode:
          if(gwet(i) < 0.5_RKIND) then
             u_ts = MAX(0._RKIND,u_ts0*(1.2_RKIND + c0*LOG10( MAX(1.0e-3_RKIND,gwet(i)))))
          else
          !Case of wet surface, no erosion:
             u_ts = 100.0_RKIND
          endif

          srce = frac_s(n)*erod(i,m,k)*dxy(i)  ! (m2)
          if(ilwi(i) == 1) then
             dsrc = ch_dust*srce*w10m(i)**2*(w10m(i)-u_ts)*dt1 ! where dsrc is in kg.
          else
             dsrc = 0._RKIND                                   ! where dsrc is in kg.
          endif
!         dsrc = cw*ch_dust(k)*srce*w10m(i)**2*(w10m(i)-u_ts)*dt1          ! where dsrc is in kg.
!         dsrc = cw*ch_dust(n,dt(1)%mn)*srce*w10m(i)**2*(w10m(i)-u_ts)*dt1 ! where dsrc is in kg.
          if(dsrc < 0._RKIND) dsrc = 0._RKIND

          !update dust mixing ratio at first model level:
          tc(i,1,n) = tc(i,1,n) + c0*dsrc/airmas(i,1)                      ! where tc is in kg/kg.
          bems(i,n) = c0*dsrc                                              ! where bems is in kg.
       enddo
    enddo
 enddo

 end subroutine source_dust

!=================================================================================================================
 end module module_gocart_dust
!=================================================================================================================
 
