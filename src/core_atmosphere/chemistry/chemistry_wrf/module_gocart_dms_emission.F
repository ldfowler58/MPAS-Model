!=================================================================================================================
 module module_gocart_dms_emission
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND
 use mpas_log

 use module_data_gocart_chem,only: ndms,smw,airmw


 implicit none
 private
 public:: gocart_dms_emission


 contains


!=================================================================================================================
 subroutine gocart_dms_emission(its,ite,kts,kte,dt,backg_dms,area,xland,tsk,u10,v10,rho,tt,u,v,zgrid,qdms)
!=================================================================================================================

!subroutine gocart_dmsemis(chem,dz8w,p8w,dms_0,ivgtyp,isltyp,xland,dx,g, &

!input arguments:
 integer,intent(in):: its,ite,kts,kte

 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in),dimension(its:ite):: backg_dms
 real(kind=RKIND),intent(in),dimension(its:ite):: area,xland,u10,v10,tsk
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho,tt,u,v,zgrid


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qdms


!local variables and arrays:
 integer,parameter:: &
    imx = 1, &
    lmx = 1, &
    nmx = 1
 integer:: i,k
 integer,dimension(imx):: ilwi

 real(kind=RKIND):: dz
 real(kind=RKIND),dimension(imx):: dxy,emsdms,dmso,tskin,w10m
 real(kind=RKIND),dimension(imx,nmx):: bems
 real(kind=RKIND),dimension(imx,lmx):: airmas
 real(kind=RKIND),dimension(imx,lmx,nmx):: tc

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_dms_emission:')

 k = kts
 do i = its,ite

    !--- no dms emission over land.
    if(xland(i) .gt. 1.5) then     
       ilwi(imx)  = 0                   
       dxy(imx)   = area(i)
       dmso(imx)  = backg_dms(i)
       tskin(imx) = min(tsk(i),40.)

       dz = zgrid(i,k+1)-zgrid(i,k)
       airmas(imx,lmx) = dz*area(i)*rho(i,k)

       w10m(imx)  = sqrt(u10(i)*u10(i) + v10(i)*v10(i))
       if(dz .lt. 12.) w10m(imx) = sqrt(u(i,k)*u(i,k) + v(i,k)*v(i,k))

!      tc(imx,lmx,nmx) = qdms(i,k) * 1.d-6
       !--- all species are in mol/mol, so no need to multiply tc by 1.e-6:
       tc(imx,lmx,ndms) = qdms(i,k)

       call srcdms(i,imx,lmx,nmx,ilwi,dt,dxy,dmso,airmas,tskin,w10m,emsdms,bems,tc)

       qdms(i,k) = tc(imx,lmx,ndms) * 1.e6
       !--- all species are in mol/mol, so no need to multiply tc by 1.e-6:
       qdms(i,k) = tc(imx,lmx,ndms)
    endif

 enddo

 call mpas_log_write('--- end subroutine gocart_dms_emissions:')

 end subroutine gocart_dms_emission

!=================================================================================================================
 subroutine srcdms(ii,imx,lmx,nmx,ilwi,dt1,area,dmso,airmas,tskin,w10m,emsdms,bems,tc)
!=================================================================================================================

! **************************************************************************
! **                                                                      **
! **  This subroutine calculates DMS emissions from the ocean.            **
! **                                                                      **
! **************************************************************************

!input arguments:
 integer,intent(in):: ii
 integer,intent(in):: imx,lmx,nmx
 integer,intent(in),dimension(imx):: ilwi

 real(kind=RKIND),intent(in):: dt1
 real(kind=RKIND),intent(in),dimension(imx):: area,airmas,dmso,tskin,w10m


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx):: emsdms
 real(kind=RKIND),intent(inout),dimension(imx,nmx):: bems
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc


!local variables and arrays:
 integer:: i,k

 real(kind=RKIND),parameter:: dms_mw = 62.00
 real(kind=RKIND),parameter:: tcmw(1)=dms_mw

 real(kind=RKIND):: akw,c,conc,dmssrc,erate,sc,scco2,sst,w10

!-----------------------------------------------------------------------------------------------------------------

! **************************************************************************
! *  ilwi = 0: water                                                      **
! *  If ilwi = 0: DMSEMS = seawaterDMS * transfer velocity.               **
! *  Otherwise,  DMSEMS = 0.0                                             **
! **************************************************************************

 k = 1

 do i = 1,imx

    !convert tskin (= sst over water) from K to degC:
    sst = tskin(i) - 273.15

    if(ilwi(i) == 0) then

       !Schmidt number for DMS (Saltzman et al., 1993)
       sc = 2674.0 - 147.12*sst + 3.726*(sst**2) - 0.038*(sst**3)

       ! ****************************************************************************
       ! *  Calculate transfer velocity in cm/hr  (AKw)                             *
       ! *                                                                          *
       ! *  Tans et al. transfer velocity (1990) for CO2 at 25oC (Erickson, 1993)   *
       ! *                                                                          *
       ! *  Tans et al. assumed AKW=0 when W10<=3. I modified it to let             *
       ! *  DMS emit at low windseeds too. Chose 3.6m/s as the threshold.           *
       ! *                                                                          *
       ! *  Schmidt number for CO2:       Sc = 600  (20oC, fresh water)             *
       ! *                                Sc = 660  (20oC, seawater)                *
       ! *                                Sc = 428  (25oC, Erickson 93)             *
       ! ****************************************************************************

       conc = dmso(i)
       w10  = w10m(i)

       !Liss and Merlivat (1986)
       scco2 = 600.0
       if(w10 <= 3.6) then
          akw = 0.17 * w10
       elseif(w10 <= 13.0) then
          akw = 2.85 * w10 - 9.65
       else
          akw = 5.90 * w10 - 49.3
       endif

       if(sc .le. 0._RKIND) call mpas_log_write('$i - sc = $r $r $r $r $r',intArgs=(/ii/), &
                                                realArgs=(/sc,sst,conc,w10,akw/))
       if(w10 <= 3.6) then
          akw = akw * ((scco2/sc) ** 0.667)
       else
          akw = akw * SQRT(scco2/sc)
       endif


       ! ****************************************************************************
       ! *  Calculate emission flux in kg/box/timestep                              *
       ! *                                                                          *
       ! *   AKw is in cm/hr:                 AKw/100/3600    -> m/sec              *
       ! *   CONC is in nmol/L (nmol/dm3):    CONC*1E-12*1000 -> kmol/m3            *
       ! *   TCMW(NDMS)       : kgDMS/kmol                                          *
       ! *   ERATE            : kgDMS/m2/timestep                                   *
       ! *   DMSSRC           : kgDMS/box/timestep                                  *
       ! ****************************************************************************

       erate  = akw/100.0/3600.0*conc*1.0e-12*1000.0*dt1*dms_mw
       dmssrc = erate * area(i)

    else   ! ilwi /= 0 (water)

       dmssrc = 0.0

    endif


    ! ****************************************************************************
    ! *  Update DMS concentration in level 1 (where emission occurs)             *
    ! ****************************************************************************

    !Convert emission from kg/box/timestep to mixing ratio/timestep:
    c = dmssrc / airmas(i) * airmw / dms_mw
    tc(i,k,ndms) = tc(i,k,ndms) + c


    !   ---------------------------------------------------------------
    !   Diagnostics:      DMS surface emission in kgS/timestep
    !   ---------------------------------------------------------------
    emsdms(i) = emsdms(i) + dmssrc * smw / dms_mw ! kgS
!   bems(i,ndms) = c * airmas(i) / airmw * smw    ! kgS
    bems(i,ndms) = dmssrc                         ! kgDMS

 enddo

 end subroutine srcdms

!=================================================================================================================
 end module module_gocart_dms_emission
!=================================================================================================================
