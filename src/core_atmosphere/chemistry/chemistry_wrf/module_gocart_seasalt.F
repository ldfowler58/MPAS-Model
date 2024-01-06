!=================================================================================================================
 module module_gocart_seasalt
!=================================================================================================================
 use mpas_kind_types,only: RKIND,RKIND
 use mpas_log

 use module_data_gocart_seas


 implicit none
 private
 public:: gocart_seasalt_emission


 contains


!=================================================================================================================
 subroutine gocart_seasalt_emission(its,ite,kts,kte,ktep1,dt,grav,qseas1,qseas2,qseas3,qseas4,xland,area, &
                             dpres,zgrid,u10,v10,rho,u,v,emis_qseas1,emis_qseas2,emis_qseas3,emis_qseas4)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte,ktep1

 real(kind=RKIND),intent(in):: dt,grav
 real(kind=RKIND),intent(in),dimension(its:ite):: area,xland
 real(kind=RKIND),intent(in),dimension(its:ite):: u10,v10
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: dpres,rho,u,v
 real(kind=RKIND),intent(in),dimension(its:ite,kts:ktep1):: zgrid


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qseas1,qseas2,qseas3,qseas4


!output arguments:
 real(kind=RKIND),intent(out),dimension(its:ite),optional:: emis_qseas1,emis_qseas2, &
                                                            emis_qseas3,emis_qseas4


!local variables and arrays:
 integer,parameter:: &
    imx = 1, &
    lmx = 1, &
    nmx = 4
 integer:: i,k,ndt1
 integer,dimension(imx):: ilwi

 real(kind=RKIND):: conver,converi,dz,gravi
 real(kind=RKIND),dimension(imx):: w10m,gwet
 real(kind=RKIND),dimension(imx):: dxy
 real(kind=RKIND),dimension(imx,nmx):: bems
 real(kind=RKIND),dimension(imx,lmx):: airmas
 real(kind=RKIND),dimension(imx,lmx,nmx):: tc

 real(kind=RKIND):: tem1,tem2,tem3,tem4

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_seasalt_emission:')

 ndt1 = ifix(dt)
 conver  = 1.d-9
 converi = 1.d9
 gravi   = 1./grav

 if(present(emis_qseas1) .and. present(emis_qseas2) .and. present(emis_qseas3) .and. present(emis_qseas4)) then
    do i = its,ite
       emis_qseas1(i) = 0._RKIND
       emis_qseas2(i) = 0._RKIND
       emis_qseas3(i) = 0._RKIND
       emis_qseas4(i) = 0._RKIND
    enddo
 endif

 k = kts
 do i = its,ite

    if(xland(i).gt.1.5 .and. zgrid(i,k).lt.1.e-3)then
       ilwi(imx) = 0
       dxy(imx)  = area(i)
       dz = zgrid(i,k+1)-zgrid(i,k)
       airmas(imx,lmx) = rho(i,k)*dz*area(i)
!      airmas(imx,lmx) = gravi*dpres(i,k)*area(i)

       tc(imx,lmx,1) = qseas1(i,k)
       tc(imx,lmx,2) = qseas2(i,k)
       tc(imx,lmx,3) = qseas3(i,k)
       tc(imx,lmx,4) = qseas4(i,k)

       w10m(imx) = sqrt(u10(i)*u10(i)+v10(i)*v10(i))
       !--- we don't trust the u10,v10 values if model layers are very thin near surface:
       if(dz.lt.12.) w10m(imx) = sqrt(u(i,kts)*u(i,kts)+v(i,kts)*v(i,kts))

       call source_ss(imx,lmx,nmx,ndt1,tc,ilwi,dxy,w10m,airmas,bems)

       qseas1(i,k) = real(tc(imx,lmx,1),RKIND)
       qseas2(i,k) = real(tc(imx,lmx,2),RKIND)
       qseas3(i,k) = real(tc(imx,lmx,3),RKIND)
       qseas4(i,k) = real(tc(imx,lmx,4),RKIND)

       emis_qseas1(i) = bems(imx,1)/(dt*airmas(imx,lmx))
       emis_qseas2(i) = bems(imx,2)/(dt*airmas(imx,lmx))
       emis_qseas3(i) = bems(imx,3)/(dt*airmas(imx,lmx))
       emis_qseas4(i) = bems(imx,4)/(dt*airmas(imx,lmx))
    endif
 enddo

 call mpas_log_write('--- end subroutine gocart_seasalt_emission:')

 end subroutine gocart_seasalt_emission

!=================================================================================================================
 subroutine source_ss(imx,lmx,nmx,ndt1,tc,ilwi,dxy,w10m,airmas,bems)

! ***************************************************************************
! *  Evaluate the source of each seasalt particles size classes  (kg/m3) 
! *  by soil emission.
! *  Input:
! *         SSALTDEN  Sea salt density                              (kg/m3)
! *         DXY       Surface of each grid cell                     (m2)
! *         NDT1      Time step                                     (s)
! *         W10m      Velocity at the anemometer level (10meters)   (m/s)
! *      
! *  Output:
! *         DSRC      Source of each sea salt bins       (kg/timestep/cell) 
! *
! *
! * Number flux density: Original formula by Monahan et al. (1986) adapted
! * by Sunling Gong (JGR 1997 (old) and GBC 2003 (new)).  The new version is
! * to better represent emission of sub-micron sea salt particles.
!
! * dFn/dr = c1*u10**c2/(r**A) * (1+c3*r**c4)*10**(c5*exp(-B**2))
! * where B = (b1 -log(r))/b2
! * see c_old, c_new, b_old, b_new below for the constants.
! * number fluxes are at 80% RH.
! *
! * To calculate the flux:
! * 1) Calculate dFn based on Monahan et al. (1986) and Gong (2003)
! * 2) Assume that wet radius r at 80% RH = dry radius r_d *frh
! * 3) Convert particles flux to mass flux :
! *    dFM/dr_d = 4/3*pi*rho_d*r_d^3 *(dr/dr_d) * dFn/dr
! *             = 4/3*pi*rho_d*r_d^3 * frh * dFn/dr
! *               where rho_p is particle density [kg/m3]
! *    The factor 1.e-18 is to convert in micro-meter r_d^3
! ****************************************************************************
 
!=================================================================================================================

!input arguments:
 integer, intent(in):: imx,lmx,nmx,ndt1
 integer,intent(in):: ilwi(imx)

 real(kind=RKIND),intent(in),dimension(imx)    :: dxy,w10m
 real(kind=RKIND),intent(in),dimension(imx,lmx):: airmas

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(imx,lmx,nmx):: tc

!output arguments:
 real(kind=RKIND),intent(out),dimension(imx,nmx):: bems

!local variables and arrays:
 logical,parameter:: old=.TRUE.,new=.FALSE.
 logical,dimension(nmx):: aerosol

 character(len=20),dimension(nmx):: tcname,tcunits

 integer:: i,n,nr,ir

 real(kind=RKIND):: c0(5), b0(2)
!real(kind=RKIND),parameter:: c_old(5)=(/1.373, 3.41, 0.057, 1.05, 1.190/)
!real(kind=RKIND),parameter:: c_new(5)=(/1.373, 3.41, 0.057, 3.45, 1.607/)
!change suggested by mc:
 real(kind=RKIND),parameter:: c_old(5)=(/1.373, 3.2, 0.057, 1.05, 1.190/)
 real(kind=RKIND),parameter:: c_new(5)=(/1.373, 3.2, 0.057, 3.45, 1.607/)
 real(kind=RKIND),parameter:: b_old(2)=(/0.380, 0.650/)
 real(kind=RKIND),parameter:: b_new(2)=(/0.433, 0.433/)
 real(kind=RKIND),parameter:: dr=5.0d-2 ! um
 real(kind=RKIND),parameter:: theta=30.0
!swelling coefficient frh (d rwet / d rd)
!real(kind=RKIND),parameter:: frh = 1.65
 real(kind=RKIND),parameter:: frh = 2.d0

 real(kind=RKIND):: rho_d,r0,r1,r,r_w,a,b,dfn,r_d,dfm,src
 real(kind=RKIND),dimension(nmx):: tcmw,ar,tcvv
 real(kind=RKIND),dimension(nmx):: ar_wetdep,kc
 real(kind=RKIND),dimension(imx,lmx,nmx):: tc1

 real(kind=RKIND),dimension(imx,lmx,nmx),target:: tcms ! tracer mass (kg; kgS for sulfur case)
 real(kind=RKIND),dimension(imx,lmx,nmx),target:: tcgm ! g/m3

!-----------------------------------------------------------------------------------------------------------------  
! sea salt specific
!-----------------------------------------------------------------------------------------------------------------  
!real(kind=RKIND),dimension(nmx):: ra,rb
!real(kind=RKIND),dimension(nmx,12):: ch_ss

!-----------------------------------------------------------------------------------------------------------------  
! emissions (input)
!-----------------------------------------------------------------------------------------------------------------  
 real(kind=RKIND),dimension(imx,nmx):: e_bb
 real(kind=RKIND),dimension(imx,2,nmx):: e_an
 real(kind=RKIND),dimension(imx,lmx,nmx):: e_ac

!-----------------------------------------------------------------------------------------------------------------  
! diagnostics (budget)
!-----------------------------------------------------------------------------------------------------------------  
!tendencies per time step and process
!real(kind=RKIND),dimension(imx,jmx,nmx),target:: bems,bdry,bstl
!real(kind=RKIND),dimension(imx,jmx,nmx),target:: bwet,bcnv

!integrated tendencies per process
!real(kind=RKIND),dimension(imx,jmx,nmx),target:: tems,tstl
!real(kind=RKIND),dimension(imx,jmx,nmx),target:: tdry,twet,tcnv

!global mass balance per time step 
 real(kind=RKIND),dimension(nmx):: tmas0,tmas1
 real(kind=RKIND),dimension(nmx):: dtems,dttrp,dtdif,dtcnv
 real(kind=RKIND),dimension(nmx):: dtwet,dtdry,dtstl
 real(kind=RKIND),dimension(nmx):: dtems2,dttrp2,dtdif2,dtcnv2
 real(kind=RKIND),dimension(nmx):: dtwet2,dtdry2,dtstl2

!detailed integrated budgets for individual emissions
 real(kind=RKIND),dimension(imx),target:: ems_tp
 real(kind=RKIND),dimension(imx,nmx),target:: ems_an,ems_bb

 real(kind=RKIND),dimension(imx,nmx),target:: ems_co
 real(kind=RKIND),dimension(imx,lmx,nmx),target:: ems_ac

!-----------------------------------------------------------------------------------------------------------------  

 do n = 1,nmx
    bems(:,n) = 0.0
    rho_d = den_seas(n)
    r0 = ra_seas(n)*frh
    r1 = rb_seas(n)*frh
    r = r0
    nr = int((r1-r0)/dr+.001)

    do ir = 1,nr
       r_w = r + dr*0.5
       r = r + dr
       if(new) then
          a  = 4.7*(1.0 + theta*r_w)**(-0.017*r_w**(-1.44))
          c0 = c_new
          b0 = b_new
       else
          a  = 3.0
          c0 = c_old
          b0 = b_old
       endif

       b = (b0(1) - log10(r_w))/b0(2)
       dfn = (c0(1)/r_w**a)*(1.0 + c0(3)*r_w**c0(4))*10**(c0(5)*exp(-(b**2)))
       r_d = r_w/frh*1.0d-6  ! um -> m
       dfm = 4.0/3.0*pi*r_d**3*rho_d*frh*dfn*dr*ndt1

       do i = 1,imx
          if(ilwi(i) == 0) then
             src = dfm*dxy(i)*w10m(i)**c0(2)
             if(src < 0.0 ) src = 0.
             tc(i,1,n) = tc(i,1,n) + src/airmas(i,1)
          else
             src = 0.0
          endif
          bems(i,n) = bems(i,n) + src
       enddo ! j
    enddo ! ir
 enddo ! n

 end subroutine source_ss

!=================================================================================================================
 end module module_gocart_seasalt
!=================================================================================================================
