!=================================================================================================================
 module parkind
!implicit none
 save

!------------------------------------------------------------------
! rrtmg kinds
! Define integer and real kinds for various types.
!
! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!
!integer kinds
!-------------
!
!integer, parameter :: kind_ib = selected_int_kind(13)  ! 8 byte integer
!integer, parameter :: kind_im = selected_int_kind(6)   ! 4 byte integer
 integer, parameter :: kind_ib = kind(1)            
 integer, parameter :: kind_im = kind(1)            
 integer, parameter :: kind_in = kind(1)                ! native integer

!
!real kinds
!----------
!
!integer, parameter :: kind_rb = selected_real_kind(12) ! 8 byte real
!integer, parameter :: kind_rm = selected_real_kind(6)  ! 4 byte real
!integer, parameter :: kind_rn = kind(1.0)              ! native real

#ifdef SINGLE_PRECISION
#define RWORDSIZE 4
#define IWORDSIZE 4
#else
#define RWORDSIZE 8
#define IWORDSIZE 4
#endif

#if 0
!Modified for WRF:
#if (RWORDSIZE == 8)
 integer, parameter :: kind_rb = selected_real_kind(12) ! 8 byte real
#endif
#if (RWORDSIZE == 4)
 integer, parameter :: kind_rb = selected_real_kind(6)  ! 4 byte real
#endif
#else
 integer, parameter :: kind_rb = kind(1.0)              ! native real
#endif

 end module parkind

!=================================================================================================================
 module ra_rrtmg_shared


 implicit none
 private
 public:: reicalc, &
          relcalc, &
          inirad,  &
          o3data,  &
          retab

 real retab(95)
 data retab /                                                   &
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,           &
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,  &
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,  &
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,  &
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,  &
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,  &
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,  &
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,  &
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,  &
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,  &
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,  &
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,  &
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,  &
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,  &
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,  &
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639/  
 save retab


 contains


!=================================================================================================================
  subroutine relcalc(ncol, pcols, pver, t, landfrac, landm, icefrac, rel, snowh)
!----------------------------------------------------------------------------------------------------------------- 
! 
! Purpose: 
! Compute cloud water size
! 
! Method: 
! analytic formula following the formulation originally developed by J. T. Kiehl
! 
! Author: Phil Rasch
! 
!-----------------------------------------------------------------------------------------------------------------
    implicit none
!------------------------------Arguments--------------------------------------------------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol
    integer, intent(in) :: pcols, pver
    real, intent(in) :: landfrac(pcols)      ! Land fraction
    real, intent(in) :: icefrac(pcols)       ! Ice fraction
    real, intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
    real, intent(in) :: landm(pcols)         ! Land fraction ramping to zero over ocean
    real, intent(in) :: t(pcols,pver)        ! Temperature

!
! Output arguments
!
    real, intent(out) :: rel(pcols,pver)      ! Liquid effective drop size (microns)
!
!---------------------------Local workspace-----------------------------------------------------------------------
!
    integer i,k           ! Lon, lev indices
    real tmelt            ! freezing temperature of fresh water (K)
    real rliqland         ! liquid drop size if over land
    real rliqocean        ! liquid drop size if over ocean
    real rliqice          ! liquid drop size if over sea ice
!
!-----------------------------------------------------------------------------------------------------------------
!
    tmelt = 273.16
    rliqocean = 14.0
    rliqice   = 14.0
    rliqland  = 8.0
    do k=1,pver
       do i=1,ncol
! jrm Reworked effective radius algorithm
          ! Start with temperature-dependent value appropriate for continental air
          ! Note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * min(1.0,max(0.0,(tmelt-t(i,k))*0.05))
          ! Modify for snow depth over land
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0,max(0.0,snowh(i)*10.))
          ! Ramp between polluted value over land to clean value over ocean.
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0,max(0.0,1.0-landm(i)))
          ! Ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k) = rel(i,k) + (rliqice-rel(i,k)) * min(1.0,max(0.0,icefrac(i)))
! end jrm
       end do
    end do
  end subroutine relcalc

!=================================================================================================================
  subroutine reicalc(ncol, pcols, pver, t, re)
!----------------------------------------------------------------------------------------------------------------- 
!
    integer, intent(in) :: ncol, pcols, pver
    real, intent(out) :: re(pcols,pver)
    real, intent(in) :: t(pcols,pver)
    real corr
    integer i
    integer k
    integer index
!
!   Tabulated values of re(T) in the temperature interval
!   180 K -- 274 K; hexagonal columns assumed:
!
!
    do k=1,pver
       do i=1,ncol
          index = int(t(i,k)-179.)
          index = min(max(index,1),94)
          corr = t(i,k) - int(t(i,k))
          re(i,k) = retab(index)*(1.-corr) &
               +retab(index+1)*corr
!         re(i,k) = amax1(amin1(re(i,k),30.),10.)
       end do
    end do
!
  end subroutine reicalc

!=================================================================================================================
 subroutine inirad(o3prof,plev,kts,kte)
 implicit none
!-----------------------------------------------------------------------------------------------------------------
 integer,intent(in):: kts,kte
 real,dimension(kts:kte+2),intent(in):: plev
 real,dimension(kts:kte+1),intent(inout):: o3prof

!local vars:
 integer:: k

!
!compute ozone mixing ratio:
!
 do k=kts,kte+1
    o3prof(k)=0.
 enddo
 call o3data(o3prof,plev,kts,kte)

 end subroutine inirad
                                                                                 
!=================================================================================================================
 subroutine o3data(o3prof,plev,kts,kte)
 implicit none
!-----------------------------------------------------------------------------------------------------------------
 integer,intent(in):: kts, kte
 real,dimension(kts:kte+2),intent(in)::  plev
 real,dimension(kts:kte+1),intent(inout):: o3prof

!local vars:
 integer:: k,jj

 real:: prlevh(kts:kte+2),ppwrkh(32),            &
        o3wrk(31),ppwrk(31),o3sum(31),ppsum(31), &
        o3win(31),ppwin(31),o3ann(31),ppann(31)

 real:: pb1,pb2,pt1,pt2

 data o3sum  /5.297e-8,5.852e-8,6.579e-8,7.505e-8,         &
    8.577e-8,9.895e-8,1.175e-7,1.399e-7,1.677e-7,2.003e-7, &
    2.571e-7,3.325e-7,4.438e-7,6.255e-7,8.168e-7,1.036e-6, &
    1.366e-6,1.855e-6,2.514e-6,3.240e-6,4.033e-6,4.854e-6, &
    5.517e-6,6.089e-6,6.689e-6,1.106e-5,1.462e-5,1.321e-5, &
    9.856e-6,5.960e-6,5.960e-6/

 data ppsum  /955.890,850.532,754.599,667.742,589.841,       &
    519.421,455.480,398.085,347.171,301.735,261.310,225.360, &
    193.419,165.490,141.032,120.125,102.689, 87.829, 75.123, &
    64.306, 55.086, 47.209, 40.535, 34.795,  29.865, 19.122, &
    9.277,  4.660,  2.421,  1.294,  0.647/                                 
                                                                             
 data o3win  /4.629e-8,4.686e-8,5.017e-8,5.613e-8,         &
    6.871e-8,8.751e-8,1.138e-7,1.516e-7,2.161e-7,3.264e-7, &
    4.968e-7,7.338e-7,1.017e-6,1.308e-6,1.625e-6,2.011e-6, &
    2.516e-6,3.130e-6,3.840e-6,4.703e-6,5.486e-6,6.289e-6, &
    6.993e-6,7.494e-6,8.197e-6,9.632e-6,1.113e-5,1.146e-5, &
    9.389e-6,6.135e-6,6.135e-6/

 data ppwin  /955.747,841.783,740.199,649.538,568.404,       &
    495.815,431.069,373.464,322.354,277.190,237.635,203.433, &
    174.070,148.949,127.408,108.915, 93.114, 79.551, 67.940, &
    58.072, 49.593, 42.318, 36.138, 30.907, 26.362, 16.423,  &
    7.583,  3.620,  1.807,  0.938,  0.469/

 do k=1,31
    ppann(k)=ppsum(k)
 enddo

 o3ann(1)=0.5*(o3sum(1)+o3win(1))
 do k=2,31
    o3ann(k)=o3win(k-1)+(o3win(k)-o3win(k-1))/(ppwin(k)-ppwin(k-1))* &
             (ppsum(k)-ppwin(k-1))
 enddo

 do k=2,31
    o3ann(k)=0.5*(o3ann(k)+o3sum(k))
 enddo

 do k=1,31
    o3wrk(k)=o3ann(k)
    ppwrk(k)=ppann(k)
 enddo

!                                                                                
!calculate half pressure levels for model and data levels.
!plev is total p at model levels, from bottom to top.
!plev is in mb.

 do k=kts,kte+2
    prlevh(k)=plev(k)
 enddo
                                                                                
 ppwrkh(1)=1100.
 do k=2,31
    ppwrkh(k)=(ppwrk(k)+ppwrk(k-1))/2.
 enddo
 ppwrkh(32)=0.
 do k=kts,kte+1
    do 25 jj=1,31
       if((-(prlevh(k)-ppwrkh(jj))).ge.0.)then
          pb1=0.
       else
          pb1=prlevh(k)-ppwrkh(jj)
       endif
       if((-(prlevh(k)-ppwrkh(jj+1))).ge.0.)then
          pb2=0.
       else
          pb2=prlevh(k)-ppwrkh(jj+1)
       endif
       if((-(prlevh(k+1)-ppwrkh(jj))).ge.0.)then
          pt1=0.
       else
          pt1=prlevh(k+1)-ppwrkh(jj)
       endif
       if((-(prlevh(k+1)-ppwrkh(jj+1))).ge.0.)then
          pt2=0.
       else
          pt2=prlevh(k+1)-ppwrkh(jj+1)
       endif
       o3prof(k)=o3prof(k)+(pb2-pb1-pt2+pt1)*o3wrk(jj)
 25 continue
    o3prof(k)=o3prof(k)/(prlevh(k)-prlevh(k+1))
 enddo
                                                                                
 end subroutine o3data

!=================================================================================================================
 end module ra_rrtmg_shared
!=================================================================================================================
