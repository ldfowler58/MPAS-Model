!=================================================================================================================
 module module_mynn_get_pblh
!=================================================================================================================
 use mpas_kind_types

 implicit none
 private
 public:: mynn_get_pblh


 contains
 

!=================================================================================================================
 subroutine mynn_get_pblh(kts,kte,zi,thetav1d,qke1d,zw1d,dz1d,landsea,kzi)
!=================================================================================================================
!notes on the pblh formulation
!
!the 1.5-theta-increase method defines pbl heights as the level at
!which the potential temperature first exceeds the minimum potential
!temperature within the boundary layer by 1.5 k. when applied to
!observed temperatures, this method has been shown to produce pbl-
!height estimates that are unbiased relative to profiler-based
!estimates (nielsen-gammon et al. 2008). however, their study did not
!include lljs. banta and pichugina (2008) show that a tke-based
!threshold is a good estimate of the pbl height in lljs. therefore,
!a hybrid definition is implemented that uses both methods, weighting
!the tke-method more during stable conditions (pblh < 400 m).
!a variable tke threshold (tkeeps) is used since no hard-wired
!value could be found to work best in all conditions.

!input arguments:
 integer,intent(in) :: kts,kte
 real(kind=RKIND),intent(in):: landsea
 real(kind=RKIND),intent(in),dimension(kts:kte)  :: thetav1d,qke1d,dz1d
 real(kind=RKIND),intent(in),dimension(kts:kte+1):: zw1d

!output arguments:
 integer,intent(out):: kzi
 real(kind=RKIND),intent(out):: zi

!local variables:
 integer:: i,j,k,kthv,ktke,kzi2

 real(kind=RKIND),parameter:: sbl_lim  = 200. !upper limit of stable bl height (m).
 real(kind=RKIND),parameter:: sbl_damp = 400. !transition length for blending (m).

 real(kind=RKIND):: pblh_tke,qtke,qtkem1,wt,maxqke,tkeeps,minthv
 real(kind=RKIND):: delt_thv !delta theta-v; dependent on land/sea point

!-----------------------------------------------------------------------------------------------------------------

!add kpbl (kzi)
!kzi2 is the tke-based part of the hybrid kpbl
 kzi = 2
 kzi2= 2

!find min thetav in the lowest 200 m agl
 k = kts+1
 kthv = 1
 minthv = 9.e9
 do while (zw1d(k) .le. 200.)
!   do k=kts+1,kte-1
    if(minthv > thetav1d(k)) then
       minthv = thetav1d(k)
       kthv = k
    endif
    k = k+1
!   if (zw1d(k) .gt. sbl_lim) exit
 enddo

!find thetav-based pblh (best for daytime).
 zi = 0.
 k = kthv+1
 if((landsea-1.5).ge.0)then
    !water
    delt_thv = 0.75
 else
    !land
    delt_thv = 1.25
 endif

 zi=0.
 k = kthv+1
 do k = kts+1,kte-1
    if(thetav1d(k) .ge. (minthv + delt_thv))then
       !kzi = max(k-1,1)
       zi = zw1d(k) - dz1d(k-1)* &
            min((thetav1d(k)-(minthv + delt_thv))/ &
            max(thetav1d(k)-thetav1d(k-1),1e-6),1.0)
       kzi = max(k-1,1) + nint((zi-zw1d(k-1))/dz1d(k-1))
    endif
    !k = k+1
    if(k .eq. kte-1) zi = zw1d(kts+1) !exit safeguard
    if(zi .ne. 0.0) exit
 enddo

!for stable boundary layers, use tke method to complement the
!thetav-based definition (when the theta-v based pblh is below ~0.5 km).
!the tanh weighting function will make the tke-based definition negligible
!when the theta-v-based definition is above ~1 km.
 ktke = 1
 maxqke = max(qke1d(kts),0.)
!use 5% of tke max (kosovic and curry, 2000; jas)
!tkeeps = maxtke/20. = maxqke/40.
 tkeeps = maxqke/40.
 tkeeps = max(tkeeps,0.02) !0.025)
 pblh_tke=0.

 k = ktke+1
 do k = kts+1,kte-1
    !qke can be negative (if ckmod == 0)... make tke non-negative.
    qtke   = max(qke1d(k)/2.,0.)   ! maximum tke
    qtkem1 = max(qke1d(k-1)/2.,0.)
    if(qtke .le. tkeeps) then
       !kzi2 = max(k-1,1)
       pblh_tke = zw1d(k) - dz1d(k-1)* &
                  min((tkeeps-qtke)/max(qtkem1-qtke, 1e-6), 1.0)
       !in case of near zero tke, set pblh = lowest level.
       pblh_tke = max(pblh_tke,zw1d(kts+1))
       kzi2 = max(k-1,1) + nint((pblh_tke-zw1d(k-1))/dz1d(k-1))
    endif
    !k = k+1
    if (k .eq. kte-1) pblh_tke = zw1d(kts+1) !exit safeguard
    if (pblh_tke .ne. 0.) exit
 enddo

!with tke advection turned on, the tke-based pblh can be very large
!in grid points with convective precipitation (> 8 km!),
!so an artificial limit is imposed to not let pblh_tke exceed the
!theta_v-based pbl height +/- 350 m.
!this has no impact on 98-99% of the domain, but is the simplest patch
!that adequately addresses these extremely large pblhs.
 pblh_tke = min(pblh_tke,zi+350.)
 pblh_tke = max(pblh_tke,max(zi-350.,10.))

 wt = .5*tanh((zi - sbl_lim)/sbl_damp) + .5
 if(maxqke <= 0.05) then
    !cold pool situation - default to theta_v-based def
 else
    !blend the two pblh types here:
    zi=pblh_tke*(1.-wt) + zi*wt
 endif

!add kpbl (kzi) for coupling to some cu schemes
!kzi = max(int(kzi2*(1.-wt) + kzi*wt),1)

!correction to the calculation of kzi by Joe Olson (2019-12-10):
 do k = kts+1, kte-1
    if(zw1d(k) >= zi) then
       kzi = k-1
       exit
    endif
 enddo

 end subroutine mynn_get_pblh

!=================================================================================================================
 end module module_mynn_get_pblh
!=================================================================================================================

