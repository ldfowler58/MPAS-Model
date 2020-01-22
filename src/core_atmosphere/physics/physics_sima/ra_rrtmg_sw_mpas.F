!=================================================================================================================
 module rrtmg_sw_mpas
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

!use parrrsw,only: nbndsw,ngptsw,naerec

 implicit none
 public:: rrtmg_sw_mpas_run,           &
!         rrtmg_sw_mpas_timestep_init, &
          rrtmg_sw_mpas_timestep_final


 contains


!   call rrtmg_sw_run &
!              (dorrsw    , ncol     , nlay      , icld      , play      , plev      , &
!               tlay      , tlev     , tsfc      , h2ovmr    , o3vmr     , co2vmr    , &
!               ch4vmr    , n2ovmr   , o2vmr     , asdir     , asdif     , aldir     , &
!               aldif     , coszen   , adjes     , dyofyr    , scon      , inflgsw   , &
!               iceflgsw  , liqflgsw , cldfmcl   , taucmcl   , ssacmcl   , asmcmcl   , &
!               fsfcmcl   , ciwpmcl  , clwpmcl   , cswpmcl   , reicmcl   , relqmcl   , &
!               resnmcl   , tauaer   , ssaaer    , asmaer    , ecaer     , swuflx    , &
!               swdflx    , swhr     , swuflxc   , swdflxc   , swhrc     , swuflxcln , &
!               swdflxcln , aer_opt  , sibvisdir , sibvisdif , sibnirdir , sibnirdif , &
!               swdkdir   , swdkdif  , swdkdirc  , calc_clean_atm_diag   , errmsg    , &
!               errflg)



!=================================================================================================================
 subroutine rrtmg_sw_mpas_run(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out)::errflg

 end subroutine rrtmg_sw_mpas_run

!=================================================================================================================
!subroutine rrtmg_sw_mpas_timestep_init &
!=================================================================================================================
!end subroutine rrtmg_sw_mpas_timestep_init

!=================================================================================================================
 subroutine rrtmg_sw_mpas_timestep_final &
                  (dorrsw,pi3d,alswvisdir,alswvisdif,alswnirdir,alswnirdif,coszr,                  &
                   gsw,swhr,swhrc,swuflx,swuflxc,swdflx,swdflxc,swuflxcln,swdflxcln,               &
                   sibvisdir,sibvisdif,sibnirdir,sibnirdif,swdkdir,swdkdirc,swdkdif,               &
                   rthratensw,swcf,swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,            &
                   swuptcln,swdntcln,swupbcln,swdnbcln,swvisdir,swvisdif,swnirdir,swnirdif,        &
                   swddir,swddni,swddif,swddirc,swddnic,swdownc,swupflx,swupflxc,swdnflx,swdnflxc, &
                   ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!input arguments:
 logical,intent(in),dimension(its:ite):: dorrsw

 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+1):: &
    swhr,       &! all-sky shortwave heating rate (K/day).
    swhrc        ! clear-sky shortwave heating rate (K/day).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    swuflx,     &! all-sky upward shortwave radiation flux (W/m2).
    swuflxc,    &! clear-sky upward shortwave radiation flux (W/m2)
    swdflx,     &! all-sky downward shortwave radiation flux (W/m2).
    swdflxc      ! clear-sky downward shortwave radiation flux (W/m2).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    swuflxcln,  &!
    swdflxcln    !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    sibvisdir,  &!
    sibvisdif,  &!
    sibnirdir,  &!
    sibnirdif    !

 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: &
    swdkdir,    &!
    swdkdirc,   &!
    swdkdif      !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    coszr        !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    pi3d         ! exner function (-).

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
    alswvisdir, &!
    alswvisdif, &!
    alswnirdir, &!
    alswnirdif   !

!output arguments:
 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    gsw          ! all-sky net surface radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swcf,       &! all-sky shortwave cloud radiative forcing (W/m2).
    swupt,      &! all-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swuptc,     &! clear-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swdnt,      &! all-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swdntc,     &! clear-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swupb,      &! all-sky shortwave surface upward radiation (W/m2).
    swupbc,     &! clear-sky shortwave surface upward radiation (W/m2).
    swdnb,      &! all-sky shortwave surface downward radiation (W/m2).
    swdnbc       ! clear-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swuptcln,   &!
    swdntcln,   &!
    swupbcln,   &!
    swdnbcln     !

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swvisdir,   &!
    swvisdif,   &!
    swnirdir,   &!
    swnirdif     !

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
      swddir,   &!
      swddni,   &!
      swddif,   &!
      swddirc,  &!
      swddnic,  &!
      swdownc    !

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    rthratensw   !

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme+1,jms:jme ),optional:: &
    swupflx,    &!
    swupflxc,   &!
    swdnflx,    &!
    swdnflxc     !

!local variables and arrays:
 integer:: i,k,j
 integer:: calc_clean_atm_diag
 integer:: sf_surface_physics

!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_sw_mpas_timestep_final:')

!--- initialization of calc_clean_atm_diag:
 calc_clean_atm_diag = 0
 if(present(swuptcln) .and. present(swdntcln) .and. present(swupbcln) .and. present(swdnbcln)) &
    calc_clean_atm_diag = 1

!--- initialization sf_surface_physics:
 sf_surface_physics = 0
 if(present(alswvisdir) .and. present(alswvisdif) .and. present(alswnirdir) .and. present(alswnirdif)) &
    sf_surface_physics = 8

!--- latitude loop:
 j_loop: do j = jts,jte

    !--- longitude loop:
    i_loop: do i = its,ite
       if(dorrsw(i)) then 
          gsw(i,j) = swdflx(i,1) - swuflx(i,1)
          do k = kts,kte
             rthratensw(i,k,j) = swhr(i,k)/86400.
             rthratensw(i,k,j) = rthratensw(i,k,j)/pi3d(i,k,j)
          enddo

          !optional arguments:
          if(present(swupt)) then
             swcf(i,j) = (swdflx(i,kte+2)-swuflx(i,kte+2)) - (swdflxc(i,kte+2)-swuflxc(i,kte+2))
             swupt(i,j)  = swuflx(i,kte+2)
             swuptc(i,j) = swuflxc(i,kte+2)
             swdnt(i,j)  = swdflx(i,kte+2)
             swdntc(i,j) = swdflxc(i,kte+2)
             swupb(i,j)  = swuflx(i,1)
             swupbc(i,j) = swuflxc(i,1)
             swdnb(i,j)  = swdflx(i,1)
             swdnbc(i,j) = swdflxc(i,1)
          endif
          if(present(swupflx)) then
             do k = kts, kte+1
                swupflx(i,k,j)  = swuflx(i,k)
                swupflxc(i,k,j) = swuflxc(i,k)
                swdnflx(i,k,j)  = swdflx(i,k)
                swdnflxc(i,k,j) = swdflxc(i,k)
             enddo
          endif
          if(calc_clean_atm_diag .gt. 0)then
             swuptcln(i,j) = swuflxcln(i,kte+2)
             swdntcln(i,j) = swdflxcln(i,kte+2)
             swupbcln(i,j) = swuflxcln(i,1)
             swdnbcln(i,j) = swdflxcln(i,1)
          endif
          if(sf_surface_physics .eq. 8) then
             swvisdir(i,j) = sibvisdir(i,1)
             swvisdif(i,j) = sibvisdif(i,1)
             swnirdir(i,j) = sibnirdir(i,1)
             swnirdif(i,j) = sibnirdif(i,1)
          endif
          if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
             present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
             swddir(i,j)   = swdkdir(i,1)
             swddni(i,j)   = swddir(i,j) / coszr(i,j)
             swddif(i,j)   = swdkdif(i,1)
             swdownc(i, j) = swdflxc(i,1)
             swddirc(i,j)  = swdkdirc(i,1)
             swddnic(i,j)  = swddirc(i,j) / coszr(i,j)
          endif
 
       else

          gsw(i,j) = 0._kind_phys
          do k = kts,kte
             rthratensw(i,k,j) = 0._kind_phys
          enddo

          !optional arguments:
          if(present(swupt)) then
             swcf(i,j)   = 0._kind_phys
             swupt(i,j)  = 0._kind_phys
             swuptc(i,j) = 0._kind_phys
             swdnt(i,j)  = 0._kind_phys
             swdntc(i,j) = 0._kind_phys
             swupb(i,j)  = 0._kind_phys
             swupbc(i,j) = 0._kind_phys
             swdnb(i,j)  = 0._kind_phys
             swdnbc(i,j) = 0._kind_phys
          endif
          if(present(swupflx)) then
             do k = kts, kte+1
                swupflx(i,k,j)  = 0._kind_phys
                swupflxc(i,k,j) = 0._kind_phys
                swdnflx(i,k,j)  = 0._kind_phys
                swdnflxc(i,k,j) = 0._kind_phys
             enddo
          endif
          if(calc_clean_atm_diag .gt. 0) then
             swuptcln(i,j) = 0._kind_phys
             swdntcln(i,j) = 0._kind_phys
             swupbcln(i,j) = 0._kind_phys
             swdnbcln(i,j) = 0._kind_phys
          endif
          if(sf_surface_physics .eq. 8) then
             swvisdir(i,j) = 0._kind_phys
             swvisdif(i,j) = 0._kind_phys
             swnirdir(i,j) = 0._kind_phys
             swnirdif(i,j) = 0._kind_phys
          endif
          if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
             present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
             swddir(i,j)   = 0._kind_phys
             swddif(i,j)   = 0._kind_phys
             swdownc(i, j) = 0._kind_phys
             swddirc(i,j)  = 0._kind_phys
             swddnic(i,j)  = 0._kind_phys
          endif

       endif !endif daytime grid columns.

    enddo i_loop !end longitude loop.

 enddo j_loop !end latitude loop.

 call mpas_log_write('--- end subroutine rrtmg_sw_mpas_timestep_final:')
 call mpas_log_write(' ')

 end subroutine rrtmg_sw_mpas_timestep_final

!=================================================================================================================
 end module rrtmg_sw_mpas
!=================================================================================================================



