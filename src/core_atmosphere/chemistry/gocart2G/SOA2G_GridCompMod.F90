! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module SOA2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: MPAS_LOG_CRIT
 use mpas_log

 use SOA2G_StateSpecs,only: SOA2G_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: Avogadro = 6.02214076e23

!--- molecular weight of dry air (grams):
 real(kind=RKIND),parameter:: fMassAir  = 28.97_RKIND


 type,public:: SOA2G_GridComp
    real(kind=RKIND):: cdt               ! chemistry timestep (secs)
    real(kind=RKIND):: ratPOM = 1._RKIND ! ratio of POM to SOA (not sure)

    contains
       procedure:: emissions_GridComp => emissions_SOA2G_GridComp
       procedure:: load_GridComp      => load_SOA2G_GridComp
 end type SOA2G_GridComp


 contains


!==================================================================================================================
 subroutine load_SOA2G_GridComp(self)
!==================================================================================================================

!--- inout arguments:
 class(SOA2G_GridComp),intent(inout) :: self

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_SOA2G_GridComp:')


 call mpas_log_write('--- end subroutine load_SOA2G_GridComp.')

 end subroutine load_SOA2G_GridComp

!==================================================================================================================
 subroutine emissions_SOA2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

 class(SOA2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(SOA2G_State),intent(inout):: self

!--- local variables and arrays:
 integer:: i,j,k,km

 real(kind=RKIND):: cdt
 real(kind=RKIND),dimension(:,:,:),allocatable:: rk_OA_OH
 real(kind=RKIND),dimension(:,:,:),allocatable:: dsoap
 real(kind=RKIND),dimension(:,:,:),allocatable:: dOAanth,dOAbiob
 real(kind=RKIND),dimension(:,:,:),allocatable:: fanth 

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_SOA2G_GridComp:')

 cdt = self_params%cdt


!--- production of the lumped SOA precursors (SOAP) as functions of CO anthropogenic and
!    biomass burning emissions:
 call SOAemission( &
         its          = its,                &
         ite          = ite,                &
         jts          = jts,                &
         jte          = jte,                &
         kts          = kts,                &
         kte          = kte,                &
         cdt          = self_params%cdt,    &
         ratPOM       = self_params%ratPOM, &
         zpbl         = self%zpbl,          &
         delp         = self%delp,          &
         delz         = self%delz,          &
         ple          = self%ple,           &
         zle          = self%zle,           &
         eocant1_src  = self%soap_anthro,   &
         biomass_src  = self%soap_biomass,  &
         biogenic_src = self%soap_biogenic, &
         soap_a       = self%soap_a,        &
         soap_bb      = self%soap_bb        &
                 )


!--- production of the lumped simple SOA due to oxidation by OH (following Kim et al. 2015):
 allocate(rk_OA_OH(its:ite,jts:jte,kts:kte))
 allocate(fanth(its:ite,jts:jte,kts:kte)   )
 allocate(dsoap(its:ite,jts:jte,kts:kte)   )
 allocate(dOAanth(its:ite,jts:jte,kts:kte) )
 allocate(dOAbiob(its:ite,jts:jte,kts:kte) )

 rk_OA_OH = 0._RKIND
 fanth    = 0._RKIND
 dsoap    = 0._RKIND
 dOAanth  = 0._RKIND
 dOAbiob  = 0._RKIND


 rk_OA_OH = 1.25d-11*Avogadro*self%soap_oh*self%airdens/fMassAir*(1.0e-6)*cdt 
 dsoap = (self%soap_a + self%soap_bb)*(1.-exp(-rk_OA_OH)) ! total loss of SOAP due to oxidation (kg/kg)

 where(dsoap .gt. 1.e-32)
    fanth        = self%soap_a/(self%soap_a + self%soap_bb)
    dOAanth      = fanth*dsoap                     ! loss of anthropogenic SOAP/production of lumped SOA (kg/kg)
    self%soap_a  = self%soap_a - fanth*dsoap       ! update anthropogenic SOA
    dOAbiob      = (1.-fanth)*dsoap                ! loss of biomass burning SOAP/production of lumped SOA (kg/kg)
    self%soap_bb = self%soap_bb - (1.-fanth)*dsoap ! update biomass biurning SOA
 end where

 self%soapa_prod  = self%airdens*dOAanth/cdt
 self%soapbb_prod = self%airdens*dOAbiob/cdt

 
 deallocate(rk_OA_OH)
 deallocate(fanth   )
 deallocate(dsoap   )
 deallocate(dOAanth )
 deallocate(dOAbiob )


 call mpas_log_write('--- end subroutine emissions_SOA2G_GridComp.')

 end subroutine emissions_SOA2G_GridComp

!==================================================================================================================
 subroutine SOAemission(its,ite,jts,jte,kts,kte,cdt,ratPOM,zpbl,delp,delz,ple,zle,eocant1_src, &
                        biomass_src,biogenic_src,soap_a,soap_bb)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

 real(kind=RKIND),intent(in):: cdt
 real(kind=RKIND),intent(in):: ratPOM
 real(kind=RKIND),intent(in),dimension(:,:),pointer:: zpbl
 real(kind=RKIND),intent(in),dimension(:,:,:),pointer:: delp,delz
 real(kind=RKIND),intent(in),dimension(:,:,:),pointer:: ple,zle

 real(kind=RKIND),intent(in),dimension(:,:),pointer:: eocant1_src,biomass_src,biogenic_src

!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(:,:,:),pointer:: soap_a,soap_bb

!--- local variables and arrays:
 integer:: i,j,k,kk
 integer,dimension(:,:),allocatable:: i100,i500,ipbl

 real(kind=RKIND):: eAnthro,eBiomass
 real(kind=RKIND):: f100,f500,fpbl,fbot,zfactor,zp
 real(kind=RKIND),dimension(:),allocatable:: zs
 real(kind=RKIND),dimension(:,:),allocatable:: p100,p500,ppbl
 real(kind=RKIND),dimension(:,:),allocatable:: srcAnthro,srcBiomass

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine SOAemission:')

 allocate(srcAnthro(its:ite,jts:jte) )
 allocate(srcBiomass(its:ite,jts:jte))

 allocate(i100(its:ite,jts:jte))
 allocate(i500(its:ite,jts:jte))
 allocate(ipbl(its:ite,kts:kte))

 allocate(zs(kts:kte+1))
 allocate(p100(its:ite,jts:jte))
 allocate(p500(its:ite,jts:jte))
 allocate(ppbl(its:ite,kts:kte))


 eAnthro  = ratPOM
 eBiomass = ratPOM


!--- finds the pressure of the 100 meters, 500 meters, and PBL height:
 do j = jts,jte
    do i = its,ite
       zp = max(zpbl(i,j),100._RKIND)
       zs(kte+1) = 0._RKIND
       do k = kte,kts,-1
          zs(k) = zs(k+1) + delz(i,j,k)
       enddo

       i100(i,j) = kte
       i500(i,j) = kte
       ipbl(i,j) = kte
       p100(i,j) = ple(i,j,kte)
       p500(i,j) = ple(i,j,kte)
       ppbl(i,j) = ple(i,j,kte)
       do k = kte+1,kts+1,-1
          if(zs(k).lt.100._RKIND .and. zs(k-1).ge.100._RKIND) then
             i100(i,j) = k-1
             zfactor = (zs(k-1)-100.)/(zs(k-1)-zs(k))
             p100(i,j) = ple(i,j,k-1)-zfactor*(ple(i,j,k-1)-ple(i,j,k))
          elseif(zs(k).lt.500._RKIND .and. zs(k-1).ge.500._RKIND) then
             i500(i,j) = k-1
             zfactor = (zs(k-1)-500.)/(zs(k-1)-zs(k))
             p500(i,j) = ple(i,j,k-1)-zfactor*(ple(i,j,k-1)-ple(i,j,k))
          elseif(zs(k).lt.zp .and. zs(k-1).ge.zp) then
             ipbl(i,j) = k-1
             zfactor = (zs(k-1)-zp)/(zs(k-1)-zs(k))
             ppbl(i,j) = ple(i,j,k-1)-zfactor*(ple(i,j,k-1)-ple(i,j,k))
          endif
       enddo

       do k = kte,kts,-1
          f100 = 0._RKIND
          zfactor = 1./(ple(i,j,kte+1)-p100(i,j))
          if(ple(i,j,k) .ge. p100(i,j)) then
             f100 = zfactor*delp(i,j,k)
          elseif(ple(i,j,k+1).ge.p100(i,j) .and. ple(i,j,k).lt.p100(i,j)) then
             f100 = zfactor*(ple(i,j,k+1)-p100(i,j))
          endif

          f500 = 0._RKIND
          zfactor = 1./(p100(i,j)-p500(i,j))
          if(ple(i,j,k+1).ge.p100(i,j) .and. ple(i,j,k).lt.p100(i,j) .and. ple(i,j,k).ge.p500(i,j)) then
             f500 = zfactor*(p100(i,j)-ple(i,j,k))
          elseif(ple(i,j,k+1).lt.p100(i,j) .and. ple(i,j,k).ge.p500(i,j)) then
             f500 = zfactor*delp(i,j,k)
          elseif(ple(i,j,k+1).ge.p500(i,j) .and. ple(i,j,k).lt.p500(i,j)) then
             f500 = zfactor*(ple(i,j,k+1)-p500(i,j))
          endif

          fpbl = 0._RKIND
          zfactor = 1./(ple(i,j,kte+1)-ppbl(i,j))
          if(ple(i,j,k) .ge. ppbl(i,j)) then
             fpbl = zfactor*delp(i,j,k)
          elseif(ple(i,j,k+1).ge.ppbl(i,j) .and. ple(i,j,k).lt.ppbl(i,j)) then
             fpbl = zfactor*(ple(i,j,k+1)-ppbl(i,j))
          endif

          srcAnthro(i,j) = f100*eAnthro*eocant1_src(i,j)
          srcBiomass(i,j) = fpbl*eBiomass*biomass_src(i,j)

          zfactor = cdt*grav/delp(i,j,k)
          soap_a(i,j,k) = soap_a(i,j,k) + zfactor*srcAnthro(i,j)
          soap_bb(i,j,k) = soap_bb(i,j,k) + zfactor*srcBiomass(i,j)
       enddo
    enddo
 enddo

 deallocate(srcAnthro )
 deallocate(srcBiomass)

 deallocate(i100)
 deallocate(i500)
 deallocate(ipbl)

 deallocate(zs)
 deallocate(p100)
 deallocate(p500)
 deallocate(ppbl)

 call mpas_log_write('--- end subroutine SOAemission.')

 end subroutine SOAemission

!==================================================================================================================
 end module SOA2G_GridCompMod
!==================================================================================================================
