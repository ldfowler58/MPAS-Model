! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module GOCART2G_GridCompMod
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_derived_types,only: MPAS_LOG_ERR
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_instance
 use GOCART2G_StateSpecs

 use CA2G_bc_StateSpecs
 use CA2G_br_StateSpecs
 use CA2G_oc_StateSpecs
 use DU2G_StateSpecs
 use NI2G_StateSpecs
 use SS2G_StateSpecs
 use SU2G_StateSpecs


 implicit none
 private


!--- types needed to define GOCART2G:
 type,public:: GOCART2G_GridComp
    real(kind=RKIND),dimension(:),allocatable:: wavelengths_profile
    real(kind=RKIND),dimension(:),allocatable:: wavelengths_vertint

    contains
       procedure:: processes_GridComp => processes_gocart2G_GridComp
 end type GOCART2G_GridComp


 contains


!==================================================================================================================
 subroutine processes_GOCART2G_GridComp(self_params,self,CA2G_bc,CA2G_br,CA2G_oc,DU2G,NI2G,SS2G,SU2G, &
                                         its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 class(CA2G_bc_State),intent(in):: CA2G_bc
 class(CA2G_br_State),intent(in):: CA2G_br
 class(CA2G_oc_State),intent(in):: CA2G_oc
 class(DU2G_State),intent(in):: DU2G
 class(NI2G_State),intent(in):: NI2G
 class(SS2G_State),intent(in):: SS2G
 class(SU2G_State),intent(in):: SU2G

 integer,intent(in):: its,ite,jts,jte,kts,kte


!--- inout arguments:
 class(GOCART2G_GridComp),intent(inout):: self_params
 class(GOCART2G_State),intent(inout):: self

!--- local variables:
 integer:: w,nw_profile,nw_vertint
 integer:: ind550,ind532
 real(kind=RKIND):: c1,c2,c3
 real(kind=RKIND),dimension(:,:),allocatable:: tau1,tau2

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_GOCART2G_GridComp:')


!--- initialization of wavelengths_profile and wavelengths_vertint:
 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)
!call mpas_log_write('--- nw_profile = $i',intArgs=(/nw_profile/))
!call mpas_log_write('--- nw_vertint = $i',intArgs=(/nw_vertint/))

 if(.not.allocated(self_params%wavelengths_profile)) allocate(self_params%wavelengths_profile(nw_profile)) 
 if(.not.allocated(self_params%wavelengths_vertint)) allocate(self_params%wavelengths_vertint(nw_vertint))

 do w = 1,nw_profile
    self_params%wavelengths_profile(w) = wavelengths_for_profile_aop_in_nm(w)
 enddo
 do w = 1,nw_vertint
    self_params%wavelengths_vertint(w) = wavelengths_for_vertically_integrated_aop_in_nm(w)
 enddo


 if(associated(self%totangstr)) then
    ind550 = 0
    do w = 1, size(self_params%wavelengths_vertint) ! find index for 550nm to compute total angstrom
       if((self_params%wavelengths_vertint(w)*1.e-9 .ge. 5.49e-7) .and. &
          (self_params%wavelengths_vertint(w)*1.e-9 .le. 5.51e-7)) then
           ind550 = w
           exit
       endif
    enddo
!   call mpas_log_write('--- ind550     = $i',intArgs=(/ind550/))

    if(ind550 == 0) then
       call mpas_log_write('550nm wavelength not in GOCART2G_instance: cannot compute TOTANGSTR', &
                           messageType=MPAS_LOG_ERR)
    endif

    self%totangstr = 0._RKIND
    if(.not.allocated(tau1)) allocate(tau1(its:ite,jts:jte))
    if(.not.allocated(tau2)) allocate(tau2(its:ite,jts:jte))
    tau1(:,:) = tiny(1.0)
    tau2(:,:) = tiny(1.0)
    c1 = -log(470./550.)
    c2 = -log(870./550.)
    c3 = -log(470./870.)
!   call mpas_log_write('--- c1 = $r',realArgs=(/c1/))
!   call mpas_log_write('--- c2 = $r',realArgs=(/c2/))
!   call mpas_log_write('--- c3 = $r',realArgs=(/c3/))
 endif


!--- initialization of total AOPs and PM 2.5 diagnostics:
 if(associated(self%totexttau)     ) self%totexttau      = 0._RKIND
 if(associated(self%totstexttau)   ) self%totstexttau    = 0._RKIND
 if(associated(self%totscatau)     ) self%totscatau      = 0._RKIND
 if(associated(self%totstscatau)   ) self%totstscatau    = 0._RKIND
 if(associated(self%totextt25)     ) self%totextt25      = 0._RKIND
 if(associated(self%totscat25)     ) self%totscat25      = 0._RKIND
 if(associated(self%totexttfm)     ) self%totexttfm      = 0._RKIND
 if(associated(self%totscatfm)     ) self%totscatfm      = 0._RKIND
 if(associated(self%totextcoef)    ) self%totextcoef     = 0._RKIND
 if(associated(self%totextcoefrh20)) self%totextcoefrh20 = 0._RKIND
 if(associated(self%totextcoefrh80)) self%totextcoefrh80 = 0._RKIND
 if(associated(self%totscacoef)    ) self%totscacoef     = 0._RKIND
 if(associated(self%totscacoefrh20)) self%totscacoefrh20 = 0._RKIND
 if(associated(self%totscacoefrh80)) self%totscacoefrh80 = 0._RKIND
 if(associated(self%totbckcoef)    ) self%totbckcoef     = 0._RKIND
 if(associated(self%totabcktoa)    ) self%totabcktoa     = 0._RKIND
 if(associated(self%totabcksfc)    ) self%totabcksfc     = 0._RKIND
 if(associated(self%pm)            ) self%pm             = 0._RKIND
 if(associated(self%pm25)          ) self%pm25           = 0._RKIND
 if(associated(self%pm_rh35)       ) self%pm_rh35        = 0._RKIND
 if(associated(self%pm25_rh35)     ) self%pm25_rh35      = 0._RKIND
 if(associated(self%pm_rh50)       ) self%pm_rh50        = 0._RKIND
 if(associated(self%pm25_rh50)     ) self%pm25_rh50      = 0._RKIND
 if(associated(self%pso4tot)       ) self%pso4tot        = 0._RKIND


!--- add AOPs from black carbon to total AOPs:
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(CA2G_bc%bcexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + CA2G_bc%bcexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(CA2G_bc%bcstexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + CA2G_bc%bcstexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(CA2G_bc%bcscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + CA2G_bc%bcscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(CA2G_bc%bcstscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + CA2G_bc%bcstscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(CA2G_bc%bcexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + CA2G_bc%bcexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(CA2G_bc%bcscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + CA2G_bc%bcscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(CA2G_bc%bcexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + CA2G_bc%bcexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(CA2G_bc%bcscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + CA2G_bc%bcscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(CA2G_bc%bcextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + CA2G_bc%bcextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(CA2G_bc%bcextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + CA2G_bc%bcextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(CA2G_bc%bcextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + CA2G_bc%bcextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(CA2G_bc%bcscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + CA2G_bc%bcscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(CA2G_bc%bcscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + CA2G_bc%bcscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(CA2G_bc%bcscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + CA2G_bc%bcscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(CA2G_bc%bcbckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + CA2G_bc%bcbckcoef(:,:,:,w)
 enddo

 if(associated(self%pm)        .and. associated(CA2G_bc%bcsmass)) self%pm        = self%pm        + CA2G_bc%bcsmass
 if(associated(self%pm25)      .and. associated(CA2G_bc%bcsmass)) self%pm25      = self%pm25      + CA2G_bc%bcsmass
 if(associated(self%pm_rh35)   .and. associated(CA2G_bc%bcsmass)) self%pm_rh35   = self%pm_rh35   + CA2G_bc%bcsmass
 if(associated(self%pm25_rh35) .and. associated(CA2G_bc%bcsmass)) self%pm25_rh35 = self%pm25_rh35 + CA2G_bc%bcsmass
 if(associated(self%pm_rh50)   .and. associated(CA2G_bc%bcsmass)) self%pm_rh50   = self%pm_rh50   + CA2G_bc%bcsmass
 if(associated(self%pm25_rh50) .and. associated(CA2G_bc%bcsmass)) self%pm25_rh50 = self%pm25_rh50 + CA2G_bc%bcsmass

 if(associated(self%totangstr) .and. associated(CA2G_bc%bcexttau) .and. associated(CA2G_bc%bcangstr)) then
    tau1 = tau1 + CA2G_bc%bcexttau(:,:,ind550)*exp(c1*CA2G_bc%bcangstr)
    tau2 = tau2 + CA2G_bc%bcexttau(:,:,ind550)*exp(c2*CA2G_bc%bcangstr)
 endif


!--- add AOPs from black carbon to total AOPs:
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(CA2G_br%brexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + CA2G_br%brexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(CA2G_br%brstexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + CA2G_br%brstexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(CA2G_br%brscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + CA2G_br%brscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(CA2G_br%brstscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + CA2G_br%brstscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(CA2G_br%brexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + CA2G_br%brexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(CA2G_br%brscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + CA2G_br%brscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(CA2G_br%brexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + CA2G_br%brexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(CA2G_br%brscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + CA2G_br%brscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(CA2G_br%brextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + CA2G_br%brextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(CA2G_br%brextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + CA2G_br%brextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(CA2G_br%brextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + CA2G_br%brextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(CA2G_br%brscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + CA2G_br%brscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(CA2G_br%brscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + CA2G_br%brscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(CA2G_br%brscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + CA2G_br%brscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(CA2G_br%brbckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + CA2G_br%brbckcoef(:,:,:,w)
 enddo

 if(associated(self%pm)        .and. associated(CA2G_br%brsmass)) self%pm        = self%pm        + CA2G_br%brsmass
 if(associated(self%pm25)      .and. associated(CA2G_br%brsmass)) self%pm25      = self%pm25      + CA2G_br%brsmass
 if(associated(self%pm_rh35)   .and. associated(CA2G_br%brsmass)) self%pm_rh35   = self%pm_rh35   + CA2G_br%brsmass
 if(associated(self%pm25_rh35) .and. associated(CA2G_br%brsmass)) self%pm25_rh35 = self%pm25_rh35 + CA2G_br%brsmass
 if(associated(self%pm_rh50)   .and. associated(CA2G_br%brsmass)) self%pm_rh50   = self%pm_rh50   + CA2G_br%brsmass
 if(associated(self%pm25_rh50) .and. associated(CA2G_br%brsmass)) self%pm25_rh50 = self%pm25_rh50 + CA2G_br%brsmass

 if(associated(self%totangstr) .and. associated(CA2G_br%brexttau) .and. associated(CA2G_br%brangstr)) then
    tau1 = tau1 + CA2G_br%brexttau(:,:,ind550)*exp(c1*CA2G_br%brangstr)
    tau2 = tau2 + CA2G_br%brexttau(:,:,ind550)*exp(c2*CA2G_br%brangstr)
 endif


!--- add AOPs from organic carbon to total AOPs:
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(CA2G_oc%ocexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + CA2G_oc%ocexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(CA2G_oc%ocstexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + CA2G_oc%ocstexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(CA2G_oc%ocscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + CA2G_oc%ocscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(CA2G_oc%ocstscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + CA2G_oc%ocstscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(CA2G_oc%ocexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + CA2G_oc%ocexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(CA2G_oc%ocscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + CA2G_oc%ocscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(CA2G_oc%ocexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + CA2G_oc%ocexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(CA2G_oc%ocscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + CA2G_oc%ocscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(CA2G_oc%ocextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + CA2G_oc%ocextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(CA2G_oc%ocextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + CA2G_oc%ocextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(CA2G_oc%ocextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + CA2G_oc%ocextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(CA2G_oc%ocscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + CA2G_oc%ocscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(CA2G_oc%ocscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + CA2G_oc%ocscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(CA2G_oc%ocscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + CA2G_oc%ocscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(CA2G_oc%ocbckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + CA2G_oc%ocbckcoef(:,:,:,w)
 enddo

 if(associated(self%pm)        .and. associated(CA2G_oc%ocsmass)) self%pm        = self%pm        + CA2G_oc%ocsmass
 if(associated(self%pm25)      .and. associated(CA2G_oc%ocsmass)) self%pm25      = self%pm25      + CA2G_oc%ocsmass
 if(associated(self%pm_rh35)   .and. associated(CA2G_oc%ocsmass)) self%pm_rh35   = self%pm_rh35   + CA2G_oc%ocsmass
 if(associated(self%pm25_rh35) .and. associated(CA2G_oc%ocsmass)) self%pm25_rh35 = self%pm25_rh35 + CA2G_oc%ocsmass
 if(associated(self%pm_rh50)   .and. associated(CA2G_oc%ocsmass)) self%pm_rh50   = self%pm_rh50   + CA2G_oc%ocsmass
 if(associated(self%pm25_rh50) .and. associated(CA2G_oc%ocsmass)) self%pm25_rh50 = self%pm25_rh50 + CA2G_oc%ocsmass

 if(associated(self%totangstr) .and. associated(CA2G_oc%ocexttau) .and. associated(CA2G_oc%ocangstr)) then
    tau1 = tau1 + CA2G_oc%ocexttau(:,:,ind550)*exp(c1*CA2G_oc%ocangstr)
    tau2 = tau2 + CA2G_oc%ocexttau(:,:,ind550)*exp(c2*CA2G_oc%ocangstr)
 endif


!--- add AOPs from mineral dust to total AOPs:
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(DU2G%duexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + DU2G%duexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(DU2G%dustexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + DU2G%dustexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(DU2G%duscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + DU2G%duscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(DU2G%dustscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + DU2G%dustscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(DU2G%duexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + DU2G%duexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(DU2G%duscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + DU2G%duscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(DU2G%duexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + DU2G%duexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(DU2G%duscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + DU2G%duscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(DU2G%duextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + DU2G%duextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(DU2G%duextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + DU2G%duextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(DU2G%duextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + DU2G%duextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(DU2G%duscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + DU2G%duscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(DU2G%duscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + DU2G%duscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(DU2G%duscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + DU2G%duscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(DU2G%dubckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + DU2G%dubckcoef(:,:,:,w)
 enddo

 if(associated(self%pm)        .and. associated(DU2G%dusmass)) self%pm        = self%pm        + DU2G%dusmass
 if(associated(self%pm25)      .and. associated(DU2G%dusmass)) self%pm25      = self%pm25      + DU2G%dusmass25
 if(associated(self%pm_rh35)   .and. associated(DU2G%dusmass)) self%pm_rh35   = self%pm_rh35   + DU2G%dusmass
 if(associated(self%pm25_rh35) .and. associated(DU2G%dusmass)) self%pm25_rh35 = self%pm25_rh35 + DU2G%dusmass25
 if(associated(self%pm_rh50)   .and. associated(DU2G%dusmass)) self%pm_rh50   = self%pm_rh50   + DU2G%dusmass
 if(associated(self%pm25_rh50) .and. associated(DU2G%dusmass)) self%pm25_rh50 = self%pm25_rh50 + DU2G%dusmass25

 if(associated(self%totangstr) .and. associated(DU2G%duexttau) .and. associated(DU2G%duangstr)) then
    tau1 = tau1 + DU2G%duexttau(:,:,ind550)*exp(c1*DU2G%duangstr)
    tau2 = tau2 + DU2G%duexttau(:,:,ind550)*exp(c2*DU2G%duangstr)
 endif


!--- add AOPs from nitrates to total AOPs (note that nitrates only support one active substance):
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(NI2G%niexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + NI2G%niexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(NI2G%nistexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + NI2G%nistexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(NI2G%niscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + NI2G%niscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(NI2G%nistscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + NI2G%nistscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(NI2G%niexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + NI2G%niexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(NI2G%niscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + NI2G%niscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(NI2G%niexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + NI2G%niexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(NI2G%niscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + NI2G%niscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(NI2G%niextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + NI2G%niextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(NI2G%niextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + NI2G%niextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(NI2G%niextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + NI2G%niextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(NI2G%niscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + NI2G%niscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(NI2G%niscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + NI2G%niscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(NI2G%niscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + NI2G%niscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(NI2G%nibckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + NI2G%nibckcoef(:,:,:,w)
 enddo

 if(associated(self%pm) .and. associated(NI2G%nismass) .and. associated(NI2G%nh4smass)) &
    self%pm = self%pm + NI2G%nismass + NI2G%nh4smass
 if(associated(self%pm25) .and. associated(NI2G%nismass25) .and. associated(NI2G%nh4smass)) &
    self%pm25 = self%pm25 + NI2G%nismass + NI2G%nh4smass
 if(associated(self%pm_rh35) .and. associated(NI2G%nismass) .and. associated(NI2G%nh4smass)) &
    self%pm_rh35   = self%pm_rh35 + 1.33*(NI2G%nismass + NI2G%nh4smass)
 if(associated(self%pm25_rh35) .and. associated(NI2G%nismass25) .and. associated(NI2G%nh4smass)) &
    self%pm25_rh35 = self%pm25_rh35 + 1.33*(NI2G%nismass25 + NI2G%nh4smass)
 if(associated(self%pm_rh50) .and. associated(NI2G%nismass) .and. associated(NI2G%nh4smass)) &
    self%pm_rh50   = self%pm_rh50   + 1.51*(NI2G%nismass + NI2G%nh4smass)
 if(associated(self%pm25_rh50) .and. associated(NI2G%nismass25) .and. associated(NI2G%nh4smass)) &
    self%pm25_rh50 = self%pm25_rh50 + 1.51*(NI2G%nismass25 + NI2G%nh4smass)

 if(associated(self%totangstr) .and. associated(NI2G%niexttau) .and. associated(NI2G%niangstr)) then
    tau1 = tau1 + NI2G%niexttau(:,:,ind550)*exp(c1*NI2G%niangstr)
    tau2 = tau2 + NI2G%niexttau(:,:,ind550)*exp(c2*NI2G%niangstr)
 endif


!--- add AOPs from sea-salt to total AOPs:
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(SS2G%ssexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + SS2G%ssexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(SS2G%ssstexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + SS2G%ssstexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(SS2G%ssscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + SS2G%ssscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(SS2G%ssstscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + SS2G%ssstscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(SS2G%ssexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + SS2G%ssexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(SS2G%ssscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + SS2G%ssscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(SS2G%ssexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + SS2G%ssexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(SS2G%ssscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + SS2G%ssscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(SS2G%ssextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + SS2G%ssextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(SS2G%ssextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + SS2G%ssextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(SS2G%ssextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + SS2G%ssextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(SS2G%ssscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + SS2G%ssscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(SS2G%ssscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + SS2G%ssscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(SS2G%ssscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + SS2G%ssscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(SS2G%ssbckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + SS2G%ssbckcoef(:,:,:,w)
 enddo

 if(associated(self%pm)        .and. associated(SS2G%sssmass)) self%pm        = self%pm        + SS2G%sssmass
 if(associated(self%pm25)      .and. associated(SS2G%sssmass)) self%pm25      = self%pm25      + SS2G%sssmass25
 if(associated(self%pm_rh35)   .and. associated(SS2G%sssmass)) self%pm_rh35   = self%pm_rh35   + 1.86*SS2G%sssmass
 if(associated(self%pm25_rh35) .and. associated(SS2G%sssmass)) self%pm25_rh35 = self%pm25_rh35 + 1.86*SS2G%sssmass25
 if(associated(self%pm_rh50)   .and. associated(SS2G%sssmass)) self%pm_rh50   = self%pm_rh50   + 2.42*SS2G%sssmass
 if(associated(self%pm25_rh50) .and. associated(SS2G%sssmass)) self%pm25_rh50 = self%pm25_rh50 + 2.42*SS2G%sssmass25

 if(associated(self%totangstr) .and. associated(SS2G%ssexttau) .and. associated(SS2G%ssangstr)) then
    tau1 = tau1 + SS2G%ssexttau(:,:,ind550)*exp(c1*SS2G%ssangstr)
    tau2 = tau2 + SS2G%ssexttau(:,:,ind550)*exp(c2*SS2G%ssangstr)
 endif


!--- add AOPs from sulfates to total AOPs:
 do w = 1,nw_vertint
    if(associated(self%totexttau) .and. associated(SU2G%suexttau)) &
       self%totexttau(:,:,w) = self%totexttau(:,:,w) + SU2G%suexttau(:,:,w)
    if(associated(self%totstexttau) .and. associated(SU2G%sustexttau)) &
       self%totstexttau(:,:,w) = self%totstexttau(:,:,w) + SU2G%sustexttau(:,:,w)
    if(associated(self%totscatau) .and. associated(SU2G%suscatau)) &
       self%totscatau(:,:,w) = self%totscatau(:,:,w) + SU2G%suscatau(:,:,w)
    if(associated(self%totstscatau) .and. associated(SU2G%sustscatau)) &
       self%totstscatau(:,:,w) = self%totstscatau(:,:,w) + SU2G%sustscatau(:,:,w)

    if(associated(self%totextt25) .and. associated(SU2G%suexttau)) &
       self%totextt25(:,:,w) = self%totextt25(:,:,w) + SU2G%suexttau(:,:,w)
    if(associated(self%totscat25) .and. associated(SU2G%suscatau)) &
       self%totscat25(:,:,w) = self%totscat25(:,:,w) + SU2G%suscatau(:,:,w)
    if(associated(self%totexttfm) .and. associated(SU2G%suexttau)) &
       self%totexttfm(:,:,w) = self%totexttfm(:,:,w) + SU2G%suexttau(:,:,w)
    if(associated(self%totscatfm) .and. associated(SU2G%suscatau)) &
       self%totscatfm(:,:,w) = self%totscatfm(:,:,w) + SU2G%suscatau(:,:,w)
 enddo

 do w = 1,nw_profile
    if(associated(self%totextcoef) .and. associated(SU2G%suextcoef)) &
       self%totextcoef(:,:,:,w) = self%totextcoef(:,:,:,w) + SU2G%suextcoef(:,:,:,w)
    if(associated(self%totextcoefrh20) .and. associated(SU2G%suextcoefrh20)) &
       self%totextcoefrh20(:,:,:,w) = self%totextcoefrh20(:,:,:,w) + SU2G%suextcoefrh20(:,:,:,w)
    if(associated(self%totextcoefrh80) .and. associated(SU2G%suextcoefrh80)) &
       self%totextcoefrh80(:,:,:,w) = self%totextcoefrh80(:,:,:,w) + SU2G%suextcoefrh80(:,:,:,w)
    if(associated(self%totscacoef) .and. associated(SU2G%suscacoef)) &
       self%totscacoef(:,:,:,w) = self%totscacoef(:,:,:,w) + SU2G%suscacoef(:,:,:,w)
    if(associated(self%totscacoefrh20) .and. associated(SU2G%suscacoefrh20)) &
       self%totscacoefrh20(:,:,:,w) = self%totscacoefrh20(:,:,:,w) + SU2G%suscacoefrh20(:,:,:,w)
    if(associated(self%totscacoefrh80) .and. associated(SU2G%suscacoefrh80)) &
       self%totscacoefrh80(:,:,:,w) = self%totscacoefrh80(:,:,:,w) + SU2G%suscacoefrh80(:,:,:,w)
    if(associated(self%totbckcoef) .and. associated(SU2G%subckcoef)) &
       self%totbckcoef(:,:,:,w) = self%totbckcoef(:,:,:,w) + SU2G%subckcoef(:,:,:,w)
 enddo

 if(associated(self%pm)        .and. associated(SU2G%so4smass)) self%pm        = self%pm        + 1.38*SU2G%so4smass
 if(associated(self%pm25)      .and. associated(SU2G%so4smass)) self%pm25      = self%pm25      + 1.38*SU2G%so4smass
 if(associated(self%pm_rh35)   .and. associated(SU2G%so4smass)) self%pm_rh35   = self%pm_rh35   + 1.33*SU2G%so4smass
 if(associated(self%pm25_rh35) .and. associated(SU2G%so4smass)) self%pm25_rh35 = self%pm25_rh35 + 1.33*SU2G%so4smass
 if(associated(self%pm_rh50)   .and. associated(SU2G%so4smass)) self%pm_rh50   = self%pm_rh50   + 1.33*SU2G%so4smass
 if(associated(self%pm25_rh50) .and. associated(SU2G%so4smass)) self%pm25_rh50 = self%pm25_rh50 + 1.33*SU2G%so4smass

 if(associated(self%totangstr) .and. associated(SU2G%suexttau) .and. associated(SU2G%suangstr)) then
    tau1 = tau1 + SU2G%suexttau(:,:,ind550)*exp(c1*SU2G%suangstr)
    tau2 = tau2 + SU2G%suexttau(:,:,ind550)*exp(c2*SU2G%suangstr)
 endif


!--- finish to compute the total Angstrom coefficient:
 if(associated(self%totangstr)) self%totangstr = log(tau1/tau2)/c3
 

!--- deallocation of wavelengths_profile and wavelengths_vertint:
 if(allocated(self_params%wavelengths_profile)) deallocate(self_params%wavelengths_profile) 
 if(allocated(self_params%wavelengths_vertint)) deallocate(self_params%wavelengths_vertint)


 call mpas_log_write('--- end subroutine processes_GOCART2G_GridComp.')

 end subroutine processes_GOCART2G_GridComp

!==================================================================================================================
 end module GOCART2G_GridCompMod
!==================================================================================================================
