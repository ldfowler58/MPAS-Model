!==================================================================================================================
 module GOCART2G_StateSpecs
 use mpas_kind_types,only: RKIND

 use GOCART2G_instance,only: wavelengths_for_profile_aop_in_nm, &
                             wavelengths_for_vertically_integrated_aop_in_nm

 implicit none
 public

!this module is the state variable specification file for AODs. it is the same as GOCART2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp.

!schema_version: 2.0.0
!component: GOCART2G


 type GOCART2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:,:),pointer:: delp           !pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),pointer:: rh2            !rel_Hum_after_moist (1)
 real(kind=RKIND),dimension(:,:,:),pointer:: airdens        !moist air density  (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),pointer:: t              !air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),pointer:: ple            !air pressure (Pa)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),pointer:: totexttau        ! total aerosol extinction aot [550 nm]
 real(kind=RKIND),dimension(:,:,:),pointer:: totstexttau      ! total aerosol extinction aot [550 nm] stratosphere
 real(kind=RKIND),dimension(:,:,:),pointer:: totscatau        ! total aerosol scattering aot [550 nm]
 real(kind=RKIND),dimension(:,:,:),pointer:: totstscatau      ! total aerosol scattering aot [550 nm] stratosphere
 real(kind=RKIND),dimension(:,:,:),pointer:: totextt25        ! total aerosol extinction aot [550 nm] - PM2.5
 real(kind=RKIND),dimension(:,:,:),pointer:: totscat25        ! total aerosol extinction aot [550 nm] - PM2.5
 real(kind=RKIND),dimension(:,:,:),pointer:: totexttfm        ! total aerosol extinction aot [550 nm] - PM1.0
 real(kind=RKIND),dimension(:,:,:),pointer:: totscatfm        ! total aerosol extinction aot [550 nm] - PM1.0
 real(kind=RKIND),dimension(:,:),pointer:: totangstr          ! total aerosol angstrom parameter [470-870 nm]
 real(kind=RKIND),dimension(:,:),pointer:: pm                 ! total reconstructed PM (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer:: pm_rh35            ! total reconstructed PM (RH=35%) (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer:: pm_rh50            ! total reconstructed PM (RH=50%) (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer:: pm25               ! total reconstructed PM2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer:: pm25_rh35          ! total reconstructed PM2.5(RH=35%) (kg m-3)
 real(kind=RKIND),dimension(:,:),pointer:: pm25_rh50          ! total reconstructed PM2.5(RH=50%) (kg m-3)
 real(kind=RKIND),dimension(:,:,:),pointer:: pso4tot          ! total sulfate produced in gocart (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totextcoef     ! total aerosol extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totextcoefrh20 ! total aerosol extinction coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totextcoefrh80 ! total aerosol extinction coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totscacoef     ! total aerosol scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totscacoefrh20 ! total aerosol scattering coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totscacoefrh80 ! total aerosol scattering coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),pointer:: totbckcoef     ! total aerosol single scattering backscatter coefficient (m-1 sr-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: totabcktoa       ! total attenuated backscatter coefficient from toa [532nm] (m-1 sr-1)
 real(kind=RKIND),dimension(:,:,:),pointer:: totabcksfc       ! total attenuated backscatter coefficient from surface [532nm] (m-1 sr-1)

!category: INTERNAL

 contains
    procedure:: gocart2G_allocate   => GOCART2G_StateSpecsInit
    procedure:: gocart2G_deallocate => GOCART2G_StateSpecsFinalize

 end type GOCART2G_State


 contains


!==================================================================================================================
 subroutine GOCART2G_StateSpecsInit(self,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(GOCART2G_State),intent(inout):: self

!--- local variables:
 integer:: nw_profile,nw_vertint

!------------------------------------------------------------------------------------------------------------------

 nw_profile = size(wavelengths_for_profile_aop_in_nm)
 nw_vertint = size(wavelengths_for_vertically_integrated_aop_in_nm)

!category: IMPORT
 if(.not.associated(self%delp)   ) allocate(self%delp(its:ite,jts:jte,kts:kte)   )
 if(.not.associated(self%rh2)    ) allocate(self%rh2(its:ite,jts:jte,kts:kte)    )
 if(.not.associated(self%airdens)) allocate(self%airdens(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%t)      ) allocate(self%t(its:ite,jts:jte,kts:kte)      )
 if(.not.associated(self%ple)    ) allocate(self%ple(its:ite,jts:jte,kts:kte+1)  )

!category: EXPORT
 if(.not.associated(self%totexttau)     ) allocate(self%totexttau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%totstexttau)   ) allocate(self%totstexttau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%totscatau)     ) allocate(self%totscatau(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%totstscatau)   ) allocate(self%totstscatau(its:ite,jts:jte,nw_vertint))
 if(.not.associated(self%totextt25)     ) allocate(self%totextt25(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%totscat25)     ) allocate(self%totscat25(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%totexttfm)     ) allocate(self%totexttfm(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%totscatfm)     ) allocate(self%totscatfm(its:ite,jts:jte,nw_vertint)  )
 if(.not.associated(self%totangstr)     ) allocate(self%totangstr(its:ite,jts:jte)      )
 if(.not.associated(self%pm)            ) allocate(self%pm(its:ite,jts:jte)             )
 if(.not.associated(self%pm_rh35)       ) allocate(self%pm_rh35(its:ite,jts:jte)        )
 if(.not.associated(self%pm_rh50)       ) allocate(self%pm_rh50(its:ite,jts:jte)        )
 if(.not.associated(self%pm25)          ) allocate(self%pm25(its:ite,jts:jte)           )
 if(.not.associated(self%pm25_rh35)     ) allocate(self%pm25_rh35(its:ite,jts:jte)      )
 if(.not.associated(self%pm25_rh50)     ) allocate(self%pm25_rh50(its:ite,jts:jte)      )
 if(.not.associated(self%pso4tot)       ) allocate(self%pso4tot(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%totextcoef)    ) allocate(self%totextcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%totextcoefrh20)) allocate(self%totextcoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%totextcoefrh80)) allocate(self%totextcoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%totscacoef)    ) allocate(self%totscacoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%totscacoefrh20)) allocate(self%totscacoefrh20(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%totscacoefrh80)) allocate(self%totscacoefrh80(its:ite,jts:jte,kts:kte,nw_profile))
 if(.not.associated(self%totbckcoef)    ) allocate(self%totbckcoef(its:ite,jts:jte,kts:kte,nw_profile)    )
 if(.not.associated(self%totabcktoa)    ) allocate(self%totabcktoa(its:ite,jts:jte,kts:kte))
 if(.not.associated(self%totabcksfc)    ) allocate(self%totabcksfc(its:ite,jts:jte,kts:kte))

!category: INTERNAL

 end subroutine GOCART2G_StateSpecsInit

!==================================================================================================================
 subroutine GOCART2G_StateSpecsFinalize(self)
!==================================================================================================================

!--- inout arguments:
 class(GOCART2G_State),intent(inout) :: self

!------------------------------------------------------------------------------------------------------------------

!category: IMPORT
 if(associated(self%delp)   ) deallocate(self%delp   )
 if(associated(self%rh2)    ) deallocate(self%rh2    )
 if(associated(self%airdens)) deallocate(self%airdens)
 if(associated(self%t)      ) deallocate(self%t      )
 if(associated(self%ple)    ) deallocate(self%ple    )

!category: EXPORT
 if(associated(self%totexttau)     ) deallocate(self%totexttau     )
 if(associated(self%totstexttau)   ) deallocate(self%totstexttau   )
 if(associated(self%totscatau)     ) deallocate(self%totscatau     )
 if(associated(self%totstscatau)   ) deallocate(self%totstscatau   )
 if(associated(self%totextt25)     ) deallocate(self%totextt25     )
 if(associated(self%totscat25)     ) deallocate(self%totscat25     )
 if(associated(self%totexttfm)     ) deallocate(self%totexttfm     )
 if(associated(self%totscatfm)     ) deallocate(self%totscatfm     )
 if(associated(self%totangstr)     ) deallocate(self%totangstr     )
 if(associated(self%pm)            ) deallocate(self%pm            )
 if(associated(self%pm_rh35)       ) deallocate(self%pm_rh35       )
 if(associated(self%pm_rh50)       ) deallocate(self%pm_rh50       )
 if(associated(self%pm25)          ) deallocate(self%pm25          )
 if(associated(self%pm25_rh35)     ) deallocate(self%pm25_rh35     )
 if(associated(self%pm25_rh50)     ) deallocate(self%pm25_rh50     )
 if(associated(self%pso4tot)       ) deallocate(self%pso4tot       )
 if(associated(self%totextcoef)    ) deallocate(self%totextcoef    )
 if(associated(self%totextcoefrh20)) deallocate(self%totextcoefrh20)
 if(associated(self%totextcoefrh80)) deallocate(self%totextcoefrh80)
 if(associated(self%totscacoef)    ) deallocate(self%totscacoef    )
 if(associated(self%totscacoefrh20)) deallocate(self%totscacoefrh20)
 if(associated(self%totscacoefrh80)) deallocate(self%totscacoefrh80)
 if(associated(self%totbckcoef)    ) deallocate(self%totbckcoef    )
 if(associated(self%totabcktoa)    ) deallocate(self%totabcktoa    )
 if(associated(self%totabcksfc)    ) deallocate(self%totabcksfc    )

!category: INTERNAL

 end subroutine GOCART2G_StateSpecsFinalize

!==================================================================================================================
 end module GOCART2G_StateSpecs
!==================================================================================================================
