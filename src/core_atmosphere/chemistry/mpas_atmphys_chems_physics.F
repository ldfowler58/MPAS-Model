!=================================================================================================================
 module mpas_atmphys_chems_physics
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_vars,only: ims,ime,jms,jme,kms,kme, &
                             its,ite,jts,jte,kts,kte

 implicit none
 private


 type,public:: chems_physics
    integer:: nchems     ! number of chemical species.
    integer:: nchems_mix !
    integer:: kdepvel    ! number of vertical levels for deposition of chemical species.
    integer:: ndepvel    ! number of chemical species with deposition.
    real(kind=RKIND),dimension(:,:,:,:),pointer:: chems_p       => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: depvels_p     => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: rchemsblten_p => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: rchemscuten_p => null()

    integer:: naers      ! number of aerosol species that describe nwfa and nifa (Thompson aerosol-aware MP).
    real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_p        => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_weight_p => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_scav_p   => null()

    real(kind=RKIND),dimension(:,:),pointer:: nifa_evapt_p => null()
    real(kind=RKIND),dimension(:,:),pointer:: nifa_scavt_p => null()
    real(kind=RKIND),dimension(:,:),pointer:: nwfa_evapt_p => null()
    real(kind=RKIND),dimension(:,:),pointer:: nwfa_scavt_p => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nifa_evap_p => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nifa_scav_p => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nwfa_evap_p => null()
    real(kind=RKIND),dimension(:,:,:),pointer:: nwfa_scav_p => null()

    contains
       procedure:: allocate_mem       => mpas_chems_allocate
       procedure:: deallocate_mem     => mpas_chems_deallocate
       procedure:: chems_to_physics   => mpas_chems_to_physics
       procedure:: chems_from_physics => mpas_chems_from_physics
 end type


 contains


!=================================================================================================================
 subroutine mpas_chems_allocate(this,mesh,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

!inout arguments:
 class(chems_physics),intent(inout):: this

!local pointers:
 integer,pointer:: num_chems
 integer,pointer:: num_aers

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chems_allocate:')
 call mpas_log_write('ims = $i   ime = $i',intArgs=(/ims,ime/))
 call mpas_log_write('jms = $i   jme = $i',intArgs=(/jms,jme/))
 call mpas_log_write('kms = $i   kme = $i',intArgs=(/kms,kme/))
 call mpas_log_write('its = $i   ite = $i',intArgs=(/its,ite/))
 call mpas_log_write('jts = $i   jte = $i',intArgs=(/jts,jte/))
 call mpas_log_write('kts = $i   kte = $i',intArgs=(/kts,kte/))

 call mpas_pool_get_dimension(state,'num_chems',num_chems)
 this%nchems     = num_chems
 this%nchems_mix = num_chems
 this%ndepvel    = num_chems
 this%kdepvel    = 1
!call mpas_log_write(' ')
 call mpas_log_write('nchems     = $i',intArgs=(/this%nchems/)    )
 call mpas_log_write('nchems_mix = $i',intArgs=(/this%nchems_mix/))
 call mpas_log_write('ndepvel    = $i',intArgs=(/this%ndepvel/)   )
 call mpas_log_write('kdepvel    = $i',intArgs=(/this%kdepvel/)   )

 if(.not.associated(this%chems_p)      ) allocate(this%chems_p(ims:ime,kms:kme,jms:jme,this%nchems)        )
 if(.not.associated(this%rchemsblten_p)) allocate(this%rchemsblten_p(ims:ime,kms:kme,jms:jme,this%nchems)  )
 if(.not.associated(this%rchemscuten_p)) allocate(this%rchemscuten_p(ims:ime,kms:kme,jms:jme,this%nchems)  )
 if(.not.associated(this%depvels_p)    ) allocate(this%depvels_p(ims:ime,this%kdepvel,jms:jme,this%ndepvel))

!--- Thompson aerosol-aware cloud microphysics scheme:
 call mpas_pool_get_dimension(mesh,'num_aerosols',num_aers)
 this%naers = num_aers
 call mpas_log_write('naers      = $i',intArgs=(/this%naers/))

 if(.not.associated(this%aerosols_p)       ) allocate(this%aerosols_p(ims:ime,kms:kme,jms:jme,this%naers)       )
 if(.not.associated(this%aerosols_weight_p)) allocate(this%aerosols_weight_p(ims:ime,kms:kme,jms:jme,this%naers))
 if(.not.associated(this%aerosols_scav_p)  ) allocate(this%aerosols_scav_p(ims:ime,kms:kme,jms:jme,this%naers)  )

 if(.not.associated(this%nifa_evapt_p)) allocate(this%nifa_evapt_p(ims:ime,jms:jme))
 if(.not.associated(this%nifa_scavt_p)) allocate(this%nifa_scavt_p(ims:ime,jms:jme))
 if(.not.associated(this%nwfa_evapt_p)) allocate(this%nwfa_evapt_p(ims:ime,jms:jme))
 if(.not.associated(this%nwfa_scavt_p)) allocate(this%nwfa_scavt_p(ims:ime,jms:jme))
 if(.not.associated(this%nifa_evap_p) ) allocate(this%nifa_evap_p(ims:ime,kms:kme,jms:jme))
 if(.not.associated(this%nifa_scav_p) ) allocate(this%nifa_scav_p(ims:ime,kms:kme,jms:jme))
 if(.not.associated(this%nwfa_evap_p) ) allocate(this%nwfa_evap_p(ims:ime,kms:kme,jms:jme))
 if(.not.associated(this%nwfa_scav_p) ) allocate(this%nwfa_scav_p(ims:ime,kms:kme,jms:jme))

 call mpas_log_write('--- end subroutine mpas_chems_allocate:')

 end subroutine mpas_chems_allocate

!=================================================================================================================
 subroutine mpas_chems_deallocate(this)
!=================================================================================================================

 class(chems_physics),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine mpas_chems_deallocate:')

 if(associated(this%chems_p)      ) deallocate(this%chems_p      )
 if(associated(this%rchemsblten_p)) deallocate(this%rchemsblten_p)
 if(associated(this%rchemscuten_p)) deallocate(this%rchemscuten_p)
 if(associated(this%depvels_p)    ) deallocate(this%depvels_p    )

 if(associated(this%aerosols_p)       ) deallocate(this%aerosols_p       )
 if(associated(this%aerosols_weight_p)) deallocate(this%aerosols_weight_p)
 if(associated(this%aerosols_scav_p)  ) deallocate(this%aerosols_scav_p  )

 if(associated(this%nifa_evapt_p)) deallocate(this%nifa_evapt_p)
 if(associated(this%nifa_scavt_p)) deallocate(this%nifa_scavt_p)
 if(associated(this%nwfa_evapt_p)) deallocate(this%nwfa_evapt_p)
 if(associated(this%nwfa_scavt_p)) deallocate(this%nwfa_scavt_p)
 if(associated(this%nifa_evap_p) ) deallocate(this%nifa_evap_p )
 if(associated(this%nifa_scav_p) ) deallocate(this%nifa_scav_p )
 if(associated(this%nwfa_evap_p) ) deallocate(this%nwfa_evap_p )
 if(associated(this%nwfa_scav_p) ) deallocate(this%nwfa_scav_p )

 call mpas_log_write('--- end subroutine mpas_chems_deallocate:')
 call mpas_log_write(' ')

 end subroutine mpas_chems_deallocate

!=================================================================================================================
 subroutine mpas_chems_to_physics(this,mesh,state,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state
 integer,intent(in):: time_lev

!inout arguments:
 class(chems_physics),intent(inout):: this

!local variables and arrays:
 integer:: i,k,j,n

!local pointers associated to this:
 integer:: nchems,kdepvel,ndepvel
 integer:: naers

 real(kind=RKIND),dimension(:,:,:),pointer  :: chems
 real(kind=RKIND),dimension(:,:,:,:),pointer:: chems_p
 real(kind=RKIND),dimension(:,:,:,:),pointer:: depvels_p
 real(kind=RKIND),dimension(:,:,:,:),pointer:: rchemsblten_p,rchemscuten_p

 integer,pointer:: index_dust1,index_dust2,index_dust3,index_dust4,index_dust5
 integer,pointer:: index_oc,index_so4,index_seas1,index_seas2
 real(kind=RKIND),dimension(:,:,:),pointer  :: aerosols,aerosols_weight
 real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_p,aerosols_weight_p,aerosols_scav_p

 real(kind=RKIND),dimension(:,:),pointer  :: nifa_evapt_p,nifa_scavt_p,nwfa_evapt_p,nwfa_scavt_p
 real(kind=RKIND),dimension(:,:,:),pointer:: nifa_evap_p,nifa_scav_p,nwfa_evap_p,nwfa_scav_p

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chems_to_physics:')

 call mpas_pool_get_array(state,'chems',chems,time_lev)

 nchems  = this%nchems
 kdepvel = this%kdepvel
 ndepvel = this%ndepvel
 chems_p       => this%chems_p
 depvels_p     => this%depvels_p
 rchemsblten_p => this%rchemsblten_p
 rchemscuten_p => this%rchemscuten_p

!call mpas_log_write('--- nchems = $i',intArgs=(/nchems/))!

 do n = 1,nchems
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             chems_p(i,k,j,n) = chems(n,k,i)
             rchemsblten_p(i,k,j,n) = 0._RKIND
             rchemscuten_p(i,k,j,n) = 0._RKIND
          enddo
       enddo
    enddo
 enddo
 do n = 1,ndepvel
    do j = jts,jte
       do k = 1,kdepvel
          do i = its,ite
             depvels_p(i,k,j,n) = 0._RKIND
          enddo
       enddo
    enddo
 enddo

!initialize the aerosols for activation and wet scavenging of the Thompson aerosol-aware cloud microphysics:
 naers = this%naers

 call mpas_pool_get_dimension(mesh,'index_dust1',index_dust1)
 call mpas_pool_get_dimension(mesh,'index_dust2',index_dust2)
 call mpas_pool_get_dimension(mesh,'index_dust3',index_dust3)
 call mpas_pool_get_dimension(mesh,'index_dust4',index_dust4)
 call mpas_pool_get_dimension(mesh,'index_dust5',index_dust5)
 call mpas_pool_get_dimension(mesh,'index_oc'   ,index_oc   )
 call mpas_pool_get_dimension(mesh,'index_so4'  ,index_so4  )
 call mpas_pool_get_dimension(mesh,'index_seas1',index_seas1)
 call mpas_pool_get_dimension(mesh,'index_seas2',index_seas2)

 call mpas_pool_get_array(mesh,'aerosols',aerosols)
 call mpas_pool_get_array(mesh,'aerosols_weight',aerosols_weight)

 call mpas_log_write('--- naers       = $i',intArgs=(/naers/))
 call mpas_log_write('--- index_dust1 = $i',intArgs=(/index_dust1/))
 call mpas_log_write('--- index_dust2 = $i',intArgs=(/index_dust2/))
 call mpas_log_write('--- index_dust3 = $i',intArgs=(/index_dust3/))
 call mpas_log_write('--- index_dust4 = $i',intArgs=(/index_dust4/))
 call mpas_log_write('--- index_dust5 = $i',intArgs=(/index_dust5/))
 call mpas_log_write('--- index_oc    = $i',intArgs=(/index_oc/)   )
 call mpas_log_write('--- index_so4   = $i',intArgs=(/index_so4/)  )
 call mpas_log_write('--- index_seas1 = $i',intArgs=(/index_seas1/))
 call mpas_log_write('--- index_seas2 = $i',intArgs=(/index_seas2/))

 aerosols_p        => this%aerosols_p
 aerosols_weight_p => this%aerosols_weight_p
 aerosols_scav_p   => this%aerosols_scav_p
 nifa_evapt_p => this%nifa_evapt_p
 nifa_scavt_p => this%nifa_scavt_p
 nwfa_evapt_p => this%nwfa_evapt_p
 nwfa_scavt_p => this%nwfa_scavt_p
 nifa_evap_p  => this%nifa_evap_p
 nifa_scav_p  => this%nifa_scav_p
 nwfa_evap_p  => this%nwfa_evap_p
 nwfa_scav_p  => this%nwfa_scav_p

 do n = 1,naers
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             aerosols_p(i,k,j,n) = aerosols(n,k,i)
             aerosols_weight_p(i,k,j,n) = aerosols_weight(n,k,i)
             aerosols_scav_p(i,k,j,n) = 0._RKIND
          enddo
       enddo
    enddo
 enddo

 do j = jts,jte
    do i = its,ite
       nifa_evapt_p(i,j) = 0._RKIND
       nifa_scavt_p(i,j) = 0._RKIND
       nwfa_evapt_p(i,j) = 0._RKIND
       nwfa_scavt_p(i,j) = 0._RKIND
    enddo
    do k = kts,kte
       do i = its,ite
          nifa_evap_p(i,k,j) = 0._RKIND
          nifa_scav_p(i,k,j) = 0._RKIND
          nwfa_evap_p(i,k,j) = 0._RKIND
          nwfa_scav_p(i,k,j) = 0._RKIND
       enddo
    enddo
 enddo

 call mpas_log_write('--- end subroutine mpas_chems_to_physics:')
!call mpas_log_write(' ')

 end subroutine mpas_chems_to_physics

!=================================================================================================================
 subroutine mpas_chems_from_physics(this,tend_chems,diag_chems)
!=================================================================================================================

!input arguments:
 class(chems_physics),intent(in):: this

!output arguments:
 type(mpas_pool_type),intent(inout):: tend_chems
 type(mpas_pool_type),intent(inout):: diag_chems

!local variables and arrays:
 integer:: i,k,j,n

!local pointers associated to this:
 integer:: nchems
 integer:: naers
 real(kind=RKIND),dimension(:,:,:,:),pointer:: rchemsblten_p,rchemscuten_p
 real(kind=RKIND),dimension(:,:),pointer    :: nifa_evapt_p,nifa_scavt_p,nwfa_evapt_p,nwfa_scavt_p
 real(kind=RKIND),dimension(:,:,:),pointer  :: nifa_evap_p,nifa_scav_p,nwfa_evap_p,nwfa_scav_p
 real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_scav_p

 integer,pointer:: num_bl_chems
 integer,pointer:: num_cu_chems
 real(kind=RKIND),dimension(:,:,:),pointer:: bl_chems,cu_chems

 real(kind=RKIND),dimension(:),pointer  :: nifa_evapt,nifa_scavt,nwfa_evapt,nwfa_scavt
 real(kind=RKIND),dimension(:,:),pointer:: nifa_evap,nifa_scav,nwfa_evap,nwfa_scav
 real(kind=RKIND),dimension(:,:,:),pointer:: aerosols_scav

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chems_from_physics:')

 call mpas_pool_get_dimension(tend_chems,'num_bl_chems',num_bl_chems)
 call mpas_pool_get_dimension(tend_chems,'num_cu_chems',num_cu_chems)
 call mpas_log_write('--- num_bl_chems = $i',intArgs=(/num_bl_chems/))
 call mpas_log_write('--- num_cu_chems = $i',intArgs=(/num_cu_chems/))

 call mpas_pool_get_array(tend_chems,'bl_chems',bl_chems)
 call mpas_pool_get_array(tend_chems,'cu_chems',cu_chems)

 call mpas_pool_get_array(diag_chems,'aerosols_scav',aerosols_scav)
 call mpas_pool_get_array(diag_chems,'nifa_evapt',nifa_evapt)
 call mpas_pool_get_array(diag_chems,'nifa_scavt',nifa_scavt)
 call mpas_pool_get_array(diag_chems,'nwfa_evapt',nwfa_evapt)
 call mpas_pool_get_array(diag_chems,'nwfa_scavt',nwfa_scavt)
 call mpas_pool_get_array(diag_chems,'nifa_evap' ,nifa_evap )
 call mpas_pool_get_array(diag_chems,'nifa_scav' ,nifa_scav )
 call mpas_pool_get_array(diag_chems,'nwfa_evap' ,nwfa_evap )
 call mpas_pool_get_array(diag_chems,'nwfa_scav' ,nwfa_scav )

 nchems = this%nchems
 rchemsblten_p => this%rchemsblten_p
 rchemscuten_p => this%rchemscuten_p

 naers = this%naers
 nifa_evapt_p => this%nifa_evapt_p
 nifa_scavt_p => this%nifa_scavt_p
 nwfa_evapt_p => this%nwfa_evapt_p
 nwfa_scavt_p => this%nwfa_scavt_p
 nifa_evap_p  => this%nifa_evap_p
 nifa_scav_p  => this%nifa_scav_p
 nwfa_evap_p  => this%nwfa_evap_p
 nwfa_scav_p  => this%nwfa_scav_p
 aerosols_scav_p => this%aerosols_scav_p

 call mpas_log_write('--- nchems       = $i',intArgs=(/nchems/))
 do n = 1,num_bl_chems
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             bl_chems(n,k,i) = rchemsblten_p(i,k,j,n)
             cu_chems(n,k,i) = rchemscuten_p(i,k,j,n)
          enddo
       enddo
    enddo
 enddo

 call mpas_log_write('--- naers        = $i',intArgs=(/naers/))
 do n = 1,naers
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             aerosols_scav(n,k,i) = aerosols_scav_p(i,k,j,n)
          enddo
       enddo
    enddo
 enddo
 do j = jts,jte
    do i = its,ite
       nifa_evapt(i) = nifa_evapt_p(i,j)
       nifa_scavt(i) = nifa_scavt_p(i,j)
       nwfa_evapt(i) = nwfa_evapt_p(i,j)
       nwfa_scavt(i) = nwfa_scavt_p(i,j)
    enddo
    do k = kts,kte
       do i = its,ite
          nifa_evap(k,i) = nifa_evap_p(i,k,j)
          nifa_scav(k,i) = nifa_scav_p(i,k,j)
          nwfa_evap(k,i) = nwfa_evap_p(i,k,j)
          nwfa_scav(k,i) = nwfa_scav_p(i,k,j)
       enddo
    enddo
 enddo

 call mpas_log_write('--- end subroutine mpas_chems_to_physics:')
!call mpas_log_write(' ')

 end subroutine mpas_chems_from_physics

!=================================================================================================================
 end module mpas_atmphys_chems_physics
!=================================================================================================================

