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

    contains
       procedure:: allocate_mem       => mpas_chems_allocate
       procedure:: deallocate_mem     => mpas_chems_deallocate
       procedure:: chems_to_physics   => mpas_chems_to_physics
       procedure:: chems_from_physics => mpas_chems_from_physics
 end type


 contains


!=================================================================================================================
 subroutine mpas_chems_allocate(this,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: state

!inout arguments:
 class(chems_physics),intent(inout):: this

!local pointers:
 integer,pointer:: num_chems

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

 call mpas_log_write('--- end subroutine mpas_chems_deallocate:')
!call mpas_log_write(' ')

 end subroutine mpas_chems_deallocate

!=================================================================================================================
 subroutine mpas_chems_to_physics(this,state,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: state
 integer,intent(in):: time_lev

!inout arguments:
 class(chems_physics),intent(inout):: this

!local variables and arrays:
 integer:: i,k,j,n

!local pointers associated to this:
 integer:: nchems,kdepvel,ndepvel

 real(kind=RKIND),dimension(:,:,:),pointer  :: chems
 real(kind=RKIND),dimension(:,:,:,:),pointer:: chems_p
 real(kind=RKIND),dimension(:,:,:,:),pointer:: depvels_p
 real(kind=RKIND),dimension(:,:,:,:),pointer:: rchemsblten_p,rchemscuten_p

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

 call mpas_log_write('--- end subroutine mpas_chems_to_physics:')
!call mpas_log_write(' ')

 end subroutine mpas_chems_to_physics

!=================================================================================================================
 subroutine mpas_chems_from_physics(this,tend_chems)
!=================================================================================================================

!input arguments:
 class(chems_physics),intent(in):: this

!output arguments:
 type(mpas_pool_type),intent(inout):: tend_chems

!local variables and arrays:
 integer:: i,k,j,n

!local pointers associated to this:
 integer:: nchems
 real(kind=RKIND),dimension(:,:,:,:),pointer:: rchemsblten_p,rchemscuten_p

 integer,pointer:: num_bl_chems
 integer,pointer:: num_cu_chems
 real(kind=RKIND),dimension(:,:,:),pointer:: bl_chems,cu_chems

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_chems_from_physics:')

 call mpas_pool_get_dimension(tend_chems,'num_bl_chems',num_bl_chems)
 call mpas_pool_get_dimension(tend_chems,'num_cu_chems',num_cu_chems)
 call mpas_log_write('--- num_bl_chems = $i',intArgs=(/num_bl_chems/))
 call mpas_log_write('--- num_cu_chems = $i',intArgs=(/num_cu_chems/))

 call mpas_pool_get_array(tend_chems,'bl_chems',bl_chems)
 call mpas_pool_get_array(tend_chems,'cu_chems',cu_chems)

 nchems  = this%nchems
 rchemsblten_p => this%rchemsblten_p
 rchemscuten_p => this%rchemscuten_p

 call mpas_log_write('--- nchems       = $i',intArgs=(/nchems/))
 do n = 1, num_bl_chems
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             bl_chems(n,k,i) = rchemsblten_p(i,k,j,n)
             cu_chems(n,k,i) = rchemscuten_p(i,k,j,n)
          enddo
       enddo
    enddo
 enddo

!do j = jts,jte
!   do i = its,its
!      do k = kts,kte
!         call mpas_log_write('$i $i $r $r $r $r $r $r',intArgs=(/i,k/), &
!                        realArgs=(/(rchemsblten_p(i,k,j,n),n=1,nchems), &
!                        (rchemscuten_p(i,k,j,n),n=1,nchems)/))
!      enddo
!      call mpas_log_write(' ')
!   enddo
!enddo
!call mpas_log_write(' ')
!do i = its,its
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r $r $r $r',intArgs=(/i,k/), &
!                realArgs=(/(bl_chems(n,k,i),n=1,num_bl_chems), &
!                           (cu_chems(n,k,i),n=1,num_cu_chems)/))
!   enddo
!   call mpas_log_write(' ')
!enddo

 call mpas_log_write('--- end subroutine mpas_chems_to_physics:')
!call mpas_log_write(' ')

 end subroutine mpas_chems_from_physics

!=================================================================================================================
 end module mpas_atmphys_chems_physics
!=================================================================================================================

