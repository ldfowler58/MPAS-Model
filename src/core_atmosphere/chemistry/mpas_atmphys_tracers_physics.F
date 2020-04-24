!=================================================================================================================
 module mpas_atmphys_tracers_physics
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_vars,only: ims,ime,jms,jme,kms,kme, &
                             its,ite,jts,jte,kts,kte

 implicit none
 private


 type,public:: tracers_physics
    integer:: ntracers ! number of passive tracers.
    real(kind=RKIND),dimension(:,:,:,:),pointer:: tracers_p       => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: rtracersblten_p => null()
    real(kind=RKIND),dimension(:,:,:,:),pointer:: rtracerscuten_p => null()

    contains
       procedure:: allocate_mem         => mpas_tracers_allocate
       procedure:: deallocate_mem       => mpas_tracers_deallocate
       procedure:: tracers_to_physics   => mpas_tracers_to_physics
       procedure:: tracers_from_physics => mpas_tracers_from_physics
 end type


 contains


!=================================================================================================================
 subroutine mpas_tracers_allocate(this,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: state

!inout arguments:
 class(tracers_physics),intent(inout):: this

!local pointers:
 integer,pointer:: num_tracers

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_tracers_allocate:')
 call mpas_log_write('ims = $i   ime = $i',intArgs=(/ims,ime/))
 call mpas_log_write('jms = $i   jme = $i',intArgs=(/jms,jme/))
 call mpas_log_write('kms = $i   kme = $i',intArgs=(/kms,kme/))
 call mpas_log_write('its = $i   ite = $i',intArgs=(/its,ite/))
 call mpas_log_write('jts = $i   jte = $i',intArgs=(/jts,jte/))
 call mpas_log_write('kts = $i   kte = $i',intArgs=(/kts,kte/))

 call mpas_pool_get_dimension(state,'num_tracers',num_tracers)
 this%ntracers = num_tracers
!call mpas_log_write(' ')
 call mpas_log_write('ntracers   = $i',intArgs=(/this%ntracers/))

 if(.not.associated(this%tracers_p)      ) allocate(this%tracers_p(ims:ime,kms:kme,jms:jme,this%ntracers)      )
 if(.not.associated(this%rtracersblten_p)) allocate(this%rtracersblten_p(ims:ime,kms:kme,jms:jme,this%ntracers))
 if(.not.associated(this%rtracerscuten_p)) allocate(this%rtracerscuten_p(ims:ime,kms:kme,jms:jme,this%ntracers))

 call mpas_log_write('--- end subroutine mpas_tracers_allocate:')

 end subroutine mpas_tracers_allocate

!=================================================================================================================
 subroutine mpas_tracers_deallocate(this)
!=================================================================================================================

 class(tracers_physics),intent(inout):: this

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine mpas_tracers_deallocate:')

 if(associated(this%tracers_p)      ) deallocate(this%tracers_p      )
 if(associated(this%rtracersblten_p)) deallocate(this%rtracersblten_p)
 if(associated(this%rtracerscuten_p)) deallocate(this%rtracerscuten_p)

 call mpas_log_write('--- end subroutine mpas_tracers_deallocate:')
 call mpas_log_write(' ')

 end subroutine mpas_tracers_deallocate

!=================================================================================================================
 subroutine mpas_tracers_to_physics(this,state,time_lev)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: state
 integer,intent(in):: time_lev

!inout arguments:
 class(tracers_physics),intent(inout):: this

!local variables and arrays:
 integer:: i,k,j,n

!local pointers associated to this:
 integer:: ntracers

 real(kind=RKIND),dimension(:,:,:),pointer  :: tracers
 real(kind=RKIND),dimension(:,:,:,:),pointer:: tracers_p
 real(kind=RKIND),dimension(:,:,:,:),pointer:: rtracersblten_p,rtracerscuten_p

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_tracers_to_physics:')

 call mpas_pool_get_array(state,'tracers',tracers,time_lev)

 ntracers = this%ntracers
 tracers_p       => this%tracers_p
 rtracersblten_p => this%rtracersblten_p
 rtracerscuten_p => this%rtracerscuten_p

!call mpas_log_write('--- ntracers = $i',intArgs=(/ntracers/))!

 do n = 1,ntracers
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             tracers_p(i,k,j,n) = tracers(n,k,i)
             rtracersblten_p(i,k,j,n) = 0._RKIND
             rtracerscuten_p(i,k,j,n) = 0._RKIND
          enddo
       enddo
    enddo
 enddo

 call mpas_log_write('--- enter subroutine mpas_tracers_to_physics:')
!call mpas_log_write(' ')

 end subroutine mpas_tracers_to_physics

!=================================================================================================================
 subroutine mpas_tracers_from_physics(this,tend_tracers)
!=================================================================================================================

!input arguments:
 class(tracers_physics),intent(in):: this

!output arguments:
 type(mpas_pool_type),intent(inout):: tend_tracers

!local variables and arrays:
 integer:: i,k,j,n

!local pointers associated to this:
 integer:: ntracers
 real(kind=RKIND),dimension(:,:,:,:),pointer:: rtracersblten_p,rtracerscuten_p

 integer,pointer:: num_bl_tracers
 integer,pointer:: num_cu_tracers
 real(kind=RKIND),dimension(:,:,:),pointer:: bl_tracers,cu_tracers

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine mpas_tracers_from_physics:')

 call mpas_pool_get_dimension(tend_tracers,'num_bl_tracers',num_bl_tracers)
 call mpas_pool_get_dimension(tend_tracers,'num_cu_tracers',num_cu_tracers)
 call mpas_log_write('--- num_bl_tracers = $i',intArgs=(/num_bl_tracers/))
 call mpas_log_write('--- num_cu_tracers = $i',intArgs=(/num_cu_tracers/))

 call mpas_pool_get_array(tend_tracers,'bl_tracers',bl_tracers)
 call mpas_pool_get_array(tend_tracers,'cu_tracers',cu_tracers)

 ntracers  = this%ntracers
 rtracersblten_p => this%rtracersblten_p
 rtracerscuten_p => this%rtracerscuten_p

 call mpas_log_write('--- ntracers       = $i',intArgs=(/ntracers/))
 do n = 1, num_bl_tracers
    do j = jts,jte
       do k = kts,kte
          do i = its,ite
             bl_tracers(n,k,i) = rtracersblten_p(i,k,j,n)
             cu_tracers(n,k,i) = rtracerscuten_p(i,k,j,n)
          enddo
       enddo
    enddo
 enddo

!do j = jts,jte
!   do i = its,its
!      do k = kts,kte
!         call mpas_log_write('$i $i $r $r $r $r',intArgs=(/i,k/), &
!                        realArgs=(/(rtracersblten_p(i,k,j,n),n=1,ntracers), &
!                        (rtracerscuten_p(i,k,j,n),n=1,ntracers)/))
!      enddo
!      call mpas_log_write(' ')
!   enddo
!enddo
!call mpas_log_write(' ')
!do i = its,its
!   do k = kts,kte
!      call mpas_log_write('$i $i $r $r $r $r',intArgs=(/i,k/), &
!                realArgs=(/(bl_tracers(n,k,i),n=1,num_bl_tracers), &
!                           (cu_tracers(n,k,i),n=1,num_cu_tracers)/))
!   enddo
!   call mpas_log_write(' ')
!enddo

 call mpas_log_write('--- end subroutine mpas_tracers_to_physics:')
!call mpas_log_write(' ')

 end subroutine mpas_tracers_from_physics

!=================================================================================================================
 end module mpas_atmphys_tracers_physics
!=================================================================================================================

