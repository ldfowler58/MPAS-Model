!=================================================================================================================
 module GOCART2G_AeroGeneric
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_log

 implicit none
 private
 public:: findKlid,     &
          setZeroKlid,  &
          setZeroKlid4d


 contains


!=================================================================================================================
!BOP
! !IROUTINE: setZeroKlid
   subroutine setZeroKlid(km, klid, int_ptr)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km   ! total model levels
   integer, intent(in) :: klid ! index for pressure level

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: int_ptr ! aerosol pointer

! !DESCRIPTION: Set values to 0 where above klid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k

!EOP
!----------------------------------------------------------------------------------
!  Begin...

    do k = 1, km
       if (k < klid) then
          int_ptr(:,:,k) = 0.0
       else if (k >= klid) then
          exit
       end if
    end do

   end subroutine setZeroKlid
!===================================================================================
!BOP
! !IROUTINE: setZeroKlid
   subroutine setZeroKlid4d (km, klid, int_ptr)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km   ! total model levels
   integer, intent(in) :: klid ! index for pressure level

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:,:), intent(inout) :: int_ptr ! aerosol pointer

! !DESCRIPTION: Set values to 0 where above klid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k, n

!EOP
!----------------------------------------------------------------------------------
!  Begin...

   do n = 1, ubound(int_ptr, 4)
      do k = 1, km
         if (k < klid) then
            int_ptr(:,:,k,n) = 0.0
         else if (k >= klid) then
            exit
         end if
      end do
   end do

   end subroutine setZeroKlid4d


!===================================================================================
!BOP
! !IROUTINE: findKlid
   subroutine findKlid (klid, plid, ple, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(inout) :: klid ! index for pressure lid
   real, intent(in)       :: plid ! pressure lid [hPa]
   real, dimension(:,:,:), intent(in) :: ple  ! air pressure [Pa]

! !OUTPUT PARAMETERS:
   integer, intent(out) :: rc ! return code; 0 - all is good
!                                            1 - bad

! !DESCRIPTION: Finds corresponding vertical index for defined pressure lid
!
! !REVISION HISTORY:
!
! 25Aug2020 E.Sherman - Written 
!
! !Local Variables
   integer :: k, kk, j, i
   real :: plid_, diff, refDiff
   real, allocatable, dimension(:) :: pres  ! pressure at each model level [Pa]

!EOP
!----------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine findKlid:')
 call mpas_log_write('--- ubound(ple,3) = $i',intArgs=(/ubound(ple,3)/))


!  Begin...
   klid = 1
   rc = 0

!  convert from hPa to Pa
   plid_ = plid*100.0

   allocate(pres(ubound(ple,3)))

!  find pressure at each model level
   do k = 1, ubound(ple,3)
      pres(k) = ple(1,1,k)
   end do
   do k = 1,ubound(ple,3)
      call mpas_log_write('$i $r',intArgs=(/k/),realArgs=(/pres(k)/))
   enddo

!  find smallest absolute difference between plid and average pressure at each model level
   refDiff = 150000.0
   do k = 1, ubound(ple,3)
      diff = abs(pres(k) - plid_)
      if (diff < refDiff) then
         klid = k
         refDiff = diff
         call mpas_log_write('$i $i $r $r $r',intArgs=(/k,klid/),realArgs=(/pres(k),diff,refDiff/))
      end if
   end do

!  Check to make sure that all pressures at (i,j) were the same
   do j = 1, ubound(ple,2)
      do i = 1, ubound(ple,1)
         if (pres(klid) /= ple(i,j,klid)) then
            rc = 1
            return
         end if
      end do
   end do

   end subroutine findKlid

!=================================================================================================================
 end module GOCART2G_AeroGeneric
!=================================================================================================================


