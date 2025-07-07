!=================================================================================================================
 module cu_ntiedtke_ctrans_post
 use ccpp_kind_types,only: kind_phys

 implicit none
 private
 public:: cu_ntiedtke_ctrans_post_init,     &
          cu_ntiedtke_ctrans_post_finalize, &
          cu_ntiedtke_ctrans_post_run


 contains


!=================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_post_init
!!\html\include cu_ntiedtke_ctrans_post_init.html
!!
 subroutine cu_ntiedtke_ctrans_post_init(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine cu_ntiedtke_ctrans_post_init

!=================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_post_finalize
!!\html\include cu_ntiedtke_ctrans_post_finalize.html
!!
 subroutine cu_ntiedtke_ctrans_post_finalize(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: &
    errmsg      ! output error message (-).

 integer,intent(out):: &
    errflg      ! output error flag (-).

!-----------------------------------------------------------------------------------------------------------------

!--- output error flag and message:
 errflg = 0
 errmsg = " "

 end subroutine cu_ntiedtke_ctrans_post_finalize

!=================================================================================================================
!>\section arg_table_cu_ntiedtke_ctrans_post_run
!!\html\include cu_ntiedtke_ctrans_post_run.html
!!
 subroutine cu_ntiedtke_ctrans_post_run(its,ite,kts,kte,nchem,ptenc,ptencf,errmsg,errflg)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte
 integer,intent(in):: nchem

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,nchem):: ptenc

!--- inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,nchem):: ptencf

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!--- local variables and arrays:
 integer:: i,ic,k,pp,zz

!-----------------------------------------------------------------------------------------------------------------

 pp = 0
 do k = kts,kte
    zz = kte-pp
    do i = its,ite
       do ic = 1,nchem
          ptencf(i,k,ic) = ptenc(i,zz,ic)
       enddo
    enddo
    pp = pp+1
 enddo

 errmsg = 'cu_ntiedtke_ctrans_post_run OK'
 errflg = 0

 end subroutine cu_ntiedtke_ctrans_post_run

!=================================================================================================================
 end module cu_ntiedtke_ctrans_post
!=================================================================================================================
