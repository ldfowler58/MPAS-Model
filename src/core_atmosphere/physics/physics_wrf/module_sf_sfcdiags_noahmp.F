! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module module_sf_sfcdiags_noahmp
 use mpas_kind_types

 implicit none
 private
 public:: sfcdiags_noahmp


 contains


!=================================================================================================================
 subroutine sfcdiags_noahmp(isice,isurban,iswater,lcz_1,lcz_2,lcz_3,lcz_4,lcz_5,lcz_6,lcz_7,lcz_8,lcz_9,lcz_10, &
                    lcz_11,cp,R_d,rcp,xice_threshold,ivgtyp,chs2,cqs2,psfc,hfx,qfx,qsfc,tsk,xice,q2mbxy,q2mvxy, &
                    t2mbxy,t2mvxy,fvegxy,q2,t2,th2,ims,ime,jms,jme,its,ite,jts,jte)
!=================================================================================================================

!input arguments:
 integer,intent(in):: ims,ime,jms,jme,its,ite,jts,jte
 integer,intent(in):: isice,isurban,iswater
 integer,intent(in):: lcz_1,lcz_2,lcz_3,lcz_4,lcz_5,lcz_6,lcz_7,lcz_8,lcz_9,lcz_10,lcz_11
 integer,intent(in),dimension(ims:ime,jms:jme):: ivgtyp

 real(kind=RKIND),intent(in):: cp,R_d,rcp
 real(kind=RKIND),intent(in):: xice_threshold
 real(kind=RKIND),intent(in),dimension(ims:ime,jms:jme):: chs2,cqs2,psfc,hfx,qfx,qsfc,tsk,xice
 real(kind=RKIND),intent(in),dimension(ims:ime,jms:jme):: fvegxy,q2mbxy,q2mvxy,t2mbxy,t2mvxy

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(ims:ime,jms:jme):: q2,t2,th2

!local variables and arrays:
 integer:: i,j

!-----------------------------------------------------------------------------------------------------------------

 do j = jts,jte
    do i = its,ite
       if(ivgtyp(i,j)==iswater .or. (ivgtyp(i,j)==isice .and. xice(i,j).ge.xice_threshold)) then
          if(cqs2(i,j) .lt. 1.e-5) then
             q2(i,j) = qsfc(i,j)
          else
             q2(i,j) = qsfc(i,j)-qfx(i,j)/(psfc(i,j)/(R_d*tsk(i,j))*cqs2(i,j))
          endif
          if(chs2(i,j) .lt. 1.e-5) then
             t2(i,j) = tsk(i,j)
          else
             t2(i,j) = tsk(i,j)-hfx(i,j)/(psfc(i,j)/(R_d*tsk(i,j))*cp*chs2(i,j))
          endif
          th2(i,j) = t2(i,j)*(1.e5/psfc(i,j))**rcp
       elseif(ivgtyp(i,j).eq.isurban .or. (ivgtyp(i,j).eq.isice .and. xice(i,j).lt.xice_threshold) .or. &
              ivgtyp(i,j).eq.lcz_1   .or. ivgtyp(i,j).eq.lcz_2  .or. ivgtyp(i,j).eq.lcz_3          .or. &
              ivgtyp(i,j).eq.lcz_4   .or. ivgtyp(i,j).eq.lcz_5  .or. ivgtyp(i,j).eq.lcz_6          .or. &
              ivgtyp(i,j).eq.lcz_7   .or. ivgtyp(i,j).eq.lcz_8  .or. ivgtyp(i,j).eq.lcz_9          .or. &
              ivgtyp(i,j).eq.lcz_10  .or. ivgtyp(i,j).eq.lcz_11) then
          q2(i,j)  = q2mbxy(i,j)
          t2(i,j)  = t2mbxy(i,j)
          th2(i,j) = t2(i,j)*(1.e5/psfc(i,j))**rcp
       else
          q2(i,j)  = fvegxy(i,j)*q2mvxy(i,j) + (1.-fvegxy(i,j))*q2mbxy(i,j)
          t2(i,j)  = fvegxy(i,j)*t2mvxy(i,j) + (1.-fvegxy(i,j))*t2mbxy(i,j)
          th2(i,j) = t2(i,j)*(1.e5/psfc(i,j))**rcp
       endif
    enddo
 enddo

 end subroutine sfcdiags_noahmp

!=================================================================================================================
 end module module_sf_sfcdiags_noahmp
!=================================================================================================================
