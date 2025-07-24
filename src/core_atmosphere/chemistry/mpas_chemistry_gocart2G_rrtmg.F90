!=================================================================================================================
 module mpas_chemistry_gocart2G_rrtmg
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_config,mpas_pool_get_dimension
 use mpas_log

 use CA2G_bc_StateSpecs,only: CA2G_bc_State
 use CA2G_br_StateSpecs,only: CA2G_br_State
 use CA2G_oc_StateSpecs,only: CA2G_oc_State
 use DU2G_StateSpecs,only: DU2G_State
 use NI2G_StateSpecs,only: NI2G_State
 use SS2G_StateSpecs,only: SS2G_State
 use SU2G_StateSpecs,only: SU2G_State
 use GOCART2G_StateSpecs,only: GOCART2G_State


 implicit none
 private
 public:: GOCART2G_rrtmg


 contains


!==================================================================================================================
 subroutine GOCART2G_rrtmg(configs,mesh,diag_physics,CA2G_bc,CA2G_br,CA2G_oc,DU2G,NI2G,SS2G,SU2G, &
                           its,ite,jts,jte,kts,kte)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(CA2G_bc_State),intent(in):: CA2G_bc
 type(CA2G_br_State),intent(in):: CA2G_br
 type(CA2G_oc_State),intent(in):: CA2G_oc
 type(DU2G_State),intent(in)   :: DU2G
 type(NI2G_State),intent(in)   :: NI2G
 type(SS2G_State),intent(in)   :: SS2G
 type(SU2G_State),intent(in)   :: SU2G

 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!--- local arguments and arrays:
 logical,pointer:: do_CA2Gbc,do_CA2Gbr,do_CA2Goc,do_NI2G,do_DU2G,do_SS2G,do_SU2G

 integer,pointer:: nbndLW,nbndSW
 integer:: i,j,k,kk,n

 real(kind=RKIND),dimension(:,:,:),pointer:: bctau_lw,bcasy_lw,bcssa_lw,bctau_sw,bcasy_sw,bcssa_sw
 real(kind=RKIND),dimension(:,:,:),pointer:: brtau_lw,brasy_lw,brssa_lw,brtau_sw,brasy_sw,brssa_sw
 real(kind=RKIND),dimension(:,:,:),pointer:: octau_lw,ocasy_lw,ocssa_lw,octau_sw,ocasy_sw,ocssa_sw
 real(kind=RKIND),dimension(:,:,:),pointer:: dutau_lw,duasy_lw,dussa_lw,dutau_sw,duasy_sw,dussa_sw
 real(kind=RKIND),dimension(:,:,:),pointer:: nitau_lw,niasy_lw,nissa_lw,nitau_sw,niasy_sw,nissa_sw
 real(kind=RKIND),dimension(:,:,:),pointer:: sstau_lw,ssasy_lw,ssssa_lw,sstau_sw,ssasy_sw,ssssa_sw
 real(kind=RKIND),dimension(:,:,:),pointer:: sutau_lw,suasy_lw,sussa_lw,sutau_sw,suasy_sw,sussa_sw
 
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaer_lw,ssaaer_lw,asyaer_lw
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaer_sw,ssaaer_sw,asyaer_sw

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine GOCART2G_rrtmg:')


 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbc',do_CA2Gbc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Gbr',do_CA2Gbr)
 call mpas_pool_get_config(configs,'config_gocart2G_do_CA2Goc',do_CA2Goc)
 call mpas_pool_get_config(configs,'config_gocart2G_do_DU2G'  ,do_DU2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_NI2G'  ,do_NI2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SS2G'  ,do_SS2G  )
 call mpas_pool_get_config(configs,'config_gocart2G_do_SU2G'  ,do_SU2G  )

 call mpas_pool_get_dimension(mesh,'nbndLW',nbndLW)
 call mpas_pool_get_dimension(mesh,'nbndSW',nbndSW)


!--- total optical properties:
 call mpas_pool_get_array(diag_physics,'tauaer_lw',tauaer_lw)
 call mpas_pool_get_array(diag_physics,'ssaaer_lw',ssaaer_lw)
 call mpas_pool_get_array(diag_physics,'asyaer_lw',asyaer_lw)

 call mpas_pool_get_array(diag_physics,'tauaer_sw',tauaer_sw)
 call mpas_pool_get_array(diag_physics,'ssaaer_sw',ssaaer_sw)
 call mpas_pool_get_array(diag_physics,'asyaer_sw',asyaer_sw)

 tauaer_lw(:,:,:) = 0._RKIND
 ssaaer_lw(:,:,:) = 0._RKIND
 asyaer_lw(:,:,:) = 0._RKIND

 tauaer_sw(:,:,:) = 0._RKIND
 ssaaer_sw(:,:,:) = 0._RKIND
 asyaer_sw(:,:,:) = 0._RKIND


!--- black carbon:
 if(do_CA2Gbc) then
    call mpas_pool_get_array(diag_physics,'bctau_lw',bctau_lw)
    call mpas_pool_get_array(diag_physics,'bcssa_lw',bcssa_lw)
    call mpas_pool_get_array(diag_physics,'bcasy_lw',bcasy_lw)
    call mpas_pool_get_array(diag_physics,'bctau_sw',bctau_sw)
    call mpas_pool_get_array(diag_physics,'bcssa_sw',bcssa_sw)
    call mpas_pool_get_array(diag_physics,'bcasy_sw',bcasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             bctau_lw(:,kk,i) = CA2G_bc%bctau_lw(i,j,k,:)
             bcasy_lw(:,kk,i) = CA2G_bc%bcasy_lw(i,j,k,:)
             bcssa_lw(:,kk,i) = CA2G_bc%bcssa_lw(i,j,k,:)

             bctau_sw(:,kk,i) = CA2G_bc%bctau_sw(i,j,k,:)
             bcasy_sw(:,kk,i) = CA2G_bc%bcasy_sw(i,j,k,:)
             bcssa_sw(:,kk,i) = CA2G_bc%bcssa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + bctau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + bcasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + bcssa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + bctau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + bcasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + bcssa_sw(:,:,:)
 endif


!--- brown carbon:
 if(do_CA2Gbr) then
    call mpas_pool_get_array(diag_physics,'brtau_lw',brtau_lw)
    call mpas_pool_get_array(diag_physics,'brssa_lw',brssa_lw)
    call mpas_pool_get_array(diag_physics,'brasy_lw',brasy_lw)
    call mpas_pool_get_array(diag_physics,'brtau_sw',brtau_sw)
    call mpas_pool_get_array(diag_physics,'brssa_sw',brssa_sw)
    call mpas_pool_get_array(diag_physics,'brasy_sw',brasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             brtau_lw(:,kk,i) = CA2G_br%brtau_lw(i,j,k,:)
             brasy_lw(:,kk,i) = CA2G_br%brasy_lw(i,j,k,:)
             brssa_lw(:,kk,i) = CA2G_br%brssa_lw(i,j,k,:)

             brtau_sw(:,kk,i) = CA2G_br%brtau_sw(i,j,k,:)
             brasy_sw(:,kk,i) = CA2G_br%brasy_sw(i,j,k,:)
             brssa_sw(:,kk,i) = CA2G_br%brssa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + brtau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + brasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + brssa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + brtau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + brasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + brssa_sw(:,:,:)
 endif


!--- organic carbon:
 if(do_CA2Goc) then
    call mpas_pool_get_array(diag_physics,'octau_lw',octau_lw)
    call mpas_pool_get_array(diag_physics,'ocssa_lw',ocssa_lw)
    call mpas_pool_get_array(diag_physics,'ocasy_lw',ocasy_lw)
    call mpas_pool_get_array(diag_physics,'octau_sw',octau_sw)
    call mpas_pool_get_array(diag_physics,'ocssa_sw',ocssa_sw)
    call mpas_pool_get_array(diag_physics,'ocasy_sw',ocasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             octau_lw(:,kk,i) = CA2G_oc%octau_lw(i,j,k,:)
             ocasy_lw(:,kk,i) = CA2G_oc%ocasy_lw(i,j,k,:)
             ocssa_lw(:,kk,i) = CA2G_oc%ocssa_lw(i,j,k,:)

             octau_sw(:,kk,i) = CA2G_oc%octau_sw(i,j,k,:)
             ocasy_sw(:,kk,i) = CA2G_oc%ocasy_sw(i,j,k,:)
             ocssa_sw(:,kk,i) = CA2G_oc%ocssa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + octau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + ocasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + ocssa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + octau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + ocasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + ocssa_sw(:,:,:)
 endif


!--- mineral dust:
 if(do_DU2G) then
    call mpas_pool_get_array(diag_physics,'dutau_lw',dutau_lw)
    call mpas_pool_get_array(diag_physics,'dussa_lw',dussa_lw)
    call mpas_pool_get_array(diag_physics,'duasy_lw',duasy_lw)
    call mpas_pool_get_array(diag_physics,'dutau_sw',dutau_sw)
    call mpas_pool_get_array(diag_physics,'dussa_sw',dussa_sw)
    call mpas_pool_get_array(diag_physics,'duasy_sw',duasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             dutau_lw(:,kk,i) = DU2G%dutau_lw(i,j,k,:)
             duasy_lw(:,kk,i) = DU2G%duasy_lw(i,j,k,:)
             dussa_lw(:,kk,i) = DU2G%dussa_lw(i,j,k,:)

             dutau_sw(:,kk,i) = DU2G%dutau_sw(i,j,k,:)
             duasy_sw(:,kk,i) = DU2G%duasy_sw(i,j,k,:)
             dussa_sw(:,kk,i) = DU2G%dussa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + dutau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + duasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + dussa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + dutau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + duasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + dussa_sw(:,:,:)
 endif


!--- nitrate:
 if(do_NI2G) then
    call mpas_pool_get_array(diag_physics,'nitau_lw',nitau_lw)
    call mpas_pool_get_array(diag_physics,'nissa_lw',nissa_lw)
    call mpas_pool_get_array(diag_physics,'niasy_lw',niasy_lw)
    call mpas_pool_get_array(diag_physics,'nitau_sw',nitau_sw)
    call mpas_pool_get_array(diag_physics,'nissa_sw',nissa_sw)
    call mpas_pool_get_array(diag_physics,'niasy_sw',niasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             nitau_lw(:,kk,i) = NI2G%nitau_lw(i,j,k,:)
             niasy_lw(:,kk,i) = NI2G%niasy_lw(i,j,k,:)
             nissa_lw(:,kk,i) = NI2G%nissa_lw(i,j,k,:)

             nitau_sw(:,kk,i) = NI2G%nitau_sw(i,j,k,:)
             niasy_sw(:,kk,i) = NI2G%niasy_sw(i,j,k,:)
             nissa_sw(:,kk,i) = NI2G%nissa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + nitau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + niasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + nissa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + nitau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + niasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + nissa_sw(:,:,:)
 endif


!--- sea salt:
 if(do_SS2G) then
    call mpas_pool_get_array(diag_physics,'sstau_lw',sstau_lw)
    call mpas_pool_get_array(diag_physics,'ssssa_lw',ssssa_lw)
    call mpas_pool_get_array(diag_physics,'ssasy_lw',ssasy_lw)
    call mpas_pool_get_array(diag_physics,'sstau_sw',sstau_sw)
    call mpas_pool_get_array(diag_physics,'ssssa_sw',ssssa_sw)
    call mpas_pool_get_array(diag_physics,'ssasy_sw',ssasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             sstau_lw(:,kk,i) = SS2G%sstau_lw(i,j,k,:)
             ssasy_lw(:,kk,i) = SS2G%ssasy_lw(i,j,k,:)
             ssssa_lw(:,kk,i) = SS2G%ssssa_lw(i,j,k,:)

             sstau_sw(:,kk,i) = SS2G%sstau_sw(i,j,k,:)
             ssasy_sw(:,kk,i) = SS2G%ssasy_sw(i,j,k,:)
             ssssa_sw(:,kk,i) = SS2G%ssssa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + sstau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + ssasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + ssssa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + sstau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + ssasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + ssssa_sw(:,:,:)
 endif


!--- sulfate:
 if(do_SU2G) then
    call mpas_pool_get_array(diag_physics,'sutau_lw',sutau_lw)
    call mpas_pool_get_array(diag_physics,'sussa_lw',sussa_lw)
    call mpas_pool_get_array(diag_physics,'suasy_lw',suasy_lw)
    call mpas_pool_get_array(diag_physics,'sutau_sw',sutau_sw)
    call mpas_pool_get_array(diag_physics,'sussa_sw',sussa_sw)
    call mpas_pool_get_array(diag_physics,'suasy_sw',suasy_sw)

    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             kk = kte+1-k
             sutau_lw(:,kk,i) = SU2G%sutau_lw(i,j,k,:)
             suasy_lw(:,kk,i) = SU2G%suasy_lw(i,j,k,:)
             sussa_lw(:,kk,i) = SU2G%sussa_lw(i,j,k,:)

             sutau_sw(:,kk,i) = SU2G%sutau_sw(i,j,k,:)
             suasy_sw(:,kk,i) = SU2G%suasy_sw(i,j,k,:)
             sussa_sw(:,kk,i) = SU2G%sussa_sw(i,j,k,:)
          enddo
       enddo
    enddo

    tauaer_lw(:,:,:) = tauaer_lw(:,:,:) + sutau_lw(:,:,:)
    asyaer_lw(:,:,:) = asyaer_lw(:,:,:) + suasy_lw(:,:,:)
    ssaaer_lw(:,:,:) = ssaaer_lw(:,:,:) + sussa_lw(:,:,:)

    tauaer_sw(:,:,:) = tauaer_sw(:,:,:) + sutau_sw(:,:,:)
    asyaer_sw(:,:,:) = asyaer_sw(:,:,:) + suasy_sw(:,:,:)
    ssaaer_sw(:,:,:) = ssaaer_sw(:,:,:) + sussa_sw(:,:,:)
 endif


!---
 do n = 1,nbndLW
    where(ssaaer_lw(n,:,:) .gt. 0._RKIND)
       asyaer_lw(n,:,:) = asyaer_lw(n,:,:)/ssaaer_lw(n,:,:)
    end where
    where(tauaer_lw(n,:,:) .gt. 0._RKIND)
       ssaaer_lw(n,:,:) = ssaaer_lw(n,:,:)/tauaer_lw(n,:,:)
    end where
 enddo
 do n = 1,nbndSW
    where(ssaaer_sw(n,:,:) .gt. 0._RKIND)
       asyaer_sw(n,:,:) = asyaer_sw(n,:,:)/ssaaer_sw(n,:,:)
    end where
    where(tauaer_sw(n,:,:) .gt. 0._RKIND)
       ssaaer_sw(n,:,:) = ssaaer_sw(n,:,:)/tauaer_sw(n,:,:)
    end where
 enddo


 call mpas_log_write('--- end subroutine GOCART2G_rrtmg.')

 end subroutine GOCART2G_rrtmg

!==================================================================================================================
 end module mpas_chemistry_gocart2G_rrtmg
!==================================================================================================================
