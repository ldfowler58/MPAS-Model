!==================================================================================================================
 module mpas_chemistry_gocart2G_diagnostics
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_array,mpas_pool_get_dimension
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
 public:: CA2G_bc_diagnostics, &
          CA2G_br_diagnostics, &
          CA2G_oc_diagnostics, &
          DU2G_diagnostics,    &
          NI2G_diagnostics,    &
          SS2G_diagnostics,    &
          SU2G_diagnostics,    &
          GOCART2G_diagnostics


 contains


!==================================================================================================================
 subroutine CA2G_bc_diagnostics(mesh,CA2G_bc,CA2G_bc_diags,CA2G_bc_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(CA2G_bc_State),intent(in):: CA2G_bc

!inout arguments:
 type(mpas_pool_type),intent(inout):: CA2G_bc_diags
 type(mpas_pool_type),intent(inout):: CA2G_bc_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: bcEMAN,bcEMBB,bcEMBF,bcEMBG,bcEM_phobic,bcEM_philic
 real(kind=RKIND),dimension(:),pointer:: bcHYPHIL,bcSD_phobic,bcSD_philic,bcDP_phobic,bcDp_philic, &
                                         bcWT_phobic,bcWT_philic,bcSV_phobic,bcSV_philic
 real(kind=RKIND),dimension(:),pointer:: bcSMASS,bcCMASS,bcFLUXU,bcFLUXV
 real(kind=RKIND),dimension(:,:),pointer:: bcMASS,bcCONC

 real(kind=RKIND),dimension(:),pointer:: bcANGSTR,bcAERIDX
 real(kind=RKIND),dimension(:,:),pointer:: bcEXTTAU,bcSTEXTTAU,bcSCATAU,bcSTSCATAU
 real(kind=RKIND),dimension(:,:,:),pointer:: bcEXTCOEF,bcEXTCOEFRH20,bcEXTCOEFRH80,bcSCACOEF,bcSCACOEFRH20, &
                                             bcSCACOEFRH80,bcBCKCOEF

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine CA2G_bc_diagnostics:')

 call mpas_pool_get_array(CA2G_bc_diags,'bcEMAN',bcEMAN)
 call mpas_pool_get_array(CA2G_bc_diags,'bcEMBB',bcEMBB)
 call mpas_pool_get_array(CA2G_bc_diags,'bcEMBF',bcEMBF)
 call mpas_pool_get_array(CA2G_bc_diags,'bcEMBG',bcEMBG)
 call mpas_pool_get_array(CA2G_bc_diags,'bcEM_phobic',bcEM_phobic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcEM_philic',bcEM_philic)

 call mpas_pool_get_array(CA2G_bc_diags,'bcHYPHIL'   ,bcHYPHIL   )
 call mpas_pool_get_array(CA2G_bc_diags,'bcSD_phobic',bcSD_phobic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcSD_philic',bcSD_philic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcDP_phobic',bcDP_phobic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcDP_philic',bcDP_philic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcWT_phobic',bcWT_phobic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcWT_philic',bcWT_philic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcSV_phobic',bcSV_phobic)
 call mpas_pool_get_array(CA2G_bc_diags,'bcSV_philic',bcSV_philic)

 call mpas_pool_get_array(CA2G_bc_diags,'bcSMASS',bcSMASS)
 call mpas_pool_get_array(CA2G_bc_diags,'bcCMASS',bcCMASS)
 call mpas_pool_get_array(CA2G_bc_diags,'bcFLUXU',bcFLUXU)
 call mpas_pool_get_array(CA2G_bc_diags,'bcFLUXV',bcFLUXV)

 call mpas_pool_get_array(CA2G_bc_diags,'bcMASS',bcMASS)
 call mpas_pool_get_array(CA2G_bc_diags,'bcCONC',bcCONC)

 do j = jts,jte
    do i = its,ite
       bcEMAN(i) = CA2G_bc%bceman(i,j)
       bcEMBB(i) = CA2G_bc%bcembb(i,j)
       bcEMBF(i) = CA2G_bc%bcembf(i,j)
       bcEMBG(i) = CA2G_bc%bcembg(i,j)
       bcEM_phobic(i) = CA2G_bc%bcem(i,j,1)
       bcEM_philic(i) = CA2G_bc%bcem(i,j,2)
    enddo

    do i = its,ite
       bcHYPHIL(i)    = CA2G_bc%bchyphil(i,j)
       bcSD_phobic(i) = CA2G_bc%bcsd(i,j,1)
       bcSD_philic(i) = CA2G_bc%bcsd(i,j,2)
       bcDP_phobic(i) = CA2G_bc%bcdp(i,j,1)
       bcDP_philic(i) = CA2G_bc%bcdp(i,j,2)
       bcWT_phobic(i) = CA2G_bc%bcwt(i,j,1)
       bcWT_philic(i) = CA2G_bc%bcwt(i,j,2)
       bcSV_phobic(i) = CA2G_bc%bcsv(i,j,1)
       bcSV_philic(i) = CA2G_bc%bcsv(i,j,2)

       bcSMASS(i) = CA2G_bc%bcsmass(i,j)
       bcCMASS(i) = CA2G_bc%bccmass(i,j)
       bcFLUXU(i) = CA2G_BC%bcfluxu(i,j)
       bcFLUXV(i) = CA2G_BC%bcfluxv(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          bcMASS(kk,i) = CA2G_bc%bcmass(i,j,k)
          bcCONC(kk,i) = CA2G_bc%bcconc(i,j,k)
       enddo
    enddo
 enddo


!--- black carbon aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(CA2G_bc_aops,'bcEXTTAU'  ,bcEXTTAU  )
 call mpas_pool_get_array(CA2G_bc_aops,'bcSTEXTTAU',bcSTEXTTAU)
 call mpas_pool_get_array(CA2G_bc_aops,'bcSCATAU'  ,bcSCATAU  )
 call mpas_pool_get_array(CA2G_bc_aops,'bcSTSCATAU',bcSTSCATAU)
 call mpas_pool_get_array(CA2G_bc_aops,'bcANGSTR'  ,bcANGSTR  )
 call mpas_pool_get_array(CA2G_bc_aops,'bcAERIDX'  ,bcAERIDX  )

 call mpas_pool_get_array(CA2G_bc_aops,'bcEXTCOEF'    ,bcEXTCOEF    )
 call mpas_pool_get_array(CA2G_bc_aops,'bcEXTCOEFRH20',bcEXTCOEFRH20)
 call mpas_pool_get_array(CA2G_bc_aops,'bcEXTCOEFRH80',bcEXTCOEFRH80)
 call mpas_pool_get_array(CA2G_bc_aops,'bcSCACOEF'    ,bcSCACOEF    )
 call mpas_pool_get_array(CA2G_bc_aops,'bcSCACOEFRH20',bcSCACOEFRH20)
 call mpas_pool_get_array(CA2G_bc_aops,'bcSCACOEFRH80',bcSCACOEFRH80)
 call mpas_pool_get_array(CA2G_bc_aops,'bcBCKCOEF'    ,bcBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       bcANGSTR(i)   = CA2G_bc%bcangstr(i,j)
       bcAERIDX(i)   = CA2G_bc%bcaeridx(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          bcEXTTAU(n,i)   = CA2G_bc%bcexttau(i,j,n)
          bcSTEXTTAU(n,i) = CA2G_bc%bcstexttau(i,j,n)
          bcSCATAU(n,i)   = CA2G_bc%bcscatau(i,j,n)
          bcSTSCATAU(n,i) = CA2G_bc%bcstscatau(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             bcEXTCOEF(kk,n,i)     = CA2G_bc%bcextcoef(i,j,k,n)
             bcEXTCOEFRH20(kk,n,i) = CA2G_bc%bcextcoefrh20(i,j,k,n)
             bcEXTCOEFRH80(kk,n,i) = CA2G_bc%bcextcoefrh80(i,j,k,n)
             bcSCACOEF(kk,n,i)     = CA2G_bc%bcscacoef(i,j,k,n)
             bcSCACOEFRH20(kk,n,i) = CA2G_bc%bcscacoefrh20(i,j,k,n)
             bcSCACOEFRH80(kk,n,i) = CA2G_bc%bcscacoefrh80(i,j,k,n)
             bcBCKCOEF(kk,n,i)     = CA2G_bc%bcbckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine CA2G_bc_diagnostics.')

 end subroutine CA2G_bc_diagnostics

!==================================================================================================================
 subroutine CA2G_br_diagnostics(mesh,CA2G_br,CA2G_br_diags,CA2G_br_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(CA2G_br_State),intent(in):: CA2G_br

!inout arguments:
 type(mpas_pool_type),intent(inout):: CA2G_br_diags
 type(mpas_pool_type),intent(inout):: CA2G_br_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: brEMAN,brEMBB,brEMBF,brEMBG,brEM_phobic,brEM_philic
 real(kind=RKIND),dimension(:),pointer:: brHYPHIL,brSD_phobic,brSD_philic,brDP_phobic,brDp_philic, &
                                         brWT_phobic,brWT_philic,brSV_phobic,brSV_philic
 real(kind=RKIND),dimension(:),pointer:: brSMASS,brCMASS,brFLUXU,brFLUXV
 real(kind=RKIND),dimension(:,:),pointer:: brMASS,brCONC

 real(kind=RKIND),dimension(:),pointer:: brPSOA,brANGSTR,brAERIDX
 real(kind=RKIND),dimension(:,:),pointer:: brEXTTAU,brSTEXTTAU,brSCATAU,brSTSCATAU
 real(kind=RKIND),dimension(:,:,:),pointer:: brEXTCOEF,brEXTCOEFRH20,brEXTCOEFRH80,brSCACOEF,brSCACOEFRH20, &
                                             brSCACOEFRH80,brBCKCOEF

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine CA2G_br_diagnostics:')

 call mpas_pool_get_array(CA2G_br_diags,'brEMAN',brEMAN)
 call mpas_pool_get_array(CA2G_br_diags,'brEMBB',brEMBB)
 call mpas_pool_get_array(CA2G_br_diags,'brEMBF',brEMBF)
 call mpas_pool_get_array(CA2G_br_diags,'brEMBG',brEMBG)
 call mpas_pool_get_array(CA2G_br_diags,'brEM_phobic',brEM_phobic)
 call mpas_pool_get_array(CA2G_br_diags,'brEM_philic',brEM_philic)

 call mpas_pool_get_array(CA2G_br_diags,'brHYPHIL'   ,brHYPHIL   )
 call mpas_pool_get_array(CA2G_br_diags,'brSD_phobic',brSD_phobic)
 call mpas_pool_get_array(CA2G_br_diags,'brSD_philic',brSD_philic)
 call mpas_pool_get_array(CA2G_br_diags,'brDP_phobic',brDP_phobic)
 call mpas_pool_get_array(CA2G_br_diags,'brDP_philic',brDP_philic)
 call mpas_pool_get_array(CA2G_br_diags,'brWT_phobic',brWT_phobic)
 call mpas_pool_get_array(CA2G_br_diags,'brWT_philic',brWT_philic)
 call mpas_pool_get_array(CA2G_br_diags,'brSV_phobic',brSV_phobic)
 call mpas_pool_get_array(CA2G_br_diags,'brSV_philic',brSV_philic)

 call mpas_pool_get_array(CA2G_br_diags,'brSMASS',brSMASS)
 call mpas_pool_get_array(CA2G_br_diags,'brCMASS',brCMASS)
 call mpas_pool_get_array(CA2G_br_diags,'brFLUXU',brFLUXU)
 call mpas_pool_get_array(CA2G_br_diags,'brFLUXV',brFLUXV)

 call mpas_pool_get_array(CA2G_br_diags,'brMASS',brMASS)
 call mpas_pool_get_array(CA2G_br_diags,'brCONC',brCONC)

 do j = jts,jte
    do i = its,ite
       brEMAN(i) = CA2G_br%breman(i,j)
       brEMBB(i) = CA2G_br%brembb(i,j)
       brEMBF(i) = CA2G_br%brembf(i,j)
       brEMBG(i) = CA2G_br%brembg(i,j)
       brEM_phobic(i) = CA2G_br%brem(i,j,1)
       brEM_philic(i) = CA2G_br%brem(i,j,2)
    enddo

    do i = its,ite
       brHYPHIL(i)    = CA2G_br%brhyphil(i,j)
       brSD_phobic(i) = CA2G_br%brsd(i,j,1)
       brSD_philic(i) = CA2G_br%brsd(i,j,2)
       brDP_phobic(i) = CA2G_br%brdp(i,j,1)
       brDP_philic(i) = CA2G_br%brdp(i,j,2)
       brWT_phobic(i) = CA2G_br%brwt(i,j,1)
       brWT_philic(i) = CA2G_br%brwt(i,j,2)
       brSV_phobic(i) = CA2G_br%brsv(i,j,1)
       brSV_philic(i) = CA2G_br%brsv(i,j,2)

       brSMASS(i) = CA2G_br%brsmass(i,j)
       brCMASS(i) = CA2G_br%brcmass(i,j)
       brFLUXU(i) = CA2G_br%brfluxu(i,j)
       brFLUXV(i) = CA2G_br%brfluxv(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          brMASS(kk,i) = CA2G_br%brmass(i,j,k)
          brCONC(kk,i) = CA2G_br%brconc(i,j,k)
       enddo
    enddo
 enddo


!--- brown carbon aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(CA2G_br_aops,'brEXTTAU'  ,brEXTTAU  )
 call mpas_pool_get_array(CA2G_br_aops,'brSTEXTTAU',brSTEXTTAU)
 call mpas_pool_get_array(CA2G_br_aops,'brSCATAU'  ,brSCATAU  )
 call mpas_pool_get_array(CA2G_br_aops,'brSTSCATAU',brSTSCATAU)
 call mpas_pool_get_array(CA2G_br_aops,'brPSOA'    ,brPSOA    )
 call mpas_pool_get_array(CA2G_br_aops,'brANGSTR'  ,brANGSTR  )
 call mpas_pool_get_array(CA2G_br_aops,'brAERIDX'  ,brAERIDX  )

 call mpas_pool_get_array(CA2G_br_aops,'brEXTCOEF'    ,brEXTCOEF    )
 call mpas_pool_get_array(CA2G_br_aops,'brEXTCOEFRH20',brEXTCOEFRH20)
 call mpas_pool_get_array(CA2G_br_aops,'brEXTCOEFRH80',brEXTCOEFRH80)
 call mpas_pool_get_array(CA2G_br_aops,'brSCACOEF'    ,brSCACOEF    )
 call mpas_pool_get_array(CA2G_br_aops,'brSCACOEFRH20',brSCACOEFRH20)
 call mpas_pool_get_array(CA2G_br_aops,'brSCACOEFRH80',brSCACOEFRH80)
 call mpas_pool_get_array(CA2G_br_aops,'brBCKCOEF'    ,brBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       brPSOA(i)     = CA2G_br%brpsoa(i,j)
       brANGSTR(i)   = CA2G_br%brangstr(i,j)
       brAERIDX(i)   = CA2G_br%braeridx(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          brEXTTAU(n,i)   = CA2G_br%brexttau(i,j,n)
          brSTEXTTAU(n,i) = CA2G_br%brstexttau(i,j,n)
          brSCATAU(n,i)   = CA2G_br%brscatau(i,j,n)
          brSTSCATAU(n,i) = CA2G_br%brstscatau(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             brEXTCOEF(kk,n,i)     = CA2G_br%brextcoef(i,j,k,n)
             brEXTCOEFRH20(kk,n,i) = CA2G_br%brextcoefrh20(i,j,k,n)
             brEXTCOEFRH80(kk,n,i) = CA2G_br%brextcoefrh80(i,j,k,n)
             brSCACOEF(kk,n,i)     = CA2G_br%brscacoef(i,j,k,n)
             brSCACOEFRH20(kk,n,i) = CA2G_br%brscacoefrh20(i,j,k,n)
             brSCACOEFRH80(kk,n,i) = CA2G_br%brscacoefrh80(i,j,k,n)
             brBCKCOEF(kk,n,i)     = CA2G_br%brbckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine CA2G_br_diagnostics.')

 end subroutine CA2G_br_diagnostics

!==================================================================================================================
 subroutine CA2G_oc_diagnostics(mesh,CA2G_oc,CA2G_oc_diags,CA2G_oc_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(CA2G_oc_State),intent(in):: CA2G_oc

!inout arguments:
 type(mpas_pool_type),intent(inout):: CA2G_oc_diags
 type(mpas_pool_type),intent(inout):: CA2G_oc_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: ocEMAN,ocEMBB,ocEMBF,ocEMBG,ocEM_phobic,ocEM_philic
 real(kind=RKIND),dimension(:),pointer:: ocHYPHIL,ocSD_phobic,ocSD_philic,ocDP_phobic,ocDp_philic, &
                                         ocWT_phobic,ocWT_philic,ocSV_phobic,ocSV_philic
 real(kind=RKIND),dimension(:),pointer:: ocSMASS,ocCMASS,ocFLUXU,ocFLUXV
 real(kind=RKIND),dimension(:,:),pointer:: ocMASS,ocCONC

 real(kind=RKIND),dimension(:),pointer:: ocPSOA,ocANGSTR,ocAERIDX
 real(kind=RKIND),dimension(:,:),pointer:: ocEXTTAU,ocSTEXTTAU,ocSCATAU,ocSTSCATAU
 real(kind=RKIND),dimension(:,:,:),pointer:: ocEXTCOEF,ocEXTCOEFRH20,ocEXTCOEFRH80,ocSCACOEF,ocSCACOEFRH20, &
                                             ocSCACOEFRH80,ocBCKCOEF

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine CA2G_oc_diagnostics:')

 call mpas_pool_get_array(CA2G_oc_diags,'ocEMAN',ocEMAN)
 call mpas_pool_get_array(CA2G_oc_diags,'ocEMBB',ocEMBB)
 call mpas_pool_get_array(CA2G_oc_diags,'ocEMBF',ocEMBF)
 call mpas_pool_get_array(CA2G_oc_diags,'ocEMBG',ocEMBG)
 call mpas_pool_get_array(CA2G_oc_diags,'ocEM_phobic',ocEM_phobic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocEM_philic',ocEM_philic)

 call mpas_pool_get_array(CA2G_oc_diags,'ocHYPHIL'   ,ocHYPHIL   )
 call mpas_pool_get_array(CA2G_oc_diags,'ocSD_phobic',ocSD_phobic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocSD_philic',ocSD_philic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocDP_phobic',ocDP_phobic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocDP_philic',ocDP_philic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocWT_phobic',ocWT_phobic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocWT_philic',ocWT_philic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocSV_phobic',ocSV_phobic)
 call mpas_pool_get_array(CA2G_oc_diags,'ocSV_philic',ocSV_philic)

 call mpas_pool_get_array(CA2G_oc_diags,'ocSMASS',ocSMASS)
 call mpas_pool_get_array(CA2G_oc_diags,'ocCMASS',ocCMASS)
 call mpas_pool_get_array(CA2G_oc_diags,'ocFLUXU',ocFLUXU)
 call mpas_pool_get_array(CA2G_oc_diags,'ocFLUXV',ocFLUXV)

 call mpas_pool_get_array(CA2G_oc_diags,'ocMASS',ocMASS)
 call mpas_pool_get_array(CA2G_oc_diags,'ocCONC',ocCONC)

 do j = jts,jte
    do i = its,ite
       ocEMAN(i) = CA2G_oc%oceman(i,j)
       ocEMBB(i) = CA2G_oc%ocembb(i,j)
       ocEMBF(i) = CA2G_oc%ocembf(i,j)
       ocEMBG(i) = CA2G_oc%ocembg(i,j)
       ocEM_phobic(i) = CA2G_oc%ocem(i,j,1)
       ocEM_philic(i) = CA2G_oc%ocem(i,j,2)
    enddo

    do i = its,ite
       ocHYPHIL(i)    = CA2G_oc%ochyphil(i,j)
       ocSD_phobic(i) = CA2G_oc%ocsd(i,j,1)
       ocSD_philic(i) = CA2G_oc%ocsd(i,j,2)
       ocDP_phobic(i) = CA2G_oc%ocdp(i,j,1)
       ocDP_philic(i) = CA2G_oc%ocdp(i,j,2)
       ocWT_phobic(i) = CA2G_oc%ocwt(i,j,1)
       ocWT_philic(i) = CA2G_oc%ocwt(i,j,2)
       ocSV_phobic(i) = CA2G_oc%ocsv(i,j,1)
       ocSV_philic(i) = CA2G_oc%ocsv(i,j,2)

       ocSMASS(i) = CA2G_oc%ocsmass(i,j)
       ocCMASS(i) = CA2G_oc%occmass(i,j)
       ocFLUXU(i) = CA2G_oc%ocfluxu(i,j)
       ocFLUXV(i) = CA2G_oc%ocfluxv(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          ocMASS(kk,i) = CA2G_oc%ocmass(i,j,k)
          ocCONC(kk,i) = CA2G_oc%occonc(i,j,k)
       enddo
    enddo
 enddo


!--- organic carbon aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(CA2G_oc_aops,'ocEXTTAU'  ,ocEXTTAU  )
 call mpas_pool_get_array(CA2G_oc_aops,'ocSTEXTTAU',ocSTEXTTAU)
 call mpas_pool_get_array(CA2G_oc_aops,'ocSCATAU'  ,ocSCATAU  )
 call mpas_pool_get_array(CA2G_oc_aops,'ocSTSCATAU',ocSTSCATAU)
 call mpas_pool_get_array(CA2G_oc_aops,'ocPSOA'    ,ocPSOA    )
 call mpas_pool_get_array(CA2G_oc_aops,'ocANGSTR'  ,ocANGSTR  )
 call mpas_pool_get_array(CA2G_oc_aops,'ocAERIDX'  ,ocAERIDX  )

 call mpas_pool_get_array(CA2G_oc_aops,'ocEXTCOEF'    ,ocEXTCOEF    )
 call mpas_pool_get_array(CA2G_oc_aops,'ocEXTCOEFRH20',ocEXTCOEFRH20)
 call mpas_pool_get_array(CA2G_oc_aops,'ocEXTCOEFRH80',ocEXTCOEFRH80)
 call mpas_pool_get_array(CA2G_oc_aops,'ocSCACOEF'    ,ocSCACOEF    )
 call mpas_pool_get_array(CA2G_oc_aops,'ocSCACOEFRH20',ocSCACOEFRH20)
 call mpas_pool_get_array(CA2G_oc_aops,'ocSCACOEFRH80',ocSCACOEFRH80)
 call mpas_pool_get_array(CA2G_oc_aops,'ocBCKCOEF'    ,ocBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       ocPSOA(i)     = CA2G_oc%ocpsoa(i,j)
       ocANGSTR(i)   = CA2G_oc%ocangstr(i,j)
       ocAERIDX(i)   = CA2G_oc%ocaeridx(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          ocEXTTAU(n,i)   = CA2G_oc%ocexttau(i,j,n)
          ocSTEXTTAU(n,i) = CA2G_oc%ocstexttau(i,j,n)
          ocSCATAU(n,i)   = CA2G_oc%ocscatau(i,j,n)
          ocSTSCATAU(n,i) = CA2G_oc%ocstscatau(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             ocEXTCOEF(kk,n,i)     = CA2G_oc%ocextcoef(i,j,k,n)
             ocEXTCOEFRH20(kk,n,i) = CA2G_oc%ocextcoefrh20(i,j,k,n)
             ocEXTCOEFRH80(kk,n,i) = CA2G_oc%ocextcoefrh80(i,j,k,n)
             ocSCACOEF(kk,n,i)     = CA2G_oc%ocscacoef(i,j,k,n)
             ocSCACOEFRH20(kk,n,i) = CA2G_oc%ocscacoefrh20(i,j,k,n)
             ocSCACOEFRH80(kk,n,i) = CA2G_oc%ocscacoefrh80(i,j,k,n)
             ocBCKCOEF(kk,n,i)     = CA2G_oc%ocbckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine CA2G_oc_diagnostics.')

 end subroutine CA2G_oc_diagnostics

!==================================================================================================================
 subroutine DU2G_diagnostics(mesh,DU2G,DU2G_diags,DU2G_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(DU2G_State),intent(in):: DU2G

!inout arguments:
 type(mpas_pool_type),intent(inout):: DU2G_diags
 type(mpas_pool_type),intent(inout):: DU2G_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: duEM_bin1,duEM_bin2,duEM_bin3,duEM_bin4,duEM_bin5
 real(kind=RKIND),dimension(:),pointer:: duSD_bin1,duSD_bin2,duSD_bin3,duSD_bin4,duSD_bin5
 real(kind=RKIND),dimension(:),pointer:: duDP_bin1,duDP_bin2,duDP_bin3,duDP_bin4,duDP_bin5
 real(kind=RKIND),dimension(:),pointer:: duWT_bin1,duWT_bin2,duWT_bin3,duWT_bin4,duWT_bin5
 real(kind=RKIND),dimension(:),pointer:: duSV_bin1,duSV_bin2,duSV_bin3,duSV_bin4,duSV_bin5
 real(kind=RKIND),dimension(:),pointer:: duSMASS,duSMASS25,duCMASS,duCMASS25,duFLUXU,duFLUXV
 real(kind=RKIND),dimension(:,:),pointer:: duMASS,duMASS25,duCONC

 real(kind=RKIND),dimension(:),pointer:: duANGSTR,duAERIDX
 real(kind=RKIND),dimension(:,:),pointer:: duEXTTAU,duSTEXTTAU,duSCATAU,duSTSCATAU
 real(kind=RKIND),dimension(:,:),pointer:: duEXTT25,duSCAT25,duEXTTFM,duSCATFM
 real(kind=RKIND),dimension(:,:,:),pointer:: duEXTCOEF,duEXTCOEFRH20,duEXTCOEFRH80,duSCACOEF,duSCACOEFRH20, &
                                             duSCACOEFRH80,duBCKCOEF

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine DU2G_diagnostics:')

 call mpas_pool_get_array(DU2G_diags,'duEM_bin1',duEM_bin1)
 call mpas_pool_get_array(DU2G_diags,'duEM_bin2',duEM_bin2)
 call mpas_pool_get_array(DU2G_diags,'duEM_bin3',duEM_bin3)
 call mpas_pool_get_array(DU2G_diags,'duEM_bin4',duEM_bin4)
 call mpas_pool_get_array(DU2G_diags,'duEM_bin5',duEM_bin5)

 call mpas_pool_get_array(DU2G_diags,'duSD_bin1',duSD_bin1)
 call mpas_pool_get_array(DU2G_diags,'duSD_bin2',duSD_bin2)
 call mpas_pool_get_array(DU2G_diags,'duSD_bin3',duSD_bin3)
 call mpas_pool_get_array(DU2G_diags,'duSD_bin4',duSD_bin4)
 call mpas_pool_get_array(DU2G_diags,'duSD_bin5',duSD_bin5)

 call mpas_pool_get_array(DU2G_diags,'duDP_bin1',duDP_bin1)
 call mpas_pool_get_array(DU2G_diags,'duDP_bin2',duDP_bin2)
 call mpas_pool_get_array(DU2G_diags,'duDP_bin3',duDP_bin3)
 call mpas_pool_get_array(DU2G_diags,'duDP_bin4',duDP_bin4)
 call mpas_pool_get_array(DU2G_diags,'duDP_bin5',duDP_bin5)

 call mpas_pool_get_array(DU2G_diags,'duWT_bin1',duWT_bin1)
 call mpas_pool_get_array(DU2G_diags,'duWT_bin2',duWT_bin2)
 call mpas_pool_get_array(DU2G_diags,'duWT_bin3',duWT_bin3)
 call mpas_pool_get_array(DU2G_diags,'duWT_bin4',duWT_bin4)
 call mpas_pool_get_array(DU2G_diags,'duWT_bin5',duWT_bin5)

 call mpas_pool_get_array(DU2G_diags,'duSV_bin1',duSV_bin1)
 call mpas_pool_get_array(DU2G_diags,'duSV_bin2',duSV_bin2)
 call mpas_pool_get_array(DU2G_diags,'duSV_bin3',duSV_bin3)
 call mpas_pool_get_array(DU2G_diags,'duSV_bin4',duSV_bin4)
 call mpas_pool_get_array(DU2G_diags,'duSV_bin5',duSV_bin5)

 call mpas_pool_get_array(DU2G_diags,'duSMASS'  ,duSMASS  )
 call mpas_pool_get_array(DU2G_diags,'duCMASS'  ,duCMASS  )
 call mpas_pool_get_array(DU2G_diags,'duSMASS25',duSMASS25)
 call mpas_pool_get_array(DU2G_diags,'duCMASS25',duCMASS25)
 call mpas_pool_get_array(DU2G_diags,'duFLUXU'  ,duFLUXU  )
 call mpas_pool_get_array(DU2G_diags,'duFLUXV'  ,duFLUXV  )

 call mpas_pool_get_array(DU2G_diags,'duMASS'  ,duMASS  )
 call mpas_pool_get_array(DU2G_diags,'duMASS25',duMASS25)
 call mpas_pool_get_array(DU2G_diags,'duCONC'  ,duCONC  )

 do j = jts,jte
    do i = its,ite
       duEM_bin1(i) = DU2G%duem(i,j,1)
       duEM_bin2(i) = DU2G%duem(i,j,2)
       duEM_bin3(i) = DU2G%duem(i,j,3)
       duEM_bin4(i) = DU2G%duem(i,j,4)
       duEM_bin5(i) = DU2G%duem(i,j,5)
       duSD_bin1(i) = DU2G%dusd(i,j,1)
       duSD_bin2(i) = DU2G%dusd(i,j,2)
       duSD_bin3(i) = DU2G%dusd(i,j,3)
       duSD_bin4(i) = DU2G%dusd(i,j,4)
       duSD_bin5(i) = DU2G%dusd(i,j,5)
       duDP_bin1(i) = DU2G%dudp(i,j,1)
       duDP_bin2(i) = DU2G%dudp(i,j,2)
       duDP_bin3(i) = DU2G%dudp(i,j,3)
       duDP_bin4(i) = DU2G%dudp(i,j,4)
       duDP_bin5(i) = DU2G%dudp(i,j,5)
       duWT_bin1(i) = DU2G%duwt(i,j,1)
       duWT_bin2(i) = DU2G%duwt(i,j,2)
       duWT_bin3(i) = DU2G%duwt(i,j,3)
       duWT_bin4(i) = DU2G%duwt(i,j,4)
       duWT_bin5(i) = DU2G%duwt(i,j,5)
       duSV_bin1(i) = DU2G%dusv(i,j,1)
       duSV_bin2(i) = DU2G%dusv(i,j,2)
       duSV_bin3(i) = DU2G%dusv(i,j,3)
       duSV_bin4(i) = DU2G%dusv(i,j,4)
       duSV_bin5(i) = DU2G%dusv(i,j,5)

       duSMASS(i)   = DU2G%dusmass(i,j)
       duCMASS(i)   = DU2G%ducmass(i,j)
       duSMASS25(i) = DU2G%dusmass25(i,j)
       duCMASS25(i) = DU2G%ducmass25(i,j)
       duFLUXU(i)   = DU2G%dufluxu(i,j)
       duFLUXV(i)   = DU2G%dufluxv(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          duMASS(kk,i)   = DU2G%dumass(i,j,k)
          duMASS25(kk,i) = DU2G%dumass25(i,j,k)
          duCONC(kk,i)   = DU2G%duconc(i,j,k)
       enddo
    enddo
 enddo


!--- dust aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(DU2G_aops,'duEXTTAU'  ,duEXTTAU  )
 call mpas_pool_get_array(DU2G_aops,'duSTEXTTAU',duSTEXTTAU)
 call mpas_pool_get_array(DU2G_aops,'duSCATAU'  ,duSCATAU  )
 call mpas_pool_get_array(DU2G_aops,'duSTSCATAU',duSTSCATAU)
 call mpas_pool_get_array(DU2G_aops,'duANGSTR'  ,duANGSTR  )
 call mpas_pool_get_array(DU2G_aops,'duAERIDX'  ,duAERIDX  )
 call mpas_pool_get_array(DU2G_aops,'duEXTT25'  ,duEXTT25  )
 call mpas_pool_get_array(DU2G_aops,'duSCAT25'  ,duSCAT25  )
 call mpas_pool_get_array(DU2G_aops,'duEXTTFM'  ,duEXTTFM  )
 call mpas_pool_get_array(DU2G_aops,'duSCATFM'  ,duSCATFM  )

 call mpas_pool_get_array(DU2G_aops,'duEXTCOEF'    ,duEXTCOEF    )
 call mpas_pool_get_array(DU2G_aops,'duEXTCOEFRH20',duEXTCOEFRH20)
 call mpas_pool_get_array(DU2G_aops,'duEXTCOEFRH80',duEXTCOEFRH80)
 call mpas_pool_get_array(DU2G_aops,'duSCACOEF'    ,duSCACOEF    )
 call mpas_pool_get_array(DU2G_aops,'duSCACOEFRH20',duSCACOEFRH20)
 call mpas_pool_get_array(DU2G_aops,'duSCACOEFRH80',duSCACOEFRH80)
 call mpas_pool_get_array(DU2G_aops,'duBCKCOEF'    ,duBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       duANGSTR(i) = DU2G%duangstr(i,j)
       duAERIDX(i) = DU2G%duaeridx(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          duEXTTAU(n,i)   = DU2G%duexttau(i,j,n)
          duSTEXTTAU(n,i) = DU2G%dustexttau(i,j,n)
          duSCATAU(n,i)   = DU2G%duscatau(i,j,n)
          duSTSCATAU(n,i) = DU2G%dustscatau(i,j,n)
          duEXTT25(n,i)   = DU2G%duextt25(i,j,n)
          duSCAT25(n,i)   = DU2G%duscat25(i,j,n)
          duEXTTFM(n,i)   = DU2G%duexttfm(i,j,n)
          duSCATFM(n,i)   = DU2G%duscatfm(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             duEXTCOEF(kk,n,i)     = DU2G%duextcoef(i,j,k,n)
             duEXTCOEFRH20(kk,n,i) = DU2G%duextcoefrh20(i,j,k,n)
             duEXTCOEFRH80(kk,n,i) = DU2G%duextcoefrh80(i,j,k,n)
             duSCACOEF(kk,n,i)     = DU2G%duscacoef(i,j,k,n)
             duSCACOEFRH20(kk,n,i) = DU2G%duscacoefrh20(i,j,k,n)
             duSCACOEFRH80(kk,n,i) = DU2G%duscacoefrh80(i,j,k,n)
             duBCKCOEF(kk,n,i)     = DU2G%dubckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine DU2G_diagnostics.')

 end subroutine DU2G_diagnostics

!==================================================================================================================
 subroutine NI2G_diagnostics(mesh,NI2G,NI2G_diags,NI2G_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(NI2G_State),intent(in):: NI2G

!inout arguments:
 type(mpas_pool_type),intent(inout):: NI2G_diags
 type(mpas_pool_type),intent(inout):: NI2G_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: niHT_bin1,niHT_bin2,niHT_bin3
 real(kind=RKIND),dimension(:),pointer:: niSD_bin1,niSD_bin2,niSD_bin3
 real(kind=RKIND),dimension(:),pointer:: niDP_bin1,niDP_bin2,niDP_bin3
 real(kind=RKIND),dimension(:),pointer:: niWT_bin1,niWT_bin2,niWT_bin3
 real(kind=RKIND),dimension(:),pointer:: niSV_bin1,niSV_bin2,niSV_bin3

 real(kind=RKIND),dimension(:),pointer:: nh3EM,nh3DP,nh3WT,nh3SV
 real(kind=RKIND),dimension(:),pointer:: nh4SD,nh4DP,nh4WT,nh4SV
 real(kind=RKIND),dimension(:),pointer:: niPNO3AQ,niPNH4AQ,niPNH3AQ
 real(kind=RKIND),dimension(:),pointer:: hno3SMASS,nh3SMASS,nh4SMASS,niSMASS,niSMASS25
 real(kind=RKIND),dimension(:),pointer:: hno3CMASS,nh3CMASS,nh4CMASS,niCMASS,niCMASS25

 real(kind=RKIND),dimension(:),pointer:: niFLUXU,niFLUXV
 real(kind=RKIND),dimension(:,:),pointer:: niCONC,niCONC25
 real(kind=RKIND),dimension(:,:),pointer:: nh3MASS,nh4MASS,niMASS,niMASS25,hno3CONC,nh3CONC,nh4CONC

 real(kind=RKIND),dimension(:),pointer:: niANGSTR
 real(kind=RKIND),dimension(:,:),pointer:: niEXTTAU,niSTEXTTAU,niSCATAU,niSTSCATAU
 real(kind=RKIND),dimension(:,:),pointer:: niEXTT25,niSCAT25,niEXTTFM,niSCATFM
 real(kind=RKIND),dimension(:,:,:),pointer:: niEXTCOEF,niEXTCOEFRH20,niEXTCOEFRH80,niSCACOEF,niSCACOEFRH20, &
                                             niSCACOEFRH80,niBCKCOEF

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine NI2G_diagnostics:')

 call mpas_pool_get_array(NI2G_diags,'niHT_bin1',niHT_bin1)
 call mpas_pool_get_array(NI2G_diags,'niHT_bin2',niHT_bin2)
 call mpas_pool_get_array(NI2G_diags,'niHT_bin3',niHT_bin3)
 call mpas_pool_get_array(NI2G_diags,'niSD_bin1',niSD_bin1)
 call mpas_pool_get_array(NI2G_diags,'niSD_bin2',niSD_bin2)
 call mpas_pool_get_array(NI2G_diags,'niSD_bin3',niSD_bin3)
 call mpas_pool_get_array(NI2G_diags,'niDP_bin1',niDP_bin1)
 call mpas_pool_get_array(NI2G_diags,'niDP_bin2',niDP_bin2)
 call mpas_pool_get_array(NI2G_diags,'niDP_bin3',niDP_bin3)
 call mpas_pool_get_array(NI2G_diags,'niWT_bin1',niWT_bin1)
 call mpas_pool_get_array(NI2G_diags,'niWT_bin2',niWT_bin2)
 call mpas_pool_get_array(NI2G_diags,'niWT_bin3',niWT_bin3)
 call mpas_pool_get_array(NI2G_diags,'niSV_bin1',niSV_bin1)
 call mpas_pool_get_array(NI2G_diags,'niSV_bin2',niSV_bin2)
 call mpas_pool_get_array(NI2G_diags,'niSV_bin3',niSV_bin3)

 call mpas_pool_get_array(NI2G_diags,'nh3EM',nh3EM)
 call mpas_pool_get_array(NI2G_diags,'nh3DP',nh3DP)
 call mpas_pool_get_array(NI2G_diags,'nh3WT',nh3WT)
 call mpas_pool_get_array(NI2G_diags,'nh3SV',nh3SV)
 call mpas_pool_get_array(NI2G_diags,'nh4SD',nh4SD)
 call mpas_pool_get_array(NI2G_diags,'nh4DP',nh4DP)
 call mpas_pool_get_array(NI2G_diags,'nh4WT',nh4WT)
 call mpas_pool_get_array(NI2G_diags,'nh4SV',nh4SV)

 call mpas_pool_get_array(NI2G_diags,'hno3SMASS',hno3SMASS)
 call mpas_pool_get_array(NI2G_diags,'hno3CMASS',hno3CMASS)
 call mpas_pool_get_array(NI2G_diags,'nh3SMASS' ,nh3SMASS )
 call mpas_pool_get_array(NI2G_diags,'nh3CMASS' ,nh3CMASS )
 call mpas_pool_get_array(NI2G_diags,'nh4SMASS' ,nh4SMASS )
 call mpas_pool_get_array(NI2G_diags,'nh4CMASS' ,nh4CMASS )
 call mpas_pool_get_array(NI2G_diags,'niSMASS'  ,niSMASS  )
 call mpas_pool_get_array(NI2G_diags,'niCMASS'  ,niCMASS  )
 call mpas_pool_get_array(NI2G_diags,'niSMASS25',niSMASS25)
 call mpas_pool_get_array(NI2G_diags,'niCMASS25',niCMASS25)

 call mpas_pool_get_array(NI2G_diags,'niPNO3AQ',niPNO3AQ)
 call mpas_pool_get_array(NI2G_diags,'niPNH4AQ',niPNH4AQ)
 call mpas_pool_get_array(NI2G_diags,'niPNH3AQ',niPNH3AQ)

 call mpas_pool_get_array(NI2G_diags,'niFLUXU',niFLUXU)
 call mpas_pool_get_array(NI2G_diags,'niFLUXV',niFLUXV)

 call mpas_pool_get_array(NI2G_diags,'nh3MASS' ,nh3MASS )
 call mpas_pool_get_array(NI2G_diags,'nh4MASS' ,nh4MASS )
 call mpas_pool_get_array(NI2G_diags,'niMASS'  ,niMASS  )
 call mpas_pool_get_array(NI2G_diags,'niMASS25',niMASS25)
 call mpas_pool_get_array(NI2G_diags,'niCONC'  ,niCONC  )
 call mpas_pool_get_array(NI2G_diags,'niCONC25',niCONC25)
 call mpas_pool_get_array(NI2G_diags,'hno3CONC',hno3CONC)
 call mpas_pool_get_array(NI2G_diags,'nh3CONC' ,nh3CONC )
 call mpas_pool_get_array(NI2G_diags,'nh4CONC' ,nh4CONC )

 do j = jts,jte
    do i = its,ite
       niHT_bin1(i) = NI2G%niht(i,j,1)
       niHT_bin2(i) = NI2G%niht(i,j,2)
       niHT_bin3(i) = NI2G%niht(i,j,3)
       niSD_bin1(i) = NI2G%nisd(i,j,1)
       niSD_bin2(i) = NI2G%nisd(i,j,2)
       niSD_bin3(i) = NI2G%nisd(i,j,3)
       niDP_bin1(i) = NI2G%nidp(i,j,1)
       niDP_bin2(i) = NI2G%nidp(i,j,2)
       niDP_bin3(i) = NI2G%nidp(i,j,3)
       niWT_bin1(i) = NI2G%niwt(i,j,1)
       niWT_bin2(i) = NI2G%niwt(i,j,2)
       niWT_bin3(i) = NI2G%niwt(i,j,3)
       niSV_bin1(i) = NI2G%nisv(i,j,1)
       niSV_bin2(i) = NI2G%nisv(i,j,2)
       niSV_bin3(i) = NI2G%nisv(i,j,3)

       nh3EM(i) = NI2G%nh3em(i,j)
       nh3DP(i) = NI2G%nh3dp(i,j)
       nh3WT(i) = NI2G%nh3wt(i,j)
       nh3SV(i) = NI2G%nh3sv(i,j)
       nh4SD(i) = NI2G%nh4sd(i,j)
       nh4DP(i) = NI2G%nh4dp(i,j)
       nh4WT(i) = NI2G%nh4wt(i,j)
       nh4SV(i) = NI2G%nh4sv(i,j)

       hno3SMASS(i) = NI2G%hno3smass(i,j)
       hno3CMASS(i) = NI2G%hno3cmass(i,j)
       nh3SMASS(i)  = NI2G%nh3smass(i,j)
       nh3CMASS(i)  = NI2G%nh3cmass(i,j)
       nh4SMASS(i)  = NI2G%nh4smass(i,j)
       nh4CMASS(i)  = NI2G%nh4cmass(i,j)
       niSMASS(i)   = NI2G%nismass(i,j)
       niCMASS(i)   = NI2G%nicmass(i,j)
       niSMASS25(i) = NI2G%nismass25(i,j)
       niCMASS25(i) = NI2G%nicmass25(i,j)

       niPNO3AQ(i) = NI2G%nipno3aq(i,j)
       niPNH4AQ(i) = NI2G%nipnh4aq(i,j)
       niPNH3AQ(i) = NI2G%nipnh3aq(i,j)

       niFLUXU(i) = NI2G%nifluxu(i,j)
       niFLUXV(i) = NI2G%nifluxv(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          nh3MASS(kk,i)  = NI2G%nh3mass(i,j,k)
          nh4MASS(kk,i)  = NI2G%nh4mass(i,j,k)
          niMASS(kk,i)   = NI2G%nimass(i,j,k)
          niMASS25(kk,i) = NI2G%nimass25(i,j,k)
          niCONC(kk,i)   = NI2G%niconc(i,j,k)
          niCONC25(kk,i) = NI2G%niconc25(i,j,k)
          hno3CONC(kk,i) = NI2G%hno3conc(i,j,k)
          nh3CONC(kk,i)  = NI2G%nh3conc(i,j,k)
          nh4CONC(kk,i)  = NI2G%nh4conc(i,j,k)
       enddo
    enddo
 enddo


!--- nitrate aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(NI2G_aops,'niEXTTAU'  ,niEXTTAU  )
 call mpas_pool_get_array(NI2G_aops,'niSTEXTTAU',niSTEXTTAU)
 call mpas_pool_get_array(NI2G_aops,'niSCATAU'  ,niSCATAU  )
 call mpas_pool_get_array(NI2G_aops,'niSTSCATAU',niSTSCATAU)
 call mpas_pool_get_array(NI2G_aops,'niANGSTR'  ,niANGSTR  )
 call mpas_pool_get_array(NI2G_aops,'niEXTT25'  ,niEXTT25  )
 call mpas_pool_get_array(NI2G_aops,'niSCAT25'  ,niSCAT25  )
 call mpas_pool_get_array(NI2G_aops,'niEXTTFM'  ,niEXTTFM  )
 call mpas_pool_get_array(NI2G_aops,'niSCATFM'  ,niSCATFM  )

 call mpas_pool_get_array(NI2G_aops,'niEXTCOEF'    ,niEXTCOEF    )
 call mpas_pool_get_array(NI2G_aops,'niEXTCOEFRH20',niEXTCOEFRH20)
 call mpas_pool_get_array(NI2G_aops,'niEXTCOEFRH80',niEXTCOEFRH80)
 call mpas_pool_get_array(NI2G_aops,'niSCACOEF'    ,niSCACOEF    )
 call mpas_pool_get_array(NI2G_aops,'niSCACOEFRH20',niSCACOEFRH20)
 call mpas_pool_get_array(NI2G_aops,'niSCACOEFRH80',niSCACOEFRH80)
 call mpas_pool_get_array(NI2G_aops,'niBCKCOEF'    ,niBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       niANGSTR(i)   = NI2G%niangstr(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          niEXTTAU(n,i)   = NI2G%niexttau(i,j,n)
          niSTEXTTAU(n,i) = NI2G%nistexttau(i,j,n)
          niSCATAU(n,i)   = NI2G%niscatau(i,j,n)
          niSTSCATAU(n,i) = NI2G%nistscatau(i,j,n)
          niEXTT25(n,i)   = NI2G%niextt25(i,j,n)
          niSCAT25(n,i)   = NI2G%niscat25(i,j,n)
          niEXTTFM(n,i)   = NI2G%niexttfm(i,j,n)
          niSCATFM(n,i)   = NI2G%niscatfm(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             niEXTCOEF(kk,n,i)     = NI2G%niextcoef(i,j,k,n)
             niEXTCOEFRH20(kk,n,i) = NI2G%niextcoefrh20(i,j,k,n)
             niEXTCOEFRH80(kk,n,i) = NI2G%niextcoefrh80(i,j,k,n)
             niSCACOEF(kk,n,i)     = NI2G%niscacoef(i,j,k,n)
             niSCACOEFRH20(kk,n,i) = NI2G%niscacoefrh20(i,j,k,n)
             niSCACOEFRH80(kk,n,i) = NI2G%niscacoefrh80(i,j,k,n)
             niBCKCOEF(kk,n,i)     = NI2G%nibckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine NI2G_diagnostics.')

 end subroutine NI2G_diagnostics

!==================================================================================================================
 subroutine SS2G_diagnostics(mesh,SS2G,SS2G_diags,SS2G_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(SS2G_State),intent(in):: SS2G

!inout arguments:
 type(mpas_pool_type),intent(inout):: SS2G_diags
 type(mpas_pool_type),intent(inout):: SS2G_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: ssEM_bin1,ssEM_bin2,ssEM_bin3,ssEM_bin4,ssEM_bin5
 real(kind=RKIND),dimension(:),pointer:: ssSD_bin1,ssSD_bin2,ssSD_bin3,ssSD_bin4,ssSD_bin5
 real(kind=RKIND),dimension(:),pointer:: ssDP_bin1,ssDP_bin2,ssDP_bin3,ssDP_bin4,ssDP_bin5
 real(kind=RKIND),dimension(:),pointer:: ssWT_bin1,ssWT_bin2,ssWT_bin3,ssWT_bin4,ssWT_bin5
 real(kind=RKIND),dimension(:),pointer:: ssSV_bin1,ssSV_bin2,ssSV_bin3,ssSV_bin4,ssSV_bin5
 real(kind=RKIND),dimension(:),pointer:: ssSMASS,ssSMASS25,ssCMASS,ssCMASS25,ssFLUXU,ssFLUXV
 real(kind=RKIND),dimension(:,:),pointer:: ssMASS,ssMASS25,ssCONC

 real(kind=RKIND),dimension(:),pointer:: ssANGSTR,ssAERIDX
 real(kind=RKIND),dimension(:,:),pointer:: ssEXTTAU,ssSTEXTTAU,ssSCATAU,ssSTSCATAU
 real(kind=RKIND),dimension(:,:),pointer:: ssEXTT25,ssSCAT25,ssEXTTFM,ssSCATFM
 real(kind=RKIND),dimension(:,:,:),pointer:: ssEXTCOEF,ssEXTCOEFRH20,ssEXTCOEFRH80,ssSCACOEF,ssSCACOEFRH20, &
                                             ssSCACOEFRH80,ssBCKCOEF

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine SS2G_diagnostics:')

 call mpas_pool_get_array(SS2G_diags,'ssEM_bin1',ssEM_bin1)
 call mpas_pool_get_array(SS2G_diags,'ssEM_bin2',ssEM_bin2)
 call mpas_pool_get_array(SS2G_diags,'ssEM_bin3',ssEM_bin3)
 call mpas_pool_get_array(SS2G_diags,'ssEM_bin4',ssEM_bin4)
 call mpas_pool_get_array(SS2G_diags,'ssEM_bin5',ssEM_bin5)

 call mpas_pool_get_array(SS2G_diags,'ssSD_bin1',ssSD_bin1)
 call mpas_pool_get_array(SS2G_diags,'ssSD_bin2',ssSD_bin2)
 call mpas_pool_get_array(SS2G_diags,'ssSD_bin3',ssSD_bin3)
 call mpas_pool_get_array(SS2G_diags,'ssSD_bin4',ssSD_bin4)
 call mpas_pool_get_array(SS2G_diags,'ssSD_bin5',ssSD_bin5)

 call mpas_pool_get_array(SS2G_diags,'ssDP_bin1',ssDP_bin1)
 call mpas_pool_get_array(SS2G_diags,'ssDP_bin2',ssDP_bin2)
 call mpas_pool_get_array(SS2G_diags,'ssDP_bin3',ssDP_bin3)
 call mpas_pool_get_array(SS2G_diags,'ssDP_bin4',ssDP_bin4)
 call mpas_pool_get_array(SS2G_diags,'ssDP_bin5',ssDP_bin5)

 call mpas_pool_get_array(SS2G_diags,'ssWT_bin1',ssWT_bin1)
 call mpas_pool_get_array(SS2G_diags,'ssWT_bin2',ssWT_bin2)
 call mpas_pool_get_array(SS2G_diags,'ssWT_bin3',ssWT_bin3)
 call mpas_pool_get_array(SS2G_diags,'ssWT_bin4',ssWT_bin4)
 call mpas_pool_get_array(SS2G_diags,'ssWT_bin5',ssWT_bin5)

 call mpas_pool_get_array(SS2G_diags,'ssSV_bin1',ssSV_bin1)
 call mpas_pool_get_array(SS2G_diags,'ssSV_bin2',ssSV_bin2)
 call mpas_pool_get_array(SS2G_diags,'ssSV_bin3',ssSV_bin3)
 call mpas_pool_get_array(SS2G_diags,'ssSV_bin4',ssSV_bin4)
 call mpas_pool_get_array(SS2G_diags,'ssSV_bin5',ssSV_bin5)

 call mpas_pool_get_array(SS2G_diags,'ssSMASS'  ,ssSMASS  )
 call mpas_pool_get_array(SS2G_diags,'ssCMASS'  ,ssCMASS  )
 call mpas_pool_get_array(SS2G_diags,'ssSMASS25',ssSMASS25)
 call mpas_pool_get_array(SS2G_diags,'ssCMASS25',ssCMASS25)
 call mpas_pool_get_array(SS2G_diags,'ssFLUXU'  ,ssFLUXU  )
 call mpas_pool_get_array(SS2G_diags,'ssFLUXV'  ,ssFLUXV  )

 call mpas_pool_get_array(SS2G_diags,'ssMASS'  ,ssMASS  )
 call mpas_pool_get_array(SS2G_diags,'ssMASS25',ssMASS25)
 call mpas_pool_get_array(SS2G_diags,'ssCONC'  ,ssCONC  )

 do j = jts,jte
    do i = its,ite
       ssEM_bin1(i) = SS2G%ssem(i,j,1)
       ssEM_bin2(i) = SS2G%ssem(i,j,2)
       ssEM_bin3(i) = SS2G%ssem(i,j,3)
       ssEM_bin4(i) = SS2G%ssem(i,j,4)
       ssEM_bin5(i) = SS2G%ssem(i,j,5)
       ssSD_bin1(i) = SS2G%sssd(i,j,1)
       ssSD_bin2(i) = SS2G%sssd(i,j,2)
       ssSD_bin3(i) = SS2G%sssd(i,j,3)
       ssSD_bin4(i) = SS2G%sssd(i,j,4)
       ssSD_bin5(i) = SS2G%sssd(i,j,5)
       ssDP_bin1(i) = SS2G%ssdp(i,j,1)
       ssDP_bin2(i) = SS2G%ssdp(i,j,2)
       ssDP_bin3(i) = SS2G%ssdp(i,j,3)
       ssDP_bin4(i) = SS2G%ssdp(i,j,4)
       ssDP_bin5(i) = SS2G%ssdp(i,j,5)
       ssWT_bin1(i) = SS2G%sswt(i,j,1)
       ssWT_bin2(i) = SS2G%sswt(i,j,2)
       ssWT_bin3(i) = SS2G%sswt(i,j,3)
       ssWT_bin4(i) = SS2G%sswt(i,j,4)
       ssWT_bin5(i) = SS2G%sswt(i,j,5)
       ssSV_bin1(i) = SS2G%sssv(i,j,1)
       ssSV_bin2(i) = SS2G%sssv(i,j,2)
       ssSV_bin3(i) = SS2G%sssv(i,j,3)
       ssSV_bin4(i) = SS2G%sssv(i,j,4)
       ssSV_bin5(i) = SS2G%sssv(i,j,5)

       ssSMASS(i)   = SS2G%sssmass(i,j)
       ssCMASS(i)   = SS2G%sscmass(i,j)
       ssSMASS25(i) = SS2G%sssmass25(i,j)
       ssCMASS25(i) = SS2G%sscmass25(i,j)
       ssFLUXU(i)   = SS2G%ssfluxu(i,j)
       ssFLUXV(i)   = SS2G%ssfluxv(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          ssMASS(kk,i)   = SS2G%ssmass(i,j,k)
          ssMASS25(kk,i) = SS2G%ssmass25(i,j,k)
          ssCONC(kk,i)   = SS2G%ssconc(i,j,k)
       enddo
    enddo
 enddo


!--- sea-salt aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(SS2G_aops,'ssEXTTAU'  ,ssEXTTAU  )
 call mpas_pool_get_array(SS2G_aops,'ssSTEXTTAU',ssSTEXTTAU)
 call mpas_pool_get_array(SS2G_aops,'ssSCATAU'  ,ssSCATAU  )
 call mpas_pool_get_array(SS2G_aops,'ssSTSCATAU',ssSTSCATAU)
 call mpas_pool_get_array(SS2G_aops,'ssANGSTR'  ,ssANGSTR  )
 call mpas_pool_get_array(SS2G_aops,'ssAERIDX'  ,ssAERIDX  )
 call mpas_pool_get_array(SS2G_aops,'ssEXTT25'  ,ssEXTT25  )
 call mpas_pool_get_array(SS2G_aops,'ssSCAT25'  ,ssSCAT25  )
 call mpas_pool_get_array(SS2G_aops,'ssEXTTFM'  ,ssEXTTFM  )
 call mpas_pool_get_array(SS2G_aops,'ssSCATFM'  ,ssSCATFM  )

 call mpas_pool_get_array(SS2G_aops,'ssEXTCOEF'    ,ssEXTCOEF    )
 call mpas_pool_get_array(SS2G_aops,'ssEXTCOEFRH20',ssEXTCOEFRH20)
 call mpas_pool_get_array(SS2G_aops,'ssEXTCOEFRH80',ssEXTCOEFRH80)
 call mpas_pool_get_array(SS2G_aops,'ssSCACOEF'    ,ssSCACOEF    )
 call mpas_pool_get_array(SS2G_aops,'ssSCACOEFRH20',ssSCACOEFRH20)
 call mpas_pool_get_array(SS2G_aops,'ssSCACOEFRH80',ssSCACOEFRH80)
 call mpas_pool_get_array(SS2G_aops,'ssBCKCOEF'    ,ssBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       ssANGSTR(i) = ss2G%ssangstr(i,j)
       ssAERIDX(i) = ss2G%ssaeridx(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          ssEXTTAU(n,i)   = SS2G%ssexttau(i,j,n)
          ssSTEXTTAU(n,i) = SS2G%ssstexttau(i,j,n)
          ssSCATAU(n,i)   = SS2G%ssscatau(i,j,n)
          ssSTSCATAU(n,i) = SS2G%ssstscatau(i,j,n)
          ssEXTT25(n,i)   = SS2G%ssextt25(i,j,n)
          ssSCAT25(n,i)   = SS2G%ssscat25(i,j,n)
          ssEXTTFM(n,i)   = SS2G%ssexttfm(i,j,n)
          ssSCATFM(n,i)   = SS2G%ssscatfm(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             ssEXTCOEF(kk,n,i)     = SS2G%ssextcoef(i,j,k,n)
             ssEXTCOEFRH20(kk,n,i) = SS2G%ssextcoefrh20(i,j,k,n)
             ssEXTCOEFRH80(kk,n,i) = SS2G%ssextcoefrh80(i,j,k,n)
             ssSCACOEF(kk,n,i)     = SS2G%ssscacoef(i,j,k,n)
             ssSCACOEFRH20(kk,n,i) = SS2G%ssscacoefrh20(i,j,k,n)
             ssSCACOEFRH80(kk,n,i) = SS2G%ssscacoefrh80(i,j,k,n)
             ssBCKCOEF(kk,n,i)     = SS2G%ssbckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine SS2G_diagnostics.')

 end subroutine SS2G_diagnostics

!==================================================================================================================
 subroutine SU2G_diagnostics(mesh,SU2G,SU2G_diags,SU2G_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(SU2G_State),intent(in):: SU2G

!inout arguments:
 type(mpas_pool_type),intent(inout):: SU2G_diags
 type(mpas_pool_type),intent(inout):: SU2G_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: suEM_dms,suEM_so2,suEM_so4,suEM_msa
 real(kind=RKIND),dimension(:),pointer:: suSD_dms,suSD_so2,suSD_so4,suSD_msa
 real(kind=RKIND),dimension(:),pointer:: suDP_dms,suDP_so2,suDP_so4,suDP_msa
 real(kind=RKIND),dimension(:),pointer:: suWT_dms,suWT_so2,suWT_so4,suWT_msa
 real(kind=RKIND),dimension(:),pointer:: suSV_dms,suSV_so2,suSV_so4,suSV_msa
 real(kind=RKIND),dimension(:),pointer:: SO4eman,SO2eman,SO2embb,SO2emvn,SO2emve
 real(kind=RKIND),dimension(:),pointer:: suPSO2,suPMSA,suPSO4,suPSO4g,suPSO4aq,suPSO4wt
 real(kind=RKIND),dimension(:),pointer:: SO2smass,SO2cmass,SO4smass,SO4cmass,DMSsmass,DMScmass,MSAsmass,MSAcmass
 real(kind=RKIND),dimension(:,:),pointer:: SO4mass,SO4sarea,SO4snum
 real(kind=RKIND),dimension(:,:),pointer:: pSO2,pMSA,pSO4,pSO4g,pSO4aq,pSO4wt

 real(kind=RKIND),dimension(:),pointer:: suANGSTR
 real(kind=RKIND),dimension(:,:),pointer:: suEXTTAU,suSTEXTTAU,suSCATAU,suSTSCATAU
 real(kind=RKIND),dimension(:,:,:),pointer:: suEXTCOEF,suEXTCOEFRH20,suEXTCOEFRH80,suSCACOEF,suSCACOEFRH20, &
                                             suSCACOEFRH80,suBCKCOEF
!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine SU2G_diagnostics:')

 call mpas_pool_get_array(SU2G_diags,'suEM_dms',suEM_dms)
 call mpas_pool_get_array(SU2G_diags,'suEM_so2',suEM_so2)
 call mpas_pool_get_array(SU2G_diags,'suEM_so4',suEM_so4)
 call mpas_pool_get_array(SU2G_diags,'suEM_msa',suEM_msa)

 call mpas_pool_get_array(SU2G_diags,'suSD_dms',suSD_dms)
 call mpas_pool_get_array(SU2G_diags,'suSD_so2',suSD_so2)
 call mpas_pool_get_array(SU2G_diags,'suSD_so4',suSD_so4)
 call mpas_pool_get_array(SU2G_diags,'suSD_msa',suSD_msa)

 call mpas_pool_get_array(SU2G_diags,'suDP_dms',suDP_dms)
 call mpas_pool_get_array(SU2G_diags,'suDP_so2',suDP_so2)
 call mpas_pool_get_array(SU2G_diags,'suDP_so4',suDP_so4)
 call mpas_pool_get_array(SU2G_diags,'suDP_msa',suDP_msa)

 call mpas_pool_get_array(SU2G_diags,'suWT_dms',suWT_dms)
 call mpas_pool_get_array(SU2G_diags,'suWT_so2',suWT_so2)
 call mpas_pool_get_array(SU2G_diags,'suWT_so4',suWT_so4)
 call mpas_pool_get_array(SU2G_diags,'suWT_msa',suWT_msa)

 call mpas_pool_get_array(SU2G_diags,'suSV_dms',suSV_dms)
 call mpas_pool_get_array(SU2G_diags,'suSV_so2',suSV_so2)
 call mpas_pool_get_array(SU2G_diags,'suSV_so4',suSV_so4)
 call mpas_pool_get_array(SU2G_diags,'suSV_msa',suSV_msa)

 call mpas_pool_get_array(SU2G_diags,'suPSO2'  ,suPSO2  )
 call mpas_pool_get_array(SU2G_diags,'suPMSA'  ,suPMSA  )
 call mpas_pool_get_array(SU2G_diags,'suPSO4'  ,suPSO4  )
 call mpas_pool_get_array(SU2G_diags,'suPSO4g' ,suPSO4g )
 call mpas_pool_get_array(SU2G_diags,'suPSO4aq',suPSO4aq)
 call mpas_pool_get_array(SU2G_diags,'suPSO4wt',suPSO4wt)

 call mpas_pool_get_array(SU2G_diags,'SO4eman',SO4eman)
 call mpas_pool_get_array(SU2G_diags,'SO2eman',SO2eman)
 call mpas_pool_get_array(SU2G_diags,'SO2embb',SO2embb)
 call mpas_pool_get_array(SU2G_diags,'SO2emvn',SO2emvn)
 call mpas_pool_get_array(SU2G_diags,'SO2emve',SO2emve)

 call mpas_pool_get_array(SU2G_diags,'SO2smass',SO2smass)
 call mpas_pool_get_array(SU2G_diags,'SO2cmass',SO2cmass)
 call mpas_pool_get_array(SU2G_diags,'SO4smass',SO4smass)
 call mpas_pool_get_array(SU2G_diags,'SO4cmass',SO4cmass)
 call mpas_pool_get_array(SU2G_diags,'DMSsmass',DMSsmass)
 call mpas_pool_get_array(SU2G_diags,'DMScmass',DMScmass)
 call mpas_pool_get_array(SU2G_diags,'MSAsmass',MSAsmass)
 call mpas_pool_get_array(SU2G_diags,'MSAcmass',MSAcmass)
 call mpas_pool_get_array(SU2G_diags,'SO4mass' ,SO4mass )
 call mpas_pool_get_array(SU2G_diags,'SO4snum' ,SO4snum )
 call mpas_pool_get_array(SU2G_diags,'SO4sarea',SO4sarea)

 call mpas_pool_get_array(SU2G_diags,'pSO2'  ,pSO2  )
 call mpas_pool_get_array(SU2G_diags,'pMSA'  ,pMSA  )
 call mpas_pool_get_array(SU2G_diags,'pSO4'  ,pSO4  )
 call mpas_pool_get_array(SU2G_diags,'pSO4g' ,pSO4g )
 call mpas_pool_get_array(SU2G_diags,'pSO4aq',pSO4aq)
 call mpas_pool_get_array(SU2G_diags,'pSO4wt',pSO4wt)

 do j = jts,jte
    do i = its,ite
       suEM_dms(i) = SU2G%suem(i,j,1)
       suEM_so2(i) = SU2G%suem(i,j,2)
       suEM_so4(i) = SU2G%suem(i,j,3)
       suEM_msa(i) = SU2G%suem(i,j,4)

       suSD_dms(i) = SU2G%susd(i,j,1)
       suSD_so2(i) = SU2G%susd(i,j,2)
       suSD_so4(i) = SU2G%susd(i,j,3)
       suSD_msa(i) = SU2G%susd(i,j,4)

       suDP_dms(i) = SU2G%sudp(i,j,1)
       suDP_so2(i) = SU2G%sudp(i,j,2)
       suDP_so4(i) = SU2G%sudp(i,j,3)
       suDP_msa(i) = SU2G%sudp(i,j,4)

       suWT_dms(i) = SU2G%suwt(i,j,1)
       suWT_so2(i) = SU2G%suwt(i,j,2)
       suWT_so4(i) = SU2G%suwt(i,j,3)
       suWT_msa(i) = SU2G%suwt(i,j,4)

       suSV_dms(i) = SU2G%susv(i,j,1)
       suSV_so2(i) = SU2G%susv(i,j,2)
       suSV_so4(i) = SU2G%susv(i,j,3)
       suSV_msa(i) = SU2G%susv(i,j,4)

       suPSO2(i)   = SU2G%supso2(i,j)
       suPMSA(i)   = SU2G%supmsa(i,j)
       suPSO4(i)   = SU2G%supso4(i,j)
       suPSO4g(i)  = SU2G%supso4g(i,j)
       suPSO4aq(i) = SU2G%supso4aq(i,j)
       suPSO4wt(i) = SU2G%supso4wt(i,j)

       SO4eman(i)  = SU2G%so4eman(i,j)
       SO2eman(i)  = SU2G%so2eman(i,j)
       SO2embb(i)  = SU2G%so2embb(i,j)
       SO2emvn(i)  = SU2G%so2emvn(i,j)
       SO2emve(i)  = SU2G%so2emve(i,j)

       SO2smass(i) = SU2G%so2smass(i,j)
       SO2cmass(i) = SU2G%so2cmass(i,j)
       SO4smass(i) = SU2G%so4smass(i,j)
       SO4cmass(i) = SU2G%so4cmass(i,j)
       DMSsmass(i) = SU2G%dmssmass(i,j)
       DMScmass(i) = SU2G%dmscmass(i,j)
       MSAsmass(i) = SU2G%msasmass(i,j)
       MSAcmass(i) = SU2G%msacmass(i,j)
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          SO4mass(kk,i)  = SU2G%so4mass(i,j,k)
          SO4snum(kk,i)  = SU2G%so4snum(i,j,k)
          SO4sarea(kk,i) = SU2G%so4sarea(i,j,k)
          pSO2(kk,i)     = SU2G%pso2(i,j,k)
          pMSA(kk,i)     = SU2G%pmsa(i,j,k)
          pSO4(kk,i)     = SU2G%pso4(i,j,k)
          pSO4g(kk,i)    = SU2G%pso4g(i,j,k)
          pSO4aq(kk,i)   = SU2G%pso4aq(i,j,k)
          pSO4wt(kk,i)   = SU2G%pso4wet(i,j,k)
       enddo
    enddo
 enddo


!--- sulfate aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(SU2G_aops,'suEXTTAU'  ,suEXTTAU  )
 call mpas_pool_get_array(SU2G_aops,'suSTEXTTAU',suSTEXTTAU)
 call mpas_pool_get_array(SU2G_aops,'suSCATAU'  ,suSCATAU  )
 call mpas_pool_get_array(SU2G_aops,'suSTSCATAU',suSTSCATAU)
 call mpas_pool_get_array(SU2G_aops,'suANGSTR'  ,suANGSTR  )

 call mpas_pool_get_array(SU2G_aops,'suEXTCOEF'    ,suEXTCOEF    )
 call mpas_pool_get_array(SU2G_aops,'suEXTCOEFRH20',suEXTCOEFRH20)
 call mpas_pool_get_array(SU2G_aops,'suEXTCOEFRH80',suEXTCOEFRH80)
 call mpas_pool_get_array(SU2G_aops,'suSCACOEF'    ,suSCACOEF    )
 call mpas_pool_get_array(SU2G_aops,'suSCACOEFRH20',suSCACOEFRH20)
 call mpas_pool_get_array(SU2G_aops,'suSCACOEFRH80',suSCACOEFRH80)
 call mpas_pool_get_array(SU2G_aops,'suBCKCOEF'    ,suBCKCOEF    )

 do j = jts,jte
    do i = its,ite
       suANGSTR(i) = SU2G%suangstr(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          suEXTTAU(n,i)   = SU2G%suexttau(i,j,n)
          suSTEXTTAU(n,i) = SU2G%sustexttau(i,j,n)
          suSCATAU(n,i)   = SU2G%suscatau(i,j,n)
          suSTSCATAU(n,i) = SU2G%sustscatau(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             suEXTCOEF(kk,n,i)     = SU2G%suextcoef(i,j,k,n)
             suEXTCOEFRH20(kk,n,i) = SU2G%suextcoefrh20(i,j,k,n)
             suEXTCOEFRH80(kk,n,i) = SU2G%suextcoefrh80(i,j,k,n)
             suSCACOEF(kk,n,i)     = SU2G%suscacoef(i,j,k,n)
             suSCACOEFRH20(kk,n,i) = SU2G%suscacoefrh20(i,j,k,n)
             suSCACOEFRH80(kk,n,i) = SU2G%suscacoefrh80(i,j,k,n)
             suBCKCOEF(kk,n,i)     = SU2G%subckcoef(i,j,k,n)
          enddo
       enddo
    enddo
 enddo

 
 call mpas_log_write('--- end subroutine SU2G_diagnostics.')

 end subroutine SU2G_diagnostics

!==================================================================================================================
 subroutine GOCART2G_diagnostics(mesh,GOCART2G,GOCART2G_diags,GOCART2G_aops,its,ite,jts,jte,kts,kte)
!==================================================================================================================

!input arguments:
 integer:: its,ite,jts,jte,kts,kte
 type(mpas_pool_type),intent(in):: mesh
 type(GOCART2G_State),intent(in):: GOCART2G

!inout arguments:
 type(mpas_pool_type),intent(inout):: GOCART2G_diags
 type(mpas_pool_type),intent(inout):: GOCART2G_aops

!local arguments and arrays:
 integer:: i,j,k,kk,n
 integer,pointer:: npAOPs,nvAOPs

 real(kind=RKIND),dimension(:),pointer:: totPM,totPMRH35,totPMRH50
 real(kind=RKIND),dimension(:),pointer:: totPM25,totPM25RH35,totPM25RH50

 real(kind=RKIND),dimension(:),pointer:: totANGSTR
 real(kind=RKIND),dimension(:,:),pointer:: totEXTTAU,totSTEXTTAU,totEXTT25,totEXTTFM
 real(kind=RKIND),dimension(:,:),pointer:: totSCATAU,totSTSCATAU,totSCAT25,totSCATFM
 real(kind=RKIND),dimension(:,:),pointer:: totABCKTOA,totABCKSFC
 real(kind=RKIND),dimension(:,:,:),pointer:: totEXTCOEF,totEXTCOEFRH20,totEXTCOEFRH80,totSCACOEF,totSCACOEFRH20, &
                                             totSCACOEFRH80,totBCKCOEF
!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine GOCART2G_diagnostics:')

 call mpas_pool_get_array(GOCART2G_diags,'totPM'      ,totPM      )
 call mpas_pool_get_array(GOCART2G_diags,'totPMRH35'  ,totPMRH35  )
 call mpas_pool_get_array(GOCART2G_diags,'totPMRH50'  ,totPMRH50  )
 call mpas_pool_get_array(GOCART2G_diags,'totPM25'    ,totPM25    )
 call mpas_pool_get_array(GOCART2G_diags,'totPM25RH35',totPM25RH35)
 call mpas_pool_get_array(GOCART2G_diags,'totPM25RH50',totPM25RH50)

 do j = jts,jte
    do i = its,ite
       totPM(i)       = GOCART2G%pm(i,j)
       totPMRH35(i)   = GOCART2G%pm_rh35(i,j)
       totPMRH50(i)   = GOCART2G%pm_rh50(i,j)
       totPM25(i)     = GOCART2G%pm(i,j)
       totPM25RH35(i) = GOCART2G%pm25_rh35(i,j)
       totPM25RH50(i) = GOCART2G%pm25_rh50(i,j)
    enddo
 enddo


!--- total aerosol optical properties:
 call mpas_pool_get_dimension(mesh,'npAOPs',npAOPs)
 call mpas_pool_get_dimension(mesh,'nvAOPs',nvAOPs)

 call mpas_pool_get_array(GOCART2G_aops,'totANGSTR'  ,totANGSTR  )

 call mpas_pool_get_array(GOCART2G_aops,'totEXTTAU'  ,totEXTTAU  )
 call mpas_pool_get_array(GOCART2G_aops,'totSTEXTTAU',totSTEXTTAU)
 call mpas_pool_get_array(GOCART2G_aops,'totEXTT25'  ,totEXTT25  )
 call mpas_pool_get_array(GOCART2G_aops,'totEXTTFM'  ,totEXTTFM  )
 call mpas_pool_get_array(GOCART2G_aops,'totSCATAU'  ,totSCATAU  )
 call mpas_pool_get_array(GOCART2G_aops,'totSTSCATAU',totSTSCATAU)
 call mpas_pool_get_array(GOCART2G_aops,'totSCAT25'  ,totSCAT25  )
 call mpas_pool_get_array(GOCART2G_aops,'totSCATFM'  ,totSCATFM  )

 call mpas_pool_get_array(GOCART2G_aops,'totEXTCOEF'    ,totEXTCOEF    )
 call mpas_pool_get_array(GOCART2G_aops,'totEXTCOEFRH20',totEXTCOEFRH20)
 call mpas_pool_get_array(GOCART2G_aops,'totEXTCOEFRH80',totEXTCOEFRH80)
 call mpas_pool_get_array(GOCART2G_aops,'totSCACOEF'    ,totSCACOEF    )
 call mpas_pool_get_array(GOCART2G_aops,'totSCACOEFRH20',totSCACOEFRH20)
 call mpas_pool_get_array(GOCART2G_aops,'totSCACOEFRH80',totSCACOEFRH80)
 call mpas_pool_get_array(GOCART2G_aops,'totBCKCOEF'    ,totBCKCOEF    )
 call mpas_pool_get_array(GOCART2G_aops,'totABCKTOA'    ,totABCKTOA    )
 call mpas_pool_get_array(GOCART2G_aops,'totABCKSFC'    ,totABCKSFC    )


 do j = jts,jte
    do i = its,ite
       totANGSTR(i) = GOCART2G%totangstr(i,j)
    enddo

    do n = 1,nvAOPs
       do i = its,ite
          totEXTTAU(n,i)   = GOCART2G%totexttau(i,j,n)
          totSTEXTTAU(n,i) = GOCART2G%totstexttau(i,j,n)
          totEXTT25(n,i)   = GOCART2G%totextt25(i,j,n)
          totEXTTFM(n,i)   = GOCART2G%totexttfm(i,j,n)
          totSCATAU(n,i)   = GOCART2G%totscatau(i,j,n)
          totSTSCATAU(n,i) = GOCART2G%totstscatau(i,j,n)
          totSCAT25(n,i)   = GOCART2G%totscat25(i,j,n)
          totSCATFM(n,i)   = GOCART2G%totscatfm(i,j,n)
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do n = 1,npAOPs
          do i = its,ite
             totEXTCOEF(kk,n,i)     = GOCART2G%totextcoef(i,j,k,n)
             totEXTCOEFRH20(kk,n,i) = GOCART2G%totextcoefrh20(i,j,k,n)
             totEXTCOEFRH80(kk,n,i) = GOCART2G%totextcoefrh80(i,j,k,n)
             totSCACOEF(kk,n,i)     = GOCART2G%totscacoef(i,j,k,n)
             totSCACOEFRH20(kk,n,i) = GOCART2G%totscacoefrh20(i,j,k,n)
             totSCACOEFRH80(kk,n,i) = GOCART2G%totscacoefrh80(i,j,k,n)
             totBCKCOEF(kk,n,i)     = GOCART2G%totbckcoef(i,j,k,n)
          enddo
       enddo
    enddo

    do k = kts,kte
       kk = kte+1-k
       do i = its,ite
          totABCKTOA(kk,i) = GOCART2G%totabcktoa(i,j,k)
          totABCKSFC(kk,i) = GOCART2G%totabcksfc(i,j,k)
       enddo
    enddo
 enddo


 call mpas_log_write('--- end subroutine GOCART2G_diagnostics.')

 end subroutine GOCART2G_diagnostics

!==================================================================================================================
 end module mpas_chemistry_gocart2G_diagnostics
!==================================================================================================================
