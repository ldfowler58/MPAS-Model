!=================================================================================================================
 module module_gocart_aerosols
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND
 use mpas_log


 implicit none
 private
 public:: gocart_carbon_hydrophobic_to_hydrophilic


 contains


!=================================================================================================================
 subroutine gocart_carbon_hydrophobic_to_hydrophilic(its,ite,kts,kte,ktep1,dt,qbc1,qbc2,qoc1,qoc2,rho,area,zgrid)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,kts,kte,ktep1

 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in),dimension(its:ite):: area
 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: rho
 real(kind=RKIND),intent(in),dimension(its:ite,kts:ktep1):: zgrid


!inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: qbc1,qbc2,qoc1,qoc2


!local variables and arrays:
 integer,parameter:: &
    imx = 1,  &
    lmx = 1,  &
    nmx = 4
 integer,parameter:: &
    nbc1 = 1, &
    noc1 = 2, &
    nbc2 = 3, &
    noc2 = 4
 integer:: i,k,n,ndt1

 real(kind=RKIND),parameter:: mw_c   = 12.
 real(kind=RKIND),parameter:: mw_dry = 28.966 ! molecular weight of dry air (g/mole) 

 real(kind=R8KIND):: dz,tt1,tt2
 real(kind=R8KIND),dimension(imx,4):: bchmlos
 real(kind=R8KIND),dimension(imx,lmx):: airmas
 real(kind=R8KIND),dimension(imx,lmx,4):: tc
 real(kind=R8KIND),dimension(imx,lmx,2):: pc2
 real(kind=R8KIND),dimension(imx,lmx,4):: chmlos

 real(kind=RKIND):: tem1,tem2,tem3,tem4

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine gocart_carbon_hydrophobic_to_hydrophilic:')

 ndt1 = ifix(dt)

 do i = its,ite
    do k = kts,kte-1
       dz = zgrid(i,k+1)-zgrid(i,k)
       dz = max(0.,dz)
       airmas(imx,lmx) = rho(i,k)*dz*area(i)

       pc2(imx,lmx,1) = 0.
       pc2(imx,lmx,2) = 0.
       do n = 1,nmx
          chmlos(imx,lmx,n) = 0.
       enddo

!      tc(imx,lmx,1) = qbc1(i,k)/mw_c*mw_dry*1.d-9
!      tc(imx,lmx,2) = qoc1(i,k)/mw_c*mw_dry*1.d-9
!      tc(imx,lmx,3) = qbc2(i,k)/mw_c*mw_dry*1.d-9
!      tc(imx,lmx,4) = qoc2(i,k)/mw_c*mw_dry*1.d-9
       tc(imx,lmx,1) = qbc1(i,k)/mw_c*mw_dry
       tc(imx,lmx,2) = qoc1(i,k)/mw_c*mw_dry
       tc(imx,lmx,3) = qbc2(i,k)/mw_c*mw_dry
       tc(imx,lmx,4) = qoc2(i,k)/mw_c*mw_dry
       tt1 = tc(imx,lmx,3)

       tem1 = qbc1(i,k)
       tem2 = qoc1(i,k)
       tem3 = qbc2(i,k)
       tem4 = qoc2(i,k)

       call chem_1(imx,lmx,nmx,nbc1,noc1,ndt1,airmas,tc,chmlos,bchmlos,pc2)
       call chem_2(imx,lmx,nmx,nbc2,noc2,ndt1,airmas,tc,pc2)

       !how much bc is aged, oc production (since no soa) is made dependent on that (?):
       tt2 = tc(imx,lmx,3)-tt1
       qbc1(i,k) = tc(imx,lmx,1)/mw_dry*mw_c
       qoc1(i,k) = tc(imx,lmx,2)/mw_dry*mw_c
       qbc2(i,k) = tc(imx,lmx,3)/mw_dry*mw_c
       qoc2(i,k) = (tc(imx,lmx,4)+8.*tt2)/mw_dry*mw_c
!      qbc1(i,k) = tc(imx,lmx,1)/mw_dry*mw_c*1.e9
!      qoc1(i,k) = tc(imx,lmx,2)/mw_dry*mw_c*1.e9
!      qbc2(i,k) = tc(imx,lmx,3)/mw_dry*mw_c*1.e9
!      qoc2(i,k) = (tc(imx,lmx,4)+8.*tt2)/mw_dry*mw_c*1.e9

       tem1 = qbc1(i,k) - tem1
       tem2 = qoc1(i,k) - tem2
       tem3 = qbc2(i,k) - tem3
       tem4 = qoc2(i,k) - tem4
    enddo
 enddo
 
 call mpas_log_write('--- end subroutine gocart_carbon_hydrophobic_to_hydrophilic:')

 end subroutine gocart_carbon_hydrophobic_to_hydrophilic

!=================================================================================================================
 subroutine chem_1(imx,lmx,nmx,nbc1,noc1,ndt1,airm,tc,chmlos,bchmlos,pc2)
!=================================================================================================================

!input arguments:
 integer, intent(in):: imx,lmx,nmx,nbc1,noc1,ndt1
 real(kind=R8KIND),intent(in):: airm(imx,lmx)

!inout arguments:
 real(kind=R8KIND),intent(inout):: tc(imx,lmx,nmx)
 real(kind=R8KIND),intent(inout):: chmlos(imx,lmx,nmx)
 real(kind=R8KIND),intent(inout):: bchmlos(imx,nmx)

!output arguments:
 real(kind=R8KIND),intent(out):: pc2(imx,lmx,2)

!local variables and arrays:
 integer:: np,n,i,l
 real(kind=R8KIND):: r1,c0,r2,rkt,c1

!-----------------------------------------------------------------------------------------------------------------

 r1 = 4.63d-6
  
 do n = 1,nmx
    if(n == nbc1 .or. n == noc1) then
       if (n == nbc1) np = 1
       if (n == noc1) np = 2

       do l = 1,lmx
          do i = 1,imx
             c0 = tc(i,l,n)
             r2 = 0.0 ! used to be loss due to dry dep
             rkt = (r1+r2)*real(ndt1)

             c1 = c0* exp(-rkt)
             c1 = max(c1,1.0d-32)
             tc(i,l,n) = c1
             pc2(i,l,np) = (c0-c1)*r1/(r1+r2)

             !diagnostics:  
             !chmlos(i,l,np) = chmlos(i,l,np) + pc2(i,l,n)*airm(i,l)
             chmlos(i,l,n) = chmlos(i,l,n) + pc2(i,l,np)*airm(i,l)
          enddo
       enddo

       do i = 1,imx
          bchmlos(i,n) = 0.
          do l = 1,lmx
             bchmlos(i,n) = bchmlos(i,n) + chmlos(i,l,n)
          enddo
       enddo
    endif
 enddo
  
 end subroutine chem_1

!=================================================================================================================
 subroutine chem_2(imx,lmx,nmx,nbc2,noc2,ndt1,airm,tc,pc2)
!=================================================================================================================

!input arguments:
 integer,intent(in):: lmx,imx,nmx,nbc2,noc2,ndt1
 real(kind=R8KIND),intent(in),dimension(imx,lmx):: airm
 real(kind=R8KIND),intent(in),dimension(imx,lmx,2):: pc2

!inout arguments:
 real(kind=R8KIND),intent(inout),dimension(imx,lmx,nmx):: tc

!local variables:
 integer :: np,n,i,l
 real(kind=R8KIND):: c0,pp,rkt,c1

!-----------------------------------------------------------------------------------------------------------------
  
 do n = 1,nmx
    if(n == nbc2 .or. n == noc2) then
       if(n == nbc2) np = 1
       if(n == noc2) np = 2

       do l = 1,lmx
          do i = 1,imx
             c0 = tc(i,l,n)
             !pp = pc2(i,l,n)
             pp = pc2(i,l,np)
             c1 = c0 + pp
                 
             c1 = max(c1,1.0d-32)
             tc(i,l,n) = c1
          enddo
       enddo
    endif
 enddo

 end subroutine chem_2

!=================================================================================================================
 end module module_gocart_aerosols
!=================================================================================================================

