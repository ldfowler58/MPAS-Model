!=================================================================================================================
 module module_cu_ntiedtke
 use mpas_log
 use mpas_kind_types,only: RKIND,StrKIND

 use cu_ntiedtke,only: cu_ntiedtke_run,           &
                       cu_ntiedtke_init,          &
                       cu_ntiedtke_timestep_init, &
                       cu_ntiedtke_timestep_final
 use cu_ntiedtke_common

 implicit none
 private
 public:: cu_ntiedtke_driver, &
          ntiedtkeinit


 contains


!=================================================================================================================
 subroutine cu_ntiedtke_driver(                          &
                 dt,itimestep,stepcu                     &
                ,raincv,pratec,qfx,hfx                   &
                ,u3d,v3d,w,t3d,qv3d,qc3d,qi3d,pi3d,rho3d &
                ,qvften,thften                           &
                ,dz8w,pcps,p8w,xland,cu_act_flag,dx      &
                ,f_qv,f_qc,f_qr,f_qi,f_qs                &
                ,grav,xlf,xls,xlv,rd,rv,cp               &
                ,rthcuten,rqvcuten,rqccuten,rqicuten     &
                ,rucuten,rvcuten                         &
                ,ids,ide,jds,jde,kds,kde                 &
                ,ims,ime,jms,jme,kms,kme                 &
                ,its,ite,jts,jte,kts,kte                 &
                ,errmsg,errflg)
!=================================================================================================================
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- qc3d        3d cloud mixing ratio (kg/kg)
!-- qi3d        3d ice mixing ratio (kg/kg)
!-- rho3d       3d air density (kg/m^3)
!-- p8w         3d hydrostatic pressure at full levels (pa)
!-- pcps        3d hydrostatic pressure at half levels (pa)
!-- pi3d        3d exner function (dimensionless)
!-- qvften      3d total advective + PBL moisture tendency (kg kg-1 s-1)
!-- thften      3d total advective + PBL + radiative temperature tendency (k s-1)
!-- rthcuten      theta tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rucuten       u wind tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rvcuten       v wind tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rqvcuten      qv tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqccuten      qc tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqicuten      qi tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rainc         accumulated total cumulus scheme precipitation (mm)
!-- raincv        cumulus scheme precipitation (mm)
!-- pratec        precipitiation rate from cumulus scheme (mm/s)
!-- dz8w        dz between full levels (m)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- hfx         upward heat flux at the surface (w/m^2) 
!-- dt          time step (s)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-----------------------------------------------------------------------------------------------------------------

!--- input arguments:
 logical,intent(in),optional:: f_qv,f_qc,f_qr,f_qi,f_qs

 integer,intent(in):: ids,ide,jds,jde,kds,kde, &
                      ims,ime,jms,jme,kms,kme, &
                      its,ite,jts,jte,kts,kte

 integer,intent(in):: itimestep,stepcu

 real(kind=RKIND),intent(in):: cp,grav,rd,rv,xlf,xls,xlv

 real(kind=kind_phys),intent(in):: dt

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: dx,hfx,qfx,xland

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz8w,   &
    pcps,   &
    p8w,    &
    pi3d,   &
    qc3d,   &
    qvften, &
    thften, &
    qi3d,   &
    qv3d,   &
    rho3d,  &
    t3d,    &
    u3d,    &
    v3d,    &
    w

!--- inout arguments:
 logical,intent(inout),dimension(ims:ime,jms:jme):: cu_act_flag

 real(kind=kind_phys),intent(inout),dimension(ims:ime,jms:jme):: raincv, pratec

 real(kind=kind_phys),intent(inout),dimension(ims:ime,kms:kme,jms:jme),optional:: &
 rqccuten,  &
 rqicuten,  &
 rqvcuten,  &
 rthcuten,  &
 rucuten,   &
 rvcuten

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!--- local variables and arrays:
 integer:: i,im,j,k,kx,kx1
 integer,dimension(its:ite)::  slimsk

 real(kind=kind_phys):: delt
 real(kind=kind_phys),dimension(its:ite):: rn
 real(kind=kind_phys),dimension(its:ite,kts:kte):: prsl,omg,ghtl
 real(kind=kind_phys),dimension(its:ite,kts:kte):: uf,vf,tf,qvf,qcf,qif
 real(kind=kind_phys),dimension(its:ite,kts:kte):: qvftenz,thftenz
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: prsi,ghti,zi

 real(kind=kind_phys),dimension(its:ite):: dx_hv,hfx_hv,qfx_hv,xland_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: dz_hv,pi_hv,prsl_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: qv_hv,qc_hv,qi_hv,rho_hv,t_hv,u_hv,v_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: qvften_hv,thften_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: prsi_hv,w_hv

 real(kind=kind_phys),dimension(its:ite):: raincv_hv,pratec_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: rthcuten_hv,rqvcuten_hv,rqccuten_hv,rqicuten_hv, &
                                                   rucuten_hv,rvcuten_hv

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine cu_ntiedtke:')

 errmsg = ' '
 errflg = 0

 call cu_ntiedtke_init(cp,rd,rv,xlv,xls,xlf,grav,errmsg,errflg)

 call mpas_log_write('--- rd  = $r',realArgs=(/rd/))
 call mpas_log_write('--- rv  = $r',realArgs=(/rv/))
 call mpas_log_write('--- cpd = $r',realArgs=(/cp/))
 call mpas_log_write('--- alv = $r',realArgs=(/xlv/))
 call mpas_log_write('--- als = $r',realArgs=(/xls/))
 call mpas_log_write('--- alf = $r',realArgs=(/xlf/))
 call mpas_log_write('--- g   = $r',realArgs=(/grav/))
 call mpas_log_write(' ')
 call mpas_log_write('--- t13     = $r',realArgs=(/t13/)    )
 call mpas_log_write('--- rcpd    = $r',realArgs=(/rcpd/)   )
 call mpas_log_write('--- tmelt   = $r',realArgs=(/tmelt/)  )
 call mpas_log_write('--- zrg     = $r',realArgs=(/zrg/)    )
 call mpas_log_write('--- c1es    = $r',realArgs=(/c1es/)   )
 call mpas_log_write('--- c2es    = $r',realArgs=(/c2es/)   )
 call mpas_log_write('--- c3les   = $r',realArgs=(/c3les/)  )
 call mpas_log_write('--- c3ies   = $r',realArgs=(/c3ies/)  )
 call mpas_log_write('--- c4les   = $r',realArgs=(/c4les/)  )
 call mpas_log_write('--- c4ies   = $r',realArgs=(/c4ies/)  )
 call mpas_log_write('--- c5les   = $r',realArgs=(/c5les/)  )
 call mpas_log_write('--- c5ies   = $r',realArgs=(/c5ies/)  )
 call mpas_log_write('--- r5alvcp = $r',realArgs=(/r5alvcp/))
 call mpas_log_write('--- r5alscp = $r',realArgs=(/r5alscp/))
 call mpas_log_write('--- ralvdcp = $r',realArgs=(/ralvdcp/))
 call mpas_log_write('--- ralsdcp = $r',realArgs=(/ralsdcp/))
 call mpas_log_write('--- ralfdcp = $r',realArgs=(/ralfdcp/))
 call mpas_log_write('--- rtwat   = $r',realArgs=(/rtwat/)  )
 call mpas_log_write('--- rtber   = $r',realArgs=(/rtber/)  )
 call mpas_log_write('--- rtice   = $r',realArgs=(/rtice/)  )
 call mpas_log_write('--- vtmpc1  = $r',realArgs=(/vtmpc1/) )
 call mpas_log_write(' ')
 call mpas_log_write('--- momtrans = $i',intArgs=(/momtrans/))
 call mpas_log_write('--- entrdd   = $r',realArgs=(/entrdd/) )
 call mpas_log_write('--- cmfcmax  = $r',realArgs=(/cmfcmax/))
 call mpas_log_write('--- cmfcmin  = $r',realArgs=(/cmfcmin/))
 call mpas_log_write('--- cmfdeps  = $r',realArgs=(/cmfdeps/))
 call mpas_log_write('--- zdnoprc  = $r',realArgs=(/zdnoprc/))
 call mpas_log_write('--- cprcon   = $r',realArgs=(/cprcon/) )
 call mpas_log_write('--- pgcoef   = $r',realArgs=(/pgcoef/) )
 call mpas_log_write(' ')
 call mpas_log_write('--- nonequil = $l',logicArgs = (/nonequil/))
 call mpas_log_write('--- lmfpen   = $l',logicArgs = (/lmfpen/)  )
 call mpas_log_write('--- lmfmid   = $l',logicArgs = (/lmfmid/)  )
 call mpas_log_write('--- lmfscv   = $l',logicArgs = (/lmfscv/)  )
 call mpas_log_write('--- lmfdd    = $l',logicArgs = (/lmfdd/)   )
 call mpas_log_write('--- lmfdudv  = $l',logicArgs = (/lmfdudv/) )
 call mpas_log_write('--- end subroutine cu_ntiedtke:')
 call mpas_log_write(' ')

 do j = jts,jte
    do i = its,ite
       cu_act_flag(i,j)=.true.
    enddo
 enddo

 do j = jts,jte

    do i = its,ite
       dx_hv(i)    = dx(i,j)
       hfx_hv(i)   = hfx(i,j)
       qfx_hv(i)   = qfx(i,j)
       xland_hv(i) = xland(i,j)
    enddo

    do k = kts,kte
       do i = its,ite
          dz_hv(i,k)   = dz8w(i,k,j)
          pi_hv(i,k)   = pi3d(i,k,j)
          prsl_hv(i,k) = pcps(i,k,j)
          qv_hv(i,k)   = qv3d(i,k,j)
          qc_hv(i,k)   = qc3d(i,k,j)
          qi_hv(i,k)   = qi3d(i,k,j)
          rho_hv(i,k)  = rho3d(i,k,j)
          t_hv(i,k)    = t3d(i,k,j)
          u_hv(i,k)    = u3d(i,k,j)
          v_hv(i,k)    = v3d(i,k,j)

          qvften_hv(i,k) = qvften(i,k,j)
          thften_hv(i,k) = thften(i,k,j)
       enddo
    enddo
    do k = kts,kte+1
       do i = its,ite
          prsi_hv(i,k) = p8w(i,k,j)
          w_hv(i,k)    = w(i,k,j)
       enddo
    enddo

    call cu_ntiedtke_timestep_init(its,ite,kts,kte,im,kx,kx1,itimestep,stepcu,dt,grav,xland_hv,dz_hv, &
                    prsl_hv,prsi_hv,t_hv,rho_hv,qv_hv,qc_hv,qi_hv,u_hv,v_hv,w_hv,qvften_hv,thften_hv, &
                    qvftenz,thftenz,slimsk,delt,prsl,ghtl,tf,qvf,qcf,qif,uf,vf,prsi,ghti,omg,         &
                    errmsg,errflg)

    call cu_ntiedtke_run(uf,vf,tf,qvf,qcf,qif,qvftenz,thftenz,ghtl,ghti,omg,prsl,prsi,qfx_hv,hfx_hv, &
                         rn,slimsk,im,kx,kx1,delt,dx_hv,errmsg,errflg)

    call cu_ntiedtke_timestep_final(its,ite,kts,kte,stepcu,dt,pi_hv,qv_hv,qc_hv,qi_hv,t_hv,u_hv,v_hv, &
                     qvf,qcf,qif,tf,uf,vf,rn,raincv_hv,pratec_hv,rthcuten_hv,rqvcuten_hv,rqccuten_hv, &
                     rqicuten_hv,rucuten_hv,rvcuten_hv,errmsg,errflg)

    do i = its,ite
       raincv(i,j) = raincv_hv(i)
       pratec(i,j) = pratec_hv(i)
    enddo

    do k = kts,kte
       do i = its,ite
          rucuten(i,k,j)  = rucuten_hv(i,k)
          rvcuten(i,k,j)  = rvcuten_hv(i,k)
          rthcuten(i,k,j) = rthcuten_hv(i,k)
          rqvcuten(i,k,j) = rqvcuten_hv(i,k)
       enddo
    enddo

    if(present(rqccuten))then
       if(f_qc) then
          do k = kts,kte
             do i = its,ite
                rqccuten(i,k,j) = rqccuten_hv(i,k)
             enddo
          enddo
       endif
    endif

    if(present(rqicuten))then
       if(f_qi) then
          do k = kts,kte
             do i = its,ite
                rqicuten(i,k,j) = rqicuten_hv(i,k)
             enddo
          enddo
       endif
    endif

 enddo

 end subroutine cu_ntiedtke_driver

!=================================================================================================================
 subroutine ntiedtkeinit(rthcuten,rqvcuten,rqccuten,rqicuten, &
                         rucuten,rvcuten,rthften,rqvften,     &
                         restart,p_qc,p_qi,p_first_scalar,    &
                         allowed_to_read,                     &
                         ids, ide, jds, jde, kds, kde,        &
                         ims, ime, jms, jme, kms, kme,        &
                         its, ite, jts, jte, kts, kte)
!=================================================================================================================

!--- input arguments:
 logical,intent(in):: allowed_to_read,restart

 integer,intent(in):: ids, ide, jds, jde, kds, kde, &
                      ims, ime, jms, jme, kms, kme, &
                      its, ite, jts, jte, kts, kte
 integer,intent(in):: p_first_scalar,p_qi,p_qc

!--- output arguments:
 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme )::  &
    rthcuten,rqvcuten,rqccuten,rqicuten,rucuten,rvcuten,rthften,rqvften

!--- local variables and arrays:
 integer:: i,j,k,itf,jtf,ktf

!-----------------------------------------------------------------------------------------------------------------

 jtf = min0(jte,jde-1)
 ktf = min0(kte,kde-1)
 itf = min0(ite,ide-1)

 if(.not.restart)then
    do j = jts,jtf
       do k = kts,ktf
          do i = its,itf
             rthcuten(i,k,j) = 0.
             rqvcuten(i,k,j) = 0.
             rucuten(i,k,j)  = 0.
             rvcuten(i,k,j)  = 0.
          enddo
       enddo
    enddo

    do j = jts,jtf
       do k = kts,ktf
          do i = its,itf
             rthften(i,k,j)=0.
             rqvften(i,k,j)=0.
          enddo
       enddo
    enddo

    if(p_qc .ge. p_first_scalar) then
       do j = jts,jtf
          do k = kts,ktf
             do i = its,itf
                rqccuten(i,k,j)=0.
             enddo
          enddo
       enddo
    endif

    if(p_qi .ge. p_first_scalar) then
       do j = jts,jtf
          do k = kts,ktf
             do i = its,itf
                rqicuten(i,k,j)=0.
             enddo
          enddo
       enddo
    endif
 endif

 end subroutine ntiedtkeinit

!=================================================================================================================
 end module module_cu_ntiedtke
!=================================================================================================================
