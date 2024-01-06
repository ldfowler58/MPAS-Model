!=================================================================================================================
 module module_data_gocart_seas
!=================================================================================================================
 use mpas_kind_types,only: RKIND,R8KIND

 implicit none
 save

 real(kind=RKIND),dimension(4,12):: ch_ss

 real(kind=RKIND),parameter:: pi=3.141592653559
 real(kind=RKIND),dimension(4),parameter:: reff_seas(4)=(/0.30e-6,1.00e-6,3.25e-6,7.50e-6/)
 real(kind=RKIND),dimension(4),parameter:: den_seas(4)=(/2.2e3,2.2e3,2.2e3,2.2e3/)

 real(kind=RKIND),dimension(4),parameter:: ra_seas(4)=(/1.e-1,5.e-1,1.5e0,5.0e0/)
 real(kind=RKIND),dimension(4),parameter:: rb_seas(4)=(/5.e-1,1.5e0,5.e0,1.e1/)

!=================================================================================================================
 end module module_data_gocart_seas
!=================================================================================================================
