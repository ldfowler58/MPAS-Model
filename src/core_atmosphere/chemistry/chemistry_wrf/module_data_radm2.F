!=================================================================================================================
 module module_data_radm2
 use mpas_kind_types,only: RKIND

 implicit none

 real(kind=RKIND),parameter:: epsilc = 1.E-16

!--- for radm solver
! .. Parameters ..
 integer,parameter:: ldiag = 18, lpred = 39, lss = 2, &
        lump = 4, naqre = 70,  nreacj = 21, nreack = 140, &
        ntroe = 7, numchem_radm = 41
 integer,parameter:: lspec = lpred + lss
 integer,dimension(1:ntroe):: itroe = (/11, 22, 10, 15, 21, 24, 28/)
!
!
!
 integer,parameter:: lso2=1
 integer,parameter:: lsulf=2
 integer,parameter:: lno2=3
 integer,parameter:: lno=4
 integer,parameter:: lo3=5
 integer,parameter:: lhno3=6
 integer,parameter:: lh2o2=7
 integer,parameter:: lald=8
 integer,parameter:: lhcho=9
 integer,parameter:: lop1=10
 integer,parameter:: lop2=11
 integer,parameter:: lpaa=12
 integer,parameter:: lora1=13
      
 integer,parameter:: lora2=14
 integer,parameter:: lnh3=15
 integer,parameter:: ln2o5=16
 integer,parameter:: lno3=17
 integer,parameter:: lpan=18
 integer,parameter:: lhc3=19
 integer,parameter:: lhc5=20
 integer,parameter:: lhc8=21
      
 integer,parameter:: leth=22
 integer,parameter:: lco=23
 integer,parameter:: lol2=24
 integer,parameter:: lolt=25
 integer,parameter:: loli=26
 integer,parameter:: ltol=27
 integer,parameter:: lxyl=28
 integer,parameter:: laco3=29
      
 integer,parameter:: ltpan=30
 integer,parameter:: lhono=31
 integer,parameter:: lhno4=32
 integer,parameter:: lket=33
 integer,parameter:: lgly=34
 integer,parameter:: lmgly=35
 integer,parameter:: ldcb=36
 integer,parameter:: lonit=37
      
 integer,parameter:: lcsl=38
 integer,parameter:: liso=39
 integer,parameter:: lho=40
 integer,parameter:: lho2=41
! parameters for timestep, integration
 integer,dimension(1:lpred):: intgrt = (/1, 1, 1, 0, 1, &
                          1, 0, 1, 1, 1, 1, 1, 1, 1, 1, &
                          0, 1, 0, 1, 1, 1, 1, 1, 1, 1, &
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                            1, 1, 1, 1 /)
!integer,dimension(1:lspec)::  qdtc  = (/0, 0, 1, 0, 1, &
!                         0, 1, 0, 1, 0, 0, 0, 1, 0, 0, &
!                         1, 1, 1, 0, 0, 0, 0, 0, 0, 0, &
!                         0, 0, 0, 1, 1, 1, 1, 0, 0, 0, &
!                                     0, 0, 0, 0, 0, 0 /)
 integer,DIMENSION(1:lspec)::  qdtc  = (/1, 1, 1, 0, 1, &
                          1, 0, 1, 1, 1, 1, 1, 1, 1, 1, &
                          0, 1, 0, 1, 1, 1, 1, 1, 1, 1, &
                          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                      1, 1, 1, 1, 0, 0 /)
! max, min values,
 integer:: itrdu
!
 real(kind=RKIND),dimension(1:lspec):: cmin =(/(1.E-16,itrdu=1,lspec)/)
!
 real(kind=RKIND),dimension(1:lspec):: cmax=(/1.,  1.,  1.,   1., .2, &
                     3., .05, .01, .01, .01, .05, .01, .05,  .05,.05, &
                     1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,   1., 1., &
                     1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,.0001, .1, &
                                         1., .001, .01, .01, .01, .01/) 

!
!
!
 integer,parameter:: lo3p=1
 integer,parameter:: lo1d=2
 integer,parameter:: ltco3=3
 integer,parameter:: lhc3p=4
 integer,parameter:: lhc5p=5
 integer,parameter:: lhc8p=6
       
 integer,parameter:: lol2p=7
 integer,parameter:: loltp=8
 integer,parameter:: lolip=9
 integer,parameter:: ltolp=10
 integer,parameter:: lxylp=11
 integer,parameter:: lethp=12
 integer,parameter:: lketp=13
 integer,parameter:: loln=14
       
 integer,parameter:: lxo2=15
 integer,parameter:: lxno2=16
 integer,parameter:: lxho=17
 integer,parameter:: lmo2=18
!
!
 integer,parameter:: lnox=1
 integer,parameter:: lhox=2
 integer,parameter:: lpao3=3
 integer,parameter:: ln2n3=4
! ..
 real(kind=RKIND),parameter:: ch4=1.7
 real(kind=RKIND),parameter:: co2=350.
 real(kind=RKIND),parameter:: n2=7.81E5
 real(kind=RKIND),parameter:: o2=2.09E5
 real(kind=RKIND),parameter:: pi=3.141592654

! ..
 real(kind=RKIND),dimension(2):: afac,bfac
 real(kind=RKIND),dimension(3):: const
 real(kind=RKIND),dimension(nreack):: eor,thafac
 real(kind=RKIND),dimension(ntroe) :: xk0300,xkf300,xmtroe,xntroe

! ..
! .. Data Statements ..
 data thafac/0.00, 6.50E-12, 1.80E-11, 3.20E-11, 2.20E-10, 2.00E-12,       &
     1.60E-12, 1.10E-14, 3.70E-12, 4*0.00, 3.30E-12, 0.00, 3.30E-19,       &
     1.40E-13, 1.70E-11, 2.50E-14, 2.50E-12, 2*0.00, 2.00E-21, 2*0.00,     &
     1.30E-12, 4.60E-11, 2*0.00, 6.95E-18, 1.37E-17, 1.59E-11, 1.73E-11,   &
     3.64E-11, 2.15E-12, 5.32E-12, 1.07E-11, 2.10E-12, 1.89E-11, 4.00E-11, &
     9.00E-12, 6.87E-12, 1.20E-11, 1.15E-11, 1.70E-11, 2.80E-11, 1.00E-11, &
     1.00E-11, 1.00E-11, 6.85E-18, 1.55E-11, 2.55E-11, 2.80E-12, 1.95E+16, &
     4.70E-12, 1.95E+16, 4.20E-12, 4.20E-12, 0.00, 4.20E-12, 0.00,         &
     4.20E-12, 0.00, 10*4.20E-12, 6.00E-13, 1.40E-12, 6.00E-13, 1.40E-12,  &
     1.40E-12, 2.20E-11, 2.00E-12, 1.00E-11, 3.23E-11, 5.81E-13, 1.20E-14, &
     1.32E-14, 7.29E-15, 1.23E-14, 14*7.70E-14, 1.90E-13, 1.40E-13,        &
     4.20E-14, 3.40E-14, 2.90E-14, 1.40E-13, 1.40E-13, 1.70E-14, 1.70E-14, &
     9.60E-13, 1.70E-14, 1.70E-14, 9.60E-13, 3.40E-13, 1.00E-13, 8.40E-14, &
     7.20E-14, 3.40E-13, 3.40E-13, 4.20E-14, 4.20E-14, 1.19E-12, 4.20E-14, &
     4.20E-14, 1.19E-12, 7.70E-14, 1.70E-14, 4.20E-14, 3.60E-16, 4.20E-12, &
     4.20E-12, 7.70E-14, 1.70E-14, 4.20E-14, 3.60E-16, 0.00, 1.70E-14,     &
     4.20E-14, 3.60E-16/
! ..
! constants for RADM2 rate coefficients
 data eor/0., -120., -110., -70., 0., 1400., 940., 500., -240., 0., 0.,          &
          0., 0., 200., 0., -530., 2500., -150., 1230., 0., 0., 0., 0., 0., 0.,  &
       -380., -230., 0., 0., 1280., 444., 540., 380., 380., -411., -504.,        &
       -549., -322., -116., 0., 0., -256., 745., 0., 0., 0., 0., 0., 0.,         &
        444., 540., -409., -181., 13543., 0., 13543., -180., -180., 0., -180.,   &
          0., -180., 0., -180., -180., -180., -180., -180., -180., -180., -180., &
       -180., -180., 2058., 1900., 2058., 1900., 1900., 0., 2923., 1895.,        &
        975., 0., 2633., 2105., 1136., 2013., -1300., -1300., -1300., -1300.,    &
      -1300., -1300., -1300., -1300., -1300., -1300., -1300., -1300.,            &
      -1300., -1300., 25* -220., -1300., -220., -220., -220., -180., -180.,      &
      -1300., -220., -220., 0., 0., -220., -220., -220./

 data xk0300/1.8E-31, 2.2E-30, 1.8E-31, 7.E-31, 2.2E-30, 2.6E-30, 3.E-31/
 data xntroe/3.2, 4.3, 3.2, 2.6, 4.3, 3.2, 3.3/
 data xkf300/4.7E-12, 1.5E-12, 4.7E-12, 1.5E-11, 1.5E-12, 2.4E-11,1.5E-12/
 data xmtroe/1.4, 0.5, 1.4, 2*.5, 1.3, 0./
 data afac/2.1E-27, 1.1E-27/
 data bfac/10900., 11200./
 data const/7.34E21, 4.4E17, 3.23E33/

!=================================================================================================================
 end module module_data_radm2
!=================================================================================================================
