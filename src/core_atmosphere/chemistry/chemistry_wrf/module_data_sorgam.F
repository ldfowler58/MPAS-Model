!=================================================================================================================
 module module_data_sorgam
 use mpas_kind_types,only: RKIND,R8KIND

 use module_data_radm2,only: epsilc

 implicit none

!
!param.inc start
 integer, parameter :: np = 8                   !bs maximum expected value of n
!integer numaer
!parameter (numaer=50)

 integer, parameter ::  maxits=100            !bs maximum number of iterations

 real(kind=RKIND), parameter :: tolf=1.e-09   !bs convergence criterion on function values
 real(kind=RKIND), parameter :: tolmin=1.e-12 !bs criterion whether superios convergence to
 real(kind=RKIND), parameter :: tolx=1.e-10   !bs convergence criterion on delta_x
 real(kind=RKIND), parameter :: stpmx=100.    !bs scaled maximum step length allowed
 real(kind=RKIND), parameter :: c303=19.83
 real(kind=RKIND), parameter :: c302=5417.4

 integer, parameter :: lcva=4
 integer, parameter :: lcvb=4
 integer, parameter :: lspcv=lcva+lcvb
 integer, parameter :: ldesn=13

!mh    ldesn is number of deposition species
!mh    true number of deposited species may be larger since there
!mh    are species which are deposited with the same rate

 integer, parameter :: laerdvc=39
 integer, parameter :: lnonaerdvc=8+lspcv
 integer, parameter :: l1ae=laerdvc+lnonaerdvc
 integer, parameter :: laero=4
 integer, parameter :: imodes=4
 integer, parameter :: aspec=1

!PARAMETER (laerdvc=34,lnonaerdvc=8+lspcv)
!LAERDVC  number of advected aerosol dynamic parameters for a given
!component species
!ia     L1AE        advected parameters+non-advected parameters
!ia     LAERO       number of aerosol component species
!ia     imodes      number of aerosol modes
!ia     ASPEC       number of gas phase comp. that are added dynamically
!ia		    currently only sulfate (=1)
!bs
!bs * BS ** BS ** BS ** BS ** BS ** BS ** BS ** BS ** BS ** BS ** BS **
!bs

 integer, parameter ::  aemiss=4
!bs *  AEMISS      # of aerosol species with emissions link to gas phase
!bs                currently ECI, ECJ, BCI, BCJ
 integer, parameter :: ldroga=11
 integer, parameter :: ldrogb=6
 integer, parameter :: ldrog=ldroga+ldrogb
!bs * LDROGA      # of anthropogenic organic aerosol precursor gases (DR
!bs * LDROGB      # of biogenic organic aerosol precursor gases (DROG)
!bs * LSPCV       # of condensable organic vapor interacting between gas
!bs               aerosol phase with SORGAM
!bs
!     param.inc stop

! //////////////////////////////////////////////////////////////////////
! FSB include file

! *** declare and set flag for organic aerosol production method
! *** Two method are available:

! *** The method of Pandis,Harley, Cass, and Seinfeld, 1992,
!     Secondary aerosol formation and transport, Atmos. Environ., 26A,
!     pp 2453-2466
!     Bowman et al. Atmospheric Environment
!     Vol 29, pp 579-589, 1995.
! *** and
! *** The method of Odum, Hoffmann, Bowman, Collins, Flagen and
!     Seinfeld, 1996, Gas/particle partitioning and secondary organic ae
!     yields, Environ. Sci, Technol, 30, pp 2580-2585.

 integer, parameter :: orgaer=2 ! 1 = Pandis et al.  1992 method is used
                                ! 2 = Pankow 1994/Odum et al. 1996 method is
! ***

! *** information about visibility variables
! number of visibility variables    
 integer, parameter :: n_ae_vis_spc=2

! index for visual range in deciview
 integer, parameter :: idcvw=1             

! index for extinction [ 1/km ]
 integer, parameter :: ibext=2

! *** set up indices for array  CBLK
! index for Accumulation mode sulfate aerosol
 integer, parameter :: vso4aj=1

! index for Aitken mode sulfate concentration
 integer, parameter :: vso4ai=2

! index for Accumulation mode aerosol ammonium
 integer, parameter :: vnh4aj=3

! index for Aitken mode ammonium concentration
 integer, parameter :: vnh4ai=4 

! index for Accumulation mode aerosol nitrate
 integer, parameter :: vno3aj=5

! index for Aitken mode nitrate concentration
 integer, parameter :: vno3ai=6

! index for Accumulation mode aerosol sodium
 integer, parameter :: vnaaj=7

! index for Aitken mode sodium concentration
 integer, parameter :: vnaai=8

! index for Accumulation mode aerosol chloride
 integer, parameter :: vclaj=9

! index for Aitken mode chloride concentration
 integer, parameter :: vclai=10

! index for Accumulation mode anthropogenic
 integer, parameter :: vorgaro1j=11

! index for Aitken mode anthropogenic org
 integer, parameter :: vorgaro1i=12 

! index for Accumulation mode anthropogenic
 integer, parameter :: vorgaro2j=13

! index for Aitken mode anthropogenic org
 integer, parameter :: vorgaro2i=14 

! index for Accumulation mode anthropogenic
 integer, parameter :: vorgalk1j=15

! index for Aitken mode anthropogenic org
 integer, parameter :: vorgalk1i=16 

! index for Accumulation mode anthropogenic
 integer, parameter :: vorgole1j=17 

! index for Aitken mode anthropogenic org
 integer, parameter :: vorgole1i=18 

! index for Accumulation mode biogenic aerosol
 integer, parameter :: vorgba1j=19 

! index for Aitken mode biogenic aerosol concentration
 integer, parameter :: vorgba1i=20 

! index for Accumulation mode biogenic aerosol
 integer, parameter :: vorgba2j=21 

! index for Aitken mode biogenic aerosol concentration
 integer, parameter :: vorgba2i=22  

! index for Accumulation mode biogenic aerosol
 integer, parameter :: vorgba3j=23 

! index for Aitken mode biogenic aerosol concentration
 integer, parameter :: vorgba3i=24 

! index for Accumulation mode biogenic aerosol
 integer, parameter :: vorgba4j=25

! index for Aitken mode biogenic aerosol concentration
 integer, parameter :: vorgba4i=26 

! index for Accumulation mode primary anthropogenic
 integer, parameter :: vorgpaj=27 

! index for Aitken mode primary anthropogenic
 integer, parameter :: vorgpai=28 

! index for Accumulation mode aerosol elemen
 integer, parameter :: vecj=29 

! index for Aitken mode elemental carbon
 integer, parameter :: veci=30     

! index for Accumulation mode primary PM2.5
 integer, parameter :: vp25aj=31 

! index for Aitken mode primary PM2.5 concentration
 integer, parameter :: vp25ai=32 

! index for coarse mode anthropogenic aerosol
 integer, parameter :: vantha=33 

! index for coarse mode marine aerosol concentration
 integer, parameter :: vseas=34 

! index for coarse mode soil-derived aerosol
 integer, parameter :: vsoila=35 

! index for Aitken mode number
 integer, parameter :: vnu0=36               

! index for accum  mode number
 integer, parameter :: vac0=37               

! index for coarse mode number
 integer, parameter :: vcorn=38 

! index for Accumulation mode aerosol water
 integer, parameter :: vh2oaj=39  

! index for Aitken mode aerosol water concentration
 integer, parameter :: vh2oai=40 

! index for Aitken mode 3'rd moment
 integer, parameter :: vnu3=41         

! index for Accumulation mode 3'rd moment
 integer, parameter :: vac3=42   

! index for coarse mode 3rd moment
 integer, parameter :: vcor3=43           

! index for sulfuric acid vapor concentration
 integer, parameter :: vsulf=44 

! index for nitric acid vapor concentration
 integer, parameter :: vhno3=45 

! index for ammonia gas concentration
 integer, parameter :: vnh3=46 

! index for HCL gas concentration
 integer, parameter :: vhcl=47 

! index for cond. vapor from aromatics
 integer, parameter :: vcvaro1=48    

! index for cond. vapor from aromatics
 integer, parameter :: vcvaro2=49    

! index for cond. vapor from anth. alkane
 integer, parameter :: vcvalk1=50 

! index for cond. vapor from anth. olefin
 integer, parameter :: vcvole1=51 

! index for cond. vapor from biogenics
 integer, parameter :: vcvapi1=52    

! index for cond. vapor from biogenics
 integer, parameter :: vcvapi2=53    

! index for cond. vapor from biogenics
 integer, parameter :: vcvlim1=54    

! index for cond. vapor from biogenics
 integer, parameter :: vcvlim2=55    


!           COMMON /CBLKINDCS/
!    &         VSO4AJ,VSO4AI,VNH4AJ,VNH4AI,VNO3AJ,VNO3AI,
!    &         VORGAJ,VORGAI, VORGPAJ,VORGPAI,
!    &         VORGBAJ,VORGBAI,VECJ,VECI,
!    &         VP25AJ,VP25AI,VANTHA,VSEAS,VSOILA,
!    &         VNU0,VAC0,VCORN,
!    &         VH2OAJ,VH2OAI,
!    &         VNU3,VAC3,VCOR3,
!    &         VSULF,VHNO3,VNH3

! *** set up species dimension and indices for sedimentation
!     velocity array VSED

! number of sedimentation velocities
 integer, parameter :: naspcssed=6          

! index for Aitken mode number
 integer, parameter :: vsnnuc=1                   

! index for Accumulation mode number
 integer, parameter :: vsnacc=2             

! index for coarse mode number
 integer, parameter :: vsncor=3                   

! index for Aitken mode mass
 integer, parameter :: vsmnuc=4                      

! index for accumulation mode mass
 integer, parameter :: vsmacc=5                

! index for coarse mass
 integer, parameter :: vsmcor=6                          

! *** set up species dimension and indices for deposition
!     velocity array VDEP

! number of deposition velocities
 integer, parameter :: naspcsdep=7             

! index for Aitken mode number
 integer, parameter :: vdnnuc=1                   

! index for accumulation mode number
 integer, parameter :: vdnacc=2 

! index for coarse mode number
 integer, parameter :: vdncor=3                   

! index for Aitken mode mass
 integer, parameter :: vdmnuc=4                     

! index for accumulation mode
 integer, parameter :: vdmacc=5                    

! index for fine mode mass (Aitken + accumulation)
 integer, parameter :: vdmfine=6 

! index for coarse mode mass
 integer, parameter :: vdmcor=7                     

! *** END AEROSTUFF.EXT
!bs
!BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * !
!BS *                                                                * !
!BS *            include file used in SORGAM routines                * !
!BS *                                                                * !
!BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * BS * !
!bs
!bs
!bs * species pointer for condensable vapor production
!bs
!bs XYL + OH                            
 integer, parameter :: pxyl=1
!bs TOL + OH                            
 integer, parameter :: ptol=2
!bs CSL + OH                            
 integer, parameter :: pcsl1=3
!bs CSL + NO                            
 integer, parameter :: pcsl2=4
!bs HC8 + OH                            
 integer, parameter :: phc8=5
!bs OLI + OH                            
 integer, parameter :: poli1=6
!bs OLI + NO                            
 integer, parameter :: poli2=7
!bs OLI + O3                            
 integer, parameter :: poli3=8
!bs OLT + OH                            
 integer, parameter :: polt1=9
!bs OLT + NO                            
 integer, parameter :: polt2=10
!bs OLT + O3                            
 integer, parameter :: polt3=11
!bs API + OH                            
 integer, parameter :: papi1=12
!bs API + NO                            
 integer, parameter :: papi2=13
!bs API + O3                            
 integer, parameter :: papi3=14
!bs LIM + OH                            
 integer, parameter :: plim1=15
!bs LIM + NO                            
 integer, parameter :: plim2=16
!bs LIM + O3                            
 integer, parameter :: plim3=17
!bs
!bs * Number of lumped condensable vapors in SORGAM
!bs
!bs      INTEGER NACV              !bs # of anth. cond. vapors
!bs      PARAMETER (NACV = 2)
!bs      INTEGER NBCV              !bs # of bio. cond. vapors
!bs      PARAMETER (NBCV = 1)
!bs      INTEGER NCV               !bs total # of cond. vapor
!bs      PARAMETER (NCV = NACV + NBCV)
!bs
!bs * species pointer for SOA species
!bs
 integer, parameter :: psoaaro1=1
 integer, parameter :: psoaaro2=2
 integer, parameter :: psoaalk1=3
 integer, parameter :: psoaole1=4
 integer, parameter :: psoaapi1=5
 integer, parameter :: psoaapi2=6
 integer, parameter :: psoalim1=7
 integer, parameter :: psoalim2=8
!bs
!bs * end of AERO_SOA.EXT *
!bs

! *** include file for aerosol routines


!....................................................................

!  CONTAINS: Fundamental constants for air quality modeling

!  DEPENDENT UPON:  none

!  REVISION HISTORY:

!    Adapted 6/92 by CJC from ROM's PI.EXT.

!    Revised 3/1/93 John McHenry to include constants needed by
!    LCM aqueous chemistry
!    Revised 9/93 by John McHenry to include additional constants
!    needed for FMEM clouds and aqueous chemistry

!    Revised 3/4/96 by Dr. Francis S. Binkowski to reflect current
!    Models3 view that MKS units should be used wherever possible,
!    and that sources be documentated. Some variables have been added
!    names changed, and values revised.

!    Revised 3/7/96 to have universal gas constant input and compute
!    gas constant is chemical form. TWOPI is now calculated rather than

!    Revised 3/13/96 to group declarations and parameter statements.

!    Revised 9/13/96 to include more physical constants.
!    Revised 12/24/96 eliminate silly EPSILON, AMISS

!    Revised 1/06/97 to eliminate most derived constants

! FSB REFERENCES:

!      CRC76,        CRC Handbook of Chemistry and Physics (76th Ed),
!                     CRC Press, 1995
!      Hobbs, P.V.   Basic Physical Chemistry for the Atmospheric Scien
!                     Cambridge Univ. Press, 206 pp, 1995.
!      Snyder, J.P., Map Projections-A Working Manual, U.S. Geological
!                     Paper 1395 U.S.GPO, Washington, DC, 1987.
!      Stull, R. B., An Introduction to Bounday Layer Meteorology, Klu
!                     Dordrecht, 1988

! Geometric Constants:

 real(kind=R8KIND), parameter :: pirs=3.14159265358979324
!real(kind=RKIND), parameter ::  PIRS = 3.141593
! Fundamental Constants: ( Source: CRC76, pp 1-1 to 1-6)

! Avogadro's Constant [ 1/mol ]
 real(kind=RKIND), parameter :: avo=6.0221367e23

! universal gas constant [ J/mol-K ]
 real(kind=RKIND), parameter :: rgasuniv=8.314510

! standard atmosphere  [ Pa ]
 real(kind=RKIND), parameter :: stdatmpa=101325.0

! Standard Temperature [ K ]
 real(kind=RKIND), parameter :: stdtemp=273.15

! Stefan-Boltzmann [ W/(m**2 K**4) ]
 real(kind=RKIND), parameter :: stfblz=5.67051e-8


! mean gravitational acceleration [ m/sec**2 ]
 real(kind=RKIND), parameter :: grav=9.80622
! FSB Non MKS qualtities:

! Molar volume at STP [ L/mol ] Non MKS units
 real(kind=RKIND), parameter :: molvol=22.41410


! Atmospheric Constants:

! FSB                     78.06%  N2, 21% O2 and 0.943% A on a mole
 real(kind=RKIND), parameter :: mwair=28.9628
                        ! fraction basis. ( Source : Hobbs, 1995) pp 69-
! mean molecular weight for dry air [ g/mol ]

! dry-air gas constant [ J / kg-K ]
 real(kind=RKIND), parameter :: rdgas=1.0E3*rgasuniv/mwair

!  3*PI
 real(kind=RKIND), parameter :: threepi=3.0*pirs

!  6/PI
 real(kind=RKIND), parameter :: f6dpi=6.0/pirs

!  1.0e9 * 6/PIRS
 real(kind=RKIND), parameter :: f6dpi9=1.0E9*f6dpi

! 1.0e-9 * 6/PIRS
 real(kind=RKIND), parameter :: f6dpim9=1.0E-9*f6dpi

!  SQRT( PI )
 real(kind=RKIND), parameter :: sqrtpi=1.7724539

!  SQRT( 2 )
 real(kind=RKIND), parameter :: sqrt2=1.4142135623731

!  ln( sqrt( 2 ) )
 real(kind=RKIND), parameter :: lgsqt2=0.34657359027997

!  1/ln( sqrt( 2 ) )
 real(kind=RKIND), parameter :: dlgsqt2=1.0/lgsqt2

!  1/3
 real(kind=RKIND), parameter :: one3=1.0/3.0

!  2/3
 real(kind=RKIND), parameter :: two3=2.0/3.0


! *** physical constants:

! Boltzmann's Constant [ J / K ]
 real(kind=RKIND), parameter :: boltz=rgasuniv/avo


! *** component densities [ kg/m**3 ] :


!  bulk density of aerosol sulfate
 real(kind=RKIND), parameter :: rhoso4=1.8e3

!  bulk density of aerosol ammonium
 real(kind=RKIND), parameter :: rhonh4=1.8e3

! bulk density of aerosol nitrate
 real(kind=RKIND), parameter :: rhono3=1.8e3

!  bulk density of aerosol water
 real(kind=RKIND), parameter :: rhoh2o=1.0e3

! bulk density for aerosol organics
 real(kind=RKIND), parameter :: rhoorg=1.0e3

! bulk density for aerosol soil dust
 real(kind=RKIND), parameter :: rhosoil=2.6e3

! bulk density for marine aerosol
 real(kind=RKIND), parameter :: rhoseas=2.2e3

! bulk density for anthropogenic aerosol
 real(kind=RKIND), parameter :: rhoanth=2.2e3

! bulk density of aerosol sodium
 real(kind=RKIND), parameter :: rhona=2.2e3

! bulk density of aerosol chloride
 real(kind=RKIND), parameter :: rhocl=2.2e3

! *** Factors for converting aerosol mass concentration [ ug m**-3] to
!      	  to 3rd moment concentration [ m**3 m^-3]

 real(kind=RKIND), parameter :: so4fac=f6dpim9/rhoso4

 real(kind=RKIND), parameter :: nh4fac=f6dpim9/rhonh4

 real(kind=RKIND), parameter :: h2ofac=f6dpim9/rhoh2o

 real(kind=RKIND), parameter :: no3fac=f6dpim9/rhono3

 real(kind=RKIND), parameter :: orgfac=f6dpim9/rhoorg

 real(kind=RKIND), parameter :: soilfac=f6dpim9/rhosoil

 real(kind=RKIND), parameter :: seasfac=f6dpim9/rhoseas

 real(kind=RKIND), parameter :: anthfac=f6dpim9/rhoanth

 real(kind=RKIND), parameter :: nafac=f6dpim9/rhona

 real(kind=RKIND), parameter :: clfac=f6dpim9/rhocl

!  starting standard surface pressure [ Pa ]  
 real(kind=RKIND), parameter :: pss0=101325.0

!  starting standard surface temperature [ K ]
 real(kind=RKIND), parameter :: tss0=288.15

!  initial sigma-G for nucleimode                 
 real(kind=RKIND), parameter :: sginin=1.70

!  initial sigma-G for accumulation mode          
 real(kind=RKIND), parameter :: sginia=2.00

! initial sigma-G for coarse mode               
 real(kind=RKIND), parameter :: sginic=2.5

!  initial mean diameter for nuclei mode [ m ]    
 real(kind=RKIND), parameter :: dginin=0.01e-6

!  initial mean diameter for accumulation mode [ m ]
 real(kind=RKIND), parameter :: dginia=0.07e-6

! initial mean diameter for coarse mode [ m ]  
 real(kind=RKIND), parameter :: dginic=1.0e-6



!................   end   AERO3box.EXT   ...............................
!///////////////////////////////////////////////////////////////////////





!     LOGICAL diagnostics
! *** Scalar variables for fixed standard deviations.

! Flag for writing diagnostics to file       
! nuclei mode exp( log^2( sigmag )/8 )  
 real(kind=RKIND) :: en1
! accumulation mode exp( log^2( sigmag )
 real(kind=RKIND) :: ea1
! coarse mode exp( log^2( sigmag )/8 )  
 real(kind=RKIND) :: ec1
! nuclei        **4                    
 real(kind=RKIND) :: esn04
! accumulation                         
 real(kind=RKIND) :: esa04
! coarse                               
 real(kind=RKIND) :: esc04
! nuclei        **5                    
 real(kind=RKIND) :: esn05

 real(kind=RKIND) :: esa05
! accumulation                         
! nuclei        **8                    
 real(kind=RKIND) :: esn08
! accumulation                         
 real(kind=RKIND) :: esa08

 real(kind=RKIND) :: esc08
! coarse                               
! nuclei        **9                    
 real(kind=RKIND) :: esn09

 real(kind=RKIND) :: esa09
! accumulation                         
! nuclei        **12                   
 real(kind=RKIND) :: esn12
! accumulation                         
 real(kind=RKIND) :: esa12

 real(kind=RKIND) :: esc12
! coarse mode                          
! nuclei        **16                   
 real(kind=RKIND) :: esn16
! accumulation                         
 real(kind=RKIND) :: esa16

 real(kind=RKIND) :: esc16
! coarse                               
! nuclei        **20                   
 real(kind=RKIND) :: esn20
! accumulation                         
 real(kind=RKIND) :: esa20

 real(kind=RKIND) :: esc20
! coarse                               
! nuclei        **25                   
 real(kind=RKIND) :: esn25

 real(kind=RKIND) :: esa25
! accumulation                         
! nuclei        **24                   
 real(kind=RKIND) :: esn24
! accumulation                         
 real(kind=RKIND) :: esa24

 real(kind=RKIND) :: esc24
! coarse                               
! nuclei        **28                   
 real(kind=RKIND) :: esn28
! accumulation                         
 real(kind=RKIND) :: esa28

 real(kind=RKIND) :: esc28
! coarse                               
! nuclei        **32                   
 real(kind=RKIND) :: esn32
! accumulation                         
 real(kind=RKIND) :: esa32

 real(kind=RKIND) :: esc32
! coarese                              
! nuclei        **36                   
 real(kind=RKIND) :: esn36
! accumulation                         
 real(kind=RKIND) :: esa36

 real(kind=RKIND) :: esc36
! coarse                               
! nuclei        **49                   
 real(kind=RKIND) :: esn49

 real(kind=RKIND) :: esa49
! accumulation                         
! nuclei        **52                   
 real(kind=RKIND) :: esn52

 real(kind=RKIND) :: esa52
! accumulation                         
! nuclei        **64                   
 real(kind=RKIND) :: esn64
! accumulation                         
 real(kind=RKIND) :: esa64

 real(kind=RKIND) :: esc64
! coarse                               

 real(kind=RKIND) :: esn100
! nuclei        **100                  
! nuclei        **(-20)                
 real(kind=RKIND) :: esnm20
! accumulation                         
 real(kind=RKIND) :: esam20

 real(kind=RKIND) :: escm20
! coarse                               
! nuclei        **(-32)                
 real(kind=RKIND) :: esnm32
! accumulation                         
 real(kind=RKIND) :: esam32

 real(kind=RKIND) :: escm32
! coarse                               
! log(sginin)                           
 real(kind=RKIND) :: xxlsgn
! log(sginia)                           
 real(kind=RKIND) :: xxlsga

 real(kind=RKIND) :: xxlsgc
! log(sginic )                          
! log(sginin ) ** 2                           
 real(kind=RKIND) :: l2sginin
! log(sginia ) ** 2                           
 real(kind=RKIND) :: l2sginia

 real(kind=RKIND) :: l2sginic


! *** set up COMMON blocks for esg's:



! log(sginic ) ** 2                           

! *** SET NUCLEATION FLAG:

 integer, parameter :: inucl=2           ! INUCL = 0, Kerminen & Wexler Mechanism
                                         ! INUCL = 1, Youngblood and Kreidenweis mech
                                         ! INUCL = 2, Kulmala et al. mechanism

! *** Set flag for sedimentation velocities:
 logical, parameter ::  icoarse=.false.  ! *** END AERO_INTERNAL.EXT

! *** Diameters and standard deviations for emissions
!     the diameters are the volume (mass) geometric mean diameters

! *** Aitken mode:
! special factor to compute mass transfer           
 real(kind=RKIND), parameter :: dgvem_i=0.03e-6 ! [ m ]
 real(kind=RKIND), parameter :: sgem_i=1.7

! *** Accumulation mode:
 real(kind=RKIND), parameter :: dgvem_j=0.3e-6 ! [ m ]
 real(kind=RKIND), parameter :: sgem_j=2.0

! *** Coarse mode
 real(kind=RKIND), parameter :: dgvem_c=6.0e-6 ! [ m ] <<< Corrected 11/19/97
 real(kind=RKIND), parameter :: sgem_c=2.2

! *** factors for getting number emissions rate from mass emissions rate
! Aitken mode                                       
 real(kind=RKIND) :: factnumn
! accumulation mode
 real(kind=RKIND) :: factnuma
! coarse mode
 real(kind=RKIND) :: factnumc

 real(kind=RKIND), parameter :: facatkn_min=0.04
 real(kind=RKIND), parameter :: facacc_min=1.0-facatkn_min
 real(kind=RKIND), parameter :: conmin=epsilc
 real(kind=RKIND):: xxm3

! [ ug/m**3 ] ! changed 1/6/98 
 real(kind=R8KIND) :: nummin_i ! factor to set minimum for Aitken mode number  
 real(kind=R8KIND) :: nummin_j ! factor to set minimum for accumulation mode nu
 real(kind=R8KIND) :: nummin_c ! factor to set minimum for coarse mode number  

!bs
!bs      real(kind=RKIND):: ALPHSULF ! Accommodation coefficient for sulfuric acid
!bs      PARAMETER ( ALPHSULF = 0.05 ) ! my be set to one in future
!bs
!bs      real(kind=RKIND):: DIFFSULF ! molecular diffusivity for sulfuric acid [ m**2
!bs      PARAMETER( DIFFSULF = 0.08E-4 ) ! may be changed in future
!bs
!bs * 23/03/99 updates of ALPHSULF and DIFFSULF adopted fro new code fro
!bs * DIFFSULF is calculated from Reid, Prausnitz, and Poling, The prope
!bs * of gases and liquids, 4th edition, McGraw-Hill, 1987, pp 587-588.
!bs * Equation (11-4.4) was used.
!bs * The value is at T = 273.16 K and P = 1.01325E05 Pa
!bs * Temperature dependence is included for DIFFSULF via DIFFCORR (see
!bs
! Accommodation coefficient for sulfuric
 real(kind=RKIND), parameter :: alphsulf=1.0

!bs updated from code of FSB         
! molecular weight for sulfuric acid [ kg/mole ] MKS 
 real(kind=RKIND), parameter :: mwh2so4=98.07354e-3

!cia corrected error 24/11/97
! molecular diffusivity for sulfuric acid [ m**2 /se
 real(kind=RKIND), parameter :: diffsulf=9.362223e-06

!bs updated from code of FSB 
!bs Accomodation coefficient for organic
 real(kind=RKIND), parameter :: alphaorg=1.0 !bs Kleeman et al. '99 propose alpha

!bs Bowman et al. '97 uses alpha = 1.
!bs mean molecular weight of organics [k
 real(kind=RKIND), parameter :: mworg=175.0e-03

!bs
!bs * DIFFORG is calculated from the same formula as DIFFSULF.
!bs * An average elemental composition of C=8, O=3, N=1, H=17 is asuumed
!bs * to calculate DIFFORG at T = 273.16K and  P = 1.01325E05 Pa.
!bs * Temepratur dependence is included below.
!bs molecular diffusivity for organics [
 real(kind=RKIND), parameter :: difforg=5.151174e-06

! *** CCONC is the factor for near-continuum condensation.
! ccofm * sqrt( ta )                    
 real(kind=RKIND), parameter :: cconc=2.0*pirs*diffsulf

!bs * factor for NC condensation for organics
! [ m**2 / sec ]       
 real(kind=RKIND), parameter :: cconc_org=2.0*pirs*difforg
! [ m**2 / sec ]    

!bs analogue to CCOFM but for organics  
 real(kind=RKIND) :: ccofm_org

! FSB  CCOFM is  the accommodation coefficient
!      times the mean molecular velocity for h2so4 without the temperatu
!      after some algebra
!bs CCOFM_ORG * sqrt(TA)                
! set to a value below                  
 real(kind=RKIND) :: ccofm

! minimum aerosol sulfate concentration          
 real(kind=RKIND), parameter :: aeroconcmin=0.0001


!*******************************************************************
!*                                                                 *
!*  start parameters and variables for aerosol-cloud interactions  *
!*                                                                 *
!*******************************************************************
!
!   maxd_atype = maximum allowable number of aerosol types
!   maxd_asize = maximum allowable number of aerosol size bins
!   maxd_acomp = maximum allowable number of chemical components
!	in each aerosol size bin
!   maxd_aphase = maximum allowable number of aerosol phases (gas, cloud, ice, rain, ...)
!
!   ntype_aer = number of aerosol types
!   nsize_aer(t) = number of aerosol size bins for aerosol type t. each bin w/ same set of components
!   nphase_aer = number of aerosol phases
!
!   msectional - if positive, moving-center sectional code is utilized,
!	and each mode is actually a section.
!   maerosolincw - if positive, both unactivated/interstitial and activated
!	aerosol species are simulated.  if zero/negative, only the
!	unactivated are simulated.
!
!   ncomp_aer(t) = number of chemical components for aerosol type t
!   ncomp_aer_nontracer(t) = number of "non-tracer" chemical components while in gchm code
!   mastercompptr_aer(c,t) = mastercomp type/i.d. for chemical component c
!	(1=sulfate, others to be defined) and aerosol type t.
!   massptr_aer(c,s,t,p) = gchm r-array index for the mixing ratio
!	(moles-x/mole-air) for chemical component c in size bin s for type t and phase p
!
!   waterptr_aer(s,t) = mixing ratio (moles-water/mole-air) for water
!       associated with aerosol size bin s and type t
!   hygroptr_aer(s,t) = gchm r-array index for the bulk hygroscopicity of the size bin and type
!   numptr_aer(s,t,p) = gchm r-array index for the number mixing ratio
!	(particles/mole-air) for aerosol size bin s, type t, and phase p
!       If zero or negative, then number is not being simulated.
!
!   mprognum_aer(s,t,p) - if positive, number mixing-ratio for size s, type t,
!       and phase p will be prognosed.  Otherwise, no.
!
!   ntot_mastercomp_aer = number of aerosol chemical components defined
!   dens_mastercomp_aer(mc) = dry density (g/cm^3) of aerosol master chemical component type c
!   mw_mastercomp_aer(mc) = molecular weight of aerosol master chemical component type mc
!   name_mastercomp_aer(mc) = name of aerosol master chemical component type mc
!   mc=mastercompptr_aer(c,t)
!   dens_aer(c,t) = dry density (g/cm^3) of aerosol chemical component type c and type t
!   mw_aer(c,t) = molecular weight of aerosol chemical component type c and type t
!   name_aer(c,t) = name of aerosol chemical component type c and type t
!
!   lptr_so4_aer(s,t,p) = gchm r-array index for the
!	mixing ratio for sulfate associated with aerosol size bin s, type t, and phase p
!   (similar for msa, oc, bc, nacl, dust)
!
!-----------------------------------------------------------------------
!
!   volumcen_sect(s,t)= volume (cm^3) at center of section m
!   volumlo_sect(s,t) = volume (cm^3) at lower boundary of section m
!   volumhi_sect(s,t) = volume (cm^3) at upper boundary of section m
!
!   dlo_sect(s,t) = diameter (cm) at lower boundary of section m
!   dhi_sect(s,t) = diameter (cm) at upper boundary of section m
!   dcen_sect(s,t) = volume arithmetic-mean diameter (cm) of section m
!	(corresponds to volumcen_sect == 0.5*(volumlo_sect + volumhi_sect)
!
!-----------------------------------------------------------------------
!   nov-04 sg ! replaced amode with aer and expanded aerosol dimension to include type and phase

 integer, parameter :: maxd_atype = 2
 integer, parameter :: maxd_asize = 2
 integer, parameter :: maxd_acomp = 19
 integer, parameter :: maxd_aphase = 2
 integer, save :: ai_phase ! interstitial phase of aerosol
 integer, save :: cw_phase ! cloud water phase of aerosol
 integer, save :: ci_phase ! cloud ice  phase of aerosol
 integer, save :: cr_phase ! rain  phase of aerosol
 integer, save :: cs_phase ! snow  phase of aerosol
 integer, save :: cg_phase ! graupel phase of aerosol

 integer, save :: ntype_aer = 0           ! number of types
 integer, save :: ntot_mastercomp_aer = 0 ! number of master components
 integer, save :: nphase_aer = 0          ! number of phases

 integer, save :: msectional , maerosolincw

 integer, dimension(maxd_atype), save ::   &
    nsize_aer,                          & ! number of size bins
    ncomp_aer,                          & ! number of chemical components
    ncomp_aer_nontracer                   !

 integer, dimension(maxd_acomp,maxd_atype), save ::   &
    mastercompptr_aer,                  & !  mastercomp index
    waterptr_aer,                       & ! index for aerosol water
    hygroptr_aer                          ! index for aerosol hygroscopicity

 integer, dimension(maxd_asize,maxd_atype,maxd_aphase), save :: &
    numptr_aer,                         & ! index for the number mixing ratio
    mprognum_aer

 integer, dimension(maxd_acomp,maxd_asize,maxd_atype,maxd_aphase), save :: &
    massptr_aer                           ! index for mixing ratio

 real(kind=RKIND), dimension(maxd_acomp), save ::   &
         dens_mastercomp_aer, mw_mastercomp_aer, hygro_mastercomp_aer

 real(kind=RKIND), dimension(maxd_acomp,maxd_atype), save :: dens_aer , mw_aer, hygro_aer

 character(len=10),dimension(maxd_acomp), save :: name_mastercomp_aer
 character(len=10),dimension(maxd_acomp,maxd_atype), save :: name_aer

 real(kind=RKIND), dimension(maxd_asize,maxd_atype), save :: &
    volumcen_sect , volumlo_sect , volumhi_sect ,   &
    dcen_sect     ,  dlo_sect    , dhi_sect     ,   &
    sigmag_aer

 integer, dimension(maxd_asize,maxd_atype,maxd_aphase), save :: &
    lptr_so4_aer     , lptr_nh4_aer     , lptr_no3_aer    , &
    lptr_orgaro1_aer , lptr_orgaro2_aer , lptr_orgalk_aer , &
    lptr_orgole_aer  , lptr_orgba1_aer  , lptr_orgba2_aer , &
    lptr_orgba3_aer  , lptr_orgba4_aer  , lptr_orgpa_aer  , &
    lptr_ec_aer      , lptr_p25_aer     , lptr_anth_aer   , &
    lptr_cl_aer      , lptr_na_aer      , lptr_seas_aer   , &
    lptr_soil_aer 

 logical, dimension(maxd_asize,maxd_atype), save :: do_cloudchem_aer


!   molecular weights (g/mol)
 real(kind=RKIND), parameter :: mw_so4_aer   = 96.066
 real(kind=RKIND), parameter :: mw_no3_aer   = 62.007
 real(kind=RKIND), parameter :: mw_nh4_aer   = 18.042
 real(kind=RKIND), parameter :: mw_oc_aer    = 250.0
 real(kind=RKIND), parameter :: mw_ec_aer    = 1.0
 real(kind=RKIND), parameter :: mw_oin_aer   = 1.0
 real(kind=RKIND), parameter :: mw_dust_aer  = 100.087
 real(kind=RKIND), parameter :: mw_seas_aer  = 58.440
 real(kind=RKIND), parameter :: mw_cl_aer    = 35.450
 real(kind=RKIND), parameter :: mw_na_aer    = 22.990
 real(kind=RKIND), parameter :: mw_water_aer = 18.016

!   dry densities (g/cm3)
 real(kind=RKIND), parameter :: dens_so4_aer  = 1.80   ! = rhoso4
 real(kind=RKIND), parameter :: dens_no3_aer  = 1.80   ! = rhono3
 real(kind=RKIND), parameter :: dens_nh4_aer  = 1.80   ! = rhonh4
 real(kind=RKIND), parameter :: dens_oc_aer   = 1.00   ! = rhoorg
 real(kind=RKIND), parameter :: dens_ec_aer   = 1.70
 real(kind=RKIND), parameter :: dens_dust_aer = 2.60   ! = rhosoil
 real(kind=RKIND), parameter :: dens_oin_aer  = 2.20   ! = rhoanth
 real(kind=RKIND), parameter :: dens_seas_aer = 2.20   ! = rhoseas
 real(kind=RKIND), parameter :: dens_cl_aer   = 2.20
 real(kind=RKIND), parameter :: dens_na_aer   = 2.20

!   water density (g/cm3)
 real(kind=RKIND), parameter :: dens_water_aer  = 1.0

!   hygroscopicity (dimensionless)
 real(kind=RKIND), parameter :: hygro_so4_aer  = 0.5
 real(kind=RKIND), parameter :: hygro_no3_aer  = 0.5
 real(kind=RKIND), parameter :: hygro_nh4_aer  = 0.5
 real(kind=RKIND), parameter :: hygro_oc_aer   = 0.14
 real(kind=RKIND), parameter :: hygro_ec_aer   = 1.e-6
 real(kind=RKIND), parameter :: hygro_oin_aer  = 0.14
 real(kind=RKIND), parameter :: hygro_dust_aer = 0.1
 real(kind=RKIND), parameter :: hygro_seas_aer = 1.16
 real(kind=RKIND), parameter :: hygro_cl_aer   = 1.16
 real(kind=RKIND), parameter :: hygro_na_aer   = 1.16

! table lookup of aerosol impaction/interception scavenging rates
 integer, parameter :: nimptblgrow_mind = -14
 integer, parameter :: nimptblgrow_maxd = 24

 real(kind=RKIND) :: dlndg_nimptblgrow

 real(kind=RKIND), dimension(4,nimptblgrow_mind:nimptblgrow_maxd,maxd_asize,maxd_atype) :: scavimptblnum
 real(kind=RKIND), dimension(4,nimptblgrow_mind:nimptblgrow_maxd,maxd_asize,maxd_atype) :: scavimptblvol


!SAM 10/08 Gaussian quadrature constants for SORGAM deposition numerical integration
 integer, parameter :: NGAUSdv = 7 ! Number of Gaussian Quadrature Points
                                   ! (constants defined in aerosols_sorgam_init).
 real(kind=RKIND), dimension(NGAUSdv) :: y_gq, wgaus

!*****************************************************************
!*                                                               *
!*  end parameters and variables for aerosol-cloud interactions  *
!*                                                               *
!*****************************************************************

!=================================================================================================================
 end module module_data_sorgam
!=================================================================================================================
