!=================================================================================================================
 module module_data_gocart_chem
 use mpas_kind_types,only: RKIND

 implicit none

!-----------------------------------------------------------------------  
! tracer info
!-----------------------------------------------------------------------
! Tracer index:
! default initialization for all sulfur and carbon species is 0 (undefined)
!     1. DMS       = Dimethyl sulfide            = CH3SCH3  
!     2. SO2       = Sulfur dioxide              = SO2               
!     3. SO4       = Sulfate                     = SO4            
!     4. MSA       = Methane sulfonic acid       = CH3SO3H             
 integer:: ndms=1, nso2=2, nso4=3, nmsa=4 
 real(kind=RKIND),parameter:: airmw    = 28.97
 real(kind=RKIND),parameter:: smw      = 32.00 
 real(kind=RKIND),parameter:: nh4_mfac = 1.375 ! increase sulf (output ond AOD only)
                                               ! to account for missing nh4
 real(kind=RKIND),parameter:: oc_mfac = 1.8    ! increase oc (output ond AOD nly)
                                               ! to account for Carbon to Organic mass

!=================================================================================================================
 end module module_data_gocart_chem
!=================================================================================================================
