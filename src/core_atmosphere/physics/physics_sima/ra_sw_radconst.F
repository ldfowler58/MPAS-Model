!=================================================================================================================
 module ra_sw_radconst
 use mpas_kind_types,only: kind_phys => RKIND

 implicit none
 private
 public:: radconst


 contains


!=================================================================================================================
 subroutine radconst(declin,solcon,julian,degrad,dpd,solcon_0)
!=================================================================================================================

!input arguments:
!integer,intent(in):: julian
 real(kind=kind_phys),intent(in):: julian
 real(kind=kind_phys),intent(in):: degrad,dpd
 real(kind=kind_phys),intent(in):: solcon_0

!output arguments:
 real(kind=kind_phys),intent(out):: declin,solcon

!local variables:
 real(kind=kind_phys):: obecl,sinob,sxlong,arg,decdeg,djul,rjul,eccfac

!-----------------------------------------------------------------------------------------------------------------

 declin=0.
 solcon=0.

!obecl : obliquity = 23.5 degree.
        
 obecl=23.5*degrad
 sinob=sin(obecl)
        
!calculate longitude of the sun from vernal equinox:        

 if(julian.ge.80.)sxlong=dpd*(julian-80.)
 if(julian.lt.80.)sxlong=dpd*(julian+285.)
 sxlong=sxlong*degrad
 arg=sinob*sin(sxlong)
 declin=asin(arg)
 decdeg=declin/degrad

!solar constant eccentricity factor (paltridge and platt 1976)

 djul=julian*360./365.
 rjul=djul*degrad
 eccfac=1.000110+0.034221*cos(rjul)+0.001280*sin(rjul)+0.000719*  &
        cos(2*rjul)+0.000077*sin(2*rjul)
 solcon=solcon_0*eccfac

 end subroutine radconst

!=================================================================================================================
 end module ra_sw_radconst
!=================================================================================================================
