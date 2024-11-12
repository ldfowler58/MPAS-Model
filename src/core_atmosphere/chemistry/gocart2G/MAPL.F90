!=================================================================================================================
 module MAPL

 implicit none
 private
 public:: MAPL_UnpackTime,MAPL_PackTime


 contains


!=================================================================================================================
 subroutine MAPL_UnpackTime(TIME,IYY,IMM,IDD)
    integer, intent (IN ) :: TIME
    integer, intent (OUT) :: IYY
    integer, intent (OUT) :: IMM
    integer, intent (OUT) :: IDD
    IYY = TIME/10000 
    IMM = mod(TIME/100,100)
    IDD = mod(TIME,100)
 end subroutine MAPL_UnpackTime

 subroutine MAPL_PackTime(TIME,IYY,IMM,IDD)
    integer, intent (OUT) :: TIME
    integer, intent (IN ) :: IYY
    integer, intent (IN ) :: IMM
    integer, intent (IN ) :: IDD
    TIME=IYY*10000+IMM*100+IDD              
 end subroutine MAPL_PackTime

!=================================================================================================================
 end module MAPL
!=================================================================================================================
