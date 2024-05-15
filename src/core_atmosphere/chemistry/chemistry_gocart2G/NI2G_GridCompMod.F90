!=================================================================================================================
 module NI2G_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process

 use NI2G_instance_NI,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                            fnum,rhFlag,pressure_lid_in_hPa
 use NI2G_StateSpecs,only: NI2G_State


 implicit none
!private
!public:: load_NI2G_GridComp


 type:: ThreadWorkspace
    logical:: first = .true.
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: NI2G_GridComp
!   logical:: first
!   logical:: recycle_HNO3 = .false.
!   real(kind=RKIND),dimension(:),allocatable:: rmedDU,rmedSS ! DU and SS radius
!   real(kind=RKIND),dimension(:),allocatable:: fnumDU,fnumSS ! DU and SS particles per kg mass
!   type(ThreadWorkspace),dimension(:),allocatable:: workspaces

    contains
       procedure:: load_GridComp => load_NI2G_GridComp
 end type NI2G_GridComp

 type wrap_
    type(NI2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_NI2G_GridComp(self)
!=================================================================================================================

!--- inout arguments:
 class(NI2G_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine load_NI2G_GridComp:')

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo

 call mpas_log_write('--- end subroutine load_NI2G_GridCOMP:')
!call mpas_log_write(' ')

 end subroutine load_NI2G_GridComp

!=================================================================================================================
 endmodule NI2G_GridCompMod
!=================================================================================================================
