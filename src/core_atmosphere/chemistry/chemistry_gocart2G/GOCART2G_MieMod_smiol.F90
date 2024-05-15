!=================================================================================================================
 module gocart2G_MieMod_smiol
 use mpas_kind_types
 use mpas_derived_types
 use mpas_log,only: mpas_log_write

 implicit none
 private
 public:: GOCART2G_Mie

 integer,parameter:: NRH_BINS = 991

 type GOCART2G_Mie
    character(len=:),allocatable :: table_name
    integer:: nch  ! number of channels in table (replacement of nlamfda)
    integer:: nrh  ! number of RH values in table
    integer:: nbin ! number of size bins in table
    integer:: nMom ! number of moments of phase function
    integer:: nPol ! number of elements of scattering phase matrix

!   c=channel, r=rh, b=bin, m=moments, p=nPol
    real,dimension(:),pointer        :: wavelengths => null() ! (c) wavelengths [m]
    real,dimension(:),pointer        :: rh          => null() ! (r) RH values   [fraction]
    real,dimension(:,:),pointer      :: reff        => null() ! (r,b) effective radius [m]
    real,dimension(:,:,:),pointer    :: bext        => null() ! (r,c,b) bext values [m2 kg-1]
    real,dimension(:,:,:),pointer    :: bsca        => null() ! (r,c,b) bsca values [m2 kg-1]
    real,dimension(:,:,:),pointer    :: bbck        => null() ! (r,c,b) bbck values [m2 kg-1]
    real,dimension(:,:,:),pointer    :: g           => null() ! (r,c,b) asymmetry parameter
!   real,dimension(:,:,:,:),pointer  :: pback       => null() ! (r,c,b,p) Backscatter phase function
    real,dimension(:,:,:),pointer    :: p11         => null() ! (r,c,b) Backscatter phase function, index 1 
    real,dimension(:,:,:),pointer    :: p22         => null() ! (r,c,b) Backscatter phase function, index 5
    real,dimension(:,:,:,:,:),pointer:: pmom        => null() ! (r,c,b,m,p) moments of phase function
    real,dimension(:,:),pointer      :: gf          => null() ! (r,b) hygroscopic growth factor
    real,dimension(:,:),pointer      :: rhop        => null() ! (r,b) wet particle density [kg m-3]
    real,dimension(:,:),pointer      :: rhod        => null() ! (r,b) wet particle density [kg m-3]
    real,dimension(:,:),pointer      :: vol         => null() ! (r,b) wet particle volume [m3 kg-1]
    real,dimension(:,:),pointer      :: area        => null() ! (r,b) wet particle cross section [m2 kg-1]
    real,dimension(:,:,:),pointer    :: refr        => null() ! (r,c,b) real part of refractive index
    real,dimension(:,:,:),pointer    :: refi        => null() ! (r,c,b) imaginary part of refractive index

    integer,dimension(NRH_BINS):: rhi ! pointer to rh LUT
    real,dimension(NRH_BINS)   :: rha ! slope on rh LUT


    contains
       procedure :: QueryByWavelength_1d
       procedure :: QueryByWavelength_2d
       procedure :: QueryByWavelength_3d
       procedure :: QueryByChannel_1d
       procedure :: QueryByChannel_2d
       procedure :: QueryByChannel_3d
       generic   :: Query => QueryByWavelength_1d, &
                             QueryByWavelength_2d, &
                             QueryByWavelength_3d, &
                             QueryByChannel_1d,    &
                             QueryByChannel_2d,    &
                             QueryByChannel_3d

       procedure:: getChannel
       procedure:: getWavelength
 end type GOCART2G_Mie

 interface GOCART2G_Mie
    module procedure GOCART2G_MieCreate
 end interface GOCART2G_Mie

 contains


!=================================================================================================================
 type(GOCART2G_Mie) function GOCART2G_MieCreate(dminfo,MieFile,wavelengths,nmom) result(self)
 use SMIOLf
#include "smiol_codes.inc"
!=================================================================================================================

!--- input arguments:
 type(dm_info),intent(in):: dminfo

 character(len=*),intent(in):: MieFile ! Mie table file name
 integer,intent(in),optional:: nmom
 real,intent(in),dimension(:), optional:: wavelengths

!--- local variables:
 integer:: nmom_,nPol
 integer:: stat
 type(SMIOLf_context),pointer :: context
 type(SMIOLf_file),pointer    :: aop_file
 type(SMIOLf_decomp),pointer  :: decomp   ! not used for non-decomposed variables


 integer (kind=I8KIND):: radius_size,rh_size,lambda_size
 real(kind=R4KIND),dimension(:),pointer    :: rh,lambda,radius,rLow,rUp
 real(kind=R4KIND),dimension(:,:),pointer  :: rEff,rMass
 real(kind=R4KIND),dimension(:,:,:),pointer:: qsca,qext,bsca,bext,g,bbck,refreal,refimag

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter function GOCART2G_MieCreate_smiol:')


!
!--- whether or not we are doing phase function:
!
 if(present(nmom)) then
    nmom_ = nmom 
 else
    nmom_ = 0
 endif


!--- set up nPol for reading backscatter phase function
 nPol = 6
      

!
!--- set up a context, within which SMIOL can do parallel I/O:
!
#ifdef MPAS_USE_MPI_F08
 stat = SMIOLf_init(dminfo%comm%mpi_val,1,dminfo%nprocs,context)
#else
 stat = SMIOLf_init(dminfo%comm,1,dminfo%nprocs,context)
#endif
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error initializing SMIOL context', messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    return
 endif


!
!--- open the netCDF file to be read:
!
 stat = SMIOLf_open_file(context,trim(MieFile),SMIOL_FILE_READ,aop_file)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error opening file opticsBands_SU.v1_3.RRTMG.nc', messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    stat = SMIOLf_finalize(context)
    return
 endif


!
!--- inquire about the size of dimensions:
!
 stat = SMIOLf_inquire_dim(aop_file,'radius',dimsize=radius_size)
 stat = SMIOLf_inquire_dim(aop_file,'rh',dimsize=rh_size)
 stat = SMIOLf_inquire_dim(aop_file,'lambda',dimsize=lambda_size)

 call mpas_log_write('radius = $i',intArgs=[int(radius_size)])
 call mpas_log_write('rh     = $i',intArgs=[int(rh_size)])
 call mpas_log_write('lambda = $i',intArgs=[int(lambda_size)])

 self%nch  = lambda_size
 self%nrh  = rh_size
 self%nbin = radius_size


!
!--- allocate arrays in netCDF file:
!
 if(.not.associated(rh)     ) allocate(rh(rh_size)        )
 if(.not.associated(lambda) ) allocate(lambda(lambda_size))
 if(.not.associated(radius) ) allocate(radius(radius_size))
 if(.not.associated(rLow)   ) allocate(rLow(radius_size)  )
 if(.not.associated(rUp)    ) allocate(rUp(radius_size)   )
 if(.not.associated(rEff)   ) allocate(rEff(radius_size,rh_size) )
 if(.not.associated(rMass)  ) allocate(rMass(radius_size,rh_size))
 if(.not.associated(qsca)   ) allocate(qsca(radius_size,rh_size,lambda_size))
 if(.not.associated(qext)   ) allocate(qext(radius_size,rh_size,lambda_size))
 if(.not.associated(bsca)   ) allocate(bsca(radius_size,rh_size,lambda_size))
 if(.not.associated(bext)   ) allocate(bext(radius_size,rh_size,lambda_size))
 if(.not.associated(g)      ) allocate(g(radius_size,rh_size,lambda_size)   )
 if(.not.associated(bbck)   ) allocate(bbck(radius_size,rh_size,lambda_size))
 if(.not.associated(refreal)) allocate(refreal(radius_size,rh_size,lambda_size))
 if(.not.associated(refimag)) allocate(refimag(radius_size,rh_size,lambda_size))


!
!--- read the variable qsca as a non-decomposed variable
!    i.e., every MPI tasks reads the full variable:
!
 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'rh',decomp,rh)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable rh',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'lambda',decomp,lambda)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable lambda',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'radius',decomp,radius)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable radius',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'rLow',decomp,rLow)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable rLow',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'rUp',decomp,rUp)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable rUp',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'rEff',decomp,rEff)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable rEff',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'rMass',decomp,rMass)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable rMass',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'qsca',decomp,qsca)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable qsca',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'qext',decomp,qext)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable qext',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'bsca',decomp,bsca)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable bsca',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'bext',decomp,bext)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable bext',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'g',decomp,g)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable g',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'bbck',decomp,bbck)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable bbck',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

!--- wet particle real part of refractive index:
 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'refreal',decomp,refreal)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable refreal',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

!--- wet particle imaginary part of refractive index (ensure positive):
 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'refimag',decomp,refimag)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable refimag',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif


 end function GOCART2G_MieCreate

!=================================================================================================================
!--- Query subroutines:
#define RANK_ 1
#include "MieQuery.H"
#undef RANK_

#define RANK_ 2
#include "MieQuery.H"
#undef RANK_

#define RANK_ 3
#include "MieQuery.H"
#undef RANK_

!=================================================================================================================
 integer function getChannel(this, wavelength, rc) result (ch)
 class (GOCART2G_Mie), intent(in) :: this
 real, intent(in) :: wavelength
 integer, optional, intent(out) :: rc
 real, parameter :: w_tol = 1.e-9
 integer :: i

 ch = -1
 do i = 1, this%nch
    if (abs(this%wavelengths(i)-wavelength) <= w_tol) then
       ch = i
       exit
    endif
 enddo

 if (present(rc)) rc = 0

 if (ch < 0) then
    !$omp critical (GetCha)
    print*, "wavelength of ",wavelength, " is an invalid value."
    !$omp end critical (GetCha)
    if (present(rc)) rc = -1
 endif

 end function getChannel

!=================================================================================================================
 real function getWavelength(this, ith_channel, rc) result (wavelength)
 class (GOCART2G_Mie), intent(in) :: this
 integer, intent(in) :: ith_channel
 integer, optional, intent(out) :: rc
 real, parameter :: w_tol = 1.e-9
 integer :: i

 if (present(rc)) rc = 0

 if (ith_channel <=0 .or. ith_channel > this%nch ) then
    !$omp critical (GetWav)
    print*, "The channel of ",ith_channel, " is an invalid channel number."
    !$omp end critical (GetWav)
    if (present(rc)) rc = -1
    wavelength = -1. ! meanlingless nagative
    return
 endif

  wavelength = this%wavelengths(ith_channel)

  end function getWavelength

!=================================================================================================================
 end module gocart2G_MieMod_smiol
!=================================================================================================================
