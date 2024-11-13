! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module gocart2G_MieMod_smiol
 use mpas_kind_types, only: R4KIND,RKIND
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

    !c=channel, r=rh, b=bin, m=moments, p=nPol
    !--- pointers available in all netCDF files:
    real,dimension(:),pointer        :: rh          => null() ! (r) RH values   [fraction]
    real,dimension(:,:),pointer      :: reff        => null() ! (r,b) effective radius [m]
    real,dimension(:,:,:),pointer    :: bext        => null() ! (r,c,b) bext values [m2 kg-1]
    real,dimension(:,:,:),pointer    :: bsca        => null() ! (r,c,b) bsca values [m2 kg-1]
    real,dimension(:,:,:),pointer    :: bbck        => null() ! (r,c,b) bbck values [m2 kg-1]
    real,dimension(:,:,:),pointer    :: g           => null() ! (r,c,b) asymmetry parameter
    real,dimension(:,:,:),pointer    :: refr        => null() ! (r,c,b) real part of refractive index
    real,dimension(:,:,:),pointer    :: refi        => null() ! (r,c,b) imaginary part of refractive index

    !--- pointers available in netCDF files with the "wavelength" option:
    real,dimension(:),pointer        :: wavelengths => null() ! (c) wavelengths [m]
    real,dimension(:,:,:,:),pointer  :: pback       => null() ! (r,c,b,m,p) backscatter phase function
    real,dimension(:,:,:,:,:),pointer:: pmom        => null() ! (r,c,b,m,p) moments of phase function

    !--- pointers (and derived pointers) that are sometimes available in netCDF with the "wavelength" option:
    real,dimension(:,:),pointer      :: gf          => null() ! (r,b) hygroscopic growth factor
    real,dimension(:,:),pointer      :: rhop        => null() ! (r,b) wet particle density [kg m-3]
    real,dimension(:,:),pointer      :: rhod        => null() ! (r,b) wet particle density [kg m-3]
    real,dimension(:,:),pointer      :: vol         => null() ! (r,b) wet particle volume [m3 kg-1]
    real,dimension(:,:),pointer      :: area        => null() ! (r,b) wet particle cross section [m2 kg-1]

    real,dimension(:,:,:),pointer    :: p11         => null() ! (r,c,b) backscatter phase function, index 1
    real,dimension(:,:,:),pointer    :: p22         => null() ! (r,c,b) backscatter phase function, index 5

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
 type(GOCART2G_Mie) function GOCART2G_MieCreate(dminfo,MieFile,wavelengths) result(self)
 use SMIOLf
#include "smiol_codes.inc"
!=================================================================================================================

!--- input arguments:
 type(dm_info),intent(in):: dminfo

 character(len=*),intent(in):: MieFile ! Mie table file name
 real,intent(in),dimension(:), optional:: wavelengths

!--- local variables:
 integer:: i,imom,ip1,ipol,j,n,nn,nMom,nPol
 integer:: stat,ndims
 type(SMIOLf_context),pointer :: context
 type(SMIOLf_file),pointer    :: aop_file
 type(SMIOLf_decomp),pointer  :: decomp   ! not used for non-decomposed variables

 real(kind=RKIND):: yerr
 real(kind=R4KIND),parameter:: missing = -999.

 integer(kind=I8KIND):: radius_size,rh_size,lambda_size,nMom_size,nPol_size
!the arrays below are available in all the opticsBands_*.nc and optics_*.nc files for BC,BR,DU,NI,OC,SS,and SU:
 real(kind=R4KIND),dimension(:),pointer:: rh,lambda,radius,rLow,rUp
 real(kind=R4KIND),dimension(:,:),pointer:: rEff,rMass
 real(kind=R4KIND),dimension(:,:,:),pointer:: qsca,qext,bsca,bext,g,bbck,refreal,refimag

!in addition to the arrays above,the arrays pback and pmom are available in all the optics_*.nc files for BC,
!BR,DU,NI,OC,SS,and SU:
 real(kind=R4KIND),dimension(:,:,:,:),pointer:: pback
 real(kind=R4KIND),dimension(:,:,:,:,:),pointer:: pmom

!in addition to the arrays above the arrays gf and rhop are available in the optics_*.nc files for BR and NI:
 real(kind=R4KIND),dimension(:,:),pointer:: gf,rhop

!extra arrays not always available in netCDF files:
 real(kind=RKIND),dimension(:,:),pointer:: rhod,vol,area

!arrays needed to call subroutine polint:
 real(kind=RKIND),dimension(:),allocatable:: lambda_r,input_r

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter function GOCART2G_MieCreate_smiol:')


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

 nMom = 0
 if(present(wavelengths)) then
    stat = SMIOLf_inquire_dim(aop_file,'nMom',dimsize=nMom_size)
    nMom = int(nMom_size)
 else
    nMom_size = 0
 endif

 nPol = 6
 if(present(wavelengths)) then
    stat = SMIOLf_inquire_dim(aop_file,'nPol',dimsize=nPol_size)
    nPol = int(nPol_size)
 else
    nPol_size = nPol
 endif

 call mpas_log_write('radius    = $i',intArgs=[int(radius_size)])
 call mpas_log_write('rh        = $i',intArgs=[int(rh_size)])
 call mpas_log_write('lambda    = $i',intArgs=[int(lambda_size)])
 call mpas_log_write('nMom      = $i',intArgs=(/nMom/))
 call mpas_log_write('nPol      = $i',intArgs=(/nPol/))
 call mpas_log_write('nMom_size = $i',intArgs=(/int(nMom_size)/))
 call mpas_log_write('nPol_size = $i',intArgs=(/int(nPol_size)/))


!
!--- allocate arrays in netCDF files:
!
 if(.not.associated(rh)     ) allocate(rh(rh_size)        )
 if(.not.associated(lambda) ) allocate(lambda(lambda_size))
 if(.not.associated(radius) ) allocate(radius(radius_size))
 if(.not.associated(rLow)   ) allocate(rLow(radius_size)  )
 if(.not.associated(rUp)    ) allocate(rUp(radius_size)   )

 if(.not.associated(rEff)   ) allocate(rEff(rh_size,radius_size) )
 if(.not.associated(rMass)  ) allocate(rMass(rh_size,radius_size))
 if(.not.associated(gf)     ) allocate(gf(rh_size,radius_size)   )
 if(.not.associated(rhop)   ) allocate(rhop(rh_size,radius_size) )
 if(.not.associated(rhod)   ) allocate(rhod(rh_size,radius_size) )
 if(.not.associated(area)   ) allocate(area(rh_size,radius_size) )
 if(.not.associated(vol)    ) allocate(vol(rh_size,radius_size)  )

 if(.not.associated(qsca)   ) allocate(qsca(lambda_size,rh_size,radius_size))
 if(.not.associated(qext)   ) allocate(qext(lambda_size,rh_size,radius_size))
 if(.not.associated(bsca)   ) allocate(bsca(lambda_size,rh_size,radius_size))
 if(.not.associated(bext)   ) allocate(bext(lambda_size,rh_size,radius_size))
 if(.not.associated(g)      ) allocate(g(lambda_size,rh_size,radius_size)   )
 if(.not.associated(bbck)   ) allocate(bbck(lambda_size,rh_size,radius_size))
 if(.not.associated(refreal)) allocate(refreal(lambda_size,rh_size,radius_size))
 if(.not.associated(refimag)) allocate(refimag(lambda_size,rh_size,radius_size))

 if(.not.associated(pback)  ) allocate(pback(lambda_size,rh_size,radius_size,nPol_size))


!
!--- read variable as a non-decomposed variable
!    i.e., every MPI tasks reads the full variable:
!
 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'rh',decomp,rh)
 do n = 1,rh_size
    call mpas_log_write('--- rh: n = $i $r',intArgs=(/n/),realArgs=(/rh(n)/))
 enddo
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable rh',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif
 

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'lambda',decomp,lambda)
 do n = 1,lambda_size
    call mpas_log_write('--- lambda: n = $i $r',intArgs=(/n/),realArgs=(/lambda(n)/))
 enddo
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

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'refreal',decomp,refreal)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable refreal',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 nullify(decomp)
 stat = SMIOLf_get_var(aop_file,'refimag',decomp,refimag)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error reading variable refimag',messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
    stat = SMIOLf_close_file(aop_file)
    stat = SMIOLf_finalize(context)
    return
 endif

 if(present(wavelengths)) then
    if(.not.associated(pmom)   ) allocate(pmom(lambda_size,rh_size,radius_size,nMom,nPol_size))

    nullify(decomp)
    stat = SMIOLf_get_var(aop_file,'pmom',decomp,pmom)
    if(stat /= SMIOL_SUCCESS) then
       call mpas_log_write('Error reading variable pmom',messageType=MPAS_LOG_ERR)
       call mpas_log_write(trim(SMIOLf_error_string(stat)),messageType=MPAS_LOG_ERR)
       stat = SMIOLf_close_file(aop_file)
       stat = SMIOLf_finalize(context)
       return
    endif
 endif


!--- AOPs not necessarily stored in all versions of netCDF files:
 if(present(wavelengths)) then
    !particle growth factor:
    nullify(decomp)
    stat = SMIOLf_inquire_var(aop_file,'growth_factor',ndims=ndims)
    if(stat /= SMIOL_SUCCESS) then
       gf(:,:) = missing
!      call mpas_log_write('--- GROWTH FACTOR GF is not available in input file',messageType=MPAS_LOG_OUT)
    else
       stat = SMIOLf_get_var(aop_file,'growth_factor',decomp,gf)
    endif

    !wet particle density:
    nullify(decomp)
    stat = SMIOLf_inquire_var(aop_file,'rhop',ndims=ndims)
    if(stat /= SMIOL_SUCCESS) then
       rhop(:,:) = missing
!      call mpas_log_write('--- WET PARTICLE DENSITY not available in input file',messageType=MPAS_LOG_OUT)
    else
       stat = SMIOLf_get_var(aop_file,'rhop',decomp,rhop)
    endif

    !dry particle density (pulled from wet particle radius):
    nullify(decomp)
    stat = SMIOLf_inquire_var(aop_file,'rhop',ndims=ndims)
    if(stat /= SMIOL_SUCCESS) then
       rhod(:,:) = missing
!      call mpas_log_write('--- DRY PARTICLE DENSITY not available in input file',messageType=MPAS_LOG_OUT)
    else
       stat = SMIOLf_get_var(aop_file,'rhop',decomp,rhod)
       do n = 1, rh_size
          rhod(n,:) = rhod(1,:)
       enddo
    endif

    !--- backscatter phase function:
    nullify(decomp)
    stat = SMIOLf_inquire_var(aop_file,'pback',ndims=ndims)
    if(stat /= SMIOL_SUCCESS) then
       pback(:,:,:,:) = 1._RKIND
!      call mpas_log_write('--- BACKSCATTER PHASE FUNCTION not available in input file',messageType=MPAS_LOG_OUT)
    else
       stat = SMIOLf_get_var(aop_file,'pback',decomp,pback)
    endif

    !--- wet particle volume [m3 kg-1]. the ratio of wet to dry volume is gf^3, hence the following
    do n = 1, rh_size
       do nn = 1, radius_size
          if(rhod(n,nn) == missing .or. gf(n,nn) == missing) then
             vol(n,nn) = missing
          else
             vol(n,nn) = gf(n,nn)**3/rhod(n,nn)
          endif
       enddo
    enddo

    !--- wet particle cross sectional area [m2 kg-1]. assume area is volume divided by (4./3.*reff)
    do n = 1, rh_size
       do nn = 1, radius_size
          if(rhod(n,nn) == missing) then
             area(n,nn) = missing
          else
             area(n,nn) = vol(n,nn)/(4./3.*rEff(n,nn))
          endif
       enddo
    enddo
 endif


!--- close netCDF file:
 stat = SMIOLf_close_file(aop_file)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    stat = SMIOLf_finalize(context)
    return
 endif

 stat = SMIOLf_finalize(context)
 if(stat /= SMIOL_SUCCESS) then
    call mpas_log_write('Error finalizing SMIOL context', messageType=MPAS_LOG_ERR)
    call mpas_log_write(trim(SMIOLf_error_string(stat)), messageType=MPAS_LOG_ERR)
    return
 endif


!--- output data to GOCART2G_Mie:
 self%nrh  = rh_size
 self%nbin = radius_size
 self%nMom = nMom
 self%nPol = nPol

 if(present(wavelengths)) then
    self%nch = size(wavelengths)
 else
    self%nch = lambda_size
 endif
 if(.not.associated(self%wavelengths)) allocate(self%wavelengths(self%nch))
 if(present(wavelengths)) then
    self%wavelengths = wavelengths
 else
    self%wavelengths = lambda
 endif

 call mpas_log_write(' ')
 call mpas_log_write('self%nbin = $i',intArgs=(/self%nbin/))
 call mpas_log_write('self%nrh  = $i',intArgs=(/self%nrh/))
 call mpas_log_write('self%nch  = $i',intArgs=(/self%nch/))
 call mpas_log_write('self%nMom = $i',intArgs=(/self%nMom/))
 call mpas_log_write('self%nPol = $i',intArgs=(/self%nPol/))

 if(.not.associated(self%rh)   ) allocate(self%rh(self%nrh))
 if(.not.associated(self%reff) ) allocate(self%reff(self%nrh,self%nbin))
 if(.not.associated(self%gf)   ) allocate(self%gf(self%nrh,self%nbin)  )
 if(.not.associated(self%rhop) ) allocate(self%rhop(self%nrh,self%nbin))
 if(.not.associated(self%rhod) ) allocate(self%rhod(self%nrh,self%nbin))
 if(.not.associated(self%vol)  ) allocate(self%vol(self%nrh,self%nbin) )
 if(.not.associated(self%area) ) allocate(self%area(self%nrh,self%nbin))
 if(.not.associated(self%bext) ) allocate(self%bext(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%bsca) ) allocate(self%bsca(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%bbck) ) allocate(self%bbck(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%refr) ) allocate(self%refr(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%refi) ) allocate(self%refi(self%nrh,self%nch,self%nbin))
 if(.not.associated(self%g)    ) allocate(self%g(self%nrh,self%nch,self%nbin)   )
 if(.not.associated(self%p11)  ) allocate(self%p11(self%nrh,self%nch,self%nbin) )
 if(.not.associated(self%p22)  ) allocate(self%p22(self%nrh,self%nch,self%nbin) )
 if(.not.associated(self%pback)) allocate(self%pback(self%nrh,self%nch,self%nbin,self%nPol))
 if(nMom_size > 0) then
    if(.not.associated(self%pmom)) allocate(self%pmom(self%nrh,self%nch,self%nbin,self%nMom,self%nPol))
 endif

 self%rh   = real(rh,kind=RKIND)   ! relative humidity (fraction).
 self%reff = real(rEff,kind=RKIND) ! effective radius of bin (m).
 self%gf   = real(gf,kind=RKIND)   ! growth factor.
 self%rhop = real(rhop,kind=RKIND) ! wet particle density (kg m^-3).
 self%rhod = real(rhod,kind=RKIND) ! dry particle density (kg m^-3).
 self%vol  = real(vol,kind=RKIND)  ! volume (m^3 kg^-1).
 self%area = real(area,kind=RKIND) ! area (m^2 kg^-1).

 if(present(wavelengths)) then
    if(.not.allocated(input_r) ) allocate(input_r(int(lambda_size)) )
    if(.not.allocated(lambda_r)) allocate(lambda_r(int(lambda_size)))
    do nn = 1,int(lambda_size)
       lambda_r(nn) = real(lambda(nn),kind=RKIND)
    enddo
    do j = 1,self%nbin
       do i = 1,self%nrh
          do n = 1,self%nch
             input_r(:) = real(bext(:,i,j),kind=RKIND)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%bext(i,n,j),yerr)
             input_r(:) = real(bsca(:,i,j),kind=RKIND)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%bsca(i,n,j),yerr)
             input_r(:) = real(bbck(:,i,j),kind=RKIND)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%bbck(i,n,j),yerr)
             input_r(:) = real(g(:,i,j),kind=RKIND)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%g(i,n,j)   ,yerr)
             input_r(:) = real(refreal(:,i,j),kind=RKIND)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%refr(i,n,j),yerr)
             input_r(:) = real(refimag(:,i,j),kind=RKIND)
             call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),self%refi(i,n,j),yerr)

             do ipol = 1,self%nPol
                input_r(:) = real(pback(:,i,j,ipol),kind=RKIND)
                call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n),pback(i,n,j,ipol),yerr)
             enddo

             if(nMom > 0) then
                do imom = 1,self%nMom
                   do ipol = 1,self%nPol
                      input_r(:) = real(pmom(:,i,j,imom,ipol),kind=RKIND)
                      call polint(lambda_r,input_r,int(lambda_size),self%wavelengths(n), &
                                  self%pmom(i,n,j,imom,ipol),yerr)
                   enddo
                enddo
             endif

!--- originally written sourcecode:
!            call polint(lambda,bext(:,i,j)   ,lambda_size,self%wavelengths(n),self%bext(i,n,j),yerr)
!            call polint(lambda,bsca(:,i,j)   ,lambda_size,self%wavelengths(n),self%bsca(i,n,j),yerr)
!            call polint(lambda,bbck(:,i,j)   ,lambda_size,self%wavelengths(n),self%bbck(i,n,j),yerr)
!            call polint(lambda,g(:,i,j)      ,lambda_size,self%wavelengths(n),self%g(i,n,j)   ,yerr)
!            call polint(lambda,refreal(:,i,j),lambda_size,self%wavelengths(n),self%refr(i,n,j),yerr)
!            call polint(lambda,refimag(:,i,j),lambda_size,self%wavelengths(n),self%refi(i,n,j),yerr)
!            do ipol = 1,self%nPol
!               call polint(lambda,pback(:,i,j,ipol),lambda_size,self%wavelengths(n),pback(i,n,j,ipol),yerr)
!            enddo

!            if(nMom > 0) then
!               do imom = 1,self%nMom
!                  do ipol = 1,self%nPol
!                     call polint(lambda,pmom(:,i,j,imom,ipol),lambda_size,self%wavelengths(n), &
!                                 self%pmom(i,n,j,imom,ipol),yerr)
!                  enddo
!               enddo
!            endif

          enddo
       enddo
    enddo
    if(allocated(lambda_r)) deallocate(lambda_r)
    if(allocated(input_r) ) deallocate(input_r )
 else
    !--- swap the order:
    self%bext  = reshape(real(bext,kind=RKIND)   ,[int(rh_size),int(lambda_size),int(radius_size)],order =[2,1,3])
    self%bsca  = reshape(real(bsca,kind=RKIND)   ,[int(rh_size),int(lambda_size),int(radius_size)],order =[2,1,3])
    self%bbck  = reshape(real(bbck,kind=RKIND)   ,[int(rh_size),int(lambda_size),int(radius_size)],order =[2,1,3])
    self%g     = reshape(real(g,kind=RKIND)      ,[int(rh_size),int(lambda_size),int(radius_size)],order =[2,1,3])
    self%refr  = reshape(real(refreal,kind=RKIND),[int(rh_size),int(lambda_size),int(radius_size)],order =[2,1,3])
    self%refi  = reshape(real(refimag,kind=RKIND),[int(rh_size),int(lambda_size),int(radius_size)],order =[2,1,3])
    self%pback = reshape(real(pback,kind=RKIND)  , &
                 [int(rh_size),int(lambda_size),int(radius_size),int(nPol_size)],order =[2,1,3,4])
    if(nMom_size > 0 ) then
       self%pmom = reshape(real(pmom,kind=RKIND),[int(rh_size),int(lambda_size),int(radius_size), &
                   int(nMom_size),int(nPol_size)],order = [2,1,3,4,5])
    endif

!--- originally written sourcecode:
!   self%bext  = reshape(bext   ,[rh_size,lambda_size,radius_size],order =[2,1,3])
!   self%bsca  = reshape(bsca   ,[rh_size,lambda_size,radius_size],order =[2,1,3])
!   self%bbck  = reshape(bbck   ,[rh_size,lambda_size,radius_size],order =[2,1,3])
!   self%g     = reshape(g      ,[rh_size,lambda_size,radius_size],order =[2,1,3])
!   self%refr  = reshape(refreal,[rh_size,lambda_size,radius_size],order =[2,1,3])
!   self%refi  = reshape(refimag,[rh_size,lambda_size,radius_size],order =[2,1,3])
!   self%pback = reshape(pback,[rh_size,lambda_size,radius_size,nPol_size],order =[2,1,3,4])
!   if(nMom_size > 0 ) then
!      self%pmom = reshape(pmom,[rh_size,lambda_size,radius_size,nMom_size,nPol_size],order = [2,1,3,4,5])
!   endif
 endif

 self%p11 = self%pback(:,:,:,1)
 self%p22 = self%pback(:,:,:,5)


!--- remapping of RH to some high resolution representation for later use. RH input is scaled 0 - 0.99.
!    we resolve the map to 0 - 0.990 in steps of 0.001 (991 total steps):
 do j = 1, NRH_BINS
    do i = self%nrh, 1, -1
       if((j-1) .ge. int(self%rh(i)*1000)) then
          ip1 = i + 1
          self%rhi(j) = i
          if(ip1 .gt. self%nrh) then
             self%rha(j) = 0.
          else
             self%rha(j) = ( (j-1)/1000.-self%rh(i)) /  (self%rh(ip1)-self%rh(i))
          endif
          exit
       endif
    enddo
 enddo


 call mpas_log_write('--- end function GOCART2G_MieCreate_smiol:')
 return


 contains


!-----------------------------------------------------------------------------------------------------------------
 subroutine polint(x,y,n,xWant,yWant,yErr)
 integer,intent(in):: n
!recall, table hard-wired single precision
 real(kind=RKIND),intent(in):: x(n),y(n)
 real(kind=RKIND),intent(inout):: xWant,yWant,yErr

!given array x(n) of independent variables and array y(n) of dependent
!variables, compute the linear interpolated result yWant at xWant and return
!with a dummy error estimate yErr.  Hacked up from Numerical Recipes Chapter 3

 character(len=255):: msg
 integer:: i, j
 real(kind=RKIND):: dx, slope

!---on out of bounds, set i to lower or upper limit:
 i = 0
 if(xWant .lt. x(1)) then
    i = 1
 endif
 if(xWant .gt. x(n)) then
    i = n
 endif

!--- if i is still zero find i less than xWant:
 if(i .eq. 0) then
    do j = 1, n
       if(xWant .ge. x(j)) i = j
    enddo
 endif

!--- slope:
 if(i .eq. n) then
    slope = 0.
 else
    slope = (y(i+1)-y(i)) / (x(i+1)-x(i))
 endif
 dx = xWant - x(i)
 yWant = y(i) + slope*dx

 yErr = 0.

 end subroutine polint
!-----------------------------------------------------------------------------------------------------------------

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
