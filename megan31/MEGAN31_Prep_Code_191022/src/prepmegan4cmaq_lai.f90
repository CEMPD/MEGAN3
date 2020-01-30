!--------------------------------------------------------------------
! Based of bio_emiss.f90_newef to only process LAI
! 
! Serena H. Chung
! 2014-05-02
!
!--------------------------------------------------------------------

   program map_megan2_emissions

   use area_mapper, only : xlong => lon, xlat => lat
   use bio_types

   implicit none

   INTEGER, parameter :: PS = 2                             ! Polar Stereographic map projection
   INTEGER, parameter :: lower = 0
   INTEGER, parameter :: upper = 1
   INTEGER :: icnt
   INTEGER :: ids, ide, jds, jde
   INTEGER :: i, j, n
   INTEGER :: ii, jj
   integer :: ncid
   integer :: cell
   integer :: map_proj
   integer :: ierr, astat
   integer :: dimid, varid
   integer :: nlon_megan, nlat_megan
   integer :: domain
   integer :: mnth, mnth_s, mnth_e, megan_month
   integer :: start_lai_mnth = 1
   integer :: end_lai_mnth   = 12
   integer :: x0,y0,ncolsin,nrowsin !shc
   integer :: xndx_megan(2)
   integer :: yndx_megan(2)
   integer, allocatable :: ix(:,:,:)                        ! index used by interpolation
   integer, allocatable :: jy(:,:,:)                        ! index used by interpolation

   real    :: missing_value
   real    :: scale_factor
   real    :: wrk_sum
   real    :: ds1, ds2
   real    :: xl, xu
   real    :: yl, yu, dy
   real    :: wrf_lon_min
   real    :: wrf_lon_max
   real    :: wrf_lat_min
   real    :: wrf_lat_max
   real    :: cen_lon
   real    :: cen_lat
   real    :: stand_lon
   real    :: truelat1
   real    :: truelat2
   real    :: dx
   real(8), allocatable :: xedge_megan(:)
   real(8), allocatable :: yedge_megan(:)
   real, allocatable :: megan_lons(:)
   real, allocatable :: megan_lats(:)

!-----allocate lai, tsa, and swdown arrays----------- 
   real, allocatable :: mlai(:,:,:)
   real, allocatable :: tmp3(:,:,:)
   real, allocatable :: ax(:,:,:)                        ! weight coef. all domain
   real, allocatable :: by(:,:,:)                        ! weight coef. all domain

   CHARACTER (LEN=132) :: varname
   CHARACTER (LEN=300) :: filespec
   CHARACTER (LEN=80)  :: message
   CHARACTER (LEN=80)  :: attribute
   CHARACTER (LEN=80)  :: units_attribute
   CHARACTER (LEN=80)  :: description_attribute
   CHARACTER (LEN=80)  :: stagger_attribute
   CHARACTER (LEN=80)  :: coor_attribute
   CHARACTER (LEN=80)  :: memord_attribute
   CHARACTER (LEN=80)  :: inpname
   CHARACTER (LEN=300)  :: outpname
   CHARACTER (LEN=300)  :: wrffile
   CHARACTER (LEN=300)  :: megan_dir
   CHARACTER (LEN=300)  :: out_dir !shc
!  CHARACTER (LEN=19)  :: Times(1)
   CHARACTER (LEN=4)   :: num  !used to calcualte LAI date
   CHARACTER (LEN=3)   :: char_mnth(12)

   logical :: has_area_map
   logical :: new_grid

   namelist /control/ start_lai_mnth, end_lai_mnth, &
                      wrffile, megan_dir, out_dir !shc added out_dir
   namelist /windowdefs/ x0, y0, ncolsin, nrowsin !shc

!---------------------------------------------------------------------
!	... include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

   char_mnth(:) = (/ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)
!   wrf_dir   = '.'
!   megan_dir = '.'

!-----------------------------------------------------------------
!     read namelist variables
!-----------------------------------------------------------------
   read(*,nml=control,iostat=ierr)
   write(*,*) start_lai_mnth, end_lai_mnth, wrffile, megan_dir
   if( ierr /= 0 ) then
     write(*,*) 'convert_emissions: failed to read namelist; error = ',ierr
     stop 'bio_emiss abort'
   endif
   !shc the following is added for trimming down the output to
   !    MCIP/CMAQ domain
   read(*,nml=windowdefs,iostat=ierr)
   write(*,*) x0, y0, ncolsin, nrowsin
   if( ierr /= 0 ) then
     write(*,*) 'convert_emissions: failed to read namelist; error = ',ierr
     stop 'bio_emiss abort'
   endif
!-----------------------------------------------------------------
!     check namelist inputs
!-----------------------------------------------------------------
!  if( domains < 1 ) then
!    write(*,*) 'convert_emissions: domains must be >= 1'
!    stop 'bio_emiss abort'
!  endif
   if( start_lai_mnth < 1 .or. start_lai_mnth > 12 ) then
     write(*,*) 'convert_emissions: start month must be in set [1,12]'
     stop 'bio_emiss abort'
   endif
   if( end_lai_mnth < 1 .or. end_lai_mnth > 12 ) then
     write(*,*) 'convert_emissions: end month must be in set [1,12]'
     stop 'bio_emiss abort'
   endif
   if( end_lai_mnth < start_lai_mnth ) then
     write(*,*) 'convert_emissions: end month must >= start_month'
      stop 'bio_emiss abort'
   endif

!-----------------------------------------------------------------
!     loop over domains
!-----------------------------------------------------------------
!domain_loop : &
!   do domain = 1,domains

!jxy      write(*,*) ' '
!jxy      write(*,*) '========================================='
!jxy      write(*,*) 'Domain = ',domain
!jxy      write(*,*) '========================================='
!jxy      write(*,*) ' '

!-----------------------------------------------------------
!     ... read wrfinput file
!-----------------------------------------------------------
      call wrf_file

!-----------------------------------------------------------
!     ... read and interpolate megan datasets
!     ... need to change missing_value, scale_factor if you 
!     ... have differet values for the datasets
!-----------------------------------------------------------


!-----read LAI----------------------------------------------
!Ling provided the fix to 4km by comment out the following line
!      mnth_e        = 1
      scale_factor  = 1.e-1
      missing_value = -1.
      do megan_month = 1, 46  !there are 46 8-day files 
        write(num,'(i4)') 1000+megan_month*8-7
         write(*,*)'num=',num
         inpname       = 'laina2008' // num(2:4) // '.nc'
         varname       = 'LAIv'
         CALL  megan2_bioemiss
      end do


!-----------------------------------------------------------
!     ... write to output file(s)
!-----------------------------------------------------------
      write(*,*) 'map_megan2_emissions: Before write_lai3'
!write to LAI3.csv file---------------------
      CALL write_lai3

!-----------------------------------------------------------
!     cleanup domain variables
!-----------------------------------------------------------

      deallocate( mlai )
      if( allocated( xlong ) ) then
        deallocate( xlong )
      endif
      if( allocated( xlat ) ) then
        deallocate( xlat )
      endif
!      do n = 1,grid_cnt
!        if( associated( grid_specs(n)%lon ) ) then
!          deallocate( grid_specs(n)%lon )
!        endif
!        if( associated( grid_specs(n)%lat ) ) then
!          deallocate( grid_specs(n)%lat )
!        endif
!        if( associated( grid_specs(n)%model_area_type ) ) then
!          do j = 1,jde
!            do i = 1,ide
!              if( associated( grid_specs(n)%model_area_type(i,j)%dcell_lon_ndx ) ) then
!                deallocate( grid_specs(n)%model_area_type(i,j)%dcell_lon_ndx )
!              endif
!              if( associated( grid_specs(n)%model_area_type(i,j)%dcell_lat_ndx ) ) then
!                deallocate( grid_specs(n)%model_area_type(i,j)%dcell_lat_ndx )
!              endif
!              if( associated( grid_specs(n)%model_area_type(i,j)%wght ) ) then
!                deallocate( grid_specs(n)%model_area_type(i,j)%wght )
!              endif
!            end do
!          end do
!          deallocate( grid_specs(n)%model_area_type )
!        endif
!      end do
!++sw
      grid_cnt = 0
!--sw
!  end do domain_loop

   CONTAINS

   subroutine wrf_file
!---------------------------------------------------------------------
!   read wrf file
!---------------------------------------------------------------------

   use area_mapper, only : proj_init

!  write(num,'(i3)') 100+domain
!  inpname = 'wrfinput_d' // num(2:3)
!  filespec = trim( wrf_dir ) // '/' // trim( inpname )
   filespec = wrffile
!---------------------------------------------------------------------
!   open wrf input file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to open ' // trim(inpname)
   call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get wrf dimesions
!---------------------------------------------------------------------
   message = 'Failed to get lon dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'west_east', dimid ), message )
   message = 'Failed to get lon dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, ide ), message )
   message = 'Failed to get lat dimension id'
   call handle_ncerr( nf_inq_dimid( ncid, 'south_north', dimid ), message )
   message = 'Failed to get lat dimension'
   call handle_ncerr( nf_inq_dimlen( ncid, dimid, jde ), message )
!---------------------------------------------------------------------
!   get wrf map projection variables
!---------------------------------------------------------------------
   message = 'Failed to get map_proj'
   call handle_ncerr( nf_get_att_int( ncid, nf_global, 'MAP_PROJ', map_proj ), message )
   if( map_proj /= PS ) then
      write(*,*) 'wrf_file: MAP_PROJ is not polar stereographic'
   else
      write(*,*) 'wrf_file: MAP_PROJ is polar stereographic'
   endif
   message = 'wrf_file: Failed to get cen_lon'
   call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LON', cen_lon ), message )
   write(*,*) 'wrf_file: CEN_LON = ',cen_lon
   message = 'wrf_file: Failed to get cen_lat'
   call handle_ncerr( nf_get_att_real( ncid, nf_global, 'CEN_LAT', cen_lat ), message )
   write(*,*) 'wrf_file: CEN_LAT = ',cen_lat
   message = 'wrf_file: Failed to get stand_lon'
   call handle_ncerr( nf_get_att_real( ncid, nf_global, 'STAND_LON', stand_lon ), message )
   write(*,*) 'wrf_file: STAND_LON = ',stand_lon
   message = 'wrf_file: Failed to get truelat1'
   call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT1', truelat1 ), message )
   write(*,*) 'wrf_file: TRUELAT1 = ',truelat1
   message = 'wrf_file: Failed to get truelat2'
   call handle_ncerr( nf_get_att_real( ncid, nf_global, 'TRUELAT2', truelat2 ), message )
   write(*,*) 'wrf_file: TRUELAT2 = ',truelat2
   message = 'wrf_file: Failed to get dx'
   call handle_ncerr( nf_get_att_real( ncid, nf_global, 'DX', dx ), message )
   write(*,*) 'wrf_file: DX = ',dx

!---------------------------------------------------------------------
!   initialize map projection
!---------------------------------------------------------------------
   call proj_init( map_proj, cen_lon, cen_lat, truelat1, truelat2, &
                   stand_lon, dx, ide, jde )

   ids = 1
   jds = 1

!  message = 'Failed to get Times id'
!  call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), message )
!  message = 'Failed to read Times'
!  call handle_ncerr( nf_get_var_text( ncid, varid, Times ), message )

   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!  write(*,*) 'wrf_file: time = ',trim(Times(1))
   write(*,*) 'wrf_file: grid dimensions'
   write(*,*) 'wrf_file: ids,ide,jds,jde'
   write(*,'(4i6)') ids,ide,jds,jde
   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!---------------------------------------------------------------------
!   close wrf file
!---------------------------------------------------------------------
   message = 'wrf_file: Failed to close ' // trim(inpname)
   call handle_ncerr( nf_close( ncid ), message )       
!---------------------------------------------------------------------
!   allocate final bioemission variables
!---------------------------------------------------------------------
   allocate( mlai(ide,jde,46),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate mlai; error = ',astat
     stop 'allocate failed'
   endif

   end subroutine wrf_file

   subroutine megan2_bioemiss
!---------------------------------------------------------------------
!   map megan dataset to wrf grid
!---------------------------------------------------------------------

   use area_mapper, only : area_interp
   use constants_module, only : rad_per_deg, earth_radius_m
   use bio_types

!---------------------------------------------------------------------
!   local variables
!---------------------------------------------------------------------
    integer :: il, iu, jl, ju
    integer :: n
    integer :: status
    real    :: wrf_lon, wrf_lat
    real    :: data_dx
    real, allocatable :: wrk_data(:,:)
    logical :: debug = .false.

    write(*,*) ' '
    write(*,*) 'Reading megan2 bio emiss file ' // trim(inpname)
!---------------------------------------------------------------------
!   open megan dataset file
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to open ' // trim(inpname)
    filespec = trim( megan_dir ) // '/' // trim(inpname)
    call handle_ncerr( nf_open( trim(filespec), nf_noclobber, ncid ), message )       
!---------------------------------------------------------------------
!   get megan dataset dimesions
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to get lon dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid ), message )
    message = 'megan2_bioemiss: Failed to get lon dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon_megan ), message )
    message = 'megan2_bioemiss: Failed to get lat dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid ), message )
    message = 'megan2_bioemiss: Failed to get lat dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlat_megan ), message )
    write(*,*) 'megan2_bioemiss:  nlon_megan, nlat_megan = ',nlon_megan,nlat_megan

!---------------------------------------------------------------------
!   allocate working variable
!---------------------------------------------------------------------
    if( allocated( megan_lons ) ) then
      deallocate( megan_lons )
    endif
    allocate( megan_lons(nlon_megan),stat=astat )
    if( astat /= 0 ) then
      write(*,*) 'megan2_bioemiss: Failed to allocate megan_lons; error = ',astat
      stop 'allocate failed'
    endif
!---------------------------------------------------------------------
!   read megan longitude variable
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to get lon variable id'
    call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), message )
    message = 'megan2_bioemiss: Failed to read lon variable'
    call handle_ncerr( nf_get_var_real( ncid, varid, megan_lons ), message )

    if( allocated( megan_lats ) ) then
      deallocate( megan_lats )
    endif
    allocate( megan_lats(nlat_megan),stat=ierr )
    if( ierr /= 0 ) then
      write(*,*) 'megan2_bioemiss: Failed to allocate megan_lats; error = ',ierr
      stop 'allocate failed'
    endif
!---------------------------------------------------------------------
!   read megan latitude variable
!---------------------------------------------------------------------
    message = 'megan2_bioemiss: Failed to get lat variable id'
    call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), message )
    message = 'megan2_bioemiss: Failed to read lat variable'
    call handle_ncerr( nf_get_var_real( ncid, varid, megan_lats ), message )

!---------------------------------------------------------------------
!   determine interpolation type; bilinear or area conserving
!---------------------------------------------------------------------
    data_dx = earth_radius_m * (megan_lats(2) - megan_lats(1)) * rad_per_deg
    has_area_map = data_dx < dx
    write(*,*) 'megan2_bioemiss: data_dx,dx,has_area_map = ',data_dx,dx,has_area_map

!-------------------------------------------------------------
!   check for match against prior datasets
!-------------------------------------------------------------
   if( grid_cnt >= grid_max ) then
     write(*,*) 'megan2_bioemiss: reached grid cache max: ',grid_max
     stop
   endif
   grid_ndx = 0
   new_grid = .true.
   do n = 1,grid_cnt
     if( grid_specs(n)%nlons /= nlon_megan .or. grid_specs(n)%nlats /= nlat_megan ) then
       cycle
     endif
     if( any( grid_specs(n)%lon(:) /= megan_lons(:) ) ) then
       cycle
     endif
     if( any( grid_specs(n)%lat(:) /= megan_lats(:) ) ) then
       cycle
     endif
     grid_ndx = n
     new_grid = .false.
     exit
   end do
!-------------------------------------------------------------
!   new data grid to cache
!-------------------------------------------------------------
   if( new_grid ) then
     grid_cnt = grid_cnt + 1
     grid_specs(grid_cnt)%nlons = nlon_megan
     grid_specs(grid_cnt)%nlats = nlat_megan
     grid_specs(grid_cnt)%has_area_map = has_area_map
     allocate( grid_specs(grid_cnt)%lon(nlon_megan),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'megan2_bioemiss: Failed to allocate megan_lats; error = ',ierr
       stop 'allocate failed'
     endif
     allocate( grid_specs(grid_cnt)%lat(nlat_megan),stat=ierr )
     if( ierr /= 0 ) then
       write(*,*) 'megan2_bioemiss: Failed to allocate megan_lats; error = ',ierr
       stop 'allocate failed'
     endif
     grid_specs(grid_cnt)%lon(:) = megan_lons(:)
     grid_specs(grid_cnt)%lat(:) = megan_lats(:)
     if( has_area_map ) then
       allocate( grid_specs(grid_cnt)%model_area_type(ide,jde),stat=astat )
       if( astat /= 0 ) then
         write(*,*) 'proj_init; failed to allocate model_area_type: error = ',astat
         stop
       endif
       grid_specs(grid_cnt)%model_area_type(:,:)%has_data = .false.
       grid_specs(grid_cnt)%model_area_type(:,:)%active_dcell_cnt = 0
       grid_specs(grid_cnt)%model_area_type(:,:)%total_dcell_cnt  = 0
       grid_specs(grid_cnt)%model_area_type(:,:)%interior_dcell_cnt = 0
       grid_specs(grid_cnt)%model_area_type(:,:)%partial_dcell_cnt  = 0
     endif
     grid_ndx = grid_cnt
     write(*,*) 'megan2_bioemiss: file ' // trim(inpname),' has a new grid'
   endif

is_area_map : &
    if( has_area_map ) then
!---------------------------------------------------------------------
!   form megan longitude edges
!---------------------------------------------------------------------
      allocate( xedge_megan(nlon_megan+1),stat=astat )
      if( astat /= 0 ) then
        write(*,*) 'megan2_bioemiss: Failed to allocate xedge_megan; error = ',astat
        stop 'allocate error'
      endif
      xedge_megan(2:nlon_megan) = .5_8*(megan_lons(1:nlon_megan-1) + megan_lons(2:nlon_megan))
      xedge_megan(1)            = megan_lons(1) - .5_8*(megan_lons(2) - megan_lons(1))
      xedge_megan(nlon_megan+1) = megan_lons(nlon_megan) + .5_8*(megan_lons(nlon_megan) - megan_lons(nlon_megan-1))
      write(*,'(''megan2_bioemiss: xcen_megan(1,2)  = '',1p,2g22.15)') megan_lons(1:2)
      write(*,'(''megan2_bioemiss: xedge_megan(1,2) = '',1p,2g22.15)') xedge_megan(1:2)
      write(*,'(''megan2_bioemiss: dx = '',1pg22.15)') int( 1./(megan_lons(2) - megan_lons(1)) )
!---------------------------------------------------------------------
!   form megan latitude edges
!---------------------------------------------------------------------
      allocate( yedge_megan(nlat_megan+1),stat=ierr )
      if( ierr /= 0 ) then
        write(*,*) 'megan2_bioemiss: Failed to allocate yedge_megan; error = ',ierr
        stop 'allocate error'
      endif

      yedge_megan(2:nlat_megan) = .5_8*(megan_lats(1:nlat_megan-1) + megan_lats(2:nlat_megan))
      yedge_megan(1)            = megan_lats(1) - .5_8*(megan_lats(2) - megan_lats(1))
      yedge_megan(nlat_megan+1) = megan_lats(nlat_megan) + .5_8*(megan_lats(nlat_megan) - megan_lats(nlat_megan-1))

      write(*,'(''megan2_bioemiss: nlon_megan,nlat_megan = '',i6,1x,i6)') nlon_megan,nlat_megan
      write(*,'(''megan2_bioemiss: ycen_megan  = '',1p,2g22.15)') megan_lats(nlat_megan-1:nlat_megan)
      write(*,'(''megan2_bioemiss: yedge_megan = '',1p,2g22.15)') yedge_megan(nlat_megan:nlat_megan+1)

      if( allocated( wrk_data ) ) then
        deallocate( wrk_data )
      endif
      allocate( wrk_data(ide,jde),stat=status )
      if( status /= 0 ) then
        write(*,*) 'megan2_bioemiss: allocate for wrk_data failed; error = ',ierr
        stop 'allocation error'
      endif
!---------------------------------------------------------------------
!   area conserving interpolation
!---------------------------------------------------------------------
      call area_interp( xedge_megan, yedge_megan, nlon_megan, nlat_megan, int(missing_value,2), &
                        wrk_data, ncid, varname, grid_ndx, new_grid )
      if( varname(:3) == 'LAI' ) then
      write(*,*)'megan_month=',megan_month
        mlai(:,:,megan_month) = scale_factor * wrk_data(:,:)
      endif

      deallocate( wrk_data )
    else is_area_map
!---------------------------------------------------------------------
!   form megan coordinate limits
!---------------------------------------------------------------------
      wrf_lon_min = minval( xlong(ids:ide,jds:jde) )
      wrf_lon_max = maxval( xlong(ids:ide,jds:jde) )
      wrf_lat_min = minval( xlat(ids:ide,jds:jde) )
      wrf_lat_max = maxval( xlat(ids:ide,jds:jde) )
      write(*,*) ' '
      write(*,'('' megan2_bioemiss: model lon limits = '',1p,2g25.16)') wrf_lon_min,wrf_lon_max
      write(*,'('' megan2_bioemiss: model lat limits = '',1p,2g25.16)') wrf_lat_min,wrf_lat_max
      write(*,*) ' '
      write(*,'('' megan2_bioemiss: megan lon limits = '',1p,2g25.16)') megan_lons(1),megan_lons(nlon_megan)
      write(*,'('' megan2_bioemiss: megan lat limits = '',1p,2g25.16)') megan_lats(1),megan_lats(nlat_megan)
      write(*,*) ' '
!---------------------------------------------------------------
!     allocate memory space to store interpolation coef.
!---------------------------------------------------------------
      ierr = 0
      allocate( ax(ids:ide,jds:jde,0:1), &
                by(ids:ide,jds:jde,0:1), stat=status )
      ierr = ierr + status
      allocate( ix(ids:ide,jds:jde,0:1), &
                jy(ids:ide,jds:jde,0:1), stat=status )
      ierr = ierr + status
      if( ierr /= 0 ) then
        write(*,*) 'allocate for ax ... jy failed; error = ',ierr
      endif
!---------------------------------------------------------------------
!   set bilinear interp variables
!---------------------------------------------------------------------
lat_loop : &
      do j = jds,jde
        do i = ids,ide
!---------------------------------------------------------------------
!   longitudes
!---------------------------------------------------------------------
          wrf_lon = xlong(i,j)
          if( wrf_lon >= megan_lons(nlon_megan) .or. &
              wrf_lon < megan_lons(1) ) then
            ix(i,j,lower) = nlon_megan
          else
            do n = 2,nlon_megan
              if( wrf_lon < megan_lons(n) ) then
                ix(i,j,lower) = min( nlon_megan-1, max(n-1,1) )
                exit
              endif
            end do
          endif
          ix(i,j,upper) = mod( ix(i,j,lower),nlon_megan ) + 1
          data_dx = megan_lons(ix(i,j,upper)) - megan_lons(ix(i,j,lower))
          if( data_dx < 0. ) then
            data_dx = 360. + data_dx
          endif
          ds1 = wrf_lon - megan_lons(ix(i,j,lower))
          if( ds1 < 0. ) then
            ds1 = 360. + ds1
          endif
          ax(i,j,lower) = ds1/data_dx
          ax(i,j,upper) = 1.0 - ax(i,j,lower)
!---------------------------------------------------------------------
!   latitudes
!---------------------------------------------------------------------
          wrf_lat = xlat(i,j)
          if( wrf_lat < megan_lats(1) ) then
            jy(i,j,0:1) = -1
            by(i,j,0:1) = 0.
          elseif( wrf_lat > megan_lats(nlat_megan) ) then
            jy(i,j,0:1) = -2
            by(i,j,0:1) = 0.
          else
            do n = 1,nlat_megan
              if( wrf_lat < megan_lats(n) ) then
                exit
              endif
            end do
            jy(i,j,lower) = min( nlat_megan-1, max(n-1,1) )
            jy(i,j,upper) = jy(i,j,lower) + 1
            by(i,j,lower) = (wrf_lat - megan_lats(jy(i,j,lower))) &
                             /(megan_lats(jy(i,j,upper)) - megan_lats(jy(i,j,lower)))
            by(i,j,upper) = 1.0 - by(i,j,lower)
          endif
        end do
      end do lat_loop
!---------------------------------------------------------------------
!   form dataset index limits
!---------------------------------------------------------------------
      write(*,*) 'megan2_bioemiss: count of points <,> data min,max lat = ',count(ix(:,:,0) == nlon_megan )
      xndx_megan(1) = minval( ix(:,:,lower) )
      xndx_megan(2) = maxval( ix(:,:,upper) )
      write(*,*) 'xndx_megan = ',xndx_megan(:)
      write(*,*) 'megan2_bioemiss: count of points < data min lat = ',count(jy(:,:,0) == -1)
      write(*,*) 'megan2_bioemiss: count of points > data max lat = ',count(jy(:,:,0) == -2)
      yndx_megan(1) = minval( jy(:,:,lower),mask=jy(:,:,lower)>0 )
      yndx_megan(2) = maxval( jy(:,:,upper),mask=jy(:,:,upper)>0 )
      write(*,*) 'yndx_megan = ',yndx_megan(:)

      if( debug ) then
      write(*,*) ' '
      write(*,*) 'megan2_bioemiss: bilinear interp diagnostics'
      write(*,*) 'megan2_bioemiss: ix'
      write(*,*) ix(ids,jds,:)
      write(*,*) 'megan2_bioemiss: ax'
      write(*,*) ax(ids,jds,:)
      write(*,*) 'megan2_bioemiss: megan lons'
      write(*,*) megan_lons(ix(ids,jds,0)),megan_lons(ix(ids,jds,1))
      write(*,*) 'megan2_bioemiss: wrf lon = ',xlong(ids,jds)
      write(*,*) 'megan2_bioemiss: jy'
      write(*,*) jy(ids,jds,:)
      write(*,*) 'megan2_bioemiss: by'
      write(*,*) by(ids,jds,:)
      write(*,*) 'megan2_bioemiss: megan lats'
      write(*,*) megan_lats(jy(ids,jds,0)),megan_lats(jy(ids,jds,1))
      write(*,*) 'megan2_bioemiss: wrf lat = ',xlat(ids,jds)
      write(*,*) ' '
      do j = jds,jde
        do i = ids,ide
          if( ix(i,j,lower) == nlon_megan ) then
      write(*,*) 'megan2_bioemiss: bilinear interp diagnostics'
      write(*,*) 'megan2_bioemiss: ix'
      write(*,*) ix(i,j,:)
      write(*,*) 'megan2_bioemiss: ax'
      write(*,*) ax(i,j,:)
      write(*,*) 'megan2_bioemiss: megan lons'
      write(*,*) megan_lons(ix(i,j,0)),megan_lons(ix(i,j,1))
      write(*,*) 'megan2_bioemiss: wrf lon = ',xlong(i,j)
      write(*,*) 'megan2_bioemiss: jy'
      write(*,*) jy(i,j,:)
      write(*,*) 'megan2_bioemiss: by'
      write(*,*) by(i,j,:)
      write(*,*) 'megan2_bioemiss: megan lats'
      write(*,*) megan_lats(jy(i,j,0)),megan_lats(jy(i,j,1))
      write(*,*) 'megan2_bioemiss: wrf lat = ',xlat(i,j)
            stop 'diagnostics'
          endif
        end do
      end do
      stop 'diagnostics'
      endif

!---------------------------------------------------------------------
!   allocate and read dataset variable
!---------------------------------------------------------------------
       if( allocated( tmp3 ) ) then
          deallocate( tmp3 )
       endif
       allocate( tmp3(xndx_megan(1):xndx_megan(2),yndx_megan(1):yndx_megan(2),mnth_s:mnth_e),stat=ierr )
                 
       if( ierr /= 0 ) then
         write(message,*) 'Failed to allocate tmp3 for lai; error = ',ierr
         stop 'bio_emiss abort'
       endif
       message = 'Failed to get variable id'
       call handle_ncerr( nf_inq_varid( ncid, trim(varname), varid ), message )
       message = 'Failed to read variable'
       if( mnth_s == mnth_e ) then
          call handle_ncerr( nf_get_vara_real( ncid, varid, &
                                               (/ xndx_megan(1),yndx_megan(1) /), &
                                               (/ xndx_megan(2)-xndx_megan(1)+1,yndx_megan(2)-yndx_megan(1)+1 /), &
                                               tmp3 ), message )
       else
          call handle_ncerr( nf_get_vara_real( ncid, varid, &
                                               (/ xndx_megan(1),yndx_megan(1),mnth_s /), &
                                               (/ xndx_megan(2)-xndx_megan(1)+1, &
                                                  yndx_megan(2)-yndx_megan(1)+1,mnth_e-mnth_s+1 /), &
                                               tmp3 ), message )
       endif

       write(*,*)  'dataset size = ',size(tmp3)
       write(*,*)  'dataset min,max values = ',minval(tmp3(:,:,:)),maxval(tmp3(:,:,:))
       write(*,*)  'dataset missing value count = ',count(tmp3(:,:,:) == missing_value )
       write(*,*)  '% valid data = ',100.* real( count(tmp3(:,:,:) /= missing_value ) ) /real(size(tmp3))

!---------------------------------------------------------------------
!   replace missing values with zero
!---------------------------------------------------------------------
       where( tmp3(:,:,:) == missing_value )
          tmp3(:,:,:) = 0.
       endwhere
!---------------------------------------------------------------------
!   set wrf bioemission variable
!---------------------------------------------------------------------
mnth_loop : &
       do mnth = mnth_s,mnth_e
         do j = jds,jde
           do i = ids,ide
             jl = jy(i,j,0)
             if( jl > 0 ) then
               il = ix(i,j,0)
               iu = ix(i,j,1)
               ju = jy(i,j,1)
               wrk_sum  = tmp3(il,jl,mnth)*ax(i,j,upper)*by(i,j,upper) &
                        + tmp3(il,ju,mnth)*ax(i,j,upper)*by(i,j,lower) &
                        + tmp3(iu,jl,mnth)*ax(i,j,lower)*by(i,j,upper) &
                        + tmp3(iu,ju,mnth)*ax(i,j,lower)*by(i,j,lower)
             else
               wrk_sum  = 0.
             endif
             if( varname(:3) == 'LAI' ) then
               mlai(i,j,megan_month) = scale_factor * wrk_sum
             endif

           end do
         end do
       end do mnth_loop
    endif is_area_map

!---------------------------------------------------------------------
!   exit
!---------------------------------------------------------------------
    if( has_area_map ) then
      deallocate( xedge_megan, yedge_megan )
    else
       deallocate( ix, jy, ax, by, tmp3 )
    endif

    write(*,*) ' Finished megan2 bio emiss dataset ',trim(inpname)
!---------------------------------------------------------------------
!   close megan dataset file
!---------------------------------------------------------------------
    message = 'Failed to close ' // trim(inpname)
    call handle_ncerr( nf_close( ncid ), message )       

    end subroutine megan2_bioemiss

   subroutine handle_ncerr( ret, mes )
!---------------------------------------------------------------------
!	... netcdf error handling routine
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!	... dummy arguments
!---------------------------------------------------------------------
   integer, intent(in) :: ret
   character(len=*), intent(in) :: mes

   if( ret /= nf_noerr ) then
      write(*,*) nf_strerror( ret )
      stop 'netcdf error'
   endif

   end subroutine handle_ncerr

      subroutine write_lai3
!---------------------------------------------------------------------
!       ... write the lai  file
!---------------------------------------------------------------------
!       ... local variables
!---------------------------------------------------------------------
      integer :: cell_id,istep,it,im
      integer :: ilat,ilon
      integer :: ilat_mcip, ilon_mcip !shc
      integer :: dsrad,dtemp
      real :: xlai(ide,jde,46)

      outpname = trim(out_dir)//'/'//'LAI3.csv' !shc added out_dir
      open(unit=99,file=outpname,status='unknown')
      write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LAI01,LAI02,LAI03,LAI04,LAI05,LAI06,LAI07,LAI08,LAI09,LAI10, &
                      LAI11,LAI12,LAI13,LAI14,LAI15,LAI16,LAI17,LAI18,LAI19,LAI20, &
                      LAI21,LAI22,LAI23,LAI24,LAI25,LAI26,LAI27,LAI28,LAI29,LAI30, &
                      LAI31,LAI32,LAI33,LAI34,LAI35,LAI36,LAI37,LAI38,LAI39,LAI40, &
                      LAI41,LAI42,LAI43,LAI44,LAI45,LAI46"
      

!-----if you have only monthly LAI, you will need to interpolate monthly LAI to 8-day LAI xlai using mlai      
!      istep=0
!      do im=1,11
!        
!         istep=istep+1
!         write(*,*)'Interpolating LAI',im,istep
!         xlai(:,:,istep)=mlai(:,:,im)-(mlai(:,:,im+1)-mlai(:,:,im))/6.0
 !        istep=istep+1
 !        write(*,*)'Interpolating LAI',im,istep
 !        xlai(:,:,istep)=mlai(:,:,im)+(mlai(:,:,im+1)-mlai(:,:,im))/3.0 
  !       istep=istep+1
 !        write(*,*)'Interpolating LAI',im,istep
 !        xlai(:,:,istep)=mlai(:,:,im)+(mlai(:,:,im+1)-mlai(:,:,im))*2.0/3.0 
  !       istep=istep+1
  !       write(*,*)'Interpolating LAI',im,istep
  !       xlai(:,:,istep)=mlai(:,:,im)+(mlai(:,:,im+1)-mlai(:,:,im))*3.0/3.0
 !
 !      end do         
 !        istep=istep+1
 !        write(*,*)'Interpolating LAI',im,istep
 !        xlai(:,:,istep)=mlai(:,:,im)+(mlai(:,:,im+1)-mlai(:,:,im))/6.0
 !        istep=istep+1
  !       write(*,*)'Interpolating LAI',im,istep
  !       xlai(:,:,istep)=mlai(:,:,im+1)+(mlai(:,:,im+1)-mlai(:,:,im))/6.0

     xlai(:,:,:)=mlai(:,:,:)

      do ilat=1,jde-jds+1
        do ilon=1,ide-ids+1
        do im=1,46 
      if(xlai(ilon,ilat,im).gt.20.or.xlai(ilon,ilat,im).lt.0.0)xlai(ilon,ilat,im)=0.0
        end do
      end do
      end do
!-----write new xlai to output


!shc     cell_id=0
!shc      do ilat=2,jde-jds
!shc       do ilon=2,ide-ids
!shc      cell_id=cell_id+1
!shc      write (99, '(3(I,","), 46(F10.4,","))')cell_id,ilon-1,ilat-1,xlai(ilon,ilat,:)
!shc       end do
!shc      end do
!shc     print*,'jiang finished LAI'

     !shc below is modified from the above to trim output to the MCIP/CMAQ domain
      cell_id=0
      do ilat=2,jde-jds
         do ilon=2,ide-ids
            if ( (ilon.gt.x0) .and. (ilon.le.(ncolsin+x0)) .and.   &
                 (ilat.gt.y0) .and. (ilat.le.(nrowsin+y0))  &
                 ) then
               cell_id=cell_id+1
               ilon_mcip = (ilon - x0 + 1) - 1
               ilat_mcip = (ilat - y0 + 1) - 1
               write (99, '(3(I,","), 48(F10.4,","))')   &
                    cell_id,ilon_mcip,ilat_mcip,XLAT(ilon,ilat),XLONG(ilon,ilat), &
                     xlai(ilon,ilat,:)
            endif
         end do
      end do
      print*,'finished writing LAI to .csv file for CMAQ'


     close(99)
      end subroutine write_lai3

end program map_megan2_emissions
