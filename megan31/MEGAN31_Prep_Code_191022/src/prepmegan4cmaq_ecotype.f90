!--------------------------------------------------------------------
! Based of prepmegan4cmaq_pft.f90 to process ecotype for MEGANv3
!
! Serena H. Chung                                                                                                                                                                                      
! 2014-05-02                                                                                                                                                                                           
! Updated by Jeremiah Johnson 03/30/17
! 
!--------------------------------------------------------------------

   program map_megan3_ecotype

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
   integer :: ietype
   integer :: mxetype   = 100
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
!----------------------------------------------------
!----allocate ecotype arrays-----------------------------
!----------------------------------------------------
   integer, allocatable :: ecotypeid(:,:,:)
   real, allocatable :: ecotypefrac(:,:,:)

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
   CHARACTER (LEN=120)  :: outpname
   CHARACTER (LEN=300)  :: wrffile
   CHARACTER (LEN=200)  :: megan_dir
   CHARACTER (LEN=200)  :: out_dir !shc
   CHARACTER (LEN=19)  :: Times(1)
   CHARACTER (LEN=4)   :: num  !used to calcualte LAI date

   logical :: new_grid

   namelist /control/ start_lai_mnth, end_lai_mnth, &
                      wrffile, megan_dir, out_dir !shc added out_dir
   namelist /windowdefs/ x0, y0, ncolsin, nrowsin !shc

!---------------------------------------------------------------------
!	... include files
!---------------------------------------------------------------------
   include 'netcdf.inc'

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
!   if( domains < 1 ) then
!     write(*,*) 'convert_emissions: domains must be >= 1'
!     stop 'bio_emiss abort'
!   endif

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


!-----read ecotype----------------------------------------------
      scale_factor  = 1.
      inpname       = 'ecotype30s.nc'
      varname       = 'ecotype'
      missing_value = -32768.
      CALL  megan2_bioemiss

!-----------------------------------------------------------
!     ... write to output file(s)
!-----------------------------------------------------------
!write to grid_ecotype.csv file---------------------
      CALL write_ecotype(ecotypeid,ecotypefrac)

!-----------------------------------------------------------
!     cleanup domain variables
!-----------------------------------------------------------
!----- include more emissions factors
      deallocate( ecotypeid )
      deallocate( ecotypefrac )
      if( allocated( xlong ) ) then
        deallocate( xlong )
      endif
      if( allocated( xlat ) ) then
        deallocate( xlat )
      endif
      do n = 1,grid_cnt
        if( associated( grid_specs(n)%lon ) ) then
          deallocate( grid_specs(n)%lon )
        endif
        if( associated( grid_specs(n)%lat ) ) then
          deallocate( grid_specs(n)%lat )
        endif
        if( associated( grid_specs(n)%model_area_type ) ) then
          do j = 1,jde
            do i = 1,ide
              if( associated( grid_specs(n)%model_area_type(i,j)%dcell_lon_ndx ) ) then
                deallocate( grid_specs(n)%model_area_type(i,j)%dcell_lon_ndx )
              endif
              if( associated( grid_specs(n)%model_area_type(i,j)%dcell_lat_ndx ) ) then
                deallocate( grid_specs(n)%model_area_type(i,j)%dcell_lat_ndx )
              endif
              if( associated( grid_specs(n)%model_area_type(i,j)%wght ) ) then
                deallocate( grid_specs(n)%model_area_type(i,j)%wght )
              endif
            end do
          end do
          deallocate( grid_specs(n)%model_area_type )
        endif
      end do
!++sw
      grid_cnt = 0
!--sw
!   end do domain_loop

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

!   message = 'Failed to get Times id'
!   call handle_ncerr( nf_inq_varid( ncid, 'Times', varid ), message )
!   message = 'Failed to read Times'
!   call handle_ncerr( nf_get_var_text( ncid, varid, Times ), message )

   write(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!   write(*,*) 'wrf_file: time = ',trim(Times(1))
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

!---------------------------

   allocate( ecotypeid(ide,jde,mxetype),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate ecotypeid; error = ',astat
     stop 'allocate failed'
   endif
   allocate( ecotypefrac(ide,jde,mxetype),stat=astat ) 
   if( astat /= 0 ) then
     write(*,*) 'wrf_file: failed to allocate ecotypefrac; error = ',astat
     stop 'allocate failed'
   endif

   end subroutine wrf_file

   subroutine megan2_bioemiss
!---------------------------------------------------------------------
!   map megan dataset to wrf grid
!---------------------------------------------------------------------

   use area_mapper, only : discrete_frac
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
    real, allocatable :: wrk_data(:,:,:,:)
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
    call handle_ncerr( nf_inq_dimid( ncid, 'x', dimid ), message )
    message = 'megan2_bioemiss: Failed to get lon dimension'
    call handle_ncerr( nf_inq_dimlen( ncid, dimid, nlon_megan ), message )
    message = 'megan2_bioemiss: Failed to get lat dimension id'
    call handle_ncerr( nf_inq_dimid( ncid, 'y', dimid ), message )
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
    call handle_ncerr( nf_inq_varid( ncid, 'x', varid ), message )
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
    call handle_ncerr( nf_inq_varid( ncid, 'y', varid ), message )
    message = 'megan2_bioemiss: Failed to read lat variable'
    call handle_ncerr( nf_get_var_real( ncid, varid, megan_lats ), message )

!---------------------------------------------------------------------
!   Assume discrete fraction methodology for ecotype
!---------------------------------------------------------------------

    write(*,*) 'megan2_bioemiss: data_dx,dx = ',data_dx,dx

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
     grid_ndx = grid_cnt
     write(*,*) 'megan2_bioemiss: file ' // trim(inpname),' has a new grid'
   endif

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
   allocate( wrk_data(ide,jde,mxetype,2),stat=status )
   if( status /= 0 ) then
      write(*,*) 'megan2_bioemiss: allocate for wrk_data failed; error = ',ierr
      stop 'allocation error'
   endif
!---------------------------------------------------------------------
!   discrete fraction calculation
!---------------------------------------------------------------------
   call discrete_frac( xedge_megan, yedge_megan, nlon_megan, nlat_megan, int(missing_value,2), &
                        wrk_data, ncid, varname, grid_ndx, new_grid )
   ecotypeid(:,:,:)   = wrk_data(:,:,:,1)
   ecotypefrac(:,:,:) = wrk_data(:,:,:,2)

   deallocate( wrk_data )

!---------------------------------------------------------------------
!   exit
!---------------------------------------------------------------------
   deallocate( xedge_megan, yedge_megan )

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

   subroutine write_ecotype(ecotypeid,ecotypefrac)
!---------------------------------------------------------------------
!       ... write the ecotype  file
!---------------------------------------------------------------------
!       ... local variables
!---------------------------------------------------------------------
   integer :: cell_id,istep,it,im
   integer :: ilat,ilon
   integer :: ilat_mcip, ilon_mcip !shc
   integer :: dsrad,dtemp
   integer :: ecotypeid(:,:,:)
   real    :: ecotypefrac(:,:,:)

   outpname = trim(out_dir)//'/'//'grid_ecotype.csv' !shc added out_dir
   open(unit=99,file=outpname,status='unknown')
   write(99,'(a)')"gridID,EcotypeID,EcotypeFrac"

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
              do ietype=1,mxetype
                if (ecotypefrac(ilon,ilat,ietype) > 0.) then
                    write (99, '( 2(I0,","),F0.4)')cell_id,&
                         ecotypeid(ilon,ilat,ietype),ecotypefrac(ilon,ilat,ietype)
                else 
                  exit
                endif
              end do
         endif
      end do
   end do

   print*,'megan_bio_emiss is done'
   close(99)

   end subroutine write_ecotype

end program map_megan3_ecotype
