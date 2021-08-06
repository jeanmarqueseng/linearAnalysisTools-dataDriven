
! Tulio R.
! qui jun 11 22:44:58 -03 2020
! Jean H. M. Ribeiro
! qua set 22 13:52:00 -03 2020

!****************************************************************************
! 
!     LIST OF SUBROUTINES :   
! 
! 
! 
!############################################################################
module mod_CGNS

  !character(len=16)  :: working_var

  character(len=132) :: CGNS_filename
  character(len=132) :: CGNS_gridname
  character(len=132) :: CGNS_solnname

  integer(kind=4)    :: index_grid, index_soln, index_file, index_mean
  integer(kind=4)    :: index_base, index_zone, index_node, index_flow, index_mode
  integer(kind=4)    :: index_field, index_coord, ifile
  
  !integer(kind=4)    :: index_hole
 
  integer(kind=4)    :: ier
  integer(kind=4)    :: location
  integer(kind=4)    :: cell_dim, phys_dim

!  integer(kind=4), allocatable :: isize(:,:)
  integer(kind=4), allocatable :: isize(:)

  character(len=16)  :: basename
  character(len=32)  :: zonename
  character(len=32)  :: solnname
  character(len=32)  :: fieldname
  character(len=5)   :: citer
  character(len=2)   :: cbin

  character(len=5)   :: cmode
  
  logical :: CGNS_file_exists

end module
!############################################################################

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine open_file_read(path, ifile, nbases)

  use cgns
  implicit none
  character(len=200), intent(in)  :: path
  integer,            intent(out) :: ifile, nbases
  integer :: ier

  call cg_open_f(trim(path),CG_MODE_READ,ifile,ier)
  if (ier.ne.CG_OK) call cg_error_exit_f

  call cg_nbases_f(ifile, nbases, ier)
  if (ier.ne.CG_OK) call cg_error_exit_f
  
  return
    
end subroutine open_file_read
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine close_file(ifile)
  
  use cgns
  implicit none

  integer, intent(in) :: ifile
  integer :: ier

  call cg_close_f(ifile,ier)
  if(ier.ne.CG_OK) call cg_error_exit_f
  
  return
    
end subroutine close_file
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!------------------------------------------------------------------------------------
subroutine create_file_CGNS(char_filename,aux_dim)

  use cgns
  use mod_CGNS
  implicit none
  
  character(len=*), intent(in) :: char_filename
  character(len=*), intent(in) :: aux_dim
  
  logical :: check
  
  CGNS_solnname = char_filename
  
  check = .false.
  
  inquire(file=CGNS_solnname,exist=CGNS_file_exists)
  if (CGNS_file_exists .eqv. .false.) then

    !open solution file
    call cg_open_f(trim(CGNS_solnname),CG_MODE_WRITE,index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    !create base in the recently opened solution file
    if (aux_dim .eq. '2D' .or. aux_dim .eq. '2d') then
      call cg_base_write_f(index_file,'Base',2,2,index_base,ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
      check = .true.
    endif

    if (aux_dim .eq. '3D' .or. aux_dim .eq. '3d') then
      call cg_base_write_f(index_file,'Base',3,3,index_base,ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
      check = .true.
    endif
    
    if (check .eqv. .false.) stop ' Tem cagada na chamada do create_file_CGNS ... '

    !close the file
    call cg_close_f(index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    write(*,'(A)') trim(char_filename) // ' already exists!'

  endif

  return

end subroutine create_file_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------

!      |||||||||||||||     |||||||||||||       ||||||||||||  |||||||||||||||
!     |||||||||||||||||   |||||||||||||||      ||||||||||||  |||||||||||||||||
!     |||||||             ||||||    ||||||        ||||||     ||||||     |||||||
!     ||||||              ||||||     ||||||       ||||||     ||||||       ||||||
!     ||||||   ||||||||   ||||||    ||||||        ||||||     ||||||       ||||||
!     ||||||   |||||||||  |||||||||||||||         ||||||     ||||||       ||||||
!     ||||||      ||||||  |||||||||||||||         ||||||     ||||||       ||||||
!     |||||||     ||||||  ||||||     ||||||       ||||||     ||||||     |||||||
!     ||||||||||||||||||  ||||||      ||||||   ||||||||||||  |||||||||||||||||
!      ||||||||||||||||   ||||||       ||||||  ||||||||||||  |||||||||||||||

!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_unstructured_grid_3D_CGNS(zone_number, nnode,nelem,cell_shape,xcoord,ycoord,zcoord,conn,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: nnode, nelem, cell_shape
  real(kind=8), dimension(:), intent(in)     :: xcoord
  real(kind=8), dimension(:), intent(in)     :: ycoord
  real(kind=8), dimension(:), intent(in)     :: zcoord
  integer(kind=4), dimension(:,:), intent(in)  :: conn
  character(len=*), intent(in) :: char_filename
  integer(kind=4)  :: m, idummy
  
  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  !allocate(isize(1,cell_dim))
  allocate(isize(3))

  m = zone_number

    write(zonename,'(a4,i4.4)') 'Zone', m

    call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
      if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

    isize(1) = nnode
    isize(2) = nelem
    isize(3) = 0
    
    call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Unstructured,index_zone,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateX', xcoord(:), index_coord, ier )
      if (ier .ne. CG_OK) call cg_error_exit_f
      print*, 'X ok...'

    call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateY', ycoord(:), index_coord, ier )
      if (ier .ne. CG_OK) call cg_error_exit_f
      print*, 'Y ok...'
      
    call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateZ', zcoord(:), index_coord, ier )
      if (ier .ne. CG_OK) call cg_error_exit_f
      print*, 'Z ok...'

    call cg_section_write_f(index_grid, index_base, index_zone, 'Elems', HEXA_8,  1, nelem, 0, conn, idummy, ier)
      print*, ier
      if (ier .ne. CG_OK) call cg_error_exit_f
      print*, 'conn ok...'

  else
  
    write(*,'(A,i0,A)') 'The grid coordinates for zone ', m, ' already exist! Skipping ...'

  endif

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)

  return

end subroutine write_unstructured_grid_3D_CGNS
!------------------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine zone_size_read(ifile, ibase, izone, idim, isize, nx, ny, nz)

  use cgns
  implicit none

  integer, intent(in)  :: ifile, ibase, izone, idim
  integer, dimension(idim,3), intent(out) :: isize
  integer,                    intent(out) :: nx, ny, nz
  character(len=200) :: zonename
  integer :: ier

  call cg_zone_read_f(ifile,ibase,izone,zonename,isize,ier)
  if(ier.ne.CG_OK) call cg_error_exit_f

  ! Allocate grid stuff and read grid coordinates
  if(idim==2)then
      nx = isize(1,1)
      ny = isize(2,1)  ! 2D problems
      nz = 1
  else
      nx = isize(1,1)
      ny = isize(2,1)  ! 3D problems
      nz = isize(3,1)
  endif
  
  return
    
end subroutine zone_size_read
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!------------------------------------------------------------------------------------

!      |||||||||||||||     |||||||||||||      |||||         ||||||     ||||
!     |||||||||||||||||   ||||||||||||||||    |||||         |||||||    ||||  
!     |||||||             ||||||    ||||||    |||||         ||||||||   ||||   
!     ||||||              |||||      |||||    |||||         |||||||||  ||||    
!     |||||||||||||       ||||        ||||    |||||         |||| ||||| ||||          
!         |||||||||||     ||||        ||||    |||||         ||||  |||||||||
!            ||||||||||   |||||       ||||    |||||         ||||   |||||||| 
!     ||        ||||||||  ||||||     |||||    |||||         ||||    |||||||
!     ||||||||||||||||||  ||||||||||||||||    ||||||||||||  ||||     ||||||  
!      ||||||||||||||||    ||||||||||||||     ||||||||||||  ||||      |||||

!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_unstructured_soln_3D_CGNS(zone_number,nnode,nelem,cell_shape,qdata,var_location,char_filename)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: nnode, nelem, cell_shape
  real(kind=8), dimension(:,:), intent(in)     :: qdata
  character(len=*), intent(in) :: char_filename
  character(len=*), intent(in) :: var_location
  integer(kind=4)  :: m, idummy
  
  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_WRITE,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_base_write_f(index_grid,'Base',3,3,index_base,ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  !allocate(isize(1,cell_dim))
  allocate(isize(3))

  m = zone_number

  write(zonename,'(a4,i4.4)') 'Zone', m

  write(*,'(A,i4.4,A)') '  -> Creating zone ', m, ' in ' // trim(char_filename)

  isize(1) = nnode
  isize(2) = nelem
  isize(3) = 0
    
  call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Unstructured,index_zone,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  
  if ( trim(var_location) .eq. "Vertex") then
    call cg_sol_write_f(index_grid,index_base,index_zone,'FlowSolution',Vertex,index_flow,ier)
      print*, ier
      if (ier .ne. CG_OK) call cg_error_exit_f
      print*, 'vertex soln ok...'
  elseif ( trim(var_location) .eq. "CellCenter") then
    call cg_sol_write_f(index_grid,index_base,index_zone,'FlowSolution',CellCenter,index_flow,ier)
      print*, ier
      if (ier .ne. CG_OK) call cg_error_exit_f
      print*, 'cellcentered soln ok...'
  else
    print*, 'Only Vertex and CellCenter soln files are possible'
  endif

  call cg_field_write_f(index_grid,index_base,index_zone,index_flow, RealDouble,'Density',qdata(:,1),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  call cg_field_write_f(index_grid,index_base,index_zone,index_flow, RealDouble,'Velocity-X',qdata(:,2),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  call cg_field_write_f(index_grid,index_base,index_zone,index_flow, RealDouble,'Velocity-Y',qdata(:,3),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  call cg_field_write_f(index_grid,index_base,index_zone,index_flow, RealDouble,'Velocity-Z',qdata(:,4),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  call cg_field_write_f(index_grid,index_base,index_zone,index_flow, RealDouble,'Pressure',qdata(:,5),index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

!  DEALLOCATE(isize)

  return

end subroutine write_unstructured_soln_3D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_link_CGNS(zone_number,char_solnname,char_pathname)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  character(len=*), intent(in) :: char_pathname, char_solnname
  integer(kind=4) :: m
  character(len=132) linkpath, linkpathg, linkpathe

  write(*,'(A)') '  -> Linking file ' // trim(char_solnname) // ' to ' // trim(char_pathname)
  inquire(file=trim(char_solnname),exist=CGNS_file_exists)
    if (CGNS_file_exists .eqv. .false.) then
      print*, trim(char_solnname)
      stop ' O arquivo CGNS não existe!!!'
    endif
  
  !open solution file
  call cg_open_f(trim(char_solnname), CG_MODE_MODIFY, ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    
    write(zonename,'(a4,i4.4)') 'Zone', zone_number

    linkpath = '/Base/' // trim(zonename)

    print*, trim(linkpath)
    call cg_gopath_f(ifile, trim(linkpath), ier)
      if (ier .ne. CG_OK) call cg_error_print_f

    linkpathg = trim(linkpath) // '/GridCoordinates'

    call cg_link_write_f('GridCoordinates', trim(char_pathname), trim(linkpathg), ier)
      if (ier .ne. CG_OK) call cg_error_print_f
      
    linkpathe = trim(linkpath) // '/Elems'

    call cg_link_write_f('Elems', trim(char_pathname), trim(linkpathe), ier)
      if (ier .ne. CG_OK) call cg_error_print_f

    write(*,'(A,i0)') 'Linking zone ', zone_number

  call cg_close_f(ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  return

end subroutine write_link_CGNS
