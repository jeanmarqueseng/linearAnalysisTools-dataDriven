
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
subroutine create_file_CGNS(char_filename,cell_dim)

  use cgns
  implicit none
  
  character(len=*), intent(in) :: char_filename
  character(len=132) :: CGNS_solnname
  integer(kind=4)    :: index_file, index_base, cell_dim, ier
  
  CGNS_solnname = char_filename
  
  !open solution file
  call cg_open_f(trim(CGNS_solnname),CG_MODE_WRITE,index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_base_write_f(index_file,'Base',cell_dim,3,index_base,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  !close the file
  call cg_close_f(index_file,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

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
subroutine write_unstructured_grid_2D_CGNS(nnode,nelem,xcoord,ycoord,zcoord,conn,char_filename)

  use cgns
  implicit none
  
  integer(kind=4), intent(in)  :: nnode, nelem
  real(kind=8), dimension(:), intent(in)     :: xcoord
  real(kind=8), dimension(:), intent(in)     :: ycoord
  real(kind=8), dimension(:), intent(in)     :: zcoord
  integer(kind=4), dimension(:,:), intent(in)  :: conn
  character(len=*), intent(in) :: char_filename

  character(len=132) :: CGNS_filename
  character(len=16)  :: basename
  character(len=32)  :: zonename
  
  integer(kind=4), dimension(:), allocatable :: isize
  integer(kind=4)    :: idummy, cell_dim, phys_dim
  integer(kind=4)    :: index_grid, index_base, index_zone, index_coord, ier

  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(3))

  zonename='Zone0001'

  call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  isize(1) = nnode
  isize(2) = nelem
  isize(3) = 0
    
  call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Unstructured,index_zone,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateX', xcoord(:), index_coord, ier )
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, '...X coordinate for grid is set to CGNS grid file...'

  call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateY', ycoord(:), index_coord, ier )
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, '...Y coordinate for grid is set to CGNS grid file...'
      
  call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateZ', zcoord(:), index_coord, ier )
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, '...Z coordinate for grid is set to CGNS grid file...'

  call cg_section_write_f(index_grid, index_base, index_zone, 'Elems', QUAD_4,  1, nelem, 0, conn, idummy, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, 'Element connectivity table is set to CGNS grid file'

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_unstructured_grid_2D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_unstructured_grid_3D_CGNS(nnode,nelem,xcoord,ycoord,zcoord,conn,char_filename)

  use cgns
  implicit none
  
  integer(kind=4), intent(in)  :: nnode, nelem
  real(kind=8), dimension(:), intent(in)     :: xcoord
  real(kind=8), dimension(:), intent(in)     :: ycoord
  real(kind=8), dimension(:), intent(in)     :: zcoord
  integer(kind=4), dimension(:,:), intent(in)  :: conn
  character(len=*), intent(in) :: char_filename

  character(len=132) :: CGNS_filename
  character(len=16)  :: basename
  character(len=32)  :: zonename
  
  integer(kind=4), dimension(:), allocatable :: isize
  integer(kind=4)    :: idummy, cell_dim, phys_dim
  integer(kind=4)    :: index_grid, index_base, index_zone, index_coord, ier

  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  index_base = 1
  call cg_base_read_f(index_grid, index_base, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(3))

  zonename='Zone0001'

  call cg_gopath_f(index_grid, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  isize(1) = nnode
  isize(2) = nelem
  isize(3) = 0
    
  call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Unstructured,index_zone,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateX', xcoord(:), index_coord, ier )
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, '...X coordinate for grid is set to CGNS grid file...'

  call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateY', ycoord(:), index_coord, ier )
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, '...Y coordinate for grid is set to CGNS grid file...'
      
  call cg_coord_write_f(index_grid, index_base, index_zone, RealDouble, 'CoordinateZ', zcoord(:), index_coord, ier )
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, '...Z coordinate for grid is set to CGNS grid file...'

  call cg_section_write_f(index_grid, index_base, index_zone, 'Elems', HEXA_8,  1, nelem, 0, conn, idummy, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    !print*, 'Element connectivity table is set to CGNS grid file'

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_unstructured_grid_3D_CGNS
!------------------------------------------------------------------------------------

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine zone_size_read(index_file, index_base, index_zone, nnodes, nelem, idummy)

  use cgns
  implicit none

  integer, intent(in)   :: index_file, index_base, index_zone
  integer, intent(out)  :: nnodes, nelem, idummy
  integer, dimension(:), allocatable :: isize 
  character(len=32) :: zonename
  integer :: ier
  
  allocate(isize(3))

  call cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
  if(ier.ne.CG_OK) call cg_error_exit_f

  ! Allocate grid stuff 
  nnodes = isize(1)
  nelem  = isize(2)  
  idummy = isize(3)
  
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
subroutine write_unstructured_solndata_3D_CGNS(nnode,nelem,var_location,char_filename)

  use cgns
  implicit none
  
  integer(kind=4), intent(in)  :: nnode, nelem
  character(len=*), intent(in) :: char_filename
  character(len=*), intent(in) :: var_location
  
  character(len=132) :: CGNS_filename
  character(len=32)  :: zonename

  integer(kind=4)    :: index_grid, index_base, index_zone
  integer(kind=4)    :: index_flow, index_field, ier
  integer(kind=4), allocatable :: isize(:)

  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_WRITE,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_base_write_f(index_grid,'Base',3,3,index_base,ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  !allocate(isize(1,cell_dim))
  allocate(isize(3))

  zonename='Zone0001'

  isize(1) = nnode
  isize(2) = nelem
  isize(3) = 0
    
  call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Unstructured,index_zone,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  
  if ( trim(var_location) .eq. "Vertex") then
    call cg_sol_write_f(index_grid,index_base,index_zone,'FlowSolution',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      !print*, '...CGNS soln saved to Nodes...'
  elseif ( trim(var_location) .eq. "CellCenter") then
    call cg_sol_write_f(index_grid,index_base,index_zone,'FlowSolution',CellCenter,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      !print*, '...CGNS soln saved to Cell Center...'
  else
    !print*, 'Only Vertex and CellCenter soln files are possible'
  endif

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_unstructured_solndata_3D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_unstructured_fielddata_3D_CGNS(qdata,ivar,char_filename)

  use cgns
  implicit none
  
  integer(kind=4), intent(in)  :: ivar
  real(kind=8), dimension(:), intent(in)     :: qdata
  character(len=*), intent(in) :: char_filename
  
  character(len=132) :: CGNS_filename
  character(len=32)  :: varname

  integer(kind=4)    :: index_grid, index_field, ier
  
  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  if (ivar .eq. 0) then
    varname = 'Density'
  elseif (ivar .eq. 1) then
    varname = 'Velocity-X'
  elseif (ivar .eq. 2) then
    varname = 'Velocity-Y'
  elseif (ivar .eq. 3) then
    varname = 'Velocity-Z'
  elseif (ivar .eq. 4) then
    varname = 'Vorticity-X'
  elseif (ivar .eq. 5) then
    varname = 'Vorticity-Y'
  elseif (ivar .eq. 6) then
    varname = 'Vorticity-Z'
  elseif (ivar .eq. 7) then
    varname = 'Pressure'
  else
    varname = 'Auxiliary'
  endif

  call cg_field_write_f(index_grid,1,1,1, RealDouble,varname,qdata,index_field,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  return

end subroutine write_unstructured_fielddata_3D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_unstructured_soln_3D_CGNS(nnode,nelem,qdata,var_location,char_filename)

  use cgns
  implicit none
  
  integer(kind=4), intent(in)  :: nnode, nelem
  real(kind=8), dimension(:,:), intent(in)     :: qdata
  character(len=*), intent(in) :: char_filename
  character(len=*), intent(in) :: var_location
  
  character(len=132) :: CGNS_filename
  character(len=32)  :: zonename

  integer(kind=4)    :: index_grid, index_base, index_zone
  integer(kind=4)    :: index_flow, index_field, ier
  integer(kind=4), allocatable :: isize(:)

  !... open CGNS file
  CGNS_filename = trim(char_filename)
  call cg_open_f(trim(CGNS_filename),CG_MODE_WRITE,index_grid,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_base_write_f(index_grid,'Base',3,3,index_base,ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  !allocate(isize(1,cell_dim))
  allocate(isize(3))

  zonename='Zone0001'

  isize(1) = nnode
  isize(2) = nelem
  isize(3) = 0
    
  call cg_zone_write_f(index_grid,index_base,trim(zonename),isize,Unstructured,index_zone,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
  
  if ( trim(var_location) .eq. "Vertex") then
    call cg_sol_write_f(index_grid,index_base,index_zone,'FlowSolution',Vertex,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      !print*, '...CGNS soln saved to Nodes...'
  elseif ( trim(var_location) .eq. "CellCenter") then
    call cg_sol_write_f(index_grid,index_base,index_zone,'FlowSolution',CellCenter,index_flow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
      !print*, '...CGNS soln saved to Cell Center...'
  else
    !print*, 'Only Vertex and CellCenter soln files are possible'
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

  deallocate(isize)

  return

end subroutine write_unstructured_soln_3D_CGNS
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine write_link_CGNS(char_solnname,char_pathname)

  use cgns
  implicit none

  character(len=*), intent(in) :: char_pathname, char_solnname

  integer(kind=4)    :: index_file, ier
  character(len=132) :: linkpath, linkpathg, linkpathe
  character(len=32)  :: zonename

  !open solution file
  call cg_open_f(trim(char_solnname), CG_MODE_MODIFY, index_file, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    
  zonename='Zone0001'

  linkpath = '/Base/' // trim(zonename)

  call cg_gopath_f(index_file, trim(linkpath), ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  linkpathg = trim(linkpath) // '/GridCoordinates'

  call cg_link_write_f('GridCoordinates', trim(char_pathname), trim(linkpathg), ier)
    if (ier .ne. CG_OK) call cg_error_print_f
      
  linkpathe = trim(linkpath) // '/Elems'

  call cg_link_write_f('Elems', trim(char_pathname), trim(linkpathe), ier)
    if (ier .ne. CG_OK) call cg_error_print_f

  call cg_close_f(index_file, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  return

end subroutine write_link_CGNS
