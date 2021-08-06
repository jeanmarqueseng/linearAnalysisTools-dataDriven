! This is a wrap for reading CGNS files.
! 
module mod_CGNS

  integer(kind=4) :: cell_dim, phys_dim
  integer(kind=4) :: ifile, ibase, izone, ifield, iflow, icoord
  integer :: ier
  
  integer(kind=4), allocatable :: isize(:,:)  
  
  character(len=8) :: zonename, basename
  
  logical :: check
  logical :: CGNS_file_exists
  
end module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! 
! 
! 
! 
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine nzones_read(ifile, ibase, numzones)

  use cgns
  implicit none
  integer, intent(in)  :: ifile, ibase
  integer, intent(out) :: numzones
  integer :: ier

  call cg_nzones_f(ifile,ibase,numzones,ier)
  if (ier.ne.CG_OK) call cg_error_exit_f

  return
  
end subroutine nzones_read
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine zonedim_read(ifile, ibase, izone, idim)

  use cgns
  implicit none
  integer, intent(in)  :: ifile, ibase, izone
  integer, intent(out) :: idim
  integer :: ier

  call cg_index_dim_f(ifile,ibase,izone,idim,ier)
  if(ier.ne.CG_OK) call cg_error_exit_f

  return

end subroutine zonedim_read
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Reads NX, NY, NZ dimensions on each zone
subroutine read_size_CGNS(CGNS_filename,zone_number,ngrid)

  use cgns
  use mod_CGNS
  implicit none
  
  character(len=*), intent(in) :: CGNS_filename
  integer(kind=4),intent(in)   :: zone_number
  integer(kind=4), intent(out) :: ngrid(3)

  inquire(file=trim(CGNS_filename),exist=CGNS_file_exists)
    if (CGNS_file_exists .eqv. .false.) then
      write(*,'(A)') trim(CGNS_filename)
      stop ' O arquivo CGNS não existe!!! ' 
    endif

  !... open CGNS file
  call cg_open_f(trim(CGNS_filename), CG_MODE_READ, ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  ibase = 1
  call cg_base_read_f(ifile, ibase, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  izone = zone_number
  call cg_zone_read_f(ifile, ibase, izone, zonename, isize, ier)

  ngrid(1) = isize(1,1)
  ngrid(2) = isize(2,1)
  if (cell_dim .gt. 2) then
    ngrid(3) = isize(3,1)
  else
    ngrid(3) = 1
  endif

  print*, 'ngrid(1) =', ngrid(1)
  print*, 'ngrid(2) =', ngrid(2)
  print*, 'ngrid(3) =', ngrid(3)

  call cg_close_f(ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine read_size_CGNS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine create_file_cgns(CGNS_solnname,aux_dim)

  use cgns
  use mod_CGNS
  implicit none
  
  character(len=*), intent(in) :: CGNS_solnname
  character(len=*), intent(in) :: aux_dim

  
  
  
  check = .false.
  
  inquire(file=trim(CGNS_solnname),exist=CGNS_file_exists)
  if (CGNS_file_exists .eqv. .false.) then

    !open solution file
    call cg_open_f(trim(CGNS_solnname),CG_MODE_WRITE, ifile,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

    !create base in the recently opened solution file
    if (aux_dim .eq. '2D' .or. aux_dim .eq. '2d') then
      call cg_base_write_f(ifile, 'Base', 2, 2, ibase, ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
      check = .true.
    endif

    if (aux_dim .eq. '3D' .or. aux_dim .eq. '3d') then
      call cg_base_write_f(ifile, 'Base', 3, 3, ibase, ier)
        if (ier .ne. CG_OK) call cg_error_exit_f
      check = .true.
    endif
    
    if (check .eqv. .false.) stop ' Tem cagada na chamada do create_file_CGNS ... '

    !close the file
    call cg_close_f(ifile, ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    write(*,'(A)') trim(CGNS_solnname) // ' already exists!'

  endif

  return

end subroutine create_file_CGNS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine write_link_CGNS(zone_number,CGNS_solnname,char_pathname)

  use cgns
  use mod_CGNS
  implicit none
  
  integer(kind=4), intent(in)  :: zone_number
  character(len=*), intent(in) :: CGNS_solnname
  character(len=*), intent(in) :: char_pathname
  integer(kind=4) :: m

  character(len=132) linkpath
  
  

  write(*,'(A)') '  -> Linking file ' // trim(CGNS_solnname) // ' to ' // trim(char_pathname)
  inquire(file=trim(CGNS_solnname),exist=CGNS_file_exists)
    if (CGNS_file_exists .eqv. .false.) then
      print*, trim(CGNS_solnname)
      stop ' O arquivo CGNS não existe!!!'
    endif
  
  !open solution file
  call cg_open_f(trim(CGNS_solnname), CG_MODE_MODIFY, ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f
    
    write(zonename,'(a4,i4.4)') 'Zone', zone_number

    linkpath = '/Base/' // trim(zonename)

    print*, trim(linkpath)
    call cg_gopath_f(ifile, trim(linkpath), ier)
      if (ier .ne. CG_OK) call cg_error_print_f

    linkpath = trim(linkpath) // '/GridCoordinates'

    call cg_link_write_f('GridCoordinates', trim(char_pathname), trim(linkpath), ier)
      if (ier .ne. CG_OK) call cg_error_print_f
      
    write(*,'(A,i0)') 'Linking zone ', zone_number

  call cg_close_f(ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  return

end subroutine write_link_CGNS

! 
! 
! 
! 
! 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_2D_coord(var,ifile, ibase, izone, ijk_min, ijk_max, nx, ny, x)

  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone
  integer, dimension(:), intent(in) :: ijk_min, ijk_max
  integer, intent(in) :: nx, ny
  real(kind=8), dimension(nx, ny), intent(out) :: x
  integer :: ier

  call cg_coord_read_f(ifile,ibase,izone,trim(var),RealDouble,ijk_min,ijk_max,x(:,:),ier)
    if(ier.ne.CG_OK) call cg_error_print_f
  
  return
    
end subroutine read_2D_coord
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_3D_coord(var,ifile, ibase, izone, ijk_min, ijk_max, nx, ny, nz, x)

  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone
  integer, dimension(:,:), intent(in) :: ijk_min, ijk_max
  integer, intent(in) :: nx, ny, nz
  real(kind=8), dimension(nx, ny, nz), intent(out) :: x
  integer :: ier
  
  call cg_coord_read_f(ifile,ibase,izone,trim(var),RealDouble,ijk_min,ijk_max,x(:,:,:),ier)
    if(ier.ne.CG_OK) call cg_error_exit_f
  
  return
    
end subroutine read_3D_coord
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine write_2D_coord(CGNS_filename, zone_number, nx, ny, xcoord, ycoord)

  use cgns
  use mod_CGNS
  implicit none

  character(len=32), intent(in) :: CGNS_filename
  
  integer(kind=4), intent(in) :: zone_number
  integer(kind=4), intent(in) :: nx, ny  
  real(kind=8), dimension(:,:), intent(in) :: xcoord, ycoord

   
    
  !... open CGNS file
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,ifile,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  ibase = 1
  call cg_base_read_f(ifile, ibase, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

    write(zonename,'(a4,i4.4)') 'Zone', zone_number

    call cg_gopath_f(ifile, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
    if (ier .eq. CG_NODE_NOT_FOUND) then

      write(*,'(A,i4.4,A)') '  -> Creating zone ', zone_number, ' in ' // trim(CGNS_filename)

      isize(1,1) = nx
      isize(2,1) = ny

      isize(1,2) = isize(1,1) - 1
      isize(2,2) = isize(2,1) - 1

      isize(1,3) = 0
      isize(2,3) = 0
      
      call cg_zone_write_f(ifile, ibase, trim(zonename), isize, Structured, izone, ier)
        if (ier .ne. CG_OK) call cg_error_exit_f

      call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateX', xcoord(:,:), icoord, ier )
        if (ier .ne. CG_OK) call cg_error_exit_f
      call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateY', ycoord(:,:), icoord, ier )
        if (ier .ne. CG_OK) call cg_error_exit_f

    else
    
      write(*,'(A)') 'The grid for ' // trim(zonename) // ' already exist! Skipping ...'

    endif

    call cg_close_f(ifile, ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)
  
  return
    
end subroutine write_2D_coord
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine write_3D_coord(CGNS_filename, zone_number, nx, ny, nz, xcoord, ycoord, zcoord)

  use cgns
  use mod_CGNS
  implicit none

  character(len=32), intent(in) :: CGNS_filename
  
  integer(kind=4), intent(in) :: zone_number
  integer(kind=4), intent(in) :: nx, ny, nz
  real(kind=8), dimension(:,:,:), intent(in) :: xcoord, ycoord, zcoord
    
  !... open CGNS file
  call cg_open_f(trim(CGNS_filename),CG_MODE_MODIFY,ifile,ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  ibase = 1
  call cg_base_read_f(ifile, ibase, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

    write(zonename,'(a4,i4.4)') 'Zone', zone_number

    call cg_gopath_f(ifile, '/Base/'//trim(zonename)//'/GridCoordinates/', ier)
    if (ier .eq. CG_NODE_NOT_FOUND) then

      write(*,'(A,i4.4,A)') '  -> Creating zone ', zone_number, ' in ' // trim(CGNS_filename)

      isize(1,1) = nx
      isize(2,1) = ny
      isize(3,1) = nz

      isize(1,2) = isize(1,1) - 1
      isize(2,2) = isize(2,1) - 1
      isize(3,2) = isize(3,1) - 1

      isize(1,3) = 0
      isize(2,3) = 0
      isize(3,3) = 0    
      
      call cg_zone_write_f(ifile, ibase, trim(zonename), isize, Structured, izone, ier)
        if (ier .ne. CG_OK) call cg_error_exit_f

      call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateX', xcoord(:,:,:), icoord, ier )
        if (ier .ne. CG_OK) call cg_error_exit_f
      call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateY', ycoord(:,:,:), icoord, ier )
        if (ier .ne. CG_OK) call cg_error_exit_f
      call cg_coord_write_f(ifile, ibase, izone, RealDouble, 'CoordinateZ', zcoord(:,:,:), icoord, ier )
        if (ier .ne. CG_OK) call cg_error_exit_f

    else
    
      write(*,'(A)') 'The grid for ' // trim(zonename) // ' already exist! Skipping ...'

    endif

    call cg_close_f(ifile, ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  DEALLOCATE(isize)
  
  return
    
end subroutine write_3D_coord
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! 
! 
! 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_2D_flow(var, ifile, ibase, izone, ijk_min, ijk_max, nx, ny, q)
  
  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone, nx, ny
  integer, dimension(:), intent(in) :: ijk_min, ijk_max
  real(kind=8), dimension(nx, ny), intent(out) :: q
  integer :: iflow, ier

  iflow = 1
  call cg_field_read_f(ifile,ibase,izone,iflow,trim(var),RealDouble,ijk_min,ijk_max,q(:,:),ier)
    if(ier.ne.CG_OK) call cg_error_exit_f
    
  return
    
end subroutine read_2D_flow
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_3D_flow(var, ifile, ibase, izone, ijk_min, ijk_max, nx, ny, nz, q)
  
  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone, nx, ny, nz
  integer, dimension(:,:), intent(in) :: ijk_min, ijk_max
  real(kind=8), dimension(nx, ny, nz), intent(out) :: q
  integer :: iflow, ier

  iflow = 1
  call cg_field_read_f(ifile,ibase,izone,iflow,trim(var),RealDouble,ijk_min, ijk_max,q(:,:,:),ier)
    if(ier.ne.CG_OK) call cg_error_exit_f
    
  return
    
end subroutine read_3D_flow
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine write_soln_2D(CGNS_solnname, zone_number, nx, ny, solution, var_name)

  use cgns
  use mod_CGNS
  implicit none
  
  character(len=*), intent(in) :: CGNS_solnname
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: nx, ny
  real(kind=8), dimension(:,:), intent(in) :: solution
  character(len=*), intent(in) :: var_name

  
  
  character(len=128) :: cdummy

  !open solution file
  call cg_open_f(trim(CGNS_solnname), CG_MODE_MODIFY, ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  ibase = 1
  call cg_base_read_f(ifile, ibase, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  write(zonename,'(a4,i4.4)') 'Zone', zone_number

  call cg_gopath_f(ifile, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', zone_number, ' in ' // trim(CGNS_solnname)

    isize(1,1) = nx
    isize(2,1) = ny

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0

    call cg_zone_write_f(ifile, ibase,trim(zonename), isize, Structured, izone, ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else
  
    izone = zone_number !só funciona quando todas as zonas tiverem sido escritas..

  endif

  call cg_zone_read_f(ifile, ibase, izone, cdummy, isize, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(ifile, '/Base/'//trim(zonename)//'/FlowSolution/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(ifile, ibase, izone, 'FlowSolution', Vertex, iflow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    iflow = 1

  endif

  call cg_field_write_f(ifile, ibase, izone, iflow, RealDouble, trim(var_name), solution(1:nx,1:ny), ifield, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_soln_2D
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine write_soln_3D(CGNS_solnname, zone_number, nx, ny, nz, solution, var_name)

  use cgns
  use mod_CGNS
  implicit none
  
  character(len=*), intent(in) :: CGNS_solnname
  integer(kind=4), intent(in)  :: zone_number
  integer(kind=4), intent(in)  :: nx, ny, nz
  real(kind=8), dimension(:,:,:), intent(in) :: solution
  character(len=*), intent(in) :: var_name

  
  
  character(len=128) :: cdummy

  !open solution file
  call cg_open_f(trim(CGNS_solnname), CG_MODE_MODIFY, ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  ibase = 1
  call cg_base_read_f(ifile, ibase, basename, cell_dim, phys_dim, ier)
    if (ier .ne. CG_OK) call cg_error_print_f
    if (trim(basename) .ne. "Base") stop ' A base nao se chama "Base"... '

  allocate(isize(cell_dim,3))

  write(zonename,'(a4,i4.4)') 'Zone', zone_number

  call cg_gopath_f(ifile, '/Base/'//trim(zonename)//'/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    write(*,'(A,i4.4,A)') '  -> Creating zone ', zone_number, ' in ' // trim(CGNS_solnname)

    isize(1,1) = nx
    isize(2,1) = ny
    isize(3,1) = nz

    isize(1,2) = isize(1,1) - 1
    isize(2,2) = isize(2,1) - 1
    isize(3,2) = isize(3,1) - 1

    isize(1,3) = 0
    isize(2,3) = 0
    isize(3,3) = 0

    call cg_zone_write_f(ifile, ibase,trim(zonename), isize, Structured, izone, ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else
  
    izone = zone_number !só funciona quando todas as zonas tiverem sido escritas..

  endif

  call cg_zone_read_f(ifile, ibase, izone, cdummy, isize, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_gopath_f(ifile, '/Base/'//trim(zonename)//'/FlowSolution/', ier)
  if (ier .eq. CG_NODE_NOT_FOUND) then

    call cg_sol_write_f(ifile, ibase, izone, 'FlowSolution', Vertex, iflow,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f

  else

    iflow = 1

  endif

  call cg_field_write_f(ifile, ibase, izone, iflow, RealDouble, trim(var_name), solution(1:nx,1:ny,1:nz), ifield, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  call cg_close_f(ifile, ier)
    if (ier .ne. CG_OK) call cg_error_exit_f

  deallocate(isize)

  return

end subroutine write_soln_3D
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! 
! 
! 
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine descriptors_read(ifile, ibase, time)

  use cgns
  implicit none

  integer, intent(in) :: ifile, ibase
  real(kind=8), intent(out) :: time
  integer :: ndescriptors, descSize, idesc, ier
  character(len=200) :: descName

  call cg_gopath_f(ifile, '/Base/', ier)
    if(ier.ne.CG_OK) call cg_error_exit_f

  ! go to another directory to read the simulation time
  call cg_gopath_f(ifile, '/Base/TimeIterValues', ier)
    if(ier.ne.CG_OK) call cg_error_exit_f

  ! read simulation time
  call cg_array_read_f(1,time,ier)
    if(ier.ne.CG_OK) call cg_error_exit_f
  
  return
    
end subroutine descriptors_read
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
