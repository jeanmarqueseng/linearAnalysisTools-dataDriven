! This is a wrap for reading CGNS files.
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
subroutine read_2D_coord(var,ifile, ibase, izone, ijk_min, ijk_max, nx, ny, x)

  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone, nx, ny
  integer, dimension(:), intent(in) :: ijk_min, ijk_max
  real(kind=8), dimension(nx, ny), intent(out) :: x
  integer :: ier

  call cg_coord_read_f(ifile,ibase,izone,trim(var),RealDouble,ijk_min,ijk_max,x(:,:),ier)
    if(ier.ne.CG_OK) call cg_error_exit_f
  
  return
    
end subroutine read_2D_coord
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine read_3D_coord(var,ifile, ibase, izone, isize, nx, ny, nz, x)

  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone, nx, ny, nz
  integer, dimension(:,:), intent(in) :: isize
  real(kind=8), dimension(nx, ny, nz), intent(out) :: x
  integer :: ier
  
  print*, 'BUGADO'
  
  return

  call cg_coord_read_f(ifile,ibase,izone,trim(var),RealDouble,1,isize(:,1),x(:,:,:),ier)

  if(ier.ne.CG_OK) call cg_error_exit_f
  
  return
    
end subroutine read_3D_coord
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
subroutine read_3D_flow(var, ifile, ibase, izone, isize, nx, ny, nz, q)
  
  use cgns
  implicit none

  character(len=16), intent(in)  :: var
  integer, intent(in) :: ifile, ibase, izone, nx, ny, nz
  integer, dimension(:,:), intent(in) :: isize
  real(kind=8), dimension(nx, ny, nz), intent(out) :: q
  integer :: iflow, ier
  
  print*, 'BUGADO'
  
  return

  iflow = 1
  call cg_field_read_f(ifile,ibase,izone,iflow,trim(var),RealDouble,1,isize(:,1),q(:,:,:),ier)
    if(ier.ne.CG_OK) call cg_error_exit_f
    
  return
    
end subroutine read_3D_flow
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
