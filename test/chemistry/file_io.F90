module file_io

#define CHECK_STATUS(status, msg) call check_status(status, msg, __FILE__, __LINE__)

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private

  public :: file_io_t

  integer, parameter :: kUnknownFileId = -9999

  type :: file_io_t
    character(len=:), allocatable :: path_
    integer :: id_ = kUnknownFileId
  contains
    procedure :: read_1D_int
    procedure :: read_1D_real
    procedure :: read_2D_int
    procedure :: read_2D_real
    procedure :: read_3D_int
    procedure :: read_3D_real
    procedure :: read_4D_int
    procedure :: read_4D_real
    generic :: read => read_1D_int, read_1D_real, read_2D_int, read_2D_real, &
                       read_3D_int, read_3D_real, read_4D_int, read_4D_real
    final :: finalize
  end type file_io_t

  interface file_io_t
    procedure :: constructor
  end interface file_io_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor(path, do_log) result(self)

    use netcdf, only : nf90_open, nf90_inquire, nf90_inquire_dimension, &
                       nf90_inquire_variable, NF90_NOWRITE, NF90_MAX_NAME

    character(len=*), intent(in) :: path
    logical, intent(in), optional :: do_log
    
    type(file_io_t), pointer :: self
    integer :: nvars, ndims, varid, dimids(10)
    integer, allocatable :: dims(:)
    character(len=NF90_MAX_NAME) :: varname
    integer :: i

    allocate(self)

    self%path_ = path
    self%id_ = 0

    CHECK_STATUS(nf90_open(self%path_, NF90_NOWRITE, self%id_), \
                 "Error opening file: "//trim(self%path_))
    if (present(do_log)) then
      if (do_log) then
        write(*,*) "Opened file: ", trim(self%path_)

        ! Get the number of variables and dimensions in the file
        CHECK_STATUS(nf90_inquire(self%id_, nVariables=nvars), \
                     "Error inquiring number of variables")
        CHECK_STATUS(nf90_inquire(self%id_, nDimensions=ndims), \
                     "Error inquiring number of dimensions")

        allocate(dims(ndims))
        do i = 1, ndims
          CHECK_STATUS(nf90_inquire_dimension(self%id_, i, len=dims(i)), \
                      "Error inquiring dimension length")
        end do

        ! Loop over all variables and print their names and dimensions
        do i = 1, nvars
          CHECK_STATUS(nf90_inquire_variable(self%id_, i, name=varname, \
                                             ndims=ndims, dimids=dimids), \
                       "Error inquiring variable")
          write(*,*) "Variable: ", trim(varname), " Dimensions: ", dims(dimids(1:ndims))
        end do
      end if
    end if

  end function constructor  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_1D_int(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    integer, dimension(:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_1D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_1D_real(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    real(r8), dimension(:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_1D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_2D_int(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    integer, dimension(:,:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_2D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_2D_real(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    real(r8), dimension(:,:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_2D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_3D_int(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    integer, dimension(:,:,:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_3D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_3D_real(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    real(r8), dimension(:,:,:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_3D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_4D_int(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    integer, dimension(:,:,:,:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_4D_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_4D_real(self, variable_name, data)

    use netcdf, only : nf90_inq_varid, nf90_get_var

    class(file_io_t), intent(in) :: self
    character(len=*), intent(in) :: variable_name
    real(r8), dimension(:,:,:,:), intent(out) :: data

    integer :: var_id

    CHECK_STATUS(nf90_inq_varid(self%id_, trim(variable_name), var_id), \
                 "Error getting variable id: '"//trim(variable_name)//"'")
    CHECK_STATUS(nf90_get_var(self%id_, var_id, data), \
                 "Error reading variable: '"//trim(variable_name)//"'")

  end subroutine read_4D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize(self)

    use netcdf, only : nf90_close

    type(file_io_t), intent(inout) :: self

    if (self%id_ .ne. kUnknownFileId) then
      CHECK_STATUS(nf90_close(self%id_), \
                   "Error closing file: "//trim(self%path_))
      self%id_ = kUnknownFileId
    end if
    if (allocated(self%path_)) deallocate(self%path_)

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_status(status, msg, file, line)

    use netcdf, only : nf90_strerror, NF90_NOERR

    integer, intent(in) :: status
    character(len=*), intent(in) :: msg
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    if (status /= NF90_NOERR) then
      write(*,*) "Error: ", trim( nf90_strerror(status) ), " - ", trim(msg)
      write(*,*) "File: ", file
      write(*,*) "Line: ", line
      stop 3
    end if

  end subroutine check_status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module file_io