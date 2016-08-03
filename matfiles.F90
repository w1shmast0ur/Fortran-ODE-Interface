!!! $Id: matfiles.F90,v 1.12 2012/11/22 11:40:26 mjlaine Exp $
!!! ------------------------------------------------------------------------
!!! mcmc library
!!! File: matfiles.F90
!!! Purpose: Making and reading Matlab V4 MAT-files
!!!
!!! Marko Laine 2001 <marko.laine@fmi.fi>
!!! Copyrights licensed under a MIT License.
!!! See the accompanying LICENSE.txt file for terms.
!!! ------------------------------------------------------------------------
!!! 
module matfiles

#if defined(__IFC) || defined(__INTEL_COMPILER)
use ifport
#endif

#if defined(__PGI)
interface
   integer function fseek(lu, offset, from)
     integer lu, offset, from
   end function fseek
end interface
#endif

  interface writemat
     module procedure writemat4_mat, writemat4_vec
  end interface

  interface readmat
     module procedure readmat4_mat, readmat4_vec
  end interface readmat

  interface addtomat
     module procedure addtomat4_mat
  end interface

  !! mat v4 header
  type Fmatrix
     sequence
     integer(kind=4) :: type 
     integer(kind=4) :: mrows
     integer(kind=4) :: ncols
     integer(kind=4) :: imagf
     integer(kind=4) :: namelen
  end type Fmatrix

contains

!!!
!!! Write a matrix with 32 bit reals into matlab Level 1.0 (V4) MAT file.
!!!
!!! Format is
!!! 5*4 = 20 first bytes:
!!!    type  = 0
!!!    mrows = number of rows
!!!    ncols = number of columns
!!!    imagf = 0
!!!    namelen = length of the name + 1 for achar(0)
!!! namelen next bytes:  name_of_x+\0 as char bytes
!!! mrow*ncol*8 next bytes: data as 8 byte floats
!!!

subroutine writemat4_mat(fname, x, xname, status)
  implicit none
  real(kind=8), intent(in) :: x(:,:)
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: xname
  integer, intent(out), optional :: status

  integer :: fstat, i, j
  type(Fmatrix) :: xtype

  if (present(status)) status = 0

#if defined(__GFORTRAN__) || defined(AIX)
  open(unit=101, file=fname, status='replace', &
       form='unformatted',access='stream',iostat=fstat)
#else
  open(unit=101, file=fname, status='replace', &
       form='binary', position='rewind', iostat=fstat)
#endif
  if (fstat /= 0) then
     write(*,*) 'Error opening file ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
     else
        stop
     end if
     return
  end if

  xtype%type    = 0 ! Little Endian, 1000 for Big Endian
  xtype%mrows   = size(x,1)
  xtype%ncols   = size(x,2)
  xtype%imagf   = 0
  xtype%namelen = len_trim(xname)+1

!!! write x to file in binary form
! write(101, iostat=fstat) xtype, trim(xname)//achar(0), x(:,:)
  write(101, iostat=fstat) xtype
  write(101, iostat=fstat) trim(xname)//achar(0)

  if (fstat /= 0) then
     write(*,*) 'Error writing to file', fname
     if (present(status)) then
        status = fstat
        close(101)
     else
        stop
     end if
     return
  end if

!  write(101, iostat=fstat) x(:,:)
!  write(101, iostat=fstat) x
! write x in loop to avoid stack overflow
  do j=1,size(x,2)
     do i=1,size(x,1)
        write(101) x(i,j)
     end do
  end do
  close(101)
end subroutine writemat4_mat


subroutine writemat4_vec(fname, x, xname, status)
  implicit none
  real(kind=8), intent(in) :: x(:)
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: xname
  integer, intent(out), optional :: status

  integer :: fstat, i
  type(Fmatrix) :: xtype

  if (present(status)) status = 0

#if defined(__GFORTRAN__) || defined(AIX)
  open(unit=101, file=fname, status='replace', &
       form='unformatted', position='rewind', access='stream', iostat=fstat)
#else
  open(unit=101, file=fname, status='replace', &
       form='binary', position='rewind', iostat=fstat)
#endif
  if (fstat /= 0) then
     write(*,*) 'Error opening file', fname
     if (present(status)) then
        status = fstat
     else
        stop
     end if
     return
  end if

  xtype%type    = 0 ! Little Endian, 1000 for Big Endian
  xtype%mrows   = size(x,1)
  xtype%ncols   = 1
  xtype%imagf   = 0
  xtype%namelen = len_trim(xname)+1

!!! write x to file in binary form
!  write(101, iostat=fstat) xtype, trim(xname)//achar(0), x(:)
  write(101, iostat=fstat) xtype, trim(xname)//achar(0)
  if (fstat /= 0) then
     write(*,*) 'Error writing to file', fname
     if (present(status)) then
        status = fstat
        close(101)
     else
        stop
     end if
     return
  end if
! write in loop to avoid stack overflow
  do i=1,size(x)
     write(101) x(i)
  end do
  close(101)
end subroutine writemat4_vec

!!!
!!! Adding columns to existing mat v4 file.
!!! 
subroutine addtomat4_mat(fname, x, status)
  implicit none
  real(kind=8), intent(in) :: x(:,:)
  character(len=*), intent(in) :: fname
  integer, intent(out), optional :: status

  integer :: fstat !, reclen
  type(Fmatrix) :: xtype

  if (present(status)) status = 0

 ! inquire(iolength=reclen) 
#if defined(__GFORTRAN__) || defined(AIX)
  open(unit=101, file=fname, status='old', &
       form='unformatted', position='rewind', access='stream', iostat=fstat)
#else
  open(unit=101, file=fname, status='old', &
       form='binary', position='rewind', access='direct', iostat=fstat)
#endif
  if (fstat /= 0) then
     write(*,*) 'Error opening file ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
        return
     else
        stop
     end if
  end if

  !! read the header
  read(101, iostat=fstat) xtype
  if (fstat /= 0) then
     write(*,*) 'Error reading file', fname
     close(101)
     if (present(status)) then
        status = fstat
        return
     else
        stop
     end if
  end if

  !! check the size
  if (xtype%mrows > 0 .and. xtype%mrows /= size(x,1) ) then
     write(*,*) 'error: xmat should have same number of rows'
     write(*,*) xtype
     write(*,*) size(x,2)
     if (present(status)) then
        status = -1
        return
     else
        stop
     end if
  end if

  !! goto end of file
  !! how to do it in fortran?  call fseek(101,0,2)?
  !! NOW: close and open with append
  close(101)
#if defined(__GFORTRAN__) || defined(AIX)
  open(unit=101, file=fname, status='old', &
       form='unformatted', position='append', access='stream', iostat=fstat)
#else
  open(unit=101, file=fname, status='old', &
       form='binary', access='append', iostat=fstat)
#endif
  if (fstat /= 0) then
     write(*,*) 'Error opening file (2)', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
        return
     else
        stop
     end if
  end if

!! append the new data
  write(101, iostat=fstat) x
  if (fstat /= 0) then
     write(*,*) 'Error writing new data to ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
        return
     else
        stop
     end if
  end if

!! goto bof
  rewind(101)
!! write new header
  xtype%mrows = size(x,1)
  xtype%ncols = xtype%ncols + size(x,2)
  write(101, iostat=fstat) xtype
  if (fstat /= 0) then
     write(*,*) 'Error writing new header ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
        return
     else
        stop
     end if
  end if

!! close file
  close(101)
end subroutine addtomat4_mat

!!!
!!! read the first matrix form V4 mat file
!!!
subroutine readmat4_mat(fname,x,status)
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=8), pointer :: x(:,:)
  integer, intent(out), optional :: status

  integer :: fstat
  type(Fmatrix) :: xtype

#if defined(__GFORTRAN__) || defined(AIX)
  open(unit=101, file=fname, status='old', &
       form='unformatted',access='stream',iostat=fstat)
#else
  open(unit=101, file=fname, status='old', &
       form='binary', position='rewind', iostat=fstat)
#endif
  if (fstat /= 0) then
     write(*,*) 'Error opening file ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
     else
        stop
     end if
     return
  end if

  read(101) xtype
  if (xtype%imagf.ne.0.or.xtype%imagf.ne.0) stop 'Not a mat file?'

#if defined(__IFC) || defined(__INTEL_COMPILER)
  fstat = fseek(lunit=101,from=1,offset=xtype%namelen)
#elif defined(__PGI)
  fstat = fseek(lu=101,from=1,offset=xtype%namelen)
#else
  call fseek(unit=101,whence=1,offset=xtype%namelen)
#endif

  allocate(x(xtype%mrows,xtype%ncols))
  read(101,iostat=fstat) x
  if (fstat /= 0) then
     write(*,*) 'Error reading the file ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
     else
        stop
     end if
     return
  end if
  close(101)
  if (present(status)) status=0
end subroutine readmat4_mat

subroutine readmat4_vec(fname,x,status)
  implicit none
  character(len=*), intent(in) :: fname
  real(kind=8), pointer :: x(:)
  integer, intent(out), optional :: status

  integer :: fstat
  type(Fmatrix) :: xtype

#if defined(__GFORTRAN__) || defined(AIX)
  open(unit=101, file=fname, status='old', &
       form='unformatted',access='stream',iostat=fstat)
#else
  open(unit=101, file=fname, status='old', &
       form='binary', position='rewind', iostat=fstat)
#endif
  if (fstat /= 0) then
     write(*,*) 'Error opening file ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
     else
        stop
     end if
     return
  end if

  read(101) xtype
  if (xtype%imagf.ne.0.or.xtype%imagf.ne.0) stop 'Not a mat file?'

#if defined(__IFC) || defined(__INTEL_COMPILER)
  fstat = fseek(lunit=101,from=1,offset=xtype%namelen)
#elif defined(__PGI)
  fstat = fseek(lu=101,from=1,offset=xtype%namelen)
#else
  call fseek(unit=101,whence=1,offset=xtype%namelen)
#endif

  allocate(x(xtype%mrows*xtype%ncols))
  read(101,iostat=fstat) x
  if (fstat /= 0) then
     write(*,*) 'Error reading the file ', trim(fname),' fstat: ',fstat
     if (present(status)) then
        status = fstat
     else
        stop
     end if
     return
  end if
  close(101)
  if (present(status)) status=0
end subroutine readmat4_vec


end module matfiles
