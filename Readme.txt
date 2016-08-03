To compile the executable, run the following script in cmd:
make ODE=%ode_file_name.F90% all
, where %ode_file_name.F90% is the name of the Fortran file, 
containing the module for implementation of ODE system, which has the following structure:
module ode
  use mcmcprec  ! floating point precision
  implicit none
  public
  ...
  contains
  function odefun(time, init, params) result(ydot)
      implicit none
      real(kind=dbl), intent(in)  :: time
	  real(kind=dbl), intent(in), target ::  init(:)
      real(kind=dbl), intent(in), target ::  params(:)
      real(kind=dbl), target :: ydot(size(init))
      ... 
  end function odefun
end module ode