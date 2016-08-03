!!!
!!!
!!! Simple ode for test purposes
!!!
!!!
module ode
  use mcmcprec  ! floating point precision
  implicit none
  public
!!! Subroutines and functions
contains
!!!
!!! Lorenz95 odefile, both full and the forecast models
!!!
function odefun(time, init, params) result(ydot)
  implicit none
  real(kind=dbl), intent(in)  :: time
  real(kind=dbl), intent(in), target ::  init(:)
  real(kind=dbl), intent(in), target ::  params(:)
  real(kind=dbl), target :: ydot(size(init))
  ydot(1) = cos(time)
end function odefun
end module ode