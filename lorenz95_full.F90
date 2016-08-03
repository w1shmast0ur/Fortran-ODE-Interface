!!!
!!!
!!! Stochastic Lorenz 95 module
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
  
  real(kind=dbl), pointer :: F(:), Fy(:), h, c, b
  real(kind=dbl), pointer :: xslow(:), xfast(:)
  real(kind=dbl), pointer :: dxfast(:), dxslow(:)
  integer :: J, K, nfast, jind, kind, nparams

  J = 8
  K = 40

  if (size(init) /= J*K+K) then
    write(*,*) 'The dimension of the initial state is wrong'
    call exit(1)
  end if

  nparams = J*K+K+3
  if (size(params) /= nparams) then
    write(*,*) 'The dimension of the parameter vector is wrong'
    call exit(1)
  end if

  nfast = J*K
  xslow => init(1:K)
  xfast => init(K+1:)
  dxslow => ydot(1:K)
  dxfast => ydot(K+1:)
  F => params(1:K)
  Fy => params(K+1:J*K+K)
  b => params(nparams-2)
  c => params(nparams-1)
  h => params(nparams)
  
  !! 1) the slow part
  dxslow(1) = -xslow(K-1)*xslow(K)+xslow(K)*xslow(2)-xslow(1) + & 
        F(1)-h*c/b*sum(xfast(1:J))
  dxslow(2) = -xslow(K)*xslow(1)+xslow(1)*xslow(3)-xslow(2) + &
        F(2)-h*c/b*sum(xfast(J+1:2*J))
  !$omp parallel do
    do kind=3,K-1
      dxslow(kind) = -xslow(kind-2)*xslow(kind-1) + &
         xslow(kind-1)*xslow(kind+1)-xslow(kind) + &
         F(kind)-h*c/b*sum(xfast(J*(kind-1)+1:J*kind))
    end do
  !$omp end parallel do
  dxslow(K) = -xslow(K-2)*xslow(K-1)+xslow(K-1)*xslow(1)-xslow(K) + &
        F(K)-h*c/b*sum(xfast(J*(K-1):J*K))
  !! 2) the fast part
  dxfast(1) = -c*b*xfast(2)*(xfast(3)-xfast(nfast))-c*xfast(1)+c/b*Fy(1)+h*c/b*xslow(1)
  !$omp parallel do
     do jind=2,J*K-2
        dxfast(jind) = -c*b*xfast(jind+1)*(xfast(jind+2)-xfast(jind-1))-c*xfast(jind) + &
             c/b*Fy(jind)+h*c/b*xslow(1+floor(dble(jind-1)/dble(J)))
     end do
  !$omp end parallel do
  dxfast(nfast-1) = -c*b*xfast(nfast)*(xfast(1)-xfast(nfast-2))-c*xfast(nfast-1)+ &
        c/b*Fy(nfast-1)+h*c/b*xslow(K)
  dxfast(nfast) = -c*b*xfast(1)*(xfast(2)-xfast(nfast-1))-c*xfast(nfast)+ & 
        c/b*Fy(nfast)+h*c/b*xslow(K)
end function odefun
end module ode