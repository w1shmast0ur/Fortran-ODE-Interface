!!!
!!! ode utils for Lorenz95 tests
!!!
module odeutils

  use mcmcprec
  implicit none
  public

contains

!!!
!!! RK45
!!!
function rk45solve(odefun,time,init,params,opt) result(y)
    real(kind=dbl), intent(in) :: init(:), time(:)
    real(kind=dbl), intent(in) ::  params(:)
    real(kind=dbl), intent(in) ::  opt(:)
    real(kind=dbl) :: y(size(time),size(init))
    interface 
        function odefun(t,y,params) result(ydot)
          import :: dbl         
          implicit none
          real(kind=dbl), intent(in) :: t
          real(kind=dbl), intent(in), target :: y(:)
          real(kind=dbl), intent(in), target :: params(:)
          real(kind=dbl) :: ydot(size(y))
        end function odefun
     end interface

    real(kind=dbl), allocatable :: k1(:),k2(:),k3(:),k4(:)
    real(kind=dbl), allocatable :: tspan(:), xx(:,:)
    real(kind=dbl) :: dt,t1,t2,t12
    integer :: i, j, ntime, neq
    
    !! the first el. of opt for this solver is inner dt-step
    dt = opt(1)
    !! use dt time stepping internally
    ntime = floor((time(size(time))-time(1)) / dt)+1
    !! avoid too small and too big intervals
    if (ntime<2 .or. ntime > 1000000) then
      write(*,*) 'ntime:',ntime
      stop 'ntime too small (or big), see lorenz95.F90'
    end if
    !! dimension of the problem
    neq = size(init)
    !! allocate auxillary arrays
    allocate(tspan(ntime),xx(ntime,neq))    
    tspan = linspace(time(1),time(size(time)),ntime)
    !! allocate coef. of rg  
    allocate(k1(neq),k2(neq),k3(neq),k4(neq))
    
    xx(1,:) = init
    do i=1,ntime-1
       t1 = tspan(i); t2 = tspan(i+1); dt = t2-t1; t12 = t1+dt/2.0
       k1 = dt * odefun(t1,  xx(i,:), params)
       k2 = dt * odefun(t12, xx(i,:) + k1/2.0, params)
       k3 = dt * odefun(t12, xx(i,:) + k2/2.0, params)
       k4 = dt * odefun(t2,  xx(i,:) + k3, params)
       xx(i+1,:) = xx(i,:) + (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
    end do  
    
    !! interpolate results to given input times
    !$omp parallel do
      do j=1,size(y,2)
         do i=1,size(y,1)
            y(i,j) = interp1(tspan,xx(:,j),time(i))
         end do
      end do
    !$omp end parallel do
    deallocate(k1,k2,k3,k4,tspan,xx)
end function rk45solve

!!!
!!! linspace ala Matlab
!!!
function linspace(a,b,n) result(y)
    implicit none
    real(kind=dbl), intent(in) :: a, b
    integer, intent(in) :: n
    real(kind=dbl), dimension(n) :: y

    integer :: i
    real(kind=dbl) :: step

    step = (b-a)/(n-1.0d0)
    y(1) = a
    y(n) = b
    do i=2,n-1
       y(i) = y(i-1) + step
    end do
end function linspace

!!!
!!! 1d table lookup by linear interpolation
!!!
function interp1(x,y,xin) result(yout)
    implicit none
    real(kind=dbl) :: x(:), y(:), xin
    real(kind=dbl) :: yout

    integer :: n, il, iu, im
    real(kind=dbl) :: x1, x2, y1, y2

    n = size(x,1)
    if (size(y,1) .ne. n) then
       stop 'wrong lenghts in interp1'
    end if

    !!  find location by bisection
    il=1
    iu=n
    do while (iu-il>1)
       im = floor((iu+il)/2.0)
       if (xin>x(im)) then
          il=im
       else
          iu=im
       end if
    end do
    x1 = x(il)
    x2 = x(iu)
    y1 = y(il)
    y2 = y(iu)
    ! linear interpolation 
    yout = y1 + (y2-y1)/(x2-x1)*(xin-x1) 
end function interp1
end module odeutils
