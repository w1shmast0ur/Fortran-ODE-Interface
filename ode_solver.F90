!!!
!!! common interface for solving ODEs
!!! accept time: points where to calculate the solution,
!!! MxN matrix of initial values, where M is dimension of state vector, N is ensemble size,
!!! PxN matrix of parameters, where P is dimension of param.vector, N is ensemble size,
!!! opt is vector containing some options for solver 
program ode_solver  
  use matfiles              ! read/write mat files
  use mcmcprec              ! floating point precision 
  use odeutils              ! rk45 etc utils
  use ode                   ! file contain the implementation of the ode system
  
  implicit none
  character(len=256) :: time_file, init_file, params_file, opt_file
  real(kind=dbl), allocatable :: x(:,:,:)
  real(kind=dbl), pointer :: x2(:,:)
  real(kind=dbl), pointer :: init_ens(:,:), params_ens(:,:)
  real(kind=dbl), pointer :: tspan(:), opt(:)
  integer :: ntime, nens, nstate, nparams, iens
    
  !! namelist file from the command line
  if (command_argument_count() < 4) then
     write(*,*) 'USAGE: <prog> time.mat init.mat params.mat opt.mat'
     call exit(1)
  end if

  !! read configuration
  call get_command_argument(1,time_file)
  call get_command_argument(2,init_file)
  call get_command_argument(3,params_file)
  call get_command_argument(4,opt_file)  

  !! read content of the files
  call readmat(time_file,tspan)
  call readmat(init_file,init_ens)  
  call readmat(params_file,params_ens)
  call readmat(opt_file,opt)
  
  !! get the dimensions
  nstate = size(init_ens,1)
  nens = size(init_ens,2)
  ntime = size(tspan,1)
  nparams = size(params_ens,1)
  
  !! allocate variables with given dimensions  
  allocate(x(ntime,nstate,nens))
  allocate(x2(ntime,nens*nstate))
  
  !! solve the given ode for ensemble
  do iens=1,nens
     x(:,:,iens) = rk45solve(odefun,tspan,init_ens(:,iens),params_ens(:,iens),opt)
  end do  

  !!  save 3d matrix as 2d for matlab -v4 mat file compatibility
  do iens=1,nens
     x2(:,(1+(iens-1)*nstate):(nstate+(iens-1)*nstate)) = x(:,:,iens)
  end do

  !! save output
  call writemat('solution.mat',x2,'solution');

  !! deallocate used variables
  deallocate(x,x2,init_ens,params_ens, tspan, opt) 

end program ode_solver
