##
## Lorenz95 f90 code Makefile
##

F90=gfortran
F90FLAGS= -O3 -mtune=native -Wall
DATA=Makefile Readme.txt
SRC=mcmcprec.F90 matfiles.F90 odeutils.F90 ode_solver.F90 $(ODE)
ZIP="c:\Program Files\7-Zip\7z.exe"

%.o : %.F90
	$(F90) $(F90FLAGS) -c $< 

odeutils.o : mcmcprec.o

ode.o : mcmcprec.o $(ODE)
	$(F90) $(F90FLAGS) -c $(ODE) -o ode.o

all: ode_solver

ode_solver: ode_solver.F90 mcmcprec.o odeutils.o matfiles.o ode.o
	$(F90) $(F90FLAGS) -o $@ $^ 

zip: 
	$(ZIP)  a -tzip ode_solver.zip $(SRC) $(DATA)

clean:
	del *.o *.mod *.exe *.mat
