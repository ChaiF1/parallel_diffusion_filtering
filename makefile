F90=gfortran -O
OBJ = gif_module.o 

phantom:phantom.o $(OBJ)
	$(F90) phantom.o $(OBJ) -o phantom

gif_module.o:   gif_module.f90
	$(F90) $(FFLAGS) -c gif_module.f90

phantom.o:phantom.f90 $(OBJ)
	$(F90) -c phantom.f90 

FC=caf
# use these flags for optimized build
FFLAGS=-O3 -march=native -fopenmp -std=f2018 -cpp -ffree-line-length-none

NP ?= 4

parallel_phantom: parallel_phantom.o $(OBJ)
	$(FC) parallel_phantom.o $(OBJ) -o parallel_phantom

parallel_phantom.o: parallel_phantom.f90 $(OBJ)
	$(FC) $(FFLAGS) -c parallel_phantom.f90

run: parallel_phantom
	mpirun -np $(NP) ./parallel_phantom


clean:
	rm -f *.o *.mod *.gif phantom parallel_phantom
