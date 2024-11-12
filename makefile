F90=gfortran -O
OBJ = gif_module.o 

phantom:phantom.o $(OBJ)
	$(F90) phantom.o $(OBJ) -o phantom

gif_module.o:   gif_module.f90
	$(F90) $(FFLAGS) -c gif_module.f90

phantom.o:phantom.f90 $(OBJ)
	$(F90) -c phantom.f90 
