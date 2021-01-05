COMPILER= gfortran

FLAGS = 

EXEC = schrodinger

SRC = $(wildcard *.f90) 

OBJ = $(SRC:.f90=.o)

$(EXEC): $(OBJ)
	$(COMPILER) $(FLAGS) -o $@ $^ $(LAPACK)

types.o: types.f90
	$(COMPILER) $(FLAGS) -c $<

read_write.o: read_write.f90 quantum.o types.o
	$(COMPILER) $(FLAGS) -c $<

linear_algebra.o: linear_algebra.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

hamiltonian.o: hamiltonian.f90 types.o
	$(COMPILER) $(FLAGS) -c $<

quantum.o: quantum.f90 hamiltonian.o linear_algebra.o types.o
	$(COMPILER) $(FLAGS) -c $<

main.o: main.f90 types.o read_write.o quantum.o
	$(COMPILER) $(FLAGS) -c $<

clean:
	rm -rf *.o *.mod

mrproper: clean
	rm -rf $(EXEC)
