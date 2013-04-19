# Set fortran compiler here
FC=ifort


%.o : %.f90; $(FC) -c $(FFLAGS) $< -o $@

OBJECTS= CPLAP.o

CPLAP:	$(OBJECTS)
	rm -f CPLAP
	$(FC) -o CPLAP $(OBJECTS) 
	rm -f CPLAP.o

clean:
	-rm -f *.o; touch *.f90

