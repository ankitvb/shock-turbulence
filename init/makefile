F77COMP     = ftn 
F90COMP     = ftn 
CCOMP       = cc 
INCLUDE     = -I/global/u1/a/ankitb/local/gsl/include
LIBS        = -L/global/u1/a/ankitb/local/gsl/lib -lgsl -lgslcblas
OPTS        = -fast 

OBJECT     = initialize_mpi.o define.o genfield.o \
             grid.o cvers.o legendre.o \
             version.o ctof.o trace.o

genfield.out: $(OBJECT)
	$(F77COMP) -o genfield.out $(OBJECT) $(LIBS) -lm  

define.o:  header define.f
	$(F77COMP) -c $(OPTS) define.f

genfield.o:  header genfield.f
	$(F77COMP) -c $(OPTS) genfield.f

genfield_st.o:  header genfield_st.f
	$(F77COMP) -c $(OPTS) genfield_st.f

grid.o: header grid.f
	$(F77COMP) -c $(OPTS) grid.f	

harmonics.o: header harmonics.f
	$(F77COMP) -c $(OPTS) harmonics.f

initialize_mpi.o: header initialize_mpi.f90
	$(F90COMP) -c $(OPTS) initialize_mpi.f90

openfile.o:  openfile.f
	$(F77COMP) -c $(OPTS) openfile.f
 
trace.o: header trace.f
	$(F77COMP) -c $(OPTS) trace.f

legendre.o: legendre.c
	$(CCOMP) -c $(OPTS) $(INCLUDE) legendre.c

version.o:  version.h version.c
	$(CCOMP) -c version.c

cvers.o:  ctof.h version.h cvers.c
	$(CCOMP) -c cvers.c

ctof.o:  ctof.h ctof.c
	$(CCOMP) -c ctof.c

#  remove object files
clean:
	rm -f *.o
