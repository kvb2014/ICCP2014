FC = gfortran
FFLAGS = -Wall -Wextra -O3 -ffast-math 
LDFlAGS = 
LIBS =  -llapack

PROGRAMS = threelevel

OBJS =
#OBJS += plot.o

all: $(PROGRAMS) 

main.o: $(OBJS)

main: $(OBJS)

%: %.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

%.o: %.f95
	$(FC) $(FFLAGS) -c $<
%.o: %.F95
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
%.o: %.F90 
	$(FC) $(FFLAGS) -c $<
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.F 
	$(FC) $(FFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD $(PROGRAMS)

veryclean: clean
	rm -f *~ $(PROGRAMS) 
