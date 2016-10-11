CXX    = g++
F90    = gfortran
FLAG   = -g -fcheck=all -fcheck=bounds -fbacktrace -fdefault-real-8 -fdefault-double-8 -O0
LIBS   = -L/usr/local/lib/ -lnleq1 -lblas -llapack
INCSAD = -I$(HOME)/Codes/cppad-2.0.0/cppad-2.0.0/include/
LINK   = -lgfortran

all: exec

exec: eigen.o rhsjac.o chem.o csp.o main.o
	$(CXX) -o exec eigen.o rhsjac.o chem.o csp.o main.o $(LIBS) $(LINK)

eigen.o : eigen.f90
	$(F90) -c $(FLAG) eigen.f90

rhsjac.o : rhsjac.cpp
	$(CXX) -c -g $< $(INCSAD)

chem.o : chem.f90 eigen.o rhsjac.o
	$(F90) -c $(FLAG) $< 

csp.o : csp.f90 chem.o eigen.o
	$(F90) -c $(FLAG) $< 

main.o : main.f90 csp.o chem.o eigen.o
	$(F90) -c $(FLAG) $< 

clean:
	rm *.o

real-clean:
	rm *~

