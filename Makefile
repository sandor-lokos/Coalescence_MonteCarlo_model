CXX           = g++
CFLAGS        = -O3 -Wall -fPIC -fno-inline -std=c++11 -c
LD            = g++
LDFLAGS       = -O3 -std=c++11
GSLFLAGS      = -lgsl -lgslcblas -lm

all: numerical_integrator.exe

numerical_integrator.o: *.cc *.h
	$(CXX) $(CFLAGS) numerical_integrator.cc

numerical_integrator.exe: numerical_integrator.o
	$(LD) $(LDFLAGS) numerical_integrator.o -o numerical_integrator.exe $(GSLFLAGS)
  
clean:
	rm -rf *.o *.exe
