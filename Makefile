################################
CXX = mpicc
DEBUG = -g
################################

CXXFLAGS = -std=c99 -Wall -ferror-limit=1 $(DEBUG)
################################

################################

################################

all: Hw2_serial pseq ppar

Hw2_serial_dep = Hw2_serial.o
Hw2_serial: $(Hw2_serial_dep)
	$(CXX) $(CXXFLAGS) $(Hw2_serial_dep) -o $@

Hw2_serial.o: Hw2_serial.c
	$(CXX) $(CXXFLAGS) -c -o $@ Hw2_serial.c

pseq_dep = pseq.o
pseq: $(pseq_dep)
	$(CXX) $(CXXFLAGS) $(pseq_dep) -o $@

pseq.o: pseq.c
	$(CXX) $(CXXFLAGS) -c -o $@ pseq.c

ppar_dep = ppar.o
ppar: $(ppar_dep)
	$(CXX) $(CXXFLAGS) $(ppar_dep) -o $@

ppar.o: ppar.c
	$(CXX) $(CXXFLAGS) -c -o $@ ppar.c


clean:
	$(RM) *.o *.out

realclean: clean
	$(RM) -r *.dSYM
	$(RM) *~
	$(RM) -r docs