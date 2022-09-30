SDSL_DIR=../sdsl/

CXX = g++
CXXFLAGS = -O3 -std=c++11 
SDSLFLAGS = -ffast-math -funroll-loops -msse4.2 -I $(SDSL_DIR)include/ -L $(SDSL_DIR)lib/ -lsdsl -ldivsufsort -ldivsufsort64 -lpthread

OBJS = io.o gapinfo.o vcf.o custom_int_vector.o saa.o 

PROGRAMS = buildSAA buildFMA locate extract

all: $(PROGRAMS)

buildSAA: buildSAA.cpp $(OBJS)
	$(CXX) -o $@ buildSAA.cpp $(OBJS) $(CXXFLAGS) $(SDSLFLAGS)

buildFMA: buildFMA.cpp gapinfo.o
	$(CXX) -o $@ buildFMA.cpp gapinfo.o $(CXXFLAGS) $(SDSLFLAGS)

locate: locate.cpp gapinfo.o
	$(CXX) -o $@ locate.cpp gapinfo.o $(CXXFLAGS) $(SDSLFLAGS)

extract: extract.cpp gapinfo.o
	$(CXX) -o $@ extract.cpp gapinfo.o $(CXXFLAGS) $(SDSLFLAGS)

io.o: io.cpp io.hpp
	$(CXX) -c io.cpp $(CXXFLAGS) $(SDSLFLAGS)

gapinfo.o: gapinfo.cpp gapinfo.hpp
	$(CXX) -c gapinfo.cpp $(CXXFLAGS) $(SDSLFLAGS)

vcf.o: vcf.cpp vcf.hpp io.o gapinfo.o
	$(CXX) -c vcf.cpp $(CXXFLAGS) $(SDSLFLAGS)

custom_int_vector.o: custom_int_vector.cpp custom_int_vector.hpp
	$(CXX) -c custom_int_vector.cpp $(CXXFLAGS) $(SDSLFLAGS)

saa.o: saa.cpp saa.hpp gapinfo.o io.o custom_int_vector.o
	$(CXX) -c saa.cpp $(CXXFLAGS) $(SDSLFLAGS)

clean:
	rm -f saa.a
	rm -f $(PROGRAMS)
	rm -f *.o
