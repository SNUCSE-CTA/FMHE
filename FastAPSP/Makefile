CC = g++

CPPFLAGS =  -fopenmp -ansi -pedantic -O3 -std=c++11
#CPPFLAGS =  -fopenmp -ansi -pedantic -std=c++11 -g -ggdb

objects = FastAPSP.o Ctree.o AuxData.o oList.o LinkedList.o

default: $(objects)
	$(CC) $(CPPFLAGS) -o FastAPSP $(objects)

clean:
	rm -f core *.exe *.o *~ FastAPSP



