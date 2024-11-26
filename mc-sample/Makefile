#.PHONY all clean

PROG = sampleMC
CC = g++-13
PATHFLAGS = -I/opt/homebrew/Cellar/fftw/3.3.10_1/include -I/opt/homebrew/Cellar/open-mpi/5.0.2/include
LIBFLAGS = -L/opt/homebrew/Cellar/fftw/3.3.10_1/lib -lfftw3 -L/opt/homebrew/Cellar/open-mpi/5.0.2/lib -lmpi
OBJS = main.o functions.o system.o utility.o
$(PROG) : $(OBJS)
	$(CC) -g3 -o $(PROG) $(OBJS) $(PATHFLAGS) $(LIBFLAGS) 

main.o : main.cpp
	$(CC) -g3 -c main.cpp $(PATHFLAGS) $(LIBFLAGS)

functions.o : functions.cpp functions.h
	$(CC) -g3 -c functions.cpp $(PATHFLAGS) $(LIBFLAGS)
	
system.o : system.cpp system.h
	$(CC) -g3 -c system.cpp $(PATHFLAGS) $(LIBFLAGS)

utility.o : utility.cpp utility.h
	$(CC) -g3 -c utility.cpp $(PATHFLAGS) $(LIBFLAGS)

clean:
	rm -f core $(PROG) $(OBJS)
