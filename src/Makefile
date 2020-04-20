CXX = g++
CXXFLAGS = -std=c++14 -g -I/usr/local/include -Wall
LDFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lm
PROG = test
OBJS = main.o
LINKER = $(CXX) -o

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

$(PROG): $(OBJS)
	$(LINKER) $(PROG) $(OBJS) $(LDFLAGS) 

main.o: main.cpp

clean:
	rm *.o
