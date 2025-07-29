# Makefile
CXX=mpicxx
CXXFLAGS=-std=c++11 -O2 -Wall

SRCS= contact.cpp graph.cpp io_utils.cpp utils.cpp main_mpi.cpp
OBJS=$(SRCS:.cpp=.o)

all: detectar_cadenas_mpi

detectar_cadenas_mpi: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) detectar_cadenas_mpi

