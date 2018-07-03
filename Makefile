#=================================================
# Hans A. Winther (2018) (hans.a.winther@gmail.com)
#=================================================

SHELL := /bin/bash

#====================================
# Executable
#====================================
EXEC = HODMatch

#====================================
# Set options, uncomment to use it
#====================================
#USE_OMP     = 1 # OpenMP for correlation function
#USE_MPI     = 1 # MPI for correlation function
#VELOCITY    = 1 # Include velocities in mock

#====================================
# Set compiler and add defines
#====================================
CC = g++-mp-6 -std=c++11 -O3
OPTIONS = -DPERIODIC -DVELOCITY
ifdef USE_MPI
CC = mpicc-openmpi-gcc6
OPTIONS += -DUSE_MPI
EXEC = HODMatch_MPI
endif
ifdef USE_OMP
OPTIONS += -DUSE_OMP -fopenmp
endif
C = -O3 $(OPTIONS)

#====================================
# GSL library
#====================================
INC = -I/opt/local/include/
LIB = -L/opt/local/lib/ -lgsl -lm

#====================================
# OBJECT FILES
#====================================
OBJS = src/simplex.o src/random_methods.o src/gsl_spline_wrapper.o src/cuter_library.o src/simplex.o src/hod.o src/rockstar_halo.o main.o

$(EXEC): $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIB)

%.o: %.cpp
	${CC} -c -o $@ $< $C $(INC)

%.o: %.c
	${CC} -c -o $@ $< $C $(INC)

clean:
	rm -rf $(TARGETS) *.o src/*.o $(EXEC)

