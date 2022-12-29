# Makefile for g++
CC = g++
SRC = main.cpp diagnose.cpp plasma.cpp particles.cpp poisson_solver_2d.cpp bfield.cpp 
OBJ := $(SRC:.cpp=.o)
CFLAGS = -fopenmp -O3 -lgsl -lgslcblas -lm -std=c++11
TARGET = pic
NTHREADS = 4

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(TARGET)
$(OBJ): %.o: %.cpp input.h
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	-rm $(OBJ) $(TARGET)
run:$(TARGET)
	-OMP_NUM_THREADS=$(NTHREADS) ./$(TARGET)

