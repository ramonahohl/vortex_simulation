CC = g++
CFLAGS = -O3 -g -fopenmp -march=native -mtune=native -std=c++11 -mavx
#CFLAGS = -O0 -g -std=c++11 -fopenmp -Wall -Wpedantic
CFLAGS+= -D ORDER=12
CFLAGS += -MP -MD
ANIM_DIR = animationcvs

all: vortex

vortex: main.o kernels.o kernel2.o tree.o Vortexsimulation.o Celllist.o
	$(CC) $(CFLAGS) $^ -o $@ 
	
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

kernels.cpp: kernels.m4
	m4 -D ORDER=12 kernels.m4 > kernels.cpp

clean:
	rm -f *.o *.d kernels.cpp vortex #$(ANIM_DIR)/mesh.*
	
ifneq "$(MAKECMDGOALS)" "clean"
-include kernels.d main.d tree.d 
endif



