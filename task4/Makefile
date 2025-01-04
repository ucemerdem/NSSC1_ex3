CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp -O3 -std=c++17 -march=native -ffast-math

SRCS = main.cpp

all: solver_static solver_static1 solver_dynamic


solver_static: $(SRCS)
	$(CXX) $(CXXFLAGS) -D SCHEDULE_MODE="static" -D OUTPUT_MODE=0 -o solver_static $(SRCS)
	

solver_static1: $(SRCS)
	$(CXX) $(CXXFLAGS) -D SCHEDULE_MODE="static,1" -D OUTPUT_MODE=1 -o solver_static1 $(SRCS)
	

solver_dynamic: $(SRCS)
	$(CXX) $(CXXFLAGS) -D SCHEDULE_MODE="dynamic" -D OUTPUT_MODE=2 -o solver_dynamic $(SRCS)
	

clean:
	rm -f solver_static solver_static1 solver_dynamic
