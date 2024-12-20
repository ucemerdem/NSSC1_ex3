CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math
#EIGEN_DIR = /usr/include/eigen3		#eigen

TARGET = cg

SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
.PHONY: all clean
