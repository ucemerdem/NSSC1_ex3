CXX = g++
CXXFLAGS = -Wall -Wextra -fopenmp


TARGET = mcint
SRCS = main.cpp

all: $(TARGET)

$(TARGET): main.o
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.o

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c main.cpp

clean:
	rm -f $(TARGET) *.o