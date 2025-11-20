CXX = g++ -fopenmp
CXXFLAGS = -std=c++17 -Wall -Wextra
LDFLAGS = -lSDL2
TARGET = raytracer
SRC = raytracer.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run