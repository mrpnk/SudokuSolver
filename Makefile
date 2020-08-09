CXX = g++
CXXFLAGS = -std=c++2a


SRC := src
OBJ := obj

SOURCES := $(wildcard $(SRC)/*.cpp)
OBJECTS := $(patsubst $(SRC)/%.cpp, $(OBJ)/%.o, $(SOURCES))

all: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o sudokuSolver
	
	
$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) $< -I $(SRC) -c -o $@
	
