# MPI compiler and runner
MPICC   ?= mpic++
MPIEXEC ?= mpiexec
NP      ?= 4   

# Paths
SRC_DIR = src
OBJ_DIR = obj

# Files
TARGET  = nbody               
SRC     = $(SRC_DIR)/main.cpp
OBJ     = $(OBJ_DIR)/main.o



CXXFLAGS  = -O2 $(GL_CFLAGS)
LDFLAGS   = $(GL_LIBS)

# Default rule
all: $(TARGET)

# Build target
$(TARGET): $(OBJ)
	$(MPICC) -o $(TARGET) $(OBJ) $(LDFLAGS)

# Compile .cpp → .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(MPICC) $(CXXFLAGS) -c $< -o $@

# Run with mpiexec
run: $(TARGET)
	$(MPIEXEC) -n $(NP) ./$(TARGET) -i ./input/nb-10000.txt -s 1 -o ./outputs/nb-out.txt

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
