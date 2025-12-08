# MPI compiler and runner
MPICC   ?= mpic++
MPIEXEC ?= mpiexec
NP      ?= 1   

# Paths
SRC_DIR = src
OBJ_DIR = obj

# Files
TARGET  = main
SRC     = $(SRC_DIR)/main.cpp
OBJ     = $(OBJ_DIR)/main.o

# OpenGL / GLFW / GLEW / GLUT libraries
# (Assuming you've installed: libglew-dev libglfw3-dev libglm-dev freeglut3-dev libglu1-mesa-dev)
GL_LIBS   = -lglfw -lGLEW -lGL -lGLU -lglut
GL_CFLAGS =    # add -I... here if you installed to a non-standard location

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
	$(MPIEXEC) -n $(NP) ./$(TARGET) -i ./inputs/nb-10000.txt -s 1

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
