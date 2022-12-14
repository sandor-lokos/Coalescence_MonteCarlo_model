INCLUDE_DIR = includes
SOURCE_DIR = sources
OBJECT_DIR = objects

CXX					 = g++
CFLAGS				= -O3 -Wall -fPIC -fno-inline -std=c++11 -lm -I$(INCLUDE_DIR)

_DEPS = model_funcs.h
DEPS = $(patsubst %,$(INCLUDE_DIR)/%,$(_DEPS))

_OBJ = numerical_integrator.o rapidity_momentum_distribution.o rapidity_momentum_distribution_GKL.o helper_functions.o
OBJ = $(patsubst %,$(OBJECT_DIR)/%,$(_OBJ))

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cc $(DEPS)
	$(shell mkdir -p $(OBJECT_DIR))
	$(CXX) -c -o $@ $< $(CFLAGS)

numerical_integrator.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)
	
clean:
	rm -rf $(OBJECT_DIR)/*.o *.exe
