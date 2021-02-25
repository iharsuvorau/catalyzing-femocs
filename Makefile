.PHONY: clean

# default values, can be overwritten through command line, e.g. "make BUILD_DIR=foobar test"
BUILD_DIR = build
FEMOCS_DIR = femocs
PARAVIEW_DIR = paraview_build
DEAL_II_DIR = $(FEMOCS_DIR)/dealii/lib/cmake/deal.II
CMAKE_PREFIX_PATH = "$(FEMOCS_DIR)/lib;$(FEMOCS_DIR)/GETELEC/lib"

all: cmake

clean:
	@rm -r ./build

test:
	@echo $(BUILD_DIR)

cmake:
	@mkdir -p $(BUILD_DIR); \
	cmake -Ddeal.II_DIR=$(DEAL_II_DIR) \
		  -DParaView_DIR=$(PARAVIEW_DIR) \
		  -DFEMOCS_DIR=$(FEMOCS_DIR) \
		  -DCMAKE_PREFIX_PATH=$(CMAKE_PREFIX_PATH) \
		  -S . -B $(BUILD_DIR);
	cmake --build $(BUILD_DIR)

