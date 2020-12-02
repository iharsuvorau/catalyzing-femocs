FEMOCS_HEADPATH=-I femocs/include -I femocs/lib -I femocs/dealii/include -I femocs/GETELEC/modules -std=c++1z 

FEMOCS_LIBPATH=-L femocs/lib -L femocs/GETELEC/lib -L femocs/dealii/lib

FEMOCS_LIB=-lfemocs -ltet -ldeal_II -lgetelec -lslatec -fopenmp -ltbb -llapack -lz -lm -lstdc++ -lgfortran 

VTKFLAGS=-I/usr/include/vtk-6.2

CMAKE_FLAGS=

.PHONY: clean

test: src/main.cpp
	mkdir -p build; cp -r femocs/in build; c++ $^ -std=c++14 ${FEMOCS_HEADPATH} ${FEMOCS_LIBPATH} ${FEMOCS_LIB} ${VTKFLAGS} -o build/$@

clean:
	rm -r ./build

cmake:
	mkdir -p build; cd build; cmake ${CMAKE_FLAGS} ..; make
