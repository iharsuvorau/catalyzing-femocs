include femocs/share/makefile.femocs

HDFLAGS=$(patsubst -I%, -Ifemocs/%, $(FEMOCS_HEADPATH))
NEW_PATHS=$(patsubst -L%, -Lfemocs/%, $(FEMOCS_LIBPATH))
VTKFLAGS=-I/usr/include/vtk-6.2
CMAKE_FLAGS=

.PHONY: clean

test: src/main.cpp
	mkdir -p build; cp -r femocs/in build; c++ $^ -std=c++14 ${HDFLAGS} ${NEW_PATHS} ${FEMOCS_LIB} ${VTKFLAGS} -o build/$@

clean:
	rm -r ./build

cmake:
	mkdir -p build; cd build; cmake ${CMAKE_FLAGS} ..; make
