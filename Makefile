include ../femocs/share/makefile.femocs

HDFLAGS=$(patsubst -I%, -I../femocs/%, $(FEMOCS_HEADPATH))
NEW_PATHS=$(patsubst -L%, -L../femocs/%, $(FEMOCS_LIBPATH))

.PHONY: clean

test: src/main.cpp
	mkdir -p build; cp -r ../femocs/in build; c++ $^ -std=c++14 ${HDFLAGS} ${NEW_PATHS} ${FEMOCS_LIB} -o build/$@

clean:
	rm -r ./build