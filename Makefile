.PHONY: clean

all: cmake

clean:
	rm -r ./build

cmake:
	mkdir -p build; cp -r femocs/in build; cd build; cmake ..; make
