all: wire2d.cpp
	$(CXX) -std=c++11 -lSDL2 -O3 -o wire2d $<

clean:
	rm wire2d
