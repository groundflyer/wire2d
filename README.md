# Wire2d

It implemets a simple realtime 2D wire dynamics program written as a testwork. It simulates one wire following the mouse pointer.

It uses Verlet integration and position-based approach to solve constraints.

## Requirements
* C++11 compatible compiler
* SDL2

## Building
Run `make`

## Usage
The are some command-line options
`-l <float>(default 200)	wire length in pixels.
-n <int>(default 20)		number of wire segments
-s <float>(default 1)		wire stiffness. Possible values are in range(0.1,1]
-g <float>(default 20)		gravity force`
